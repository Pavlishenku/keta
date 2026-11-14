#include "keta/soilzones.h"

#include "keta/bef.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {

// 8-point Gauss-Legendre on [0,1], exact through degree 15. See the header for
// why 8 points (a degree-6 integrand only needs 4).
constexpr int kG = 8;
constexpr double kGx[kG] = {0.019855071751231850, 0.101666761293186630, 0.237233795041835507,
                            0.408282678752175098, 0.591717321247824902, 0.762766204958164493,
                            0.898333238706813370, 0.980144928248768150};
constexpr double kGw[kG] = {0.050614268145188130, 0.111190517226687235, 0.156853322938943644,
                            0.181341891689180991, 0.181341891689180991, 0.156853322938943644,
                            0.111190517226687235, 0.050614268145188130};

// Shape functions of one bed component at abscissa x — IDENTICAL to the ones
// element.cpp integrates in closed form for the full-bar bed.
struct BedShape {
  int n = 2;
  int dof[4]{};
  double N[4]{};
};

BedShape bed_shape(int comp, double x, double L) {
  BedShape s;
  double xi = x / L;
  if (comp == kBedFx || comp == kBedMx) {
    s.n = 2;
    s.dof[0] = comp == kBedFx ? 0 : 3;
    s.dof[1] = comp == kBedFx ? 6 : 9;
    s.N[0] = 1 - xi;
    s.N[1] = xi;
    return s;
  }
  // Hermite cubics. l = +L in the x-y plane, -L in the x-z plane where the
  // rotation DOF is theta_y = -dw/dx (same sgn trick as element.cpp).
  bool y_plane = comp == kBedFy;
  double l = y_plane ? L : -L;
  s.n = 4;
  s.dof[0] = y_plane ? 1 : 2;
  s.dof[1] = y_plane ? 5 : 4;
  s.dof[2] = y_plane ? 7 : 8;
  s.dof[3] = y_plane ? 11 : 10;
  double x2 = xi * xi, x3 = x2 * xi;
  s.N[0] = 1 - 3 * x2 + 2 * x3;
  s.N[1] = l * (xi - 2 * x2 + x3);
  s.N[2] = 3 * x2 - 2 * x3;
  s.N[3] = l * (-x2 + x3);
  return s;
}

// u_comp(x) from the element-end local displacements.
double eval_u(int comp, const Vec12& u, double x, double L) {
  BedShape s = bed_shape(comp, x, L);
  double v = 0;
  for (int a = 0; a < s.n; ++a) v += s.N[a] * u(s.dof[a]);
  return v;
}

// INTEGRAL of u_comp over [x0, x1] (Gauss, exact for the cubic field).
double integrate_u(int comp, const Vec12& u, double x0, double x1, double L) {
  double len = x1 - x0;
  if (len <= 0) return 0;
  double acc = 0;
  for (int g = 0; g < kG; ++g) acc += kGw[g] * len * eval_u(comp, u, x0 + kGx[g] * len, L);
  return acc;
}

// Consistent nodal vector of a CONSTANT line intensity q of component `comp`
// over [x0, x1]: f += INTEGRAL N^T q dx. Same shapes as the bed matrix, so the
// capped reaction is work-conjugate to the bed it replaces.
void add_constant_line_load(int comp, double q, double x0, double x1, double L, Vec12& f) {
  double len = x1 - x0;
  if (len <= 0 || q == 0) return;
  for (int g = 0; g < kG; ++g) {
    double w = kGw[g] * len * q;
    BedShape s = bed_shape(comp, x0 + kGx[g] * len, L);
    for (int a = 0; a < s.n; ++a) f(s.dof[a]) += w * s.N[a];
  }
}

std::array<double, kBedComps> zone_moduli(const SoilZone& z) {
  std::array<double, kBedComps> k{};
  k[kBedFx] = z.kfx;
  k[kBedMx] = z.kmx;
  k[kBedFy] = z.kfy;
  k[kBedFz] = z.kfz;
  return k;
}

std::array<double, kBedComps> plastic_thresholds(const SoilPlastic& p) {
  std::array<double, kBedComps> f{};
  f[kBedFx] = p.ffx;
  f[kBedMx] = p.fmx;
  f[kBedFy] = p.ffy;
  f[kBedFz] = p.ffz;
  return f;
}

// Scatter a 12x12 nodal-frame element matrix into free-DOF triplets.
void scatter12(const Model& model, const DofMap& d, const Bar& b, const Mat12& k,
               std::vector<Triplet>& trips) {
  int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
  int gdof[12];
  for (int c = 0; c < 6; ++c) {
    gdof[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
    gdof[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
  }
  for (int a = 0; a < 12; ++a) {
    if (gdof[a] < 0) continue;
    for (int bb = 0; bb < 12; ++bb) {
      if (gdof[bb] < 0) continue;
      double v = k(a, bb);
      if (v != 0.0) trips.emplace_back(gdof[a], gdof[bb], v);
    }
  }
}

}  // namespace

// ---------------------------------------------------------------------------
// Zone resolution and bed matrices
// ---------------------------------------------------------------------------

std::vector<BedSegment> resolve_zones(const std::vector<SoilZone>& zones, double L) {
  require(!zones.empty(), "soil zones: empty zone list");
  require(finite(L) && L > 0, "soil zones: non-positive bar length");
  std::vector<BedSegment> out;
  out.reserve(zones.size());
  double prev = 0.0;
  int n = static_cast<int>(zones.size());
  for (int i = 0; i < n; ++i) {
    const SoilZone& z = zones[static_cast<std::size_t>(i)];
    for (double k : {z.kfx, z.kfy, z.kfz, z.kmx})
      require(finite(k) && k >= 0, "soil zones: negative or non-finite modulus in zone " +
                                       std::to_string(i + 1));
    // The LAST zone always runs to the bar end, whatever XL it carries
    // (ST1 manual p. 61 and p. 66).
    double end = (i == n - 1) ? L : zones[static_cast<std::size_t>(i)].xl_end;
    if (i < n - 1) {
      require(finite(end), "soil zones: non-finite XL in zone " + std::to_string(i + 1));
      require(end > prev + 1e-12, "soil zones: XL of zone " + std::to_string(i + 1) +
                                      " is not strictly increasing");
      require(end < L - 1e-12, "soil zones: XL of zone " + std::to_string(i + 1) +
                                   " reaches or exceeds the bar length");
    }
    BedSegment s;
    s.x0 = prev;
    s.x1 = end;
    s.k = zone_moduli(z);
    s.zone = i;
    out.push_back(s);
    prev = end;
  }
  return out;
}

Mat12 bed_matrix(const std::vector<BedSegment>& segs, double L) {
  require(finite(L) && L > 0, "bed matrix: non-positive bar length");
  Mat12 k = Mat12::Zero();
  for (const BedSegment& s : segs) {
    double len = s.length();
    if (len <= 0) continue;
    for (int comp = 0; comp < kBedComps; ++comp) {
      double kk = s.k[static_cast<std::size_t>(comp)];
      if (kk <= 0) continue;
      for (int g = 0; g < kG; ++g) {
        double w = kGw[g] * len * kk;
        BedShape sh = bed_shape(comp, s.x0 + kGx[g] * len, L);
        for (int a = 0; a < sh.n; ++a)
          for (int b = 0; b < sh.n; ++b)
            k(sh.dof[a], sh.dof[b]) += w * sh.N[a] * sh.N[b];
      }
    }
  }
  return k;
}

Mat12 zoned_bed_matrix(const std::vector<SoilZone>& zones, double L) {
  return bed_matrix(resolve_zones(zones, L), L);
}

Mat12 zoned_bed_matrix_nodal(const BarElement& el, const std::vector<SoilZone>& zones) {
  Mat12 k = zoned_bed_matrix(zones, el.L());
  return el.T().transpose() * k * el.T();
}

void validate_soil_zones(const Model& model, const SoilZoneMap& soil) {
  model.index();
  for (const auto& [id, bs] : soil) {
    const Bar& b = model.bar(id);  // throws on a dangling id
    std::string where = "soil zones on bar " + std::to_string(id) + ": ";
    require(b.kind == BarKind::Frame, where + "only frame bars carry a Winkler bed");
    require(!b.soil.any(), where + "the bar already carries a uniform Bar::soil bed "
                                   "(a zoned bed would double it)");
    require(!bs.zones.empty(), where + "empty zone list");
    // Length-independent checks; the "< L" bound is enforced by resolve_zones
    // against the true element length (which the eccentricity may change).
    double prev = 0;
    for (std::size_t i = 0; i + 1 < bs.zones.size(); ++i) {
      double xl = bs.zones[i].xl_end;
      require(finite(xl), where + "non-finite XL in zone " + std::to_string(i + 1));
      require(xl > prev + 1e-12,
              where + "XL of zone " + std::to_string(i + 1) + " is not strictly increasing");
      prev = xl;
    }
    for (const SoilZone& z : bs.zones)
      for (double k : {z.kfx, z.kfy, z.kfz, z.kmx})
        require(finite(k) && k >= 0, where + "negative or non-finite modulus");
    if (bs.plastic) {
      const SoilPlastic& p = *bs.plastic;
      for (double f : {p.ffx, p.ffy, p.ffz, p.fmx})
        require(finite(f) && f >= 0, where + "negative or non-finite plastic threshold");
    }
  }
}

// ---------------------------------------------------------------------------
// Pseudo-plastic driver
// ---------------------------------------------------------------------------

namespace {

// Per-bar working state of the active set.
struct BarBed {
  std::size_t e = 0;  // element / bar index
  Id id = 0;
  double L = 1;
  std::vector<BedSegment> subs;    // active-set sub-segments (zones when linear)
  std::vector<BedSegState> state;  // parallel to subs
  std::optional<SoilPlastic> plastic;
  bool nonlinear = false;
};

// Split each zone segment into `sub` equal parts.
std::vector<BedSegment> subdivide(const std::vector<BedSegment>& zones, int sub) {
  if (sub <= 1) return zones;
  std::vector<BedSegment> out;
  out.reserve(zones.size() * static_cast<std::size_t>(sub));
  for (const BedSegment& z : zones) {
    double len = z.length();
    for (int i = 0; i < sub; ++i) {
      BedSegment s = z;
      s.x0 = z.x0 + len * i / sub;
      s.x1 = i == sub - 1 ? z.x1 : z.x0 + len * (i + 1) / sub;
      out.push_back(s);
    }
  }
  return out;
}

// Moduli actually assembled this sweep: a non-elastic component contributes no
// stiffness.
std::vector<BedSegment> active_segments(const BarBed& bb) {
  std::vector<BedSegment> out = bb.subs;
  for (std::size_t i = 0; i < out.size(); ++i)
    for (int c = 0; c < kBedComps; ++c)
      if (bb.state[i].state[static_cast<std::size_t>(c)] != BedState::Elastic)
        out[i].k[static_cast<std::size_t>(c)] = 0.0;
  return out;
}

// Local force vector of every CAPPED sub-segment: a constant line intensity of
// magnitude ff opposing the motion. Uncondensed, exactly like the bed matrix.
Vec12 capped_load_local(const BarBed& bb) {
  Vec12 f = Vec12::Zero();
  for (std::size_t i = 0; i < bb.subs.size(); ++i) {
    for (int c = 0; c < kBedComps; ++c) {
      if (bb.state[i].state[static_cast<std::size_t>(c)] != BedState::Capped) continue;
      // r = reported reaction intensity (keta sign, follows the displacement);
      // the force ON THE BAR is -r.
      double r = bb.state[i].r_mean[static_cast<std::size_t>(c)];
      add_constant_line_load(c, -r, bb.subs[i].x0, bb.subs[i].x1, bb.L, f);
    }
  }
  return f;
}

// Local end displacements the bed sees: T * u_nodal, WITHOUT release
// re-injection — the bed is added after condensation, so it is work-conjugate
// to the raw transformed nodal displacements (mirrors element.cpp exactly).
Vec12 bed_local_disp(const Model& model, const BarElement& el, const DofMap& d,
                     const VecX& u_full) {
  const Bar& b = el.bar();
  int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
  Vec12 u_nodal;
  for (int c = 0; c < 6; ++c) {
    u_nodal(c) = u_full(d.full(ni, c));
    u_nodal(6 + c) = u_full(d.full(nj, c));
  }
  return el.T() * u_nodal;
}

VecX full_displacements(const DofMap& d, const CaseLoads& cl, const VecX& u_free) {
  VecX u_full = VecX::Zero(kDofPerNode * d.n_nodes);
  for (int f = 0; f < kDofPerNode * d.n_nodes; ++f) {
    int fi = d.free_ix[static_cast<std::size_t>(f)];
    if (fi >= 0) u_full(f) = u_free(fi);
  }
  for (const auto& [pf, pv] : cl.prescribed) u_full(pf) = pv;
  return u_full;
}

// Decide the state of every sub-segment from the current displacement field.
// Returns true when at least one state changed.
// `mean_u(c, i)`: the mean displacement of component c over sub-segment i.
template <class MeanU>
bool update_states(BarBed& bb, const MeanU& mean_u, double tol) {
  bool changed = false;
  std::array<double, kBedComps> ff{};
  bool bilateral = true, pos = false;
  if (bb.plastic) {
    ff = plastic_thresholds(*bb.plastic);
    bilateral = bb.plastic->bilateral();
    pos = bb.plastic->pos;
  }
  for (std::size_t i = 0; i < bb.subs.size(); ++i) {
    const BedSegment& s = bb.subs[i];
    BedSegState& st = bb.state[i];
    st.x0 = s.x0;
    st.x1 = s.x1;
    double len = s.length();
    for (int c = 0; c < kBedComps; ++c) {
      std::size_t ci = static_cast<std::size_t>(c);
      double k = s.k[ci];
      if (k <= 0 || len <= 0) {
        st.p_mean[ci] = st.r_mean[ci] = 0;
        st.state[ci] = BedState::Elastic;  // nothing to switch off
        continue;
      }
      // Mean reported reaction intensity over the sub-segment (keta sign).
      double p = k * mean_u(c, i);
      st.p_mean[ci] = p;
      // ST1 "pression" carries the sign of k.u_local: Exemple11 (a beam loaded
      // downwards on a PLASTIQUE NEG bed) keeps its bed where it moves down.
      double p_st1 = p;
      BedState want;
      double scale = std::max(std::fabs(p), 1.0);
      if (!bilateral && (pos ? p_st1 < -tol * scale : p_st1 > tol * scale)) {
        want = BedState::Inactive;
      } else if (ff[ci] > 0 && std::fabs(p) > ff[ci] * (1.0 + tol)) {
        want = BedState::Capped;
      } else {
        want = BedState::Elastic;
      }
      // Hysteresis: a capped/inactive segment only returns to elastic once it
      // is clearly back inside the admissible domain (kills 2-cycles).
      if (st.state[ci] == BedState::Capped && want == BedState::Elastic && ff[ci] > 0 &&
          std::fabs(p) > ff[ci] * (1.0 - tol))
        want = BedState::Capped;
      if (st.state[ci] == BedState::Inactive && want != BedState::Inactive && !bilateral &&
          std::fabs(p_st1) <= tol * scale)
        want = BedState::Inactive;
      if (want != st.state[ci]) {
        st.state[ci] = want;
        changed = true;
      }
      double r_old = st.r_mean[ci];
      st.r_mean[ci] = st.state[ci] == BedState::Inactive ? 0.0
                      : st.state[ci] == BedState::Capped
                          ? (p >= 0 ? ff[ci] : -ff[ci])
                          : p;
      // A capped segment whose pressure changed SIGN keeps its state but flips
      // the constant force it feeds back — that is a change of the active set.
      if (st.state[ci] == BedState::Capped && r_old * st.r_mean[ci] < 0) changed = true;
    }
  }
  return changed;
}

// Station recovery of the actual reaction intensity + its resultant, at the
// stations recover_case gave the bar (`xs`: the study sections when set, not
// necessarily the uniform grid).
void recover_soil(const BarElement& el, const BarBed& bb, const Vec12& u,
                  const std::vector<double>& xs, PseBarSoil& out) {
  int n = static_cast<int>(xs.size());
  double L = bb.L;
  out.bedded = true;
  out.soil_p.assign(static_cast<std::size_t>(n), Vec3::Zero());
  out.soil_m.assign(static_cast<std::size_t>(n), 0.0);
  if (bb.subs.empty()) return;
  auto seg_at = [&](double x) -> std::size_t {
    for (std::size_t i = 0; i < bb.subs.size(); ++i)
      if (x < bb.subs[i].x1 || i + 1 == bb.subs.size()) return i;
    return bb.subs.size() - 1;
  };
  // Per zone and component: the mean intensity of the actual reaction (the
  // capped force or k·u over each sub-segment) and the mean displacement.
  struct ZoneSecant {
    double r = 0, um = 0, len = 0;
    bool capped = false;
  };
  int nz = 0;
  for (const BedSegment& sg : bb.subs) nz = std::max(nz, sg.zone + 1);
  std::vector<std::array<ZoneSecant, kBedComps>> secant(static_cast<std::size_t>(nz));
  {
    static const double gx[5] = {-0.9061798459386640, -0.5384693101056831, 0.0,
                                 0.5384693101056831, 0.9061798459386640};
    static const double gw[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                                 0.4786286704993665, 0.2369268850561891};
    for (std::size_t i = 0; i < bb.subs.size(); ++i) {
      const BedSegment& sg = bb.subs[i];
      for (int c = 0; c < kBedComps; ++c) {
        const std::size_t ci = static_cast<std::size_t>(c);
        double iu = 0;  // integral of u over the sub-segment
        for (int g = 0; g < 5; ++g)
          iu += 0.5 * sg.length() * gw[g] *
                eval_u(c, u, 0.5 * (sg.x0 + sg.x1) + 0.5 * sg.length() * gx[g], L);
        ZoneSecant& zs = secant[static_cast<std::size_t>(sg.zone)][ci];
        const BedState bs = bb.state[i].state[ci];
        zs.r += bs == BedState::Inactive ? 0.0
                : bs == BedState::Capped ? bb.state[i].r_mean[ci] * sg.length()
                                         : sg.k[ci] * iu;
        zs.um += iu;
        zs.len += sg.length();
        zs.capped = zs.capped || bs == BedState::Capped;
      }
    }
    for (auto& z : secant)
      for (ZoneSecant& zs : z)
        if (zs.len > 0) {
          zs.r /= zs.len;
          zs.um /= zs.len;
        }
  }
  for (int j = 0; j < n; ++j) {
    double x = std::clamp(xs[static_cast<std::size_t>(j)], 0.0, L);
    std::size_t i = seg_at(x);
    const BedSegment& s = bb.subs[i];
    const BedSegState& st = bb.state[i];
    double val[kBedComps]{};
    for (int c = 0; c < kBedComps; ++c) {
      std::size_t ci = static_cast<std::size_t>(c);
      double k = s.k[ci];
      if (k <= 0) continue;
      const BedState state =
          secant[static_cast<std::size_t>(s.zone)][ci].capped ? BedState::Capped : st.state[ci];
      switch (state) {
        case BedState::Inactive: val[c] = 0; break;
        case BedState::Capped: {
          // ST1's pseudo-plastic bed caps the MEAN pressure of the bar and
          // spreads it as the displacement, a secant spring (p.121): a zone
          // with a capped part has its mean reaction spread as u(x)/mean(u)
          // all along (Exemple11: -287.80 / -312.42 on bar 35).
          const ZoneSecant& zs = secant[static_cast<std::size_t>(s.zone)][ci];
          const double ux = eval_u(c, u, x, L);
          val[c] = std::fabs(zs.um) > 1e-300 ? zs.r * ux / zs.um : st.r_mean[ci];
          break;
        }
        case BedState::Elastic: val[c] = k * eval_u(c, u, x, L); break;
      }
    }
    out.soil_p[static_cast<std::size_t>(j)] = Vec3(val[kBedFx], val[kBedFy], val[kBedFz]);
    out.soil_m[static_cast<std::size_t>(j)] = val[kBedMx];
  }
  // Resultant: exact integral of the piecewise field.
  double tot[kBedComps]{};
  for (std::size_t i = 0; i < bb.subs.size(); ++i) {
    const BedSegment& s = bb.subs[i];
    const BedSegState& st = bb.state[i];
    for (int c = 0; c < kBedComps; ++c) {
      std::size_t ci = static_cast<std::size_t>(c);
      double k = s.k[ci];
      if (k <= 0) continue;
      switch (st.state[ci]) {
        case BedState::Inactive: break;
        case BedState::Capped: tot[c] += st.r_mean[ci] * s.length(); break;
        case BedState::Elastic: tot[c] += k * integrate_u(c, u, s.x0, s.x1, L); break;
      }
    }
  }
  out.total_force_local = Vec3(tot[kBedFx], tot[kBedFy], tot[kBedFz]);
  out.total_torque_local = tot[kBedMx];
  out.total_force_global = el.R().transpose() * out.total_force_local;
  out.segments = bb.state;
}

// The member stations of a bedded bar WITH the soil reaction inside it: the
// actual reaction of the converged state (-k u where elastic, the capped
// plateau, nothing where inactive) is applied to the bar as line loads on
// fine linear pieces, its fixed-end forces taken off the end forces — the
// member then stands in equilibrium under its loads AND its soil, and a free
// end carries no shear (as in ST1: Ty = 0 at the free ends, Ty = 0
// all along under self weight).
void recover_bedded_stations(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                             const DofMap& d, const VecX& u_full, const BarBed& bb,
                             const Vec12& u, const std::vector<double>& xs, CaseResult& res) {
  const BarElement& el = asm_.elements[bb.e];
  const Bar& b = model.bars[bb.e];
  if (b.kind != BarKind::Frame || bb.subs.empty() || bb.e >= res.bars.size()) return;
  BarResult& br = res.bars[bb.e];
  const double L = bb.L;
  const Material& mat = model.material(b.material);
  BarSpanLoads sp = cl.span[bb.e];
  Vec12 f0_soil = Vec12::Zero();
  constexpr int kPieces = 8;
  for (std::size_t i = 0; i < bb.subs.size(); ++i) {
    const BedSegment& sg = bb.subs[i];
    const BedSegState& st = bb.state[i];
    auto reaction = [&](double x) {  // the soil's line force / torque on the bar
      double v[kBedComps]{};
      for (int c = 0; c < kBedComps; ++c) {
        const std::size_t ci = static_cast<std::size_t>(c);
        if (sg.k[ci] <= 0) continue;
        switch (st.state[ci]) {
          case BedState::Inactive: break;
          case BedState::Capped: v[c] = -st.r_mean[ci]; break;
          case BedState::Elastic: v[c] = -sg.k[ci] * eval_u(c, u, x, L); break;
        }
      }
      return std::make_pair(Vec3(v[kBedFx], v[kBedFy], v[kBedFz]), v[kBedMx]);
    };
    for (int k = 0; k < kPieces; ++k) {
      const double xa = sg.x0 + sg.length() * k / kPieces;
      const double xb = sg.x0 + sg.length() * (k + 1) / kPieces;
      if (xb - xa <= 0) continue;
      const auto [qa, ta] = reaction(xa);
      const auto [qb, tb] = reaction(xb);
      Load l = Load::linear(b.id, qa, qb, LoadFrame::Local, xa, xb);
      l.M = Vec3(ta, 0, 0);
      l.M2 = Vec3(tb, 0, 0);
      f0_soil += el.f0_local(l, mat);
      DistSeg ds;
      ds.a = xa;
      ds.b = xb;
      ds.qa = qa;
      ds.qb = qb;
      ds.ma = Vec3(ta, 0, 0);
      ds.mb = Vec3(tb, 0, 0);
      sp.segs.push_back(ds);
    }
  }
  const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
  Vec12 u_nodal;
  for (int c = 0; c < 6; ++c) {
    u_nodal(c) = u_full(d.full(ni, c));
    u_nodal(6 + c) = u_full(d.full(nj, c));
  }
  const Vec12 u_local = el.local_displacements(u_nodal, cl.f0[bb.e]);
  br.end_forces_local -= f0_soil;
  std::vector<double> x = xs;
  recover_bar_stations(model, el, sp, br.end_forces_local, u_local, x, br);
}

}  // namespace

PseResults solve_static_pse(const Model& model, const SoilZoneMap& soil, const StaticOptions& opt,
                            const PseOptions& pse) {
  require(pse.sub_per_zone >= 1, "solve_static_pse: sub_per_zone must be >= 1");
  require(pse.max_iterations >= 1, "solve_static_pse: max_iterations must be >= 1");
  require(finite(pse.tol) && pse.tol >= 0, "solve_static_pse: bad tolerance");
  model.validate();
  validate_soil_zones(model, soil);
  for (const Support& s : model.supports)
    for (int c = 0; c < 6; ++c)
      require(s.oneside[static_cast<std::size_t>(c)] == 0,
              "solve_static_pse: one-sided (uplift) supports are not supported by the "
              "pseudo-plastic soil driver");

  Assembled asm_ = assemble_elements(model);
  const DofMap& d = asm_.dofs;
  SpMat K_base = assemble_stiffness(model, asm_);

  // Per-bar working beds, in model bar order (deterministic).
  std::vector<BarBed> beds;
  bool nonlinear_any = false;
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    auto it = soil.find(model.bars[e].id);
    if (it == soil.end()) continue;
    BarBed bb;
    bb.e = e;
    bb.id = model.bars[e].id;
    bb.L = asm_.elements[e].L();
    bb.plastic = it->second.plastic;
    bb.nonlinear = bb.plastic && bb.plastic->nonlinear();
    std::vector<BedSegment> zones = resolve_zones(it->second.zones, bb.L);
    bb.subs = bb.nonlinear ? subdivide(zones, pse.sub_per_zone) : zones;
    bb.state.assign(bb.subs.size(), BedSegState{});
    for (std::size_t i = 0; i < bb.subs.size(); ++i) {
      bb.state[i].x0 = bb.subs[i].x0;
      bb.state[i].x1 = bb.subs[i].x1;
    }
    nonlinear_any = nonlinear_any || bb.nonlinear;
    beds.push_back(std::move(bb));
  }

  // Assemble K for the current active set.
  auto assemble_K = [&]() {
    std::vector<Triplet> trips;
    trips.reserve(beds.size() * 100 + 1);
    for (const BarBed& bb : beds) {
      const BarElement& el = asm_.elements[bb.e];
      Mat12 kb = bed_matrix(active_segments(bb), bb.L);
      scatter12(model, d, model.bars[bb.e], el.T().transpose() * kb * el.T(), trips);
    }
    SpMat Kb(d.n_free, d.n_free);
    Kb.setFromTriplets(trips.begin(), trips.end());
    return SpMat(K_base + Kb);
  };

  // Fold prescribed displacements through the BED part of K as well
  // (assembly.cpp's fold_prescribed only knows about ke(), which is bed-free
  // for a zoned bar).
  auto fold_bed_prescribed = [&](CaseLoads& cl) {
    if (cl.prescribed.empty()) return;
    for (const BarBed& bb : beds) {
      const BarElement& el = asm_.elements[bb.e];
      const Bar& b = model.bars[bb.e];
      Mat12 kb = el.T().transpose() * bed_matrix(active_segments(bb), bb.L) * el.T();
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      int gdof[12];
      for (int c = 0; c < 6; ++c) {
        gdof[c] = d.full(ni, c);
        gdof[6 + c] = d.full(nj, c);
      }
      for (const auto& [pf, pv] : cl.prescribed)
        for (int a = 0; a < 12; ++a) {
          if (gdof[a] != pf) continue;
          for (int r = 0; r < 12; ++r) {
            int fr = d.free_ix[static_cast<std::size_t>(gdof[r])];
            if (fr >= 0) cl.rhs(fr) -= kb(r, a) * pv;
          }
        }
    }
  };

  // Add the capped-segment line forces to the RHS (they are NOT recorded as
  // span loads: the bed never enters member-force recovery — keta invariant).
  auto add_capped_rhs = [&](VecX& rhs, std::vector<Vec12>* nodal_out) {
    for (std::size_t i = 0; i < beds.size(); ++i) {
      const BarBed& bb = beds[i];
      const BarElement& el = asm_.elements[bb.e];
      const Bar& b = model.bars[bb.e];
      Vec12 f = el.T().transpose() * capped_load_local(bb);
      if (nodal_out) (*nodal_out)[i] = f;
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      for (int c = 0; c < 6; ++c) {
        int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
        if (fi >= 0) rhs(fi) += f(c);
        int fj = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
        if (fj >= 0) rhs(fj) += f(6 + c);
      }
    }
  };

  // Reactions: recover_case builds R = K u - F from ke() only, which excludes
  // the zoned bed and the capped line forces. Correct both terms here.
  auto correct_reactions = [&](CaseResult& r, const VecX& u_full,
                               const std::vector<Vec12>& capped_nodal) {
    VecX corr = VecX::Zero(kDofPerNode * d.n_nodes);
    for (std::size_t i = 0; i < beds.size(); ++i) {
      const BarBed& bb = beds[i];
      const BarElement& el = asm_.elements[bb.e];
      const Bar& b = model.bars[bb.e];
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      Vec12 u_nodal;
      for (int c = 0; c < 6; ++c) {
        u_nodal(c) = u_full(d.full(ni, c));
        u_nodal(6 + c) = u_full(d.full(nj, c));
      }
      Mat12 kb = el.T().transpose() * bed_matrix(active_segments(bb), bb.L) * el.T();
      Vec12 f = kb * u_nodal - capped_nodal[i];
      for (int c = 0; c < 6; ++c) {
        corr(d.full(ni, c)) += f(c);
        corr(d.full(nj, c)) += f(6 + c);
      }
    }
    for (Reaction& rea : r.reactions) {
      int ni = model.node_index(rea.node);
      for (int c = 0; c < 6; ++c) rea.R(c) += corr(d.full(ni, c));
    }
  };

  PseResults out;
  FactorCache cache;

  // The exact path: every bedded bar takes the exact beam on elastic
  // foundation (bef.h) with the moduli of its current states (0 where capped
  // or lifted), a capped segment's constant reaction is a line load, and the
  // mean pressure of a sub-segment is the integral of the exact deflection
  // (8-point Gauss on pieces cut at the loads and no longer than 1/beta) —
  // the transfer solution the note RDM §6 itself uses, instead of the
  // Hermite cubic bed. Bars with releases or a CARA VAR kernel keep the
  // Hermite bed (install_bef rejects them).
  bool exact_ok = pse.exact;
  for (const BarBed& bb : beds) {
    const Bar& b = model.bars[bb.e];
    for (bool r : b.releases) exact_ok = exact_ok && !r;
    exact_ok = exact_ok && !model.tapered.count(b.id) && b.kind == BarKind::Frame;
  }
  auto solve_case_exact = [&](const LoadCase& lc) {
    PseCaseResult pcr;
    pcr.soil.resize(model.bars.size());
    for (std::size_t e = 0; e < model.bars.size(); ++e) pcr.soil[e].bar = model.bars[e].id;
    for (BarBed& bb : beds)
      for (BedSegState& st : bb.state) {
        st.state.fill(BedState::Elastic);
        st.r_mean.fill(0.0);
      }
    // Gauss abscissae of each sub-segment: (x, weight / length, sub index)
    struct Gp {
      double x, w;
      std::size_t i;
    };
    std::map<Id, std::vector<Gp>> gauss;
    std::map<Id, std::vector<double>> gx_of;
    static const double g8x[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                  -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                  0.7966664774136267,  0.9602898564975363};
    static const double g8w[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                  0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                  0.2223810344533745, 0.1012285362903763};
    for (const BarBed& bb : beds) {
      const Bar& b = model.bars[bb.e];
      const SectionProps& sp = model.section(b.section).props;
      const Material& mat = model.material(b.material);
      double beta = 0;
      for (const BedSegment& sg : bb.subs) {
        if (sg.k[kBedFy] > 0 && sp.Iz > 0) beta = std::max(beta, std::pow(sg.k[kBedFy] / (4 * mat.E * sp.Iz), 0.25));
        if (sg.k[kBedFz] > 0 && sp.Iy > 0) beta = std::max(beta, std::pow(sg.k[kBedFz] / (4 * mat.E * sp.Iy), 0.25));
        if (sg.k[kBedFx] > 0 && sp.A > 0) beta = std::max(beta, std::sqrt(sg.k[kBedFx] / (mat.E * sp.A)));
        if (sg.k[kBedMx] > 0 && sp.J > 0) beta = std::max(beta, std::sqrt(sg.k[kBedMx] / (mat.G() * sp.J)));
      }
      std::vector<double> cuts;
      for (const Load& l : lc.loads) {
        if (l.bar != b.id) continue;
        const double sc = l.relative ? bb.L : 1.0;
        cuts.push_back(l.a * sc);
        if (l.b >= 0) cuts.push_back(l.b * sc);
      }
      std::vector<Gp>& g = gauss[b.id];
      for (std::size_t i = 0; i < bb.subs.size(); ++i) {
        const BedSegment& sg = bb.subs[i];
        if (sg.length() <= 0) continue;
        std::vector<double> c{sg.x0};
        for (double x : cuts)
          if (x > sg.x0 + 1e-12 && x < sg.x1 - 1e-12) c.push_back(x);
        std::sort(c.begin(), c.end());
        c.push_back(sg.x1);
        for (std::size_t q = 0; q + 1 < c.size(); ++q) {
          const double len = c[q + 1] - c[q];
          const int m = std::max(1, static_cast<int>(std::ceil(len * beta)));
          for (int j = 0; j < m; ++j) {
            const double a = c[q] + len * j / m, h = len / m;
            for (int k = 0; k < 8; ++k)
              g.push_back({a + 0.5 * h * (1 + g8x[k]), 0.5 * h * g8w[k] / sg.length(), i});
          }
        }
      }
      std::sort(g.begin(), g.end(), [](const Gp& a, const Gp& b) { return a.x < b.x; });
      for (const Gp& p : g) gx_of[b.id].push_back(p.x);
    }
    auto build = [&]() {
      Model mm = model;
      mm.cases.clear();
      mm.combos.clear();
      LoadCase lc2 = lc;
      for (const BarBed& bb : beds) {
        std::vector<BefZone> zs;
        for (std::size_t i = 0; i < bb.subs.size(); ++i) {
          const BedSegment& sg = bb.subs[i];
          const BedSegState& st = bb.state[i];
          auto kk = [&](int c) {
            return st.state[static_cast<std::size_t>(c)] == BedState::Elastic ? sg.k[static_cast<std::size_t>(c)] : 0.0;
          };
          // contiguous sub-segments of equal moduli are ONE zone of the
          // exact block (a short piece only degrades its conditioning: the
          // bed is ~(beta l)^4 of the bending stiffness of a piece)
          const BefZone z{sg.x1, kk(kBedFx), kk(kBedFy), kk(kBedFz), kk(kBedMx)};
          if (!zs.empty() && zs.back().kx == z.kx && zs.back().ky == z.ky && zs.back().kz == z.kz &&
              zs.back().kt == z.kt)
            zs.back().x1 = z.x1;
          else
            zs.push_back(z);
          // a capped segment: its constant reaction, the soil's force on the bar
          Vec3 q = Vec3::Zero();
          double t = 0;
          for (int c = 0; c < kBedComps; ++c)
            if (st.state[static_cast<std::size_t>(c)] == BedState::Capped) {
              const double v = -st.r_mean[static_cast<std::size_t>(c)];
              if (c == kBedMx) t = v;
              else q(c == kBedFx ? 0 : c == kBedFy ? 1 : 2) = v;
            }
          if (q != Vec3::Zero() || t != 0) {
            Load l = Load::linear(bb.id, q, q, LoadFrame::Local, sg.x0, sg.x1);
            l.M = Vec3(t, 0, 0);
            l.M2 = Vec3(t, 0, 0);
            lc2.loads.push_back(l);
          }
        }
        install_bef(mm, bb.id, zs);
      }
      mm.cases.push_back(lc2);
      return mm;
    };
    auto disp = [](const Station& st, int c) {
      return c == kBedFx ? st.defl.x() : c == kBedFy ? st.defl.y() : c == kBedFz ? st.defl.z() : st.rot.x();
    };
    int it = 0;
    bool converged = false;
    for (; it < pse.max_iterations; ++it) {
      Model mm = build();
      StaticOptions so = opt;
      so.stations_abs = &gx_of;
      const StaticResults r = solve_static(mm, so);
      bool changed = false;
      for (BarBed& bb : beds) {
        const BarResult* br = nullptr;
        for (const BarResult& b : r.cases[0].bars)
          if (b.bar == bb.id) br = &b;
        if (!br) continue;
        const std::vector<Gp>& g = gauss[bb.id];
        std::vector<std::array<double, kBedComps>> mean(bb.subs.size(), std::array<double, kBedComps>{});
        for (std::size_t k = 0; k < g.size() && k < br->stations.size(); ++k)
          for (int c = 0; c < kBedComps; ++c)
            mean[g[k].i][static_cast<std::size_t>(c)] += g[k].w * disp(br->stations[k], c);
        changed = update_states(bb, [&](int c, std::size_t i) { return mean[i][static_cast<std::size_t>(c)]; },
                                pse.tol) ||
                  changed;
      }
      if (!changed) {
        converged = true;
        break;
      }
    }
    require(converged, "solve_static_pse: pseudo-plastic soil active set did not converge for "
                       "case '" + lc.name + "' in " + std::to_string(pse.max_iterations) +
                           " iterations");
    pcr.iterations = it + 1;
    const StaticResults rf = solve_static(build(), opt);
    pcr.result = rf.cases[0];
    pcr.result.non_superposable = nonlinear_any;
    for (const BarBed& bb : beds) {
      PseBarSoil& ps = pcr.soil[bb.e];
      ps.bedded = true;
      ps.segments = bb.state;
      const BarResult* br = nullptr;
      for (const BarResult& b : pcr.result.bars)
        if (b.bar == bb.id) br = &b;
      Vec3 tot = Vec3::Zero();
      double tt = 0;
      for (std::size_t i = 0; i < bb.subs.size(); ++i) {
        const double len = bb.subs[i].length();
        for (int c = 0; c < kBedComps; ++c) {
          const BedState sc = bb.state[i].state[static_cast<std::size_t>(c)];
          const double r = sc == BedState::Inactive ? 0.0
                           : sc == BedState::Capped ? bb.state[i].r_mean[static_cast<std::size_t>(c)]
                                                    : bb.state[i].p_mean[static_cast<std::size_t>(c)];
          if (c == kBedMx) tt += r * len;
          else tot(c == kBedFx ? 0 : c == kBedFy ? 1 : 2) += r * len;
          if (sc != BedState::Elastic) pcr.plastic = true;
        }
      }
      ps.total_force_local = tot;
      ps.total_torque_local = tt;
      ps.total_force_global = asm_.elements[bb.e].R().transpose() * tot;
      if (!br) continue;
      // ST1's edition of a zone with a capped part: its mean reaction spread
      // as the displacement, a secant spring (p.121, Exemple11 bar 35), as
      // recover_soil does on the Hermite path
      int nz = 0;
      for (const BedSegment& sg : bb.subs) nz = std::max(nz, sg.zone + 1);
      struct Sec {
        double r = 0, um = 0, len = 0;
        bool capped = false;
      };
      std::vector<std::array<Sec, kBedComps>> sec(static_cast<std::size_t>(nz));
      for (std::size_t i = 0; i < bb.subs.size(); ++i) {
        const BedSegment& sg = bb.subs[i];
        for (int c = 0; c < kBedComps; ++c) {
          const std::size_t ci = static_cast<std::size_t>(c);
          if (sg.k[ci] <= 0) continue;
          const BedState sc = bb.state[i].state[ci];
          Sec& z = sec[static_cast<std::size_t>(sg.zone)][ci];
          z.r += (sc == BedState::Inactive ? 0.0
                  : sc == BedState::Capped ? bb.state[i].r_mean[ci]
                                           : bb.state[i].p_mean[ci]) *
                 sg.length();
          z.um += bb.state[i].p_mean[ci] / sg.k[ci] * sg.length();
          z.len += sg.length();
          z.capped = z.capped || sc == BedState::Capped;
        }
      }
      for (const Station& st : br->stations) {
        std::size_t i = 0;
        while (i + 1 < bb.subs.size() && st.x > bb.subs[i].x1 + 1e-12) ++i;
        double v[kBedComps]{};
        for (int c = 0; c < kBedComps; ++c) {
          const std::size_t ci = static_cast<std::size_t>(c);
          const BedState sc = bb.state[i].state[ci];
          const Sec& z = sec[static_cast<std::size_t>(bb.subs[i].zone)][ci];
          if (bb.subs[i].k[ci] <= 0) continue;
          if (z.capped)
            v[c] = std::fabs(z.um) > 1e-300 ? (z.r / z.len) * disp(st, c) / (z.um / z.len)
                                             : bb.state[i].r_mean[ci];
          else
            v[c] = sc == BedState::Inactive ? 0.0 : bb.subs[i].k[ci] * disp(st, c);
        }
        ps.soil_p.push_back(Vec3(v[kBedFx], v[kBedFy], v[kBedFz]));
        ps.soil_m.push_back(v[kBedMx]);
      }
      // the bar's own pressure field is the actual reaction (the exact block
      // only knows the moduli of the current states: 0 where capped)
      for (BarResult& b : pcr.result.bars)
        if (b.bar == bb.id) {
          b.soil_p = ps.soil_p;
          b.soil_m = ps.soil_m;
        }
    }
    return pcr;
  };

  for (const LoadCase& lc : model.cases) {
    if (exact_ok) {
      out.cases.push_back(solve_case_exact(lc));
      continue;
    }
    PseCaseResult pcr;
    pcr.soil.resize(model.bars.size());
    for (std::size_t e = 0; e < model.bars.size(); ++e) pcr.soil[e].bar = model.bars[e].id;

    // Reset the active set for every case (pseudo-plastic soil is history-free
    // here: each case is solved from the virgin elastic state).
    for (BarBed& bb : beds)
      for (BedSegState& st : bb.state) st.state.fill(BedState::Elastic);

    // Member-load lowering is active-set independent: do it once per case.
    CaseLoads base = lower_case(model, asm_, lc);
    fold_prescribed(model, asm_, base);

    VecX u_free;
    VecX u_full;
    CaseLoads cl;
    std::vector<Vec12> capped_nodal(beds.size(), Vec12::Zero());
    int it = 0;
    bool converged = false;
    for (; it < pse.max_iterations; ++it) {
      SpMat K = assemble_K();
      cache.factorize(K);
      cl = base;
      fold_bed_prescribed(cl);
      for (Vec12& f : capped_nodal) f.setZero();
      add_capped_rhs(cl.rhs, &capped_nodal);
      u_free = cache.solve(cl.rhs);
      u_full = full_displacements(d, cl, u_free);

      bool changed = false;
      for (BarBed& bb : beds) {
        Vec12 u = bed_local_disp(model, asm_.elements[bb.e], d, u_full);
        changed = update_states(
                      bb,
                      [&](int c, std::size_t i) {
                        const BedSegment& sg = bb.subs[i];
                        return integrate_u(c, u, sg.x0, sg.x1, bb.L) / sg.length();
                      },
                      pse.tol) ||
                  changed;
      }
      if (!changed) {
        converged = true;
        break;
      }
    }
    require(converged, "solve_static_pse: pseudo-plastic soil active set did not converge for "
                       "case '" + lc.name + "' in " + std::to_string(pse.max_iterations) +
                           " iterations");
    pcr.iterations = it + 1;

    // The converged state may have been reached with a stale K if the LAST
    // sweep changed nothing; assemble_K/solve above already used it, so the
    // displacements are consistent with the reported states.
    pcr.result = recover_case(model, asm_, cl, u_free, opt, lc.id, lc.name);
    correct_reactions(pcr.result, u_full, capped_nodal);
    pcr.result.non_superposable = nonlinear_any;

    for (const BarBed& bb : beds) {
      Vec12 u = bed_local_disp(model, asm_.elements[bb.e], d, u_full);
      std::vector<double> xs;
      for (const BarResult& br : pcr.result.bars)
        if (br.bar == pcr.soil[bb.e].bar) {
          for (const Station& st : br.stations) xs.push_back(st.x);
          break;
        }
      if (xs.empty()) {
        const int n = std::max(1, opt.stations_per_bar);
        for (int j = 0; j < n; ++j) xs.push_back(n == 1 ? bb.L / 2 : bb.L * j / (n - 1));
      }
      recover_soil(asm_.elements[bb.e], bb, u, xs, pcr.soil[bb.e]);
      recover_bedded_stations(model, asm_, cl, d, u_full, bb, u, xs, pcr.result);
      for (const BedSegState& st : bb.state)
        for (int c = 0; c < kBedComps; ++c)
          if (st.state[static_cast<std::size_t>(c)] != BedState::Elastic) pcr.plastic = true;
    }
    out.cases.push_back(std::move(pcr));
  }

  if (!model.combos.empty()) {
    require(!nonlinear_any,
            "solve_static_pse: pseudo-plastic / one-sided soil results are non-superposable, "
            "combinations are not allowed");
    for (const Combination& cb : model.combos) {
      std::vector<const CaseResult*> srcs;
      std::vector<double> facs;
      for (const auto& [cid, f] : cb.terms) {
        for (std::size_t i = 0; i < model.cases.size(); ++i)
          if (model.cases[i].id == cid) srcs.push_back(&out.cases[i].result);
        facs.push_back(f);
      }
      require(srcs.size() == cb.terms.size(), "combination '" + cb.name + "': missing case");
      out.combos.push_back(combine_results(srcs, facs, cb.id, cb.name));
    }
  }
  return out;
}

}  // namespace keta

#include "keta/firechain.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <tuple>
#include <utility>

namespace keta {

// ---------------------------------------------------------------------------
// EN 1992-1-2 §3.2.2.1
// ---------------------------------------------------------------------------

double ConcreteFireLaw::sigma(double eps) const {
  if (eps >= 0) return tension ? e_origin() * eps : 0.0;
  double x = -eps;  // compression magnitude
  if (x <= ec1) {
    double u = x / ec1;
    return -3.0 * x * fc / (ec1 * (2.0 + u * u * u));
  }
  if (x < ecu1) return -fc * (ecu1 - x) / (ecu1 - ec1);  // linear descending branch
  return 0.0;                                            // crushed
}

double ConcreteFireLaw::tangent(double eps) const {
  if (eps > 0) return tension ? e_origin() : 0.0;
  double x = -eps;
  if (x <= ec1) {
    double u = x / ec1;
    double d = 2.0 + u * u * u;
    return 6.0 * fc * (1.0 - u * u * u) / (ec1 * d * d);
  }
  if (x < ecu1) return -fc / (ecu1 - ec1);
  return 0.0;
}

double ConcreteFireLaw::secant(double eps) const {
  if (std::fabs(eps) < 1e-14) return tension ? e_origin() : e_origin();
  return sigma(eps) / eps;
}

// ---------------------------------------------------------------------------
// Reinforcement
// ---------------------------------------------------------------------------

double RebarSteel::bar_area() const {
  if (area > 0) return area;
  require(diameter > 0, "firechain: ACIER needs a DIAMETRE or a SECTION");
  return kPi * diameter * diameter / 4.0;
}

std::vector<Rebar> resolve_rebars(const FireSection& section,
                                  const std::vector<RebarLayer>& layers) {
  std::vector<Rebar> out;
  for (const RebarLayer& l : layers) {
    require(l.count >= 1, "firechain: ARMATURE NOMBRE must be >= 1");
    double a = l.steel.bar_area();
    require(a > 0 && finite(a), "firechain: non-finite rebar area");
    if (l.placement == RebarLayer::Placement::Local) {
      require(finite(l.y) && finite(l.z), "firechain: non-finite ARMATURE LOCAL position");
      for (int i = 0; i < l.count; ++i) out.push_back({l.y, l.z, a, l.steel});
      continue;
    }
    require(section.kind != FireSection::Kind::Perimeter,
            "firechain: ARMATURE ENROBAGE SUP/INF is only defined for DALLE / RECTANGLE "
            "sections — use LOCAL Y Z on a PERIMETRE section (ST1 manual p223)");
    require(l.cover > 0 && l.cover < section.thickness,
            "firechain: ENROBAGE " + std::to_string(l.cover) +
                " m falls outside the section thickness " + std::to_string(section.thickness) +
                " m");
    // SUP / INF are the VY / WY faces (manual p220-223): the thickness runs
    // along the section's local Y, the width along Z.
    double y = (l.placement == RebarLayer::Placement::CoverSup)
                   ? section.thickness / 2 - l.cover
                   : -section.thickness / 2 + l.cover;
    // Spread the layer evenly across the width.
    for (int i = 0; i < l.count; ++i) {
      double z = section.width * ((i + 0.5) / l.count - 0.5);
      out.push_back({y, z, a, l.steel});
    }
  }
  return out;
}

// ---------------------------------------------------------------------------
// Fibre discretisation
// ---------------------------------------------------------------------------

double HotFibers::sigma(const Fibre& f, double eps_e) const {
  return f.is_concrete ? concrete[static_cast<std::size_t>(f.law)].sigma(eps_e)
                       : steel[static_cast<std::size_t>(f.law)].sigma(eps_e);
}

double HotFibers::tangent(const Fibre& f, double eps_e) const {
  return f.is_concrete ? concrete[static_cast<std::size_t>(f.law)].tangent(eps_e)
                       : steel[static_cast<std::size_t>(f.law)].tangent(eps_e);
}

double HotFibers::secant(const Fibre& f, double eps_e) const {
  if (f.is_concrete) return concrete[static_cast<std::size_t>(f.law)].secant(eps_e);
  const UniaxialLaw& l = steel[static_cast<std::size_t>(f.law)];
  if (std::fabs(eps_e) < 1e-14) return l.E;
  return l.sigma(eps_e) / eps_e;
}

double HotFibers::max_rebar_temp() const {
  double t = 0;
  for (const Fibre& f : fibres)
    if (!f.is_concrete) t = std::max(t, f.theta);
  return t;
}

HotFibers build_hot_fibers(const FireBarSpec& spec, const HeatResult& heat, int snap) {
  const std::vector<double>& T = heat.field(snap);
  const HeatGrid& g = heat.grid;
  // heat2d meshes a DALLE / RECTANGLE with its thickness along the grid z
  // (SUP = +z) and its width along the grid y; the section frame has the
  // thickness along local Y (SUP = VY, manual p220-223) and the width along
  // Z (VZ = the grid's +y face), so the grid axes are relabelled: Y = z_grid,
  // Z = y_grid. A PERIMETRE is meshed in its own (Y, Z) and needs nothing.
  const bool turned = spec.section.kind != FireSection::Kind::Perimeter;
  auto to_sec = [&](double gy, double gz) {
    return turned ? std::pair<double, double>{gz, gy} : std::pair<double, double>{gy, gz};
  };
  const auto& to_grid = to_sec;  // the relabelling is its own inverse
  HotFibers out;
  out.fibres.reserve(g.cells.size() + 8);

  // One concrete law per 1 C temperature bin (the EN tables are interpolated
  // linearly on 100 C steps, so 1 C bins are far below the data resolution and
  // keep the law table small).
  std::map<int, int> cbin;
  const FireConcrete& mc = spec.material;
  for (std::size_t k = 0; k < g.cells.size(); ++k) {
    double th = T[k];
    int key = static_cast<int>(std::lround(th));
    auto it = cbin.find(key);
    int law;
    if (it == cbin.end()) {
      ConcreteFireLaw cl;
      cl.fc = mc.fck * concrete_kc(th, mc.thermal.aggregate);
      cl.ec1 = concrete_eps_c1(th);
      cl.ecu1 = concrete_eps_cu1(th);
      cl.tension = mc.tension;
      // fc -> 0 at 1200 C: keep the law well-formed (zero strength, zero moduli).
      cl.fc = std::max(cl.fc, 0.0);
      law = static_cast<int>(out.concrete.size());
      out.concrete.push_back(cl);
      cbin.emplace(key, law);
    } else {
      law = it->second;
    }
    HotFibers::Fibre f;
    std::tie(f.y, f.z) = to_sec(g.cells[k].y, g.cells[k].z);
    f.area = g.cells[k].area;
    f.theta = th;
    f.eps_th = concrete_thermal_strain(th, mc.thermal.aggregate);
    f.law = law;
    f.is_concrete = true;
    out.fibres.push_back(f);
  }

  for (const Rebar& r : resolve_rebars(spec.section, spec.rebars)) {
    // Manual p223: the bar diameter selects the concrete cell whose temperature
    // the rebar takes — use the hottest cell it overlaps.
    double th = -1e300;
    double rad = 0.5 * std::max(r.steel.diameter, 0.0);
    const auto [ry, rz] = to_grid(r.y, r.z);
    // A cell counts when the bar's disc enters it: the nearest point of the
    // cell to the bar centre is closer than phi/2. (Sampling points spaced by
    // the cell pitch around the centre took two cells on one side and one on
    // the other whenever the centre fell on a cell boundary — 0.15 m on a
    // 0.01 m MAILLE — and heated the bars of a symmetric section unequally.)
    for (std::size_t c = 0; c < g.cells.size(); ++c) {
      double dy = std::max(std::fabs(g.cells[c].y - ry) - 0.5 * g.hy, 0.0);
      double dz = std::max(std::fabs(g.cells[c].z - rz) - 0.5 * g.hz, 0.0);
      if (dy * dy + dz * dz < rad * rad * (1 - 1e-9)) th = std::max(th, T[c]);
    }
    if (th < -1e299) {
      int c = g.locate(ry, rz);
      require(c >= 0, "firechain: rebar at (" + std::to_string(r.y) + ", " +
                          std::to_string(r.z) + ") falls outside the meshed section");
      th = T[static_cast<std::size_t>(c)];
    }
    HotFibers::Fibre f;
    f.y = r.y;
    f.z = r.z;
    f.area = r.area;
    f.theta = th;
    f.eps_th = rebar_thermal_strain(th);
    f.law = static_cast<int>(out.steel.size());
    f.is_concrete = false;
    out.steel.push_back(hot_rebar_law(r.steel.E, r.steel.fyk, th, r.steel.cls, r.steel.forming));
    out.fibres.push_back(f);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Hot section state
// ---------------------------------------------------------------------------

HotSectionState hot_section_state(const HotFibers& fib, const Vec3& e, double yref, double zref,
                                  bool secant, bool centroid_var, bool ignore_thermal) {
  HotSectionState s;
  s.secant = secant;
  double S0 = 0, Sy = 0, Sz = 0, Syy = 0, Szz = 0;
  double Se = 0, Sey = 0, Sez = 0;
  double area = 0, t_area = 0, t_max = 0;
  for (const HotFibers::Fibre& f : fib.fibres) {
    double eps = e(0) + e(1) * (f.z - zref) - e(2) * (f.y - yref);
    double eth = ignore_thermal ? 0.0 : f.eps_th;
    double eps_e = eps - eth;
    double E = secant ? fib.secant(f, eps_e) : fib.tangent(f, eps_e);
    // A negative modulus (softening branch of the concrete law) would destroy the
    // positive definiteness of the frame stiffness: clamp it to zero and let the
    // load shed to the fibres that still carry — the documented secant-approach
    // limitation.
    if (!(E > 0) || !finite(E)) E = 0.0;
    double EA = E * f.area;
    S0 += EA;
    Sy += EA * f.y;
    Sz += EA * f.z;
    Syy += EA * f.y * f.y;
    Szz += EA * f.z * f.z;
    Se += EA * eth;
    Sey += EA * eth * f.y;
    Sez += EA * eth * f.z;
    if (f.is_concrete) {
      area += f.area;
      t_area += f.area * f.theta;
      t_max = std::max(t_max, f.theta);
    }
  }
  s.theta_max = t_max;
  s.theta_mean = area > 0 ? t_area / area : 0.0;
  s.theta_rebar_max = fib.max_rebar_temp();
  require(S0 > 0, "firechain: the heated section has lost all axial rigidity (EA = 0) — "
                  "total section failure");
  s.EA = S0;
  s.yc = centroid_var ? Sy / S0 : 0.0;
  s.zc = centroid_var ? Sz / S0 : 0.0;
  s.EIy = Szz - S0 * s.zc * s.zc;
  s.EIz = Syy - S0 * s.yc * s.yc;
  require(s.EIy > 0 && s.EIz > 0,
          "firechain: the heated section has lost all bending rigidity (EI = 0) — "
          "total section failure");
  s.N_th = Se;
  s.M_th_y = Sez - Se * s.zc;
  s.M_th_z = -(Sey - Se * s.yc);
  s.eps_th = s.N_th / s.EA;
  s.kappa_th_y = s.M_th_y / s.EIy;
  s.kappa_th_z = s.M_th_z / s.EIz;
  return s;
}

namespace {

// Fibre resultants (N, My, Mz) of a strain plane about (yref, zref).
Vec3 fibre_resultants(const HotFibers& fib, const Vec3& e, double yref, double zref,
                      bool ignore_thermal) {
  Vec3 r = Vec3::Zero();
  for (const HotFibers::Fibre& f : fib.fibres) {
    double dy = f.y - yref, dz = f.z - zref;
    double eps = e(0) + e(1) * dz - e(2) * dy;
    double eps_e = eps - (ignore_thermal ? 0.0 : f.eps_th);
    double s = fib.sigma(f, eps_e) * f.area;
    r(0) += s;
    r(1) += s * dz;
    r(2) -= s * dy;
  }
  return r;
}

Mat3 fibre_tangent(const HotFibers& fib, const Vec3& e, double yref, double zref,
                   bool ignore_thermal) {
  Mat3 k = Mat3::Zero();
  for (const HotFibers::Fibre& f : fib.fibres) {
    double dy = f.y - yref, dz = f.z - zref;
    double eps = e(0) + e(1) * dz - e(2) * dy;
    double eps_e = eps - (ignore_thermal ? 0.0 : f.eps_th);
    double Et = fib.tangent(f, eps_e);
    if (!(Et > 0) || !finite(Et)) Et = 0.0;
    double EA = Et * f.area;
    k(0, 0) += EA;
    k(0, 1) += EA * dz;
    k(0, 2) -= EA * dy;
    k(1, 1) += EA * dz * dz;
    k(1, 2) -= EA * dz * dy;
    k(2, 2) += EA * dy * dy;
  }
  k(1, 0) = k(0, 1);
  k(2, 0) = k(0, 2);
  k(2, 1) = k(1, 2);
  return k;
}

}  // namespace

Vec3 section_newton(const HotFibers& fib, const Vec3& target, const Vec3& start, double yref,
                    double zref, bool ignore_thermal, int max_iter, double tol) {
  Vec3 e = start;
  double scale = std::max({std::fabs(target(0)), std::fabs(target(1)), std::fabs(target(2)), 1.0});
  for (int it = 0; it < max_iter; ++it) {
    Vec3 r = fibre_resultants(fib, e, yref, zref, ignore_thermal) - target;
    if (r.lpNorm<Eigen::Infinity>() <= tol * scale) break;
    Mat3 k = fibre_tangent(fib, e, yref, zref, ignore_thermal);
    Eigen::FullPivLU<Mat3> lu(k);
    if (!lu.isInvertible()) break;  // fully cracked / crushed: keep the last plane
    Vec3 de = lu.solve(-r);
    if (!de.allFinite()) break;
    // Damped step: the fibre laws have plateaus, an undamped Newton overshoots.
    double a = 1.0;
    double r0 = r.norm();
    for (int ls = 0; ls < 12; ++ls) {
      Vec3 trial = e + a * de;
      Vec3 rt = fibre_resultants(fib, trial, yref, zref, ignore_thermal) - target;
      if (rt.allFinite() && rt.norm() < r0) {
        e = trial;
        break;
      }
      a *= 0.5;
      if (ls == 11) return e;
    }
  }
  return e;
}

// ---------------------------------------------------------------------------
// INCENDIE driver
// ---------------------------------------------------------------------------

namespace {

// Trapezoidal mean of a station field over the bar (manual p213: the fire
// analysis uses the MEAN deformation of the whole bar, no sub-discretisation).
double mean_station(const std::vector<Station>& st, double Station::*field) {
  require(!st.empty(), "firechain: bar has no recovery station");
  if (st.size() == 1) return st[0].*field;
  double num = 0, den = 0;
  for (std::size_t i = 1; i < st.size(); ++i) {
    double dx = st[i].x - st[i - 1].x;
    num += 0.5 * (st[i].*field + st[i - 1].*field) * dx;
    den += dx;
  }
  return den > 0 ? num / den : st[0].*field;
}

}  // namespace

FireOutcome solve_fire(const Model& model, const FireProblem& p, const FireAnalysisOptions& opt) {
  require(p.option != ModelOption::Grill,
          "firechain: fire analysis is not available for ST1 OPTION GRILL — PLANE and "
          "SPATIALE only (ST1 manual p213)");
  require(!p.bars.empty(), "firechain: INCENDIE needs at least one fire-exposed bar");
  require(opt.max_iter >= 1, "firechain: ANALYSE FEU MAX_ITER must be >= 1");
  require(opt.period >= 1, "firechain: ANALYSE FEU PERIODE must be >= 1");
  require(opt.tol > 0 && finite(opt.tol), "firechain: convergence tolerance must be > 0");
  model.validate();

  const std::size_t nf = p.bars.size();
  // the mean axial force of each bar at the previous iteration (N.e coupling)
  std::vector<double> n_prev(nf, 0.0);
  FireOutcome out;
  out.fields.reserve(nf);
  out.states.resize(nf);

  HeatOptions ho;
  ho.duration = p.duration;
  ho.t_init = p.t_init;
  ho.dt = p.dt;

  std::vector<int> bar_slot(nf), sec_slot(nf), mat_slot(nf);
  std::vector<HotFibers> fibres(nf);
  std::vector<SectionProps> cold_props(nf);
  std::vector<double> cold_ea(nf), cold_rho_a(nf), e_ref(nf);

  Model hot = model;
  Id next_sec = 0, next_mat = 0;
  for (const Section& s : hot.sections) next_sec = std::max(next_sec, s.id);
  for (const Material& m : hot.materials) next_mat = std::max(next_mat, m.id);

  for (std::size_t k = 0; k < nf; ++k) {
    const FireBarSpec& spec = p.bars[k];
    const Bar& b = model.bar(spec.bar);
    require(b.kind == BarKind::Frame,
            "firechain: bar " + std::to_string(spec.bar) +
                " is not a frame bar — fire analysis models bending of beam elements");
    require(!b.soil.any(), "firechain: bar " + std::to_string(spec.bar) +
                               " has a Winkler bed; the hot-section substitution would "
                               "silently rescale it — not supported");
    // Manual p218 defines DALLE for PLANE and RECTANGLE for SPATIALE; the ST1
    // layer lowers a widthful DALLE onto its equivalent RECTANGLE (lateral
    // faces adiabatic), so both kinds are accepted in PLANE.
    if (p.option == ModelOption::Plane)
      require(spec.section.kind == FireSection::Kind::Slab ||
                  spec.section.kind == FireSection::Kind::Rectangle,
              "firechain: ST1 OPTION PLANE defines the DALLE fire section (or its "
              "RECTANGLE equivalent) — manual p218");

    out.fields.push_back(solve_heat2d(spec.section, spec.material.thermal, p.exposures, ho));
    fibres[k] = build_hot_fibers(spec, out.fields.back());

    int bi = model.bar_index(spec.bar);
    bar_slot[k] = bi;
    const Material& m0 = model.material(b.material);
    const SectionProps& sp = model.section(b.section).props;
    require(sp.A > 0 && sp.Iy > 0 && sp.Iz > 0,
            "firechain: bar " + std::to_string(spec.bar) + " has a degenerate cold section");
    cold_props[k] = sp;
    e_ref[k] = m0.E;
    require(e_ref[k] > 0, "firechain: fire bar material needs E > 0");
    cold_ea[k] = m0.E * sp.A;
    cold_rho_a[k] = m0.rho * sp.A;

    Section hs;
    hs.id = ++next_sec;
    hs.name = "hot:" + std::to_string(spec.bar);
    hs.shape = ShapeGeneric{};
    hs.props = sp;
    sec_slot[k] = static_cast<int>(hot.sections.size());
    hot.sections.push_back(hs);

    Material hm = m0;
    hm.id = ++next_mat;
    hm.name = m0.name + " (hot)";
    mat_slot[k] = static_cast<int>(hot.materials.size());
    hot.materials.push_back(hm);

    hot.bars[static_cast<std::size_t>(bi)].section = hs.id;
    hot.bars[static_cast<std::size_t>(bi)].material = hm.id;
  }

  const LoadCase& preload = model.load_case(p.preload_case);
  hot.validate();

  StaticOptions sopt;
  sopt.stations_per_bar = opt.stations_per_bar;

  std::vector<Vec3> plane(nf, Vec3::Zero());
  std::vector<double> yref(nf, 0.0), zref(nf, 0.0);
  VecX u_prev;
  bool converged = false;
  // Stall guard: the plain secant fixed point on a cracking section can lock
  // into a 2-cycle (exemple9 at the true ST1 self-weight: the displacements
  // alternate with a relative change frozen at 2.3e-3 from iteration 15 on).
  // When the residual stops falling — not below 0.99 x its value two sweeps
  // earlier — the strain-plane update is under-relaxed, halving the factor at
  // every further stall (floor 1/16). A run that converges without stalling
  // never relaxes: its numbers are unchanged.
  // When the stall outlives the floor (25 further stalls at 1/16), the
  // section-level tangent refinement (PERIODE) is what keeps the iteration
  // cycling: a support section near its plastic moment, where the refined
  // plane and the secant pseudo-load disagree (e.g. a protected slab at 60
  // min). The iteration
  // then continues on the plain secant fixed point, relaxation reset.
  std::vector<double> res_hist;
  double omega = 1.0;
  bool refine = opt.period > 1;
  int floor_stalls = 0;

  for (int iter = 0; iter < opt.max_iter; ++iter) {
    out.iterations = iter + 1;
    // --- hot secant rigidities from the current strain planes ---
    for (std::size_t k = 0; k < nf; ++k) {
      HotSectionState st = hot_section_state(fibres[k], plane[k], yref[k], zref[k],
                                             /*secant=*/true, opt.centroid_var, p.local);
      st.bar = p.bars[k].bar;
      out.states[k] = st;

      SectionProps sp = cold_props[k];
      double ratio = st.EA / cold_ea[k];
      sp.A = st.EA / e_ref[k];
      sp.Iy = st.EIy / e_ref[k];
      sp.Iz = st.EIz / e_ref[k];
      // Shear and torsion rigidities are not integrated over the fibre field
      // (that needs the warping/Mindlin problems of annex B5 p384-388, which keta
      // does not port): they are scaled by the axial degradation ratio.
      sp.Avy = cold_props[k].Avy * ratio;
      sp.Avz = cold_props[k].Avz * ratio;
      sp.J = cold_props[k].J * ratio;
      hot.sections[static_cast<std::size_t>(sec_slot[k])].props = sp;
      // Keep self-weight identical to the cold model despite the area substitution.
      hot.materials[static_cast<std::size_t>(mat_slot[k])].rho = cold_rho_a[k] / sp.A;
      if (opt.fold_centroid_ecc) {
        Bar& b = hot.bars[static_cast<std::size_t>(bar_slot[k])];
        Mat3 R = bar_rotation(model.node(b.ni).x, model.node(b.nj).x, b.roll);
        Vec3 d = R.transpose() * Vec3(0.0, st.yc, st.zc);  // local (0, yc, zc) -> global
        b.ecc[0] = model.bars[static_cast<std::size_t>(bar_slot[k])].ecc[0] + d;
        b.ecc[1] = model.bars[static_cast<std::size_t>(bar_slot[k])].ecc[1] + d;
      }
    }

    // --- fire load case: preload + thermal equivalent strains ---
    LoadCase lc;
    lc.id = p.preload_case;
    lc.name = "INCENDIE";
    lc.loads = preload.loads;
    for (std::size_t k = 0; k < nf; ++k) {
      const HotSectionState& st = out.states[k];
      double ky = 0, kz = 0;
      Load l;
      l.kind = Load::Kind::Prestrain;
      l.bar = p.bars[k].bar;
      if (!p.local) {
        l.strain = st.eps_th;
        ky = st.kappa_th_y;
        kz = st.kappa_th_z;
      }
      if (!opt.fold_centroid_ecc) {
        // Hot eccentricity (as ST1): the bar stays on the node line
        // with its hot stiffnesses about the hot centroid (yc, zc); the axial
        // force of the node line acts at (-yc, -zc) from it, a moment
        // (-zc N, yc N) handed to the element as an imposed curvature, with N
        // from the previous iteration (fire column: reaction 0.0392 / ST1
        // 0.0352, 0.0538 without).
        ky += -st.zc * n_prev[k] / st.EIy;
        kz += st.yc * n_prev[k] / st.EIz;
      }
      if (p.local && ky == 0 && kz == 0) continue;
      // The section frame strains eps = e0 + ky z - kz y; a keta Prestrain
      // curvature has the opposite sense on both axes (kappa_z > 0 lengthens
      // the +y fibre, kappa_y > 0 the -z fibre, cf. BarThermal).
      l.kappa_y = -ky;
      l.kappa_z = -kz;
      lc.loads.push_back(l);
    }

    StaticContext ctx = StaticContext::build(hot);
    CaseResult r = ctx.solve_case(lc, sopt);
    out.result = r;

    VecX u(static_cast<int>(r.node_disp.size()) * 6);
    for (std::size_t i = 0; i < r.node_disp.size(); ++i)
      u.segment<6>(static_cast<int>(i) * 6) = r.node_disp[i];
    require(u.allFinite(), "firechain: non-finite displacement field");
    if (u_prev.size() == u.size()) {
      double num = (u - u_prev).lpNorm<Eigen::Infinity>();
      double den = std::max(u.lpNorm<Eigen::Infinity>(), 1e-12);
      out.residual = num / den;
      // under an under-relaxed strain-plane update the displacement change
      // shrinks by omega: the convergence test takes the full update
      if (out.residual / omega < opt.tol) {
        converged = true;
        break;
      }
      res_hist.push_back(out.residual);
      const std::size_t n = res_hist.size();
      if (n >= 3 && res_hist[n - 1] >= 0.99 * res_hist[n - 3]) {
        if (omega == 1.0 / 16 && refine && ++floor_stalls >= 25) {
          refine = false;
          omega = 1.0;
          res_hist.clear();
        } else {
          omega = std::max(omega / 2, 1.0 / 16);
        }
      }
    }
    u_prev = u;

    // --- new strain plane from the mean member forces of the bar ---
    for (std::size_t k = 0; k < nf; ++k) {
      const BarResult& br = r.bars[static_cast<std::size_t>(bar_slot[k])];
      double N = mean_station(br.stations, &Station::N);
      double My = mean_station(br.stations, &Station::My);
      double Mz = mean_station(br.stations, &Station::Mz);
      n_prev[k] = N;
      const HotSectionState& st = out.states[k];
      yref[k] = st.yc;
      zref[k] = st.zc;
      if (!opt.fold_centroid_ecc) {
        // the node-line member moments transferred to the hot centroid:
        // Mz_c = Mz_0 + yc N, My_c = My_0 - zc N
        Mz += st.yc * N;
        My -= st.zc * N;
      }
      Vec3 e((N + st.N_th) / st.EA, (My + st.M_th_y) / st.EIy, (Mz + st.M_th_z) / st.EIz);
      // Tangent refinement of the strain plane at frozen member forces
      // (ST1 PERIODE, manual p226): period - 1 section-level Newton passes.
      if (refine)
        e = section_newton(fibres[k], Vec3(N, My, Mz), e, yref[k], zref[k], p.local,
                           opt.period - 1);
      require(e.allFinite(), "firechain: non-finite section strain plane");
      plane[k] = omega == 1.0 ? e : Vec3(plane[k] + omega * (e - plane[k]));
    }
  }

  require(converged,
          "firechain: the fire analysis did not converge in " + std::to_string(out.iterations) +
              " iterations (last relative displacement change " + std::to_string(out.residual) +
              ") — a non-converging fire analysis IS the collapse signal, the structure has "
              "no equilibrium under this load at this fire duration");

  if (opt.keep_hot) {  // ST1 CARA GARDER (manual p227)
    out.hot_model = hot;
    out.hot_kept = true;
  }
  return out;
}

}  // namespace keta

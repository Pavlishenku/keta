#include "keta/phasing.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <set>

#include "keta/bef.h"
#include "keta/cable.h"
#include "keta/composite.h"
#include "keta/creep.h"
#include "keta/element.h"
#include "keta/hereditary.h"
#include "keta/matlaws.h"
#include "keta/shell_tendon.h"

namespace keta {

namespace {
// Bar loads with the same support — bar, kind, frame, abscissae — summed into
// one: the equivalent loads of many tendons share their spans (the same
// uniform subdivision of each bar), and every station of a recovery runs
// over the loads of its bar. Exact by superposition (round-off aside).
struct MergeKey {
  Id bar;
  int kind, frame;
  bool relative;
  double a, b;
  bool operator==(const MergeKey& o) const {
    return bar == o.bar && kind == o.kind && frame == o.frame && relative == o.relative && a == o.a &&
           b == o.b;  // -0.0 == 0.0, as the ordered key compared them
  }
};
struct MergeKeyHash {
  std::size_t operator()(const MergeKey& k) const {
    std::size_t h = std::hash<Id>()(k.bar);
    auto mix = [&h](std::size_t v) { h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2); };
    mix(static_cast<std::size_t>(k.kind * 8 + k.frame * 2 + (k.relative ? 1 : 0)));
    mix(std::hash<double>()(k.a));  // std::hash<double>: 0.0 and -0.0 alike
    mix(std::hash<double>()(k.b));
    return h;
  }
};
std::vector<Load> merge_bar_loads(const std::vector<Load>& in) {
  std::vector<Load> out;
  out.reserve(in.size());
  // hashed (the creep and loss loads of a long history run into the tens of
  // thousands per increment); the merge order is the order of first
  // occurrence, the sums in input order — as with the ordered map before
  std::unordered_map<MergeKey, std::size_t, MergeKeyHash> at;
  at.reserve(in.size());
  for (const Load& l : in) {
    const bool mergeable = l.kind == Load::Kind::BarLinear || l.kind == Load::Kind::BarUdl ||
                           l.kind == Load::Kind::BarPoint;
    if (!mergeable) {
      out.push_back(l);
      continue;
    }
    const MergeKey key{l.bar, static_cast<int>(l.kind), static_cast<int>(l.frame), l.relative, l.a,
                       l.kind == Load::Kind::BarPoint ? 0.0 : l.b};
    auto it = at.find(key);
    if (it == at.end()) {
      at.emplace(key, out.size());
      out.push_back(l);
      continue;
    }
    Load& m = out[it->second];
    m.F += l.F;
    m.M += l.M;
    m.F2 += l.F2;
    m.M2 += l.M2;
  }
  return out;
}

// ---------------- rate-type creep memory (CreepEngine::Kelvin) ----------------
//
// Boltzmann superposition (McHenry 1943; Bazant & Wu 1974): the delayed strain
// of an aging linear-viscoelastic bar under a staircase stress history is the
// SUM over loading events of the step responses J(t, t_event) * dsigma_event.
// Each event therefore gets its OWN Kelvin chain, fitted at the age of that
// event, carrying the stress increment of that event; the total creep strain
// of an interval is the sum of every event chain's advance.

// A generalized-stress increment (N/A, My/Iy, Mz/Iz — three independent scalar
// laws, the elastic section response being s = E * e component-wise) applied
// at one date, carried by one chain fitted at that date.
struct CreepEvent {
  std::shared_ptr<const KelvinChain> chain;
  std::array<KelvinState, 3> state;
  double scale = 1.0;  // non-linear creep factor of the event (ConcreteLaw::nonlinear_fcm)
};

// Creep memory of ONE station of a bar (the stress varies along the bar, so
// each station carries its own history).
// A loading event of the exact superposition (CreepEngine::Exact): the
// generalized stress step, its age and the modulus phi is referred to there.
struct ExactEvent {
  double t = 0, e_ref = 1, scale = 1;
  std::array<double, 3> ds{0.0, 0.0, 0.0};
};
struct StationCreep {
  std::vector<CreepEvent> events;
  std::vector<ExactEvent> exact;
  // Total generalized stress the events carry; the difference between the
  // current stress and `carried` at an interval start is the new increment.
  std::array<double, 3> carried{0.0, 0.0, 0.0};
};

// ---------------- ST1's incremental method (CreepEngine::St1) ----------------
// State of one section fibre (see keta/phasing.h, CreepEngine::St1): the
// stress it carries, the amplitude S = sum dsigma_j A(t_j) / E_ref(t_j) of its
// steps, its creep strain e, and the branch it runs on.
struct St1Fibre {
  double sigma = 0, S = 0, e = 0;
  double smax = 0;  // largest |stress| seen: the scale of a significant step
  enum class Mode { Virgin, Loading, Recovery } mode = Mode::Virgin;
  double te = 0;                       // Loading: equivalent duration [days]
  double e1 = 0, target = 0;           // Recovery: start strain, target
  double t_ref = -1;                   // first loading age: the time function's
  bool reversal = false;               // reversal jump: recover towards S next
};

// The time-dependent law of an aging concrete: the generic `concrete_law`
// when set, else the EC2 laws of `concrete` referenced to `Material::E`.
bool aging(const Material& m) { return m.concrete_law || m.concrete; }
double law_phi(const Material& m, double t, double t0) {
  return m.concrete_law ? m.concrete_law->phi(t, t0) : ec2_phi(t, t0, *m.concrete);
}
double law_E_ref(const Material& m, double t0) {
  return m.concrete_law ? m.concrete_law->E_ref(t0) : m.E;
}
double law_shrinkage(const Material& m, double t) {
  return m.concrete_law ? m.concrete_law->shrinkage(t) : ec2_shrinkage(t, *m.concrete);
}
// J(t, t0) = 1/E + phi(t, t0)/E_ref(t0): the instantaneous term is the model's
// own constant modulus, as the bar stiffness uses it.
// Shrinkage increment of an interval [t1, t2] of a member cast at `cast`:
// the law counted from the age of 0.5 day on (S(a) = 0 below), on the TRUE
// ages — the ages clamped for the creep coefficient (a1 >= 0.5, a2 > a1)
// would count S(0.5) once more on every sub-step ending before 0.5 day.
double shrinkage_increment(const Material& m, double t1, double t2, double cast, double shift = 0.0) {
  auto S = [&](double age) { return age > 0.5 ? law_shrinkage(m, age - shift) : 0.0; };
  return S(t2 - cast) - S(t1 - cast);
}
ComplianceFn law_compliance(const Material& m) {
  if (!m.concrete_law) return ec2_compliance(*m.concrete, m.E);
  std::shared_ptr<const ConcreteLaw> law = m.concrete_law;
  double E = m.E;
  return [law, E](double t, double t0) { return 1.0 / E + law->phi(t, t0) / law->E_ref(t0); };
}

// Event threshold: a stress-component change below 1e-8 of max(|new|, |carried|)
// is folded into the LATEST chain instead of spawning one. 1e-8 is invisible
// next to the >= 0.5 % chain fit error, yet far above increment-solve roundoff
// (~1e-15 relative), so pure interval subdivision never spawns a chain and the
// telescoping property of the exponential algorithm is preserved.
constexpr double kKelvinEventRelTol = 1e-8;

// Shared NNLS fits: one per (material id, quantized loading age). All bars of
// one material loaded at one phase date share a single fit. The age is
// quantized on a log grid of 4096 steps per decade (0.056 % resolution — far
// below the fit error, so snapping the fit age to the grid representative is
// harmless) and the fit is made AT the representative age, which makes the
// cache deterministic regardless of bar order.
struct KelvinFitCache {
  std::map<std::pair<Id, long long>, std::shared_ptr<const KelvinChain>> chains;
  int fits = 0;

  std::shared_ptr<const KelvinChain> get(const Material& mat, double age, int n_units) {
    long long q = std::llround(std::log10(age) * 4096.0);
    auto key = std::make_pair(mat.id, q);
    auto it = chains.find(key);
    if (it != chains.end()) return it->second;
    double age_rep = std::pow(10.0, static_cast<double>(q) / 4096.0);
    KelvinFitOptions ko;
    ko.n_units = n_units;
    // keep the retardation window pinned at 1e-2 .. 1e5 days whatever the
    // unit count (8 units -> the 1-decade default spacing, bit-identical to
    // the historical grid; 16 -> half-decade, fit < 0.1 %)
    ko.decades = 7.0 / std::max(1, n_units - 1);
    auto c = std::make_shared<KelvinChain>(
        fit_kelvin_chain(law_compliance(mat), age_rep, ko));
    ++fits;
    chains.emplace(key, c);
    return c;
  }
};

struct TendonState {
  double frac_applied = 0;  // sum of jack fractions
  double frac_eff = 0;      // after elastic/creep-coupled losses
  bool live = false;
  double t_tension = 0;  // date of the last tensioning action [days]
  double sigma_pi = 0;   // effective tendon stress just after it [Pa] (relaxation base)
  double relax_cum = 0;  // relaxation cumulated since then [Pa] (BPEL equivalent time)
  bool on_bed = false;   // pretension strand held by its bed (TensionAction::Mode::Bed)
  bool unbonded = false; // tensioned NON INJECTE and not grouted yet (TensionAction::unbonded)
  // Phase::inject_partial: the grouted intervals [s0, s1] of an unbonded
  // tendon (bonded there, sliding elsewhere)
  std::vector<std::pair<double, double>> grouted;
  // PhasingOptions::local_losses: force lost at each point of the tendon's
  // grid (TendonGrid), and the BPEL relaxation cumulated there
  std::vector<double> L, relax_k;
};

// The stations of a tendon's chain bars along its arc (between its anchors):
// where the local losses are evaluated and kept.
struct TendonGrid {
  std::vector<double> s, shape, ey, ez;  // arc abscissa, jack-profile force, eccentricities
  std::vector<Id> bar;
  std::vector<std::size_t> st;           // station index in the bar's record
  // per point, filled at the first loss pass: the bar, its dense index and
  // the section at the station (the hot loss loop then does no lookup)
  bool cached = false;
  std::vector<const Bar*> bp;
  std::vector<int> bix;
  std::vector<double> x;
  std::vector<SectionProps> props;
};

// Mean force of a tendon's friction / draw-in profile over its anchored
// length: the reference that turns a jack fraction into a mean steel stress.
double profile_mean_force(const ForceProfile& fp, const Tendon& t) {
  double area = 0, len = 0;
  for (std::size_t i = 1; i < fp.s.size(); ++i) {
    if (fp.P[i - 1] <= 0 || fp.P[i] <= 0) continue;  // outside the anchors
    double ds = fp.s[i] - fp.s[i - 1];
    area += 0.5 * (fp.P[i - 1] + fp.P[i]) * ds;
    len += ds;
  }
  return len > 0 ? area / len : t.jack_force;
}

// Zero-shaped full-layout result.
CaseResult zero_result(const Model& model, int stations,
                       const std::map<Id, std::vector<double>>* stations_abs,
                       const std::vector<Support>& more) {
  CaseResult r;
  r.node_disp.assign(model.nodes.size(), Vec6::Zero());
  r.shells.resize(model.shells.size());
  for (std::size_t i = 0; i < model.shells.size(); ++i) {
    r.shells[i].shell = model.shells[i].id;
    r.shells[i].slab = model.shells[i].slab;
  }
  // one record per supported node, in order of first occurrence
  std::unordered_set<Id> seen;
  seen.reserve(model.supports.size() + more.size());
  std::vector<Support> all = model.supports;
  all.insert(all.end(), more.begin(), more.end());
  for (const Support& s : all) {
    if (s.implicit) continue;  // modelling restraint: no reaction reported
    if (!seen.insert(s.node).second) continue;
    r.reactions.push_back({s.node, Vec6::Zero()});
  }
  r.bars.resize(model.bars.size());
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    r.bars[i].bar = model.bars[i].id;
    if (stations_abs) {  // the static solve's own abscissae (StaticOptions::stations_abs)
      auto it = stations_abs->find(model.bars[i].id);
      if (it != stations_abs->end() && !it->second.empty()) {
        for (double x : it->second) {
          Station st;
          st.x = x;
          r.bars[i].stations.push_back(st);
        }
        continue;
      }
    }
    int n = std::max(1, stations);
    r.bars[i].stations.assign(static_cast<std::size_t>(n), Station{});
    double L = model.bar_chord_length(model.bars[i]);
    for (int k = 0; k < n; ++k)
      r.bars[i].stations[static_cast<std::size_t>(k)].x = n == 1 ? L / 2 : L * k / (n - 1);
  }
  r.links.resize(model.links.size());
  for (std::size_t i = 0; i < model.links.size(); ++i) r.links[i].link = model.links[i].id;
  r.non_superposable = true;
  return r;
}

void accumulate(CaseResult& acc, const Model& model, const CaseResult& inc,
                const Model& work) {
  // node displacements: same node set
  for (std::size_t i = 0; i < acc.node_disp.size(); ++i) acc.node_disp[i] += inc.node_disp[i];
  // reactions: match by node id (indexed: the support count of a long deck
  // squared, at every increment, is not negligible)
  {
    std::unordered_map<Id, std::size_t> ix;
    ix.reserve(inc.reactions.size());
    for (std::size_t q = 0; q < inc.reactions.size(); ++q) ix.emplace(inc.reactions[q].node, q);
    for (auto& ra : acc.reactions) {
      const auto hit = ix.find(ra.node);
      if (hit != ix.end()) ra.R += inc.reactions[hit->second].R;
    }
  }
  // bars: match by id (ids are unique) through an index of the increment —
  // a quadratic search here dominated the long creep histories
  std::unordered_map<Id, std::size_t> inc_ix;
  inc_ix.reserve(inc.bars.size());
  for (std::size_t q = 0; q < inc.bars.size(); ++q) inc_ix.emplace(inc.bars[q].bar, q);
  for (auto& ba : acc.bars) {
    const auto hit = inc_ix.find(ba.bar);
    if (hit == inc_ix.end()) continue;
    const BarResult& bi = inc.bars[hit->second];
    ba.end_forces_local += bi.end_forces_local;
    for (std::size_t k = 0; k < ba.stations.size() && k < bi.stations.size(); ++k) {
      Station& a = ba.stations[k];
      const Station& b = bi.stations[k];
      a.N += b.N;
      a.Vy += b.Vy;
      a.Vz += b.Vz;
      a.T += b.T;
      a.My += b.My;
      a.Mz += b.Mz;
      a.defl += b.defl;
      a.rot += b.rot;
    }
    // soil pressure of a bedded bar (keta/bef.h)
    if (!bi.soil_p.empty()) {
      ba.soil_p.resize(ba.stations.size(), Vec3::Zero());
      ba.soil_m.resize(ba.stations.size(), 0.0);
      for (std::size_t k = 0; k < ba.stations.size() && k < bi.soil_p.size(); ++k) {
        ba.soil_p[k] += bi.soil_p[k];
        if (k < bi.soil_m.size()) ba.soil_m[k] += bi.soil_m[k];
      }
    }
  }
  // elastic links: match by id
  for (auto& la : acc.links)
    for (const LinkResult& li : inc.links)
      if (li.link == la.link) {
        la.deformation += li.deformation;
        la.force += li.force;
      }
  // Facet shells: ACTIVE FROM THE FIRST PHASE by design (a phase plan names
  // BARS; slabs are cast with the initial layout — activating shells per phase
  // is out of scope and would need its own stress-free-cast bookkeeping).
  // Every working model carries the full shell set, so the records are
  // index-parallel: superpose the eight tensor components, then RECOMPUTE the
  // principals (never combine principals linearly).
  for (std::size_t s = 0; s < acc.shells.size(); ++s) {
    ShellResultants& a = acc.shells[s].r;
    const ShellResultants& b = inc.shells[s].r;
    a.nx += b.nx;
    a.ny += b.ny;
    a.nxy += b.nxy;
    a.mx += b.mx;
    a.my += b.my;
    a.mxy += b.mxy;
    a.qx += b.qx;
    a.qy += b.qy;
    a.recompute_principals();
  }
  (void)model;
  (void)work;
}

// every superposable field of a result times f (what accumulate adds)
void scale_result(CaseResult& r, double f) {
  for (Vec6& v : r.node_disp) v *= f;
  for (Reaction& q : r.reactions) q.R *= f;
  for (BarResult& b : r.bars) {
    b.end_forces_local *= f;
    for (Station& st : b.stations) {
      st.N *= f;
      st.Vy *= f;
      st.Vz *= f;
      st.T *= f;
      st.My *= f;
      st.Mz *= f;
      st.defl *= f;
      st.rot *= f;
    }
    for (Vec3& p : b.soil_p) p *= f;
    for (double& m : b.soil_m) m *= f;
  }
  for (LinkResult& l : r.links) {
    l.deformation *= f;
    l.force *= f;
  }
  for (ShellResult& sh : r.shells) {
    ShellResultants& a = sh.r;
    for (double* x : {&a.nx, &a.ny, &a.nxy, &a.mx, &a.my, &a.mxy, &a.qx, &a.qy}) *x *= f;
    a.recompute_principals();
  }
}

// per tendon: the loads of a unit loss at each node of its grid
struct LossBasis {
  bool built = false, ok = false;
  std::vector<double> s;
  std::vector<std::vector<Load>> loads;
  // per node, the loads its unit loss actually touches (index, F, M, F2, M2)
  struct Nz {
    std::size_t q;
    Vec3 F, M, F2, M2;
  };
  std::vector<std::vector<Nz>> nz;
  std::vector<Load> zero;  // loads.front() with zero forces: the sum's start
  std::vector<int> bix;    // dense index of each load's bar
};

// Pseudo-plastic soil (PhasingOptions::pse): per bar, the committed state
// of each force component of its bed (kBedFx / kBedFy / kBedFz; kBedMx stays
// elastic), the sign of a capped reaction, and per station the local
// displacement since the bar's activation and the soil pressure it carries.
struct PseBar {
  const BarSoil* soil = nullptr;
  std::array<BedState, kBedComps> st{};
  std::array<double, kBedComps> sign{};
  std::vector<std::array<double, kBedComps>> u, p;
};
// the modulus of component c at abscissa x (the zone holding x)
double pse_k(const BarSoil& bs, int c, double x) {
  const SoilZone* z = &bs.zones.back();
  for (const SoilZone& zz : bs.zones)
    if (&zz == &bs.zones.back() || x <= zz.xl_end + 1e-9) {
      z = &zz;
      break;
    }
  return c == kBedFx ? z->kfx : c == kBedMx ? z->kmx : c == kBedFy ? z->kfy : z->kfz;
}
double pse_ff(const SoilPlastic& pl, int c) {
  return c == kBedFx ? pl.ffx : c == kBedMx ? pl.fmx : c == kBedFy ? pl.ffy : pl.ffz;
}
double station_u(const Station& s, int c) {
  return c == kBedFx ? s.defl.x() : c == kBedMx ? s.rot.x() : c == kBedFy ? s.defl.y()
                                                                       : s.defl.z();
}
// the soil pressure of component c at station q under state `st`, `sg` and
// displacement `uq`
double pse_p(const PseBar& pb, int c, BedState st, double sg, double x, double uq) {
  if (st == BedState::Inactive) return 0.0;
  if (st == BedState::Capped) return sg * pse_ff(*pb.soil->plastic, c);
  return pse_k(*pb.soil, c, x) * uq;
}

// Two tendons with the same trace, chain, anchors and loss data — everything
// the equivalent loads, the force profile and the losses are computed from.
bool same_segment(const PolySegment& a, const PolySegment& b) {
  return a.s0 == b.s0 && a.s1 == b.s1 && a.cy == b.cy && a.cz == b.cz && a.exact == b.exact &&
         a.p[0] == b.p[0] && a.p[1] == b.p[1] && a.p[2] == b.p[2] && a.p[3] == b.p[3] &&
         a.C == b.C && a.U == b.U && a.D == b.D && a.R == b.R && a.th == b.th;
}
bool same_losses(const TendonLosses& a, const TendonLosses& b) {
  return a.mu == b.mu && a.k == b.k && a.slip == b.slip && a.Ep == b.Ep && a.Ap == b.Ap &&
         a.rho1000 == b.rho1000 && a.fpk == b.fpk && a.relaxation_class == b.relaxation_class &&
         a.relaxation_bpel == b.relaxation_bpel && a.relax_teb_days == b.relax_teb_days &&
         a.thermal_loss == b.thermal_loss && a.phi_abs == b.phi_abs && a.non_simul == b.non_simul;
}
bool same_tendon(const Tendon& a, const Tendon& b) {
  if (a.id != b.id || a.bars != b.bars || a.flips != b.flips || a.profile.size() != b.profile.size() ||
      a.jack_force != b.jack_force || a.jack_both_ends != b.jack_both_ends ||
      !same_losses(a.losses, b.losses) || a.force_profile != b.force_profile ||
      a.s_begin != b.s_begin || a.s_end != b.s_end || a.exact_angle != b.exact_angle ||
      a.jack_from_end != b.jack_from_end || a.release_end_first != b.release_end_first ||
      a.free_runs != b.free_runs)
    return false;
  for (std::size_t k = 0; k < a.profile.size(); ++k)
    if (!same_segment(a.profile[k], b.profile[k])) return false;
  return true;
}
bool same_props(const SectionProps& a, const SectionProps& b) {
  return a.A == b.A && a.Avy == b.Avy && a.Avz == b.Avz && a.Iy == b.Iy && a.Iz == b.Iz &&
         a.J == b.J && a.Wely == b.Wely && a.Welz == b.Welz && a.Wply == b.Wply &&
         a.Wplz == b.Wplz && a.Ip == b.Ip;
}
bool same_section(const Section& a, const Section& b) {
  if (!same_props(a.props, b.props) || a.stress_points.size() != b.stress_points.size()) return false;
  for (std::size_t k = 0; k < a.stress_points.size(); ++k) {
    const StressPoint &p = a.stress_points[k], &q = b.stress_points[k];
    if (p.y != q.y || p.z != q.z || p.tau_factor_y != q.tau_factor_y || p.tau_factor_z != q.tau_factor_z)
      return false;
  }
  return true;
}
bool same_concrete(const std::optional<ConcreteTime>& a, const std::optional<ConcreteTime>& b) {
  if (a.has_value() != b.has_value()) return false;
  if (!a) return true;
  return a->fck == b->fck && a->rh == b->rh && a->h0 == b->h0 && a->cement_s == b->cement_s &&
         a->t_shrink0 == b->t_shrink0 && a->autogenous == b->autogenous;
}
bool same_material(const Material& a, const Material& b) {
  return a.E == b.E && a.nu == b.nu && a.rho == b.rho && a.alpha == b.alpha && a.fy == b.fy &&
         a.damping == b.damping && same_concrete(a.concrete, b.concrete) &&
         (a.concrete_law != nullptr) == (b.concrete_law != nullptr);
}
bool same_bar(const Bar& a, const Bar& b) {
  if (a.id != b.id || a.ni != b.ni || a.nj != b.nj || a.kind != b.kind || a.roll != b.roll ||
      a.releases != b.releases || a.ecc[0] != b.ecc[0] || a.ecc[1] != b.ecc[1] ||
      a.soil.kfx != b.soil.kfx || a.soil.kfy != b.soil.kfy || a.soil.kfz != b.soil.kfz ||
      a.soil.kmx != b.soil.kmx || a.cable.has_value() != b.cable.has_value())
    return false;
  if (a.cable)
    return a.cable->pretension == b.cable->pretension &&
           a.cable->unstressed_length == b.cable->unstressed_length && a.cable->ea == b.cable->ea &&
           a.cable->w == b.cable->w && a.cable->end_tension == b.cable->end_tension;
  return true;
}

// KETA_TIMING: where a phasing spends its time
struct PhTimer {
  double solve = 0, creep = 0, loss = 0, t_merge = 0, t_lower = 0, t_back = 0, t_recover = 0;
  long nsolve = 0, nk1 = 0, nk1it = 0, ngrid = 0;
  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  void print() const {
    if (!std::getenv("KETA_TIMING")) return;
    std::fprintf(stderr,
                 "PHASAGE: %.2f s dont %ld resolutions %.2f s, sous-pas de fluage %.2f s, pertes %.2f s "
                 "(%ld appels, %ld iterations, %ld points de grille)\n",
                 std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), nsolve,
                 solve, creep, loss, nk1, nk1it, ngrid);
    std::fprintf(stderr, "PHASAGE resolutions: fusion %.2f s, abaissement %.2f s, descente %.2f s, recouvrement %.2f s\n",
                 t_merge, t_lower, t_back, t_recover);
  }
};

// The stations of the layout a tendon's local-loss grid runs over (the chain
// bars' records of `layout`).
TendonGrid build_grid(const Model& model, const Tendon& t, const ForceProfile& fp,
                      const CaseResult& layout) {
  TendonGrid g;
  const double Lchain = tendon_chain_length(model, t);
  const double sb = t.s_begin, se = t.s_end < 0 ? Lchain : t.s_end;
  double s0 = 0;
  for (std::size_t bi = 0; bi < t.bars.size(); ++bi) {
    const Bar& b = model.bar(t.bars[bi]);
    const double Lb = model.bar_chord_length(b);
    const bool flip = bi < t.flips.size() && t.flips[bi];
    const BarResult& br = layout.bars[static_cast<std::size_t>(model.bar_index(b.id))];
    const std::size_t n = br.stations.size();
    for (std::size_t q = 0; q < n; ++q) {
      const std::size_t k = flip ? n - 1 - q : q;  // ascending along the chain
      const double s = s0 + (flip ? Lb - br.stations[k].x : br.stations[k].x);
      if (s < sb - 1e-9 || s > se + 1e-9) continue;
      g.s.push_back(s);
      g.shape.push_back(fp.at(s));
      g.ey.push_back(tendon_ey(t, s));
      g.ez.push_back(tendon_ez(t, s));
      g.bar.push_back(b.id);
      g.st.push_back(k);
    }
    s0 += Lb;
  }
  return g;
}

}  // namespace

bool equal_support(const Support& a, const Support& b) {
  if (a.node != b.node || a.fixed != b.fixed || a.spring != b.spring || a.implicit != b.implicit ||
      a.oneside != b.oneside || a.cap != b.cap || a.rot != b.rot || a.decol_axis != b.decol_axis ||
      a.decol_sign != b.decol_sign || a.spring6.has_value() != b.spring6.has_value())
    return false;
  return !a.spring6 || *a.spring6 == *b.spring6;
}

FactorCache::Stats PhasedFactorCache::stats() const {
  FactorCache::Stats agg;
  for (const auto& [key, fc] : slots_) {
    agg.full += fc.stats().full;
    agg.symbolic_reused += fc.stats().symbolic_reused;
    agg.factor_reused += fc.stats().factor_reused;
  }
  return agg;
}

// ============================================================ the session --
struct PhasingSession::Impl {
  Model model;
  PhasingOptions opt;
  // The working models copy this one: the model without its load cases,
  // combinations, tendons (the plan lowers them itself through its
  // TensionActions; a tendon running into bars cast later would not even
  // validate on a substructure) and bars (each phase's active ones).
  Model model_lite;
  PhTimer ph_timer;
  // phi(t, t0) of a material is a pure function, evaluated for every
  // station and every past stress event of the superposition: memoised
  struct PhiKey {
    const Material* m;
    double t, t0;
    bool operator==(const PhiKey& o) const { return m == o.m && t == o.t && t0 == o.t0; }
  };
  struct PhiHash {
    std::size_t operator()(const PhiKey& k) const {
      std::size_t h = std::hash<const void*>()(k.m);
      h ^= std::hash<double>()(k.t) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      h ^= std::hash<double>()(k.t0) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<PhiKey, double, PhiHash> phi_memo;
  double phi_of(const Material& m, double t, double t0) {
    auto [it, fresh] = phi_memo.try_emplace(PhiKey{&m, t, t0}, 0.0);
    if (fresh) it->second = law_phi(m, t, t0);
    return it->second;
  }

  // ---- the state carried from phase to phase ----
  std::set<Id> active_bars;
  // the same set by dense bar index (hot loops: one load at a time), kept in
  // step with every insert / erase of active_bars
  std::vector<char> active_ix;
  std::vector<Support> active_supports;
  // degree of each node in the active structure: active bars at the node +
  // non-implicit active supports there (Phase::activate_rel asks whether a
  // node already belongs to the structure)
  std::unordered_map<Id, int> struct_deg;
  // every support the plans run so far add (the reaction slots of master)
  std::vector<Support> plan_supports;
  std::unordered_set<Id> reaction_nodes;  // the nodes master.reactions records
  // Full-layout master state (reactions of every support the plans add too).
  CaseResult master;
  std::vector<Vec12> bar_nodal_acc;  // global frame, by dense bar index
  // node-to-node links: the active ones, and each elastic link's accumulated
  // nodal forces (global, spring-back on deactivation)
  std::set<Id> active_links, active_rigid;
  std::map<Id, Vec12> link_nodal_acc;
  // composite bars (keta/composite.h): per bar, per part — active, casting
  // date, the accumulated resultants (N, My, Mz about the part's own
  // centroid) per station of master, the exact creep memory per station
  struct PartState {
    bool active = false;
    double cast = 0;
    std::vector<Vec3> res;
    std::vector<StationCreep> mem;
  };
  std::map<Id, std::vector<PartState>> parts;
  std::set<std::pair<Id, int>> deferred_parts;
  // The creep step of the parts about to be solved (one sub-step of the
  // interval), per composite bar. Step-by-step method with the stress
  // increment of the step applied at its start (Bazant 1972, first order in
  // the log-time steps, unconditionally stable): over [a1, a2] the strain of
  // a part is
  //     d_eps = d_sigma / E + d_sigma phi(a2, a1) / E_ref(a1) + d_hist,
  // d_hist the creep over the step of every earlier stress event (and the
  // free shrinkage), so the part answers the step with the effective modulus
  //     E_eff = 1 / (1 / E + phi(a2, a1) / E_ref(a1))
  // to the strain beyond d_hist; its d_sigma is then an event at a1. (The
  // explicit variant — d_sigma crept from the next step on — oscillates
  // once phi(a2, a1) of a part almost fully restrained by the others, the
  // slab's own bending, exceeds 2.)
  struct CompPending {
    std::vector<std::vector<std::array<double, 3>>> part;  // d_hist: part -> station
    std::vector<std::array<double, 3>> star;               // member imposed strains: station
    std::vector<double> E;                                 // moduli of the step (E_eff of the aging parts)
    std::vector<char> creeps;                              // part takes an event at a1
    std::vector<double> a1, e_ref;
  };
  std::map<Id, CompPending> comp_pending;

  void init_parts(Id b, double cast) {
    const auto it = model.composite.find(b);
    if (it == model.composite.end()) return;
    const std::size_t n = master.bars[static_cast<std::size_t>(model.bar_index(b))].stations.size();
    std::vector<PartState>& ps = parts[b];
    ps.assign(it->second.parts.size(), PartState{});
    for (std::size_t i = 0; i < ps.size(); ++i) {
      ps[i].active = !deferred_parts.count({b, static_cast<int>(i)});
      ps[i].cast = cast;
      ps[i].res.assign(n, Vec3::Zero());
    }
  }
  // the moduli of a composite bar's parts (0 = not cast); those of the
  // pending creep step while one is being solved
  std::vector<double> part_E(Id b) const {
    const auto pend = comp_pending.find(b);
    if (pend != comp_pending.end()) return pend->second.E;
    const std::vector<SectionPart>& q = model.composite.at(b).parts;
    const std::vector<PartState>& ps = parts.at(b);
    std::vector<double> E(q.size(), 0.0);
    for (std::size_t i = 0; i < q.size(); ++i)
      if (ps[i].active) E[i] = model.material(q[i].material).E;
    return E;
  }
  // One solved increment on the composite bars: the part resultants from
  // its generalized strains (the imposed ones of a creep step included, at
  // the step's moduli), the creep step's stress increments recorded as
  // events at its start, then the bar's record moved from the centroid of the
  // increment's member to the node line (results superpose across phases).
  void composite_increment(CaseResult& inc) {
    for (auto& [bid, ps] : parts) {
      if (!active_bars.count(bid)) continue;
      BarResult* br = nullptr;
      for (BarResult& x : inc.bars)
        if (x.bar == bid) br = &x;
      if (!br) continue;
      const std::vector<SectionPart>& q = model.composite.at(bid).parts;
      const CompositeState cs = composite_state(model, bid, part_E(bid));
      const auto pend = comp_pending.find(bid);
      const bool creep = pend != comp_pending.end();
      const double EA = cs.E_ref * cs.props.A, EIy = cs.E_ref * cs.props.Iy, EIz = cs.E_ref * cs.props.Iz;
      const std::size_t n = br->stations.size();
      for (std::size_t k = 0; k < n; ++k) {
        Station& st = br->stations[k];
        std::array<double, 3> tot{st.N / EA, EIy > 0 ? st.My / EIy : 0.0, EIz > 0 ? st.Mz / EIz : 0.0};
        if (creep)
          for (std::size_t c = 0; c < 3; ++c) tot[c] += pend->second.star[k][c];
        for (std::size_t i = 0; i < q.size(); ++i) {
          if (!ps[i].active || k >= ps[i].res.size()) continue;
          std::array<double, 3> d{0.0, 0.0, 0.0};
          if (creep) d = pend->second.part[i][k];
          const double Ei = cs.E[i], z = q[i].z - cs.zc;
          const Vec3 dr(Ei * q[i].props.A * (tot[0] + tot[1] * z - d[0]), Ei * q[i].props.Iy * (tot[1] - d[1]),
                        Ei * q[i].props.Iz * (tot[2] - d[2]));
          ps[i].res[k] += dr;
          // the step's stress increment: an event at its start, carried
          if (creep && pend->second.creeps[i] && (along || k == n / 2)) {
            const SectionProps& sp = q[i].props;
            StationCreep& sc = ps[i].mem[along ? k : 0];
            ExactEvent ev;
            ev.t = pend->second.a1[i];
            ev.e_ref = pend->second.e_ref[i];
            ev.ds = {dr(0) / sp.A, sp.Iy > 0 ? dr(1) / sp.Iy : 0.0, sp.Iz > 0 ? dr(2) / sp.Iz : 0.0};
            sc.exact.push_back(ev);
            for (std::size_t c = 0; c < 3; ++c) sc.carried[c] += ev.ds[c];
          }
        }
        // to the node line: r = (0, 0, zc) local
        st.My += cs.zc * st.N;
        st.T -= cs.zc * st.Vy;
        st.defl.x() -= st.rot.y() * cs.zc;
        st.defl.y() += st.rot.x() * cs.zc;
      }
      Vec12& f = br->end_forces_local;
      for (int e : {0, 6}) {
        f(e + 4) += cs.zc * f(e);
        f(e + 3) -= cs.zc * f(e + 1);
      }
    }
    comp_pending.clear();
  }

  // The creep step [t1, t2] of the parts of composite bar `bid` (see
  // CompPending): the stress events since the last step (loads applied at
  // the phase date) recorded at a1, d_hist per part and station, the step
  // moduli, and the member's imposed strains — the restraint force of d_hist
  // about the centroid of the step's member over its rigidities. Appends the
  // member's loads; the step goes to `out`.
  void composite_creep(Id bid, double t1, double t2, std::vector<Load>& loads,
                       std::map<Id, CompPending>& out) {
    const std::vector<SectionPart>& q = model.composite.at(bid).parts;
    std::vector<PartState>& ps = parts.at(bid);
    const BarResult& br = master.bars[static_cast<std::size_t>(model.bar_index(bid))];
    const std::size_t n = br.stations.size();
    if (n == 0) return;
    CompPending cp;
    cp.part.assign(q.size(), std::vector<std::array<double, 3>>(n, {0.0, 0.0, 0.0}));
    cp.star.assign(n, {0.0, 0.0, 0.0});
    cp.E = part_E(bid);
    cp.creeps.assign(q.size(), 0);
    cp.a1.assign(q.size(), 0.0);
    cp.e_ref.assign(q.size(), 1.0);
    bool any = false;
    const std::size_t kb = along ? 0 : n / 2, ke = along ? n : n / 2 + 1;
    for (std::size_t i = 0; i < q.size(); ++i) {
      const Material& mat = model.material(q[i].material);
      if (!ps[i].active || !aging(mat)) continue;
      const double cast_t = ps[i].cast;
      // the laws start at the age of 0.5 day: a sub-step that ends before
      // has nothing (its clamped ages would count phi(0.5+, 0.5) again)
      if (t2 - cast_t <= 0.5) continue;
      require(opt.creep_engine == CreepEngine::Exact,
              "phasing: the creep of composite parts runs on the exact engine (CreepEngine::Exact)");
      require(!(mat.concrete_law && mat.concrete_law->nonlinear_fcm),
              "phasing: non-linear creep is not supported in composite parts");
      const double a1 = std::max(0.5, t1 - cast_t), a2 = std::max(a1 + 1e-6, t2 - cast_t);
      const double deps_sh = shrinkage_increment(mat, t1, t2, cast_t);
      const double e_ref = law_E_ref(mat, a1);
      cp.creeps[i] = 1;
      cp.a1[i] = a1;
      cp.e_ref[i] = e_ref;
      cp.E[i] = 1.0 / (1.0 / cp.E[i] + phi_of(mat, a2, a1) / e_ref);
      ps[i].mem.resize(along ? n : 1);
      const SectionProps& sp = q[i].props;
      for (std::size_t k = kb; k < ke; ++k) {
        const Vec3& r = ps[i].res[k];
        const std::array<double, 3> s{r(0) / sp.A, sp.Iy > 0 ? r(1) / sp.Iy : 0.0, sp.Iz > 0 ? r(2) / sp.Iz : 0.0};
        StationCreep& sc = ps[i].mem[along ? k : 0];
        // the loads since the last step: an event at a1
        ExactEvent ev;
        ev.t = a1;
        ev.e_ref = e_ref;
        bool step = false;
        for (std::size_t c = 0; c < 3; ++c) {
          ev.ds[c] = s[c] - sc.carried[c];
          step = step || ev.ds[c] != 0;
        }
        if (step) sc.exact.push_back(ev);
        sc.carried = s;
        std::array<double, 3>& de = cp.part[i][k];
        for (const ExactEvent& x : sc.exact) {
          const double dph = phi_of(mat, a2, x.t) - (a1 > x.t ? phi_of(mat, a1, x.t) : 0.0);
          for (std::size_t c = 0; c < 3; ++c) de[c] += x.ds[c] * dph / x.e_ref;
        }
        de[0] += deps_sh;
      }
      if (!along)
        for (std::size_t k = 0; k < n; ++k) cp.part[i][k] = cp.part[i][n / 2];
      any = true;
    }
    if (!any) return;
    const CompositeState cs = composite_state(model, bid, cp.E);
    for (std::size_t k = 0; k < n; ++k) {
      double N = 0, My = 0, Mz = 0;
      for (std::size_t i = 0; i < q.size(); ++i) {
        if (!ps[i].active) continue;
        const std::array<double, 3>& d = cp.part[i][k];
        const double Ei = cs.E[i];
        N += Ei * q[i].props.A * d[0];
        My += Ei * (q[i].props.A * d[0] * (q[i].z - cs.zc) + q[i].props.Iy * d[1]);
        Mz += Ei * q[i].props.Iz * d[2];
      }
      cp.star[k] = {N / (cs.E_ref * cs.props.A), cs.props.Iy > 0 ? My / (cs.E_ref * cs.props.Iy) : 0.0,
                    cs.props.Iz > 0 ? Mz / (cs.E_ref * cs.props.Iz) : 0.0};
    }
    // lowered as the existing bar creep: the curvatures negated (an imposed
    // kappa = +M/EI bows the beam against M)
    auto lowered = [](std::array<double, 3> v) {
      v[1] = -v[1];
      v[2] = -v[2];
      return v;
    };
    auto nonzero = [](const std::array<double, 3>& v) { return v[0] != 0 || v[1] != 0 || v[2] != 0; };
    if (!along || n == 1) {
      const std::array<double, 3> de = lowered(cp.star[n / 2]);
      if (nonzero(de)) loads.push_back(Load::prestrain(bid, de[0], de[1], de[2]));
    } else {
      const double g = 0.5 / std::sqrt(3.0);
      for (std::size_t k = 0; k + 1 < n; ++k) {
        const double xa = br.stations[k].x, xb = br.stations[k + 1].x, dx = xb - xa;
        if (dx <= 0) continue;
        const std::array<double, 3> da = lowered(cp.star[k]), db = lowered(cp.star[k + 1]);
        for (double t : {0.5 - g, 0.5 + g}) {
          std::array<double, 3> m{0, 0, 0};
          for (std::size_t c = 0; c < 3; ++c) m[c] = 0.5 * dx * (da[c] + t * (db[c] - da[c]));
          if (nonzero(m)) loads.push_back(Load::prestrain_point(bid, xa + t * dx, m[0], m[1], m[2]));
        }
      }
    }
    out[bid] = std::move(cp);
  }
  // ---- facet shells (active from the first phase, cast at t = 0) ----
  // The creep of a facet of aging concrete: Boltzmann superposition of its
  // stress events, each kept as the element displacement whose strain field
  // is the event's elastic strain at the material's E (membrane, curvature,
  // transverse shear: the whole field), so the creep strain over [a1, a2]
  //     sum_j d_j E (phi(a2, t_j) - phi(a1, t_j)) / E_ref(t_j)
  // is imposed EXACTLY as such a displacement d* (Load::ShellImposed: f =
  // Ke d*, recovered from d - d*), the free shrinkage as the in-plane field
  // eps_cs (x, y). The step is the implicit one of the section parts: the
  // facets answer it with E_eff = 1 / (1 / E + phi(a2, a1) / E_ref(a1)) and
  // its stress increment is an event at a1.
  struct ShellEvent {
    double t = 0, e_ref = 1;
    Vec18 d = Vec18::Zero();
  };
  std::vector<std::vector<ShellEvent>> shell_mem;  // model.shells order
  std::vector<ShellElem> shell_els;                // the model's facets
  std::vector<Vec18> shell_dstar;                  // the imposed field of the increment being solved
  std::map<Id, double> shell_E;                    // material -> modulus of its facets in the step
  std::map<Id, std::pair<double, double>> shell_step;  // material -> (a1, E_ref(a1)) of the step
  double cur_time = 0;                             // the date of the increments
  // shell tendons (Model::shell_tendons): jacked by TensionAction::shell,
  // bonded — every recorded increment strains them (the mean concrete strain
  // along the tendon, shell_tendon_mean_strain), the loss loads lose again,
  // iterated; relaxation EN 1992-1-1 by equivalent time on the current stress
  struct ShellTendonState {
    double frac_applied = 0, frac_eff = 0;
    bool live = false;
    double t_tension = 0, sigma_pi = 0, relax_cum = 0;
  };
  std::map<Id, ShellTendonState> sh_t;
  std::map<Id, double> sh_mean;  // mean force of the P0 profile
  bool in_shell_loss = false;
  bool shell_aging(std::size_t i) const { return aging(model.material(model.shells[i].material)); }
  const ShellTendon& shell_tendon(Id id) const {
    for (const ShellTendon& t : model.shell_tendons)
      if (t.id == id) return t;
    fail("phasing: no shell tendon " + std::to_string(id));
  }
  // an increment's events (recorded increments only)
  void shell_increment(const CaseResult& inc) {
    if (inc.node_disp.empty()) return;
    for (std::size_t i = 0; i < shell_mem.size(); ++i) {
      if (!shell_aging(i)) continue;
      const ShellTri& sh = model.shells[i];
      const Material& mat = model.material(sh.material);
      Vec18 d;
      for (int c = 0; c < 3; ++c)
        d.segment<6>(6 * c) = inc.node_disp[static_cast<std::size_t>(shell_els[i].nodes[static_cast<std::size_t>(c)])];
      if (!shell_dstar.empty()) d -= shell_dstar[i];
      const auto e = shell_E.find(sh.material);
      const double Euse = e != shell_E.end() ? e->second : mat.E;
      ShellEvent ev;
      const auto st = shell_step.find(sh.material);
      if (st != shell_step.end()) {
        ev.t = st->second.first;
        ev.e_ref = st->second.second;
      } else {
        ev.t = std::max(0.5, cur_time);
        ev.e_ref = law_E_ref(mat, ev.t);
      }
      ev.d = (Euse / mat.E) * d;
      if (ev.d.squaredNorm() == 0) continue;
      std::vector<ShellEvent>& m = shell_mem[i];
      if (!m.empty() && m.back().t == ev.t && m.back().e_ref == ev.e_ref)
        m.back().d += ev.d;
      else
        m.push_back(ev);
    }
  }
  // the loss loads of the live shell tendons for an increment (bonded:
  // Ep Ap times the mean concrete strain along the tendon)
  std::vector<Load> shell_losses(const CaseResult& inc) {
    std::vector<Load> out;
    if (inc.node_disp.empty()) return out;
    for (const ShellTendon& t : model.shell_tendons) {
      ShellTendonState& st = sh_t[t.id];
      if (!st.live || t.Ap <= 0 || t.P0 <= 0) continue;
      const double dP = t.Ep * t.Ap * shell_tendon_mean_strain(model, shell_els, t, inc.node_disp);
      if (std::fabs(dP) <= 1e-13 * t.P0) continue;
      const double dfrac = dP / sh_mean.at(t.id);
      st.frac_eff += dfrac;
      out.push_back(Load::shell_prestress(t.id, dfrac));
    }
    return out;
  }

  std::unordered_map<Id, Vec6> reaction_acc;
  std::unordered_map<Id, TendonState> tstate;
  bool along = false;  // see PhasingOptions::creep_along_bar
  bool local = false;  // local losses: the grid of every tendon on the stations of the layout
  // tendon id -> force a jack fraction stands for in the relaxation base: the
  // mean of its friction / draw-in profile (along the bar), else P0
  std::unordered_map<Id, double> mean_force;
  // the friction / draw-in profile of each tendon, computed once (a pure
  // function of the trace and the chain; ~10^5 trace evaluations a tendon)
  std::unordered_map<Id, ForceProfile> force_profile;
  std::unordered_map<Id, TendonGrid> tgrid;
  std::unordered_map<Id, LossBasis> loss_basis;
  // A tendon's equivalent loads are proportional to its force fraction (the
  // trace and the force profile are fixed): those of the full force, once
  // per tendon, scaled for each uniform loss of the passes and intervals.
  std::unordered_map<Id, std::vector<Load>> unit_loads;
  // per-station, per-event creep memory (Kelvin / exact engines)
  std::unordered_map<Id, std::vector<StationCreep>> bar_creep;
  // per-station, per-fibre state (St1 engine)
  std::unordered_map<Id, std::vector<std::vector<St1Fibre>>> bar_st1;
  KelvinFitCache fit_cache;
  std::unordered_map<Id, double> cast_time;  // bar id -> casting date [days]
  std::map<Id, PseBar> pse;
  std::vector<PhaseResult> results;
  // PhasingOptions::origin_ledgers: the cumulated state split by origin
  // (0 = S(G), 1 = S(Po), 2 = S(dP)), every increment accumulated into the
  // ledger of `origin_now` as well as into master, and every direct write of
  // master mirrored, so the three always sum to master; `led_invalid` names
  // the first event whose increment mixed the origins.
  bool led_on = false;
  std::array<CaseResult, 3> led;
  int origin_now = 0;
  std::string led_invalid;
  void for_states(const std::function<void(CaseResult&)>& f) {
    f(master);
    if (led_on)
      for (CaseResult& l : led) f(l);
  }
  void invalidate_ledgers(const std::string& why) {
    if (led_on && led_invalid.empty()) led_invalid = why;
  }
  CaseResult isostatic_state();
  // catenary cables (Phase::cable_actions): the state of each cable
  // bar (L0 is its only state), and the event the next increment carries
  struct CableSt {
    bool active = false, seen = false;
    double L0 = 0, EA = 0, w = 0;
  };
  std::map<Id, CableSt> cab;
  std::map<Id, double> cab_next;  // L0 after the next increment
  // a force target of the next increment, met at its EQUILIBRIUM (the
  // anchors move with the re-tension): kind 0 horizontal, 1 end-j tension
  std::map<Id, std::pair<int, double>> cab_target;
  std::set<Id> cab_new, cab_off;  // activated / ruptured by the next increment
  bool cab_live() const {
    for (const auto& [b, c] : cab)
      if (c.active) return true;
    return !cab_new.empty();
  }
  Vec3 node_now(Id n) const {
    return model.node(n).x + master.node_disp[static_cast<std::size_t>(model.node_index(n))].head<3>();
  }
  CableSt& cable_init(Id bid) {
    CableSt& c = cab[bid];
    if (!c.seen) {
      const Bar& b = model.bar(bid);
      require(b.cable.has_value(), "phasing: cable bar " + std::to_string(bid) + " without CableProps");
      const auto [EA, w] = cable_ea_w(model, b);
      c.EA = EA;
      c.w = w;
      // the CableProps routes of resolve_cable_l0, at the anchors AS THEY
      // STAND when the cable is installed (the structure may have moved)
      const CableProps& cp = *b.cable;
      const Vec3 d = node_now(b.nj) - node_now(b.ni);
      const double dxh = std::hypot(d.x(), d.y());
      if (cp.unstressed_length > 0) c.L0 = cp.unstressed_length;
      else if (cp.end_tension > 0) c.L0 = catenary_l0_for_end_tension(dxh, d.z(), EA, w, cp.end_tension);
      else {
        require(cp.pretension > 0, "phasing: cable " + std::to_string(bid) + ": give unstressed_length or a positive pretension");
        c.L0 = catenary_l0_for_tension(dxh, d.z(), EA, w, cp.pretension);
      }
      c.seen = true;
    }
    return c;
  }
  CaseResult rupture_peak;
  bool rupture_peak_set = false;
  // second-order regime (Phase::regime): the displacement of each
  // bar's nodes at its stress-free activation, the injected imperfection
  std::map<Id, Vec12> bar_ref;
  std::map<Id, Vec6> imperf;
  bool so_now = false;
  Vec12 bar_d(const Bar& b) const {  // displacement since activation + imperfection
    Vec12 d = Vec12::Zero();
    if (master.node_disp.empty()) return d;
    d.head<6>() = master.node_disp[static_cast<std::size_t>(model.node_index(b.ni))];
    d.tail<6>() = master.node_disp[static_cast<std::size_t>(model.node_index(b.nj))];
    const auto r = bar_ref.find(b.id);
    if (r != bar_ref.end()) d -= r->second;
    const auto i1 = imperf.find(b.ni), i2 = imperf.find(b.nj);
    if (i1 != imperf.end()) d.head<6>() += i1->second;
    if (i2 != imperf.end()) d.tail<6>() += i2->second;
    return d;
  }
  double prev_time = 0;
  std::size_t phases_run = 0;
  PhasedFactorCache own_cache;
  PhasedFactorCache* cache_ = nullptr;

  Impl(const Model& m, const PhasingOptions& o, PhasedFactorCache* cache)
      : model(m), opt(o), cache_(cache ? cache : &own_cache) {
    model.validate();
    model_lite = model;
    model_lite.cases.clear();
    model_lite.combos.clear();
    model_lite.tendons.clear();
    model_lite.bars.clear();
    along = opt.creep_along_bar;
    local = opt.local_losses && along;
    active_ix.assign(model.bars.size(), 0);
    bar_nodal_acc.assign(model.bars.size(), Vec12::Zero());
    // the pseudo-plastic beds (PhasingOptions::pse): the bars carrying the
    // exact bed with a non-linear law, elastic until an increment says
    if (opt.pse)
      for (const auto& [bid, bs] : *opt.pse)
        if (bs.plastic && bs.plastic->nonlinear() && !bs.zones.empty() && model.tapered.count(bid))
          pse[bid].soil = &bs;
  }

  // The pristine state: `from_start` active and cast at t = 0, the model's
  // supports, an empty master.
  void start(const std::set<Id>& from_start) {
    for (Id b : from_start) {
      require(model.bar_index(b) >= 0, "phasing: unknown bar " + std::to_string(b));
      if (model.bar(b).kind == BarKind::Cable) continue;  // a cable enters by name only
      activate_bar(b);
      cast_time[b] = 0.0;
    }
    for (const ElasticLink& l : model.links) active_links.insert(l.id);
    for (const RigidLink& l : model.rigid_links) active_rigid.insert(l.id);
    active_supports = model.supports;
    for (const Support& s : active_supports)
      if (!s.implicit) ++struct_deg[s.node];
    master = zero_result(model, opt.stations_per_bar, opt.stations_abs, plan_supports);
    for (const auto& re : master.reactions) {
      reaction_acc[re.node] = Vec6::Zero();
      reaction_nodes.insert(re.node);
    }
    led_on = opt.origin_ledgers;
    if (led_on) led.fill(master);
    for (Id b : from_start) init_parts(b, 0.0);
    init_shells();
    for (const Tendon& t : model.tendons)
      for (Id b : t.bars)
        require(!model.composite.count(b), "phasing: tendon " + std::to_string(t.id) +
                                               " runs in composite bar " + std::to_string(b) +
                                               " — tendons in section parts are not supported");
    init_tendons();
  }

  void init_shells() {
    shell_els = shell_elems(model);
    shell_mem.resize(model.shells.size());
    for (std::size_t i = 0; i < model.shells.size(); ++i)
      if (shell_aging(i))
        require(!model.shells[i].enrich.ortho && !model.shells[i].enrich.laminate && model.shells[i].enrich.kbed == 0,
                "phasing: shell " + std::to_string(model.shells[i].id) +
                    ": the creep of an orthotropic facet or of a facet on a bed is not supported");
    for (const ShellTendon& t : model.shell_tendons) {
      sh_t.try_emplace(t.id);
      sh_mean[t.id] = shell_tendon_mean_force(t);
    }
  }

  void activate_bar(Id b) {
    active_bars.insert(b);
    active_ix[static_cast<std::size_t>(model.bar_index(b))] = 1;
    const Bar& bar = model.bar(b);
    if (!master.node_disp.empty()) {  // stress-free in the current geometry
      Vec12 r;
      r.head<6>() = master.node_disp[static_cast<std::size_t>(model.node_index(bar.ni))];
      r.tail<6>() = master.node_disp[static_cast<std::size_t>(model.node_index(bar.nj))];
      bar_ref[b] = r;
    }
    ++struct_deg[bar.ni];
    ++struct_deg[bar.nj];
  }
  void deactivate_bar(Id b) {
    active_bars.erase(b);
    bar_ref.erase(b);
    active_ix[static_cast<std::size_t>(model.bar_index(b))] = 0;
    const Bar& bar = model.bar(b);
    --struct_deg[bar.ni];
    --struct_deg[bar.nj];
  }
  void add_support(const Support& s) {
    active_supports.push_back(s);
    if (!s.implicit) ++struct_deg[s.node];
  }

  // The tendons of the model: their states (kept where they exist), the
  // reference force of a jack fraction, their local-loss grids on the
  // stations of master. The profiles are computed once (inherited by a fork).
  void init_tendons() {
    for (const Tendon& t : model.tendons) {
      tstate.try_emplace(t.id);
      mean_force[t.id] = along ? profile_mean_force(profile_of(t), t) : t.jack_force;
    }
    if (local)
      for (const Tendon& t : model.tendons) {
        TendonGrid g = build_grid(model, t, profile_of(t), master);
        const TendonState& ts = tstate.at(t.id);
        require(ts.L.empty() || ts.L.size() == g.s.size(),
                "phasing: tendon " + std::to_string(t.id) + " — its loss grid changed");
        tgrid[t.id] = std::move(g);
      }
  }
  const ForceProfile& profile_of(const Tendon& t) {
    auto it = force_profile.find(t.id);
    if (it == force_profile.end()) it = force_profile.emplace(t.id, tendon_force_profile(model, t)).first;
    return it->second;
  }
  // a loss field on a tendon's grid as (s, value) breakpoints, strictly
  // increasing (the two records of a node shared by two bars averaged)
  std::vector<std::pair<double, double>> grid_profile(Id tid, const std::vector<double>& v) const;
  // An unbonded tendon (NON INJECTE) slides in its duct against Coulomb
  // friction: see the definition.
  void smear(Id tid, std::vector<double>& dL);
  void tendon_loads_scaled(const Tendon& t, double fraction, std::vector<Load>& out,
                           const std::set<Id>* only = nullptr);
  // a result's record of bar `bid`: results list their bars in model order
  // (binary search on the model index), else a scan
  const BarResult* record_of(const CaseResult& res, Id bid) const {
    const int ix = model.bar_index(bid);
    std::size_t lo = 0, hi = res.bars.size();
    while (lo < hi) {
      const std::size_t mid = (lo + hi) / 2;
      if (model.bar_index(res.bars[mid].bar) < ix)
        lo = mid + 1;
      else
        hi = mid;
    }
    if (lo < res.bars.size() && res.bars[lo].bar == bid) return &res.bars[lo];
    for (const auto& bb : res.bars)
      if (bb.bar == bid) return &bb;
    return nullptr;
  }

  // The reaction slots of the supports a plan adds, appended to master (one
  // per node, in order of first occurrence — the layout zero_result builds).
  void add_plan_supports(const PhasePlan& plan) {
    for (const Phase& ph : plan.phases)
      for (const Support& s : ph.add_supports) {
        plan_supports.push_back(s);
        if (s.implicit || !reaction_nodes.insert(s.node).second) continue;
        for_states([&](CaseResult& c) { c.reactions.push_back({s.node, Vec6::Zero()}); });
        reaction_acc.emplace(s.node, Vec6::Zero());
      }
  }

  void run(const PhasePlan& plan) {
    hooks_check(opt.hooks);
    require(!plan.phases.empty(), "phasing: empty plan");
    ph_timer = PhTimer{};
    add_plan_supports(plan);
    if (phases_run == 0) prev_time = plan.phases.front().time;
    for (std::size_t i = 0; i < plan.phases.size(); ++i) step(plan.phases[i], i, plan.phases.size());
    ph_timer.ngrid = 0;
    for (const auto& kv : tgrid) ph_timer.ngrid += static_cast<long>(kv.second.s.size());
    ph_timer.print();
    hooks_report(opt.hooks, 1.0, "done");
  }

  // One phase on the current state.
  void step(const Phase& ph, std::size_t run_i, std::size_t run_n);
  // Whether the interval before `ph` has anything time-dependent to do: the
  // MATERIALS decide (an active aging bar not cast this phase), and so do
  // the live tendons carrying relaxation data — never a flag of the plan,
  // so a chain's prefix integrates the same laws whatever the chain adds
  // later (a continuation that brings the first concrete to an all-steel
  // prefix resumes from it).
  bool time_dependent(const Phase& ph) const {
    const std::unordered_set<Id> cast_now(ph.activate_bars.begin(), ph.activate_bars.end());
    for (const Bar& b : model.bars)
      if (active_ix[static_cast<std::size_t>(model.bar_index(b.id))] && !cast_now.count(b.id) &&
          !model.composite.count(b.id) && aging(model.material(b.material)))
        return true;
    for (const auto& [bid, ps] : parts) {
      if (!active_bars.count(bid) || cast_now.count(bid)) continue;
      const std::vector<SectionPart>& q = model.composite.at(bid).parts;
      for (std::size_t i = 0; i < ps.size(); ++i)
        if (ps[i].active && aging(model.material(q[i].material))) return true;
    }
    for (const Tendon& t : model.tendons) {
      const TendonState& ts = tstate.at(t.id);
      if (ts.live && t.losses.Ap > 0 && t.jack_force > 0 && t.losses.fpk > 0 && t.losses.rho1000 > 0)
        return true;
    }
    for (std::size_t i = 0; i < model.shells.size(); ++i)
      if (shell_aging(i)) return true;
    for (const ShellTendon& t : model.shell_tendons) {
      const auto it = sh_t.find(t.id);
      if (it != sh_t.end() && it->second.live && t.Ap > 0 && t.fpk > 0 && t.rho1000 > 0) return true;
    }
    return false;
  }

  std::string fork_reject(const Model& m2, const PhasingOptions& o2) const;
  std::unique_ptr<Impl> fork(const Model& m2, const PhasingOptions& o2) const;
};

std::vector<std::pair<double, double>> PhasingSession::Impl::grid_profile(
  Id tid, const std::vector<double>& v) const {
  const TendonGrid& g = tgrid.at(tid);
  std::vector<std::pair<double, double>> out;
  int cnt = 0;
  for (std::size_t k = 0; k < g.s.size(); ++k) {
    if (!out.empty() && g.s[k] <= out.back().first + 1e-9) {
      out.back().second = (out.back().second * cnt + v[k]) / (cnt + 1);
      ++cnt;
      continue;
    }
    out.push_back({g.s[k], v[k]});
    cnt = 1;
  }
  return out;
}

void PhasingSession::Impl::smear(Id tid, std::vector<double>& dL) {
  TendonState& ts = tstate[tid];
  if (!ts.unbonded) return;
  const TendonGrid& g = tgrid.at(tid);
  const Tendon& t = model.tendon(tid);
  // the anchored grid points, merged where they share an abscissa
  std::vector<std::vector<std::size_t>> mem;
  std::vector<double> ns, T, No;
  for (std::size_t k = 0; k < g.s.size(); ++k) {
    if (g.shape[k] <= 1e-12 * t.jack_force) continue;
    const double nold = ts.frac_eff * g.shape[k] - (k < ts.L.size() ? ts.L[k] : 0.0);
    if (!ns.empty() && g.s[k] <= ns.back() + 1e-9) {
      const double c = static_cast<double>(mem.back().size());
      mem.back().push_back(k);
      No.back() = (No.back() * c + nold) / (c + 1);
      T.back() = (T.back() * c + nold - dL[k]) / (c + 1);
      continue;
    }
    mem.push_back({k});
    ns.push_back(g.s[k]);
    No.push_back(nold);
    T.push_back(nold - dL[k]);
  }
  const std::size_t m = ns.size();
  if (m < 2) return;
  // grouted nodes (Phase::inject_partial) keep their bonded value, and a
  // segment touching one constrains nothing (status 3): the grouted zone
  // ANCHORS the free stretches next to it (the force passes to the bond over
  // a short transfer length), each free stretch then slides between its own
  // anchors — its force set by the compatibility with its own duct
  std::vector<char> gr(m, 0);
  for (std::size_t j = 0; j < m; ++j)
    for (const auto& [a, b] : ts.grouted) gr[j] = gr[j] || (ns[j] >= a - 1e-9 && ns[j] <= b + 1e-9);
  // Coulomb along the duct, exactly: dN/ds = -+ N (mu (dtheta/ds + k) +
  // phi) where the tendon slides, so a segment of friction exponent
  // a = mu (dtheta + k h) + phi h bounds the force ratio of its ends,
  //     e^-a N_j <= N_j+1 <= e^a N_j,
  // linear constraints of the projection (the linearised N_old a bound
  // was first order in a). dtheta is the angle the tangent sweeps (16
  // samples a segment: an S curve or a slope break inside it counts).
  std::vector<double> w(m, 0.0), c(m - 1, 1.0);
  auto tangent = [&](double sx) {
    Vec3 v(1.0, tendon_dey(t, sx), tendon_dez(t, sx));
    return Vec3(v / v.norm());
  };
  for (std::size_t k = 0; k + 1 < m; ++k) {
    const double h = ns[k + 1] - ns[k];
    w[k] += 0.5 * h;
    w[k + 1] += 0.5 * h;
    double dth = 0;
    Vec3 prev = tangent(ns[k] + 1e-9 * h);
    for (int q = 1; q <= 16; ++q) {
      const Vec3 cur = tangent(q == 16 ? ns[k + 1] - 1e-9 * h : ns[k] + h * q / 16.0);
      dth += std::atan2(prev.cross(cur).norm(), prev.dot(cur));
      prev = cur;
    }
    bool free_run = false;
    for (const auto& [f0, f1] : t.free_runs) free_run = free_run || (ns[k] >= f0 - 1e-9 && ns[k + 1] <= f1 + 1e-9);
    const double rub = free_run ? 0.0 : h;
    c[k] = std::exp(t.losses.mu * (dth + t.losses.k * rub) + t.losses.phi_abs * rub);
  }
  // a free stretch runs from its anchor (the grouted boundary) on: the half
  // of the boundary segment beside the grouted node belongs to it too, with
  // the trial INCREMENT of the grouted end (its strain; not its force, which
  // the bond holds) — the compatibility over the stretch's whole length
  for (std::size_t k = 0; k + 1 < m; ++k) {
    if (gr[k] == gr[k + 1]) continue;
    const std::size_t g = gr[k] ? k : k + 1, f = gr[k] ? k + 1 : k;
    const double h = 0.5 * (ns[k + 1] - ns[k]);
    const double dlf = No[f] - T[f], dlg = No[g] - T[g];
    T[f] = No[f] - (w[f] * dlf + h * dlg) / (w[f] + h);
    w[f] += h;
  }
  for (double& x : w) x = std::max(x, 1e-12);
  // segment status: 0 sticking, +1 sliding with N_j+1 = c N_j (the force
  // rises along s), -1 with N_j = c N_j+1, 2 frictionless (N_j+1 = N_j,
  // either sign). Active segments chain nodes into groups N_i = rho_i N_g;
  // each group's N_g minimises sum w (rho N_g - T)^2, and the multipliers
  // follow node by node from stationarity (>= 0 on a sliding segment).
  std::vector<int> st(m - 1, 0);
  const double tiny = 1e-12 * std::max(1.0, t.jack_force);
  for (std::size_t k = 0; k + 1 < m; ++k) {
    if (gr[k] || gr[k + 1]) st[k] = 3;
    else if (c[k] - 1.0 <= 1e-14) st[k] = 2;
    else if (T[k + 1] > c[k] * T[k]) st[k] = 1;
    else if (T[k] > c[k] * T[k + 1]) st[k] = -1;
  }
  std::vector<double> N(m, 0.0), u(m - 1, 0.0), rho(m, 1.0);
  bool settled = false;
  for (int it = 0; it < 50 + 20 * static_cast<int>(m); ++it) {
    std::fill(u.begin(), u.end(), 0.0);
    for (std::size_t j0 = 0; j0 < m;) {
      std::size_t j1 = j0;  // the group [j0, j1]
      rho[j0] = 1.0;
      while (j1 + 1 < m && st[j1] != 0 && st[j1] != 3) {
        rho[j1 + 1] = rho[j1] * (st[j1] == 1 ? c[j1] : st[j1] == -1 ? 1.0 / c[j1] : 1.0);
        ++j1;
      }
      double num = 0, den = 0;
      for (std::size_t j = j0; j <= j1; ++j) {
        num += w[j] * rho[j] * T[j];
        den += w[j] * rho[j] * rho[j];
      }
      const double Ng = num / den;
      for (std::size_t j = j0; j <= j1; ++j) N[j] = rho[j] * Ng;
      // stationarity w_j (N_j - T_j) + u_j-1 dg_j-1/dN_j + u_j dg_j/dN_j = 0
      // with g = N_k+1 - c N_k (+1), N_k - c N_k+1 (-1), N_k+1 - N_k (2)
      // (a grouted node is a group of its own: N = T)
      double carry = 0;  // u_j-1 dg_j-1/dN_j
      for (std::size_t j = j0; j < j1; ++j) {
        const double dgj = st[j] == 1 ? -c[j] : st[j] == -1 ? 1.0 : -1.0;  // d g_j / d N_j
        u[j] = -(w[j] * (N[j] - T[j]) + carry) / dgj;
        const double dgn = st[j] == 1 ? 1.0 : st[j] == -1 ? -c[j] : 1.0;  // d g_j / d N_j+1
        carry = u[j] * dgn;
      }
      j0 = j1 + 1;
    }
    bool changed = false;
    if (it < 50) {  // all the changes at once (fast when it settles)
      for (std::size_t k = 0; k + 1 < m; ++k) {
        if (st[k] == 2 || st[k] == 3) continue;
        if (st[k] != 0 && u[k] < -tiny * 1e-6) {
          st[k] = 0;
          changed = true;
        } else if (st[k] == 0 && N[k + 1] > c[k] * N[k] + tiny) {
          st[k] = 1;
          changed = true;
        } else if (st[k] == 0 && N[k] > c[k] * N[k + 1] + tiny) {
          st[k] = -1;
          changed = true;
        }
      }
    } else {
      // a cycle: one change at a time, the worst one (the most negative
      // multiplier, else the most violated segment) — the classical
      // active-set step of a strictly convex QP
      std::size_t kb = m;
      int to = 0;
      double worst = 0;
      for (std::size_t k = 0; k + 1 < m; ++k)
        if (st[k] != 0 && st[k] != 2 && st[k] != 3 && u[k] < -tiny * 1e-6 && -u[k] > worst) {
          worst = -u[k];
          kb = k;
          to = 0;
        }
      if (kb == m)
        for (std::size_t k = 0; k + 1 < m; ++k) {
          if (st[k] != 0) continue;  // (3: unconstrained, never activated)
          const double vp = N[k + 1] - c[k] * N[k], vm = N[k] - c[k] * N[k + 1];
          if (vp > tiny && vp > worst) {
            worst = vp;
            kb = k;
            to = 1;
          }
          if (vm > tiny && vm > worst) {
            worst = vm;
            kb = k;
            to = -1;
          }
        }
      if (kb < m) {
        st[kb] = to;
        changed = true;
      }
    }
    if (!changed) {
      settled = true;
      break;
    }
  }
  require(settled, "phasing: NON INJECTE tendon " + std::to_string(tid) +
                       " — the Coulomb slip states do not settle");
  for (std::size_t j = 0; j < m; ++j)
    for (std::size_t k : mem[j])
      dL[k] = ts.frac_eff * g.shape[k] - (k < ts.L.size() ? ts.L[k] : 0.0) - N[j];
}

void PhasingSession::Impl::tendon_loads_scaled(const Tendon& t, double fraction,
                                             std::vector<Load>& out, const std::set<Id>* only) {
  require(fraction >= -1.0 && fraction <= 1.0 + 1e-12, "tendon: fraction out of range");
  auto it = unit_loads.find(t.id);
  if (it == unit_loads.end()) it = unit_loads.emplace(t.id, tendon_loads(model, t, 1.0)).first;
  for (const Load& u : it->second) {
    if (only && !only->count(u.bar)) continue;
    Load l = u;
    l.F *= fraction;
    l.M *= fraction;
    l.F2 *= fraction;
    l.M2 *= fraction;
    out.push_back(l);
  }
}

void PhasingSession::Impl::step(const Phase& ph, std::size_t run_i, std::size_t run_n) {
  hooks_check(opt.hooks);
  hooks_report(opt.hooks, static_cast<double>(run_i) / static_cast<double>(run_n),
               "phase " + std::to_string(run_i + 1) + "/" + std::to_string(run_n));
  auto since = [](std::chrono::steady_clock::time_point a) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - a).count();
  };
  PhasedFactorCache& cache = *cache_;
  const std::unordered_set<Id> cast_now(ph.activate_bars.begin(), ph.activate_bars.end());
  CaseResult incremental =
      zero_result(model, opt.stations_per_bar, opt.stations_abs, plan_supports);
  // the time-dependent effects of the interval, decided by the materials
  // and the live tendons (time_dependent), not by the plan
  const bool do_creep = phases_run > 0 && ph.time > prev_time && time_dependent(ph);
  double creep_t1 = prev_time, creep_t2 = ph.time;
  cur_time = ph.time;
  so_now = ph.regime == Phase::Regime::SecondOrder;
  for (const auto& [n, v] : ph.imperfection) {
    require(model.node_index(n) >= 0, "phasing: imperfection at unknown node " + std::to_string(n));
    require(v.allFinite(), "phasing: non-finite imperfection at node " + std::to_string(n));
    auto it = imperf.try_emplace(n, Vec6::Zero()).first;
    it->second += v;
  }

  // 0. relaxation of the strands held by their bed over the interval: the
  //    bed takes it, the structure sees nothing (t_eb of a heat cure added
  //    to the clock, EN 1992-1-1 §10.3.2)
  if (phases_run > 0 && ph.time > prev_time)
    for (const Tendon& t : model.tendons) {
      TendonState& ts = tstate[t.id];
      if (!ts.on_bed || t.losses.Ap <= 0 || t.jack_force <= 0) continue;
      if (t.losses.fpk <= 0 || t.losses.rho1000 <= 0 || ts.sigma_pi <= 0) continue;
      double dsig = 0;
      if (t.losses.relaxation_bpel) {
        const double sigma_now = ts.frac_eff * t.jack_force / t.losses.Ap;
        // heat cure (BPEL 99 annexe 6 §4.1): the loss is that of t + t_eb
        // from the tensioning — the first interval runs t_eb longer
        const double teb = ts.relax_cum == 0 ? t.losses.relax_teb_days : 0.0;
        dsig = bpel_relaxation_step(ph.time - prev_time + teb, sigma_now, ts.relax_cum,
                                    t.losses.fpk, t.losses.rho1000);
        ts.relax_cum += dsig;
      } else {
        // cumulated r(t + t_eb) since the tensioning (none at t = 0)
        const double teb = t.losses.relax_teb_days;
        auto r = [&](double t_abs) {
          const double dt = t_abs - ts.t_tension;
          return dt > 0 ? ec2_relaxation_ratio((dt + teb) * 24.0, ts.sigma_pi, t.losses.fpk,
                                               t.losses.relaxation_class, t.losses.rho1000)
                        : 0.0;
        };
        dsig = ts.sigma_pi * std::max(0.0, r(ph.time) - r(prev_time));
      }
      ts.frac_eff -= dsig * t.losses.Ap / t.jack_force;
    }

  // 1. activations (stress-free) — supports first, as a deck orders them
  for (const Support& s : ph.add_supports) add_support(s);
  std::vector<Id> cable_activations;  // catenary cables: their own coupled increments (below)
  for (Id b : ph.activate_bars) {
    if (model.bar(b).kind == BarKind::Cable) {
      cable_activations.push_back(b);
      continue;
    }
    require(!active_bars.count(b), "phasing: bar " + std::to_string(b) + " already active");
    if (ph.activate_rel.count(b)) {
      // Prolongation of the deformed structure: a node not yet in it takes
      // the rigid continuation of the bar's other end, so do the stations.
      // (an implicit support is the out-of-plane restraint of EVERY node
      // of a PLANE / GRILL model: it does not put the node in the structure)
      // — the degree of the node in the active structure, kept by
      // activate_bar / deactivate_bar / add_support (O(1) a query)
      auto in_structure = [&](Id node) {
        const auto it = struct_deg.find(node);
        return it != struct_deg.end() && it->second > 0;
      };
      const Bar& bar = model.bar(b);
      bool in_i = in_structure(bar.ni), in_j = in_structure(bar.nj);
      if (in_i != in_j) {
        Id from = in_i ? bar.ni : bar.nj, to = in_i ? bar.nj : bar.ni;
        // linear in the state: applied to master and to each origin ledger
        for_states([&](CaseResult& cs) {
          const Vec6 ua = cs.node_disp[static_cast<std::size_t>(model.node_index(from))];
          const Vec3 xa = model.node(from).x;
          const Vec3 tr = ua.head<3>(), th = ua.tail<3>();
          auto carried = [&](const Vec3& p) { return Vec3(tr + th.cross(p - xa)); };
          Vec6& ub = cs.node_disp[static_cast<std::size_t>(model.node_index(to))];
          ub.head<3>() = carried(model.node(to).x);
          ub.tail<3>() = th;
          const Vec3 xi = model.node(bar.ni).x, xj = model.node(bar.nj).x;
          const Mat3 R = bar_rotation(xi, xj, bar.roll);
          const double L = (xj - xi).norm();
          for (Station& stn : cs.bars[static_cast<std::size_t>(model.bar_index(b))].stations) {
            stn.defl += R * carried(xi + (L > 0 ? stn.x / L : 0.0) * (xj - xi));
            stn.rot += R * th;
          }
        });
        activate_bar(b);
        auto age = ph.activate_age.find(b);
        cast_time[b] = ph.time - (age != ph.activate_age.end() ? age->second : 0.0);
        init_parts(b, cast_time[b]);
        continue;
      }
    }
    if (opt.born_deformed) {
      // A bar cast between nodes that have already moved is born
      // stress-free and STRAIGHT on their displaced positions: its stations
      // start from the chord of the end translations (rotations: the
      // chord's, twist interpolated), so the printed station displacements
      // are totals, continuous with the nodes in translation, with a kink
      // at a node the structure had already rotated (ST1, Fontenoy bar 4
      // at node 301: Dz equal, Ry -0.00103 against the node's +0.00046).
      const Bar& bar = model.bar(b);
      const Vec3 xi = model.node(bar.ni).x, xj = model.node(bar.nj).x;
      const Mat3 R = bar_rotation(xi, xj, bar.roll);
      const double L = (xj - xi).norm();
      for_states([&](CaseResult& cs) {
        const Vec6 gi = cs.node_disp[static_cast<std::size_t>(model.node_index(bar.ni))];
        const Vec6 gj = cs.node_disp[static_cast<std::size_t>(model.node_index(bar.nj))];
        if (L > 0 && (gi.squaredNorm() > 0 || gj.squaredNorm() > 0)) {
          const Vec3 ui = R * gi.head<3>(), uj = R * gj.head<3>();
          const double rxi = (R * gi.tail<3>())(0), rxj = (R * gj.tail<3>())(0);
          // theta_z = v', theta_y = -w' of the chord
          const Vec3 chord_rot(0.0, -(uj(2) - ui(2)) / L, (uj(1) - ui(1)) / L);
          for (Station& stn : cs.bars[static_cast<std::size_t>(model.bar_index(b))].stations) {
            const double t = std::clamp(stn.x / L, 0.0, 1.0);
            stn.defl += ui + t * (uj - ui);
            stn.rot += chord_rot + Vec3(rxi + t * (rxj - rxi), 0.0, 0.0);
          }
        }
      });
    }
    activate_bar(b);
    // The concrete age at activation is independent of the activation
    // date (ST1 AGE, p.187): a bar activated at age a was cast a days ago.
    auto age = ph.activate_age.find(b);
    cast_time[b] = ph.time - (age != ph.activate_age.end() ? age->second : 0.0);
    init_parts(b, cast_time[b]);
  }

  for (Id l : ph.activate_links) {
    bool known = false;
    for (const ElasticLink& x : model.links) known = known || x.id == l;
    require(known, "phasing: unknown elastic link " + std::to_string(l));
    require(active_links.insert(l).second, "phasing: elastic link " + std::to_string(l) + " already active");
  }
  for (Id l : ph.activate_rigid_links) {
    bool known = false;
    for (const RigidLink& x : model.rigid_links) known = known || x.id == l;
    require(known, "phasing: unknown rigid link " + std::to_string(l));
    require(active_rigid.insert(l).second, "phasing: rigid link " + std::to_string(l) + " already active");
  }

  // The modulus of an active bar for this phase's increments: E at its age
  // (PhasingOptions::modulus_at), else its material's.
  // (the phase date and the casting dates are fixed within a phase: each
  // bar's modulus is evaluated once)
  std::unordered_map<Id, double> e_now_cache;
  auto E_now = [&](const Bar& b) {
    const Material& mat = model.material(b.material);
    if (!opt.modulus_at || !aging(mat)) return mat.E;
    auto it = e_now_cache.find(b.id);
    if (it != e_now_cache.end()) return it->second;
    const double age = cast_time.count(b.id) ? ph.time - cast_time[b.id] : 28.0;
    const double E = opt.modulus_at(b.id, std::max(0.5, age));
    e_now_cache.emplace(b.id, E);
    return E;
  };

  // The section at abscissa x of a bar: a CARA VAR bar's own law at x / L
  // (Model::tapered props), else its constant section.
  auto props_at = [&](const Bar& b, double x) {
    const TaperedBar* tb = model.tapered_of(b.id);
    if (tb && tb->props) {
      const double L = model.bar_chord_length(b);
      return tb->props(L > 0 ? std::clamp(x / L, 0.0, 1.0) : 0.0);
    }
    return model.section(b.section).props;
  };

  // Working model of the CURRENT active substructure.
  auto build_work = [&]() {
    Model w = model_lite;  // the model without its cases, tendons and bars
    Id next_mat = 0;
    for (const Material& mm : model.materials) next_mat = std::max(next_mat, mm.id);
    for (const Bar& b : model.bars) {
      if (!active_bars.count(b.id)) continue;
      w.bars.push_back(b);
      const Material& mat = model.material(b.material);
      const double E = E_now(b);
      if (E != mat.E) {  // its own material at the modulus of its age
        Material aged = mat;
        aged.id = ++next_mat;
        aged.E = E;
        w.materials.push_back(aged);
        w.bars.back().material = aged.id;
      }
    }
    w.supports = active_supports;
    w.links.clear();
    for (const ElasticLink& l : model.links)
      if (active_links.count(l.id)) w.links.push_back(l);
    w.rigid_links.clear();
    for (const RigidLink& l : model.rigid_links)
      if (active_rigid.count(l.id)) w.rigid_links.push_back(l);
    for (const auto& [bid, ps] : parts)
      if (active_bars.count(bid)) lower_composite_bar(w, bid, composite_state(model, bid, part_E(bid)));
    w.composite.clear();
    // the facets of a creep step at their effective modulus (own materials)
    for (const auto& [mid, E] : shell_E) {
      Material m2 = model.material(mid);
      m2.id = ++next_mat;
      m2.E = E;
      w.materials.push_back(m2);
      for (ShellTri& sh : w.shells)
        if (sh.material == mid) sh.material = m2.id;
    }
    require(!w.bars.empty() || !w.shells.empty(), "phasing: no active bars");
    return w;
  };

  // The work structure assembled and factorised once per state of the
  // substructure (active bars, their moduli, the supports): an increment
  // then costs its load lowering, one back-substitution and the recovery.
  // Unilateral supports keep the general path (an active set per solve).
  struct Prepared {
    std::vector<std::pair<Id, double>> sig;
    std::vector<Support> sups;
    std::set<Id> links, rigid;
    std::vector<double> psig;
    Model w;
    std::unique_ptr<Assembled> asm_;
    SpMat K;  // extended stiffness (kept when rigid links reduce the solve)
    FactorCache* fac = nullptr;
    bool general = false;  // unilateral supports: solve_static_cached
  };
  std::unique_ptr<Prepared> prep;
  auto prepared = [&]() -> Prepared& {
    std::vector<std::pair<Id, double>> sig;
    sig.reserve(active_bars.size());
    for (const Bar& b : model.bars)
      if (active_bars.count(b.id)) sig.emplace_back(b.id, E_now(b));
    std::vector<double> psig;  // the parts' moduli (0 = not cast; E_eff in a creep step)
    for (const auto& [bid, ps] : parts)
      if (active_bars.count(bid))
        for (double e : part_E(bid)) psig.push_back(e);
    for (const auto& [mid, E] : shell_E) {
      psig.push_back(static_cast<double>(mid));
      psig.push_back(E);
    }
    if (prep && prep->sig == sig && prep->sups.size() == active_supports.size() &&
        prep->links == active_links && prep->rigid == active_rigid && prep->psig == psig) {
      bool same = true;
      for (std::size_t i = 0; i < active_supports.size() && same; ++i)
        same = equal_support(prep->sups[i], active_supports[i]);
      if (same) return *prep;
    }
    auto P = std::make_unique<Prepared>();
    P->sig = std::move(sig);
    P->sups = active_supports;
    P->links = active_links;
    P->rigid = active_rigid;
    P->psig = psig;
    P->w = build_work();
    P->w.validate();
    for (const Support& sup : P->w.supports)
      for (int c = 0; c < 6; ++c)
        if (sup.oneside[static_cast<std::size_t>(c)] != 0) P->general = true;
    for (const Support& sup : P->w.supports)
      if (sup.decol_axis >= 0) P->general = true;
    PhasedFactorCache::Key key;
    key.first.assign(active_bars.begin(), active_bars.end());  // std::set: sorted
    for (const Support& sup : active_supports) key.second.push_back(sup.node);
    std::sort(key.second.begin(), key.second.end());
    P->fac = &cache.slot(key);
    if (!P->general) {
      P->asm_ = std::make_unique<Assembled>(assemble_elements(P->w, true));
      SpMat K = assemble_stiffness(P->w, *P->asm_);
      P->fac->set_reduction(P->asm_->red.active ? P->asm_->red.T : SpMat());
      P->fac->factorize(K);
      if (P->asm_->red.active) P->K = std::move(K);
    }
    prep = std::move(P);
    return *prep;
  };
  // one increment on the prepared structure: its lowered loads and result
  auto solve_prepared = [&](Prepared& P, const std::vector<Load>& loads, const StaticOptions& so,
                            CaseLoads& cl) {
    auto t0 = std::chrono::steady_clock::now();
    LoadCase lc{1, "incr", merge_bar_loads(loads)};
    ph_timer.t_merge += since(t0);
    t0 = std::chrono::steady_clock::now();
    cl = lower_case(P.w, *P.asm_, lc, opt.threads);
    fold_prescribed(P.w, *P.asm_, cl);
    if (P.asm_->red.active) fold_link_offsets(P.w, *P.asm_, P.K, cl);
    ph_timer.t_lower += since(t0);
    t0 = std::chrono::steady_clock::now();
    const VecX u = P.fac->solve(cl.rhs);
    ph_timer.t_back += since(t0);
    t0 = std::chrono::steady_clock::now();
    CaseResult res = recover_case(P.w, *P.asm_, cl, u, so, 1, "incr");
    ph_timer.t_recover += since(t0);
    return res;
  };

  // The coupled cable + frame increment: the frame linear about
  // the accumulated state, every active cable its exact catenary at the
  // true anchor positions x = X + u_acc + du,
  //     K du = dF + sum_c [f_c(x + du; L0') - f_c(x; L0)]   (f = force ON the nodes),
  // solved by modified Newton on K + K_t (K_t the cables' tangent at du = 0,
  // factorised once): (K + K_t) du_k+1 = dF + df(du_k) + K_t du_k, until the
  // true residual is at round-off. The converged cable force change is then a
  // nodal load of the increment, so the recovery (bar forces, reactions at a
  // supported anchor) is the linear one. A cable activated by this increment
  // has no baseline, a ruptured one no final force.
  auto solve_cable_increment = [&](Prepared& P, const std::vector<Load>& loads, const StaticOptions& so,
                                    CaseLoads& cf) {
    const DofMap& dm = P.asm_->dofs;
    struct Cab {
      Id bid, n1, n2;
      Vec3 x1, x2;
      double L0b, L0a, EA, w;
      bool before, after;
      std::array<int, 3> d1, d2;
      Vec3 fb1 = Vec3::Zero(), fb2 = Vec3::Zero();
      Mat3 k3 = Mat3::Zero();
    };
    std::vector<Cab> cs;
    for (const auto& [bid, st] : cab) {
      const bool before = st.active, after = (st.active || cab_new.count(bid)) && !cab_off.count(bid);
      if (!before && !after) continue;
      const Bar& b = model.bar(bid);
      Cab c{bid, b.ni, b.nj, node_now(b.ni), node_now(b.nj), st.L0, st.L0, st.EA, st.w, before, after, {}, {}};
      const auto it = cab_next.find(bid);
      if (it != cab_next.end()) c.L0a = it->second;
      for (int k = 0; k < 3; ++k) {
        c.d1[static_cast<std::size_t>(k)] = dm.free_of(model.node_index(b.ni), k);
        c.d2[static_cast<std::size_t>(k)] = dm.free_of(model.node_index(b.nj), k);
      }
      if (before) {
        const CableEval e = eval_cable(c.x1, c.x2, c.L0b, c.EA, c.w);
        c.fb1 = e.f1;
        c.fb2 = e.f2;
      }
      if (after) c.k3 = eval_cable(c.x1, c.x2, c.L0a, c.EA, c.w).k3;
      cs.push_back(c);
    }
    const SpMat K = assemble_stiffness(P.w, *P.asm_);
    const int n = static_cast<int>(K.rows());
    // K + K_t: the pair blocks [k3 -k3; -k3 k3] on the free translations
    std::vector<Triplet> tr;
    tr.reserve(static_cast<std::size_t>(K.nonZeros()) + 36 * cs.size());
    for (int k = 0; k < K.outerSize(); ++k)
      for (SpMat::InnerIterator it(K, k); it; ++it) tr.emplace_back(static_cast<int>(it.row()), static_cast<int>(it.col()), it.value());
    for (const Cab& c : cs)
      for (int a = 0; a < 3; ++a)
        for (int b2 = 0; b2 < 3; ++b2) {
          const double v = c.k3(a, b2);
          if (v == 0) continue;
          const int i1 = c.d1[static_cast<std::size_t>(a)], i2 = c.d2[static_cast<std::size_t>(a)];
          const int j1 = c.d1[static_cast<std::size_t>(b2)], j2 = c.d2[static_cast<std::size_t>(b2)];
          if (i1 >= 0 && j1 >= 0) tr.emplace_back(i1, j1, v);
          if (i2 >= 0 && j2 >= 0) tr.emplace_back(i2, j2, v);
          if (i1 >= 0 && j2 >= 0) tr.emplace_back(i1, j2, -v);
          if (i2 >= 0 && j1 >= 0) tr.emplace_back(i2, j1, -v);
        }
    SpMat M(n, n);
    M.setFromTriplets(tr.begin(), tr.end());
    FactorCache fc;
    fc.factorize(M);
    LoadCase lc{1, "incr", merge_bar_loads(loads)};
    CaseLoads cl = lower_case(P.w, *P.asm_, lc, opt.threads);
    fold_prescribed(P.w, *P.asm_, cl);
    auto trans = [&](const std::array<int, 3>& d, const VecX& u) {
      Vec3 v = Vec3::Zero();
      for (int k = 0; k < 3; ++k)
        if (d[static_cast<std::size_t>(k)] >= 0) v(k) = u(d[static_cast<std::size_t>(k)]);
      return v;
    };
    // df at du: the cable force change on the nodes, per cable (node 1, node 2)
    std::vector<std::pair<Vec3, Vec3>> df(cs.size());
    auto eval_df = [&](const VecX& u) {
      for (std::size_t i = 0; i < cs.size(); ++i) {
        const Cab& c = cs[i];
        Vec3 f1 = -c.fb1, f2 = -c.fb2;
        if (c.after) {
          const CableEval e = eval_cable(c.x1 + trans(c.d1, u), c.x2 + trans(c.d2, u), c.L0a, c.EA, c.w);
          f1 += e.f1;
          f2 += e.f2;
        }
        df[i] = {f1, f2};
      }
    };
    auto scatter = [&](VecX& v, const std::array<int, 3>& d, const Vec3& f) {
      for (int k = 0; k < 3; ++k)
        if (d[static_cast<std::size_t>(k)] >= 0) v(d[static_cast<std::size_t>(k)]) += f(k);
    };
    VecX du = VecX::Zero(n);
    auto newton = [&]() {
    bool ok = false;
    for (int it = 0; it < 100 && !ok; ++it) {
      eval_df(du);
      VecX rhs = cl.rhs, res = cl.rhs - K * du;
      for (std::size_t i = 0; i < cs.size(); ++i) {
        scatter(rhs, cs[i].d1, df[i].first);
        scatter(rhs, cs[i].d2, df[i].second);
        scatter(res, cs[i].d1, df[i].first);
        scatter(res, cs[i].d2, df[i].second);
        const Vec3 rel = trans(cs[i].d2, du) - trans(cs[i].d1, du);  // K_t du on the pair
        if (cs[i].after) {
          scatter(rhs, cs[i].d1, -cs[i].k3 * rel);
          scatter(rhs, cs[i].d2, cs[i].k3 * rel);
        }
      }
      double scale = std::max(cl.rhs.size() ? cl.rhs.lpNorm<Eigen::Infinity>() : 0.0, 1e-30);
      for (const auto& [a, b] : df) scale = std::max({scale, a.lpNorm<Eigen::Infinity>(), b.lpNorm<Eigen::Infinity>()});
      if (it > 0 && (n == 0 || res.lpNorm<Eigen::Infinity>() <= 1e-11 * scale)) {
        ok = true;
        break;
      }
      if (n > 0) du = fc.solve(rhs);
    }
    require(ok, "phasing: the coupled cable + frame increment does not converge");
    };
    newton();
    // force targets: L0 such that the target holds at the converged anchors —
    // a secant on g(L0) = L0 - l0_target(anchors(L0)) per targeted cable
    // (its first step the fixed point), each evaluation one warm-started
    // coupled solve on the same factorisation
    if (!cab_target.empty()) {
      std::vector<std::size_t> ti;
      for (std::size_t i = 0; i < cs.size(); ++i)
        if (cab_target.count(cs[i].bid)) ti.push_back(i);
      std::vector<double> gp(ti.size(), 0), lp(ti.size(), 0);
      bool done = false;
      for (int k = 0; k < 60 && !done; ++k) {
        done = true;
        for (std::size_t q = 0; q < ti.size(); ++q) {
          Cab& c = cs[ti[q]];
          const auto [kind, val] = cab_target.at(c.bid);
          const Vec3 d = (c.x2 + trans(c.d2, du)) - (c.x1 + trans(c.d1, du));
          const double dxh = std::hypot(d.x(), d.y());
          const double l0s = kind == 0 ? catenary_l0_for_tension(dxh, d.z(), c.EA, c.w, val)
                                       : catenary_l0_for_end_tension(dxh, d.z(), c.EA, c.w, val);
          const double g = c.L0a - l0s;
          if (std::abs(g) <= 1e-13 * c.L0a) continue;
          done = false;
          double next = l0s;
          if (k > 0 && g != gp[q]) next = c.L0a - g * (c.L0a - lp[q]) / (g - gp[q]);
          gp[q] = g;
          lp[q] = c.L0a;
          c.L0a = next;
        }
        if (!done) newton();
      }
      require(done, "phasing: the cable force target is not reached at equilibrium");
      for (std::size_t q : ti) cab_next[cs[q].bid] = cs[q].L0a;
    }
    // the converged cable force change as nodal loads of the increment
    std::vector<Load> all = loads;
    for (std::size_t i = 0; i < cs.size(); ++i) {
      all.push_back(Load::nodal(cs[i].n1, df[i].first));
      all.push_back(Load::nodal(cs[i].n2, df[i].second));
    }
    LoadCase lf{1, "incr", merge_bar_loads(all)};
    cf = lower_case(P.w, *P.asm_, lf, opt.threads);
    fold_prescribed(P.w, *P.asm_, cf);
    return recover_case(P.w, *P.asm_, cf, du, so, 1, "incr");
  };

  // The second-order increment (Phase::regime), about the
  // accumulated state: [K + Kg(N1)] du = dF - sum_e [kg_e(N1) - kg_e(N0)] d_e,
  // N0 the accumulated axial force of each bar, N1 = N0 + N(du) (fixed
  // point, linear recovery of du as solve_second_order), d_e the bar's
  // displacement since its activation plus the imperfection. Bar stations
  // are the linear recovery of du (the solve_second_order convention: the
  // P-Delta moments come through the nodal displacements); the reactions
  // gain the geometric terms sum_e kg_e(N1)(d_e + du_e) - kg_e(N0) d_e.
  auto solve_so_increment = [&](Prepared& P, const std::vector<Load>& loads, const StaticOptions& so,
                                CaseLoads& cl) {
    const Model& w = P.w;
    const Assembled& as = *P.asm_;
    const DofMap& dm = as.dofs;
    LoadCase lc{1, "incr", merge_bar_loads(loads)};
    cl = lower_case(w, as, lc, opt.threads);
    require(cl.prescribed.empty(), "phasing: imposed displacements in a second-order phase are not supported");
    const std::size_t ne = as.elements.size();
    const int n = dm.n_free;
    std::vector<double> N0(ne, 0.0), N1(ne, 0.0);
    std::vector<Vec12> d(ne, Vec12::Zero());
    std::vector<std::array<int, 12>> g(ne);
    for (std::size_t e = 0; e < ne; ++e) {
      const Bar& b = w.bars[e];
      N0[e] = -master.bars[static_cast<std::size_t>(model.bar_index(b.id))].end_forces_local(0);
      d[e] = bar_d(model.bar(b.id));
      const int ni = w.node_index(b.ni), nj = w.node_index(b.nj);
      for (int c = 0; c < 6; ++c) {
        g[e][static_cast<std::size_t>(c)] = dm.free_ix[static_cast<std::size_t>(dm.full(ni, c))];
        g[e][static_cast<std::size_t>(6 + c)] = dm.free_ix[static_cast<std::size_t>(dm.full(nj, c))];
      }
    }
    auto gather = [&](std::size_t e, const VecX& u) {
      Vec12 v = Vec12::Zero();
      for (int a = 0; a < 12; ++a)
        if (g[e][static_cast<std::size_t>(a)] >= 0) v(a) = u(g[e][static_cast<std::size_t>(a)]);
      return v;
    };
    const SpMat K = assemble_stiffness(w, as);
    VecX du = n > 0 ? P.fac->solve(cl.rhs) : VecX();
    StaticOptions s2;
    s2.stations_per_bar = 2;
    bool ok = n == 0;
    for (int it = 0; it < 60 && !ok; ++it) {
      const CaseResult tmp = recover_case(w, as, cl, du, s2, 0, "");
      for (std::size_t e = 0; e < ne; ++e) N1[e] = N0[e] - tmp.bars[e].end_forces_local(0);
      SpMat Kt = (K + assemble_geometric(w, as, N1)).eval();
      FactorCache ft;
      ft.factorize(Kt);
      require(ft.inertia().n_neg == 0, "phasing: second-order phase beyond the critical load of the current state");
      VecX rhs = cl.rhs;
      for (std::size_t e = 0; e < ne; ++e) {
        if (N1[e] == N0[e] || d[e].isZero(0)) continue;
        const Vec12 f = (as.elements[e].kg(N1[e]) - as.elements[e].kg(N0[e])) * d[e];
        for (int a = 0; a < 12; ++a)
          if (g[e][static_cast<std::size_t>(a)] >= 0) rhs(g[e][static_cast<std::size_t>(a)]) -= f(a);
      }
      const VecX du1 = ft.solve(rhs);
      const double ch = (du1 - du).lpNorm<Eigen::Infinity>(), un = du1.lpNorm<Eigen::Infinity>();
      du = du1;
      ok = ch <= 1e-12 * std::max(un, 1e-300);
    }
    require(ok, "phasing: the second-order increment does not converge");
    CaseResult res = recover_case(w, as, cl, du, so, 1, "incr");
    // the geometric reaction terms (final N1 of the converged du)
    {
      const CaseResult tmp = recover_case(w, as, cl, du, s2, 0, "");
      for (std::size_t e = 0; e < ne; ++e) N1[e] = N0[e] - tmp.bars[e].end_forces_local(0);
    }
    std::unordered_map<Id, Vec6> fg;
    for (std::size_t e = 0; e < ne; ++e) {
      if (N1[e] == 0.0 && N0[e] == 0.0) continue;
      const Vec12 de = gather(e, du);
      const Vec12 f = as.elements[e].kg(N1[e]) * (d[e] + de) - as.elements[e].kg(N0[e]) * d[e];
      fg.try_emplace(w.bars[e].ni, Vec6::Zero()).first->second += f.head<6>();
      fg.try_emplace(w.bars[e].nj, Vec6::Zero()).first->second += f.tail<6>();
    }
    for (Reaction& re : res.reactions) {
      const auto it = fg.find(re.node);
      if (it != fg.end()) re.R += it->second;
    }
    res.non_superposable = true;
    return res;
  };

  // Solve one increment given loads; accumulate into master + incremental.
  auto solve_increment_raw = [&](const std::vector<Load>& loads) {
    const bool cable_event = !cab_next.empty() || !cab_new.empty() || !cab_off.empty();
    if (loads.empty() && !cable_event) return CaseResult{};
    const auto ts0 = std::chrono::steady_clock::now();
    struct Acc {
      double& acc;
      long& n;
      std::chrono::steady_clock::time_point t;
      decltype(since)& f;
      ~Acc() {
        acc += f(t);
        ++n;
      }
    } acc_guard{ph_timer.solve, ph_timer.nsolve, ts0, since};
    StaticOptions so;
    so.stations_per_bar = opt.stations_per_bar;
    so.bar_threads = opt.threads;
    so.stations_abs = opt.stations_abs;
    StaticResults r;
    Model w;
    const Model* W = &w;            // the work structure of this increment
    const Assembled* A = nullptr;   // its assembly, when prepared
    CaseLoads cl_prep;
    // Pseudo-plastic soil: iterate on the bed states of the increment.
    std::map<Id, PseBar> trial;
    for (auto& [bid, pb] : pse)
      if (active_bars.count(bid)) {
        const std::size_t n =
            master.bars[static_cast<std::size_t>(model.bar_index(bid))].stations.size();
        pb.u.resize(n, {});
        pb.p.resize(n, {});
        trial[bid] = pb;
      }
    if (so_now) {
      require(trial.empty(), "phasing: nonlinear soil springs in a second-order phase are not supported");
      require(!cab_live(), "phasing: catenary cables in a second-order phase are not supported");
      require(model.shells.empty(), "phasing: facet shells in a second-order phase are not supported");
      Prepared& P = prepared();
      require(!P.general && !P.asm_->red.active,
              "phasing: unilateral supports or rigid links in a second-order phase are not supported");
      r.cases.push_back(solve_so_increment(P, loads, so, cl_prep));
      W = &P.w;
      A = P.asm_.get();
      invalidate_ledgers("second-order increment (phase " + ph.name + ")");
    } else if (trial.empty() && cab_live()) {
      Prepared& P = prepared();
      require(!P.general && !P.asm_->red.active,
              "phasing: catenary cables with unilateral supports or rigid links are not supported");
      r.cases.push_back(solve_cable_increment(P, loads, so, cl_prep));
      W = &P.w;
      A = P.asm_.get();
    } else if (trial.empty()) {
      Prepared& P = prepared();
      if (!P.general) {
        r.cases.push_back(solve_prepared(P, loads, so, cl_prep));
        W = &P.w;
        A = P.asm_.get();
      } else {
        // This substructure's factor slot (see PhasedFactorCache).
        w = P.w;
        w.cases.push_back({1, "incr", merge_bar_loads(loads)});
        r = solve_static_cached(w, so, *P.fac);
      }
    } else {
      require(!cab_live(), "phasing: catenary cables with nonlinear soil springs are not supported");
      w = build_work();
      w.cases.push_back({1, "incr", merge_bar_loads(loads)});
    }
    for (int it = 0; !trial.empty(); ++it) {
      w = build_work();
      std::vector<Load> lw = loads;
      for (const auto& [bid, tb] : trial) {
        const PseBar& pb = pse.at(bid);
        std::vector<BefZone> zs;
        for (const SoilZone& z : tb.soil->zones) {
          BefZone bz{z.xl_end, z.kfx, z.kfy, z.kfz, z.kmx};
          if (tb.st[kBedFx] != BedState::Elastic) bz.kx = 0;
          if (tb.st[kBedFy] != BedState::Elastic) bz.ky = 0;
          if (tb.st[kBedFz] != BedState::Elastic) bz.kz = 0;
          zs.push_back(bz);
        }
        install_bef(w, bid, zs);
        // The soil force that changes with the state at the start of the
        // increment, released as a line load (p_prev - p_new on the bar).
        const auto& sts = master.bars[static_cast<std::size_t>(model.bar_index(bid))].stations;
        std::vector<Vec3> q(sts.size(), Vec3::Zero());
        bool any = false;
        for (int c : {kBedFx, kBedFy, kBedFz}) {
          if (tb.st[c] == pb.st[c] && (tb.st[c] != BedState::Capped || tb.sign[c] == pb.sign[c]))
            continue;
          for (std::size_t k = 0; k < sts.size(); ++k) {
            const double d =
                pb.p[k][c] - pse_p(tb, c, tb.st[c], tb.sign[c], sts[k].x, pb.u[k][c]);
            q[k](c == kBedFx ? 0 : c == kBedFy ? 1 : 2) += d;
            any = any || d != 0.0;
          }
        }
        if (any)
          for (std::size_t k = 0; k + 1 < sts.size(); ++k)
            if (sts[k + 1].x > sts[k].x)
              lw.push_back(Load::linear(bid, q[k], q[k + 1], LoadFrame::Local, sts[k].x,
                                        sts[k + 1].x));
      }
      w.cases.push_back({1, "incr", lw});
      r = solve_static(w, so);
      // the states the total displacement asks for (ST1: the mean pressure);
      // an active set still moving after 100 passes is an error, not a
      // silently accepted state
      require(it < 100, "phasing: the pseudo-plastic soil states do not converge (EXEC PHASAGE PSE)");
      bool changed = false;
      for (auto& [bid, tb] : trial) {
        const PseBar& pb = pse.at(bid);
        const BarResult* br = record_of(r.cases[0], bid);
        if (!br) continue;
        const SoilPlastic& pl = *tb.soil->plastic;
        const std::size_t n = std::min(br->stations.size(), pb.u.size());
        for (int c : {kBedFx, kBedFy, kBedFz}) {
          double sum = 0, len = 0;
          for (std::size_t k = 0; k + 1 < n; ++k) {
            const double x0 = br->stations[k].x, x1 = br->stations[k + 1].x;
            const double p0 = pse_k(*tb.soil, c, x0) * (pb.u[k][c] + station_u(br->stations[k], c));
            const double p1 =
                pse_k(*tb.soil, c, x1) * (pb.u[k + 1][c] + station_u(br->stations[k + 1], c));
            sum += 0.5 * (x1 - x0) * (p0 + p1);
            len += x1 - x0;
          }
          if (len <= 0) continue;
          const double pm = sum / len, ff = pse_ff(pl, c);
          const double tol = 1e-9 * std::max(std::fabs(pm), 1.0);
          BedState want = BedState::Elastic;
          double sg = 0;
          if (!pl.bilateral() && (pl.pos ? pm < -tol : pm > tol)) {
            want = BedState::Inactive;
          } else if (ff > 0 && std::fabs(pm) > ff) {
            want = BedState::Capped;
            sg = pm >= 0 ? 1.0 : -1.0;
          }
          // One state at a time, through the elastic one: the trial pressure
          // of a lifted (or capped) bar comes from a structure without its
          // soil stiffness and overshoots — jumping straight to the opposite
          // state would make the active set swing between the two.
          if (tb.st[c] != BedState::Elastic && want != BedState::Elastic && want != tb.st[c]) {
            want = BedState::Elastic;
            sg = 0;
          }
          if (want != tb.st[c] || sg != tb.sign[c]) {
            tb.st[c] = want;
            tb.sign[c] = sg;
            changed = true;
          }
        }
      }
      if (!changed) break;
    }
    CaseResult& inc = r.cases[0];
    // commit the soil: displacement since activation and the pressure of the
    // converged states
    for (auto& [bid, tb] : trial) {
      PseBar& pb = pse.at(bid);
      pb.st = tb.st;
      pb.sign = tb.sign;
      BarResult& mb = master.bars[static_cast<std::size_t>(model.bar_index(bid))];
      const BarResult* br = record_of(inc, bid);
      for (std::size_t k = 0; k < pb.u.size(); ++k)
        for (int c = 0; c < kBedComps; ++c) {
          if (br && k < br->stations.size()) pb.u[k][c] += station_u(br->stations[k], c);
          pb.p[k][c] = c == kBedMx ? pb.p[k][c] + (br && k < br->soil_m.size() ? br->soil_m[k] : 0.0)
                                   : pse_p(pb, c, pb.st[c], pb.sign[c], mb.stations[k].x, pb.u[k][c]);
        }
    }
    if (!parts.empty()) composite_increment(inc);
    accumulate(master, model, inc, *W);
    if (led_on) accumulate(led[static_cast<std::size_t>(origin_now)], model, inc, *W);
    // the cable event of this increment is done: commit the cable states
    for (Id b : cab_new) cab[b].active = true;
    for (Id b : cab_off) cab[b].active = false;
    for (const auto& [b, L0] : cab_next) cab[b].L0 = L0;
    cab_new.clear();
    cab_off.clear();
    cab_next.clear();
    cab_target.clear();
    accumulate(incremental, model, inc, *W);
    shell_increment(inc);
    shell_dstar.clear();  // the imposed field belongs to this increment only
    for (const auto& [bid, tb] : trial) {
      (void)tb;
      const PseBar& pb = pse.at(bid);
      BarResult& mb = master.bars[static_cast<std::size_t>(model.bar_index(bid))];
      mb.soil_p.assign(pb.p.size(), Vec3::Zero());
      mb.soil_m.assign(pb.p.size(), 0.0);
      for (std::size_t k = 0; k < pb.p.size(); ++k) {
        mb.soil_p[k] = Vec3(pb.p[k][kBedFx], pb.p[k][kBedFy], pb.p[k][kBedFz]);
        mb.soil_m[k] = pb.p[k][kBedMx];
      }
    }
    // bar nodal force accumulation for future spring-back
    Assembled asm_own;
    CaseLoads cl_own;
    if (!A) {
      asm_own = assemble_elements(w, true);
      cl_own = lower_case(w, asm_own, w.cases[0]);
    }
    const Assembled& asm_ = A ? *A : asm_own;
    const CaseLoads& cl = A ? cl_prep : cl_own;
    for (std::size_t e = 0; e < W->bars.size(); ++e) {
      int ix = model.bar_index(W->bars[e].id);
      int ni = W->node_index(W->bars[e].ni), nj = W->node_index(W->bars[e].nj);
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        ue(c) = inc.node_disp[static_cast<std::size_t>(ni)](c);
        ue(6 + c) = inc.node_disp[static_cast<std::size_t>(nj)](c);
      }
      Vec12 f_eq = asm_.elements[e].T().transpose() * asm_.elements[e].condense_f0(cl.f0[e]);
      bar_nodal_acc[static_cast<std::size_t>(ix)] += asm_.elements[e].ke() * ue - f_eq;
    }
    for (const LinkElem& le : asm_.links) {
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        ue(c) = inc.node_disp[static_cast<std::size_t>(le.ni)](c);
        ue(6 + c) = inc.node_disp[static_cast<std::size_t>(le.nj)](c);
      }
      auto it = link_nodal_acc.find(le.id);
      if (it == link_nodal_acc.end()) it = link_nodal_acc.emplace(le.id, Vec12::Zero()).first;
      it->second += le.K * ue;
    }
    // NOTE: `reaction_acc[node] += R` would DEFAULT-CONSTRUCT the Vec6 on a
    // node absent from the map, and a default-constructed Eigen matrix is
    // UNINITIALISED — the running sum would start from heap garbage. Nodes
    // supported only from a later phase (`Phase::add_supports`) are not in
    // the pre-seeded set, so this really happens: a support added by a phase
    // and removed by a later one sprang back with a garbage force.
    for (const auto& re : inc.reactions) {
      auto it = reaction_acc.find(re.node);
      if (it == reaction_acc.end()) it = reaction_acc.emplace(re.node, Vec6::Zero()).first;
      it->second += re.R;
    }
    return inc;
  };

  // A recorded increment, then the bonded losses of the live shell tendons
  // it strains: the loss loads lose again in every one of them, iterated
  // (each a recorded increment: their creep follows).
  auto solve_increment = [&](const std::vector<Load>& loads) {
    CaseResult inc = solve_increment_raw(loads);
    if (in_shell_loss || model.shell_tendons.empty()) return inc;
    in_shell_loss = true;
    std::vector<Load> ll = shell_losses(inc);
    for (int it = 0; !ll.empty(); ++it) {
      require(it < 200, "phasing: the shell tendon losses do not converge (coupling >= 1)");
      ll = shell_losses(solve_increment_raw(ll));
    }
    in_shell_loss = false;
    return inc;
  };

  // Length-weighted mean, over the ACTIVE part of a tendon between its
  // anchors, of f(bar, its record in `res`, station index, e_y, e_z) sampled
  // at the stations of `res` (trapezoid rule along each bar). Returns false
  // when the tendon has no active anchored station. Mid-bar scheme: the
  // mid station of every active chain bar, weighted by the bar length.
  // the tendon's eccentricities at the stations (fixed geometry), per
  // (tendon, chain bar, station abscissa)
  std::unordered_map<Id, std::vector<std::vector<std::array<double, 3>>>> ecc_cache;
  auto ecc_at = [&](const Tendon& t, std::size_t bi, std::size_t k, double x, double s) {
    auto& per_bar = ecc_cache[t.id];
    if (per_bar.size() != t.bars.size()) per_bar.assign(t.bars.size(), {});
    auto& v = per_bar[bi];
    if (v.size() <= k) v.resize(k + 1, {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
    if (!(v[k][0] == x)) v[k] = {x, tendon_ey(t, s), tendon_ez(t, s)};
    return std::make_pair(v[k][1], v[k][2]);
  };
  auto tendon_mean = [&](const Tendon& t, const CaseResult& res, auto&& f, double& mean) {
    const double Lchain = tendon_chain_length(model, t);
    const double sb = t.s_begin, se = t.s_end < 0 ? Lchain : t.s_end;
    double sum = 0, wsum = 0, s0 = 0;
    for (std::size_t bi = 0; bi < t.bars.size(); ++bi) {
      Id bid = t.bars[bi];
      const Bar& b = model.bar(bid);
      double Lb = model.bar_chord_length(b);
      const BarResult* br = active_bars.count(bid) ? record_of(res, bid) : nullptr;
      if (br && !br->stations.empty() && !along) {
        double smid = s0 + Lb / 2;
        sum += Lb * f(b, *br, br->stations.size() / 2, tendon_ey(t, smid), tendon_ez(t, smid));
        wsum += Lb;
      } else if (br && !br->stations.empty()) {
        const bool flip = bi < t.flips.size() && t.flips[bi];
        const std::size_t n = br->stations.size();
        for (std::size_t k = 0; k < n; ++k) {
          double x = br->stations[k].x;
          double s = s0 + (flip ? Lb - x : x);
          if (s < sb - 1e-9 || s > se + 1e-9) continue;
          double xa = k > 0 ? br->stations[k - 1].x : x;
          double xb = k + 1 < n ? br->stations[k + 1].x : x;
          double w = n == 1 ? Lb : 0.5 * (xb - xa);
          const auto [ey, ez] = ecc_at(t, bi, k, x, s);
          sum += w * f(b, *br, k, ey, ez);
          wsum += w;
        }
      }
      s0 += Lb;
    }
    if (wsum <= 0) return false;
    mean = sum / wsum;
    return true;
  };

  // The loads of a loss field increment dL (force lost at each grid point)
  // on the active bars: the tendon with that field as its authored profile.
  // The equivalent loads of a loss field are LINEAR in its grid values (the
  // positions depend on the geometry only): the loads of a unit loss at
  // each grid node are computed once per tendon and combined — the same
  // loads as `tendon_loads` on the field, without re-sampling the trace at
  // every pass of the losses.
  auto local_loss_loads = [&](const Tendon& t, const std::vector<double>& dL, double sign) {
    std::vector<Load> lf;
    bool any = false;
    for (double d : dL) any = any || std::fabs(d) > 1e-12 * std::max(1.0, t.jack_force);
    if (!any) return lf;
    const std::vector<std::pair<double, double>> prof = grid_profile(t.id, dL);
    auto direct = [&]() {
      Tendon tc = t;
      tc.force_profile = prof;
      std::vector<Load> o;
      for (const Load& l : tendon_loads(model, tc, sign))
        if (active_bars.count(l.bar)) o.push_back(l);
      return o;
    };
    LossBasis& lb = loss_basis[t.id];
    if (!lb.built) {
      lb.built = true;
      lb.ok = true;
      for (std::size_t j = 0; j < prof.size() && lb.ok; ++j) {
        Tendon tc = t;
        tc.force_profile = prof;
        for (std::size_t q = 0; q < prof.size(); ++q) tc.force_profile[q].second = q == j ? 1.0 : 0.0;
        lb.s.push_back(prof[j].first);
        // the unit loss at node j is a hat between its neighbours
        const double lo = j > 0 ? prof[j - 1].first : -std::numeric_limits<double>::infinity();
        const double hi = j + 1 < prof.size() ? prof[j + 1].first : std::numeric_limits<double>::infinity();
        lb.loads.push_back(tendon_loads(model, tc, 1.0, 8, lo, hi));
        lb.ok = lb.loads.back().size() == lb.loads.front().size();
        for (std::size_t q = 0; q < lb.loads.back().size() && lb.ok; ++q)
          lb.ok = lb.loads.back()[q].kind == lb.loads.front()[q].kind &&
                  lb.loads.back()[q].bar == lb.loads.front()[q].bar;
      }
      // a unit loss at a node only loads the spans around it: keep those
      if (lb.ok) {
        lb.nz.resize(lb.loads.size());
        for (std::size_t j = 0; j < lb.loads.size(); ++j)
          for (std::size_t q = 0; q < lb.loads[j].size(); ++q) {
            const Load& l = lb.loads[j][q];
            if (l.F != Vec3::Zero() || l.M != Vec3::Zero() || l.F2 != Vec3::Zero() || l.M2 != Vec3::Zero())
              lb.nz[j].push_back({q, l.F, l.M, l.F2, l.M2});
          }
      }
    }
    if (!lb.ok || lb.s.size() != prof.size()) return direct();
    if (lb.zero.empty()) {
      lb.zero = lb.loads.front();
      for (Load& l : lb.zero) {
        l.F = l.M = l.F2 = l.M2 = Vec3::Zero();
        lb.bix.push_back(model.bar_index(l.bar));
      }
    }
    std::vector<Load> sum = lb.zero;
    for (std::size_t j = 0; j < prof.size(); ++j) {
      if (std::fabs(prof[j].first - lb.s[j]) > 1e-9) return direct();
      const double c = sign * prof[j].second;
      if (c == 0) continue;
      for (const LossBasis::Nz& z : lb.nz[j]) {
        Load& o = sum[z.q];
        o.F += c * z.F;
        o.M += c * z.M;
        o.F2 += c * z.F2;
        o.M2 += c * z.M2;
      }
    }
    lf.reserve(sum.size());
    for (std::size_t q = 0; q < sum.size(); ++q)
      if (active_ix[static_cast<std::size_t>(lb.bix[q])]) lf.push_back(sum[q]);
    return lf;
  };
  // the bar records of a result by bar id (active bars only)
  auto bar_records = [&](const CaseResult& res) {
    std::unordered_map<Id, const BarResult*> m;
    for (const BarResult& br : res.bars)
      if (active_bars.count(br.bar) && !br.stations.empty()) m[br.bar] = &br;
    return m;
  };

  // K=1 elastic redistribution on previously-live tendons after an increment:
  // each loses Ep Ap times the mean concrete strain change at tendon level
  // (PhasingOptions::local_losses: the strain at each point of its grid).
  // `scale` (when given) weighs each tendon's loss — the non-simultaneous
  // jacking share of a tendon's own shortening; absent tendons take none.
  // One pass of the bonded-tendon (K = 1) losses of an increment: each live
  // tendon follows the concrete strain at its level; returns the loss loads.
  auto k1_pass = [&](const CaseResult& inc, const std::map<Id, double>* scale) {
    const auto tk0 = std::chrono::steady_clock::now();
    struct AccL {
      double& acc;
      std::chrono::steady_clock::time_point t;
      decltype(since)& f;
      ~AccL() { acc += f(t); }
    } accl{ph_timer.loss, tk0, since};
    std::vector<Load> out;
    if (inc.node_disp.empty()) return out;
    std::unordered_map<Id, const BarResult*> recs;  // built at the first local tendon
    bool recs_built = false;
    std::vector<const BarResult*> recv;  // the same, by dense bar index
    for (const Tendon& t : model.tendons) {
      TendonState& st = tstate[t.id];
      if (!st.live || t.losses.Ap <= 0 || t.jack_force <= 0) continue;
      double w = 1.0;
      if (scale) {
        auto it = scale->find(t.id);
        if (it == scale->end()) continue;
        w = it->second;
      }
      double eps = 0;
      auto strain = [&](const Bar& b, const BarResult& br, std::size_t k, double ey, double ez) {
        const Station& stn = br.stations[k];
        const SectionProps sp = props_at(b, stn.x);
        double dsig = stn.N / sp.A;
        if (sp.Iy > 0) dsig += stn.My * ez / sp.Iy;
        if (sp.Iz > 0) dsig -= stn.Mz * ey / sp.Iz;
        return dsig / E_now(b);
      };
      if (local) {
        TendonGrid& g = tgrid.at(t.id);
        if (!recs_built) {
          recs = bar_records(inc);
          recv.assign(model.bars.size(), nullptr);
          for (const auto& [bid, br] : recs) recv[static_cast<std::size_t>(model.bar_index(bid))] = br;
          recs_built = true;
        }
        if (!g.cached) {
          g.cached = true;
          for (std::size_t k = 0; k < g.s.size(); ++k) {
            const Bar& b = model.bar(g.bar[k]);
            const int ix = model.bar_index(b.id);
            const BarResult& mb = master.bars[static_cast<std::size_t>(ix)];
            const double x = g.st[k] < mb.stations.size() ? mb.stations[g.st[k]].x : 0.0;
            g.bp.push_back(&b);
            g.bix.push_back(ix);
            g.x.push_back(x);
            g.props.push_back(props_at(b, x));
          }
        }
        // strain() point by point, its lookups hoisted: the section of the
        // grid point (the same station abscissa), E per bar run
        std::vector<double> dL(g.s.size(), 0.0);
        const Bar* eb = nullptr;
        double E = 0;
        for (std::size_t k = 0; k < g.s.size(); ++k) {
          if (g.shape[k] <= 1e-12 * t.jack_force) continue;  // not bonded there
          const BarResult* r = recv[static_cast<std::size_t>(g.bix[k])];
          if (!r || g.st[k] >= r->stations.size()) continue;
          const Station& stn = r->stations[g.st[k]];
          double sv;
          if (stn.x == g.x[k]) {
            if (g.bp[k] != eb) {
              eb = g.bp[k];
              E = E_now(*eb);
            }
            const SectionProps& sp = g.props[k];
            double dsig = stn.N / sp.A;
            if (sp.Iy > 0) dsig += stn.My * g.ez[k] / sp.Iy;
            if (sp.Iz > 0) dsig -= stn.Mz * g.ey[k] / sp.Iz;
            sv = dsig / E;
          } else {
            sv = strain(*g.bp[k], *r, g.st[k], g.ey[k], g.ez[k]);
          }
          dL[k] = -w * t.losses.Ep * t.losses.Ap * sv;
        }
        double big = 0;
        for (double v : dL) big = std::max(big, std::fabs(v));
        if (big < 1e-9 * std::max(1.0, t.jack_force)) continue;
        smear(t.id, dL);
        st.L.resize(g.s.size(), 0.0);
        for (std::size_t k = 0; k < dL.size(); ++k) st.L[k] += dL[k];
        {
          const std::vector<Load> ll = local_loss_loads(t, dL, -1.0);
          out.insert(out.end(), ll.begin(), ll.end());
        }
        continue;
      }
      if (!tendon_mean(t, inc, strain, eps)) continue;
      double dp = -w * t.losses.Ep * t.losses.Ap * eps;  // compression shortens the tendon
      if (std::fabs(dp) < 1e-9 * std::max(1.0, t.jack_force)) continue;
      double dfrac = -dp / t.jack_force;
      st.frac_eff += dfrac;
      tendon_loads_scaled(t, dfrac, out, &active_bars);
    }
    return out;
  };
  // The exact bonded-tendon losses of an increment: the losses load the
  // structure, whose response loses again in every tendon (coupling
  // included), iterated to convergence (the factor 1 / (1 + sum rho) of the
  // section formulas, without their isolated-section hypothesis);
  // PhasingOptions::exact_losses false keeps the single pass.
  // The structure is linear between two recorded states (no pseudo-plastic
  // soil in play): the passes read only the internal forces of their loss
  // increments, solved on the cached factor without recording, and the
  // state takes the summed losses once — the same total, one recovery.
  auto pse_live = [&]() {
    for (const auto& kv : pse)
      if (active_bars.count(kv.first)) return true;
    return false;
  };
  auto solve_forces = [&](const std::vector<Load>& loads) {
    const auto ts0 = std::chrono::steady_clock::now();
    struct Acc {
      double& acc;
      long& n;
      std::chrono::steady_clock::time_point t;
      decltype(since)& f;
      ~Acc() {
        acc += f(t);
        ++n;
      }
    } acc_guard{ph_timer.solve, ph_timer.nsolve, ts0, since};
    StaticOptions so;
    so.stations_per_bar = opt.stations_per_bar;
    so.bar_threads = opt.threads;
    so.stations_abs = opt.stations_abs;
    so.forces_only = true;
    Prepared& P = prepared();
    if (!P.general) {
      CaseLoads cl;
      return solve_prepared(P, loads, so, cl);
    }
    Model w = P.w;
    w.cases.push_back({1, "incr", merge_bar_loads(loads)});
    StaticResults r = solve_static_cached(w, so, *P.fac);
    return std::move(r.cases[0]);
  };
  // the elastic (K = 1) losses an increment causes: instantaneous (S(Po))
  // unless the increment is deferred (origin 2: the re-losses of a deferred
  // loss, the shortening of the tendons by the structural creep)
  auto k1_losses = [&](const CaseResult& inc, const std::map<Id, double>* scale = nullptr, int origin = 1) {
    const int saved = origin_now;
    origin_now = origin;
    struct Restore {
      int& o;
      int v;
      ~Restore() { o = v; }
    } restore{origin_now, saved};
    std::vector<Load> ll = k1_pass(inc, scale);
    ++ph_timer.nk1;
    // (with live cables the inner passes must be the coupled increments too)
    if (opt.exact_losses && !pse_live() && !cab_live() && !so_now) {
      std::vector<Load> total;
      for (int it = 0; !ll.empty(); ++it) {
        ++ph_timer.nk1it;
        require(it < 200, "phasing: the tendon losses do not converge (coupling >= 1)");
        total.insert(total.end(), ll.begin(), ll.end());
        ll = k1_pass(solve_forces(ll), scale);
      }
      if (!total.empty()) solve_increment(total);
      return;
    }
    // the same tendons, with the same weights, lose again at every pass
    for (int it = 0; !ll.empty(); ++it) {
      ++ph_timer.nk1it;
      require(it < 200, "phasing: the tendon losses do not converge (coupling >= 1)");
      const CaseResult r = solve_increment(ll);
      if (!opt.exact_losses) break;
      ll = k1_pass(r, scale);
    }
  };

  // 2. creep + shrinkage over [previous phase date, this phase date] —
  //    applied FIRST, on the substructure as it stood entering this phase.
  //    Mid-bar scheme: one stress per bar, a uniform imposed strain. Along
  //    the bar (PhasingOptions::creep_along_bar): the creep strain is
  //    computed at EVERY station from the stress there (a span-long bar of a
  //    continuous deck changes sign along its length) and lowered as the
  //    imposed-strain field it is — two concentrated dislocations at the
  //    Gauss points of each segment (exact for the piecewise-linear field),
  //    on every bar (CARA VAR and soil kernels take dislocations exactly).
  // The exact engine follows the hereditary integral through the interval:
  // log-spaced sub-steps (PhasingOptions::creep_substeps_per_decade from
  // 0.1 day after the interval start), the stress — redistribution of the
  // previous sub-steps included — re-evaluated at each, so the result
  // converges on the Volterra solution; ST1's method keeps one step.
  std::vector<double> creep_ends;
  if (do_creep) {
    const int m = opt.creep_engine == CreepEngine::Exact ? opt.creep_substeps_per_decade : 0;
    const double span = creep_t2 - creep_t1, h0 = 0.1;
    if (m > 0 && span > h0)
      for (int k = 0;; ++k) {
        const double e = creep_t1 + h0 * std::pow(10.0, static_cast<double>(k) / m);
        if (e >= creep_t2 - 1e-9 * span) break;
        creep_ends.push_back(e);
      }
    creep_ends.push_back(creep_t2);
  }
  for (std::size_t sub = 0; sub < creep_ends.size(); ++sub) {
    struct SubT {
      double& acc;
      std::chrono::steady_clock::time_point t = std::chrono::steady_clock::now();
      ~SubT() { acc += std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }
    } sub_t{ph_timer.creep};
    // superposition: per (material, interval ages), the creep increment of
    // each past event date — shared by the bars cast together
    std::map<std::tuple<const Material*, double, double>, std::unordered_map<double, double>>
        dph_groups;
    double t1 = sub == 0 ? creep_t1 : creep_ends[sub - 1], t2 = creep_ends[sub];
    std::vector<Load> loads;
    // Per-bar interval quantities kept for the tendon deferred-loss pass:
    // per-station creep strains (axial, kappa_y, kappa_z), the free
    // shrinkage increment and the interval creep coefficient phi(a2, a1)
    // used in the eq. 5.46 interaction denominator.
    struct BarDe {
      std::vector<std::array<double, 3>> de;
      double sh = 0, dphi = 0;
    };
    std::unordered_map<Id, BarDe> bar_de;
    std::map<Id, CompPending> comp_free;  // composite bars of this interval
    // the facets: their creep and shrinkage over the step, their modulus
    std::vector<Vec18> sh_dstar;
    for (std::size_t i = 0; i < shell_mem.size(); ++i) {
      if (!shell_aging(i) || t2 <= 0.5) continue;
      const ShellTri& sh = model.shells[i];
      const Material& mat = model.material(sh.material);
      const double a1 = std::max(0.5, t1), a2 = std::max(a1 + 1e-6, t2);
      const double e_ref = law_E_ref(mat, a1);
      shell_step[sh.material] = {a1, e_ref};
      shell_E[sh.material] = 1.0 / (1.0 / mat.E + phi_of(mat, a2, a1) / e_ref);
      Vec18 ds = Vec18::Zero();
      for (const ShellEvent& ev : shell_mem[i]) {
        const double dph = phi_of(mat, a2, ev.t) - (a1 > ev.t ? phi_of(mat, a1, ev.t) : 0.0);
        ds += ev.d * (mat.E * dph / ev.e_ref);
      }
      const double esh = shrinkage_increment(mat, t1, t2, 0.0);
      if (esh != 0) {
        // u = eps_cs (x, y) in the facet plane about its centroid: the
        // uniform in-plane strain, no curvature, no shear, no spin
        const ShellElem& se = shell_els[i];
        Vec3 xc = Vec3::Zero();
        for (int c = 0; c < 3; ++c) xc += model.nodes[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(c)])].x / 3.0;
        const Vec3 n = se.el.normal();
        for (int c = 0; c < 3; ++c) {
          const Vec3 r = model.nodes[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(c)])].x - xc;
          ds.segment<3>(6 * c) += esh * (r - n * n.dot(r));
        }
      }
      if (ds.squaredNorm() == 0) continue;
      if (sh_dstar.empty()) sh_dstar.assign(shell_mem.size(), Vec18::Zero());
      sh_dstar[i] = ds;
      loads.push_back(Load::shell_imposed(sh.id, ds));
    }
    // the shell tendons' relaxation over the step (EC2, equivalent time on
    // the current stress): a loss load with the creep increment
    for (const ShellTendon& t : model.shell_tendons) {
      ShellTendonState& st = sh_t[t.id];
      if (!st.live || t.Ap <= 0 || t.fpk <= 0 || t.rho1000 <= 0) continue;
      const double sig = st.frac_eff * sh_mean.at(t.id) / t.Ap;
      const double dt = std::max(0.0, t2 - std::max(t1, st.t_tension));
      if (sig <= 0 || dt <= 0) continue;
      const double dsig = ec2_relaxation_step(24.0 * dt, sig, st.relax_cum, t.fpk, t.relaxation_class, t.rho1000);
      st.relax_cum += dsig;
      const double dfrac = -dsig * t.Ap / sh_mean.at(t.id);
      st.frac_eff += dfrac;
      if (dfrac != 0) loads.push_back(Load::shell_prestress(t.id, dfrac));
    }
    std::set<Id> bed_bars;  // bars carrying a strand still on its bed
    for (const Tendon& t : model.tendons)
      if (tstate[t.id].on_bed) bed_bars.insert(t.bars.begin(), t.bars.end());
    for (const Bar& b : model.bars) {
      if (!active_bars.count(b.id)) continue;
      if (cast_now.count(b.id)) continue;  // cast this phase: no creep over the preceding interval
      if (model.composite.count(b.id)) {   // its parts creep (keta/composite.h)
        composite_creep(b.id, t1, t2, loads, comp_free);
        continue;
      }
      const Material& mat = model.material(b.material);
      if (!aging(mat)) continue;
      double cast_t = cast_time.count(b.id) ? cast_time[b.id] : t1;
      // not cast yet (ST1 `ACTIVER BARRE AGE` < 0): no creep, no shrinkage
      if (t2 - cast_t <= 0) continue;
      double a1 = std::max(0.5, t1 - cast_t);
      double a2 = std::max(a1 + 1e-6, t2 - cast_t);
      double deps_sh = shrinkage_increment(mat, t1, t2, cast_t);
      // bonded to a strand still held by its bed: the concrete slides on the
      // bed and the strand restrains it by Ep Ap / (Ec Ac) only (about 3 %),
      // so it shrinks freely (as in ST1: Dx = free shrinkage at the
      // real age, the heat-cure shift t_eb coming with the release),
      // unloaded, hence without creep
      if (bed_bars.count(b.id)) {
        const double sh = mat.concrete_law ? mat.concrete_law->age_shift : 0.0;
        const double d = shrinkage_increment(mat, t1, t2, cast_t, sh);
        if (d != 0) loads.push_back(Load::prestrain(b.id, d));
        continue;
      }
      const BarResult& br = master.bars[static_cast<std::size_t>(model.bar_index(b.id))];
      if (br.stations.empty()) continue;
      const std::size_t n = br.stations.size();
      const SectionProps sp_bar = model.section(b.section).props;
      double dphi = phi_of(mat, a2, a1);
      double e_ref = law_E_ref(mat, a1);
      BarDe& bd = bar_de[b.id];
      bd.de.assign(n, {0.0, 0.0, 0.0});
      bd.sh = deps_sh;
      bd.dphi = dphi;
      // this bar's engine: St1 only where it applies (st1_creep_bars), the
      // superposition elsewhere
      CreepEngine engine = opt.creep_engine;
      if (engine == CreepEngine::St1 && opt.st1_creep_bars && !opt.st1_creep_bars->count(b.id))
        engine = CreepEngine::Exact;
      std::vector<StationCreep>& memory = bar_creep[b.id];
      if (engine == CreepEngine::Kelvin || engine == CreepEngine::Exact) memory.resize(along ? n : 1);
      // stations evaluated: every one along the bar, else the mid station
      const std::size_t kb = along ? 0 : n / 2, ke = along ? n : n / 2 + 1;
      // Section fibres (y, z) for the fibre-wise creep options of the law
      // (ConcreteLaw::creep_tension / nonlinear_fcm): the CARA fibres, else
      // the centroid. A fibre's stress is s0 + s1 z - s2 y, its creep strain
      // d0 + d1 z - d2 y.
      const Section& sec = model.section(b.section);
      std::vector<std::array<double, 2>> fibres;
      for (const StressPoint& p : sec.stress_points) fibres.push_back({p.y, p.z});
      if (fibres.empty()) fibres.push_back({0.0, 0.0});
      const ConcreteLaw* claw = mat.concrete_law.get();
      auto fibre_value = [](const std::array<double, 3>& v, const std::array<double, 2>& f) {
        return v[0] + v[1] * f[1] - v[2] * f[0];
      };
      // St1 engine: the distinct fibres, the law's final value A(t0) and its
      // time function f(tau) = phi(t_ref + tau, t_ref) / A(t_ref), inverted
      // by bisection on log(tau) (f increases from 0 to 1)
      // ST1's fibres (annex B3.1): 2 in a plane section, 4 in space, at the
      // radii of gyration — equal shares of the area that carry N and M
      // exactly (the CARA fibres VY/WY crept too fast)
      std::vector<std::array<double, 2>> st1_fibres;
      {
        const double iy = sp_bar.A > 0 && sp_bar.Iz > 0 ? std::sqrt(sp_bar.Iz / sp_bar.A) : 0.0;
        const double iz = sp_bar.A > 0 && sp_bar.Iy > 0 ? std::sqrt(sp_bar.Iy / sp_bar.A) : 0.0;
        const int lay = opt.st1_fibre_layout;
        for (double y : lay == 1 ? std::vector<double>{0.0} : std::vector<double>{iy, -iy})
          for (double z : lay == 0 ? std::vector<double>{0.0} : std::vector<double>{iz, -iz})
            st1_fibres.push_back({y, z});
      }
      constexpr double kSt1Inf = 1e15;  // days: the final value of phi (BPEL
                                        // converges as 1/sqrt: f(1e9) = 0.9993)
      auto st1_amp = [&](double t0) { return phi_of(mat, t0 + kSt1Inf, t0); };
      auto st1_f = [&](double t_ref, double tau) {
        if (tau <= 0) return 0.0;
        const double a = st1_amp(t_ref);
        return a > 0 ? phi_of(mat, t_ref + tau, t_ref) / a : 0.0;
      };
      auto st1_inv = [&](double t_ref, double r) {
        double lo = -8, hi = 14;  // log10 of the duration
        for (int it = 0; it < 80; ++it) {
          const double mid = 0.5 * (lo + hi);
          (st1_f(t_ref, std::pow(10.0, mid)) < r ? lo : hi) = mid;
        }
        return std::pow(10.0, 0.5 * (lo + hi));
      };
      // tension creep (BPEL family) or none (EC2 family)
      const bool tension_creep = opt.st1_tension_creep_bars && opt.st1_tension_creep_bars->count(b.id);
      std::vector<std::vector<St1Fibre>>& st1_mem = bar_st1[b.id];
      if (engine == CreepEngine::St1) st1_mem.resize(along ? n : 1);
      // the most compressed fibre's stress (<= 0), for the non-linear factor
      auto min_fibre = [&](const std::array<double, 3>& v) {
        double m = 0;
        for (const auto& f : fibres) m = std::min(m, fibre_value(v, f));
        return m;
      };
      // superposition: per event date, its creep increment over the interval
      std::unordered_map<double, double>& dph_group = dph_groups[std::make_tuple(&mat, a1, a2)];
      for (std::size_t k = kb; k < ke; ++k) {
        const Station& stn = br.stations[k];
        // Generalized stresses (N/A, My/Iy, Mz/Iz) and their elastic
        // strains e = s/E: axial strain and the two curvatures.
        const SectionProps sp = props_at(b, stn.x);  // CARA VAR: the section at x
        std::array<double, 3> s{stn.N / sp.A, sp.Iy > 0 ? stn.My / sp.Iy : 0.0,
                                sp.Iz > 0 ? stn.Mz / sp.Iz : 0.0};
        std::array<double, 3>& de = bd.de[k];
        if (engine == CreepEngine::Exact) {
          StationCreep& sc = memory[along ? k : 0];
          ExactEvent ev;
          ev.t = a1;
          ev.e_ref = e_ref;
          bool any = false;
          for (std::size_t c = 0; c < 3; ++c) {
            ev.ds[c] = s[c] - sc.carried[c];
            any = any || ev.ds[c] != 0;
          }
          if (any) {
            if (claw && claw->nonlinear_fcm)
              ev.scale = nonlinear_creep_factor(min_fibre(s), claw->nonlinear_fcm(a1));
            sc.exact.push_back(ev);
          }
          sc.carried = s;
          // the event dates are the bar's, shared by its stations: their
          // creep increments are computed once per group (material,
          // interval) — hashed on the exact date, O(1) an event
          for (const ExactEvent& x : sc.exact) {
            auto [it, fresh] = dph_group.try_emplace(x.t, 0.0);
            if (fresh) it->second = phi_of(mat, a2, x.t) - (a1 > x.t ? phi_of(mat, a1, x.t) : 0.0);
            const double dph = it->second;
            for (std::size_t c = 0; c < 3; ++c) de[c] += x.scale * x.ds[c] * dph / x.e_ref;
          }
        } else if (engine == CreepEngine::Kelvin) {
          StationCreep& sc = memory[along ? k : 0];
          // Stress increment since the last event: above the threshold it
          // gets its own chain fitted at THIS age (Boltzmann superposition);
          // below, it is folded into the latest chain (sub-fit-error drift).
          std::array<double, 3> ds{0, 0, 0};
          bool event = false;
          for (std::size_t c = 0; c < 3; ++c) {
            ds[c] = s[c] - sc.carried[c];
            double scale = std::max(std::fabs(s[c]), std::fabs(sc.carried[c]));
            if (std::fabs(ds[c]) > kKelvinEventRelTol * scale) event = true;
          }
          if (event) {
            CreepEvent ev;
            ev.chain = fit_cache.get(mat, a1, opt.kelvin_units);
            if (claw && claw->nonlinear_fcm)
              ev.scale = nonlinear_creep_factor(min_fibre(s), claw->nonlinear_fcm(a1));
            for (std::size_t c = 0; c < 3; ++c) {
              ev.state[c] = kelvin_state(*ev.chain, a1);
              ev.state[c].sigma = ds[c];  // stress step at the event date
            }
            sc.events.push_back(std::move(ev));
          } else if (!sc.events.empty()) {
            for (std::size_t c = 0; c < 3; ++c) sc.events.back().state[c].sigma += ds[c];
          }
          sc.carried = s;
          for (CreepEvent& ev : sc.events)
            for (std::size_t c = 0; c < 3; ++c)
              de[c] += ev.scale * kelvin_creep(ev.state[c], *ev.chain, a2);
        } else if (engine == CreepEngine::St1) {
          // ST1's incremental method, fibre by fibre (CreepEngine::St1)
          std::vector<St1Fibre>& fs = st1_mem[along ? k : 0];
          fs.resize(st1_fibres.size());
          const double scale = claw && claw->nonlinear_fcm
                                   ? nonlinear_creep_factor(min_fibre(s), claw->nonlinear_fcm(a1))
                                   : 1.0;
          Eigen::Matrix3d A = 1e-12 * Eigen::Matrix3d::Identity();
          Eigen::Vector3d rhs = Eigen::Vector3d::Zero();
          // Plane sections (annex B3.1): the part of the fibre creep strains
          // a plane cannot follow (4 fibres, 3 section strains: SPATIALE)
          // is restrained, a self-equilibrated stress -E (e - plane fit)
          // added to each fibre (nil with 2 fibres)
          std::vector<double> eig(st1_fibres.size(), 0.0);
          if (st1_fibres.size() > 3) {
            Eigen::Matrix3d M = 1e-12 * Eigen::Matrix3d::Identity();
            Eigen::Vector3d v = Eigen::Vector3d::Zero();
            for (std::size_t j = 0; j < st1_fibres.size(); ++j) {
              const Eigen::Vector3d rr(1.0, st1_fibres[j][1], -st1_fibres[j][0]);
              M += rr * rr.transpose();
              v += rr * fs[j].e;
            }
            const Eigen::Vector3d p = M.ldlt().solve(v);
            const double En = E_now(b);
            for (std::size_t j = 0; j < st1_fibres.size(); ++j) {
              const Eigen::Vector3d rr(1.0, st1_fibres[j][1], -st1_fibres[j][0]);
              eig[j] = -En * (fs[j].e - rr.dot(p));
            }
          }
          for (std::size_t j = 0; j < st1_fibres.size(); ++j) {
            St1Fibre& F = fs[j];
            const double sig = fibre_value(s, st1_fibres[j]) + eig[j];
            const double ds = sig - F.sigma;
            double dej = 0;
            if (!(claw && !claw->creep_tension && sig > 0)) {
              F.smax = std::max(F.smax, std::fabs(sig));
              // a continuous recovery towards `target`: the equivalent
              // duration on the recovery curve, else a new one
              auto recover = [&](double target) {
                const double u = target != F.e1 ? (F.e - F.e1) / (target - F.e1) : -1.0;
                if (F.mode == St1Fibre::Mode::Recovery && u > 0 && u < 1) {
                  F.te = st1_inv(F.t_ref, u);
                } else {
                  F.mode = St1Fibre::Mode::Recovery;
                  F.e1 = F.e;
                  F.te = 0;
                }
                F.target = target;
              };
              if (F.reversal) {
                F.reversal = false;
                F.mode = St1Fibre::Mode::Recovery;
                F.e1 = F.e;
                F.te = 0;
                F.target = F.S;
              }
              if (std::fabs(ds) > 1e-6 * F.smax) {
                if (F.t_ref < 0) F.t_ref = a1;
                F.sigma = sig;
                if (!tension_creep && sig > 1e-6 * F.smax) {
                  // EC2 family: no creep in tension. A virgin fibre keeps
                  // no amplitude; a crept one recovers towards
                  // S while S keeps the sign of its creep,
                  // towards 0 once S changes sign; with
                  // 0 < r < 1 (no data) it is held
                  if (F.e == 0) {
                    F.S = 0;
                    F.mode = St1Fibre::Mode::Virgin;
                  } else {
                    F.S += scale * ds * st1_amp(a1) / law_E_ref(mat, a1);
                    const double r = F.S != 0 ? F.e / F.S : 0.0;
                    if (r >= 1) recover(F.S);
                    else if (r <= 0) recover(0.0);
                    else F.mode = St1Fibre::Mode::Virgin;
                  }
                } else {
                  F.S += scale * ds * st1_amp(a1) / law_E_ref(mat, a1);
                  const double r = F.S != 0 ? F.e / F.S : 0.0;
                  if (F.e > 0 && sig < 0 && sig - ds > 0 && F.S < F.e) {
                    // a tension creep under compression: a loading
                    // of S from the equivalent duration f(te) = e/(2e - S);
                    // if S stays positive, a recovery towards
                    // S from the next step. A compression creep going into
                    // tension recovers towards S
                    F.mode = St1Fibre::Mode::Loading;
                    F.te = st1_inv(F.t_ref, F.e / (2 * F.e - F.S));
                    F.reversal = F.S * F.e > 0;
                  } else if (F.mode == St1Fibre::Mode::Recovery &&
                             (F.e - F.e1) * (F.S - F.e1) > 0 &&
                             std::fabs(F.e - F.e1) < std::fabs(F.S - F.e1)) {
                    // a recovery towards S goes on while e lies between
                    // e1 and S, though e may change sign
                    recover(F.S);
                  } else if (F.e == 0 || (r > 0 && r < 1)) {
                    F.mode = St1Fibre::Mode::Loading;
                    F.te = F.e == 0 ? 0.0 : st1_inv(F.t_ref, r);
                  } else {
                    recover(F.S);
                  }
                }
              }
              double e_new = F.e;
              if (F.mode == St1Fibre::Mode::Loading) {
                F.te += a2 - a1;
                e_new = F.S * st1_f(F.t_ref, F.te);
              } else if (F.mode == St1Fibre::Mode::Recovery) {
                F.te += a2 - a1;
                e_new = F.e1 + (F.target - F.e1) * st1_f(F.t_ref, F.te);
              }
              dej = e_new - F.e;
              F.e = e_new;
            }
            const Eigen::Vector3d rr(1.0, st1_fibres[j][1], -st1_fibres[j][0]);
            A += rr * rr.transpose();
            rhs += rr * dej;
          }
          const Eigen::Vector3d d = A.ldlt().solve(rhs);
          for (std::size_t c = 0; c < 3; ++c) de[c] = d[static_cast<int>(c)];
        } else {
          const double f = claw && claw->nonlinear_fcm
                               ? nonlinear_creep_factor(min_fibre(s), claw->nonlinear_fcm(a1))
                               : 1.0;
          for (std::size_t c = 0; c < 3; ++c) de[c] = f * dphi * s[c] / e_ref;
        }
        // Compression-only creep (ST1 TRACTION NON): a fibre in tension takes
        // no increment; the section's increment is refitted (least squares,
        // anchored on the unrestricted one for what the fibres leave free).
        // (The St1 engine applies it fibre by fibre above.)
        if (claw && !claw->creep_tension && engine != CreepEngine::St1) {
          bool any_tension = false;
          for (const auto& fb : fibres) any_tension |= fibre_value(s, fb) > 0;
          if (any_tension) {
            Eigen::Matrix3d A = 1e-9 * Eigen::Matrix3d::Identity();
            Eigen::Vector3d rhs = 1e-9 * Eigen::Vector3d(de[0], de[1], de[2]);
            for (const auto& fb : fibres) {
              Eigen::Vector3d r(1.0, fb[1], -fb[0]);
              const double target = fibre_value(s, fb) > 0 ? 0.0 : fibre_value(de, fb);
              A += r * r.transpose();
              rhs += r * target;
            }
            const Eigen::Vector3d d = A.ldlt().solve(rhs);
            for (std::size_t c = 0; c < 3; ++c) de[c] = d[static_cast<int>(c)];
          }
        }
      }
      // `de` holds the creep strains in the sign of the stresses (fibre
      // strain d0 + d1 ez - d2 ey, as the tendon pass reads them); the
      // imposed curvature of a Load is the opposite measure — kappa = +M/EI
      // bows a beam AGAINST its elastic deflection under M — so the
      // curvatures are lowered negated: creep then adds to the elastic
      // deflection and the redistributions take their true sign (Fontenoy
      // bar 23, creep from day 202 to 250: ST1 x1.30 on the deflection).
      auto lowered = [](std::array<double, 3> v) {
        v[1] = -v[1];
        v[2] = -v[2];
        return v;
      };
      if (!along) {
        // mid-bar scheme: one uniform imposed strain, shrinkage included
        const std::array<double, 3> de = lowered(bd.de[n / 2]);
        Load l = Load::prestrain(b.id, de[0] + deps_sh, de[1], de[2]);
        if (l.strain != 0 || l.kappa_y != 0 || l.kappa_z != 0) loads.push_back(l);
        continue;
      }
      if (deps_sh != 0) loads.push_back(Load::prestrain(b.id, deps_sh));
      auto nonzero = [](const std::array<double, 3>& v) {
        return v[0] != 0 || v[1] != 0 || v[2] != 0;
      };
      if (n == 1) {
        const std::array<double, 3> de = lowered(bd.de[0]);
        if (nonzero(de)) loads.push_back(Load::prestrain(b.id, de[0], de[1], de[2]));
      } else {
        // the linear field of each segment as two dislocations at its Gauss
        // points (weights dx/2): the nodal loads (integrand of degree 3)
        // and the station displacements and rotations (degree <= 2) of the
        // distributed field, exactly
        const double g = 0.5 / std::sqrt(3.0);
        for (std::size_t k = 0; k + 1 < n; ++k) {
          double xa = br.stations[k].x, xb = br.stations[k + 1].x, dx = xb - xa;
          if (dx <= 0) continue;
          const std::array<double, 3> da = lowered(bd.de[k]), db = lowered(bd.de[k + 1]);
          for (double t : {0.5 - g, 0.5 + g}) {
            std::array<double, 3> m{0, 0, 0};
            for (std::size_t c = 0; c < 3; ++c) m[c] = 0.5 * dx * (da[c] + t * (db[c] - da[c]));
            if (nonzero(m)) loads.push_back(Load::prestrain_point(b.id, xa + t * dx, m[0], m[1], m[2]));
          }
        }
      }
    }

    // Deferred prestress losses over the interval — EN 1992-1-1 eq. 5.46
    // applied INCREMENTALLY: every live tendon loses
    //     dsigma_p = (Ep |d_eps_cs| + c_r dsigma_pr + Ep |d_eps_creep|) / den
    //     den      = 1 + (Ep Ap / (E Ac)) (1 + Ac zcp^2 / Ic) (1 + 0.8 dphi)
    // where d_eps_creep is the interval's concrete creep strain AT TENDON
    // LEVEL from the creep engine above (Kelvin or AAEM — the same strains
    // that load the structure), d_eps_cs the interval's free shrinkage, and
    // dsigma_pr the relaxation of the interval: the EC2 §3.3.2 curve
    // evaluated at the time since the last tensioning and DIFFERENCED over
    // [t1, t2] (c_r = 0.8, eq. 5.46), or the BPEL 91 law by the equivalent-
    // time method on the CURRENT tendon stress, which already carries the
    // reduction 0.8 stands for (c_r = 1, ST1). The tendon stress is its jack
    // fraction times P0, or the mean of its friction / draw-in profile along
    // the bar. All section terms are length-weighted means over the
    // tendon's active chain — mid-bar values, or the stations of its
    // anchored length along the bar (biaxial zcp^2/Ic generalization:
    // Ac (ez^2/Iy + ey^2/Iz)); the interaction denominator uses the
    // interval's own phi(a2, a1). The loss is applied to the cable force
    // BEFORE the K = 1 redistribution of the concrete creep increment below;
    // the loss increment itself takes NO K = 1 pass — its elastic + creep
    // feedback on the tendon is exactly what the 5.46 denominator carries
    // (cross-tendon coupling within the interval is second-order and picked
    // up by the next interval's state). A tendon with Ap = 0 (force-only) is
    // skipped; with fpk = 0 the relaxation term alone is disabled. An inert
    // configuration (no creep, no shrinkage, no relaxation) yields an exact
    // zero loss.
    {
      std::vector<Load> loss_loads;
      std::unordered_map<Id, const BarResult*> recs;
      if (local) recs = bar_records(master);
      for (const Tendon& t : model.tendons) {
        TendonState& ts = tstate[t.id];
        if (!ts.live || t.losses.Ap <= 0 || t.jack_force <= 0) continue;
        if (local) {
          // eq. 5.46 at each point of the grid, from the creep, shrinkage
          // and stress there; the relaxation of the local stress
          const TendonGrid& g = tgrid.at(t.id);
          ts.L.resize(g.s.size(), 0.0);
          ts.relax_k.resize(g.s.size(), 0.0);
          std::vector<double> dL(g.s.size(), 0.0);
          const double mf = mean_force[t.id];
          double r1 = 0, r2 = 0;
          const bool ec2_relax = t.losses.fpk > 0 && t.losses.rho1000 > 0 && ts.sigma_pi > 0 &&
                                 !t.losses.relaxation_bpel;
          if (ec2_relax) {
            const double teb = t.losses.relax_teb_days;
            const double h1 = t1 > ts.t_tension ? (t1 - ts.t_tension + teb) * 24.0 : 0.0;
            const double h2 =
                std::max(h1, t2 > ts.t_tension ? (t2 - ts.t_tension + teb) * 24.0 : 0.0);
            auto r = [&](double h) {
              return h > 0 ? ec2_relaxation_ratio(h, ts.sigma_pi, t.losses.fpk,
                                                  t.losses.relaxation_class, t.losses.rho1000)
                           : 0.0;
            };
            r1 = r(h1);
            r2 = r(h2);
          }
          for (std::size_t k = 0; k < g.s.size(); ++k) {
            if (g.shape[k] <= 1e-12 * t.jack_force) continue;
            auto rr = recs.find(g.bar[k]);
            auto bd = bar_de.find(g.bar[k]);
            if (rr == recs.end()) continue;
            const Bar& b = model.bar(g.bar[k]);
            const Material& mat = model.material(b.material);
            // an aging bar absent from the interval (cast this phase, not
            // cast yet, on a bed) carries no loss here; a bar whose
            // material does not age (steel) carries the relaxation alone
            const bool crept = bd != bar_de.end() && g.st[k] < bd->second.de.size();
            if (!crept && aging(mat)) continue;
            const SectionProps sp = props_at(b, rr->second->stations[g.st[k]].x);
            const std::array<double, 3> d =
                crept ? bd->second.de[g.st[k]] : std::array<double, 3>{0.0, 0.0, 0.0};
            const double ec = d[0] + d[1] * g.ez[k] - d[2] * g.ey[k];
            const double esh = crept ? bd->second.sh : 0.0;
            const double dphi_k = crept ? bd->second.dphi : 0.0;
            const double rho = (t.losses.Ep * t.losses.Ap / (mat.E * sp.A)) *
                               (1.0 + (sp.Iy > 0 ? sp.A * g.ez[k] * g.ez[k] / sp.Iy : 0.0) +
                                (sp.Iz > 0 ? sp.A * g.ey[k] * g.ey[k] / sp.Iz : 0.0));
            const double sig = (ts.frac_eff * g.shape[k] - ts.L[k]) / t.losses.Ap;
            double dsig_pr = 0, c_r = 0.8;
            if (opt.exact_losses && ec2_relax && sig > 0) {
              // EC2 relaxation on the current local stress, equivalent time
              // (annex D); t_eb (heat cure) on the first interval
              const double dt = std::max(0.0, t2 - std::max(t1, ts.t_tension)) +
                                (ts.relax_k[k] == 0 ? t.losses.relax_teb_days : 0.0);
              dsig_pr = ec2_relaxation_step(24.0 * dt, sig, ts.relax_k[k], t.losses.fpk,
                                            t.losses.relaxation_class, t.losses.rho1000);
              ts.relax_k[k] += dsig_pr;
              c_r = 1.0;
            } else if (t.losses.fpk > 0 && t.losses.rho1000 > 0 && sig > 0 &&
                t.losses.relaxation_bpel) {
              const double dt = t2 - std::max(t1, ts.t_tension);
              dsig_pr = bpel_relaxation_step(std::max(0.0, dt), sig, ts.relax_k[k],
                                             t.losses.fpk, t.losses.rho1000);
              ts.relax_k[k] += dsig_pr;
              c_r = 1.0;
            } else if (ec2_relax) {
              // the stress the clock started from, at this point
              const double spi = mf > 0 ? ts.sigma_pi * g.shape[k] / mf : ts.sigma_pi;
              dsig_pr = spi * std::max(0.0, r2 - r1);
            }
            const double num = -t.losses.Ep * (ec + esh) + c_r * dsig_pr;
            // exact: the loss's own feedback through the structure (K = 1
            // iterated) and its creep (the exact engine's sub-steps)
            const double den = opt.exact_losses ? 1.0 : 1.0 + rho * (1.0 + 0.8 * dphi_k);
            dL[k] = (num / den) * t.losses.Ap;
          }
          smear(t.id, dL);
          for (std::size_t k = 0; k < dL.size(); ++k) ts.L[k] += dL[k];
          for (const Load& l : local_loss_loads(t, dL, -1.0)) loss_loads.push_back(l);
          continue;
        }
        double eps_c = 0, eps_sh = 0, rho = 0, dphi_m = 0;
        auto creep_at = [&](const Bar& b, const BarResult&, std::size_t k, double ey, double ez) {
          auto it = bar_de.find(b.id);
          if (it == bar_de.end()) return 0.0;
          const auto& d = it->second.de[k];
          return d[0] + d[1] * ez - d[2] * ey;
        };
        auto shrink_at = [&](const Bar& b, const BarResult&, std::size_t, double, double) {
          auto it = bar_de.find(b.id);
          return it == bar_de.end() ? 0.0 : it->second.sh;
        };
        auto dphi_at = [&](const Bar& b, const BarResult&, std::size_t, double, double) {
          auto it = bar_de.find(b.id);
          return it == bar_de.end() ? 0.0 : it->second.dphi;
        };
        auto rho_at = [&](const Bar& b, const BarResult& br, std::size_t k, double ey, double ez) {
          const Material& mat = model.material(b.material);
          const SectionProps sp = props_at(b, br.stations[k].x);
          return (t.losses.Ep * t.losses.Ap / (mat.E * sp.A)) *
                 (1.0 + (sp.Iy > 0 ? sp.A * ez * ez / sp.Iy : 0.0) +
                  (sp.Iz > 0 ? sp.A * ey * ey / sp.Iz : 0.0));
        };
        if (!tendon_mean(t, master, creep_at, eps_c)) continue;
        tendon_mean(t, master, shrink_at, eps_sh);
        tendon_mean(t, master, dphi_at, dphi_m);
        tendon_mean(t, master, rho_at, rho);
        double sigma_now = ts.frac_eff * mean_force[t.id] / t.losses.Ap;
        double dsig_pr = 0, c_r = 0.8;
        if (opt.exact_losses && t.losses.fpk > 0 && t.losses.rho1000 > 0 && sigma_now > 0 &&
            !t.losses.relaxation_bpel && ts.sigma_pi > 0) {
          const double dt = std::max(0.0, t2 - std::max(t1, ts.t_tension)) +
                            (ts.relax_cum == 0 ? t.losses.relax_teb_days : 0.0);
          dsig_pr = ec2_relaxation_step(24.0 * dt, sigma_now, ts.relax_cum, t.losses.fpk,
                                        t.losses.relaxation_class, t.losses.rho1000);
          ts.relax_cum += dsig_pr;
          c_r = 1.0;
        } else if (t.losses.fpk > 0 && t.losses.rho1000 > 0 && sigma_now > 0 && t.losses.relaxation_bpel) {
          double dt = t2 - std::max(t1, ts.t_tension);
          dsig_pr = bpel_relaxation_step(std::max(0.0, dt), sigma_now, ts.relax_cum,
                                         t.losses.fpk, t.losses.rho1000);
          ts.relax_cum += dsig_pr;
          c_r = 1.0;
        } else if (t.losses.fpk > 0 && t.losses.rho1000 > 0 && ts.sigma_pi > 0) {
          // t_eb (heat cure) shifts the clock from the tensioning on
          const double teb = t.losses.relax_teb_days;
          double h1 = t1 > ts.t_tension ? (t1 - ts.t_tension + teb) * 24.0 : 0.0;
          double h2 = std::max(h1, t2 > ts.t_tension ? (t2 - ts.t_tension + teb) * 24.0 : 0.0);
          double r1 = h1 > 0 ? ec2_relaxation_ratio(h1, ts.sigma_pi, t.losses.fpk,
                                                    t.losses.relaxation_class, t.losses.rho1000)
                             : 0.0;
          double r2 = h2 > 0 ? ec2_relaxation_ratio(h2, ts.sigma_pi, t.losses.fpk,
                                                    t.losses.relaxation_class, t.losses.rho1000)
                             : 0.0;
          dsig_pr = ts.sigma_pi * std::max(0.0, r2 - r1);
        }
        // creep + shrinkage shorten the concrete (negative strain at tendon
        // level) -> positive loss; algebraic throughout, so a net extension
        // (creep recovery) is a gain, symmetrically.
        double num = -t.losses.Ep * (eps_c + eps_sh) + c_r * dsig_pr;
        double den = opt.exact_losses ? 1.0 : 1.0 + rho * (1.0 + 0.8 * dphi_m);
        double dp = (num / den) * t.losses.Ap;
        if (std::fabs(dp) < 1e-9 * std::max(1.0, t.jack_force)) continue;
        double dfrac = -dp / t.jack_force;
        ts.frac_eff += dfrac;
        tendon_loads_scaled(t, dfrac, loss_loads, &active_bars);
      }
      if (!loss_loads.empty()) {
        origin_now = 2;
        const CaseResult r = solve_increment(loss_loads);
        if (opt.exact_losses) k1_losses(r, nullptr, 2);  // the losses lose again, iterated
      }
    }

    if (!loads.empty()) {
      comp_pending = std::move(comp_free);  // the parts' strains of THIS increment
      shell_dstar = std::move(sh_dstar);    // the facets' imposed fields of THIS increment
      origin_now = 0;  // structural creep and shrinkage: with the permanent effects
      CaseResult inc = solve_increment(loads);
      k1_losses(inc, nullptr, 2);  // the tendons shortened by creep: a deferred loss
    }
    comp_pending.clear();
    shell_dstar.clear();
    shell_E.clear();
    shell_step.clear();
  }

  // 2b. section parts cast in this phase (stress-free, after the creep of
  //     the interval: a part cast now neither crept nor restrained it)
  for (const auto& [bid, pi] : ph.activate_parts) {
    require(model.composite.count(bid), "phasing: bar " + std::to_string(bid) + " is not composite");
    require(active_bars.count(bid), "phasing: casting a part of inactive bar " + std::to_string(bid));
    std::vector<PartState>& ps = parts.at(bid);
    require(pi >= 0 && static_cast<std::size_t>(pi) < ps.size(),
            "phasing: bar " + std::to_string(bid) + " has no part " + std::to_string(pi));
    PartState& p = ps[static_cast<std::size_t>(pi)];
    require(!p.active, "phasing: part " + std::to_string(pi) + " of bar " + std::to_string(bid) + " already cast");
    p.active = true;
    p.cast = ph.time;
    std::fill(p.res.begin(), p.res.end(), Vec3::Zero());
    p.mem.clear();
  }

  // 3. deactivations (spring-back, one increment)
  {
    std::vector<Load> rel;
    for (Id bid : ph.deactivate_bars) {
      require(active_bars.count(bid), "phasing: deactivating inactive bar " + std::to_string(bid));
      const Bar& b = model.bar(bid);
      require(!b.soil.any(), "phasing: deactivating a soiled bar is not supported");
      std::size_t ix = static_cast<std::size_t>(model.bar_index(bid));
      const Vec12& f = bar_nodal_acc[ix];
      // the removed bar was exerting -f on its nodes; releasing it applies +f
      rel.push_back(Load::nodal(b.ni, f.segment<3>(0), f.segment<3>(3)));
      rel.push_back(Load::nodal(b.nj, f.segment<3>(6), f.segment<3>(9)));
      deactivate_bar(bid);
      // The bar's accumulated state is released with it: wipe its creep
      // memory so a later re-activation starts stress-free AND memory-free.
      bar_creep.erase(bid);
      bar_st1.erase(bid);
      parts.erase(bid);
    }
    for (Id lid : ph.deactivate_links) {
      require(active_links.count(lid), "phasing: deactivating inactive elastic link " + std::to_string(lid));
      const ElasticLink* lk = nullptr;
      for (const ElasticLink& x : model.links)
        if (x.id == lid) lk = &x;
      auto it = link_nodal_acc.find(lid);
      if (it != link_nodal_acc.end()) {
        const Vec12& f = it->second;
        rel.push_back(Load::nodal(lk->ni, f.segment<3>(0), f.segment<3>(3)));
        rel.push_back(Load::nodal(lk->nj, f.segment<3>(6), f.segment<3>(9)));
      }
      active_links.erase(lid);
    }
    for (Id nid : ph.remove_support_nodes) {
      bool found = false;
      for (auto it = active_supports.begin(); it != active_supports.end();)
        if (it->node == nid && !it->implicit) {
          it = active_supports.erase(it);
          --struct_deg[nid];
          found = true;
        } else {
          ++it;
        }
      require(found, "phasing: removing unknown support at node " + std::to_string(nid));
      Vec6 R = reaction_acc.count(nid) ? reaction_acc[nid] : Vec6::Zero();
      rel.push_back(Load::nodal(nid, -R.segment<3>(0), -R.segment<3>(3)));
      reaction_acc[nid].setZero();
    }
    if (!rel.empty()) {
      // zero the removed bars' accumulated state AFTER building the release
      for (Id lid : ph.deactivate_links) {
        link_nodal_acc.erase(lid);
        for_states([&](CaseResult& cs) {
          for (auto& lr : cs.links)
            if (lr.link == lid) {
              lr.deformation.setZero();
              lr.force.setZero();
            }
        });
      }
      // a removed support reacts no more: its recorded reaction is released
      // with the support (the release load carries it into the structure)
      for (Id nid : ph.remove_support_nodes)
        for_states([&](CaseResult& cs) {
          for (Reaction& re : cs.reactions)
            if (re.node == nid) re.R.setZero();
        });
      for (Id bid : ph.deactivate_bars) {
        std::size_t ix = static_cast<std::size_t>(model.bar_index(bid));
        bar_nodal_acc[ix].setZero();
        // wipe master's bar record (it no longer carries force)
        for_states([&](CaseResult& cs) {
          BarResult& bb = cs.bars[ix];
          bb.end_forces_local.setZero();
          for (auto& st : bb.stations) st.N = st.Vy = st.Vz = st.T = st.My = st.Mz = 0;
        });
      }
      // the release is a MIXED-origin locked-in state: one increment that
      // cannot be split, kept in S(G) (the ledgers still sum to master)
      invalidate_ledgers("phase " + ph.name + " : liberation de barres, liaisons ou appuis");
      origin_now = 0;
      CaseResult inc = solve_increment(rel);
      k1_losses(inc);
    }
  }

  // 3c. section changes and phase hinges
  for (const auto& [bid, sid] : ph.modify_sections) {
    require(active_bars.count(bid), "phasing: changing the section of inactive bar " + std::to_string(bid));
    (void)model.section(sid);  // throws on an unknown section
    Bar& b = model.bars[static_cast<std::size_t>(model.bar_index(bid))];
    const Material& mt = model.material(b.material);
    require(!mt.concrete && !mt.concrete_law,
            "phasing: bar " + std::to_string(bid) +
                " creeps: a section change would make its locked-in stress creep as the new section's — cast the added "
                "concrete as a composite part (Phase::activate_parts)");
    require(!model.tapered_of(bid) && !model.composite.count(bid),
            "phasing: section change of a tapered or composite bar " + std::to_string(bid));
    if (b.section == sid) continue;
    b.section = sid;
    // the tendon grids cache the section at their points
    for (auto& [tid, g] : tgrid)
      if (std::find(g.bar.begin(), g.bar.end(), bid) != g.bar.end()) {
        g.cached = false;
        g.bp.clear();
        g.bix.clear();
        g.x.clear();
        g.props.clear();
      }
  }
  for (const Phase::HingeAction& h : ph.place_hinges) {
    require(active_bars.count(h.bar), "phasing: hinge on inactive bar " + std::to_string(h.bar));
    require(h.end == 0 || h.end == 1, "phasing: hinge end 0 (i) or 1 (j)");
    Bar& b = model.bars[static_cast<std::size_t>(model.bar_index(h.bar))];
    const Vec12& f = master.bars[static_cast<std::size_t>(model.bar_index(h.bar))].end_forces_local;
    const double scale = std::max(1.0, f.cwiseAbs().maxCoeff());
    bool any = false;
    for (int c = 0; c < 6; ++c) {
      if (!h.dofs[static_cast<std::size_t>(c)]) continue;
      const std::size_t k = static_cast<std::size_t>(6 * h.end + c);
      require(!b.releases[k], "phasing: bar " + std::to_string(h.bar) + " already released there");
      require(std::abs(f(static_cast<int>(k))) <= 1e-9 * scale,
              "phasing: hinge on bar " + std::to_string(h.bar) +
                  " where a force is locked in (release it by a deactivation, not a hinge)");
      b.releases[k] = true;
      any = true;
    }
    require(any, "phasing: hinge without a released component");
  }
  for (const Phase::HingeAction& h : ph.close_hinges) {
    require(active_bars.count(h.bar), "phasing: clavage of inactive bar " + std::to_string(h.bar));
    require(h.end == 0 || h.end == 1, "phasing: hinge end 0 (i) or 1 (j)");
    Bar& b = model.bars[static_cast<std::size_t>(model.bar_index(h.bar))];
    bool any = false;
    for (int c = 0; c < 6; ++c) {
      if (!h.dofs[static_cast<std::size_t>(c)]) continue;
      const std::size_t k = static_cast<std::size_t>(6 * h.end + c);
      require(b.releases[k], "phasing: clavage of bar " + std::to_string(h.bar) + " where nothing is released");
      b.releases[k] = false;
      any = true;
    }
    require(any, "phasing: clavage without a component");
    if (h.jump[0] != 0 || h.jump[1] != 0 || h.jump[2] != 0 || h.jump[3] != 0) {
      Load d;
      d.kind = Load::Kind::PrestrainPoint;
      d.bar = h.bar;
      d.a = h.end == 0 ? 0.0 : model.bar_chord_length(b);
      d.strain = h.jump[0];
      d.dthx = h.jump[1];
      d.kappa_y = h.jump[2];
      d.kappa_z = h.jump[3];
      origin_now = 0;  // an external imposed action
      CaseResult inc = solve_increment({d});
      k1_losses(inc);
    }
  }

  // 3d. catenary cables: each event one coupled cable + frame increment
  for (Id b : cable_activations) {
    CableSt& c = cable_init(b);
    require(!c.active, "phasing: cable " + std::to_string(b) + " already active");
    const Bar& bar = model.bar(b);
    for (Id n : {bar.ni, bar.nj}) {
      const auto it = struct_deg.find(n);
      require(it != struct_deg.end() && it->second > 0,
              "phasing: cable " + std::to_string(b) + " anchored at node " + std::to_string(n) + " outside the active structure");
    }
    cab_new = {b};
    cab_next[b] = c.L0;
    origin_now = 0;
    CaseResult inc = solve_increment({});
    k1_losses(inc);
  }
  for (const Phase::CableAction& a : ph.cable_actions) {
    require(model.bar_index(a.bar) >= 0 && model.bar(a.bar).kind == BarKind::Cable,
            "phasing: cable action on bar " + std::to_string(a.bar) + " which is not a cable");
    CableSt& c = cable_init(a.bar);
    require(c.active, "phasing: cable action on inactive cable " + std::to_string(a.bar));
    const Bar& bar = model.bar(a.bar);
    const Vec3 d = node_now(bar.nj) - node_now(bar.ni);
    const double dxh = std::hypot(d.x(), d.y());
    using K = Phase::CableAction::Kind;
    switch (a.kind) {
      case K::Length:
        require(a.value > 0, "phasing: cable length > 0");
        cab_next[a.bar] = a.value;
        break;
      // the force targets: first guess at the anchors as they stand, met at
      // the equilibrium of the increment (solve_cable_increment)
      case K::Horizontal:
        cab_next[a.bar] = catenary_l0_for_tension(dxh, d.z(), c.EA, c.w, a.value);
        cab_target[a.bar] = {0, a.value};
        break;
      case K::EndTension:
        cab_next[a.bar] = catenary_l0_for_end_tension(dxh, d.z(), c.EA, c.w, a.value);
        cab_target[a.bar] = {1, a.value};
        break;
      case K::Thermal:
        require(a.value > -1, "phasing: cable thermal strain > -1");
        cab_next[a.bar] = c.L0 * (1 + a.value);
        break;
      case K::Rupture:
        require(a.value >= 1 || a.value == 0, "phasing: cable rupture dynamic factor >= 1 (0 = 1)");
        cab_off = {a.bar};
        break;
    }
    // ledgers: a re-tension is a tensioning, S(Po);
    // a thermal change or a rupture an action on the structure, S(G)
    origin_now = (a.kind == K::Length || a.kind == K::Horizontal || a.kind == K::EndTension) ? 1 : 0;
    CaseResult inc = solve_increment({});
    k1_losses(inc);
    if (a.kind == K::Rupture && a.value > 1) {
      // the accidental peak: (phi - 1) x the static release, on top of the state
      const double f = a.value - 1;
      CaseResult pk = inc;
      scale_result(pk, f);
      if (!rupture_peak_set) rupture_peak = std::move(pk);
      else accumulate(rupture_peak, model, pk, model);
      rupture_peak_set = true;
    }
  }

  // grouting (INJECTER CABLE): the tendon is bonded from now on
  for (Id tid : ph.inject) {
    tstate[tid].unbonded = false;
    tstate[tid].grouted.clear();
  }
  // partial grouting (INJECTER CABLE PARTIELLEMENT): bonded on the
  // interval, sliding elsewhere
  for (const Phase::PartialGrout& pg : ph.inject_partial) {
    TendonState& ts = tstate[pg.tendon];
    require(ts.unbonded, "phasing: partial grouting of tendon " + std::to_string(pg.tendon) +
                             " which is not an unbonded tensioned tendon");
    require(local, "phasing: partial grouting needs PhasingOptions::local_losses (the tendon's loss grid)");
    require(pg.s1 >= pg.s0, "phasing: partial grouting interval s0 <= s1");
    ts.grouted.push_back({pg.s0, pg.s1});
  }

  // 3b. support modifications (MODIFIER APPUI) then jacking (VERINER APPUI)
  {
    std::vector<Load> rel;
    for (const Support& ns : ph.modify_supports) {
      auto it = std::find_if(active_supports.begin(), active_supports.end(),
                             [&](const Support& s) { return s.node == ns.node && !s.implicit; });
      require(it != active_supports.end(),
              "phasing: modifying an inactive support at node " + std::to_string(ns.node));
      auto ra = reaction_acc.find(ns.node);
      Vec6 R = ra != reaction_acc.end() ? ra->second : Vec6::Zero();
      Vec6 freed = Vec6::Zero();
      for (int c = 0; c < 6; ++c)
        if (it->fixed[static_cast<std::size_t>(c)] && !ns.fixed[static_cast<std::size_t>(c)]) {
          freed(c) = R(c);
          if (ra != reaction_acc.end()) ra->second(c) = 0;
          for_states([&](CaseResult& cs) {
            for (Reaction& re : cs.reactions)  // a free component reacts no more
              if (re.node == ns.node) re.R(c) = 0;
          });
        }
      *it = ns;
      if (freed.squaredNorm() > 0)
        rel.push_back(Load::nodal(ns.node, -freed.segment<3>(0), -freed.segment<3>(3)));
    }
    if (!rel.empty()) {
      invalidate_ledgers("phase " + ph.name + " : modification d'appui");
      origin_now = 0;
      CaseResult inc = solve_increment(rel);
      k1_losses(inc);
    }
    std::vector<Load> jack;
    for (Id nid : ph.jack_support_nodes) {
      auto it = std::find_if(active_supports.begin(), active_supports.end(),
                             [&](const Support& s) { return s.node == nid && !s.implicit; });
      require(it != active_supports.end(),
              "phasing: jacking an inactive support at node " + std::to_string(nid));
      require(it->rot == 0.0, "phasing: jacking a ROT support is not supported");
      const Vec6 u = master.node_disp[static_cast<std::size_t>(model.node_index(nid))];
      std::array<std::optional<double>, 6> d{};
      bool any = false;
      for (std::size_t c = 0; c < 6; ++c)
        if (it->fixed[c] && u(static_cast<int>(c)) != 0) {
          d[c] = -u(static_cast<int>(c));
          any = true;
        }
      if (any) jack.push_back(Load::support_disp(nid, d));
    }
    if (!jack.empty()) {
      invalidate_ledgers("phase " + ph.name + " : verinage d'appui");
      origin_now = 0;
      CaseResult inc = solve_increment(jack);
      k1_losses(inc);
    }
  }

  // 4. case loads
  {
    std::vector<Load> loads;
    for (Id cid : ph.cases) {
      const LoadCase& lc = model.load_case(cid);
      // a span load on a bar not (yet / any more) active loads nothing
      for (const Load& l : lc.loads) {
        using K = Load::Kind;
        bool on_bar = l.kind == K::BarUdl || l.kind == K::BarLinear || l.kind == K::BarParabolic || l.kind == K::BarPoint ||
                      l.kind == K::BarThermal || l.kind == K::Prestrain ||
                      l.kind == K::PrestrainPoint || l.kind == K::PrestrainLinear || l.kind == K::PrestrainParabolic;
        if (!on_bar || active_bars.count(l.bar)) loads.push_back(l);
      }
    }
    if (!loads.empty()) {
      origin_now = 0;
      CaseResult inc = solve_increment(loads);
      k1_losses(inc);
    }
  }

  // 5. tensioning / detensioning, group after group (TensionAction::group):
  //    simultaneous within a group (K = 0 among the acted-on tendons — the
  //    jack compensates its own action) while the OTHER live tendons, those
  //    of the earlier groups included, take the K = 1 elastic pass
  std::vector<int> groups;  // in order of first occurrence
  {
    std::unordered_set<int> seen_groups;
    for (const TensionAction& ta : ph.tension)
      if (seen_groups.insert(ta.group).second) groups.push_back(ta.group);
  }
  for (int grp : groups) {
    std::vector<Load> loads;
    // (tendon id, live-after-this-phase): acted-on tendons are excluded
    // from the K = 1 pass of their own increment, then set to this state.
    std::vector<std::pair<Id, bool>> acted;
    std::map<Id, std::pair<double, int>> units;  // NON SIMUL tendons: unit fraction, n
    // released pretension strands: bonded, so their own shortening in the
    // release increment is a loss (they take the K = 1 pass below)
    std::vector<Id> released;
    std::vector<std::pair<Id, bool>> sh_acted;  // shell tendons acted on (K = 0 among the group)
    for (const TensionAction& ta : ph.tension) {
      if (ta.group != grp) continue;
      if (ta.shell) {
        // a SHELL tendon (bonded, jacked): the share of its P0 profile; a
        // target below the applied fraction releases the same share of the
        // current effective force (the spring-back of tensioning)
        const ShellTendon& t = shell_tendon(ta.tendon);
        require(ta.mode == TensionAction::Mode::Jack, "phasing: shell tendon " + std::to_string(t.id) +
                                                          ": jacking only (no bed, no release)");
        ShellTendonState& st = sh_t[t.id];
        double dfe = 0;
        if (ta.fraction < st.frac_applied - 1e-12) {
          const double target = std::max(0.0, ta.fraction);
          dfe = -(st.frac_applied - target) / st.frac_applied * st.frac_eff;
          st.frac_applied = target;
        } else {
          dfe = ta.fraction - st.frac_applied;
          if (dfe <= 0) continue;
          st.frac_applied = ta.fraction;
        }
        st.frac_eff += dfe;
        loads.push_back(Load::shell_prestress(t.id, dfe));
        st.t_tension = ph.time;
        st.sigma_pi = t.Ap > 0 ? st.frac_eff * sh_mean.at(t.id) / t.Ap : 0.0;
        st.relax_cum = 0;
        sh_acted.push_back({t.id, st.frac_applied > 0});
        continue;
      }
      const Tendon& t = model.tendon(ta.tendon);
      TendonState& st = tstate[t.id];
      if (ta.mode == TensionAction::Mode::Bed) {
        // tensioned against the bed: the strand carries P0, the structure
        // nothing; its relaxation clock starts now
        require(!st.live && !st.on_bed,
                "phasing: TENDRE_BANC on tendon " + std::to_string(t.id) + " already tensioned");
        st.on_bed = true;
        st.frac_applied = 1.0;
        st.frac_eff = 1.0;
        st.t_tension = ph.time;
        st.sigma_pi = t.losses.Ap > 0 ? t.jack_force / t.losses.Ap : 0.0;
        st.relax_cum = 0;
        continue;
      }
      for (Id bid : t.bars)
        require(active_bars.count(bid),
                "phasing: tensioning tendon on inactive bar " + std::to_string(bid));
      if (ta.mode == TensionAction::Mode::Release) {
        // released from the bed: the strand force, less the thermal loss of
        // a heat cure (EN 1992-1-1 §10.5.2), is transmitted by bond along
        // the tendon's force profile; the relaxation clock runs on
        require(st.on_bed, "phasing: RELACHER_BANC on tendon " + std::to_string(t.id) +
                               " not tensioned on a bed");
        st.on_bed = false;
        if (t.losses.Ap > 0 && t.jack_force > 0)
          st.frac_eff -= t.losses.thermal_loss * t.losses.Ap / t.jack_force;
        tendon_loads_scaled(t, st.frac_eff, loads);
        released.push_back(t.id);
        acted.push_back({t.id, true});
        continue;
      }
      if (ta.fraction < st.frac_applied - 1e-12) {
        // DETENSION down to the target fraction (a negative target = full
        // release, the ST1 `DETENDRE` marker). The jack lets go of the
        // released share of the CURRENT EFFECTIVE force: the structure
        // receives the exact opposite of the equivalent loads it carries for
        // that share — the spring-back symmetric to tensioning. In an inert
        // (creep-free) history a tension -> full-detension cycle therefore
        // cancels exactly. A partial detension keeps the tendon live and
        // restarts its relaxation clock at the new effective stress; a full
        // detension takes it out of service (it may be re-tensioned later,
        // which is a fresh jacking).
        require(st.frac_applied > 0, "phasing: detensioning tendon " + std::to_string(t.id) +
                                         " that was never tensioned");
        double target = std::max(0.0, ta.fraction);
        double share = (st.frac_applied - target) / st.frac_applied;
        double dfrac_eff = -share * st.frac_eff;
        tendon_loads_scaled(t, dfrac_eff, loads);
        if (local && !st.L.empty()) {  // the same share of the local losses comes back
          std::vector<double> back(st.L.size());
          for (std::size_t k = 0; k < back.size(); ++k) back[k] = share * st.L[k];
          for (const Load& l : local_loss_loads(t, back, 1.0)) loads.push_back(l);
          for (double& v : st.L) v *= 1.0 - share;
          for (double& v : st.relax_k) v *= 1.0 - share;
        }
        st.frac_applied = target;
        st.frac_eff += dfrac_eff;
        bool live_after = target > 0;
        st.t_tension = ph.time;
        st.sigma_pi = (live_after && t.losses.Ap > 0)
                          ? st.frac_eff * mean_force[t.id] / t.losses.Ap
                          : 0.0;
        st.relax_cum = 0;
        acted.push_back({t.id, live_after});
        continue;
      }
      double dfrac = ta.fraction - st.frac_applied;
      if (dfrac <= 0) continue;
      st.unbonded = ta.unbonded;
      if (opt.exact_losses && t.losses.non_simul > 1 && t.losses.Ap > 0) {
        // NON SIMUL n, exact: its units jacked one after the other (below)
        units[t.id] = {dfrac / t.losses.non_simul, t.losses.non_simul};
      } else {
        tendon_loads_scaled(t, dfrac, loads);
      }
      st.frac_applied = ta.fraction;
      st.frac_eff += dfrac;
      // relaxation base: mean effective stress just after this jacking, and
      // the date its relaxation clock starts from
      st.t_tension = ph.time;
      st.sigma_pi = t.losses.Ap > 0 ? st.frac_eff * mean_force[t.id] / t.losses.Ap : 0.0;
      st.relax_cum = 0;
      acted.push_back({t.id, true});
    }
    if (!units.empty()) {
      // NON SIMUL n (annex B1.1.3), exact: unit j of n is jacked while the
      // j - 1 units before it are anchored. They follow the concrete at the
      // tendon's level with (j - 1) / n of its area — they lose, and their
      // stiffness takes part in what the next unit shortens (the
      // (n - 1) / (2 n) of ST1 is the limit Ep Ap << E A of this sequence).
      // The first units go with the phase's other jackings; the tendons
      // acted on simultaneously are held by their jacks meanwhile.
      for (const auto& [id, u] : units)
        tendon_loads_scaled(model.tendon(id), u.first, loads);
    }
    if (!loads.empty()) {
      for (const auto& a : sh_acted) sh_t[a.first].live = false;
      origin_now = 1;
      CaseResult inc = solve_increment(loads);
      // K = 1 elastic losses for the live tendons NOT acted on in this phase
      // (a released pretension strand is bonded: it takes its own)
      for (const auto& a : acted) tstate[a.first].live = false;
      for (Id r : released) tstate[r].live = true;
      k1_losses(inc);
      int jmax = 0;
      for (const auto& [id, u] : units) jmax = std::max(jmax, u.second);
      for (int j = 2; j <= jmax; ++j) {
        std::map<Id, double> w;  // who loses: the live tendons, the anchored units
        for (const Tendon& tt : model.tendons)
          if (tstate[tt.id].live) w[tt.id] = 1.0;
        std::vector<Load> lj;
        for (const auto& [id, u] : units) {
          tstate[id].live = true;
          w[id] = std::min(j - 1, u.second) / static_cast<double>(u.second);
          if (j <= u.second)
            tendon_loads_scaled(model.tendon(id), u.first, lj);
        }
        k1_losses(solve_increment(lj), &w);
      }
      for (const auto& a : acted) tstate[a.first].live = a.second;
      for (const auto& a : sh_acted) sh_t[a.first].live = a.second;
      // NON SIMUL n (annex B1.1.3): the units jacked one after the other —
      // each the shortening of those jacked after it, (n - 1) / (2 n) of
      // the tendon's own elastic loss on average
      std::map<Id, double> ns;
      for (const auto& a : acted) {
        const int n = model.tendon(a.first).losses.non_simul;
        if (a.second && n > 1 && !units.count(a.first)) ns[a.first] = (n - 1.0) / (2.0 * n);
      }
      if (!ns.empty()) k1_losses(inc, &ns);
    }
  }

  // record phase
  PhaseResult pr;
  pr.name = ph.name;
  pr.time = ph.time;
  pr.cumulated = master;
  pr.incremental = incremental;
  for (const auto& [bid, ps] : parts)
    for (std::size_t i = 0; i < ps.size(); ++i)
      if (ps[i].active && active_bars.count(bid))
        pr.parts.push_back({bid, static_cast<int>(i), ps[i].res});
  for (const Tendon& t : model.tendons) {
    const TendonState& st = tstate.at(t.id);
    if (st.frac_applied > 0)
      pr.tendons.push_back({t.id, st.frac_applied * t.jack_force,
                            (st.frac_applied - st.frac_eff) * t.jack_force,
                            st.frac_eff * t.jack_force, st.on_bed,
                            local && !st.L.empty() ? grid_profile(t.id, st.L)
                                                   : std::vector<std::pair<double, double>>{}});
  }
  for (const ShellTendon& t : model.shell_tendons) {
    const ShellTendonState& st = sh_t.at(t.id);
    if (st.frac_applied > 0 || st.frac_eff != 0)
      pr.shell_tendons.push_back({t.id, st.frac_applied * t.P0, (st.frac_applied - st.frac_eff) * t.P0,
                                  st.frac_eff * t.P0, false, {}});
  }
  for (const auto& [bid, c] : cab) {
    if (!c.seen) continue;
    PhaseResult::CableState cst;
    cst.bar = bid;
    cst.active = c.active;
    cst.L0 = c.L0;
    if (c.active) {
      const Bar& b = model.bar(bid);
      const Vec3 d = node_now(b.nj) - node_now(b.ni);
      const Catenary2D k = solve_catenary(std::hypot(d.x(), d.y()), d.z(), c.L0, c.EA, c.w);
      cst.H = k.H;
      cst.T1 = k.T1();
      cst.T2 = k.T2();
    }
    pr.cables.push_back(cst);
  }
  if (rupture_peak_set) {
    pr.cable_rupture_peak = std::move(rupture_peak);
    rupture_peak = CaseResult{};
    rupture_peak_set = false;
  }
  if (led_on) {
    pr.origins.assign(led.begin(), led.end());
    pr.isostatic = isostatic_state();
    pr.origins_invalid = led_invalid;
  }
  if (opt.keep_structures) pr.structure = std::make_shared<const Model>(build_work());
  results.push_back(std::move(pr));
  prev_time = ph.time;
  ++phases_run;
}

// The primary (isostatic) prestress of the tendons acting on the structure at
// the stations of their chain — the section resultant of the tendon force cut
// there: N = -P c, My = -P c e_z, Mz = +P c e_y, Vz = -P c e_z', Vy = -P c
// e_y' (the transverse components of the force; keta's sign convention), c =
// cos(alpha) for an exact-angle tendon, 1 otherwise; P = the effective force
// (the profile fraction, minus the local loss field when kept). Under
// friction dMy/dx is NOT Vz: the drag acts at the tendon, a distributed couple
// P' e on the section. O(points log breakpoints) per tendon.
CaseResult PhasingSession::Impl::isostatic_state() {
  CaseResult iso = zero_result(model, opt.stations_per_bar, opt.stations_abs, plan_supports);
  for (const Tendon& t : model.tendons) {
    const TendonState& st = tstate.at(t.id);
    if (st.frac_applied <= 0 || st.on_bed || st.frac_eff == 0) continue;
    const ForceProfile& fp = profile_of(t);
    TendonGrid own;
    const TendonGrid* g = nullptr;
    if (local) {
      const auto it = tgrid.find(t.id);
      if (it != tgrid.end()) g = &it->second;
    }
    if (!g) {
      own = build_grid(model, t, fp, iso);
      g = &own;
    }
    std::unordered_map<Id, double> sense;  // +1 along the bar, -1 flipped
    for (std::size_t i = 0; i < t.bars.size(); ++i) sense[t.bars[i]] = i < t.flips.size() && t.flips[i] ? -1.0 : 1.0;
    const std::size_t n = g->s.size();
    // each point read from INSIDE its own bar: at a bar end shared with the
    // next bar the trace may kink and the force drop (a concentrated
    // deviation), so the slope and the force are the one-sided values
    const double eps = 1e-9 * std::max(1.0, tendon_chain_length(model, t));
    std::vector<double> se(n);
    for (std::size_t k = 0; k < n; ++k) {
      const bool first = k == 0 || g->bar[k - 1] != g->bar[k];
      const bool last = k + 1 == n || g->bar[k + 1] != g->bar[k];
      se[k] = g->s[k] + (first && !last ? eps : 0.0) - (last && !first ? eps : 0.0);
    }
    const bool use_L = local && st.L.size() == n;
    for (std::size_t k = 0; k < n; ++k) {
      const double s = se[k];
      const double dey = tendon_dey(t, s), dez = tendon_dez(t, s);
      const double c = t.exact_angle ? 1.0 / std::sqrt(1.0 + dey * dey + dez * dez) : 1.0;
      const double P = st.frac_eff * fp.at(s) - (use_L ? st.L[k] : 0.0);
      const double Pc = P * c, sg = sense.at(g->bar[k]);
      Station& o = iso.bars[static_cast<std::size_t>(model.bar_index(g->bar[k]))].stations[g->st[k]];
      o.N += -Pc;
      o.My += -Pc * g->ez[k];
      o.Mz += Pc * g->ey[k];
      o.Vz += -sg * Pc * dez;
      o.Vy += -sg * Pc * dey;
    }
  }
  return iso;
}

// ----------------------------------------------------------------- fork --
std::string PhasingSession::Impl::fork_reject(const Model& m2, const PhasingOptions& o2) const {
  // nodes: every node of the footprint, with its coordinates, in the same
  // relative order (the free-DOF numbering of a substructure follows it)
  {
    std::unordered_map<Id, std::size_t> at;
    at.reserve(m2.nodes.size());
    for (std::size_t k = 0; k < m2.nodes.size(); ++k) at.emplace(m2.nodes[k].id, k);
    std::size_t last = 0;
    for (std::size_t k = 0; k < model.nodes.size(); ++k) {
      const Node& n = model.nodes[k];
      const auto it = at.find(n.id);
      if (it == at.end()) return "noeud " + std::to_string(n.id) + " absent du modele repris";
      if (m2.nodes[it->second].x != n.x) return "noeud " + std::to_string(n.id) + " deplace";
      if (k > 0 && it->second <= last) return "ordre des noeuds change (noeud " + std::to_string(n.id) + ")";
      last = it->second;
    }
  }
  // bars: the same, with their section, material constants and kernels
  {
    std::unordered_map<Id, std::size_t> at;
    at.reserve(m2.bars.size());
    for (std::size_t k = 0; k < m2.bars.size(); ++k) at.emplace(m2.bars[k].id, k);
    std::size_t last = 0;
    for (std::size_t k = 0; k < model.bars.size(); ++k) {
      const Bar& b = model.bars[k];
      const auto it = at.find(b.id);
      if (it == at.end()) return "barre " + std::to_string(b.id) + " absente du modele repris";
      const Bar& c = m2.bars[it->second];
      if (!same_bar(b, c)) return "barre " + std::to_string(b.id) + " redefinie";
      if (!same_section(model.section(b.section), m2.section(c.section)))
        return "barre " + std::to_string(b.id) + " : caracteristiques changees";
      if (!same_material(model.material(b.material), m2.material(c.material)))
        return "barre " + std::to_string(b.id) + " : materiau change";
      if ((model.tapered_of(b.id) != nullptr) != (m2.tapered_of(b.id) != nullptr))
        return "barre " + std::to_string(b.id) + " : noyau exact change";
      if (k > 0 && it->second <= last) return "ordre des barres change (barre " + std::to_string(b.id) + ")";
      last = it->second;
      // the bar's creep engine and the sets the options point at
      auto in = [&](const std::set<Id>* s, bool dflt) { return s ? s->count(b.id) > 0 : dflt; };
      if (in(opt.st1_creep_bars, true) != in(o2.st1_creep_bars, true) ||
          in(opt.st1_tension_creep_bars, false) != in(o2.st1_tension_creep_bars, false))
        return "barre " + std::to_string(b.id) + " : moteur de fluage change";
    }
  }
  // the model's own supports (the active ones carry over as they stand)
  if (model.supports.size() != m2.supports.size()) return "appuis du modele changes";
  for (std::size_t k = 0; k < model.supports.size(); ++k)
    if (!equal_support(model.supports[k], m2.supports[k])) return "appuis du modele changes";
  {
    bool same = model.links.size() == m2.links.size() && model.rigid_links.size() == m2.rigid_links.size();
    for (std::size_t k = 0; same && k < model.links.size(); ++k) {
      const ElasticLink &a = model.links[k], &b = m2.links[k];
      same = a.id == b.id && a.ni == b.ni && a.nj == b.nj && a.k == b.k && a.s == b.s && a.roll == b.roll &&
             a.axes.has_value() == b.axes.has_value() && (!a.axes || *a.axes == *b.axes);
    }
    for (std::size_t k = 0; same && k < model.rigid_links.size(); ++k) {
      const RigidLink &a = model.rigid_links[k], &b = m2.rigid_links[k];
      same = a.id == b.id && a.master == b.master && a.slaves == b.slaves && a.dofs == b.dofs;
    }
    if (!same) return "liaisons changees";
    same = model.composite.size() == m2.composite.size();
    for (auto a = model.composite.begin(), b = m2.composite.begin(); same && a != model.composite.end(); ++a, ++b) {
      same = a->first == b->first && a->second.parts.size() == b->second.parts.size();
      for (std::size_t i = 0; same && i < a->second.parts.size(); ++i) {
        const SectionPart &x = a->second.parts[i], &y = b->second.parts[i];
        same = x.material == y.material && x.z == y.z && same_props(x.props, y.props);
      }
    }
    if (!same) return "barres mixtes changees";
  }
  if (model.shells.size() != m2.shells.size()) return "coques changees";
  for (std::size_t k = 0; k < model.shells.size(); ++k)
    if (model.shells[k].id != m2.shells[k].id) return "coques changees";
  for (const ShellTendon& a : model.shell_tendons)
    for (const ShellTendon& b : m2.shell_tendons)
      if (a.id == b.id && (a.pts != b.pts || a.P0 != b.P0 || a.mu != b.mu || a.k != b.k || a.jack != b.jack ||
                           a.slab != b.slab || a.Ap != b.Ap || a.Ep != b.Ep || a.fpk != b.fpk ||
                           a.rho1000 != b.rho1000 || a.relaxation_class != b.relaxation_class))
        return "cable de coque " + std::to_string(a.id) + " change";
  // the options the increments depend on
  if (opt.stations_per_bar != o2.stations_per_bar || opt.creep_engine != o2.creep_engine ||
      opt.kelvin_units != o2.kelvin_units || opt.creep_along_bar != o2.creep_along_bar ||
      opt.creep_substeps_per_decade != o2.creep_substeps_per_decade ||
      opt.exact_losses != o2.exact_losses || opt.born_deformed != o2.born_deformed ||
      opt.local_losses != o2.local_losses || opt.st1_fibre_layout != o2.st1_fibre_layout ||
      (opt.modulus_at != nullptr) != (o2.modulus_at != nullptr))
    return "options du phasage changees";
  // the station layout of every bar of the footprint
  const CaseResult layout = zero_result(m2, o2.stations_per_bar, o2.stations_abs, plan_supports);
  for (const BarResult& br : master.bars) {
    const BarResult& c = layout.bars[static_cast<std::size_t>(m2.bar_index(br.bar))];
    if (c.stations.size() != br.stations.size())
      return "barre " + std::to_string(br.bar) + " : sections d'etude changees";
    for (std::size_t q = 0; q < c.stations.size(); ++q)
      if (c.stations[q].x != br.stations[q].x)
        return "barre " + std::to_string(br.bar) + " : sections d'etude changees";
  }
  // the tendons already acted on: same tendon, same loss grid
  for (const Tendon& t : model.tendons) {
    const TendonState& ts = tstate.at(t.id);
    const bool acted = ts.frac_applied != 0 || ts.frac_eff != 0 || ts.live || ts.on_bed ||
                       ts.unbonded || !ts.L.empty();
    if (!acted) continue;
    const Tendon* c = nullptr;
    for (const Tendon& x : m2.tendons)
      if (x.id == t.id) c = &x;
    if (!c) return "cable " + std::to_string(t.id) + " absent du modele repris";
    if (!same_tendon(t, *c)) return "cable " + std::to_string(t.id) + " redefini";
    if (local && !ts.L.empty()) {
      const TendonGrid g = build_grid(m2, *c, tendon_force_profile(m2, *c), layout);
      if (g.s.size() != ts.L.size()) return "cable " + std::to_string(t.id) + " : grille de pertes changee";
    }
  }
  // the pseudo-plastic beds with a state
  for (const auto& [bid, pb] : pse) {
    (void)pb;
    if (!o2.pse || !o2.pse->count(bid)) return "barre " + std::to_string(bid) + " : sol plastique absent";
  }
  return {};
}

std::unique_ptr<PhasingSession::Impl> PhasingSession::Impl::fork(const Model& m2,
                                                                 const PhasingOptions& o2) const {
  const std::string why = fork_reject(m2, o2);
  require(why.empty(), "phasing: reprise impossible — " + why);
  auto c = std::make_unique<Impl>(m2, o2, nullptr);
  // the state, re-keyed on the new model where it was dense
  c->plan_supports = plan_supports;
  // a full-layout state re-keyed on the new model (master and each ledger)
  auto rekey = [&](const CaseResult& src) {
    CaseResult out = zero_result(c->model, c->opt.stations_per_bar, c->opt.stations_abs, c->plan_supports);
    for (std::size_t k = 0; k < model.nodes.size(); ++k)
      out.node_disp[static_cast<std::size_t>(c->model.node_index(model.nodes[k].id))] = src.node_disp[k];
    std::unordered_map<Id, const Reaction*> at;
    for (const Reaction& re : src.reactions) at.emplace(re.node, &re);
    for (Reaction& re : out.reactions) {
      const auto it = at.find(re.node);
      if (it != at.end()) re.R = it->second->R;
    }
    for (std::size_t k = 0; k < model.bars.size(); ++k)
      out.bars[static_cast<std::size_t>(c->model.bar_index(model.bars[k].id))] = src.bars[k];
    out.shells = src.shells;
    out.links = src.links;
    out.non_superposable = src.non_superposable;
    return out;
  };
  c->master = rekey(master);
  for (const Reaction& re : c->master.reactions) c->reaction_nodes.insert(re.node);
  for (std::size_t k = 0; k < model.bars.size(); ++k)
    c->bar_nodal_acc[static_cast<std::size_t>(c->model.bar_index(model.bars[k].id))] = bar_nodal_acc[k];
  c->led_on = led_on && c->opt.origin_ledgers;
  if (c->led_on)
    for (std::size_t o = 0; o < 3; ++o) c->led[o] = rekey(led[o]);
  c->led_invalid = led_invalid;
  c->cab = cab;
  c->bar_ref = bar_ref;
  c->imperf = imperf;
  c->active_links = active_links;
  c->active_rigid = active_rigid;
  c->link_nodal_acc = link_nodal_acc;
  c->parts = parts;
  c->deferred_parts = deferred_parts;
  c->reaction_acc = reaction_acc;
  for (Id b : active_bars) c->activate_bar(b);
  c->active_supports = active_supports;
  for (const Support& s : active_supports)
    if (!s.implicit) ++c->struct_deg[s.node];
  c->cast_time = cast_time;
  c->bar_creep = bar_creep;
  c->bar_st1 = bar_st1;
  c->prev_time = prev_time;
  c->phases_run = phases_run;
  for (const auto& [bid, pb] : pse) {
    PseBar& q = c->pse[bid];
    q = pb;
    q.soil = &c->opt.pse->at(bid);
  }
  // the tendons: their states, and the derived data of the unchanged ones
  // (profiles, unit loads, loss bases — pure functions of the tendon and
  // its chain, the costly part of a phasing's set-up)
  {
    std::unordered_map<Id, const Tendon*> mine;
    for (const Tendon& t : model.tendons) mine.emplace(t.id, &t);
    for (const Tendon& t : c->model.tendons) {
      const auto it = mine.find(t.id);
      if (it == mine.end() || !same_tendon(*it->second, t)) continue;
      c->tstate[t.id] = tstate.at(t.id);
      if (auto f = force_profile.find(t.id); f != force_profile.end()) c->force_profile.emplace(t.id, f->second);
      if (auto u = unit_loads.find(t.id); u != unit_loads.end()) c->unit_loads.emplace(t.id, u->second);
      if (auto l = loss_basis.find(t.id); l != loss_basis.end()) c->loss_basis.emplace(t.id, l->second);
    }
  }
  c->init_tendons();
  // the facets' creep memory and the shell tendons (same facets, same
  // tendons: fork_reject compares the shells)
  c->init_shells();
  c->shell_mem = shell_mem;
  for (const auto& [id, st] : sh_t)
    if (c->sh_t.count(id)) c->sh_t[id] = st;
  return c;
}

// ------------------------------------------------------------- the API --
PhasingSession::PhasingSession(const Model& model, const PhasingOptions& opt,
                               const std::set<Id>& active_from_start, PhasedFactorCache* cache)
    : impl_(std::make_unique<Impl>(model, opt, cache)) {
  impl_->start(active_from_start);
}
PhasingSession::PhasingSession(std::unique_ptr<Impl> impl) : impl_(std::move(impl)) {}
PhasingSession::~PhasingSession() = default;
PhasingSession::PhasingSession(PhasingSession&&) noexcept = default;
PhasingSession& PhasingSession::operator=(PhasingSession&&) noexcept = default;

void PhasingSession::run(const PhasePlan& plan) { impl_->run(plan); }
void PhasingSession::defer_parts(const std::set<std::pair<Id, int>>& later) {
  require(impl_->phases_run == 0, "phasing: defer_parts before the first run only");
  impl_->deferred_parts = later;
  for (const auto& [b, i] : later) {
    auto it = impl_->parts.find(b);
    if (it != impl_->parts.end() && i >= 0 && static_cast<std::size_t>(i) < it->second.size())
      it->second[static_cast<std::size_t>(i)].active = false;
  }
}
void PhasingSession::defer_links(const std::set<Id>& elastic, const std::set<Id>& rigid) {
  require(impl_->phases_run == 0, "phasing: defer_links before the first run only");
  for (Id l : elastic) impl_->active_links.erase(l);
  for (Id l : rigid) impl_->active_rigid.erase(l);
}
const std::vector<PhaseResult>& PhasingSession::results() const { return impl_->results; }
std::vector<PhaseResult> PhasingSession::take_results() { return std::move(impl_->results); }
std::size_t PhasingSession::phases_run() const { return impl_->phases_run; }
double PhasingSession::time() const { return impl_->prev_time; }
std::unique_ptr<PhasingSession> PhasingSession::fork(const Model& model, const PhasingOptions& opt) const {
  return std::unique_ptr<PhasingSession>(new PhasingSession(impl_->fork(model, opt)));
}
std::string PhasingSession::fork_reject(const Model& model, const PhasingOptions& opt) const {
  return impl_->fork_reject(model, opt);
}
const ForceProfile& PhasingSession::force_profile(Id tendon) {
  return impl_->profile_of(impl_->model.tendon(tendon));
}
int PhasingSession::kelvin_fits() const { return impl_->fit_cache.fits; }
FactorCache::Stats PhasingSession::cache_stats() const { return impl_->cache_->stats(); }
const Model& PhasingSession::model() const { return impl_->model; }

PhasedResults solve_phased(const Model& model, const PhasePlan& plan, const PhasingOptions& opt) {
  require(model.continua.empty(), "solve_phased: continuum elements are not supported (linear statics, modal and dynamics take them)");
  PhasedFactorCache cache;
  return solve_phased_cached(model, plan, opt, cache);
}

PhasedResults solve_phased_cached(const Model& model, const PhasePlan& plan,
                                  const PhasingOptions& opt, PhasedFactorCache& cache) {
  require(model.continua.empty(), "solve_phased_cached: continuum elements are not supported (linear statics, modal and dynamics take them)");
  hooks_check(opt.hooks);
  model.validate();
  require(!plan.phases.empty(), "phasing: empty plan");
  // Bars never named in any activate list are active from the start.
  std::set<Id> ever_activated;
  for (const Phase& ph : plan.phases)
    for (Id b : ph.activate_bars) ever_activated.insert(b);
  std::set<Id> from_start;
  for (const Bar& b : model.bars)
    if (!ever_activated.count(b.id)) from_start.insert(b.id);
  PhasingSession session(model, opt, from_start, &cache);
  std::set<Id> later_links, later_rigid;
  for (const Phase& ph : plan.phases) {
    later_links.insert(ph.activate_links.begin(), ph.activate_links.end());
    later_rigid.insert(ph.activate_rigid_links.begin(), ph.activate_rigid_links.end());
  }
  session.defer_links(later_links, later_rigid);
  std::set<std::pair<Id, int>> later_parts;
  for (const Phase& ph : plan.phases) later_parts.insert(ph.activate_parts.begin(), ph.activate_parts.end());
  session.defer_parts(later_parts);
  session.run(plan);
  PhasedResults out;
  out.phases = session.take_results();
  out.kelvin_fits = session.kelvin_fits();
  out.cache_stats = cache.stats();
  return out;
}

// ------------------------------------------------- weighted states (ETAT PONDERE)
namespace {
// o (op) s, scalar by scalar, over the fields the ledgers accumulate (same layout)
template <class F>
void zip_state(CaseResult& o, const CaseResult& s, F f) {
  require(o.node_disp.size() == s.node_disp.size() && o.reactions.size() == s.reactions.size() &&
              o.bars.size() == s.bars.size() && o.links.size() == s.links.size(),
          "phasing: origin ledgers of another layout");
  for (std::size_t i = 0; i < o.node_disp.size(); ++i)
    for (int c = 0; c < 6; ++c) f(o.node_disp[i](c), s.node_disp[i](c));
  for (std::size_t i = 0; i < o.reactions.size(); ++i)
    for (int c = 0; c < 6; ++c) f(o.reactions[i].R(c), s.reactions[i].R(c));
  for (std::size_t i = 0; i < o.bars.size(); ++i) {
    BarResult& a = o.bars[i];
    const BarResult& b = s.bars[i];
    for (int c = 0; c < 12; ++c) f(a.end_forces_local(c), b.end_forces_local(c));
    require(a.stations.size() == b.stations.size(), "phasing: origin ledgers of another layout");
    for (std::size_t k = 0; k < a.stations.size(); ++k) {
      Station& x = a.stations[k];
      const Station& y = b.stations[k];
      f(x.N, y.N);
      f(x.Vy, y.Vy);
      f(x.Vz, y.Vz);
      f(x.T, y.T);
      f(x.My, y.My);
      f(x.Mz, y.Mz);
      for (int c = 0; c < 3; ++c) {
        f(x.defl(c), y.defl(c));
        f(x.rot(c), y.rot(c));
      }
    }
  }
  for (std::size_t i = 0; i < o.links.size(); ++i)
    for (int c = 0; c < 6; ++c) {
      f(o.links[i].force(c), s.links[i].force(c));
      f(o.links[i].deformation(c), s.links[i].deformation(c));
    }
}
void check_ledgers(const PhaseResult& ph) {
  require(ph.origins.size() == 3, "phasing: phase " + ph.name + " run without PhasingOptions::origin_ledgers");
  require(ph.origins_invalid.empty(), "phasing: origins no longer separable — " + ph.origins_invalid);
}
}  // namespace

CaseResult weighted_state(const PhaseResult& phase, double w_g, double w_po, double w_dp) {
  check_ledgers(phase);
  require(std::isfinite(w_g) && std::isfinite(w_po) && std::isfinite(w_dp), "phasing: weights must be finite");
  CaseResult out = phase.cumulated;
  const double w[3] = {w_g, w_po, w_dp};
  for (std::size_t o = 0; o < 3; ++o)
    if (w[o] != 1.0) {
      const double a = w[o] - 1.0;
      zip_state(out, phase.origins[o], [a](double& x, double y) { x += a * y; });
    }
  return out;
}

std::pair<CaseResult, CaseResult> weighted_envelope(const PhaseResult& phase, const WeightingCoefficients& c) {
  check_ledgers(phase);
  CaseResult lo = phase.cumulated, hi = phase.cumulated;
  auto term = [&](const CaseResult& S, std::pair<double, double> k) {
    const double a = k.first - 1.0, b = k.second - 1.0;
    if (a == 0 && b == 0) return;
    zip_state(lo, S, [a, b](double& x, double y) { x += std::min(a * y, b * y); });
    zip_state(hi, S, [a, b](double& x, double y) { x += std::max(a * y, b * y); });
  };
  term(phase.origins[0], c.kg);
  if (c.total) {
    CaseResult P = phase.origins[1];
    zip_state(P, phase.origins[2], [](double& x, double y) { x += y; });
    term(P, c.kpt);
  } else {
    term(phase.origins[1], c.kpo);
    term(phase.origins[2], c.kpe);
  }
  return {std::move(lo), std::move(hi)};
}

CaseResult hyperstatic_prestress(const PhaseResult& phase) {
  require(phase.origins.size() == 3, "phasing: phase " + phase.name + " run without PhasingOptions::origin_ledgers");
  CaseResult h = phase.origins[1];
  zip_state(h, phase.origins[2], [](double& x, double y) { x += y; });
  for (std::size_t i = 0; i < h.bars.size(); ++i) {
    BarResult& a = h.bars[i];
    const BarResult& b = phase.isostatic.bars[i];
    for (std::size_t k = 0; k < a.stations.size() && k < b.stations.size(); ++k) {
      a.stations[k].N -= b.stations[k].N;
      a.stations[k].Vy -= b.stations[k].Vy;
      a.stations[k].Vz -= b.stations[k].Vz;
      a.stations[k].My -= b.stations[k].My;
      a.stations[k].Mz -= b.stations[k].Mz;
    }
  }
  return h;
}

}  // namespace keta


// Moving loads: influence lines by a factor-once ADJOINT unit-load sweep (one
// solve K w = g for the response gradient g, then eta_k = w.r_k per load
// position — Muller-Breslau by Maxwell-Betti; see unit_load_sweep),
// lean scalar extraction (displacement / reaction / bar force at a station),
// sign-split influence-line algebra (exact zero-crossing trapezoids), convoy
// patterns (axles + uniform segments) slid for extremes, and UDL placement on
// the sign-matched zones only.
#pragma once
#include "keta/statics.h"

namespace keta {

// Travel path: a contiguous chain of bars; flips handled automatically.
struct PathStation {
  Id bar = 0;
  double x_local = 0;  // abscissa on the bar (from its node i)
  double s = 0;        // curvilinear abscissa along the path
};
struct LoadPath {
  std::vector<PathStation> stations;
  double length = 0;
};
// n_per_bar stations per bar (>= 2, ends included). Throws on a broken chain.
LoadPath build_path(const Model& model, const std::vector<Id>& bars, int n_per_bar = 11);

struct ResponseKind {
  // BarDisplacement: the displacement / section rotation of a bar SECTION
  // (a study section of ST1), GLOBAL components 0..5 = Dx Dy Dz Rx Ry Rz,
  // read from the same station recovery as the static results.
  enum class What { Displacement, Reaction, BarForce, BarDisplacement };
  What what = What::Displacement;
  Id node = 0;   // Displacement / Reaction
  int comp = 0;  // 0..5 (dof) or station component 0..5 = N,Vy,Vz,T,My,Mz
  Id bar = 0;    // BarForce / BarDisplacement
  double x = 0;  // abscissa on the target bar [m]
  static ResponseKind displacement(Id node, int comp);
  static ResponseKind reaction(Id node, int comp);
  static ResponseKind bar_force(Id bar, double x, int comp);
  static ResponseKind bar_displacement(Id bar, double x, int comp);
};

struct InfluenceLine {
  std::vector<double> s;    // path abscissa
  std::vector<double> eta;  // response per unit load
  // Slopes d eta / ds just left / right of each sample (empty: a linear
  // line). Within a prismatic bar the line is the deflected shape of the dual
  // problem — a cubic between samples: with the slopes the line is the exact
  // cubic Hermite interpolant (right slope at s_i, left slope at s_i+1), and
  // `at`, the integrals and the sign split are exact for it.
  std::vector<double> dl, dr;
  bool cubic() const { return dl.size() == s.size() && dr.size() == s.size(); }
  // the line on interval i (s_i <= x <= s_i+1)
  double on(std::size_t i, double x) const;
  double at(double x) const;                       // interpolation (cubic when slopes are set)
  double integral(double a, double b) const;       // exact trapezoid, clipped
  double integral_pos() const;                     // over eta > 0 zones (exact split)
  double integral_neg() const;
  std::pair<double, double> point_extreme() const; // (max, min) single unit load
};

// Prefix-integral tables of an influence line: any window integral (and the
// sign-split areas) in O(log n) instead of an O(n) walk. Built once in O(n);
// a window is the difference of two clipped prefix evaluations, so the result
// equals InfluenceLine::integral up to round-off. Borrows the line: it must
// outlive the prefix.
class IlPrefix {
 public:
  explicit IlPrefix(const InfluenceLine& il);
  explicit IlPrefix(const InfluenceLine&&) = delete;  // would dangle
  double integral(double a, double b) const;  // exact trapezoid, clipped
  // prefix_at(q) when the interval of q is known: s[lo] <= q < s[lo + 1]
  double prefix_in(std::size_t lo, double q) const;
  double total() const { return cum_.empty() ? 0.0 : cum_.back(); }
 private:
  double prefix_at(double q) const;           // integral over [s.front(), clamp(q)]
  const InfluenceLine& il_;
  std::vector<double> cum_;                   // cum_[i] = integral over [s0, s_i]
};

// Unit load of direction `dir` (e.g. (0,0,-1)) swept along the path.
InfluenceLine influence_line(const Model& model, const LoadPath& path, const Vec3& dir,
                             const ResponseKind& response);
// Same sweep against a caller-owned factor cache: an unchanged stiffness
// (same model, any path / target / direction) reuses the factorization
// (FactorCache's exact-compare ladder), so only the multi-RHS back-substitution
// is paid. Results are bit-identical to the cache-free overload.
// `hooks` (keta/hooks.h, nullable): see unit_load_sweep.
InfluenceLine influence_line(const Model& model, const LoadPath& path, const Vec3& dir,
                             const ResponseKind& response, FactorCache& cache,
                             const SolveHooks* hooks = nullptr);

// Lean scalar extraction of ONE response from a solved free-DOF state `u` and
// the lowered loads `cl` that produced it. Shared by the influence LINE sweep
// (influence.h) and the influence SURFACE sweep (influence_surface.h) — the
// unit-load position only enters through `cl`, so the same extractor serves a
// 1-D path and a 2-D node grid. Throws on a non-finite ordinate.
// Note: the reaction branch accumulates the stiffness blocks incident to the
// node — bars (bar_ke_assembled: Winkler bed and rigid-arm eccentricity
// included), facet shells and elastic links — then subtracts the load applied
// at the node.
double extract_response(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                        const VecX& u, const ResponseKind& response);
// Several BarDisplacement responses of ONE bar read from a single station
// recovery (the sections asked added to extract_response's grid): the
// recovery is exact between the load cuts, so each value equals
// extract_response's to round-off.
void extract_bar_displacements(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                               const VecX& u, Id bar, const std::vector<const ResponseKind*>& rs,
                               std::vector<double>& out);

// Full DOF indices `extract_response` reads for this response (the support of
// its gradient with respect to the displacements).
std::vector<int> response_dofs(const Model& model, const Assembled& asm_, const ResponseKind& r);

// The unit-load sweep behind influence_line and influence_surface, by the
// adjoint (Maxwell-Betti) route. The response is affine in the free
// displacements, value = g.u + h(load), with g exact from unit probes on the
// few DOFs the extractor reads. K being symmetric, g.K^-1.r = (K^-1.g).r: ONE
// solve K w = g replaces one back-substitution per load position, and each
// ordinate is eta_k = w.r_k + h_k with r_k (<= 12 non-zeros) lowered into one
// reused CaseLoads. h_k is the load's own local term: minus the applied nodal
// force for a reaction at the loaded node, the span record for a bar force on
// the loaded bar, zero otherwise. Cost: one factorization (reused through
// `fac` when K is unchanged), one solve, O(1) per position (a reaction target
// adds its incident elements once). Mathematically identical to the direct
// multi-RHS sweep (equal up to round-off). `unit` holds NodalForce or BarPoint
// loads only; throws on a non-finite ordinate. `hooks` (keta/hooks.h,
// nullable) is polled before the factorization, before the adjoint solve and
// every 64 positions ("ordinates k/n"); it never changes the numbers.
std::vector<double> unit_load_sweep(const Model& model, const std::vector<Load>& unit,
                                    const ResponseKind& response, FactorCache& fac,
                                    const SolveHooks* hooks = nullptr);

// ---------------- convoys ----------------

struct ConvoyPattern {
  // Point axles at offsets from the moving reference (offset, load>0 factor).
  std::vector<std::pair<double, double>> axles;
  // Uniform segments (start offset, length, TOTAL load).
  std::vector<std::array<double, 3>> segments;
  ConvoyPattern mirrored() const;  // travel in the other direction
};

struct Extremes {
  double max = 0, min = 0;
  double r_max = 0, r_min = 0;  // reference positions producing them
};
// Slide the convoy along the line (both keeping loads on-deck and overhanging;
// off-deck axles contribute zero). steps = sweep resolution (>= 1: steps + 1
// reference positions, both ends included). Segment loads integrate through an
// IlPrefix, so a sweep costs O(steps * (axles + segments) * log n).
Extremes convoy_extreme(const InfluenceLine& il, const ConvoyPattern& convoy, int steps = 500);
// UDL free to load ANY number of zones: max = q * integral_pos, min = q * integral_neg.
Extremes udl_extreme(const InfluenceLine& il, double q);

}  // namespace keta

namespace keta {
// Moving-convoy nodal force history: the convoy traverses the path at speed v;
// at each time step every axle inside the deck is lowered as a BarPoint load
// and its equivalent nodal forces accumulated on the free DOFs. Feed the
// result to time_integrate / time_history_modal for convoy dynamics
// (resonance at v_cr = D f1 for regular axle spacings D).
std::function<VecX(int)> moving_convoy_force(const Model& model, const LoadPath& path,
                                             const ConvoyPattern& convoy, const Vec3& dir,
                                             double speed, double dt);
// The same convoy at constant acceleration: reference r(t) = r0 + v0 t + a t^2 / 2,
// a braking convoy (a < 0) stopping at t* = v0 / |a| and staying there (it
// does not reverse). The loads are the vertical (dir) ones only: the
// longitudinal inertia force m a of the convoy is not applied.
double convoy_position_accel(double r0, double v0, double a, double t);
std::function<VecX(int)> moving_convoy_force_accel(const Model& model, const LoadPath& path,
                                                   const ConvoyPattern& convoy, const Vec3& dir, double r0,
                                                   double v0, double a, double dt);
}  // namespace keta

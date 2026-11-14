// Linear statics: multi-case one-factorization pipeline, prescribed support
// displacements by direct elimination, one-sided (uplift) supports by
// active-set iteration, exact station recovery (equilibrium integration),
// deflection recovery (double curvature integration anchored on the end
// displacements), reactions, soil pressure, section stress points,
// combinations by superposition. Per-case recovery can be spread over
// std::threads (StaticOptions::threads, opt-in, default sequential) with
// slot-indexed writes — bit-deterministic by construction.
#pragma once
#include <array>
#include <map>

#include "keta/assembly.h"
#include "keta/hooks.h"
#include "keta/linsolve.h"

namespace keta {

struct Station {
  double x = 0;   // abscissa along the bar [m]
  double N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0;
  Vec3 defl = Vec3::Zero();  // local deflections (u, v, w)
  // Local section rotations (theta_x, theta_y, theta_z): end rotation plus the
  // integrated curvature (theta_y = -dw/dx convention), linear twist.
  Vec3 rot = Vec3::Zero();
};

struct StressAtPoint {
  std::string point;      // stress-point name
  double sigma = 0;       // normal stress N/A + My z/Iy - Mz y/Iz [Pa]
  double tau = 0;         // shear estimate from Vy/Avy, Vz/Avz at this point
  double von_mises = 0;
};

struct BarResult {
  Id bar = 0;
  Vec12 end_forces_local = Vec12::Zero();  // member end forces (soil excluded)
  std::vector<Station> stations;
  // Soil pressure at stations (only when the bar has a Winkler bed):
  std::vector<Vec3> soil_p;   // px, py, pz [N/m]
  std::vector<double> soil_m; // mx [N·m/m]
};

struct Reaction {
  Id node = 0;
  Vec6 R = Vec6::Zero();
};

// Elastic node-to-node link (model.h ElasticLink), link frame: the relative
// displacement (d, dtheta) of node j with respect to node i at the spring
// point, and the spring force k (d, dtheta) (axial component positive in
// tension when the link frame x runs from i to j); node j receives -force.
struct LinkResult {
  Id link = 0;
  Vec6 deformation = Vec6::Zero();
  Vec6 force = Vec6::Zero();
};

struct CaseResult {
  Id id = 0;
  std::string name;
  std::vector<Vec6> node_disp;     // per model node, global frame
  std::vector<Reaction> reactions; // at supported nodes (declaration order)
  std::vector<BarResult> bars;
  // Per-facet-shell centroid resultants, model.shells order (empty in a
  // bar-only model). Appended field — every existing consumer is untouched.
  std::vector<ShellResult> shells;
  bool non_superposable = false;
  std::vector<LinkResult> links;  // model.links order (appended)
  // General constraints (model.constraints order): the multiplier mu of
  // each, its force on the structure being -mu c over its DOFs (the support
  // reaction convention: K u - F = mu c). Skew supports: the reaction in the
  // support frame (blocked local components) and in global components.
  std::vector<double> constraint_forces;
  struct SkewReaction {
    Id node = 0;
    Vec6 local = Vec6::Zero(), global = Vec6::Zero();
  };
  std::vector<SkewReaction> skew_reactions;
};

struct StaticResults {
  std::vector<CaseResult> cases;
  std::vector<CaseResult> combos;
  std::vector<AutoFixedDof> auto_fixed;
};

struct StaticOptions {
  int stations_per_bar = 11;
  // Custom study sections (ST1 ETUDE, see etude.h): per bar id, the sorted
  // abscissas [m] where stations are recovered. Bars absent from the map keep
  // the uniform `stations_per_bar` grid; nullptr (the default) leaves the whole
  // recovery on the uniform path. Not owned — must outlive the solve.
  const std::map<Id, std::vector<double>>* stations_abs = nullptr;
  // Multi-case PARALLEL recovery (appended). After the single factorization
  // and the (already vectorized, sequential) multi-RHS solve, the per-case
  // station/reaction/deflection recovery is spread over `threads`
  // std::threads, each case written into its own pre-allocated slot indexed
  // by case order — NO shared mutation, so the results are BIT-DETERMINISTIC
  // BY CONSTRUCTION whatever the scheduling (an indexed write order has no
  // effect on values). 1 = sequential, the DEFAULT (historical behavior);
  // 0 = auto (hardware_concurrency); n > 1 = that many threads. A typed error
  // raised inside a worker is captured and rethrown for the LOWEST failing
  // case index (deterministic). The uplift (active-set) path stays sequential
  // regardless — its cases re-factorize per iteration.
  int threads = 1;
  // Cancellation + progress (keta/hooks.h); not owned, must outlive the
  // solve. Null (the default): no polling, results unchanged.
  const SolveHooks* hooks = nullptr;
  // Internal forces only at the stations (N, V, T, M): the deflection
  // recovery (double integration over the load cuts) is skipped, `defl` and
  // `rot` stay zero. For drivers that read only the forces of an
  // intermediate solve (the bonded-tendon loss iterations of a phasing).
  bool forces_only = false;
  // Station recovery of ONE case spread over the bars (keta/parallel.h):
  // 1 = sequential (default), 0 = hardware concurrency. Each bar writes its
  // own record and the internal-force sum for the reactions runs afterwards
  // in bar order: results are bit-identical for any thread count.
  int bar_threads = 1;
  // A uniform Bar::soil bed is solved by the exact member block (bef.h) on
  // the bars it applies to (exact_bed_applies); false keeps the element's
  // Hermite bed everywhere — the discretization the nonlinear, modal and
  // stability drivers share, to compare them with the linear statics.
  bool exact_soil_beds = true;
};

// Full pipeline over all cases + combinations of the model.
StaticResults solve_static(const Model& model, const StaticOptions& opt = {});
StaticResults solve_static_cached(const Model& model, const StaticOptions& opt, FactorCache& cache);

// Single-case solve on a prebuilt assembly (used by phasing/influence/nonlinear).
struct StaticContext {
  const Model* model = nullptr;
  Assembled asm_;
  SpMat K;
  FactorCache cache;

  // `links_supported`: the caller works on the reduced space of the rigid
  // links (solve_case / solve_displacements do; ctx.K stays EXTENDED, so a
  // caller using K directly must reduce it with asm_.red). A model with links
  // is rejected otherwise.
  static StaticContext build(const Model& model, bool links_supported = false);
  // Solve one lowered case and recover a full CaseResult.
  CaseResult solve_case(const LoadCase& lc, const StaticOptions& opt);
  // Displacement-only solve (fast path for sweeps).
  VecX solve_displacements(const CaseLoads& cl);
};

// Recover a CaseResult from a free-DOF displacement vector + lowered loads.
CaseResult recover_case(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                        const VecX& u_free, const StaticOptions& opt, Id id,
                        const std::string& name);

// Linear superposition of case results (fieldwise).
// One bar's station recovery at the abscissas `xs` [m, ascending] — the exact
// routine behind every CaseResult station (internal forces by equilibrium,
// deflections and section rotations by curvature integration anchored on the
// end displacements). `f_end` / `u_local` are the member end forces and local
// end displacements. Shared with the influence-line extractor (section
// displacements, influence.h) so both read one definition.
void recover_bar_stations(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                          const Vec12& f_end, const Vec12& u_local, const std::vector<double>& xs,
                          BarResult& out);

// ---- exact constant-N beam-column member (second order, exact stability) ----
// The member of `solve_second_order` with `exact_stability`: Euler-Bernoulli
// bending under a CONSTANT axial force N (tension +), the closed-form
// solution of EI y'''' - N y'' = q whose end stiffness is the stability-function
// block `BarElement::k_exact_local`. Per bending plane the deflection y
// relative to end i solves y'' - (N/EI) y = h(x), h the first-order curvature
// (end forces + span loads + imposed curvature), integrated EXACTLY piece by
// piece between the load cuts (h is a cubic there) with the transfer
// functions cosh/sinh (tension), cos/sin (compression) and their particular
// responses (Taylor series for |N x^2 / EI| <= 4).
//   * beam_column_f0: the fixed-end forces f0(N) of the bar's span loads
//     (bending components; the axial and torsional ones are `f0_linear`'s).
//   * recover_beam_column: end forces k_exact(N) u - f0(N) (releases
//     re-injected with the exact block), stations with M = M1 + N w (w the
//     exact deflected shape, M1 the first-order moment of the exact end
//     forces) and the exact deflections / section rotations.
// Typed rejections: tapered bars, Winkler bed, concentrated dislocations,
// AXIAL span loads (N would vary along the member), and tension with
// L sqrt(N/EI) > 18 (the transfer from end i amplifies round-off by e^(kL):
// split the member).
// The transfer functions behind it: F_m(s) = sum_j alpha^j s^(2j+m)/(2j+m)!,
// m = 0..6 (F_0 = cosh/cos, F_1 = sinh/k or sin/k, F_(m+2) the response from
// rest of y'' - alpha y = s^m/m!; F_m' = F_(m-1)). Shared with warping.h.
std::array<double, 7> beam_transfer(double alpha, double s);
Vec12 beam_column_f0(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                     const Vec12& f0_linear, double N);
void recover_beam_column(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                         const Vec12& f0_exact, double N, const Vec12& u_nodal,
                         const StaticOptions& opt, BarResult& br);

CaseResult combine_results(const std::vector<const CaseResult*>& sources,
                           const std::vector<double>& factors, Id id, const std::string& name);

// Per-bar axial force at the origin station (tension-positive) — feeds Kg.
std::vector<double> axial_field(const CaseResult& r);

// Per-shell membrane resultants (nx, ny, nxy) [N/m] of a case, model.shells
// order (tension-positive) — feeds the shell part of assemble_geometric.
std::vector<Vec3> shell_membrane_field(const CaseResult& r);

// Section stress recovery at a station (closed-form stress points).
std::vector<StressAtPoint> station_stresses(const Model& model, const Bar& bar,
                                            const Station& st);

}  // namespace keta

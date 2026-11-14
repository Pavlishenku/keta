// Geometric nonlinearity:
//  * solve_second_order  — P-Delta fixed point on the axial field: the tangent
//    [K + Kg(N)] sits inside the factorization so the 1/(1 - P/Pcr)
//    amplification is captured after one tangent solve; LDLT inertia guard
//    rejects beyond-critical states (never iterates on an indefinite operator).
//  * solve_corotational  — large-displacement/large-rotation statics. Element
//    strain energy is the EXACT corotational energy (chord frame + mean-twist
//    orientation, log-map deformational rotations, cancellation-free stretch)
//    of the UNCHANGED linear Timoshenko local stiffness; nodal rotations are
//    total rotation vectors (additive updates, re-canonicalized past 0.9 pi so
//    a tip can roll through a full circle). Internal force and tangent are the
//    exact gradient/Hessian of that energy by forward-mode automatic
//    differentiation to second order (keta/jet.h): exact to round-off, so the
//    tangent is the true consistent tangent (quadratic Newton) and converged
//    equilibria carry no derivative error (the former finite differences left
//    ~1e-9 relative noise in the forces).
//    Newton + backtracking line search + load substepping. Facet SHELLS
//    enter as element-independent corotational triangles (EICR: rigid frame
//    of the current facet, deformational translations / log-map rotations,
//    energy 1/2 d^T K_l d of the linear Allman + MITC3 facet), exact gradient
//    and Hessian by AD (Jet<18>) — large rotations, small strains.
//    RIGID ARMS (Bar::ecc) are exact for finite rotations: the flexible end
//    is x_node + R(theta) a, the arm moment (R a) x F of the span loads'
//    end forces joins the node; END RELEASES through the release-condensed
//    local block (its null space keeps the rigid modes); TAPERED members
//    (CARA VAR) through their exact local block and equivalent loads.
//    WINKLER beds: ground springs along the undeformed flexible axis on the
//    absolute displacement of the current axis (chord + Hermite deflection
//    of the deformational rotations) and on the end twists e1_0 . theta —
//    an energy that linearises exactly to the Hermite bed of linear statics.
//    CABLE bars: the exact elastic catenary (keta/cable.h) at the current
//    positions of its ends — forces and their analytic tangent from
//    eval_cable, unstrained length from resolve_cable_l0, its weight in
//    full from the start (not scaled by the load factor): solve_corotational
//    only, the path-following drivers reject them. End tensions in cable_T.
//    Custom foundation blocks (bef.h): typed rejection. Newton stops at the tolerance, or within the elements'
//    round-off floor (64 eps ||K|| |X|) once the last correction is
//    negligible.
//  * arc_length_trace    — Crisfield cylindrical arc-length continuation with
//    deterministic root selection, step cutting, and LDLT inertia stamps
//    (limit/bifurcation crossings flagged).
//  * ruin_load           — certified bracket on the first failure multiplier
//    (Newton stall / indefinite converged tangent / singular factorization).
#pragma once
#include <map>
#include <functional>
#include "keta/hooks.h"
#include "keta/statics.h"

namespace keta {

struct SecondOrderOptions {
  int max_iterations = 30;
  double tol_rel_du = 1e-12;
  double tol_rel_res = 1e-10;
  int stations_per_bar = 11;
  const SolveHooks* hooks = nullptr;  // polled every P-Delta iteration (keta/hooks.h), not owned
  // Use the EXACT stability-function member tangent (ST1-style closed-form
  // beam-column) instead of the consistent cubic Kg: one element per member
  // then yields the exact second-order response — no meshing needed.
  bool exact_stability = false;
};
struct SecondOrderResults {
  CaseResult result;
  int iterations = 0;
};
SecondOrderResults solve_second_order(const Model& model, Id case_id,
                                      const SecondOrderOptions& opt = {});

// How the corotational element measures its deformational rotations from
// Rbar = E^T R E0: Log = the rotation vector log(Rbar) (keta's default, the
// standard EICR measure of Battini & Pacoste), Sin = axial(Rbar - Rbar^T)/2,
// i.e. sin(theta) a (the sine form, guarded at |s| <= 0.7 ~ 44 deg). Both
// agree to third order in the local rotation and converge to the same answer
// under mesh refinement; on one element per member at tens of degrees they
// differ by a few percent. The parametric nodes (Arc Length, Nonlinear) use
// Sin.
enum class CorotRotationMeasure { Log, Sin };

struct CorotOptions {
  CorotRotationMeasure measure = CorotRotationMeasure::Log;
  int substeps = 1;
  int max_iterations = 40;
  double tol_rel = 1e-8;
  int max_backtracks = 20;
  const SolveHooks* hooks = nullptr;  // polled every Newton step (keta/hooks.h), not owned
};
struct CorotResults {
  std::vector<Vec6> node_disp;     // translations + total rotation vectors
  std::vector<Vec12> bar_forces;   // local member end forces (deformational)
  VecX u_free;                     // raw free-DOF solution
  int iterations = 0;              // Newton steps over all substeps
  bool non_superposable = true;
  // Support reactions at the converged state (appended): element forces
  // minus external actions at each supported node's blocked / sprung DOFs,
  // global frame, supports in declaration order.
  std::vector<Reaction> reactions;
  // Cable bars (exact catenaries, appended): the end tensions (end i, end j)
  // and the horizontal tension in the deformed state; zero for other bars.
  std::vector<std::array<double, 2>> cable_T;
  std::vector<double> cable_H;
  // Skew supports (appended): reactions in the support frame and globally.
  std::vector<CaseResult::SkewReaction> skew_reactions;
};
// Loads: nodal forces/moments AND bar span loads (UDL / linear / point /
// self-weight), single case. Span loads are lowered as equivalent nodal
// forces on the CURRENT corotated frame at every residual evaluation:
// Global-frame intensities (and self-weight) keep their fixed spatial
// direction — conservative dead load — while Local-frame intensities FOLLOW
// the corotated element frame (documented choice: local = follower).
// Equivalent-load moments are spin-conjugate (J(theta)^T). The load
// stiffness — d(external force)/du of the nodal moments and span loads — is
// in the tangent (exact, by AD), which is then non-symmetric and solved by
// sparse LU; Newton starts each substep from the secant prediction of the
// previous increment. The lowering uses the undeformed length: accuracy is mesh-bounded like
// the consistent mass. Projected (PRO) bar loads and every other load kind
// are typed rejects.
// LINKS (also in arc_length_trace and ruin_load): a rigid link is the exact
// finite-rotation constraint u_s = u_m + (R(theta_m) - I) r, theta_s =
// theta_m, solved on the independent DOFs with the constraint's second
// variation in the tangent; an elastic link is the objective energy of its
// relative displacement read in the mean frame of its two nodes (arms turning
// with them). Chains, partial rotation masks and link laws are refused.
CorotResults solve_corotational(const Model& model, Id case_id, const CorotOptions& opt = {});

struct ArcLengthOptions {
  CorotRotationMeasure measure = CorotRotationMeasure::Log;
  double dl0 = 0.1;          // initial arc radius (displacement norm)
  int max_steps = 100;
  int max_iterations = 20;
  int desired_iterations = 5;  // 0 = fixed step; else dl *= sqrt(desired / used) in [1/2, 2], dl in [dl0/16, 4 dl0]
  int max_cuts = 8;
  double tol_rel = 1e-8;
  double stop_lambda = 0;      // 0 = no lambda stop
  double stop_disp_inf = 0;    // 0 = no displacement stop
  const SolveHooks* hooks = nullptr;  // polled every continuation step (keta/hooks.h), not owned
};
struct ArcPoint {
  double lambda = 0;
  VecX u;
  int n_neg = 0;      // tangent inertia at the converged point
  bool crossing = false;
};
struct ArcTrace {
  std::vector<ArcPoint> points;
};
ArcTrace arc_length_trace(const Model& model, Id case_id, const ArcLengthOptions& opt = {});

struct RuinOptions {
  double lambda_init = 1.0;
  double growth = 2.0;
  double lambda_max = 1e6;
  double tol_rel = 1e-6;
  CorotOptions corot;  // its `hooks` are polled at every probe's Newton steps;
                       // "bracketing probe k" / "bisection probe k" reported
};
struct RuinResult {
  double lambda_lower = 0, lambda_upper = 0;
  // "stall" (Newton diverged) | "indefinite" | "singular" | "snap" (the
  // tangent is indefinite midway between two definite equilibria)
  std::string oracle;
};
RuinResult ruin_load(const Model& model, Id case_id, const RuinOptions& opt = {});

// Verification seam: one bar's corotational strain energy, internal force (the
// exact energy gradient) and consistent tangent (the exact energy Hessian) at
// its 12 element DOFs q = [u_i, theta_i, u_j, theta_j] (global translations and
// total rotation vectors). Forward-mode AD to second order (keta/jet.h): exact
// to round-off, the same evaluation solve_corotational / arc_length_trace use.
struct CorotElementState {
  double energy = 0;
  Vec12 force = Vec12::Zero();
  Mat12 tangent = Mat12::Zero();
};
CorotElementState corot_element_state(const Model& model, Id bar, const Vec12& q,
                                      CorotRotationMeasure measure = CorotRotationMeasure::Log);

// Verification seam: the assembled corotational residual r(u) = F(u) - f(u)
// at load factor lam on the free DOFs, and the consistent tangent the Newton
// driver factorizes, -dr/du = K_int - lam K_load (non-symmetric once nodal
// moments or span loads make F depend on u).
struct CorotTangentProbe {
  VecX residual;
  MatX tangent;
};
// `allow_cables` (appended): exact catenary cables enter the residual and the
// tangent as in solve_corotational (the path-following drivers refuse them).
CorotTangentProbe corot_tangent_probe(const Model& model, Id case_id, const VecX& u, double lam,
                                      CorotRotationMeasure measure = CorotRotationMeasure::Log,
                                      bool allow_cables = false);

// ---------------- staged construction in large displacements ----------------
// Stages of a corotational structure: stage k activates bars (frames,
// trusses, catenary cables) and ADDS the loads of its case (they stay
// applied), optionally resetting cable unstrained lengths; each stage is
// solved to equilibrium from the previous state. A bar activated on a node
// already in place and a NEW node is placed stress-free in the deformed
// geometry: the new node continues its neighbour rigidly (the same total
// rotation, the bar's undeformed vector turned by it) — the element's
// reference is then exactly its current state, so no strain is
// introduced, whatever the rotations. A bar joining two nodes already in
// place (a closure) is activated with its own undeformed geometry between
// them as they stand: the gap strains it, as a closure segment cast to its
// design length would be. Rigid links, constraints and skew supports are
// rejected.
struct CorotStage {
  std::vector<Id> activate;                       // bars activated at this stage
  Id case_id = 0;                                 // its loads, added (0: none)
  std::vector<std::pair<Id, double>> cable_l0;    // (cable bar, new unstrained length)
};
struct CorotStagedResults {
  std::vector<std::vector<Vec6>> node_disp;  // per stage, per model node (translations + rotation vectors)
  int iterations = 0;
};
CorotStagedResults solve_corotational_staged(const Model& model, const std::vector<CorotStage>& stages,
                                             const CorotOptions& opt = {});

// ---------------- large displacements in time ----------------
// Newmark (average acceleration by default) on the corotational system:
//     M a + C v + f_int(u) = lam(t) F(u),
// f_int, its consistent tangent and the (follower) loads of the large-
// displacement solver (solve_corotational: frames, trusses, shells, links,
// and the exact elastic catenaries, their tangent rebuilt at every
// iteration, their mass lumped at the ends by assemble_mass), M the model's
// mass (lumped by default:
// translational masses, so the inertia is exact for finite rotations; the
// rotary inertia of a consistent mass is kept constant on the total rotation
// vectors — exact to first order in the rotation increments), C = a M + b K0
// (Rayleigh, b on the initial stiffness). Newton on every step with the
// effective tangent K_t - lam K_load + M / (beta dt^2) + gamma / (beta dt) C.
// Starts at rest (or from CorotDynOptions::initial_displacement), the
// initial acceleration from M a0 = lam(0) F - f_int(u0) - C v0 on the DOFs
// carrying mass (0 on the massless ones).
struct CorotDynOptions {
  double dt = 0.01;
  int steps = 100;
  double gamma = 0.5, beta = 0.25;
  double tol_rel = 1e-9;
  int max_iterations = 40;
  bool lumped_mass = true;
  double rayleigh_a = 0, rayleigh_b = 0;
  std::vector<std::pair<Id, int>> observe;  // (node, dof)
  CorotRotationMeasure measure = CorotRotationMeasure::Log;
  // Appended: initial velocities (node -> translations, rotation rates);
  // the motion then starts from them (the initial acceleration from the
  // equilibrium at t = 0, as at rest)
  std::map<Id, Vec6> initial_velocity;
  // Appended: GENERALIZED-alpha (Chung & Hulbert 1993) when 0 <= rho_inf
  // <= 1 (the spectral radius at infinite frequency; 1: no dissipation):
  // alpha_m = (2 rho - 1) / (rho + 1), alpha_f = rho / (rho + 1), gamma = 1/2
  // - alpha_m + alpha_f, beta = (1 - alpha_m + alpha_f)^2 / 4 (gamma, beta
  // above then ignored); the forces at the generalized midpoint blended
  // from both ends (Kuhl & Crisfield): r = (1 - af) (F - f)(n+1) + af (F -
  // f)(n) - M ((1 - am) a(n+1) + am a(n)) - C ((1 - af) v(n+1) + af v(n)).
  // Second order, with the high (stiff axial) modes damped — a usual scheme
  // for large-rotation dynamics. -1: Newmark as above.
  double rho_inf = -1;
  // --- appended ---
  // Initial displacements per model node (translations and total rotation
  // vectors, the layout of CorotResults::node_disp; empty = the reference
  // geometry): the motion starts from that state — typically the converged
  // static equilibrium of solve_corotational, so a cable structure starts
  // hanging in its sag instead of falling into it. The Rayleigh stiffness
  // term takes the tangent of this initial state.
  std::vector<Vec6> initial_displacement;
};
struct CorotDynResult {
  std::vector<double> t;
  std::vector<std::vector<double>> u;  // per step: observed DOFs (translations, total rotation vectors)
  std::vector<Vec6> node_disp;         // final state
  int iterations = 0;
  // --- appended ---
  std::vector<Vec6> node_vel;          // final velocities (translations, rotation-vector rates)
};
CorotDynResult corot_newmark(const Model& model, Id case_id, const std::function<double(double)>& lam,
                             const CorotDynOptions& opt);

// Rotation helpers (exposed for direct use).

Mat3 rot_exp(const Vec3& theta);
Vec3 rot_log(const Mat3& R);

}  // namespace keta

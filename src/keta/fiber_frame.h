// Force-based (flexibility) nonlinear fiber beam-column — SPATIAL element with
// biaxial bending, optional corotational geometry and three continuation
// drivers (load / displacement / arc-length control) — plus the original
// PLANE X-Z API kept as a compatibility overlay on the same formulation.
//
// THE ELEMENT (Spacone, Ciampi & Filippou, "Mixed formulation of nonlinear
// beam finite element", Computers & Structures 58(1), 1996, 71-83; Neuenhofer &
// Filippou, "Evaluation of nonlinear frame finite-element models", J. Struct.
// Eng. 123(7), 1997, 958-966).
//
// A displacement-based element interpolates the DISPLACEMENT field (cubic
// Hermite), so its curvature field is linear and cannot follow the sharp
// curvature localisation of a plastic hinge without meshing. A force-based
// element interpolates the FORCE field instead, and for a prismatic member
// under end forces the internal forces are known EXACTLY from equilibrium —
// so one element per member carries the exact moment diagram at any level of
// material nonlinearity, and refinement is needed only for the integration of
// the deformations.
//
// SPATIAL basic (simply-supported) system, 6 modes:
//
//     v = (elongation, ty_i, ty_j, tz_i, tz_j, phi)   [m, rad x5]
//     q = (N, My_i, My_j, Mz_i, Mz_j, T)              [N, N.m x5]
//     ty^def = theta_y + (w_j - w_i)/L        (local z deflection w, ry = -w')
//     tz^def = theta_z - (v_j - v_i)/L        (local y deflection v, rz = +v')
//     phi    = rx_j - rx_i                    (Saint-Venant twist)
//
// with the exact force interpolation (no span load on a fiber bar — typed
// reject), section resultants in the `fiber.h` convention
// (eps = eps0 + kappa_y z - kappa_z y):
//
//     N(x)  = q0                    (constant)
//     My(x) = (1 - x/L) q1 - (x/L) q2         (linear)
//     Mz(x) = (1 - x/L) q3 - (x/L) q4         (linear)
//
// The (a, b) pair is contragredient — a6^t kb a6 reproduces keta's linear
// Euler-Bernoulli local stiffness bit for bit (elastic case to 1e-9), which
// fixes every sign. TORSION is ELASTIC and DECOUPLED: T = (G J / L) phi with J
// from the bar's `Section` — Saint-Venant only. Nonlinear/warping torsion is
// OUTSIDE the formulation domain of a point-fiber section (a fiber carries a
// normal stress only; restrained warping needs a sectorial coordinate field),
// as is nonlinear SHEAR (no shear law on a fiber): the element is
// Euler-Bernoulli + elastic Saint-Venant by construction.
//
// ELEMENT STATE DETERMINATION. Given a target basic deformation v: iterate on
// the 5 fiber forces q (torsion is closed-form):
//   1. s(x_k) = b(x_k) q at the Gauss-Lobatto points (5 by default: Lobatto
//      puts a point exactly on each END, where the hinges of a frame form);
//   2. section state determination d(x_k) = s^-1(x_k) — the 3x3 fiber Newton
//      (`HystFiberSection::solve_state`, biaxial target (N, My, Mz)) at the
//      COMMITTED material states;
//   3. compatibility residual  s_res = v - SUM_k w_k b_k^t d_k;
//   4. f_b = SUM_k w_k b_k^t f_s(x_k) b_k,  k_b = f_b^-1;
//   5. q += k_b s_res, repeat (with a fallback bisection towards the last
//      cleanly-evaluated q when a trial q leaves the section capacity).
// Because step 2 inverts the section exactly, the per-section residual
// deformations of the original algorithm vanish and only the element residual
// survives (Neuenhofer & Filippou 1997 §3). f_s is the inverse of the section
// tangent regularised as k_s + `tangent_floor` k_s0 — a numerical device on
// the TANGENT only; resisting forces always come from the laws.
//
// CONVERGENCE POLICY (element AND section). Both iterations converge on a
// scaled residual to their tolerance, and both fall back on a GRACEFUL FLOOR:
// when the residual stops improving (6 non-improving passes at the element,
// 10 at the section) the best state so far is accepted PROVIDED it already
// sits far below any structural relevance (1e3 x `tol_element` for the
// element, 1e-5 relative for the section). This is not a loosened tolerance:
// a deeply plastified section has a near-flat response, so the last decades
// of a force->deformation inversion crawl at a linear rate and stagnate at a
// solver-precision floor that no amount of iteration clears. States that are
// genuinely non-equilibrable stagnate orders of magnitude above the floor and
// still raise `NonConvergence` — which is exactly what a bracketing of the
// limit load relies on. The global Newton judges equilibrium in FORCE
// space independently, so an accepted floor state is still checked there.
//
// MATERIAL HISTORY. Sections are `HystFiberSection` (fiber_hyst.h): each
// Lobatto point of each element owns one committed `HystState` per fiber.
// Every evaluation inside a step is a TRIAL from the committed state; the
// driver COMMITS once per accepted step (backtracked / cut steps are reverted
// by simply not committing). Memoryless sections wrap `UniaxialLaw` and make
// the plane API a strict special case.
//
// GEOMETRIC NONLINEARITY (opt-in `corotational`): corotational formulation of
// the BASIC system (Filippou & Fenves, "Methods of analysis for earthquake-
// resistant structures", ch. 6 in Bozorgnia & Bertero 2004; de Souza, "Force-
// based finite element for large displacement inelastic analysis of frames",
// PhD, UC Berkeley 2000; Crisfield, "Non-linear Finite Element Analysis of
// Solids and Structures" vol. 1 §7). The element deformations are measured in
// the rotated CHORD frame: elongation from current node positions
// (cancellation-safe (ln^2 - L0^2)/(ln + L0)), the chord triad built from the
// chord direction + the mean of the two rotated nodal y axes, and the
// deformational rotations extracted by the rotation log map
// Theta_i = log(E^t R(theta_i) E0). Nodal rotations are total rotation
// vectors with additive updates (moderate rotations; same convention as
// `nonlinear.cpp`). Internal force f = B^t q with B = dv/du and tangent
// K = B^t k_b B + sum_a q_a d2v_a/du2, both evaluated by high-order central
// finite differences of the exact deformation extraction (the house pattern
// of the corotational element: FD error ~1e-10 on equilibrium, tangent
// accuracy affects convergence rate only). The P-delta INSIDE one element
// (bowing between the nodes) is NOT in the force interpolation, which stays
// linear in the chord frame: subdivide compressed members into 2-4 elements
// (standard force-based practice — the corotational chords then track the
// member curve and the assembly captures the internal P-delta; the Euler
// load is reached to ~1 % with 4 elements).
//
// DRIVERS. All three share the same assembled Newton core (full 6-DOF nodes),
// commit-on-acceptance and per-step trace:
//   * Load control — equal increments to lambda = 1 (limit point is a typed
//     `NonConvergence`, bracketing stays possible);
//   * Displacement control — one controlled DOF follows a piecewise-linear
//     target history (cyclic pushover), lambda is the unknown: the classic
//     bordered solve du = du_r + dlam du_F with dlam fixing the controlled
//     component. Limit points of the load are CROSSED; post-peak and cyclic
//     hysteresis are traced;
//   * Arc-length — Crisfield cylindrical continuation, the `nonlinear.cpp`
//     motif (predictor sign from the previous increment, deterministic root
//     selection maximizing dot(Du_new, Du), step cutting, adaptive radius).
//
// DOCUMENTED LIMITS (typed rejects, never silent):
//   * no span loads, releases, eccentricity, Winkler bed, taper, truss/cable
//     kind on a FIBER bar (elastic bars keep all of them in the geometrically
//     linear mode through `BarElement`);
//   * corotational mode: nodal loads only, and every bar (elastic too) must be
//     free of eccentricity / Winkler bed / taper (elastic bars then run
//     through the same corotational basic system with their exact basic
//     stiffness extracted from the condensed local matrix);
//   * no shear deformation in a fiber bar; elastic comparisons need Avy = Avz
//     = 0 sections to match keta's Timoshenko element bit for bit;
//   * torsion elastic Saint-Venant, warping and NL torsion out of domain (see
//     above);
//   * only plain fixed / diagonal-spring supports; no prescribed support
//     displacements;
//   * dense solve (Eigen `FullPivLU`): a small-frame driver.
//
// MODELLING NOTE — never a reject, but it decides whether a run is traceable.
// A force-based element INVERTS the section (force -> deformation), so a law
// with an exactly FLAT plateau makes that inversion ill-posed once a section
// saturates: d(kappa)/dM is unbounded and no finite-precision iteration can
// place the deformation. Load control still works (the limit load is reached
// and reported as a typed `NonConvergence`, which is what a bracketing of
// the limit load uses), but tracing a POST-PEAK branch needs a small hardening —
// `HystLaw::steel_kinematic(E, fy, b)` with b of order 1e-3..1e-2, the
// standard Steel01 practice. That hardening is far below the P-Delta
// softening slope, so the limit point and the descending branch survive.
//
// The PLANE API (`solve_fiber_frame`, options/results below) is the original
// v1 surface: same solver, plane gates kept (X-Z models, plane loads), results
// projected on (ux, uz, ry). It is a wrapper over the spatial core — one
// formulation, two surfaces.
#pragma once
#include <map>

#include "keta/fiber.h"
#include "keta/drying.h"
#include "keta/fiber_hyst.h"
#include "keta/statics.h"

namespace keta {

// Typed non-convergence: the load step could not be equilibrated (limit load
// reached, or a section driven outside its capacity polytope).
class NonConvergence : public Error {
 public:
  using Error::Error;
};

// ---------------- plane (X-Z) compatibility surface ----------------

struct FiberFrameOptions {
  // Bar id -> fiber section. Bars absent from the map stay linear elastic.
  std::map<Id, FiberSection> sections;
  int n_ip = 5;                  // Gauss-Lobatto points per fiber element (3..6)
  int load_steps = 10;           // equal load increments
  int max_global_iter = 40;
  int max_element_iter = 30;
  int max_backtracks = 12;
  double tol_force = 1e-9;       // ||residual|| / ||applied|| at a converged step
  double tol_element = 1e-11;    // relative energy norm of the element residual
                                 // (the section Newton runs one decade tighter)
  double tangent_floor = 1e-6;   // section tangent regularisation (see header)
};

// State of one fiber element at a Gauss-Lobatto point (plane view).
struct FiberStation {
  double x = 0;        // abscissa along the bar [m]
  double weight = 0;   // Lobatto weight * L [m]
  double N = 0;        // section axial force [N] (tension +)
  double My = 0;       // section moment INT sigma z dA [N.m]
  double eps0 = 0;     // section axial strain at the reference axis [-]
  double kappa = 0;    // section curvature kappa_y [1/m]
};

struct FiberBarState {
  Id bar = 0;
  bool fiber = false;
  Vec3 q = Vec3::Zero();  // basic forces (N, M_i, M_j) [N, N.m, N.m]
  std::vector<FiberStation> stations;
};

struct FiberFrameResult {
  // Per model node, in node order: (ux, uz, ry) — the plane DOFs [m, m, rad].
  std::vector<Vec3> node_disp;
  std::vector<FiberBarState> bars;  // parallel to model.bars
  int steps = 0;                    // load steps performed
  int iterations = 0;               // total global Newton iterations
  double lambda = 1.0;              // load factor reached (always 1 on success)
};

// Solve one load case of a plane frame whose bars may carry fiber sections.
// Throws `NonConvergence` when a step cannot be equilibrated, `Error` on a
// model that violates one of the documented limits above.
FiberFrameResult solve_fiber_frame(const Model& model, Id case_id,
                                   const FiberFrameOptions& opt);

// Gauss-Lobatto abscissae/weights on [0, 1] for n = 3..6 (exposed for direct use).
void gauss_lobatto01(int n, std::vector<double>& x, std::vector<double>& w);

// ---------------- spatial driver ----------------

struct FiberFrame3DOptions {
  // Bar id -> hysteretic fiber section (wrap a memoryless `FiberSection` with
  // `HystFiberSection::wrap`). Bars absent from the map stay linear elastic.
  std::map<Id, HystFiberSection> sections;
  int n_ip = 5;               // Gauss-Lobatto points per fiber element (3..6)
  bool corotational = false;  // corotational basic system (P-Delta and beyond)

  enum class Control { Load, Displacement, ArcLength };
  Control control = Control::Load;

  // Load control: equal increments to lambda = 1.
  int load_steps = 10;

  // Displacement control: DOF `ctrl_comp` (0..5: ux..rz) of node `ctrl_node`
  // follows the piecewise-linear history 0 -> targets[0] -> targets[1] -> ...
  // with `steps_per_leg` increments per leg; lambda scales the case loads.
  Id ctrl_node = 0;
  int ctrl_comp = 2;
  std::vector<double> ctrl_targets;
  int steps_per_leg = 10;

  // Arc length (Crisfield cylindrical, `nonlinear.cpp` motif).
  double dl0 = 0.01;          // initial arc radius (displacement norm) [m-ish]
  int max_steps = 60;
  int max_arc_iterations = 20;
  int max_cuts = 8;
  double stop_lambda = 0;     // 0 = no lambda stop
  double stop_disp_inf = 0;   // 0 = no displacement stop

  // Shared Newton / element knobs (same semantics as the plane options).
  int max_global_iter = 40;
  int max_element_iter = 30;
  int max_backtracks = 12;
  double tol_force = 1e-9;
  double tol_element = 1e-11; // element energy metric; section runs 10x tighter
  double tangent_floor = 1e-6;

  // --- appended: STAGED construction (material nonlinearity across phases) ---
  // Non-empty: the case_id argument is ignored and the stages run in order
  // under load control. Stage k activates bars and ADDS the loads of its case
  // (they stay applied) in `load_steps` increments, from the committed
  // fiber states of the previous stages (history carried: a fiber that
  // yielded stays yielded). A node is in place once it carries an active bar
  // or has its three translations held; a new node is placed by rigid continuation of its
  // neighbour (small rotations: u_j = u_i + theta_i x (X_j - X_i), theta_j =
  // theta_i), and EVERY activated bar is stress-free at activation — its
  // deformation is counted from its state then (a closure cast in place on
  // the deformed structure). A SelfWeight load in a stage case weighs the
  // bars activated by that stage only. Nodal loads on nodes not in place
  // are rejected.
  // Appended: CREEP between phases — after its loads, a stage may last
  // `duration` days in `time_steps` steps under constant loads: the fibers
  // whose law creeps (HystLaw::creep_a, keta/fiber_hyst.h: a Kelvin chain,
  // linear or EN 1992-1-1 §3.1.4(4) nonlinear) take their exact
  // exponential-algorithm step, solved implicitly with equilibrium (Newton
  // on the consistent tangent); the loading itself is instantaneous. Their
  // SHRINKAGE (HystLaw::shrink, the age counted from the bar's activation)
  // is imposed at the end of each time step.
  struct Stage {
    std::vector<Id> activate;
    Id case_id = 0;  // 0: no load
    int load_steps = 1;
    double duration = 0;  // [days] of creep after the loads
    int time_steps = 1;
    // (appended) the age [days] of the activated bars' concrete at their
    // activation (precast: > 0); their shrinkage counts from then
    double activation_age = 0;
  };
  std::vector<Stage> stages;
  // Appended: DRYING of the sections (bar id -> its rectangular fiber grid
  // and the diffusion parameters, keta/drying.h): from its activation the
  // bar's humidity field diffuses over each time step of the stages, and
  // each concrete fiber of the grid (the first ny x nz fibers of the
  // section; rebars do not shrink) takes the shrinkage of its own humidity
  // — the drying gradient's eigenstresses, the curvature of an asymmetric
  // exposure. Not together with a law's shrinkage table on those fibers.
  std::map<Id, DryingSection> drying;
};

// State of one fiber element at a Gauss-Lobatto point (spatial view).
struct FiberStation3D {
  double x = 0, weight = 0;             // abscissa / Lobatto weight * L [m]
  double N = 0, My = 0, Mz = 0;         // section resultants (fiber convention)
  double eps0 = 0, kappa_y = 0, kappa_z = 0;
};

struct FiberBarState3D {
  Id bar = 0;
  bool fiber = false;
  Vec6 q = Vec6::Zero();  // basic forces (N, My_i, My_j, Mz_i, Mz_j, T)
  std::vector<FiberStation3D> stations;
};

// One accepted equilibrium point of the continuation.
struct FiberTracePoint {
  double lambda = 0;
  double u_ctrl = 0;  // controlled-DOF value (Displacement control, else 0)
  std::vector<Vec6> node_disp;
};

struct FiberFrame3DResult {
  std::vector<Vec6> node_disp;          // final state, per model node
  std::vector<FiberBarState3D> bars;    // parallel to model.bars
  std::vector<FiberTracePoint> trace;   // every accepted step, in order
  int steps = 0;
  int iterations = 0;
  double lambda = 0;                    // load factor at the final state
  std::vector<std::size_t> stage_end;   // staged run: trace index closing each stage
  std::vector<double> trace_time;       // staged run: the time of each trace point [days from the start]
  std::map<Id, std::vector<double>> humidity;  // drying bars: the final humidity per grid cell
};

FiberFrame3DResult solve_fiber_frame3d(const Model& model, Id case_id,
                                       const FiberFrame3DOptions& opt);

}  // namespace keta

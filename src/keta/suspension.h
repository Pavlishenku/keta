// Suspension bridges: main-cable FORM-FINDING under hanger loads (vertical OR
// inclined), model assembly of the found geometry (optionally with back-stay
// side spans and their saddle regime), tower-flexibility coupling, and the
// live-load response driver.
//
// ---------------------------------------------------------------------------
// Physics
// ---------------------------------------------------------------------------
// The main cable is the FUNICULAR of the loads the hangers deliver. Between two
// consecutive hangers the cable is an elastic catenary of constant horizontal
// tension; at a hanger node the internal tension components jump by the hanger
// force VECTOR. Writing (H_k, V_k) for the horizontal / vertical components of
// the tension carried by segment k (the force the right part of the cable
// exerts on the left part), node equilibrium at hanger k gives
//
//     V_start(k+1) = V_end(k) + P_k                                        (2)
//     H_{k+1}      = H_k - P_k * dx_k / (-dz_k)                            (2')
//
// with P_k >= 0 the DOWNWARD load carried by the hanger and (dx_k, dz_k) the
// in-plane vector from the cable node to the hanger's deck attachment point
// (dz_k < 0: the deck hangs below). A vertical hanger has dx_k = 0 and (2')
// collapses to the classical "H = const along the whole cable"; an inclined
// hanger (Severn-type systems; Gimsing & Georgakis, *Cable Supported Bridges*,
// 3rd ed., Wiley 2012, ch. 2 & 5) trades horizontal tension across the node
// because the hanger force is ALIGNED WITH ITS AXIS — that is its definition.
// The hanger axial force is then P_k * |(dx_k, dz_k)| / (-dz_k).
//
// In the weightless vertical-hanger limit the funicular is the classical
// "bending-moment diagram of the equivalent simply supported beam divided by
// H": for a UDL q over a span L with mid-sag f this is H = q L^2 / (8 f),
// reproduced to 0.5 %.
//
// Segment kinematics use the EXACT elastic catenary (Irvine, *Cable
// Structures*, MIT Press 1981, §1.3-1.4 — the same closed forms as `cable.h`),
// rewritten here in a CANCELLATION-FREE form so that the w -> 0 chord limit is
// exact instead of losing ~log10(V/(w L0)) digits to the difference of two
// asinh's:
//
//     lx(s) = (s B / H) * sinhc^-1(A) + H s / EA ,   B = T1 - V1(V1+V2)/(T1+T2)
//     lz(s) = s (V1+V2)/(T1+T2) + (s/EA)(V1 + w s / 2)                    (3)
//     A     = w s B / H^2 ,   sinhc^-1(A) = asinh(A)/A  (-> 1 as A -> 0)
//
// with V2 = V1 + w s, T1 = hypot(H, V1), T2 = hypot(H, V2). Setting w = 0 gives
// lx = H s / T1 + H s / EA, the elastic straight chord, exactly.
//
// ---------------------------------------------------------------------------
// Form-finding algorithm (documented choice)
// ---------------------------------------------------------------------------
// Two unknowns only: (H, V0) — the tension components of the FIRST segment at
// the left saddle. Everything else follows by a forward MARCH:
//
//   for each segment k, solve lx(s) = dx_k for its unstressed length L0_k
//     (1-D, strictly increasing: dlx/ds = H/T2 + H/EA > 0, safeguarded Newton),
//   accumulate z by lz(L0_k), then apply the node jumps (2) + (2').
//
// INCLINED HANGERS, the robustness choice: the CABLE attachment abscissas are
// imposed (`hanger_x`) and the deck attachment points are known (`hanger_deck`)
// — so when the march reaches hanger k it already knows the cable node position
// (x imposed, z just marched) and the jump direction is pure geometry. No inner
// iteration exists at all; the only implicit dependence (the marched z feeds
// the jump) is closed by the SAME outer 2x2 Newton that closes the geometry,
// which therefore still sees a smooth 2-D map. The alternative — imposing deck
// abscissas and solving each cable abscissa by an inner fixed point — was
// rejected: it nests an iteration inside every residual evaluation and its
// convergence degrades exactly in the strongly inclined configurations the
// feature exists for. A hanger whose deck point is NOT below the cable node
// would be in compression: typed error, never a sign flip.
//
// The two residuals closed on by the 2x2 Newton (central-difference Jacobian,
// step-limited so H stays positive) are
//
//   r1 = z(right saddle) - z_right                                (geometry)
//   r2 = sag(H, V0) - sag_target       or      H - H_target        (the target)
//
// This is the classic "two-level" cable iteration flattened into one Newton:
// the inner level (segment lengths) is solved to machine precision inside every
// residual evaluation. Deterministic — no random restarts, fixed iteration
// budget, every failure a typed `Error`. `target_H` targets the FIRST segment
// (equal to every other one when all hangers are vertical).
//
// `sag` is measured at the HORIZONTAL MIDPOINT of the two saddles, from the
// saddle chord down to the cable (positive = the cable hangs below the chord).
// Unequal saddle heights are fully supported.
//
// ---------------------------------------------------------------------------
// Back stays and saddle regimes (side spans)
// ---------------------------------------------------------------------------
// The cable may continue past a saddle to a back anchorage (`BackstaySpec`).
// Two saddle regimes bound the physics (Gimsing & Georgakis, op. cit., ch. 5 —
// tower / cable interaction):
//
//   * CLAMPED — the cable is fixed to the saddle. The two sides carry
//     independent horizontal tensions; the difference H_main - H_back descends
//     into the tower as shear/moment (that is what the tower-flexibility
//     coupling below equilibrates). The back-stay is erected to a
//     horizontal-tension target `back_H` (default: matched to the adjacent
//     main-span segment, i.e. a balanced saddle).
//   * FRICTIONLESS — the cable can slide: the TENSION MAGNITUDE is equal on
//     the two sides of the saddle (rope-on-pulley equilibrium). The back-stay
//     unstressed length is solved (safeguarded 1-D Newton/bisection on its
//     horizontal tension, residual <= 1e-12 relative) so that
//     |T_main(saddle)| = |T_back(saddle)|.
//
// The regime is a FORM-FINDING (dead-state) condition: the assembled model
// carries both cables with their found L0 fixed, so live load may develop a
// tension difference across a frictionless saddle (no slip under live —
// documented domain edge, see below).
//
// ---------------------------------------------------------------------------
// Tower flexibility (fixed-point coupling)
// ---------------------------------------------------------------------------
// `build_suspension` alone puts the saddles on the given nodes and form-finds
// with those positions fixed. `build_suspension_on_towers` couples the
// form-finding with modelled towers: the saddle nodes are the TOPS of tower
// bars already in the model, and the driver iterates
//
//   form-find with the saddles at the current (displaced) positions
//   -> assemble on a trial copy -> solve the permanent case
//   -> read the tower-top displacements -> move the saddles -> repeat
//
// until the saddle position update is below `saddle_tol` (default 1e-8 m,
// typed non-convergence guard). WHY THE FIXED POINT CONVERGES: each trial
// solve is itself an equilibrium of tower + cables, so the update map
// contracts like k_cable / (k_cable + k_tower) — always < 1, but close to 1
// when a taut back-stay ties the top elastically to ground (measured ~0.8).
// The update is therefore ACCELERATED by Aitken / Irons-Tuck dynamic
// relaxation (Irons & Tuck, IJNME 1969 — the standard fixed-point accelerator
// of partitioned couplings), which restores convergence in a handful of
// evaluations; theta is clamped and deterministic. The hanger cable-clamp
// abscissas are re-expressed in the displaced chord frame each iteration so
// the attachment points stay put in SPACE (a tower lean must not drag the
// clamps sideways off their deck nodes). The cable bars connect to the
// UNDEFORMED tower-top node while their L0 comes from the DISPLACED-saddle
// form-finding: the FEM solve then reproduces exactly the displaced funicular
// — that consistency is the fixed point.
//
// ---------------------------------------------------------------------------
// Assembly and the acid test
// ---------------------------------------------------------------------------
// `build_suspension` writes the found geometry into a `Model`: cable nodes at
// the funicular vertices, `BarKind::Cable` bars carrying the found L0 in
// `CableProps::unstressed_length` (back stays included), and one hanger bar
// per hanger down to an existing deck node. Hangers are `BarKind::Truss` by
// default — a straight axial link, which is exactly the form-finding's model
// of a hanger; a `Cable` hanger is rejected (`eval_cable` has no vertical
// chord, and a catenary hanger would not be the assumed straight axial link).
//
// HANGER LENGTHS. A hanger carrying T is T L0 / EA longer than it was made:
// laid unstressed at its geometric length, it lets the deck drop by that
// stretch under the dead load. The builder reports the unstrained length of
// every hanger (`SuspensionModel::hanger_l0`, `hanger_unstrained_length`) and,
// given the case that holds the permanent state (`SuspensionSpec::
// hanger_fit_case`), writes the lack of fit there as a Prestrain load on each
// hanger bar — the kernel's statement of a bar fabricated short, the one
// `tune_cables` uses: the hangers then carry their form-found tension at the
// form-found geometry and the deck hangs at its drawn level — exactly when
// the hanger loads of the form are the reactions the deck hands them (a deck
// hinged at every hanger); a continuous deck whose reactions differ keeps
// the deflection of that redistribution. A live case
// combined with that permanent case (`suspension_live_response`) keeps it.
// Each hanger has its own length: grouping them by fabrication lengths is not
// done here.
//
// The acid test of a form-finding is that the DEAD-LOAD solve reproduces it:
// `solve_cable_static` on the assembled model under the very load case whose
// tributary reactions fed the form-finding must return the same H, the same
// hanger tensions, and a deck bending moment that is NEGLIGIBLE compared with
// the simple-beam q L^2 / 8 — the deck is hung, not spanned. The dead-load
// solve reproduces the form to 0.20 % with vertical hangers, 0.4 % at 20
// degrees of hanger inclination, and the hanger tensions to 1 %.
//
// ---------------------------------------------------------------------------
// Units and design domain
// ---------------------------------------------------------------------------
// SI throughout: m, N, Pa, kg, rad. Z is up.
//   * ONE CABLE PLANE PER CALL (design choice, not a deferral): the saddle
//     chord defines a vertical plane; every hanger deck point, anchorage and
//     saddle sits on it. A real bridge has two (or more) cable planes: call
//     the builder once per plane with that plane's tributary loads.
//   * hangers may be vertical or inclined in the cable plane; a hanger whose
//     deck point is not strictly below its cable node is a typed error
//     (compression), as is one out of the plane;
//   * saddle regimes cover the two FRICTION EXTREMES: clamped (no slip) and
//     frictionless (tension equality). A finite-friction saddle
//     (|ln(T1/T2)| <= mu * wrap angle) lies between the two and is out of the
//     modelled domain — bracket it with the two extremes;
//   * the frictionless regime is enforced at the DEAD state; the assembled
//     model does not re-slide under live load;
//   * `solve_cable_static` (the driver used for the response) floors a zero
//     cable weight internally, so a rho = 0 main cable form-found here with the
//     exact w = 0 branch will drift by that floor.
#pragma once
#include "keta/cable.h"
#include "keta/cable_tuning.h"

namespace keta {

// ---------------- form-finding ----------------

// Deck-side attachment point of a hanger, CABLE-PLANE coordinates:
// x = horizontal offset from the left saddle along the saddle chord's plan
// direction [m], z = ABSOLUTE elevation [m] (same datum as the saddles).
struct HangerAnchor {
  double x = 0, z = 0;
};

struct FormFindSpec {
  Vec3 left = Vec3::Zero();          // left tower saddle point [m]
  Vec3 right = Vec3::Zero();         // right tower saddle point [m]
  std::vector<double> hanger_x;      // CABLE attachment abscissa from `left` [m], strictly increasing
  std::vector<double> hanger_load;   // DOWNWARD load carried by each hanger [N], >= 0
                                     // (the VERTICAL component for an inclined hanger)
  double EA = 0;                     // main-cable axial rigidity [N] (> 0)
  double w = 0;                      // main-cable weight per unit unstressed length [N/m] (>= 0)
  std::optional<double> target_sag;  // mid-span sag below the saddle chord [m] (> 0)
  std::optional<double> target_H;    // ... or the first-segment horizontal tension [N] (> 0). Exactly one.
  int max_iterations = 80;
  double tol_rel = 1e-12;
  // Inclined hangers: deck-side attachment point per hanger (size n), driving
  // the jump direction (2'). EMPTY = every hanger vertical (bit-identical to
  // the classical constant-H march).
  std::vector<HangerAnchor> hanger_deck;
};

struct FormFound {
  double H = 0;                        // horizontal tension of the FIRST segment [N]
                                       // (constant along the cable when all hangers are vertical)
  std::vector<Vec3> nodes;             // saddle_L, hanger nodes..., saddle_R (n + 2)
  std::vector<double> segment_l0;      // unstressed length per segment [m] (n + 1)
  std::vector<double> segment_v1;      // vertical tension component at each segment start [N]
  std::vector<double> segment_tmax;    // max tension per segment [N] (n + 1)
  std::vector<double> hanger_tension;  // hanger AXIAL force at dead load [N] (n)
  double sag = 0;                      // realised mid-span sag [m]
  Vec3 saddle_force_left = Vec3::Zero();   // force the cable exerts ON the left saddle [N]
  Vec3 saddle_force_right = Vec3::Zero();  // ... and on the right saddle
  double total_l0 = 0;                 // sum of the segment unstressed lengths [m]
  double EA = 0, w = 0;                // the cable properties actually used
  int iterations = 0;
  std::vector<double> segment_h;       // horizontal tension per segment [N] (n + 1)
  std::vector<Vec3> hanger_force;      // force each hanger exerts ON its cable node [N] (n)
};

// Solve the funicular. Throws a typed `Error` on a bad spec, an unreachable
// target (slack / unstable configuration), a hanger in compression, or
// non-convergence.
FormFound form_find_cable(const FormFindSpec& spec);

// Cable elevation of a form-found shape at a horizontal offset from the left
// saddle (exact catenary interpolation INSIDE the segment, not a chord).
double form_found_elevation(const FormFound& ff, double x);

// ---------------- back stays (side spans) ----------------

// Saddle regime of a back-stay: see the header comment. The two friction
// extremes; a finite-friction saddle is out of the modelled domain.
enum class SaddleRegime { Clamped, Frictionless };

struct BackstaySpec {
  Vec3 anchor = Vec3::Zero();     // back-anchorage point [m], in the cable plane,
                                  // on the FAR side of the saddle (away from the span)
  SaddleRegime regime = SaddleRegime::Clamped;
  std::optional<double> back_H;   // Clamped only: erection horizontal-tension target [N];
                                  // empty = matched to the adjacent main-span segment.
                                  // Setting it in the Frictionless regime is a typed error
                                  // (the tension-equality condition determines it).
  Id anchor_node = 0;             // existing anchorage node; 0 = create one, fully fixed
};

// Solved back-stay: one elastic catenary from the saddle to the anchorage.
struct BackstayFound {
  double H = 0;                        // horizontal tension [N]
  double l0 = 0;                       // unstressed length [m]
  double T_saddle = 0, T_anchor = 0;   // tension at the saddle / anchorage end [N]
  Vec3 saddle_force = Vec3::Zero();    // force the back-stay exerts ON the saddle node [N]
  Vec3 anchor_force = Vec3::Zero();    // ... and on the anchorage
};

// ---------------- model assembly ----------------

struct SuspensionSpec {
  FormFindSpec form;                  // EA / w may be left at 0: taken from the cable section.
                                      // `form.hanger_deck` is DERIVED from `hanger_deck_nodes`
                                      // by the builder (any user-filled value is overwritten).
  Id cable_material = 0, cable_section = 0;
  Id hanger_material = 0, hanger_section = 0;
  BarKind hanger_kind = BarKind::Truss;  // Truss (default) or Frame; Cable is rejected
  std::vector<Id> hanger_deck_nodes;  // existing deck node under each hanger (size = hangers)
  Id left_node = 0, right_node = 0;   // existing saddle nodes; 0 = create one
  bool fix_saddles = true;            // fully fix a saddle node this call creates
  std::vector<Id> deck_bars;          // the deck chain — reporting only
  double plane_tol = 1e-6;            // relative out-of-plane tolerance of the
                                      // "hanger / anchorage sits in the cable plane" check
  double saddle_pos_tol = 1e-6;       // relative tolerance of the "existing saddle node is at
                                      // the saddle point" check (the tower coupling loosens it:
                                      // the node sits at the UNDEFORMED position while the
                                      // saddle point is the displaced one)
  std::optional<BackstaySpec> back_left, back_right;  // side spans (none by default)
  // --- appended: hangers at their unstrained length ---
  // The case that receives the lack of fit of every hanger (a Prestrain load
  // -(L - L0) / L0 on its bar, see hanger_unstrained_length): the hangers
  // then carry their form-found tension AT the form-found geometry, and the
  // deck hangs at its drawn level instead of dropping by the stretch T L / EA
  // of hangers laid unstressed. 0 = none (build_suspension: the hangers are
  // bars of their geometric length; build_suspension_on_towers and solve_uelm
  // take their dead case).
  Id hanger_fit_case = 0;
  // Weight per unit length that the fit case applies to a hanger [N/m] (rho A g
  // under a gravity load; 0 = a weightless link). The form's hanger load is
  // the load at the CABLE node, so it then counts the hanger itself.
  double hanger_weight = 0;
};

struct SuspensionModel {
  FormFound form;
  std::vector<Id> cable_nodes;  // saddle_L, hanger nodes..., saddle_R
  std::vector<Id> cable_bars;   // main-cable segments, left to right
  std::vector<Id> hanger_bars;  // one per hanger, cable node -> deck node
  std::vector<Id> deck_nodes;   // echoed from the spec
  std::vector<Id> deck_bars;    // echoed from the spec
  // Back stays (empty / nullopt without side spans).
  std::optional<BackstayFound> back_left, back_right;
  std::vector<Id> backstay_bars;  // left then right (present ones only)
  std::vector<Id> anchor_nodes;   // matching anchorage nodes
  // --- appended ---
  std::vector<double> hanger_l0;  // unstrained length per hanger [m] (hanger_unstrained_length
                                  // at the form-found tension, the spec's hanger_weight)
};

// Unstrained length of a straight hanger of stretched length `length` [m] and
// axial rigidity EA [N] whose tension is `tension` [N] at one end and varies
// by `w_axial` [N/m] per unit unstrained length towards the other end (Hooke
// on the unstrained length, as the elastic catenary of keta/cable.h):
//
//     length = L0 + (tension L0 + w_axial L0^2 / 2) / EA
//     L0 = 2 length / (a + sqrt(a^2 + 2 w_axial length / EA)),  a = 1 + tension / EA
//
// From the FOOT of a vertical hanger of weight w per unit length, w_axial = w
// (the tension grows upwards); from its TOP, w_axial = -w. An inclined hanger
// is a straight link whose axis carries w drop / length (the sag of its own
// weight is not modelled). Typed error: length, EA not > 0, a tension that
// turns negative along the hanger.
double hanger_unstrained_length(double length, double EA, double tension, double w_axial = 0);

// Form-find, then write the cable + hangers (+ back stays) into `model`. The
// deck chain must already exist; every `hanger_deck_nodes[i]` must sit in the
// cable plane, strictly BELOW the hanger's cable node (vertical or inclined —
// the hanger direction is taken from the geometry). Throws a typed `Error`
// otherwise.
SuspensionModel build_suspension(Model& model, const SuspensionSpec& spec);

// ---------------- tower-flexibility coupling ----------------

struct SuspensionFlexOptions {
  int max_iterations = 40;
  double saddle_tol = 1e-8;  // fixed-point convergence: max saddle position update [m]
};

struct SuspensionFlexResult {
  SuspensionModel sm;                        // assembled (in the caller's model) at the
                                             // CONVERGED saddle positions
  Vec3 saddle_disp_left = Vec3::Zero();      // converged tower-top displacement under the
  Vec3 saddle_disp_right = Vec3::Zero();     // permanent case [m]
  double saddle_gap = 0;                     // last fixed-point update [m] (<= saddle_tol)
  int iterations = 0;
};

// Fixed-point coupling of the form-finding with modelled towers (see the
// header comment). `spec.left_node` / `spec.right_node` MUST name existing
// tower-top nodes (held by tower bars, not by supports) and `spec.form.left` /
// `spec.form.right` their undeformed positions; `dead_case` is the permanent
// load case (typically the deck dead load — the cables add their own weight).
// Trial assemblies run on private copies; the caller's `model` receives ONLY
// the final converged assembly. Typed `Error` on non-convergence.
SuspensionFlexResult build_suspension_on_towers(Model& model, const SuspensionSpec& spec,
                                                Id dead_case,
                                                const SuspensionFlexOptions& opt = {});

// ---------------- response ----------------

struct SuspensionResponse {
  double cable_H = 0;                   // horizontal tension of the first segment [N]
  double cable_Tmax = 0;                // largest tension anywhere in the main cable [N]
  std::vector<double> segment_H;        // per main-cable segment [N]
  std::vector<double> hanger_tension;   // per hanger [N], tension positive
  double deck_my_max = 0;               // extreme station My over the deck chain [N·m]
  double deck_my_min = 0;
  double deck_dz_min = 0;               // most negative vertical deck node displacement [m]
  double deck_dz_max = 0;
  CableStaticResults raw;               // the underlying coupled cable+frame solve
  std::vector<double> backstay_H;       // per back-stay bar (sm.backstay_bars order) [N]
  std::vector<double> backstay_Tmax;
};

// One coupled cable + frame solve of `case_id`, reduced to the suspension
// quantities. Results are NOT superposable (the cable is nonlinear).
SuspensionResponse suspension_response(const Model& model, const SuspensionModel& sm, Id case_id,
                                       const CableStaticOptions& opt = {});

struct SuspensionDelta {
  SuspensionResponse dead;              // dead case alone
  SuspensionResponse total;             // dead + live in ONE case (cables never superpose)
  std::vector<double> hanger_delta;     // total - dead per hanger [N]
  double cable_H_delta = 0;             // total - dead [N]
  double deck_my_max_delta = 0, deck_my_min_delta = 0;
};

// Live-load response: solves the dead case and a case holding dead + live
// loads TOGETHER (a nonlinear cable system cannot be superposed), and reports
// the deltas. The combined case is built in a private copy of the model, so
// `model` is untouched and every id stays valid.
SuspensionDelta suspension_live_response(const Model& model, const SuspensionModel& sm,
                                         Id dead_case, Id live_case,
                                         const CableStaticOptions& opt = {});

// ---------------- UELM two-stage initial state ----------------
// Stage 1: form_find_cable + build_suspension (the funicular of the
// tributary hanger loads, segment L0 inverted analytically). Stage 2: the
// coupled cable + frame equilibrium of `dead_case` with those L0 frozen —
// the deck's stiffness, the hangers' stretch and the cable weight move the
// cable away from the stage-1 profile. The correction: the segment L0 are
// re-solved until the equilibrium sits on the stage-1 profile (z of every
// hanger node of the main cable, x of the sag node: a square set), by the
// Gauss-Newton of tune_cables on the measured deviation (a unit-gain fixed
// point on the stage-1 inputs crawls when the deck side is stiff). The model
// receives the stage-1 assembly with the corrected unstrained lengths.
struct UelmOptions {
  CableTuningOptions tuning;  // tol, iterations, statics of the correction
};
struct UelmResult {
  SuspensionModel sm;                // the stage-1 assembly (form, node and bar ids)
  std::vector<double> l0;            // corrected unstrained length per main-cable segment [m]
  double deviation_stage1 = 0;       // max |z, x deviation| of the stage-2 equilibrium at stage-1 L0 [m]
  double deviation = 0;              // ... after the correction [m]
  int iterations = 0;                // Gauss-Newton steps of the correction
  CaseResult state;                  // the corrected equilibrium of dead_case
};
UelmResult solve_uelm(Model& model, const SuspensionSpec& spec, Id dead_case, const UelmOptions& opt = {});

// ---------------- deck camber ----------------
// Fabrication camber: the vertical offsets of chosen deck nodes at
// fabrication such that the deformed profile under the dead case meets the
// target elevations. Unit-gain fixed point c <- c - (z_deformed - z_target),
// the structure re-fabricated (node z = model z + c) and re-solved every
// pass — the coupled cable statics when the model carries catenary cables,
// the linear statics otherwise; the bars at a cambered node follow the new
// geometry (fabrication, not displacement). For a frame the deflection
// depends on the camber only through the geometry (second order in the
// slopes): two or three passes; for a hung deck the hanger / cable coupling
// slows it as the hangers stiffen (typed error past max_iterations, with the
// deviation trace). A target on a node whose vertical DOF is blocked is a
// typed error (no camber can move it).
struct CamberOptions {
  int max_iterations = 30;
  double tol_rel = 1e-9;       // max deviation / bounding-box diagonal of the model
  CableStaticOptions statics;  // of the cable-carrying passes
};
struct CamberResults {
  std::vector<double> camber;      // per target: fabrication offset (+ = fabricated high) [m]
  std::vector<double> deviation;   // per target: deformed z - target [m]
  double max_deviation = 0;
  double span_scale = 0;           // the bounding-box diagonal the tolerance is taken against
  int iterations = 0;              // camber updates
  std::vector<double> trace;       // max deviation after each solve (iterations + 1)
  Model cambered;                  // the model at its fabrication geometry
  CaseResult state;                // the dead-case equilibrium of `cambered`
};
// targets: (node id, target ABSOLUTE elevation z after the dead load) [m].
CamberResults solve_deck_camber(const Model& model, Id dead_case, const std::vector<std::pair<Id, double>>& targets,
                                const CamberOptions& opt = {});

}  // namespace keta

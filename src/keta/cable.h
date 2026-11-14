// Exact elastic catenary cable (Irvine): closed-form projection integrals,
// analytic 2x2 flexibility, Newton on (H, V1), 3-D lift with the out-of-plane
// pendulum stiffness, unstrained-length resolution from a length or a
// horizontal-tension target, and the coupled cable + linear-frame static
// driver (Newton with backtracking line search on the true residual).
//
// 2-D kernel state: H > 0 horizontal tension, V1 vertical end-1 tension
// component, unstrained length L0, weight w per unit UNSTRAINED length,
// axial rigidity EA. Anchor forces on the CABLE: F1 = (H, V1) at end 1,
// F2 = (-H, -V2) at end 2, V2 = V1 + w L0; the cable pulls the structure
// with -F1 / -F2. Identity: F1 + F2 = (0, w L0) (the weight hangs).
#pragma once
#include "keta/statics.h"

namespace keta {

struct Catenary2D {
  double H = 0, V1 = 0;   // solved state
  double L0 = 0, EA = 0, w = 0;
  double V2() const { return V1 + w * L0; }
  double T1() const { return std::hypot(H, V1); }
  double T2() const { return std::hypot(H, V2()); }
};

// Solve (H, V1) so the elastic catenary spans (dx > 0, dz). Throws on
// non-convergence or invalid input.
Catenary2D solve_catenary(double dx, double dz, double L0, double EA, double w);

// Projection integrals and their Jacobian at a state (exposed for direct use).
void catenary_projections(const Catenary2D& c, double& lx, double& lz);
Mat3 catenary_flexibility22(const Catenary2D& c);  // top-left 2x2 used

// Unstrained length such that the solved horizontal tension equals H_target.
double catenary_l0_for_tension(double dx, double dz, double EA, double w, double H_target);
// Unstrained length such that the tension MAGNITUDE at end 2 (s = L0) equals
// T_target, on the TAUT branch: the end tension first falls with slack then
// rises again with self-weight, so the bracket starts just inside a 4.5 %
// strain (0.9 x the 5 % elastic-validity cap) and grows by 5 % steps;
// a target below the tension minimum (no taut root) or above the cap throws.
double catenary_l0_for_end_tension(double dx, double dz, double EA, double w, double T_target);
// The same with the tension measured at end 1 (s = 0, end = 1) or end 2
// (end = 2, the function above).
double catenary_l0_for_end_tension(double dx, double dz, double EA, double w, double T_target, int end);



// 3-D cable element evaluation at given anchor positions.
struct CableEval {
  Vec3 f1 = Vec3::Zero(), f2 = Vec3::Zero();  // forces the cable exerts ON the nodes
  Mat3 k3 = Mat3::Zero();                     // translational pair tangent (node1 diagonal block)
  double H = 0, T_max = 0;
};
CableEval eval_cable(const Vec3& p1, const Vec3& p2, double L0, double EA, double w);

// Coupled cable + frame statics. Cables carry their own weight always; the
// load case supplies the other actions. Results are NOT superposable.
struct CableStaticOptions {
  int max_iterations = 60;
  double tol_rel = 1e-10;
  int substeps = 1;
  int stations_per_bar = 11;
  // Cooperative cancel + progress (keta/hooks.h): polled every Newton step,
  // "substep k/n" reported per load substep. Null: unhooked.
  const SolveHooks* hooks = nullptr;
};
struct CableStaticResults {
  CaseResult result;                       // frame recovery (cables excluded)
  std::vector<double> cable_H;             // per model bar (0 for non-cables)
  std::vector<double> cable_Tmax;
  int iterations = 0;                      // Newton steps over all substeps
};
CableStaticResults solve_cable_static(const Model& model, Id case_id,
                                      const CableStaticOptions& opt = {});

// Resolve every cable's unstrained length at the undeformed anchors, from
// CableProps in priority order: explicit unstressed_length, end_tension
// (tension magnitude at end j), pretension (horizontal tension).
std::vector<double> resolve_cable_l0(const Model& model);

// A cable bar's axial rigidity and weight per unstrained length: the
// CableProps overrides (ea > 0, w set) else E*A and rho*A*g of the bar.
std::pair<double, double> cable_ea_w(const Model& model, const Bar& bar);

// ---------------- Ernst equivalent stiffness ----------------
// Ernst's equivalent modulus of a cable of specific weight gamma = w / A
// [N/m^3], horizontal projection l [m] and tension stress sigma [Pa]:
//   E_eq = E / (1 + (gamma l)^2 E / (12 sigma^3)).
double ernst_modulus(double E, double gamma, double l, double sigma);
// The linear chord stiffness that replaces an elastic catenary of chord
// (dx > 0 horizontal, dz vertical) under horizontal tension H: the elastic
// bar in series with the sag spring,
//   1 / k = L0 / EA + (w dx)^2 Lc / (12 T^3),  T = H Lc / dx (mid-chord tension),
// (the classical Ernst form, l = dx, sigma = T / A; the elastic term on the
// unstrained length L0 of the stretched cable). Its error on the exact
// catenary chord stiffness is about chi^2 / 5 with the sag parameter
// chi = w Lc cos^2(alpha) / H (cos alpha = dx / Lc); `note` is set beyond
// kErnstSagNote (error above ~2e-3: use the exact catenary).
// The sag term must not be written w^2 Lc^3 / (12 H^3) (horizontal tension,
// chord length): exact for a level chord, it is a factor 1 / cos^5(alpha)
// too compliant on an inclined one (-18 % stiffness at 30 deg, -34 % at
// 60 deg, chi = 0.08).
inline constexpr double kErnstSagNote = 0.1;
struct ErnstStiffness {
  double k = 0;          // combined chord stiffness [N/m]
  double k_elastic = 0;  // EA / L0
  double k_sag = 0;      // 12 T^3 / ((w dx)^2 Lc) (infinite for w = 0)
  double chi = 0;        // sag parameter w Lc cos^2(alpha) / H
  bool note = false;     // chi > kErnstSagNote
};
ErnstStiffness ernst_chord_stiffness(double EA, double L0, double dx, double dz, double w, double H);

// The model with every catenary cable replaced by the truss of its Ernst
// stiffness at a cable-static state (the linear substitute the linear
// analyses — statics, modes, spectra — can take): per cable bar a material of
// modulus k Lc_model / A and of density keeping the cable mass w L0 / g;
// `axial` holds the tension of the state per bar (cables: T at mid-chord;
// frame bars: N of the state), the prestress of a modal analysis about it
// (with it, the string stiffness T / Lc of the stays). The trusses lie on
// the model's geometry: a state whose anchors moved by u turns the chords by
// ~u / Lc, an error of that order (the exact linearization about any state
// is tangent_of_cable_state, keta/tangent.h).
struct ErnstSubstitute {
  Model model;
  std::vector<double> axial;
  std::vector<ErnstStiffness> cables;  // per model bar (default for non-cables)
};
ErnstSubstitute ernst_substitute(const Model& model, const CableStaticResults& state);

}  // namespace keta

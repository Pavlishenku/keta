// History-carrying (hysteretic) uniaxial laws + the fiber section that drives
// them — the state layer the memoryless `fiber.h` kernel deliberately lacks.
//
// STATE PATTERN (the standard one of nonlinear FE codes): a law is a PURE
// PARAMETER SET; the mutable material memory lives in a separate `HystState`
// POD owned by the caller (one per fiber per integration point). Within a
// global step every evaluation is a TRIAL: `response(state, eps)` maps the
// COMMITTED state + a trial strain to (stress, tangent) without mutating
// anything, so a Newton iteration or a rejected/backtracked step needs no
// undo — reverting IS simply not committing. When a step is accepted the
// caller calls `update(state, eps)` (or `HystFiberSection::commit`) once to
// roll the committed state forward. Both laws here are "branch-explicit":
// stress at a trial strain depends only on (committed state, trial strain),
// never on the sub-step path, so the pattern is exact, not an approximation.
//
// LAW 1 — BILINEAR STEEL, PURE KINEMATIC HARDENING ("Steel01"-class, i.e. the
// Menegotto-Pinto ancestor with sharp corners). Rate-independent 1-D
// plasticity with a linear kinematic back-stress (Simo & Hughes,
// "Computational Inelasticity", box 1.4):
//     sigma  = E (eps - eps_p)
//     f      = |sigma - alpha| - fy        (yield function)
//     alpha  = H eps_p,  H = b E / (1 - b) (so the elastoplastic tangent is
//                                           E H / (E + H) = b E exactly)
// Return mapping: sigma_trial = E (eps - eps_p_committed); if f_trial > 0,
// dgamma = f_trial / (E + H), eps_p += dgamma * sign(sigma_trial - alpha).
// Unloading is ELASTIC at slope E from any point; the pure kinematic rule
// translates the yield surface without inflating it, which yields the
// Bauschinger effect: after a plastic excursion the reversed-yield stress is
// lowered by exactly the hardened amount (elastic range always 2 fy wide).
// The stabilized symmetric-cycle loop is a parallelogram of closed-form area
//     W_cycle = 4 fy (1 - b) (eps_a - fy/E)      per unit volume
// (sides of slope E and bE).
//
// LAW 2 — CONCRETE WITH KARSAN-JIRSA UNLOADING/RELOADING. Envelope = the
// existing EC2 parabola-rectangle (EN 1992-1-1 §3.1.7), zero tension, crushed
// to zero beyond eps_cu2 — bit-compatible with `UniaxialLaw::ConcretePR` on
// virgin loading. Unloading from the envelope point (eps_m, sigma_env(eps_m))
// goes LINEARLY to the residual plastic strain eps_p on the strain axis, and
// reloading retraces the same line back to the envelope (the single
// unload/reload line of the classic "Concrete01" simplification of the
// Karsan-Jirsa focal model — no inner-loop stiffness degradation). The
// residual strain is Karsan & Jirsa's regression (I.D. Karsan, J.O. Jirsa,
// "Behavior of concrete under compressive loadings", ASCE J. Struct. Div.
// 95(ST12), 1969), in eps_c2 units with eta = |eps_m| / eps_c2:
//     eps_p / eps_c2 = 0.145 eta^2 + 0.127 eta            (eta <= 2)
//     eps_p / eps_c2 = 0.707 (eta - 2) + 0.834            (eta >  2)
// (the linear extension is C1-continuous at eta = 2 and keeps eps_p < |eps_m|).
// State: eps_min, the most compressive committed strain (<= 0). Tension side
// of the residual: zero stress, zero tangent (open crack), like the envelope.
//
// Units: SI (Pa, strains dimensionless). Sign: tension-positive strains and
// stresses, matching `fiber.h`.
#pragma once
#include "keta/fiber.h"

namespace keta {

// Material memory of ONE fiber. A plain POD: copy = save, assign = restore.
struct HystState {
  double eps_p = 0;    // SteelKinematic: committed plastic strain
  double eps_min = 0;  // ConcreteKJ: most compressive committed strain (<= 0)
  // --- appended (the laws of G10) ---
  double alpha = 0;    // SteelMixed: back stress
  double ebar = 0;     // cumulated (effective) plastic strain
  double r = 0;        // damage threshold reached (0: the initial r0)
  // --- appended: CREEP (a law with creep_a non-empty) ---
  double cr_k[6] = {0, 0, 0, 0, 0, 0};  // Kelvin unit strains (their sum: the creep strain)
  double cr_s = 0;                       // the committed driving stress g(sigma) sigma
  double cr_dt = 0;                      // duration of the step being solved [days] (set by the driver)
  double age = 0;                        // concrete age at the start of that step [days] (set by the driver)
  double sh = 0;                         // committed shrinkage strain (imposed)
  double age0 = 0;                       // age at the activation (the shrinkage counts from it)
  // the shrinkage at the end of the step given by the driver (drying of the
  // section, keta/drying.h) instead of the law's table
  double sh_drv = 0;
  int sh_mode = 0;  // 1: sh_drv
};

struct HystLaw {
  enum class Kind { Memoryless, SteelKinematic, ConcreteKJ, SteelMixed, AsymVonMises, DamageSimoJu, PlasticDamageOliver };
  Kind kind = Kind::Memoryless;
  UniaxialLaw base;  // Memoryless: the wrapped path-independent law
  // SteelKinematic:
  double E = 210e9;
  double fy = 355e6;
  double b = 0.0;  // hardening ratio, 0 <= b < 1 (plastic tangent = b E)
  // ConcreteKJ (compression negative, like ConcretePR):
  double fcd = 20e6;
  double eps_c2 = 0.002, eps_cu2 = 0.0035, n_exp = 2.0;

  // --- appended: mixed hardening, asymmetric von Mises, damage laws ---
  // SteelMixed: plastic modulus H (slope of the yield stress against the
  // plastic strain; elasto-plastic tangent E H / (E + H)), share `iso` of
  // ISOTROPIC hardening (1: isotropic, 0: kinematic = SteelKinematic with
  // b = H / (E + H)): f = |sigma - alpha| - (fy + iso H ebar), alpha_dot =
  // (1 - iso) H eps_p_dot; exact return mapping (linear hardening).
  double H = 0, iso = 1;
  // AsymVonMises — a "modified von Mises", uniaxial: distinct tension
  // and compression yield stresses with their own isotropic hardenings on
  // one effective plastic strain ebar, sigma_t = ft0 + Ht ebar, sigma_c =
  // fc0 + Hc ebar (fc0 > 0, a magnitude); the flow of the pressure-dependent
  // surface, d eps_p = d ebar (sigma + gamma/2) / sqrt(sigma^2 + gamma^2/2),
  // gamma = sigma_c - sigma_t, integrated by
  // backward Euler over each step.
  double ft0 = 0, Ht = 0, fc0 = 0, Hc = 0;
  // DamageSimoJu: elastic-damage, r = max over the history of
  // sqrt(sigma^2 / E) (nominal stress), r >= r0,
  // d = 1 - r0 (1 - A) / r - A exp(B (r0 - r)), sigma = (1 - d) E eps.
  // PlasticDamageOliver: d = 1 - (r0 / r) exp(A (1 - r / r0)), r = beta
  // sqrt(sigma^2 / E), beta = 1 in tension, 1/n in compression; plus
  // isotropic plasticity of the nominal stress |sigma| <= fy + h ebar
  // (h = H): eps = sigma / ((1 - d) E) + eps_p.
  double r0 = 0, A = 0, B = 0, n_ratio = 1;
  // the stress-driven strain of these laws from the committed state (they
  // are formulated in the nominal stress; eps(sigma) is strictly
  // increasing, so `response` inverts it exactly — safeguarded Newton)
  double strain_at(const HystState& st, double sig) const;
  // the state after a step to the stress sig (stress-driven commit)
  HystState update_stress(const HystState& st, double sig) const;
  double damage(const HystState& st) const;  // the committed damage d

  static HystLaw steel_mixed(double E, double fy, double H, double iso);
  static HystLaw asym_von_mises(double E, double ft0, double Ht, double fc0, double Hc);
  static HystLaw damage_simo_ju(double E, double r0, double A, double B);
  static HystLaw plastic_damage_oliver(double E, double fy0, double h, double r0, double A, double n);

  // Trial response: committed state + trial strain -> (stress, tangent).
  // Never mutates; safe to call any number of times per step.
  void response(const HystState& st, double eps, double& sig, double& Et) const;
  // Committed state after accepting the step at strain `eps`.
  HystState update(const HystState& st, double eps) const;

  // --- appended: CREEP of the fiber (non-aging Kelvin chain, the rate form
  // of a creep compliance J(t, t0) - 1/E: fit_kelvin_chain of
  // keta/hereditary.h, or given). Units a_i [1/Pa] (strain per unit
  // stress), tau_i [days], at most 6. The creep strain sum e_i enters the
  // law as an imposed strain: sigma = law(eps - sum e_i). Each unit
  // e_i' = (a_i s(sigma) - e_i) / tau_i is integrated EXACTLY over a step of
  // HystState::cr_dt days for a driving stress s linear in the step (the
  // exponential algorithm: e_i <- e_i E + a_i s0 (1 - E) + a_i (s1 - s0) (1 -
  // (1 - E) / x), E = exp(-x), x = dt / tau_i), the end stress solved
  // implicitly per fiber (monotone scalar equation, Newton to round-off) and
  // its consistent tangent Et / (1 + Et B s'(sigma)). s(sigma) = g sigma with
  // EN 1992-1-1 §3.1.4(4) NONLINEAR creep when creep_fck > 0: g = exp(1.5
  // (|sigma| / fck - 0.45)) for a compression beyond 0.45 fck (1 otherwise).
  // Under a constant stress the creep strain is g(sigma) sigma sum a_i (1 -
  // exp(-t / tau_i)) exactly.
  std::vector<double> creep_a, creep_tau;
  double creep_fck = 0;
  bool creeps() const { return !creep_a.empty() || !shrink.empty(); }
  // --- appended: SHRINKAGE, an imposed strain eps_sh(age) - eps_sh(age0)
  // of the fiber (age in days; age0 its age at its bar's activation, the
  // stage's activation_age; NEGATIVE = shortening): the table (age, eps_sh),
  // piecewise linear, flat beyond its ends (e.g. sampled from
  // ec2_shrinkage, keta/creep.h). Taken at the end of each time step,
  // with the creep of the same step.
  std::vector<std::pair<double, double>> shrink;
  double shrink_at(double age) const;
  void response_mech(const HystState& st, double eps, double& sig, double& Et) const;  // without creep
  HystState update_mech(const HystState& st, double eps) const;
  // Upper bound of the tangent over all states (regularisation reference).
  double sup_tangent() const;

  // Karsan-Jirsa residual plastic strain (returned >= 0, strain units) for a
  // most-compressive strain eps_min <= 0. Exposed for direct evaluation.
  double kj_residual(double eps_min) const;

  static HystLaw wrap(const UniaxialLaw& l);
  static HystLaw steel_kinematic(double E, double fy, double b);
  static HystLaw concrete_kj(double fcd, double eps_c2 = 0.002, double eps_cu2 = 0.0035,
                             double n = 2.0);
};

// Fiber section over history-carrying laws. Same fiber geometry and strain
// plane convention as `FiberSection` (eps = eps0 + kappa_y z - kappa_z y,
// resultants (N, My, Mz)); the caller owns one `HystState` per fiber
// (`virgin_states()` builds the zeroed vector) and passes it to every query,
// so one section object can serve many integration points.
struct HystFiberSection {
  std::vector<Fiber> fibers;
  std::vector<HystLaw> laws;

  std::vector<HystState> virgin_states() const {
    return std::vector<HystState>(fibers.size());
  }
  Vec3 resultants(const Vec3& e, const std::vector<HystState>& st) const;
  Mat3 tangent(const Vec3& e, const std::vector<HystState>& st) const;
  Mat3 sup_stiffness() const;  // sum of sup-tangent fiber contributions (k0)
  // State determination at COMMITTED state `st`: find e with
  // resultants(e, st) = target. Scaled-residual Newton with sup-tangent
  // fallback and backtracking, warm-startable from `e_start` (a cracked/cyclic
  // section solved from the previous strain plane converges where a cold
  // start can stall on a zero-stiffness branch). Unlike
  // `FiberSection::solve_state` it follows the fibres past the strain bounds
  // of their laws: a history may crush or break them. Throws on unreachable
  // targets.
  Vec3 solve_state(const Vec3& target, const std::vector<HystState>& st,
                   const Vec3& e_start = Vec3::Zero(), int max_iter = 100,
                   double tol = 1e-8) const;
  // Roll the committed states forward to the accepted strain plane `e`.
  void commit(const Vec3& e, std::vector<HystState>& st) const;

  static HystFiberSection wrap(const FiberSection& s);  // memoryless clone
  static HystFiberSection rect(double b, double h, int ny, int nz, const HystLaw& law);
  void add_rebar(double y, double z, double A, const HystLaw& law);
};

}  // namespace keta

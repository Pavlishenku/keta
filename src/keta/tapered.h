// Tapered (variable-inertia) bar kernel — ST1 `CARA VAR` — plus the `CARA SPECIAL`
// force-capped element helpers.
//
// ============================ 1. CARA VAR ====================================
//
// ST1 v24 manual, §2.10 "Caractéristiques R.D.M. des barres" (pages 59-68) knows
// three variation laws for a bar's section characteristics:
//
//   * `CARA VAR LIN Y` / `CARA VAR LIN Z` (p. 59 syntax, p. 62-64 rules) — the bar
//     behaves like a RECTANGULAR section one of whose dimensions varies LINEARLY
//     along the bar; the OR (x = 0) and EX (x = L) values of every characteristic
//     are given and each one is interpolated with the power that this dimension
//     carries in the rectangle formulas.  Writing h(x) for the linearly varying
//     dimension and b for the constant one, the manual states verbatim:
//       - `LIN Y`, OPTION PLANE (p. 62): "L'inertie IZ varie en degré 3 en
//         b×h(x)³/12" and, for `LIN Z`: "L'inertie IZ varie linéairement en
//         b(x)×h³/12".
//       - `LIN Y`, OPTION GRILL (p. 63): "L'inertie IY varie linéairement en
//         b(x)×h³/12"; `LIN Z`: "L'inertie IY varie en degré 3 en b×h(x)³/12".
//       - `LIN Y`, OPTION SPATIALE (p. 64): "L'inertie IY varie linéairement …
//         et l'inertie IZ varie en degré 3".
//     i.e. `LIN Y` makes the dimension along local y vary, so IZ = b·h_y³/12 goes
//     as h_y³ while IY = h_y·b³/12 and SX = b·h_y go as h_y; `LIN Z` swaps the two
//     inertias.  keta implements exactly that, as
//         C(t) = [ (1-t)·C_OR^(1/p) + t·C_EX^(1/p) ]^p ,   t = x/L,
//     which is "C is the p-th power of a quantity varying linearly from C_OR^(1/p)
//     to C_EX^(1/p)" — the only reading of "varie en degré 3 de 1.0 à 2.0" that is
//     consistent with the rectangle rule.  Exponents used (see `TaperLaw`):
//         SX, SY, SZ : p = 1  (all area-like, proportional to h)
//         LIN Y      : IY p = 1, IZ p = 3        LIN Z : IY p = 3, IZ p = 1
//         IX         : p = 1  (keta choice — see "Deviations" below)
//     The p = 3 branch is evaluated through `cbrt`, so a p = 1 characteristic is
//     bit-exactly linear and OR == EX reproduces the constant section exactly.
//
//   * `CARA VAR PARA` (p. 59 syntax, p. 61 example) — "définition des poutres dont
//     l'INVERSE des caractéristiques est interpolable par une parabole", from the
//     three control stations OR (x = 0), MI (x = L/2), EX (x = L).  Hence
//         1/C(t) = (1-t)(1-2t)/C_OR + 4t(1-t)/C_MI + t(2t-1)/C_EX ,
//     the quadratic Lagrange basis at t = 0, ½, 1.  A control set whose parabola
//     crosses zero inside the bar is a typed rejection, never a silent infinity.
//
//   IPOL (polar inertia for rotary mass) is NOT part of `SectionProps` and is
//   therefore excluded from the interpolation, as in ST1.  A characteristic left
//   undefined (0 at every control station) stays 0 — for SY/SZ that is ST1's
//   documented way of switching shear deformation off (p. 61, p. 66: "ne pas
//   définir de section réduite de cisaillement revient à ne pas prendre en compte
//   l'effet du tranchant").  A characteristic that is 0 at one station and > 0 at
//   another is a typed rejection.
//
// ---- Exact stiffness by flexibility integration -----------------------------
//
// ST1's doctrine (and keta's) is to use the exact R.D.M. solution rather than an
// FE interpolation.  For a variable section there is no closed form, but the
// FLEXIBILITY of the statically determinate cantilever base system (end i
// clamped, end j free) is an exact quadrature of the internal-force fields:
//
//   F = ∫₀^L Uᵀ(x) D(x) U(x) dx ,   D = diag[1/EA, 1/GA_vy, 1/GA_vz, 1/GJ, 1/EIy, 1/EIz]
//
// where the rows of U(x) are the internal resultants (n, v_y, v_z, t, m_y, m_z) at
// the cut x produced by the six UNIT end actions at j, reduced at (x, 0, 0):
//
//   n = F_x ,  v_y = F_y ,  v_z = F_z ,  t = M_x ,
//   m_y = M_y − (L−x)·F_z ,  m_z = M_z + (L−x)·F_y .
//
// F is symmetric positive definite; K_jj = F⁻¹ is the 6×6 end stiffness of the
// clamped-free member, and the full 12×12 follows by RIGID EQUILIBRIUM TRANSFER
//   K = Bᵀ K_jj B ,  B = [ −Gᵀ  I ] ,  G = [[I, 0], [S(r), I]] ,  r = (L, 0, 0)
// (p_i = −G p_j is force + moment equilibrium of the member; d_j^rigid = Gᵀ d_i is
// the conjugate kinematics).  For a constant section this reproduces the prismatic
// Timoshenko 12×12 of `element.cpp` bit-for-bit (tested to 1e-10), shear included:
// the shear terms enter as ∫ v_i v_j /(G A_v) dx exactly like φ = 12EI/(GA_vL²) does.
//
// ---- Equivalent (fixed-end) load vectors ------------------------------------
//
// `tapered_f0` returns f0 with keta's element convention: f0 is the WORK-EQUIVALENT
// NODAL LOAD vector (it is added to the global load vector), and the member end
// forces are `k·u − f0`, so the clamped-clamped fixed-end forces are −f0.  It is
// built by the unit-load theorem on the same base system:
//   d⁰ = ∫ Uᵀ D s⁰ dx   (s⁰ = internal resultants of the SPAN LOAD alone on [x, L]),
//   p_j = −K_jj d⁰ ,  p_i = −G p_j − q_res  (q_res = the span load reduced at node i),
//   f0 = −[p_i ; p_j].
// The end-i part is obtained from equilibrium, so the total force and the total
// moment carried by f0 match the applied load to machine precision by construction.
// The quadrature is SPLIT at every load breakpoint (segment ends, point-load
// abscissae) so the piecewise-smooth s⁰ fields are integrated accurately.
//
// ---- Station recovery -------------------------------------------------------
//
// `tapered_stations` mirrors `statics.cpp::recover_stations`: internal forces are
// EXACT (equilibrium reduction of the end-i member forces plus the span loads on
// [0, x], sign convention identical to `keta::Station`), and deflections come from
// a double integration of the curvature M(x)/(E I(x)) (+ the shear angle
// V(x)/(G A_v(x)) and, axially, N(x)/(E A(x))) anchored on the end displacements.
//
// ---- Quadrature -------------------------------------------------------------
//
// Gauss-Legendre, 16 points by default, nodes/weights generated by Newton
// iteration on the Legendre polynomials (deterministic, no tables).  Rationale:
//   * `VAR PARA` — 1/I(x) IS the parabola, so the integrands ξ^k/(E I(x)) are
//     polynomials of degree ≤ 4: 16 points are exact (as are 3).
//   * `VAR LIN` — the integrands are ξ^k/(a + bξ)^p, analytic on [0, L] with the
//     nearest pole outside it.  Gauss-Legendre then converges geometrically as
//     ρ^(−2n) with ρ the Bernstein-ellipse parameter of that pole; measured on the
//     12x12 stiffness, 16 points reach 6e-16 relative for a 2:1 taper and 4e-7 for
//     an extreme 10:1 taper — hence 16 as the default and an explicit `gauss`
//     argument for the pathological cases.
//   * a constant section is degree ≤ 2 and therefore exact.
//
// ============================ 2. CARA SPECIAL ================================
//
// ST1 manual p. 210, §10.13 "Éléments spéciaux": "La formulation est celle d'une
// barre classique RDM, à la différence que l'on peut définir un seuil limite en
// effort, pour les différentes composantes d'une barre (effort normal,
// cisaillement, flexion et torsion)" and "Si l'on définit une valeur nulle selon
// une composante d'effort […] le comportement de la composante reste élastique,
// sinon la composante est limitée par l'effort non nul défini".
// This module provides the PURE per-element pieces (demand, ratio, clamp, secant
// scaling); the iteration itself belongs to the runtime driver.
//
// ============================ Units & limitations ============================
//
// SI everywhere: m, N, Pa, kg, rad.  All matrices/vectors are in the ELEMENT LOCAL
// frame with keta's DOF order [ux uy uz rx ry rz]_i then _j; N > 0 in tension.
//
// Deviations from / gaps in the manual, all documented rather than silent:
//   * p. 64 prints, for `CARA VAR LIN Z` in OPTION SPATIALE, "l'inertie IZ varie en
//     degré 3" — inconsistent with p. 62 (OPTION PLANE, `LIN Z`: "IZ varie
//     linéairement") and with the rectangle rule itself.  keta follows the
//     rectangle rule and p. 62: for `LIN Z`, IZ is linear.  Reading p. 64 literally
//     would make both inertias cubic, which no rectangle can do.
//   * The manual never states how IX (torsion) varies.  A rectangle's St-Venant
//     constant is not a power of h (β(b/h)·b³h), so no exponent is exact; keta
//     interpolates IX with p = 1, the thin-rectangle limit, and documents it here.
//     IX must be > 0 along the bar: the flexibility method inverts ∫dx/(GJ), so a
//     zero torsion constant (OPTION PLANE, where IX is simply not given) is a typed
//     rejection — the caller supplies a value.  SX, IY, IZ must likewise be > 0.
//   * Wel/Wpl are not ST1 characteristics (ST1 carries VY/WY/VZ/WZ instead).  They
//     are interpolated with the consistent rectangle exponents (Wel = I/v, v ∝ h in
//     the varying direction) and take part in NO stiffness computation.
//   * Span loads: UDL, LIN (linearly varying) and CON (concentrated force/moment)
//     only.  SELF-WEIGHT of a `VAR LIN` bar lowers EXACTLY to a `Linear` span load
//     (SX is linear under both LIN laws); under `VAR PARA` it does not, and the
//     caller must subdivide.  Thermal / imposed-strain loads (uniform or linear
//     along the bar) lower exactly through `tapered_f0_strain`; a concentrated
//     dislocation (DEFOR CON) is not lowered.
//   * Releases, rigid-arm eccentricity and the Winkler bed are NOT folded in: this
//     module returns the RAW member matrices and the caller (the ST1 layer)
//     condenses / transforms them, exactly as `BarElement` does for prismatic bars.
//   * Deflection recovery is a Simpson double integration over the stations, so it
//     is O(h⁴) in the station spacing rather than exact (with a variable EI there
//     is no polynomial closed form); the internal FORCES are exact at every station.
#pragma once
#include "keta/model.h"

namespace keta {

using Mat6 = Eigen::Matrix<double, 6, 6>;

// Default Gauss-Legendre order for the flexibility / load integrals (see header doc).
inline constexpr int kTaperedGauss = 16;

// ---------------------------------------------------------------------------
// Variation law
// ---------------------------------------------------------------------------

enum class TaperLaw {
  LinY,  // CARA VAR LIN Y — rectangle whose dimension along local y varies linearly
  LinZ,  // CARA VAR LIN Z — rectangle whose dimension along local z varies linearly
  Para,  // CARA VAR PARA  — the INVERSE of each characteristic is a parabola OR/MI/EX
};

// One bar's variation descriptor. `mid` is read only by TaperLaw::Para.
// The `SectionProps` fields map onto the ST1 keywords as
//   A <- SX, Avy <- SY, Avz <- SZ, J <- IX, Iy <- IY, Iz <- IZ.
struct TaperedVar {
  TaperLaw law = TaperLaw::LinY;
  SectionProps orig;  // OR — x = 0
  SectionProps mid;   // MI — x = L/2 (Para only)
  SectionProps ext;   // EX — x = L
};

// Interpolated characteristics at the NORMALIZED abscissa t = x/L in [0, 1].
// Typed rejection on a mixed 0 / non-0 control set or a sign-changing PARA parabola.
SectionProps tapered_props_at(const TaperedVar& var, double t);

// Section compliances at t: [1/EA, 1/GA_vy, 1/GA_vz, 1/GJ, 1/EIy, 1/EIz].
// A zero shear area yields 0 (ST1: shear deformation not taken into account).
Vec6 tapered_compliance(const TaperedVar& var, const Material& mat, double t);

// ---------------------------------------------------------------------------
// Stiffness
// ---------------------------------------------------------------------------

// 6x6 end flexibility of the cantilever base system (end i clamped, end j free),
// in the local frame, DOF order [ux uy uz rx ry rz] at end j.
Mat6 tapered_flexibility(const TaperedVar& var, const Material& mat, double L,
                         int gauss = kTaperedGauss);

// K_jj = flexibility^-1 (6x6 stiffness of the clamped-free member at end j).
Mat6 tapered_end_stiffness(const TaperedVar& var, const Material& mat, double L,
                           int gauss = kTaperedGauss);

// EXACT 12x12 local stiffness of the tapered member (no releases, no eccentricity,
// no soil — the caller folds those in).
Mat12 tapered_ke(const TaperedVar& var, const Material& mat, double L,
                 int gauss = kTaperedGauss);

// ---------------------------------------------------------------------------
// Span loads
// ---------------------------------------------------------------------------

// A span load on the tapered bar, LOCAL frame, ABSOLUTE abscissae [m]
// (the ST1 layer resolves REL/ABS/XL and rotates GLO/PRO intensities to local).
struct TaperedSpanLoad {
  enum class Kind {
    Udl,     // ST1 CHARG BARRE UNI  — constant q over [a, b]
    Linear,  // ST1 CHARG BARRE LIN  — q varies linearly qa -> qb over [a, b]
    Point,   // ST1 CHARG BARRE CON  — concentrated force `qa` + moment `ma` at a
  };
  Kind kind = Kind::Udl;
  Vec3 qa = Vec3::Zero();  // Udl/Linear: intensity at a [N/m]; Point: force at a [N]
  Vec3 qb = Vec3::Zero();  // Linear: intensity at b [N/m]
  Vec3 ma = Vec3::Zero();  // Point: concentrated moment at a [N·m]
  double a = 0.0;          // start abscissa [m]
  double b = -1.0;         // end abscissa [m]; < 0 = full length (Udl/Linear only)

  static TaperedSpanLoad udl(const Vec3& q, double a = 0, double b = -1);
  static TaperedSpanLoad linear(const Vec3& qa, const Vec3& qb, double a = 0, double b = -1);
  static TaperedSpanLoad point(double a, const Vec3& F, const Vec3& M = Vec3::Zero());
};

// Reduction at the point (x, 0, 0) of the span loads acting on [s0, s1]:
// F = the resultant force, M = the resultant moment about that point.
// `open_at_s1` excludes a concentrated load sitting exactly at s1 (the value just
// before the jump) — used by station recovery, as `statics.cpp` does.
struct SpanResultant {
  Vec3 F = Vec3::Zero();
  Vec3 M = Vec3::Zero();
};
SpanResultant tapered_reduce_loads(const std::vector<TaperedSpanLoad>& loads, double L, double s0,
                                   double s1, double x, bool open_at_s1 = false);

// Work-equivalent NODAL load vector f0 (local frame, 12): the member end forces
// are k·u − f0, so the clamped-clamped fixed-end forces are −f0.
Vec12 tapered_f0(const TaperedVar& var, const Material& mat, double L,
                 const std::vector<TaperedSpanLoad>& loads, int gauss = kTaperedGauss);

// Convenience wrappers for the three ST1 span-load families.
Vec12 tapered_f0_udl(const TaperedVar& var, const Material& mat, double L, const Vec3& q,
                     double a = 0, double b = -1, int gauss = kTaperedGauss);
Vec12 tapered_f0_linear(const TaperedVar& var, const Material& mat, double L, const Vec3& qa,
                        const Vec3& qb, double a = 0, double b = -1, int gauss = kTaperedGauss);
// Imposed generalized strains (thermal, ST1 DEFOR) varying LINEARLY from `ea`
// at x = 0 to `eb` at x = L, each (eps0 axial, kappa_y, kappa_z) with keta's
// Prestrain convention: the exact work-equivalent vector of the tapered bar
// (unit-load theorem on the base system, no applied force, so f0 is
// self-equilibrated). A constant section reproduces BarElement::f0_local.
Vec12 tapered_f0_strain(const TaperedVar& var, const Material& mat, double L, const Vec3& ea,
                        const Vec3& eb, int gauss = kTaperedGauss);
Vec12 tapered_f0_point(const TaperedVar& var, const Material& mat, double L, double a,
                       const Vec3& F, const Vec3& M = Vec3::Zero(), int gauss = kTaperedGauss);

// ---------------------------------------------------------------------------
// Station recovery
// ---------------------------------------------------------------------------

// Same fields and same sign convention as `keta::Station` (statics.h).
struct TaperedStation {
  double x = 0;  // abscissa along the bar [m]
  double N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0;
  Vec3 defl = Vec3::Zero();  // local deflections (u, v, w) [m]
};

// `f_end_local` = member end forces (k·u − f0, local frame, 12);
// `u_local` = local end displacements (12). Internal forces are exact; deflections
// are a Simpson double integration anchored on the end displacements.
std::vector<TaperedStation> tapered_stations(const TaperedVar& var, const Material& mat, double L,
                                             const std::vector<TaperedSpanLoad>& loads,
                                             const Vec12& f_end_local, const Vec12& u_local,
                                             int n_stations);

// ---------------------------------------------------------------------------
// CARA SPECIAL — force-capped elastic element (manual p. 210)
// ---------------------------------------------------------------------------

struct SpecialLimits {
  // Threshold per component [N] for fx/fy/fz, [N·m] for mx/my/mz.
  // 0 = UNLIMITED: "si l'on définit une valeur nulle selon une composante d'effort
  // […] le comportement de la composante reste élastique" (p. 210).
  std::array<double, 6> cap{};
  bool any() const;
  // Finiteness / non-negativity sweep; typed rejection otherwise.
  void validate() const;
};

// Demand per component = max |end force| over the two ends of the member.
std::array<double, 6> special_demands(const Vec12& f_local);

// Secant ratio r_c = min(1, cap_c / demand_c) in (0, 1]; 1 when the component is
// elastic (cap 0, or demand below the cap, or zero demand).
std::array<double, 6> special_ratios(const SpecialLimits& lim, const Vec12& f_local);

// Member end forces clamped onto the caps: both end entries of a component are
// scaled by r_c.  For N and T (equal and opposite at the two ends) this preserves
// member equilibrium exactly; for the shear/bending components it is the SECANT
// reduction of that component only — the global equilibrium is restored by the
// runtime's iteration on `special_secant_stiffness`, never by an unbalanced return.
Vec12 special_clamp(const SpecialLimits& lim, const Vec12& f_local);

// Per-DOF secant scale s_d = sqrt(r_c(d)), so that k_sec = S k S (S = diag(s)) is a
// CONGRUENCE of the elastic stiffness: symmetric, positive semi-definite, and it
// scales the diagonal term of a capped component exactly by r_c (the off-diagonal
// coupling between a capped shear and a capped moment takes the geometric mean —
// the standard diagonal-secant compromise, documented, never silent).
std::array<double, 12> special_dof_scales(const SpecialLimits& lim, const Vec12& f_local);
Mat12 special_secant_stiffness(const Mat12& k_local, const std::array<double, 12>& scales);

// Convergence measure for the runtime iteration: max over capped components of
// demand/cap − 1, floored at 0 (0 = every component within its threshold).
double special_overshoot(const SpecialLimits& lim, const Vec12& f_local);

// Installs the exact tapered kernel of bar `bar` on `model.tapered`: the
// 12x12 member stiffness, the section law along the bar for the station
// recovery, and the equivalent nodal loads of the span loads (uniform, linear,
// point, thermal, imposed strains, self-weight; a concentrated dislocation in
// local axes). Call it once the bar's nodes, eccentricities, roll and material
// are final. Typed errors: unknown bar, unsupported load kind at lowering.
void install_tapered(Model& model, Id bar, const TaperedVar& var);

}  // namespace keta

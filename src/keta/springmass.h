// Full elastic supports (ST1 `APPUI EL`), all-or-nothing uplift (`DECOL`) and
// modal masses (`MASSE PROPRE`) — builders + assembly helpers.
//
// This module owns three ST1 features that the rest of the kernel only *hooks*:
//
//  1. FULL ELASTIC SUPPORTS — ST1 manual v24 p. 47-51 (`APPUI ... EL`).
//     ST1 gives a support stiffness as the LOWER TRIANGLE of a symmetric
//     matrix, read ROW BY ROW (p. 49 for the syntax, p. 50 for the worked
//     example, p. 51 for the effort = matrix x displacement statement):
//
//         SPATIALE, 21 terms, DOF order (dx dy dz rx ry rz) — keta's own order:
//             m1
//             m2  m3
//             m4  m5  m6
//             m7  m8  m9  m10
//             m11 m12 m13 m14 m15
//             m16 m17 m18 m19 m20 m21
//
//         PLANE, 6 terms, DOF order (dx dy rz) -> keta DOFs (0, 1, 5):
//             m1
//             m2  m3
//             m4  m5  m6
//
//     `EL DI` gives only the diagonal (3 terms PLANE, 6 terms SPATIALE).
//     `ROT alpha` (p. 49, radians, about global Z, trigonometric sense) says
//     the matrix is expressed in the frame rotated by alpha; see
//     `rotate_support_stiffness` for the derivation of K_global.
//
//     The matrix is stored VERBATIM (support frame) in `Support::spring6` with
//     the angle in `Support::rot`; assembly asks `support_stiffness_global()`
//     for the rotated global-frame matrix. The per-DOF `Support::spring`
//     diagonal stays an independent, additive path (both are summed into K).
//
//     Reactions follow the existing diagonal-spring convention exactly: an
//     elastic support DOF stays FREE, so the reported reaction R = K_elem u - F
//     equals -K_spring u by nodal equilibrium — nothing is added at the
//     reaction site (see `recover_case`). A DOF that is BOTH rigidly blocked
//     and elastically supported is eliminated, so its spring force is not part
//     of the reaction (pre-existing convention of the diagonal path).
//
//  2. DECOL, ALL-OR-NOTHING (p. 47 and p. 49). The uplift criterion bears on
//     ONE global DOF (X, Y or Z); when it is met, *every* DOF of that support
//     is released. `POS` (default) = the support stays glued while its reaction
//     on that axis is >= 0. Statics only (ST1 forbids it with moving loads);
//     the active-set loop lives in `statics.cpp`.
//
//     This is a support-level switch, orthogonal to the pre-existing per-DOF
//     `Support::oneside` (which releases only the offending DOF). Combining
//     both on the same support is a typed error.
//
//  3. MASSE PROPRE (p. 69-71). Nodal mass matrices (full lower triangle or
//     diagonal, same ordering as the support stiffness above, `SUP` = add to
//     what is already there) and per-bar line masses `UNI m1` / `LIN m1 m2` /
//     `PARABOL m1 m2 m3` in kg/m, optionally attached with an eccentricity
//     (`EXC OR/EX`, global frame), plus `MASSE NUL` which cancels a bar's
//     automatic rho*A mass. ST1's automatic bar mass is RO/G_DYN * A; keta's
//     `Material::rho` is already a mass density, so it is simply rho*A.
//
//     Line masses are TRANSLATIONAL ONLY (ST1 p. 69: the added line masses are
//     masses, not inertias — no rotary inertia term, no torsional term). They
//     are integrated with the SAME shape functions as the element consistent
//     mass (linear axial, cubic Hermite transverse with the theta_y = -dw/dx
//     flip), so a constant mu reproduces `BarElement::me()`'s translational
//     blocks bit-for-bit; the varying-mu integral int mu(x) N^t N dx is taken
//     by 5-point Gauss-Legendre (exact through degree 9, the integrand peaks at
//     degree 2+3+3 = 8 for the parabolic law). Line masses stay CONSISTENT even
//     when the caller asks `assemble_mass` for lumped masses — ST1's
//     `MASSE` / `MASSE DI` switch bears on the automatic bar masses only.
//
//     Known limitations, all documented rather than silent:
//       * `DECOL` is honoured by the linear-statics driver only (as in ST1);
//         every other driver sees a bilateral support, exactly like the
//         pre-existing per-DOF `Support::oneside`.
//       * A full 6x6 support is a typed reject in the corotational driver
//         (its residual carries springs as a diagonal vector).
//       * The new `Support` / `NodalMass` / `BarMass` data is not part of the
//         JSON round trip yet (io.cpp).
//
// SI units everywhere: N/m, N·m/rad, N/rad for the stiffness terms, kg,
// kg·m^2 for the masses, kg/m for the line masses, m for the offsets, rad for
// the angles.
#pragma once
#include "keta/assembly.h"

namespace keta {

// ---------------- APPUI EL: builders ----------------

// Symmetric 6x6 from ST1's row-by-row lower triangle (manual p. 49-51).
// `plane` = true: 6 terms on the PLANE DOFs (dx, dy, rz) = keta DOFs 0, 1, 5.
// `plane` = false: 21 terms on the full SPATIALE order (dx dy dz rx ry rz).
// Throws on a wrong term count or a non-finite term.
Mat6 sym_from_lower_triangle(const std::vector<double>& values, bool plane);

// Symmetric 6x6 from ST1's diagonal form `... DI d1 ... dj` (p. 49, p. 71).
// `plane` = true: 3 terms on (dx, dy, rz); false: 6 terms on the full order.
Mat6 sym_from_diagonal(const std::vector<double>& values, bool plane);

// ST1 `APPUI ... EL m1 ... m21` (or 6 terms in PLANE). Alias of
// `sym_from_lower_triangle` that additionally rejects a negative diagonal
// term — a support stiffness is positive semi-definite by construction.
Mat6 spring_from_lower_triangle(const std::vector<double>& values, bool plane);
// ST1 `APPUI ... EL DI d1 ... dj`.
Mat6 spring_from_diagonal(const std::vector<double>& values, bool plane);

// ST1 `MASSE PROPRE NOEUD ... MASSE m1 ...` / `... MASSE DI d1 ...` (p. 71).
Mat6 mass_from_lower_triangle(const std::vector<double>& values, bool plane);
Mat6 mass_from_diagonal(const std::vector<double>& values, bool plane);

// Frame map for ST1 `ROT alpha` (p. 49): rows are the support-frame axes
// expressed in the global frame, i.e. v_support = Rz(alpha) * v_global, with
// the support frame obtained by rotating the global frame by +alpha about Z in
// the trigonometric sense.
Mat3 rot_z_frame(double alpha);

// K_global for a support stiffness expressed in the ROT-alpha frame.
// Derivation: with R = blkdiag(Rz, Rz) mapping global -> support frame,
//   u_s = R u_g  and  f_s = K_s u_s.
// R is orthogonal, so forces map the same way (f_s = R f_g) and
//   f_g = R^t f_s = R^t K_s R u_g   =>   K_global = R^t K_s R.
// Sanity check, an X-only spring at alpha = pi/2 becomes a Y-only spring:
// R^t e_x = e_y, so K_g = k (R^t e_x)(R^t e_x)^t = k e_y e_y^t.
Mat6 rotate_support_stiffness(const Mat6& k_support, double alpha);

// The support's full 6x6 in the GLOBAL frame (rotated + symmetrized).
// Zero when the support carries no `spring6`.
Mat6 support_stiffness_global(const Support& s);

// Any support carrying a full 6x6? (drivers that cannot honour it reject.)
bool has_full_spring(const Model& model);
// Any support carrying a DECOL switch? (statics takes the active-set path.)
bool has_decol(const Model& model);

// Scatter every full 6x6 support stiffness into a triplet list, skipping
// blocked / auto-fixed DOFs. Mirrors the per-DOF diagonal spring path.
void scatter_full_springs(const Assembled& asm_, std::vector<Triplet>& trips);

// ---------------- MASSE PROPRE ----------------

// Line mass law mu(xi), xi = x / L in [0, 1] (manual p. 69).
double line_mass_at(const BarMass& bm, double xi);
// Exact analytic total added mass int_0^L mu(x) dx [kg].
double line_mass_total(const BarMass& bm, double L);

// Consistent 12x12 line-mass matrix int mu(x) N(x)^t N(x) dx in the ELEMENT
// LOCAL frame, translational DOFs only (no rotary/torsional inertia).
// Shape functions: linear axial, cubic Hermite transverse — the exact ones
// behind `BarElement::me()`, including the theta_y = -dw/dx sign flip.
// A constant mu equals me()'s translational blocks to round-off.
Mat12 line_mass_local(const BarMass& bm, double L);

// The line mass lumped onto the two bar ends as translational point masses,
// int mu (1 - xi) dx at i and int mu xi dx at j, in all three directions —
// ST1's form on a MASSE NUL bar (a SUP UNI mass on such a portal beam
// feeds the sway mode but adds no bending mode of the beam).
Mat12 line_mass_lumped_local(const BarMass& bm, double L);

// Same, folded to the model NODAL frame: T^t M T with
// T = rotation * rigid-arm(bar.ecc + mass EXC offsets) — keta's own
// eccentricity transform, so the mass line can hang off the bar axis.
// `lumped` picks line_mass_lumped_local over the consistent matrix.
Mat12 line_mass_nodal(const BarElement& el, const BarMass& bm, bool lumped = false);

// Does this bar carry a `MASSE NUL` entry (automatic rho*A mass cancelled)?
bool bar_mass_is_null(const Model& model, Id bar);
// True when `bar` carries an automatic-mass entry (BarMassKind::Auto).
bool bar_mass_is_auto(const Model& model, Id bar);

// Lumped (diagonal) or consistent mass for this bar: its Auto entry's own form
// when it names one, the model-wide `lumped` otherwise.
bool bar_mass_lumped(const Model& model, Id bar, bool lumped);

// `MASSE DI` (manual p.69): the COMPLETE diagonal form, translations AND
// rotations, against the default diagonal form in translations only.
bool bar_mass_full_diagonal(const Model& model, Id bar);

// Effective nodal mass per node index, folding `NodalMass::add` in declaration
// order (add = false replaces the running matrix, add = true accumulates —
// ST1's `SUP`, p. 69). Only nodes carrying at least one entry are returned.
std::vector<std::pair<int, Mat6>> effective_nodal_masses(const Model& model);

// Scatter helpers used by `assemble_mass`.
void scatter_bar_line_masses(const Model& model, const Assembled& asm_,
                             std::vector<Triplet>& trips);
void scatter_nodal_masses(const Model& model, const Assembled& asm_, std::vector<Triplet>& trips);

}  // namespace keta

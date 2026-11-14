// Material-nonlinear facet shells: layered von Mises (J2) plane-stress
// plasticity through the thickness, and the load-controlled Newton driver of
// a shell model (geometrically linear, material-nonlinear analysis with
// plastic steel facets).
//
// MATERIAL. Each thickness point carries the plane-stress state (sx, sy,
// txy), the plastic strain and the equivalent plastic strain alpha. Yield
// f = 1/2 s' P s - 1/3 R(alpha)^2, R = fy + H alpha (H = 0: perfect
// plasticity), P the plane-stress deviatoric projector. Return mapping:
// the closed-form-in-the-eigenbasis algorithm of Simo & Taylor (1986; Simo
// & Hughes, "Computational Inelasticity", §3.4): C and P share the
// eigenvectors (sx + sy, sx - sy, txy), so the backward-Euler stress is
//   s(dg) = (C^-1 + dg P)^-1 C^-1 s_trial
// in closed form and the consistency condition is ONE scalar equation in dg,
// solved by Newton to round-off. The algorithmic (consistent) tangent
//   C_ep = Xi - (Xi P s)(Xi P s)' / (s' P Xi P s + beta)
// makes the global Newton quadratic.
//
// SECTION. Kirchhoff-Love layering on the facet kinematics: eps(z) = eps0 +
// z kappa with the Allman membrane strain at each of the three mid-edge
// points of the membrane rule and the constant MITC3 bending curvature;
// TRANSVERSE SHEAR in the criterion too (J2 on sx, sy, txy, txz, tyz with
// sz = 0): the MITC3 shear strain gbar at the same points, distributed
// through the thickness as f(z) gbar, f = 5/4 (1 - 4 z^2 / t^2) — the
// profile whose energy gives kappa_s = 5/6, so the elastic limit is the
// MITC3 element exactly. Thickness integrated by a composite 4-point
// Gauss-Lobatto rule (exact to degree 5 per layer: z^2 bending and the
// quartic shear energy exact) whose end points are the faces (surface
// stresses are read, not extrapolated). The drilling stabilisation stays
// elastic (a numerical device, never a stress). With an elastic material the
// tangent is ke() exactly.
//
// DRIVER. Equal load increments lambda = k / steps of the case's loads; per
// increment a full Newton with the consistent tangent, converged on the
// residual (relative to the external load) and the increment; the history is
// committed per accepted increment. A load beyond the limit load raises a
// typed NonConvergence (no equilibrium: the section's moment cannot reach
// it). BARS, springs and full 6x6 elastic supports in the same model are the
// linear part of the system: each bar with its exact linear kernel (releases,
// rigid arms, Winkler bed, tapered members — bar_ke_assembled, as linear
// statics), span loads through the case's lowering, end forces k u - lam f0
// recovered per increment. Cables are a typed rejection.
//
// LARGE DISPLACEMENTS (opt.large_displacements): the facet's deformational
// vector d(q) of the EICR kinematics (corot_kin.h: rigid frame of the current
// triangle, centroid-relative translations, log-map rotations) replaces
// T u; the layered response f_l(d), K_l(d) is pulled back exactly,
//   f = J^T f_l,   K = J^T K_l J + d2(f_l . d)/dq2   (J = dd/dq),
// J by first-order and the geometric term by second-order forward AD. Bars
// are corotational elastic members (energy formulation, exact AD tangent).
// Loads: nodal forces and surface pressures as DEAD loads (fixed direction,
// per undeformed area), nodal moments as spatially fixed moments (J(theta)^T
// M, its load stiffness in the tangent). Elastic bars keep rigid arms (exact
// finite-rotation arm kinematics), releases, taper and Winkler beds (ground
// springs on the current axis, linearising to the Hermite bed), and their
// span loads are re-lowered on the current corotated frame (equivalent loads
// of solve_corotational, load stiffness by AD). Rotation vectors are re-canonicalised past 0.9 pi. Newton stops
// at the tolerance (relative to the load, span loads included) or within
// the elements' round-off floor once the last correction is negligible.
#pragma once
#include <map>
#include <string>
#include <vector>

#include "keta/fiber_hyst.h"
#include "keta/plastic_hinge.h"
#include "keta/statics.h"

namespace keta {

// HILL'S orthotropic yield criterion (appended): yield stresses relative to
// the reference fy (R = fy + H alpha scales them all) in the material axes
// 1, 2 (in the facet plane) and 3 (the normal):
//   F (s22 - s33)^2 + G (s33 - s11)^2 + H (s11 - s22)^2 + 2 L s23^2
//     + 2 M s13^2 + 2 N s12^2 = R^2,
//   F = (1/r22^2 + 1/r33^2 - 1/r11^2) / 2 (cyclic), N = 1 / (2 r12^2) ...
// with s33 = 0 (plane stress — r33 still enters through F, G, H). All ratios
// 1 and the shears 1/sqrt(3): von Mises. The plane-stress form may be
// degenerate (semi-definite: an open yield cylinder), never indefinite.
// Material axis 1 = `axis1` (a
// global direction) projected on the facet plane; zero: the facet local x.
// The return is the general backward Euler of an associative quadratic
// criterion, s(dg) = (C^-1 + dg P)^-1 C^-1 s_trial, the consistency a scalar
// equation in dg solved by safeguarded Newton to round-off, and the
// consistent tangent of the J2 algorithm (it holds for any constant P).
struct HillYield {
  double r11 = 1, r22 = 1, r33 = 1;
  double r12 = 0.5773502691896258, r13 = 0.5773502691896258, r23 = 0.5773502691896258;
  Vec3 axis1 = Vec3::Zero();
};

struct ShellPlasticOptions {
  int steps = 20;
  int max_iterations = 50;
  double tol = 1e-10;
  int layers = 40;        // composite 4-point Gauss-Lobatto: 3 layers + 1 thickness points
  double hardening = 0;   // isotropic linear hardening modulus H [Pa]
  // transverse shear in the yield criterion (J2 on (sx, sy, txy, txz, tyz),
  // sz = 0) with the parabolic shear-strain profile through the thickness;
  // false keeps the transverse shear elastic (MITC3, kappa = 5/6)
  bool plastic_transverse_shear = true;
  // LARGE DISPLACEMENTS: facets as element-independent corotational
  // triangles (the EICR kinematics of solve_corotational) carrying the
  // layered plasticity in their corotated frame; bars as corotational
  // elastic members. Small strains, arbitrary rotations.
  bool large_displacements = false;
  // NONLINEAR BARS, in the bar basic system v = (elongation, ty_i, ty_j,
  // tz_i, tz_j, twist) — B_lin u in linear kinematics, the corotational
  // chord-frame deformations (exact AD Jacobian and geometric Hessian) in
  // large displacements:
  //   * force-based FIBER bars (bar id -> hysteretic section, the element of
  //     solve_fiber_frame3d, fiber_points Gauss-Lobatto points). In LINEAR
  //     kinematics: span loads by the exact particular solution of the basic
  //     system (section forces b(x) q + s_p(x); its end reactions join the
  //     nodal forces) on a composite Lobatto rule cut at the load
  //     discontinuities (a point-load abscissa carries two points, the left
  //     and the right limit); RELEASES through a constraint basis q = Cq Q
  //     (axial / moment / torsion: q_c = 0; SHEAR: constant moment in that
  //     plane), the particular field moved off the released actions
  //     (release_statics); span TORQUES (their particular torque's elastic
  //     twist int T_p / GJ) and IMPOSED STRAINS (thermal, prestrain, linear,
  //     dislocations: the total strain in the compatibility, the mechanical
  //     one in the section); rigid arms (T); tapered members (a section per point,
  //     fiber_sections_along); a WINKLER bed (Bar::soil kfx, kfy, kfz) coupled
  //     in the element: the deflection integrated exactly from the section
  //     strains (piecewise-polynomial through the points), its pressure a
  //     span load of the same basic-system statics, compatibility and bed
  //     law solved together by Newton with the exact (non-symmetric) tangent
  //     — a collocation converging spectrally with fiber_points to the beam
  //     on elastic foundation (bef.h), not its closed form; a TORSIONAL bed
  //     (kmx) through the exact torsion block of the beam on foundation
  //     (bef.h: torsion is elastic and uncoupled from the fibers). In LARGE
  //     displacements the same, with the
  //     exact arm kinematics and the span loads on the current frame: local
  //     intensities w(u) (dead: E^T g; follower: constant), particular field
  //     sum S w, basic-system reactions pushed back like equivalent loads;
  //     the tangent carries dq/dw through the converged compliance to the
  //     particular field; a Winkler bed coupled in the corotated frame
  //     (ground springs on the absolute axis displacement, undeformed local
  //     axes), the twist bed the exact torsion block minus its GJ part on
  //     phi = e1_0 . theta;
  //   * concentrated HINGES (HingeSpec: bilinear moment or N-M polygon with
  //     kinematic hardening), the hinge forces from a backward-Euler
  //     multi-surface return (closest point on the polygon, exact for linear
  //     facets) and the algorithmic tangent. In LINEAR kinematics a hinged
  //     bar is the 12-DOF condensation of hinge_pushover on the bar's FULL
  //     local kernel — releases, rigid arms (T), Winkler bed (its global
  //     part), tapered member (exact kernel), span loads (f0 split between
  //     the nodes and the internal hinge DOFs) — the internal DOFs solved per
  //     element to round-off; in large displacements the hinges are series
  //     springs on the basic components (rotation at the hinged end; the
  //     elongation for the axial part of an N-M hinge) of the bar's elastic
  //     basic system (release-condensed or tapered block, rigid arms in the
  //     kinematics), span loads as the fixed-end deformation shift
  //     q = kb (v - v0(w) - P e) + qn(w) plus the pushed basic-system
  //     reactions — qn, outside the range of kb, the part of the fixed-end
  //     forces a released member carries by equilibrium alone (it passes
  //     through the hinges). Shear releases accepted (kb singular along the
  //     sliding mode).
  // A bar carries either.
  std::map<Id, HystFiberSection> fiber_sections;
  // a section per Gauss-Lobatto point (fiber_points of them, ascending x) —
  // required for a tapered fiber bar, allowed for any
  std::map<Id, std::vector<HystFiberSection>> fiber_sections_along;
  int fiber_points = 5;
  std::vector<HingeSpec> hinges;
  // Appended: Hill's criterion for the shells of these materials (material
  // id -> yield ratios and axes); the others stay von Mises
  std::map<Id, HillYield> hill;
  // per shell id: its own material axis 1 (a global direction projected on
  // the facet), overriding HillYield::axis1 — meridians of a dome, ...
  std::map<Id, Vec3> hill_axis1;
  // Appended: a load increment without equilibrium (beyond the limit load)
  // is halved from the last converged state up to limit_cuts times, then the
  // analysis stops there and returns the converged steps (limit_reached)
  // instead of raising NonConvergence
  bool stop_at_limit = false;
  int limit_cuts = 6;
};

struct ShellLayerState {
  Vec3 sigma = Vec3::Zero();     // (sx, sy, txy) local frame [Pa]
  Vec3 eps_p = Vec3::Zero();     // plastic strain (engineering shear)
  Eigen::Vector2d tau = Eigen::Vector2d::Zero();      // transverse shear stresses (txz, tyz) [Pa]
  Eigen::Vector2d gamma_p = Eigen::Vector2d::Zero();  // their plastic (engineering) strains
  double alpha = 0;              // equivalent plastic strain
};

struct ShellPlasticStep {
  double lambda = 0;
  std::vector<Vec6> node_disp;  // per model node, global
  std::vector<Vec12> bar_end_forces;  // per bar, local end forces (k u - f0; nonlinear bars: a6^T q)
  std::vector<Vec6> bar_basic_forces;  // per bar, basic forces q (N, My_i, My_j, Mz_i, Mz_j, T); zero for elastic bars
  int iterations = 0;
};

struct ShellPlasticResult {
  std::vector<ShellPlasticStep> steps;
  // final state: per shell, per in-plane point (the 3 mid-edge points; the 7
  // Dunavant points of a facet whose thickness varies), per thickness point
  // (bottom z = -t/2 to top z = +t/2)
  std::vector<std::vector<std::vector<ShellLayerState>>> state;
  std::vector<double> z;  // unit thickness coordinates of the points (-1/2 .. 1/2)
  // (an orthotropic facet or a laminate integrates each ply on its own:
  // its points are per_shell_z[shell], z / t from the nodes' plane)
  std::vector<std::vector<double>> per_shell_z;
  bool limit_reached = false;  // stop_at_limit: the last increment found no equilibrium
  std::string limit_message;
};

ShellPlasticResult solve_shell_plastic(const Model& model, Id case_id, const ShellPlasticOptions& opt = {});

// Plane-stress J2 return mapping (exposed for direct use): from the committed
// state and the total strain, the updated state and the consistent tangent.
ShellLayerState j2_plane_stress(const ShellLayerState& committed, const Vec3& eps, double E, double nu,
                                double fy, double H, Mat3* tangent);
// The same with the transverse shears: eps = (ex, ey, gxy, gxz, gyz), sz = 0.
// C and the J2 projector still share the eigenvectors (sx + sy, sx - sy and
// the three shears, which all scale by 1 / (1 + 2 G dg)), so the return is
// the same closed form with A2 = (sx - sy)^2 / 2 + 2 (txy^2 + txz^2 + tyz^2).
using Mat5 = Eigen::Matrix<double, 5, 5>;
using Vec5 = Eigen::Matrix<double, 5, 1>;
ShellLayerState j2_shell(const ShellLayerState& committed, const Vec5& eps, double E, double nu, double fy,
                         double H, Mat5* tangent);
// Hill: the projector P (1/2 s' P s = R^2 / 3) in the facet's local frame for
// material axis 1 at `angle` [rad] from the local x — (sx, sy, txy, txz, tyz)
// — and the general return with any symmetric P (isotropic elasticity)
Mat5 hill_projector(const HillYield& h, double angle);
ShellLayerState hill_shell(const ShellLayerState& committed, const Vec5& eps, double E, double nu, double fy,
                           double H, const Mat5& P, Mat5* tangent);
// The same with ANY elastic matrix C (orthotropic facets, laminate plies:
// the local-axes C of the ply; the strain and stress orders as above) and
// the J2 projector of (sx, sy, txy, txz, tyz) for a von Mises ply
ShellLayerState general_shell(const ShellLayerState& committed, const Vec5& eps, const Mat5& C, double fy, double H,
                              const Mat5& P, Mat5* tangent);
Mat5 j2_projector5();

}  // namespace keta

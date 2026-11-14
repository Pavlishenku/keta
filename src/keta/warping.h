// Non-uniform (Vlasov) torsion of thin-walled bars: a SEVENTH degree of
// freedom per node, the warping psi = d(theta_x)/dx, carried by a dedicated
// driver over an ordinary keta Model (the 12-DOF bar of the rest of keta stays
// Saint-Venant). Scope: doubly-symmetric sections (shear centre = centroid,
// no Wagner asymmetry beta), straight members, every bar meeting a node with
// a free warping DOF collinear with the others there (psi is invariant to the
// bar direction, so one scalar per node expresses warping continuity).
//
// Element (local DOFs: the 12 of keta's bar, then psi_i, psi_j):
//   * axial + bending: keta's member block (Timoshenko in first order; the
//     exact stability functions under N in second order, Euler-Bernoulli);
//   * torsion: the EXACT solution of E Iw theta'''' - G_t theta'' = m_x,
//     G_t = G J (first order) or G J + N r0^2 (second order, r0^2 =
//     (Iy + Iz)/A the polar radius about the shear centre — tension
//     stiffens), homogeneous basis {1, x, F2, F3} of the transfer functions
//     (cosh/sinh, cos/sin, Taylor series near G_t = 0), so ONE element per
//     member is exact in pure torsion; fixed-end torques/bimoments of
//     distributed and point torques exact in the same basis.
//   * flexural-torsional coupling (second order and buckling): the
//     second-order energy of the thin-walled bar (Trahair, Flexural-Torsional
//     Buckling of Structures, 1993, ch. 3; doubly-symmetric):
//         U2 = int [ N/2 (v'^2 + w'^2 + r0^2 theta'^2) + My theta v'' + Mz theta w'' ] dx
//              + 1/2 sum (F . e) theta^2          (load at offset e from the shear centre)
//     with keta's moment signs (+My tends the +z fibre, +Mz compresses +y),
//     consistent cubic Hermite interpolation of v, w, theta, 5-point Gauss
//     (exact for the quadratic moment of a uniform load). This coupling is a
//     DISCRETIZATION (converges as h^2..h^4 under meshing) — the uncoupled
//     blocks are exact.
// Results: 7 displacements per node; local member end forces (12 + the two
// bimoments, B = -E Iw theta''); the torque split at the
// ends: primary M_T1 = G J psi, Wagner M_T3 = N r0^2 psi (second order),
// secondary M_T2 = M_T - M_T1 - M_T3.
//
// Typed rejections: bars with releases, rigid-arm eccentricity, Winkler bed,
// CARA VAR, cables/trusses, settlements, shells, non-collinear bars at a node
// with free warping, |G_t| L^2 / (E Iw) > 400 on one element (split it),
// Iw = 0 on a bar whose node warping DOF is free and touched by no other bar
// with Iw > 0 (no stiffness: the DOF is then dropped — Saint-Venant).
#pragma once
#include <map>
#include <set>

#include "keta/statics.h"

namespace keta {

struct WarpingSpec {
  std::map<Id, double> Iw;           // warping constant C_M per SECTION id [m^6] (absent = 0)
  std::set<Id> warping_fixed;        // nodes whose warping psi is restrained (FIX ... MB)
  // Span loads of these bars act at the offset (ey, ez) [m] from the shear
  // centre in the section's local axes: linear torque
  // ey qz - ez qy, and the load-height term (q . e) theta^2 / 2.
  std::map<Id, Eigen::Vector2d> load_offset;
};

enum class WarpingTheory { First, Second };

struct WarpingOptions {
  WarpingTheory theory = WarpingTheory::First;
  int max_iterations = 50;
  double tol_rel = 1e-10;
};

struct WarpingBarResult {
  Id bar = 0;
  Eigen::Matrix<double, 14, 1> end_forces = Eigen::Matrix<double, 14, 1>::Zero();  // local, on the member
  double N = 0;                      // axial force (tension +)
  // torque split at end i (index 0) and end j (index 1), internal-force
  // convention: M_T = -f(3) at i, +f(9) at j; bimoment M_omega = -E Iw theta''
  std::array<double, 2> MT{}, MT1{}, MT2{}, MT3{}, Mw{};
};

struct WarpingResults {
  std::vector<Eigen::Matrix<double, 7, 1>> node_disp;  // ux uy uz rx ry rz psi (psi per unit length)
  std::vector<WarpingBarResult> bars;
  std::vector<std::pair<Id, Eigen::Matrix<double, 7, 1>>> reactions;  // supported nodes
  int iterations = 0;
};

WarpingResults solve_warping(const Model& model, Id case_id, const WarpingSpec& spec,
                             const WarpingOptions& opt = {});

struct WarpingBucklingResults {
  std::vector<double> load_factors;  // ascending positive
  std::vector<std::vector<Eigen::Matrix<double, 7, 1>>> modes;  // unit crest
};
// Linear flexural-torsional buckling about the first-order state of the case:
// (K + lambda Kg) x = 0, Kg from N, My, Mz and the load heights.
WarpingBucklingResults warping_buckling(const Model& model, Id case_id, const WarpingSpec& spec,
                                        int n_modes = 3);

}  // namespace keta

// Third-order theory WITH WARPING: the
// geometrically exact corotational bar of keta (corot_kin.h: exact finite
// rotations, chord + mean-director frame, rigid arms) given the SEVENTH
// degree of freedom per node of the Vlasov theory, the warping psi = the
// rate of twist (keta/warping.h). The element energy is
//     W = 1/2 d^T Kbar d   (axial, bending: the 6-DOF corotational element)
//       + 1/2 t^T K_tw t,  t = (theta_1, psi_1, theta_2, psi_2),
// theta_k the corotational twists of the ends (rigid twist of the frame
// removed), K_tw = E Iw int N''^2 + G J int N'^2 on the cubic Hermite twist
// (Vlasov torsion), its gradient and Hessian by forward automatic
// differentiation (exact derivatives). psi is a scalar per node (invariant to
// the bar direction; collinear bars at a node with free warping, as the
// linear driver). Load stepping and Newton on the total residual.
//
// Scope: frame bars (rigid arms allowed — a load height is an arm), nodal
// forces (dead loads); warping restrained at WarpingSpec::warping_fixed and
// at the nodes no bar with Iw > 0 reaches. Typed rejections: releases,
// Winkler beds, CARA VAR, cables, nodal moments, span loads, settlements,
// shells, links. A linear local element in a corotational frame: the
// flexural-torsional coupling of the third-order theory comes from the
// frame's rotations and converges under mesh refinement (no Wagner term:
// valid without large axial force).
#pragma once
#include "keta/nonlinear.h"
#include "keta/warping.h"

namespace keta {

struct CorotWarpingResults {
  std::vector<Eigen::Matrix<double, 7, 1>> node_disp;  // translations, rotation vector, psi
  std::vector<Vec12> bar_forces;                       // local deformational end forces (as CorotResults)
  std::vector<Eigen::Vector2d> bimoments;              // per bar: B at the ends, B = -E Iw theta''
  int iterations = 0;
};
CorotWarpingResults solve_corotational_warping(const Model& model, Id case_id, const WarpingSpec& spec,
                                               const CorotOptions& opt = {});

}  // namespace keta

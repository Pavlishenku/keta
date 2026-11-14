// Modal analysis of a structure spinning at a constant rate about a fixed
// axis (rotating blades, rotors): the centrifugal prestress and the spin
// softening (a "load correction stiffness"), in the rotating
// frame, small vibrations about the steady state:
//   (K + Kg(N_c) - Omega^2 M_perp) phi = omega^2 M phi,
// * the centrifugal loads are EXACT: rho A Omega^2 P (x - p) along each bar,
//   linear in the abscissa (P = I - e e^t the projection onto the plane of
//   rotation, p a point of the axis), m Omega^2 P (x - p) on point masses;
//   the steady state is the linear static solution under them;
// * Kg integrates the axial force N(x) of the steady state along each bar
//   (quadratic in x: N_i minus the integral of the axial load), with the
//   bar's own interpolation (Hermite transverse, linear truss), by Gauss
//   quadrature exact for it — no constant-N-per-bar approximation;
// * M_perp = int rho A S^t P S dx: the mass seen by the centrifugal field of
//   the displacement (S the interpolation of the bar's consistent mass:
//   linear axial, Hermite transverse; with a lumped mass, half the bar at
//   each end), point masses m P; Coriolis
//   (gyroscopic) coupling neglected.
// Rejected: shells, end releases, rigid arms, tapered bars, rigid links /
// constraints, a spin beyond divergence (K_eff not positive definite).
#pragma once
#include "keta/eigsolve.h"

namespace keta {

struct Spin {
  Vec3 point = Vec3::Zero();      // a point of the axis
  Vec3 axis = Vec3(1, 0, 0);      // its direction (any length)
  double omega = 0;               // [rad/s]
  bool stiffening = true;         // the centrifugal prestress Kg
  bool softening = true;          // -Omega^2 M_perp
};

struct RotatingModal {
  ModalResults modal;
  CaseResult steady;              // the centrifugal steady state
};

RotatingModal solve_modal_rotating(const Model& model, const Spin& spin, const ModalOptions& opt = {});

}  // namespace keta

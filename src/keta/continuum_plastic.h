// Plasticity of continuum elements (keta/continuum.h): von Mises (J2) with
// linear isotropic hardening, small strains, plane strain, axisymmetric and
// solid elements (plane stress is a typed rejection for a plastic material:
// its return mapping is another algorithm).
//
// Every Gauss point carries the full 3-D state (plastic strain, equivalent
// plastic strain): the elastic predictor sigma = D (eps - eps_p), and if
// ||s|| > sqrt(2/3) (sigma_y + H ebar) the radial return (Simo & Hughes,
// box 3.2) dgamma = f / (2G + 2H/3), exact for linear hardening, with the
// algorithmic tangent K 1x1 + 2G beta P - 2G gbar n x n (beta = 1 - 2G
// dgamma / ||s_tr||, gbar = 1 / (1 + H / 3G) - (1 - beta)): Newton converges
// quadratically on every load step. Elastic materials (not listed) stay
// elastic. Loads: the case times the load factors given (each solved from
// the committed state of the previous one, in substeps).
#pragma once
#include <map>

#include "keta/continuum.h"

namespace keta {

struct J2Material {
  double sigma_y = 0;
  double H = 0;  // isotropic hardening modulus (slope of the yield stress against ebar)
  // Appended: power-law creep, time hardening, eps_c = A q^n t^m (q the von
  // Mises stress, t the time from the start of the analysis; A = 0: none).
  // Deviatoric flow along the J2 normal, integrated implicitly (backward
  // Euler on the stress: d eps_c = A q1^n (t1^m - t0^m)) together with the
  // plastic return, so the step stays unconditionally stable; first order
  // in the time step under a varying stress, exact under a constant one.
  // sigma_y = 0 with A > 0: viscous only (no yield surface).
  double creep_A = 0, creep_n = 1, creep_m = 1;
  // Appended: thermally activated creep (Norton-Arrhenius): the creep strain
  // is multiplied by exp(-creep_Q / T), T the temperature at the Gauss point
  // interpolated from ContinuumPlasticOptions::temperature (absolute scale,
  // T > 0; constant in time). 0: no temperature dependence.
  double creep_Q = 0;
};
// FRICTIONAL (soil) plasticity (appended), tension positive, small strains,
// perfect or linear cohesion hardening c = c0 + H ebar:
//  * Drucker-Prager f = sqrt(J2) + eta p - xi c, plastic potential with
//    eta_bar (dilatancy psi); eta, xi from phi by the outer cone (through
//    the compressive meridian of Mohr-Coulomb), the inner one (tensile
//    meridian) or the plane-strain match (eta = 3 tan phi / sqrt(9 + 12
//    tan^2 phi), xi = 3 / sqrt(...)): the closed-form return of the smooth
//    cone dgamma = f / (G + K eta eta_bar + xi^2 H), or of the apex;
//  * Mohr-Coulomb (s1 >= s2 >= s3) (s1 - s3) + (s1 + s3) sin phi = 2 c cos
//    phi, potential with psi: the principal-space return of de Souza Neto,
//    Peric & Owen (2008, box 8.4) — main plane, right / left edge (two
//    planes) or apex — each a LINEAR system A dgamma = f_trial for linear
//    hardening (exact), ebar = 2 cos phi sum dgamma (cot phi c at the apex).
// Consistent tangents: principal-space D_p = E - E N A^-1 M^T E (planes),
// the apex's, plus the spin terms (s_i - s_j) / (e_i - e_j) of the
// eigenvector rotation; the smooth DP tangent in closed form. psi != phi:
// non-symmetric tangent (sparse LU).
struct FrictionalMaterial {
  enum class Kind { DruckerPrager, MohrCoulomb };
  enum class DpFit { Outer, Inner, PlaneStrain };
  Kind kind = Kind::MohrCoulomb;
  double c = 0, phi = 0, psi = 0;  // cohesion [Pa], friction and dilatancy angles [rad]
  double H = 0;                    // cohesion hardening modulus (against ebar)
  DpFit fit = DpFit::Outer;
};
// One material-point update (exposed for direct use): from the total strain
// (xx, yy, zz, gxy, gyz, gzx) and the committed plastic strain / ebar, the
// stress, the consistent tangent, the new plastic strain and ebar; mode 0
// elastic, 1 smooth / main plane, 2 edge, 3 apex.
struct FrictionalPoint {
  Vec6 sigma = Vec6::Zero(), ep = Vec6::Zero();
  Mat6 D = Mat6::Zero();
  double ebar = 0;
  int mode = 0;
};
FrictionalPoint frictional_update(const FrictionalMaterial& m, double E, double nu, const Vec6& eps, const Vec6& ep,
                                  double ebar);

struct ContinuumPlasticOptions {
  std::vector<double> lambdas{1.0};   // load factors of the case, in order
  int substeps = 4;                   // increments from one factor to the next
  int max_iter = 60;
  double tol = 1e-10;                 // |r| / |f|
  std::map<Id, J2Material> plastic;   // material id -> J2 law
  // Appended: total Lagrangian large displacements (Green-Lagrange strain,
  // St Venant-Kirchhoff elasticity: small strains, large rotations) for
  // plane stress, plane strain and solid elements — elastic materials only;
  // the loads keep their direction (dead loads)
  bool large_displacements = false;
  // Appended: the time of each load factor (same size as lambdas; the
  // analysis starts at t = 0, substeps interpolate linearly) — required
  // when a material creeps; a repeated factor with a later time is a hold
  std::vector<double> times;
  // Appended: frictional materials (material id -> law); a material is in
  // `plastic` or here, not both
  std::map<Id, FrictionalMaterial> frictional;
  // a load step without equilibrium is halved from the last converged
  // state up to limit_cuts times, then the analysis returns its converged
  // factors (limit_reached) instead of raising
  bool stop_at_limit = false;
  int limit_cuts = 6;
  // Appended: nodal temperatures (node -> T) for materials with creep_Q > 0
  std::map<Id, double> temperature;
};
struct ContinuumPlasticResult {
  std::vector<std::vector<Vec6>> node_disp;  // per lambda, per model node
  // final state, per element, per stiffness Gauss point: stress (xx, yy, zz,
  // xy, yz, zx) — axisymmetric (rr, zz, tt, rz, 0, 0) — and the equivalent
  // plastic strain
  std::vector<std::vector<std::array<double, 7>>> gauss;
  int iterations = 0;
  // Appended: per converged load factor, the support reactions (per model
  // node, global: the internal forces on the blocked translations; loads
  // applied on blocked DOFs are not subtracted) and the factor itself;
  // prescribed support displacements of the case scale with the factor
  std::vector<std::vector<Vec6>> reactions;
  std::vector<double> lambdas_done;
  bool limit_reached = false;
};
ContinuumPlasticResult solve_continuum_plastic(const Model& model, Id case_id, const ContinuumPlasticOptions& opt);

}  // namespace keta

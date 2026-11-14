// Equivalent static loads and concomitant effects of a spectral analysis.
//
// Load-response correlation (Kasperski): for loads F of covariance C_F and a
// target response r = I . F, the load F_e = g C_F I / sqrt(I^T C_F I) gives
// r(F_e) = g sigma_r EXACTLY (I^T C_F I = sigma_r^2) and is the most probable
// load pattern producing it. In modal space (q of covariance C_q, r = psi . q)
// q_e = g C_q psi / sqrt(psi^T C_q psi) and the static nodal loads that
// produce q_e are the inertia loads F = M Phi Omega^2 q_e (K^-1 M phi_i
// omega_i^2 = phi_i): background, resonant and cross-modal parts together, the
// target's peak recovered by a linear static solve. A single mode: the
// inertia load g sigma_qi omega_i^2 M phi_i.
//
// Concomitance (Capra & Davidovici, Calcul dynamique des structures en zone
// sismique, 1980): effects E of mean m, covariance Sigma, peak factor g — the
// equiprobability ellipse(oid) (E - m)^T Sigma^-1 (E - m) = g^2 is replaced
// by the vertices of a circumscribing polygon / polyhedron. Two effects: the
// octagon whose 8 sides are the tangents x/s_a +- y/s_b = +-sqrt(2(1 +- rho))
// cut by the sigma box |E_a| <= g s_a; three effects: that octagon in each
// coordinate plane (24 vertices), the third coordinate the conditional mean
// Sigma_c,ab Sigma_ab^-1 (E_a, E_b). A rank-deficient Sigma reduces in its
// principal frame (rank 2: the tangent octagon of the ellipse, rank 1: m +-
// g sigma, rank 0: m).
#pragma once
#include <vector>

#include "keta/wind_buffeting.h"

namespace keta {

Eigen::VectorXd lrc_loads(const Eigen::MatrixXd& C_F, const Eigen::VectorXd& I, double g);
Eigen::VectorXd lrc_modal(const Eigen::MatrixXd& C_q, const Eigen::VectorXd& psi, double g);
// nodal inertia loads M Phi Omega^2 q of modal amplitudes q, as a load case
LoadCase modal_inertia_case(const Model& model, const ModalResults& modal, const Eigen::VectorXd& q, Id case_id,
                            const std::string& name);

// vertices of the circumscribing polygon (k = 2) / polyhedron (k = 3)
std::vector<Eigen::VectorXd> equiprobable_polyhedron(const Eigen::VectorXd& m, const Eigen::MatrixXd& Sigma, double g);

// The nonlinear wind loop, at the second-order level (it does not run on a
// saved tangent state of keta/tangent.h): permanent + mean wind -> second-order static -> axial
// field N -> modes about N (ModalOptions::axial_state) -> spectral buffeting
// -> peak equivalent loads (modal LRC on the target) -> second-order static
// under permanent + mean + peak -> new N, ... until the target's sigma and
// extreme change by less than tol. Without axial force anywhere it is the
// plain spectral analysis after one pass.
struct WindLoopResult {
  double sigma = 0, peak_factor = 0;
  double extreme = 0;              // target under permanent + mean + LRC peak (second order)
  int iterations = 0;
  bool converged = false;
  std::vector<double> axial;       // per bar, the last state
  std::vector<double> frequencies_hz;
};
WindLoopResult wind_nonlinear_loop(const Model& model, Id permanent_case, const BuffetingSpec& spec,
                                   const BuffetingRequest& target, int n_modes = 6, int max_iterations = 20,
                                   double tol = 1e-6);

}  // namespace keta

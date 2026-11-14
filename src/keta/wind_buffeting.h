// Multi-component spectral buffeting.
//
// Strip forces per unit length, quasi-steady (Davenport / Scanlan), q =
// rho U^2 / 2, coefficients and slopes at the mean incidence:
//   D = q B [2 C_D u / U + (C_D' - C_L) w / U]
//   L = q B [2 C_L u / U + (C_L' + C_D) w / U]
//   M = q B^2 [2 C_M u / U + C_M' w / U]
// each force row times its admittance chi_r (amplitude; chi_r chi_s on the
// cross terms). Drag along the wind, lift vertical, moment about a = z x d
// (positive nose-up for the wind d). The modal forces' cross-spectral matrix
// S_Q(n) is built per velocity component (u, w uncorrelated) from the
// cross-spectral velocity matrix, O(m N^2) per frequency.
//
// Response (unit modal masses): S_q = H S_Q H^*, H(omega) = [Omega^2 + K_ae
// - omega^2 + i omega (2 zeta Omega + C_ae)]^-1 — the cross-modal (CQC)
// terms come out of the spectra, not of a combination rule. A response r =
// psi . q: sigma^2 = integral of psi^T S_q psi dn on a log grid refined
// around every modal frequency (the resonant peaks, width zeta n_i, are
// resolved), BACKGROUND sigma_B^2 = the same with H(0) (quasi-static),
// RESONANT sigma_R^2 = sigma^2 - sigma_B^2; up-crossing rate nu = sqrt(int
// n^2 S / int S), Davenport peak factor sqrt(2 ln nu T) + 0.5772 / sqrt(2 ln
// nu T).
//
// Aeroelastic (quasi-steady, coupled): the motion of the section
// changes the relative wind, u -> -x_D', w -> -x_L', so the self-excited
// force is the buffeting map with those velocities — a modal damping matrix
// C_ae (along-wind drag damping rho U B C_D included) — and the rotation
// theta changes the incidence, an aerodynamic stiffness K_ae = -q B (C_D',
// C_L', B C_M') theta. Without aero terms the result is the buffeting one.
#pragma once
#include <functional>
#include <vector>

#include "keta/eigsolve.h"
#include "keta/wind_turbulence.h"

namespace keta {

// Coefficients as functions of the incidence alpha [rad]: a table
// linear between the rows, the slopes those of the segment (the mean of the
// two segments at a row), held flat outside (slopes 0 there).
struct AeroTable {
  std::vector<double> alpha, cd, cl, cm;
  static AeroTable from_section(const AeroSection& s);  // linear: C(alpha) = C0 + C' alpha
  struct Value {
    double cd = 0, cl = 0, cm = 0, dcd = 0, dcl = 0, dcm = 0;
  };
  Value at(double a) const;
};

// ---- stage 1: modal response to a modal force spectrum ----

struct ModalSystem {
  std::vector<double> omega, zeta;  // per mode
  Eigen::MatrixXd C_ae, K_ae;       // aeroelastic modal matrices (empty = none)
};
struct ResponseStats {
  double sigma = 0, sigma_b = 0, sigma_r = 0, nu = 0, peak_factor = 0;
};
// the frequency grid: log over [n_min, n_max] plus `per_mode` points over
// n_i (1 +- 25 zeta_i) for every mode
std::vector<double> buffeting_grid(const ModalSystem& sys, double n_min, double n_max, int n_log, int per_mode);
// S_Q(n): m x m real symmetric modal force cross-PSD; psi: one row per response
std::vector<ResponseStats> modal_response(const ModalSystem& sys, const std::vector<double>& grid,
                                          const std::function<Eigen::MatrixXd(double)>& S_Q,
                                          const std::vector<Eigen::VectorXd>& psi, double duration = 600);
// the modal displacement covariance matrix integral of S_q dn (m x m)
Eigen::MatrixXd modal_covariance(const ModalSystem& sys, const std::vector<double>& grid,
                                 const std::function<Eigen::MatrixXd(double)>& S_Q);

// ---- stage 2: strips of a model ----

struct BuffetingSpec {
  TurbulenceModel turb;
  Vec3 wind_dir = Vec3(1, 0, 0);  // horizontal
  double B = 1;                   // reference width [m]
  AeroTable aero;                 // coefficients vs incidence
  double alpha0 = 0;              // mean incidence [rad]
  std::array<Admittance, 3> chi;  // per force row D, L, M
  bool aeroelastic = false;
  DampingSpec damping;            // structural
  int n_log = 300, per_mode = 120;
  double duration = 600;
  int threads = 0;
};
struct Strip {
  WindPoint wp;
  double L = 0;
  std::vector<Vec3> phi;  // per mode: (drag, lift, rotation about a) at the midpoint
};
std::vector<Strip> buffeting_strips(const Model& model, const ModalResults& modal, const BuffetingSpec& spec);
ModalSystem buffeting_system(const ModalResults& modal, const std::vector<Strip>& strips, const BuffetingSpec& spec);
Eigen::MatrixXd modal_force_psd(const std::vector<Strip>& strips, const BuffetingSpec& spec, double n);

struct BuffetingRequest {
  Id node = 0;
  int comp = 0;  // 0..5
};
struct MultiBuffetingResults {
  std::vector<ResponseStats> stats;  // per request
  Eigen::MatrixXd modal_cov;         // the CQC covariance of the modal coordinates
  ModalSystem system;
};
MultiBuffetingResults spectral_buffeting_multi(const Model& model, const ModalResults& modal, const BuffetingSpec& spec,
                                               const std::vector<BuffetingRequest>& requests);

// The mean wind load case of the buffeting model: q(z) = rho U(z)^2 / 2 (the
// MEAN pressure, not the peak one of mean_wind_case), D, L, M of the table at
// the mean incidence, per bar at its mid height.
LoadCase mean_buffeting_case(const Model& model, const BuffetingSpec& spec, Id case_id, const std::string& name);

}  // namespace keta

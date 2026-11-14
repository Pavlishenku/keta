// Turbulent wind field synthesis and time-domain response.
//
// Spectral representation with Deodatis' double-indexed frequencies: for N
// points and a velocity component, H(n) the Cholesky factor of the N x N
// cross-spectral matrix (S = H H^*), the series
//   v_j(t) = sum_{m <= j} sum_l sqrt(2 dn) |H_jm(n_ml)| cos(2 pi n_ml t - theta_jm + phi_ml),
//   n_ml = (l - 1) dn + m dn / N,   phi_ml uniform on [0, 2 pi) (seeded),
// evaluated by one FFT per point over the period T = N / dn. Every frequency
// is distinct, so over one period the cosines are orthogonal: the time
// variance is sum_ml dn |H_jm|^2 and the time covariance of two points
// sum_ml dn Re(H_jm conj(H_km)) EXACTLY — the discrete target spectrum, not
// a statistical estimate.
// The cut-off n_max drops the spectral tail above it.
//
// Time-domain modal response to force histories: each mode (or the coupled
// aeroelastic modal system) integrated by Newmark's average acceleration,
// unconditionally stable, second order; the first `transient` seconds left out
// of the statistics.
#pragma once
#include <cstdint>
#include <vector>

#include "keta/wind_buffeting.h"

namespace keta {

// in-place radix-2 FFT (size a power of two), sign -1 forward, +1 inverse
// (unscaled)
void fft(std::vector<std::complex<double>>& a, int sign);

struct WindSeries {
  double dt = 0;
  int nt = 0;
  // series[j][c]: point j, component c (0 u, 1 w), nt samples over one period
  std::vector<std::array<std::vector<double>, 2>> series;
  // the exact one-period variance / covariance the synthesis carries:
  // cov[c](j, k) = sum_ml dn Re(H_jm conj(H_km))
  std::array<Eigen::MatrixXd, 2> cov;
};
WindSeries synthesize_wind(const TurbulenceModel& m, const std::vector<WindPoint>& pts, double n_max, int n_freq,
                           std::uint64_t seed);

// Newmark (average acceleration) on the modal system: q'' + (2 zeta Omega +
// C_ae) q' + (Omega^2 + K_ae) q = Q(t), Q given at dt (m rows, nt columns)
Eigen::MatrixXd modal_time_response(const ModalSystem& sys, const Eigen::MatrixXd& Q, double dt);

// The time-domain buffeting of a model: the wind synthesized at the strip
// midpoints (u, w), the quasi-steady strip forces of wind_buffeting.h (no
// admittance: a frequency-domain filter), the modal forces, Newmark on the
// modal system (aeroelastic terms when spec.aeroelastic), the requested
// nodal responses over `periods` periods of the synthesis, the first one
// discarded as transient. Statistics over the kept periods.
struct TemporalResponse {
  double sigma = 0, max = 0, min = 0;
  std::vector<double> series;  // kept samples
};
struct TemporalBuffeting {
  double dt = 0;
  std::vector<TemporalResponse> responses;  // per request
};
TemporalBuffeting temporal_buffeting(const Model& model, const ModalResults& modal, const BuffetingSpec& spec,
                                     const std::vector<BuffetingRequest>& requests, double n_max, int n_freq,
                                     std::uint64_t seed, int periods = 2);

}  // namespace keta

#include "keta/wind_turbulence.h"

#include <algorithm>
#include <cmath>

namespace keta {

double reduced_spectrum(SpectrumForm form, int a, double f) {
  require(a >= 0 && a <= 2, "spectrum: component 0 (u), 1 (v) or 2 (w)");
  require(std::isfinite(f) && f >= 0, "spectrum: reduced frequency >= 0");
  if (form == SpectrumForm::EN1991) return 6.8 * f / std::pow(1 + 10.2 * f, 5.0 / 3.0);
  if (a == 0) return 4 * f / std::pow(1 + 70.8 * f * f, 5.0 / 6.0);
  return 4 * f * (1 + 755.2 * f * f) / std::pow(1 + 283.2 * f * f, 11.0 / 6.0);
}

double auto_spectrum(SpectrumForm form, int a, double n, double U, double sigma, double L) {
  require(n > 0 && U > 0 && sigma >= 0 && L > 0, "spectrum: n > 0, U > 0, sigma >= 0, L > 0");
  return sigma * sigma * reduced_spectrum(form, a, n * L / U) / n;
}

double TurbulenceModel::length_scale_u(double z) const {
  require(profile.z0 > 0, "turbulence: z0 > 0");
  const double alpha = 0.67 + 0.05 * std::log(profile.z0);
  return 300 * std::pow(std::max(z, profile.zmin) / 200, alpha);
}

WindPoint wind_point(const TurbulenceModel& m, double y, double z) {
  require(std::isfinite(y) && std::isfinite(z), "wind point: finite position");
  WindPoint p;
  p.y = y;
  p.z = z;
  p.U = m.profile.vm(z);
  require(p.U > 0, "wind point: the mean speed must be positive");
  const double su = m.profile.iv(z) * p.U;
  p.sigma = {su, m.rv * su, m.rw * su};
  const double Lu = m.length_scale_u(z);
  p.L = {Lu, m.lv * Lu, m.lw * Lu};
  return p;
}

double davenport_coherence(double n, double dy, double dz, double Uj, double Uk, double Cy, double Cz) {
  require(n >= 0 && Uj > 0 && Uk > 0 && Cy >= 0 && Cz >= 0, "coherence: n >= 0, U > 0, C >= 0");
  return std::exp(-2 * n * std::hypot(Cy * dy, Cz * dz) / (Uj + Uk));
}

Eigen::MatrixXd cross_spectral_matrix(const TurbulenceModel& m, const std::vector<WindPoint>& pts, double n) {
  const std::size_t N = pts.size();
  Eigen::MatrixXd S = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(3 * N), static_cast<Eigen::Index>(3 * N));
  std::vector<std::array<double, 3>> sq(N);
  for (std::size_t j = 0; j < N; ++j)
    for (int a = 0; a < 3; ++a)
      sq[j][a] = pts[j].sigma[a] > 0 ? std::sqrt(auto_spectrum(m.form, a, n, pts[j].U, pts[j].sigma[a], pts[j].L[a])) : 0;
  for (std::size_t j = 0; j < N; ++j)
    for (std::size_t k = j; k < N; ++k)
      for (int a = 0; a < 3; ++a) {
        const double c = k == j ? 1.0
                                : davenport_coherence(n, pts[j].y - pts[k].y, pts[j].z - pts[k].z, pts[j].U, pts[k].U,
                                                      m.C[a][0], m.C[a][1]);
        const double v = sq[j][a] * sq[k][a] * c;
        S(static_cast<Eigen::Index>(3 * j + a), static_cast<Eigen::Index>(3 * k + a)) = v;
        S(static_cast<Eigen::Index>(3 * k + a), static_cast<Eigen::Index>(3 * j + a)) = v;
      }
  return S;
}

// ---------------- admittance ----------------

std::complex<double> theodorsen_exact(double K) {
  require(K > 0, "Theodorsen: K > 0");
  const std::complex<double> H0(std::cyl_bessel_j(0.0, K), -std::cyl_neumann(0.0, K));
  const std::complex<double> H1(std::cyl_bessel_j(1.0, K), -std::cyl_neumann(1.0, K));
  return H1 / (H1 + std::complex<double>(0, 1) * H0);
}

std::complex<double> sears_function(double K) {
  require(K > 0, "Sears: K > 0");
  const double J0 = std::cyl_bessel_j(0.0, K), J1 = std::cyl_bessel_j(1.0, K);
  const std::complex<double> i(0, 1);
  return (J0 - i * J1) * theodorsen_exact(K) + i * J1;
}

double Admittance::chi2(double K) const {
  require(K >= 0, "admittance: K >= 0");
  switch (kind) {
    case Kind::QuasiSteady: return 1;
    case Kind::Sears: return K == 0 ? 1.0 : std::norm(sears_function(K));
    case Kind::Liepmann: return 1 / (1 + 2 * kPi * K);
    case Kind::Table: {
      require(!table.empty(), "admittance: empty table");
      if (K <= table.front().first) return table.front().second;
      if (K >= table.back().first) return table.back().second;
      const auto it = std::upper_bound(table.begin(), table.end(), K, [](double x, const auto& p) { return x < p.first; });
      const auto& [k1, v1] = *it;
      const auto& [k0, v0] = *(it - 1);
      return v0 + (K - k0) / (k1 - k0) * (v1 - v0);
    }
  }
  return 1;
}

double Admittance::chi2_at(double n, double B, double U) const {
  require(n >= 0 && B > 0 && U > 0, "admittance: n >= 0, B > 0, U > 0");
  return chi2(kPi * n * B / U);  // K = omega (B/2) / U
}

}  // namespace keta

// Atmospheric turbulence for buffeting.
//
// Spectra (one-sided, per Hz), every form NORMALISED: the integral of S_a
// over n is sigma_a^2 exactly — the property a spectrum must have and that
// the u*-normalised Kaimal shapes (200 f / (1 + 50 f)^(5/3) and the like)
// do not: used with sigma^2 they carry 6 to 18 times the variance:
//   EN 1991-1-4 Annex B (Kaimal type, on the integral scale L_a of the
//   component): n S / sigma^2 = 6.8 f / (1 + 10.2 f)^(5/3), f = n L_a / U;
//   von Karman: u 4 f / (1 + 70.8 f^2)^(5/6), v and w (transverse)
//   4 f (1 + 755.2 f^2) / (1 + 283.2 f^2)^(11/6).
// Davenport coherence per component a: exp(-n sqrt((C_ay dy)^2 + (C_az
// dz)^2) / U_jk), U_jk the mean of the two speeds; the three components
// uncorrelated with each other (the u-w Reynolds stress neglected — the
// usual buffeting assumption).
// Variation with height: U(z), I_u(z) from the EN 1991-1-4 profile, I_v =
// r_v I_u, I_w = r_w I_u; L_u(z) = 300 (z / 200)^alpha, alpha = 0.67 + 0.05
// ln z0, below z_min held at L(z_min) (EN 1991-1-4 (B.1)); L_v, L_w as
// ratios of L_u. EN gives u only: the v / w ratios are the user's (defaults
// 0.75 / 0.5 for the intensities, 1/4 and 1/12 for the length scales, ESDU-
// type orders of magnitude).
//
// Aerodynamic admittance chi^2(K), K = omega B / (2 U) the reduced frequency
// on the half-chord: quasi-steady (1), Sears' exact function |S(K)|^2 with
// S = [J0 - i J1] C(K) + i J1 and Theodorsen's C(K) = H1(K) / (H1(K) + i
// H0(K)) (Hankel functions of the second kind), Liepmann's 1 / (1 + 2 pi K),
// or a table (K, chi^2) linear between, held flat outside.
#pragma once
#include <array>
#include <complex>
#include <vector>

#include "keta/wind.h"

namespace keta {

enum class SpectrumForm { EN1991, VonKarman };

// n S_a(n) / sigma_a^2 at reduced frequency f = n L_a / U (component a = 0 u, 1 v, 2 w)
double reduced_spectrum(SpectrumForm form, int component, double f);
// S_a(n) [ (m/s)^2 / Hz ]
double auto_spectrum(SpectrumForm form, int component, double n, double U, double sigma, double L);

struct TurbulenceModel {
  WindProfile profile;
  SpectrumForm form = SpectrumForm::EN1991;
  double rv = 0.75, rw = 0.5;              // I_v / I_u, I_w / I_u
  double lv = 0.25, lw = 1.0 / 12.0;       // L_v / L_u, L_w / L_u
  // Davenport decay per component (rows u, v, w) and direction (cols y, z)
  std::array<std::array<double, 2>, 3> C = {{{10, 10}, {6.5, 6.5}, {3, 3}}};
  double length_scale_u(double z) const;  // EN 1991-1-4 (B.1)
};

struct WindPoint {
  double y = 0, z = 0;               // lateral / vertical position (the wind blows along x)
  double U = 0;                      // mean speed
  std::array<double, 3> sigma{};     // sigma_u, sigma_v, sigma_w
  std::array<double, 3> L{};         // integral scales
};
WindPoint wind_point(const TurbulenceModel& m, double y, double z);

double davenport_coherence(double n, double dy, double dz, double Uj, double Uk, double Cy, double Cz);

// The 3N x 3N cross-spectral matrix of (u, v, w) at the points, row/col
// 3 j + a; symmetric, the component blocks uncoupled. O(N^2).
Eigen::MatrixXd cross_spectral_matrix(const TurbulenceModel& m, const std::vector<WindPoint>& pts, double n);

// ---------------- aerodynamic admittance ----------------

std::complex<double> theodorsen_exact(double K);  // C(K), K = omega b / U (b the half-chord)
std::complex<double> sears_function(double K);
struct Admittance {
  enum class Kind { QuasiSteady, Sears, Liepmann, Table } kind = Kind::QuasiSteady;
  std::vector<std::pair<double, double>> table;  // (K, chi^2), K increasing
  double chi2(double K) const;
  // at frequency n [Hz] for a deck of width B in a mean speed U: K = pi n B / U
  double chi2_at(double n, double B, double U) const;
};

}  // namespace keta

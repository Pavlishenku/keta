// Multimodal flutter on the modes of a model.
//
// Each strip s of length L_s moves as x_s = (h_s, alpha_s) = T_s q (h
// DOWNWARD = -lift direction, alpha nose-up), q the modal coordinates (unit
// modal masses). Its self-excited forces in harmonic motion at omega are the
// impedance Z_s(omega) x_s with Z = -omega^2 Ma + i omega Ca + Ka added to the
// left-hand side, from
//   * Theodorsen's flat plate (half-chord b, elastic axis a half-chords aft of
//     mid-chord): exactly the matrices of keta's 2-DOF flutter_speed (Jones'
//     C(k)), so two modes of a sectional model reproduce it;
//   * Scanlan's flutter derivatives tabulated against K = B omega / U (L and
//     h positive downward, M and alpha nose-up):
//       L = q [K H1* h'/U + K H2* B a'/U + K^2 H3* a + K^2 H4* h/B] B
//       M = q [K A1* h'/U + K A2* B a'/U + K^2 A3* a + K^2 A4* h/B] B^2
//     i.e. Z = -q K^2 [[H4* + i H1*, B (H3* + i H2*)], [B (A4* + i A1*),
//     B^2 (A3* + i A2*)]], q = rho U^2 / 2.
// The modal impedance sum_s L_s T_s^T Z_s T_s enters the p-k iteration of
// each branch (aero frozen at the current omega, the eigenvalue of the
// branch followed until Im p = omega); the flutter speed is the first zero
// of the least branch damping on the sweep, bisected.
#pragma once
#include <complex>
#include <vector>

#include "keta/wind.h"

namespace keta {

struct ScanlanDerivatives {
  // rows (K, H1..H4, A1..A4), K increasing; linear in K, flat outside
  std::vector<std::array<double, 9>> rows;
  std::array<double, 8> at(double K) const;
};
// the flat-plate derivatives of Theodorsen's theory (the Jones C(k) of
// flutter_speed) identified from its impedance: B = 2 b, about the axis a
ScanlanDerivatives scanlan_from_theodorsen(double a, const std::vector<double>& Ks);

struct FlutterStrip {
  double L = 1;
  Eigen::Matrix<double, 2, Eigen::Dynamic> T;  // (h down, alpha) per mode
  enum class Kind { Theodorsen, Scanlan } kind = Kind::Theodorsen;
  double b = 1, a = 0;  // Theodorsen: half-chord, axis position
  double B = 2;         // Scanlan: deck width
  ScanlanDerivatives der;
};
struct MultimodalFlutter {
  std::vector<double> omega, zeta;  // modal
  std::vector<FlutterStrip> strips;
  double rho = 1.225;
};
FlutterResult multimodal_flutter(const MultimodalFlutter& f, double u_min, double u_max, int n_sweep = 60);
// the branches at one speed: (frequency, damping ratio) per mode
std::vector<std::pair<double, double>> multimodal_branches(const MultimodalFlutter& f, double U);

// strips of a model: drag ignored, h = -(lift-direction translation), alpha
// = rotation about z x d, per bar at its midpoint; T scaled by the modes
std::vector<FlutterStrip> flutter_strips(const Model& model, const ModalResults& modal, const Vec3& wind_dir,
                                         const FlutterStrip& prototype);

}  // namespace keta

#include "keta/traffic_transverse.h"

#include <cmath>

namespace keta {

std::vector<double> courbon_coefficients(const std::vector<double>& gy, const std::vector<double>& weights, double y) {
  const std::size_t n = gy.size();
  require(n >= 2, "Courbon: at least two girders");
  require(weights.empty() || weights.size() == n, "Courbon: one weight per girder (or none)");
  auto w = [&](std::size_t i) { return weights.empty() ? 1.0 : weights[i]; };
  double W = 0, yc = 0;
  for (std::size_t i = 0; i < n; ++i) {
    require(std::isfinite(gy[i]) && w(i) > 0, "Courbon: finite positions, positive stiffnesses");
    W += w(i);
    yc += w(i) * gy[i];
  }
  yc /= W;
  double J = 0;
  for (std::size_t i = 0; i < n; ++i) J += w(i) * (gy[i] - yc) * (gy[i] - yc);
  require(J > 0, "Courbon: the girders must not all sit at the same position");
  std::vector<double> K(n);
  for (std::size_t i = 0; i < n; ++i) K[i] = w(i) * (1 / W + (y - yc) * (gy[i] - yc) / J);
  return K;
}

InfluenceSurface separable_surface(const InfluenceLine& il, const std::vector<double>& ky,
                                   const std::vector<double>& k) {
  require(il.s.size() >= 2, "separable surface: the line needs 2 samples");
  require(ky.size() >= 2 && ky.size() == k.size(), "separable surface: a table of at least 2 (y, K) pairs");
  for (std::size_t j = 1; j < ky.size(); ++j) require(ky[j] > ky[j - 1], "separable surface: y strictly increasing");
  for (std::size_t i = 1; i < il.s.size(); ++i)
    require(il.s[i] > il.s[i - 1], "separable surface: line abscissae strictly increasing");
  InfluenceSurface s;
  s.nx = static_cast<int>(il.s.size());
  s.ny = static_cast<int>(ky.size());
  s.x = il.s;
  s.y = ky;
  s.eta.resize(static_cast<std::size_t>(s.nx) * static_cast<std::size_t>(s.ny));
  for (int j = 0; j < s.ny; ++j)
    for (int i = 0; i < s.nx; ++i)
      s.eta[static_cast<std::size_t>(j) * static_cast<std::size_t>(s.nx) + static_cast<std::size_t>(i)] =
          il.eta[static_cast<std::size_t>(i)] * k[static_cast<std::size_t>(j)];
  return s;
}

}  // namespace keta

#include "keta/hinge_law.h"

namespace keta {

void HingeLawFamily::validate() const {
  require(N.size() >= 1 && N.size() == laws.size(), "hinge law: one law per axial force");
  for (std::size_t k = 0; k < N.size(); ++k) {
    require(std::isfinite(N[k]), "hinge law: non-finite axial force");
    if (k > 0) require(N[k] > N[k - 1], "hinge law: axial forces must increase");
    require(laws[k].pts.size() >= 2 && laws[k].pts.size() == laws[0].pts.size(),
            "hinge law: every law needs the same number (>= 2) of points");
    for (std::size_t i = 1; i < laws[k].pts.size(); ++i)
      require(laws[k].pts[i].first > laws[k].pts[i - 1].first, "hinge law: rotations must increase");
  }
}

LinkLaw HingeLawFamily::at(double n) const {
  validate();
  require(n >= N.front() - 1e-12 * (1 + std::fabs(N.front())) && n <= N.back() + 1e-12 * (1 + std::fabs(N.back())),
          "hinge law: axial force outside the tabulated range");
  if (N.size() == 1) return laws[0];
  std::size_t k = 1;
  while (k + 1 < N.size() && n > N[k]) ++k;
  const double t = std::min(1.0, std::max(0.0, (n - N[k - 1]) / (N[k] - N[k - 1])));
  LinkLaw out;
  for (std::size_t i = 0; i < laws[k].pts.size(); ++i) {
    const auto& a = laws[k - 1].pts[i];
    const auto& b = laws[k].pts[i];
    out.pts.emplace_back(a.first + t * (b.first - a.first), a.second + t * (b.second - a.second));
  }
  for (std::size_t i = 1; i < out.pts.size(); ++i)
    require(out.pts[i].first > out.pts[i - 1].first, "hinge law: the interpolated rotations do not increase");
  return out;
}

}  // namespace keta

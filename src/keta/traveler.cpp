#include "keta/traveler.h"

#include <cmath>

namespace keta {

std::vector<double> traveler_reactions(const FormTraveler& t) {
  const std::size_t n = t.supports.size();
  require(n >= 3, "traveler: at least 3 support points");
  require(t.gravity >= 0 && t.gravity <= 2, "traveler: gravity axis 0, 1 or 2");
  require(std::isfinite(t.weight) && t.weight > 0 && std::isfinite(t.cg_x) && std::isfinite(t.cg_y),
          "traveler: finite positive weight and finite centre of gravity");
  double mx = 0, my = 0;
  for (const TravelerSupport& s : t.supports) {
    require(std::isfinite(s.x) && std::isfinite(s.y), "traveler: finite support coordinates");
    mx += s.x;
    my += s.y;
  }
  mx /= double(n);
  my /= double(n);
  double sxx = 0, syy = 0, sxy = 0;
  for (const TravelerSupport& s : t.supports) {
    const double x = s.x - mx, y = s.y - my;
    sxx += x * x;
    syy += y * y;
    sxy += x * y;
  }
  const double det = sxx * syy - sxy * sxy;
  require(det > 1e-12 * std::max(sxx * syy, 1e-300), "traveler: support points collinear or coincident");
  // F_i = W / n + b (x_i - mx) + c (y_i - my): sum F = W; the moment balances
  // [Sxx Sxy; Sxy Syy] [b c] = W (cg - mean)
  const double rx = t.weight * (t.cg_x - mx), ry = t.weight * (t.cg_y - my);
  const double b = (rx * syy - ry * sxy) / det, c = (ry * sxx - rx * sxy) / det;
  std::vector<double> F(n);
  for (std::size_t i = 0; i < n; ++i) F[i] = t.weight / double(n) + b * (t.supports[i].x - mx) + c * (t.supports[i].y - my);
  return F;
}

TravelerStability traveler_stability(const FormTraveler& t) {
  TravelerStability s;
  s.reactions = traveler_reactions(t);
  for (std::size_t i = 0; i < t.supports.size(); ++i)
    if (!t.supports[i].anchor && s.reactions[i] < 0) s.uplifting.push_back(t.supports[i].node);
  return s;
}

std::vector<Load> traveler_loads(const FormTraveler& t) {
  const std::vector<double> F = traveler_reactions(t);
  std::vector<Load> out;
  out.reserve(F.size());
  for (std::size_t i = 0; i < F.size(); ++i) {
    Vec3 f = Vec3::Zero();
    f(t.gravity) = -F[i];
    out.push_back(Load::nodal(t.supports[i].node, f));
  }
  return out;
}

std::vector<Load> traveler_move_loads(const FormTraveler& from, const FormTraveler& to) {
  std::vector<Load> out = traveler_loads(from);
  for (Load& l : out) l.F = -l.F;
  const std::vector<Load> add = traveler_loads(to);
  out.insert(out.end(), add.begin(), add.end());
  return out;
}

}  // namespace keta

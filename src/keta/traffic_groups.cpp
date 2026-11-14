#include "keta/traffic_groups.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
double effect(const InfluenceSurface& s, const Wheels& w, double x, double y, double sgnx) {
  double v = 0;
  for (const auto& q : w) v += q[2] * s.at(x + sgnx * q[0], y + q[1]);
  return v;
}
std::vector<double> uniq(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), v.end());
  return v;
}
}  // namespace

std::vector<VehiclePlacement> vehicle_surface_placements(const InfluenceSurface& s, const Wheels& w, const VehicleBand& b) {
  require(s.nx >= 2 && s.ny >= 2, "vehicle: an empty surface");
  require(!w.empty(), "vehicle: no wheel");
  double dymin = 1e300, dymax = -1e300;
  for (const auto& q : w) {
    require(std::isfinite(q[0]) && std::isfinite(q[1]) && std::isfinite(q[2]), "vehicle: non-finite wheel");
    dymin = std::min(dymin, q[1]);
    dymax = std::max(dymax, q[1]);
  }
  const double ylo = std::max(b.y0 - dymin, s.y.front() - dymax), yhi = std::min(b.y1 - dymax, s.y.back() - dymin);
  require(yhi >= ylo - 1e-12, "vehicle: the band is narrower than the vehicle");
  const double sg = b.minimise ? -1.0 : 1.0;
  std::vector<VehiclePlacement> out;
  for (int dir = 0; dir < (b.both_directions ? 2 : 1); ++dir) {
    const double sx = dir == 0 ? 1.0 : -1.0;
    std::vector<double> X, Y{ylo, yhi};
    for (const auto& q : w) {
      for (double xg : s.x) X.push_back(xg - sx * q[0]);
      for (double yg : s.y) {
        const double y = yg - q[1];
        if (y > ylo && y < yhi) Y.push_back(y);
      }
    }
    X = uniq(X);
    Y = uniq(Y);
    for (double x : X)
      for (double y : Y) out.push_back({sg * effect(s, w, x, y, sx), x, y, dir == 1});
  }
  std::sort(out.begin(), out.end(), [](const VehiclePlacement& a, const VehiclePlacement& c) { return a.value > c.value; });
  for (auto& p : out) p.value *= sg;
  return out;
}

VehiclePlacement vehicle_surface_extreme(const InfluenceSurface& s, const Wheels& w, const VehicleBand& b) {
  const std::vector<VehiclePlacement> p = vehicle_surface_placements(s, w, b);
  return p.front();
}

Wheels lm2_wheels(double beta_q) {
  require(beta_q > 0, "LM2: beta_Q > 0");
  return {{0, -1.0, 0.5 * beta_q * 400e3}, {0, 1.0, 0.5 * beta_q * 400e3}};
}

std::vector<std::pair<double, double>> annex_a_axle_lines(const std::string& name) {
  // EN 1991-2:2003 Annex A, tables A.1 / A.2 (e = 1.50 m within a group,
  // 12 m between the two bogies of the double vehicles)
  auto train = [](std::vector<double> loads) {
    std::vector<std::pair<double, double>> l;
    for (std::size_t i = 0; i < loads.size(); ++i) l.emplace_back(1.5 * static_cast<double>(i), loads[i] * 1e3);
    return l;
  };
  auto bogies = [](std::size_t n1, std::size_t n2, double p) {
    std::vector<std::pair<double, double>> l;
    for (std::size_t i = 0; i < n1; ++i) l.emplace_back(1.5 * static_cast<double>(i), p * 1e3);
    const double s2 = 1.5 * static_cast<double>(n1 - 1) + 12.0;
    for (std::size_t i = 0; i < n2; ++i) l.emplace_back(s2 + 1.5 * static_cast<double>(i), p * 1e3);
    return l;
  };
  if (name == "600/150") return train(std::vector<double>(4, 150));
  if (name == "900/150") return train(std::vector<double>(6, 150));
  if (name == "1200/150") return train(std::vector<double>(8, 150));
  if (name == "1200/200") return train(std::vector<double>(6, 200));
  if (name == "1500/150") return train(std::vector<double>(10, 150));
  if (name == "1500/200") {
    std::vector<double> v(7, 200);
    v.push_back(100);
    return train(v);
  }
  if (name == "1800/150") return train(std::vector<double>(12, 150));
  if (name == "1800/200") return train(std::vector<double>(9, 200));
  if (name == "2400/200") return train(std::vector<double>(12, 200));
  if (name == "2400/240") return train(std::vector<double>(10, 240));
  if (name == "2400/200/200") return bogies(6, 6, 200);
  if (name == "3000/200") return train(std::vector<double>(15, 200));
  if (name == "3000/240") {
    std::vector<double> v(12, 240);
    v.push_back(120);
    return train(v);
  }
  if (name == "3000/200/200") return bogies(8, 7, 200);
  if (name == "3600/200") return train(std::vector<double>(18, 200));
  if (name == "3600/240") return train(std::vector<double>(15, 240));
  if (name == "3600/240/240") return bogies(8, 7, 240);
  throw Error("Annex A: unknown special vehicle " + name);
}

Wheels axle_lines_to_wheels(const std::vector<std::pair<double, double>>& lines,
                            const std::vector<std::pair<double, double>>& layout) {
  require(!layout.empty(), "axle lines: a transverse layout is required");
  double sum = 0;
  for (const auto& [dy, f] : layout) sum += f;
  require(std::fabs(sum - 1) < 1e-9, "axle lines: the layout fractions must sum to 1");
  Wheels w;
  for (const auto& [dx, P] : lines)
    for (const auto& [dy, f] : layout) w.push_back({dx, dy, f * P});
  return w;
}

GroupResult group_1a(const InfluenceSurface& s, const Lm1SurfaceOptions& lm1, const std::vector<Footway>& footways,
                     double q_footway) {
  GroupResult g;
  g.lm1 = lm1_surface_extreme(s, lm1);
  const int sign = lm1.minimise ? -1 : 1;
  const double xa = std::max(lm1.x0, s.x.front()), xb = std::min(lm1.x1, s.x.back());
  for (const Footway& f : footways) {
    require(f.y1 > f.y0, "group 1a: a footway needs y0 < y1");
    g.footways += q_footway * surface_integral_where(s, s, sign, xa, xb, f.y0, f.y1);
  }
  g.value = g.lm1.value + g.footways;
  return g;
}

GroupResult group_4(const InfluenceSurface& s, double y0, double y1, double q, bool minimise) {
  require(y1 > y0, "group 4: y0 < y1");
  GroupResult g;
  g.value = q * surface_integral_where(s, s, minimise ? -1 : 1, s.x.front(), s.x.back(), y0, y1);
  return g;
}

GroupResult group_5(const InfluenceSurface& s, const Wheels& vehicle, const VehicleBand& band, const Lm1SurfaceOptions& lm1,
                    double clear) {
  require(clear >= 0, "group 5: the clearance must be >= 0");
  Lm1SurfaceOptions f = lm1;  // frequent values (EN 1990 A2, psi1)
  for (double& a : f.factors.aQ) a *= 0.75;
  f.factors.aq1 *= 0.40;
  f.factors.aq_other *= 0.40;
  f.minimise = band.minimise;
  VehicleBand vb = band;
  const std::vector<VehiclePlacement> P = vehicle_surface_placements(s, vehicle, vb);
  const double sg = band.minimise ? -1.0 : 1.0;
  const Lm1SurfaceResult free_lm1 = lm1_surface_extreme(s, f);  // an upper bound of every restricted LM1
  double dxmin = 1e300, dxmax = -1e300, dymin = 1e300, dymax = -1e300;
  for (const auto& q : vehicle) {
    dxmin = std::min(dxmin, q[0]);
    dxmax = std::max(dxmax, q[0]);
    dymin = std::min(dymin, q[1]);
    dymax = std::max(dymax, q[1]);
  }
  GroupResult best;
  best.value = -1e300;
  bool found = false;
  for (const VehiclePlacement& p : P) {
    if (found && sg * (p.value + free_lm1.value) <= sg * best.value) break;  // nothing better below
    const double sx = p.reversed ? -1.0 : 1.0;
    const double xf = p.x + std::min(sx * dxmin, sx * dxmax), xr = p.x + std::max(sx * dxmin, sx * dxmax);
    Lm1SurfaceOptions o = f;
    o.exclude.push_back({xf - clear, xr + clear, p.y + dymin, p.y + dymax});
    const Lm1SurfaceResult l = lm1_surface_extreme(s, o);
    ++best.vehicle_positions_tried;
    const double tot = p.value + l.value;
    if (!found || sg * tot > sg * best.value) {
      const int tried = best.vehicle_positions_tried;
      best.value = tot;
      best.lm1 = l;
      best.vehicle = p;
      best.vehicle_positions_tried = tried;
      found = true;
    }
  }
  return best;
}

double vehicle_value(const InfluenceSurface& target, const Wheels& w, const VehiclePlacement& p) {
  return effect(target, w, p.x, p.y, p.reversed ? -1.0 : 1.0);
}

namespace {
Lm1SurfaceOptions frequent(Lm1SurfaceOptions f, bool minimise) {
  for (double& a : f.factors.aQ) a *= 0.75;
  f.factors.aq1 *= 0.40;
  f.factors.aq_other *= 0.40;
  f.minimise = minimise;
  return f;
}
std::array<double, 4> vehicle_box(const Wheels& w, const VehiclePlacement& p, double clear) {
  double dxmin = 1e300, dxmax = -1e300, dymin = 1e300, dymax = -1e300;
  for (const auto& q : w) {
    dxmin = std::min(dxmin, q[0]);
    dxmax = std::max(dxmax, q[0]);
    dymin = std::min(dymin, q[1]);
    dymax = std::max(dymax, q[1]);
  }
  const double sx = p.reversed ? -1.0 : 1.0;
  const double xf = p.x + std::min(sx * dxmin, sx * dxmax), xr = p.x + std::max(sx * dxmin, sx * dxmax);
  return {xf - clear, xr + clear, p.y + dymin, p.y + dymax};
}
}  // namespace

GroupResult group_extreme(const InfluenceSurface& driver, const GroupSpec& g) {
  switch (g.group) {
    case LoadGroup::Gr1a: return group_1a(driver, g.lm1, g.footways, g.q_footway);
    case LoadGroup::Gr4: return group_4(driver, g.y0, g.y1, g.q_crowd, g.lm1.minimise);
    case LoadGroup::Gr5: return group_5(driver, g.vehicle, g.band, g.lm1, g.clear);
  }
  return {};
}

double group_concomitant(const InfluenceSurface& target, const InfluenceSurface& driver, const GroupResult& r,
                         const GroupSpec& g) {
  switch (g.group) {
    case LoadGroup::Gr1a: {
      double v = lm1_surface_value(target, driver, r.lm1, g.lm1);
      const int sign = g.lm1.minimise ? -1 : 1;
      const double xa = std::max(g.lm1.x0, driver.x.front()), xb = std::min(g.lm1.x1, driver.x.back());
      for (const Footway& f : g.footways) v += g.q_footway * surface_integral_where(target, driver, sign, xa, xb, f.y0, f.y1);
      return v;
    }
    case LoadGroup::Gr4:
      return g.q_crowd * surface_integral_where(target, driver, g.lm1.minimise ? -1 : 1, driver.x.front(), driver.x.back(), g.y0, g.y1);
    case LoadGroup::Gr5: {
      Lm1SurfaceOptions f = frequent(g.lm1, g.band.minimise);
      f.exclude.push_back(vehicle_box(g.vehicle, r.vehicle, g.clear));
      return vehicle_value(target, g.vehicle, r.vehicle) + lm1_surface_value(target, driver, r.lm1, f);
    }
  }
  return 0;
}

}  // namespace keta

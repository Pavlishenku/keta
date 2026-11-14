#include "keta/rail_actions.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
void check_len(double L) { require(std::isfinite(L) && L >= 0, "rail action: loaded length >= 0"); }
void check_alpha(double a) { require(std::isfinite(a) && a > 0, "rail action: alpha > 0"); }
}  // namespace

double rail_phi2(double L) {
  require(std::isfinite(L) && L > 0.04, "Phi2: L_phi > 0.04 m");
  return std::clamp(1.44 / (std::sqrt(L) - 0.2) + 0.82, 1.0, 1.67);
}
double rail_phi3(double L) {
  require(std::isfinite(L) && L > 0.04, "Phi3: L_phi > 0.04 m");
  return std::clamp(2.16 / (std::sqrt(L) - 0.2) + 0.73, 1.0, 2.0);
}

double rail_traction(double L, double alpha) {
  check_len(L);
  check_alpha(alpha);
  return std::min(33e3 * L, 1000e3) * alpha;
}

double rail_braking(double L, RailModel model, double alpha) {
  check_len(L);
  check_alpha(alpha);
  if (model == RailModel::SW2) return 35e3 * L;  // no cap, no alpha (SW/2 is not classified)
  return std::min(20e3 * L, 6000e3) * alpha;
}

double rail_nosing(double alpha) {
  check_alpha(alpha);
  return 100e3 * std::max(alpha, 1.0);
}

double rail_centrifugal_factor(double v, double Lf) {
  require(std::isfinite(v) && v >= 0, "centrifugal: V >= 0 [km/h]");
  require(std::isfinite(Lf) && Lf > 0, "centrifugal: L_f > 0");
  if (v <= 120 || Lf <= 2.88) return 1.0;
  v = std::min(v, 300.0);
  Lf = std::min(Lf, 150.0);  // EN 1991-2 6.5.1(8): L_f > 150 m is taken as 150 m
  return 1 - (v - 120) / 1000 * (814 / v + 1.75) * (1 - std::sqrt(2.88 / Lf));
}

double rail_centrifugal_force(double v, double r, double q, double Lf) {
  require(r > 0, "centrifugal: radius > 0 (infinite for a straight track)");
  require(std::isfinite(q), "centrifugal: finite vertical load");
  if (std::isinf(r)) return 0;
  return v * v / (127 * r) * rail_centrifugal_factor(v, Lf) * q;
}

EccentricRail eccentric_rail(double total, double ratio, double gauge) {
  require(std::isfinite(total), "eccentric rail load: finite total");
  require(ratio >= 1 && ratio <= 1.25 + 1e-12, "eccentric rail load: 1 <= q2/q1 <= 1.25 (EN 1991-2 §6.3.5)");
  require(gauge > 0, "eccentric rail load: gauge > 0");
  EccentricRail e;
  e.q1 = total / (1 + ratio);
  e.q2 = total - e.q1;
  e.gauge = gauge;
  return e;
}

std::vector<LoadPatch> ballast_patches(double Q, const BallastOptions& o) {
  require(std::isfinite(Q), "ballast: finite load");
  require(o.depth >= 0 && o.sleeper_width > 0 && o.sleeper_length > 0 && o.block_width > 0 && o.gauge > 0 &&
              o.spacing > 0,
          "ballast: positive dimensions, depth >= 0");
  const double ax = o.sleeper_width + o.depth / 2;  // 4:1 on each side: + 2 h / 4
  std::vector<std::pair<double, double>> along;     // (dx, share)
  if (o.three_sleepers)
    along = {{-o.spacing, 0.25}, {0.0, 0.5}, {o.spacing, 0.25}};
  else
    along = {{0.0, 1.0}};
  std::vector<LoadPatch> out;
  for (const auto& [dx, f] : along) {
    if (o.sleeper == Sleeper::Monobloc) {
      out.push_back({dx, 0.0, f * Q, ax, o.sleeper_length + o.depth / 2});
    } else {
      const double ay = o.block_width + o.depth / 2;
      out.push_back({dx, -o.gauge / 2, 0.5 * f * Q, ax, ay});
      out.push_back({dx, o.gauge / 2, 0.5 * f * Q, ax, ay});
    }
  }
  return out;
}

double patches_value(const InfluenceSurface& s, const std::vector<LoadPatch>& ps, double rx, double ry) {
  double v = 0;
  for (const LoadPatch& p : ps) {
    const double cx = rx + p.dx, cy = ry + p.dy;
    v += p.load / (p.ax * p.ay) *
         s.integral_box(cx - p.ax / 2, cx + p.ax / 2, cy - p.ay / 2, cy + p.ay / 2);
  }
  return v;
}

}  // namespace keta

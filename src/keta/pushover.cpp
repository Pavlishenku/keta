// Performance-point procedures — see pushover.h.
#include "keta/pushover.h"

#include <algorithm>
#include <cmath>

#include "keta/core.h"

namespace keta {

N2Result n2_performance_point(double Sdy, double Sy, const std::function<double(double)>& Se, double TC) {
  require(Sdy > 0 && Sy > 0, "n2: yield point must be positive");
  N2Result r;
  r.T = 2 * kPi * std::sqrt(Sdy / Sy);
  r.Se = Se(r.T);
  r.Sde = r.Se * r.T * r.T / (4 * kPi * kPi);
  r.Ry = r.Se / Sy;
  if (r.Ry <= 1) {
    r.mu = r.Ry;  // elastic: displacement ratio to yield
    r.Sdp = r.Sde;
    r.Sp = r.Se;
    return r;
  }
  r.mu = r.T >= TC ? r.Ry : (r.Ry - 1) * TC / r.T + 1;
  r.Sdp = r.mu * r.Sde / r.Ry;
  r.Sp = r.Se / r.Ry;
  return r;
}

namespace {
struct Curve {
  std::vector<CapacityPoint> p;  // with the origin
  double sa(double d) const {
    for (std::size_t i = 1; i < p.size(); ++i)
      if (d <= p[i].Sd) return p[i - 1].Sa + (p[i].Sa - p[i - 1].Sa) * (d - p[i - 1].Sd) / (p[i].Sd - p[i - 1].Sd);
    return p.back().Sa;
  }
  double area(double d) const {  // integral of Sa dSd from 0 to d
    double a = 0;
    for (std::size_t i = 1; i < p.size(); ++i) {
      const double x1 = std::min(d, p[i].Sd);
      if (x1 <= p[i - 1].Sd) break;
      a += 0.5 * (p[i - 1].Sa + sa(x1)) * (x1 - p[i - 1].Sd);
    }
    return a;
  }
};
}  // namespace

Atc40Result atc40_performance_point(const std::vector<CapacityPoint>& capacity, double CA, double CV, double g,
                                    Atc40Behaviour type) {
  require(capacity.size() >= 2, "atc40: capacity curve needs two points or more");
  Curve c;
  c.p.push_back({0, 0});
  for (const CapacityPoint& q : capacity) {
    require(q.Sd > c.p.back().Sd, "atc40: capacity Sd must increase");
    c.p.push_back(q);
  }
  const double k0 = c.p[1].Sa / c.p[1].Sd;  // initial slope
  // state at a trial displacement dpi
  auto state = [&](double dpi, Atc40Result& r) {
    const double api = c.sa(dpi);
    // equal area: bilinear O-(dy, k0 dy)-(dpi, api); area = 0.5 k0 dy^2 +
    // 0.5 (k0 dy + api)(dpi - dy) = A  ->  0.5 dy (api - k0 dpi) ... linear in dy
    const double A = c.area(dpi);
    double dy = (2 * A - api * dpi) / (k0 * dpi - api);
    if (!(k0 * dpi > api)) dy = dpi;  // still elastic
    dy = std::min(std::max(dy, 0.0), dpi);
    const double ay = k0 * dy;
    const double x = (ay * dpi - dy * api) / (api * dpi);
    const double beta0 = 63.7 * x;
    double kappa = 0.33, sra_min = 0.56, srv_min = 0.67;
    if (type == Atc40Behaviour::A) {
      kappa = beta0 <= 16.25 ? 1.0 : 1.13 - 0.51 * x;
      sra_min = 0.33;
      srv_min = 0.50;
    } else if (type == Atc40Behaviour::B) {
      kappa = beta0 <= 25.0 ? 0.67 : 0.845 - 0.446 * x;
      sra_min = 0.44;
      srv_min = 0.56;
    }
    r.beta_eff = kappa * beta0 + 5.0;
    r.SRA = std::max(sra_min, (3.21 - 0.68 * std::log(r.beta_eff)) / 2.12);
    r.SRV = std::max(srv_min, (2.31 - 0.41 * std::log(r.beta_eff)) / 1.65);
    r.Sdy = dy;
    r.Say = ay;
    r.Sdp = dpi;
    r.Sap = api;
    // reduced demand acceleration on the radial line through (dpi, api):
    // T = 2 pi sqrt(d / a); Sa_demand(T) = min(2.5 CA SRA, CV SRV / T) g
    const double T = 2 * kPi * std::sqrt(dpi / api);
    const double demand = std::min(2.5 * CA * r.SRA, CV * r.SRV / T) * g;
    return api - demand;
  };
  Atc40Result r;
  double lo = 1e-9 * c.p[1].Sd, hi = c.p.back().Sd;
  double flo = state(lo, r), fhi = state(hi, r);
  require(flo < 0 && fhi >= 0, "atc40: no performance point on the capacity curve (demand exceeds capacity)");
  for (int i = 0; i < 300 && hi - lo > 1e-15 * hi; ++i) {
    const double m = 0.5 * (lo + hi);
    (state(m, r) < 0 ? lo : hi) = m;
  }
  state(0.5 * (lo + hi), r);
  return r;
}

}  // namespace keta

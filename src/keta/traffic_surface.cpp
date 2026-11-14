#include "keta/traffic_surface.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
struct Pt {
  double x, y, v, d;  // position, value, driver
};
using Poly = std::vector<Pt>;

// keep the part of the polygon where f(p) >= 0 (f linear along edges)
template <class F>
Poly clip(const Poly& in, F f) {
  Poly out;
  const std::size_t n = in.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Pt& a = in[i];
    const Pt& b = in[(i + 1) % n];
    const double fa = f(a), fb = f(b);
    if (fa >= 0) out.push_back(a);
    if ((fa >= 0) != (fb >= 0)) {
      const double t = fa / (fa - fb);
      out.push_back({a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), a.v + t * (b.v - a.v), a.d + t * (b.d - a.d)});
    }
  }
  return out;
}

double integrate(const Poly& p) {  // value linear over the polygon: fan of triangles
  double s = 0;
  for (std::size_t i = 1; i + 1 < p.size(); ++i) {
    const Pt &a = p[0], &b = p[i], &c = p[i + 1];
    const double area = 0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
    s += area * (a.v + b.v + c.v) / 3;
  }
  return s;
}

double rank_Q(int r, const Lm1Factors& f) {
  static const double Qk[3] = {300e3, 200e3, 100e3};
  return r <= 3 ? f.aQ[r - 1] * Qk[r - 1] : 0.0;
}
double rank_q(int r, const Lm1Factors& f) { return r == 1 ? f.aq1 * 9.0e3 : f.aq_other * 2.5e3; }

// the tandem's effect (unit axle load, two wheels of 1/2) at first axle x,
// lane centre yc
double tandem(const InfluenceSurface& s, double x, double yc) {
  return 0.5 * (s.at(x, yc - 1) + s.at(x, yc + 1) + s.at(x + 1.2, yc - 1) + s.at(x + 1.2, yc + 1));
}
}  // namespace

double surface_integral_where(const InfluenceSurface& value, const InfluenceSurface& driver, int sign, double x0,
                              double x1, double y0, double y1) {
  require(value.nx == driver.nx && value.ny == driver.ny && value.x == driver.x && value.y == driver.y,
          "surface integral: the two surfaces must share their grid");
  if (x1 <= x0 || y1 <= y0) return 0;
  double total = 0;
  for (int j = 0; j + 1 < value.ny; ++j) {
    const double ya = value.y[static_cast<std::size_t>(j)], yb = value.y[static_cast<std::size_t>(j) + 1];
    if (yb <= y0 || ya >= y1) continue;
    for (int i = 0; i + 1 < value.nx; ++i) {
      const double xa = value.x[static_cast<std::size_t>(i)], xb = value.x[static_cast<std::size_t>(i) + 1];
      if (xb <= x0 || xa >= x1) continue;
      const Pt a{xa, ya, value.node_eta(i, j), driver.node_eta(i, j)};
      const Pt b{xb, ya, value.node_eta(i + 1, j), driver.node_eta(i + 1, j)};
      const Pt c{xa, yb, value.node_eta(i, j + 1), driver.node_eta(i, j + 1)};
      const Pt d{xb, yb, value.node_eta(i + 1, j + 1), driver.node_eta(i + 1, j + 1)};
      const Poly tris[4] = {{a, b, d}, {a, d, c}, {a, b, c}, {b, d, c}};
      double cell = 0;
      for (const Poly& t : tris) {
        Poly p = clip(t, [&](const Pt& q) { return q.x - x0; });
        p = clip(p, [&](const Pt& q) { return x1 - q.x; });
        p = clip(p, [&](const Pt& q) { return q.y - y0; });
        p = clip(p, [&](const Pt& q) { return y1 - q.y; });
        if (sign != 0) p = clip(p, [&](const Pt& q) { return sign * q.d; });
        if (p.size() >= 3) cell += integrate(p);
      }
      total += 0.5 * cell;  // the two diagonal splittings averaged
    }
  }
  return total;
}

Lm1SurfaceResult lm1_surface_extreme(const InfluenceSurface& s0, const Lm1SurfaceOptions& o) {
  require(s0.nx >= 2 && s0.ny >= 2, "LM1 surface: an empty surface");
  require(o.y1 > o.y0 && o.lane_width > 0 && o.step_y > 0, "LM1 surface: carriageway y0 < y1, lane width and step > 0");
  InfluenceSurface s = s0;
  const double sg = o.minimise ? -1.0 : 1.0;
  for (double& e : s.eta) e *= sg;
  const double xa = std::max(o.x0, s.x.front()), xb = std::min(o.x1, s.x.back());
  require(xb - xa > 1.2, "LM1 surface: the loaded length must exceed the tandem");
  const double W = o.y1 - o.y0, w = o.lane_width;
  const int n = o.n_lanes > 0 ? o.n_lanes : en_lane_division(W).n_lanes;
  require(n >= 1 && n * w <= W + 1e-9, "LM1 surface: the lanes do not fit the carriageway");
  const double w_lane = o.n_lanes > 0 ? w : (en_lane_division(W).lane_width > 0 ? en_lane_division(W).lane_width : w);
  // lane edges searched: the step, and every breakpoint of the effects
  std::vector<double> Y;
  for (double y = o.y0; y <= o.y1 - w_lane + 1e-12; y += o.step_y) Y.push_back(y);
  Y.push_back(o.y1 - w_lane);
  for (double yg : s.y)
    for (double c : {yg, yg - w_lane, yg - w_lane / 2 + 1, yg - w_lane / 2 - 1}) Y.push_back(c);
  std::sort(Y.begin(), Y.end());
  Y.erase(std::remove_if(Y.begin(), Y.end(), [&](double y) { return y < o.y0 - 1e-12 || y > o.y1 - w_lane + 1e-12; }),
          Y.end());
  Y.erase(std::unique(Y.begin(), Y.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), Y.end());
  // per lane position: the tandem's best place and value, the strip's favourable integral
  std::vector<double> X;
  for (double xg : s.x)
    for (double c : {xg, xg - 1.2}) X.push_back(c);
  X.push_back(xa);
  X.push_back(xb - 1.2);
  for (const auto& bx : o.exclude) {
    require(bx[1] > bx[0] && bx[3] > bx[2], "LM1 surface: an exclusion box needs x0 < x1, y0 < y1");
    for (double c : {bx[0] - 1.2, bx[1]}) X.push_back(c);  // the admissible intervals' ends
    for (double c : {bx[2] - w_lane, bx[3], bx[2] - w_lane / 2 - 1, bx[3] - w_lane / 2 + 1}) Y.push_back(c);
  }
  std::sort(Y.begin(), Y.end());
  Y.erase(std::remove_if(Y.begin(), Y.end(), [&](double y) { return y < o.y0 - 1e-12 || y > o.y1 - w_lane + 1e-12; }),
          Y.end());
  Y.erase(std::unique(Y.begin(), Y.end(), [](double a, double b) { return std::fabs(a - b) < 1e-9; }), Y.end());
  X.erase(std::remove_if(X.begin(), X.end(), [&](double x) { return x < xa - 1e-12 || x > xb - 1.2 + 1e-12; }), X.end());
  // a tandem (axles x, x + 1.2, wheels yc -+ 1) touching an exclusion box
  auto blocked = [&](double x, double yc) {
    for (const auto& bx : o.exclude)
      if (x < bx[1] && x + 1.2 > bx[0] && yc - 1 < bx[3] && yc + 1 > bx[2]) return true;
    return false;
  };
  auto excluded_integral = [&](double ya, double yb) {
    double e = 0;
    for (const auto& bx : o.exclude)
      e += surface_integral_where(s, s, 1, std::max(xa, bx[0]), std::min(xb, bx[1]), std::max(ya, bx[2]), std::min(yb, bx[3]));
    return e;
  };
  const std::size_t m = Y.size();
  std::vector<double> T(m, -1e300), TX(m, xa), U(m);
  for (std::size_t k = 0; k < m; ++k) {
    const double yc = Y[k] + w_lane / 2;
    for (double x : X) {
      if (blocked(x, yc)) continue;
      const double t = tandem(s, x, yc);
      if (t > T[k]) {
        T[k] = t;
        TX[k] = x;
      }
    }
    U[k] = surface_integral_where(s, s, 1, xa, xb, Y[k], Y[k] + w_lane) - excluded_integral(Y[k], Y[k] + w_lane);
  }
  const double Utot = surface_integral_where(s, s, 1, xa, xb, o.y0, o.y1) - excluded_integral(o.y0, o.y1);
  const double qrem = o.factors.aq_other * 2.5e3;
  const int nr = std::min(n, 3);
  auto gain = [&](int r, std::size_t k) {
    return rank_Q(r, o.factors) * std::max(T[k], 0.0) + (rank_q(r, o.factors) - qrem) * U[k];
  };
  // the other n - nr lanes must still fit in the gaps
  auto fits = [&](const std::vector<double>& starts) {
    std::vector<double> e = starts;
    std::sort(e.begin(), e.end());
    int room = 0;
    double from = o.y0;
    for (double y : e) {
      room += static_cast<int>(std::floor((y - from) / w_lane + 1e-9));
      from = y + w_lane;
    }
    room += static_cast<int>(std::floor((o.y1 - from) / w_lane + 1e-9));
    return room >= n - nr;
  };
  double best = -1e300;
  std::vector<std::size_t> pick;
  const double gap = w_lane - 1e-9;
  for (std::size_t a = 0; a < m; ++a) {
    const double ga = gain(1, a);
    if (nr == 1) {
      if (ga > best && fits({Y[a]})) {
        best = ga;
        pick = {a};
      }
      continue;
    }
    for (std::size_t b = 0; b < m; ++b) {
      if (std::fabs(Y[a] - Y[b]) < gap) continue;
      const double gb = ga + gain(2, b);
      if (nr == 2) {
        if (gb > best && fits({Y[a], Y[b]})) {
          best = gb;
          pick = {a, b};
        }
        continue;
      }
      for (std::size_t c = 0; c < m; ++c) {
        if (std::fabs(Y[c] - Y[a]) < gap || std::fabs(Y[c] - Y[b]) < gap) continue;
        const double gc = gb + gain(3, c);
        if (gc > best && fits({Y[a], Y[b], Y[c]})) {
          best = gc;
          pick = {a, b, c};
        }
      }
    }
  }
  require(!pick.empty(), "LM1 surface: no admissible lane arrangement");
  Lm1SurfaceResult r;
  double usum = 0;
  for (int i = 0; i < nr; ++i) {
    const std::size_t k = pick[static_cast<std::size_t>(i)];
    Lm1Lane l;
    l.rank = i + 1;
    l.y_start = Y[k];
    l.tandem = T[k] > 0;
    l.ts_x = TX[k];
    l.ts = sg * rank_Q(l.rank, o.factors) * std::max(T[k], 0.0);
    l.udl = sg * rank_q(l.rank, o.factors) * U[k];
    usum += U[k];
    r.lanes.push_back(l);
  }
  r.remaining = sg * qrem * (Utot - usum);
  r.value = sg * (qrem * Utot + best);
  return r;
}

double lm1_surface_value(const InfluenceSurface& target, const InfluenceSurface& driver, const Lm1SurfaceResult& r,
                         const Lm1SurfaceOptions& o) {
  const int sg = o.minimise ? -1 : 1;
  const double xa = std::max(o.x0, driver.x.front()), xb = std::min(o.x1, driver.x.back());
  const double W = o.y1 - o.y0;
  const double w_lane = o.n_lanes > 0 ? o.lane_width
                                      : (en_lane_division(W).lane_width > 0 ? en_lane_division(W).lane_width : o.lane_width);
  const double qrem = o.factors.aq_other * 2.5e3;
  // the UDL zones: where the driver is unfavourable, outside the exclusion boxes
  auto zone = [&](double ya, double yb) {
    double z = surface_integral_where(target, driver, sg, xa, xb, ya, yb);
    for (const auto& bx : o.exclude)
      z -= surface_integral_where(target, driver, sg, std::max(xa, bx[0]), std::min(xb, bx[1]), std::max(ya, bx[2]),
                                  std::min(yb, bx[3]));
    return z;
  };
  double v = qrem * zone(o.y0, o.y1);
  for (const Lm1Lane& l : r.lanes) {
    if (l.tandem) v += rank_Q(l.rank, o.factors) * tandem(target, l.ts_x, l.y_start + w_lane / 2);
    v += (rank_q(l.rank, o.factors) - qrem) * zone(l.y_start, l.y_start + w_lane);
  }
  return v;
}

}  // namespace keta

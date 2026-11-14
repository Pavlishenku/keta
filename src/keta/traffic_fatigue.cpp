#include "keta/traffic_fatigue.h"

#include <algorithm>
#include <cmath>

namespace keta {

std::vector<double> turning_points(const std::vector<double>& series) {
  std::vector<double> s;
  s.reserve(series.size());
  for (double v : series) {
    require(std::isfinite(v), "rainflow: non-finite value in the series");
    if (s.empty() || s.back() != v) s.push_back(v);
  }
  if (s.size() <= 2) return s;
  std::vector<double> tp;
  tp.reserve(s.size());
  tp.push_back(s.front());
  for (std::size_t i = 1; i + 1 < s.size(); ++i)
    if ((s[i] - s[i - 1] > 0) != (s[i + 1] - s[i] > 0)) tp.push_back(s[i]);
  tp.push_back(s.back());
  return tp;
}

std::vector<RainflowCycle> rainflow(const std::vector<double>& series) {
  const std::vector<double> tp = turning_points(series);
  std::vector<double> st;  // the stack is st[base..]
  st.reserve(tp.size());
  std::size_t base = 0;
  std::vector<RainflowCycle> out;
  for (double p : tp) {
    st.push_back(p);
    while (st.size() - base >= 3) {
      const std::size_t n = st.size();
      const double x = std::abs(st[n - 1] - st[n - 2]), y = std::abs(st[n - 2] - st[n - 3]);
      if (x < y) break;
      const double a = st[n - 3], b = st[n - 2];
      if (n - base == 3) {
        // Y holds the start point: a half cycle, the start moves on
        out.push_back({y, 0.5 * (a + b), 0.5});
        ++base;
      } else {
        // a closed cycle: its two points leave, the newest point drops down
        out.push_back({y, 0.5 * (a + b), 1.0});
        st[n - 3] = st[n - 1];
        st.resize(n - 2);
      }
    }
  }
  for (std::size_t i = base; i + 1 < st.size(); ++i)
    out.push_back({std::abs(st[i + 1] - st[i]), 0.5 * (st[i] + st[i + 1]), 0.5});
  return out;
}

double EffectHistory::max() const { return effect.empty() ? 0.0 : *std::max_element(effect.begin(), effect.end()); }
double EffectHistory::min() const { return effect.empty() ? 0.0 : *std::min_element(effect.begin(), effect.end()); }

double convoy_effect_at(const InfluenceLine& il, const IlPrefix& prefix, const ConvoyPattern& convoy, double r) {
  double E = 0;
  for (const auto& [a, p] : convoy.axles) E += p * il.at(r + a);
  for (const auto& seg : convoy.segments)
    if (seg[1] > 0) E += seg[2] / seg[1] * prefix.integral(r + seg[0], r + seg[0] + seg[1]);
  return E;
}

namespace {
// the convoy_extreme sweep bounds: every position where it touches the line
std::pair<double, double> bounds(const InfluenceLine& il, const ConvoyPattern& c) {
  require(il.s.size() >= 2, "convoy history: the influence line needs 2 samples");
  require(!c.axles.empty() || !c.segments.empty(), "convoy history: empty convoy");
  double a_min = 0, a_max = 0;
  for (const auto& [a, p] : c.axles) {
    a_min = std::min(a_min, a);
    a_max = std::max(a_max, a);
  }
  for (const auto& seg : c.segments) {
    a_min = std::min(a_min, seg[0]);
    a_max = std::max(a_max, seg[0] + seg[1]);
  }
  return {il.s.front() - a_max, il.s.back() - a_min};
}
}  // namespace

EffectHistory convoy_history(const InfluenceLine& il, const ConvoyPattern& convoy, int steps) {
  require(steps >= 1, "convoy history: steps >= 1");
  const auto [r0, r1] = bounds(il, convoy);
  const IlPrefix prefix(il);
  EffectHistory h;
  h.r.reserve(static_cast<std::size_t>(steps) + 1);
  h.effect.reserve(static_cast<std::size_t>(steps) + 1);
  for (int k = 0; k <= steps; ++k) {
    const double r = r0 + (r1 - r0) * k / steps;
    h.r.push_back(r);
    h.effect.push_back(convoy_effect_at(il, prefix, convoy, r));
  }
  return h;
}

EffectHistory convoy_history_exact(const InfluenceLine& il, const ConvoyPattern& convoy, int per_interval) {
  const auto [r0, r1] = bounds(il, convoy);
  if (per_interval < 0) per_interval = (il.cubic() || !convoy.segments.empty()) ? 8 : 0;
  std::vector<double> br;
  br.reserve(il.s.size() * (convoy.axles.size() + 2 * convoy.segments.size()) + 2);
  br.push_back(r0);
  br.push_back(r1);
  for (double s : il.s) {
    for (const auto& [a, p] : convoy.axles) br.push_back(s - a);
    for (const auto& seg : convoy.segments) {
      br.push_back(s - seg[0]);
      br.push_back(s - seg[0] - seg[1]);
    }
  }
  std::sort(br.begin(), br.end());
  const double tol = 1e-12 * std::max(1.0, r1 - r0);
  std::vector<double> pts;
  pts.reserve(br.size() * static_cast<std::size_t>(per_interval + 1));
  for (double r : br) {
    if (r < r0 - tol || r > r1 + tol) continue;
    if (!pts.empty() && r - pts.back() <= tol) continue;
    if (!pts.empty() && per_interval > 0) {
      const double a = pts.back();
      for (int k = 1; k <= per_interval; ++k) pts.push_back(a + (r - a) * k / (per_interval + 1));
    }
    pts.push_back(r);
  }
  const IlPrefix prefix(il);
  EffectHistory h;
  h.r = std::move(pts);
  h.effect.reserve(h.r.size());
  for (double r : h.r) h.effect.push_back(convoy_effect_at(il, prefix, convoy, r));
  return h;
}

// ---------------- EN 1991-2 §4.6 ----------------

namespace {
constexpr double kN = 1e3;
// Tables 4.6 / 4.7: the five lorry silhouettes (axle spacings), shared
const std::vector<double> kSpacing[5] = {{4.50}, {4.20, 1.30}, {3.20, 5.20, 1.30, 1.30}, {3.40, 6.00, 1.80},
                                         {4.80, 3.60, 4.40, 1.30}};
const std::vector<double> kFlm2[5] = {{90, 190}, {80, 140, 140}, {90, 180, 120, 120, 120}, {90, 190, 140, 140},
                                      {90, 180, 120, 110, 110}};
const std::vector<double> kFlm4[5] = {{70, 130}, {70, 120, 120}, {70, 150, 90, 90, 90}, {70, 140, 90, 90},
                                      {70, 130, 90, 80, 80}};
// Table 4.7 shares [%]: long distance, medium distance, local
constexpr double kShare[5][3] = {{20, 40, 80}, {5, 10, 5}, {50, 30, 5}, {15, 15, 5}, {10, 5, 5}};

ConvoyPattern lorry(int i, const std::vector<double>* loads) {
  require(i >= 1 && i <= 5, "fatigue lorry: index 1..5");
  const std::vector<double>& sp = kSpacing[i - 1];
  const std::vector<double>& q = loads[i - 1];
  ConvoyPattern c;
  double off = 0;
  for (std::size_t k = 0; k < q.size(); ++k) {
    if (k > 0) off += sp[k - 1];
    c.axles.push_back({off, q[k] * kN});
  }
  return c;
}
}  // namespace

ConvoyPattern flm1_tandem(int lane) {
  require(lane >= 1, "FLM1: lane rank >= 1");
  const double Q = lane == 1 ? 300 : lane == 2 ? 200 : lane == 3 ? 100 : 0;
  ConvoyPattern c;
  c.axles = {{0.0, 0.7 * Q * kN}, {1.2, 0.7 * Q * kN}};
  return c;
}

double flm1_udl(int lane, double lane_width) {
  require(lane >= 1 && lane_width > 0, "FLM1: lane rank >= 1 and a positive width");
  return 0.3 * (lane == 1 ? 9.0 : 2.5) * kN * lane_width;
}

ConvoyPattern flm2_lorry(int i) { return lorry(i, kFlm2); }
ConvoyPattern flm4_lorry(int i) { return lorry(i, kFlm4); }

double flm4_share(int i, FatigueTraffic t) {
  require(i >= 1 && i <= 5, "FLM4: lorry index 1..5");
  return kShare[i - 1][static_cast<int>(t)];
}

ConvoyPattern flm3_vehicle() {
  ConvoyPattern c;
  for (double a : {0.0, 1.2, 7.2, 8.4}) c.axles.push_back({a, 120 * kN});
  return c;
}

ConvoyPattern flm3_double(double gap) {
  require(gap >= 40.0, "FLM3: the follower is at least 40 m behind (EN 1991-2 §4.6.4(3))");
  ConvoyPattern c = flm3_vehicle();
  for (double a : {0.0, 1.2, 7.2, 8.4}) c.axles.push_back({a + gap, 36 * kN});
  return c;
}

}  // namespace keta

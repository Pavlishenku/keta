#include "keta/traffic.h"

namespace keta {

LaneDivision en_lane_division(double w) {
  require(w > 0, "lane division: width must be positive");
  LaneDivision d;
  if (w < 5.4) {
    d.n_lanes = 1;
    d.lane_width = 3.0;
    d.remaining = std::max(0.0, w - 3.0);
  } else if (w < 6.0) {
    d.n_lanes = 2;
    d.lane_width = w / 2;
    d.remaining = 0;
  } else {
    d.n_lanes = static_cast<int>(w / 3.0);
    d.lane_width = 3.0;
    d.remaining = w - 3.0 * d.n_lanes;
  }
  return d;
}

Lm1Factors Lm1Factors::class1() { return {}; }
Lm1Factors Lm1Factors::class2() {
  // French NA 2nd class: TS 0.9/0.8/0.8... use the published alpha set
  Lm1Factors f;
  f.aQ[0] = 0.9;
  f.aQ[1] = 0.8;
  f.aQ[2] = 0.8;
  f.aq1 = 0.7;
  f.aq_other = 1.0;
  return f;
}

namespace {
constexpr double kQk[3] = {300e3, 200e3, 100e3};  // TS axle loads per lane rank [N]
constexpr double kQ1k = 9.0e3;                    // lane-1 UDL [N/m^2]
constexpr double kQok = 2.5e3;                    // other lanes UDL [N/m^2]
}  // namespace

ConvoyPattern lm1_tandem(int rank, const Lm1Factors& f) {
  require(rank >= 1 && rank <= 3, "LM1 TS exists on lanes 1..3 only");
  double Q = f.aQ[rank - 1] * kQk[rank - 1];
  ConvoyPattern c;
  c.axles = {{0.0, Q}, {1.2, Q}};
  return c;
}

double lm1_udl(int rank, double lane_width, const Lm1Factors& f) {
  double q = rank == 1 ? f.aq1 * kQ1k : f.aq_other * kQok;
  return q * lane_width;
}

ConvoyPattern lm2_axle(double beta_q, bool french_na) {
  ConvoyPattern c;
  c.axles = {{0.0, beta_q * (french_na ? 320e3 : 400e3)}};
  return c;
}

double lm4_crowd(double lane_width) { return 5.0e3 * lane_width; }

Extremes lm1_lane_extreme(const InfluenceLine& il, int rank, double lane_width,
                          const Lm1Factors& f, int steps) {
  Extremes ts = rank <= 3 ? convoy_extreme(il, lm1_tandem(rank, f), steps) : Extremes{};
  Extremes udl = udl_extreme(il, lm1_udl(rank, lane_width, f));
  Extremes ex;
  ex.max = ts.max + udl.max;
  ex.min = ts.min + udl.min;
  ex.r_max = ts.r_max;
  ex.r_min = ts.r_min;
  return ex;
}

// ---------------- Fascicule 61 ----------------

double fasc61_a_of_l(double l) {
  require(l > 0, "A(l): loaded length must be positive");
  // 0.23 + 36/(l+12) t/m^2 -> N/m^2 (1 t/m^2 = 9810 N/m^2 g-standard 9.81;
  // use 9.80665 kN base: 1 t = 9806.65 N)
  return (0.23 + 36.0 / (l + 12.0)) * 9806.65;
}

double fasc61_a1(int bridge_class, int n_loaded) {
  require(bridge_class >= 1 && bridge_class <= 3, "a1: class 1..3");
  require(n_loaded >= 1, "a1: at least one loaded lane");
  // Fascicule 61 titre II, art. 4.2.2: class x loaded-lane count (>= 5 as 5);
  // 2nd- and 3rd-class bridges are loaded on two lanes at most.
  static const double t1[5] = {1.0, 1.0, 0.9, 0.75, 0.7};
  static const double t2[2] = {1.0, 0.9};
  static const double t3[2] = {0.9, 0.8};
  switch (bridge_class) {
    case 1: return t1[std::min(n_loaded, 5) - 1];
    case 2:
      require(n_loaded <= 2, "a1: 2nd-class bridges carry at most 2 loaded lanes");
      return t2[n_loaded - 1];
    default:
      require(n_loaded <= 2, "a1: 3rd-class bridges carry at most 2 loaded lanes");
      return t3[n_loaded - 1];
  }
}

double fasc61_a2(int bridge_class, double lane_width) {
  static const double v0[3] = {3.5, 3.0, 2.75};
  require(bridge_class >= 1 && bridge_class <= 3, "a2: class 1..3");
  require(lane_width > 0, "a2: lane width must be positive");
  return v0[bridge_class - 1] / lane_width;
}

double fasc61_a_lane_load(double l, int bridge_class, double lane_width) {
  double a1 = fasc61_a1(bridge_class, 1);
  double aA = a1 * fasc61_a_of_l(l);
  // the floor (400 - 0.2 l) kg/m^2 = (400 - 0.2 l) * 9.80665 N/m^2
  const double floor_pa = (400.0 - 0.2 * l) * 9.80665;
  double A2 = std::max(aA, floor_pa);
  return A2 * fasc61_a2(bridge_class, lane_width) * lane_width;
}

ConvoyPattern fasc61_bc() {
  // One file of two 30 t trucks: axles (front to rear) 6 t at 0, 12 t at 4.5 m,
  // 12 t at 6.0 m; second truck 4.5 m behind the first's rear axle.
  ConvoyPattern c;
  double t = 9806.65;
  auto truck = [&](double x0) {
    c.axles.emplace_back(x0 + 0.0, 6 * t);
    c.axles.emplace_back(x0 + 4.5, 12 * t);
    c.axles.emplace_back(x0 + 6.0, 12 * t);
  };
  truck(0.0);
  truck(6.0 + 4.5);
  return c;
}

ConvoyPattern fasc61_bt() {
  // Two 16 t-axle tandems (1.35 m spacing) — one file.
  ConvoyPattern c;
  double t = 9806.65;
  c.axles = {{0.0, 16 * t}, {1.35, 16 * t}};
  return c;
}

ConvoyPattern fasc61_mc120() {
  // Mc120: 110 t on two tracks of 6.10 m -> uniform segment.
  ConvoyPattern c;
  c.segments.push_back({0.0, 6.10, 110 * 9806.65});
  return c;
}

double fasc61_delta(double L, double G, double S) {
  require(L > 0, "delta: element length must be positive");
  double d = 1.0 + 0.4 / (1.0 + 0.2 * L);
  if (S > 0) d += 0.6 / (1.0 + 4.0 * G / S);
  return d;
}

}  // namespace keta

namespace keta {

ConvoyPattern rail_lm71(double alpha) {
  require(alpha >= 0.7 && alpha <= 1.5, "LM71: alpha out of the standard range");
  ConvoyPattern c;
  for (double x : {0.0, 1.6, 3.2, 4.8}) c.axles.emplace_back(x, alpha * 250e3);
  // 80 kN/m tails starting 0.8 m outside the end axles (finite 20 m lengths;
  // the sign-matched UDL beyond is the caller's udl_extreme business).
  c.segments.push_back({-0.8 - 20.0, 20.0, alpha * 80e3 * 20.0});
  c.segments.push_back({4.8 + 0.8, 20.0, alpha * 80e3 * 20.0});
  return c;
}

ConvoyPattern rail_sw0(double alpha) {
  ConvoyPattern c;
  c.segments.push_back({0.0, 15.0, alpha * 133e3 * 15.0});
  c.segments.push_back({15.0 + 5.3, 15.0, alpha * 133e3 * 15.0});
  return c;
}

ConvoyPattern rail_sw2() {
  ConvoyPattern c;
  c.segments.push_back({0.0, 25.0, 150e3 * 25.0});
  c.segments.push_back({25.0 + 7.0, 25.0, 150e3 * 25.0});
  return c;
}

double rail_unloaded_train() { return 10e3; }

ConvoyPattern rail_hslm_b(int n, double d) {
  require(n > 0 && d > 0, "HSLM-B: bad axle data");
  ConvoyPattern c;
  for (int i = 0; i < n; ++i) c.axles.emplace_back(i * d, 170e3);
  return c;
}

}  // namespace keta

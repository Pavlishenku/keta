#include "keta/rail_fatigue.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
constexpr double kNperT = 10e3;  // Annex D: 10 kN per tonne of axle mass [N/t]

struct Car {
  double axle_t;
  int n_axles;
  double intra, bogie_gap;
  std::vector<double> offsets() const {
    std::vector<double> o;
    const int half = n_axles / 2;
    double x = 0;
    for (int k = 0; k < half; ++k, x += intra) o.push_back(x);
    double xb = o.back() + bogie_gap;
    for (int k = 0; k < half; ++k, xb += intra) o.push_back(xb);
    return o;
  }
};
struct Regular {
  Car loco;
  bool rear_loco;
  Car wagon;
  int n_wagons;
  double d_lw, d_ww;  // couplings: last loco axle -> first wagon axle, wagon -> wagon
};
struct Row {
  int number;
  double mass_t, v_kmh;
  bool regular;
  Regular reg;
  std::vector<std::pair<double, double>> axles_t;  // explicit (offset, tonnes)
};

// Axle files of the trains whose vehicles differ along the train (offset from
// the first axle [m], tonnes), read on figures D.4, D.6, D.9 and D.10
const std::vector<std::pair<double, double>> kT4 = {
    {0.0, 17},   {3.0, 17},   {14.0, 17},  {17.0, 17},  {20.3, 17},  {23.3, 17},  {39.0, 17},  {42.0, 17},
    {57.7, 17},  {60.7, 17},  {76.4, 17},  {79.4, 17},  {95.1, 17},  {98.1, 17},  {113.8, 17}, {116.8, 17},
    {132.5, 17}, {135.5, 17}, {151.2, 17}, {154.2, 17}, {169.9, 17}, {172.9, 17}, {188.6, 17}, {191.6, 17},
    {207.3, 17}, {210.3, 17}, {213.6, 17}, {216.6, 17}, {227.6, 17}, {230.6, 17}};
const std::vector<std::pair<double, double>> kT6 = {
    {0.0, 22.5},   {2.1, 22.5},   {4.2, 22.5},   {8.6, 22.5},   {10.7, 22.5},  {12.8, 22.5},  {16.7, 7},
    {23.2, 7},     {27.0, 7},     {33.5, 7},     {37.2, 22.5},  {39.0, 22.5},  {51.8, 22.5},  {53.6, 22.5},
    {57.3, 7},     {63.8, 7},     {67.3, 22.5},  {69.1, 22.5},  {77.1, 22.5},  {78.9, 22.5},  {82.1, 22.5},
    {83.9, 22.5},  {91.9, 22.5},  {93.7, 22.5},  {97.2, 7},     {103.7, 7},    {107.4, 22.5}, {109.2, 22.5},
    {122.0, 22.5}, {123.8, 22.5}, {127.4, 22.5}, {129.2, 22.5}, {142.0, 22.5}, {143.8, 22.5}, {147.4, 22.5},
    {149.2, 22.5}, {162.0, 22.5}, {163.8, 22.5}, {167.5, 7},    {174.0, 7},    {177.8, 7},    {184.3, 7},
    {187.8, 22.5}, {189.6, 22.5}, {197.6, 22.5}, {199.4, 22.5}, {202.6, 22.5}, {204.4, 22.5}, {212.4, 22.5},
    {214.2, 22.5}, {217.7, 7},    {224.2, 7},    {227.9, 22.5}, {229.7, 22.5}, {242.5, 22.5}, {244.3, 22.5},
    {247.7, 22.5}, {249.5, 22.5}, {257.5, 22.5}, {259.3, 22.5}, {262.8, 7},    {269.3, 7},    {273.1, 7},
    {279.6, 7},    {283.1, 22.5}, {284.9, 22.5}, {292.9, 22.5}, {294.7, 22.5}, {297.9, 22.5}, {299.7, 22.5},
    {307.7, 22.5}, {309.5, 22.5}, {312.9, 22.5}, {314.7, 22.5}, {327.5, 22.5}, {329.3, 22.5}};
const std::vector<std::pair<double, double>> kT9 = {
    {0.0, 13},   {2.5, 13},   {16.5, 13},  {19.0, 13},  {23.3, 11},  {25.8, 11},  {37.3, 11},  {39.8, 11},
    {44.1, 13},  {46.6, 13},  {60.6, 13},  {63.1, 13},  {67.4, 13},  {69.9, 13},  {83.9, 13},  {86.4, 13},
    {90.7, 11},  {93.2, 11},  {104.7, 11}, {107.2, 11}, {111.5, 13}, {114.0, 13}, {128.0, 13}, {130.5, 13}};
// train 10 (figure D.10): cars A B B A A B B A of 16.2 m, bogie axles at 1.75
// and 4.15 m from each car end; A carries 15, 10, 10, 15 t, B four axles of 10 t
const std::vector<std::pair<double, double>> kT10 = {
    {0.0, 15},   {2.4, 10},   {10.3, 10},  {12.7, 15},  {16.2, 10},  {18.6, 10},  {26.5, 10},  {28.9, 10},
    {32.4, 10},  {34.8, 10},  {42.7, 10},  {45.1, 10},  {48.6, 15},  {51.0, 10},  {58.9, 10},  {61.3, 15},
    {64.8, 15},  {67.2, 10},  {75.1, 10},  {77.5, 15},  {81.0, 10},  {83.4, 10},  {91.3, 10},  {93.7, 10},
    {97.2, 10},  {99.6, 10},  {107.5, 10}, {109.9, 10}, {113.4, 15}, {115.8, 10}, {123.7, 10}, {126.1, 15}};

const Row& row(int n) {
  static const std::vector<Row> rows = {
      {1, 663, 200, true, {{22.5, 6, 2.2, 6.9}, false, {11, 4, 2.6, 11.5}, 12, 3.2, 3.6}, {}},
      {2, 530, 160, true, {{22.5, 4, 3.3, 6.7}, false, {11, 4, 2.5, 16.5}, 10, 3.9, 5.0}, {}},
      {3, 940, 250, true, {{20, 4, 3.0, 8.46}, true, {15, 4, 2.5, 16.5}, 13, 4.45, 4.90}, {}},
      {4, 510, 250, false, {}, kT4},
      {5, 2160, 80, true, {{22.5, 6, 2.1, 4.4}, false, {22.5, 6, 1.8, 5.7}, 15, 4.0, 4.0}, {}},
      {6, 1431, 100, false, {}, kT6},
      {7, 1035, 120, true, {{22.5, 6, 2.2, 6.9}, false, {22.5, 4, 1.8, 11.0}, 10, 3.0, 3.2}, {}},
      {8, 1035, 100, true, {{22.5, 6, 2.2, 6.9}, false, {22.5, 2, 0.0, 5.5}, 20, 3.5, 4.2}, {}},
      {9, 296, 120, false, {}, kT9},
      {10, 360, 120, false, {}, kT10},
      {11, 1135, 120, true, {{22.5, 6, 2.2, 6.9}, false, {25, 4, 2.0, 11.0}, 10, 2.9, 3.0}, {}},
      {12, 1135, 100, true, {{22.5, 6, 2.2, 6.9}, false, {25, 2, 0.0, 5.5}, 20, 3.5, 4.2}, {}},
  };
  require(n >= 1 && n <= 12, "Annex D: service train number 1..12");
  return rows[static_cast<std::size_t>(n - 1)];
}
}  // namespace

double ServiceTrain::total_load() const {
  double t = 0;
  for (const auto& [a, p] : pattern.axles) t += p;
  return t;
}
double ServiceTrain::length() const {
  if (pattern.axles.empty()) return 0;
  double lo = pattern.axles.front().first, hi = lo;
  for (const auto& [a, p] : pattern.axles) {
    lo = std::min(lo, a);
    hi = std::max(hi, a);
  }
  return hi - lo;
}

ServiceTrain annex_d_train(int number) {
  const Row& r = row(number);
  ServiceTrain t;
  t.number = r.number;
  t.mass_t = r.mass_t;
  t.v_max_kmh = r.v_kmh;
  auto& ax = t.pattern.axles;
  if (!r.regular) {
    for (const auto& [o, m] : r.axles_t) ax.push_back({o, m * kNperT});
    return t;
  }
  const Regular& g = r.reg;
  for (double o : g.loco.offsets()) ax.push_back({o, g.loco.axle_t * kNperT});
  const std::vector<double> wo = g.wagon.offsets();
  double x = ax.back().first + g.d_lw;
  for (int w = 0; w < g.n_wagons; ++w) {
    for (double o : wo) ax.push_back({x + o, g.wagon.axle_t * kNperT});
    x = ax.back().first + g.d_ww;
  }
  if (g.rear_loco) {
    x = ax.back().first + g.d_lw;
    for (double o : g.loco.offsets()) ax.push_back({x + o, g.loco.axle_t * kNperT});
  }
  return t;
}

std::vector<RailMixEntry> annex_d_mix(RailFatigueMix mix) {
  std::vector<std::pair<int, int>> rows;
  switch (mix) {
    case RailFatigueMix::Standard: rows = {{1, 12}, {2, 12}, {3, 5}, {4, 5}, {5, 7}, {6, 12}, {7, 8}, {8, 6}}; break;
    case RailFatigueMix::Heavy250: rows = {{5, 6}, {6, 13}, {11, 16}, {12, 16}}; break;
    case RailFatigueMix::Light: rows = {{1, 10}, {2, 5}, {5, 2}, {9, 190}}; break;
  }
  std::vector<RailMixEntry> out;
  for (const auto& [n, d] : rows) out.push_back({n, d, d * row(n).mass_t * 365.0 / 1e6});
  return out;
}

double annex_d_mix_tonnage(RailFatigueMix mix) {
  double s = 0;
  for (const RailMixEntry& e : annex_d_mix(mix)) s += e.tonnage_mt;
  return s;
}

double rail_n0_lower(double L) {
  require(L >= 4 && L <= 100, "EN 1991-2 Figure 6.10: 4 <= L_phi <= 100 m");
  return L <= 20 ? 80 / L : 23.58 * std::pow(L, -0.592);
}

std::pair<double, double> rail_real_train_increments(double v, double L, double n0) {
  require(std::isfinite(v) && v >= 0 && L > 0 && n0 > 0, "Annex C: v >= 0, L_phi > 0, n0 > 0");
  const double K = v / (2 * L * n0);
  const double p1 = K < 0.76 ? K / (1 - K + std::pow(K, 4)) : 1.325;
  const double alpha = std::min(v / 22.0, 1.0);
  const double p2 = std::max(
      0.0, alpha / 100 * (56 * std::exp(-(L / 10) * (L / 10)) + 50 * (L * n0 / 80 - 1) * std::exp(-(L / 20) * (L / 20))));
  return {p1, p2};
}

double rail_fatigue_dynamic_factor(double v, double L, double n0) {
  const auto [p1, p2] = rail_real_train_increments(v, L, n0);
  return 1 + 0.5 * (p1 + 0.5 * p2);
}

}  // namespace keta

#include "keta/traffic_multi.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace keta {

MobileVehicle MobileVehicle::from_margins(std::vector<std::pair<double, double>> axles, double occ_back,
                                          double occ_front, double rem_back, double rem_front) {
  require(!axles.empty(), "vehicle: at least one axle");
  require(occ_back >= 0 && occ_front >= 0 && rem_back >= 0 && rem_front >= 0, "vehicle: margins >= 0");
  double lo = axles.front().first, hi = lo;
  for (const auto& [a, p] : axles) {
    lo = std::min(lo, a);
    hi = std::max(hi, a);
  }
  MobileVehicle v;
  v.axles = std::move(axles);
  v.occ_lo = lo - occ_back;
  v.occ_hi = hi + occ_front;
  v.rem_lo = lo - rem_back;
  v.rem_hi = hi + rem_front;
  return v;
}

// ---------------- sign-split prefix ----------------

IlSignPrefix::IlSignPrefix(const InfluenceLine& il) : il_(il) {
  const std::size_t n = il.s.size();
  require(n >= 2, "influence line: 2 samples at least");
  x_.reserve(2 * n);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    x_.push_back(il.s[i]);
    const double a = il.s[i], b = il.s[i + 1], h = b - a;
    if (h <= 0) continue;
    if (!il.cubic()) {
      const double ea = il.eta[i], eb = il.eta[i + 1];
      if ((ea < 0 && eb > 0) || (ea > 0 && eb < 0)) x_.push_back(a + h * ea / (ea - eb));
      continue;
    }
    // cubic piece: sign changes on 8 sub-intervals, each root by bisection
    constexpr int m = 8;
    double xl = a, fl = il.on(i, a);
    for (int k = 1; k <= m; ++k) {
      const double xr = a + h * k / m, fr = il.on(i, xr);
      if ((fl < 0 && fr > 0) || (fl > 0 && fr < 0)) {
        double lo = xl, hi = xr, flo = fl;
        for (int it = 0; it < 80 && hi - lo > 1e-15 * std::max(1.0, std::abs(hi)); ++it) {
          const double mid = 0.5 * (lo + hi), fm = il.on(i, mid);
          if ((fm < 0) == (flo < 0) && fm != 0) {
            lo = mid;
            flo = fm;
          } else {
            hi = mid;
          }
        }
        x_.push_back(0.5 * (lo + hi));
      }
      xl = xr;
      fl = fr;
    }
  }
  x_.push_back(il.s.back());
  cpos_.assign(x_.size(), 0.0);
  cneg_.assign(x_.size(), 0.0);
  for (std::size_t k = 0; k + 1 < x_.size(); ++k) {
    const double I = x_[k + 1] > x_[k] ? il.integral(x_[k], x_[k + 1]) : 0.0;
    cpos_[k + 1] = cpos_[k] + std::max(I, 0.0);
    cneg_[k + 1] = cneg_[k] + std::min(I, 0.0);
  }
}

std::pair<double, double> IlSignPrefix::integrals(double a, double b) const {
  a = std::max(a, x_.front());
  b = std::min(b, x_.back());
  if (b <= a) return {0.0, 0.0};
  // pieces of constant sign: [x_k, x_k+1]
  auto piece = [&](double x) {
    std::size_t k = static_cast<std::size_t>(std::upper_bound(x_.begin(), x_.end(), x) - x_.begin());
    return k == 0 ? 0 : std::min(k - 1, x_.size() - 2);
  };
  const std::size_t ka = piece(a), kb = piece(b);
  auto split = [](double I) { return std::make_pair(std::max(I, 0.0), std::min(I, 0.0)); };
  if (ka == kb) return split(il_.integral(a, b));
  const auto [pa, na] = split(il_.integral(a, x_[ka + 1]));
  const auto [pb, nb] = split(il_.integral(x_[kb], b));
  return {pa + (cpos_[kb] - cpos_[ka + 1]) + pb, na + (cneg_[kb] - cneg_[ka + 1]) + nb};
}

// ---------------- the exact placement ----------------

namespace {
struct Place {
  int v;
  double lo, hi, r, score;
};
// the best sum of one score per vehicle over non-overlapping occupancies
// (touching allowed), by DP over (vehicles placed, rightmost end)
double best_placement(const std::vector<Place>& pl, int k, std::vector<double>& r_out) {
  const std::size_t P = pl.size();
  const std::size_t M = std::size_t(1) << k;
  constexpr double NEG = -std::numeric_limits<double>::infinity();
  std::vector<std::size_t> by_lo(P), by_hi(P);
  std::iota(by_lo.begin(), by_lo.end(), 0);
  std::iota(by_hi.begin(), by_hi.end(), 0);
  std::sort(by_lo.begin(), by_lo.end(), [&](std::size_t a, std::size_t b) {
    return pl[a].lo != pl[b].lo ? pl[a].lo < pl[b].lo : pl[a].hi < pl[b].hi;
  });
  std::sort(by_hi.begin(), by_hi.end(), [&](std::size_t a, std::size_t b) { return pl[a].hi < pl[b].hi; });
  std::vector<double> g(P * M, NEG);       // best total ending with placement p, vehicles `mask` placed
  std::vector<int> pred(P * M, -1);        // the placement before p in that configuration
  std::vector<double> avail(M, NEG);       // best total of `mask` among the placements ended so far
  std::vector<int> avail_arg(M, -1);
  avail[0] = 0;
  std::size_t h = 0;
  constexpr double tol = 1e-9;
  for (std::size_t p : by_lo) {
    while (h < P && pl[by_hi[h]].hi <= pl[p].lo + tol) {
      const std::size_t q = by_hi[h++];
      for (std::size_t mask = 0; mask < M; ++mask)
        if (g[q * M + mask] > avail[mask]) {
          avail[mask] = g[q * M + mask];
          avail_arg[mask] = static_cast<int>(q);
        }
    }
    const std::size_t bit = std::size_t(1) << pl[p].v;
    for (std::size_t mask = 0; mask < M; ++mask) {
      if (mask & bit || avail[mask] == NEG) continue;
      g[p * M + (mask | bit)] = avail[mask] + pl[p].score;
      pred[p * M + (mask | bit)] = avail_arg[mask];
    }
  }
  const std::size_t full = M - 1;
  double best = NEG;
  std::size_t bp = P;
  for (std::size_t p = 0; p < P; ++p)
    if (g[p * M + full] > best) {
      best = g[p * M + full];
      bp = p;
    }
  require(bp < P, "multi-vehicle: the deck cannot hold every vehicle without overlapping occupancies");
  r_out.assign(static_cast<std::size_t>(k), 0.0);
  std::size_t mask = full;
  for (int p = static_cast<int>(bp); p >= 0;) {
    const Place& x = pl[static_cast<std::size_t>(p)];
    r_out[static_cast<std::size_t>(x.v)] = x.r;
    const int prev = pred[static_cast<std::size_t>(p) * M + mask];
    mask &= ~(std::size_t(1) << x.v);
    p = prev;
  }
  return best;
}
}  // namespace

MultiVehicleResult multi_vehicle_extreme(const InfluenceLine& il, const std::vector<MobileVehicle>& vehicles,
                                         double udl, int steps) {
  require(il.s.size() >= 2, "multi-vehicle: the influence line needs 2 samples");
  require(steps >= 1, "multi-vehicle: steps >= 1");
  const int k = static_cast<int>(vehicles.size());
  require(k >= 1 && k <= 10, "multi-vehicle: 1 to 10 vehicles (the exact search is O(P k 2^k))");
  for (const MobileVehicle& v : vehicles) {
    require(!v.axles.empty() && v.occ_lo <= v.occ_hi && v.rem_lo <= v.rem_hi, "multi-vehicle: bad vehicle");
    require(v.rem_lo >= v.occ_lo - 1e-12 && v.rem_hi <= v.occ_hi + 1e-12,
            "multi-vehicle: the UDL removal interval must lie inside the occupancy");
  }
  const double s0 = il.s.front(), s1 = il.s.back();
  const IlSignPrefix sp(il);
  const auto [Apos, Aneg] = sp.integrals(s0, s1);
  std::vector<Place> pmax, pmin;
  pmax.reserve(static_cast<std::size_t>(k) * static_cast<std::size_t>(steps + 1));
  pmin.reserve(pmax.capacity());
  for (int v = 0; v < k; ++v) {
    const MobileVehicle& mv = vehicles[static_cast<std::size_t>(v)];
    for (int g = 0; g <= steps; ++g) {
      const double r = s0 + (s1 - s0) * g / steps;
      double E = 0;
      for (const auto& [a, p] : mv.axles) E += p * il.at(r + a);
      double rp = 0, rn = 0;
      if (udl != 0) std::tie(rp, rn) = sp.integrals(r + mv.rem_lo, r + mv.rem_hi);
      pmax.push_back({v, r + mv.occ_lo, r + mv.occ_hi, r, E - udl * rp});
      pmin.push_back({v, r + mv.occ_lo, r + mv.occ_hi, r, -(E - udl * rn)});
    }
  }
  MultiVehicleResult out;
  out.extremes.max = udl * Apos + best_placement(pmax, k, out.max_r);
  out.extremes.min = udl * Aneg - best_placement(pmin, k, out.min_r);
  out.extremes.r_max = out.max_r.front();
  out.extremes.r_min = out.min_r.front();
  return out;
}

}  // namespace keta

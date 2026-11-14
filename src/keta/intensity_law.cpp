#include "keta/intensity_law.h"

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>

#include "keta/traffic_multi.h"

namespace keta {

IntensityLaw IntensityLaw::constant(double c) {
  require(std::isfinite(c), "intensity law: finite constant");
  IntensityLaw l;
  l.shape_ = [c](double) { return c; };
  return l;
}

IntensityLaw IntensityLaw::tabulated(std::vector<std::pair<double, double>> pts) {
  require(!pts.empty(), "intensity law: empty table");
  for (std::size_t i = 0; i < pts.size(); ++i) {
    require(std::isfinite(pts[i].first) && std::isfinite(pts[i].second) && pts[i].first >= 0,
            "intensity law: finite samples, L >= 0");
    if (i) require(pts[i].first > pts[i - 1].first, "intensity law: L strictly increasing");
  }
  IntensityLaw l;
  l.shape_ = [p = std::move(pts)](double L) {
    if (L <= p.front().first) return p.front().second;
    if (L >= p.back().first) return p.back().second;
    const auto it = std::upper_bound(p.begin(), p.end(), L, [](double x, const auto& q) { return x < q.first; });
    const auto& [x1, y1] = *it;
    const auto& [x0, y0] = *(it - 1);
    return y0 + (L - x0) / (x1 - x0) * (y1 - y0);
  };
  return l;
}

IntensityLaw IntensityLaw::sqrt_law(double a) {
  IntensityLaw l;
  l.shape_ = [a](double L) { return a * std::sqrt(L); };
  return l;
}

IntensityLaw IntensityLaw::power(double b, double p) {
  IntensityLaw l;
  l.shape_ = [b, p](double L) {
    require(L > 0, "intensity law b / L^p: L > 0 (bound it with l_min)");
    return b / std::pow(L, p);
  };
  return l;
}

IntensityLaw IntensityLaw::log_law(double b, double c) {
  IntensityLaw l;
  l.shape_ = [b, c](double L) {
    require(L > 0, "intensity law b + c ln L: L > 0 (bound it with l_min)");
    return b + c * std::log(L);
  };
  return l;
}

IntensityLaw IntensityLaw::custom(std::function<double(double)> f) {
  require(static_cast<bool>(f), "intensity law: empty function");
  IntensityLaw l;
  l.shape_ = std::move(f);
  return l;
}

IntensityLaw& IntensityLaw::bounds(std::optional<double> lmin, std::optional<double> lmax,
                                   std::optional<double> total) {
  require(!lmin || (std::isfinite(*lmin) && *lmin >= 0), "intensity law: l_min >= 0");
  require(!lmax || (std::isfinite(*lmax) && *lmax >= 0), "intensity law: l_max >= 0");
  require(!lmin || !lmax || *lmin <= *lmax, "intensity law: l_min <= l_max");
  require(!total || (std::isfinite(*total) && *total >= 0), "intensity law: max_total >= 0");
  l_min_ = lmin;
  l_max_ = lmax;
  max_total_ = total;
  return *this;
}

double IntensityLaw::at(double L) const {
  require(std::isfinite(L) && L >= 0, "intensity law: loaded length >= 0");
  if (l_min_) L = std::max(L, *l_min_);
  if (l_max_) L = std::min(L, *l_max_);
  double v = shape_(L);
  require(std::isfinite(v), "intensity law: non-finite value");
  if (max_total_ && L > 0 && std::abs(v) * L > *max_total_) v = std::copysign(*max_total_ / L, v);
  return v;
}

std::vector<LoadedZone> sign_zones(const InfluenceLine& il, int sign) {
  require(sign == 1 || sign == -1, "sign zones: sign +1 or -1");
  const IlSignPrefix sp(il);
  // the samples with the zero crossings inserted (as IlSignPrefix finds
  // them); a zone is a maximal run of pieces of the sign
  std::vector<LoadedZone> out;
  const std::size_t n = il.s.size();
  auto val = [&](double x) { return il.at(x); };
  std::vector<double> pts;
  pts.reserve(2 * n);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    pts.push_back(il.s[i]);
    const double a = il.s[i], b = il.s[i + 1];
    if (b <= a) continue;
    const int m = il.cubic() ? 8 : 1;
    double xl = a, fl = il.on(i, a);
    for (int k = 1; k <= m; ++k) {
      const double xr = a + (b - a) * k / m, fr = il.on(i, xr);
      if ((fl < 0 && fr > 0) || (fl > 0 && fr < 0)) {
        double lo = xl, hi = xr, flo = fl;
        for (int it = 0; it < 80; ++it) {
          const double mid = 0.5 * (lo + hi), fm = il.on(i, mid);
          if ((fm < 0) == (flo < 0) && fm != 0) {
            lo = mid;
            flo = fm;
          } else {
            hi = mid;
          }
        }
        pts.push_back(0.5 * (lo + hi));
      }
      xl = xr;
      fl = fr;
    }
  }
  pts.push_back(il.s.back());
  for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
    if (pts[k + 1] <= pts[k]) continue;
    const double mid = val(0.5 * (pts[k] + pts[k + 1]));
    if (!(sign > 0 ? mid > 0 : mid < 0)) continue;
    if (!out.empty() && out.back().b == pts[k])
      out.back().b = pts[k + 1];
    else
      out.push_back({pts[k], pts[k + 1], 0});
  }
  for (LoadedZone& z : out) {
    const auto [p, q] = sp.integrals(z.a, z.b);
    z.area = sign > 0 ? p : q;
  }
  return out;
}

IntensityExtreme intensity_udl_extreme(const InfluenceLine& il, const IntensityLaw& law, int sign) {
  const std::vector<LoadedZone> zs = sign_zones(il, sign);
  const int n = static_cast<int>(zs.size());
  require(n <= 24, "intensity extreme: more than 24 zones of one sign (2^n subsets)");
  IntensityExtreme best;
  // Gray code: subset g(k) = k ^ (k >> 1); step k flips the zone of the
  // lowest set bit of k
  double L = 0, A = 0;
  std::uint32_t mask = 0, best_mask = 0;
  const std::uint32_t N = std::uint32_t(1) << n;
  for (std::uint32_t k = 1; k < N; ++k) {
    const int z = std::countr_zero(k);
    const std::uint32_t bit = std::uint32_t(1) << z;
    const double s = (mask & bit) ? -1.0 : 1.0;
    mask ^= bit;
    L += s * (zs[static_cast<std::size_t>(z)].b - zs[static_cast<std::size_t>(z)].a);
    A += s * zs[static_cast<std::size_t>(z)].area;
    const double E = law.at(std::max(L, 0.0)) * A;
    if (sign > 0 ? E > best.value : E < best.value) {
      best.value = E;
      best.loaded_length = L;
      best_mask = mask;
    }
  }
  // the winner re-summed from its zones (no Gray-walk round-off)
  if (best_mask) {
    double Lb = 0, Ab = 0;
    for (int z = 0; z < n; ++z)
      if (best_mask & (std::uint32_t(1) << z)) {
        best.zones.push_back(z);
        Lb += zs[static_cast<std::size_t>(z)].b - zs[static_cast<std::size_t>(z)].a;
        Ab += zs[static_cast<std::size_t>(z)].area;
      }
    best.loaded_length = Lb;
    best.value = law.at(Lb) * Ab;
  }
  return best;
}

}  // namespace keta

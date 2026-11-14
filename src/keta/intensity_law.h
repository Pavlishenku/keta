// Length-dependent traffic intensity laws and the exact choice of the loaded
// zones.
//
// A law v(L) of the loaded length L: constant, tabulated (linear, held flat
// outside), a sqrt(L), b L^-p, b + c ln L, or a caller's function; then two
// optional bounds — the length clamped to [l_min, l_max] before the shape
// (l_min = l_max: an imposed length) and the total v(L) L capped (v reduced
// to max_total / L beyond). Read as a UDL density or a scale of point loads,
// the caller's choice.
//
// When the density depends on the loaded length, loading every favourable
// zone is not the worst case: E = v(L(S)) A(S) over the subsets S of the
// sign-matched zones of the influence line (lengths and areas). The exact
// maximum is found by enumerating the subsets in Gray-code order — one zone
// enters or leaves per step, so L(S) and A(S) are updated in O(1): O(2^n)
// for n zones (n <= 24, a typed error beyond). keta's SURCH does the same
// for ST1's density tables (repli au-delà de 12 zones).
#pragma once
#include <functional>
#include <optional>
#include <vector>

#include "keta/influence.h"

namespace keta {

class IntensityLaw {
 public:
  static IntensityLaw constant(double c);
  static IntensityLaw tabulated(std::vector<std::pair<double, double>> points);  // (L, v), L strictly increasing
  static IntensityLaw sqrt_law(double a);           // a sqrt(L)
  static IntensityLaw power(double b, double p);    // b / L^p (L > 0)
  static IntensityLaw log_law(double b, double c);  // b + c ln L (L > 0)
  static IntensityLaw custom(std::function<double(double)> f);
  // bounds: clamp of L, cap of v L (nullopt = none)
  IntensityLaw& bounds(std::optional<double> l_min, std::optional<double> l_max, std::optional<double> max_total);
  double at(double loaded_length) const;

 private:
  std::function<double(double)> shape_;
  std::optional<double> l_min_, l_max_, max_total_;
};

struct LoadedZone {
  double a = 0, b = 0;  // abscissae of the zone
  double area = 0;      // integral of eta over it (signed)
};
// the maximal zones of one sign (+1 / -1) of the line (zero crossings exact
// for a linear line, bisection on cubic pieces)
std::vector<LoadedZone> sign_zones(const InfluenceLine& il, int sign);

struct IntensityExtreme {
  double value = 0;            // v(L(S)) A(S) at the best subset (0: nothing loaded)
  double loaded_length = 0;
  std::vector<int> zones;      // indices into sign_zones(il, sign)
};
// the most unfavourable subset of the zones of `sign` for a UDL of density
// law v(L) per unit width (value = max of v A for sign +1, min for -1)
IntensityExtreme intensity_udl_extreme(const InfluenceLine& il, const IntensityLaw& law, int sign);

}  // namespace keta

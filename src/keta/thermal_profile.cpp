// Non-linear temperature profile decomposition — see thermal_profile.h.
#include "keta/thermal_profile.h"

#include <algorithm>
#include <cmath>

#include "keta/core.h"

namespace keta {

double profile_temperature(const std::vector<TempPoint>& p, double z) {
  require(!p.empty(), "thermal profile: no point");
  if (z <= p.front().z) return p.front().dT;
  if (z >= p.back().z) return p.back().dT;
  const auto it = std::upper_bound(p.begin(), p.end(), z, [](double v, const TempPoint& q) { return v < q.z; });
  const TempPoint& b = *it;
  const TempPoint& a = *(it - 1);
  return a.dT + (b.dT - a.dT) * (z - a.z) / (b.z - a.z);
}

EquivalentTemperature equivalent_temperature(const std::vector<ThermalStrip>& strips,
                                             const std::vector<TempPoint>& profile) {
  require(!strips.empty(), "thermal profile: no strip");
  for (std::size_t i = 1; i < profile.size(); ++i)
    require(profile[i].z > profile[i - 1].z, "thermal profile: depths must increase");
  static const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
  static const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
  // integrate f(z) E b(z) over every strip, cut at the profile points
  auto integrate = [&](auto&& f) {
    double s = 0;
    for (const ThermalStrip& st : strips) {
      require(st.z1 > st.z0, "thermal profile: empty strip");
      std::vector<double> cuts = {st.z0, st.z1};
      for (const TempPoint& q : profile)
        if (q.z > st.z0 && q.z < st.z1) cuts.push_back(q.z);
      std::sort(cuts.begin(), cuts.end());
      for (std::size_t k = 1; k < cuts.size(); ++k) {
        const double a = cuts[k - 1], b = cuts[k], m = 0.5 * (a + b), hh = 0.5 * (b - a);
        for (int g = 0; g < 3; ++g) {
          const double z = m + hh * gx[g];
          const double w = st.b0 + (st.b1 - st.b0) * (z - st.z0) / (st.z1 - st.z0);
          s += gw[g] * hh * st.E * w * f(z);
        }
      }
    }
    return s;
  };
  EquivalentTemperature e;
  double ztop = strips.front().z0, zbot = strips.front().z1;
  for (const ThermalStrip& st : strips) {
    ztop = std::min(ztop, st.z0);
    zbot = std::max(zbot, st.z1);
  }
  e.h = zbot - ztop;
  e.EA = integrate([](double) { return 1.0; });
  e.zbar = integrate([](double z) { return z; }) / e.EA;
  e.EI = integrate([&](double z) { return (z - e.zbar) * (z - e.zbar); });
  e.dT_eq = integrate([&](double z) { return profile_temperature(profile, z); }) / e.EA;
  e.dTz_eq = e.h * integrate([&](double z) { return profile_temperature(profile, z) * (z - e.zbar); }) / e.EI;
  return e;
}

double thermal_eigenstress(const EquivalentTemperature& eq, const std::vector<TempPoint>& profile,
                           double z, double E, double alpha_t) {
  return -E * alpha_t * (profile_temperature(profile, z) - eq.dT_eq - eq.dTz_eq * (z - eq.zbar) / eq.h);
}

}  // namespace keta

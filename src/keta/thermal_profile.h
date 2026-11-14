// Non-linear temperature profile through the depth of a cross-section
// (EN 1991-1-5 §6.1.4.2 / Annex B profiles; Keogh & O'Brien): its
// decomposition into the equivalent uniform part, the equivalent linear part
// and the self-equilibrating eigenstresses, on a section made of horizontal
// trapezoidal strips of possibly different moduli (composite sections are
// weighted by E):
//   dT_eq   = INT E b dT dz / INT E b dz
//   dTz,eq  = h INT E b dT (z - zbar) dz / INT E b (z - zbar)^2 dz
//   sigma_e(z) = -E alpha [dT(z) - dT_eq - dTz,eq (z - zbar) / h]
// z measured downward from the top fibre, h the total depth, zbar the
// E-weighted centroid. dT is piecewise linear between the profile points
// and b linear in each strip, so every integrand is a polynomial of degree
// <= 3 between the breakpoints: 3-point Gauss there is exact.
#pragma once
#include <vector>

namespace keta {

struct ThermalStrip {
  double z0 = 0, z1 = 0;  // depth range [m], z down from the top fibre
  double b0 = 0, b1 = 0;  // widths at z0 and z1 [m] (linear between)
  double E = 1;           // modulus (only ratios matter)
};
struct TempPoint {
  double z = 0, dT = 0;  // depth [m], temperature change [K]
};
struct EquivalentTemperature {
  double dT_eq = 0, dTz_eq = 0;  // [K]
  double zbar = 0, h = 0;        // E-weighted centroid depth, total depth [m]
  double EA = 0, EI = 0;         // E-weighted area and second moment (units of E)
};
EquivalentTemperature equivalent_temperature(const std::vector<ThermalStrip>& strips,
                                             const std::vector<TempPoint>& profile);
// Self-equilibrating eigenstress at depth z (strip modulus E there, alpha_t).
double thermal_eigenstress(const EquivalentTemperature& eq, const std::vector<TempPoint>& profile,
                           double z, double E, double alpha_t);
// The profile value at z (linear interpolation, constant beyond the ends).
double profile_temperature(const std::vector<TempPoint>& profile, double z);

}  // namespace keta

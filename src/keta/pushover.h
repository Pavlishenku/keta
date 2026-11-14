// Performance point of an equivalent SDOF capacity curve (ADRS format: Sd
// [m], Sa [m/s2]):
//   * N2 method, EN 1998-1 Annex B (Fajfar): the idealised elasto-plastic
//     system of yield point (Sdy, Sy), T* = 2 pi sqrt(Sdy / Sy), Ry = Se(T*) /
//     Sy, mu = Ry (T* >= TC) or (Ry - 1) TC / T* + 1 (T* < TC), Sdp = mu Sde /
//     Ry, Sp = Se / Ry (elastic when Ry <= 1: PP = PE);
//   * ATC-40 capacity-spectrum method, procedure A (ATC-40 1996 §8.2.2.1):
//     at a trial point (dpi, api) of the multilinear capacity curve, the
//     bilinear idealisation with the initial slope and EQUAL AREA up to dpi
//     gives (dy, ay); beta0 = 63.7 (ay dpi - dy api) / (api dpi) [%], kappa
//     of the structural behaviour type (ATC-40 Table 8-1: A 1.0 up to beta0 =
//     16.25 then 1.13 - 0.51 x; B 0.67 up to 25 then 0.845 - 0.446 x; C 0.33;
//     x = (ay dpi - dy api) / (api dpi)), beta_eff = kappa beta0 + 5;
//     SRA = (3.21 - 0.68 ln beta_eff) / 2.12, SRV = (2.31 - 0.41 ln
//     beta_eff) / 1.65, not below Table 8-2 (A 0.33 / 0.50, B 0.44 / 0.56,
//     C 0.56 / 0.67); the UBC demand (CA, CV) reduced by them. The
//     performance point is the ROOT (bisection to round-off) of capacity minus
//     reduced demand along the capacity curve — the converged procedure, not
//     the 5 % graphical tolerance.
#pragma once
#include <functional>
#include <vector>

namespace keta {

struct N2Result {
  double T = 0, Se = 0, Sde = 0, Ry = 0, mu = 0, Sdp = 0, Sp = 0;
};
// Se(T) the elastic acceleration spectrum [m/s2], TC its corner period.
N2Result n2_performance_point(double Sdy, double Sy, const std::function<double(double)>& Se, double TC);

struct CapacityPoint {
  double Sd = 0, Sa = 0;
};
struct Atc40Result {
  double beta_eff = 0, SRA = 0, SRV = 0;  // beta_eff in %
  double Sdy = 0, Say = 0, Sdp = 0, Sap = 0;
};
// UBC 1997 5 % spectrum in g: CA, CV; capacity: points after the origin,
// increasing Sd. g the gravity used for the ADRS conversion.
enum class Atc40Behaviour { A, B, C };
Atc40Result atc40_performance_point(const std::vector<CapacityPoint>& capacity, double CA, double CV, double g,
                                    Atc40Behaviour type = Atc40Behaviour::A);

}  // namespace keta

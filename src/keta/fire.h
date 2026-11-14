// Structural fire engineering:
//  * exposure curves (ISO 834 standard, external, hydrocarbon HC, French HCM,
//    user table) with their convection coefficients;
//  * EN 1993-1-2 steel member heating by the §4.2.5 incremental method
//    (unprotected with shadow factor, protected §4.2.5.2), net heat flux =
//    convection + radiation (eps_m = 0.7, sigma = 5.67e-8, +273 offset);
//  * EN 1993-1-2 Table 3.1 reduction factors ky/kp/kE(theta) and the critical
//    temperature 39.19 ln(1/(0.9674 mu0^3.833) - 1) + 482;
//  * EN 1992-1-2 Table 3.1 concrete kc(theta) for SILICEOUS and CALCAREOUS
//    aggregates with eps_c1(theta), eps_cu1(theta);
//  * EN 1992-1-2 Tables 3.2a/3.2b reinforcing-steel reduction factors
//    ks/kp/kE(theta) for class N (hot rolled / cold worked) and class X;
//  * EN 1992-1-2 §3.3.1 free thermal elongation of concrete and reinforcement;
//  * hot fiber laws: steel four-branch elliptic, concrete fire parabola with
//    linear descending branch — for M-kappa capacity at temperature.
#pragma once
#include "keta/fiber.h"

namespace keta {

// EN 1992-1-2 aggregate family (ST1 `GRANULAT SILICE / CALCAIRE`, manual p214).
// Drives kc(theta) (Table 3.1) and the free thermal elongation (§3.3.1).
enum class Aggregate { Siliceous, Calcareous };

// EN 1992-1-2 §3.2.3 reinforcing-steel class (ST1 `ACIER ... CLASSE <N, X>`)
// and forming route (ST1 `LAMINE <CHAUD, FROID>`; ST1 default is cold worked,
// manual p223).
enum class SteelClass { N, X };
enum class SteelForming { HotRolled, ColdWorked };

struct ExposureCurve {
  enum class Kind { Iso834, External, Hydrocarbon, Hcm, Table };
  Kind kind = Kind::Iso834;
  std::vector<std::pair<double, double>> table;  // (minutes, deg C)
  double alpha_c_user = 25.0;

  double gas_temp(double t_minutes) const;  // deg C (starts at 20)
  double alpha_c() const;                   // convection coefficient W/m^2K
};

// Net heat flux [W/m^2] on a surface at theta_m from gas at theta_g (deg C).
double net_heat_flux(double theta_g, double theta_m, double alpha_c, double eps_m = 0.7);

// EN 1993-1-2 §4.2.5.1 unprotected steel heating: section factor Am/V [1/m],
// shadow factor ksh; returns member temperature history at dt steps [s].
std::vector<double> steel_temperature_unprotected(const ExposureCurve& fire, double am_v,
                                                  double ksh, double t_end_minutes,
                                                  double dt_seconds = 5.0);
// §4.2.5.2 protected: dp thickness [m], lambda_p [W/mK], rho_p [kg/m^3], cp [J/kgK].
std::vector<double> steel_temperature_protected(const ExposureCurve& fire, double ap_v,
                                                double dp, double lambda_p, double rho_p,
                                                double cp_p, double t_end_minutes,
                                                double dt_seconds = 5.0);

// EN 1993-1-2 Table 3.1 (carbon steel): interpolated.
double steel_ky(double theta_c);  // effective yield
double steel_kp(double theta_c);  // proportional limit
double steel_ke(double theta_c);  // Young's modulus
// EN 1993-1-2 §4.2.4 critical temperature from the utilization mu0.
double steel_critical_temperature(double mu0);

// EN 1992-1-2 Table 3.1 concrete. `agg` selects the siliceous or calcareous
// kc(theta) column; eps_c1 / eps_cu1 are shared by both columns in Table 3.1.
double concrete_kc(double theta_c, Aggregate agg = Aggregate::Siliceous);
double concrete_eps_c1(double theta_c);
double concrete_eps_cu1(double theta_c);

// EN 1992-1-2 Tables 3.2a (class N) / 3.2b (class X): reinforcing steel.
//   ks = f_sy,theta / f_yk   (effective yield)
//   kp = f_sp,theta / f_yk   (proportional limit)
//   kE = E_s,theta / E_s     (slope of the linear branch)
// Table 3.2a splits hot rolled / cold worked; Table 3.2b (class X) is a single
// column, so `forming` is ignored for class X.
double rebar_ks(double theta_c, SteelClass cls = SteelClass::N,
                SteelForming forming = SteelForming::ColdWorked);
double rebar_kp(double theta_c, SteelClass cls = SteelClass::N,
                SteelForming forming = SteelForming::ColdWorked);
double rebar_kE(double theta_c, SteelClass cls = SteelClass::N,
                SteelForming forming = SteelForming::ColdWorked);

// EN 1992-1-2 §3.3.1 free thermal elongation (dl/l, dimensionless), theta in C.
//   concrete siliceous : -1.8e-4 + 9e-6 th + 2.3e-11 th^3   (20..700 C), 14e-3 above
//   concrete calcareous: -1.2e-4 + 6e-6 th + 1.4e-11 th^3   (20..805 C), 12e-3 above
//   reinforcement      : -2.416e-4 + 1.2e-5 th + 0.4e-8 th^2 (20..750 C),
//                        11e-3 (750..860 C), -6.2e-3 + 2e-5 th (860..1200 C)
double concrete_thermal_strain(double theta_c, Aggregate agg = Aggregate::Siliceous);
double rebar_thermal_strain(double theta_c);

// Hot laws at a uniform member temperature.
UniaxialLaw hot_steel_law(double E20, double fy20, double theta_c);
UniaxialLaw hot_concrete_law(double fcd20, double theta_c,
                             Aggregate agg = Aggregate::Siliceous);
// Reinforcing bar hot law (EN 1992-1-2 §3.2.3, same four-branch elliptic shape
// as structural steel but with the Table 3.2a/3.2b factors).
UniaxialLaw hot_rebar_law(double E20, double fyk, double theta_c,
                          SteelClass cls = SteelClass::N,
                          SteelForming forming = SteelForming::ColdWorked);

// Fire resistance of a steel fiber section at uniform temperature: plastic
// moment ratio M_fi / M_20 (M-kappa capacity at large curvature).
double fire_moment_ratio(const FiberSection& cold, double E20, double fy20, double theta_c);

}  // namespace keta

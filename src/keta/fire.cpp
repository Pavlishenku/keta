#include "keta/fire.h"

namespace keta {

double ExposureCurve::gas_temp(double t) const {
  switch (kind) {
    case Kind::Iso834: return 20.0 + 345.0 * std::log10(8.0 * t + 1.0);
    case Kind::External:
      return 20.0 + 660.0 * (1.0 - 0.687 * std::exp(-0.32 * t) - 0.313 * std::exp(-3.8 * t));
    case Kind::Hydrocarbon:
      return 20.0 + 1080.0 * (1.0 - 0.325 * std::exp(-0.167 * t) - 0.675 * std::exp(-2.5 * t));
    case Kind::Hcm:
      return 20.0 + 1280.0 * (1.0 - 0.325 * std::exp(-0.167 * t) - 0.675 * std::exp(-2.5 * t));
    case Kind::Table: {
      require(!table.empty(), "fire: empty gas-temperature table");
      if (t <= table.front().first) return table.front().second;
      if (t >= table.back().first) return table.back().second;
      for (std::size_t i = 1; i < table.size(); ++i)
        if (t <= table[i].first) {
          double f = (t - table[i - 1].first) / (table[i].first - table[i - 1].first);
          return table[i - 1].second + f * (table[i].second - table[i - 1].second);
        }
      return table.back().second;
    }
  }
  return 20.0;
}

double ExposureCurve::alpha_c() const {
  switch (kind) {
    case Kind::Hydrocarbon:
    case Kind::Hcm: return 50.0;
    case Kind::Table: return alpha_c_user;
    default: return 25.0;
  }
}

double net_heat_flux(double tg, double tm, double ac, double eps_m) {
  constexpr double sigma_sb = 5.67e-8;
  double Tg = tg + 273.0, Tm = tm + 273.0;  // the standards' literal offset
  return ac * (tg - tm) + sigma_sb * eps_m * 1.0 * (Tg * Tg * Tg * Tg - Tm * Tm * Tm * Tm);
}

namespace {
// EN 1993-1-2 §3.4.1.2 steel specific heat [J/kgK], theta in deg C.
double steel_cp(double th) {
  if (th < 600)
    return 425.0 + 0.773 * th - 1.69e-3 * th * th + 2.22e-6 * th * th * th;
  if (th < 735) return 666.0 + 13002.0 / (738.0 - th);
  if (th < 900) return 545.0 + 17820.0 / (th - 731.0);
  return 650.0;
}
constexpr double kSteelRho = 7850.0;
}  // namespace

std::vector<double> steel_temperature_unprotected(const ExposureCurve& fire, double am_v,
                                                  double ksh, double t_end_min, double dt) {
  require(am_v > 0 && dt > 0, "fire: bad section factor or step");
  std::vector<double> out;
  double th = 20.0;
  double t = 0, t_end = t_end_min * 60.0;
  out.push_back(th);
  while (t < t_end) {
    double tg = fire.gas_temp((t + dt) / 60.0);
    double h = net_heat_flux(tg, th, fire.alpha_c());
    th += ksh * am_v / (steel_cp(th) * kSteelRho) * h * dt;
    t += dt;
    out.push_back(th);
  }
  return out;
}

std::vector<double> steel_temperature_protected(const ExposureCurve& fire, double ap_v, double dp,
                                                double lp, double rp, double cp, double t_end_min,
                                                double dt) {
  require(ap_v > 0 && dp > 0 && lp > 0, "fire: bad protection data");
  std::vector<double> out;
  double th = 20.0;
  double t = 0, t_end = t_end_min * 60.0;
  out.push_back(th);
  double tg_prev = fire.gas_temp(0.0);
  while (t < t_end) {
    double tg = fire.gas_temp((t + dt) / 60.0);
    double phi = (cp * rp / (steel_cp(th) * kSteelRho)) * dp * ap_v;
    double dth = (lp * ap_v / (dp * steel_cp(th) * kSteelRho)) * (tg - th) / (1.0 + phi / 3.0) * dt -
                 (std::exp(phi / 10.0) - 1.0) * (tg - tg_prev);
    if (dth < 0 && tg > th) dth = 0;  // EN rule: no cooling while the gas heats
    th += dth;
    tg_prev = tg;
    t += dt;
    out.push_back(th);
  }
  return out;
}

namespace {
struct KRow {
  double th, ky, kp, ke;
};
// EN 1993-1-2 Table 3.1.
const KRow kSteel[] = {
    {20, 1.000, 1.000, 1.000},  {100, 1.000, 1.000, 1.000}, {200, 1.000, 0.807, 0.900},
    {300, 1.000, 0.613, 0.800}, {400, 1.000, 0.420, 0.700}, {500, 0.780, 0.360, 0.600},
    {600, 0.470, 0.180, 0.310}, {700, 0.230, 0.075, 0.130}, {800, 0.110, 0.050, 0.090},
    {900, 0.060, 0.0375, 0.0675}, {1000, 0.040, 0.0250, 0.0450}, {1100, 0.020, 0.0125, 0.0225},
    {1200, 0.000, 0.000, 0.000}};

double interp_steel(double th, double KRow::*field) {
  if (th <= kSteel[0].th) return kSteel[0].*field;
  for (std::size_t i = 1; i < sizeof(kSteel) / sizeof(KRow); ++i)
    if (th <= kSteel[i].th) {
      double f = (th - kSteel[i - 1].th) / (kSteel[i].th - kSteel[i - 1].th);
      return kSteel[i - 1].*field + f * (kSteel[i].*field - kSteel[i - 1].*field);
    }
  return 0.0;
}

struct CRow {
  double th, kc, kc_cal, ec1, ecu1;
};
// EN 1992-1-2 Table 3.1: kc(theta) for siliceous and calcareous aggregates
// (the eps_c1 / eps_cu1 columns are common to both aggregate families).
const CRow kConc[] = {{20, 1.00, 1.00, 0.0025, 0.0200},   {100, 1.00, 1.00, 0.0040, 0.0225},
                      {200, 0.95, 0.97, 0.0055, 0.0250},  {300, 0.85, 0.91, 0.0070, 0.0275},
                      {400, 0.75, 0.85, 0.0100, 0.0300},  {500, 0.60, 0.74, 0.0150, 0.0325},
                      {600, 0.45, 0.60, 0.0250, 0.0350},  {700, 0.30, 0.43, 0.0250, 0.0375},
                      {800, 0.15, 0.27, 0.0250, 0.0400},  {900, 0.08, 0.15, 0.0250, 0.0425},
                      {1000, 0.04, 0.06, 0.0250, 0.0450}, {1100, 0.01, 0.02, 0.0250, 0.0475},
                      {1200, 0.00, 0.00, 0.0250, 0.0500}};

double interp_conc(double th, double CRow::*field) {
  if (th <= kConc[0].th) return kConc[0].*field;
  for (std::size_t i = 1; i < sizeof(kConc) / sizeof(CRow); ++i)
    if (th <= kConc[i].th) {
      double f = (th - kConc[i - 1].th) / (kConc[i].th - kConc[i - 1].th);
      return kConc[i - 1].*field + f * (kConc[i].*field - kConc[i - 1].*field);
    }
  return kConc[sizeof(kConc) / sizeof(CRow) - 1].*field;
}

// EN 1992-1-2 Tables 3.2a / 3.2b — reinforcing steel reduction factors.
// Class N hot rolled reproduces EN 1993-1-2 Table 3.1 column for column
// (rebar_ks(th, N, HotRolled) == steel_ky(th)).
struct RRow {
  double th;
  double ks_nh, kp_nh, ke_nh;  // class N, hot rolled
  double ks_nc, kp_nc, ke_nc;  // class N, cold worked
  double ks_x, kp_x, ke_x;     // class X (single column, Table 3.2b)
};
const RRow kRebar[] = {
    {20, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.00, 1.00, 1.00},
    {100, 1.000, 1.000, 1.000, 1.000, 0.960, 1.000, 1.00, 1.00, 1.00},
    {200, 1.000, 0.810, 0.900, 1.000, 0.920, 0.870, 1.00, 0.87, 0.95},
    {300, 1.000, 0.610, 0.800, 1.000, 0.810, 0.720, 1.00, 0.74, 0.90},
    {400, 1.000, 0.420, 0.700, 0.940, 0.630, 0.560, 0.90, 0.70, 0.75},
    {500, 0.780, 0.360, 0.600, 0.670, 0.440, 0.400, 0.70, 0.51, 0.60},
    {600, 0.470, 0.180, 0.310, 0.400, 0.260, 0.240, 0.47, 0.18, 0.31},
    {700, 0.230, 0.075, 0.130, 0.120, 0.080, 0.080, 0.23, 0.07, 0.13},
    {800, 0.110, 0.050, 0.090, 0.110, 0.060, 0.060, 0.11, 0.05, 0.09},
    {900, 0.060, 0.0375, 0.0675, 0.080, 0.050, 0.050, 0.06, 0.04, 0.07},
    {1000, 0.040, 0.0250, 0.0450, 0.050, 0.030, 0.030, 0.04, 0.02, 0.04},
    {1100, 0.020, 0.0125, 0.0225, 0.030, 0.020, 0.020, 0.02, 0.01, 0.02},
    {1200, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.00, 0.00, 0.00}};

double interp_rebar(double th, double RRow::*field) {
  if (th <= kRebar[0].th) return kRebar[0].*field;
  for (std::size_t i = 1; i < sizeof(kRebar) / sizeof(RRow); ++i)
    if (th <= kRebar[i].th) {
      double f = (th - kRebar[i - 1].th) / (kRebar[i].th - kRebar[i - 1].th);
      return kRebar[i - 1].*field + f * (kRebar[i].*field - kRebar[i - 1].*field);
    }
  return 0.0;
}

double RRow::*rebar_field(SteelClass cls, SteelForming forming, int which) {
  if (cls == SteelClass::X)
    return which == 0 ? &RRow::ks_x : (which == 1 ? &RRow::kp_x : &RRow::ke_x);
  if (forming == SteelForming::HotRolled)
    return which == 0 ? &RRow::ks_nh : (which == 1 ? &RRow::kp_nh : &RRow::ke_nh);
  return which == 0 ? &RRow::ks_nc : (which == 1 ? &RRow::kp_nc : &RRow::ke_nc);
}
}  // namespace

double steel_ky(double th) { return interp_steel(th, &KRow::ky); }
double steel_kp(double th) { return interp_steel(th, &KRow::kp); }
double steel_ke(double th) { return interp_steel(th, &KRow::ke); }

double rebar_ks(double th, SteelClass c, SteelForming f) {
  return interp_rebar(th, rebar_field(c, f, 0));
}
double rebar_kp(double th, SteelClass c, SteelForming f) {
  return interp_rebar(th, rebar_field(c, f, 1));
}
double rebar_kE(double th, SteelClass c, SteelForming f) {
  return interp_rebar(th, rebar_field(c, f, 2));
}

double steel_critical_temperature(double mu0) {
  require(mu0 >= 0.013 && mu0 <= 1.0, "critical temperature: mu0 in [0.013, 1]");
  return 39.19 * std::log(1.0 / (0.9674 * std::pow(mu0, 3.833)) - 1.0) + 482.0;
}

double concrete_kc(double th, Aggregate agg) {
  return interp_conc(th, agg == Aggregate::Calcareous ? &CRow::kc_cal : &CRow::kc);
}
double concrete_eps_c1(double th) { return interp_conc(th, &CRow::ec1); }
double concrete_eps_cu1(double th) { return interp_conc(th, &CRow::ecu1); }

double concrete_thermal_strain(double th, Aggregate agg) {
  double t = std::min(std::max(th, 20.0), 1200.0);
  if (agg == Aggregate::Calcareous) {
    if (t <= 805.0) return -1.2e-4 + 6e-6 * t + 1.4e-11 * t * t * t;
    return 12e-3;
  }
  if (t <= 700.0) return -1.8e-4 + 9e-6 * t + 2.3e-11 * t * t * t;
  return 14e-3;
}

double rebar_thermal_strain(double th) {
  double t = std::min(std::max(th, 20.0), 1200.0);
  if (t <= 750.0) return -2.416e-4 + 1.2e-5 * t + 0.4e-8 * t * t;
  if (t <= 860.0) return 11e-3;
  return -6.2e-3 + 2e-5 * t;
}

UniaxialLaw hot_steel_law(double E20, double fy20, double th) {
  return UniaxialLaw::steel_elliptic(E20 * steel_ke(th), fy20 * steel_kp(th),
                                     fy20 * steel_ky(th));
}

UniaxialLaw hot_rebar_law(double E20, double fyk, double th, SteelClass cls,
                          SteelForming forming) {
  // At 20 C the proportional limit meets the yield plateau (kp = ks = 1): the
  // elliptic branch then degenerates continuously to the bilinear law, so the
  // only guard needed is fp <= fy.
  double fy = fyk * rebar_ks(th, cls, forming);
  double fp = std::min(fyk * rebar_kp(th, cls, forming), fy);
  return UniaxialLaw::steel_elliptic(E20 * rebar_kE(th, cls, forming), fp, fy);
}

UniaxialLaw hot_concrete_law(double fcd20, double th, Aggregate agg) {
  UniaxialLaw l = UniaxialLaw::concrete_pr(fcd20 * concrete_kc(th, agg), concrete_eps_c1(th),
                                           concrete_eps_cu1(th), 2.0);
  return l;
}

double fire_moment_ratio(const FiberSection& cold, double E20, double fy20, double theta_c) {
  require(cold.tendons.empty(), "fire_moment_ratio: prestressing tendons are not degraded by the fire laws");
  for (const Fiber& f : cold.fibers)
    require(f.eps0 == 0, "fire_moment_ratio: the initial strain of a fibre is not carried through the fire laws");
  FiberSection hot = cold;
  for (auto& law : hot.laws)
    if (law.kind == UniaxialLaw::Kind::SteelBilinear || law.kind == UniaxialLaw::Kind::Elastic)
      law = hot_steel_law(E20, fy20, theta_c);
  // capacity at a large curvature (fully plastic) at N = 0
  double h = 0;
  for (const Fiber& f : cold.fibers) h = std::max(h, std::fabs(f.z));
  double kappa_pl = 10.0 * 0.02 / std::max(h, 1e-6);  // ~10x yield strain at the extreme fiber
  auto mk_cold = cold.m_kappa(0.0, {kappa_pl});
  auto mk_hot = hot.m_kappa(0.0, {kappa_pl});
  require(std::fabs(mk_cold[0].moment) > 0, "fire: zero cold capacity");
  return mk_hot[0].moment / mk_cold[0].moment;
}

}  // namespace keta

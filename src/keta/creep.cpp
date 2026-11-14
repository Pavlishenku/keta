#include "keta/creep.h"

#include <cmath>

namespace keta {

double ec2_t0_adjusted(double t0, double cement_s) {
  // alpha exponent from cement class: s = 0.38 -> -1 (S), 0.25 -> 0 (N), 0.20 -> +1 (R)
  double alpha = cement_s <= 0.21 ? 1.0 : (cement_s >= 0.3 ? -1.0 : 0.0);
  double t0a = t0 * std::pow(9.0 / (2.0 + std::pow(t0, 1.2)) + 1.0, alpha);
  return std::max(0.5, t0a);
}

double ec2_beta_h(const ConcreteTime& c) {
  double fcm = c.fck / 1e6 + 8.0;
  double h0mm = c.h0 * 1000.0;
  double bh;
  if (fcm <= 35.0) {
    bh = 1.5 * (1.0 + std::pow(0.012 * c.rh, 18)) * h0mm + 250.0;
    bh = std::min(bh, 1500.0);
  } else {
    double a3 = std::sqrt(35.0 / fcm);
    bh = 1.5 * (1.0 + std::pow(0.012 * c.rh, 18)) * h0mm + 250.0 * a3;
    bh = std::min(bh, 1500.0 * a3);
  }
  return bh;
}

double ec2_phi(double t, double t0, const ConcreteTime& c) {
  require(t >= t0 && t0 > 0, "phi(t,t0): need t >= t0 > 0");
  double fcm = c.fck / 1e6 + 8.0;
  double h0mm = c.h0 * 1000.0;
  double phi_rh;
  if (fcm <= 35.0) {
    phi_rh = 1.0 + (1.0 - c.rh / 100.0) / (0.1 * std::cbrt(h0mm));
  } else {
    double a1 = std::pow(35.0 / fcm, 0.7), a2 = std::pow(35.0 / fcm, 0.2);
    phi_rh = (1.0 + (1.0 - c.rh / 100.0) / (0.1 * std::cbrt(h0mm)) * a1) * a2;
  }
  double beta_fcm = 16.8 / std::sqrt(fcm);
  double t0a = ec2_t0_adjusted(t0, c.cement_s);
  double beta_t0 = 1.0 / (0.1 + std::pow(t0a, 0.2));
  double bh = ec2_beta_h(c);
  double beta_c = std::pow((t - t0) / (bh + (t - t0)), 0.3);
  return phi_rh * beta_fcm * beta_t0 * beta_c;
}

double ec2_shrinkage(double t, const ConcreteTime& c) {
  double fck = c.fck / 1e6, fcm = fck + 8.0;
  double h0mm = c.h0 * 1000.0;
  // drying
  double ads1, ads2;
  if (c.cement_s <= 0.21) {         // R
    ads1 = 6.0;
    ads2 = 0.11;
  } else if (c.cement_s >= 0.3) {   // S
    ads1 = 3.0;
    ads2 = 0.13;
  } else {                          // N
    ads1 = 4.0;
    ads2 = 0.12;
  }
  double beta_rh = 1.55 * (1.0 - std::pow(c.rh / 100.0, 3));
  double eps_cd0 = 0.85 * (220.0 + 110.0 * ads1) * std::exp(-ads2 * fcm / 10.0) * 1e-6 * beta_rh;
  double kh = h0mm <= 100 ? 1.0
              : h0mm <= 200 ? 1.0 + (0.85 - 1.0) * (h0mm - 100) / 100
              : h0mm <= 300 ? 0.85 + (0.75 - 0.85) * (h0mm - 200) / 100
              : h0mm <= 500 ? 0.75 + (0.70 - 0.75) * (h0mm - 300) / 200
                            : 0.70;
  double ts = c.t_shrink0;
  double beta_ds = t > ts ? (t - ts) / ((t - ts) + 0.04 * std::sqrt(h0mm * h0mm * h0mm)) : 0.0;
  double eps_cd = beta_ds * kh * eps_cd0;
  // autogenous
  double eps_ca_inf = 2.5 * std::max(0.0, fck - 10.0) * 1e-6;
  double beta_as = 1.0 - std::exp(-0.2 * std::sqrt(t));
  double eps_ca = c.autogenous ? beta_as * eps_ca_inf : 0.0;
  return -(eps_cd + eps_ca);
}

double ec2_ecm(double fck) {
  double fcm = fck / 1e6 + 8.0;
  return 22000.0 * std::pow(fcm / 10.0, 0.3) * 1e6;
}

double aaem_modulus(double E, double dphi, double chi) { return E / (1.0 + chi * dphi); }

double ec2_relaxation_ratio(double t_hours, double sigma_pi, double fpk, int relaxation_class,
                            double rho1000) {
  require(t_hours >= 0 && fpk > 0, "relaxation: bad input");
  double mu = std::clamp(sigma_pi / fpk, 0.0, 1.0);
  double k1, k3;
  switch (relaxation_class) {
    case 1: k1 = 5.39; k3 = 6.7; break;
    case 2: k1 = 0.66; k3 = 9.1; break;
    case 3: k1 = 1.98; k3 = 8.0; break;
    default: fail("relaxation: class must be 1, 2 or 3");
  }
  return k1 * rho1000 * std::exp(k3 * mu) * std::pow(t_hours / 1000.0, 0.75 * (1.0 - mu)) * 1e-5;
}

double bpel_relaxation(double t_days, double sigma, double fprg, double rho1000) {
  require(t_days >= 0 && fprg > 0, "relaxation: bad input");
  if (sigma <= 0 || t_days == 0) return 0.0;
  double mu = std::min(sigma / fprg, 1.0);
  double k1 = rho1000 <= 2.5 ? 6e-3 : 8e-3;
  double k2 = rho1000 < 2.5 ? 1.10 : (rho1000 < 8.0 ? 1.25 : 1.50);
  return k1 * rho1000 * std::pow(24.0 * t_days / 1000.0, 0.75 * (1.0 - mu)) *
         std::exp((10.0 * mu - 7.5) / k2) * sigma;
}

double ip2_relaxation(double t_hours, double sigma_pi, double rg, double rho1000, double rho3000) {
  require(t_hours >= 0 && rg > 0 && rho1000 >= 0 && rho3000 >= 0, "IP2 relaxation: bad input");
  if (sigma_pi <= 0 || t_hours == 0) return 0.0;
  const double f = std::max(0.0, (sigma_pi - 0.55 * rg) / (0.25 * rg));
  return std::max(2.4 * rho1000, rho3000 + 2.5) / 100.0 * f * sigma_pi * std::pow(t_hours / 1e5, 0.19);
}

double bpel_relaxation_step(double dt_days, double sigma, double cum, double fprg,
                            double rho1000) {
  require(dt_days >= 0, "relaxation: negative time step");
  double s0 = sigma + cum;
  if (s0 <= 0 || dt_days == 0) return 0.0;
  // equivalent time: the loss grows monotonically with t (mu < 1). It is a
  // power of t, L(t) = L(t_ref) (t / t_ref)^p with p = 3/4 (1 - mu):
  // inverted in closed form (exact), the bisection kept for mu -> 1
  double t_eq = 0;
  const double p = 0.75 * (1.0 - std::min(s0 / fprg, 1.0));
  if (cum > 0 && p > 1e-9) {
    const double tref = 1000.0 / 24.0, Lref = bpel_relaxation(tref, s0, fprg, rho1000);
    if (Lref > 0) {
      t_eq = tref * std::pow(cum / Lref, 1.0 / p);
      if (std::isfinite(t_eq)) return std::max(0.0, bpel_relaxation(t_eq + dt_days, s0, fprg, rho1000) - cum);
      t_eq = 0;
    }
  }
  if (cum > 0) {
    double hi = 1.0;
    while (bpel_relaxation(hi, s0, fprg, rho1000) < cum && hi < 1e12) hi *= 2;
    double lo = 0;
    for (int it = 0; it < 200 && hi - lo > 1e-12 * hi; ++it) {
      double mid = 0.5 * (lo + hi);
      (bpel_relaxation(mid, s0, fprg, rho1000) < cum ? lo : hi) = mid;
    }
    t_eq = 0.5 * (lo + hi);
  }
  return std::max(0.0, bpel_relaxation(t_eq + dt_days, s0, fprg, rho1000) - cum);
}

double ec2_relaxation_step(double dt_hours, double sigma, double cum, double fpk,
                           int relaxation_class, double rho1000) {
  require(dt_hours >= 0, "relaxation: negative time step");
  const double s0 = sigma + cum;  // the stress before the relaxation so far
  if (s0 <= 0 || dt_hours == 0 || fpk <= 0) return 0.0;
  auto loss = [&](double th) {
    return th > 0 ? ec2_relaxation_ratio(th, s0, fpk, relaxation_class, rho1000) * s0 : 0.0;
  };
  // equivalent time on the current stress (EN 1992-1-1 annex D): the loss
  // grows monotonically with t, as a power of it (p = 3/4 (1 - mu)):
  // inverted in closed form, the bisection kept for mu -> 1
  double t_eq = 0;
  const double p = 0.75 * (1.0 - std::clamp(s0 / fpk, 0.0, 1.0));
  if (cum > 0 && p > 1e-9) {
    const double Lref = loss(1000.0);
    if (Lref > 0) {
      t_eq = 1000.0 * std::pow(cum / Lref, 1.0 / p);
      if (std::isfinite(t_eq)) return std::max(0.0, loss(t_eq + dt_hours) - cum);
      t_eq = 0;
    }
  }
  if (cum > 0) {
    double hi = 1.0;
    while (loss(hi) < cum && hi < 1e12) hi *= 2;
    double lo = 0;
    for (int it = 0; it < 200 && hi - lo > 1e-12 * hi; ++it) {
      const double mid = 0.5 * (lo + hi);
      (loss(mid) < cum ? lo : hi) = mid;
    }
    t_eq = 0.5 * (lo + hi);
  }
  return std::max(0.0, loss(t_eq + dt_hours) - cum);
}

double ec2_deferred_loss(double eps_cs, double Ep, double d_sigma_pr, double Ecm, double phi,
                         double sigma_c, double Ap, double Ac, double Ic, double zcp) {
  double num = std::fabs(eps_cs) * Ep + 0.8 * d_sigma_pr + (Ep / Ecm) * phi * sigma_c;
  double den = 1.0 + (Ep * Ap / (Ecm * Ac)) * (1.0 + Ac * zcp * zcp / Ic) * (1.0 + 0.8 * phi);
  return num / den;
}

}  // namespace keta

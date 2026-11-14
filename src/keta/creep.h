// Time-dependent concrete material laws, EN 1992-1-1 Annex B (+ §3.1.4):
// creep coefficient phi(t, t0) with the cement-class age adjustment and the
// high-strength alpha corrections, drying + autogenous shrinkage, and the
// EC2 §3.3.2 prestressing-steel relaxation classes. Ages in DAYS, strengths
// in Pa (converted internally to MPa), RH in %, h0 in m.
#pragma once
#include "keta/model.h"

namespace keta {

// Creep coefficient phi(t, t0) — EN 1992-1-1 eq. B.1..B.8c.
double ec2_phi(double t, double t0, const ConcreteTime& c);
// beta_H (helper, exposed for direct use).
double ec2_beta_h(const ConcreteTime& c);
// Cement-class adjusted loading age t0_adj (B.9); cement_s: 0.20 R, 0.25 N, 0.38 S.
double ec2_t0_adjusted(double t0, double cement_s);

// Total shrinkage eps_cs(t) = eps_cd (drying, from t_shrink0) + eps_ca
// (autogenous). NEGATIVE = shortening.
double ec2_shrinkage(double t, const ConcreteTime& c);

// Secant modulus Ecm from fck [Pa]: 22000 (fcm/10)^0.3 MPa, fcm = fck + 8 MPa.
double ec2_ecm(double fck);

// Age-adjusted effective modulus for the interval creep increment (AAEM with
// aging coefficient chi, default 0.8): E_eff = E / (1 + chi * dphi).
double aaem_modulus(double E, double dphi, double chi = 0.8);

// Prestressing-steel relaxation loss ratio d_sigma_pr / sigma_pi at time t [h]
// (EC2 eq. 3.28-3.30). relaxation_class: 1 (ordinary), 2 (low relaxation),
// 3 (hot-rolled bars). mu = sigma_pi / fpk clamped to [0, 1].
double ec2_relaxation_ratio(double t_hours, double sigma_pi, double fpk, int relaxation_class,
                            double rho1000);

// BPEL 91 annex 2 §3 pure relaxation loss after t days of a tendon tensioned
// to `sigma` (same unit as sigma), as ST1 states it (manual annex B4.3.2):
//   k1 rho1000 (24 t / 1000)^(3/4 (1 - mu)) exp((10 mu - 7.5) / k2) sigma,
//   mu = sigma / fprg, k1 = 6e-3 (rho1000 <= 2.5) else 8e-3,
//   k2 = 1.10 (rho1000 < 2.5), 1.25 (< 8), 1.50.
double bpel_relaxation(double t_days, double sigma, double fprg, double rho1000);
// IP2 art. 9.2 pure relaxation:
// p(h) = p (h / 100000)^0.19, p the greater of 2.4 rho1000/100 and (rho3000 +
// 2.5)/100, times (sigma_pi - 0.55 Rg)/(0.25 Rg) sigma_pi (none below 0.55 Rg);
// h in HOURS, rho in %, same unit as sigma_pi.
double ip2_relaxation(double t_hours, double sigma_pi, double rg, double rho1000, double rho3000);
// Increment over `dt_days` of the cumulated relaxation `cum` of a tendon now
// at `sigma`, by the equivalent-time method (BPEL annex 2 §4): t* solves
// bpel_relaxation(t*, sigma + cum) = cum, and the cumulated loss becomes
// bpel_relaxation(t* + dt, sigma + cum). Telescopes over subdivisions of a
// constant-strain interval.
double bpel_relaxation_step(double dt_days, double sigma, double cum, double fprg,
                            double rho1000);
// EC2 relaxation loss over dt_hours on the current stress `sigma`, `cum`
// lost so far: the equivalent time on the curve of sigma + cum (EN 1992-1-1
// annex D), as bpel_relaxation_step.
double ec2_relaxation_step(double dt_hours, double sigma, double cum, double fpk,
                           int relaxation_class, double rho1000);

// Combined deferred loss (EN 1992-1-1 eq. 5.46 style):
// d_sigma = (eps_cs Ep + 0.8 d_sigma_pr + (Ep/Ecm) phi sigma_c) /
//           (1 + (Ep Ap / Ecm Ac)(1 + Ac zcp^2 / Ic)(1 + 0.8 phi))
// sigma_c compression-POSITIVE at tendon level.
double ec2_deferred_loss(double eps_cs, double Ep, double d_sigma_pr, double Ecm, double phi,
                         double sigma_c, double Ap, double Ac, double Ic, double zcp);

}  // namespace keta

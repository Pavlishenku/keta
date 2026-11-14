// EN 1992-1-1 ULS design formulas — see ec2_design.h.
#include "keta/ec2_design.h"

#include <algorithm>
#include <cmath>

namespace keta::ec2 {

namespace {
constexpr double MPa = 1e6;
}

double shear_k(double d) { return std::min(2.0, 1.0 + std::sqrt(0.2 / d)); }

namespace {
double vrd_c(const VrdcInput& in, double CRdc, double k1, double vmin) {
  require(in.bw > 0 && in.d > 0 && in.fck > 0, "ec2 shear: bw, d, fck must be positive");
  const double k = shear_k(in.d);
  const double rho = std::min(0.02, in.Asl / (in.bw * in.d));
  double scp = in.sigma_cp;
  if (in.fcd > 0) scp = std::min(scp, 0.2 * in.fcd);
  const double v = CRdc * k * std::cbrt(100 * rho * in.fck / MPa) * MPa + k1 * scp;
  return std::max(v, vmin + k1 * scp) * in.bw * in.d;
}
}  // namespace

double vrd_c_en(const VrdcInput& in) {
  const double k = shear_k(in.d);
  const double vmin = 0.035 * std::pow(k, 1.5) * std::sqrt(in.fck / MPa) * MPa;
  return vrd_c(in, 0.18 / in.gamma_c, 0.15, vmin);
}

double vrd_c_din(const VrdcInput& in) {
  require(in.d <= 0.6, "ec2 shear (DIN NA): v_min of eq. 6.3aDE is for d <= 600 mm");
  const double k = shear_k(in.d);
  const double vmin = 0.0525 / in.gamma_c * std::pow(k, 1.5) * std::sqrt(in.fck / MPa) * MPa;
  return vrd_c(in, 0.15 / in.gamma_c, 0.12, vmin);
}

double nu_en(double fck) { return 0.6 * (1 - fck / (250 * MPa)); }
double nu1_din(double fck) { return 0.75 * std::min(1.0, 1.1 - fck / (500 * MPa)); }

double z_din(double d, double c_vl) {
  return std::min(0.9 * d, std::max(d - c_vl - 0.03, d - 2 * c_vl));
}

double vrd_cc_din(double fck, double sigma_cd, double fcd, double bw, double z) {
  return 0.5 * 0.48 * std::cbrt(fck / MPa) * (1 - 1.2 * sigma_cd / fcd) * bw * z * MPa;
}

double cot_theta_din(double sigma_cd, double fcd, double VRdcc, double VEd, double cot_max) {
  const double c = (1.2 + 1.4 * sigma_cd / fcd) / (1 - VRdcc / VEd);
  if (!(VEd > VRdcc)) return cot_max;  // the denominator vanishes or turns negative: the upper bound
  return std::min(cot_max, std::max(1.0, c));
}

double asf_flange(double vEd, double hf, double cot_theta_f, double fyd) { return vEd * hf / (cot_theta_f * fyd); }
double torsion_shear_flow(double TEd, double Ak) { return TEd / (2 * Ak); }

double asw_s(double VEd, double z, double fywd, double cot_theta) { return VEd / (z * fywd * cot_theta); }

double vrd_max(double bw, double z, double nu1, double fcd, double cot_theta, double alpha_cw) {
  return alpha_cw * bw * z * nu1 * fcd / (cot_theta + 1.0 / cot_theta);
}

double asw_min_s_en(double fck, double fyk, double bw, double alpha) {
  return 0.08 * std::sqrt(fck / MPa) / (fyk / MPa) * bw * std::sin(alpha);
}

TorsionBox torsion_box(double b, double h, double tef) {
  require(tef > 0 && tef < std::min(b, h), "ec2 torsion: t_ef out of range");
  return {tef, (b - tef) * (h - tef), 2 * (b + h - 2 * tef)};
}

double asw_s_torsion(double TEd, const TorsionBox& t, double fyd, double cot_theta) {
  return TEd / cot_theta / (fyd * 2 * t.Ak);
}

double asl_torsion(double TEd, const TorsionBox& t, double fyd, double cot_theta) {
  return TEd * t.uk * cot_theta / (fyd * 2 * t.Ak);
}

double trd_max(double nu, double fcd, const TorsionBox& t, double cot_theta, double alpha_cw) {
  const double th = std::atan(1.0 / cot_theta);
  return 2 * nu * alpha_cw * fcd * t.Ak * t.tef * std::sin(th) * std::cos(th);
}


double steel_stress(const SteelDesign& s, double eps) {
  const double a = std::fabs(eps), eyd = s.fyd / s.Es;
  require(a <= s.limit() * (1 + 1e-12), "ec2: steel strain beyond the design limit");
  const double sig = a <= eyd ? s.Es * a
                     : a >= s.eps_ud ? s.ftd
                                     : s.fyd + (a - eyd) * (s.ftd - s.fyd) / (s.eps_ud - eyd);
  return eps < 0 ? -sig : sig;
}

namespace {
struct PR {
  double n, ec2, ecu2;
};
PR pr_params(double fck) {
  const double f = fck / MPa;
  if (f <= 50) return {2.0, 2.0e-3, 3.5e-3};
  const double r = (90 - f) / 100;
  return {1.4 + 23.4 * std::pow(r, 4), (2.0 + 0.085 * std::pow(f - 50, 0.53)) * 1e-3,
          (2.6 + 35 * std::pow(r, 4)) * 1e-3};
}
// I0 = int_0^e s(t) dt, I1 = int_0^e t s(t) dt, s = sigma / fcd
void pr_integrals(const PR& p, double e, double& I0, double& I1) {
  const double ep = std::min(e, p.ec2), u1 = 1 - ep / p.ec2, n = p.n;
  I0 = ep - p.ec2 / (n + 1) * (1 - std::pow(u1, n + 1));
  I1 = ep * ep / 2 - p.ec2 * p.ec2 * ((1 - std::pow(u1, n + 1)) / (n + 1) - (1 - std::pow(u1, n + 2)) / (n + 2));
  if (e > p.ec2) {
    I0 += e - p.ec2;
    I1 += (e * e - p.ec2 * p.ec2) / 2;
  }
}
}  // namespace

Block parabola_rectangle_block(double fck, double eps_c) {
  const PR p = pr_params(fck);
  require(eps_c > 0 && eps_c <= p.ecu2 * (1 + 1e-12), "ec2: top strain outside (0, eps_cu2]");
  double I0, I1;
  pr_integrals(p, eps_c, I0, I1);
  return {I0 / eps_c, 1 - I1 / (eps_c * I0)};
}

RectDesign design_rect_bending(double b, double d, double d2, double fck, double fcd,
                               const SteelDesign& st, double MEds, double NEd, double xi_lim) {
  const PR p = pr_params(fck);
  RectDesign r;
  r.mu = MEds / (b * d * d * fcd);
  const Block bu = parabola_rectangle_block(fck, p.ecu2);
  auto mu3 = [&](double xi) { return bu.alpha * xi * (1 - bu.ka * xi); };
  r.mu_lim = mu3(xi_lim);
  const double esu = st.limit();
  const double xi23 = p.ecu2 / (p.ecu2 + esu);
  double xi, ec;
  if (r.mu <= mu3(xi23)) {
    // domain 2: eps_s1 = the steel limit, bisection on eps_c (mu increases with it)
    double lo = 0, hi = p.ecu2;
    for (int it = 0; it < 200 && hi - lo > 1e-18; ++it) {
      const double m = 0.5 * (lo + hi), x = m / (m + esu);
      const Block bl = parabola_rectangle_block(fck, m);
      (bl.alpha * x * (1 - bl.ka * x) < r.mu ? lo : hi) = m;
    }
    ec = 0.5 * (lo + hi);
    xi = ec / (ec + esu);
    const Block bl = parabola_rectangle_block(fck, ec);
    r.zeta = 1 - bl.ka * xi;
    r.eps_s1 = esu;
  } else {
    ec = p.ecu2;
    const double mu_c = std::min(r.mu, r.mu_lim);
    xi = (1 - std::sqrt(1 - 4 * bu.ka * mu_c / bu.alpha)) / (2 * bu.ka);
    xi = std::min(xi, xi_lim);
    r.zeta = 1 - bu.ka * xi;
    r.eps_s1 = ec * (1 - xi) / xi;
  }
  r.xi = xi;
  r.eps_c = ec;
  r.sigma_s1 = steel_stress(st, r.eps_s1);
  if (r.mu <= r.mu_lim) {
    r.As1 = (MEds / (r.zeta * d) + NEd) / r.sigma_s1;
    return r;
  }
  const double Mlim = r.mu_lim * b * d * d * fcd, dM = MEds - Mlim;
  r.eps_s2 = -ec * (xi * d - d2) / (xi * d);
  r.sigma_s2 = steel_stress(st, r.eps_s2);
  r.As1 = (Mlim / (r.zeta * d) + dM / (d - d2) + NEd) / r.sigma_s1;
  r.As2 = dM / ((d - d2) * std::fabs(r.sigma_s2));
  return r;
}


double rho_interface(double vEdi, double c, double fctd, double mu, double sigma_n, double fyd,
                     double alpha, bool din) {
  const double k = din ? 1.2 : 1.0;
  return std::max(0.0, (vEdi - c * fctd - mu * sigma_n) / (fyd * (k * mu * std::sin(alpha) + std::cos(alpha))));
}

double vrdi_max(double nu, double fcd) { return 0.5 * nu * fcd; }


double kc_rect(double sigma_c, double h, double fct_eff) {
  const double hs = std::min(h, 1.0);
  const double k1 = sigma_c >= 0 ? 1.5 : 2 * hs / (3 * h);
  return std::min(1.0, std::max(0.0, 0.4 * (1 - sigma_c / (k1 * (h / hs) * fct_eff))));
}

namespace {
double lerp_h(double h, double a, double b) {
  if (h <= 0.3) return a;
  if (h >= 0.8) return b;
  return a + (b - a) * (h - 0.3) / 0.5;
}
}  // namespace
double k_selfeq_din(double h) { return lerp_h(h, 0.8, 0.5); }

double sigma_s_tab72de(double wk, double phi_star) {
  return std::sqrt(wk * 1e3 * 3.48e6 / (phi_star * 1e3)) * MPa;
}
double phi_star_min_reinf(double phi_s, double fct_eff) { return phi_s * 2.9 * MPa / fct_eff; }
double phi_star_load(double phi_s, double sigma_s, double As, double h, double d, double b, double fct_eff) {
  const double f = std::max(sigma_s * As / (4 * (h - d) * b * 2.9 * MPa), fct_eff / (2.9 * MPa));
  return phi_s / f;
}
double crack_width_din(double sigma_s, double As, double Ac_eff, double phi, double fct_eff, double kt,
                       double alpha_e, double Es) {
  const double rho = As / Ac_eff;
  const double de = std::max((sigma_s - kt * fct_eff / rho * (1 + alpha_e * rho)) / Es, 0.6 * sigma_s / Es);
  const double sr = std::min(phi / (3.6 * rho), sigma_s * phi / (3.6 * fct_eff));
  return sr * de;
}

double as_min(double kc, double k, double fct_eff, double Act, double sigma_s) {
  return kc * k * fct_eff * Act / sigma_s;
}


double fcd_fat(double k1, double beta_cc_t0, double fcd, double fck) {
  return k1 * beta_cc_t0 * fcd * (1 - fck / (250 * MPa));
}
double dsigma_rsk_bars(bool din) { return (din ? 175.0 : 162.5) * MPa; }
double dsigma_s_equ(double lambda_s, double sigma_max, double sigma_min) {
  return lambda_s * (sigma_max - sigma_min);
}


double imperfection_theta(double l, double m, double theta0) {
  const double ah = std::min(1.0, std::max(2.0 / 3.0, 2.0 / std::sqrt(l)));
  return theta0 * ah * std::sqrt(0.5 * (1 + 1 / m));
}
double imperfection_ei(double theta_i, double l0) { return theta_i * l0 / 2; }
double nominal_curvature(double Kr, double Kphi, double eps_yd, double d) {
  return Kr * Kphi * eps_yd / (0.45 * d);
}
double second_order_e2(double K1, double inv_r, double l0, double c) { return K1 * inv_r * l0 * l0 / c; }


double punching_u1_rect(double a, double b, double d) { return 2 * (a + b) + 4 * kPi * d; }

namespace {
double vrdc_punch(const PunchingInput& in, double CRdc, double vmin, double rho_cap) {
  const double k = shear_k(in.d);
  const double rho = std::min(rho_cap, std::sqrt(in.rho_ly * in.rho_lz));
  const double v = CRdc * k * std::cbrt(100 * rho * in.fck / MPa) * MPa + 0.1 * in.sigma_cp;
  return std::max(v, vmin + 0.1 * in.sigma_cp);
}
}  // namespace

double vrd_c_punching_en(const PunchingInput& in) {
  const double k = shear_k(in.d);
  return vrdc_punch(in, 0.18 / in.gamma_c, 0.035 * std::pow(k, 1.5) * std::sqrt(in.fck / MPa) * MPa, 0.02);
}

double vrd_c_punching_din(const PunchingInput& in) {
  const double k = shear_k(in.d);
  double c = 0.18 / in.gamma_c;
  if (in.u0 > 0 && in.u0 / in.d < 4) c *= 0.1 * in.u0 / in.d + 0.6;
  require(in.d <= 0.6 || in.d > 0.8, "ec2 punching (DIN NA): v_min interpolation 600 < d <= 800 mm not here");
  const double cv = in.d <= 0.6 ? 0.0525 : 0.0375;
  const double vmin = cv / in.gamma_c * std::pow(k, 1.5) * std::sqrt(in.fck / MPa) * MPa;
  const double cap = in.fyd > 0 ? std::min(0.02, 0.5 * in.fcd / in.fyd) : 0.02;
  return vrdc_punch(in, c, vmin, cap);
}

double beta_punching_din(PunchingPosition p) {
  switch (p) {
    case PunchingPosition::Interior: return 1.10;
    case PunchingPosition::Edge: return 1.40;
    case PunchingPosition::Corner: return 1.50;
    case PunchingPosition::WallEnd: return 1.35;
    case PunchingPosition::WallCorner: return 1.20;
  }
  return 1.5;
}

double ved_punching(double beta, double VEd, double u, double d) { return beta * VEd / (u * d); }
double vrd_max_punching_din(double vrdc) { return 1.4 * vrdc; }

double sandwich_panel_force(double n, double m, double z) { return n / 2 + m / z; }

MembraneSteel membrane_reinforcement(double Nx, double Ny, double Nxy) {
  MembraneSteel r;
  const double a = std::fabs(Nxy);
  if (Nx >= -a && Ny >= -a) {
    r.Zx = Nx + a;
    r.Zy = Ny + a;
    r.Fc = 2 * a;
  } else if (Nx < -a) {
    r.Zx = 0;
    r.Zy = Ny + Nxy * Nxy / std::fabs(Nx);
    r.Fc = std::fabs(Nx) + Nxy * Nxy / std::fabs(Nx);
  } else {
    r.Zy = 0;
    r.Zx = Nx + Nxy * Nxy / std::fabs(Ny);
    r.Fc = std::fabs(Ny) + Nxy * Nxy / std::fabs(Ny);
  }
  r.Zx = std::max(0.0, r.Zx);
  r.Zy = std::max(0.0, r.Zy);
  return r;
}
}  // namespace keta::ec2

namespace keta::sia262 {
double tau_cd(double fck, double gamma_c, double eta_t) { return 0.3 * eta_t * std::sqrt(fck / 1e6) / gamma_c * 1e6; }
double perimeter_u(double a, double b, double d) { return 2 * a + 2 * b + kPi * d; }
double area_inside_u(double a, double b, double d) { return a * b + (a + b) * d + kPi * d * d / 4; }
double ke(double e_u, double b_u) { return 1.0 / (1.0 + e_u / b_u); }
double psi_loa3(double rs, double d, double fsd, double Es, double msd, double mRd) {
  return 1.2 * rs / d * fsd / Es * std::pow(msd / mRd, 1.5);
}
double kr(double psi, double d, double Dmax) {
  const double kg = 48.0 / (16.0 + Dmax * 1e3);
  return 1.0 / (0.45 + 0.18 * psi * d * 1e3 * kg);
}
double vrd_c(double kr, double tau_cd, double d, double u_red) { return kr * tau_cd * d * u_red; }
}  // namespace keta::sia262

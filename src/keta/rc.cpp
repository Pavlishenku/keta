#include "keta/rc.h"

#include <algorithm>
#include <cmath>

namespace keta {

// ---------------- geometry ----------------

void RcSection::validate() const {
  require(finite(b) && b > 0 && finite(h) && h > 0, "rc: b and h must be positive");
  require(finite(d) && d > 0 && d < h, "rc: d must lie inside the section");
  require(finite(As) && As >= 0 && finite(As2) && As2 >= 0, "rc: steel areas must be >= 0");
  require(As2 == 0 || (d2 > 0 && d2 < h), "rc: d2 must lie inside the section");
  require(finite(bf) && finite(hf) && hf >= 0 && bf >= 0, "rc: non-finite flange");
  if (hf > 0 && bf > b) require(hf < h, "rc: the flange must be thinner than the section");
}

double RcSection::width(double z) const {
  if (hf > 0 && bf > b && z < hf) return bf;
  return b;
}

double RcSection::conc_area(double t) const {
  double tt = std::clamp(t, 0.0, h);
  if (hf > 0 && bf > b) {
    double a = bf * std::min(tt, hf);
    if (tt > hf) a += b * (tt - hf);
    return a;
  }
  return b * tt;
}

double RcSection::conc_first_moment(double t) const {
  double tt = std::clamp(t, 0.0, h);
  auto strip = [](double w, double z1, double z2) { return w * 0.5 * (z2 * z2 - z1 * z1); };
  if (hf > 0 && bf > b) {
    double s = strip(bf, 0.0, std::min(tt, hf));
    if (tt > hf) s += strip(b, hf, tt);
    return s;
  }
  return strip(b, 0.0, tt);
}

double RcSection::conc_second_moment(double t) const {
  double tt = std::clamp(t, 0.0, h);
  auto strip = [](double w, double z1, double z2) {
    return w * (z2 * z2 * z2 - z1 * z1 * z1) / 3.0;
  };
  if (hf > 0 && bf > b) {
    double s = strip(bf, 0.0, std::min(tt, hf));
    if (tt > hf) s += strip(b, hf, tt);
    return s;
  }
  return strip(b, 0.0, tt);
}

double RcSection::centroid() const {
  double a = area();
  require(a > 0, "rc: zero concrete area");
  return conc_first_moment(h) / a;
}

// ---------------- material helpers ----------------

double ec2_fctm(double fck) {
  require(finite(fck) && fck > 0, "rc: fck must be positive");
  double f = fck / 1e6;
  if (f <= 50.0) return 0.30 * std::pow(f, 2.0 / 3.0) * 1e6;
  return 2.12 * std::log(1.0 + (f + 8.0) / 10.0) * 1e6;
}

double ec2_eps_cu3(double fck) {
  double f = fck / 1e6;
  if (f <= 50.0) return 0.0035;
  return (2.6 + 35.0 * std::pow((90.0 - f) / 100.0, 4)) * 1e-3;
}

double ec2_block_lambda(double fck) {
  double f = fck / 1e6;
  return f <= 50.0 ? 0.8 : 0.8 - (f - 50.0) / 400.0;
}

double ec2_block_eta(double fck) {
  double f = fck / 1e6;
  return f <= 50.0 ? 1.0 : 1.0 - (f - 50.0) / 200.0;
}

double rc_fctm(const RcConcrete& c) { return c.fctm > 0 ? c.fctm : ec2_fctm(c.fck); }
double rc_ecm(const RcConcrete& c) { return c.Ecm > 0 ? c.Ecm : ec2_ecm(c.fck); }

double rc_ec_eff(const RcConcrete& c, double phi) {
  require(finite(phi) && phi >= 0, "rc: creep coefficient must be >= 0");
  return rc_ecm(c) / (1.0 + phi);
}

double rc_alpha_e(const RcMaterials& m, double phi) {
  double ec = rc_ec_eff(m.c, phi);
  require(ec > 0 && m.s.Es > 0, "rc: non-positive modulus");
  return m.s.Es / ec;
}

// ---------------- state I ----------------

RcStateI rc_uncracked(const RcSection& s, const RcMaterials& m, double N, double M, double phi) {
  s.validate();
  require(finite(N) && finite(M), "rc: non-finite action");
  RcStateI r;
  r.alpha_e = rc_alpha_e(m, phi);
  double n = r.alpha_e - 1.0;  // steel replaces concrete
  double Ac = s.area(), Sc = s.conc_first_moment(s.h), Ic0 = s.conc_second_moment(s.h);
  r.A = Ac + n * (s.As + s.As2);
  double S = Sc + n * (s.As * s.d + s.As2 * s.d2);
  r.x = S / r.A;
  double I0 = Ic0 + n * (s.As * s.d * s.d + s.As2 * s.d2 * s.d2);  // about the top fibre
  r.I = I0 - r.A * r.x * r.x;                                      // parallel axis
  require(r.I > 0, "rc: degenerate uncracked section");

  double zg = s.centroid();
  double M_top = M + N * zg;             // internal moment about the top fibre
  double MI = M_top - N * r.x;           // about the transformed centroid
  auto sig = [&](double z) { return N / r.A + MI * (z - r.x) / r.I; };
  r.sigma_c_top = sig(0.0);
  r.sigma_c_bot = sig(s.h);
  r.sigma_s = r.alpha_e * sig(s.d);
  r.sigma_s2 = s.As2 > 0 ? r.alpha_e * sig(s.d2) : 0.0;
  r.kappa = MI / (rc_ec_eff(m.c, phi) * r.I);
  return r;
}

double rc_cracking_moment(const RcSection& s, const RcMaterials& m, double N, double phi,
                          double fct_eff) {
  RcStateI r = rc_uncracked(s, m, N, 0.0, phi);
  double fct = fct_eff > 0 ? fct_eff : rc_fctm(m.c);
  double zg = s.centroid();
  // sigma_c(h) = N/A + (M + N zg - N x)(h - x)/I = fct
  double lever = s.h - r.x;
  require(lever > 0, "rc: degenerate cracking-moment geometry");
  return (fct - N / r.A) * r.I / lever - N * (zg - r.x);
}

// ---------------- state II ----------------

namespace {

// Transformed cracked-section resultants for a UNIT stress gradient
// sigma_c(z) = -(x - z): axial force (tension +) and moment about the top
// fibre. `alpha_e` on the tension side, `alpha_e - 1` on the compression side.
void cracked_hat(const RcSection& s, double alpha_e, double x, double& nhat, double& mhat) {
  double t = std::clamp(x, 0.0, s.h);
  // concrete block: -INT_0^t w (x - z) dz  and  -INT_0^t w (x - z) z dz
  double A0 = s.conc_area(t), A1 = s.conc_first_moment(t), A2 = s.conc_second_moment(t);
  nhat = -(x * A0 - A1);
  mhat = -(x * A1 - A2);
  auto steel = [&](double As, double ds) {
    if (As <= 0) return;
    double n = (ds < x) ? (alpha_e - 1.0) : alpha_e;
    double f = -n * As * (x - ds);
    nhat += f;
    mhat += f * ds;
  };
  steel(s.As, s.d);
  steel(s.As2, s.d2);
}

}  // namespace

RcStateII rc_cracked(const RcSection& s, const RcMaterials& m, double N, double M, double phi) {
  s.validate();
  require(finite(N) && finite(M), "rc: non-finite action");
  RcStateII r;
  r.alpha_e = rc_alpha_e(m, phi);
  double zg = s.centroid();
  double M_top = M + N * zg;

  // f(x) = N * mhat(x) - M_top * nhat(x) = 0  (the resultant of the stress
  // block must act at the eccentricity of the applied action).
  auto f = [&](double x) {
    double nh, mh;
    cracked_hat(s, r.alpha_e, x, nh, mh);
    return N * mh - M_top * nh;
  };
  const int nscan = 800;
  double x_lo = 1e-6 * s.h, x_hi = 4.0 * s.h;
  double a = x_lo, fa = f(a), found_lo = 0, found_hi = 0;
  bool bracket = false;
  for (int i = 1; i <= nscan; ++i) {
    double x = x_lo + (x_hi - x_lo) * i / nscan;
    double fx = f(x);
    if ((fa <= 0 && fx >= 0) || (fa >= 0 && fx <= 0)) {
      found_lo = a;
      found_hi = x;
      bracket = true;
      break;
    }
    a = x;
    fa = fx;
  }
  require(bracket, "rc: no cracked neutral axis for this (N, M) — section fully in tension?");
  {
    double fl = f(found_lo);
    for (int i = 0; i < 120; ++i) {
      double mid = 0.5 * (found_lo + found_hi);
      double fm = f(mid);
      if ((fl <= 0) == (fm <= 0)) {
        found_lo = mid;
        fl = fm;
      } else {
        found_hi = mid;
      }
    }
  }
  r.x = 0.5 * (found_lo + found_hi);

  double nh, mh;
  cracked_hat(s, r.alpha_e, r.x, nh, mh);
  // Both equilibrium equations (N = lambda nhat, M_top = lambda mhat) hold at
  // the root; take the least-squares scale so neither degenerate case
  // (pure bending -> nhat = 0, pure compression -> mhat = 0) divides by zero.
  double den = nh * nh + mh * mh;
  require(den > 0, "rc: degenerate cracked equilibrium");
  double lambda = (N * nh + M_top * mh) / den;  // stress gradient [Pa/m]
  require(finite(lambda), "rc: degenerate cracked equilibrium");

  r.sigma_c_top = -lambda * r.x;
  r.sigma_s = lambda * r.alpha_e * (s.d - r.x);
  r.sigma_s2 = s.As2 > 0 ? -lambda * r.alpha_e * (r.x - s.d2) : 0.0;
  r.kappa = lambda / rc_ec_eff(m.c, phi);

  // cracked transformed inertia about the neutral axis
  double t = std::clamp(r.x, 0.0, s.h);
  double A0 = s.conc_area(t), A1 = s.conc_first_moment(t), A2 = s.conc_second_moment(t);
  r.I = A2 - 2 * r.x * A1 + r.x * r.x * A0;
  auto add = [&](double As, double ds) {
    if (As <= 0) return;
    double n = (ds < r.x) ? (r.alpha_e - 1.0) : r.alpha_e;
    r.I += n * As * (ds - r.x) * (ds - r.x);
  };
  add(s.As, s.d);
  add(s.As2, s.d2);
  r.lever_arm = (s.As > 0 && std::fabs(r.sigma_s) > 0) ? M / (s.As * r.sigma_s) : 0.0;
  return r;
}

// ---------------- crack width ----------------

RcCrackResult rc_crack_width(const RcSection& s, const RcMaterials& m, const RcCrackParams& p,
                             double N, double M, double phi) {
  require(p.cover >= 0 && p.bar_dia > 0, "rc: bad crack parameters");
  RcStateII st = rc_cracked(s, m, N, M, phi);
  RcCrackResult r;
  r.x = st.x;
  r.sigma_s = st.sigma_s;
  double fct = p.fct_eff > 0 ? p.fct_eff : rc_fctm(m.c);

  // §7.3.2(3): h_c,ef = min(2.5 (h - d), (h - x)/3, h/2)
  r.hc_ef = std::min({2.5 * (s.h - s.d), (s.h - st.x) / 3.0, s.h / 2.0});
  require(r.hc_ef > 0, "rc: degenerate effective tension height");
  r.Ac_eff = s.b * r.hc_ef;  // the web width carries the tension zone
  require(s.As > 0, "rc: crack width needs tension reinforcement");
  r.rho_p_eff = s.As / r.Ac_eff;  // eq. 7.10 (no prestress: xi1 = 0)

  // eq. 7.9, floored at 0.6 sigma_s / Es
  double alpha_e = st.alpha_e;
  double num = r.sigma_s - p.kt * fct * (1.0 + alpha_e * r.rho_p_eff) / r.rho_p_eff;
  r.d_eps = std::max(num / m.s.Es, 0.6 * r.sigma_s / m.s.Es);

  // eq. 7.11, or eq. 7.14 when the bars are farther apart than 5(c + phi/2)
  double s_lim = 5.0 * (p.cover + p.bar_dia / 2.0);
  r.wide_spacing = (p.spacing > s_lim);
  r.s_r_max = r.wide_spacing
                  ? 1.3 * (s.h - st.x)
                  : p.k3 * p.cover + p.k1 * p.k2 * p.k4 * p.bar_dia / r.rho_p_eff;
  r.wk = r.s_r_max * r.d_eps;  // eq. 7.8
  return r;
}

// ---------------- tension stiffening ----------------

RcCurvature rc_curvature(const RcSection& s, const RcMaterials& m, double N, double M, double phi,
                         double beta) {
  require(finite(beta) && beta > 0 && beta <= 1.0, "rc: beta must lie in (0, 1]");
  RcCurvature r;
  r.M_cr = rc_cracking_moment(s, m, N, phi);
  r.kappa_I = rc_uncracked(s, m, N, M, phi).kappa;
  if (std::fabs(M) <= std::fabs(r.M_cr) || s.As <= 0) {
    r.zeta = 0.0;
    r.kappa_II = 0.0;
    r.kappa = r.kappa_I;
    return r;
  }
  r.kappa_II = rc_cracked(s, m, N, M, phi).kappa;
  double ratio = r.M_cr / M;                    // = sigma_sr / sigma_s (eq. 7.19)
  r.zeta = 1.0 - beta * ratio * ratio;
  r.zeta = std::clamp(r.zeta, 0.0, 1.0);
  r.kappa = r.zeta * r.kappa_II + (1.0 - r.zeta) * r.kappa_I;  // eq. 7.18
  return r;
}

double rc_deflection_ss(const RcSection& s, const RcMaterials& m, double L,
                        const std::function<double(double)>& M_of_x, double N, double phi,
                        double beta, int n_per_half) {
  require(finite(L) && L > 0, "rc: span must be positive");
  require(n_per_half >= 2, "rc: need at least 2 intervals per half span");
  int n = n_per_half + (n_per_half % 2);  // Simpson needs an even interval count
  auto integrand = [&](double x) {
    double mx = (x <= L / 2) ? x / 2.0 : (L - x) / 2.0;  // unit midspan load
    return rc_curvature(s, m, N, M_of_x(x), phi, beta).kappa * mx;
  };
  auto simpson = [&](double x0, double x1) {
    double hh = (x1 - x0) / n, sum = integrand(x0) + integrand(x1);
    for (int i = 1; i < n; ++i)
      sum += (i % 2 ? 4.0 : 2.0) * integrand(x0 + hh * i);
    return sum * hh / 3.0;
  };
  return simpson(0.0, L / 2) + simpson(L / 2, L);
}

double rc_shrinkage_curvature(const RcSection& s, const RcMaterials& m, double eps_cs, double phi,
                              bool cracked) {
  const double ae = rc_alpha_e(m, phi);
  if (!cracked) {
    const RcStateI st = rc_uncracked(s, m, 0.0, 0.0, phi);
    double S = s.As * (s.d - st.x);
    if (s.As2 > 0) S -= s.As2 * (st.x - s.d2);
    return -eps_cs * ae * S / st.I;
  }
  // pure state II: neutral axis and inertia of the cracked section (they do
  // not depend on M for N = 0)
  const RcStateII st = rc_cracked(s, m, 0.0, 1.0, phi);
  double S = s.As * (s.d - st.x);
  if (s.As2 > 0) S -= s.As2 * (st.x - s.d2);
  return -eps_cs * ae * S / st.I;
}

double rc_deflection_ss_cs(const RcSection& s, const RcMaterials& m, double L,
                           const std::function<double(double)>& M_of_x, double phi, double eps_cs,
                           double beta, double phi_cr) {
  require(finite(L) && L > 0, "rc: span must be positive");
  const double Mcr = std::fabs(rc_cracking_moment(s, m, 0.0, phi_cr < 0 ? phi : phi_cr));
  const double kI_cs = rc_shrinkage_curvature(s, m, eps_cs, phi, false);
  const double kII_cs = s.As > 0 ? rc_shrinkage_curvature(s, m, eps_cs, phi, true) : 0.0;
  auto f = [&](double x) {
    const double mx = (x <= L / 2) ? x / 2.0 : (L - x) / 2.0;
    const double M = M_of_x(x);
    const double kI = rc_uncracked(s, m, 0.0, M, phi).kappa + kI_cs;
    if (std::fabs(M) <= Mcr || s.As <= 0) return kI * mx;
    const double zeta = std::clamp(1.0 - beta * (Mcr / M) * (Mcr / M), 0.0, 1.0);
    const double kII = rc_cracked(s, m, 0.0, M, phi).kappa + kII_cs;
    return (zeta * kII + (1 - zeta) * kI) * mx;
  };
  // cuts: mid-span (kink of m) and the cracking points
  std::vector<double> cuts = {0.0, L / 2, L};
  const int ns = 256;
  auto g = [&](double x) { return std::fabs(M_of_x(x)) - Mcr; };
  for (int i = 0; i < ns; ++i) {
    double a = L * i / ns, b = L * (i + 1) / ns;
    if ((g(a) > 0) == (g(b) > 0)) continue;
    for (int it = 0; it < 200 && b - a > 1e-15 * L; ++it) {
      const double c = 0.5 * (a + b);
      ((g(c) > 0) == (g(a) > 0) ? a : b) = c;
    }
    cuts.push_back(0.5 * (a + b));
  }
  std::sort(cuts.begin(), cuts.end());
  static const double gx[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                               0.9061798459386640};
  static const double gw[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                               0.4786286704993665, 0.2369268850561891};
  auto gauss = [&](double a, double b) {
    double s5 = 0;
    for (int k = 0; k < 5; ++k) s5 += gw[k] * f(0.5 * (a + b) + 0.5 * (b - a) * gx[k]);
    return 0.5 * (b - a) * s5;
  };
  std::function<double(double, double, double, int)> adapt = [&](double a, double b, double whole, int depth) {
    const double c = 0.5 * (a + b), l = gauss(a, c), r = gauss(c, b);
    if (depth >= 30 || std::fabs(l + r - whole) <= 1e-13 * std::max(std::fabs(l + r), 1e-300)) return l + r;
    return adapt(a, c, l, depth + 1) + adapt(c, b, r, depth + 1);
  };
  double total = 0;
  for (std::size_t i = 1; i < cuts.size(); ++i)
    if (cuts[i] > cuts[i - 1]) total += adapt(cuts[i - 1], cuts[i], gauss(cuts[i - 1], cuts[i]), 0);
  return total;
}

// ---------------- ULS bending ----------------

RcUls rc_mrd(const RcSection& s, const RcMaterials& m, double N) {
  s.validate();
  require(finite(N), "rc: non-finite axial force");
  RcUls r;
  r.fcd = m.c.alpha_cc * m.c.fck / m.c.gamma_c;
  r.fyd = m.s.fyk / m.s.gamma_s;
  r.lambda = ec2_block_lambda(m.c.fck);
  r.eta = ec2_block_eta(m.c.fck);
  r.eps_cu3 = ec2_eps_cu3(m.c.fck);
  require(r.fcd > 0 && r.fyd > 0, "rc: non-positive design strengths");

  // Forces, tension positive, as functions of the neutral axis depth x.
  auto sig_steel = [&](double eps) { return std::clamp(m.s.Es * eps, -r.fyd, r.fyd); };
  auto forces = [&](double x, double& Fc, double& Fs, double& Fs2, double& zc, double& eps_s,
                    double& ss, double& ss2) {
    double blk = std::min(r.lambda * x, s.h);
    double A = s.conc_area(blk), S = s.conc_first_moment(blk);
    Fc = -r.eta * r.fcd * A;
    zc = A > 0 ? S / A : 0.0;
    eps_s = x > 0 ? r.eps_cu3 * (s.d - x) / x : 0.0;
    ss = sig_steel(eps_s);
    Fs = s.As * ss;
    double eps2 = x > 0 ? r.eps_cu3 * (s.d2 - x) / x : 0.0;
    ss2 = s.As2 > 0 ? sig_steel(eps2) : 0.0;
    Fs2 = s.As2 * ss2;
  };
  auto resid = [&](double x) {
    double Fc, Fs, Fs2, zc, es, ss, ss2;
    forces(x, Fc, Fs, Fs2, zc, es, ss, ss2);
    return Fc + Fs + Fs2 - N;
  };

  // resid decreases with x (more concrete in compression): bracket and bisect.
  double lo = 1e-6 * s.h, hi = 4.0 * s.h;
  require(resid(lo) >= 0 && resid(hi) <= 0,
          "rc: the section cannot equilibrate this axial force");
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (lo + hi);
    (resid(mid) >= 0 ? lo : hi) = mid;
  }
  r.x = 0.5 * (lo + hi);

  double Fc, Fs, Fs2, zc, es, ss, ss2;
  forces(r.x, Fc, Fs, Fs2, zc, es, ss, ss2);
  r.block = std::min(r.lambda * r.x, s.h);
  r.Fc = Fc;
  r.eps_s = es;
  r.sigma_s = ss;
  r.sigma_s2 = ss2;
  r.steel_yields = std::fabs(ss) >= r.fyd * (1.0 - 1e-9);
  // M about the gross concrete centroid: M = SUM F_i z_i - N z_g
  double zg = s.centroid();
  r.M_Rd = Fc * zc + Fs * s.d + Fs2 * s.d2 - N * zg;
  return r;
}

}  // namespace keta

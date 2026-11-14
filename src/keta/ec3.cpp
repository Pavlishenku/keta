// EN 1993-1-1 / EN 1993-1-5 checks — see ec3.h.
#include "keta/ec3.h"

#include <algorithm>
#include <cmath>
#include <string>

namespace keta::ec3 {

double epsilon(double fy) {
  require(fy > 0, "ec3: fy must be positive");
  return std::sqrt(235e6 / fy);
}

int class_internal(double c_t, double eps, double alpha, double psi) {
  require(alpha > 0 && alpha <= 1, "ec3: alpha must lie in (0, 1]");
  const double l1 = alpha > 0.5 ? 396 * eps / (13 * alpha - 1) : 36 * eps / alpha;
  const double l2 = alpha > 0.5 ? 456 * eps / (13 * alpha - 1) : 41.5 * eps / alpha;
  const double l3 = psi > -1 ? 42 * eps / (0.67 + 0.33 * psi) : 62 * eps * (1 - psi) * std::sqrt(-psi);
  if (c_t <= l1) return 1;
  if (c_t <= l2) return 2;
  if (c_t <= l3) return 3;
  return 4;
}

int class_outstand(double c_t, double eps, double alpha, double psi, bool tip_in_compression) {
  require(alpha > 0 && alpha <= 1, "ec3: alpha must lie in (0, 1]");
  double l1, l2, l3;
  if (alpha == 1.0 && psi == 1.0) {  // uniform compression
    l1 = 9 * eps;
    l2 = 10 * eps;
    l3 = 14 * eps;
  } else {
    const double d = tip_in_compression ? alpha : alpha * std::sqrt(alpha);
    l1 = 9 * eps / d;
    l2 = 10 * eps / d;
    l3 = 21 * eps * std::sqrt(ksigma_outstand(psi, tip_in_compression));
  }
  if (c_t <= l1) return 1;
  if (c_t <= l2) return 2;
  if (c_t <= l3) return 3;
  return 4;
}

double c_flange(const ISection& s) { return s.b / 2 - s.tw / 2 - s.r; }
double c_web(const ISection& s) { return s.h - 2 * s.tf - 2 * s.r; }

double ksigma_internal(double psi) {
  require(psi >= -3 && psi <= 1, "ec3: psi outside [-3, 1] (EN 1993-1-5 Table 4.1)");
  if (psi == 1) return 4.0;
  if (psi > 0) return 8.2 / (1.05 + psi);
  if (psi == 0) return 7.81;
  if (psi > -1) return 7.81 - 6.29 * psi + 9.78 * psi * psi;
  if (psi == -1) return 23.9;
  return 5.98 * (1 - psi) * (1 - psi);
}

double ksigma_outstand(double psi, bool max_compression_at_free_edge) {
  if (max_compression_at_free_edge) {
    require(psi >= -1 && psi <= 1, "ec3: psi outside [-1, 1] (EN 1993-1-5 Table 4.2)");
    if (psi >= 0) return 0.578 / (psi + 0.34);
    return 1.70 - 5 * psi + 17.1 * psi * psi;
  }
  require(psi > -3 && psi <= 1, "ec3: psi outside (-3, 1] (EN 1993-1-5 Table 4.2)");
  if (psi == 1) return 0.43;
  return 0.57 - 0.21 * psi + 0.07 * psi * psi;
}

double lambda_p(double b_t, double eps, double ksigma) {
  return b_t / (28.4 * eps * std::sqrt(ksigma));
}

double rho_internal(double lp, double psi) {
  if (lp <= 0.5 + std::sqrt(0.085 - 0.055 * psi)) return 1.0;
  return std::min(1.0, (lp - 0.055 * (3 + psi)) / (lp * lp));
}

double rho_outstand(double lp) {
  if (lp <= 0.748) return 1.0;
  return std::min(1.0, (lp - 0.188) / (lp * lp));
}

double a_c(double t, double b1, double b2, double I_sl) {
  return 4.33 * std::pow(I_sl * b1 * b1 * b2 * b2 / (t * t * t * (b1 + b2)), 0.25);
}

double sigma_cr_sl(double E, double nu, double t, double a, double b1, double b2, double A_sl,
                   double I_sl) {
  const double b = b1 + b2;
  require(a <= a_c(t, b1, b2, I_sl),
          "ec3: a > a_c — the second formula of EN 1993-1-5 A.2.2(1) is not implemented");
  return kPi * kPi * E * I_sl / (A_sl * a * a) +
         E * t * t * t * b * a * a / (4 * kPi * kPi * (1 - nu * nu) * A_sl * b1 * b1 * b2 * b2);
}

void flat_stiffener(double t, double bp, double hs, double ts, double& A_sl, double& I_sl) {
  // plate strip centred on z = 0, stiffener on the face z = t/2
  const double Ap = bp * t, As = hs * ts, zs = t / 2 + hs / 2;
  A_sl = Ap + As;
  const double zg = As * zs / A_sl;
  I_sl = bp * t * t * t / 12 + Ap * zg * zg + ts * hs * hs * hs / 12 + As * (zs - zg) * (zs - zg);
}

RectProps rect_props(const std::vector<Rect>& add, const std::vector<Rect>& remove) {
  RectProps p;
  double Sy = 0, Sz = 0, Iyy = 0, Izz = 0;  // about the origin
  auto acc = [&](const Rect& r, double sg) {
    const double b = r.y1 - r.y0, h = r.z1 - r.z0, a = sg * b * h;
    const double yc = 0.5 * (r.y0 + r.y1), zc = 0.5 * (r.z0 + r.z1);
    p.A += a;
    Sy += a * zc;
    Sz += a * yc;
    Iyy += sg * (b * h * h * h / 12) + a * zc * zc;
    Izz += sg * (h * b * b * b / 12) + a * yc * yc;
  };
  for (const Rect& r : add) acc(r, 1.0);
  for (const Rect& r : remove) acc(r, -1.0);
  require(p.A > 0, "ec3: section of zero area");
  p.zg = Sy / p.A;
  p.yg = Sz / p.A;
  p.Iy = Iyy - p.A * p.zg * p.zg;
  p.Iz = Izz - p.A * p.yg * p.yg;
  return p;
}

EffectiveSection effective_uniform_compression(const std::vector<Rect>& gross,
                                               const std::vector<InternalPlate>& plates,
                                               double fy) {
  const double eps = epsilon(fy);
  EffectiveSection out;
  const RectProps g = rect_props(gross);
  out.A = g.A;
  out.yg = g.yg;
  out.zg = g.zg;
  std::vector<Rect> removed;
  for (const InternalPlate& p : plates) {
    const double lp = lambda_p(p.b / p.t, eps, ksigma_internal(1.0));
    const double rho = rho_internal(lp, 1.0);
    const double w = (1 - rho) * p.b;  // the central strip (be1 = be2 = 0.5 rho b, Table 4.1)
    out.lambda_p.push_back(lp);
    out.rho.push_back(rho);
    out.removed.push_back(w);
    if (w <= 0) continue;
    Rect r;
    if (p.along_y) {
      r = {p.yc - w / 2, p.yc + w / 2, p.zc - p.t / 2, p.zc + p.t / 2};
    } else {
      r = {p.yc - p.t / 2, p.yc + p.t / 2, p.zc - w / 2, p.zc + w / 2};
    }
    removed.push_back(r);
  }
  const RectProps e = rect_props(gross, removed);
  out.Aeff = e.A;
  out.yg_eff = e.yg;
  out.zg_eff = e.zg;
  out.Iy_eff = e.Iy;
  out.Iz_eff = e.Iz;
  return out;
}

double vpl_rd(double Av, double fy, double gM0) { return Av * fy / (std::sqrt(3.0) * gM0); }

double av_z_rolled_i(const ISection& s, double A, double eta) {
  const double hw = s.h - 2 * s.tf;
  return std::max(A - 2 * s.b * s.tf + (s.tw + 2 * s.r) * s.tf, eta * hw * s.tw);
}

double av_y_flanges(const ISection& s) { return 2 * s.b * s.tf; }

double rho_shear(double VEd, double VplRd) {
  const double v = std::fabs(VEd) / VplRd;
  return v <= 0.5 ? 0.0 : (2 * v - 1) * (2 * v - 1);
}

double my_v_rd_i(double Wply, double rho, double hw, double tw, double fy, double gM0) {
  const double Aw = hw * tw;
  return (Wply - rho * Aw * Aw / (4 * tw)) * fy / gM0;
}

double mz_v_rd_i(double Wplz, double rho, const ISection& s, double fy, double gM0) {
  return (Wplz - rho * 2 * s.tf * s.b * s.b / 4) * fy / gM0;
}

double nv_rd(double A, double Av, double rho, double fy, double gM0) { return (A - rho * Av) * fy / gM0; }

double mn_y_rd_i(double M, double n, double a) {
  a = std::min(a, 0.5);
  return std::min(M, M * (1 - n) / (1 - 0.5 * a));
}

bool axial_negligible_y(double NEd, double N, double hw, double tw, double rho, double fy,
                        double gM0) {
  return NEd <= 0.25 * N && NEd <= 0.5 * hw * tw * (1 - rho) * fy / gM0;
}

double imperfection(Curve c) {
  switch (c) {
    case Curve::a0: return 0.13;
    case Curve::a: return 0.21;
    case Curve::b: return 0.34;
    case Curve::c: return 0.49;
    case Curve::d: return 0.76;
  }
  return 0.76;
}

Curve curve_flexural_i(bool rolled, double h, double b, double tf, bool about_y) {
  if (rolled) {
    if (h / b > 1.2) {
      require(tf <= 0.1, "ec3: Table 6.2 row h/b > 1.2 limited to tf <= 100 mm");
      if (tf <= 0.04) return about_y ? Curve::a : Curve::b;
      return about_y ? Curve::b : Curve::c;
    }
    if (tf <= 0.1) return about_y ? Curve::b : Curve::c;
    return Curve::d;
  }
  if (tf <= 0.04) return about_y ? Curve::b : Curve::c;
  return about_y ? Curve::c : Curve::d;
}

Curve curve_lt_general_i(bool rolled, double h, double b) {
  if (rolled) return h / b <= 2 ? Curve::a : Curve::b;
  return h / b <= 2 ? Curve::c : Curve::d;
}

Curve curve_lt_rolled_i(bool rolled, double h, double b) {
  if (rolled) return h / b <= 2 ? Curve::b : Curve::c;
  return h / b <= 2 ? Curve::c : Curve::d;
}

double chi(double lambda, double alpha) {
  const double phi = 0.5 * (1 + alpha * (lambda - 0.2) + lambda * lambda);
  return std::min(1.0, 1.0 / (phi + std::sqrt(phi * phi - lambda * lambda)));
}

double chi_lt_rolled(double lambda, double alpha, double lambda0, double beta) {
  const double phi = 0.5 * (1 + alpha * (lambda - lambda0) + beta * lambda * lambda);
  const double c = 1.0 / (phi + std::sqrt(phi * phi - beta * lambda * lambda));
  return std::min({c, 1.0, 1.0 / (lambda * lambda)});
}

double f_modification(double kc, double lambda_lt) {
  const double d = lambda_lt - 0.8;
  return std::min(1.0, 1 - 0.5 * (1 - kc) * (1 - 2 * d * d));
}

double ncr(double EI, double Lcr) { return kPi * kPi * EI / (Lcr * Lcr); }

double ncr_torsional(double A, double Iy, double Iz, double G, double IT, double E, double Iw,
                     double LT) {
  return A / (Iy + Iz) * (G * IT + kPi * kPi * E * Iw / (LT * LT));
}

double mcr_din18800(double zeta, double Ncrz, double Iw, double IT, double Iz, double L, double zp) {
  const double c2 = (Iw + 0.039 * L * L * IT) / Iz;
  return zeta * Ncrz * (std::sqrt(c2 + 0.25 * zp * zp) + 0.5 * zp);
}

double cm0_linear(double psi, double NEd, double Ncr) {
  return 0.79 + 0.21 * psi + 0.36 * (psi - 0.33) * NEd / Ncr;
}

double cm0_transverse(double EI, double L, double delta, double M, double NEd, double Ncr) {
  return 1 + (kPi * kPi * EI * std::fabs(delta) / (L * L * std::fabs(M)) - 1) * NEd / Ncr;
}

AnnexAResult annex_a(const AnnexAInput& in) {
  AnnexAResult r;
  const double N = in.NEd;
  require(N >= 0 && N < in.Ncry && N < in.Ncrz, "ec3 annex A: NEd must lie below Ncr,y and Ncr,z");
  r.mu_y = (1 - N / in.Ncry) / (1 - in.chiy * N / in.Ncry);
  r.mu_z = (1 - N / in.Ncrz) / (1 - in.chiz * N / in.Ncrz);
  r.wy = std::min(in.Wply / in.Wely, 1.5);
  r.wz = std::min(in.Wplz / in.Welz, 1.5);
  const double NRk = in.A * in.fy;
  r.npl = N / (NRk / in.gM1);
  r.aLT = std::max(0.0, 1 - in.IT / in.Iy);
  const double Mply = in.Wply * in.fy, Mplz = in.Wplz * in.fy;  // Mpl,Rk
  r.lambda0 = std::sqrt(Mply / in.Mcr0);
  const double prodT = (1 - N / in.Ncrz) * (1 - N / in.NcrT);
  r.lambda0_lim = 0.2 * std::sqrt(in.C1) * std::pow(prodT, 0.25);
  r.eps_y = N > 0 ? in.MyEd / N * in.A / in.Wely : 0.0;
  if (r.lambda0 <= r.lambda0_lim) {
    r.Cmy = in.Cmy0;
    r.Cmz = in.Cmz0;
    r.CmLT = 1.0;
  } else {
    const double se = std::sqrt(r.eps_y) * r.aLT;
    r.Cmy = in.Cmy0 + (1 - in.Cmy0) * se / (1 + se);
    r.Cmz = in.Cmz0;
    r.CmLT = std::max(1.0, r.Cmy * r.Cmy * r.aLT / std::sqrt(prodT));
  }
  const double lmax = std::max(in.lambday, in.lambdaz);
  const double my = in.MyEd / (in.chiLT * Mply / in.gM1), mz = in.MzEd / (Mplz / in.gM1);
  const double lz4 = std::pow(in.lambdaz, 4);
  r.bLT = 0.5 * r.aLT * r.lambda0 * r.lambda0 * my * mz;
  r.cLT = 10 * r.aLT * r.lambda0 * r.lambda0 / (5 + lz4) * in.MyEd / (r.Cmy * in.chiLT * Mply / in.gM1);
  r.dLT = 2 * r.aLT * r.lambda0 / (0.1 + lz4) * in.MyEd / (r.Cmy * in.chiLT * Mply / in.gM1) * in.MzEd /
          (r.Cmz * Mplz / in.gM1);
  r.eLT = 1.7 * r.aLT * r.lambda0 / (0.1 + lz4) * in.MyEd / (r.Cmy * in.chiLT * Mply / in.gM1);
  const double wy = r.wy, wz = r.wz, n = r.npl;
  r.Cyy = std::max(in.Wely / in.Wply,
                   1 + (wy - 1) * ((2 - 1.6 / wy * r.Cmy * r.Cmy * lmax - 1.6 / wy * r.Cmy * r.Cmy * lmax * lmax) * n -
                                   r.bLT));
  r.Cyz = std::max(0.6 * std::sqrt(wz / wy) * in.Welz / in.Wplz,
                   1 + (wz - 1) * ((2 - 14 * r.Cmz * r.Cmz * lmax * lmax / std::pow(wz, 5)) * n - r.cLT));
  r.Czy = std::max(0.6 * std::sqrt(wy / wz) * in.Wely / in.Wply,
                   1 + (wy - 1) * ((2 - 14 * r.Cmy * r.Cmy * lmax * lmax / std::pow(wy, 5)) * n - r.dLT));
  // Table A.1 as printed: e_LT sits INSIDE the bracket multiplied by n_pl
  // (b_LT, c_LT, d_LT are subtracted after it) — ECCS 119 agrees.
  r.Czz = std::max(in.Welz / in.Wplz,
                   1 + (wz - 1) * (2 - 1.6 / wz * r.Cmz * r.Cmz * lmax - 1.6 / wz * r.Cmz * r.Cmz * lmax * lmax -
                                   r.eLT) * n);
  r.kyy = r.Cmy * r.CmLT * r.mu_y / (1 - N / in.Ncry) / r.Cyy;
  r.kyz = r.Cmz * r.mu_y / (1 - N / in.Ncrz) / r.Cyz * 0.6 * std::sqrt(wz / wy);
  r.kzy = r.Cmy * r.CmLT * r.mu_z / (1 - N / in.Ncry) / r.Czy * 0.6 * std::sqrt(wy / wz);
  r.kzz = r.Cmz * r.mu_z / (1 - N / in.Ncrz) / r.Czz;
  r.eq661 = N / (in.chiy * NRk / in.gM1) + r.kyy * in.MyEd / (in.chiLT * Mply / in.gM1) +
            r.kyz * in.MzEd / (Mplz / in.gM1);
  r.eq662 = N / (in.chiz * NRk / in.gM1) + r.kzy * in.MyEd / (in.chiLT * Mply / in.gM1) +
            r.kzz * in.MzEd / (Mplz / in.gM1);
  return r;
}

double cm_b3_linear(double psi) { return std::max(0.4, 0.6 + 0.4 * psi); }
double cm_b3_uniform_load_ah(double alpha_h) {
  require(alpha_h >= 0 && alpha_h <= 1, "ec3: Table B.3 alpha_h in [0, 1] here");
  return 0.95 + 0.05 * alpha_h;
}
double cm_b3_uniform_load_as(double alpha_s) {
  require(alpha_s >= 0 && alpha_s <= 1, "ec3: Table B.3 alpha_s in [0, 1] here");
  return std::max(0.4, 0.2 + 0.8 * alpha_s);
}

double kyy_b(const AnnexBInput& in) {
  const double n = in.NEd / (in.chiy * in.NRk / in.gM1);
  if (in.class34) return std::min(in.Cmy * (1 + 0.6 * in.lambday * n), in.Cmy * (1 + 0.6 * n));
  return std::min(in.Cmy * (1 + (in.lambday - 0.2) * n), in.Cmy * (1 + 0.8 * n));
}

double kzy_b2(const AnnexBInput& in) {
  const double n = in.NEd / (in.chiz * in.NRk / in.gM1);
  const double c = in.class34 ? 0.05 : 0.1;
  if (!in.class34 && in.lambdaz < 0.4)
    return std::min(0.6 + in.lambdaz, 1 - c * in.lambdaz * n / (in.CmLT - 0.25));
  return std::max(1 - c * in.lambdaz * n / (in.CmLT - 0.25), 1 - c * n / (in.CmLT - 0.25));
}

}  // namespace keta::ec3

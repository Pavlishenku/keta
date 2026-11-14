// Layered J2 plasticity of facet shells — see shell_plastic.h.
#include "keta/shell_plastic.h"

#include <cmath>
#include <limits>
#include <unordered_set>

#include <Eigen/SparseLU>

#include "keta/assembly.h"
#include "keta/bef.h"
#include "keta/corot_kin.h"
#include "keta/fiber_elem.h"
#include "keta/fiber_frame.h"  // NonConvergence
#include "keta/linsolve.h"

namespace keta {

namespace {
Mat3 plane_stress_C(double E, double nu) {
  Mat3 C;
  const double c = E / (1 - nu * nu);
  C << c, c * nu, 0, c * nu, c, 0, 0, 0, c * (1 - nu) / 2;
  return C;
}
Mat3 projector() {
  Mat3 P;
  P << 2, -1, 0, -1, 2, 0, 0, 0, 6;
  return P / 3.0;
}
}  // namespace

ShellLayerState j2_shell(const ShellLayerState& n, const Vec5& eps, double E, double nu, double fy, double H,
                         Mat5* tangent) {
  const Mat3 Cp = plane_stress_C(E, nu);
  const double G = E / (2 * (1 + nu)), a = E / (3 * (1 - nu));
  Vec5 epl;
  epl << n.eps_p, n.gamma_p;
  const Vec5 ee = eps - epl;
  Vec5 str;
  str.head<3>() = Cp * ee.head<3>();
  str.tail<2>() = G * ee.tail<2>();
  const double A1 = (str(0) + str(1)) * (str(0) + str(1)) / 6;
  const double A2 = 0.5 * (str(0) - str(1)) * (str(0) - str(1)) + 2 * (str(2) * str(2) + str(3) * str(3) + str(4) * str(4));
  auto fbar2 = [&](double dg) {
    return A1 / ((1 + a * dg) * (1 + a * dg)) + A2 / ((1 + 2 * G * dg) * (1 + 2 * G * dg));
  };
  Mat5 C = Mat5::Zero();
  C.topLeftCorner<3, 3>() = Cp;
  C(3, 3) = C(4, 4) = G;
  const double R0 = fy + H * n.alpha;
  ShellLayerState out = n;
  if (0.5 * fbar2(0) - R0 * R0 / 3 <= 0) {  // elastic
    out.sigma = str.head<3>();
    out.tau = str.tail<2>();
    if (tangent) *tangent = C;
    return out;
  }
  // consistency: phi(dg) = fbar2 / 2 - R(alpha_n + dg sqrt(2/3) fbar)^2 / 3,
  // phi(0) > 0, phi decreasing (Simo & Taylor): Newton safeguarded by a
  // bracket [lo, hi] grown until phi(hi) < 0, bisection when Newton leaves it
  const double s23 = std::sqrt(2.0 / 3.0);
  auto phi_of = [&](double dg, double* dphi) {
    const double f2 = fbar2(dg), f = std::sqrt(f2);
    const double df2 = -2 * a * A1 / std::pow(1 + a * dg, 3) - 4 * G * A2 / std::pow(1 + 2 * G * dg, 3);
    const double R = fy + H * (n.alpha + dg * s23 * f);
    const double dR = H * s23 * (f + dg * df2 / (2 * f));
    if (dphi) *dphi = 0.5 * df2 - 2.0 / 3.0 * R * dR;
    return 0.5 * f2 - R * R / 3;
  };
  double lo = 0, hi = 1 / (2 * G);
  for (int it = 0; it < 200 && phi_of(hi, nullptr) > 0; ++it) hi *= 2;
  require(phi_of(hi, nullptr) <= 0, "j2 shell: no plastic multiplier bracket");
  double dg = 0;
  for (int it = 0; it < 300; ++it) {
    double dphi = 0;
    const double phi = phi_of(dg, &dphi);
    const double R = fy + H * n.alpha;
    if (std::fabs(phi) <= 1e-15 * R * R) break;  // on the yield surface to round-off
    (phi > 0 ? lo : hi) = dg;
    double next = dg - phi / dphi;
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    if (std::fabs(next - dg) <= 1e-16 * std::fabs(next) || hi - lo <= 1e-16 * hi) {
      dg = next;
      break;
    }
    dg = next;
  }
  const double sum = (str(0) + str(1)) / (1 + a * dg), sh = 1 / (1 + 2 * G * dg);
  const double dif = (str(0) - str(1)) * sh;
  Vec5 sig;
  sig << 0.5 * (sum + dif), 0.5 * (sum - dif), str(2) * sh, str(3) * sh, str(4) * sh;
  out.sigma = sig.head<3>();
  out.tau = sig.tail<2>();
  Vec5 ee_new;
  ee_new.head<3>() = Cp.inverse() * sig.head<3>();
  ee_new.tail<2>() = sig.tail<2>() / G;
  const Vec5 ep_new = eps - ee_new;
  out.eps_p = ep_new.head<3>();
  out.gamma_p = ep_new.tail<2>();
  Mat5 P = Mat5::Zero();
  P.topLeftCorner<3, 3>() = projector();
  P(3, 3) = P(4, 4) = 2.0;
  const double sPs = sig.dot(P * sig);
  out.alpha = n.alpha + dg * s23 * std::sqrt(sPs);
  if (tangent) {
    Mat5 Cinv = Mat5::Zero();
    Cinv.topLeftCorner<3, 3>() = Cp.inverse();
    Cinv(3, 3) = Cinv(4, 4) = 1 / G;
    const Mat5 Xi = (Cinv + dg * P).inverse();
    const Vec5 XPs = Xi * (P * sig);
    const double beta = 2.0 / 3.0 * H * sPs / (1 - 2.0 / 3.0 * H * dg);
    *tangent = Xi - XPs * XPs.transpose() / ((P * sig).dot(XPs) + beta);
  }
  return out;
}

Mat5 hill_projector(const HillYield& h, double angle) {
  require(h.r11 > 0 && h.r22 > 0 && h.r33 > 0 && h.r12 > 0 && h.r13 > 0 && h.r23 > 0,
          "Hill criterion: positive yield ratios");
  const double i1 = 1 / (h.r11 * h.r11), i2 = 1 / (h.r22 * h.r22), i3 = 1 / (h.r33 * h.r33);
  const double F = 0.5 * (i2 + i3 - i1), G = 0.5 * (i3 + i1 - i2), Hh = 0.5 * (i1 + i2 - i3);
  const double L = 0.5 / (h.r23 * h.r23), M = 0.5 / (h.r13 * h.r13), N = 0.5 / (h.r12 * h.r12);
  // s' A s = R^2 with s33 = 0; P = 2/3 A
  Mat5 A = Mat5::Zero();
  A(0, 0) = G + Hh;
  A(1, 1) = F + Hh;
  A(0, 1) = A(1, 0) = -Hh;
  A(2, 2) = 2 * N;
  A(3, 3) = 2 * M;
  A(4, 4) = 2 * L;
  // the plane-stress quadratic form must be positive SEMI-definite: a
  // degenerate locus (an open cylinder, e.g. r33 = r11 > r22: only 0.5 s11 -
  // s22 counts) is Hill's own and is accepted; indefinite is not
  const double sc = std::max({std::fabs(A(0, 0)), std::fabs(A(1, 1)), std::fabs(Hh)});
  require(A(0, 0) >= -1e-12 * sc && A(1, 1) >= -1e-12 * sc && A(0, 0) * A(1, 1) - Hh * Hh >= -1e-12 * sc * sc &&
              A(0, 0) + A(1, 1) > 0,
          "Hill criterion: the plane-stress quadratic form is indefinite (inconsistent yield ratios)");
  // local stresses -> material axes
  const double c = std::cos(angle), sn = std::sin(angle);
  Mat5 T = Mat5::Zero();
  T(0, 0) = c * c;
  T(0, 1) = sn * sn;
  T(0, 2) = 2 * c * sn;
  T(1, 0) = sn * sn;
  T(1, 1) = c * c;
  T(1, 2) = -2 * c * sn;
  T(2, 0) = -c * sn;
  T(2, 1) = c * sn;
  T(2, 2) = c * c - sn * sn;
  T(3, 3) = c;  // t13 = c txz + s tyz
  T(3, 4) = sn;
  T(4, 3) = -sn;  // t23 = -s txz + c tyz
  T(4, 4) = c;
  return 2.0 / 3.0 * T.transpose() * A * T;
}

ShellLayerState hill_shell(const ShellLayerState& n, const Vec5& eps, double E, double nu, double fy, double H,
                           const Mat5& P, Mat5* tangent) {
  const double G = E / (2 * (1 + nu));
  Mat5 C = Mat5::Zero();
  C.topLeftCorner<3, 3>() = plane_stress_C(E, nu);
  C(3, 3) = C(4, 4) = G;
  return general_shell(n, eps, C, fy, H, P, tangent);
}

Mat5 j2_projector5() {
  Mat5 P = Mat5::Zero();
  P.topLeftCorner<3, 3>() = projector();
  P(3, 3) = P(4, 4) = 2.0;
  return P;
}

ShellLayerState general_shell(const ShellLayerState& n, const Vec5& eps, const Mat5& C, double fy, double H,
                              const Mat5& P, Mat5* tangent) {
  const Mat5 Cinv = C.inverse();
  const double G = std::max(C(3, 3), C(2, 2));  // a modulus scale for the multiplier bracket
  Vec5 epl;
  epl << n.eps_p, n.gamma_p;
  const Vec5 str = C * (eps - epl);
  const Vec5 ctr = eps - epl;  // C^-1 s_trial
  const double R0 = fy + H * n.alpha;
  ShellLayerState out = n;
  if (0.5 * str.dot(P * str) - R0 * R0 / 3 <= 0) {  // elastic
    out.sigma = str.head<3>();
    out.tau = str.tail<2>();
    if (tangent) *tangent = C;
    return out;
  }
  const double s23 = std::sqrt(2.0 / 3.0);
  Mat5 Xi;
  Vec5 sig;
  auto phi_of = [&](double dg, double* dphi) {
    Xi = (Cinv + dg * P).inverse();
    sig = Xi * ctr;
    const Vec5 Ps = P * sig;
    const double f2 = std::max(sig.dot(Ps), 0.0), f = std::sqrt(f2);
    const double df2 = -2 * Ps.dot(Xi * Ps);
    const double R = fy + H * (n.alpha + dg * s23 * f);
    const double dR = H * s23 * (f + (f > 0 ? dg * df2 / (2 * f) : 0.0));
    if (dphi) *dphi = 0.5 * df2 - 2.0 / 3.0 * R * dR;
    return 0.5 * f2 - R * R / 3;
  };
  double lo = 0, hi = 1 / (2 * G);
  for (int it = 0; it < 200 && phi_of(hi, nullptr) > 0; ++it) hi *= 2;
  require(phi_of(hi, nullptr) <= 0, "Hill shell: no plastic multiplier bracket");
  double dg = 0;
  for (int it = 0; it < 300; ++it) {
    double dphi = 0;
    const double phi = phi_of(dg, &dphi);
    if (std::fabs(phi) <= 1e-15 * R0 * R0) break;
    (phi > 0 ? lo : hi) = dg;
    double next = dphi < 0 ? dg - phi / dphi : 0.5 * (lo + hi);
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    if (std::fabs(next - dg) <= 1e-16 * std::fabs(next) || hi - lo <= 1e-16 * hi) {
      dg = next;
      break;
    }
    dg = next;
  }
  phi_of(dg, nullptr);  // Xi, sig at the root
  out.sigma = sig.head<3>();
  out.tau = sig.tail<2>();
  const Vec5 ep_new = eps - Cinv * sig;
  out.eps_p = ep_new.head<3>();
  out.gamma_p = ep_new.tail<2>();
  const double sPs = sig.dot(P * sig);
  out.alpha = n.alpha + dg * s23 * std::sqrt(sPs);
  if (tangent) {
    const Vec5 XPs = Xi * (P * sig);
    const double beta = 2.0 / 3.0 * H * sPs / (1 - 2.0 / 3.0 * H * dg);
    *tangent = Xi - XPs * XPs.transpose() / ((P * sig).dot(XPs) + beta);
  }
  return out;
}

ShellLayerState j2_plane_stress(const ShellLayerState& n, const Vec3& eps, double E, double nu, double fy,
                                double H, Mat3* tangent) {
  // pure plane stress: the transverse shear strains held at their plastic
  // values (zero shear stress); the shear blocks decouple at zero shear
  Vec5 e5;
  e5 << eps, n.gamma_p;
  Mat5 t5;
  const ShellLayerState s = j2_shell(n, e5, E, nu, fy, H, tangent ? &t5 : nullptr);
  if (tangent) *tangent = t5.topLeftCorner<3, 3>();
  return s;
}

namespace {
using ShellState = std::vector<std::vector<ShellLayerState>>;

// composite 4-point Gauss-Lobatto through the thickness (unit coordinates):
// exact to degree 5 per layer, so the elastic z^2 bending and the quartic
// parabolic-shear energy are integrated exactly; the faces are points
struct ThicknessRule {
  std::vector<double> zeta, w;
  std::vector<int> ply;  // the ply of each point (ShellNlOps::plies), -1: isotropic
};
// a laminate / orthotropic facet: each ply integrated on its own (its faces
// are points of both neighbours: the material jumps there), with about
// `layers` sub-layers over the whole thickness
ThicknessRule ply_rule(const ShellNlOps& o, int layers) {
  ThicknessRule r;
  const double xi[4] = {-1.0, -1.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0), 1.0};
  const double wi[4] = {1.0 / 6, 5.0 / 6, 5.0 / 6, 1.0 / 6};
  for (std::size_t q = 0; q < o.plies.size(); ++q) {
    const double z0 = o.plies[q].z0, z1 = o.plies[q].z1;
    const int ns = std::max(1, static_cast<int>(std::lround(layers * (z1 - z0))));
    const double h = (z1 - z0) / ns;
    r.zeta.push_back(z0);
    r.w.push_back(0.0);
    r.ply.push_back(static_cast<int>(q));
    for (int l = 0; l < ns; ++l) {
      const double zm = z0 + (l + 0.5) * h;
      r.w.back() += wi[0] * h / 2;
      for (int k = 1; k < 4; ++k) {
        r.zeta.push_back(zm + xi[k] * h / 2);
        r.w.push_back(wi[k] * h / 2);
        r.ply.push_back(static_cast<int>(q));
      }
    }
  }
  return r;
}
ThicknessRule lobatto4(int layers) {
  ThicknessRule r;
  const double xi[4] = {-1.0, -1.0 / std::sqrt(5.0), 1.0 / std::sqrt(5.0), 1.0};
  const double wi[4] = {1.0 / 6, 5.0 / 6, 5.0 / 6, 1.0 / 6};
  const double h = 1.0 / layers;
  r.zeta.push_back(-0.5);
  r.w.push_back(0.0);
  for (int l = 0; l < layers; ++l) {
    const double zm = -0.5 + (l + 0.5) * h;
    r.w.back() += wi[0] * h / 2;
    for (int k = 1; k < 4; ++k) {
      r.zeta.push_back(zm + xi[k] * h / 2);
      r.w.push_back(wi[k] * h / 2);
    }
  }
  return r;
}

// Local response of one layered facet: d the 18 local element DOFs (element
// order, local frame), f the local internal force, K the local tangent.
// Generalized strains per mid-edge point: membrane e0 (Allman), curvature
// kappa (constant), transverse shear gbar (MITC3); fibre strains e0 + z kappa
// and f(z) gbar, f = 5/4 (1 - 4 z^2 / t^2) (energy-consistent with kappa_s =
// 5/6: the elastic limit is ke() exactly).
void shell_local_response(const ShellNlOps& o, const Vec18& d, const ShellState& committed, ShellState& trial,
                          double E, double nu, double fy, double H, const ThicknessRule& rule, bool plastic_shear,
                          Vec18& f, Mat18* K, const Mat5* Phill = nullptr) {
  if (o.offset || o.k_bed.any()) {
    // enriched facet: the layered mid-surface at te d, pulled back to the
    // node DOFs (the offset is a rigid arm), plus the elastic bed
    ShellNlOps om = o;
    om.offset = false;
    om.k_bed.setZero();
    Vec18 fm;
    Mat18 km;
    shell_local_response(om, o.te * d, committed, trial, E, nu, fy, H, rule, plastic_shear, fm, K ? &km : nullptr, Phill);
    f = o.te.transpose() * fm + o.k_bed * d;
    if (K) *K = o.te.transpose() * km * o.te + o.k_bed;
    return;
  }
  Eigen::Matrix<double, 9, 1> um, ub;
  for (int i = 0; i < 9; ++i) {
    um(i) = d(o.mem_idx[static_cast<std::size_t>(i)]);
    ub(i) = d(o.bend_idx[static_cast<std::size_t>(i)]);
  }
  // generalized-strain operator in the local 18-DOF space: rows (e0, kappa, gbar)
  auto Bgen = [&](int p) {
    Eigen::Matrix<double, 8, 18> B = Eigen::Matrix<double, 8, 18>::Zero();
    for (int i = 0; i < 9; ++i) {
      const int mi = o.mem_idx[static_cast<std::size_t>(i)], bi = o.bend_idx[static_cast<std::size_t>(i)];
      B.block<3, 1>(0, mi) = o.bm[static_cast<std::size_t>(p)].col(i);
      B.block<3, 1>(3, bi) = o.bb.col(i);
      B.block<2, 1>(6, bi) = o.bs[static_cast<std::size_t>(p)].col(i);
    }
    return B;
  };
  f.setZero();
  if (K) K->setZero();
  // elastic parts: drilling stabilisation; transverse shear when not plastic
  {
    const Eigen::Matrix<double, 9, 1> fm = o.k_drill * um;
    for (int i = 0; i < 9; ++i) f(o.mem_idx[static_cast<std::size_t>(i)]) += fm(i);
    if (K)
      for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
          (*K)(o.mem_idx[static_cast<std::size_t>(i)], o.mem_idx[static_cast<std::size_t>(j)]) += o.k_drill(i, j);
    if (!plastic_shear) {
      const Eigen::Matrix<double, 9, 1> fs = o.k_shear * ub;
      for (int i = 0; i < 9; ++i) f(o.bend_idx[static_cast<std::size_t>(i)]) += fs(i);
      if (K)
        for (int i = 0; i < 9; ++i)
          for (int j = 0; j < 9; ++j)
            (*K)(o.bend_idx[static_cast<std::size_t>(i)], o.bend_idx[static_cast<std::size_t>(j)]) += o.k_shear(i, j);
    }
  }
  const std::size_t nz = rule.zeta.size();
  for (int p = 0; p < static_cast<int>(o.w.size()); ++p) {
    const double t = o.tp[static_cast<std::size_t>(p)], wA = o.w[static_cast<std::size_t>(p)];
    const Eigen::Matrix<double, 8, 18> B = Bgen(p);
    const Eigen::Matrix<double, 8, 1> Eg = B * d;
    Eigen::Matrix<double, 8, 1> S = Eigen::Matrix<double, 8, 1>::Zero();
    Eigen::Matrix<double, 8, 8> D = Eigen::Matrix<double, 8, 8>::Zero();
    for (std::size_t k = 0; k < nz; ++k) {
      const double zt = rule.zeta[k], z = zt * t, w = rule.w[k] * t;
      const double prof = plastic_shear ? 1.25 * (1 - 4 * zt * zt) : 0.0;
      const ShellLayerState& c = committed[static_cast<std::size_t>(p)][k];
      Vec5 e5;
      e5.head<3>() = Eg.head<3>() + z * Eg.segment<3>(3);
      if (plastic_shear)
        e5.tail<2>() = prof * Eg.tail<2>();
      else
        e5.tail<2>() = c.gamma_p;  // no transverse shear stress in the criterion
      Mat5 Ct;
      const int ply = rule.ply.empty() ? -1 : rule.ply[k];
      const ShellLayerState sl =
          ply >= 0 ? general_shell(c, e5, o.plies[static_cast<std::size_t>(ply)].C, fy, H, Phill ? *Phill : j2_projector5(),
                                   K ? &Ct : nullptr)
          : Phill  ? hill_shell(c, e5, E, nu, fy, H, *Phill, K ? &Ct : nullptr)
                   : j2_shell(c, e5, E, nu, fy, H, K ? &Ct : nullptr);
      trial[static_cast<std::size_t>(p)][k] = sl;
      // S(z): E8 -> eps5
      Eigen::Matrix<double, 5, 8> Sz = Eigen::Matrix<double, 5, 8>::Zero();
      Sz.block<3, 3>(0, 0).setIdentity();
      Sz.block<3, 3>(0, 3) = z * Mat3::Identity();
      Sz.block<2, 2>(3, 6) = prof * Eigen::Matrix2d::Identity();
      Vec5 sig;
      sig << sl.sigma, sl.tau;
      S += w * Sz.transpose() * sig;
      if (K) D += w * Sz.transpose() * Ct * Sz;
    }
    f += wA * B.transpose() * S;
    if (K) *K += wA * B.transpose() * D * B;
  }
}
}  // namespace

namespace {
// A concentrated hinge in the bar basic system: m components (comp, sign)
// of v, elastic Ke, kinematic hardening H, facets n . (s - beta) <= r.
struct BasicHinge {
  int m = 1;
  int comp[2] = {0, 0};
  double sign[2] = {1, 1};
  Eigen::Matrix2d Ke = Eigen::Matrix2d::Zero(), H = Eigen::Matrix2d::Zero();
  std::vector<Eigen::Vector2d> n;
  std::vector<double> r;
  Eigen::Vector2d ep = Eigen::Vector2d::Zero(), beta = Eigen::Vector2d::Zero();      // committed
  Eigen::Vector2d ep_t = Eigen::Vector2d::Zero(), beta_t = Eigen::Vector2d::Zero();  // trial
};

// Backward-Euler multi-surface return on a polygon with linear kinematic
// hardening (Simo & Hughes, Computational Inelasticity, sec. 2.6 / 5.2):
// active set iterated until every multiplier is >= 0 and every facet
// admissible; exact for linear facets. s, algorithmic tangent C.
void hinge_return(BasicHinge& h, const Eigen::Vector2d& e, Eigen::Vector2d& s, Eigen::Matrix2d& C) {
  const int m = h.m;
  const Eigen::MatrixXd Ke = h.Ke.topLeftCorner(m, m), Hh = h.H.topLeftCorner(m, m);
  const Eigen::VectorXd str = Ke * (e - h.ep).head(m), xi = str - h.beta.head(m);
  std::vector<int> act;
  Eigen::VectorXd dl;
  Eigen::MatrixXd N, Ginv;
  auto f_at = [&](int j, const Eigen::VectorXd& sv, const Eigen::VectorXd& bv) {
    return h.n[static_cast<std::size_t>(j)].head(m).dot(sv - bv) - h.r[static_cast<std::size_t>(j)];
  };
  Eigen::VectorXd sv = str, bv = h.beta.head(m), dep = Eigen::VectorXd::Zero(m);
  for (std::size_t pass = 0; pass < 4 * h.n.size() + 8; ++pass) {
    if (!act.empty()) {
      N.resize(static_cast<Eigen::Index>(act.size()), m);
      Eigen::VectorXd ftr(static_cast<Eigen::Index>(act.size()));
      for (std::size_t k = 0; k < act.size(); ++k) {
        N.row(static_cast<Eigen::Index>(k)) = h.n[static_cast<std::size_t>(act[k])].head(m).transpose();
        ftr(static_cast<Eigen::Index>(k)) = f_at(act[k], str, h.beta.head(m));
      }
      const Eigen::MatrixXd G = N * (Ke + Hh) * N.transpose();
      const Eigen::FullPivLU<Eigen::MatrixXd> lu(G);
      require(lu.isInvertible(), "hinge return: dependent active facets");
      Ginv = lu.inverse();
      dl = Ginv * ftr;
      Eigen::Index kmin;
      if (dl.minCoeff(&kmin) < 0) {  // drop the most negative multiplier
        act.erase(act.begin() + kmin);
        continue;
      }
      dep = N.transpose() * dl;
      sv = str - Ke * dep;
      bv = h.beta.head(m) + Hh * dep;
    } else {
      sv = str;
      bv = h.beta.head(m);
      dep.setZero();
    }
    int worst = -1;
    double fw = 0;
    for (int j = 0; j < static_cast<int>(h.n.size()); ++j) {
      if (std::find(act.begin(), act.end(), j) != act.end()) continue;
      const double fj = f_at(j, sv, bv);
      if (fj > 1e-12 * h.r[static_cast<std::size_t>(j)] && fj > fw) {
        fw = fj;
        worst = j;
      }
    }
    if (worst < 0) break;
    require(static_cast<int>(act.size()) < m, "hinge return: more active facets than deformations");
    act.push_back(worst);
  }
  s.setZero();
  s.head(m) = sv;
  h.ep_t = h.ep;
  h.ep_t.head(m) += dep;
  h.beta_t.setZero();
  h.beta_t.head(m) = bv;
  C.setZero();
  C.topLeftCorner(m, m) = act.empty() ? Ke : Eigen::MatrixXd(Ke - Ke * N.transpose() * Ginv * N * Ke);
}

// A nonlinear bar in the basic system: fiber element, or elastic basic
// stiffness with series hinges (internal hinge deformations solved per call).
struct NlBar {
  std::size_t e = 0;
  bool fiber = false;
  fiberel::FiberElement3D fe;
  Mat6 kb = Mat6::Zero();       // elastic basic stiffness (hinged bars)
  std::vector<BasicHinge> hs;
  Eigen::VectorXd ehinge;       // warm start of the hinge deformations
  fiberel::Mat612 Blin = fiberel::Mat612::Zero();
  corot::CorotElem ce;          // large displacements
  int gdof[12];
  double L = 0;
  // ---- 12-DOF hinged bar (linear kinematics): the bar's full local kernel
  // (release-condensed, tapered), rigid arms through T, the Winkler bed as
  // its separate global part, span loads f0c (release-condensed, per unit
  // load factor); hinges on local DOFs with internal DOFs t solved per call
  bool mode12 = false;
  // fiber bars: span-load particular field (unit), basic reactions, the
  // release-condensed equivalent load the case vector carries
  std::vector<Vec3> sp_unit;
  Vec12 p_ss = Vec12::Zero();
  bool span = false;
  double GJ = 0;                    // Saint-Venant torsion stiffness of the bar
  // imposed strains / dislocations (per unit load factor) and the twist of
  // the span torques (tor_int / GJ)
  std::vector<Vec3> e0_unit;
  Vec6 vshift_unit = Vec6::Zero();
  // torsion OUTSIDE the fiber element (torsional bed, linear kinematics):
  // the exact beam-on-foundation torsion block (bef.h) on (rx_i, rx_j),
  // release-condensed, and its fixed-end torques per unit load factor
  bool tor_ext = false;
  Eigen::Matrix2d ktor = Eigen::Matrix2d::Zero();
  Eigen::Vector2d f0tor = Eigen::Vector2d::Zero();
  Vec12 tor_cur = Vec12::Zero();    // its current local end torques
  // ---- span loads in LARGE displacements, per load: its corotational record
  // (local intensities w(q): E^T g dead, constant follower), the basic-system
  // reactions per unit local intensity (pushed back like equivalent loads),
  // and the particular fields per unit intensity (fiber) or the basic
  // deformation shift per unit intensity (hinged: elastic basic system,
  // q = kb (v - v0) with kb v0 = -q0, q0 the fixed-end basic forces)
  struct LdLoad {
    corot::CorotBarLoad w, ss;
    std::array<std::vector<Vec3>, 9> S;
    Eigen::Matrix<double, 6, 9> V0 = Eigen::Matrix<double, 6, 9>::Zero();
    Eigen::Matrix<double, 6, 9> Qn = Eigen::Matrix<double, 6, 9>::Zero();  // hinged: determinate part
  };
  std::vector<LdLoad> ldl;
  std::array<bool, 6> rel6{};       // released basic components
  Vec12 pss_cur = Vec12::Zero();    // current basic-system reactions (local, scaled)
  Mat12 kloc = Mat12::Zero(), Tn = Mat12::Identity(), ksoil = Mat12::Zero();
  Vec12 f0c = Vec12::Zero();
  std::vector<int> q12;         // hinged local DOFs, hinge order
  Eigen::VectorXd t12;          // internal DOFs (warm start)
  // global internal force (with the bar's own equivalent load added back, as
  // the case vector carries it) and tangent at nodal displacements ug
  void response12(const Vec12& ug, double lam, Vec12& f, Mat12& K, Vec12& p_out) {
    const int M = static_cast<int>(q12.size());
    const Vec12 ul = Tn * ug;
    if (t12.size() != M) {
      t12.resize(M);
      for (int a = 0; a < M; ++a) t12(a) = ul(q12[static_cast<std::size_t>(a)]);
    }
    Eigen::VectorXd S(M);
    Eigen::MatrixXd Cb = Eigen::MatrixXd::Zero(M, M);
    auto eval = [&](const Eigen::VectorXd& t, Eigen::VectorXd& r) {
      Vec12 v = ul;
      for (int a = 0; a < M; ++a) v(q12[static_cast<std::size_t>(a)]) = t(a);
      const Vec12 pb = kloc * v - lam * f0c;
      int o = 0;
      for (BasicHinge& h : hs) {
        Eigen::Vector2d eh = Eigen::Vector2d::Zero(), sh;
        Eigen::Matrix2d Ch;
        for (int a = 0; a < h.m; ++a) eh(a) = ul(h.comp[a]) - t(o + a);
        hinge_return(h, eh, sh, Ch);
        S.segment(o, h.m) = sh.head(h.m);
        Cb.block(o, o, h.m, h.m) = Ch.topLeftCorner(h.m, h.m);
        o += h.m;
      }
      r.resize(M);
      for (int a = 0; a < M; ++a) r(a) = pb(q12[static_cast<std::size_t>(a)]) - S(a);
      return pb;
    };
    Eigen::MatrixXd kqq(M, M);
    for (int a = 0; a < M; ++a)
      for (int b = 0; b < M; ++b) kqq(a, b) = kloc(q12[static_cast<std::size_t>(a)], q12[static_cast<std::size_t>(b)]);
    Eigen::VectorXd r;
    Vec12 pb = eval(t12, r);
    for (int it = 0;; ++it) {
      double sc = 1e-300;
      for (int a = 0; a < M; ++a) {
        double row = 0;
        for (int b = 0; b < 12; ++b) row += std::fabs(kloc(q12[static_cast<std::size_t>(a)], b)) * std::fabs(ul(b));
        sc = std::max(sc, row + std::fabs(lam * f0c(q12[static_cast<std::size_t>(a)])) + std::fabs(S(a)));
      }
      if (r.lpNorm<Eigen::Infinity>() <= 64 * std::numeric_limits<double>::epsilon() * sc) break;
      require(it < 200, "shell plasticity: hinged bar state determination did not converge");
      const Eigen::MatrixXd J = kqq + Cb;  // dr/dt
      const Eigen::FullPivLU<Eigen::MatrixXd> lu(J);
      require(lu.isInvertible(), "shell plasticity: singular hinged bar (mechanism within the bar)");
      const Eigen::VectorXd dt = -lu.solve(r);
      double a = 1.0;
      Eigen::VectorXd rt;
      Vec12 pt;
      for (int back = 0; back < 40; ++back, a *= 0.5) {
        pt = eval(t12 + a * dt, rt);
        if (rt.norm() < r.norm()) break;
      }
      t12 += a * dt;
      r = rt;
      pb = pt;
      if ((a * dt).norm() <= 16 * std::numeric_limits<double>::epsilon() * t12.norm()) break;
    }
    // node forces: the bar's at the free DOFs, the hinge forces at the hinged ones
    p_out = pb;
    for (int a = 0; a < M; ++a) p_out(q12[static_cast<std::size_t>(a)]) = S(a);
    f = Tn.transpose() * (p_out + lam * f0c) + ksoil * ug;
    // tangent: the augmented system with the algorithmic hinge stiffness, condensed
    MatX Ka = MatX::Zero(12 + M, 12 + M);
    std::vector<int> map(12);
    for (int a = 0; a < 12; ++a) map[static_cast<std::size_t>(a)] = a;
    for (int a = 0; a < M; ++a) map[static_cast<std::size_t>(q12[static_cast<std::size_t>(a)])] = 12 + a;
    for (int a = 0; a < 12; ++a)
      for (int b = 0; b < 12; ++b) Ka(map[static_cast<std::size_t>(a)], map[static_cast<std::size_t>(b)]) += kloc(a, b);
    for (int a = 0; a < M; ++a)
      for (int b = 0; b < M; ++b) {
        const int qa = q12[static_cast<std::size_t>(a)], qb = q12[static_cast<std::size_t>(b)];
        Ka(qa, qb) += Cb(a, b);
        Ka(12 + a, 12 + b) += Cb(a, b);
        Ka(qa, 12 + b) -= Cb(a, b);
        Ka(12 + a, qb) -= Cb(a, b);
      }
    const MatX Kii = Ka.bottomRightCorner(M, M);
    const Eigen::FullPivLU<MatX> lk(Kii);
    require(lk.isInvertible(), "shell plasticity: singular hinged bar condensation (mechanism within the bar)");
    const Mat12 Kc = Ka.topLeftCorner(12, 12) - Ka.topRightCorner(12, M) * lk.solve(Ka.bottomLeftCorner(M, 12));
    K = Tn.transpose() * Kc * Tn + ksoil;
  }
  // basic response q(v), tangent dq/dv
  // (hinged bars: qn, a statically determined basic force outside the
  // range of kb — the part of the span loads' fixed-end forces a released
  // member carries by equilibrium alone — added to q; Gn = dq/dqn)
  void response(const Vec6& v, Vec6& q, Mat6& Kq, const Vec6& qn = Vec6::Zero(), Mat6* Gn = nullptr) {
    if (fiber) {
      fe.determine(v);
      q = fe.q6();
      Kq = fe.kb;
      return;
    }
    int M = 0;
    for (const BasicHinge& h : hs) M += h.m;
    if (M == 0) {
      q = kb * v + qn;
      Kq = kb;
      if (Gn) *Gn = Mat6::Identity();
      return;
    }
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(6, M);
    int off = 0;
    for (const BasicHinge& h : hs) {
      for (int a = 0; a < h.m; ++a) P(h.comp[a], off + a) = h.sign[a];
      off += h.m;
    }
    if (ehinge.size() != M) ehinge = Eigen::VectorXd::Zero(M);
    Eigen::VectorXd S(M);
    Eigen::MatrixXd Cb = Eigen::MatrixXd::Zero(M, M);
    const Eigen::MatrixXd PtK = P.transpose() * kb, PtKP = PtK * P;
    // r(e) = P^T kb (v - P e) - S(e): the bar force equals the hinge force;
    // Newton on the piecewise-linear r with a residual-decrease safeguard
    // (halving) so a kink of the hinge law cannot make it cycle
    auto eval = [&](const Eigen::VectorXd& eh_all, Eigen::VectorXd& r, Eigen::MatrixXd* J) {
      int o = 0;
      for (BasicHinge& h : hs) {
        Eigen::Vector2d eh = Eigen::Vector2d::Zero(), sh;
        Eigen::Matrix2d Ch;
        eh.head(h.m) = eh_all.segment(o, h.m);
        hinge_return(h, eh, sh, Ch);
        S.segment(o, h.m) = sh.head(h.m);
        Cb.block(o, o, h.m, h.m) = Ch.topLeftCorner(h.m, h.m);
        o += h.m;
      }
      r = PtK * v + P.transpose() * qn - PtKP * eh_all - S;
      if (J) *J = PtKP + Cb;
    };
    Eigen::VectorXd r;
    Eigen::MatrixXd J;
    eval(ehinge, r, &J);
    for (int it = 0;; ++it) {
      // round-off bound of r: |P^T kb| |v| + |P^T kb P| |e| + |S| (absolute
      // values entrywise — the products cancel inside)
      const double sc = (PtK.cwiseAbs() * v.cwiseAbs() + P.cwiseAbs().transpose() * qn.cwiseAbs() +
                         PtKP.cwiseAbs() * ehinge.cwiseAbs() + S.cwiseAbs())
                            .lpNorm<Eigen::Infinity>() + 1e-300;
      if (r.lpNorm<Eigen::Infinity>() <= 64 * std::numeric_limits<double>::epsilon() * sc) break;
      require(it < 200, "shell plasticity: hinge state determination did not converge");
      const Eigen::FullPivLU<Eigen::MatrixXd> lu(J);
      require(lu.isInvertible(), "shell plasticity: singular hinged bar (mechanism within the bar)");
      const Eigen::VectorXd de = lu.solve(r);
      double a = 1.0;
      Eigen::VectorXd rt;
      Eigen::MatrixXd Jt;
      for (int back = 0; back < 40; ++back, a *= 0.5) {
        eval(ehinge + a * de, rt, &Jt);
        if (rt.norm() < r.norm()) break;
      }
      ehinge += a * de;
      r = rt;
      J = Jt;
      // converged in deformation: the hinge forces carry the round-off of
      // Ke (e - ep) (stiff elastic hinges), below which r cannot drop
      if ((a * de).norm() <= 16 * std::numeric_limits<double>::epsilon() * ehinge.norm()) break;
    }
    const Eigen::FullPivLU<Eigen::MatrixXd> lu(J);
    require(lu.isInvertible(), "shell plasticity: singular hinged bar (mechanism within the bar)");
    q = kb * (v - P * ehinge) + qn;
    Kq = kb - kb * P * lu.solve(PtK);
    if (Gn) *Gn = Mat6::Identity() - kb * P * lu.solve(Eigen::MatrixXd(P.transpose()));
  }
  void commit() {
    if (fiber) fe.commit_states();
    for (BasicHinge& h : hs) {
      h.ep = h.ep_t;
      h.beta = h.beta_t;
    }
  }
};
}  // namespace

ShellPlasticResult solve_shell_plastic(const Model& model, Id case_id, const ShellPlasticOptions& opt) {
  require(model.continua.empty(), "solve_shell_plastic: continuum elements are not supported (linear statics, modal and dynamics take them)");
  require(opt.steps >= 1 && opt.layers >= 1, "shell plasticity: steps and layers must be >= 1");
  // node-to-node links (small displacements only): the elastic links join
  // the linear part, the rigid links reduce every correction exactly
  // (u = T u_r: the tangent factored as T^t K T, the residual measured on the
  // reduced space — a slave's own residual is its constraint force)
  const Assembled asm_ = assemble_elements(model, !opt.large_displacements);
  const DofMap& d = asm_.dofs;
  const bool red = asm_.red.active;
  // the linear part: bars (their exact kernel — releases, rigid arms,
  // Winkler bed, tapered members — as linear statics assembles it), the
  // per-DOF springs and the full 6x6 elastic supports. Its internal force is
  // K_lin u at every iterate, its tangent K_lin.
  std::vector<Triplet> lin_trips;
  auto bar_gdofs = [&](const Bar& b, int g[12]) {
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (int c = 0; c < 6; ++c) {
      g[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
      g[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
    }
  };
  const bool ld = opt.large_displacements;
  std::vector<corot::CorotElem> cbars;
  const CaseLoads cl_bars = lower_case(model, asm_, model.load_case(case_id));  // span records
  // nonlinear bars (fibers / hinges), in the basic system
  std::map<Id, std::vector<const HingeSpec*>> hinges_of;
  for (const HingeSpec& hs : opt.hinges) hinges_of[hs.bar].push_back(&hs);
  for (const auto& [bid, v] : hinges_of) (void)model.bar_index(bid);  // (ids checked)
  for (const auto& [bid, sec] : opt.fiber_sections) {
    (void)model.bar_index(bid);
    require(!hinges_of.count(bid), "shell plasticity: a bar carries either a fiber section or hinges");
  }
  for (const auto& [bid, v] : opt.fiber_sections_along) {
    (void)model.bar_index(bid);
    require(!hinges_of.count(bid) && !opt.fiber_sections.count(bid),
            "shell plasticity: a bar carries one fiber description, or hinges");
    require(static_cast<int>(v.size()) == opt.fiber_points,
            "shell plasticity: fiber_sections_along needs one section per fiber point");
  }
  std::vector<NlBar> nlbars;
  bool nonsym = false;  // a non-symmetric element tangent (fiber bar on a bed): LU
  std::vector<int> nl_of(model.bars.size(), -1), cbar_of(model.bars.size(), -1);
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    const bool fib = opt.fiber_sections.count(b.id) > 0 || opt.fiber_sections_along.count(b.id) > 0;
    const bool hin = hinges_of.count(b.id) > 0;
    if (!fib && !hin) continue;
    require(b.kind == BarKind::Frame, "shell plasticity: a nonlinear bar must be a Frame bar");
    const bool m12 = hin && !ld;  // 12-DOF hinged bar: keeps arms, bed, taper, releases, span loads
    if (fib) {
      // force-based: releases as zero basic forces (shear releases have no
      // basic counterpart), rigid arms through T (linear) or the exact arm
      // kinematics (large displacements), taper through per-point sections;
      // the Winkler bed coupled in the element (linear kinematics)
      // (large displacements: kmx with span torques or a twist release
      // would need the exact torsion block in a corotated form)
      require(!ld || b.soil.kmx == 0 || (!b.releases[3] && !b.releases[9]),
              "shell plasticity (large displacements): a torsional bed on a fiber bar with a twist release is not supported");
      require(b.soil.kmx == 0 || model.tapered_of(b.id) == nullptr,
              "shell plasticity: a torsional bed on a tapered fiber bar is not supported");
      require(model.tapered_of(b.id) == nullptr || opt.fiber_sections_along.count(b.id),
              "shell plasticity: a tapered fiber bar needs fiber_sections_along (a section per point)");
    } else if (!m12) {
      // hinged bar in large displacements: series hinges on the elastic basic
      // system of the bar's local block (release-condensed or tapered)
    }
    const BarElement& el = asm_.elements[e];
    NlBar nb;
    nb.e = e;
    nb.L = el.L();
    Mat12 T = Mat12::Zero();
    for (int blk = 0; blk < 4; ++blk) T.block<3, 3>(3 * blk, 3 * blk) = el.R();
    nb.Blin = fiberel::basic_transform6(nb.L) * T;
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (int c = 0; c < 6; ++c) {
      nb.gdof[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
      nb.gdof[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
    }
    if (ld) {
      nb.ce = corot::make_corot_elem(model, asm_, e, CorotRotationMeasure::Log);
      nb.ce.bed_in_energy = false;  // (the bar adds its bed itself)
      // the kinematic round-off of the corotational deformations
      const double size = std::max({nb.ce.X1.cwiseAbs().maxCoeff(), nb.ce.X2.cwiseAbs().maxCoeff(), nb.ce.l0});
      // (well below the global force floor 64 eps |K| |X|, so an element
      // accepted at its floor never holds the global residual above it)
      nb.fe.v_noise_len = 4 * std::numeric_limits<double>::epsilon() * size;
      nb.fe.v_noise_ang = nb.fe.v_noise_len / nb.ce.l0;
    }
    {
      static const int RELB6[12] = {0, -1, -1, 5, 1, 3, 0, -1, -1, 5, 2, 4};  // local DOF -> basic component
      for (int c = 0; c < 12; ++c)
        if (b.releases[static_cast<std::size_t>(c)] && RELB6[c] >= 0) nb.rel6[static_cast<std::size_t>(RELB6[c])] = true;
    }
    if (m12) {
      nb.mode12 = true;
      const TaperedBar* tb = model.tapered_of(b.id);
      nb.kloc = tb ? tb->ke : el.k_local();
      nb.Tn = el.T();
      nb.ksoil = bar_ke_assembled(model, el) - el.T().transpose() * nb.kloc * el.T();
      for (const HingeSpec* hs : hinges_of.at(b.id)) {
        require(hs->end == 0 || hs->end == 1, "shell plasticity: hinge end must be 0 or 1");
        require(hs->axis >= 0 && hs->axis <= 2, "shell plasticity: hinge axis must be 0, 1 or 2");
        BasicHinge h;
        const int qrot = 6 * hs->end + 3 + hs->axis;
        if (hs->nm) {
          const NMHinge& nm = *hs->nm;
          require(nm.kN > 0 && nm.kM > 0 && !nm.surface.facets.empty(), "shell plasticity: invalid N-M hinge");
          h.m = 2;
          h.comp[0] = 6 * hs->end;  // spring coordinates: local DOFs
          h.comp[1] = qrot;
          const double sgnN = hs->end == 1 ? 1.0 : -1.0;  // N = sgnN s_x (tension positive)
          const Eigen::Matrix2d D = Eigen::Vector2d(sgnN, 1.0).asDiagonal();
          h.Ke << nm.kN, 0, 0, nm.kM;
          h.H = D * nm.H * D;
          for (const Vec3& f : nm.surface.facets) {
            h.n.push_back(D * Eigen::Vector2d(f(0), f(1)));
            h.r.push_back(f(2));
          }
        } else {
          const BilinearHinge& L = hs->law;
          require(L.k1 > 0 && L.My > 0 && L.k2 >= 0 && L.k2 < L.k1, "shell plasticity: invalid bilinear hinge");
          h.m = 1;
          h.comp[0] = qrot;
          h.Ke(0, 0) = L.k1;
          h.H(0, 0) = L.k2 == 0 ? 0.0 : L.k1 * L.k2 / (L.k1 - L.k2);
          h.n = {Eigen::Vector2d(1, 0), Eigen::Vector2d(-1, 0)};
          h.r = {L.My, L.My};
        }
        for (int a = 0; a < h.m; ++a) {
          require(!b.releases[static_cast<std::size_t>(h.comp[a])], "shell plasticity: a hinge on a released DOF");
          require(std::find(nb.q12.begin(), nb.q12.end(), h.comp[a]) == nb.q12.end(),
                  "shell plasticity: two hinges on one DOF");
          nb.q12.push_back(h.comp[a]);
        }
        nb.hs.push_back(h);
      }
    } else if (fib) {
      nb.fiber = true;
      if (!ld) nb.Blin = fiberel::basic_transform6(nb.L) * el.T();  // rigid arms in T
      nb.Tn = el.T();
      // axial / moment / torsion releases (q_c = 0) and SHEAR releases
      // (constant moment in that plane): the constraint basis of the element
      nb.fe.set_releases(b.releases);
      double GJ = model.material(b.material).G() * model.section(b.section).props.J;
      if (const TaperedBar* tb = model.tapered_of(b.id)) {
        // torsion of a tapered member: GJ/L from its exact kernel
        GJ = tb->ke(9, 9) * nb.L;
      }
      {
        // integration breakpoints at the span loads' DISCONTINUITIES: point
        // loads, and segment ends where the distributed intensity jumps (a
        // continuous chain of segments — the 64 pieces of a tapered bar's
        // self-weight — has none: the moment field stays smooth to its
        // second derivative there)
        const BarSpanLoads& sl = cl_bars.span[e];
        auto jump_at = [&](double x) {
          Vec3 qin = Vec3::Zero(), qout = Vec3::Zero(), min = Vec3::Zero(), mout = Vec3::Zero();
        require_linear_span(sl.segs, "shell plasticity (bars)");
          for (const DistSeg& g : sl.segs) {
            if (g.b == x) {
              qin += g.qb;
              min += g.mb;
            }
            if (g.a == x) {
              qout += g.qa;
              mout += g.ma;
            }
          }
          const double sc = std::max({qin.norm(), qout.norm(), min.norm(), mout.norm(), 1e-300});
          return (qin - qout).norm() > 1e-12 * sc || (min - mout).norm() > 1e-12 * sc;
        };
        for (const DistSeg& g : sl.segs) {
          if (jump_at(g.a)) nb.fe.breaks.push_back(g.a);
          if (jump_at(g.b)) nb.fe.breaks.push_back(g.b);
        }
        for (const PointAt& pt : sl.points) nb.fe.breaks.push_back(pt.a);
      }
      if (opt.fiber_sections_along.count(b.id)) {
        require(nb.fe.breaks.empty() || std::all_of(nb.fe.breaks.begin(), nb.fe.breaks.end(), [&](double c) {
                  return c <= 1e-12 * nb.L || c >= nb.L * (1 - 1e-12);
                }),
                "shell plasticity: a fiber bar with sections along takes no interior span-load breakpoint");
        const std::vector<HystFiberSection>& along = opt.fiber_sections_along.at(b.id);
        std::vector<const HystFiberSection*> ptrs;
        for (const HystFiberSection& hs : along) ptrs.push_back(&hs);
        nb.fe.set_sections(ptrs);
        nb.fe.init(along.front(), nb.L, GJ, opt.fiber_points, 1e-6, 30, 1e-13);
      } else {
        nb.fe.init(opt.fiber_sections.at(b.id), nb.L, GJ, opt.fiber_points, 1e-6, 30, 1e-13);
      }
      nb.GJ = GJ;
      if (ld) {
        // the fiber bar's translational bed is coupled in the element; the
        // twist bed acts on phi = e1_0 . theta as the EXACT torsion block of
        // the beam on foundation minus its GJ part (carried by the basic
        // twist): the linear limit is the linear-kinematics bar's exact block
        nb.ce.kbed.setZero();
        nb.ce.Kphi.setZero();
        if (b.soil.kmx != 0) {
          BefBar bb;
          bb.L = nb.L;
          bb.EA = bb.EIy = bb.EIz = 1.0;
          bb.GJ = GJ;
          bb.zones = {{nb.L, 0, 0, 0, b.soil.kmx}};
          const Mat12 kt = bef_ke(bb);
          nb.ce.Kphi << kt(3, 3), kt(3, 9), kt(9, 3), kt(9, 9);
          nb.ce.Kphi -= GJ / nb.L * (Eigen::Matrix2d() << 1, -1, -1, 1).finished();
        }
      }
      if (b.soil.kfx != 0 || b.soil.kfy != 0 || b.soil.kfz != 0) {
        nb.fe.init_bed(Vec3(b.soil.kfx, b.soil.kfy, b.soil.kfz));
        nonsym = true;
      }
      if (b.soil.kmx != 0 && !ld) {
        // torsional bed: GJ th'' - kmx th = -m, elastic and uncoupled from
        // the fibers — its EXACT block (bef.h transfer functions)
        nb.tor_ext = nb.fe.tor_external = true;
        nb.fe.kb(5, 5) = 0;
        BefBar bb;
        bb.L = nb.L;
        bb.EA = bb.EIy = bb.EIz = 1.0;  // (bending / axial unused)
        bb.GJ = GJ;
        bb.zones = {{nb.L, 0, 0, 0, b.soil.kmx}};
        const Mat12 kt = bef_ke(bb);
        nb.ktor << kt(3, 3), kt(3, 9), kt(9, 3), kt(9, 9);
        BarSpanLoads tor_only;
        require_linear_span(cl_bars.span[e].segs, "shell plasticity (bars)");
        for (const DistSeg& g : cl_bars.span[e].segs) {
          DistSeg t{g.a, g.b, Vec3::Zero(), Vec3::Zero()};
          t.ma = Vec3(g.ma.x(), 0, 0);
          t.mb = Vec3(g.mb.x(), 0, 0);
          if (t.ma.x() != 0 || t.mb.x() != 0) tor_only.segs.push_back(t);
        }
        for (const PointAt& pt : cl_bars.span[e].points)
          if (pt.M.x() != 0) tor_only.points.push_back({pt.a, Vec3::Zero(), Vec3(pt.M.x(), 0, 0)});
        for (const KinkAt& kk : cl_bars.span[e].kinks)
          require(kk.dthx == 0, "shell plasticity: a twist dislocation on a fiber bar on a torsional bed is not supported");
        const Vec12 ft = bef_f0(bb, tor_only);
        nb.f0tor << ft(3), ft(9);
        // release of the twist at one end: static condensation
        const bool ri = b.releases[3], rj = b.releases[9];
        if (ri || rj) {
          const int r = ri ? 0 : 1, k = 1 - r;
          Eigen::Matrix2d kc = Eigen::Matrix2d::Zero();
          Eigen::Vector2d fc = Eigen::Vector2d::Zero();
          if (!(ri && rj)) {
            kc(k, k) = nb.ktor(k, k) - nb.ktor(k, r) * nb.ktor(r, k) / nb.ktor(r, r);
            fc(k) = nb.f0tor(k) - nb.ktor(k, r) / nb.ktor(r, r) * nb.f0tor(r);
          }
          nb.ktor = kc;
          nb.f0tor = fc;
        }
        nonsym = nonsym || false;
      }
    } else {
      {
        const TaperedBar* tb = model.tapered_of(b.id);
        nb.kb = fiberel::elastic_basic6(tb ? tb->ke : el.k_local(), nb.L);
      }
      for (const HingeSpec* hs : hinges_of.at(b.id)) {
        BasicHinge h;
        require(hs->end == 0 || hs->end == 1, "shell plasticity: hinge end must be 0 or 1");
        int rc = 0;
        double rsign = 1;
        if (hs->axis == 1)
          rc = hs->end == 0 ? 1 : 2;
        else if (hs->axis == 2)
          rc = hs->end == 0 ? 3 : 4;
        else {
          rc = 5;  // twist = rx_j - rx_i: a hinge at end j reduces it
          rsign = hs->end == 0 ? 1 : -1;
        }
        if (hs->nm) {
          const NMHinge& nm = *hs->nm;
          require(nm.kN > 0 && nm.kM > 0 && !nm.surface.facets.empty(), "shell plasticity: invalid N-M hinge");
          h.m = 2;
          h.comp[0] = 0;  // elongation, conjugate to N (tension positive)
          h.comp[1] = rc;
          h.sign[1] = rsign;
          h.Ke << nm.kN, 0, 0, nm.kM;
          h.H = nm.H;
          for (const Vec3& f : nm.surface.facets) {
            h.n.emplace_back(f(0), f(1));
            h.r.push_back(f(2));
          }
        } else {
          const BilinearHinge& L = hs->law;
          require(L.k1 > 0 && L.My > 0 && L.k2 >= 0 && L.k2 < L.k1, "shell plasticity: invalid bilinear hinge");
          h.m = 1;
          h.comp[0] = rc;
          h.sign[0] = rsign;
          h.Ke(0, 0) = L.k1;
          h.H(0, 0) = L.k2 == 0 ? 0.0 : L.k1 * L.k2 / (L.k1 - L.k2);
          h.n = {Eigen::Vector2d(1, 0), Eigen::Vector2d(-1, 0)};
          h.r = {L.My, L.My};
        }
        for (int a = 0; a < h.m; ++a)
          require(!nb.rel6[static_cast<std::size_t>(h.comp[a])], "shell plasticity: a hinge on a released DOF");
        nb.hs.push_back(h);
      }
    }
    nl_of[e] = static_cast<int>(nlbars.size());
    nlbars.push_back(std::move(nb));
  }
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    require(b.kind != BarKind::Cable, "shell plasticity: cables belong to solve_cable_static");
    if (nl_of[e] >= 0) continue;
    if (ld) {
      // rigid arms, releases, taper and the Winkler bed (ground springs on
      // the current axis) are in the corotational element
      cbar_of[e] = static_cast<int>(cbars.size());
      cbars.push_back(corot::make_corot_elem(model, asm_, e, CorotRotationMeasure::Log));
      continue;
    }
    const Mat12 k = bar_ke_assembled(model, asm_.elements[e]);
    int g[12];
    bar_gdofs(b, g);
    for (int a = 0; a < 12; ++a)
      for (int c = 0; c < 12; ++c)
        if (g[a] >= 0 && g[c] >= 0 && k(a, c) != 0.0) lin_trips.emplace_back(g[a], g[c], k(a, c));
  }
  for (std::size_t f = 0; f < asm_.spring.size(); ++f)
    if (asm_.spring[f] != 0.0 && d.free_ix[f] >= 0) lin_trips.emplace_back(d.free_ix[f], d.free_ix[f], asm_.spring[f]);
  for (const auto& [ni, K6] : asm_.spring6)
    for (int a = 0; a < 6; ++a)
      for (int c = 0; c < 6; ++c) {
        const int fa = d.free_ix[static_cast<std::size_t>(d.full(ni, a))], fc = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
        if (fa >= 0 && fc >= 0 && K6(a, c) != 0.0) lin_trips.emplace_back(fa, fc, K6(a, c));
      }
  for (const LinkElem& le : asm_.links)
    for (int a = 0; a < 12; ++a) {
      const int fa = d.free_of(a < 6 ? le.ni : le.nj, a % 6);
      if (fa < 0) continue;
      for (int c = 0; c < 12; ++c) {
        const int fc = d.free_of(c < 6 ? le.ni : le.nj, c % 6);
        if (fc >= 0 && le.K(a, c) != 0.0) lin_trips.emplace_back(fa, fc, le.K(a, c));
      }
    }
  SpMat K_lin(d.n_free, d.n_free);
  K_lin.setFromTriplets(lin_trips.begin(), lin_trips.end());
  const CaseLoads cl = lower_case(model, asm_, model.load_case(case_id));
  require(cl.prescribed.empty(), "shell plasticity: prescribed displacements are not supported");
  for (const ShellRecovery& sr : cl.shell_recovery)
    require(sr.eps0.isZero() && sr.kappa0.isZero() && sr.dstar.isZero(),
            "shell plasticity: thermal loads and imposed shell fields are not supported");
  const VecX F = cl.rhs;
  require(F.norm() > 0, "shell plasticity: the load case is empty");
  for (NlBar& nb : nlbars) {
    if (nb.mode12) {
      const BarSpanLoads& sp = cl.span[nb.e];
      require(sp.kinks.empty() && sp.eps0 == 0 && sp.kappa_y == 0 && sp.kappa_z == 0 && sp.deps0 == 0 &&
                  sp.dkappa_y == 0 && sp.dkappa_z == 0 && sp.kappa_x == 0 && sp.dkappa_x == 0 && !has_parabolic_strain(sp),
              "shell plasticity: imposed strains / dislocations on a hinged bar are not supported");
      nb.f0c = asm_.elements[nb.e].condense_f0(cl.f0[nb.e]);
      continue;
    }
    if (nb.fiber) {
      // imposed strains / dislocations of the case (both kinematics: local)
      const BarSpanLoads& rec = cl.span[nb.e];
      Vec6 vk;
      fiberel::imposed_fields(rec, nb.L, nb.fe.xg, nb.e0_unit, vk);
      require(rec.kinks.empty() || !nb.fe.has_bed(),
              "shell plasticity: a dislocation on a fiber bar on a Winkler bed is not supported");
      nb.vshift_unit = vk;
      if (nb.tor_ext) nb.vshift_unit(5) = 0;  // (the twist lives in the exact torsion block)
    }
    if (nb.fiber && !ld && cl.f0[nb.e].cwiseAbs().maxCoeff() > 0.0) {
      fiberel::SpanParticular spp = fiberel::span_particular(cl.span[nb.e], nb.L, nb.fe.xg, &nb.fe.incl);
      if (nb.tor_ext) {  // torques: the exact torsion block
        spp.p_ss(3) = spp.p_ss(9) = 0;
        spp.tor_int = 0;
      }
      fiberel::release_statics(spp, model.bars[nb.e].releases, nb.L, nb.fe.xg);
      nb.sp_unit = spp.sp;
      nb.p_ss = spp.p_ss;
      if (!nb.tor_ext) nb.vshift_unit(5) += spp.tor_int / nb.GJ;
      nb.f0c = asm_.elements[nb.e].condense_f0(cl.f0[nb.e]);
      nb.span = true;
      continue;
    }
    if (!ld) require(cl.f0[nb.e].cwiseAbs().maxCoeff() == 0.0, "shell plasticity: span loads on a nonlinear bar are not supported");
  }
  // large displacements: dead nodal forces + spatially fixed nodal moments
  VecX Fforce = F;
  std::vector<std::pair<int, Vec3>> moments;
  std::vector<corot::CorotBarLoad> cbar_loads;  // span loads of the corotational elastic bars
  VecX bar_full = VecX::Zero(cl.nodal_full.size());  // the bars' equivalent loads, all DOFs
  if (ld) {
    // span loads on bars leave the case vector: they are re-lowered on the
    // CURRENT corotated frame (elastic bars: the equivalent loads of
    // solve_corotational; nonlinear bars: basic-system reactions pushed back
    // + the particular field / deformation shift inside the element)
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      if (cl.f0[e].cwiseAbs().maxCoeff() == 0.0) continue;
      const BarSpanLoads& sp = cl.span[e];
      // (a fiber bar carries them in its own frame: imposed_fields)
      require((nl_of[e] >= 0 && nlbars[static_cast<std::size_t>(nl_of[e])].fiber) ||
                  (sp.kinks.empty() && sp.eps0 == 0 && sp.kappa_y == 0 && sp.kappa_z == 0 && sp.deps0 == 0 &&
                   sp.dkappa_y == 0 && sp.dkappa_z == 0 && sp.kappa_x == 0 && sp.dkappa_x == 0 && !has_parabolic_strain(sp)),
              "shell plasticity (large displacements): imposed strains / dislocations on elastic or hinged bars are not supported");
      const BarElement& el = asm_.elements[e];
      const Vec12 fg = el.T().transpose() * el.condense_f0(cl.f0[e]);
      int g[12];
      bar_gdofs(model.bars[e], g);
      const int nn[2] = {model.node_index(model.bars[e].ni), model.node_index(model.bars[e].nj)};
      for (int a = 0; a < 12; ++a) {
        if (g[a] >= 0) Fforce(g[a]) -= fg(a);
        bar_full(d.full(nn[a / 6], a % 6)) += fg(a);  // nodal_full carries them too
      }
    }
    // sw: a self-weight (l = the gravity as a global UDL record): its weight
    // profile rho A(x) — exact for a tapered bar's variable section
    auto add_bar_load = [&](std::size_t e, const Load& l, bool sw) {
      require(l.frame != LoadFrame::GlobalProjected,
              "shell plasticity (large displacements): projected (PRO) bar loads are not supported");
      if (nl_of[e] < 0) {
        const std::size_t ce = static_cast<std::size_t>(cbar_of[e]);
        cbar_loads.push_back(sw ? corot::make_corot_self_weight(model, asm_, e, ce, l.F)
                                : corot::make_corot_bar_load(model, asm_, e, ce, l));
        return;
      }
      NlBar& nb = nlbars[static_cast<std::size_t>(nl_of[e])];
      NlBar::LdLoad L;
      L.w = sw ? corot::make_corot_self_weight(model, asm_, e, 0, l.F) : corot::make_corot_bar_load(model, asm_, e, 0, l);
      L.ss = L.w;
      // the weight per unit length along the bar (self-weight)
      const Bar& bw = model.bars[e];
      const TaperedBar* tbw = model.tapered_of(bw.id);
      const double rho = model.material(bw.material).rho;
      const double A0 = model.section(bw.section).props.A;
      auto weight = [&](double x) { return rho * (tbw && tbw->props ? tbw->props(x / nb.L).A : A0); };
      const double Lb = nb.L;
      const fiberel::Mat612 a6 = fiberel::basic_transform6(Lb);
      const Mat6 aat = a6 * a6.transpose();
      const Eigen::CompleteOrthogonalDecomposition<Mat6> kb_cod(nb.kb);
      for (int k = 0; k < 9; ++k) {
        // a local intensity component that stays zero at every state (a
        // dead load's zero vector, a follower's zero component, the unused
        // F2 of a uniform / point load) needs no field
        const Vec3& gk = k < 3 ? l.F : (k < 6 ? l.M : l.F2);
        if ((k >= 6 && l.kind != Load::Kind::BarLinear) || (L.w.follower ? gk(k % 3) == 0.0 : gk.isZero())) continue;
        // the unit local intensity component k as a span record
        Vec3 F = Vec3::Zero(), M = Vec3::Zero(), F2 = Vec3::Zero();
        (k < 3 ? F : k < 6 ? M : F2)(k % 3) = 1.0;
        const std::vector<double> xs = nb.fiber ? nb.fe.xg : std::vector<double>{};
        const std::vector<char>* inc = nb.fiber ? &nb.fe.incl : nullptr;
        fiberel::SpanParticular spp;
        if (sw) {
          // the weight profile rho A(x) e_c, integrated by Gauss-Legendre
          // (12 points per span: exact for a polynomial A of degree <= 22)
          std::vector<double> gx, gw;
          fiberel::FiberElement3D::gauss_legendre01(12, gx, gw);
          auto loads_upto = [&](double x, double x0, Vec3& Fq, Vec3& Mq, bool = false) {
            const double c = std::min(x, Lb);
            double f0 = 0, f1 = 0;
            for (std::size_t i = 0; i < gx.size() && c > 0; ++i) {
              const double t = gx[i] * c;
              f0 += gw[i] * c * weight(t);
              f1 += gw[i] * c * (t - x0) * weight(t);
            }
            Fq = f0 * F;
            Mq = f1 * Vec3(0, -F.z(), F.y());
          };
          spp = fiberel::span_particular_core(Lb, xs, inc, loads_upto);
        } else {
          BarSpanLoads rec;
          const double a = std::clamp(l.relative ? l.a * Lb : l.a, 0.0, Lb);
          if (l.kind == Load::Kind::BarPoint) {
            rec.points.push_back({a, F, M});
          } else {
            const double bb = std::clamp(l.b < 0 ? Lb : (l.relative ? l.b * Lb : l.b), 0.0, Lb);
            DistSeg sg{a, bb, F, l.kind == Load::Kind::BarUdl ? F : F2};
            sg.ma = sg.mb = M;
            rec.segs.push_back(sg);
          }
          spp = fiberel::span_particular(rec, Lb, xs, inc);
        }
        // fiber bars: the particular field moved off the released actions
        // (their sections see the total); hinged bars keep the simply
        // supported basic system, the released statics going into q0
        if (nb.fiber) fiberel::release_statics(spp, model.bars[e].releases, Lb, nb.fe.xg);
        L.ss.basis.col(k) = spp.p_ss;
        if (nb.fiber) {
          require(spp.tor_int == 0 || model.bars[e].soil.kmx == 0,
                  "shell plasticity (large displacements): span torques on a fiber bar on a torsional bed are not supported");
          L.S[static_cast<std::size_t>(k)] = spp.sp;
          L.V0(5, k) = spp.tor_int / nb.GJ;  // elastic twist of the particular torque
        } else {
          // fixed-end basic forces: a6^T q0 = -f0 - p_ss (exact: the
          // right-hand side is a self-equilibrated basic-system field)
          const Vec6 q0 = aat.ldlt().solve(a6 * (-L.w.basis.col(k) - spp.p_ss));
          // q0 = q0r + q0n: q0r in the range of kb, the elastic shift
          // kb v0 = -q0r (least norm: a release — moment or SHEAR, whose null
          // direction is not a coordinate axis — leaves free directions);
          // q0n in its null space, the part the released member carries by
          // equilibrium alone (it passes through the hinges as is)
          const Vec6 v0 = -kb_cod.solve(q0);
          const Vec6 q0r = -nb.kb * v0;
          L.V0.col(k) = v0;
          L.Qn.col(k) = q0 - q0r;
        }
      }
      nb.ldl.push_back(std::move(L));
    };
    for (const Load& l : model.load_case(case_id).loads) {
      switch (l.kind) {
        case Load::Kind::BarUdl:
        case Load::Kind::BarLinear:
        case Load::Kind::BarPoint:
          add_bar_load(static_cast<std::size_t>(model.bar_index(l.bar)), l, false);
          break;
        case Load::Kind::BarParabolic:
          fail("shell plasticity: parabolic bar loads are not supported here (linear span loads only)");
        case Load::Kind::SelfWeight:
          for (std::size_t e = 0; e < model.bars.size(); ++e)
            add_bar_load(e, Load::udl(model.bars[e].id, l.g, LoadFrame::Global), true);
          break;
        default:
          break;
      }
    }
    for (int n = 0; n < d.n_nodes; ++n) {
      Vec3 M = Vec3::Zero();
      for (int c = 0; c < 3; ++c) {
        M(c) = cl.nodal_full(d.full(n, 3 + c)) - bar_full(d.full(n, 3 + c));
        const int fi = d.free_ix[static_cast<std::size_t>(d.full(n, 3 + c))];
        if (fi >= 0) Fforce(fi) = 0.0;
      }
      if (M.squaredNorm() > 0) moments.emplace_back(n, M);
    }
  }
  auto rot_dofs = [&](int n, int fx[3]) {
    for (int c = 0; c < 3; ++c) fx[c] = d.free_ix[static_cast<std::size_t>(d.full(n, 3 + c))];
  };
  // external force at u, lambda, and its derivative (load stiffness) per unit lambda
  auto external = [&](const VecX& u, double lam, std::vector<Triplet>* Kload) {
    VecX Fx = lam * Fforce;
    for (const corot::CorotBarLoad& bl : cbar_loads) {
      const corot::CorotElem& ce = cbars[bl.ei];
      int g[12];
      std::array<ad::Grad<12>, 12> x;
      for (int c = 0; c < 6; ++c) {
        g[c] = d.free_ix[static_cast<std::size_t>(d.full(ce.ni, c))];
        g[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(ce.nj, c))];
      }
      for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(g[i] >= 0 ? u(g[i]) : 0.0, i);
      const std::array<ad::Grad<12>, 12> f = corot::bar_load_vector(ce, bl, x);
      for (int r = 0; r < 12; ++r) {
        if (g[r] < 0) continue;
        Fx(g[r]) += lam * f[static_cast<std::size_t>(r)].v;
        if (Kload)
          for (int c = 0; c < 12; ++c)
            if (g[c] >= 0) Kload->emplace_back(g[r], g[c], lam * f[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)]);
      }
    }
    for (const auto& [n, M] : moments) {
      int fx[3];
      rot_dofs(n, fx);
      corot::V3t<ad::Grad<3>> th;
      for (int c = 0; c < 3; ++c) th[static_cast<std::size_t>(c)] = ad::Grad<3>::seed(fx[c] >= 0 ? u(fx[c]) : 0.0, c);
      const corot::V3t<double> m{M(0), M(1), M(2)};
      const corot::V3t<ad::Grad<3>> g = corot::jt_times(th, m);
      for (int r = 0; r < 3; ++r) {
        if (fx[r] < 0) continue;
        Fx(fx[r]) += lam * g[static_cast<std::size_t>(r)].v;
        if (Kload)
          for (int c = 0; c < 3; ++c)
            if (fx[c] >= 0) Kload->emplace_back(fx[r], fx[c], lam * g[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)]);
      }
    }
    return Fx;
  };

  struct Elem {
    ShellNlOps ops;
    corot::CorotShell cs;
    int gdof[18];
    double E, nu, fy;
    bool hill = false;
    Mat5 P = Mat5::Zero();
    ThicknessRule rule;  // its points through the thickness (per ply for a laminate)
  };
  double force_floor = 0;  // round-off floor of the corotational elements
  // corotational bars: chord and frame are differences of positions
  // (cancellation ~ eps |X|): stretch noise eps |X|, angle noise eps |X| / l0
  auto bar_floor = [&](const corot::CorotElem& ce) {
    const double size = std::max({ce.X1.cwiseAbs().maxCoeff(), ce.X2.cwiseAbs().maxCoeff(), ce.l0});
    const double krot = ce.Kbar.bottomRightCorner<6, 6>().cwiseAbs().rowwise().sum().maxCoeff();
    force_floor = std::max(force_floor, 64 * std::numeric_limits<double>::epsilon() * size *
                                            std::max({std::fabs(ce.Kbar(0, 0)), krot / ce.l0, krot / (ce.l0 * ce.l0)}));
  };
  if (ld) {
    for (const corot::CorotElem& ce : cbars) bar_floor(ce);
    for (const NlBar& nb : nlbars) bar_floor(nb.ce);
  }
  double kmax = 0;         // largest element stiffness row sum (round-off floor of K u)
  for (const ShellElem& se : asm_.shells) kmax = std::max(kmax, se.el.ke().cwiseAbs().rowwise().sum().maxCoeff());
  for (std::size_t e = 0; e < model.bars.size(); ++e)
    if (model.bars[e].kind != BarKind::Cable)
      kmax = std::max(kmax, asm_.elements[e].ke().cwiseAbs().rowwise().sum().maxCoeff());
  for (const auto& [mid, h] : opt.hill) (void)model.material(mid);
  std::unordered_set<Id> hill_shells;  // the shells of a Hill material
  if (!opt.hill_axis1.empty())
    for (const ShellTri& sh : model.shells)
      if (opt.hill.count(sh.material)) hill_shells.insert(sh.id);
  for (const auto& [sid, a] : opt.hill_axis1) {
    const bool found = hill_shells.count(sid) > 0;
    require(found, "shell plasticity: hill_axis1 names shell " + std::to_string(sid) + ", not a Hill shell");
  }
  std::vector<Elem> els;
  for (const ShellElem& se : asm_.shells) {
    Elem el;
    el.ops = se.el.nl_ops();
    for (int k = 0; k < 3; ++k)
      for (int c = 0; c < 6; ++c)
        el.gdof[6 * k + c] = d.free_ix[static_cast<std::size_t>(d.full(se.nodes[static_cast<std::size_t>(k)], c))];
    const Material& mat = model.material(model.shells[se.index].material);
    require(mat.fy > 0, "shell plasticity: the shell material needs fy > 0");
    el.E = mat.E;
    el.nu = mat.nu;
    el.fy = mat.fy;
    if (auto hi = opt.hill.find(model.shells[se.index].material); hi != opt.hill.end()) {
      // material axis 1: axis1 projected on the facet plane, its angle from the local x
      const Mat3 R = el.ops.T.block<3, 3>(0, 0);  // rows: local axes in global
      double angle = 0;
      const auto ax = opt.hill_axis1.find(model.shells[se.index].id);
      const Vec3 axis = ax != opt.hill_axis1.end() ? ax->second : hi->second.axis1;
      if (axis.norm() > 0) {
        const Vec3 a = R * axis;
        require(std::hypot(a(0), a(1)) > 1e-6 * a.norm(), "Hill criterion: axis1 is normal to a facet");
        angle = std::atan2(a(1), a(0));
      }
      el.hill = true;
      el.P = hill_projector(hi->second, angle);
    }
    el.rule = el.ops.plies.empty() ? lobatto4(opt.layers) : ply_rule(el.ops, opt.layers);
    require(!el.ops.laminate || !opt.plastic_transverse_shear,
            "shell plasticity: a laminate takes the transverse shear elastic (plastic_transverse_shear = false: "
            "the parabolic shear profile is not the plies' own)");
    if (ld) {
      el.cs = corot::make_corot_shell(se);
      double size = 0;
      for (const Vec3& x : el.cs.X) size = std::max(size, x.cwiseAbs().maxCoeff());
      force_floor = std::max(force_floor, 64 * std::numeric_limits<double>::epsilon() *
                                              el.cs.Kl.cwiseAbs().rowwise().sum().maxCoeff() * size);
    }
    els.push_back(el);
  }
  const ThicknessRule rule = lobatto4(opt.layers);
  std::vector<ShellState> committed(els.size());
  for (std::size_t e = 0; e < els.size(); ++e)
    committed[e].assign(els[e].ops.w.size(), std::vector<ShellLayerState>(els[e].rule.zeta.size()));
  std::vector<ShellState> trial = committed;

  // internal force and tangent at u (trial states from the committed ones)
  std::vector<Vec6>* out_q = nullptr;  // basic forces of the nonlinear bars, when requested
  std::vector<Vec12>* out_end = nullptr;  // local end forces of the 12-DOF hinged bars
  double lam_cur = 0;                   // load factor of the current evaluation (span loads)
  double elem_noise = 0;  // force noise of the fiber elements' section solves (last evaluation)
  auto evaluate = [&](const VecX& u, VecX& fint, SpMat* K) {
    elem_noise = 0;
    fint = VecX::Zero(d.n_free);
    std::vector<Triplet> trips;
    if (K) trips.reserve(els.size() * 324 + 16);
    for (std::size_t e = 0; e < els.size(); ++e) {
      const Elem& el = els[e];
      const ShellNlOps& o = el.ops;
      Vec18 ug;
      for (int a = 0; a < 18; ++a) ug(a) = el.gdof[a] >= 0 ? u(el.gdof[a]) : 0.0;
      Vec18 fl;
      Mat18 kl;
      Vec18 fg;
      Mat18 kg;
      if (!ld) {
        shell_local_response(o, o.T * ug, committed[e], trial[e], el.E, el.nu, el.fy, opt.hardening, el.rule,
                             opt.plastic_transverse_shear, fl, K ? &kl : nullptr, el.hill ? &el.P : nullptr);
        fg = o.T.transpose() * fl;
        if (K) kg = o.T.transpose() * kl * o.T;
      } else {
        // d(q) and J = dd/dq (first-order AD), the local response at d, then
        // the pull-back J^T f, J^T K J + d2(f . d)/dq2 (second-order AD)
        std::array<ad::Grad<18>, 18> qg;
        for (int i = 0; i < 18; ++i) qg[static_cast<std::size_t>(i)] = ad::Grad<18>::seed(ug(i), i);
        const std::array<ad::Grad<18>, 18> dg = el.cs.deform_t(qg);
        Vec18 dl;
        Mat18 J;
        for (int i = 0; i < 18; ++i) {
          dl(i) = dg[static_cast<std::size_t>(i)].v;
          for (int j = 0; j < 18; ++j) J(i, j) = dg[static_cast<std::size_t>(i)].g[static_cast<std::size_t>(j)];
        }
        shell_local_response(o, dl, committed[e], trial[e], el.E, el.nu, el.fy, opt.hardening, el.rule,
                             opt.plastic_transverse_shear, fl, K ? &kl : nullptr, el.hill ? &el.P : nullptr);
        fg = J.transpose() * fl;
        if (K) {
          std::array<ad::Jet<18>, 18> qj;
          for (int i = 0; i < 18; ++i) qj[static_cast<std::size_t>(i)] = ad::Jet<18>::seed(ug(i), i);
          const std::array<ad::Jet<18>, 18> dj = el.cs.deform_t(qj);
          ad::Jet<18> phi = dj[0] * fl(0);
          for (int i = 1; i < 18; ++i) phi = phi + dj[static_cast<std::size_t>(i)] * fl(i);
          kg = J.transpose() * kl * J;
          for (int i = 0; i < 18; ++i)
            for (int j = 0; j < 18; ++j) kg(i, j) += phi.h[static_cast<std::size_t>(18 * i + j)];
          kg = 0.5 * (kg + kg.transpose()).eval();
        }
      }
      for (int a = 0; a < 18; ++a)
        if (el.gdof[a] >= 0) fint(el.gdof[a]) += fg(a);
      if (K) {
        for (int a = 0; a < 18; ++a)
          for (int b = 0; b < 18; ++b)
            if (el.gdof[a] >= 0 && el.gdof[b] >= 0 && kg(a, b) != 0.0) trips.emplace_back(el.gdof[a], el.gdof[b], kg(a, b));
      }
    }
    // corotational elastic bars (energy formulation, exact AD tangent)
    for (const corot::CorotElem& ce : cbars) {
      Eigen::Matrix<double, 12, 1> q = Eigen::Matrix<double, 12, 1>::Zero();
      int g[12];
      for (int c = 0; c < 6; ++c) {
        g[c] = d.free_ix[static_cast<std::size_t>(d.full(ce.ni, c))];
        g[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(ce.nj, c))];
      }
      for (int a = 0; a < 12; ++a)
        if (g[a] >= 0) q(a) = u(g[a]);
      const auto fe = ce.grad(q);
      for (int a = 0; a < 12; ++a)
        if (g[a] >= 0) fint(g[a]) += fe(a);
      if (K) {
        const auto ke = ce.hess(q);
        for (int a = 0; a < 12; ++a)
          for (int b = 0; b < 12; ++b)
            if (g[a] >= 0 && g[b] >= 0) trips.emplace_back(g[a], g[b], 0.5 * (ke(a, b) + ke(b, a)));
      }
    }
    // nonlinear bars: basic response pulled back through B (linear) or the
    // corotational Jacobian (exact AD, plus the geometric Hessian of q . v)
    for (NlBar& nb : nlbars) {
      Vec12 ue;
      for (int a = 0; a < 12; ++a) ue(a) = nb.gdof[a] >= 0 ? u(nb.gdof[a]) : 0.0;
      if (nb.mode12) {
        Vec12 fe, p;
        Mat12 ke;
        nb.response12(ue, lam_cur, fe, ke, p);
        for (int a = 0; a < 12; ++a) {
          if (nb.gdof[a] < 0) continue;
          fint(nb.gdof[a]) += fe(a);
          if (K)
            for (int b = 0; b < 12; ++b)
              if (nb.gdof[b] >= 0) trips.emplace_back(nb.gdof[a], nb.gdof[b], 0.5 * (ke(a, b) + ke(b, a)));
        }
        if (out_q) (*out_q)[nb.e] << p(6), p(4), p(10), p(5), p(11), p(9);
        if (out_end) (*out_end)[nb.e] = p;
        continue;
      }
      if (nb.fiber) {
        // imposed strains and basic deformation shifts at this load level
        for (std::size_t k = 0; k < nb.e0_unit.size(); ++k) nb.fe.e0[k] = lam_cur * nb.e0_unit[k];
        nb.fe.vshift = lam_cur * nb.vshift_unit;
      }
      // the exact torsion block of a fiber bar on a torsional bed (linear)
      auto add_tor = [&](Vec12& fe_g, Mat12* ke_g) {
        if (!nb.tor_ext) return;
        const Vec12 ul = nb.Tn * ue;
        const Eigen::Vector2d t = nb.ktor * Eigen::Vector2d(ul(3), ul(9)) - lam_cur * nb.f0tor;
        nb.tor_cur.setZero();
        nb.tor_cur(3) = t(0);
        nb.tor_cur(9) = t(1);
        fe_g += nb.Tn.transpose() * nb.tor_cur;
        if (ke_g) {
          Mat12 kl = Mat12::Zero();
          kl(3, 3) = nb.ktor(0, 0);
          kl(3, 9) = nb.ktor(0, 1);
          kl(9, 3) = nb.ktor(1, 0);
          kl(9, 9) = nb.ktor(1, 1);
          *ke_g += nb.Tn.transpose() * kl * nb.Tn;
        }
      };
      if (!ld && nb.fiber && nb.fe.has_bed()) {
        // fiber bar on a Winkler bed (linear kinematics): the coupled element
        if (nb.span)
          for (std::size_t k = 0; k < nb.sp_unit.size(); ++k) nb.fe.sp[k] = lam_cur * nb.sp_unit[k];
        Vec12 fl;
        Mat12 Kl;
        nb.fe.determine_bed(nb.Tn * ue, fl, Kl);
        const Vec6 qb = nb.fe.q6();
        {
          // force noise: the section solves' floor through kb, the bed
          // pressures' round-off through the reactions
          Vec6 vs = Vec6::Zero();
          vs.head<5>() = nb.fe.section_noise();
          const Vec12 nl = fiberel::basic_transform6(nb.L).transpose().cwiseAbs() * (nb.fe.kb * vs).cwiseAbs() +
                           nb.fe.Rop.cwiseAbs() * nb.fe.bed_pnoise;
          elem_noise = std::max(elem_noise, (nb.Tn.transpose().cwiseAbs() * nl).sum());
        }
        nb.pss_cur = fl - fiberel::basic_transform6(nb.L).transpose() * qb;  // bed reactions
        Vec12 fe = nb.Tn.transpose() * fl;
        if (nb.span) {
          fe += lam_cur * nb.Tn.transpose() * (nb.p_ss + nb.f0c);  // the case vector carries T^t f0c
          nb.pss_cur += lam_cur * nb.p_ss;
        }
        Mat12 ke = nb.Tn.transpose() * Kl * nb.Tn;
        add_tor(fe, &ke);
        for (int a = 0; a < 12; ++a)
          if (nb.gdof[a] >= 0) fint(nb.gdof[a]) += fe(a);
        if (K) {
          for (int a = 0; a < 12; ++a)
            for (int b = 0; b < 12; ++b)
              if (nb.gdof[a] >= 0 && nb.gdof[b] >= 0 && ke(a, b) != 0.0) trips.emplace_back(nb.gdof[a], nb.gdof[b], ke(a, b));
        }
        if (out_q) (*out_q)[nb.e] = qb;
        continue;
      }
      Vec6 v, q;
      Mat6 Kq;
      fiberel::Mat612 B;
      static const int MAPV[6][2] = {{0, 0}, {2, 0}, {5, 0}, {3, 0}, {6, 0}, {4, 1}};  // v_k = d[i] (- d[1] for twist)
      if (!ld) {
        B = nb.Blin;
        v = B * ue;
      } else {
        std::array<ad::Grad<12>, 12> x;
        for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(ue(i), i);
        const std::array<ad::Grad<12>, 7> dd = nb.ce.deform_t(x);
        for (int k = 0; k < 6; ++k) {
          const auto& di = dd[static_cast<std::size_t>(MAPV[k][0])];
          v(k) = di.v - (MAPV[k][1] ? dd[1].v : 0.0);
          for (int j = 0; j < 12; ++j)
            B(k, j) = di.g[static_cast<std::size_t>(j)] - (MAPV[k][1] ? dd[1].g[static_cast<std::size_t>(j)] : 0.0);
        }
      }
      if (nb.span)
        for (std::size_t k = 0; k < nb.sp_unit.size(); ++k) nb.fe.sp[k] = lam_cur * nb.sp_unit[k];
      // large displacements, span loads: local intensities w_l(u) on the
      // current corotated frame and their Jacobians (exact AD)
      const std::size_t nl = ld ? nb.ldl.size() : 0;
      Vec6 qn = Vec6::Zero();
      std::vector<Eigen::Matrix<double, 9, 1>> wl(nl);
      std::vector<Eigen::Matrix<double, 9, 12>> Wl(nl);
      std::array<ad::Grad<12>, 12> xg12;
      if (nl > 0) {
        for (int i = 0; i < 12; ++i) xg12[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(ue(i), i);
        const corot::M3t<ad::Grad<12>> R1 = corot::exp_t<ad::Grad<12>>({xg12[3], xg12[4], xg12[5]});
        const corot::M3t<ad::Grad<12>> R2 = corot::exp_t<ad::Grad<12>>({xg12[9], xg12[10], xg12[11]});
        const corot::M3t<ad::Grad<12>> E = nb.ce.frame_t(xg12, R1, R2);
        nb.pss_cur.setZero();
        if (nb.fiber)
          for (Vec3& z : nb.fe.sp) z.setZero();
        for (std::size_t l = 0; l < nl; ++l) {
          const NlBar::LdLoad& L = nb.ldl[l];
          const std::array<ad::Grad<12>, 9> w = corot::local_intensities(nb.ce, L.w, xg12, E);
          for (int c = 0; c < 9; ++c) {
            wl[l](c) = w[static_cast<std::size_t>(c)].v;
            for (int j = 0; j < 12; ++j) Wl[l](c, j) = w[static_cast<std::size_t>(c)].g[static_cast<std::size_t>(j)];
          }
          nb.pss_cur += lam_cur * (L.ss.basis * wl[l]);
          if (nb.fiber) {
            for (int c = 0; c < 9; ++c)
              if (wl[l](c) != 0.0 && !L.S[static_cast<std::size_t>(c)].empty())
                for (std::size_t k = 0; k < nb.fe.sp.size(); ++k) nb.fe.sp[k] += lam_cur * wl[l](c) * L.S[static_cast<std::size_t>(c)][k];
            nb.fe.vshift(5) += lam_cur * L.V0.row(5).dot(wl[l]);  // the twist of the span torques
          } else {
            v -= lam_cur * (L.V0 * wl[l]);  // q = kb (v - v0 - P e) + qn
            qn += lam_cur * (L.Qn * wl[l]);
          }
        }
      }
      if (ld && nb.fiber && nb.fe.has_bed()) {
        // FIBER BAR ON A WINKLER BED, large displacements: the coupled element
        // in the corotated frame — the bed reacts to the ABSOLUTE axis
        // displacement c + E delta (chord + deflection from the strains), in
        // the undeformed local axes (Ks = E0 diag(k) E0^T); local pressures P
        // with M = E^T Ks E, g0 = E^T Ks c; the reactions Rop P pushed back
        const corot::CorotElem& ce = nb.ce;
        std::array<ad::Grad<12>, 12> xu;
        for (int i = 0; i < 12; ++i) xu[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(ue(i), i);
        const corot::M3t<ad::Grad<12>> R1 = corot::exp_t<ad::Grad<12>>({xu[3], xu[4], xu[5]});
        const corot::M3t<ad::Grad<12>> R2 = corot::exp_t<ad::Grad<12>>({xu[9], xu[10], xu[11]});
        const corot::M3t<ad::Grad<12>> Eg = ce.frame_t(xu, R1, R2);
        const corot::V3t<ad::Grad<12>> x21 = ce.chord_t(xu, R1, R2);
        const corot::V3t<ad::Grad<12>> arm1 = ce.arm_t(0, R1);
        Mat3 Ed;
        for (int i = 0; i < 9; ++i) Ed(i / 3, i % 3) = Eg[static_cast<std::size_t>(i)].v;
        const Mat3 Ks = ce.E0 * nb.fe.kbed.asDiagonal() * ce.E0.transpose();
        const std::size_t np = nb.fe.xg.size();
        const int n3 = static_cast<int>(3 * np);
        const int nfq = static_cast<int>(nb.fe.Cq.cols());
        // absolute chord displacement at the points, generic
        std::vector<corot::V3t<ad::Grad<12>>> cg(np);
        Eigen::VectorXd g0(n3);
        const Vec3 X21 = ce.X2 - ce.X1;
        const Mat3 EtK = Ed.transpose() * Ks;
        for (std::size_t k = 0; k < np; ++k) {
          const double xi = nb.fe.xg[k] / nb.L;
          Vec3 cd;
          for (int c = 0; c < 3; ++c) {
            cg[k][static_cast<std::size_t>(c)] = xu[static_cast<std::size_t>(c)] + arm1[static_cast<std::size_t>(c)] -
                                                 ad::Grad<12>::constant(ce.a1(c)) + x21[static_cast<std::size_t>(c)] * xi -
                                                 ad::Grad<12>::constant(X21(c) * xi);
            cd(c) = cg[k][static_cast<std::size_t>(c)].v;
          }
          g0.segment<3>(static_cast<Eigen::Index>(3 * k)) = EtK * cd;
        }
        nb.fe.solve_bed(v, Ed.transpose() * Ks * Ed, g0);
        const Vec6 qb = nb.fe.q6();
        {
          Vec6 vs = Vec6::Zero();
          vs.head<5>() = nb.fe.section_noise();
          const Vec12 nl = B.transpose().cwiseAbs() * (nb.fe.kb * vs).cwiseAbs() +
                           (nb.fe.Rop.cwiseAbs() * nb.fe.bed_pnoise);
          elem_noise = std::max(elem_noise, nl.sum());
        }
        // dr2/du at fixed strains: h(u) = E^T Ks (c + E delta)
        Eigen::MatrixXd dh(n3, 12);
        for (std::size_t k = 0; k < np; ++k) {
          const Vec3 dl = nb.fe.bed_delta.segment<3>(static_cast<Eigen::Index>(3 * k));
          corot::V3t<ad::Grad<12>> ua;
          for (int c = 0; c < 3; ++c) {
            ua[static_cast<std::size_t>(c)] = cg[k][static_cast<std::size_t>(c)];
            for (int j = 0; j < 3; ++j) ua[static_cast<std::size_t>(c)] = ua[static_cast<std::size_t>(c)] + Eg[static_cast<std::size_t>(3 * c + j)] * dl(j);
          }
          for (int a = 0; a < 3; ++a) {
            // (E^T Ks ua)_a = sum_c E(c, a) (Ks ua)_c
            ad::Grad<12> h = ad::Grad<12>::constant(0.0);
            for (int c = 0; c < 3; ++c) {
              ad::Grad<12> ksu = ad::Grad<12>::constant(0.0);
              for (int j = 0; j < 3; ++j)
                if (Ks(c, j) != 0.0) ksu = ksu + ua[static_cast<std::size_t>(j)] * Ks(c, j);
              h = h + Eg[static_cast<std::size_t>(3 * c + a)] * ksu;
            }
            for (int j = 0; j < 12; ++j) dh(static_cast<Eigen::Index>(3 * k + a), j) = h.g[static_cast<std::size_t>(j)];
          }
        }
        // the bed reactions pushed back, and the push operator Pu = d push / d floc
        const Vec12 fl = nb.fe.Rop * nb.fe.Pbed;
        std::array<ad::Grad<12>, 12> flg;
        for (int i = 0; i < 12; ++i) flg[static_cast<std::size_t>(i)] = ad::Grad<12>::constant(fl(i));
        const std::array<ad::Grad<12>, 12> pg = corot::push_local(ce, xu, R1, R2, Eg, flg);
        Mat12 Pu;
        {
          const std::array<double, 12> qd = corot::CorotElem::arr(ue);
          const corot::M3t<double> R1d = corot::exp_t<double>({qd[3], qd[4], qd[5]}), R2d = corot::exp_t<double>({qd[9], qd[10], qd[11]});
          const corot::M3t<double> Edd = ce.frame_t(qd, R1d, R2d);
          for (int j = 0; j < 12; ++j) {
            std::array<double, 12> unit{};
            unit[static_cast<std::size_t>(j)] = 1.0;
            const std::array<double, 12> col = corot::push_local(ce, qd, R1d, R2d, Edd, unit);
            for (int i = 0; i < 12; ++i) Pu(i, j) = col[static_cast<std::size_t>(i)];
          }
        }
        Vec12 fe = B.transpose() * qb;
        Mat12 kk = Mat12::Zero();
        for (int r = 0; r < 12; ++r) {
          fe(r) += pg[static_cast<std::size_t>(r)].v;
          for (int c = 0; c < 12; ++c) kk(r, c) += pg[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)];
        }
        nb.pss_cur = fl;
        // tangent: (dQ, dP)/du from the coupled system
        Eigen::MatrixXd rhs(nfq + n3, 12);
        rhs.topRows(nfq) = nb.fe.Cq.transpose() * B.topRows<5>();
        rhs.bottomRows(n3) = dh;
        const Eigen::MatrixXd dz = nb.fe.bed_sensitivity(rhs);
        Eigen::Matrix<double, 6, 12> dq6 = Eigen::Matrix<double, 6, 12>::Zero();
        dq6.topRows<5>() = nb.fe.Cq * dz.topRows(nfq);
        dq6.row(5) = nb.fe.torsion_stiffness() * B.row(5);
        kk += B.transpose() * dq6 + Pu * (nb.fe.Rop * dz.bottomRows(n3));
        // span loads on the current frame: their reactions pushed, their
        // particular fields' sensitivity through the coupled system
        for (std::size_t l = 0; l < nl; ++l) {
          const NlBar::LdLoad& Ld = nb.ldl[l];
          const std::array<ad::Grad<12>, 12> f = corot::bar_load_vector(ce, Ld.ss, xg12);
          for (int r = 0; r < 12; ++r) {
            fe(r) += lam_cur * f[static_cast<std::size_t>(r)].v;
            for (int c = 0; c < 12; ++c) kk(r, c) += lam_cur * f[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)];
          }
          nb.pss_cur += lam_cur * (Ld.ss.basis * wl[l]);
          if (!K) continue;
          Eigen::MatrixXd dsp = Eigen::MatrixXd::Zero(n3, 9);
          for (int c = 0; c < 9; ++c)
            if (!Ld.S[static_cast<std::size_t>(c)].empty())
              for (std::size_t k = 0; k < np; ++k) dsp.block<3, 1>(static_cast<Eigen::Index>(3 * k), c) = Ld.S[static_cast<std::size_t>(c)][k];
          const Eigen::MatrixXd dzw = nb.fe.bed_sensitivity(nb.fe.bed_dr_dsp(dsp));
          Eigen::Matrix<double, 6, 9> dqw = Eigen::Matrix<double, 6, 9>::Zero();
          dqw.topRows<5>() = nb.fe.Cq * dzw.topRows(nfq);
          dqw.row(5) = -nb.fe.torsion_stiffness() * Ld.V0.row(5);
          kk += lam_cur * ((B.transpose() * dqw + Pu * (nb.fe.Rop * dzw.bottomRows(n3))) * Wl[l]);
        }
        if (nb.ce.has_bed()) {  // its twist bed
          fe += nb.ce.bed_grad(ue);
          kk += nb.ce.bed_hess(ue);
        }
        for (int a = 0; a < 12; ++a)
          if (nb.gdof[a] >= 0) fint(nb.gdof[a]) += fe(a);
        if (K) {
          // + d2(q . v)/du2 at frozen q
          std::array<ad::Jet<12>, 12> xj;
          for (int i = 0; i < 12; ++i) xj[static_cast<std::size_t>(i)] = ad::Jet<12>::seed(ue(i), i);
          const std::array<ad::Jet<12>, 7> dd = ce.deform_t(xj);
          ad::Jet<12> g = dd[0] * qb(0);
          for (int k = 1; k < 6; ++k) {
            g = g + dd[static_cast<std::size_t>(MAPV[k][0])] * qb(k);
            if (MAPV[k][1]) g = g - dd[1] * qb(k);
          }
          for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j) kk(i, j) += g.h[static_cast<std::size_t>(12 * i + j)];
          for (int a = 0; a < 12; ++a)
            for (int b = 0; b < 12; ++b)
              if (nb.gdof[a] >= 0 && nb.gdof[b] >= 0 && kk(a, b) != 0.0) trips.emplace_back(nb.gdof[a], nb.gdof[b], kk(a, b));
        }
        if (out_q) (*out_q)[nb.e] = qb;
        continue;
      }
      Mat6 Gn = Mat6::Identity();
      nb.response(v, q, Kq, qn, nl > 0 && !nb.fiber ? &Gn : nullptr);
      if (nb.fiber) {
        Vec6 vs = Vec6::Zero();
        vs.head<5>() = nb.fe.section_noise();
        elem_noise = std::max(elem_noise, (B.transpose() * (Kq * vs)).cwiseAbs().sum());
      }
      Vec12 fe = B.transpose() * q;
      if (nb.span) fe += lam_cur * nb.Tn.transpose() * (nb.p_ss + nb.f0c);  // the case vector carries T^t f0c
      Mat12 kspan = Mat12::Zero();  // non-symmetric span-load tangent (large displacements)
      add_tor(fe, &kspan);
      if (ld && nb.ce.has_bed()) {
        // hinged bar on a Winkler bed: the ground springs on the current
        // axis (a fiber bar: only its twist bed here)
        const Eigen::Matrix<double, 12, 1> q12 = ue;
        fe += nb.ce.bed_grad(q12);
        if (K) kspan += nb.ce.bed_hess(q12);
      }
      for (std::size_t l = 0; l < nl; ++l) {
        const NlBar::LdLoad& L = nb.ldl[l];
        // basic-system reactions pushed back like equivalent loads (exact AD)
        const std::array<ad::Grad<12>, 12> f = corot::bar_load_vector(nb.ce, L.ss, xg12);
        for (int r = 0; r < 12; ++r) {
          fe(r) += lam_cur * f[static_cast<std::size_t>(r)].v;
          for (int c = 0; c < 12; ++c) kspan(r, c) += lam_cur * f[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)];
        }
        if (!K) continue;
        // dq/du through the intensities: a basic deformation shift V w(u) —
        // V0 (hinged, constant) or the converged compliance to the
        // particular fields G(S_c) (fiber): dq = -Kq V W du
        Eigen::Matrix<double, 6, 9> V = L.V0;
        if (nb.fiber) {
          V.topRows<5>().setZero();  // (row 5: the torsion twist per unit intensity)
          for (int c = 0; c < 9; ++c)
            if (!L.S[static_cast<std::size_t>(c)].empty()) V.col(c).head<5>() = nb.fe.sp_compliance(L.S[static_cast<std::size_t>(c)]);
        }
        kspan -= lam_cur * (B.transpose() * (Kq * (V * Wl[l])));
        if (!nb.fiber) kspan += lam_cur * (B.transpose() * (Gn * (L.Qn * Wl[l])));
      }
      for (int a = 0; a < 12; ++a)
        if (nb.gdof[a] >= 0) fint(nb.gdof[a]) += fe(a);
      if (K) {
        for (int a = 0; a < 12; ++a)
          for (int b = 0; b < 12; ++b)
            if (nb.gdof[a] >= 0 && nb.gdof[b] >= 0 && kspan(a, b) != 0.0) trips.emplace_back(nb.gdof[a], nb.gdof[b], kspan(a, b));
        Mat12 ke = B.transpose() * Kq * B;
        if (ld) {  // d2(q . v)/du2 at frozen q
          std::array<ad::Jet<12>, 12> x;
          for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Jet<12>::seed(ue(i), i);
          const std::array<ad::Jet<12>, 7> dd = nb.ce.deform_t(x);
          ad::Jet<12> g = dd[0] * q(0);
          for (int k = 1; k < 6; ++k) {
            g = g + dd[static_cast<std::size_t>(MAPV[k][0])] * q(k);
            if (MAPV[k][1]) g = g - dd[1] * q(k);
          }
          for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j) ke(i, j) += g.h[static_cast<std::size_t>(12 * i + j)];
        }
        for (int a = 0; a < 12; ++a)
          for (int b = 0; b < 12; ++b)
            if (nb.gdof[a] >= 0 && nb.gdof[b] >= 0) trips.emplace_back(nb.gdof[a], nb.gdof[b], 0.5 * (ke(a, b) + ke(b, a)));
      }
      if (out_q) (*out_q)[nb.e] = q;
    }
    fint += K_lin * u;
    if (K) {
      K->resize(d.n_free, d.n_free);
      K->setFromTriplets(trips.begin(), trips.end());
      *K += K_lin;
    }
  };

  ShellPlasticResult res;
  VecX u = VecX::Zero(d.n_free);
  FactorCache fc;
  if (red) fc.set_reduction(asm_.red.T);
  std::vector<Vec6> basic_q(model.bars.size(), Vec6::Zero());
  std::vector<Vec12> end12(model.bars.size(), Vec12::Zero());
  out_q = &basic_q;
  out_end = &end12;
  auto recanonicalize = [&](VecX& v) {
    for (int n = 0; n < d.n_nodes; ++n) {
      int fx[3];
      rot_dofs(n, fx);
      Vec3 th = Vec3::Zero();
      for (int c = 0; c < 3; ++c)
        if (fx[c] >= 0) th(c) = v(fx[c]);
      const double t = th.norm();
      if (t > 0.9 * kPi) {
        const Vec3 th2 = (1.0 - 2 * kPi / t) * th;
        for (int c = 0; c < 3; ++c)
          if (fx[c] >= 0) v(fx[c]) = th2(c);
      }
    }
  };
  // load increments: 1 / steps; with stop_at_limit a failed increment is
  // halved (limit_cuts times) from the last converged state, then the
  // analysis ends there
  double lam_prev = 0, dl = 1.0 / opt.steps;
  int cuts = 0, st = 0;
  while (lam_prev < 1.0 - 1e-12) {
    const double lam = cuts == 0 ? static_cast<double>(st + 1) / opt.steps : std::min(1.0, lam_prev + dl);
    lam_cur = lam;
    int it = 0;
    double du_last = std::numeric_limits<double>::infinity();
    double rn_last = std::numeric_limits<double>::infinity();  // residual the last correction answered
    const VecX u_start = u;
    try {
    for (;; ++it) {
      if (it >= opt.max_iterations)
        throw NonConvergence("shell plasticity: no equilibrium at load factor " + std::to_string(lam) +
                             " (limit load exceeded or too large a step)");
      VecX fint;
      SpMat K;
      evaluate(u, fint, &K);
      std::vector<Triplet> kload;
      const VecX Fl = external(u, lam, ld ? &kload : nullptr);
      const VecX r = Fl - fint;
      const double rn = red ? VecX(asm_.red.reduce(r)).norm() : r.norm();
      // (the case vector F still holds the bars' equivalent span loads, which
      // large displacements carry inside the nonlinear bars instead)
      const double scale = std::max({Fl.norm(), lam * F.norm(), 1e-300});
      // round-off floor of K u: the stiff membrane / drilling terms of thin
      // facets (and long members) carry absolute errors ~ eps ||K_e|| |u|
      // that no iteration removes; below it the residual is noise
      // — accepted only once Newton has converged in displacement (the last
      // correction negligible): near a mechanism |u| grows and the
      // corrections do not shrink, so a collapse is never mistaken for it
      const double ufloor = 64 * std::numeric_limits<double>::epsilon() * kmax * u.lpNorm<Eigen::Infinity>();
      if (rn <= opt.tol * scale) break;
      // (and the fiber elements' section-solve noise, kb times its
      // deformation floor pulled back to the nodes)
      const bool noise = rn <= std::max({force_floor, ufloor, elem_noise * std::sqrt(static_cast<double>(d.n_free))});
      // (the correction the noise itself induces, (|du| / |r|)_last x floor,
      // bounds the attainable accuracy: once below it — and far below |u|,
      // which a mechanism's drift is not — the state is as converged as the
      // floors allow)
      const double floor_r = std::max({force_floor, ufloor, elem_noise * std::sqrt(static_cast<double>(d.n_free))});
      const double du_noise = std::isfinite(rn_last) && rn_last > 0 ? 4 * du_last / rn_last * floor_r : 0.0;
      if (it > 0 && noise && (du_last <= 1e-8 * u.norm() || (du_last <= du_noise && du_last <= 1e-2 * u.norm()))) break;
      rn_last = rn;
      VecX du;
      if (!ld && !nonsym) {
        try {
          fc.factorize(K);
        } catch (const Error&) {
          throw NonConvergence("shell plasticity: singular tangent at load factor " + std::to_string(lam) +
                               " (mechanism: the limit load is reached)");
        }
        du = fc.solve(r);
      } else {
        SpMat KL(d.n_free, d.n_free);
        KL.setFromTriplets(kload.begin(), kload.end());
        SpMat Kt = K - KL;  // d(f_int - f_ext)/du (non-symmetric: load stiffness, bedded fibers)
        if (red) Kt = asm_.red.reduce(Kt);
        Kt.makeCompressed();
        Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
        lu.compute(Kt);
        if (lu.info() != Eigen::Success)
          throw NonConvergence("shell plasticity: singular tangent at load factor " + std::to_string(lam));
        du = red ? asm_.red.expand(lu.solve(asm_.red.reduce(r))) : VecX(lu.solve(r));
      }
      require(du.allFinite(), "shell plasticity: non-finite increment");
      if (!ld || noise) {
        // (below the round-off floor there is no descent to measure: the
        // full, noise-sized correction is taken and its size decides the stop)
        u += du;
        du_last = du.norm();
        if (ld) recanonicalize(u);
        continue;
      }
      // line search (Crisfield vol. 1 sec. 9.3, as solve_corotational): accept
      // once the directional derivative of the potential along du has dropped
      // below 0.8 of its start value, or the inf-norm residual has dropped —
      // the linearised predictor of a slender member carries a large but
      // harmless axial residual that a norm test would reject
      const double g0 = std::fabs(du.dot(r)), r0 = r.lpNorm<Eigen::Infinity>();
      double sfac = 1.0;
      bool ok = false;
      for (int back = 0; back < 20 && !ok; ++back, sfac *= 0.5) {
        VecX ut = u + sfac * du;
        try {
          VecX ft;
          evaluate(ut, ft, nullptr);
          const VecX rt = external(ut, lam, nullptr) - ft;
          if (rt.allFinite() && (std::fabs(du.dot(rt)) <= 0.8 * g0 || rt.lpNorm<Eigen::Infinity>() < r0)) {
            u = ut;
            du_last = sfac * du.norm();
            ok = true;
          }
        } catch (const Error&) {
        }
      }
      if (!ok) throw NonConvergence("shell plasticity: line search failed at load factor " + std::to_string(lam));
      recanonicalize(u);
    }
    } catch (const NonConvergence& e) {
      if (!opt.stop_at_limit) throw;
      u = u_start;
      if (cuts < opt.limit_cuts) {
        ++cuts;
        dl = (cuts == 1 ? 1.0 / opt.steps : dl) / 2;
        continue;
      }
      res.limit_reached = true;
      res.limit_message = e.what();
      break;
    }
    ++st;
    lam_prev = lam;
    committed = trial;
    for (NlBar& nb : nlbars) nb.commit();
    ShellPlasticStep s;
    s.lambda = lam;
    s.iterations = it;
    s.node_disp.assign(model.nodes.size(), Vec6::Zero());
    for (std::size_t n = 0; n < model.nodes.size(); ++n)
      for (int c = 0; c < 6; ++c) {
        const int fi = d.free_ix[static_cast<std::size_t>(d.full(static_cast<int>(n), c))];
        if (fi >= 0) s.node_disp[n](c) = u(fi);
      }
    // bar end forces (local frame), per model bar: elastic k u - lam f0; the
    // corotational elastic bars from Kbar d in the corotated frame; nonlinear
    // bars a6^T q of their converged basic forces
    s.bar_end_forces.assign(model.bars.size(), Vec12::Zero());
    s.bar_basic_forces = basic_q;
    std::size_t ci = 0;
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      if (nl_of[e] >= 0) {
        const NlBar& nb = nlbars[static_cast<std::size_t>(nl_of[e])];
        s.bar_end_forces[e] = nb.mode12 ? end12[e]
                                         : Vec12(fiberel::basic_transform6(nb.L).transpose() * basic_q[e] +
                                                 ((ld || nb.fe.has_bed()) ? nb.pss_cur : Vec12(lam * nb.p_ss)) +
                                                 nb.tor_cur);
        continue;
      }
      if (ld) {
        const corot::CorotElem& ce = cbars[ci++];
        Eigen::Matrix<double, 12, 1> q;
        q << s.node_disp[static_cast<std::size_t>(ce.ni)], s.node_disp[static_cast<std::size_t>(ce.nj)];
        const Eigen::Matrix<double, 7, 1> fl = ce.Kbar * ce.deform(q);
        Vec12 f = Vec12::Zero();
        f(6) = fl(0);
        f(0) = -fl(0);
        f.segment<3>(3) = fl.segment<3>(1);
        f.segment<3>(9) = fl.segment<3>(4);
        const double ln = ce.chord(q).norm();
        f(1) = -(fl(3) + fl(6)) / ln;
        f(7) = -f(1);
        f(2) = (fl(2) + fl(5)) / ln;
        f(8) = -f(2);
        s.bar_end_forces[e] = f;
        continue;
      }
      const Bar& b = model.bars[e];
      const BarElement& el = asm_.elements[e];
      Vec12 un;
      un << s.node_disp[static_cast<std::size_t>(model.node_index(b.ni))], s.node_disp[static_cast<std::size_t>(model.node_index(b.nj))];
      const Vec12 f0 = lam * cl.f0[e];
      const Vec12 ul = el.local_displacements(un, f0);
      const TaperedBar* tb = model.tapered_of(b.id);  // as recover_case
      s.bar_end_forces[e] = tb ? Vec12(tb->ke * ul - f0) : el.end_forces_local(ul, f0);
    }
    res.steps.push_back(std::move(s));
  }
  res.state = committed;
  res.z = rule.zeta;
  for (const Elem& el : els) res.per_shell_z.push_back(el.rule.zeta);
  return res;
}

}  // namespace keta

#include "keta/wind_flutter.h"

#include <algorithm>
#include <cmath>

namespace keta {

using C = std::complex<double>;

std::array<double, 8> ScanlanDerivatives::at(double K) const {
  require(!rows.empty(), "Scanlan derivatives: empty table");
  std::array<double, 8> v{};
  if (K <= rows.front()[0] || rows.size() == 1) {
    for (int i = 0; i < 8; ++i) v[i] = rows.front()[i + 1];
    return v;
  }
  if (K >= rows.back()[0]) {
    for (int i = 0; i < 8; ++i) v[i] = rows.back()[i + 1];
    return v;
  }
  const auto it = std::upper_bound(rows.begin(), rows.end(), K, [](double x, const auto& r) { return x < r[0]; });
  const auto& r1 = *it;
  const auto& r0 = *(it - 1);
  const double t = (K - r0[0]) / (r1[0] - r0[0]);
  for (int i = 0; i < 8; ++i) v[i] = r0[i + 1] + t * (r1[i + 1] - r0[i + 1]);
  return v;
}

namespace {
// Theodorsen's flat-plate matrices of flutter_speed (h down, alpha nose-up,
// row 0 = +L_up, row 1 = -M on the left-hand side), per unit span
void theodorsen_mats(double rho, double b, double a, double U, double om, Eigen::Matrix2d& Ma, Eigen::Matrix2cd& Ca,
                     Eigen::Matrix2cd& Ka) {
  const double k = std::max(1e-4, om * b / U);
  const C Ck = theodorsen_c(k);
  const double pb2 = kPi * rho * b * b;
  const C two = 2.0 * kPi * rho * U * b * Ck;
  Ma << pb2, -pb2 * b * a, -pb2 * b * a, pb2 * b * b * (0.125 + a * a);
  Ca(0, 0) = two;
  Ca(0, 1) = pb2 * U + two * b * (0.5 - a);
  Ca(1, 0) = -two * b * (a + 0.5);
  Ca(1, 1) = pb2 * U * b * (0.5 - a) - two * b * (a + 0.5) * b * (0.5 - a);
  Ka(0, 0) = 0;
  Ka(0, 1) = two * U;
  Ka(1, 0) = 0;
  Ka(1, 1) = -two * U * b * (a + 0.5);
}
// Scanlan impedance Z = Ka + i omega Ca at K = B omega / U
Eigen::Matrix2cd scanlan_z(const ScanlanDerivatives& d, double rho, double B, double U, double om) {
  const double K = B * om / U, q = 0.5 * rho * U * U;
  const std::array<double, 8> h = d.at(K);  // H1..H4, A1..A4
  Eigen::Matrix2cd Z;
  Z(0, 0) = -q * K * K * C(h[3], h[0]);
  Z(0, 1) = -q * K * K * B * C(h[2], h[1]);
  Z(1, 0) = -q * K * K * B * C(h[7], h[4]);
  Z(1, 1) = -q * K * K * B * B * C(h[6], h[5]);
  return Z;
}
}  // namespace

ScanlanDerivatives scanlan_from_theodorsen(double a, const std::vector<double>& Ks) {
  ScanlanDerivatives d;
  const double rho = 1, U = 1, b = 1, B = 2, q = 0.5;
  for (double K : Ks) {
    require(K > 0, "Scanlan from Theodorsen: K > 0");
    const double om = K * U / B;
    Eigen::Matrix2d Ma;
    Eigen::Matrix2cd Ca, Ka;
    theodorsen_mats(rho, b, a, U, om, Ma, Ca, Ka);
    const Eigen::Matrix2cd Z = -om * om * Ma.cast<C>() + C(0, om) * Ca + Ka;
    const double s = q * K * K;
    const C h4h1 = -Z(0, 0) / s, h3h2 = -Z(0, 1) / (s * B), a4a1 = -Z(1, 0) / (s * B), a3a2 = -Z(1, 1) / (s * B * B);
    d.rows.push_back({K, h4h1.imag(), h3h2.imag(), h3h2.real(), h4h1.real(), a4a1.imag(), a3a2.imag(), a3a2.real(),
                      a4a1.real()});
  }
  return d;
}

namespace {
// the modal matrices M, C, K of the p-k problem at (U, omega)
void modal_mats(const MultimodalFlutter& f, double U, double om, Eigen::MatrixXcd& M, Eigen::MatrixXcd& Cm,
                Eigen::MatrixXcd& Km) {
  const Eigen::Index m = static_cast<Eigen::Index>(f.omega.size());
  M = Eigen::MatrixXcd::Identity(m, m);
  Cm = Eigen::MatrixXcd::Zero(m, m);
  Km = Eigen::MatrixXcd::Zero(m, m);
  for (Eigen::Index i = 0; i < m; ++i) {
    const double w = f.omega[static_cast<std::size_t>(i)];
    Cm(i, i) = 2 * f.zeta[static_cast<std::size_t>(i)] * w;
    Km(i, i) = w * w;
  }
  for (const FlutterStrip& s : f.strips) {
    require(s.T.cols() == m, "flutter strip: one column per mode");
    const Eigen::MatrixXcd T = s.T.cast<C>();
    if (s.kind == FlutterStrip::Kind::Theodorsen) {
      Eigen::Matrix2d Ma;
      Eigen::Matrix2cd Ca, Ka;
      theodorsen_mats(f.rho, s.b, s.a, U, om, Ma, Ca, Ka);
      M += s.L * T.transpose() * Ma.cast<C>() * T;
      Cm += s.L * T.transpose() * Ca * T;
      Km += s.L * T.transpose() * Ka * T;
    } else {
      const Eigen::Matrix2cd Z = scanlan_z(s.der, f.rho, s.B, U, om);
      const Eigen::Matrix2cd Ka = Z.real().cast<C>(), Ca = (Z.imag() / om).cast<C>();
      Cm += s.L * T.transpose() * Ca * T;
      Km += s.L * T.transpose() * Ka * T;
    }
  }
}
}  // namespace

std::vector<std::pair<double, double>> multimodal_branches(const MultimodalFlutter& f, double U) {
  const Eigen::Index m = static_cast<Eigen::Index>(f.omega.size());
  require(m >= 1 && f.zeta.size() == f.omega.size() && !f.strips.empty() && U > 0, "multimodal flutter: modes, strips, U > 0");
  std::vector<std::pair<double, double>> out;
  for (Eigen::Index j = 0; j < m; ++j) {
    const double w0 = f.omega[static_cast<std::size_t>(j)];
    double om = w0;
    C p(0, om);
    for (int it = 0; it < 100; ++it) {
      Eigen::MatrixXcd M, Cm, Km;
      modal_mats(f, U, om, M, Cm, Km);
      const Eigen::MatrixXcd Mi = M.inverse();
      Eigen::MatrixXcd A = Eigen::MatrixXcd::Zero(2 * m, 2 * m);
      A.topRightCorner(m, m).setIdentity();
      A.bottomLeftCorner(m, m) = -Mi * Km;
      A.bottomRightCorner(m, m) = -Mi * Cm;
      const Eigen::ComplexEigenSolver<Eigen::MatrixXcd> es(A, false);
      double best = 1e300;
      C pick = p;
      for (Eigen::Index i = 0; i < 2 * m; ++i) {
        const C pi = es.eigenvalues()(i);
        if (pi.imag() <= 0) continue;
        const double dist = std::abs(pi.imag() - om);
        if (dist < best) {
          best = dist;
          pick = pi;
        }
      }
      p = pick;
      const double om_new = pick.imag();
      if (std::abs(om_new - om) < 1e-10 * w0) {
        om = om_new;
        break;
      }
      om = om_new;
    }
    out.push_back({om, -p.real() / std::abs(p)});
  }
  return out;
}

FlutterResult multimodal_flutter(const MultimodalFlutter& f, double u_min, double u_max, int n_sweep) {
  require(u_min > 0 && u_max > u_min && n_sweep >= 1, "multimodal flutter: 0 < u_min < u_max");
  auto least = [&](double U) {
    double z = 1e300, om = 0;
    for (const auto& [w, zz] : multimodal_branches(f, U))
      if (zz < z) {
        z = zz;
        om = w;
      }
    return std::make_pair(z, om);
  };
  FlutterResult out;
  double u_prev = u_min;
  double z_prev = least(u_min).first;
  for (int i = 1; i <= n_sweep; ++i) {
    const double U = u_min + (u_max - u_min) * i / n_sweep;
    auto [z, om] = least(U);
    if (z_prev >= 0 && z < 0) {
      double lo = u_prev, hi = U;
      for (int it = 0; it < 60; ++it) {
        const double mid = 0.5 * (lo + hi);
        const auto [zm, omm] = least(mid);
        (zm >= 0 ? lo : hi) = mid;
        om = omm;
      }
      out.U_flutter = 0.5 * (lo + hi);
      out.omega_flutter = om;
      out.found = true;
      return out;
    }
    u_prev = U;
    z_prev = z;
  }
  return out;
}

std::vector<FlutterStrip> flutter_strips(const Model& model, const ModalResults& modal, const Vec3& wind_dir,
                                         const FlutterStrip& proto) {
  Vec3 d = wind_dir;
  d.z() = 0;
  require(d.norm() > 1e-9, "flutter strips: horizontal wind direction");
  d.normalize();
  const Vec3 z(0, 0, 1), ax = z.cross(d);
  const Assembled asm_ = assemble_elements(model);
  const Eigen::Index m = static_cast<Eigen::Index>(modal.modes_free.size());
  std::vector<FlutterStrip> out;
  for (const Bar& b : model.bars) {
    if (b.kind == BarKind::Cable) continue;
    FlutterStrip s = proto;
    s.L = (model.node(b.nj).x - model.node(b.ni).x).norm();
    s.T.resize(2, m);
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (Eigen::Index i = 0; i < m; ++i) {
      Vec6 v = Vec6::Zero();
      for (int c = 0; c < 6; ++c) {
        const int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        const int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nj, c))];
        const VecX& ph = modal.modes_free[static_cast<std::size_t>(i)];
        v(c) = 0.5 * ((f1 >= 0 ? ph(f1) : 0.0) + (f2 >= 0 ? ph(f2) : 0.0));
      }
      s.T(0, i) = -v.head<3>().dot(z);  // h downward
      s.T(1, i) = v.tail<3>().dot(ax);  // alpha nose-up
    }
    out.push_back(std::move(s));
  }
  return out;
}

}  // namespace keta

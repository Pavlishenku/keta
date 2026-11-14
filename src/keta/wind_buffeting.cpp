#include "keta/wind_buffeting.h"

#include <algorithm>
#include <cmath>
#include <complex>

namespace keta {

using CMat = Eigen::MatrixXcd;

// ---------------- coefficients vs incidence ----------------

AeroTable AeroTable::from_section(const AeroSection& s) {
  AeroTable t;
  t.alpha = {-1.0, 1.0};
  t.cd = {s.cd, s.cd};  // AeroSection carries no drag slope
  t.cl = {s.cl - s.dcl, s.cl + s.dcl};
  t.cm = {s.cm - s.dcm, s.cm + s.dcm};
  return t;
}

AeroTable::Value AeroTable::at(double a) const {
  const std::size_t n = alpha.size();
  require(n >= 1 && cd.size() == n && cl.size() == n && cm.size() == n, "aero table: same length columns, >= 1 row");
  for (std::size_t i = 1; i < n; ++i) require(alpha[i] > alpha[i - 1], "aero table: incidence strictly increasing");
  Value v;
  if (n == 1 || a <= alpha.front() || a >= alpha.back()) {
    const std::size_t k = (n == 1 || a <= alpha.front()) ? 0 : n - 1;
    v.cd = cd[k];
    v.cl = cl[k];
    v.cm = cm[k];
    return v;  // held flat outside: no slope
  }
  const std::size_t k = static_cast<std::size_t>(std::upper_bound(alpha.begin(), alpha.end(), a) - alpha.begin()) - 1;
  const double h = alpha[k + 1] - alpha[k], t = (a - alpha[k]) / h;
  auto lin = [&](const std::vector<double>& c) { return c[k] + t * (c[k + 1] - c[k]); };
  auto slope = [&](const std::vector<double>& c, std::size_t i) { return (c[i + 1] - c[i]) / (alpha[i + 1] - alpha[i]); };
  auto dc = [&](const std::vector<double>& c) {
    // on a row: the mean of the two segments (one segment at the table ends)
    if (t == 0 && k > 0) return 0.5 * (slope(c, k - 1) + slope(c, k));
    return slope(c, k);
  };
  v.cd = lin(cd);
  v.cl = lin(cl);
  v.cm = lin(cm);
  v.dcd = dc(cd);
  v.dcl = dc(cl);
  v.dcm = dc(cm);
  return v;
}

// ---------------- stage 1 ----------------

std::vector<double> buffeting_grid(const ModalSystem& sys, double n_min, double n_max, int n_log, int per_mode) {
  require(n_min > 0 && n_max > n_min && n_log >= 2 && per_mode >= 0, "buffeting grid: 0 < n_min < n_max, n_log >= 2");
  std::vector<double> g;
  for (int k = 0; k <= n_log; ++k) g.push_back(n_min * std::pow(n_max / n_min, double(k) / n_log));
  for (std::size_t i = 0; i < sys.omega.size(); ++i) {
    const double ni = sys.omega[i] / (2 * kPi), z = std::max(sys.zeta[i], 1e-4);
    for (int k = 0; k < per_mode; ++k) {
      const double n = ni * (1 + 25 * z * (2.0 * k / std::max(per_mode - 1, 1) - 1));
      if (n > n_min && n < n_max) g.push_back(n);
    }
  }
  std::sort(g.begin(), g.end());
  g.erase(std::unique(g.begin(), g.end()), g.end());
  return g;
}

namespace {
void check_system(const ModalSystem& s) {
  const std::size_t m = s.omega.size();
  require(m >= 1 && s.zeta.size() == m, "modal system: one damping ratio per mode");
  require(s.C_ae.size() == 0 || (s.C_ae.rows() == Eigen::Index(m) && s.C_ae.cols() == Eigen::Index(m)),
          "modal system: C_ae m x m");
  require(s.K_ae.size() == 0 || (s.K_ae.rows() == Eigen::Index(m) && s.K_ae.cols() == Eigen::Index(m)),
          "modal system: K_ae m x m");
}
// H(omega) and H(0) of the modal system
CMat frf(const ModalSystem& s, double om) {
  const Eigen::Index m = static_cast<Eigen::Index>(s.omega.size());
  CMat D = CMat::Zero(m, m);
  for (Eigen::Index i = 0; i < m; ++i) {
    const double w = s.omega[static_cast<std::size_t>(i)];
    D(i, i) = std::complex<double>(w * w - om * om, 2 * s.zeta[static_cast<std::size_t>(i)] * w * om);
  }
  if (s.K_ae.size()) D += s.K_ae.cast<std::complex<double>>();
  if (s.C_ae.size()) D += std::complex<double>(0, om) * s.C_ae.cast<std::complex<double>>();
  return D.partialPivLu().inverse();
}
Eigen::MatrixXd static_frf(const ModalSystem& s) {
  const Eigen::Index m = static_cast<Eigen::Index>(s.omega.size());
  Eigen::MatrixXd D = Eigen::MatrixXd::Zero(m, m);
  for (Eigen::Index i = 0; i < m; ++i) D(i, i) = s.omega[static_cast<std::size_t>(i)] * s.omega[static_cast<std::size_t>(i)];
  if (s.K_ae.size()) D += s.K_ae;
  return D.partialPivLu().inverse();
}
}  // namespace

static void check_stability(const ModalSystem& s) {
  if (!s.C_ae.size() && !s.K_ae.size()) return;
  // the state matrix [0 I; -(Omega^2 + K_ae) -(2 zeta Omega + C_ae)]: every
  // eigenvalue in the open left half-plane, else the spectral response does
  // not exist (galloping, flutter or divergence)
  const Eigen::Index m = static_cast<Eigen::Index>(s.omega.size());
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2 * m, 2 * m);
  A.topRightCorner(m, m).setIdentity();
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(m, m), C = K;
  for (Eigen::Index i = 0; i < m; ++i) {
    K(i, i) = s.omega[static_cast<std::size_t>(i)] * s.omega[static_cast<std::size_t>(i)];
    C(i, i) = 2 * s.zeta[static_cast<std::size_t>(i)] * s.omega[static_cast<std::size_t>(i)];
  }
  if (s.K_ae.size()) K += s.K_ae;
  if (s.C_ae.size()) C += s.C_ae;
  A.bottomLeftCorner(m, m) = -K;
  A.bottomRightCorner(m, m) = -C;
  const Eigen::EigenSolver<Eigen::MatrixXd> es(A, false);
  require(es.eigenvalues().real().maxCoeff() < 0,
          "buffeting: the aeroelastic system is unstable (negative total damping or divergence) at this speed");
}

std::vector<ResponseStats> modal_response(const ModalSystem& sys, const std::vector<double>& grid,
                                          const std::function<Eigen::MatrixXd(double)>& S_Q,
                                          const std::vector<Eigen::VectorXd>& psi, double T) {
  check_system(sys);
  check_stability(sys);
  require(grid.size() >= 2, "modal response: a grid of 2 frequencies at least");
  const std::size_t R = psi.size();
  std::vector<double> var(R, 0), varb(R, 0), varn(R, 0), ps(R, 0), pb(R, 0), pn(R, 0);
  const Eigen::MatrixXd H0 = static_frf(sys);
  for (std::size_t g = 0; g < grid.size(); ++g) {
    const double n = grid[g];
    const Eigen::MatrixXd SQ = S_Q(n);
    const CMat H = frf(sys, 2 * kPi * n);
    const CMat Sq = H * SQ.cast<std::complex<double>>() * H.adjoint();
    const Eigen::MatrixXd Sb = H0 * SQ * H0.transpose();
    for (std::size_t r = 0; r < R; ++r) {
      const double s = (psi[r].transpose().cast<std::complex<double>>() * Sq * psi[r].cast<std::complex<double>>())(0, 0).real();
      const double sb = psi[r].dot(Sb * psi[r]);
      if (g > 0) {
        const double dn = n - grid[g - 1];
        var[r] += 0.5 * (ps[r] + s) * dn;
        varb[r] += 0.5 * (pb[r] + sb) * dn;
        varn[r] += 0.5 * (pn[r] + s * n * n) * dn;
      }
      ps[r] = s;
      pb[r] = sb;
      pn[r] = s * n * n;
    }
  }
  std::vector<ResponseStats> out(R);
  for (std::size_t r = 0; r < R; ++r) {
    ResponseStats& o = out[r];
    o.sigma = std::sqrt(std::max(var[r], 0.0));
    o.sigma_b = std::sqrt(std::max(varb[r], 0.0));
    o.sigma_r = std::sqrt(std::max(var[r] - varb[r], 0.0));
    o.nu = var[r] > 0 ? std::sqrt(varn[r] / var[r]) : 0;
    const double x = o.nu * T;
    if (x > 1) {
      const double l = std::sqrt(2 * std::log(x));
      o.peak_factor = l + 0.5772156649 / l;
    }
  }
  return out;
}

Eigen::MatrixXd modal_covariance(const ModalSystem& sys, const std::vector<double>& grid,
                                 const std::function<Eigen::MatrixXd(double)>& S_Q) {
  check_system(sys);
  check_stability(sys);
  const Eigen::Index m = static_cast<Eigen::Index>(sys.omega.size());
  Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(m, m), prev = cov;
  for (std::size_t g = 0; g < grid.size(); ++g) {
    const CMat H = frf(sys, 2 * kPi * grid[g]);
    const Eigen::MatrixXd Sq = (H * S_Q(grid[g]).cast<std::complex<double>>() * H.adjoint()).real();
    if (g > 0) cov += 0.5 * (prev + Sq) * (grid[g] - grid[g - 1]);
    prev = Sq;
  }
  return cov;
}

// ---------------- stage 2 ----------------

namespace {
struct Axes {
  Vec3 d, z, a;
};
Axes axes_of(const Vec3& wind_dir) {
  Vec3 d = wind_dir;
  d.z() = 0;
  require(d.norm() > 1e-9, "buffeting: the wind direction needs a horizontal component");
  d.normalize();
  const Vec3 z(0, 0, 1);
  return {d, z, z.cross(d)};
}
// the 3 x 2 quasi-steady map (u, w) -> (D, L, M) per unit length at speed U
Eigen::Matrix<double, 3, 2> qs_map(const BuffetingSpec& sp, double U) {
  const AeroTable::Value c = sp.aero.at(sp.alpha0);
  const double q = 0.5 * sp.turb.profile.rho * U * U, B = sp.B;
  Eigen::Matrix<double, 3, 2> A;
  A << q * B * 2 * c.cd / U, q * B * (c.dcd - c.cl) / U, q * B * 2 * c.cl / U, q * B * (c.dcl + c.cd) / U,
      q * B * B * 2 * c.cm / U, q * B * B * c.dcm / U;
  return A;
}
}  // namespace

std::vector<Strip> buffeting_strips(const Model& model, const ModalResults& modal, const BuffetingSpec& spec) {
  const Axes ax = axes_of(spec.wind_dir);
  const Assembled asm_ = assemble_elements(model);
  std::vector<Strip> out;
  for (const Bar& b : model.bars) {
    if (b.kind == BarKind::Cable) continue;
    const Vec3 x1 = model.node(b.ni).x, x2 = model.node(b.nj).x, mid = 0.5 * (x1 + x2);
    Strip s;
    s.L = (x2 - x1).norm();
    s.wp = wind_point(spec.turb, mid.dot(ax.a), mid.z());
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (std::size_t i = 0; i < modal.modes_free.size(); ++i) {
      Vec6 v = Vec6::Zero();
      for (int c = 0; c < 6; ++c) {
        const int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        const int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nj, c))];
        v(c) = 0.5 * ((f1 >= 0 ? modal.modes_free[i](f1) : 0.0) + (f2 >= 0 ? modal.modes_free[i](f2) : 0.0));
      }
      s.phi.push_back(Vec3(v.head<3>().dot(ax.d), v.head<3>().dot(ax.z), v.tail<3>().dot(ax.a)));
    }
    out.push_back(std::move(s));
  }
  require(!out.empty(), "buffeting: no bars");
  return out;
}

ModalSystem buffeting_system(const ModalResults& modal, const std::vector<Strip>& strips, const BuffetingSpec& spec) {
  ModalSystem s;
  s.omega = modal.omegas;
  for (std::size_t i = 0; i < modal.omegas.size(); ++i) {
    s.zeta.push_back(spec.damping.zeta_for(modal.omegas[i], i));
    require(s.zeta.back() >= 0, "buffeting: structural damping >= 0");
  }
  if (!spec.aeroelastic) return s;
  const Eigen::Index m = static_cast<Eigen::Index>(modal.omegas.size());
  s.C_ae = Eigen::MatrixXd::Zero(m, m);
  s.K_ae = Eigen::MatrixXd::Zero(m, m);
  const AeroTable::Value c = spec.aero.at(spec.alpha0);
  for (const Strip& st : strips) {
    const Eigen::Matrix<double, 3, 2> A = qs_map(spec, st.wp.U);
    Eigen::Matrix3d Cs = Eigen::Matrix3d::Zero(), Ks = Eigen::Matrix3d::Zero();
    Cs.col(0) = A.col(0);  // u -> -x_D'
    Cs.col(1) = A.col(1);  // w -> -x_L'
    const double q = 0.5 * spec.turb.profile.rho * st.wp.U * st.wp.U;
    Ks.col(2) = -q * spec.B * Vec3(c.dcd, c.dcl, spec.B * c.dcm);  // incidence = theta
    for (Eigen::Index i = 0; i < m; ++i)
      for (Eigen::Index j = 0; j < m; ++j) {
        const Vec3 &pi = st.phi[static_cast<std::size_t>(i)], &pj = st.phi[static_cast<std::size_t>(j)];
        s.C_ae(i, j) += st.L * pi.dot(Cs * pj);
        s.K_ae(i, j) += st.L * pi.dot(Ks * pj);
      }
  }
  return s;
}

Eigen::MatrixXd modal_force_psd(const std::vector<Strip>& strips, const BuffetingSpec& spec, double n) {
  const std::size_t N = strips.size(), m = strips.front().phi.size();
  // G_c (m x N): modal placement of the strip forces per velocity component
  Eigen::MatrixXd Gu(static_cast<Eigen::Index>(m), static_cast<Eigen::Index>(N)), Gw = Gu;
  Eigen::VectorXd su(static_cast<Eigen::Index>(N)), sw(static_cast<Eigen::Index>(N));
  for (std::size_t s = 0; s < N; ++s) {
    const Strip& st = strips[s];
    Eigen::Matrix<double, 3, 2> A = qs_map(spec, st.wp.U);
    for (int r = 0; r < 3; ++r) A.row(r) *= std::sqrt(spec.chi[static_cast<std::size_t>(r)].chi2_at(n, spec.B, st.wp.U));
    for (std::size_t i = 0; i < m; ++i) {
      const Vec3& p = st.phi[i];
      Gu(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(s)) = st.L * p.dot(A.col(0));
      Gw(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(s)) = st.L * p.dot(A.col(1));
    }
    const WindPoint& w = st.wp;
    su(static_cast<Eigen::Index>(s)) = w.sigma[0] > 0 ? std::sqrt(auto_spectrum(spec.turb.form, 0, n, w.U, w.sigma[0], w.L[0])) : 0;
    sw(static_cast<Eigen::Index>(s)) = w.sigma[2] > 0 ? std::sqrt(auto_spectrum(spec.turb.form, 2, n, w.U, w.sigma[2], w.L[2])) : 0;
  }
  Eigen::MatrixXd Su(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N)), Sw = Su;
  for (std::size_t s = 0; s < N; ++s)
    for (std::size_t t = s; t < N; ++t) {
      const WindPoint &a = strips[s].wp, &b = strips[t].wp;
      const double cu = s == t ? 1.0 : davenport_coherence(n, a.y - b.y, a.z - b.z, a.U, b.U, spec.turb.C[0][0], spec.turb.C[0][1]);
      const double cw = s == t ? 1.0 : davenport_coherence(n, a.y - b.y, a.z - b.z, a.U, b.U, spec.turb.C[2][0], spec.turb.C[2][1]);
      const Eigen::Index i = static_cast<Eigen::Index>(s), j = static_cast<Eigen::Index>(t);
      Su(i, j) = Su(j, i) = su(i) * su(j) * cu;
      Sw(i, j) = Sw(j, i) = sw(i) * sw(j) * cw;
    }
  return Gu * Su * Gu.transpose() + Gw * Sw * Gw.transpose();
}

MultiBuffetingResults spectral_buffeting_multi(const Model& model, const ModalResults& modal, const BuffetingSpec& spec,
                                               const std::vector<BuffetingRequest>& requests) {
  require(!modal.omegas.empty(), "buffeting: no modes");
  const std::vector<Strip> strips = buffeting_strips(model, modal, spec);
  MultiBuffetingResults out;
  out.system = buffeting_system(modal, strips, spec);
  const double n_min = 0.02 * modal.frequencies_hz.front(), n_max = 5 * modal.frequencies_hz.back();
  const std::vector<double> grid = buffeting_grid(out.system, n_min, n_max, spec.n_log, spec.per_mode);
  const Assembled asm_ = assemble_elements(model);
  std::vector<Eigen::VectorXd> psi;
  for (const BuffetingRequest& r : requests) {
    require(r.comp >= 0 && r.comp < 6 && model.node_index(r.node) >= 0, "buffeting: request node / component");
    const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(model.node_index(r.node), r.comp))];
    Eigen::VectorXd p = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(modal.omegas.size()));
    if (f >= 0)
      for (std::size_t i = 0; i < modal.modes_free.size(); ++i) p(static_cast<Eigen::Index>(i)) = modal.modes_free[i](f);
    psi.push_back(p);
  }
  const auto SQ = [&](double n) { return modal_force_psd(strips, spec, n); };
  out.stats = modal_response(out.system, grid, SQ, psi, spec.duration);
  out.modal_cov = modal_covariance(out.system, grid, SQ);
  return out;
}

LoadCase mean_buffeting_case(const Model& model, const BuffetingSpec& spec, Id case_id, const std::string& name) {
  const Axes ax = axes_of(spec.wind_dir);
  const AeroTable::Value c = spec.aero.at(spec.alpha0);
  LoadCase lc;
  lc.id = case_id;
  lc.name = name;
  for (const Bar& b : model.bars) {
    if (b.kind == BarKind::Cable) continue;
    const Vec3 x1 = model.node(b.ni).x, x2 = model.node(b.nj).x;
    auto q = [&](double z) {
      const double U = spec.turb.profile.vm(z);
      return 0.5 * spec.turb.profile.rho * U * U;
    };
    const double q1 = q(x1.z()), q2 = q(x2.z()), L = (x2 - x1).norm();
    lc.loads.push_back(Load::linear(b.id, q1 * spec.B * (c.cd * ax.d + c.cl * ax.z),
                                    q2 * spec.B * (c.cd * ax.d + c.cl * ax.z), LoadFrame::Global));
    if (c.cm != 0) {
      const double m1 = q1 * spec.B * spec.B * c.cm, m2 = q2 * spec.B * spec.B * c.cm;
      lc.loads.push_back(Load::point(b.id, 0.0, Vec3::Zero(), (L * (m1 / 3 + m2 / 6)) * ax.a));
      lc.loads.push_back(Load::point(b.id, L, Vec3::Zero(), (L * (m1 / 6 + m2 / 3)) * ax.a));
    }
  }
  return lc;
}

}  // namespace keta

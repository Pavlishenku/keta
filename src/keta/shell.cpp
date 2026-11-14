#include "keta/shell.h"

#include <cmath>

namespace keta {

// ===========================================================================
// Local frame + small helpers
// ===========================================================================

Mat3 plane_frame(const Vec3& normal, const std::optional<Vec3>& local_x) {
  require(finite(normal), "plane frame: non-finite normal");
  double nl = normal.norm();
  require(nl > 0, "plane frame: zero-length normal");
  Vec3 zhat = normal / nl;

  // Seed for local x: the optional override, else global X with a global-Y
  // fallback when the projection collapses.
  auto project = [&](const Vec3& a) -> std::optional<Vec3> {
    if (!finite(a)) return std::nullopt;
    double al = a.norm();
    if (!(al > 0)) return std::nullopt;
    Vec3 au = a / al;
    Vec3 pa = au - zhat * au.dot(zhat);
    if (pa.norm() < kShellNearParallelTol) return std::nullopt;
    return Vec3(pa / pa.norm());
  };
  std::optional<Vec3> xhat;
  if (local_x) xhat = project(*local_x);  // defensive fallback if degenerate
  if (!xhat) xhat = project(Vec3(1, 0, 0));
  if (!xhat) xhat = project(Vec3(0, 1, 0));
  require(xhat.has_value(), "plane frame: cannot build a local frame");
  Vec3 yhat = zhat.cross(*xhat);

  Mat3 r;
  r.row(0) = *xhat;
  r.row(1) = yhat;
  r.row(2) = zhat;
  return r;
}

Mat3 facet_frame(const std::array<Vec3, 3>& p, const std::optional<Vec3>& local_x) {
  for (const Vec3& v : p) require(finite(v), "shell facet: non-finite node coordinate");
  Vec3 v1 = p[1] - p[0];
  Vec3 v2 = p[2] - p[0];
  Vec3 n = v1.cross(v2);
  double two_area = n.norm();
  double char_len = std::max({v1.norm(), v2.norm(), (p[2] - p[1]).norm()});
  require(finite(two_area) && two_area > kShellAreaDegenTol * char_len * char_len,
          "shell facet: degenerate triangle (collinear or coincident nodes)");
  return plane_frame(n, local_x);
}

void principal2(double a, double b, double c, double& p1, double& p2, double& theta) {
  double avg = 0.5 * (a + b);
  double d = 0.5 * (a - b);
  double radius = std::sqrt(d * d + c * c);
  p1 = avg + radius;
  p2 = avg - radius;
  theta = 0.5 * std::atan2(2.0 * c, a - b);
}

void ShellResultants::recompute_principals() {
  principal2(mx, my, mxy, m1, m2, m_theta);
  principal2(nx, ny, nxy, n1, n2, n_theta);
}

namespace {

// Local DOF index of engineering DOF `which` (0..5) at node `node` (0..2).
constexpr int dof_of(int node, int which) { return 6 * node + which; }

// Membrane block [u, v, rz] per node (Allman), in [u1,v1,t1,...,u3,v3,t3] order.
constexpr std::array<int, 9> kMemIdx{dof_of(0, 0), dof_of(0, 1), dof_of(0, 5),
                                     dof_of(1, 0), dof_of(1, 1), dof_of(1, 5),
                                     dof_of(2, 0), dof_of(2, 1), dof_of(2, 5)};
// Bending/shear block [w, rx, ry] per node.
constexpr std::array<int, 9> kBendIdx{dof_of(0, 2), dof_of(0, 3), dof_of(0, 4),
                                      dof_of(1, 2), dof_of(1, 3), dof_of(1, 4),
                                      dof_of(2, 2), dof_of(2, 3), dof_of(2, 4)};

// s * [[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]] — the plane-stress / plate shape.
Mat3 plane_stress(double s, double nu) {
  Mat3 d = Mat3::Zero();
  d(0, 0) = s;
  d(0, 1) = s * nu;
  d(1, 0) = s * nu;
  d(1, 1) = s;
  d(2, 2) = s * (1.0 - nu) / 2.0;
  return d;
}

// Cartesian shape-function derivatives of the linear triangle.
void dn_cartesian(const std::array<std::array<double, 2>, 3>& xy, double dndx[3],
                  double dndy[3]) {
  double x21 = xy[1][0] - xy[0][0], y21 = xy[1][1] - xy[0][1];
  double x31 = xy[2][0] - xy[0][0], y31 = xy[2][1] - xy[0][1];
  double det = x21 * y31 - x31 * y21;  // = 2A
  const double dndr[3] = {-1.0, 1.0, 0.0};
  const double dnds[3] = {-1.0, 0.0, 1.0};
  for (int i = 0; i < 3; ++i) {
    dndx[i] = (y31 * dndr[i] - y21 * dnds[i]) / det;
    dndy[i] = (-x31 * dndr[i] + x21 * dnds[i]) / det;
  }
}

// ---------------------------------------------------------------------------
// Allman membrane = constrained LST (Cook 1986)
// ---------------------------------------------------------------------------

// LST (6-node quadratic triangle) membrane B (3x12) at area coordinates
// (L1, L2, L3), DOF order [u1,v1,...,u3,v3, u4,v4,...,u6,v6] with midside 4 on
// edge 1-2, 5 on edge 2-3, 6 on edge 3-1. Strain [ex, ey, gxy].
Eigen::Matrix<double, 3, 12> lst_b(const std::array<std::array<double, 2>, 3>& xy, double l1,
                                   double l2, double l3) {
  double x1 = xy[0][0], y1 = xy[0][1];
  double x2 = xy[1][0], y2 = xy[1][1];
  double x3 = xy[2][0], y3 = xy[2][1];
  double two_a = x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3 + x1 * y2 - x2 * y1;  // signed 2A
  // dLi/dx = bi / 2A, dLi/dy = ci / 2A.
  double bl[3] = {(y2 - y3) / two_a, (y3 - y1) / two_a, (y1 - y2) / two_a};
  double cl[3] = {(x3 - x2) / two_a, (x1 - x3) / two_a, (x2 - x1) / two_a};
  double L[3] = {l1, l2, l3};
  double dndx[6], dndy[6];
  for (int i = 0; i < 3; ++i) {
    dndx[i] = (4.0 * L[i] - 1.0) * bl[i];
    dndy[i] = (4.0 * L[i] - 1.0) * cl[i];
  }
  // N4 = 4 L1 L2, N5 = 4 L2 L3, N6 = 4 L3 L1.
  const int ea[3] = {0, 1, 2}, eb[3] = {1, 2, 0};
  for (int m = 0; m < 3; ++m) {
    dndx[3 + m] = 4.0 * (L[eb[m]] * bl[ea[m]] + L[ea[m]] * bl[eb[m]]);
    dndy[3 + m] = 4.0 * (L[eb[m]] * cl[ea[m]] + L[ea[m]] * cl[eb[m]]);
  }
  Eigen::Matrix<double, 3, 12> b = Eigen::Matrix<double, 3, 12>::Zero();
  for (int n = 0; n < 6; ++n) {
    b(0, 2 * n) = dndx[n];
    b(1, 2 * n + 1) = dndy[n];
    b(2, 2 * n) = dndy[n];
    b(2, 2 * n + 1) = dndx[n];
  }
  return b;
}

// Allman constraint (Cook 1986): the LST midside displacements slaved to the
// corner [u, v, theta_z] DOFs. Edge i->j with midside m:
//   u_m = (u_i + u_j)/2 + (y_j - y_i)(theta_j - theta_i)/8
//   v_m = (v_i + v_j)/2 - (x_j - x_i)(theta_j - theta_i)/8
// (theta = the continuum spin convention; verified against the exact quadratic
// edge-normal displacement of pure in-plane bending.)
Eigen::Matrix<double, 12, 9> allman_constraint(const std::array<std::array<double, 2>, 3>& xy) {
  Eigen::Matrix<double, 12, 9> t = Eigen::Matrix<double, 12, 9>::Zero();
  for (int n = 0; n < 3; ++n) {
    t(2 * n, 3 * n) = 1.0;
    t(2 * n + 1, 3 * n + 1) = 1.0;
  }
  const int ea[3] = {0, 1, 2}, eb[3] = {1, 2, 0};
  for (int m = 0; m < 3; ++m) {
    int i = ea[m], j = eb[m];
    double xji = xy[static_cast<std::size_t>(j)][0] - xy[static_cast<std::size_t>(i)][0];
    double yji = xy[static_cast<std::size_t>(j)][1] - xy[static_cast<std::size_t>(i)][1];
    int ru = 6 + 2 * m, rv = 7 + 2 * m;
    t(ru, 3 * i) = 0.5;
    t(ru, 3 * j) = 0.5;
    t(ru, 3 * i + 2) = -yji / 8.0;
    t(ru, 3 * j + 2) = yji / 8.0;
    t(rv, 3 * i + 1) = 0.5;
    t(rv, 3 * j + 1) = 0.5;
    t(rv, 3 * i + 2) = xji / 8.0;
    t(rv, 3 * j + 2) = -xji / 8.0;
  }
  return t;
}

// Bending curvature B (3x9), DOF per node [w, rx, ry]; curvature
// [kx, ky, kxy] with kx = d(ry)/dx, ky = -d(rx)/dy, kxy = d(ry)/dy - d(rx)/dx.
Eigen::Matrix<double, 3, 9> bending_b(const std::array<std::array<double, 2>, 3>& xy) {
  double dndx[3], dndy[3];
  dn_cartesian(xy, dndx, dndy);
  Eigen::Matrix<double, 3, 9> b = Eigen::Matrix<double, 3, 9>::Zero();
  for (int i = 0; i < 3; ++i) {
    int tx = 3 * i + 1, ty = 3 * i + 2;
    b(0, ty) = dndx[i];
    b(1, tx) = -dndy[i];
    b(2, ty) = dndy[i];
    b(2, tx) = -dndx[i];
  }
  return b;
}

// MITC3 tying data: the displacement-interpolated covariant shears
//   gamma_r = (w2-w1) + sum Ni (ry_i * x21 - rx_i * y21)
//   gamma_s = (w3-w1) + sum Ni (ry_i * x31 - rx_i * y31)
// evaluated at the tying points A(1/2,0), B(0,1/2), C(1/2,1/2), plus the
// coupling row c = gamma_r^C - gamma_s^C + gamma_s^B - gamma_r^A.
using Row9 = Eigen::Matrix<double, 1, 9>;

Row9 gamma_r_row(const double nsf[3], double x21, double y21) {
  Row9 g = Row9::Zero();
  g(0) = -1.0;  // w1
  g(3) = 1.0;   // w2
  for (int i = 0; i < 3; ++i) {
    g(3 * i + 1) += -nsf[i] * y21;  // rx_i
    g(3 * i + 2) += nsf[i] * x21;   // ry_i
  }
  return g;
}

Row9 gamma_s_row(const double nsf[3], double x31, double y31) {
  Row9 g = Row9::Zero();
  g(0) = -1.0;  // w1
  g(6) = 1.0;   // w3
  for (int i = 0; i < 3; ++i) {
    g(3 * i + 1) += -nsf[i] * y31;
    g(3 * i + 2) += nsf[i] * x31;
  }
  return g;
}

struct ShearSetup {
  Row9 gr_a, c, gs_b;
  double det = 0, x21 = 0, y21 = 0, x31 = 0, y31 = 0;
};

ShearSetup shear_setup(const std::array<std::array<double, 2>, 3>& xy) {
  ShearSetup s;
  s.x21 = xy[1][0] - xy[0][0];
  s.y21 = xy[1][1] - xy[0][1];
  s.x31 = xy[2][0] - xy[0][0];
  s.y31 = xy[2][1] - xy[0][1];
  s.det = s.x21 * s.y31 - s.x31 * s.y21;
  const double n_a[3] = {0.5, 0.5, 0.0};
  const double n_b[3] = {0.5, 0.0, 0.5};
  const double n_c[3] = {0.0, 0.5, 0.5};
  s.gr_a = gamma_r_row(n_a, s.x21, s.y21);
  s.gs_b = gamma_s_row(n_b, s.x31, s.y31);
  Row9 gr_c = gamma_r_row(n_c, s.x21, s.y21);
  Row9 gs_c = gamma_s_row(n_c, s.x31, s.y31);
  s.c = gr_c - gs_c + s.gs_b - s.gr_a;
  return s;
}

// Cartesian shear B (2x9) at the natural point (r,s) from the assumed covariant
// field gamma_rt = gamma_r^A + c*s, gamma_st = gamma_s^B - c*r, mapped through
// J^{-1}.
Eigen::Matrix<double, 2, 9> shear_b_at(double r, double s, const ShearSetup& su) {
  Eigen::Matrix<double, 2, 9> b;
  for (int i = 0; i < 9; ++i) {
    double grt = su.gr_a(i) + su.c(i) * s;
    double gst = su.gs_b(i) - su.c(i) * r;
    b(0, i) = (su.y31 * grt - su.y21 * gst) / su.det;
    b(1, i) = (-su.x31 * grt + su.x21 * gst) / su.det;
  }
  return b;
}

// Block-diagonal 18x18 transform T (R repeated over the six triples), so that
// d_local = T * d_global and K_global = T' K_local T.
Mat18 transform18(const Mat3& r) {
  Mat18 t = Mat18::Zero();
  for (int blk = 0; blk < 6; ++blk) t.block<3, 3>(3 * blk, 3 * blk) = r;
  return t;
}

}  // namespace

// ===========================================================================
// ShellElement
// ===========================================================================

ShellElement::ShellElement(const std::array<Vec3, 3>& coords, double t, double E, double nu,
                           const std::optional<Vec3>& local_x, const ShellEnrich& enrich)
    : coords_(coords), t_(t), en_(enrich) {
  for (const Vec3& p : coords) require(finite(p), "shell element: non-finite node coordinate");
  require(finite(t) && t > 0, "shell element: thickness must be finite and positive");
  require(finite(E) && E > 0, "shell element: E must be finite and positive");
  require(finite(nu) && nu > -1.0 && nu < 0.5, "shell element: nu out of range (-1, 0.5)");

  r_ = facet_frame(coords, local_x);  // also performs the degeneracy check
  Vec3 v1 = coords[1] - coords[0];
  Vec3 v2 = coords[2] - coords[0];
  area_ = 0.5 * v1.cross(v2).norm();

  Vec3 xhat = r_.row(0), yhat = r_.row(1);
  std::array<std::array<double, 2>, 3> xy{
      std::array<double, 2>{0.0, 0.0},
      std::array<double, 2>{v1.dot(xhat), v1.dot(yhat)},
      std::array<double, 2>{v2.dot(xhat), v2.dot(yhat)}};

  double G = E / (2.0 * (1.0 + nu));
  xy_ = xy;
  G_ = G;
  gd_ = G;
  gs_ = G * Eigen::Matrix2d::Identity();
  dm_ = plane_stress(E / (1.0 - nu * nu), nu);
  db_ = plane_stress(E * t * t * t / (12.0 * (1.0 - nu * nu)), nu);
  ks_coef_ = kShellShearCorrection * G * t;

  Eigen::Matrix<double, 12, 9> ta = allman_constraint(xy);
  bm_ = lst_b(xy, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0) * ta;  // centroid, for recovery
  bb_ = bending_b(xy);
  ShearSetup su = shear_setup(xy);

  Mat18 kl = Mat18::Zero();
  // Membrane: Allman (constrained LST), 3-point mid-edge rule — exact for the
  // quadratic B'DB of the linear-in-(r,s) B_Allman.
  Eigen::Matrix<double, 9, 9> km = Eigen::Matrix<double, 9, 9>::Zero();
  const double mgp[3][3] = {{0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}};
  for (const auto& gp : mgp) {
    Eigen::Matrix<double, 3, 9> bmg = lst_b(xy, gp[0], gp[1], gp[2]) * ta;
    km += (t * area_ / 3.0) * (bmg.transpose() * dm_ * bmg);
  }
  // Spurious-mode stabilisation: penalty on (theta_bar - omega_CST), the ONLY
  // fictitious stiffness left (see shell.h). g d = mean drilling rotation
  // minus the continuum spin of the corner (linear) field.
  {
    double dndx[3], dndy[3];
    dn_cartesian(xy, dndx, dndy);
    Eigen::Matrix<double, 1, 9> g = Eigen::Matrix<double, 1, 9>::Zero();
    for (int i = 0; i < 3; ++i) {
      g(3 * i) = 0.5 * dndy[i];       // -omega: -1/2 dv/dx + 1/2 du/dy -> u part
      g(3 * i + 1) = -0.5 * dndx[i];  // v part
      g(3 * i + 2) = 1.0 / 3.0;       // theta_bar
    }
    km += (kShellDrillStab * G * t * area_) * (g.transpose() * g);
  }
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 9; ++j) kl(kMemIdx[static_cast<std::size_t>(i)],
                                   kMemIdx[static_cast<std::size_t>(j)]) += km(i, j);

  // Bending: K_b = A * Bb' Db Bb; shear: K_s = sum_gp (A/3) Bs' (ks I) Bs
  // at the three mid-edge points (exact for the quadratic Bs'Bs).
  Eigen::Matrix<double, 9, 9> kpb = area_ * (bb_.transpose() * db_ * bb_);
  const double gps[3][2] = {{0.5, 0.0}, {0.0, 0.5}, {0.5, 0.5}};
  for (const auto& gp : gps) {
    Eigen::Matrix<double, 2, 9> bs = shear_b_at(gp[0], gp[1], su);
    kpb += (area_ / 3.0) * ks_coef_ * (bs.transpose() * bs);
  }
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 9; ++j) kl(kBendIdx[static_cast<std::size_t>(i)],
                                   kBendIdx[static_cast<std::size_t>(j)]) += kpb(i, j);

  tc_ = t;
  tcorner_ = {t, t, t};
  if (!en_.plain()) {
    // ---- the enriched facet: everything integrated by the 7-point
    // Dunavant rule (exact to degree 5: t linear times the quadratic B'DB
    // of the membrane and shear, t^3 of the bending, the mass) ----
    for (int c = 0; c < 3; ++c) {
      const double tv = en_.tn[static_cast<std::size_t>(c)];
      require(finite(tv) && tv >= 0, "shell element: corner thickness must be >= 0");
      tcorner_[static_cast<std::size_t>(c)] = tv > 0 ? tv : t;
    }
    tc_ = (tcorner_[0] + tcorner_[1] + tcorner_[2]) / 3;
    require(finite(en_.ecc) && finite(en_.kbed) && en_.kbed >= 0, "shell element: invalid offset or bed");
    Mat3 Q = dm_;
    gs_ = G * Eigen::Matrix2d::Identity();
    double Gd = G;
    if (en_.ortho) {
      const ShellOrtho& o = *en_.ortho;
      require(o.Ex > 0 && o.Ey > 0 && o.Gxy > 0 && o.Gxz > 0 && o.Gyz > 0 && finite(o.nu_xy) && finite(o.angle),
              "shell element: orthotropic moduli must be positive");
      const double nyx = o.nu_xy * o.Ey / o.Ex, den = 1 - o.nu_xy * nyx;
      require(den > 0, "shell element: orthotropic Poisson ratios not admissible");
      Mat3 Qm;
      Qm << o.Ex / den, o.nu_xy * o.Ey / den, 0, o.nu_xy * o.Ey / den, o.Ey / den, 0, 0, 0, o.Gxy;
      // material axes at `angle`: Q_local = T' Qm T, T the strain rotation
      // (engineering shear) from the local to the material axes
      const double c = std::cos(o.angle), sn = std::sin(o.angle);
      Mat3 Te;
      Te << c * c, sn * sn, c * sn, sn * sn, c * c, -c * sn, -2 * c * sn, 2 * c * sn, c * c - sn * sn;
      Q = Te.transpose() * Qm * Te;
      Eigen::Matrix2d R2;
      R2 << c, sn, -sn, c;
      gs_ = R2.transpose() * Eigen::Vector2d(o.Gxz, o.Gyz).asDiagonal() * R2;
      Gd = o.Gxy;
    }
    dm_ = Q;
    db_ = (tc_ * tc_ * tc_ / 12.0) * Q;  // centroid value (recovery); the stiffness integrates t^3
    if (en_.laminate) {
      // classical laminate theory about the facet mid-plane (see model.h)
      const ShellLaminate& lm = *en_.laminate;
      require(!en_.ortho, "shell element: a laminate replaces the orthotropic material (not both)");
      require(en_.tn[0] == 0 && en_.tn[1] == 0 && en_.tn[2] == 0, "shell element: a laminate takes a uniform thickness");
      require(!lm.plies.empty(), "shell element: a laminate needs plies");
      double tsum = 0;
      for (const ShellPly& p : lm.plies) {
        require(p.t > 0 && p.E1 > 0 && p.E2 > 0 && p.G12 > 0 && p.G13 > 0 && p.G23 > 0 && finite(p.nu12) && finite(p.angle),
                "shell element: ply thickness and moduli must be positive");
        tsum += p.t;
      }
      require(std::fabs(tsum - t) <= 1e-9 * t, "shell element: the plies must add up to the facet thickness");
      require(finite(lm.z_bottom), "shell element: non-finite laminate bottom");
      double z0 = lm.centred ? -tsum / 2 : lm.z_bottom;
      lam_A_.setZero();
      lam_B_.setZero();
      lam_D_.setZero();
      Eigen::Matrix2d Hs = Eigen::Matrix2d::Zero();
      for (const ShellPly& p : lm.plies) {
        const double n21 = p.nu12 * p.E2 / p.E1, den = 1 - p.nu12 * n21;
        require(den > 0, "shell element: ply Poisson ratios not admissible");
        Mat3 Qm;
        Qm << p.E1 / den, p.nu12 * p.E2 / den, 0, p.nu12 * p.E2 / den, p.E2 / den, 0, 0, 0, p.G12;
        const double c = std::cos(p.angle), sn = std::sin(p.angle);
        Mat3 Te;
        Te << c * c, sn * sn, c * sn, sn * sn, c * c, -c * sn, -2 * c * sn, 2 * c * sn, c * c - sn * sn;
        const Mat3 Qb = Te.transpose() * Qm * Te;
        Eigen::Matrix2d R2;
        R2 << c, sn, -sn, c;
        const double z1 = z0 + p.t;
        lam_A_ += (z1 - z0) * Qb;
        lam_B_ += ((z1 * z1 - z0 * z0) / 2) * Qb;
        lam_D_ += ((z1 * z1 * z1 - z0 * z0 * z0) / 3) * Qb;
        Hs += p.t * (R2.transpose() * Eigen::Vector2d(p.G13, p.G23).asDiagonal() * R2);
        z0 = z1;
      }
      lam_ = true;
      gs_ = Hs / tsum;
      Gd = 0;  // the drilling stabilisation: the plies' mean in-plane shear modulus
      for (const ShellPly& p : lm.plies) Gd += p.G12 * p.t / tsum;
      dm_ = lam_A_ / tsum;
      db_ = lam_D_;
    }
    static const double dw[7] = {0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                                 0.125939180544827, 0.125939180544827, 0.125939180544827};
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, a2 = 0.797426985353087, b2 = 0.101286507323456;
    const double dl[7][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1},
                             {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    Eigen::Matrix<double, 9, 9> km = Eigen::Matrix<double, 9, 9>::Zero(), kb = km, kmb = km;
    const ShearSetup su2 = shear_setup(xy);
    for (int g = 0; g < 7; ++g) {
      const double tg = dl[g][0] * tcorner_[0] + dl[g][1] * tcorner_[1] + dl[g][2] * tcorner_[2];
      const double w = dw[g] * area_;
      const Eigen::Matrix<double, 3, 9> bmg = lst_b(xy, dl[g][0], dl[g][1], dl[g][2]) * ta;
      if (lam_) {
        km += w * (bmg.transpose() * lam_A_ * bmg);
        kb += w * (bb_.transpose() * lam_D_ * bb_);
        kmb += w * (bmg.transpose() * lam_B_ * bb_);
      } else {
        km += (w * tg) * (bmg.transpose() * Q * bmg);
        kb += (w * tg * tg * tg / 12.0) * (bb_.transpose() * Q * bb_);
      }
      const Eigen::Matrix<double, 2, 9> bs = shear_b_at(dl[g][1], dl[g][2], su2);
      kb += (w * kShellShearCorrection * tg) * (bs.transpose() * gs_ * bs);
    }
    {
      double dndx[3], dndy[3];
      dn_cartesian(xy, dndx, dndy);
      Eigen::Matrix<double, 1, 9> gd = Eigen::Matrix<double, 1, 9>::Zero();
      for (int i = 0; i < 3; ++i) {
        gd(3 * i) = 0.5 * dndy[i];
        gd(3 * i + 1) = -0.5 * dndx[i];
        gd(3 * i + 2) = 1.0 / 3.0;
      }
      km += (kShellDrillStab * Gd * tc_ * area_) * (gd.transpose() * gd);
      gd_ = Gd;
    }
    kl.setZero();
    for (int i = 0; i < 9; ++i)
      for (int j = 0; j < 9; ++j) {
        kl(kMemIdx[static_cast<std::size_t>(i)], kMemIdx[static_cast<std::size_t>(j)]) += km(i, j);
        kl(kBendIdx[static_cast<std::size_t>(i)], kBendIdx[static_cast<std::size_t>(j)]) += kb(i, j);
        // laminate membrane-bending coupling (B)
        kl(kMemIdx[static_cast<std::size_t>(i)], kBendIdx[static_cast<std::size_t>(j)]) += kmb(i, j);
        kl(kBendIdx[static_cast<std::size_t>(j)], kMemIdx[static_cast<std::size_t>(i)]) += kmb(i, j);
      }
    ks_coef_ = kShellShearCorrection * tc_;  // times gs_ in recovery
    // the mid-surface offset e along the normal: u_m = u + e theta_y,
    // v_m = v - e theta_x (local), a rigid arm
    if (en_.ecc != 0)
      for (int nd = 0; nd < 3; ++nd) {
        te_(dof_of(nd, 0), dof_of(nd, 4)) = en_.ecc;
        te_(dof_of(nd, 1), dof_of(nd, 3)) = -en_.ecc;
      }
    kl = te_.transpose() * kl * te_;
    // Winkler bed on w: k (A / 12) (1 + delta_ij), exact for the linear w
    if (en_.kbed > 0)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) kl(dof_of(a, 2), dof_of(b, 2)) += en_.kbed * area_ * (a == b ? 2.0 : 1.0) / 12.0;
  }
  Mat18 T = transform18(r_);
  ke_ = T.transpose() * kl * T;
  ke_ = 0.5 * (ke_ + ke_.transpose().eval());  // kill round-off asymmetry

  bsc_ = shear_b_at(1.0 / 3.0, 1.0 / 3.0, su);  // centroid, for recovery
  {
    double dndx[3], dndy[3];
    dn_cartesian(xy, dndx, dndy);
    for (int i = 0; i < 3; ++i) {
      gw_(0, i) = dndx[i];
      gw_(1, i) = dndy[i];
    }
  }
}

Mat18 ShellElement::local_transform() const { return transform18(r_); }

Vec3 ShellElement::strain_at(const Vec18& d_global, const std::array<double, 3>& l, double z) const {
  Vec18 d = transform18(r_) * d_global;
  if (en_.ecc != 0) d = te_ * d;
  Eigen::Matrix<double, 9, 1> um, ub;
  for (int i = 0; i < 9; ++i) {
    um(i) = d(kMemIdx[static_cast<std::size_t>(i)]);
    ub(i) = d(kBendIdx[static_cast<std::size_t>(i)]);
  }
  const Eigen::Matrix<double, 3, 9> bm = lst_b(xy_, l[0], l[1], l[2]) * allman_constraint(xy_);
  return bm * um + z * (bb_ * ub);
}

ShellNlOps ShellElement::nl_ops() const {
  ShellNlOps o;
  o.area = area_;
  o.t = tc_;
  o.mem_idx = kMemIdx;
  o.bend_idx = kBendIdx;
  o.T = transform18(r_);
  const Eigen::Matrix<double, 12, 9> ta = allman_constraint(xy_);
  const ShearSetup su = shear_setup(xy_);
  const bool varies = tcorner_[0] != tcorner_[1] || tcorner_[1] != tcorner_[2];
  if (!varies) {
    // the mid-edge points; in MITC3 natural coordinates (r, s) = (L2, L3)
    const double mgp[3][3] = {{0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}};
    for (const auto& gp : mgp) {
      o.bm.push_back(lst_b(xy_, gp[0], gp[1], gp[2]) * ta);
      o.bs.push_back(shear_b_at(gp[1], gp[2], su));
      o.w.push_back(area_ / 3.0);
      o.tp.push_back(tc_);
    }
  } else {
    static const double dw[7] = {0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                                 0.125939180544827, 0.125939180544827, 0.125939180544827};
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, a2 = 0.797426985353087, b2 = 0.101286507323456;
    const double dl[7][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1},
                             {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    for (int g = 0; g < 7; ++g) {
      o.bm.push_back(lst_b(xy_, dl[g][0], dl[g][1], dl[g][2]) * ta);
      o.bs.push_back(shear_b_at(dl[g][1], dl[g][2], su));
      o.w.push_back(dw[g] * area_);
      o.tp.push_back(dl[g][0] * tcorner_[0] + dl[g][1] * tcorner_[1] + dl[g][2] * tcorner_[2]);
    }
  }
  if (!en_.plain()) {
    // the plies' elastic matrices (local axes), as the stiffness above
    auto ply_C = [](double E1, double E2, double nu12, double G12, double G13, double G23, double angle) {
      const double n21 = nu12 * E2 / E1, den = 1 - nu12 * n21;
      Mat3 Qm;
      Qm << E1 / den, nu12 * E2 / den, 0, nu12 * E2 / den, E2 / den, 0, 0, 0, G12;
      const double c = std::cos(angle), sn = std::sin(angle);
      Mat3 Te;
      Te << c * c, sn * sn, c * sn, sn * sn, c * c, -c * sn, -2 * c * sn, 2 * c * sn, c * c - sn * sn;
      Eigen::Matrix2d R2;
      R2 << c, sn, -sn, c;
      Eigen::Matrix<double, 5, 5> C = Eigen::Matrix<double, 5, 5>::Zero();
      C.topLeftCorner<3, 3>() = Te.transpose() * Qm * Te;
      C.bottomRightCorner<2, 2>() = R2.transpose() * Eigen::Vector2d(G13, G23).asDiagonal() * R2;
      return C;
    };
    if (en_.laminate) {
      double z0 = en_.laminate->centred ? -t_ / 2 : en_.laminate->z_bottom;
      for (const ShellPly& p : en_.laminate->plies) {
        o.plies.push_back({z0 / t_, (z0 + p.t) / t_, ply_C(p.E1, p.E2, p.nu12, p.G12, p.G13, p.G23, p.angle)});
        z0 += p.t;
      }
      o.laminate = true;
    } else if (en_.ortho) {
      const ShellOrtho& w = *en_.ortho;
      o.plies.push_back({-0.5, 0.5, ply_C(w.Ex, w.Ey, w.nu_xy, w.Gxy, w.Gxz, w.Gyz, w.angle)});
    }
    o.te = te_;
    o.offset = en_.ecc != 0;
    if (en_.kbed > 0)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) o.k_bed(dof_of(a, 2), dof_of(b, 2)) = en_.kbed * area_ * (a == b ? 2.0 : 1.0) / 12.0;
  }
  o.bb = bb_;
  double dndx[3], dndy[3];
  dn_cartesian(xy_, dndx, dndy);
  Eigen::Matrix<double, 1, 9> g = Eigen::Matrix<double, 1, 9>::Zero();
  for (int i = 0; i < 3; ++i) {
    g(3 * i) = 0.5 * dndy[i];
    g(3 * i + 1) = -0.5 * dndx[i];
    g(3 * i + 2) = 1.0 / 3.0;
  }
  o.k_drill = (kShellDrillStab * gd_ * tc_ * area_) * (g.transpose() * g);
  // the elastic transverse shear at the same points (exact for the
  // quadratic Bs'Bs times the linear t)
  o.k_shear.setZero();
  for (std::size_t p = 0; p < o.w.size(); ++p)
    o.k_shear += (o.w[p] * kShellShearCorrection * o.tp[p]) * (o.bs[p].transpose() * gs_ * o.bs[p]);
  return o;
}

Mat18 ShellElement::me(double rho, bool lumped) const {
  require(finite(rho) && rho >= 0, "shell element: rho must be finite and >= 0");
  Mat18 ml = Mat18::Zero();
  if (!en_.plain() && !lumped) {
    // int rho t Ni Nj and rho t^3/12 Ni Nj with the linear t: the monomial
    // integrals of the area coordinates (int La Lb Lc = 2A a!b!c!/(a+b+c+2)!)
    // by the 7-point rule (exact to degree 5)
    static const double dw[7] = {0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                                 0.125939180544827, 0.125939180544827, 0.125939180544827};
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, a2 = 0.797426985353087, b2 = 0.101286507323456;
    const double dl[7][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1},
                             {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    for (int g = 0; g < 7; ++g) {
      const double tg = dl[g][0] * tcorner_[0] + dl[g][1] * tcorner_[1] + dl[g][2] * tcorner_[2];
      const double w = dw[g] * area_ * rho;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
          const double nn = dl[g][a] * dl[g][b] * w;
          for (int c = 0; c < 3; ++c) ml(dof_of(a, c), dof_of(b, c)) += nn * tg;
          ml(dof_of(a, 3), dof_of(b, 3)) += nn * tg * tg * tg / 12.0;
          ml(dof_of(a, 4), dof_of(b, 4)) += nn * tg * tg * tg / 12.0;
        }
    }
    ml = te_.transpose() * ml * te_;
    Mat18 T = transform18(r_);
    Mat18 m = T.transpose() * ml * T;
    return 0.5 * (m + m.transpose().eval());
  }
  if (!en_.plain()) {
    // lumped: the corner tributaries of the linear thickness
    const double* tc = tcorner_.data();
    for (int n = 0; n < 3; ++n) {
      const double tt = (2 * tc[n] + tc[(n + 1) % 3] + tc[(n + 2) % 3]) / 4;
      for (int c = 0; c < 3; ++c) ml(dof_of(n, c), dof_of(n, c)) += rho * tt * area_ / 3.0;
      ml(dof_of(n, 3), dof_of(n, 3)) += rho * tt * tt * tt / 12.0 * area_ / 3.0;
      ml(dof_of(n, 4), dof_of(n, 4)) += rho * tt * tt * tt / 12.0 * area_ / 3.0;
    }
    ml = te_.transpose() * ml * te_;
    Mat18 T = transform18(r_);
    return T.transpose() * ml * T;
  }
  double mt = rho * t_ * area_;                 // translational mass
  double mr = rho * (t_ * t_ * t_ / 12.0) * area_;  // rotary inertia (rx, ry)
  if (lumped) {
    for (int n = 0; n < 3; ++n) {
      for (int c = 0; c < 3; ++c) ml(dof_of(n, c), dof_of(n, c)) += mt / 3.0;
      ml(dof_of(n, 3), dof_of(n, 3)) += mr / 3.0;
      ml(dof_of(n, 4), dof_of(n, 4)) += mr / 3.0;
    }
  } else {
    // Exact consistent matrix of a linear triangle: (A/12)(1 + delta_ij)
    // per scalar field => coefficients m/12 * [[2,1,1],[1,2,1],[1,1,2]].
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        double w = (a == b ? 2.0 : 1.0) / 12.0;
        for (int c = 0; c < 3; ++c) ml(dof_of(a, c), dof_of(b, c)) += mt * w;
        ml(dof_of(a, 3), dof_of(b, 3)) += mr * w;
        ml(dof_of(a, 4), dof_of(b, 4)) += mr * w;
      }
  }
  Mat18 T = transform18(r_);
  Mat18 m = T.transpose() * ml * T;
  return 0.5 * (m + m.transpose().eval());
}

Mat18 ShellElement::kg(double nx, double ny, double nxy) const {
  require(finite(nx) && finite(ny) && finite(nxy), "shell element: non-finite membrane state");
  // Kg_w = A * grad' [n] grad — grad(w) is CONSTANT (linear w), 1-point exact.
  Eigen::Matrix<double, 2, 2> n;
  n << nx, nxy, nxy, ny;
  Eigen::Matrix<double, 3, 3> kw = area_ * (gw_.transpose() * n * gw_);
  Mat18 kl = Mat18::Zero();
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) kl(dof_of(a, 2), dof_of(b, 2)) += kw(a, b);
  if (en_.ecc != 0) {
    // the offset is a rigid arm r = e n from each node to the mid-surface:
    // the in-plane force P the mid-surface hands to the arm (minus the
    // element's membrane nodal force of the prestress, int Bm' n over the
    // facet) adds the second variation of -P . (theta x (theta x r)) / 2,
    // (P . r) I - (P r' + r P') / 2, on the node's rotations — the term of
    // the rigid links (eigsolve.cpp); with P in the plane and r normal to it
    // it couples the drilling rotation with the bending ones
    const Eigen::Matrix<double, 12, 9> ta = allman_constraint(xy_);
    const double mgp[3][3] = {{0.5, 0.5, 0.0}, {0.0, 0.5, 0.5}, {0.5, 0.0, 0.5}};
    Eigen::Matrix<double, 9, 1> fm = Eigen::Matrix<double, 9, 1>::Zero();
    for (const auto& gp : mgp)
      fm += (area_ / 3.0) * (lst_b(xy_, gp[0], gp[1], gp[2]) * ta).transpose() * Vec3(nx, ny, nxy);
    const Vec3 r(0, 0, en_.ecc);
    for (int k = 0; k < 3; ++k) {
      const Vec3 P(-fm(3 * k), -fm(3 * k + 1), 0.0);
      const Mat3 B = P.dot(r) * Mat3::Identity() - 0.5 * (P * r.transpose() + r * P.transpose());
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) kl(dof_of(k, 3 + a), dof_of(k, 3 + b)) += B(a, b);
    }
  }
  Mat18 T = transform18(r_);
  Mat18 g = T.transpose() * kl * T;
  return 0.5 * (g + g.transpose().eval());
}

Vec18 ShellElement::pressure_load_global(const Vec3& q) const {
  Vec3 f = q * (area_ / 3.0);
  Vec18 out = Vec18::Zero();
  for (int node = 0; node < 3; ++node) out.segment<3>(6 * node) = f;
  if (en_.ecc != 0) {
    // applied on the offset mid-surface: f_nodes = Te' f_mid (local)
    const Mat18 T = transform18(r_);
    out = T.transpose() * (te_.transpose() * (T * out));
  }
  return out;
}

Vec18 ShellElement::pressure_load_normal(double p) const {
  return pressure_load_global(normal() * p);
}

Vec18 ShellElement::hydrostatic_load(double gamma, const Vec3& x0, const Vec3& down, double side) const {
  const Vec3 dn = down.normalized();
  struct P {
    Eigen::Vector3d l;  // barycentric
    double d;           // depth
  };
  std::vector<P> poly;
  for (int i = 0; i < 3; ++i) poly.push_back({Eigen::Vector3d::Unit(i), (coords_[static_cast<std::size_t>(i)] - x0).dot(dn)});
  // clip to the wet part d >= 0 (Sutherland-Hodgman, d linear)
  std::vector<P> wet;
  for (std::size_t k = 0; k < poly.size(); ++k) {
    const P& a = poly[k];
    const P& b = poly[(k + 1) % poly.size()];
    if (a.d >= 0) wet.push_back(a);
    if ((a.d >= 0) != (b.d >= 0)) {
      const double t = a.d / (a.d - b.d);
      wet.push_back({a.l + t * (b.l - a.l), 0.0});
    }
  }
  Eigen::Vector3d I = Eigen::Vector3d::Zero();  // int d N_i dA
  for (std::size_t k = 1; k + 1 < wet.size(); ++k) {
    const P* v[3] = {&wet[0], &wet[k], &wet[k + 1]};
    const Eigen::Vector2d e1(v[1]->l(1) - v[0]->l(1), v[1]->l(2) - v[0]->l(2));
    const Eigen::Vector2d e2(v[2]->l(1) - v[0]->l(1), v[2]->l(2) - v[0]->l(2));
    const double As = area_ * std::fabs(e1(0) * e2(1) - e1(1) * e2(0));
    const double sd = v[0]->d + v[1]->d + v[2]->d;
    for (int i = 0; i < 3; ++i) {
      double sdl = 0, sl = 0;
      for (const P* q : v) {
        sdl += q->d * q->l(i);
        sl += q->l(i);
      }
      I(i) += As / 12.0 * (sdl + sd * sl);  // int f g = A/12 (sum f g + sum f sum g), f, g linear
    }
  }
  Vec18 out = Vec18::Zero();
  for (int i = 0; i < 3; ++i) out.segment<3>(6 * i) = (gamma * side * I(i)) * normal();
  if (en_.ecc != 0) {
    const Mat18 T = transform18(r_);
    out = T.transpose() * (te_.transpose() * (T * out));
  }
  return out;
}

Vec18 ShellElement::pressure_load_projected(const Vec3& q) const {
  // Plan-area convention: the true-area force is q * A_plan = q * A * |n.z|.
  return pressure_load_global(q * std::fabs(normal().z()));
}

Vec18 ShellElement::self_weight_load(double rho, const Vec3& g) const {
  if (en_.plain()) return pressure_load_global(g * (rho * t_));
  // int rho t Ni dA = rho A (2 ti + tj + tk) / 12 for the linear t
  Vec18 out = Vec18::Zero();
  const double* tc = tcorner_.data();
  for (int n = 0; n < 3; ++n)
    out.segment<3>(6 * n) = g * (rho * area_ * (2 * tc[n] + tc[(n + 1) % 3] + tc[(n + 2) % 3]) / 12.0);
  if (en_.ecc != 0) {
    const Mat18 T = transform18(r_);
    out = T.transpose() * (te_.transpose() * (T * out));
  }
  return out;
}

Vec18 ShellElement::thermal_load(const Vec3& eps0, const Vec3& kappa0) const {
  require(!lam_, "shell element: a thermal load on a laminate is not supported (per-ply dilations)");
  Vec18 f_local = Vec18::Zero();
  if (!en_.plain()) {
    // the enriched facet: f_m = int t Bm' Q eps0, f_b = int t^3/12 Bb' Q
    // kappa0 (the 7-point rule, exact), carried to the nodes through the
    // offset (f = Te' f_mid)
    const Eigen::Matrix<double, 12, 9> ta = allman_constraint(xy_);
    static const double dw[7] = {0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                                 0.125939180544827, 0.125939180544827, 0.125939180544827};
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, a2 = 0.797426985353087, b2 = 0.101286507323456;
    const double dl[7][3] = {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {a1, b1, b1}, {b1, a1, b1}, {b1, b1, a1},
                             {a2, b2, b2}, {b2, a2, b2}, {b2, b2, a2}};
    Eigen::Matrix<double, 9, 1> fm = Eigen::Matrix<double, 9, 1>::Zero();
    double t3 = 0;
    for (int g = 0; g < 7; ++g) {
      const double tg = dl[g][0] * tcorner_[0] + dl[g][1] * tcorner_[1] + dl[g][2] * tcorner_[2];
      const Eigen::Matrix<double, 3, 9> bmg = lst_b(xy_, dl[g][0], dl[g][1], dl[g][2]) * ta;
      fm += (dw[g] * area_ * tg) * (bmg.transpose() * (dm_ * eps0));
      t3 += dw[g] * area_ * tg * tg * tg / 12.0;
    }
    const Eigen::Matrix<double, 9, 1> fb = t3 * (bb_.transpose() * (dm_ * kappa0));
    for (int i = 0; i < 9; ++i) {
      f_local(kMemIdx[static_cast<std::size_t>(i)]) += fm(i);
      f_local(kBendIdx[static_cast<std::size_t>(i)]) += fb(i);
    }
    f_local = te_.transpose() * f_local;
    return transform18(r_).transpose() * f_local;
  }
  // Membrane f_m = (t*A) Bm' (Dm eps0) on the local [u, v, rz] DOFs. B_Allman
  // is linear in (r, s), so its exact area integral is A times the centroid
  // value — the same bm_ recovery uses.
  Eigen::Matrix<double, 9, 1> fm = (t_ * area_) * (bm_.transpose() * (dm_ * eps0));
  for (int i = 0; i < 9; ++i) f_local(kMemIdx[static_cast<std::size_t>(i)]) += fm(i);
  // Bending f_b = A Bb' (Db kappa0) on the local [w, rx, ry] DOFs.
  Eigen::Matrix<double, 9, 1> fb = area_ * (bb_.transpose() * (db_ * kappa0));
  for (int i = 0; i < 9; ++i) f_local(kBendIdx[static_cast<std::size_t>(i)]) += fb(i);
  return transform18(r_).transpose() * f_local;
}

ShellResultants ShellElement::recover(const Vec18& d_global, const Vec3& eps0,
                                      const Vec3& kappa0) const {
  Vec18 d_local = transform18(r_) * d_global;
  if (en_.ecc != 0) d_local = te_ * d_local;  // the mid-surface's motion
  Eigen::Matrix<double, 9, 1> dm_vec;
  for (int i = 0; i < 9; ++i) dm_vec(i) = d_local(kMemIdx[static_cast<std::size_t>(i)]);
  Eigen::Matrix<double, 9, 1> db_vec;
  for (int i = 0; i < 9; ++i) db_vec(i) = d_local(kBendIdx[static_cast<std::size_t>(i)]);

  ShellResultants out;
  if (lam_) {
    // the laminate: N = A e + B k, M = B e + D k (no thermal state)
    const Vec3 e = bm_ * dm_vec, k = bb_ * db_vec;
    const Vec3 n = lam_A_ * e + lam_B_ * k, m = lam_B_ * e + lam_D_ * k;
    out.nx = n(0);
    out.ny = n(1);
    out.nxy = n(2);
    out.mx = m(0);
    out.my = m(1);
    out.mxy = m(2);
    const Eigen::Matrix<double, 2, 1> gamma = gs_ * (bsc_ * db_vec);
    out.qx = ks_coef_ * gamma(0);
    out.qy = ks_coef_ * gamma(1);
    out.recompute_principals();
    return out;
  }
  // Membrane: n = t * Dm * (eps - eps0), eps = Allman strain at the centroid.
  Vec3 n = dm_ * (Vec3(bm_ * dm_vec) - eps0);
  out.nx = tc_ * n(0);
  out.ny = tc_ * n(1);
  out.nxy = tc_ * n(2);
  // Bending: m = Db * (kappa - kappa0).
  Vec3 m = db_ * (Vec3(bb_ * db_vec) - kappa0);
  out.mx = m(0);
  out.my = m(1);
  out.mxy = m(2);
  // Transverse shear at the centroid (a uniform temperature induces none).
  Eigen::Matrix<double, 2, 1> gamma = bsc_ * db_vec;
  if (!en_.plain()) gamma = gs_ * gamma;  // ks_coef_ = k t there, the moduli in gs_
  out.qx = ks_coef_ * gamma(0);
  out.qy = ks_coef_ * gamma(1);
  out.recompute_principals();
  return out;
}

// ===========================================================================
// Model-side plumbing
// ===========================================================================

std::vector<ShellElem> shell_elems(const Model& model) {
  std::vector<ShellElem> out;
  out.reserve(model.shells.size());
  // Per-slab local-x override: a ShellTri inherits its owning Slab's axis so the
  // recovered mx/my/mxy resolve along the authored direction. A standalone shell
  // (slab id matching no Slab record) uses the default rule.
  std::unordered_map<Id, std::optional<Vec3>> slab_axis;
  for (const Slab& s : model.slabs) slab_axis.emplace(s.id, s.local_x);

  // Build the kernel, re-raising any typed failure with the shell id attached.
  auto build = [](const std::string& at, const std::array<Vec3, 3>& coords, double t,
                  const Material& mat, const std::optional<Vec3>& axis, const ShellEnrich& en) -> ShellElement {
    try {
      return ShellElement(coords, t, mat.E, mat.nu, axis, en);
    } catch (const Error& err) {
      fail(at + ": " + err.what());
    }
  };

  for (std::size_t i = 0; i < model.shells.size(); ++i) {
    const ShellTri& sh = model.shells[i];
    const std::string at = "shell " + std::to_string(sh.id);
    std::array<int, 3> nodes{};
    std::array<Vec3, 3> coords;
    for (int k = 0; k < 3; ++k) {
      Id nid = sh.n[static_cast<std::size_t>(k)];
      nodes[static_cast<std::size_t>(k)] = model.node_index(nid);
      coords[static_cast<std::size_t>(k)] = model.node(nid).x;
    }
    const Material& mat = model.material(sh.material);
    std::optional<Vec3> axis;
    auto it = slab_axis.find(sh.slab);
    if (it != slab_axis.end()) axis = it->second;
    out.push_back(ShellElem{i, sh.id, sh.slab, nodes, mat.rho, mat.alpha,
                            build(at, coords, sh.t, mat, axis, sh.enrich)});
  }
  return out;
}

std::vector<Id> nodes_on_segment(const Model& model, const Vec3& a, const Vec3& b, double tol) {
  require(finite(a) && finite(b), "line support: non-finite segment endpoint");
  require(finite(tol) && tol >= 0, "line support: non-finite / negative tolerance");
  Vec3 d = b - a;
  double len2 = d.squaredNorm();
  require(len2 > 1e-24, "line support: zero-length segment (a line support must have two "
                        "distinct endpoints)");
  std::vector<Id> out;
  for (const Node& n : model.nodes) {
    double s = (n.x - a).dot(d) / len2;
    s = std::min(std::max(s, 0.0), 1.0);
    if ((n.x - (a + s * d)).norm() <= tol) out.push_back(n.id);
  }
  return out;
}

std::vector<Id> add_line_support(Model& model, const Vec3& a, const Vec3& b,
                                 const std::array<bool, 6>& fixed, double tol) {
  std::vector<Id> ids = nodes_on_segment(model, a, b, tol);
  for (Id id : ids) model.add_support(id, fixed);
  return ids;
}

}  // namespace keta

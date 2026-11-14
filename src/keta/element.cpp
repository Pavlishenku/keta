#include "keta/element.h"

#include <cmath>

namespace keta {

namespace {

Mat3 skew(const Vec3& v) {
  Mat3 s;
  s << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
  return s;
}

// 5-point Gauss-Legendre on [0,1] — exact through degree 9, covers every
// polynomial shape-function x load-shape product used here.
constexpr int kGauss = 5;
constexpr double kGx[kGauss] = {0.046910077030668, 0.230765344947158, 0.5,
                                0.769234655052842, 0.953089922969332};
constexpr double kGw[kGauss] = {0.118463442528095, 0.239314335249683, 0.284444444444444,
                                0.239314335249683, 0.118463442528095};

}  // namespace

ShearPhi shear_phi(const Material& m, const SectionProps& s, double L) {
  ShearPhi p;
  double G = m.G();
  if (s.Avy > 0) p.phi_y = 12.0 * m.E * s.Iz / (G * s.Avy * L * L);
  if (s.Avz > 0) p.phi_z = 12.0 * m.E * s.Iy / (G * s.Avz * L * L);
  return p;
}

Mat3 bar_rotation(const Vec3& xi, const Vec3& xj, double roll) {
  Vec3 dx = xj - xi;
  double L = dx.norm();
  require(L > 1e-12, "element: zero-length bar");
  Vec3 ex = dx / L;
  Vec3 ey0, ez0;
  if (std::fabs(ex.z()) > 0.99999) {  // vertical bar: local y = global Y
    ey0 = Vec3(0, 1, 0);
    ez0 = ex.cross(ey0).normalized();
    ey0 = ez0.cross(ex);
  } else {
    ey0 = Vec3(0, 0, 1).cross(ex).normalized();  // horizontal, perpendicular to bar
    ez0 = ex.cross(ey0);                          // tilted towards +Z
  }
  double c = std::cos(roll), s = std::sin(roll);
  Vec3 ey = c * ey0 + s * ez0;
  Vec3 ez = -s * ey0 + c * ez0;
  Mat3 R;
  R.row(0) = ex;
  R.row(1) = ey;
  R.row(2) = ez;
  return R;
}

BarElement::BarElement(const Model& model, const Bar& bar)
    : model_(&model),
      bar_(&bar),
      mat_(&model.material(bar.material)),
      sec_(&model.section(bar.section).props) {
  const Vec3& xi = model.node(bar.ni).x;
  const Vec3& xj = model.node(bar.nj).x;
  // Element axis runs between the OFFSET ends when eccentricity is present;
  // the rigid arms connect the node line to the element line.
  Vec3 pi = xi + bar.ecc[0];
  Vec3 pj = xj + bar.ecc[1];
  L_ = (pj - pi).norm();
  require(L_ > 1e-12, "bar " + std::to_string(bar.id) + ": zero length after eccentricity");
  R_ = bar_rotation(pi, pj, bar.roll);
  build_local_stiffness();
  build_condensation();
  build_transform();
  ke_nodal_ = T_.transpose() * (k_local_cond_ + k_soil_local_) * T_;
}

void BarElement::build_local_stiffness() {
  Mat12& k = k_local_raw_;
  k.setZero();
  double L = L_, E = mat_->E, G = mat_->G();
  double A = sec_->A;

  if (bar_->kind == BarKind::Truss || bar_->kind == BarKind::Cable) {
    double eal = E * A / L;
    k(0, 0) = k(6, 6) = eal;
    k(0, 6) = k(6, 0) = -eal;
    return;
  }

  double Iy = sec_->Iy, Iz = sec_->Iz, J = sec_->J;
  ShearPhi phi = shear_phi(*mat_, *sec_, L);

  double eal = E * A / L, gjl = G * J / L;
  k(0, 0) = k(6, 6) = eal;
  k(0, 6) = k(6, 0) = -eal;
  k(3, 3) = k(9, 9) = gjl;
  k(3, 9) = k(9, 3) = -gjl;

  auto bend = [&](double I, double ph, int iv, int ir, int jv, int jr, double sgn) {
    // Plane bending block with shear parameter ph. sgn = +1 for the x-y plane
    // (v, theta_z), -1 for the x-z plane (w, theta_y) where theta_y = -dw/dx.
    double d = 1.0 + ph;
    double a = 12 * E * I / (L * L * L * d);
    double b = 6 * E * I / (L * L * d) * sgn;
    double c = (4 + ph) * E * I / (L * d);
    double e = (2 - ph) * E * I / (L * d);
    k(iv, iv) += a;
    k(iv, ir) += b;
    k(ir, iv) += b;
    k(ir, ir) += c;
    k(jv, jv) += a;
    k(jv, jr) -= b;
    k(jr, jv) -= b;
    k(jr, jr) += c;
    k(iv, jv) -= a;
    k(jv, iv) -= a;
    k(iv, jr) += b;
    k(jr, iv) += b;
    k(ir, jv) -= b;
    k(jv, ir) -= b;
    k(ir, jr) += e;
    k(jr, ir) += e;
  };
  bend(Iz, phi.phi_y, 1, 5, 7, 11, +1.0);  // x-y plane: uy, rz
  bend(Iy, phi.phi_z, 2, 4, 8, 10, -1.0);  // x-z plane: uz, ry

  // Winkler bed lives in a SEPARATE matrix: it is added to the nodal stiffness
  // after release condensation (a hinge frees the member moment but the soil
  // still reacts to the physical displacement field) and it never enters
  // member-force recovery — soil pressure is recovered separately as p = k u.
  const SoilProps& soil = bar_->soil;
  k_soil_local_.setZero();
  if (soil.any()) {
    Mat12& ks = k_soil_local_;
    auto add_axial_like = [&](double kf, int i0, int j0) {
      double m = kf * L;
      ks(i0, i0) += m / 3;
      ks(j0, j0) += m / 3;
      ks(i0, j0) += m / 6;
      ks(j0, i0) += m / 6;
    };
    auto add_hermite = [&](double kf, int iv, int ir, int jv, int jr, double sgn) {
      double m = kf * L;
      double l = L * sgn;
      ks(iv, iv) += m * 13 / 35;
      ks(iv, ir) += m * 11 * l / 210;
      ks(ir, iv) += m * 11 * l / 210;
      ks(ir, ir) += m * l * l / 105;
      ks(jv, jv) += m * 13 / 35;
      ks(jv, jr) -= m * 11 * l / 210;
      ks(jr, jv) -= m * 11 * l / 210;
      ks(jr, jr) += m * l * l / 105;
      ks(iv, jv) += m * 9 / 70;
      ks(jv, iv) += m * 9 / 70;
      ks(iv, jr) -= m * 13 * l / 420;
      ks(jr, iv) -= m * 13 * l / 420;
      ks(ir, jv) += m * 13 * l / 420;
      ks(jv, ir) += m * 13 * l / 420;
      ks(ir, jr) -= m * l * l / 140;
      ks(jr, ir) -= m * l * l / 140;
    };
    if (soil.kfx > 0) add_axial_like(soil.kfx, 0, 6);
    if (soil.kmx > 0) add_axial_like(soil.kmx, 3, 9);
    if (soil.kfy > 0) add_hermite(soil.kfy, 1, 5, 7, 11, +1.0);
    if (soil.kfz > 0) add_hermite(soil.kfz, 2, 4, 8, 10, -1.0);
  }
}

void BarElement::build_condensation() {
  rel_.clear();
  keep_.clear();
  for (int i = 0; i < 12; ++i)
    (bar_->releases[static_cast<std::size_t>(i)] ? rel_ : keep_).push_back(i);
  has_rel_ = !rel_.empty();
  if (!has_rel_) {
    k_local_cond_ = k_local_raw_;
    return;
  }
  require(rel_.size() < 12, "bar " + std::to_string(bar_->id) + ": all DOFs released");
  int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
  MatX krr(nr, nr), krk(nr, nk), kkk(nk, nk);
  for (int a = 0; a < nr; ++a) {
    for (int b = 0; b < nr; ++b) krr(a, b) = k_local_raw_(rel_[a], rel_[b]);
    for (int b = 0; b < nk; ++b) krk(a, b) = k_local_raw_(rel_[a], keep_[b]);
  }
  for (int a = 0; a < nk; ++a)
    for (int b = 0; b < nk; ++b) kkk(a, b) = k_local_raw_(keep_[a], keep_[b]);
  Eigen::FullPivLU<MatX> lu(krr);
  require(lu.isInvertible(),
          "bar " + std::to_string(bar_->id) + ": release pattern makes the bar a mechanism");
  krr_inv_ = lu.inverse();
  krr_inv_krk_ = krr_inv_ * krk;
  MatX kc = kkk - krk.transpose() * krr_inv_krk_;
  k_local_cond_.setZero();
  for (int a = 0; a < nk; ++a)
    for (int b = 0; b < nk; ++b) k_local_cond_(keep_[a], keep_[b]) = kc(a, b);
}

void BarElement::build_transform() {
  Mat12 rot = Mat12::Zero();
  for (int blk = 0; blk < 4; ++blk) rot.block<3, 3>(3 * blk, 3 * blk) = R_;
  Mat12 ecc = Mat12::Identity();
  if (bar_->has_ecc()) {
    // u_end = u_node + theta x e  ->  [I  -S(e); 0  I] per end (global frame).
    ecc.block<3, 3>(0, 3) = -skew(bar_->ecc[0]);
    ecc.block<3, 3>(6, 9) = -skew(bar_->ecc[1]);
  }
  T_ = rot * ecc;
}

Mat12 BarElement::me(bool lumped) const {
  Mat12 m = Mat12::Zero();
  double rho = mat_->rho, A = sec_->A, L = L_;
  double Ip = sec_->Ip > 0 ? sec_->Ip : sec_->Iy + sec_->Iz;  // torsional inertia
  double mass = rho * A * L;
  if (lumped) {
    for (int i : {0, 1, 2, 6, 7, 8}) m(i, i) = mass / 2;
    m(3, 3) = m(9, 9) = rho * Ip * L / 2;
  } else if (bar_->kind != BarKind::Frame) {
    // truss / cable: the consistent mass of the LINEAR interpolation on the
    // three translations, rho A L / 6 [2 1; 1 2] (the bar's own shape: no
    // bending field), the torsional rotary inertia lumped as before
    for (int c = 0; c < 3; ++c) {
      m(c, c) = m(6 + c, 6 + c) = mass / 3;
      m(c, 6 + c) = m(6 + c, c) = mass / 6;
    }
    m(3, 3) = m(9, 9) = rho * Ip * L / 2;
  } else {
    auto axial_like = [&](double mm, int i0, int j0) {
      m(i0, i0) += mm / 3;
      m(j0, j0) += mm / 3;
      m(i0, j0) += mm / 6;
      m(j0, i0) += mm / 6;
    };
    auto hermite = [&](double mm, int iv, int ir, int jv, int jr, double sgn) {
      double l = L * sgn;
      m(iv, iv) += mm * 13 / 35;
      m(iv, ir) += mm * 11 * l / 210;
      m(ir, iv) += mm * 11 * l / 210;
      m(ir, ir) += mm * l * l / 105;
      m(jv, jv) += mm * 13 / 35;
      m(jv, jr) -= mm * 11 * l / 210;
      m(jr, jv) -= mm * 11 * l / 210;
      m(jr, jr) += mm * l * l / 105;
      m(iv, jv) += mm * 9 / 70;
      m(jv, iv) += mm * 9 / 70;
      m(iv, jr) -= mm * 13 * l / 420;
      m(jr, iv) -= mm * 13 * l / 420;
      m(ir, jv) += mm * 13 * l / 420;
      m(jv, ir) += mm * 13 * l / 420;
      m(ir, jr) -= mm * l * l / 140;
      m(jr, ir) -= mm * l * l / 140;
    };
    axial_like(mass, 0, 6);
    axial_like(rho * Ip * L, 3, 9);
    hermite(mass, 1, 5, 7, 11, +1.0);
    hermite(mass, 2, 4, 8, 10, -1.0);
    // Linear-theta rotary inertia of the cross-section about y and z.
    axial_like(rho * sec_->Iz * L, 5, 11);
    axial_like(rho * sec_->Iy * L, 4, 10);
  }
  // Nodal frame: T^t m T (rigid arms carry inertia transfer as well).
  return T_.transpose() * m * T_;
}

Mat12 BarElement::kg(double N) const {
  Mat12 g = Mat12::Zero();
  double L = L_;
  if (bar_->kind != BarKind::Frame) {
    // Truss/cable string geometric stiffness on the transverse translations.
    double c = N / L;
    for (int d : {1, 2}) {
      g(d, d) += c;
      g(d + 6, d + 6) += c;
      g(d, d + 6) -= c;
      g(d + 6, d) -= c;
    }
  } else {
    auto plane = [&](int iv, int ir, int jv, int jr, double sgn) {
      double c = N / L;
      double l = L * sgn;
      g(iv, iv) += c * 6 / 5;
      g(iv, ir) += c * l / 10;
      g(ir, iv) += c * l / 10;
      g(ir, ir) += c * 2 * l * l / 15;
      g(jv, jv) += c * 6 / 5;
      g(jv, jr) -= c * l / 10;
      g(jr, jv) -= c * l / 10;
      g(jr, jr) += c * 2 * l * l / 15;
      g(iv, jv) -= c * 6 / 5;
      g(jv, iv) -= c * 6 / 5;
      g(iv, jr) += c * l / 10;
      g(jr, iv) += c * l / 10;
      g(ir, jv) -= c * l / 10;
      g(jv, ir) -= c * l / 10;
      g(ir, jr) -= c * l * l / 30;
      g(jr, ir) -= c * l * l / 30;
    };
    plane(1, 5, 7, 11, +1.0);
    plane(2, 4, 8, 10, -1.0);
    // Wagner/torsional term for spatial members: N * Ip/A coupling on twist.
    if (sec_->A > 0) {
      double c = N * (sec_->Iy + sec_->Iz) / (sec_->A * L);
      g(3, 3) += c;
      g(9, 9) += c;
      g(3, 9) -= c;
      g(9, 3) -= c;
    }
  }
  if (has_rel_) {
    // Condense the geometric stiffness with the same release pattern (using the
    // ELASTIC condensation operator — standard practical treatment).
    int nk = static_cast<int>(keep_.size());
    int nr = static_cast<int>(rel_.size());
    MatX grr(nr, nr), grk(nr, nk), gkk(nk, nk);
    for (int a = 0; a < nr; ++a) {
      for (int b = 0; b < nr; ++b) grr(a, b) = g(rel_[a], rel_[b]);
      for (int b = 0; b < nk; ++b) grk(a, b) = g(rel_[a], keep_[b]);
    }
    for (int a = 0; a < nk; ++a)
      for (int b = 0; b < nk; ++b) gkk(a, b) = g(keep_[a], keep_[b]);
    MatX gc = gkk - grk.transpose() * krr_inv_krk_ - krr_inv_krk_.transpose() * grk +
              krr_inv_krk_.transpose() * grr * krr_inv_krk_;
    g.setZero();
    for (int a = 0; a < nk; ++a)
      for (int b = 0; b < nk; ++b) g(keep_[a], keep_[b]) = gc(a, b);
  }
  return T_.transpose() * g * T_;
}

// ---------------- span loads ----------------

namespace {

// Timoshenko interpolation (IIE) for one bending plane on [0,1].
// v(x) = H1 v_i + H2 t_i + H3 v_j + H4 t_j with t = rotation DOF (already
// sign-adjusted by the caller for the x-z plane).
struct BendShapes {
  double h1, h2, h3, h4;    // deflection interpolation
  double d1, d2, d3, d4;    // SECTION-ROTATION interpolation (IIE): the work
                            // of a couple is done through the rotation, not
                            // dv/dx (they differ by the shear strain)
};

BendShapes bend_shapes(double xi, double L, double phi) {
  BendShapes s;
  double d = 1.0 + phi;
  double x2 = xi * xi, x3 = x2 * xi;
  s.h1 = (1 - 3 * x2 + 2 * x3 + phi * (1 - xi)) / d;
  s.h2 = L * (xi - 2 * x2 + x3 + phi / 2 * (xi - x2)) / d;
  s.h3 = (3 * x2 - 2 * x3 + phi * xi) / d;
  s.h4 = L * (-x2 + x3 + phi / 2 * (-xi + x2)) / d;
  s.d1 = 6 * (x2 - xi) / (d * L);
  s.d2 = (1 - 4 * xi + 3 * x2 + phi * (1 - xi)) / d;
  s.d3 = -s.d1;
  s.d4 = (3 * x2 - 2 * xi + phi * xi) / d;
  return s;
}

}  // namespace

Vec12 BarElement::f0_local(const Load& load, const Material& mat) const {
  Vec12 f = Vec12::Zero();
  double L = L_;
  ShearPhi phi = bar_->kind == BarKind::Frame ? shear_phi(*mat_, *sec_, L) : ShearPhi{};

  auto to_local_force = [&](const Vec3& Fg, LoadFrame frame) -> Vec3 {
    if (frame == LoadFrame::Local) return Fg;
    return R_ * Fg;
  };
  auto projected_scale = [&](LoadFrame frame) -> Vec3 {
    if (frame != LoadFrame::GlobalProjected) return Vec3::Ones();
    // Per-component: intensity given per unit length of the bar's projection
    // perpendicular to that global axis (ST1 PRO semantics).
    Vec3 ex = R_.row(0);
    return Vec3(std::sqrt(std::max(0.0, 1 - ex.x() * ex.x())),
                std::sqrt(std::max(0.0, 1 - ex.y() * ex.y())),
                std::sqrt(std::max(0.0, 1 - ex.z() * ex.z())));
  };
  auto range = [&](double& a, double& b) {
    a = load.a;
    b = load.b < 0 ? (load.relative ? 1.0 : L) : load.b;
    if (load.relative) {
      a *= L;
      b *= L;
    }
    require(a >= -1e-9 && b <= L * (1 + 1e-9) && b >= a,
            "load: abscissa range outside bar length");
    a = std::clamp(a, 0.0, L);
    b = std::clamp(b, 0.0, L);
  };

  // Accumulate a transverse distributed load q(x) (local) over [a,b] by Gauss
  // quadrature of the shape-function work integrals: linear qa -> qb plus the
  // parabolic excess qc at mid-range (5 points: exact to degree 9, the cubic
  // shapes times a quadratic load are degree 5).
  auto add_distributed = [&](double a, double b, const Vec3& qa, const Vec3& qb,
                             const Vec3& qc = Vec3::Zero()) {
    if (b - a < 1e-14) return;
    double len = b - a;
    for (int gp = 0; gp < kGauss; ++gp) {
      double x = a + kGx[gp] * len;
      double w = kGw[gp] * len;
      double t = (x - a) / len;
      Vec3 q = qa + (qb - qa) * t + 4.0 * t * (1.0 - t) * qc;
      double xi = x / L;
      // Axial (linear shapes)
      f(0) += w * q.x() * (1 - xi);
      f(6) += w * q.x() * xi;
      if (bar_->kind == BarKind::Frame) {
        BendShapes sy = bend_shapes(xi, L, phi.phi_y);
        f(1) += w * q.y() * sy.h1;
        f(5) += w * q.y() * sy.h2;
        f(7) += w * q.y() * sy.h3;
        f(11) += w * q.y() * sy.h4;
        BendShapes sz = bend_shapes(xi, L, phi.phi_z);
        f(2) += w * q.z() * sz.h1;
        f(4) -= w * q.z() * sz.h2;  // theta_y = -dw/dx
        f(8) += w * q.z() * sz.h3;
        f(10) -= w * q.z() * sz.h4;
      } else {
        f(1) += w * q.y() * (1 - xi);
        f(7) += w * q.y() * xi;
        f(2) += w * q.z() * (1 - xi);
        f(8) += w * q.z() * xi;
      }
    }
  };

  // A distributed moment field m(x) (local, linear ma -> mb over [a,b]): the
  // work of a point moment integrated by Gauss quadrature — torsion through
  // the linear twist shapes, bending through the slopes (theta_y = -dw/dx).
  auto add_distributed_moment = [&](double a, double b, const Vec3& ma, const Vec3& mb) {
    if (b - a < 1e-14) return;
    const double len = b - a;
    for (int gp = 0; gp < kGauss; ++gp) {
      const double x = a + kGx[gp] * len, w = kGw[gp] * len;
      const Vec3 M = ma + (mb - ma) * ((x - a) / len);
      const double xi = x / L;
      f(3) += w * M.x() * (1 - xi);
      f(9) += w * M.x() * xi;
      if (bar_->kind == BarKind::Frame) {
        BendShapes sy = bend_shapes(xi, L, phi.phi_y);
        f(1) += w * M.z() * sy.d1;
        f(5) += w * M.z() * sy.d2;
        f(7) += w * M.z() * sy.d3;
        f(11) += w * M.z() * sy.d4;
        BendShapes sz = bend_shapes(xi, L, phi.phi_z);
        f(2) -= w * M.y() * sz.d1;
        f(4) += w * M.y() * sz.d2;
        f(8) -= w * M.y() * sz.d3;
        f(10) += w * M.y() * sz.d4;
      }
    }
  };

  switch (load.kind) {
    case Load::Kind::BarUdl:
    case Load::Kind::BarLinear: {
      double a, b;
      range(a, b);
      Vec3 s = projected_scale(load.frame);
      Vec3 qa_g = load.F.cwiseProduct(s);
      Vec3 qb_g = (load.kind == Load::Kind::BarUdl ? load.F : load.F2).cwiseProduct(s);
      Vec3 qa = to_local_force(qa_g, load.frame);
      Vec3 qb = to_local_force(qb_g, load.frame);
      add_distributed(a, b, qa, qb);
      const Vec3 mb_g = load.kind == Load::Kind::BarUdl ? load.M : load.M2;
      if (load.M.squaredNorm() > 0 || mb_g.squaredNorm() > 0)
        add_distributed_moment(a, b, to_local_force(load.M, load.frame),
                               to_local_force(mb_g, load.frame));
      break;
    }
    case Load::Kind::BarParabolic: {
      double a, b;
      range(a, b);
      const Vec3 s = projected_scale(load.frame);
      const Vec3 qa = to_local_force(load.F.cwiseProduct(s), load.frame);
      const Vec3 qb = to_local_force(load.F2.cwiseProduct(s), load.frame);
      const Vec3 qmid = to_local_force(load.Fmid.cwiseProduct(s), load.frame);
      add_distributed(a, b, qa, qb, qmid - 0.5 * (qa + qb));
      break;
    }
    case Load::Kind::BarPoint: {
      double a = load.relative ? load.a * L : load.a;
      require(a >= -1e-9 && a <= L * (1 + 1e-9), "point load outside bar");
      a = std::clamp(a, 0.0, L);
      double xi = a / L;
      Vec3 F = to_local_force(load.F, load.frame);
      Vec3 M = to_local_force(load.M, load.frame);
      f(0) += F.x() * (1 - xi);
      f(6) += F.x() * xi;
      f(3) += M.x() * (1 - xi);
      f(9) += M.x() * xi;
      if (bar_->kind == BarKind::Frame) {
        BendShapes sy = bend_shapes(xi, L, phi.phi_y);
        f(1) += F.y() * sy.h1;
        f(5) += F.y() * sy.h2;
        f(7) += F.y() * sy.h3;
        f(11) += F.y() * sy.h4;
        BendShapes sz = bend_shapes(xi, L, phi.phi_z);
        f(2) += F.z() * sz.h1;
        f(4) -= F.z() * sz.h2;
        f(8) += F.z() * sz.h3;
        f(10) -= F.z() * sz.h4;
        // Point moment about local z does work through the section rotation
        // theta_z; about local y through theta_y (sign as -dw/dx).
        f(1) += M.z() * sy.d1;
        f(5) += M.z() * sy.d2;
        f(7) += M.z() * sy.d3;
        f(11) += M.z() * sy.d4;
        f(2) -= M.y() * sz.d1;
        f(4) += M.y() * sz.d2;
        f(8) -= M.y() * sz.d3;
        f(10) += M.y() * sz.d4;
      } else {
        f(1) += F.y() * (1 - xi);
        f(7) += F.y() * xi;
        f(2) += F.z() * (1 - xi);
        f(8) += F.z() * xi;
      }
      break;
    }
    case Load::Kind::BarThermal: {
      double ea = mat.E * sec_->A * mat.alpha * load.dt;
      f(0) -= ea;
      f(6) += ea;
      if (bar_->kind == BarKind::Frame) {
        // Gradient gy (dT/dy) bends in the x-y plane about z; gz about y.
        // f_eq = int B^t E eps_th dV = -E I alpha g int N'' dx = +M at end i, -M at end j.
        // The hotter face lengthens and the bar bows TOWARD it: kappa_z =
        // +alpha gy, but kappa_y = -alpha gz, because keta's rotation about y
        // is theta_y = -dw/dx (a top-hot simply supported beam must rise;
        // the + sign made it sag — caught against ST1's GZ references).
        double mz = mat.E * sec_->Iz * mat.alpha * load.gy;
        double my = -mat.E * sec_->Iy * mat.alpha * load.gz;
        f(5) += mz;
        f(11) -= mz;
        f(4) += my;
        f(10) -= my;
      }
      break;
    }
    case Load::Kind::Prestrain: {
      double ea = mat.E * sec_->A * load.strain;
      f(0) -= ea;
      f(6) += ea;
      if (bar_->kind == BarKind::Frame) {
        // imposed twist rate: the clamped shaft holds T = -GJ kappa_x
        const double gj = mat.G() * sec_->J * load.kappa_x;
        f(3) -= gj;
        f(9) += gj;
        double mz = mat.E * sec_->Iz * load.kappa_z;
        double my = mat.E * sec_->Iy * load.kappa_y;
        f(5) += mz;
        f(11) -= mz;
        f(4) += my;
        f(10) -= my;
      }
      break;
    }
    case Load::Kind::PrestrainPoint: {
      // Concentrated dislocation at x0 (ST1 DEFOR CON): axial gap `strain`,
      // transverse jumps dvy / dvz [m], twist jump dthx and the kinks
      // kappa_y / kappa_z [rad] (keta's curvature sign: the rotation jumps by
      // -kappa). Free cantilever from end i: the stress-free tip displacement
      // it produces at j, and f0 = k_member [0; tip] — the forces holding the
      // element's ends where the dislocation leaves them free (exact for the
      // exact member stiffness, Timoshenko included). A truss bar keeps only
      // the axial gap (its member stiffness has no bending).
      const Load d = load.frame == LoadFrame::Global ? dislocation_local(load, R_) : load;
      double a = d.relative ? d.a * L : d.a;
      require(a >= -1e-9 && a <= L * (1 + 1e-9), "prestrain dislocation outside bar");
      a = std::clamp(a, 0.0, L);
      const double thz = -d.kappa_z, thy = -d.kappa_y;
      Vec6 tip;
      tip << d.strain,                    // u
          d.dvy + thz * (L - a),          // v: v' = theta_z
          d.dvz - thy * (L - a),          // w: w' = -theta_y
          d.dthx, thy, thz;
      f += k_local_raw_.block<12, 6>(0, 6) * tip;
      break;
    }
    case Load::Kind::PrestrainLinear: {
      // Linearly varying imposed strain/curvatures over the full bar
      // (ST1 DEFOR LIN): exact work integrals of the linear fields against the
      // exact interpolation, in closed form. With k(xi) = ka + dk*xi:
      //   int B_b k dx over the four rotation shapes gives
      //   [dk/(L(1+phi)), -ka - dk*phi/(2(1+phi)), -dk/(L(1+phi)),
      //    ka + dk*(1 + phi/2)/(1+phi)]  (reduces to the UNI couples at dk=0).
      double e_mean = 0.5 * (load.strain + load.strain2);
      double ea = mat.E * sec_->A * e_mean;
      f(0) -= ea;
      f(6) += ea;
      if (bar_->kind == BarKind::Frame) {
        // a linear twist rate: the clamped shaft holds the constant torque
        // -GJ times its mean (int theta' = 0 over the bar)
        const double gj = mat.G() * sec_->J * 0.5 * (load.kappa_x + load.kappa_x2);
        f(3) -= gj;
        f(9) += gj;
        auto lower_plane = [&](double EI, double ka, double kb, double phi_p, int iv, int it,
                               int jv, int jt, double sv) {
          double dk = kb - ka;
          double d = 1.0 + phi_p;
          double i1 = dk / (L * d);                 // int g1 k dx
          double i2 = -ka - dk * phi_p / (2 * d);   // int g2 k dx
          double i4 = ka + dk * (1 + phi_p / 2) / d;  // int g4 k dx
          f(iv) -= sv * EI * i1;
          f(it) -= EI * i2;
          f(jv) -= sv * EI * (-i1);
          f(jt) -= EI * i4;
        };
        lower_plane(mat.E * sec_->Iz, load.kappa_z, load.kappa_z2, phi.phi_y, 1, 5, 7, 11, 1.0);
        lower_plane(mat.E * sec_->Iy, load.kappa_y, load.kappa_y2, phi.phi_z, 2, 4, 8, 10, -1.0);
      }
      break;
    }
    case Load::Kind::PrestrainParabolic: {
      // Parabolic imposed fields over the full bar: k(xi) = ka + dk xi + c2 xi^2
      // with dk = kb - ka + 4c, c2 = -4c, c = k_mid - (ka + kb)/2. The
      // curvature shapes g of the exact (Timoshenko) element are linear in
      // xi, so int g xi^2 dx = int g xi dx - (1/6) int g dx: the moments of
      // PrestrainLinear give the quadratic term in closed form,
      //   g1: (0, 1/(L d), 1/(L d)), g2: (-1, -phi/(2d), 1/6 - phi/(2d)),
      //   g4: (1, (1 + phi/2)/d, (1 + phi/2)/d - 1/6)   (d = 1 + phi),
      // the (int g, int g xi, int g xi^2) of each shape; g3 = -g1.
      auto mean = [](double a, double m, double b) {  // (1/L) int f dx = (a + 4m + b)/6
        return (a + 4.0 * m + b) / 6.0;
      };
      const double ea = mat.E * sec_->A * mean(load.strain, load.strain_m, load.strain2);
      f(0) -= ea;
      f(6) += ea;
      if (bar_->kind == BarKind::Frame) {
        // the clamped shaft holds the constant torque -GJ times the mean rate
        const double gj = mat.G() * sec_->J * mean(load.kappa_x, load.kappa_x_m, load.kappa_x2);
        f(3) -= gj;
        f(9) += gj;
        auto lower_plane = [&](double EI, double ka, double km, double kb, double phi_p, int iv, int it,
                               int jv, int jt, double sv) {
          const double c = km - 0.5 * (ka + kb);
          const double dk = kb - ka + 4.0 * c, c2 = -4.0 * c;
          const double d = 1.0 + phi_p;
          const double i1 = (dk + c2) / (L * d);
          const double i2 = -ka - dk * phi_p / (2 * d) + c2 * (1.0 / 6.0 - phi_p / (2 * d));
          const double i4 = ka + dk * (1 + phi_p / 2) / d + c2 * ((1 + phi_p / 2) / d - 1.0 / 6.0);
          f(iv) -= sv * EI * i1;
          f(it) -= EI * i2;
          f(jv) -= sv * EI * (-i1);
          f(jt) -= EI * i4;
        };
        lower_plane(mat.E * sec_->Iz, load.kappa_z, load.kappa_z_m, load.kappa_z2, phi.phi_y, 1, 5, 7, 11, 1.0);
        lower_plane(mat.E * sec_->Iy, load.kappa_y, load.kappa_y_m, load.kappa_y2, phi.phi_z, 2, 4, 8, 10, -1.0);
      }
      break;
    }
    case Load::Kind::SelfWeight: {
      Vec3 qg = mat.rho * sec_->A * load.g;  // N/m, global
      Vec3 ql = R_ * qg;
      add_distributed(0, L, ql, ql);
      break;
    }
    default: fail("f0_local: unsupported load kind for a bar element");
  }
  return f;
}

Vec12 BarElement::condense_f0(const Vec12& f0) const {
  if (!has_rel_) return f0;
  int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
  VecX fr(nr);
  for (int a = 0; a < nr; ++a) fr(a) = f0(rel_[a]);
  VecX corr = krr_inv_krk_.transpose() * fr;  // k_kr k_rr^-1 f_r
  Vec12 out = Vec12::Zero();
  for (int a = 0; a < nk; ++a) out(keep_[a]) = f0(keep_[a]) - corr(a);
  return out;
}

Vec12 BarElement::equivalent_nodal(const Load& load, const Material& mat) const {
  Vec12 f0 = f0_local(load, mat);
  Vec12 fc = condense_f0(f0);
  return T_.transpose() * fc;
}

Vec12 BarElement::local_displacements(const Vec12& u_nodal, const Vec12& f0_local_total) const {
  Vec12 u = T_ * u_nodal;
  if (has_rel_) {
    // Released DOFs carry no member force: u_r = k_rr^-1 (f0_r - k_rk u_k).
    int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
    VecX rhs(nr);
    for (int a = 0; a < nr; ++a) {
      double v = f0_local_total(rel_[a]);
      for (int b = 0; b < nk; ++b) v -= k_local_raw_(rel_[a], keep_[b]) * u(keep_[b]);
      rhs(a) = v;
    }
    VecX ur = krr_inv_ * rhs;
    for (int a = 0; a < nr; ++a) u(rel_[a]) = ur(a);
  }
  return u;
}

Vec12 BarElement::end_forces_local(const Vec12& u_local, const Vec12& f0_local_total) const {
  return k_local_raw_ * u_local - f0_local_total;
}

}  // namespace keta

namespace keta {

StabilityFunctions stability_functions(double N, double EI, double L) {
  StabilityFunctions f;
  require(EI > 0 && L > 0, "stability functions: bad EI or L");
  double u2 = -N * L * L / EI;  // >0 compression
  f.u2 = u2;
  double u = std::sqrt(std::fabs(u2));
  if (u < 1.0) {
    // The closed forms cancel catastrophically for small u: D = 2 - 2 cos u
    // - u sin u ~ u^4/12 (s = 4.000985 instead of 3.9999998 at u = 1.2e-3,
    // so the P-Delta fixed point of short bars never settled). Their Taylor
    // series in the SIGNED q = u2 (the hyperbolic forms are the same
    // functions continued to q < 0), coefficients from the exact rational
    // series division, sixteen terms: the truncation is below 1e-20 for
    // |q| < 1 (the coefficients fall by ~1/20 a term: radius (4.4934)^2,
    // set by sin u = u cos u).
    static constexpr double ks[] = {4, -0.13333333333333333, -0.001746031746031746,
        -3.7037037037037037e-05, -8.7439016010444587e-07, -2.1461489715457969e-08,
        -5.3563706247001782e-10, -1.3471819416419479e-11, -3.4007314847583161e-13,
        -8.5997439884052179e-15, -2.1765627192905307e-16, -5.5111003240982866e-18,
        -1.3957061776974721e-19, -3.5350286089192456e-21, -8.9539155913818665e-23,
        -2.2680017527226832e-24};
    static constexpr double kc[] = {0.5, 0.025000000000000001, 0.0013095238095238095,
        6.6468253968253967e-05, 3.3227685013399297e-06, 1.6511872166634071e-07,
        8.1874668751085992e-09, 4.0566847221800868e-10, 2.0094560802570245e-11,
        9.9528316883373345e-13, 4.929484385150959e-14, 2.4414723873963281e-15,
        1.2092068792119806e-16, 5.9889251306090383e-18, 2.9661764554823286e-19,
        1.4690785695343879e-20};
    double ps = 0, pc = 0;
    for (int k = 15; k >= 0; --k) {
      ps = ps * u2 + ks[k];
      pc = pc * u2 + kc[k];
    }
    f.s = ps;
    f.c = pc;
    return f;
  }
  if (u2 > 0) {
    // compression (trigonometric); reject at/inside the first singularity
    double D = 2.0 - 2.0 * std::cos(u) - u * std::sin(u);
    require(std::fabs(D) > 1e-14 && u < 2.0 * kPi,
            "stability functions: axial force at or beyond the member's own buckling range");
    double denom = std::sin(u) - u * std::cos(u);
    f.s = u * denom / D;
    f.c = (u - std::sin(u)) / denom;
  } else {
    // tension (hyperbolic) — always regular
    double D = 2.0 - 2.0 * std::cosh(u) + u * std::sinh(u);
    double denom = u * std::cosh(u) - std::sinh(u);
    f.s = u * denom / D;
    f.c = (std::sinh(u) - u) / denom;
  }
  return f;
}

Mat12 BarElement::condense_matrix(const Mat12& k) const {
  if (!has_rel_) return k;
  int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
  MatX krr(nr, nr), krk(nr, nk), kkk(nk, nk);
  for (int a = 0; a < nr; ++a) {
    for (int b = 0; b < nr; ++b) krr(a, b) = k(rel_[a], rel_[b]);
    for (int b = 0; b < nk; ++b) krk(a, b) = k(rel_[a], keep_[b]);
  }
  for (int a = 0; a < nk; ++a)
    for (int b = 0; b < nk; ++b) kkk(a, b) = k(keep_[a], keep_[b]);
  Eigen::FullPivLU<MatX> lu(krr);
  require(lu.isInvertible(), "condense: release pattern makes the bar a mechanism");
  MatX kc = kkk - krk.transpose() * lu.solve(krk);
  Mat12 out = Mat12::Zero();
  for (int a = 0; a < nk; ++a)
    for (int b = 0; b < nk; ++b) out(keep_[a], keep_[b]) = kc(a, b);
  return out;
}

Mat12 BarElement::k_exact_local(double N) const {
  if (bar_->kind != BarKind::Frame) return k_local_raw_;
  double L = L_, E = mat_->E, G = mat_->G();
  double A = sec_->A, Iy = sec_->Iy, Iz = sec_->Iz, J = sec_->J;

  Mat12 k = Mat12::Zero();
  // Uniform torsion with Wagner's term N r0^2 (r0^2 = (Iy + Iz) / A, the
  // same term as the consistent Kg): linear in N, so exact for a member
  // without warping; without it the torsional buckling N_T = GJ A / (Iy + Iz)
  // is missed.
  const double wagner = A > 0 ? N * (Iy + Iz) / A : 0.0;
  double eal = E * A / L, gjl = (G * J + wagner) / L;
  k(0, 0) = k(6, 6) = eal;
  k(0, 6) = k(6, 0) = -eal;
  k(3, 3) = k(9, 9) = gjl;
  k(3, 9) = k(9, 3) = -gjl;

  auto bend_exact = [&](double I, int iv, int ir, int jv, int jr, double sgn) {
    StabilityFunctions f = stability_functions(N, E * I, L);
    double EI_L = E * I / L;
    double kth_ii = EI_L * f.s;
    double kth_ij = EI_L * f.s * f.c;
    double kvt = E * I / (L * L) * f.s * (1.0 + f.c) * sgn;
    double kvv = E * I / (L * L * L) * (2.0 * f.s * (1.0 + f.c) - f.u2);
    k(iv, iv) += kvv;
    k(jv, jv) += kvv;
    k(iv, jv) -= kvv;
    k(jv, iv) -= kvv;
    k(ir, ir) += kth_ii;
    k(jr, jr) += kth_ii;
    k(ir, jr) += kth_ij;
    k(jr, ir) += kth_ij;
    k(iv, ir) += kvt;
    k(ir, iv) += kvt;
    k(iv, jr) += kvt;
    k(jr, iv) += kvt;
    k(jv, ir) -= kvt;
    k(ir, jv) -= kvt;
    k(jv, jr) -= kvt;
    k(jr, jv) -= kvt;
  };
  bend_exact(Iz, 1, 5, 7, 11, +1.0);
  bend_exact(Iy, 2, 4, 8, 10, -1.0);
  return k;
}

Mat12 BarElement::ke_exact(double N) const {
  if (bar_->kind != BarKind::Frame) {
    // truss/cable: linear axial + string geometric (already exact for a straight member)
    Mat12 k = k_local_cond_;
    return T_.transpose() * (condense_matrix(k) + k_soil_local_) * T_ +
           (N != 0 ? kg(N) : Mat12::Zero());
  }
  // The stability functions are exact down to N = 0 (Taylor series below
  // u = 1): no consistent-Kg patch near zero axial force.
  return T_.transpose() * (condense_matrix(k_exact_local(N)) + k_soil_local_) * T_;
}

Vec12 BarElement::condense_f0_with(const Mat12& k, const Vec12& f0) const {
  if (!has_rel_) return f0;
  int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
  MatX krr(nr, nr), krk(nr, nk);
  VecX fr(nr);
  for (int a = 0; a < nr; ++a) {
    fr(a) = f0(rel_[a]);
    for (int b = 0; b < nr; ++b) krr(a, b) = k(rel_[a], rel_[b]);
    for (int b = 0; b < nk; ++b) krk(a, b) = k(rel_[a], keep_[b]);
  }
  Eigen::FullPivLU<MatX> lu(krr);
  require(lu.isInvertible(), "condense: release pattern makes the bar a mechanism");
  VecX corr = krk.transpose() * lu.solve(fr);
  Vec12 out = Vec12::Zero();
  for (int a = 0; a < nk; ++a) out(keep_[a]) = f0(keep_[a]) - corr(a);
  return out;
}

Vec12 BarElement::local_displacements_with(const Mat12& k, const Vec12& u_nodal,
                                           const Vec12& f0) const {
  Vec12 u = T_ * u_nodal;
  if (!has_rel_) return u;
  int nr = static_cast<int>(rel_.size()), nk = static_cast<int>(keep_.size());
  MatX krr(nr, nr);
  VecX rhs(nr);
  for (int a = 0; a < nr; ++a) {
    double v = f0(rel_[a]);
    for (int b = 0; b < nk; ++b) v -= k(rel_[a], keep_[b]) * u(keep_[b]);
    rhs(a) = v;
    for (int b = 0; b < nr; ++b) krr(a, b) = k(rel_[a], rel_[b]);
  }
  Eigen::FullPivLU<MatX> lu(krr);
  require(lu.isInvertible(), "condense: release pattern makes the bar a mechanism");
  VecX ur = lu.solve(rhs);
  for (int a = 0; a < nr; ++a) u(rel_[a]) = ur(a);
  return u;
}

Load dislocation_local(const Load& l, const Mat3& R) {
  Load d = l;
  const Vec3 du = R * Vec3(l.strain, l.dvy, l.dvz);
  const Vec3 dth = R * Vec3(l.dthx, -l.kappa_y, -l.kappa_z);  // rotation jumps
  d.strain = du.x();
  d.dvy = du.y();
  d.dvz = du.z();
  d.dthx = dth.x();
  d.kappa_y = -dth.y();
  d.kappa_z = -dth.z();
  d.frame = LoadFrame::Local;
  return d;
}

}  // namespace keta

// Corotational kinematics shared by the geometrically nonlinear drivers
// (nonlinear.cpp: solve_corotational / arc_length_trace / ruin_load, and
// shell_plastic.cpp in large displacements): the exponential / logarithm
// maps written generically on the scalar type (double, ad::Grad, ad::Jet),
// the corotational bar (CorotElem, chord frame + mean director) and the
// element-independent corotational facet (CorotShell, EICR). Internal
// header: implementation detail of those drivers, not a public API.
#pragma once
#include <array>
#include <cmath>
#include <type_traits>

#include "keta/assembly.h"
#include "keta/jet.h"
#include "keta/nonlinear.h"

namespace keta::corot {

// ---------------- corotational element (energy formulation) ----------------
//
// The element energy W(q) = 1/2 d(q)^T Kbar d(q) is written ONCE, generically
// on the scalar type, and evaluated on double (value, recovery), on
// ad::Grad<12> (internal force = exact gradient, one pass) and on ad::Jet<12>
// (consistent tangent = exact Hessian, one pass). The derivatives are thus the
// analytic derivatives of the corotational energy to round-off — no finite
// difference step, no truncation error (the former 4th-order FD gradient
// carried ~1e-9 relative noise, which put a floor under every corrector
// tolerance). Smoothness at the rest state: the exponential and logarithm maps
// switch to their Taylor series in theta^2 / sin^2 below 1e-8, so no square
// root is ever taken at zero (exact to O(theta^6) — far below round-off — with
// exact first and second derivatives at theta = 0).

template <class T>
using V3t = std::array<T, 3>;
template <class T>
using M3t = std::array<T, 9>;  // row-major

template <class T>
T dot3(const V3t<T>& a, const V3t<T>& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
template <class T>
V3t<T> cross3(const V3t<T>& a, const V3t<T>& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Rotation vector -> matrix: R = I + A S + B (th th^T - t^2 I),
// A = sin t / t, B = (1 - cos t) / t^2.
template <class T>
M3t<T> exp_t(const V3t<T>& th) {
  using std::cos;
  using std::sin;
  using std::sqrt;
  const T t2 = dot3(th, th);
  T A = t2, B = t2;
  if (ad::val(t2) < 1e-8) {
    A = 1.0 - t2 / 6.0 + t2 * t2 / 120.0;
    B = 0.5 - t2 / 24.0 + t2 * t2 / 720.0;
  } else {
    const T t = sqrt(t2);
    A = sin(t) / t;
    B = (1.0 - cos(t)) / t2;
  }
  M3t<T> R;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) R[static_cast<std::size_t>(3 * i + j)] = B * th[static_cast<std::size_t>(i)] * th[static_cast<std::size_t>(j)];
  for (int i = 0; i < 3; ++i) R[static_cast<std::size_t>(4 * i)] = R[static_cast<std::size_t>(4 * i)] + (1.0 - B * t2);
  R[1] = R[1] - A * th[2];
  R[2] = R[2] + A * th[1];
  R[3] = R[3] + A * th[2];
  R[5] = R[5] - A * th[0];
  R[6] = R[6] - A * th[1];
  R[7] = R[7] + A * th[0];
  return R;
}

// Rotation matrix -> vector: with ax = axial(R - R^T) = 2 sin t a and
// cos t = (tr R - 1)/2, log R = (t / sin t) ax / 2, t = atan2(sin t, cos t)
// (well conditioned everywhere, unlike acos near the identity); t / sin t by
// its series in sin^2 t near the identity.
template <class T>
V3t<T> log_t(const M3t<T>& R) {
  using std::atan2;
  using std::sqrt;
  const V3t<T> ax{R[7] - R[5], R[2] - R[6], R[3] - R[1]};
  const T s2 = dot3(ax, ax) * 0.25;
  const T c = (R[0] + R[4] + R[8] - 1.0) * 0.5;
  T f = s2;
  if (ad::val(c) > 0 && ad::val(s2) < 1e-8) {
    f = 1.0 + s2 / 6.0 + s2 * s2 * (3.0 / 40.0);
  } else {
    const T sn = sqrt(s2);
    require(ad::val(sn) > 1e-12, "corot: deformational rotation near pi");
    f = atan2(sn, c) / sn;
  }
  return {f * ax[0] * 0.5, f * ax[1] * 0.5, f * ax[2] * 0.5};
}

struct CorotElem {
  std::size_t bar_ix;
  int ni, nj;
  Mat3 E0;          // columns = undeformed local axes
  // undeformed FLEXIBLE end positions (node + rigid arm) and the arms
  // (global, undeformed): a rigid arm turns with its node, so the current
  // flexible end is x_node + R(theta) a — exact for finite rotations, the
  // flexible ends' rotations being the nodes' own
  Vec3 X1, X2;
  Vec3 a1 = Vec3::Zero(), a2 = Vec3::Zero();
  double l0;

  // current chord x2 - x1 between the flexible ends
  template <class T>
  V3t<T> chord_t(const std::array<T, 12>& q, const M3t<T>& R1, const M3t<T>& R2) const {
    V3t<T> x21;
    for (int c = 0; c < 3; ++c) {
      T v = (q[static_cast<std::size_t>(6 + c)] - q[static_cast<std::size_t>(c)]) + (X2(c) - X1(c));
      for (int k = 0; k < 3; ++k) {
        if (a2(k) != 0.0) v = v + R2[static_cast<std::size_t>(3 * c + k)] * a2(k);
        if (a1(k) != 0.0) v = v - R1[static_cast<std::size_t>(3 * c + k)] * a1(k);
      }
      x21[static_cast<std::size_t>(c)] = v - (a2(c) - a1(c));
    }
    return x21;
  }
  // current arm R(theta) a of end k (0, 1)
  template <class T>
  V3t<T> arm_t(int k, const M3t<T>& R) const {
    const Vec3& a = k == 0 ? a1 : a2;
    V3t<T> r;
    for (int c = 0; c < 3; ++c)
      r[static_cast<std::size_t>(c)] = R[static_cast<std::size_t>(3 * c)] * a(0) + R[static_cast<std::size_t>(3 * c + 1)] * a(1) +
                                       R[static_cast<std::size_t>(3 * c + 2)] * a(2);
    return r;
  }

  // ---- Winkler bed in large displacements: ground springs ----
  // The bed is attached to the GROUND along the undeformed flexible axis and
  // reacts to the absolute displacement of the axis, in the undeformed local
  // axes (moduli kx, ky, kz), and to the end twists phi = e1_0 . theta about
  // the original axis (2x2 Kphi). The axis of the corotational element is its
  // current chord plus the Hermite deflection of its deformational end
  // rotations, so the energy
  //   W = 1/2 int (kx ux^2 + ky uy^2 + kz uz^2) dx + 1/2 phi^T Kphi phi
  // linearises EXACTLY to keta's consistent Hermite bed (NtN) of linear
  // statics; integrated by 4-point Gauss (exact for that limit), its
  // gradient and Hessian by the same AD as the element.
  Vec3 kbed = Vec3::Zero();
  Eigen::Matrix2d Kphi = Eigen::Matrix2d::Zero();
  bool bed_in_energy = true;  // false: the caller adds bed_grad / bed_hess
  bool has_bed() const { return kbed.squaredNorm() > 0 || Kphi.squaredNorm() > 0; }

  template <class T>
  T bed_energy_t(const std::array<T, 12>& q) const {
    auto cst = [](double x) {
      if constexpr (std::is_same_v<T, double>)
        return x;
      else
        return T::constant(x);
    };
    T w = cst(0.0);
    if (kbed.squaredNorm() > 0) {
      const M3t<T> R1 = exp_t<T>({q[3], q[4], q[5]}), R2 = exp_t<T>({q[9], q[10], q[11]});
      const M3t<T> E = frame_t(q, R1, R2);
      const V3t<T> x21 = chord_t(q, R1, R2);
      const std::array<T, 7> d = deform_t(q);
      V3t<T> u1;  // displacement of the flexible end i
      const V3t<T> r1a = arm_t(0, R1);
      for (int c = 0; c < 3; ++c) u1[static_cast<std::size_t>(c)] = q[static_cast<std::size_t>(c)] + (r1a[static_cast<std::size_t>(c)] - a1(c));
      static const double gx[4] = {0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970263};
      static const double gw[4] = {0.1739274225687269, 0.3260725774312731, 0.3260725774312731, 0.1739274225687269};
      const Vec3 X21 = X2 - X1;
      for (int g = 0; g < 4; ++g) {
        const double xi = gx[g];
        const double N2 = xi * (1 - xi) * (1 - xi), N4 = -xi * xi * (1 - xi);
        const T vl = (d[3] * N2 + d[6] * N4) * l0;      // rz = v'
        const T wl = (d[2] * N2 + d[5] * N4) * (-l0);   // ry = -w'
        V3t<T> ua;  // absolute displacement of the axis point
        for (int c = 0; c < 3; ++c)
          ua[static_cast<std::size_t>(c)] = u1[static_cast<std::size_t>(c)] + x21[static_cast<std::size_t>(c)] * xi - X21(c) * xi +
                                            E[static_cast<std::size_t>(3 * c + 1)] * vl + E[static_cast<std::size_t>(3 * c + 2)] * wl;
        for (int a = 0; a < 3; ++a) {
          if (kbed(a) == 0.0) continue;
          // component on the undeformed local axis a (E0 column a)
          const T ul = ua[0] * E0(0, a) + ua[1] * E0(1, a) + ua[2] * E0(2, a);
          w = w + ul * ul * (0.5 * kbed(a) * gw[g] * l0);
        }
      }
    }
    if (Kphi.squaredNorm() > 0) {
      const T p1 = q[3] * E0(0, 0) + q[4] * E0(1, 0) + q[5] * E0(2, 0);
      const T p2 = q[9] * E0(0, 0) + q[10] * E0(1, 0) + q[11] * E0(2, 0);
      w = w + p1 * p1 * (0.5 * Kphi(0, 0)) + p1 * p2 * Kphi(0, 1) + p2 * p2 * (0.5 * Kphi(1, 1));
    }
    return w;
  }
  Eigen::Matrix<double, 12, 1> bed_grad(const Eigen::Matrix<double, 12, 1>& q) const {
    std::array<ad::Grad<12>, 12> x;
    for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(q(i), i);
    const ad::Grad<12> w = bed_energy_t(x);
    Eigen::Matrix<double, 12, 1> g;
    for (int i = 0; i < 12; ++i) g(i) = w.g[static_cast<std::size_t>(i)];
    return g;
  }
  Eigen::Matrix<double, 12, 12> bed_hess(const Eigen::Matrix<double, 12, 1>& q) const {
    std::array<ad::Jet<12>, 12> x;
    for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Jet<12>::seed(q(i), i);
    const ad::Jet<12> w = bed_energy_t(x);
    Eigen::Matrix<double, 12, 12> K;
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) K(i, j) = w.h[static_cast<std::size_t>(12 * i + j)];
    return K;
  }
  Eigen::Matrix<double, 7, 7> Kbar;  // condensed local stiffness at [ux_j, rot_i, rot_j]
  CorotRotationMeasure measure = CorotRotationMeasure::Log;

  // Current corotated frame (columns = local axes): chord e1 + mean-director
  // construction — the SAME frame the deformational kinematics uses, so load
  // transfers and strain measures agree.
  template <class T>
  M3t<T> frame_t(const std::array<T, 12>& q, const M3t<T>& R1, const M3t<T>& R2) const {
    using std::sqrt;
    V3t<T> p;
    const V3t<T> x21 = chord_t(q, R1, R2);
    const T ln = sqrt(dot3(x21, x21));
    require(ad::val(ln) > 1e-9 * l0, "corot: element collapsed");
    const V3t<T> e1{x21[0] / ln, x21[1] / ln, x21[2] / ln};
    const Vec3 t0 = E0.col(1);
    for (int i = 0; i < 3; ++i) {
      T a = R1[static_cast<std::size_t>(3 * i)] * t0(0) + R1[static_cast<std::size_t>(3 * i + 1)] * t0(1) + R1[static_cast<std::size_t>(3 * i + 2)] * t0(2);
      T b = R2[static_cast<std::size_t>(3 * i)] * t0(0) + R2[static_cast<std::size_t>(3 * i + 1)] * t0(1) + R2[static_cast<std::size_t>(3 * i + 2)] * t0(2);
      p[static_cast<std::size_t>(i)] = (a + b) * 0.5;
    }
    V3t<T> e3 = cross3(e1, p);
    const T n3 = sqrt(dot3(e3, e3));
    require(ad::val(n3) > 1e-9, "corot: degenerate frame (director parallel to chord)");
    for (auto& x : e3) x = x / n3;
    const V3t<T> e2 = cross3(e3, e1);
    M3t<T> E;
    for (int i = 0; i < 3; ++i) {
      E[static_cast<std::size_t>(3 * i)] = e1[static_cast<std::size_t>(i)];
      E[static_cast<std::size_t>(3 * i + 1)] = e2[static_cast<std::size_t>(i)];
      E[static_cast<std::size_t>(3 * i + 2)] = e3[static_cast<std::size_t>(i)];
    }
    return E;
  }

  // Deformational vector d(q) — the exact corotational kinematics:
  // [stretch ln - l0 (cancellation-free), log(E^T R1 E0), log(E^T R2 E0)].
  template <class T>
  std::array<T, 7> deform_t(const std::array<T, 12>& q) const {
    using std::sqrt;
    const M3t<T> R1 = exp_t<T>({q[3], q[4], q[5]}), R2 = exp_t<T>({q[9], q[10], q[11]});
    const M3t<T> E = frame_t(q, R1, R2);
    const V3t<T> x21 = chord_t(q, R1, R2);
    const T sq = dot3(x21, x21);
    std::array<T, 7> d;
    d[0] = (sq - l0 * l0) / (sqrt(sq) + l0);
    auto rel = [&](const M3t<T>& R) {  // E^T R E0
      M3t<T> A, B;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          A[static_cast<std::size_t>(3 * i + j)] = E[static_cast<std::size_t>(i)] * R[static_cast<std::size_t>(j)] + E[static_cast<std::size_t>(3 + i)] * R[static_cast<std::size_t>(3 + j)] + E[static_cast<std::size_t>(6 + i)] * R[static_cast<std::size_t>(6 + j)];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          B[static_cast<std::size_t>(3 * i + j)] = A[static_cast<std::size_t>(3 * i)] * E0(0, j) + A[static_cast<std::size_t>(3 * i + 1)] * E0(1, j) + A[static_cast<std::size_t>(3 * i + 2)] * E0(2, j);
      if (measure == CorotRotationMeasure::Sin) {
        // Sin-form extraction s = axial(Rbar - Rbar^T)/2 = sin(t) a, with the
        // small-local-rotation guard |s| <= 0.7.
        const V3t<T> sv{(B[7] - B[5]) * 0.5, (B[2] - B[6]) * 0.5, (B[3] - B[1]) * 0.5};
        require(ad::val(dot3(sv, sv)) <= 0.49,
                "corot: local rotation beyond the sin-form guard |s| <= 0.7 (refine the mesh or add substeps)");
        return sv;
      }
      return log_t(B);
    };
    const V3t<T> r1 = rel(R1), r2 = rel(R2);
    for (int k = 0; k < 3; ++k) {
      d[static_cast<std::size_t>(1 + k)] = r1[static_cast<std::size_t>(k)];
      d[static_cast<std::size_t>(4 + k)] = r2[static_cast<std::size_t>(k)];
    }
    return d;
  }

  template <class T>
  T energy_t(const std::array<T, 12>& q) const {
    const std::array<T, 7> d = deform_t(q);
    T w = d[0] * (0.5 * Kbar(0, 0)) * d[0];
    for (int i = 0; i < 7; ++i)
      for (int j = 0; j < 7; ++j)
        if (i != 0 || j != 0) w = w + d[static_cast<std::size_t>(i)] * (0.5 * Kbar(i, j)) * d[static_cast<std::size_t>(j)];
    if (bed_in_energy && has_bed()) w = w + bed_energy_t(q);
    return w;
  }

  static std::array<double, 12> arr(const Eigen::Matrix<double, 12, 1>& q) {
    std::array<double, 12> a;
    for (int i = 0; i < 12; ++i) a[static_cast<std::size_t>(i)] = q(i);
    return a;
  }

  Mat3 frame(const Eigen::Matrix<double, 12, 1>& q) const {
    const std::array<double, 12> a = arr(q);
    const M3t<double> E = frame_t(a, exp_t<double>({a[3], a[4], a[5]}), exp_t<double>({a[9], a[10], a[11]}));
    Mat3 M;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) M(i, j) = E[static_cast<std::size_t>(3 * i + j)];
    return M;
  }

  Eigen::Matrix<double, 7, 1> deform(const Eigen::Matrix<double, 12, 1>& q) const {
    const std::array<double, 7> d = deform_t(arr(q));
    Eigen::Matrix<double, 7, 1> v;
    for (int i = 0; i < 7; ++i) v(i) = d[static_cast<std::size_t>(i)];
    return v;
  }

  double energy(const Eigen::Matrix<double, 12, 1>& q) const { return energy_t(arr(q)); }

  // current chord between the flexible ends (double)
  Vec3 chord(const Eigen::Matrix<double, 12, 1>& q) const {
    const std::array<double, 12> a = arr(q);
    const V3t<double> x = chord_t(a, exp_t<double>({a[3], a[4], a[5]}), exp_t<double>({a[9], a[10], a[11]}));
    return Vec3(x[0], x[1], x[2]);
  }

  // Internal force: the exact energy gradient, one first-order AD pass.
  Eigen::Matrix<double, 12, 1> grad(const Eigen::Matrix<double, 12, 1>& q) const {
    std::array<ad::Grad<12>, 12> x;
    for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(q(i), i);
    const ad::Grad<12> w = energy_t(x);
    Eigen::Matrix<double, 12, 1> g;
    for (int i = 0; i < 12; ++i) g(i) = w.g[static_cast<std::size_t>(i)];
    return g;
  }

  // Consistent tangent: the exact energy Hessian, one second-order AD pass.
  Eigen::Matrix<double, 12, 12> hess(const Eigen::Matrix<double, 12, 1>& q) const {
    std::array<ad::Jet<12>, 12> x;
    for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Jet<12>::seed(q(i), i);
    const ad::Jet<12> w = energy_t(x);
    Eigen::Matrix<double, 12, 12> K;
    for (int i = 0; i < 12; ++i)
      for (int j = 0; j < 12; ++j) K(i, j) = w.h[static_cast<std::size_t>(12 * i + j)];
    return K;
  }
};

// ---------------- corotational facet shell (EICR, energy formulation) ----------------
//
// Element-independent corotational triangle (Rankin & Nour-Omid, Comput.
// Struct. 30 (1988) 257-267; Battini & Pacoste, CMAME 191 (2002) 1755-1789):
// the rigid motion of the facet is the frame of its CURRENT triangle —
// A(x) = [a1, a3 x a1, a3], a1 along side 1-2, a3 the unit normal — carried
// onto the undeformed local axes, E = A(x) A(X)^T E0; the deformational
// vector, in that frame, holds the centroid-relative translations
// E^T (x_i - x_c) - E0^T (X_i - X_c) and the rotations log(E^T R_i E0). The
// energy W = 1/2 d^T K_l d uses the UNCHANGED linear facet stiffness
// (Allman membrane + MITC3 plate, local frame), whose null space holds the
// rigid modes, so W linearises to 1/2 u^T ke u. Written once on the scalar
// type: internal force and consistent tangent are its exact gradient and
// Hessian by AD (ad::Grad<18>, ad::Jet<18>). Small strains, arbitrary
// rotations (finite rotations of the triangle are exact; the strain measure
// is the linear one in the corotated frame, converging with refinement).
struct CorotShell {
  std::array<int, 3> nodes{};
  std::array<Vec3, 3> X;
  Mat3 E0;   // columns = undeformed local axes
  Mat3 A0;   // columns = auxiliary frame of the undeformed triangle
  Mat18 Kl;  // local linear stiffness

  template <class T>
  static M3t<T> aux_frame(const std::array<V3t<T>, 3>& x) {
    using std::sqrt;
    V3t<T> e12, e13;
    for (int c = 0; c < 3; ++c) {
      e12[static_cast<std::size_t>(c)] = x[1][static_cast<std::size_t>(c)] - x[0][static_cast<std::size_t>(c)];
      e13[static_cast<std::size_t>(c)] = x[2][static_cast<std::size_t>(c)] - x[0][static_cast<std::size_t>(c)];
    }
    const T l1 = sqrt(dot3(e12, e12));
    V3t<T> a3 = cross3(e12, e13);
    const T n3 = sqrt(dot3(a3, a3));
    require(ad::val(n3) > 1e-14 * ad::val(l1) * ad::val(l1), "corot shell: degenerate facet");
    V3t<T> a1;
    for (int c = 0; c < 3; ++c) {
      a1[static_cast<std::size_t>(c)] = e12[static_cast<std::size_t>(c)] / l1;
      a3[static_cast<std::size_t>(c)] = a3[static_cast<std::size_t>(c)] / n3;
    }
    const V3t<T> a2 = cross3(a3, a1);
    M3t<T> A;
    for (int i = 0; i < 3; ++i) {
      A[static_cast<std::size_t>(3 * i)] = a1[static_cast<std::size_t>(i)];
      A[static_cast<std::size_t>(3 * i + 1)] = a2[static_cast<std::size_t>(i)];
      A[static_cast<std::size_t>(3 * i + 2)] = a3[static_cast<std::size_t>(i)];
    }
    return A;
  }

  template <class T>
  std::array<T, 18> deform_t(const std::array<T, 18>& q) const {
    std::array<V3t<T>, 3> x;
    for (int k = 0; k < 3; ++k)
      for (int c = 0; c < 3; ++c)
        x[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)] = q[static_cast<std::size_t>(6 * k + c)] + X[static_cast<std::size_t>(k)](c);
    const M3t<T> A = aux_frame(x);
    // E = A (A0^T E0): a constant 3x3 on the right
    const Mat3 C = A0.transpose() * E0;
    M3t<T> E;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        E[static_cast<std::size_t>(3 * i + j)] = A[static_cast<std::size_t>(3 * i)] * C(0, j) + A[static_cast<std::size_t>(3 * i + 1)] * C(1, j) + A[static_cast<std::size_t>(3 * i + 2)] * C(2, j);
    V3t<T> xc;
    for (int c = 0; c < 3; ++c)
      xc[static_cast<std::size_t>(c)] = (x[0][static_cast<std::size_t>(c)] + x[1][static_cast<std::size_t>(c)] + x[2][static_cast<std::size_t>(c)]) / 3.0;
    const Vec3 Xc = (X[0] + X[1] + X[2]) / 3.0;
    std::array<T, 18> d;
    for (int k = 0; k < 3; ++k) {
      const Vec3 L0 = E0.transpose() * (X[static_cast<std::size_t>(k)] - Xc);
      for (int i = 0; i < 3; ++i) {
        T v = E[static_cast<std::size_t>(i)] * (x[static_cast<std::size_t>(k)][0] - xc[0]) +
              E[static_cast<std::size_t>(3 + i)] * (x[static_cast<std::size_t>(k)][1] - xc[1]) +
              E[static_cast<std::size_t>(6 + i)] * (x[static_cast<std::size_t>(k)][2] - xc[2]);
        d[static_cast<std::size_t>(6 * k + i)] = v - L0(i);
      }
      const M3t<T> R = exp_t<T>({q[static_cast<std::size_t>(6 * k + 3)], q[static_cast<std::size_t>(6 * k + 4)], q[static_cast<std::size_t>(6 * k + 5)]});
      M3t<T> Am, B;  // E^T R, then (E^T R) E0
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          Am[static_cast<std::size_t>(3 * i + j)] = E[static_cast<std::size_t>(i)] * R[static_cast<std::size_t>(j)] + E[static_cast<std::size_t>(3 + i)] * R[static_cast<std::size_t>(3 + j)] + E[static_cast<std::size_t>(6 + i)] * R[static_cast<std::size_t>(6 + j)];
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          B[static_cast<std::size_t>(3 * i + j)] = Am[static_cast<std::size_t>(3 * i)] * E0(0, j) + Am[static_cast<std::size_t>(3 * i + 1)] * E0(1, j) + Am[static_cast<std::size_t>(3 * i + 2)] * E0(2, j);
      const V3t<T> th = log_t(B);
      for (int i = 0; i < 3; ++i) d[static_cast<std::size_t>(6 * k + 3 + i)] = th[static_cast<std::size_t>(i)];
    }
    return d;
  }

  template <class T>
  T energy_t(const std::array<T, 18>& q) const {
    const std::array<T, 18> d = deform_t(q);
    auto row = [&](int i) {  // (K_l d)_i
      T r = d[0] * Kl(i, 0);
      for (int j = 1; j < 18; ++j)
        if (Kl(i, j) != 0.0) r = r + d[static_cast<std::size_t>(j)] * Kl(i, j);
      return r;
    };
    T w = row(0) * (0.5 * d[0]);
    for (int i = 1; i < 18; ++i) w = w + row(i) * (0.5 * d[static_cast<std::size_t>(i)]);
    return w;
  }

  Eigen::Matrix<double, 18, 1> grad(const Eigen::Matrix<double, 18, 1>& q) const {
    std::array<ad::Grad<18>, 18> x;
    for (int i = 0; i < 18; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<18>::seed(q(i), i);
    const ad::Grad<18> w = energy_t(x);
    Eigen::Matrix<double, 18, 1> g;
    for (int i = 0; i < 18; ++i) g(i) = w.g[static_cast<std::size_t>(i)];
    return g;
  }
  Eigen::Matrix<double, 18, 18> hess(const Eigen::Matrix<double, 18, 1>& q) const {
    std::array<ad::Jet<18>, 18> x;
    for (int i = 0; i < 18; ++i) x[static_cast<std::size_t>(i)] = ad::Jet<18>::seed(q(i), i);
    const ad::Jet<18> w = energy_t(x);
    Eigen::Matrix<double, 18, 18> K;
    for (int i = 0; i < 18; ++i)
      for (int j = 0; j < 18; ++j) K(i, j) = w.h[static_cast<std::size_t>(18 * i + j)];
    return K;
  }
};

inline CorotShell make_corot_shell(const ShellElem& se) {
  CorotShell cs;
  cs.nodes = se.nodes;
  cs.X = se.el.coords();
  cs.E0 = se.el.R().transpose();
  std::array<V3t<double>, 3> x;
  for (int k = 0; k < 3; ++k)
    for (int c = 0; c < 3; ++c) x[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)] = cs.X[static_cast<std::size_t>(k)](c);
  const M3t<double> A = CorotShell::aux_frame(x);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) cs.A0(i, j) = A[static_cast<std::size_t>(3 * i + j)];
  // the local facet is linear (offset, varying thickness, orthotropy ride on
  // the corotated frame); a Winkler bed is fixed in space, not in that frame
  require(se.el.enrichment().kbed == 0, "large displacements: a shell's Winkler bed is not corotational "
                                        "(use spring supports or a bed of links)");
  const Mat18 T = se.el.local_transform();
  cs.Kl = T * se.el.ke() * T.transpose();
  return cs;
}

inline CorotElem make_corot_elem(const Model& model, const Assembled& asm_, std::size_t e,
                          CorotRotationMeasure measure) {
  const Bar& b = model.bars[e];
  CorotElem ce;
  ce.measure = measure;
  ce.bar_ix = e;
  ce.ni = model.node_index(b.ni);
  ce.nj = model.node_index(b.nj);
  ce.E0 = asm_.elements[e].R().transpose();  // axes of the flexible part (offset ends)
  ce.a1 = b.ecc[0];
  ce.a2 = b.ecc[1];
  ce.X1 = model.node(b.ni).x + ce.a1;
  ce.X2 = model.node(b.nj).x + ce.a2;
  ce.l0 = asm_.elements[e].L();
  // Kbar: the local block restricted to the deformational DOFs (ux_j and the
  // end rotations, the corotated chord carrying the rest) — exact for any
  // local kernel whose null space holds the rigid modes: the release-condensed
  // prismatic block (releases keep that null space) or a tapered member's
  // exact block (CARA VAR). A custom block with a foundation (bef.h) has no
  // rigid modes and is rejected.
  static const int LRED[7] = {6, 3, 4, 5, 9, 10, 11};
  const TaperedBar* tb = model.tapered_of(b.id);
  require(tb == nullptr || !tb->f0_span,
          "corotational: a custom member block (beam on elastic foundation) is not supported");
  const Mat12& kl = tb ? tb->ke : asm_.elements[e].k_local();
  // the Winkler bed (as keta's Hermite bed: positive moduli)
  const SoilProps& so = b.soil;
  ce.kbed = Vec3(std::max(0.0, so.kfx), std::max(0.0, so.kfy), std::max(0.0, so.kfz));
  if (so.kmx > 0) ce.Kphi << so.kmx * ce.l0 / 3, so.kmx * ce.l0 / 6, so.kmx * ce.l0 / 6, so.kmx * ce.l0 / 3;
  for (int a = 0; a < 7; ++a)
    for (int bb = 0; bb < 7; ++bb) ce.Kbar(a, bb) = kl(LRED[a], LRED[bb]);
  return ce;
}

// Generic left Jacobian of the exponential map, J = I + c2 S + c3 S^2 with
// c2 = (1 - cos t)/t^2, c3 = (t - sin t)/t^3 (Taylor series below t^2 = 1e-8).
template <class T>
M3t<T> left_jacobian_t(const V3t<T>& th) {
  using std::cos;
  using std::sin;
  using std::sqrt;
  const T t2 = dot3(th, th);
  T c2 = t2, c3 = t2;
  if (ad::val(t2) < 1e-8) {
    c2 = 0.5 - t2 / 24.0 + t2 * t2 / 720.0;
    c3 = 1.0 / 6.0 - t2 / 120.0 + t2 * t2 / 5040.0;
  } else {
    const T t = sqrt(t2);
    c2 = (1.0 - cos(t)) / t2;
    c3 = (t - sin(t)) / (t2 * t);
  }
  M3t<T> J;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) J[static_cast<std::size_t>(3 * i + j)] = c3 * th[static_cast<std::size_t>(i)] * th[static_cast<std::size_t>(j)];
  for (int i = 0; i < 3; ++i) J[static_cast<std::size_t>(4 * i)] = J[static_cast<std::size_t>(4 * i)] + (1.0 - c3 * t2);
  J[1] = J[1] - c2 * th[2];
  J[2] = J[2] + c2 * th[1];
  J[3] = J[3] + c2 * th[2];
  J[5] = J[5] - c2 * th[0];
  J[6] = J[6] - c2 * th[1];
  J[7] = J[7] + c2 * th[0];
  return J;
}

// J(theta)^T m: the rotation-vector generalized force of a spatial moment m.
template <class T, class M>
V3t<T> jt_times(const V3t<T>& th, const V3t<M>& m) {
  const M3t<T> J = left_jacobian_t(th);
  V3t<T> r;
  for (int j = 0; j < 3; ++j)
    r[static_cast<std::size_t>(j)] = J[static_cast<std::size_t>(j)] * m[0] + J[static_cast<std::size_t>(3 + j)] * m[1] + J[static_cast<std::size_t>(6 + j)] * m[2];
  return r;
}

// One span load carried on a corotational element. The equivalent nodal
// forces are re-lowered on the CURRENT corotated frame at every residual
// evaluation ("moments d'encastrement transportes"):
//   * LoadFrame::Global (and SelfWeight): the load vector KEEPS its global
//     direction — it is expressed in the current frame (q_loc = E^T q_glob),
//     lowered with the element's exact static shape functions on the
//     UNDEFORMED length, then pushed back through E. Conservative dead load
//     (gravity!), exact resultant at any rotation (E E^T = I).
//   * LoadFrame::Local: a FOLLOWER load — the local intensities are constant
//     in the corotated frame, so the lowered vector is precomputed once and
//     only the push-back rotates. Documented choice: local = follower.
// Equivalent-load moments are conjugate to the SPIN, hence the J(theta)^T
// transfer, exactly like the fixed nodal moments. The load-stiffness term
// (derivative of these forces w.r.t. u) is IN the tangent: the lowering is
// written generically (bar_load_vector below) and differentiated by
// AD, as is J(theta)^T M for the nodal moments (CorotSystem::load_stiffness).
// The lowering itself is the standard corotational equivalent-nodal-load
// treatment (e.g. Crisfield, Non-linear Finite Element Analysis of Solids and
// Structures, vol. 2 ch. 17): exact up to the element discretization, like
// the consistent mass.
struct CorotBarLoad {
  std::size_t ei = 0;     // index into CorotSystem::elems
  Load load;              // the span load (frame Global or Local)
  bool follower = false;  // LoadFrame::Local
  Vec12 f0_follower = Vec12::Zero();  // condensed local FEF (follower loads)
  // Condensed local FEF per UNIT local intensity component, columns
  // [F_x, F_y, F_z, M_x, M_y, M_z, F2_x, F2_y, F2_z]: the lowering is linear in
  // the intensities, so floc = basis * v_local — a closed form the AD jets can
  // differentiate (the load stiffness). Equal to the direct lowering up to
  // round-off.
  Eigen::Matrix<double, 12, 9> basis = Eigen::Matrix<double, 12, 9>::Zero();
};

// The corotational record of one span load on bar e (ei: its element index
// in the caller's list): local fixed-end forces per unit local intensity —
// the tapered member's exact ones (CARA VAR) or the prismatic shape-function
// integrals — release-condensed. A linear load's distributed moment must be
// constant (M2 = M: the basis carries one moment intensity).
inline CorotBarLoad make_corot_bar_load(const Model& model, const Assembled& asm_, std::size_t e, std::size_t ei,
                                        const Load& l) {
  CorotBarLoad bl;
  bl.ei = ei;
  bl.load = l;
  bl.follower = l.frame == LoadFrame::Local;
  const Bar& b = model.bars[e];
  require(l.kind != Load::Kind::BarLinear || l.M2 == l.M,
          "corotational: a linear bar load with a varying distributed moment is not supported");
  const TaperedBar* tb = model.tapered_of(b.id);
  require(tb == nullptr || static_cast<bool>(tb->f0),
          "corotational: span load on a tapered bar without its equivalent loads (f0)");
  auto fef = [&](const Load& ll) {
    return asm_.elements[e].condense_f0(tb ? tb->f0(ll) : asm_.elements[e].f0_local(ll, model.material(b.material)));
  };
  if (bl.follower) bl.f0_follower = fef(l);
  for (int k = 0; k < 9; ++k) {
    Load unit = l;
    unit.frame = LoadFrame::Local;
    unit.F = unit.M = unit.F2 = Vec3::Zero();
    (k < 3 ? unit.F : k < 6 ? unit.M : unit.F2)(k % 3) = 1.0;
    unit.M2 = unit.M;
    bl.basis.col(k) = fef(unit);
  }
  return bl;
}

// The corotational record of a bar's SELF-WEIGHT rho A(x) g: a dead load
// whose local intensity is E^T g, its fixed-end forces per unit local
// component those of the self-weight load itself (rho A(x) e_c — the
// tapered member's exact variable weight, CARA VAR, or the prismatic one),
// release-condensed.
inline CorotBarLoad make_corot_self_weight(const Model& model, const Assembled& asm_, std::size_t e, std::size_t ei,
                                           const Vec3& g) {
  CorotBarLoad bl;
  bl.ei = ei;
  const Bar& b = model.bars[e];
  bl.load = Load::udl(b.id, g, LoadFrame::Global);
  bl.follower = false;
  const TaperedBar* tb = model.tapered_of(b.id);
  require(tb == nullptr || static_cast<bool>(tb->f0),
          "corotational: self-weight of a tapered bar without its equivalent loads (f0)");
  const Mat3& R = asm_.elements[e].R();
  for (int c = 0; c < 3; ++c) {
    // the gravity vector whose local components are e_c
    const Load sw = Load::self_weight(Vec3(R.row(c).transpose()));
    bl.basis.col(c) =
        asm_.elements[e].condense_f0(tb ? tb->f0(sw) : asm_.elements[e].f0_local(sw, model.material(b.material)));
  }
  return bl;
}

// The 9 LOCAL intensities [F, M, F2] of a span load at element DOFs q:
// E(q)^T g for a dead (global) load, constant for a follower (local) one.
template <class T>
std::array<T, 9> local_intensities(const CorotElem& ce, const CorotBarLoad& bl, const std::array<T, 12>& q,
                                   const M3t<T>& E) {
  const Vec3* glob[3] = {&bl.load.F, &bl.load.M, &bl.load.F2};
  std::array<T, 9> v;
  for (int k = 0; k < 9; ++k) {
    const Vec3& g = *glob[k / 3];
    const int col = k % 3;
    if (bl.follower) {
      if constexpr (std::is_same_v<T, double>)
        v[static_cast<std::size_t>(k)] = g(col);
      else
        v[static_cast<std::size_t>(k)] = T::constant(g(col));
    } else {
      v[static_cast<std::size_t>(k)] = E[static_cast<std::size_t>(col)] * g(0) + E[static_cast<std::size_t>(3 + col)] * g(1) + E[static_cast<std::size_t>(6 + col)] * g(2);
    }
  }
  (void)ce;
  (void)q;
  return v;
}

// Local end forces floc (corotated frame) pushed to the element DOFs:
// forces through E, the flexible ends' moments plus the CURRENT arm moment
// (R a) x F, as rotation-vector generalized forces J(theta)^T M.
template <class T>
std::array<T, 12> push_local(const CorotElem& ce, const std::array<T, 12>& q, const M3t<T>& R1, const M3t<T>& R2,
                             const M3t<T>& E, const std::array<T, 12>& floc) {
  auto cst = [](double x) {
    if constexpr (std::is_same_v<T, double>)
      return x;
    else
      return T::constant(x);
  };
  std::array<T, 12> out;
  for (int a = 0; a < 2; ++a) {
    V3t<T> mt, ft3;
    for (int i = 0; i < 3; ++i) {
      T ft = cst(0.0), m = cst(0.0);
      for (int c = 0; c < 3; ++c) {
        ft = ft + E[static_cast<std::size_t>(3 * i + c)] * floc[static_cast<std::size_t>(6 * a + c)];
        m = m + E[static_cast<std::size_t>(3 * i + c)] * floc[static_cast<std::size_t>(6 * a + 3 + c)];
      }
      out[static_cast<std::size_t>(6 * a + i)] = ft;
      mt[static_cast<std::size_t>(i)] = m;
      ft3[static_cast<std::size_t>(i)] = ft;
    }
    // rigid arm: the flexible end's force moves to the node with the moment
    // of the CURRENT arm, (R a) x F
    if ((a == 0 ? ce.a1 : ce.a2).squaredNorm() > 0) {
      const V3t<T> arm = ce.arm_t(a, a == 0 ? R1 : R2);
      const V3t<T> am = cross3(arm, ft3);
      for (int i = 0; i < 3; ++i) mt[static_cast<std::size_t>(i)] = mt[static_cast<std::size_t>(i)] + am[static_cast<std::size_t>(i)];
    }
    const V3t<T> th{q[static_cast<std::size_t>(6 * a + 3)], q[static_cast<std::size_t>(6 * a + 4)], q[static_cast<std::size_t>(6 * a + 5)]};
    const V3t<T> gm = jt_times(th, mt);
    for (int i = 0; i < 3; ++i) out[static_cast<std::size_t>(6 * a + 3 + i)] = gm[static_cast<std::size_t>(i)];
  }
  return out;
}

// Equivalent nodal forces of one span load at element DOFs q (global
// forces; moments as rotation-vector generalized forces J(theta)^T M),
// generic on the scalar type so its Jacobian — the LOAD STIFFNESS — comes
// out of one AD pass. Global intensities keep their spatial direction (the
// local intensities are E^T g); follower (Local) ones are constant in the
// corotated frame.
template <class T>
std::array<T, 12> bar_load_vector(const CorotElem& ce, const CorotBarLoad& bl, const std::array<T, 12>& q) {
  auto cst = [](double x) {
    if constexpr (std::is_same_v<T, double>)
      return x;
    else
      return T::constant(x);
  };
  const M3t<T> R1 = exp_t<T>({q[3], q[4], q[5]}), R2 = exp_t<T>({q[9], q[10], q[11]});
  const M3t<T> E = ce.frame_t(q, R1, R2);
  const std::array<T, 9> v = local_intensities(ce, bl, q, E);
  std::array<T, 12> floc;
  for (int r = 0; r < 12; ++r) {
    T acc = cst(0.0);
    for (int k = 0; k < 9; ++k)
      if (bl.basis(r, k) != 0.0) acc = acc + v[static_cast<std::size_t>(k)] * bl.basis(r, k);
    floc[static_cast<std::size_t>(r)] = acc;
  }
  return push_local(ce, q, R1, R2, E, floc);
}


}  // namespace keta::corot

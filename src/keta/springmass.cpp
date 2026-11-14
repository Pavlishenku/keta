#include "keta/springmass.h"

#include <cmath>

namespace keta {

namespace {

// PLANE DOF order (dx, dy, rz) -> keta DOF indices (ux uy uz rx ry rz).
// ST1 manual v24 p. 47 ("Option PLANE  dx, dy, rz") and p. 51 (the plane
// effort/displacement vectors).
constexpr int kPlaneDof[3] = {0, 1, 5};

// 5-point Gauss-Legendre on [0, 1] — exact through degree 9. The line-mass
// integrand mu(xi) * N_a(xi) * N_b(xi) peaks at degree 2 + 3 + 3 = 8
// (parabolic law x two cubic Hermite functions), so the quadrature is EXACT
// for every law of the vocabulary. Same nodes/weights as element.cpp.
constexpr int kGauss = 5;
constexpr double kGx[kGauss] = {0.046910077030668, 0.230765344947158, 0.5,
                                0.769234655052842, 0.953089922969332};
constexpr double kGw[kGauss] = {0.118463442528095, 0.239314335249683, 0.284444444444444,
                                0.239314335249683, 0.118463442528095};

Mat3 skew(const Vec3& v) {
  Mat3 s;
  s << 0, -v.z(), v.y(), v.z(), 0, -v.x(), -v.y(), v.x(), 0;
  return s;
}

}  // namespace

// ---------------- ST1 lower-triangle / diagonal readers ----------------

Mat6 sym_from_lower_triangle(const std::vector<double>& values, bool plane) {
  const int n = plane ? 3 : 6;
  const std::size_t need = static_cast<std::size_t>(n * (n + 1) / 2);
  require(values.size() == need,
          "ST1 lower triangle: expected " + std::to_string(need) + " terms (" +
              (plane ? "PLANE" : "SPATIALE") + "), got " + std::to_string(values.size()));
  Mat6 k = Mat6::Zero();
  std::size_t t = 0;
  // Row by row: m1 | m2 m3 | m4 m5 m6 | ... (manual p. 50 worked example).
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= i; ++j) {
      double v = values[t++];
      require(finite(v), "ST1 lower triangle: non-finite term " + std::to_string(t));
      int a = plane ? kPlaneDof[i] : i;
      int b = plane ? kPlaneDof[j] : j;
      k(a, b) = v;
      k(b, a) = v;
    }
  }
  return k;
}

Mat6 sym_from_diagonal(const std::vector<double>& values, bool plane) {
  const int n = plane ? 3 : 6;
  require(values.size() == static_cast<std::size_t>(n),
          "ST1 DI: expected " + std::to_string(n) + " terms (" + (plane ? "PLANE" : "SPATIALE") +
              "), got " + std::to_string(values.size()));
  Mat6 k = Mat6::Zero();
  for (int i = 0; i < n; ++i) {
    require(finite(values[static_cast<std::size_t>(i)]),
            "ST1 DI: non-finite term " + std::to_string(i + 1));
    int a = plane ? kPlaneDof[i] : i;
    k(a, a) = values[static_cast<std::size_t>(i)];
  }
  return k;
}

namespace {

Mat6 check_psd_diagonal(Mat6 k, const char* what) {
  for (int c = 0; c < 6; ++c)
    require(k(c, c) >= 0, std::string(what) + ": negative diagonal term on DOF " + std::to_string(c));
  return k;
}

}  // namespace

Mat6 spring_from_lower_triangle(const std::vector<double>& values, bool plane) {
  return check_psd_diagonal(sym_from_lower_triangle(values, plane), "APPUI EL");
}
Mat6 spring_from_diagonal(const std::vector<double>& values, bool plane) {
  return check_psd_diagonal(sym_from_diagonal(values, plane), "APPUI EL DI");
}
Mat6 mass_from_lower_triangle(const std::vector<double>& values, bool plane) {
  return check_psd_diagonal(sym_from_lower_triangle(values, plane), "MASSE");
}
Mat6 mass_from_diagonal(const std::vector<double>& values, bool plane) {
  return check_psd_diagonal(sym_from_diagonal(values, plane), "MASSE DI");
}

// ---------------- ROT alpha ----------------

Mat3 rot_z_frame(double alpha) {
  double c = std::cos(alpha), s = std::sin(alpha);
  Mat3 r;
  // Rows = support-frame axes in global components: v_support = R v_global.
  r << c, s, 0, -s, c, 0, 0, 0, 1;
  return r;
}

Mat6 rotate_support_stiffness(const Mat6& k_support, double alpha) {
  if (alpha == 0.0) return k_support;
  Mat3 rz = rot_z_frame(alpha);
  Mat6 r = Mat6::Zero();
  r.block<3, 3>(0, 0) = rz;  // translations
  r.block<3, 3>(3, 3) = rz;  // rotations (same axial rotation)
  return r.transpose() * k_support * r;
}

Mat6 support_stiffness_global(const Support& s) {
  if (!s.spring6) return Mat6::Zero();
  Mat6 k = rotate_support_stiffness(*s.spring6, s.rot);
  return 0.5 * (k + k.transpose());  // symmetrize defensively
}

bool has_full_spring(const Model& model) {
  for (const Support& s : model.supports)
    if (s.spring6 && s.spring6->squaredNorm() > 0) return true;
  return false;
}

bool has_decol(const Model& model) {
  for (const Support& s : model.supports)
    if (s.decol_axis >= 0) return true;
  return false;
}

void scatter_full_springs(const Assembled& asm_, std::vector<Triplet>& trips) {
  const DofMap& d = asm_.dofs;
  for (const auto& [ni, k6] : asm_.spring6) {
    int g[6];
    for (int c = 0; c < 6; ++c) g[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
    for (int a = 0; a < 6; ++a) {
      if (g[a] < 0) continue;  // blocked or auto-fixed: eliminated like the diagonal path
      for (int b = 0; b < 6; ++b) {
        if (g[b] < 0) continue;
        double v = k6(a, b);
        if (v != 0.0) trips.emplace_back(g[a], g[b], v);
      }
    }
  }
}

// ---------------- MASSE PROPRE: line masses ----------------

double line_mass_at(const BarMass& bm, double xi) {
  switch (bm.kind) {
    case BarMassKind::Uni:
      return bm.m1;
    case BarMassKind::Lin:
      return bm.m1 + (bm.m2 - bm.m1) * xi;
    case BarMassKind::Parabol:
      // Quadratic Lagrange through (0, m1), (1/2, m2), (1, m3).
      return bm.m1 * (1 - 3 * xi + 2 * xi * xi) + bm.m2 * (4 * xi - 4 * xi * xi) +
             bm.m3 * (-xi + 2 * xi * xi);
    case BarMassKind::Auto:
    case BarMassKind::Nul:
      break;
  }
  return 0.0;
}

double line_mass_total(const BarMass& bm, double L) {
  switch (bm.kind) {
    case BarMassKind::Uni:
      return bm.m1 * L;
    case BarMassKind::Lin:
      return 0.5 * (bm.m1 + bm.m2) * L;
    case BarMassKind::Parabol:
      return L * (bm.m1 + 4 * bm.m2 + bm.m3) / 6.0;  // Simpson, exact for a quadratic
    case BarMassKind::Auto:
    case BarMassKind::Nul:
      break;
  }
  return 0.0;
}

Mat12 line_mass_local(const BarMass& bm, double L) {
  Mat12 m = Mat12::Zero();
  if (bm.kind == BarMassKind::Auto || bm.kind == BarMassKind::Nul) return m;
  require(L > 0, "line mass: non-positive bar length");

  // Shape functions, xi = x/L in [0, 1]:
  //   axial    u(xi) = (1 - xi) u_i + xi u_j                       (DOFs 0, 6)
  //   bending  v(xi) = N1 v_i + N2 r_i + N3 v_j + N4 r_j with
  //            N1 = 1 - 3xi^2 + 2xi^3,  N2 = l (xi - 2xi^2 + xi^3),
  //            N3 = 3xi^2 - 2xi^3,      N4 = l (-xi^2 + xi^3),
  //   l = +L for the x-y plane (v, theta_z) and l = -L for the x-z plane
  //   (w, theta_y = -dw/dx) — the same flip the consistent mass uses.
  auto accumulate = [&](const int* dof, const double* n, int cnt, double w) {
    for (int a = 0; a < cnt; ++a)
      for (int b = 0; b < cnt; ++b) m(dof[a], dof[b]) += w * n[a] * n[b];
  };

  static constexpr int kAx[2] = {0, 6};
  static constexpr int kXY[4] = {1, 5, 7, 11};
  static constexpr int kXZ[4] = {2, 4, 8, 10};

  for (int q = 0; q < kGauss; ++q) {
    double xi = kGx[q];
    double w = kGw[q] * L * line_mass_at(bm, xi);  // dx = L dxi
    if (w == 0.0) continue;
    double na[2] = {1 - xi, xi};
    accumulate(kAx, na, 2, w);
    double h1 = 1 - 3 * xi * xi + 2 * xi * xi * xi;
    double h3 = 3 * xi * xi - 2 * xi * xi * xi;
    double h2u = xi - 2 * xi * xi + xi * xi * xi;
    double h4u = -xi * xi + xi * xi * xi;
    double nxy[4] = {h1, L * h2u, h3, L * h4u};
    accumulate(kXY, nxy, 4, w);
    double nxz[4] = {h1, -L * h2u, h3, -L * h4u};
    accumulate(kXZ, nxz, 4, w);
  }
  return m;
}

Mat12 line_mass_lumped_local(const BarMass& bm, double L) {
  Mat12 m = Mat12::Zero();
  if (bm.kind == BarMassKind::Auto || bm.kind == BarMassKind::Nul) return m;
  require(L > 0, "line mass: non-positive bar length");
  double mi = 0, mj = 0;
  for (int q = 0; q < kGauss; ++q) {
    double xi = kGx[q];
    double w = kGw[q] * L * line_mass_at(bm, xi);
    mi += w * (1 - xi);
    mj += w * xi;
  }
  for (int c = 0; c < 3; ++c) {
    m(c, c) = mi;
    m(6 + c, 6 + c) = mj;
  }
  return m;
}

Mat12 line_mass_nodal(const BarElement& el, const BarMass& bm, bool lumped) {
  Mat12 m = lumped ? line_mass_lumped_local(bm, el.L()) : line_mass_local(bm, el.L());
  // Nodal-frame transform: rotation * rigid arm, exactly like BarElement's own
  // T (element.cpp build_transform). The mass line hangs at the bar's own
  // eccentricity PLUS the ST1 EXC offsets of the mass itself.
  Mat12 rot = Mat12::Zero();
  for (int blk = 0; blk < 4; ++blk) rot.block<3, 3>(3 * blk, 3 * blk) = el.R();
  Mat12 ecc = Mat12::Identity();
  Vec3 ei = el.bar().ecc[0] + bm.ecc_or;
  Vec3 ej = el.bar().ecc[1] + bm.ecc_ex;
  ecc.block<3, 3>(0, 3) = -skew(ei);  // u_end = u_node + theta x e
  ecc.block<3, 3>(6, 9) = -skew(ej);
  Mat12 t = rot * ecc;
  return t.transpose() * m * t;
}

bool bar_mass_is_null(const Model& model, Id bar) {
  for (const BarMass& bm : model.bar_masses)
    if (bm.bar == bar && bm.kind == BarMassKind::Nul) return true;
  return false;
}

bool bar_mass_is_auto(const Model& model, Id bar) {
  for (const BarMass& bm : model.bar_masses)
    if (bm.bar == bar && bm.kind == BarMassKind::Auto) return true;
  return false;
}

bool bar_mass_lumped(const Model& model, Id bar, bool lumped) {
  for (const BarMass& bm : model.bar_masses)
    if (bm.bar == bar && bm.kind == BarMassKind::Auto && bm.form != BarMass::Form::ModelWide)
      return bm.form == BarMass::Form::Diagonal;
  return lumped;
}

bool bar_mass_full_diagonal(const Model& model, Id bar) {
  for (const BarMass& bm : model.bar_masses)
    if (bm.bar == bar && bm.kind == BarMassKind::Auto) return bm.form == BarMass::Form::Diagonal;
  return false;
}

std::vector<std::pair<int, Mat6>> effective_nodal_masses(const Model& model) {
  std::vector<std::pair<int, Mat6>> out;
  for (const NodalMass& nm : model.nodal_masses) {
    int ni = model.node_index(nm.node);
    auto it = out.begin();
    for (; it != out.end(); ++it)
      if (it->first == ni) break;
    if (it == out.end()) {
      out.emplace_back(ni, nm.m);
    } else if (nm.add) {
      it->second += nm.m;  // ST1 SUP
    } else {
      it->second = nm.m;   // no SUP: replaces the earlier declaration
    }
  }
  return out;
}

namespace {

void scatter6(const DofMap& d, int ni, const Mat6& k, std::vector<Triplet>& trips) {
  int g[6];
  for (int c = 0; c < 6; ++c) g[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
  for (int a = 0; a < 6; ++a) {
    if (g[a] < 0) continue;
    for (int b = 0; b < 6; ++b) {
      if (g[b] < 0) continue;
      if (k(a, b) != 0.0) trips.emplace_back(g[a], g[b], k(a, b));
    }
  }
}

}  // namespace

void scatter_bar_line_masses(const Model& model, const Assembled& asm_,
                             std::vector<Triplet>& trips) {
  const DofMap& d = asm_.dofs;
  for (const BarMass& bm : model.bar_masses) {
    if (bm.kind == BarMassKind::Auto || bm.kind == BarMassKind::Nul) continue;
    std::size_t e = static_cast<std::size_t>(model.bar_index(bm.bar));
    const Bar& b = model.bars[e];
    require(b.kind != BarKind::Cable, "bar " + std::to_string(b.id) + ": line mass on a cable");
    // On a MASSE NUL bar the added line mass is lumped onto the ends as point
    // masses (ST1's SUP UNI on a NUL portal beam feeds the sway mode
    // exactly and adds no bending mode of the beam); elsewhere it keeps the
    // consistent form (Exemple7 and spatial modes match ST1 with it, not
    // lumped).
    Mat12 m = line_mass_nodal(asm_.elements[e], bm, bar_mass_is_null(model, bm.bar));
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    int g[12];
    for (int c = 0; c < 6; ++c) {
      g[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
      g[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
    }
    for (int a = 0; a < 12; ++a) {
      if (g[a] < 0) continue;
      for (int bb = 0; bb < 12; ++bb) {
        if (g[bb] < 0) continue;
        if (m(a, bb) != 0.0) trips.emplace_back(g[a], g[bb], m(a, bb));
      }
    }
  }
}

void scatter_nodal_masses(const Model& model, const Assembled& asm_, std::vector<Triplet>& trips) {
  for (const auto& [ni, m] : effective_nodal_masses(model)) scatter6(asm_.dofs, ni, m, trips);
}

}  // namespace keta

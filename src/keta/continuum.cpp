#include "keta/continuum.h"
#include "keta/continuum_loads.h"
#include "keta/conduction.h"
#include "keta/continuum_plastic.h"

#include <functional>

#include <Eigen/Eigenvalues>

#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>

#include "keta/assembly.h"

#include <cmath>

namespace keta {

namespace {

// a first-order dual number in up to 3 variables (exact polynomial derivatives)
struct D3 {
  double v = 0;
  double d[3] = {0, 0, 0};
};
D3 cst(double a) {
  D3 r;
  r.v = a;
  return r;
}
D3 var(double a, int k) {
  D3 r;
  r.v = a;
  r.d[k] = 1;
  return r;
}
D3 operator+(const D3& a, const D3& b) {
  D3 r;
  r.v = a.v + b.v;
  for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] + b.d[k];
  return r;
}
D3 operator-(const D3& a, const D3& b) {
  D3 r;
  r.v = a.v - b.v;
  for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] - b.d[k];
  return r;
}
D3 operator*(const D3& a, const D3& b) {
  D3 r;
  r.v = a.v * b.v;
  for (int k = 0; k < 3; ++k) r.d[k] = a.d[k] * b.v + a.v * b.d[k];
  return r;
}
D3 operator*(double s, const D3& a) { return cst(s) * a; }
D3 operator+(double s, const D3& a) { return cst(s) + a; }
D3 operator-(double s, const D3& a) { return cst(s) - a; }

std::vector<D3> shapes(int dim, int n, const double* x) {
  std::vector<D3> N;
  const D3 xi = var(x[0], 0), eta = var(x[1], 1), ze = var(dim == 3 ? x[2] : 0.0, 2);
  if (dim == 2 && (n == 3 || n == 6)) {
    const D3 L1 = 1.0 - xi - eta, L2 = xi, L3 = eta;
    if (n == 3) return {L1, L2, L3};
    return {L1 * (2.0 * L1 - cst(1)), L2 * (2.0 * L2 - cst(1)), L3 * (2.0 * L3 - cst(1)), 4.0 * L1 * L2, 4.0 * L2 * L3,
            4.0 * L3 * L1};
  }
  if (dim == 2 && (n == 4 || n == 8)) {
    const double cx[4] = {-1, 1, 1, -1}, cy[4] = {-1, -1, 1, 1};
    for (int i = 0; i < 4; ++i) {
      const D3 a = 1.0 + cx[i] * xi, b = 1.0 + cy[i] * eta;
      if (n == 4)
        N.push_back(0.25 * a * b);
      else
        N.push_back(0.25 * a * b * (cx[i] * xi + cy[i] * eta - cst(1)));
    }
    if (n == 8) {
      N.push_back(0.5 * (1.0 - xi * xi) * (1.0 - eta));
      N.push_back(0.5 * (1.0 + xi) * (1.0 - eta * eta));
      N.push_back(0.5 * (1.0 - xi * xi) * (1.0 + eta));
      N.push_back(0.5 * (1.0 - xi) * (1.0 - eta * eta));
    }
    return N;
  }
  if (dim == 3 && n == 4) return {1.0 - xi - eta - ze, xi, eta, ze};
  if (dim == 3 && (n == 8 || n == 20)) {
    const double cx[8] = {-1, 1, 1, -1, -1, 1, 1, -1}, cy[8] = {-1, -1, 1, 1, -1, -1, 1, 1},
                 cz[8] = {-1, -1, -1, -1, 1, 1, 1, 1};
    for (int i = 0; i < 8; ++i) {
      const D3 a = 1.0 + cx[i] * xi, b = 1.0 + cy[i] * eta, c = 1.0 + cz[i] * ze;
      if (n == 8)
        N.push_back(0.125 * a * b * c);
      else
        N.push_back(0.125 * a * b * c * (cx[i] * xi + cy[i] * eta + cz[i] * ze - cst(2)));
    }
    if (n == 20) {
      // mid-edges: bottom 01 12 23 30, top 45 56 67 74, verticals 04 15 26 37
      const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
      for (const auto& ed : e) {
        const double mx = (cx[ed[0]] + cx[ed[1]]) / 2, my = (cy[ed[0]] + cy[ed[1]]) / 2, mz = (cz[ed[0]] + cz[ed[1]]) / 2;
        const D3 fx = mx == 0 ? 1.0 - xi * xi : 1.0 + mx * xi;
        const D3 fy = my == 0 ? 1.0 - eta * eta : 1.0 + my * eta;
        const D3 fz = mz == 0 ? 1.0 - ze * ze : 1.0 + mz * ze;
        N.push_back(0.25 * fx * fy * fz);
      }
    }
    return N;
  }
  throw Error("continuum: unsupported element (" + std::to_string(dim) + "-D, " + std::to_string(n) + " nodes)");
}

struct GP {
  double x[3];
  double w;
};
std::vector<GP> gauss_line(int p) {
  std::vector<GP> g;
  std::vector<double> x, w;
  if (p == 1) x = {0}, w = {2};
  if (p == 2) x = {-1 / std::sqrt(3.0), 1 / std::sqrt(3.0)}, w = {1, 1};
  if (p == 3) x = {-std::sqrt(0.6), 0, std::sqrt(0.6)}, w = {5.0 / 9, 8.0 / 9, 5.0 / 9};
  if (p == 4) {
    const double a = std::sqrt(3.0 / 7 - 2.0 / 7 * std::sqrt(1.2)), b = std::sqrt(3.0 / 7 + 2.0 / 7 * std::sqrt(1.2));
    const double wa = (18 + std::sqrt(30.0)) / 36, wb = (18 - std::sqrt(30.0)) / 36;
    x = {-b, -a, a, b}, w = {wb, wa, wa, wb};
  }
  for (std::size_t i = 0; i < x.size(); ++i) g.push_back({{x[i], 0, 0}, w[i]});
  return g;
}
// a rule exact to `deg` on the reference element
std::vector<GP> rule(int dim, int n, int deg) {
  std::vector<GP> r;
  if (dim == 2 && (n == 3 || n == 6)) {
    if (deg <= 1) return {{{1.0 / 3, 1.0 / 3, 0}, 0.5}};
    if (deg <= 2) return {{{1.0 / 6, 1.0 / 6, 0}, 1.0 / 6}, {{2.0 / 3, 1.0 / 6, 0}, 1.0 / 6}, {{1.0 / 6, 2.0 / 3, 0}, 1.0 / 6}};
    // Dunavant degree 5 (7 points)
    const double a1 = 0.059715871789770, b1 = 0.470142064105115, a2 = 0.797426985353087, b2 = 0.101286507323456;
    const double w0 = 0.225 / 2, w1 = 0.132394152788506 / 2, w2 = 0.125939180544827 / 2;
    return {{{1.0 / 3, 1.0 / 3, 0}, w0}, {{b1, b1, 0}, w1}, {{a1, b1, 0}, w1}, {{b1, a1, 0}, w1},
            {{b2, b2, 0}, w2},           {{a2, b2, 0}, w2}, {{b2, a2, 0}, w2}};
  }
  if (dim == 3 && n == 4) {
    if (deg <= 1) return {{{0.25, 0.25, 0.25}, 1.0 / 6}};
    const double a = 0.5854101966249685, b = 0.1381966011250105;
    return {{{b, b, b}, 1.0 / 24}, {{a, b, b}, 1.0 / 24}, {{b, a, b}, 1.0 / 24}, {{b, b, a}, 1.0 / 24}};
  }
  const int p = std::min(4, (deg + 2) / 2);
  const std::vector<GP> l = gauss_line(p);
  if (dim == 2) {
    for (const GP& a : l)
      for (const GP& b : l) r.push_back({{a.x[0], b.x[0], 0}, a.w * b.w});
    return r;
  }
  for (const GP& a : l)
    for (const GP& b : l)
      for (const GP& c : l) r.push_back({{a.x[0], b.x[0], c.x[0]}, a.w * b.w * c.w});
  return r;
}
int order_of(int dim, int n) {  // polynomial order of the shape functions
  return (n == 3 || (dim == 3 && n == 4) || (dim == 2 && n == 4) || (dim == 3 && n == 8)) ? 1 : 2;
}
int dim_of(Continuum::Kind k) { return k == Continuum::Kind::Solid ? 3 : 2; }

MatX elasticity(Continuum::Kind k, double E, double nu) {
  if (k == Continuum::Kind::PlaneStress) {
    MatX D(3, 3);
    const double c = E / (1 - nu * nu);
    D << c, c * nu, 0, c * nu, c, 0, 0, 0, c * (1 - nu) / 2;
    return D;
  }
  const double c = E / ((1 + nu) * (1 - 2 * nu));
  if (k == Continuum::Kind::PlaneStrain) {
    MatX D(3, 3);
    D << c * (1 - nu), c * nu, 0, c * nu, c * (1 - nu), 0, 0, 0, c * (1 - 2 * nu) / 2;
    return D;
  }
  if (k == Continuum::Kind::Axisymmetric) {  // (er, ez, etheta, grz)
    MatX D = MatX::Zero(4, 4);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) D(i, j) = c * (i == j ? 1 - nu : nu);
    D(3, 3) = c * (1 - 2 * nu) / 2;
    return D;
  }
  MatX D = MatX::Zero(6, 6);  // (ex, ey, ez, gxy, gyz, gzx)
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) D(i, j) = c * (i == j ? 1 - nu : nu);
  for (int i = 3; i < 6; ++i) D(i, i) = c * (1 - 2 * nu) / 2;
  return D;
}

struct Geo {
  int dim = 2, n = 0;
  MatX X;  // n x dim
};
Geo geometry(const Model& model, const Continuum& c) {
  Geo g;
  g.dim = dim_of(c.kind);
  g.n = static_cast<int>(c.n.size());
  g.X.resize(g.n, g.dim);
  for (int i = 0; i < g.n; ++i) {
    const Vec3& p = model.node(c.n[static_cast<std::size_t>(i)]).x;
    for (int k = 0; k < g.dim; ++k) g.X(i, k) = p(k);
    if (g.dim == 2)
      require(std::fabs(p.z() - model.node(c.n[0]).x.z()) <= 1e-9 * (1 + p.norm()),
              "continuum " + std::to_string(c.id) + ": a plane / axisymmetric element lies in a plane z = const");
  }
  return g;
}

// B (strain-displacement) at a point, its Jacobian determinant and r
struct Pt {
  MatX B;
  VecX N;
  double detJ = 0, r = 0;
};
Pt point(const Continuum& c, const Geo& g, const double* xi) {
  VecX N;
  MatX dN;
  continuum_shape(g.dim, g.n, xi, N, dN);
  const MatX J = dN.transpose() * g.X;  // dim x dim, J(a, k) = d x_k / d xi_a
  Pt p;
  p.detJ = J.determinant();
  require(p.detJ > 0, "continuum " + std::to_string(c.id) + ": non-positive Jacobian (node order or distortion)");
  const MatX dX = dN * J.inverse().transpose();  // n x dim: d N_i / d x_k
  p.N = N;
  const int nd = g.dim;
  if (c.kind == Continuum::Kind::Solid) {
    p.B = MatX::Zero(6, 3 * g.n);
    for (int i = 0; i < g.n; ++i) {
      const double a = dX(i, 0), b = dX(i, 1), d = dX(i, 2);
      p.B(0, 3 * i) = a;
      p.B(1, 3 * i + 1) = b;
      p.B(2, 3 * i + 2) = d;
      p.B(3, 3 * i) = b, p.B(3, 3 * i + 1) = a;
      p.B(4, 3 * i + 1) = d, p.B(4, 3 * i + 2) = b;
      p.B(5, 3 * i) = d, p.B(5, 3 * i + 2) = a;
    }
    return p;
  }
  if (c.kind == Continuum::Kind::Axisymmetric) {
    p.r = N.dot(g.X.col(0));
    require(p.r > 0, "continuum " + std::to_string(c.id) + ": an axisymmetric element needs r = x > 0 inside it");
    p.B = MatX::Zero(4, 2 * g.n);
    for (int i = 0; i < g.n; ++i) {
      p.B(0, 2 * i) = dX(i, 0);
      p.B(1, 2 * i + 1) = dX(i, 1);
      p.B(2, 2 * i) = N(i) / p.r;
      p.B(3, 2 * i) = dX(i, 1), p.B(3, 2 * i + 1) = dX(i, 0);
    }
    return p;
  }
  (void)nd;
  p.B = MatX::Zero(3, 2 * g.n);
  for (int i = 0; i < g.n; ++i) {
    p.B(0, 2 * i) = dX(i, 0);
    p.B(1, 2 * i + 1) = dX(i, 1);
    p.B(2, 2 * i) = dX(i, 1), p.B(2, 2 * i + 1) = dX(i, 0);
  }
  return p;
}
}  // namespace

void continuum_shape(int dim, int n, const double* xi, VecX& N, MatX& dN) {
  const std::vector<D3> s = shapes(dim, n, xi);
  N.resize(n);
  dN.resize(n, dim);
  for (int i = 0; i < n; ++i) {
    N(i) = s[static_cast<std::size_t>(i)].v;
    for (int k = 0; k < dim; ++k) dN(i, k) = s[static_cast<std::size_t>(i)].d[k];
  }
}

std::vector<std::array<double, 3>> continuum_node_naturals(int dim, int n) {
  if (dim == 2 && n == 3) return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
  if (dim == 2 && n == 6) return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}};
  if (dim == 2 && n == 4) return {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}};
  if (dim == 2 && n == 8) return {{-1, -1, 0}, {1, -1, 0}, {1, 1, 0}, {-1, 1, 0}, {0, -1, 0}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}};
  if (dim == 3 && n == 4) return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  std::vector<std::array<double, 3>> c = {{-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
                                          {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1}};
  if (n == 20) {
    const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    const auto corners = c;
    for (const auto& ed : e)
      c.push_back({(corners[ed[0]][0] + corners[ed[1]][0]) / 2, (corners[ed[0]][1] + corners[ed[1]][1]) / 2,
                   (corners[ed[0]][2] + corners[ed[1]][2]) / 2});
  }
  return c;
}

ContinuumKernel continuum_kernel(const Model& model, const Continuum& c) {
  const Material& mat = model.material(c.material);
  const Geo g = geometry(model, c);
  ContinuumKernel k;
  k.dim = g.dim;
  for (Id id : c.n) k.nodes.push_back(model.node_index(id));
  const MatX D = elasticity(c.kind, mat.E, mat.nu);
  const int nd = g.dim * g.n, o = order_of(g.dim, g.n);
  const bool axi = c.kind == Continuum::Kind::Axisymmetric;
  const double tt = (c.kind == Continuum::Kind::PlaneStress || c.kind == Continuum::Kind::PlaneStrain) ? c.t : 1.0;
  k.K = MatX::Zero(nd, nd);
  k.M = MatX::Zero(nd, nd);
  // stiffness rules (see the header): simplices exact on affine elements
  // (axisymmetric one order more); Q4 / H8 2 x 2 (x 2); Q8 2 x 2 (the
  // reduced rule);
  // H20 3 x 3 x 3 (2 x 2 x 2 leaves spurious modes in a single-element
  // column)
  int sdeg = 2 * (o - 1) + (axi ? 2 : 0);
  if (g.dim == 2 && (g.n == 4 || g.n == 8)) sdeg = 2;   // 2 x 2
  if (g.dim == 3 && g.n == 8) sdeg = 2;                 // 2 x 2 x 2
  if (g.dim == 3 && g.n == 20) sdeg = 4;                // 3 x 3 x 3
  for (const GP& q : rule(g.dim, g.n, sdeg)) {
    const Pt p = point(c, g, q.x);
    const double w = q.w * p.detJ * tt * (axi ? p.r : 1.0);
    k.K += w * (p.B.transpose() * D * p.B);
  }
  if (mat.rho != 0)
    for (const GP& q : rule(g.dim, g.n, 2 * o + (axi ? 1 : 0))) {
      const Pt p = point(c, g, q.x);
      const double w = q.w * p.detJ * tt * (axi ? p.r : 1.0) * mat.rho;
      for (int i = 0; i < g.n; ++i)
        for (int j = 0; j < g.n; ++j)
          for (int a = 0; a < g.dim; ++a) k.M(g.dim * i + a, g.dim * j + a) += w * p.N(i) * p.N(j);
    }
  return k;
}

VecX continuum_stress(const Model& model, const Continuum& c, const std::vector<Vec6>& node_disp,
                      const std::array<double, 3>& xi) {
  const Material& mat = model.material(c.material);
  const Geo g = geometry(model, c);
  const Pt p = point(c, g, xi.data());
  VecX u(g.dim * g.n);
  for (int i = 0; i < g.n; ++i)
    for (int a = 0; a < g.dim; ++a)
      u(g.dim * i + a) = node_disp[static_cast<std::size_t>(model.node_index(c.n[static_cast<std::size_t>(i)]))](a);
  const VecX s = elasticity(c.kind, mat.E, mat.nu) * (p.B * u);
  if (c.kind == Continuum::Kind::PlaneStress || c.kind == Continuum::Kind::PlaneStrain) {
    VecX out(4);
    out << s(0), s(1), s(2), c.kind == Continuum::Kind::PlaneStrain ? mat.nu * (s(0) + s(1)) : 0.0;
    return out;
  }
  return s;
}

namespace {
std::vector<Load> face_load_impl(const Model& model, Id element, int face, const Vec3& traction, const Vec3& tb, bool linear,
                                 double pressure = 0, bool is_pressure = false);
}
std::vector<Load> continuum_face_load(const Model& model, Id element, int face, const Vec3& traction) {
  return face_load_impl(model, element, face, traction, traction, false);
}
std::vector<Load> continuum_edge_load_linear(const Model& model, Id element, int edge, const Vec3& ta, const Vec3& tb) {
  return face_load_impl(model, element, edge, ta, tb, true);
}
namespace {
std::vector<Load> face_load_impl(const Model& model, Id element, int face, const Vec3& traction, const Vec3& tb, bool linear,
                                 double pressure, bool is_pressure) {
  const Continuum* cp = nullptr;
  for (const Continuum& c : model.continua)
    if (c.id == element) cp = &c;
  require(cp != nullptr, "continuum face load: no element " + std::to_string(element));
  const Continuum& c = *cp;
  const int dim = dim_of(c.kind), n = static_cast<int>(c.n.size());
  std::vector<double> f(static_cast<std::size_t>(n) * 3, 0.0);
  const auto nat = continuum_node_naturals(dim, n);
  if (dim == 2) {
    const int nc = (n == 3 || n == 6) ? 3 : 4;
    require(face >= 0 && face < nc, "continuum face load: edge out of range");
    const int a = face, b = (face + 1) % nc;
    const int deg = order_of(2, n) * 2 + (c.kind == Continuum::Kind::Axisymmetric ? 1 : 0) + (linear ? 1 : 0);
    for (const GP& q : gauss_line(std::min(4, (deg + 2) / 2))) {
      const double s = (q.x[0] + 1) / 2;  // along the edge a -> b
      double xi[2];
      for (int k = 0; k < 2; ++k) xi[k] = nat[a][k] + s * (nat[b][k] - nat[a][k]);
      VecX N;
      MatX dN;
      continuum_shape(2, n, xi, N, dN);
      Eigen::Vector2d t = Eigen::Vector2d::Zero();
      double r = 0;
      for (int i = 0; i < n; ++i) {
        const Vec3& x = model.node(c.n[static_cast<std::size_t>(i)]).x;
        Eigen::Vector2d g(dN(i, 0) * (nat[b][0] - nat[a][0]) + dN(i, 1) * (nat[b][1] - nat[a][1]), 0);
        t += g(0) * x.head<2>();
        r += N(i) * x.x();
      }
      const double dl = t.norm() * q.w / 2;  // |dx/ds| ds, ds = dxi/2
      const double w = dl * (c.kind == Continuum::Kind::Axisymmetric ? r : (c.kind == Continuum::Kind::Solid ? 1.0 : c.t));
      Vec3 tr = traction + s * (tb - traction);
      if (is_pressure) {  // the LOCAL normal of the (possibly curved) edge, pushing inwards
        Vec3 nl(t(1), -t(0), 0);
        nl /= nl.norm();
        Vec3 x = Vec3::Zero(), cen = Vec3::Zero();
        for (int i = 0; i < n; ++i) {
          x += N(i) * model.node(c.n[static_cast<std::size_t>(i)]).x;
          cen += model.node(c.n[static_cast<std::size_t>(i)]).x / n;
        }
        if (nl.dot(x - cen) < 0) nl = -nl;
        tr = -pressure * nl;
      }
      for (int i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) f[static_cast<std::size_t>(3 * i + k)] += w * N(i) * tr(k);
    }
  } else {
    require(!linear, "continuum edge load: a linear traction is for 2-D elements");
    require(n == 8 || n == 20, "continuum face load: hexahedra only (H8, H20)");
    require(face >= 0 && face < 6, "continuum face load: face out of range");
    // face: fixed natural coordinate and its value
    const int fix[6] = {2, 2, 1, 0, 1, 0};
    const double val[6] = {-1, 1, -1, 1, 1, -1};
    const int a = fix[face];
    const int u = a == 0 ? 1 : 0, v = a == 2 ? 1 : 2;
    for (const GP& q1 : gauss_line(3))
      for (const GP& q2 : gauss_line(3)) {
        double xi[3];
        xi[a] = val[face];
        xi[u] = q1.x[0];
        xi[v] = q2.x[0];
        VecX N;
        MatX dN;
        continuum_shape(3, n, xi, N, dN);
        Vec3 tu = Vec3::Zero(), tv = Vec3::Zero();
        for (int i = 0; i < n; ++i) {
          const Vec3& x = model.node(c.n[static_cast<std::size_t>(i)]).x;
          tu += dN(i, u) * x;
          tv += dN(i, v) * x;
        }
        const double w = tu.cross(tv).norm() * q1.w * q2.w;
        Vec3 tr = traction;
        if (is_pressure) {  // the local normal of the (possibly curved) face
          Vec3 nl = tu.cross(tv).normalized();
          Vec3 x = Vec3::Zero(), cen = Vec3::Zero();
          for (int i = 0; i < n; ++i) {
            x += N(i) * model.node(c.n[static_cast<std::size_t>(i)]).x;
            cen += model.node(c.n[static_cast<std::size_t>(i)]).x / n;
          }
          if (nl.dot(x - cen) < 0) nl = -nl;
          tr = -pressure * nl;
        }
        for (int i = 0; i < n; ++i)
          for (int k = 0; k < 3; ++k) f[static_cast<std::size_t>(3 * i + k)] += w * N(i) * tr(k);
      }
  }
  std::vector<Load> out;
  for (int i = 0; i < n; ++i) {
    const Vec3 F(f[static_cast<std::size_t>(3 * i)], f[static_cast<std::size_t>(3 * i + 1)], f[static_cast<std::size_t>(3 * i + 2)]);
    if (F.norm() != 0) out.push_back(Load::nodal(c.n[static_cast<std::size_t>(i)], F));
  }
  return out;
}
}  // namespace

std::vector<Load> continuum_face_pressure(const Model& model, Id element, int face, double p) {
  // a uniform pressure along the LOCAL inward normal of the edge / face,
  // exact on curved (quadratic) sides too
  return face_load_impl(model, element, face, Vec3::Zero(), Vec3::Zero(), false, p, true);
}

// ---------------- J2 plasticity ----------------

namespace {
using Mat6d = Eigen::Matrix<double, 6, 6>;
using Vec6d = Eigen::Matrix<double, 6, 1>;
struct GPState {
  Vec6d ep = Vec6d::Zero();  // plastic strain (engineering shears)
  double ebar = 0;
};
// the 3-D strain (xx, yy, zz, gxy, gyz, gzx) from the element's strain vector
MatX expand(Continuum::Kind k) {
  if (k == Continuum::Kind::Solid) return MatX::Identity(6, 6);
  MatX S = MatX::Zero(6, k == Continuum::Kind::Axisymmetric ? 4 : 3);
  if (k == Continuum::Kind::Axisymmetric) {  // (er, ez, et, grz) -> (xx = rr, yy = zz, zz = tt, xy = rz)
    S(0, 0) = S(1, 1) = S(2, 2) = S(3, 3) = 1;
  } else {  // plane strain (ex, ey, gxy)
    S(0, 0) = S(1, 1) = S(3, 2) = 1;
  }
  return S;
}
// the J2 return: stress and algorithmic tangent at the total strain eps
// the scalar root of a strictly increasing g on [lo, hi] (g(lo) <= 0 <= g(hi)):
// Newton safeguarded by bisection, to round-off
template <class G>
double increasing_root(G g, double lo, double hi) {
  double x = hi;
  for (int it = 0; it < 200; ++it) {
    double gx, dg;
    g(x, gx, dg);
    if (gx == 0) return x;
    if (gx > 0) hi = x;
    else lo = x;
    double xn = dg > 0 ? x - gx / dg : 0.5 * (lo + hi);
    if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
    if (std::abs(xn - x) <= 1e-15 * std::max(std::abs(x), 1e-300) || hi - lo <= 1e-15 * std::max(std::abs(hi), 1e-300)) return xn;
    x = xn;
  }
  return x;
}
void j2_update(const Vec6d& eps, const GPState& st, double E, double nu, const J2Material* m, Vec6d& sig, Mat6d& D,
               GPState& out, double t0 = 0, double t1 = 0, double cfac = 1) {
  const double G = E / (2 * (1 + nu)), K = E / (3 * (1 - 2 * nu));
  Vec6d one;
  one << 1, 1, 1, 0, 0, 0;
  Mat6d P = Mat6d::Zero();  // deviatoric projector, engineering strain -> stress-like
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) P(i, j) = (i == j ? 1.0 : 0.0) - 1.0 / 3;
  for (int i = 3; i < 6; ++i) P(i, i) = 0.5;
  const Mat6d De = K * one * one.transpose() + 2 * G * P;
  const Vec6d ee = eps - st.ep;
  sig = De * ee;
  D = De;
  out = st;
  if (!m) return;
  const double p = (sig(0) + sig(1) + sig(2)) / 3;
  Vec6d s = sig;
  for (int i = 0; i < 3; ++i) s(i) -= p;
  const double sn = std::sqrt(s(0) * s(0) + s(1) * s(1) + s(2) * s(2) + 2 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5)));
  const double c = m->creep_A > 0 && t1 > t0 ? cfac * m->creep_A * (std::pow(t1, m->creep_m) - std::pow(t0, m->creep_m)) : 0.0;
  if (c > 0) {
    // creep (and plasticity): the von Mises stress q = q_tr - 3G (d ep + d ec)
    // with d ec = c q^n, and q = R0 + H d ep while yielding
    const double G3 = 3 * G, qtr = std::sqrt(1.5) * sn, n_ = m->creep_n, R0 = m->sigma_y + m->H * st.ebar;
    if (qtr <= 0) return;
    double q = increasing_root(
        [&](double x, double& g, double& dg) {
          g = x - qtr + G3 * c * std::pow(x, n_);
          dg = 1 + G3 * c * n_ * std::pow(x, n_ - 1);
        },
        0.0, qtr);
    double ct = 1 / (1 + G3 * c * n_ * std::pow(q, n_ - 1)), dep = 0;
    if (m->sigma_y > 0 && q > R0) {
      if (m->H > 0) {
        q = increasing_root(
            [&](double x, double& g, double& dg) {
              g = x - qtr + G3 * (x - R0) / m->H + G3 * c * std::pow(x, n_);
              dg = 1 + G3 / m->H + G3 * c * n_ * std::pow(x, n_ - 1);
            },
            R0, q);
        dep = (q - R0) / m->H;
        ct = 1 / (1 + G3 / m->H + G3 * c * n_ * std::pow(q, n_ - 1));
      } else {
        q = R0;
        dep = (qtr - q) / G3 - c * std::pow(q, n_);
        ct = 0;
      }
    }
    const Vec6d n = s / sn;
    const double dg = std::sqrt(1.5) * (qtr - q) / G3;
    sig -= 2 * G * dg * n;
    Vec6d nstrain = n;
    for (int i = 3; i < 6; ++i) nstrain(i) *= 2;
    out.ep = st.ep + dg * nstrain;
    out.ebar = st.ebar + dep;
    const double beta = q / qtr;
    D = K * one * one.transpose() + 2 * G * beta * P + 2 * G * (ct - beta) * n * n.transpose();
    return;
  }
  if (!(m->sigma_y > 0)) return;  // viscous only, no time elapsed
  const double R = std::sqrt(2.0 / 3) * (m->sigma_y + m->H * st.ebar);
  const double f = sn - R;
  if (f <= 0) return;
  const Vec6d n = s / sn;
  const double dg = f / (2 * G + 2.0 / 3 * m->H);
  sig -= 2 * G * dg * n;  // the shear components of n are stress-like: 2G n gives the stress correction
  Vec6d nstrain = n;
  for (int i = 3; i < 6; ++i) nstrain(i) *= 2;
  out.ep = st.ep + dg * nstrain;
  out.ebar = st.ebar + std::sqrt(2.0 / 3) * dg;
  const double beta = 1 - 2 * G * dg / sn, gbar = 1 / (1 + m->H / (3 * G)) - (1 - beta);
  D = K * one * one.transpose() + 2 * G * beta * P - 2 * G * gbar * n * n.transpose();
}
}  // namespace

// ---------------- scalar conduction (keta/conduction.h) ----------------

namespace {
// the Gauss points of an element side: shape values of the WHOLE element
// and the weight (length / area times thickness, or r for axisymmetry)
void face_points(const Model& model, const Continuum& c, int face, const std::function<void(const VecX&, double)>& cb) {
  const int dim = dim_of(c.kind), n = static_cast<int>(c.n.size());
  const auto nat = continuum_node_naturals(dim, n);
  if (dim == 2) {
    const int nc = (n == 3 || n == 6) ? 3 : 4;
    require(face >= 0 && face < nc, "conduction: edge out of range");
    const int a = face, b = (face + 1) % nc;
    for (const GP& q : gauss_line(4)) {
      const double s = (q.x[0] + 1) / 2;
      double xi[2];
      for (int k = 0; k < 2; ++k) xi[k] = nat[a][k] + s * (nat[b][k] - nat[a][k]);
      VecX N;
      MatX dN;
      continuum_shape(2, n, xi, N, dN);
      Eigen::Vector2d t = Eigen::Vector2d::Zero();
      double r = 0;
      for (int i = 0; i < n; ++i) {
        const Vec3& x = model.node(c.n[static_cast<std::size_t>(i)]).x;
        t += (dN(i, 0) * (nat[b][0] - nat[a][0]) + dN(i, 1) * (nat[b][1] - nat[a][1])) * x.head<2>();
        r += N(i) * x.x();
      }
      const double dl = t.norm() * q.w / 2;
      cb(N, dl * (c.kind == Continuum::Kind::Axisymmetric ? r : c.t));
    }
    return;
  }
  require(n == 8 || n == 20, "conduction: faces of hexahedra only (H8, H20)");
  require(face >= 0 && face < 6, "conduction: face out of range");
  const int fix[6] = {2, 2, 1, 0, 1, 0};
  const double val[6] = {-1, 1, -1, 1, 1, -1};
  const int a = fix[face], u = a == 0 ? 1 : 0, v = a == 2 ? 1 : 2;
  for (const GP& q1 : gauss_line(4))
    for (const GP& q2 : gauss_line(4)) {
      double xi[3];
      xi[a] = val[face];
      xi[u] = q1.x[0];
      xi[v] = q2.x[0];
      VecX N;
      MatX dN;
      continuum_shape(3, n, xi, N, dN);
      Vec3 tu = Vec3::Zero(), tv = Vec3::Zero();
      for (int i = 0; i < n; ++i) {
        const Vec3& x = model.node(c.n[static_cast<std::size_t>(i)]).x;
        tu += dN(i, u) * x;
        tv += dN(i, v) * x;
      }
      cb(N, tu.cross(tv).norm() * q1.w * q2.w);
    }
}
const Continuum& continuum_by_id(const Model& model, Id element) {
  for (const Continuum& c : model.continua)
    if (c.id == element) return c;
  throw Error("continuum: no element " + std::to_string(element));
}
}  // namespace

ConductionResult solve_conduction(const Model& model, const ConductionProblem& p, const ConductionOptions& opt) {
  require(!model.continua.empty(), "conduction: no continuum element");
  // unknowns: the nodes of the elements (and of the links)
  std::map<Id, int> idx;
  for (const Continuum& c : model.continua)
    for (Id id : c.n) idx.emplace(id, 0);
  for (const auto& l : p.links) {
    (void)model.node(l.a);
    idx.emplace(l.a, 0);
    if (l.b != 0) idx.emplace(l.b, 0);
  }
  std::vector<Id> ids;
  for (auto& [id, k] : idx) {
    k = static_cast<int>(ids.size());
    ids.push_back(id);
  }
  const int n = static_cast<int>(ids.size());
  for (const auto& [id, v] : p.fixed) require(idx.count(id), "conduction: fixed node " + std::to_string(id) + " is not on an element");
  // element data at the Gauss points
  struct GPt {
    VecX N;
    MatX dX;
    double w = 0;
  };
  struct El {
    std::vector<int> dof;
    std::vector<GPt> gk, gm;  // conduction and capacity rules
    const ConductionMaterial* m = nullptr;
  };
  std::vector<El> els;
  for (const Continuum& c : model.continua) {
    El e;
    auto it = p.materials.find(c.material);
    require(it != p.materials.end(), "conduction: no conduction data for material " + std::to_string(c.material));
    e.m = &it->second;
    require(e.m->kx > 0 && e.m->ky > 0 && e.m->kz > 0 && e.m->capacity >= 0, "conduction: k > 0, capacity >= 0");
    for (Id id : c.n) e.dof.push_back(idx.at(id));
    const Geo g = geometry(model, c);
    const int o = order_of(g.dim, g.n);
    const bool axi = c.kind == Continuum::Kind::Axisymmetric;
    const double tt = g.dim == 2 && !axi ? c.t : 1.0;
    auto fill = [&](int deg, std::vector<GPt>& out) {
      for (const GP& q : rule(g.dim, g.n, deg)) {
        const Pt pt = point(c, g, q.x);
        VecX N;
        MatX dN;
        continuum_shape(g.dim, g.n, q.x, N, dN);
        const MatX J = dN.transpose() * g.X;
        out.push_back({N, dN * J.inverse().transpose(), q.w * pt.detJ * tt * (axi ? pt.r : 1.0)});
      }
    };
    fill(2 * o + (axi ? 1 : 0) + 1, e.gk);  // full integration (k linear in u)
    fill(2 * o + (axi ? 1 : 0), e.gm);
    els.push_back(std::move(e));
  }
  // faces: precomputed (N_i, w) per point
  struct FacePts {
    std::vector<int> dof;
    std::vector<std::pair<VecX, double>> pts;
    double q, h, ue;
  };
  std::vector<FacePts> faces;
  for (const auto& f : p.faces) {
    const Continuum& c = continuum_by_id(model, f.element);
    require(f.h >= 0, "conduction: convection h >= 0");
    FacePts fp;
    for (Id id : c.n) fp.dof.push_back(idx.at(id));
    face_points(model, c, f.face, [&](const VecX& N, double w) { fp.pts.emplace_back(N, w); });
    fp.q = f.q;
    fp.h = f.h;
    fp.ue = f.u_env;
    faces.push_back(std::move(fp));
  }
  std::vector<char> is_fixed(static_cast<std::size_t>(n), 0);
  VecX ufix = VecX::Zero(n);
  for (const auto& [id, v] : p.fixed) {
    is_fixed[static_cast<std::size_t>(idx.at(id))] = 1;
    ufix(idx.at(id)) = v;
  }
  // the residual R(u) = internal - external (all nodes) and its tangent
  auto assemble = [&](const VecX& u, VecX& R, SpMat* J) {
    R = VecX::Zero(n);
    std::vector<Triplet> t;
    for (const El& e : els) {
      const int ne = static_cast<int>(e.dof.size()), dim = static_cast<int>(e.gk[0].dX.cols());
      VecX ue(ne);
      for (int i = 0; i < ne; ++i) ue(i) = u(e.dof[static_cast<std::size_t>(i)]);
      VecX re = VecX::Zero(ne);
      MatX ke = MatX::Zero(ne, ne);
      for (const GPt& g : e.gk) {
        const double ug = g.N.dot(ue);
        const VecX grad = g.dX.transpose() * ue;  // dim
        VecX kd(dim);
        kd(0) = e.m->kx + e.m->k1 * ug;
        kd(1) = e.m->ky + e.m->k1 * ug;
        if (dim == 3) kd(2) = e.m->kz + e.m->k1 * ug;
        const VecX flux = kd.cwiseProduct(grad);
        re += g.w * g.dX * flux;
        if (J) {
          ke += g.w * g.dX * kd.asDiagonal() * g.dX.transpose();
          if (e.m->k1 != 0) ke += g.w * e.m->k1 * (g.dX * grad) * g.N.transpose();
        }
      }
      for (int i = 0; i < ne; ++i) {
        R(e.dof[static_cast<std::size_t>(i)]) += re(i);
        if (J)
          for (int j = 0; j < ne; ++j) t.emplace_back(e.dof[static_cast<std::size_t>(i)], e.dof[static_cast<std::size_t>(j)], ke(i, j));
      }
    }
    for (const FacePts& f : faces)
      for (const auto& [N, w] : f.pts) {
        double uf = 0;
        for (std::size_t i = 0; i < f.dof.size(); ++i) uf += N(static_cast<Eigen::Index>(i)) * u(f.dof[i]);
        for (std::size_t i = 0; i < f.dof.size(); ++i) {
          const double Ni = N(static_cast<Eigen::Index>(i));
          R(f.dof[i]) += w * Ni * (f.h * (uf - f.ue) - f.q);
          if (J && f.h != 0)
            for (std::size_t j = 0; j < f.dof.size(); ++j) t.emplace_back(f.dof[i], f.dof[j], w * f.h * Ni * N(static_cast<Eigen::Index>(j)));
        }
      }
    for (const auto& [id, s] : p.sources) R(idx.at(id)) -= s;
    for (const auto& l : p.links) {
      const int a = idx.at(l.a);
      const double ub = l.b != 0 ? u(idx.at(l.b)) : l.u_env;
      const double q = l.k * (u(a) - ub);
      R(a) += q;
      if (J) t.emplace_back(a, a, l.k);
      if (l.b != 0) {
        const int b = idx.at(l.b);
        R(b) -= q;
        if (J) {
          t.emplace_back(a, b, -l.k);
          t.emplace_back(b, a, -l.k);
          t.emplace_back(b, b, l.k);
        }
      }
    }
    if (J) {
      J->resize(n, n);
      J->setFromTriplets(t.begin(), t.end());
    }
  };
  // the consistent capacity matrix
  SpMat C(n, n);
  {
    std::vector<Triplet> t;
    for (const El& e : els) {
      if (e.m->capacity == 0) continue;
      for (const GPt& g : e.gm)
        for (std::size_t i = 0; i < e.dof.size(); ++i)
          for (std::size_t j = 0; j < e.dof.size(); ++j)
            t.emplace_back(e.dof[i], e.dof[j], g.w * e.m->capacity * g.N(static_cast<Eigen::Index>(i)) * g.N(static_cast<Eigen::Index>(j)));
    }
    C.setFromTriplets(t.begin(), t.end());
  }
  // free-node map
  std::vector<int> fr(static_cast<std::size_t>(n), -1);
  int nf = 0;
  for (int i = 0; i < n; ++i)
    if (!is_fixed[static_cast<std::size_t>(i)]) fr[static_cast<std::size_t>(i)] = nf++;
  ConductionResult out;
  // Newton on the free nodes for R(u) + C (u - u_old) / dt = 0 (dt = 0: steady)
  auto newton = [&](VecX& u, const VecX& u_old, double dt) {
    double scale = 0;
    for (int it = 0;; ++it) {
      require(it < opt.max_iter, "conduction: no convergence (nonlinear conductivity: reduce the steps)");
      VecX R;
      SpMat J;
      assemble(u, R, &J);
      if (dt > 0) {
        R += C * (u - u_old) / dt;
        J += C / dt;
      }
      VecX Rf(nf);
      for (int i = 0; i < n; ++i)
        if (fr[static_cast<std::size_t>(i)] >= 0) Rf(fr[static_cast<std::size_t>(i)]) = R(i);
      if (it == 0) {
        VecX R0;
        assemble(VecX::Zero(n), R0, nullptr);
        scale = std::max({R.cwiseAbs().maxCoeff(), R0.cwiseAbs().maxCoeff(), 1e-300});
      }
      if (nf == 0 || Rf.cwiseAbs().maxCoeff() <= opt.tol * scale) return;
      std::vector<Triplet> t;
      for (int k = 0; k < J.outerSize(); ++k)
        for (SpMat::InnerIterator i2(J, k); i2; ++i2) {
          const int a = fr[static_cast<std::size_t>(i2.row())], b = fr[static_cast<std::size_t>(i2.col())];
          if (a >= 0 && b >= 0) t.emplace_back(a, b, i2.value());
        }
      SpMat Jf(nf, nf);
      Jf.setFromTriplets(t.begin(), t.end());
      Jf.makeCompressed();
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(Jf);
      require(lu.info() == Eigen::Success, "conduction: singular system (a region without any fixed value or exchange?)");
      const VecX du = lu.solve(Rf);
      require(du.allFinite(), "conduction: non-finite Newton step");
      for (int i = 0; i < n; ++i)
        if (fr[static_cast<std::size_t>(i)] >= 0) u(i) -= du(fr[static_cast<std::size_t>(i)]);
      ++out.iterations;
      if (du.cwiseAbs().maxCoeff() <= 1e-15 * std::max(1.0, u.cwiseAbs().maxCoeff())) return;
    }
  };
  auto record = [&](const VecX& u, const VecX& u_old, double dt) {
    std::map<Id, double> uu, rr;
    for (int i = 0; i < n; ++i) uu[ids[static_cast<std::size_t>(i)]] = u(i);
    VecX R;
    assemble(u, R, nullptr);
    if (dt > 0) R += C * (u - u_old) / dt;
    for (int i = 0; i < n; ++i)
      if (is_fixed[static_cast<std::size_t>(i)]) rr[ids[static_cast<std::size_t>(i)]] = R(i);
    out.u.push_back(uu);
    out.reactions.push_back(rr);
  };
  VecX u = VecX::Constant(n, p.u_initial);
  for (const auto& [id, v] : p.initial) u(idx.at(id)) = v;
  for (int i = 0; i < n; ++i)
    if (is_fixed[static_cast<std::size_t>(i)]) u(i) = ufix(i);
  if (opt.times.empty()) {
    newton(u, u, 0);
    record(u, u, 0);
    return out;
  }
  require(opt.substeps >= 1, "conduction: substeps >= 1");
  double t0 = 0;
  for (double t1 : opt.times) {
    require(t1 > t0, "conduction: increasing output times");
    const double dt = (t1 - t0) / opt.substeps;
    VecX uo = u;
    for (int k = 0; k < opt.substeps; ++k) {
      uo = u;
      newton(u, uo, dt);
    }
    record(u, uo, dt);
    t0 = t1;
  }
  return out;
}

std::vector<Load> continuum_thermal_load(const Model& model, Id element, const std::map<Id, double>& T, double T0) {
  const Continuum& c = continuum_by_id(model, element);
  const Material& mat = model.material(c.material);
  require(mat.alpha != 0, "continuum thermal load: the material needs alpha");
  const Geo g = geometry(model, c);
  const bool axi = c.kind == Continuum::Kind::Axisymmetric;
  const double tt = (c.kind == Continuum::Kind::PlaneStress || c.kind == Continuum::Kind::PlaneStrain) ? c.t : 1.0;
  const MatX D = elasticity(c.kind, mat.E, mat.nu);
  VecX Te(g.n);
  for (int i = 0; i < g.n; ++i) {
    auto it = T.find(c.n[static_cast<std::size_t>(i)]);
    require(it != T.end(), "continuum thermal load: no temperature at node " + std::to_string(c.n[static_cast<std::size_t>(i)]));
    Te(i) = it->second - T0;
  }
  const int o = order_of(g.dim, g.n);
  VecX f = VecX::Zero(g.dim * g.n);
  for (const GP& q : rule(g.dim, g.n, 2 * o + (axi ? 1 : 0))) {
    const Pt pt = point(c, g, q.x);
    const double dT = pt.N.dot(Te), w = q.w * pt.detJ * tt * (axi ? pt.r : 1.0);
    VecX e0 = VecX::Zero(D.rows());
    const double a = mat.alpha * dT;
    if (c.kind == Continuum::Kind::PlaneStress) e0 << a, a, 0;
    else if (c.kind == Continuum::Kind::PlaneStrain) e0 << (1 + mat.nu) * a, (1 + mat.nu) * a, 0;
    else if (axi) e0 << a, a, a, 0;
    else e0 << a, a, a, 0, 0, 0;
    f += w * pt.B.transpose() * (D * e0);
  }
  std::vector<Load> out;
  for (int i = 0; i < g.n; ++i) {
    Vec3 F = Vec3::Zero();
    for (int k = 0; k < g.dim; ++k) F(k) = f(g.dim * i + k);
    if (F.norm() != 0) out.push_back(Load::nodal(c.n[static_cast<std::size_t>(i)], F));
  }
  return out;
}

// ---------------- frictional plasticity (Drucker-Prager, Mohr-Coulomb) ----------------

namespace {
// Voigt (xx, yy, zz, xy, yz, zx) of a symmetric 3x3 tensor (stress-like)
Vec6 voigt(const Eigen::Matrix3d& T) {
  Vec6 v;
  v << T(0, 0), T(1, 1), T(2, 2), T(0, 1), T(1, 2), T(2, 0);
  return v;
}
// the symmetric tensor of an engineering-shear strain vector
Eigen::Matrix3d strain_tensor(const Vec6& e) {
  Eigen::Matrix3d T;
  T << e(0), e(3) / 2, e(5) / 2, e(3) / 2, e(1), e(4) / 2, e(5) / 2, e(4) / 2, e(2);
  return T;
}
Vec6 strain_voigt(const Eigen::Matrix3d& T) {
  Vec6 v;
  v << T(0, 0), T(1, 1), T(2, 2), 2 * T(0, 1), 2 * T(1, 2), 2 * T(2, 0);
  return v;
}
Mat6 iso_D(double K, double G) {
  Vec6 one;
  one << 1, 1, 1, 0, 0, 0;
  Mat6 P = Mat6::Zero();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) P(i, j) = (i == j ? 1.0 : 0.0) - 1.0 / 3;
  for (int i = 3; i < 6; ++i) P(i, i) = 0.5;
  return K * one * one.transpose() + 2 * G * P;
}
void dp_coefficients(FrictionalMaterial::DpFit fit, double ang, double& eta, double& xi) {
  const double s = std::sin(ang), c = std::cos(ang), r3 = std::sqrt(3.0);
  if (fit == FrictionalMaterial::DpFit::Outer) {
    eta = 6 * s / (r3 * (3 - s));
    xi = 6 * c / (r3 * (3 - s));
  } else if (fit == FrictionalMaterial::DpFit::Inner) {
    eta = 6 * s / (r3 * (3 + s));
    xi = 6 * c / (r3 * (3 + s));
  } else {
    const double t = std::tan(ang), q = std::sqrt(9 + 12 * t * t);
    eta = 3 * t / q;
    xi = 3 / q;
  }
}
}  // namespace

FrictionalPoint frictional_update(const FrictionalMaterial& m, double E, double nu, const Vec6& eps, const Vec6& ep,
                                  double ebar) {
  const double G = E / (2 * (1 + nu)), K = E / (3 * (1 - 2 * nu));
  FrictionalPoint out;
  out.ep = ep;
  out.ebar = ebar;
  const Mat6 De = iso_D(K, G);
  const Vec6 ee = eps - ep;
  const double cn = m.c + m.H * ebar;
  Vec6 one;
  one << 1, 1, 1, 0, 0, 0;
  if (m.kind == FrictionalMaterial::Kind::DruckerPrager) {
    double eta, xi, etab, xib;
    dp_coefficients(m.fit, m.phi, eta, xi);
    dp_coefficients(m.fit, m.psi, etab, xib);
    const Vec6 st = De * ee;
    const double p = (st(0) + st(1) + st(2)) / 3;
    Vec6 s = st;
    for (int i = 0; i < 3; ++i) s(i) -= p;
    const double sn = std::sqrt(s(0) * s(0) + s(1) * s(1) + s(2) * s(2) + 2 * (s(3) * s(3) + s(4) * s(4) + s(5) * s(5)));
    const double J2s = sn / std::sqrt(2.0);  // sqrt(J2)
    const double f = J2s + eta * p - xi * cn;
    out.sigma = st;
    out.D = De;
    if (f <= 1e-14 * std::max({std::fabs(p), J2s, xi * cn, 1e-300})) return out;
    const double A = 1 / (G + K * eta * etab + xi * xi * m.H);
    const double dg = f * A;
    if (J2s - G * dg >= 0 && sn > 0) {  // smooth cone
      const Vec6 n = s / sn;
      out.sigma = (1 - G * dg / J2s) * s + (p - K * etab * dg) * one;
      out.ebar = ebar + xi * dg;
      const double beta = G * dg / J2s;
      Mat6 P = Mat6::Zero();
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) P(i, j) = (i == j ? 1.0 : 0.0) - 1.0 / 3;
      for (int i = 3; i < 6; ++i) P(i, i) = 0.5;
      out.D = 2 * G * (1 - beta) * P + 2 * G * (beta - G * A) * n * n.transpose() -
              std::sqrt(2.0) * G * A * K * (eta * n * one.transpose() + etab * one * n.transpose()) +
              K * (1 - K * eta * etab * A) * one * one.transpose();
      out.mode = 1;
    } else {  // apex: s = 0, eta p = xi c(ebar + xi / eta_bar dev)
      require(etab > 0 && eta > 0, "Drucker-Prager: the apex is reached with no dilatancy (psi = 0) or phi = 0");
      const double al = xi / etab, be = xi / eta;
      const double dev = (p - be * cn) / (K + al * be * m.H);
      out.sigma = (p - K * dev) * one;
      out.ebar = ebar + al * dev;
      out.D = K * (1 - K / (K + al * be * m.H)) * one * one.transpose();
      out.mode = 3;
    }
    Vec6 sd = out.sigma;
    const double pn = (sd(0) + sd(1) + sd(2)) / 3;
    for (int i = 0; i < 3; ++i) sd(i) -= pn;
    Vec6 een = sd / (2 * G);
    for (int i = 3; i < 6; ++i) een(i) *= 2;  // engineering shear
    for (int i = 0; i < 3; ++i) een(i) += pn / (3 * K);
    out.ep = eps - een;
    return out;
  }
  // ---- Mohr-Coulomb in principal space ----
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(strain_tensor(ee));
  Eigen::Vector3d e;
  Eigen::Matrix3d Q;  // columns: principal directions, descending strains
  for (int i = 0; i < 3; ++i) {
    e(i) = es.eigenvalues()(2 - i);
    Q.col(i) = es.eigenvectors().col(2 - i);
  }
  Eigen::Matrix3d E3 = Eigen::Matrix3d::Constant(K - 2 * G / 3);
  E3.diagonal().array() += 2 * G;
  const Eigen::Vector3d str = E3 * e;
  const double sf = std::sin(m.phi), cf = std::cos(m.phi), sp = std::sin(m.psi);
  auto plane = [&](int k, Eigen::Vector3d& mv, Eigen::Vector3d& nv) {
    // 0 main (1, 3), 1 right edge (1, 2), 2 left edge (2, 3)
    const int a = k == 2 ? 1 : 0, b = k == 1 ? 1 : 2;
    mv.setZero();
    nv.setZero();
    mv(a) = 1 + sf;
    mv(b) = -(1 - sf);
    nv(a) = 1 + sp;
    nv(b) = -(1 - sp);
  };
  Eigen::Vector3d mm, nn;
  plane(0, mm, nn);
  const double scale = std::max({std::fabs(str(0)), std::fabs(str(2)), cn, 1e-300});
  const double f0 = mm.dot(str) - 2 * cn * cf;
  Eigen::Vector3d s = str;
  Eigen::Matrix3d Dp = E3;
  int mode = 0;
  if (f0 > 1e-14 * scale) {
    auto try_planes = [&](const std::vector<int>& ks) {
      const int n = static_cast<int>(ks.size());
      Eigen::MatrixXd M(3, n), N(3, n), A(n, n);
      Eigen::VectorXd f(n);
      for (int j = 0; j < n; ++j) {
        Eigen::Vector3d a, b;
        plane(ks[static_cast<std::size_t>(j)], a, b);
        M.col(j) = a;
        N.col(j) = b;
        f(j) = a.dot(str) - 2 * cn * cf;
      }
      // dc/d(dgamma) = 2 cos phi H: each plane's -2 cos phi dc
      A = M.transpose() * E3 * N + Eigen::MatrixXd::Constant(n, n, 4 * cf * cf * m.H);
      const Eigen::VectorXd dg = A.lu().solve(f);
      for (int j = 0; j < n; ++j)
        if (dg(j) < -1e-12 * std::fabs(dg.maxCoeff()) - 1e-300) return false;
      const Eigen::Vector3d sn = str - E3 * N * dg;
      const double tol = 1e-10 * scale;
      if (sn(0) < sn(1) - tol || sn(1) < sn(2) - tol) return false;
      s = sn;
      Dp = E3 - E3 * N * A.lu().solve(M.transpose() * E3);
      out.ebar = ebar + 2 * cf * dg.sum();
      return true;
    };
    if (try_planes({0})) {
      mode = 1;
    } else if (try_planes({0, 1}) || try_planes({0, 2})) {
      mode = 2;
    } else {
      require(sf > 0 && sp > 0, "Mohr-Coulomb: the apex is reached with phi or psi = 0");
      const double ptr = str.sum() / 3, cot = cf / sf, al = cf / sp;
      const double dev = (ptr - cn * cot) / (K + m.H * al * cot);
      s = Eigen::Vector3d::Constant(ptr - K * dev);
      out.ebar = ebar + al * dev;
      Dp = Eigen::Matrix3d::Constant(K * (1 - K / (K + m.H * al * cot)));
      mode = 3;
    }
  }
  out.mode = mode;
  // back to the global frame
  Eigen::Matrix3d S = Q * s.asDiagonal() * Q.transpose();
  out.sigma = voigt(S);
  const Eigen::Vector3d een = E3.inverse() * s;
  out.ep = eps - strain_voigt(Q * een.asDiagonal() * Q.transpose());
  if (mode == 0) {
    out.D = De;
    return out;
  }
  std::array<Vec6, 3> vi;
  for (int i = 0; i < 3; ++i) vi[static_cast<std::size_t>(i)] = voigt(Q.col(i) * Q.col(i).transpose());
  out.D.setZero();
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) out.D += Dp(i, j) * vi[static_cast<std::size_t>(i)] * vi[static_cast<std::size_t>(j)].transpose();
  const double emax = e.cwiseAbs().maxCoeff();
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j) {
      const double de = e(i) - e(j);
      const double cij = std::fabs(de) > 1e-10 * std::max(emax, 1e-300) ? (s(i) - s(j)) / de : Dp(i, i) - Dp(i, j);
      const Vec6 v = voigt(0.5 * (Q.col(i) * Q.col(j).transpose() + Q.col(j) * Q.col(i).transpose()));
      out.D += 2 * cij * v * v.transpose();
    }
  return out;
}

ContinuumPlasticResult solve_continuum_plastic(const Model& model, Id case_id, const ContinuumPlasticOptions& opt) {
  model.validate();
  require(model.bars.empty() && model.shells.empty() && model.links.empty(),
          "continuum plasticity: a model of continuum elements only (bars, shells, links: not yet)");
  require(!model.continua.empty(), "continuum plasticity: no continuum element");
  require(!opt.lambdas.empty() && opt.substeps >= 1, "continuum plasticity: load factors and substeps >= 1");
  for (const auto& [id, m] : opt.plastic) {
    (void)model.material(id);
    require(m.sigma_y >= 0 && m.H >= 0 && (m.sigma_y > 0 || m.creep_A > 0),
            "continuum plasticity: sigma_y > 0 (or a creep law), H >= 0");
    require(m.creep_A >= 0 && m.creep_n > 0 && m.creep_m > 0 && m.creep_Q >= 0, "continuum creep: A >= 0, n > 0, m > 0, Q >= 0");
    if (m.creep_A > 0)
      require(opt.times.size() == opt.lambdas.size(), "continuum creep: one time per load factor (options.times)");
  }
  require(opt.times.empty() || opt.times.size() == opt.lambdas.size(), "continuum plasticity: one time per load factor");
  for (std::size_t k = 0; k < opt.times.size(); ++k)
    require(opt.times[k] >= (k ? opt.times[k - 1] : 0.0), "continuum plasticity: times start at >= 0 and never decrease");
  Assembled as = assemble_elements(model, true);
  require(!as.red.active, "continuum plasticity: rigid links / constraints are not supported");
  const DofMap& d = as.dofs;
  const CaseLoads cl = lower_case(model, as, model.load_case(case_id));
  // prescribed support displacements: scaled by the load factor
  std::vector<double> presc(static_cast<std::size_t>(kDofPerNode * d.n_nodes), 0.0);
  for (const auto& [full, val] : cl.prescribed) presc[static_cast<std::size_t>(full)] = val;
  const VecX f = cl.rhs;
  for (const auto& [id, fm] : opt.frictional) {
    (void)model.material(id);
    require(!opt.plastic.count(id), "continuum plasticity: a material is von Mises or frictional, not both");
    require(fm.c >= 0 && fm.phi >= 0 && fm.phi < kPi / 2 && fm.psi >= 0 && fm.psi <= fm.phi && fm.H >= 0 &&
                (fm.c > 0 || fm.phi > 0),
            "continuum plasticity: frictional law needs c >= 0, 0 <= psi <= phi < 90 deg, H >= 0");
  }
  const int nf = d.n_free;
  // Gauss points of every element
  struct GPt {
    MatX B;   // element strain-displacement
    MatX dX;  // n x dim: d N_i / d X_k (total Lagrangian)
    double w = 0;
    double cfac = 1;  // thermally activated creep: exp(-Q / T) at the point
  };
  struct El {
    const Continuum* c = nullptr;
    std::vector<int> dof;  // free index per element DOF (-1 blocked)
    std::vector<GPt> gp;
    MatX S;
    double E = 0, nu = 0;
    const J2Material* m = nullptr;
    const FrictionalMaterial* fm = nullptr;
    std::vector<int> full;  // full DOF index per element DOF (prescribed values)
  };
  bool unsym = false;
  std::vector<El> els;
  for (const Continuum& c : model.continua) {
    El e;
    e.c = &c;
    const Material& mat = model.material(c.material);
    e.E = mat.E;
    e.nu = mat.nu;
    auto it = opt.plastic.find(c.material);
    if (it != opt.plastic.end()) {
      e.m = &it->second;
      require(c.kind != Continuum::Kind::PlaneStress, "continuum plasticity: plane stress is not supported for a plastic material");
    }
    if (auto fi = opt.frictional.find(c.material); fi != opt.frictional.end()) {
      e.fm = &fi->second;
      unsym = unsym || fi->second.psi != fi->second.phi;
      require(c.kind != Continuum::Kind::PlaneStress, "continuum plasticity: plane stress is not supported for a frictional material");
      require(!opt.large_displacements, "continuum plasticity: frictional laws in small displacements only");
    }
    require(c.kind != Continuum::Kind::PlaneStress || !e.m, "continuum plasticity: plane stress");
    const Geo g = geometry(model, c);
    for (Id id : c.n)
      for (int a = 0; a < g.dim; ++a) {
        e.dof.push_back(d.free_of(model.node_index(id), a));
        e.full.push_back(d.full(model.node_index(id), a));
      }
    const int o = order_of(g.dim, g.n);
    const bool axi = c.kind == Continuum::Kind::Axisymmetric;
    int sdeg = 2 * (o - 1) + (axi ? 2 : 0);
    if (g.dim == 2 && (g.n == 4 || g.n == 8)) sdeg = 2;
    if (g.dim == 3 && g.n == 8) sdeg = 2;
    if (g.dim == 3 && g.n == 20) sdeg = 4;
    const double tt = c.kind == Continuum::Kind::PlaneStrain ? c.t : 1.0;
    const double tps = c.kind == Continuum::Kind::PlaneStress ? c.t : tt;
    for (const GP& q : rule(g.dim, g.n, sdeg)) {
      const Pt p = point(c, g, q.x);
      VecX N;
      MatX dN;
      continuum_shape(g.dim, g.n, q.x, N, dN);
      const MatX J = dN.transpose() * g.X;
      double cfac = 1;
      if (e.m && e.m->creep_Q > 0) {
        double T = 0;
        for (int i = 0; i < g.n; ++i) {
          const auto ti = opt.temperature.find(c.n[static_cast<std::size_t>(i)]);
          require(ti != opt.temperature.end(), "continuum creep: creep_Q needs the temperature of every node of the element");
          T += N(i) * ti->second;
        }
        require(T > 0, "continuum creep: the temperature must be positive (absolute scale)");
        cfac = std::exp(-e.m->creep_Q / T);
      }
      e.gp.push_back({p.B, dN * J.inverse().transpose(), q.w * p.detJ * tps * (axi ? p.r : 1.0), cfac});
    }
    if (opt.large_displacements)
      require(!e.m && c.kind != Continuum::Kind::Axisymmetric,
              "continuum large displacements: elastic plane stress, plane strain or solid elements only");
    e.S = expand(c.kind);
    if (c.kind == Continuum::Kind::PlaneStress) {
      // elastic plane stress: the 3-D law is not used — handled by D below
    }
    els.push_back(std::move(e));
  }
  std::vector<std::vector<GPState>> state(els.size()), trial(els.size());
  for (std::size_t k = 0; k < els.size(); ++k) state[k].assign(els[k].gp.size(), GPState{});
  ContinuumPlasticResult out;
  VecX u = VecX::Zero(nf);
  double tA = 0, tB = 0;  // the time step of the current increment
  double lam_now = 0;     // load factor of the evaluation (prescribed displacements)
  VecX fblocked = VecX::Zero(kDofPerNode * d.n_nodes);  // internal forces on the blocked DOFs
  auto element_u = [&](const El& e, const VecX& uu) {
    const int nd = static_cast<int>(e.dof.size());
    VecX ue(nd);
    for (int i = 0; i < nd; ++i) {
      const int fi = e.dof[static_cast<std::size_t>(i)];
      ue(i) = fi >= 0 ? uu(fi) : lam_now * presc[static_cast<std::size_t>(e.full[static_cast<std::size_t>(i)])];
    }
    return ue;
  };
  auto assemble = [&](const VecX& uu, VecX& fint, SpMat& Kt) {
    fint = VecX::Zero(nf);
    fblocked.setZero();
    std::vector<Triplet> t;
    for (std::size_t k = 0; k < els.size(); ++k) {
      El& e = els[k];
      const int nd = static_cast<int>(e.dof.size());
      const VecX ue = element_u(e, uu);
      VecX fe = VecX::Zero(nd);
      MatX ke = MatX::Zero(nd, nd);
      trial[k].resize(e.gp.size());
      for (std::size_t q = 0; q < e.gp.size(); ++q) {
        const GPt& g = e.gp[q];
        if (opt.large_displacements) {
          // total Lagrangian: F = I + H, E = (H + H' + H'H) / 2, S = D E
          const int dim = static_cast<int>(g.dX.cols()), nn = static_cast<int>(g.dX.rows());
          MatX H = MatX::Zero(dim, dim);
          for (int i = 0; i < nn; ++i)
            for (int a = 0; a < dim; ++a)
              for (int k2 = 0; k2 < dim; ++k2) H(a, k2) += ue(dim * i + a) * g.dX(i, k2);
          const MatX F = MatX::Identity(dim, dim) + H;
          const MatX Eg = 0.5 * (H + H.transpose() + H.transpose() * H);
          const MatX D = elasticity(e.c->kind, e.E, e.nu);
          const int ns = dim == 2 ? 3 : 6;
          VecX Ev(ns);
          if (dim == 2) Ev << Eg(0, 0), Eg(1, 1), 2 * Eg(0, 1);
          else Ev << Eg(0, 0), Eg(1, 1), Eg(2, 2), 2 * Eg(0, 1), 2 * Eg(1, 2), 2 * Eg(2, 0);
          const VecX Sv = D * Ev;
          MatX St(dim, dim);
          if (dim == 2) St << Sv(0), Sv(2), Sv(2), Sv(1);
          else St << Sv(0), Sv(3), Sv(5), Sv(3), Sv(1), Sv(4), Sv(5), Sv(4), Sv(2);
          // B_nl: dE = B_nl du
          MatX Bn = MatX::Zero(ns, dim * nn);
          for (int i = 0; i < nn; ++i)
            for (int a = 0; a < dim; ++a) {
              // dE_kl = (F_ak dN_l + F_al dN_k) / 2 per unit u_ia
              auto dE = [&](int k2, int l2) { return 0.5 * (F(a, k2) * g.dX(i, l2) + F(a, l2) * g.dX(i, k2)); };
              const int col = dim * i + a;
              if (dim == 2) {
                Bn(0, col) = dE(0, 0);
                Bn(1, col) = dE(1, 1);
                Bn(2, col) = 2 * dE(0, 1);
              } else {
                Bn(0, col) = dE(0, 0);
                Bn(1, col) = dE(1, 1);
                Bn(2, col) = dE(2, 2);
                Bn(3, col) = 2 * dE(0, 1);
                Bn(4, col) = 2 * dE(1, 2);
                Bn(5, col) = 2 * dE(2, 0);
              }
            }
          fe += g.w * Bn.transpose() * Sv;
          ke += g.w * Bn.transpose() * D * Bn;
          for (int i = 0; i < nn; ++i)
            for (int j = 0; j < nn; ++j) {
              const double gij = g.dX.row(i) * St * g.dX.row(j).transpose();
              for (int a = 0; a < dim; ++a) ke(dim * i + a, dim * j + a) += g.w * gij;
            }
          continue;
        }
        if (e.c->kind == Continuum::Kind::PlaneStress) {
          const MatX D = elasticity(e.c->kind, e.E, e.nu);
          fe += g.w * g.B.transpose() * (D * (g.B * ue));
          ke += g.w * g.B.transpose() * D * g.B;
          continue;
        }
        const Vec6d eps = e.S * (g.B * ue);
        Vec6d sig;
        Mat6d D;
        if (e.fm) {
          const FrictionalPoint fp = frictional_update(*e.fm, e.E, e.nu, eps, state[k][q].ep, state[k][q].ebar);
          sig = fp.sigma;
          D = fp.D;
          trial[k][q].ep = fp.ep;
          trial[k][q].ebar = fp.ebar;
        } else {
          j2_update(eps, state[k][q], e.E, e.nu, e.m, sig, D, trial[k][q], tA, tB, e.gp[q].cfac);
        }
        fe += g.w * g.B.transpose() * (e.S.transpose() * sig);
        ke += g.w * g.B.transpose() * (e.S.transpose() * D * e.S) * g.B;
      }
      for (int i = 0; i < nd; ++i) {
        const int fi = e.dof[static_cast<std::size_t>(i)];
        if (fi < 0) {
          fblocked(e.full[static_cast<std::size_t>(i)]) += fe(i);
          continue;
        }
        fint(fi) += fe(i);
        for (int j = 0; j < nd; ++j) {
          const int fj = e.dof[static_cast<std::size_t>(j)];
          if (fj >= 0 && ke(i, j) != 0) t.emplace_back(fi, fj, ke(i, j));
        }
      }
    }
    Kt.resize(nf, nf);
    Kt.setFromTriplets(t.begin(), t.end());
  };
  // one load increment to factor l over the time step [tA, tB]: Newton from
  // the committed state; false when it finds no equilibrium
  std::string fail_msg;
  auto increment = [&](double l) {
    lam_now = l;
    const VecX F = l * f;
    double last = 1e300;
    for (int it = 0;; ++it) {
      if (it >= opt.max_iter) {
        fail_msg = "continuum plasticity: no convergence at load factor " + std::to_string(l) +
                   " (beyond the limit load, or reduce the increments)";
        return false;
      }
      VecX fint;
      SpMat Kt;
      assemble(u, fint, Kt);
      const VecX r = F - fint;
      // the residual's tolerance, or a Newton correction at the round-off
      // of the displacements (stiff materials: |K u| carries eps |K| |u|)
      const double scale = std::max(F.norm(), fblocked.norm());
      if (r.norm() <= opt.tol * std::max(scale, 1e-30) || last <= 1e-13 * std::max(u.norm(), 1e-300)) return true;
      ++out.iterations;
      VecX du;
      if (unsym) {
        Kt.makeCompressed();
        Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
        lu.compute(Kt);
        if (lu.info() != Eigen::Success) {
          fail_msg = "continuum plasticity: singular tangent at load factor " + std::to_string(l);
          return false;
        }
        du = lu.solve(r);
      } else {
        Eigen::SimplicialLDLT<SpMat> ld(Kt);
        if (ld.info() != Eigen::Success) {
          fail_msg = "continuum plasticity: singular tangent at load factor " + std::to_string(l);
          return false;
        }
        du = ld.solve(r);
      }
      if (!du.allFinite()) {
        fail_msg = "continuum plasticity: singular tangent at load factor " + std::to_string(l);
        return false;
      }
      u += du;
      last = du.norm();
    }
  };
  double lam = 0, time = 0;
  bool stopped = false;
  for (std::size_t li = 0; li < opt.lambdas.size() && !stopped; ++li) {
    const double target = opt.lambdas[li];
    const double from = lam, tfrom = time, tto = opt.times.empty() ? 0.0 : opt.times[li];
    // progress phi in [0, 1] over the leg: substeps equal increments, halved
    // from the last converged state on a failure (stop_at_limit)
    int k = 0, cuts = 0;
    double phi = 0, dphi = 1.0 / opt.substeps;
    while (phi < 1.0 - 1e-15) {
      const double phin = cuts == 0 ? static_cast<double>(k + 1) / opt.substeps : std::min(1.0, phi + dphi);
      const double l = cuts == 0 ? from + (target - from) * (k + 1) / opt.substeps : from + (target - from) * phin;
      tA = cuts == 0 ? tfrom + (tto - tfrom) * k / opt.substeps : tfrom + (tto - tfrom) * phi;
      tB = cuts == 0 ? tfrom + (tto - tfrom) * (k + 1) / opt.substeps : tfrom + (tto - tfrom) * phin;
      const VecX u0 = u;
      if (!increment(l)) {
        require(opt.stop_at_limit, fail_msg);
        u = u0;
        lam_now = from + (target - from) * phi;
        if (cuts < opt.limit_cuts) {
          dphi = (cuts == 0 ? 1.0 / opt.substeps : dphi) / 2;
          ++cuts;
          continue;
        }
        out.limit_reached = true;
        stopped = true;
        break;
      }
      state = trial;  // commit
      phi = phin;
      ++k;
    }
    lam = from + (target - from) * phi;
    if (!stopped || phi > 0) {
      if (stopped) {
        // re-evaluate the converged state for its reactions
        lam_now = lam;
        VecX fint;
        SpMat Kt;
        assemble(u, fint, Kt);
      }
    }
    time = tfrom + (tto - tfrom) * phi;
    if (stopped && phi == 0) break;
    std::vector<Vec6> snap(model.nodes.size(), Vec6::Zero()), reac(model.nodes.size(), Vec6::Zero());
    for (std::size_t ni = 0; ni < model.nodes.size(); ++ni)
      for (int c = 0; c < 6; ++c) {
        const int fi = d.free_of(static_cast<int>(ni), c);
        const std::size_t full = static_cast<std::size_t>(d.full(static_cast<int>(ni), c));
        if (fi >= 0) snap[ni](c) = u(fi);
        else {
          snap[ni](c) = lam * presc[full];
          reac[ni](c) = fblocked(static_cast<Eigen::Index>(full));
        }
      }
    out.node_disp.push_back(snap);
    out.reactions.push_back(reac);
    out.lambdas_done.push_back(lam);
  }
  // final Gauss-point stresses
  for (std::size_t k = 0; k < els.size(); ++k) {
    El& e = els[k];
    const VecX ue = element_u(e, u);
    std::vector<std::array<double, 7>> g7;
    for (std::size_t q = 0; q < e.gp.size(); ++q) {
      std::array<double, 7> a{};
      if (e.c->kind != Continuum::Kind::PlaneStress) {
        GPState st = state[k][q];
        const Vec6d eps = e.S * (e.gp[q].B * ue);
        Vec6d sig;
        Mat6d D;
        GPState o2;
        GPState zero_inc = st;
        j2_update(eps, zero_inc, e.E, e.nu, nullptr, sig, D, o2);  // committed: the elastic law on eps - ep
        for (int i = 0; i < 6; ++i) a[static_cast<std::size_t>(i)] = sig(i);
        a[6] = st.ebar;
      }
      g7.push_back(a);
    }
    out.gauss.push_back(g7);
  }
  return out;
}

}  // namespace keta

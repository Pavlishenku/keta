#include "keta/influence_surface.h"

#include <algorithm>
#include <cmath>

#include "keta/shell.h"

namespace keta {

Id NodeGrid::at(int i, int j) const {
  require(i >= 0 && i < nx && j >= 0 && j < ny, "influence surface: grid index out of range");
  return nodes[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
               static_cast<std::size_t>(i)];
}

NodeGrid build_node_grid(const Model& model, const std::vector<Id>& nodes, int nx, int ny,
                         double tol) {
  require(nx >= 2 && ny >= 2, "influence surface: a node grid needs nx >= 2 and ny >= 2");
  require(nodes.size() == static_cast<std::size_t>(nx) * static_cast<std::size_t>(ny),
          "influence surface: the node list must hold exactly nx * ny ids (row-major)");
  require(tol > 0, "influence surface: grid tolerance must be positive");
  NodeGrid g;
  g.nx = nx;
  g.ny = ny;
  g.nodes = nodes;

  auto pos = [&](int i, int j) {
    return model
        .node(nodes[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
                    static_cast<std::size_t>(i)])
        .x;
  };
  g.x.resize(static_cast<std::size_t>(nx));
  g.y.resize(static_cast<std::size_t>(ny));
  for (int i = 0; i < nx; ++i) g.x[static_cast<std::size_t>(i)] = pos(i, 0).x();
  for (int j = 0; j < ny; ++j) g.y[static_cast<std::size_t>(j)] = pos(0, j).y();

  double ext = std::max(g.x.back() - g.x.front(), g.y.back() - g.y.front());
  require(ext > 0, "influence surface: degenerate node grid (zero extent)");
  double eps = tol * ext;
  for (int i = 1; i < nx; ++i)
    require(g.x[static_cast<std::size_t>(i)] > g.x[static_cast<std::size_t>(i) - 1] + eps,
            "influence surface: grid X lines must be strictly increasing");
  for (int j = 1; j < ny; ++j)
    require(g.y[static_cast<std::size_t>(j)] > g.y[static_cast<std::size_t>(j) - 1] + eps,
            "influence surface: grid Y lines must be strictly increasing");
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i) {
      Vec3 p = pos(i, j);
      require(std::fabs(p.x() - g.x[static_cast<std::size_t>(i)]) <= eps &&
                  std::fabs(p.y() - g.y[static_cast<std::size_t>(j)]) <= eps,
              "influence surface: node grid is not plan-rectangular at (" + std::to_string(i) +
                  ", " + std::to_string(j) + ")");
    }
  std::vector<Id> sorted = nodes;
  std::sort(sorted.begin(), sorted.end());
  require(std::adjacent_find(sorted.begin(), sorted.end()) == sorted.end(),
          "influence surface: the same node appears twice in the grid");
  return g;
}

// ---------------------------------------------------------------------------
// The sweep
// ---------------------------------------------------------------------------

InfluenceSurface influence_surface(const Model& model, const NodeGrid& grid, const Vec3& dir,
                                   const ResponseKind& response) {
  FactorCache fac;
  return influence_surface(model, grid, dir, response, fac);
}

InfluenceSurface influence_surface(const Model& model, const NodeGrid& grid, const Vec3& dir,
                                   const ResponseKind& response, FactorCache& fac,
                                   const SolveHooks* hooks) {
  require(finite(dir) && dir.squaredNorm() > 0, "influence surface: bad load direction");
  std::vector<Load> unit;
  unit.reserve(grid.nodes.size());
  for (Id n : grid.nodes) unit.push_back(Load::nodal(n, dir));
  InfluenceSurface s;
  s.nx = grid.nx;
  s.ny = grid.ny;
  s.x = grid.x;
  s.y = grid.y;
  s.eta = unit_load_sweep(model, unit, response, fac, hooks);
  return s;
}

// ---------------------------------------------------------------------------
// Interpolation and integrals
// ---------------------------------------------------------------------------

double InfluenceSurface::node_eta(int i, int j) const {
  require(i >= 0 && i < nx && j >= 0 && j < ny, "influence surface: index out of range");
  return eta[static_cast<std::size_t>(j) * static_cast<std::size_t>(nx) +
             static_cast<std::size_t>(i)];
}

namespace {
// Index of the cell holding `v` in a strictly increasing gridline vector,
// clamped to the last cell (so the upper edge belongs to it).
int cell_of(const std::vector<double>& g, double v) {
  std::size_t lo = 0, hi = g.size() - 1;
  while (hi - lo > 1) {
    std::size_t mid = (lo + hi) / 2;
    if (g[mid] <= v)
      lo = mid;
    else
      hi = mid;
  }
  return static_cast<int>(lo);
}
}  // namespace

double InfluenceSurface::at(double px, double py) const {
  if (nx < 2 || ny < 2) return 0.0;
  if (!(px >= x.front() && px <= x.back() && py >= y.front() && py <= y.back())) return 0.0;
  int i = cell_of(x, px), j = cell_of(y, py);
  double hx = x[static_cast<std::size_t>(i) + 1] - x[static_cast<std::size_t>(i)];
  double hy = y[static_cast<std::size_t>(j) + 1] - y[static_cast<std::size_t>(j)];
  double tx = (px - x[static_cast<std::size_t>(i)]) / hx;
  double ty = (py - y[static_cast<std::size_t>(j)]) / hy;
  return (1 - tx) * (1 - ty) * node_eta(i, j) + tx * (1 - ty) * node_eta(i + 1, j) +
         (1 - tx) * ty * node_eta(i, j + 1) + tx * ty * node_eta(i + 1, j + 1);
}

double InfluenceSurface::integral_box(double x0, double x1, double y0, double y1) const {
  if (nx < 2 || ny < 2 || x1 <= x0 || y1 <= y0) return 0.0;
  x0 = std::max(x0, x.front());
  x1 = std::min(x1, x.back());
  y0 = std::max(y0, y.front());
  y1 = std::min(y1, y.back());
  if (x1 <= x0 || y1 <= y0) return 0.0;
  double sum = 0;
  // Exact separable integration of the bilinear interpolant: on a cell the
  // field is N_a(x) N_b(y) e_ab, so the integral over any sub-rectangle is
  // (int N_a dx)(int N_b dy) e_ab.
  for (int j = 0; j + 1 < ny; ++j) {
    double ya = y[static_cast<std::size_t>(j)], yb = y[static_cast<std::size_t>(j) + 1];
    double c0 = std::max(y0, ya), c1 = std::min(y1, yb);
    if (c1 <= c0) continue;
    double hy = yb - ya;
    // int (yb - t)/hy dt and int (t - ya)/hy dt over [c0, c1]
    double Iy1 = ((c1 - ya) * (c1 - ya) - (c0 - ya) * (c0 - ya)) / (2 * hy);
    double Iy0 = (c1 - c0) - Iy1;
    for (int i = 0; i + 1 < nx; ++i) {
      double xa = x[static_cast<std::size_t>(i)], xb = x[static_cast<std::size_t>(i) + 1];
      double d0 = std::max(x0, xa), d1 = std::min(x1, xb);
      if (d1 <= d0) continue;
      double hx = xb - xa;
      double Ix1 = ((d1 - xa) * (d1 - xa) - (d0 - xa) * (d0 - xa)) / (2 * hx);
      double Ix0 = (d1 - d0) - Ix1;
      sum += Ix0 * Iy0 * node_eta(i, j) + Ix1 * Iy0 * node_eta(i + 1, j) +
             Ix0 * Iy1 * node_eta(i, j + 1) + Ix1 * Iy1 * node_eta(i + 1, j + 1);
    }
  }
  return sum;
}

double InfluenceSurface::integral() const {
  if (nx < 2 || ny < 2) return 0.0;
  return integral_box(x.front(), x.back(), y.front(), y.back());
}

namespace {

struct PV {  // a polygon corner carrying the interpolated field value
  double px = 0, py = 0, v = 0;
};

// Integral of the POSITIVE (or negative) part of the linear function defined by
// three corner values over the triangle. Exact: the zero set of a linear
// function is a straight line, so a Sutherland-Hodgman clip against v >= 0
// followed by a fan integration (int over a triangle = A * mean of corners) is
// closed form.
double tri_signed_part(const PV& a, const PV& b, const PV& c, bool positive) {
  PV in[3] = {a, b, c};
  std::vector<PV> poly(in, in + 3);
  std::vector<PV> clipped;
  clipped.reserve(4);
  auto keep = [&](const PV& p) { return positive ? p.v >= 0 : p.v <= 0; };
  for (std::size_t k = 0; k < poly.size(); ++k) {
    const PV& p = poly[k];
    const PV& q = poly[(k + 1) % poly.size()];
    bool kp = keep(p), kq = keep(q);
    if (kp) clipped.push_back(p);
    if (kp != kq) {
      double t = p.v / (p.v - q.v);  // exact zero crossing of the linear field
      clipped.push_back({p.px + t * (q.px - p.px), p.py + t * (q.py - p.py), 0.0});
    }
  }
  if (clipped.size() < 3) return 0.0;
  double sum = 0;
  for (std::size_t k = 1; k + 1 < clipped.size(); ++k) {
    const PV& p0 = clipped[0];
    const PV& p1 = clipped[k];
    const PV& p2 = clipped[k + 1];
    double area = 0.5 * std::fabs((p1.px - p0.px) * (p2.py - p0.py) -
                                  (p2.px - p0.px) * (p1.py - p0.py));
    sum += area * (p0.v + p1.v + p2.v) / 3.0;
  }
  return sum;
}

}  // namespace

namespace {
double signed_surface_integral(const InfluenceSurface& s, bool positive) {
  if (s.nx < 2 || s.ny < 2) return 0.0;
  double sum = 0;
  for (int j = 0; j + 1 < s.ny; ++j)
    for (int i = 0; i + 1 < s.nx; ++i) {
      double xa = s.x[static_cast<std::size_t>(i)], xb = s.x[static_cast<std::size_t>(i) + 1];
      double ya = s.y[static_cast<std::size_t>(j)], yb = s.y[static_cast<std::size_t>(j) + 1];
      PV p00{xa, ya, s.node_eta(i, j)};
      PV p10{xb, ya, s.node_eta(i + 1, j)};
      PV p01{xa, yb, s.node_eta(i, j + 1)};
      PV p11{xb, yb, s.node_eta(i + 1, j + 1)};
      // Both diagonals, averaged: exact on a sign-uniform cell (it then equals
      // the bilinear cell integral), O(h^2) on a cell the zero contour crosses.
      double dA = tri_signed_part(p00, p10, p11, positive) +
                  tri_signed_part(p00, p11, p01, positive);
      double dB = tri_signed_part(p00, p10, p01, positive) +
                  tri_signed_part(p10, p11, p01, positive);
      sum += 0.5 * (dA + dB);
    }
  return sum;
}
}  // namespace

double InfluenceSurface::integral_pos() const { return signed_surface_integral(*this, true); }
double InfluenceSurface::integral_neg() const { return signed_surface_integral(*this, false); }

std::pair<double, double> InfluenceSurface::point_extreme() const {
  double mx = 0, mn = 0;
  for (double e : eta) {
    mx = std::max(mx, e);
    mn = std::min(mn, e);
  }
  return {mx, mn};
}

// ---------------------------------------------------------------------------
// Surface algebra
// ---------------------------------------------------------------------------

namespace {
void record(SurfaceExtremes& ex, bool& first, double E, double rx, double ry) {
  if (first || E > ex.max) {
    ex.max = E;
    ex.rx_max = rx;
    ex.ry_max = ry;
  }
  if (first || E < ex.min) {
    ex.min = E;
    ex.rx_min = rx;
    ex.ry_min = ry;
  }
  first = false;
}
}  // namespace

SurfaceExtremes patch_extreme(const InfluenceSurface& s, double a, double b, double p,
                              int steps_x, int steps_y) {
  require(a > 0 && b > 0 && finite(a) && finite(b), "influence surface: patch needs a, b > 0");
  require(finite(p), "influence surface: non-finite patch pressure");
  require(steps_x >= 1 && steps_y >= 1, "influence surface: patch sweep needs >= 1 step");
  require(s.nx >= 2 && s.ny >= 2, "influence surface: empty surface");
  double x0 = s.x.front() - a / 2, x1 = s.x.back() + a / 2;
  double y0 = s.y.front() - b / 2, y1 = s.y.back() + b / 2;
  SurfaceExtremes ex;
  bool first = true;
  for (int jy = 0; jy <= steps_y; ++jy) {
    double ry = y0 + (y1 - y0) * jy / steps_y;
    for (int ix = 0; ix <= steps_x; ++ix) {
      double rx = x0 + (x1 - x0) * ix / steps_x;
      double E = p * s.integral_box(rx - a / 2, rx + a / 2, ry - b / 2, ry + b / 2);
      record(ex, first, E, rx, ry);
    }
  }
  return ex;
}

SurfaceExtremes udl_surface_extreme(const InfluenceSurface& s, double p) {
  require(finite(p) && p >= 0, "influence surface: UDL pressure must be finite and >= 0");
  SurfaceExtremes ex;
  ex.max = p * s.integral_pos();
  ex.min = p * s.integral_neg();
  return ex;
}

WheelSet WheelSet::mirrored_x() const {
  WheelSet m;
  m.wheels.reserve(wheels.size());
  for (const auto& wl : wheels) m.wheels.push_back({-wl[0], wl[1], wl[2]});
  return m;
}

double wheel_set_value(const InfluenceSurface& s, const WheelSet& ws, double rx, double ry) {
  double E = 0;
  for (const auto& wl : ws.wheels) E += wl[2] * s.at(rx + wl[0], ry + wl[1]);
  return E;
}

SurfaceExtremes wheel_set_extreme(const InfluenceSurface& s, const WheelSet& ws, int steps_x,
                                  int steps_y) {
  require(!ws.wheels.empty(), "influence surface: empty wheel set");
  require(steps_x >= 1 && steps_y >= 1, "influence surface: wheel sweep needs >= 1 step");
  require(s.nx >= 2 && s.ny >= 2, "influence surface: empty surface");
  double dx_min = 0, dx_max = 0, dy_min = 0, dy_max = 0;
  for (const auto& wl : ws.wheels) {
    require(finite(wl[0]) && finite(wl[1]) && finite(wl[2]),
            "influence surface: non-finite wheel");
    dx_min = std::min(dx_min, wl[0]);
    dx_max = std::max(dx_max, wl[0]);
    dy_min = std::min(dy_min, wl[1]);
    dy_max = std::max(dy_max, wl[1]);
  }
  // Reference range so that every wheel can reach every point of the grid.
  double x0 = s.x.front() - dx_max, x1 = s.x.back() - dx_min;
  double y0 = s.y.front() - dy_max, y1 = s.y.back() - dy_min;
  SurfaceExtremes ex;
  bool first = true;
  for (int jy = 0; jy <= steps_y; ++jy) {
    double ry = y0 + (y1 - y0) * jy / steps_y;
    for (int ix = 0; ix <= steps_x; ++ix) {
      double rx = x0 + (x1 - x0) * ix / steps_x;
      record(ex, first, wheel_set_value(s, ws, rx, ry), rx, ry);
    }
  }
  return ex;
}

}  // namespace keta

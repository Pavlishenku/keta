#include "keta/influence_scatter.h"

#include <algorithm>
#include <cmath>
#include <map>

namespace keta {

namespace {
using mesh2::Vec2;
double cross(const Vec2& o, const Vec2& a, const Vec2& b) {
  return (a.x() - o.x()) * (b.y() - o.y()) - (a.y() - o.y()) * (b.x() - o.x());
}
// convex hull, CCW, collinear boundary points KEPT (they are nodes too)
std::vector<int> hull(const std::vector<Vec2>& p) {
  std::vector<int> ix(p.size());
  for (std::size_t i = 0; i < p.size(); ++i) ix[i] = static_cast<int>(i);
  std::sort(ix.begin(), ix.end(), [&](int a, int b) {
    return p[a].x() != p[b].x() ? p[a].x() < p[b].x() : p[a].y() < p[b].y();
  });
  std::vector<int> h;
  for (int pass = 0; pass < 2; ++pass) {
    const std::size_t base = h.size();
    for (int i : ix) {
      while (h.size() >= base + 2 && cross(p[h[h.size() - 2]], p[h.back()], p[i]) < 0) h.pop_back();
      h.push_back(i);
    }
    h.pop_back();  // the last point starts the other chain
    std::reverse(ix.begin(), ix.end());
  }
  return h;
}
}  // namespace

ScatteredSurface::ScatteredSurface(const std::vector<Vec2>& pts, const std::vector<double>& values) {
  require(pts.size() >= 3 && pts.size() == values.size(), "scattered surface: at least 3 points, one value each");
  for (std::size_t i = 0; i < pts.size(); ++i)
    require(pts[i].allFinite() && std::isfinite(values[i]), "scattered surface: finite points and values");
  const std::vector<int> h = hull(pts);
  std::vector<char> on_hull(pts.size(), 0);
  mesh2::MeshInput in;
  for (int i : h) {
    require(!on_hull[static_cast<std::size_t>(i)], "scattered surface: collinear or duplicate points");
    on_hull[static_cast<std::size_t>(i)] = 1;
    in.outline.push_back(pts[static_cast<std::size_t>(i)]);
  }
  for (std::size_t i = 0; i < pts.size(); ++i)
    if (!on_hull[i]) in.points.push_back(pts[i]);
  in.refine = false;
  mesh_ = mesh2::triangulate(in);
  require(mesh_.nodes.size() == pts.size(), "scattered surface: the triangulation added or merged nodes");
  // values in the mesh's node order (the nodes are the input points verbatim)
  std::map<std::pair<double, double>, double> val;
  for (std::size_t i = 0; i < pts.size(); ++i) val[{pts[i].x(), pts[i].y()}] = values[i];
  eta_.reserve(mesh_.nodes.size());
  for (const Vec2& n : mesh_.nodes) {
    const auto it = val.find({n.x(), n.y()});
    require(it != val.end(), "scattered surface: a mesh node is not an input point");
    eta_.push_back(it->second);
  }
  // bucket grid: about one triangle per cell
  lo_ = hi_ = pts[0];
  for (const Vec2& q : pts) {
    lo_ = lo_.cwiseMin(q);
    hi_ = hi_.cwiseMax(q);
  }
  const double w = std::max(hi_.x() - lo_.x(), 1e-300), hh = std::max(hi_.y() - lo_.y(), 1e-300);
  const double nT = static_cast<double>(mesh_.tris.size());
  gx_ = std::max(1, static_cast<int>(std::ceil(std::sqrt(nT * w / hh))));
  gy_ = std::max(1, static_cast<int>(std::ceil(nT / gx_)));
  cw_ = w / gx_;
  ch_ = hh / gy_;
  cells_.assign(static_cast<std::size_t>(gx_ * gy_), {});
  auto cx = [&](double x) { return std::clamp(static_cast<int>(std::floor((x - lo_.x()) / cw_)), 0, gx_ - 1); };
  auto cy = [&](double y) { return std::clamp(static_cast<int>(std::floor((y - lo_.y()) / ch_)), 0, gy_ - 1); };
  for (std::size_t t = 0; t < mesh_.tris.size(); ++t) {
    Vec2 a = mesh_.nodes[static_cast<std::size_t>(mesh_.tris[t][0])], b = a;
    for (int k = 1; k < 3; ++k) {
      a = a.cwiseMin(mesh_.nodes[static_cast<std::size_t>(mesh_.tris[t][k])]);
      b = b.cwiseMax(mesh_.nodes[static_cast<std::size_t>(mesh_.tris[t][k])]);
    }
    for (int j = cy(a.y()); j <= cy(b.y()); ++j)
      for (int i = cx(a.x()); i <= cx(b.x()); ++i)
        cells_[static_cast<std::size_t>(j * gx_ + i)].push_back(static_cast<int>(t));
  }
}

double ScatteredSurface::at(double x, double y) const {
  const double tol = 1e-12 * std::max(hi_.x() - lo_.x(), hi_.y() - lo_.y());
  if (x < lo_.x() - tol || x > hi_.x() + tol || y < lo_.y() - tol || y > hi_.y() + tol) return 0;
  const int i = std::clamp(static_cast<int>(std::floor((x - lo_.x()) / cw_)), 0, gx_ - 1);
  const int j = std::clamp(static_cast<int>(std::floor((y - lo_.y()) / ch_)), 0, gy_ - 1);
  const Vec2 p(x, y);
  for (int t : cells_[static_cast<std::size_t>(j * gx_ + i)]) {
    const auto& tr = mesh_.tris[static_cast<std::size_t>(t)];
    const Vec2 &a = mesh_.nodes[static_cast<std::size_t>(tr[0])], &b = mesh_.nodes[static_cast<std::size_t>(tr[1])],
               &c = mesh_.nodes[static_cast<std::size_t>(tr[2])];
    const double A = cross(a, b, c);
    if (A <= 0) continue;
    const double l0 = cross(p, b, c) / A, l1 = cross(a, p, c) / A, l2 = 1 - l0 - l1;
    const double e = -1e-12;
    if (l0 >= e && l1 >= e && l2 >= e)
      return l0 * eta_[static_cast<std::size_t>(tr[0])] + l1 * eta_[static_cast<std::size_t>(tr[1])] +
             l2 * eta_[static_cast<std::size_t>(tr[2])];
  }
  return 0;
}

double ScatteredSurface::integral() const {
  double s = 0;
  for (const auto& tr : mesh_.tris) {
    const double A = 0.5 * cross(mesh_.nodes[static_cast<std::size_t>(tr[0])], mesh_.nodes[static_cast<std::size_t>(tr[1])],
                                 mesh_.nodes[static_cast<std::size_t>(tr[2])]);
    s += A * (eta_[static_cast<std::size_t>(tr[0])] + eta_[static_cast<std::size_t>(tr[1])] +
              eta_[static_cast<std::size_t>(tr[2])]) / 3;
  }
  return s;
}

// integral of max(sign eta, 0): each triangle clipped by the straight zero
// line of its linear field, the kept polygon fanned into triangles
double ScatteredSurface::clipped(int sign) const {
  double s = 0;
  for (const auto& tr : mesh_.tris) {
    std::vector<std::pair<Vec2, double>> poly, out;
    for (int k = 0; k < 3; ++k)
      poly.push_back({mesh_.nodes[static_cast<std::size_t>(tr[k])], sign * eta_[static_cast<std::size_t>(tr[k])]});
    for (std::size_t k = 0; k < 3; ++k) {
      const auto& [P, fp] = poly[k];
      const auto& [Q, fq] = poly[(k + 1) % 3];
      if (fp >= 0) out.push_back(poly[k]);
      if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) {
        const double t = fp / (fp - fq);
        out.push_back({P + t * (Q - P), 0.0});
      }
    }
    for (std::size_t k = 1; k + 1 < out.size(); ++k) {
      const double A = 0.5 * cross(out[0].first, out[k].first, out[k + 1].first);
      s += A * (out[0].second + out[k].second + out[k + 1].second) / 3;
    }
  }
  return sign * s;
}
double ScatteredSurface::integral_pos() const { return clipped(1); }
double ScatteredSurface::integral_neg() const { return clipped(-1); }

ScatteredSurface influence_surface_scattered(const Model& model, const std::vector<Id>& nodes, const Vec3& dir,
                                             const ResponseKind& response) {
  require(dir.allFinite() && dir.squaredNorm() > 0, "scattered surface: bad load direction");
  std::vector<Load> unit;
  std::vector<Vec2> pts;
  unit.reserve(nodes.size());
  pts.reserve(nodes.size());
  for (Id n : nodes) {
    unit.push_back(Load::nodal(n, dir));
    const Vec3& x = model.node(n).x;
    pts.emplace_back(x.x(), x.y());
  }
  FactorCache fac;
  return ScatteredSurface(pts, unit_load_sweep(model, unit, response, fac));
}

double wheel_set_value(const ScatteredSurface& s, const WheelSet& ws, double rx, double ry) {
  double v = 0;
  for (const auto& w : ws.wheels) v += w[2] * s.at(rx + w[0], ry + w[1]);
  return v;
}

SurfaceExtremes wheel_set_extreme(const ScatteredSurface& s, const WheelSet& ws, int nx, int ny) {
  require(nx >= 1 && ny >= 1 && !ws.wheels.empty(), "scattered wheel sweep: steps >= 1, wheels");
  double dx0 = 0, dx1 = 0, dy0 = 0, dy1 = 0;
  for (const auto& w : ws.wheels) {
    dx0 = std::min(dx0, w[0]);
    dx1 = std::max(dx1, w[0]);
    dy0 = std::min(dy0, w[1]);
    dy1 = std::max(dy1, w[1]);
  }
  const double x0 = s.xmin() - dx1, x1 = s.xmax() - dx0, y0 = s.ymin() - dy1, y1 = s.ymax() - dy0;
  SurfaceExtremes e;
  bool first = true;
  for (int j = 0; j <= ny; ++j)
    for (int i = 0; i <= nx; ++i) {
      const double rx = x0 + (x1 - x0) * i / nx, ry = y0 + (y1 - y0) * j / ny;
      const double v = wheel_set_value(s, ws, rx, ry);
      if (first || v > e.max) {
        e.max = v;
        e.rx_max = rx;
        e.ry_max = ry;
      }
      if (first || v < e.min) {
        e.min = v;
        e.rx_min = rx;
        e.ry_min = ry;
      }
      first = false;
    }
  return e;
}

}  // namespace keta

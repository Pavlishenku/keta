#include "keta/section_exact.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <sstream>

namespace keta {

// ===========================================================================
// Small geometry helpers
// ===========================================================================

double polygon_area(const std::vector<Vec2>& ring) {
  double a = 0;
  std::size_t n = ring.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& p = ring[i];
    const Vec2& q = ring[(i + 1) % n];
    a += p.x() * q.y() - q.x() * p.y();
  }
  return a / 2.0;
}

bool point_in_polygon(const std::vector<Vec2>& ring, const Vec2& p) {
  // Crossing-number test; boundary membership is undefined (callers use
  // centroids / interior points, never boundary points).
  bool in = false;
  std::size_t n = ring.size();
  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    double yi = ring[i].x(), zi = ring[i].y();
    double yj = ring[j].x(), zj = ring[j].y();
    if ((zi > p.y()) != (zj > p.y())) {
      double t = (p.y() - zi) / (zj - zi);
      if (p.x() < yi + t * (yj - yi)) in = !in;
    }
  }
  return in;
}

namespace {

double seg_distance(const Vec2& p, const Vec2& a, const Vec2& b) {
  Vec2 d = b - a;
  double L2 = d.squaredNorm();
  if (!(L2 > 0)) return (p - a).norm();
  double t = std::min(1.0, std::max(0.0, (p - a).dot(d) / L2));
  return (p - (a + t * d)).norm();
}

// Domain = inside the outline and outside every hole. `constraints` are internal
// rings (material-region boundaries) that the mesh must also conform to: they do
// not change what is inside, only where element edges have to fall.
struct Domain {
  std::vector<Vec2> outline;                 // CCW
  std::vector<std::vector<Vec2>> holes;      // CCW
  std::vector<std::vector<Vec2>> constraints;
  bool contains(const Vec2& p) const {
    if (!point_in_polygon(outline, p)) return false;
    for (const auto& h : holes)
      if (point_in_polygon(h, p)) return false;
    return true;
  }
  double area() const {
    double a = polygon_area(outline);
    for (const auto& h : holes) a -= polygon_area(h);
    return a;
  }
  double min_distance(const Vec2& p) const {
    double d = std::numeric_limits<double>::infinity();
    auto ring = [&](const std::vector<Vec2>& r) {
      std::size_t n = r.size();
      for (std::size_t i = 0; i < n; ++i) d = std::min(d, seg_distance(p, r[i], r[(i + 1) % n]));
    };
    ring(outline);
    for (const auto& h : holes) ring(h);
    for (const auto& c : constraints) ring(c);
    return d;
  }
};

std::vector<Vec2> normalized_ring(const std::vector<Vec2>& in, const char* what) {
  require(in.size() >= 3, std::string("section_exact: ") + what + " needs >= 3 vertices");
  std::vector<Vec2> r;
  r.reserve(in.size());
  for (const Vec2& p : in) {
    require(std::isfinite(p.x()) && std::isfinite(p.y()),
            std::string("section_exact: non-finite vertex in ") + what);
    // Drop an exactly repeated closing vertex.
    if (!r.empty() && (p - r.back()).squaredNorm() == 0) continue;
    r.push_back(p);
  }
  if (r.size() > 3 && (r.front() - r.back()).squaredNorm() == 0) r.pop_back();
  require(r.size() >= 3, std::string("section_exact: degenerate ") + what);
  if (polygon_area(r) < 0) std::reverse(r.begin(), r.end());
  require(polygon_area(r) > 0, std::string("section_exact: zero-area ") + what);
  return r;
}

Domain domain_of(const ExactInput& in) {
  Domain d;
  d.outline = normalized_ring(in.outline, "outline");
  for (const auto& h : in.holes) d.holes.push_back(normalized_ring(h, "hole"));
  // A material region's boundary must fall on element edges, else an element
  // straddles two materials and silently takes one of them.
  for (const ExactRegion& r : in.regions)
    if (r.outline.size() >= 3) d.constraints.push_back(normalized_ring(r.outline, "region"));
  double a = d.area();
  require(a > 0, "section_exact: holes consume the whole outline");
  return d;
}

// Point cloud with coincidence merging: a region boundary that runs along the
// outline must reuse the outline's samples, not duplicate them.
struct PointSet {
  std::vector<Vec2> pts;
  std::map<std::pair<long long, long long>, int> grid;
  double q = 1e-12;
  int add(const Vec2& p) {
    long long ky = std::llround(p.x() / q), kz = std::llround(p.y() / q);
    for (long long dy = -1; dy <= 1; ++dy)
      for (long long dz = -1; dz <= 1; ++dz) {
        auto it = grid.find({ky + dy, kz + dz});
        if (it != grid.end() && (pts[static_cast<std::size_t>(it->second)] - p).norm() <= 2 * q)
          return it->second;
      }
    int id = static_cast<int>(pts.size());
    pts.push_back(p);
    grid.emplace(std::make_pair(ky, kz), id);
    return id;
  }
};

// ===========================================================================
// Bowyer-Watson triangulator of the section (see the header)
// ===========================================================================

struct BwTri {
  int a = 0, b = 0, c = 0;
  Vec2 cc = Vec2::Zero();
  double r2 = 0;
  bool alive = true;
};

double tri_area2(const Vec2& a, const Vec2& b, const Vec2& c) {
  return (b.x() - a.x()) * (c.y() - a.y()) - (c.x() - a.x()) * (b.y() - a.y());
}

void set_circum(BwTri& t, const std::vector<Vec2>& p) {
  const Vec2& A = p[t.a];
  const Vec2& B = p[t.b];
  const Vec2& C = p[t.c];
  double d = 2.0 * ((B.x() - A.x()) * (C.y() - A.y()) - (C.x() - A.x()) * (B.y() - A.y()));
  require(std::fabs(d) > 0, "section_exact: degenerate triangle while meshing");
  double bx = B.x() - A.x(), by = B.y() - A.y();
  double cx = C.x() - A.x(), cy = C.y() - A.y();
  double b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
  double ux = (cy * b2 - by * c2) / d;
  double uy = (bx * c2 - cx * b2) / d;
  t.cc = A + Vec2(ux, uy);
  t.r2 = ux * ux + uy * uy;
}

// Bowyer-Watson incremental Delaunay over `pts` (deterministic insertion order).
// The cavity is grown by flood fill from the triangle containing the point, so
// a roundoff-inclusive in-circle test can never disconnect it.
std::vector<std::array<int, 3>> delaunay(std::vector<Vec2>& pts) {
  std::size_t n0 = pts.size();
  require(n0 >= 3, "section_exact: too few mesh points");
  Vec2 lo = pts[0], hi = pts[0];
  for (const Vec2& p : pts) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  Vec2 mid = 0.5 * (lo + hi);
  double R = 20.0 * std::max(hi.x() - lo.x(), hi.y() - lo.y());
  require(R > 0, "section_exact: degenerate point cloud");
  pts.push_back(mid + Vec2(0.0, 3.0 * R));
  pts.push_back(mid + Vec2(-3.0 * R, -1.5 * R));
  pts.push_back(mid + Vec2(3.0 * R, -1.5 * R));

  std::vector<BwTri> tris;
  {
    BwTri t;
    t.a = static_cast<int>(n0);
    t.b = static_cast<int>(n0 + 1);
    t.c = static_cast<int>(n0 + 2);
    if (tri_area2(pts[t.a], pts[t.b], pts[t.c]) < 0) std::swap(t.b, t.c);
    set_circum(t, pts);
    tris.push_back(t);
  }

  // Directed-edge ownership: a valid triangulation gives each directed edge to
  // exactly one triangle, so it doubles as the adjacency structure.
  std::map<std::pair<int, int>, int> owner;
  auto claim = [&](int ti) {
    const BwTri& t = tris[static_cast<std::size_t>(ti)];
    owner[{t.a, t.b}] = ti;
    owner[{t.b, t.c}] = ti;
    owner[{t.c, t.a}] = ti;
  };
  auto release = [&](int ti) {
    const BwTri& t = tris[static_cast<std::size_t>(ti)];
    int e[3][2] = {{t.a, t.b}, {t.b, t.c}, {t.c, t.a}};
    for (auto& ed : e) {
      auto it = owner.find({ed[0], ed[1]});
      if (it != owner.end() && it->second == ti) owner.erase(it);
    }
  };
  claim(0);

  std::vector<int> cavity, stack;
  std::vector<int> stamp(1, 0);
  int generation = 0;
  for (std::size_t ip = 0; ip < n0; ++ip) {
    const Vec2& p = pts[ip];
    // Seed: the triangle containing p. Scanned newest-first — consecutive input
    // points are neighbours in the sampling order, so the walk is short.
    int seed = -1;
    for (std::size_t k = tris.size(); k-- > 0;) {
      if (!tris[k].alive) continue;
      const BwTri& t = tris[k];
      if (tri_area2(pts[t.a], pts[t.b], p) >= 0 && tri_area2(pts[t.b], pts[t.c], p) >= 0 &&
          tri_area2(pts[t.c], pts[t.a], p) >= 0) {
        seed = static_cast<int>(k);
        break;
      }
    }
    require(seed >= 0, "section_exact: mesh point outside the working hull");

    stamp.resize(tris.size(), -1);
    ++generation;
    cavity.clear();
    stack.assign(1, seed);
    stamp[static_cast<std::size_t>(seed)] = generation;
    cavity.push_back(seed);
    while (!stack.empty()) {
      int ti = stack.back();
      stack.pop_back();
      const BwTri t = tris[static_cast<std::size_t>(ti)];
      int e[3][2] = {{t.a, t.b}, {t.b, t.c}, {t.c, t.a}};
      for (auto& ed : e) {
        auto it = owner.find({ed[1], ed[0]});  // the neighbour across this edge
        if (it == owner.end()) continue;
        int j = it->second;
        if (stamp[static_cast<std::size_t>(j)] == generation) continue;
        const BwTri& u = tris[static_cast<std::size_t>(j)];
        if ((p - u.cc).squaredNorm() < u.r2 * (1.0 + 1e-12)) {
          stamp[static_cast<std::size_t>(j)] = generation;
          cavity.push_back(j);
          stack.push_back(j);
        }
      }
    }

    // Cavity boundary = directed edges whose reverse is not in the cavity.
    std::map<std::pair<int, int>, int> count;
    for (int ti : cavity) {
      const BwTri& t = tris[static_cast<std::size_t>(ti)];
      int e[3][2] = {{t.a, t.b}, {t.b, t.c}, {t.c, t.a}};
      for (auto& ed : e) count[{ed[0], ed[1]}] += 1;
    }
    std::vector<std::pair<int, int>> hull;
    for (const auto& kv : count)
      if (count.find({kv.first.second, kv.first.first}) == count.end()) hull.push_back(kv.first);
    require(hull.size() >= 3, "section_exact: degenerate Delaunay cavity");
    for (int ti : cavity) {
      release(ti);
      tris[static_cast<std::size_t>(ti)].alive = false;
    }
    for (const auto& ed : hull) {
      BwTri t;
      t.a = ed.first;
      t.b = ed.second;
      t.c = static_cast<int>(ip);
      if (!(tri_area2(pts[t.a], pts[t.b], pts[t.c]) > 0)) continue;  // collinear sliver
      set_circum(t, pts);
      tris.push_back(t);
      stamp.push_back(-1);
      claim(static_cast<int>(tris.size()) - 1);
    }
  }

  std::vector<std::array<int, 3>> out;
  for (const BwTri& t : tris) {
    if (!t.alive) continue;
    if (t.a >= static_cast<int>(n0) || t.b >= static_cast<int>(n0) ||
        t.c >= static_cast<int>(n0))
      continue;
    out.push_back({t.a, t.b, t.c});
  }
  pts.resize(n0);
  return out;
}

// One smart-Laplacian pass: an interior node moves to the centroid of its
// neighbours only when every incident triangle keeps a positive area and the
// worst incident quality does not get worse.
void smooth_interior(std::vector<Vec2>& pts, const std::vector<std::array<int, 3>>& tris,
                     const std::vector<char>& fixed) {
  std::vector<std::vector<int>> inc(pts.size());
  for (std::size_t i = 0; i < tris.size(); ++i)
    for (int k = 0; k < 3; ++k) inc[tris[i][k]].push_back(static_cast<int>(i));
  auto quality = [&](int ti, const std::vector<Vec2>& q) {
    const Vec2& a = q[tris[ti][0]];
    const Vec2& b = q[tris[ti][1]];
    const Vec2& c = q[tris[ti][2]];
    double ar = tri_area2(a, b, c) / 2.0;
    double l2 = (b - a).squaredNorm() + (c - b).squaredNorm() + (a - c).squaredNorm();
    return l2 > 0 ? 4.0 * std::sqrt(3.0) * ar / l2 : -1.0;
  };
  for (std::size_t i = 0; i < pts.size(); ++i) {
    if (fixed[i] || inc[i].empty()) continue;
    Vec2 acc = Vec2::Zero();
    int cnt = 0;
    for (int ti : inc[i])
      for (int k = 0; k < 3; ++k)
        if (tris[ti][k] != static_cast<int>(i)) {
          acc += pts[tris[ti][k]];
          ++cnt;
        }
    if (cnt == 0) continue;
    double before = 1e300;
    for (int ti : inc[i]) before = std::min(before, quality(ti, pts));
    Vec2 old = pts[i];
    pts[i] = acc / static_cast<double>(cnt);
    double after = 1e300;
    for (int ti : inc[i]) after = std::min(after, quality(ti, pts));
    if (!(after > before)) pts[i] = old;
  }
}

// Sample a ring into sub-edges of length <= h, keeping every input vertex.
void sample_ring(const std::vector<Vec2>& ring, double h, PointSet& ps,
                 std::vector<std::pair<int, int>>& edges) {
  std::size_t n = ring.size();
  int first = -1, prev = -1;
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = ring[i];
    const Vec2& b = ring[(i + 1) % n];
    double L = (b - a).norm();
    int k = std::max(1, static_cast<int>(std::ceil(L / h - 1e-9)));
    for (int j = 0; j < k; ++j) {
      double t = static_cast<double>(j) / k;
      int id = ps.add(a + t * (b - a));
      if (first < 0) first = id;
      if (prev >= 0 && prev != id) edges.push_back({prev, id});
      prev = id;
    }
  }
  if (prev >= 0 && first >= 0 && prev != first) edges.push_back({prev, first});
}

Tri6Mesh assign_materials(Tri6Mesh m, const ExactInput& in) {
  m.mat.assign(m.tris.size(), 0);
  if (in.regions.empty()) return m;
  for (std::size_t i = 0; i < m.tris.size(); ++i) {
    Vec2 g = (m.nodes[m.tris[i][0]] + m.nodes[m.tris[i][1]] + m.nodes[m.tris[i][2]]) / 3.0;
    for (const ExactRegion& r : in.regions) {
      if (r.outline.size() < 3) continue;
      if (point_in_polygon(r.outline, g)) m.mat[i] = r.material;
    }
  }
  return m;
}

// Promote a linear triangle list to Tri6 (straight midside nodes, deduplicated).
Tri6Mesh to_tri6(const std::vector<Vec2>& corners, const std::vector<std::array<int, 3>>& tris,
                 const std::vector<char>& on_boundary) {
  Tri6Mesh m;
  m.nodes = corners;
  std::map<std::pair<int, int>, int> mid;
  auto midside = [&](int a, int b) {
    std::pair<int, int> key(std::min(a, b), std::max(a, b));
    auto it = mid.find(key);
    if (it != mid.end()) return it->second;
    int id = static_cast<int>(m.nodes.size());
    m.nodes.push_back(0.5 * (corners[a] + corners[b]));
    mid.emplace(key, id);
    return id;
  };
  m.tris.reserve(tris.size());
  for (const auto& t : tris) {
    std::array<int, 6> e{t[0], t[1], t[2], 0, 0, 0};
    e[3] = midside(t[0], t[1]);
    e[4] = midside(t[1], t[2]);
    e[5] = midside(t[2], t[0]);
    m.tris.push_back(e);
  }
  for (std::size_t i = 0; i < on_boundary.size(); ++i)
    if (on_boundary[i]) m.boundary_nodes.push_back(static_cast<int>(i));
  return m;
}

Tri6Mesh mesh_once(const ExactInput& in, const Domain& dom, double h) {
  Vec2 lo = dom.outline[0], hi = dom.outline[0];
  for (const Vec2& p : dom.outline) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  PointSet ps;
  ps.q = std::max(1e-14, 1e-9 * std::max(hi.x() - lo.x(), hi.y() - lo.y()));
  std::vector<std::pair<int, int>> bedges;
  sample_ring(dom.outline, h, ps, bedges);
  for (const auto& hole : dom.holes) sample_ring(hole, h, ps, bedges);
  int n_outer = static_cast<int>(ps.pts.size());
  for (const auto& c : dom.constraints) sample_ring(c, h, ps, bedges);
  std::vector<Vec2>& pts = ps.pts;
  int n_boundary = static_cast<int>(pts.size());
  // Triangular lattice, kept 0.65 h clear of every boundary segment so that
  // every boundary sub-edge stays a Gabriel edge (see the header).
  double dz = h * std::sqrt(3.0) / 2.0;
  int rows = static_cast<int>(std::floor((hi.y() - lo.y()) / dz)) + 1;
  for (int r = 0; r <= rows; ++r) {
    double z = lo.y() + r * dz;
    double y0 = lo.x() + ((r % 2) ? 0.5 * h : 0.0);
    int cols = static_cast<int>(std::floor((hi.x() - y0) / h)) + 1;
    for (int c = 0; c <= cols; ++c) {
      Vec2 p(y0 + c * h, z);
      if (!dom.contains(p)) continue;
      if (dom.min_distance(p) < 0.65 * h) continue;
      pts.push_back(p);
    }
  }
  require(static_cast<int>(pts.size()) <= in.max_nodes,
          "section_exact: mesh point budget exceeded (raise max_nodes or target_h)");

  std::vector<std::array<int, 3>> raw = delaunay(pts);
  std::vector<std::array<int, 3>> tris;
  tris.reserve(raw.size());
  for (const auto& t : raw) {
    Vec2 g = (pts[t[0]] + pts[t[1]] + pts[t[2]]) / 3.0;
    if (!dom.contains(g)) continue;
    std::array<int, 3> u = t;
    if (tri_area2(pts[u[0]], pts[u[1]], pts[u[2]]) < 0) std::swap(u[1], u[2]);
    tris.push_back(u);
  }
  require(!tris.empty(), "section_exact: empty mesh after domain filtering");

  // Verification 1: every boundary sub-edge was recovered.
  std::set<std::pair<int, int>> have;
  for (const auto& t : tris) {
    int e[3][2] = {{t[0], t[1]}, {t[1], t[2]}, {t[2], t[0]}};
    for (auto& ed : e) have.insert({std::min(ed[0], ed[1]), std::max(ed[0], ed[1])});
  }
  for (const auto& ed : bedges) {
    std::pair<int, int> key(std::min(ed.first, ed.second), std::max(ed.first, ed.second));
    if (have.find(key) == have.end()) return Tri6Mesh{};  // caller refines and retries
  }
  // Verification 2: the meshed area matches the exact Green integral.
  double area = 0;
  for (const auto& t : tris) area += tri_area2(pts[t[0]], pts[t[1]], pts[t[2]]) / 2.0;
  if (std::fabs(area - dom.area()) > 1e-9 * dom.area()) return Tri6Mesh{};

  std::vector<char> fixed(pts.size(), 0);
  for (int i = 0; i < n_boundary; ++i) fixed[i] = 1;  // outline, holes and interfaces
  if (in.smooth) smooth_interior(pts, tris, fixed);

  std::vector<char> on_boundary(pts.size(), 0);
  for (int i = 0; i < n_outer; ++i) on_boundary[i] = 1;  // domain boundary only
  Tri6Mesh m = to_tri6(pts, tris, on_boundary);
  return assign_materials(std::move(m), in);
}

}  // namespace

double Tri6Mesh::area() const {
  double a = 0;
  for (const auto& t : tris)
    a += tri_area2(nodes[t[0]], nodes[t[1]], nodes[t[2]]) / 2.0;
  return a;
}

Tri6Mesh mesh_section(const ExactInput& in) {
  Domain dom = domain_of(in);
  Vec2 lo = dom.outline[0], hi = dom.outline[0];
  for (const Vec2& p : dom.outline) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  double h = in.target_h;
  if (!(h > 0)) h = std::max(hi.x() - lo.x(), hi.y() - lo.y()) / 20.0;
  require(std::isfinite(h) && h > 0, "section_exact: invalid target_h");
  // Up to three halvings: a target size that does not resolve the thinnest wall
  // loses the Gabriel property and the boundary check above rejects the mesh.
  for (int attempt = 0; attempt < 4; ++attempt) {
    Tri6Mesh m = mesh_once(in, dom, h);
    if (!m.tris.empty()) return m;
    h *= 0.5;
  }
  fail("section_exact: could not build a boundary-conforming mesh — target_h does not "
       "resolve the thinnest wall (give a smaller ExactInput::target_h)");
}

Tri6Mesh tri6_upgrade(const mesh2::TriMesh& linear, const ExactInput& in) {
  require(!linear.tris.empty(), "section_exact: empty mesh2 triangulation");
  std::vector<Vec2> corners(linear.nodes.begin(), linear.nodes.end());
  std::vector<std::array<int, 3>> tris;
  tris.reserve(linear.tris.size());
  for (const auto& t : linear.tris) {
    std::array<int, 3> u{t[0], t[1], t[2]};
    if (tri_area2(corners[u[0]], corners[u[1]], corners[u[2]]) < 0) std::swap(u[1], u[2]);
    tris.push_back(u);
  }
  std::vector<char> on_boundary(corners.size(), 0);
  for (const auto& e : linear.boundary_edges) {
    on_boundary[e[0]] = 1;
    on_boundary[e[1]] = 1;
  }
  return assign_materials(to_tri6(corners, tris, on_boundary), in);
}

// ===========================================================================
// Tri6 warping / shear finite elements
// ===========================================================================

namespace {

struct GaussPt {
  double w = 0;
  double n[6] = {};
  double b[2][6] = {};  // derivatives wrt the natural coordinates (eta, xi)
};

void shape_at(double eta, double xi, double n[6], double b[2][6]) {
  double zeta = 1.0 - eta - xi;
  n[0] = eta * (2 * eta - 1);
  n[1] = xi * (2 * xi - 1);
  n[2] = zeta * (2 * zeta - 1);
  n[3] = 4 * eta * xi;
  n[4] = 4 * xi * zeta;
  n[5] = 4 * eta * zeta;
  double d0[6] = {4 * eta - 1, 0, 1 - 4 * zeta, 4 * xi, -4 * xi, 4 * (zeta - eta)};
  double d1[6] = {0, 4 * xi - 1, 1 - 4 * zeta, 4 * eta, 4 * (zeta - xi), -4 * eta};
  for (int i = 0; i < 6; ++i) {
    b[0][i] = d0[i];
    b[1][i] = d1[i];
  }
}

// Symmetric 4-point rule on the triangle, exact to cubic. The centroid weight
// is negative, which is correct for this rule.
const std::array<GaussPt, 4>& gauss_tri6() {
  static const std::array<GaussPt, 4> g = [] {
    const double pts[4][3] = {{1.0 / 3, 1.0 / 3, -27.0 / 48},
                              {0.6, 0.2, 25.0 / 48},
                              {0.2, 0.6, 25.0 / 48},
                              {0.2, 0.2, 25.0 / 48}};
    std::array<GaussPt, 4> out{};
    for (int k = 0; k < 4; ++k) {
      out[k].w = pts[k][2];
      shape_at(pts[k][0], pts[k][1], out[k].n, out[k].b);
    }
    return out;
  }();
  return g;
}

// Per-element geometry at one evaluation point: physical gradients of the six
// shape functions and the quadrature Jacobian j = detJ/2 (the element area for
// a straight-sided Tri6).
struct ElemPoint {
  double j = 0;
  double gy[6] = {};  // d/dy
  double gz[6] = {};  // d/dz
  double y = 0, z = 0;
};

ElemPoint elem_point(const Vec2 p[6], const double n[6], const double b[2][6]) {
  double dyde = 0, dzde = 0, dydx = 0, dzdx = 0;
  for (int a = 0; a < 6; ++a) {
    dyde += b[0][a] * p[a].x();
    dzde += b[0][a] * p[a].y();
    dydx += b[1][a] * p[a].x();
    dzdx += b[1][a] * p[a].y();
  }
  double detj = dyde * dzdx - dydx * dzde;
  require(detj > 0, "section_exact: non-positive Jacobian (element not counter-clockwise)");
  double inv = 1.0 / detj;
  double j00 = dzdx * inv, j01 = -dydx * inv, j10 = -dzde * inv, j11 = dyde * inv;
  ElemPoint ep;
  ep.j = 0.5 * detj;
  for (int a = 0; a < 6; ++a) {
    ep.gy[a] = j00 * b[0][a] + j10 * b[1][a];
    ep.gz[a] = j01 * b[0][a] + j11 * b[1][a];
    ep.y += n[a] * p[a].x();
    ep.z += n[a] * p[a].y();
  }
  return ep;
}

// Number of connected components of the element graph. Saint-Venant torsion
// couples the whole section through one warping field, so a disjoint section is
// rank-deficient no matter which DOF is pinned: detect it deterministically.
int components(const Tri6Mesh& m) {
  std::vector<int> parent(m.nodes.size());
  for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
  std::function<int(int)> find = [&](int x) {
    while (parent[x] != x) {
      parent[x] = parent[parent[x]];
      x = parent[x];
    }
    return x;
  };
  for (const auto& t : m.tris) {
    int r0 = find(t[0]);
    for (int k = 1; k < 6; ++k) {
      int r = find(t[k]);
      if (r != r0) parent[r] = r0;
    }
  }
  int c = 0;
  for (std::size_t i = 0; i < parent.size(); ++i)
    if (find(static_cast<int>(i)) == static_cast<int>(i)) ++c;
  return c;
}

}  // namespace

// ---------------------------------------------------------------------------

ExactSection analyze_exact(const Tri6Mesh& mesh_in, const std::vector<ExactMaterial>& mats_in) {
  ExactSection s;
  s.mesh = mesh_in;
  s.materials = mats_in;
  if (s.materials.empty()) s.materials.push_back(ExactMaterial{});
  require(!s.mesh.tris.empty(), "section_exact: empty mesh");
  require(s.mesh.mat.size() == s.mesh.tris.size(), "section_exact: material tags missing");
  for (const ExactMaterial& m : s.materials) {
    require(m.E > 0 && std::isfinite(m.E), "section_exact: material E must be > 0");
    require(m.nu > -1.0 && m.nu < 0.5, "section_exact: material nu must lie in (-1, 0.5)");
  }
  for (int mi : s.mesh.mat)
    require(mi >= 0 && mi < static_cast<int>(s.materials.size()),
            "section_exact: element material index out of range");
  require(components(s.mesh) == 1,
          "section_exact: the section falls into several disjoint regions — Saint-Venant "
          "torsion is not defined component by component; analyse each region separately");

  const std::size_t ne = s.mesh.tris.size();
  const std::size_t nn = s.mesh.nodes.size();
  s.E_ref = s.materials[0].E;

  std::vector<double> n_rel(ne), nu_e(ne);
  for (std::size_t e = 0; e < ne; ++e) {
    n_rel[e] = s.materials[s.mesh.mat[e]].E / s.E_ref;
    nu_e[e] = s.materials[s.mesh.mat[e]].nu;
  }

  const auto& gauss = gauss_tri6();
  auto gather = [&](std::size_t e, Vec2 p[6]) {
    for (int k = 0; k < 6; ++k) p[k] = s.mesh.nodes[s.mesh.tris[e][k]];
  };

  // --- pass 0: transformed area and first moments -> elastic centroid --------
  double EA = 0, Sy = 0, Sz = 0, area = 0, nu_w = 0;
  for (std::size_t e = 0; e < ne; ++e) {
    Vec2 p[6];
    gather(e, p);
    for (const GaussPt& g : gauss) {
      ElemPoint ep = elem_point(p, g.n, g.b);
      double w = g.w * ep.j;
      area += w;
      EA += w * n_rel[e];
      Sz += w * n_rel[e] * ep.y;  // integral n y dA
      Sy += w * n_rel[e] * ep.z;  // integral n z dA
      nu_w += w * n_rel[e] * nu_e[e];
    }
  }
  require(EA > 0, "section_exact: non-positive transformed area");
  s.A = area;
  s.EA_ref = EA;
  s.yc = Sz / EA;
  s.zc = Sy / EA;
  s.nu_ref = nu_w / EA;

  // Shift to the elastic centroid: every integral below is centroidal.
  for (Vec2& p : s.mesh.nodes) p -= Vec2(s.yc, s.zc);

  // --- pass 1: centroidal second moments ------------------------------------
  for (std::size_t e = 0; e < ne; ++e) {
    Vec2 p[6];
    gather(e, p);
    for (const GaussPt& g : gauss) {
      ElemPoint ep = elem_point(p, g.n, g.b);
      double w = g.w * ep.j * n_rel[e];
      s.Iy += w * ep.z * ep.z;
      s.Iz += w * ep.y * ep.y;
      s.Iyz += w * ep.y * ep.z;
    }
  }
  {
    double avg = 0.5 * (s.Iy + s.Iz);
    double rad = std::sqrt(0.25 * (s.Iy - s.Iz) * (s.Iy - s.Iz) + s.Iyz * s.Iyz);
    s.I1 = avg + rad;
    s.I2 = avg - rad;
    s.phi = std::atan2(s.Iy - s.I1, s.Iyz);
  }
  for (int b : s.mesh.boundary_nodes) {
    s.ymax = std::max(s.ymax, std::fabs(s.mesh.nodes[b].x()));
    s.zmax = std::max(s.zmax, std::fabs(s.mesh.nodes[b].y()));
  }

  // --- assemble K and the three loads ---------------------------------------
  std::vector<Triplet> trip;
  trip.reserve(ne * 36);
  std::vector<double> f_tors(nn, 0.0), f_psi_y(nn, 0.0), f_psi_z(nn, 0.0);
  for (std::size_t e = 0; e < ne; ++e) {
    Vec2 p[6];
    gather(e, p);
    double ke[6][6] = {};
    double lt[6] = {}, ly[6] = {}, lz[6] = {};
    for (const GaussPt& g : gauss) {
      ElemPoint ep = elem_point(p, g.n, g.b);
      double w = g.w * ep.j * n_rel[e];
      for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) ke[i][j] += w * (ep.gy[i] * ep.gy[j] + ep.gz[i] * ep.gz[j]);
      // Torsion: B^T (z, -y).
      for (int i = 0; i < 6; ++i) lt[i] += w * (ep.gy[i] * ep.z - ep.gz[i] * ep.y);
      // Transverse shear (Gruttmann & Wagner / Pilkey right-hand sides).
      double r = ep.y * ep.y - ep.z * ep.z;
      double q = 2.0 * ep.y * ep.z;
      double dy = s.Iy * r - s.Iyz * q, dz = s.Iyz * r + s.Iy * q;
      double hy = -s.Iyz * r + s.Iz * q, hz = -s.Iz * r - s.Iyz * q;
      double src_y = s.Iy * ep.y - s.Iyz * ep.z;
      double src_z = s.Iz * ep.z - s.Iyz * ep.y;
      // The nu/2 d term is the per-material compatibility correction; the
      // source is the equilibrium term div(tau) = n src V / ip2, whose factor
      // 2 (1 + nu_ref) must match delta_s. A per-element 2 (1 + nu_e) here
      // (as in sectionproperties) leaves a net load sum_e 2 (1 + nu_e)
      // integral n src != 0 when nu varies, so the Neumann problem is
      // inconsistent and the shear areas depend on the pinned node (on ring
      // order and mirror images, by up to ~20 %). Identical for a common nu.
      double half_nu = 0.5 * nu_e[e], two = 2.0 * (1.0 + s.nu_ref);
      for (int i = 0; i < 6; ++i) {
        ly[i] += w * (half_nu * (ep.gy[i] * dy + ep.gz[i] * dz) + two * g.n[i] * src_y);
        lz[i] += w * (half_nu * (ep.gy[i] * hy + ep.gz[i] * hz) + two * g.n[i] * src_z);
      }
    }
    for (int i = 0; i < 6; ++i) {
      int gi = s.mesh.tris[e][i];
      f_tors[gi] += lt[i];
      f_psi_y[gi] += ly[i];
      f_psi_z[gi] += lz[i];
      for (int j = 0; j < 6; ++j) {
        int gj = s.mesh.tris[e][j];
        if (gi == 0 || gj == 0) continue;  // pinned DOF
        trip.emplace_back(gi, gj, ke[i][j]);
      }
    }
  }
  trip.emplace_back(0, 0, 1.0);
  SpMat K(static_cast<int>(nn), static_cast<int>(nn));
  K.setFromTriplets(trip.begin(), trip.end());

  Eigen::SimplicialLDLT<SpMat> ldlt;
  ldlt.compute(K);
  require(ldlt.info() == Eigen::Success,
          "section_exact: the warping system is not positive definite (degenerate mesh)");
  MatX rhs(static_cast<int>(nn), 3);
  for (std::size_t i = 0; i < nn; ++i) {
    rhs(static_cast<int>(i), 0) = f_tors[i];
    rhs(static_cast<int>(i), 1) = f_psi_y[i];
    rhs(static_cast<int>(i), 2) = f_psi_z[i];
  }
  rhs.row(0).setZero();
  MatX sol = ldlt.solve(rhs);
  require(sol.allFinite(), "section_exact: warping solve produced non-finite values");
  s.omega.resize(nn);
  s.psi_y.resize(nn);
  s.psi_z.resize(nn);
  for (std::size_t i = 0; i < nn; ++i) {
    s.omega[i] = sol(static_cast<int>(i), 0);
    s.psi_y[i] = sol(static_cast<int>(i), 1);
    s.psi_z[i] = sol(static_cast<int>(i), 2);
  }

  // --- recovery -------------------------------------------------------------
  double ip2 = s.Iy * s.Iz - s.Iyz * s.Iyz;
  s.delta_s = 2.0 * (1.0 + s.nu_ref) * ip2;
  require(std::fabs(s.delta_s) > 1e-14 * std::max(std::fabs(ip2), 1.0),
          "section_exact: degenerate section (vanishing shear denominator)");

  double q_om = 0, i_om = 0, i_yom = 0, i_zom = 0;
  double sc_y = 0, sc_z = 0;
  double kap_y = 0, kap_z = 0, kap_yz = 0;
  for (std::size_t e = 0; e < ne; ++e) {
    Vec2 p[6];
    gather(e, p);
    double om[6], py[6], pz[6];
    for (int k = 0; k < 6; ++k) {
      om[k] = s.omega[s.mesh.tris[e][k]];
      py[k] = s.psi_y[s.mesh.tris[e][k]];
      pz[k] = s.psi_z[s.mesh.tris[e][k]];
    }
    for (const GaussPt& g : gauss) {
      ElemPoint ep = elem_point(p, g.n, g.b);
      double w = g.w * ep.j * n_rel[e];
      double om_g = 0;
      double gy_py = 0, gz_py = 0, gy_pz = 0, gz_pz = 0;
      for (int a = 0; a < 6; ++a) {
        om_g += g.n[a] * om[a];
        gy_py += ep.gy[a] * py[a];
        gz_py += ep.gz[a] * py[a];
        gy_pz += ep.gy[a] * pz[a];
        gz_pz += ep.gz[a] * pz[a];
      }
      double r2 = ep.y * ep.y + ep.z * ep.z;
      q_om += w * om_g;
      i_om += w * om_g * om_g;
      i_yom += w * ep.y * om_g;
      i_zom += w * ep.z * om_g;
      sc_y += w * (s.Iz * ep.y + s.Iyz * ep.z) * r2;
      sc_z += w * (s.Iy * ep.z + s.Iyz * ep.y) * r2;
      double r = ep.y * ep.y - ep.z * ep.z, q = 2.0 * ep.y * ep.z;
      double dy = s.Iy * r - s.Iyz * q, dz = s.Iyz * r + s.Iy * q;
      double hy = -s.Iyz * r + s.Iz * q, hz = -s.Iz * r - s.Iyz * q;
      double hn = 0.5 * nu_e[e];
      double vy0 = gy_py - hn * dy, vy1 = gz_py - hn * dz;
      double vz0 = gy_pz - hn * hy, vz1 = gz_pz - hn * hz;
      kap_y += w * (vy0 * vy0 + vy1 * vy1);
      kap_z += w * (vz0 * vz0 + vz1 * vz1);
      kap_yz += w * (vy0 * vz0 + vy1 * vz1);
    }
  }

  double warp_energy = 0, ft_psi_y = 0, ft_psi_z = 0;
  for (std::size_t i = 0; i < nn; ++i) {
    warp_energy += s.omega[i] * f_tors[i];
    ft_psi_y += f_tors[i] * s.psi_y[i];
    ft_psi_z += f_tors[i] * s.psi_z[i];
  }
  s.J = s.Iy + s.Iz - warp_energy;
  require(s.J > 0, "section_exact: non-positive torsion constant (degenerate mesh)");

  s.y_sc = (0.5 * s.nu_ref * sc_y - ft_psi_z) / s.delta_s;
  s.z_sc = (0.5 * s.nu_ref * sc_z + ft_psi_y) / s.delta_s;
  s.y_st = (s.Iyz * i_zom - s.Iz * i_yom) / ip2;
  s.z_st = (s.Iy * i_zom - s.Iyz * i_yom) / ip2;
  s.Iw = i_om - q_om * q_om / EA - s.z_sc * i_yom + s.y_sc * i_zom;
  s.omega_mean = q_om / EA;

  s.Avy = kap_y > 0 ? s.delta_s * s.delta_s / kap_y : 0.0;
  s.Avz = kap_z > 0 ? s.delta_s * s.delta_s / kap_z : 0.0;
  double kap_scale = std::max(std::fabs(kap_y), std::fabs(kap_z));
  s.Avyz = std::fabs(kap_yz) > 1e-10 * kap_scale ? s.delta_s * s.delta_s / kap_yz : 0.0;
  return s;
}

ExactSection analyze_exact(const ExactInput& in) {
  Tri6Mesh m = mesh_section(in);
  std::vector<ExactMaterial> mats = in.materials;
  if (mats.empty()) mats.push_back(ExactMaterial{});
  return analyze_exact(m, mats);
}

// ===========================================================================
// Evaluators
// ===========================================================================

int ExactSection::locate(double y, double z) const {
  Vec2 p(y, z);
  for (std::size_t e = 0; e < mesh.tris.size(); ++e) {
    const Vec2& a = mesh.nodes[mesh.tris[e][0]];
    const Vec2& b = mesh.nodes[mesh.tris[e][1]];
    const Vec2& c = mesh.nodes[mesh.tris[e][2]];
    double a2 = tri_area2(a, b, c);
    if (!(a2 > 0)) continue;
    double l0 = tri_area2(b, c, p) / a2;
    double l1 = tri_area2(c, a, p) / a2;
    double l2 = tri_area2(a, b, p) / a2;
    if (l0 >= -1e-9 && l1 >= -1e-9 && l2 >= -1e-9) return static_cast<int>(e);
  }
  return -1;
}

namespace {

// Natural coordinates of a point inside element `e` (straight-sided Tri6, so the
// barycentric coordinates of the corner triangle are the natural ones).
void natural_of(const ExactSection& s, int e, double y, double z, double& eta, double& xi) {
  const Vec2& a = s.mesh.nodes[s.mesh.tris[e][0]];
  const Vec2& b = s.mesh.nodes[s.mesh.tris[e][1]];
  const Vec2& c = s.mesh.nodes[s.mesh.tris[e][2]];
  double a2 = tri_area2(a, b, c);
  eta = tri_area2(b, c, Vec2(y, z)) / a2;
  xi = tri_area2(c, a, Vec2(y, z)) / a2;
}

// Vlasov sectorial coordinate about the shear centre, normalised to zero mean.
// The finite-element field satisfies grad(omega_fe) ~ (z, -y), i.e. it is the
// NEGATIVE of the classical sectorial area d(omega) = (y - y_p) dz - (z - z_p) dy
// used by the thin-wall module, hence the sign flip here. `om_fe` must already
// carry the zero-mean shift (omega - omega_mean).
double omega_sc_of(const ExactSection& s, double om_fe, double y, double z) {
  return -(om_fe + s.y_sc * z - s.z_sc * y);
}

}  // namespace

double ExactSection::omega_sc(int node) const {
  require(node >= 0 && node < static_cast<int>(omega.size()), "section_exact: node out of range");
  const Vec2& p = mesh.nodes[node];
  return omega_sc_of(*this, omega[node] - omega_mean, p.x(), p.y());
}

namespace {

// Natural coordinates of the six Tri6 nodes, matching `shape_at`.
const double kNodeNat[6][2] = {{1, 0}, {0, 1}, {0, 0}, {0.5, 0.5}, {0, 0.5}, {0.5, 0}};

// The action-independent part of the stress kernel at one point of a KNOWN
// element: position, material, and the interpolated warping / shear
// potentials and their gradients.
struct StressBasis {
  double y = 0, z = 0, n_rel = 1, nu_e = 0;
  double om = 0, gy_om = 0, gz_om = 0;
  double gy_py = 0, gz_py = 0, gy_pz = 0, gz_pz = 0;
};

StressBasis stress_basis(const ExactSection& s, int e, double eta, double xi) {
  double n[6], b[2][6];
  shape_at(eta, xi, n, b);
  Vec2 p[6];
  for (int k = 0; k < 6; ++k) p[k] = s.mesh.nodes[s.mesh.tris[e][k]];
  ElemPoint ep = elem_point(p, n, b);
  StressBasis sb;
  sb.y = ep.y;
  sb.z = ep.z;
  sb.n_rel = s.materials[s.mesh.mat[e]].E / s.E_ref;
  sb.nu_e = s.materials[s.mesh.mat[e]].nu;
  for (int k = 0; k < 6; ++k) {
    int gi = s.mesh.tris[e][k];
    sb.om += n[k] * s.omega[gi];
    sb.gy_om += ep.gy[k] * s.omega[gi];
    sb.gz_om += ep.gz[k] * s.omega[gi];
    sb.gy_py += ep.gy[k] * s.psi_y[gi];
    sb.gz_py += ep.gz[k] * s.psi_y[gi];
    sb.gy_pz += ep.gy[k] * s.psi_z[gi];
    sb.gz_pz += ep.gz[k] * s.psi_z[gi];
  }
  return sb;
}

// The stress kernel on a prepared basis (one basis serves any number of
// actions: the unit-field sweep evaluates six on it).
StressState stress_of_basis(const ExactSection& s, const StressBasis& sb, const SectionActions& a) {
  const double y = sb.y, z = sb.z, n_rel = sb.n_rel, nu_e = sb.nu_e;
  const double om = sb.om, gy_om = sb.gy_om, gz_om = sb.gz_om;
  const double gy_py = sb.gy_py, gz_py = sb.gz_py, gy_pz = sb.gy_pz, gz_pz = sb.gz_pz;

  StressState st;
  double det = s.Iy * s.Iz - s.Iyz * s.Iyz;
  st.sig = a.N / s.EA_ref +
           (a.My * (s.Iz * z - s.Iyz * y) - a.Mz * (s.Iy * y - s.Iyz * z)) / det;
  if (a.B != 0.0) {
    require(s.Iw > 0, "section_exact: bimoment requested on a section with Iw = 0");
    st.sig += a.B * omega_sc_of(s, om - s.omega_mean, y, z) / s.Iw;
  }
  st.sig *= n_rel;

  double ty = 0, tz = 0;
  if (a.T != 0.0) {
    ty += a.T / s.J * (gy_om - z);
    tz += a.T / s.J * (gz_om + y);
  }
  double r = y * y - z * z, q = 2.0 * y * z;
  double dy = s.Iy * r - s.Iyz * q, dz = s.Iyz * r + s.Iy * q;
  double hy = -s.Iyz * r + s.Iz * q, hz = -s.Iz * r - s.Iyz * q;
  double hn = 0.5 * nu_e;
  if (a.Vy != 0.0) {
    ty += a.Vy / s.delta_s * (gy_py - hn * dy);
    tz += a.Vy / s.delta_s * (gz_py - hn * dz);
  }
  if (a.Vz != 0.0) {
    ty += a.Vz / s.delta_s * (gy_pz - hn * hy);
    tz += a.Vz / s.delta_s * (gz_pz - hn * hz);
  }
  st.tau_y = n_rel * ty;
  st.tau_z = n_rel * tz;
  st.tau = std::hypot(st.tau_y, st.tau_z);
  st.vm = std::sqrt(st.sig * st.sig + 3.0 * st.tau * st.tau);
  return st;
}

// The stress kernel evaluated inside a KNOWN element, so a sweep over the whole
// mesh never pays for a point location.
StressState stress_in_element(const ExactSection& s, int e, double eta, double xi,
                              const SectionActions& a) {
  return stress_of_basis(s, stress_basis(s, e, eta, xi), a);
}

}  // namespace

StressState ExactSection::stress_at(double y, double z, const SectionActions& a) const {
  int e = locate(y, z);
  require(e >= 0, "section_exact: point outside the section");
  double eta, xi;
  natural_of(*this, e, y, z, eta, xi);
  return stress_in_element(*this, e, eta, xi, a);
}

double ExactSection::sigma_at(double y, double z, const SectionActions& a) const {
  return stress_at(y, z, a).sig;
}

double ExactSection::von_mises_max(const SectionActions& a) const {
  double best = 0;
  for (std::size_t e = 0; e < mesh.tris.size(); ++e)
    for (const auto& nat : kNodeNat)
      best = std::max(best,
                      stress_in_element(*this, static_cast<int>(e), nat[0], nat[1], a).vm);
  return best;
}

std::array<std::vector<StressState>, 6> ExactSection::nodal_unit_stresses() const {
  std::array<SectionActions, 6> unit;
  unit[0].N = 1;
  unit[1].Vy = 1;
  unit[2].Vz = 1;
  unit[3].My = 1;
  unit[4].Mz = 1;
  unit[5].T = 1;
  std::array<std::vector<StressState>, 6> out;
  for (auto& f : out) f.assign(mesh.nodes.size(), StressState{});
  std::vector<int> count(mesh.nodes.size(), 0);
  for (std::size_t e = 0; e < mesh.tris.size(); ++e)
    for (int k = 0; k < 6; ++k) {
      const StressBasis sb = stress_basis(*this, static_cast<int>(e), kNodeNat[k][0], kNodeNat[k][1]);
      const std::size_t node = static_cast<std::size_t>(mesh.tris[e][static_cast<std::size_t>(k)]);
      for (std::size_t i = 0; i < 6; ++i) {
        const StressState st = stress_of_basis(*this, sb, unit[i]);
        StressState& o = out[i][node];
        o.sig += st.sig;
        o.tau_y += st.tau_y;
        o.tau_z += st.tau_z;
      }
      ++count[node];
    }
  for (auto& f : out)
    for (std::size_t n = 0; n < f.size(); ++n) {
      if (count[n] == 0) continue;
      StressState& o = f[n];
      o.sig /= count[n];
      o.tau_y /= count[n];
      o.tau_z /= count[n];
      o.tau = std::hypot(o.tau_y, o.tau_z);
      o.vm = std::sqrt(o.sig * o.sig + 3 * o.tau * o.tau);
    }
  return out;
}

std::vector<StressState> ExactSection::nodal_stresses(const SectionActions& a) const {
  std::vector<StressState> out(mesh.nodes.size());
  std::vector<int> count(mesh.nodes.size(), 0);
  for (std::size_t e = 0; e < mesh.tris.size(); ++e)
    for (int k = 0; k < 6; ++k) {
      const StressState s = stress_in_element(*this, static_cast<int>(e), kNodeNat[k][0], kNodeNat[k][1], a);
      StressState& o = out[static_cast<std::size_t>(mesh.tris[e][static_cast<std::size_t>(k)])];
      o.sig += s.sig;
      o.tau_y += s.tau_y;
      o.tau_z += s.tau_z;
      ++count[static_cast<std::size_t>(mesh.tris[e][static_cast<std::size_t>(k)])];
    }
  for (std::size_t i = 0; i < out.size(); ++i) {
    if (count[i] == 0) continue;
    StressState& o = out[i];
    o.sig /= count[i];
    o.tau_y /= count[i];
    o.tau_z /= count[i];
    o.tau = std::hypot(o.tau_y, o.tau_z);
    o.vm = std::sqrt(o.sig * o.sig + 3 * o.tau * o.tau);
  }
  return out;
}

SigmaRange ExactSection::sigma_extremes(const SectionActions& a) const {
  require(!mesh.boundary_nodes.empty(), "section_exact: mesh carries no boundary nodes");
  double det = Iy * Iz - Iyz * Iyz;
  SigmaRange out;
  out.min.value = 1e300;
  out.max.value = -1e300;
  for (int b : mesh.boundary_nodes) {
    double y = mesh.nodes[b].x(), z = mesh.nodes[b].y();
    double sig = a.N / EA_ref + (a.My * (Iz * z - Iyz * y) - a.Mz * (Iy * y - Iyz * z)) / det;
    if (sig < out.min.value) out.min = {sig, y, z};
    if (sig > out.max.value) out.max = {sig, y, z};
  }
  return out;
}

SectionProps ExactSection::props() const {
  SectionProps p;
  p.A = EA_ref;
  p.Iy = Iy;
  p.Iz = Iz;
  p.J = J;
  p.Avy = Avy;
  p.Avz = Avz;
  p.Wely = zmax > 0 ? Iy / zmax : 0;
  p.Welz = ymax > 0 ? Iz / ymax : 0;
  return p;
}

// ===========================================================================
// Wire-in to Section
// ===========================================================================

void attach_exact(Section& s, std::shared_ptr<const ExactSection> ex) {
  require(ex != nullptr, "section_exact: null exact result");
  s.exact = ex;
  s.props.J = ex->J;
  s.props.Avy = ex->Avy;
  s.props.Avz = ex->Avz;
}

void attach_exact(Section& s, const ExactInput& in) {
  attach_exact(s, std::make_shared<const ExactSection>(analyze_exact(in)));
}

std::shared_ptr<const ExactSection> exact_of(const Section& s) { return s.exact; }

StressState exact_stress_at(const Section& s, double y, double z, const SectionActions& a) {
  require(s.exact != nullptr,
          "section_exact: section '" + s.name + "' carries no exact backend (call attach_exact)");
  return s.exact->stress_at(y, z, a);
}

SigmaRange exact_sigma_extremes(const Section& s, const SectionActions& a) {
  require(s.exact != nullptr,
          "section_exact: section '" + s.name + "' carries no exact backend (call attach_exact)");
  return s.exact->sigma_extremes(a);
}

// ===========================================================================
// DXF import
// ===========================================================================

namespace {

struct DxfPair {
  int code = 0;
  std::string value;
};

std::string trim(const std::string& s) {
  std::size_t a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  std::size_t b = s.find_last_not_of(" \t\r\n");
  return s.substr(a, b - a + 1);
}

double dxf_number(const std::string& v, const char* what) {
  try {
    std::size_t used = 0;
    double d = std::stod(v, &used);
    require(std::isfinite(d), std::string("dxf: non-finite ") + what);
    return d;
  } catch (const Error&) {
    throw;
  } catch (const std::exception&) {
    fail(std::string("dxf: cannot read ") + what + " from '" + v + "'");
  }
}

struct Loop {
  std::vector<Vec2> pts;
  bool closed = false;
};

// Sagitta-bounded angular step for polygonising an arc of radius r: the largest
// dphi with r (1 - cos(dphi/2)) <= tol (tol <= 0 -> relative default 1e-3 r).
double arc_step(double r, double sagitta_tol) {
  require(std::isfinite(r) && r > 0, "dxf: arc with non-positive radius");
  double tol = sagitta_tol > 0 ? sagitta_tol : 1e-3 * r;
  tol = std::min(tol, r);  // a sagitta beyond r is geometrically meaningless
  return 2.0 * std::acos(std::max(-1.0, 1.0 - tol / r));
}

// Expand LWPOLYLINE / POLYLINE bulges (group 42, bulge = tan(theta/4), theta
// the signed included angle, positive counter-clockwise) into intermediate arc
// points under the sagitta tolerance. The original vertices are kept verbatim
// (exact junction nodes, no duplicate at the chord/arc seam); only strictly
// interior points are inserted.
std::vector<Vec2> expand_bulges(const std::vector<Vec2>& v, const std::vector<double>& bulge,
                                bool closed, double sagitta_tol) {
  std::vector<Vec2> out;
  out.reserve(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) {
    out.push_back(v[i]);
    double b = bulge[i];
    if (b == 0.0) continue;
    require(std::isfinite(b), "dxf: non-finite bulge");
    bool last = i + 1 == v.size();
    require(!last || closed,
            "dxf: bulge on the last vertex of an OPEN polyline has no following vertex");
    const Vec2& p0 = v[i];
    const Vec2& p1 = v[last ? 0 : i + 1];
    double theta = 4.0 * std::atan(b);  // signed included angle, |theta| < 2 pi
    double d = (p1 - p0).norm();
    require(d > 0, "dxf: bulge on a zero-length segment");
    double R = d / (2.0 * std::sin(std::fabs(theta) / 2.0));
    // Centre: on the chord's mid-normal at apothem R - s (s = |b| d / 2 the
    // sagitta) — LEFT of the travel direction for a CCW (b > 0) arc, so the
    // arc itself bulges right of the chord (the convex side).
    Vec2 t = (p1 - p0) / d;
    Vec2 nrm(-t.y(), t.x());  // left normal
    double s = std::fabs(b) * d / 2.0;
    Vec2 c = 0.5 * (p0 + p1) + (b > 0 ? 1.0 : -1.0) * (R - s) * nrm;
    double phi0 = std::atan2(p0.y() - c.y(), p0.x() - c.x());
    int n = std::max(1, static_cast<int>(std::ceil(std::fabs(theta) / arc_step(R, sagitta_tol))));
    for (int k = 1; k < n; ++k) {
      double phi = phi0 + theta * k / n;
      out.push_back(c + R * Vec2(std::cos(phi), std::sin(phi)));
    }
  }
  return out;
}

}  // namespace

DxfSection import_dxf(const std::string& text, const DxfOptions& opt) {
  double scale = opt.scale;
  require(std::isfinite(scale) && scale > 0, "dxf: scale must be > 0");
  require(std::isfinite(opt.arc_sagitta), "dxf: arc_sagitta must be finite");
  // --- tokenise into (group code, value) pairs ---
  std::vector<DxfPair> pairs;
  {
    std::istringstream in(text);
    std::string a, b;
    while (std::getline(in, a)) {
      if (!std::getline(in, b)) break;
      std::string ca = trim(a);
      if (ca.empty()) continue;
      DxfPair p;
      p.code = static_cast<int>(dxf_number(ca, "group code"));
      p.value = trim(b);
      pairs.push_back(p);
    }
  }
  require(!pairs.empty(), "dxf: empty or non-ASCII DXF stream");

  std::vector<Loop> loops;                      // closed polylines
  std::vector<std::pair<Vec2, Vec2>> segments;  // loose LINE entities / open chains

  enum class Ent { Other, LwPoly, Poly, Vertex, Line, Arc };
  Ent kind = Ent::Other;
  bool in_entities = false, expect_section_name = false;
  std::vector<Vec2> pts;       // vertices of the entity being read
  std::vector<double> bulges;  // group-42 bulge per vertex (parallel to pts)
  int flags = 0;               // code 70 of the polyline being read
  bool poly_active = false;
  bool have_x = false, have_x1 = false;
  double x = 0, x1 = 0, y1 = 0;
  double arc_r = 0, arc_a0 = 0, arc_a1 = 0;  // ARC: radius, start/end angle [deg]
  bool have_r = false, have_a0 = false, have_a1 = false;

  // Turn the accumulated vertices into a closed loop, or into loose segments.
  auto close_vertices = [&](std::vector<Vec2> v, std::vector<double> bg, int f) {
    if (v.size() < 2) return;
    bool closed = (f & 1) != 0;
    if (!closed && v.size() >= 3 &&
        (v.front() - v.back()).norm() <= 1e-9 * (1.0 + v.front().norm())) {
      require(bg.back() == 0.0,
              "dxf: bulge on the duplicated closing vertex of a polyline");
      v.pop_back();
      bg.pop_back();
      closed = true;
    }
    bool any_bulge = false;
    for (double b : bg) any_bulge = any_bulge || b != 0.0;
    if (any_bulge) v = expand_bulges(v, bg, closed, opt.arc_sagitta);
    if (closed) {
      require(v.size() >= 3, "dxf: closed polyline with fewer than 3 vertices");
      Loop lp;
      lp.pts = std::move(v);
      lp.closed = true;
      loops.push_back(std::move(lp));
    } else {
      for (std::size_t i = 1; i < v.size(); ++i) segments.push_back({v[i - 1], v[i]});
    }
  };
  // Flush a self-contained entity. POLYLINE / VERTEX keep accumulating (their
  // vertex list and closed flag are only consumed by `finish_polyline`).
  auto finish = [&]() {
    if (kind == Ent::LwPoly) {
      close_vertices(pts, bulges, flags);
      pts.clear();
      bulges.clear();
      flags = 0;
    } else if (kind == Ent::Line) {
      require(!pts.empty() && have_x1, "dxf: LINE without both end points");
      segments.push_back({pts.front(), Vec2(x1, y1)});
      pts.clear();
      bulges.clear();
      flags = 0;
    } else if (kind == Ent::Arc) {
      // ARC: centre (10/20), radius (40), start/end angles 50/51 in degrees,
      // ALWAYS counter-clockwise from 50 to 51 (DXF convention). Polygonised
      // under the sagitta tolerance into loose chords, end points exact.
      require(!pts.empty() && have_r && have_a0 && have_a1,
              "dxf: ARC without centre, radius or angles");
      Vec2 c = pts.front();
      double a0 = arc_a0 * kPi / 180.0, a1 = arc_a1 * kPi / 180.0;
      double sweep = a1 - a0;
      while (sweep <= 0) sweep += 2.0 * kPi;  // CCW; 50 == 51 -> full turn
      int n = std::max(1, static_cast<int>(std::ceil(sweep / arc_step(arc_r, opt.arc_sagitta))));
      Vec2 prev = c + arc_r * Vec2(std::cos(a0), std::sin(a0));
      for (int k = 1; k <= n; ++k) {
        double phi = k == n ? a1 : a0 + sweep * k / n;  // exact junction at 51
        Vec2 p = c + arc_r * Vec2(std::cos(phi), std::sin(phi));
        segments.push_back({prev, p});
        prev = p;
      }
      pts.clear();
      bulges.clear();
      have_r = have_a0 = have_a1 = false;
    }
    have_x = false;
    have_x1 = false;
  };
  auto finish_polyline = [&]() {
    if (!poly_active) return;
    close_vertices(pts, bulges, flags);
    poly_active = false;
    pts.clear();
    bulges.clear();
    flags = 0;
  };

  for (const DxfPair& p : pairs) {
    if (p.code == 0) {
      finish();
      const std::string& e = p.value;
      if (e == "SECTION") {
        expect_section_name = true;
        kind = Ent::Other;
        continue;
      }
      if (e == "ENDSEC" || e == "EOF") {
        finish_polyline();
        in_entities = false;
        kind = Ent::Other;
        continue;
      }
      // the entities of the other sections (the BLOCKS of the dimensions a CAD
      // package writes: LINEs whose points are never read) are not geometry
      if (!in_entities) {
        kind = Ent::Other;
        continue;
      }
      if (e == "VERTEX") {
        kind = Ent::Vertex;  // keeps feeding the open POLYLINE's vertex list
        continue;
      }
      if (e == "SEQEND") {
        finish_polyline();
        kind = Ent::Other;
        continue;
      }
      finish_polyline();  // an entity that is not a VERTEX ends any open POLYLINE
      if (e == "LWPOLYLINE") {
        kind = Ent::LwPoly;
      } else if (e == "POLYLINE") {
        kind = Ent::Poly;
        poly_active = true;
      } else if (e == "LINE") {
        kind = Ent::Line;
      } else if (e == "ARC") {
        kind = Ent::Arc;
      } else {
        kind = Ent::Other;
      }
      continue;
    }
    if (p.code == 2 && expect_section_name) {
      expect_section_name = false;
      in_entities = (p.value == "ENTITIES");
      continue;
    }
    if (!in_entities) continue;
    if (kind == Ent::Other) continue;
    switch (p.code) {
      case 70:
        // The POLYLINE header carries the closed flag; a VERTEX's own code 70 is
        // a vertex flag and must not overwrite it.
        if (kind != Ent::Vertex) flags = static_cast<int>(dxf_number(p.value, "flags"));
        break;
      case 42: {
        // Vertex bulge (LWPOLYLINE / POLYLINE VERTEX): the group follows the
        // vertex's own 10/20 pair and describes the arc STARTING at it.
        double b = dxf_number(p.value, "bulge");
        if (b != 0.0) {
          require((kind == Ent::LwPoly || kind == Ent::Vertex) && !bulges.empty(),
                  "dxf: bulge outside a polyline vertex");
          bulges.back() = b;
        }
        break;
      }
      case 40:
        if (kind == Ent::Arc) {
          arc_r = dxf_number(p.value, "arc radius");
          have_r = true;
        }
        break;
      case 50:
        if (kind == Ent::Arc) {
          arc_a0 = dxf_number(p.value, "arc start angle");
          have_a0 = true;
        }
        break;
      case 51:
        if (kind == Ent::Arc) {
          arc_a1 = dxf_number(p.value, "arc end angle");
          have_a1 = true;
        }
        break;
      case 10:
        x = dxf_number(p.value, "x");
        have_x = true;
        break;
      case 20:
        require(have_x, "dxf: y ordinate without a preceding x");
        pts.push_back(Vec2(x, dxf_number(p.value, "y")));
        bulges.push_back(0.0);
        have_x = false;
        break;
      case 11:
        x1 = dxf_number(p.value, "x2");
        have_x1 = true;
        break;
      case 21:
        require(have_x1, "dxf: LINE end y without x");
        y1 = dxf_number(p.value, "y2");
        break;
      default: break;
    }
  }
  finish();
  finish_polyline();

  // --- chain loose segments into closed loops ---
  if (!segments.empty()) {
    double ext = 0;
    for (const auto& s : segments) ext = std::max({ext, s.first.norm(), s.second.norm()});
    double tol = std::max(1e-9, 1e-7 * ext);
    std::vector<char> used(segments.size(), 0);
    for (std::size_t i = 0; i < segments.size(); ++i) {
      if (used[i]) continue;
      used[i] = 1;
      Loop lp;
      lp.pts.push_back(segments[i].first);
      Vec2 head = segments[i].second;
      Vec2 start = segments[i].first;
      while ((head - start).norm() > tol) {
        bool found = false;
        for (std::size_t j = 0; j < segments.size() && !found; ++j) {
          if (used[j]) continue;
          if ((segments[j].first - head).norm() <= tol) {
            used[j] = 1;
            lp.pts.push_back(head);
            head = segments[j].second;
            found = true;
          } else if ((segments[j].second - head).norm() <= tol) {
            used[j] = 1;
            lp.pts.push_back(head);
            head = segments[j].first;
            found = true;
          }
        }
        require(found, "dxf: open loop — a chain of LINE entities does not close");
      }
      lp.pts.push_back(head);
      require(lp.pts.size() >= 3, "dxf: degenerate loop (fewer than 3 vertices)");
      lp.closed = true;
      loops.push_back(lp);
    }
  }

  require(!loops.empty(), "dxf: no closed LWPOLYLINE / POLYLINE / LINE loop in ENTITIES");
  // --- classify: largest loop is the outline, the rest are holes ---
  std::size_t big = 0;
  double best = 0;
  for (std::size_t i = 0; i < loops.size(); ++i) {
    double a = std::fabs(polygon_area(loops[i].pts));
    if (a > best) {
      best = a;
      big = i;
    }
  }
  require(best > 0, "dxf: all loops have zero area");
  DxfSection out;
  auto ccw = [&](std::vector<Vec2> r) {
    for (Vec2& p : r) p *= scale;
    if (polygon_area(r) < 0) std::reverse(r.begin(), r.end());
    return r;
  };
  out.outline = ccw(loops[big].pts);
  for (std::size_t i = 0; i < loops.size(); ++i)
    if (i != big) out.holes.push_back(ccw(loops[i].pts));
  return out;
}

DxfSection import_dxf(const std::string& text, double scale) {
  DxfOptions opt;
  opt.scale = scale;
  return import_dxf(text, opt);
}

DxfSection load_dxf(const std::string& path, const DxfOptions& opt) {
  std::ifstream f(path, std::ios::binary);
  require(f.good(), "dxf: cannot open '" + path + "'");
  std::ostringstream ss;
  ss << f.rdbuf();
  return import_dxf(ss.str(), opt);
}

DxfSection load_dxf(const std::string& path, double scale) {
  DxfOptions opt;
  opt.scale = scale;
  return load_dxf(path, opt);
}

}  // namespace keta

namespace keta {

namespace {
std::vector<Vec2> ring_rect(double y0, double z0, double y1, double z1) {
  return {{y0, z0}, {y1, z0}, {y1, z1}, {y0, z1}};
}
}  // namespace

bool energy_shear_areas(const Shape& shape, double nu, double& Avy, double& Avz) {
  require(std::isfinite(nu) && nu > -1 && nu < 0.5, "section: Poisson ratio outside (-1, 0.5)");
  ExactInput in;
  in.materials = {{"m", 1.0, nu}};
  double tmin = 0, ext = 0;
  if (const auto* r = std::get_if<ShapeRect>(&shape)) {
    // Saint-Venant's flexure of the rectangle in closed form. Shear V along
    // z (height h, width b): tau = V/2I (h^2/4 - z^2) e_z + curl psi, with
    // lap psi = c y, psi = 0 on the boundary, c = nu/(1+nu) V/I; the energy
    // integral gives
    //   A / (k V^2) U = 6/5 + A c'^2 [h b^5/720 - b^6/(8 pi^5) sum tanh(n pi h/b)/n^5],
    // c' = c/V, the sum read as zeta(5) less an exponentially small tail.
    // Shear along y: b and h exchanged.
    auto kz = [&](double b, double h) {
      const double A = b * h, I = b * h * h * h / 12, c = nu / (1 + nu) / I;
      const double zeta5 = 1.0369277551433699263;
      double tail = 0;
      for (int k = 1; k < 2000; ++k) {
        const double t = 1.0 - std::tanh(k * kPi * h / b);
        const double term = t / std::pow(k, 5);
        tail += term;
        if (term < 1e-19 * zeta5) break;
      }
      const double S = zeta5 - tail;
      const double U = 1.2 + A * c * c * (h * std::pow(b, 5) / 720 - std::pow(b, 6) / (8 * std::pow(kPi, 5)) * S);
      return A / U;
    };
    Avz = kz(r->b, r->h);
    Avy = kz(r->h, r->b);
    return true;
  } else if (const auto* c = std::get_if<ShapeCircle>(&shape)) {
    // Saint-Venant flexure energy of a solid circle, closed form
    const double A = kPi * c->d * c->d / 4;
    Avy = Avz = A * 6 * (1 + nu) * (1 + nu) / (7 + 14 * nu + 8 * nu * nu);
    return true;
  } else if (const auto* t = std::get_if<ShapeTube>(&shape)) {
    // and of a thick tube, m = r_in / r_out:
    //   k = 6 (1+nu)^2 (1+m^2)^2 / [(7+14nu+8nu^2)(1+m^2)^2 + 4(5+10nu+4nu^2) m^2]
    const double ro = t->d / 2, ri = ro - t->t, m = ri / ro, q = (1 + m * m) * (1 + m * m);
    const double A = kPi * (ro * ro - ri * ri);
    Avy = Avz = A * 6 * (1 + nu) * (1 + nu) * q /
                ((7 + 14 * nu + 8 * nu * nu) * q + 4 * (5 + 10 * nu + 4 * nu * nu) * m * m);
    return true;
  } else if (const auto* s = std::get_if<ShapeI>(&shape)) {
    const double h = s->h, b = s->b, tw = s->tw, tf = s->tf, rr = s->r;
    // outline with the four root fillets (quarter circles, 16 chords each)
    std::vector<Vec2> o;
    auto fillet = [&](double cy, double cz, double a0, double a1) {
      for (int k = 0; k <= 16; ++k) {
        const double a = a0 + (a1 - a0) * k / 16;
        o.push_back({cy + rr * std::cos(a), cz + rr * std::sin(a)});
      }
    };
    const double yw = tw / 2, zf = h / 2 - tf, q = 0.5 * kPi;
    o.push_back({-b / 2, -h / 2});
    o.push_back({b / 2, -h / 2});
    o.push_back({b / 2, -zf});
    if (rr > 0) fillet(yw + rr, -zf + rr, 3 * q, 2 * q);
    else o.push_back({yw, -zf});
    if (rr > 0) fillet(yw + rr, zf - rr, 2 * q, q);
    else o.push_back({yw, zf});
    o.push_back({b / 2, zf});
    o.push_back({b / 2, h / 2});
    o.push_back({-b / 2, h / 2});
    o.push_back({-b / 2, zf});
    if (rr > 0) fillet(-yw - rr, zf - rr, q, 0.0);
    else o.push_back({-yw, zf});
    if (rr > 0) fillet(-yw - rr, -zf + rr, 0.0, -q);
    else o.push_back({-yw, -zf});
    o.push_back({-b / 2, -zf});
    in.outline = o;
    tmin = std::min(tw, tf);
    ext = std::max(h, b);
  } else if (const auto* x = std::get_if<ShapeBox>(&shape)) {
    in.outline = ring_rect(-x->b / 2, -x->h / 2, x->b / 2, x->h / 2);
    in.holes = {ring_rect(-x->b / 2 + x->tw, -x->h / 2 + x->tf, x->b / 2 - x->tw, x->h / 2 - x->tf)};
    tmin = std::min(x->tw, x->tf);
    ext = std::max(x->h, x->b);
  } else if (const auto* x = std::get_if<ShapeT>(&shape)) {
    const double hw = x->h - x->tf;
    in.outline = {{-x->tw / 2, 0},     {x->tw / 2, 0},     {x->tw / 2, hw},
                  {x->b / 2, hw},      {x->b / 2, x->h},   {-x->b / 2, x->h},
                  {-x->b / 2, hw},     {-x->tw / 2, hw}};
    tmin = std::min(x->tw, x->tf);
    ext = std::max(x->h, x->b);
  } else if (const auto* x = std::get_if<ShapePolygon>(&shape)) {
    double ymin = 1e300, ymax = -1e300, zmin = 1e300, zmax = -1e300;
    for (const Vec3& v : x->outer) {
      in.outline.push_back({v.x(), v.y()});
      ymin = std::min(ymin, v.x());
      ymax = std::max(ymax, v.x());
      zmin = std::min(zmin, v.y());
      zmax = std::max(zmax, v.y());
    }
    for (const auto& hl : x->holes) {
      std::vector<Vec2> r;
      for (const Vec3& v : hl) r.push_back({v.x(), v.y()});
      in.holes.push_back(r);
    }
    ext = std::max(ymax - ymin, zmax - zmin);
    tmin = 0;  // walls unknown: the mesher checks the resolution
  } else {
    return false;
  }
  in.target_h = tmin > 0 ? std::min(ext / 40, tmin / 4) : ext / 40;
  // no closed form for these shapes: the warping solve, once per geometry
  // and Poisson ratio
  static std::mutex mu;
  static std::map<std::vector<double>, std::pair<double, double>> cache;
  std::vector<double> key{nu, static_cast<double>(shape.index())};
  for (const Vec2& v : in.outline) key.insert(key.end(), {v.x(), v.y()});
  for (const auto& hl : in.holes)
    for (const Vec2& v : hl) key.insert(key.end(), {v.x(), v.y()});
  {
    std::lock_guard<std::mutex> lk(mu);
    auto it = cache.find(key);
    if (it != cache.end()) {
      Avy = it->second.first;
      Avz = it->second.second;
      return true;
    }
  }
  const ExactSection s = analyze_exact(in);
  Avy = s.Avy;
  Avz = s.Avz;
  std::lock_guard<std::mutex> lk(mu);
  cache[key] = {Avy, Avz};
  return true;
}

}  // namespace keta

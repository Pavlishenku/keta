// 2-D constrained Delaunay triangulation with Ruppert-style refinement.
//
// Pipeline (every stage is deterministic; no container is iterated in an
// address- or hash-dependent order, so the output is byte-reproducible):
//
//   0. VALIDATION.  Finiteness of every coordinate; |signed area| of the
//      outline > 0; simplicity of the outline and of each hole (no proper
//      crossing, no fold-back spike); holes inside the outline and pairwise
//      disjoint; internal segments not crossing a boundary; extra `points`
//      inside the domain; no two input vertices closer than
//      kDupRel * bbox_diagonal.  Every violation is a keta::fail.
//
//   1. DELAUNAY.  Incremental Bowyer-Watson (A. Bowyer / D. F. Watson, Computer
//      Journal 24(2), 1981) inside a super-triangle of circumradius
//      5 * bbox_diagonal.  The insertion order is the input order (outline
//      vertices, then hole vertices hole by hole, then new segment endpoints,
//      then `points`), so the vertex numbering of the result is a pure function
//      of the input.  Point location is a deterministic visibility walk seeded
//      at the previously created triangle, with a linear-scan fallback if the
//      walk stalls.  The cavity grows across non-constrained edges whose
//      neighbour fails the in-circle test AND whose two far edges keep the
//      inserted point strictly to their left ("star-shape guard"): that is what
//      keeps the cavity valid in the presence of constraints and forbids the
//      creation of inverted or zero-area triangles.
//
//   2. CONSTRAINT RECOVERY.  Anglada-style edge flipping (M. V. Anglada, "An
//      improved incremental algorithm for constructing restricted Delaunay
//      triangulations", Computers & Graphics 21(2), 1997): the edges crossed by
//      the missing constraint are enumerated by a straight walk, and the first
//      strictly-convex crossing quadrilateral is flipped; the walk is redone
//      after each flip (O(k^2) per constraint, k = crossed triangles).  Three
//      special cases are handled before any flip:
//        * an existing vertex lying ON the constraint -> the constraint is
//          split there and both halves are recovered recursively;
//        * a crossed edge that is itself constrained -> the intersection point
//          is inserted and BOTH constraints are split there;
//        * no flippable crossing edge at all -> the constraint is split at its
//          midpoint and both halves are recovered (Shewchuk's fallback).
//      Split points are always evaluated as a + t*(b-a) on the ORIGINAL segment
//      endpoints, so a chain of splits stays collinear with the input segment
//      to ~1e-16 * |ab| instead of accumulating midpoint round-off.
//
//   3. DOMAIN CLASSIFICATION.  Breadth-first flood from the triangles touching
//      the super-triangle (outside by construction), toggling the inside flag
//      each time a BOUNDARY constraint (outline or hole edge) is crossed and
//      keeping it across internal constraints and free edges.  Parity, not
//      reachability, so arbitrarily nested holes work.
//
//   4. RUPPERT REFINEMENT (J. Ruppert, "A Delaunay Refinement Algorithm for
//      Quality 2-Dimensional Mesh Generation", J. Algorithms 18(3), 1995;
//      J. R. Shewchuk, "Delaunay Refinement Algorithms for Triangular Mesh
//      Generation", Comput. Geom. 22(1-3), 2002, pp. 21-74 -- the canonical
//      reference for the encroachment rules and termination bounds used here).
//      a. every constraint sub-segment longer than h is pre-split into
//         ceil(L/h) equal pieces, so the boundary discretisation honours the
//         target edge length independently of the interior;
//      b. loop: split every ENCROACHED sub-segment first (a sub-segment is
//         encroached when the apex of an adjacent IN-DOMAIN triangle lies
//         strictly inside its diametral circle, i.e. (a-p).(b-p) < 0), then
//         take the first bad triangle and insert its circumcenter -- unless the
//         circumcenter itself encroaches upon sub-segments (those are split and
//         the circumcenter discarded) or is not visible from the bad triangle
//         (the blocking sub-segment is split).
//
// QUALITY CRITERION.  A triangle of the domain is "bad" when either
//      * SIZE:   R > h / sqrt(3), i.e. its circumradius exceeds that of the
//        equilateral triangle of side h.  Since every edge of a triangle is at
//        most 2R, this bounds every edge by 2h/sqrt(3) ~ 1.155 h while
//        accepting an exactly equilateral triangle of side h; or
//      * SHAPE:  R / shortest_edge > 1 / (2 sin(min_angle)), the classical
//        radius-edge ratio.  R/e = 1/(2 sin(theta_min)) exactly, so this is the
//        minimum-angle policy written without a trigonometric evaluation.
//
// TERMINATION AND LIMITS.  Ruppert's proof of termination holds for a radius-
// edge bound of sqrt(2), i.e. min_angle <= 20.7 deg.  Between 20.7 and 34 deg
// termination is empirical (it is what every practical mesher ships); above
// 34 deg it is known to diverge, so min_angle_deg is required to lie in
// [20, 34].  Two further, unavoidable exceptions:
//      * INPUT ANGLES.  Two constraints meeting at an angle smaller than
//        min_angle_deg can never be honoured -- no triangle in that corner can
//        have a larger smallest angle.  Ruppert's algorithm may then cycle on
//        shrinking sub-segments; the size budget stops it.  The angle policy is
//        therefore guaranteed only for triangles that are NOT incident to an
//        input vertex whose input angle is itself below the policy.
//      * SIZE BUDGET.  Refinement stops when the vertex count exceeds
//        200 * the expected node count (domain_area / (sqrt(3)/4 h^2) / 2,
//        floored at 16).  Reaching the budget, or leaving a bad triangle or an
//        encroached sub-segment behind, sets refinement_complete = false, which
//        is a keta::fail unless the caller passed allow_partial.
//
// TOLERANCE POLICY (predicates).  orient2d and incircle are evaluated in double
// with the classical static error filter (|det| > c * sum of the magnitudes of
// the products).  When the filter fails, the same expression is re-evaluated in
// `long double` and accepted when it exceeds 8 * eps_ld times the same
// magnitude sum; otherwise the predicate returns an exact 0 (collinear /
// cocircular).  On x86-64 `long double` carries 64 mantissa bits (11 more than
// double) and on aarch64 113 bits, so the ladder buys 3 to 18 decimal digits;
// where `long double == double` it degenerates gracefully to "declare
// degenerate", never to a wrong sign.  Every consumer treats 0 as a hard
// degeneracy: a cocircular quadrilateral is never flipped, a collinear triple
// never grows the cavity, so a misclassified near-degeneracy costs mesh
// quality, never mesh validity.  Collinear runs, cocircular quads and
// grid-aligned point sets are handled.
//   Two length tolerances, both relative to the bounding-box diagonal:
//      * kDupRel = 1e-12 -- two input vertices closer than this are a typed
//        failure (contract); a COMPUTED point closer than this to an existing
//        vertex is snapped to it instead of inserted;
//      * the same value serves as the "on the boundary" band of the
//        containment checks.
//
// KNOWN LIMITATIONS.
//   * The validation pass is O(n^2) in the number of boundary edges (with a
//     bounding-box pre-filter) and the encroachment rescan is O(S) per inserted
//     vertex.  Both are fine for the mesh sizes the kernel produces
//     (10^3..10^5 triangles) but neither is asymptotically optimal.
//   * locate() rebuilds a vertex->triangle incidence table on every call
//     (O(N+T)) because the frozen TriMesh carries no adjacency; the walk itself
//     is O(sqrt(T)).  A persistent locator belongs to a future header revision.
//   * Curved boundaries must be pre-discretised by the caller: the mesher never
//     moves an input vertex and never adds one on a curve.
//   * Dead triangle slots are never recycled, so peak memory is ~5x the final
//     triangle count.
#include "keta/mesh2d.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <deque>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>

namespace keta::mesh2 {
namespace {

// ---------------------------------------------------------------------------
// Geometric predicates (see the tolerance policy in the module header).
// ---------------------------------------------------------------------------

constexpr double kOrientFilter = 3.3e-16;    // ~3 eps, the classical static filter
constexpr double kInCircleFilter = 1.1e-14;  // ~48 eps, Shewchuk's iccerrboundA order
constexpr double kDupRel = 1e-12;            // duplicate / on-boundary band, relative to bbox

inline long double ld_bound() { return 8.0L * std::numeric_limits<long double>::epsilon(); }
inline long double labs_(long double v) { return v < 0 ? -v : v; }

// > 0 if c is strictly left of a->b, < 0 if right, 0 if collinear.
double orient2d(const Vec2& a, const Vec2& b, const Vec2& c) {
  const double acx = a.x() - c.x(), acy = a.y() - c.y();
  const double bcx = b.x() - c.x(), bcy = b.y() - c.y();
  const double dl = acx * bcy, dr = acy * bcx;
  const double det = dl - dr;
  const double mag = std::fabs(dl) + std::fabs(dr);
  if (mag == 0.0) return 0.0;
  if (std::fabs(det) > kOrientFilter * mag) return det;
  using L = long double;
  const L lacx = (L)a.x() - (L)c.x(), lacy = (L)a.y() - (L)c.y();
  const L lbcx = (L)b.x() - (L)c.x(), lbcy = (L)b.y() - (L)c.y();
  const L ldl = lacx * lbcy, ldr = lacy * lbcx;
  const L ldet = ldl - ldr;
  const L lmag = labs_(ldl) + labs_(ldr);
  if (labs_(ldet) > ld_bound() * lmag) return ldet > 0 ? 1.0 : -1.0;
  return 0.0;
}

// > 0 if d is strictly inside the circumcircle of the CCW triangle abc.
double incircle(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
  const double adx = a.x() - d.x(), ady = a.y() - d.y();
  const double bdx = b.x() - d.x(), bdy = b.y() - d.y();
  const double cdx = c.x() - d.x(), cdy = c.y() - d.y();
  const double bc = bdx * cdy, cb = cdx * bdy;
  const double ca = cdx * ady, ac = adx * cdy;
  const double ab = adx * bdy, ba = bdx * ady;
  const double al = adx * adx + ady * ady;
  const double bl = bdx * bdx + bdy * bdy;
  const double cl = cdx * cdx + cdy * cdy;
  const double det = al * (bc - cb) + bl * (ca - ac) + cl * (ab - ba);
  const double perm = (std::fabs(bc) + std::fabs(cb)) * al + (std::fabs(ca) + std::fabs(ac)) * bl +
                      (std::fabs(ab) + std::fabs(ba)) * cl;
  if (perm == 0.0) return 0.0;
  if (std::fabs(det) > kInCircleFilter * perm) return det;
  using L = long double;
  const L ladx = (L)a.x() - (L)d.x(), lady = (L)a.y() - (L)d.y();
  const L lbdx = (L)b.x() - (L)d.x(), lbdy = (L)b.y() - (L)d.y();
  const L lcdx = (L)c.x() - (L)d.x(), lcdy = (L)c.y() - (L)d.y();
  const L lbc = lbdx * lcdy, lcb = lcdx * lbdy;
  const L lca = lcdx * lady, lac = ladx * lcdy;
  const L lab = ladx * lbdy, lba = lbdx * lady;
  const L lal = ladx * ladx + lady * lady;
  const L lbl = lbdx * lbdx + lbdy * lbdy;
  const L lcl = lcdx * lcdx + lcdy * lcdy;
  const L ldet = lal * (lbc - lcb) + lbl * (lca - lac) + lcl * (lab - lba);
  const L lperm = (labs_(lbc) + labs_(lcb)) * lal + (labs_(lca) + labs_(lac)) * lbl +
                  (labs_(lab) + labs_(lba)) * lcl;
  if (labs_(ldet) > ld_bound() * lperm) return ldet > 0 ? 1.0 : -1.0;
  return 0.0;
}

Vec2 circumcenter(const Vec2& a, const Vec2& b, const Vec2& c) {
  const double bx = b.x() - a.x(), by = b.y() - a.y();
  const double cx = c.x() - a.x(), cy = c.y() - a.y();
  const double d = 2.0 * (bx * cy - by * cx);
  const double bl = bx * bx + by * by, cl = cx * cx + cy * cy;
  return Vec2(a.x() + (cy * bl - by * cl) / d, a.y() + (bx * cl - cx * bl) / d);
}

// ---------------------------------------------------------------------------
// Polygon helpers used by the validation pass.
// ---------------------------------------------------------------------------

double poly_signed_area(const std::vector<Vec2>& p) {
  double s = 0;
  const std::size_t n = p.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = p[i];
    const Vec2& b = p[(i + 1) % n];
    s += a.x() * b.y() - b.x() * a.y();
  }
  return 0.5 * s;
}

double dist_point_seg(const Vec2& a, const Vec2& b, const Vec2& q) {
  const Vec2 ab = b - a;
  const double l2 = ab.squaredNorm();
  if (l2 == 0.0) return (q - a).norm();
  double t = (q - a).dot(ab) / l2;
  t = std::min(1.0, std::max(0.0, t));
  return (q - (a + ab * t)).norm();
}

// +1 strictly inside, 0 within `tol` of the boundary, -1 strictly outside.
int point_vs_poly(const std::vector<Vec2>& poly, const Vec2& q, double tol) {
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i)
    if (dist_point_seg(poly[i], poly[(i + 1) % n], q) <= tol) return 0;
  bool in = false;
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = poly[i];
    const Vec2& b = poly[(i + 1) % n];
    if ((a.y() > q.y()) != (b.y() > q.y())) {
      const double x = a.x() + (q.y() - a.y()) / (b.y() - a.y()) * (b.x() - a.x());
      if (q.x() < x) in = !in;
    }
  }
  return in ? 1 : -1;
}

bool strictly_between(const Vec2& a, const Vec2& b, const Vec2& q) {
  if (orient2d(a, b, q) != 0.0) return false;
  return (q - a).dot(b - a) > 0.0 && (q - b).dot(a - b) > 0.0;
}

// True when the two closed segments share more than a common endpoint: a proper
// crossing, an endpoint strictly inside the other, or a duplicated edge.
bool segments_conflict(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
  if (std::min(a.x(), b.x()) > std::max(c.x(), d.x())) return false;
  if (std::max(a.x(), b.x()) < std::min(c.x(), d.x())) return false;
  if (std::min(a.y(), b.y()) > std::max(c.y(), d.y())) return false;
  if (std::max(a.y(), b.y()) < std::min(c.y(), d.y())) return false;
  const double o1 = orient2d(a, b, c), o2 = orient2d(a, b, d);
  const double o3 = orient2d(c, d, a), o4 = orient2d(c, d, b);
  if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0))
    return true;
  if (strictly_between(a, b, c) || strictly_between(a, b, d)) return true;
  if (strictly_between(c, d, a) || strictly_between(c, d, b)) return true;
  if (o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0) {
    const bool same = ((a - c).squaredNorm() == 0 && (b - d).squaredNorm() == 0) ||
                      ((a - d).squaredNorm() == 0 && (b - c).squaredNorm() == 0);
    if (same) return true;
  }
  return false;
}

// Weaker relation used for "does this internal segment leave the domain?": a
// PROPER crossing or a collinear overlap of positive length.  Merely touching
// (an endpoint sitting on the boundary, which is the normal way to anchor a rib
// or a line support) is legal and must not be rejected.
bool segments_properly_cross(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
  if (std::min(a.x(), b.x()) > std::max(c.x(), d.x())) return false;
  if (std::max(a.x(), b.x()) < std::min(c.x(), d.x())) return false;
  if (std::min(a.y(), b.y()) > std::max(c.y(), d.y())) return false;
  if (std::max(a.y(), b.y()) < std::min(c.y(), d.y())) return false;
  const double o1 = orient2d(a, b, c), o2 = orient2d(a, b, d);
  const double o3 = orient2d(c, d, a), o4 = orient2d(c, d, b);
  if (o1 != 0 && o2 != 0 && o3 != 0 && o4 != 0 && (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0))
    return true;
  if (o1 == 0 && o2 == 0) {  // collinear: reject an overlap of positive length
    const Vec2 dir = b - a;
    const double ta = 0.0, tb = dir.dot(dir);
    const double tc = (c - a).dot(dir), td = (d - a).dot(dir);
    const double l = std::max(std::min(ta, tb), std::min(tc, td));
    const double u = std::min(std::max(ta, tb), std::max(tc, td));
    if (u - l > 0) return true;
  }
  return false;
}

// ---------------------------------------------------------------------------
// Triangulation data structures.
//
// Triangles are CCW.  nb[i] is the triangle across the edge OPPOSITE v[i], that
// is the edge (v[(i+1)%3], v[(i+2)%3]).  con[i] carries that edge's flags:
// bit 0 = constrained, bit 1 = domain boundary (outline or hole edge).
// ---------------------------------------------------------------------------

constexpr unsigned char kConstrained = 1;
constexpr unsigned char kBoundary = 2;

struct Tri {
  std::array<int, 3> v{{-1, -1, -1}};
  std::array<int, 3> nb{{-1, -1, -1}};
  std::array<unsigned char, 3> con{{0, 0, 0}};
  bool alive = true;
  bool inside = false;
};

// One original (unsplit) constraint: a polygon edge or an input segment.
struct OrigSeg {
  Vec2 a, b;
  bool boundary = false;
};

// One sub-segment of an original constraint, covering [t0, t1] of it.
struct SubSeg {
  int a = -1, b = -1;
  int orig = -1;
  double t0 = 0, t1 = 1;
  bool alive = true;
};

struct Mesher {
  std::vector<Vec2> P;    // 0,1,2 = super-triangle vertices
  std::vector<int> vtri;  // vertex -> one live incident triangle
  std::vector<Tri> T;
  std::vector<OrigSeg> orig;
  std::vector<SubSeg> segs;

  std::vector<int> stamp;  // cavity marker, compared against `gen`
  int gen = 0;

  double eps = 0;        // absolute length tolerance = kDupRel * bbox diagonal
  double h = 0;          // target edge length
  double max_ratio = 1;  // circumradius / shortest-edge bound
  double size_r = 0;     // circumradius bound = h / sqrt(3)
  int vert_budget = 0;

  // -- small accessors ------------------------------------------------------
  int local(int t, int u) const {
    for (int i = 0; i < 3; ++i)
      if (T[t].v[i] == u) return i;
    return -1;
  }
  int mate(int t, int i) const {
    const int n = T[t].nb[i];
    if (n < 0) return -1;
    for (int j = 0; j < 3; ++j)
      if (T[n].nb[j] == t) return j;
    return -1;
  }
  const Vec2& pt(int t, int i) const { return P[static_cast<std::size_t>(T[t].v[i])]; }

  int add_tri(int a, int b, int c) {
    const int idx = static_cast<int>(T.size());
    Tri t;
    t.v = {{a, b, c}};
    T.push_back(t);
    stamp.push_back(0);
    vtri[static_cast<std::size_t>(a)] = idx;
    vtri[static_cast<std::size_t>(b)] = idx;
    vtri[static_cast<std::size_t>(c)] = idx;
    return idx;
  }

  void link(int t, int i, int n, int j) {
    if (t >= 0) T[t].nb[i] = n;
    if (n >= 0 && j >= 0) T[n].nb[j] = t;
  }

  int find_edge_index(int t, int u, int w) const {
    for (int i = 0; i < 3; ++i) {
      const int a = T[t].v[(i + 1) % 3], b = T[t].v[(i + 2) % 3];
      if ((a == u && b == w) || (a == w && b == u)) return i;
    }
    return -1;
  }
  int local_apex(int t, int u, int w) const {
    for (int i = 0; i < 3; ++i)
      if (T[t].v[i] != u && T[t].v[i] != w) return i;
    return -1;
  }

  // Triangles incident to vertex u, as (triangle, local index of u).
  void ring(int u, std::vector<std::array<int, 2>>& out) const {
    out.clear();
    const int t0 = vtri[static_cast<std::size_t>(u)];
    if (t0 < 0) return;
    const int k0 = local(t0, u);
    if (k0 < 0) return;
    int t = t0, k = k0;
    for (int guard = 0; guard < (1 << 22); ++guard) {
      out.push_back({{t, k}});
      const int nx = T[t].nb[(k + 2) % 3];
      if (nx == t0) return;
      if (nx < 0) {
        // Open fan: sweep the other way from the seed.
        int t2 = T[t0].nb[(k0 + 1) % 3];
        for (int g2 = 0; t2 >= 0 && g2 < (1 << 22); ++g2) {
          const int k2 = local(t2, u);
          if (k2 < 0) return;
          out.push_back({{t2, k2}});
          t2 = T[t2].nb[(k2 + 1) % 3];
        }
        return;
      }
      t = nx;
      k = local(t, u);
      if (k < 0) return;
    }
  }

  // (triangle, edge index) of the edge {a,b}; {-1,-1} when it does not exist.
  std::array<int, 2> find_edge(int a, int b) const {
    std::vector<std::array<int, 2>> r;
    ring(a, r);
    for (const auto& tk : r) {
      const int t = tk[0], k = tk[1];
      if (T[t].v[(k + 1) % 3] == b) return {{t, (k + 2) % 3}};
      if (T[t].v[(k + 2) % 3] == b) return {{t, (k + 1) % 3}};
    }
    return {{-1, -1}};
  }

  void mark_edge(int a, int b, unsigned char mask) {
    const auto e = find_edge(a, b);
    require(e[0] >= 0, "mesh2d: internal error, constrained edge vanished");
    T[e[0]].con[e[1]] |= mask;
    const int n = T[e[0]].nb[e[1]];
    if (n >= 0) T[n].con[mate(e[0], e[1])] |= mask;
  }

  // -- point classification / location -------------------------------------
  // -1 outside, 0 strictly inside, 1+i on edge i, 10+k at vertex k.
  int classify(int t, const Vec2& p) const {
    int zeros = 0, zi = -1, zj = -1;
    for (int i = 0; i < 3; ++i) {
      const double o = orient2d(pt(t, (i + 1) % 3), pt(t, (i + 2) % 3), p);
      if (o < 0) return -1;
      if (o == 0) {
        ++zeros;
        if (zi < 0)
          zi = i;
        else
          zj = i;
      }
    }
    if (zeros == 0) return 0;
    if (zeros == 1) return 1 + zi;
    if (zeros == 2) return 10 + (3 - zi - zj);
    return -1;  // fully degenerate triangle
  }

  // Deterministic visibility walk. Returns the containing triangle, or -1.
  // With `respect_con`, a constrained edge blocking the way is reported through
  // `blocked` and the walk returns -1.
  int walk(int start, const Vec2& target, bool respect_con, std::array<int, 2>* blocked) const {
    if (blocked) *blocked = {{-1, -1}};
    int t = (start >= 0 && T[start].alive) ? start : -1;
    if (t < 0)
      for (std::size_t k = 0; k < T.size(); ++k)
        if (T[k].alive) {
          t = static_cast<int>(k);
          break;
        }
    if (t < 0) return -1;
    int prev = -1;
    const int limit = 8 * static_cast<int>(T.size()) + 64;
    for (int it = 0; it < limit; ++it) {
      bool moved = false;
      for (int k = 0; k < 3; ++k) {
        const int i = (it + k) % 3;
        if (orient2d(pt(t, (i + 1) % 3), pt(t, (i + 2) % 3), target) >= 0) continue;
        if (respect_con && (T[t].con[i] & kConstrained)) {
          if (blocked) *blocked = {{t, i}};
          return -1;
        }
        const int nx = T[t].nb[i];
        if (nx < 0) return -1;
        if (nx == prev) continue;
        prev = t;
        t = nx;
        moved = true;
        break;
      }
      if (!moved) {
        if (classify(t, target) >= 0) return t;
        break;
      }
    }
    if (respect_con) return -1;
    for (std::size_t k = 0; k < T.size(); ++k)
      if (T[k].alive && classify(static_cast<int>(k), target) >= 0) return static_cast<int>(k);
    return -1;
  }

  // -- Bowyer-Watson insertion ---------------------------------------------
  // `seed_edge` >= 0 declares that p lies on edge `seed_edge` of `t0`.
  int insert_vertex(const Vec2& p, int t0, int seed_edge) {
    const int np = static_cast<int>(P.size());
    P.push_back(p);
    vtri.push_back(-1);

    ++gen;
    stamp.resize(T.size(), 0);
    std::vector<int> cav;
    std::deque<int> q;
    auto push = [&](int t) {
      stamp[static_cast<std::size_t>(t)] = gen;
      cav.push_back(t);
      q.push_back(t);
    };
    push(t0);
    if (seed_edge >= 0) {
      const int n = T[t0].nb[seed_edge];
      if (n >= 0) push(n);
    }
    while (!q.empty()) {
      const int t = q.front();
      q.pop_front();
      for (int i = 0; i < 3; ++i) {
        const int n = T[t].nb[i];
        if (n < 0 || stamp[static_cast<std::size_t>(n)] == gen) continue;
        if (T[t].con[i] & kConstrained) continue;
        if (incircle(pt(n, 0), pt(n, 1), pt(n, 2), p) <= 0) continue;
        const int j = mate(t, i);
        if (j < 0) continue;
        bool ok = true;
        for (int k = 1; k <= 2 && ok; ++k) {
          const int e = (j + k) % 3;
          if (orient2d(pt(n, (e + 1) % 3), pt(n, (e + 2) % 3), p) <= 0) ok = false;
        }
        if (!ok) continue;
        push(n);
      }
    }

    struct RingEdge {
      int u, w, outer, outer_i;
      unsigned char con;
      bool inside;
    };
    std::vector<RingEdge> ring_e;
    for (const int t : cav)
      for (int i = 0; i < 3; ++i) {
        const int n = T[t].nb[i];
        if (n >= 0 && stamp[static_cast<std::size_t>(n)] == gen) continue;
        ring_e.push_back({T[t].v[(i + 1) % 3], T[t].v[(i + 2) % 3], n, n >= 0 ? mate(t, i) : -1,
                          T[t].con[i], T[t].inside});
      }
    for (const int t : cav) T[t].alive = false;

    std::vector<int> nt(ring_e.size(), -1);
    for (std::size_t r = 0; r < ring_e.size(); ++r) {
      const RingEdge e = ring_e[r];
      const int t = add_tri(e.u, e.w, np);
      T[t].inside = e.inside;
      T[t].con[2] = e.con;
      T[t].nb[2] = e.outer;
      if (e.outer >= 0) T[e.outer].nb[e.outer_i] = t;
      nt[r] = t;
    }
    for (std::size_t r = 0; r < ring_e.size(); ++r) {
      int a = -1, b = -1;
      for (std::size_t s = 0; s < ring_e.size(); ++s) {
        if (ring_e[s].u == ring_e[r].w) a = nt[s];
        if (ring_e[s].w == ring_e[r].u) b = nt[s];
      }
      require(a >= 0 && b >= 0, "mesh2d: internal error, insertion cavity is not a simple ring");
      T[nt[r]].nb[0] = a;
      T[nt[r]].nb[1] = b;
    }
    return np;
  }

  // Insert p, reusing an existing vertex when p falls on (or within eps of) one.
  int insert_point(const Vec2& p, int hint) {
    const int t = walk(hint, p, false, nullptr);
    require(t >= 0, "mesh2d: internal error, point location failed");
    const int c = classify(t, p);
    require(c >= 0, "mesh2d: internal error, located triangle does not contain the point");
    if (c >= 10) return T[t].v[c - 10];
    for (int i = 0; i < 3; ++i)
      if ((pt(t, i) - p).norm() <= eps) return T[t].v[i];
    return insert_vertex(p, t, c == 0 ? -1 : c - 1);
  }

  // -- edge flip ------------------------------------------------------------
  // Flip edge `i` of triangle `t`; the caller guarantees it is flippable and
  // unconstrained.  Triangle slots are reused so neighbour back-pointers of the
  // two untouched sides stay valid.
  void flip(int t, int i) {
    const int t2 = T[t].nb[i];
    const int j = mate(t, i);
    const int p = T[t].v[i], u = T[t].v[(i + 1) % 3], w = T[t].v[(i + 2) % 3];
    const int qv = T[t2].v[j];
    const int A = T[t].nb[(i + 1) % 3];   // across (w,p)
    const int B = T[t].nb[(i + 2) % 3];   // across (p,u)
    const int C = T[t2].nb[(j + 1) % 3];  // across (u,q)
    const int D = T[t2].nb[(j + 2) % 3];  // across (q,w)
    const unsigned char cA = T[t].con[(i + 1) % 3], cB = T[t].con[(i + 2) % 3];
    const unsigned char cC = T[t2].con[(j + 1) % 3], cD = T[t2].con[(j + 2) % 3];
    const int Ai = mate(t, (i + 1) % 3);
    const int Bi = mate(t, (i + 2) % 3);
    const int Ci = mate(t2, (j + 1) % 3);
    const int Di = mate(t2, (j + 2) % 3);
    const bool ins = T[t].inside;

    T[t].v = {{p, u, qv}};
    T[t].con = {{cC, 0, cB}};
    T[t2].v = {{p, qv, w}};
    T[t2].con = {{cD, cA, 0}};
    T[t].inside = ins;
    T[t2].inside = ins;
    link(t, 0, C, Ci);
    link(t, 1, t2, 2);
    link(t, 2, B, Bi);
    link(t2, 0, D, Di);
    link(t2, 1, A, Ai);
    vtri[static_cast<std::size_t>(p)] = t;
    vtri[static_cast<std::size_t>(u)] = t;
    vtri[static_cast<std::size_t>(qv)] = t;
    vtri[static_cast<std::size_t>(w)] = t2;
  }

  // -- constraints ----------------------------------------------------------
  double param_on(int oi, const Vec2& p) const {
    const Vec2 d = orig[static_cast<std::size_t>(oi)].b - orig[static_cast<std::size_t>(oi)].a;
    const double l2 = d.squaredNorm();
    return l2 > 0 ? (p - orig[static_cast<std::size_t>(oi)].a).dot(d) / l2 : 0.0;
  }

  unsigned char mask_of(int oi) const {
    return static_cast<unsigned char>(
        kConstrained | (orig[static_cast<std::size_t>(oi)].boundary ? kBoundary : 0));
  }

  void record_seg(int a, int b, int oi, double t0, double t1) {
    segs.push_back({a, b, oi, t0, t1, true});
    mark_edge(a, b, mask_of(oi));
  }

  int seg_index(int a, int b) const {
    for (std::size_t i = 0; i < segs.size(); ++i) {
      if (!segs[i].alive) continue;
      if ((segs[i].a == a && segs[i].b == b) || (segs[i].a == b && segs[i].b == a))
        return static_cast<int>(i);
    }
    return -1;
  }

  static Vec2 line_intersection(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d) {
    const Vec2 r = b - a, s = d - c;
    const double den = r.x() * s.y() - r.y() * s.x();
    if (den == 0.0) return 0.5 * (c + d);
    const double t = ((c - a).x() * s.y() - (c - a).y() * s.x()) / den;
    return a + r * t;
  }

  // Split sub-segment `si` at parameter `t` of its ORIGINAL segment; returns
  // the new vertex, or -1 when the split point would collapse onto one of the
  // endpoints (the floating-point floor of segment splitting, reached only in
  // the small-input-angle cycle -- a soft refusal, never a corrupt mesh).
  int split_subseg_at(int si, double t) {
    const SubSeg s = segs[static_cast<std::size_t>(si)];
    const OrigSeg o = orig[static_cast<std::size_t>(s.orig)];
    double tm = t;
    if (!(tm > s.t0 && tm < s.t1)) tm = 0.5 * (s.t0 + s.t1);
    const Vec2 m = o.a + (o.b - o.a) * tm;
    if ((m - P[static_cast<std::size_t>(s.a)]).squaredNorm() == 0 ||
        (m - P[static_cast<std::size_t>(s.b)]).squaredNorm() == 0)
      return -1;
    const auto e = find_edge(s.a, s.b);
    require(e[0] >= 0, "mesh2d: internal error, sub-segment edge missing");
    const int nv = insert_vertex(m, e[0], e[1]);
    segs[static_cast<std::size_t>(si)].alive = false;
    const unsigned char mask = mask_of(s.orig);
    segs.push_back({s.a, nv, s.orig, s.t0, tm, true});
    mark_edge(s.a, nv, mask);
    segs.push_back({nv, s.b, s.orig, tm, s.t1, true});
    mark_edge(nv, s.b, mask);
    return nv;
  }

  // Edges crossed by segment (a,b) starting at vertex a, as (u, w, triangle).
  // Sets *split_v when an existing vertex lies strictly on the segment.
  bool crossings(int a, int b, std::vector<std::array<int, 3>>& out, int* split_v) const {
    out.clear();
    *split_v = -1;
    const Vec2 A = P[static_cast<std::size_t>(a)];
    const Vec2 B = P[static_cast<std::size_t>(b)];
    std::vector<std::array<int, 2>> r;
    ring(a, r);
    int t = -1, u = -1, w = -1;
    for (const auto& tk : r) {
      const int tt = tk[0], k = tk[1];
      const int uu = T[tt].v[(k + 1) % 3], ww = T[tt].v[(k + 2) % 3];
      const Vec2& U = P[static_cast<std::size_t>(uu)];
      const Vec2& W = P[static_cast<std::size_t>(ww)];
      const double ou = orient2d(A, B, U), ow = orient2d(A, B, W);
      if (ou == 0 && strictly_between(A, B, U)) {
        *split_v = uu;
        return false;
      }
      if (ow == 0 && strictly_between(A, B, W)) {
        *split_v = ww;
        return false;
      }
      if (ou < 0 && ow > 0) {
        t = tt;
        u = uu;
        w = ww;
        break;
      }
    }
    if (t < 0) return false;
    const int limit = 4 * static_cast<int>(T.size()) + 64;
    for (int it = 0; it < limit; ++it) {
      const int i = find_edge_index(t, u, w);
      if (i < 0) return false;
      out.push_back({{u, w, t}});
      const int n = T[t].nb[i];
      if (n < 0) return false;
      const int j = local_apex(n, u, w);
      if (j < 0) return false;
      const int x = T[n].v[j];
      if (x == b) return true;
      const double d = orient2d(A, B, P[static_cast<std::size_t>(x)]);
      if (d == 0) {
        if (strictly_between(A, B, P[static_cast<std::size_t>(x)])) *split_v = x;
        return false;
      }
      t = n;
      if (d < 0)
        u = x;
      else
        w = x;
    }
    return false;
  }

  void insert_constraint(int a, int b, int oi, double t0, double t1, int depth) {
    require(depth < 8192, "mesh2d: constraint recovery recursed too deep");
    if (a == b) return;
    std::vector<std::array<int, 3>> cross;
    for (int guard = 0; guard < (1 << 22); ++guard) {
      if (find_edge(a, b)[0] >= 0) {
        record_seg(a, b, oi, t0, t1);
        return;
      }
      int split_v = -1;
      crossings(a, b, cross, &split_v);
      if (split_v >= 0) {
        const double tm = param_on(oi, P[static_cast<std::size_t>(split_v)]);
        insert_constraint(a, split_v, oi, t0, tm, depth + 1);
        insert_constraint(split_v, b, oi, tm, t1, depth + 1);
        return;
      }

      // A crossed edge that is itself constrained: cut both at the crossing.
      bool cut = false;
      for (const auto& c : cross) {
        const int t = c[2];
        if (!T[t].alive) continue;
        const int i = find_edge_index(t, c[0], c[1]);
        if (i < 0 || !(T[t].con[i] & kConstrained)) continue;
        const int si = seg_index(c[0], c[1]);
        require(si >= 0, "mesh2d: internal error, constrained edge has no sub-segment");
        const Vec2 x = line_intersection(P[static_cast<std::size_t>(a)],
                                         P[static_cast<std::size_t>(b)],
                                         P[static_cast<std::size_t>(c[0])],
                                         P[static_cast<std::size_t>(c[1])]);
        const int so = segs[static_cast<std::size_t>(si)].orig;
        const int nv = split_subseg_at(si, param_on(so, x));
        require(nv >= 0, "mesh2d: crossing constraints intersect at an existing vertex");
        const double tm = param_on(oi, P[static_cast<std::size_t>(nv)]);
        insert_constraint(a, nv, oi, t0, tm, depth + 1);
        insert_constraint(nv, b, oi, tm, t1, depth + 1);
        cut = true;
        break;
      }
      if (cut) return;

      bool flipped = false;
      for (const auto& c : cross) {
        const int t = c[2];
        if (!T[t].alive) continue;
        const int i = find_edge_index(t, c[0], c[1]);
        if (i < 0) continue;
        const int n = T[t].nb[i];
        if (n < 0) continue;
        const int j = mate(t, i);
        if (j < 0) continue;
        const Vec2 p = pt(t, i);
        const Vec2 U = pt(t, (i + 1) % 3);
        const Vec2 W = pt(t, (i + 2) % 3);
        const Vec2 Q = pt(n, j);
        if (orient2d(p, U, Q) <= 0 || orient2d(Q, W, p) <= 0) continue;
        flip(t, i);
        flipped = true;
        break;
      }
      if (flipped) continue;

      // Nothing recoverable by flipping: split the constraint (Shewchuk).
      const double tm = 0.5 * (t0 + t1);
      const Vec2 m = orig[static_cast<std::size_t>(oi)].a +
                     (orig[static_cast<std::size_t>(oi)].b - orig[static_cast<std::size_t>(oi)].a) * tm;
      const int nv = insert_point(m, vtri[static_cast<std::size_t>(a)]);
      require(nv != a && nv != b, "mesh2d: constraint midpoint collapsed onto an endpoint");
      insert_constraint(a, nv, oi, t0, tm, depth + 1);
      insert_constraint(nv, b, oi, tm, t1, depth + 1);
      return;
    }
    fail("mesh2d: constraint recovery did not converge");
  }

  // -- domain classification ------------------------------------------------
  void flood_domain() {
    std::vector<char> seen(T.size(), 0);
    std::deque<int> q;
    for (std::size_t i = 0; i < T.size(); ++i) {
      if (!T[i].alive) continue;
      bool super = false;
      for (int k = 0; k < 3; ++k)
        if (T[i].v[k] < 3) super = true;
      if (!super) continue;
      T[i].inside = false;
      seen[i] = 1;
      q.push_back(static_cast<int>(i));
    }
    while (!q.empty()) {
      const int t = q.front();
      q.pop_front();
      for (int i = 0; i < 3; ++i) {
        const int n = T[t].nb[i];
        if (n < 0 || seen[static_cast<std::size_t>(n)]) continue;
        T[n].inside = (T[t].con[i] & kBoundary) ? !T[t].inside : T[t].inside;
        seen[static_cast<std::size_t>(n)] = 1;
        q.push_back(n);
      }
    }
  }

  // -- quality --------------------------------------------------------------
  double circumradius(int t) const {
    const Vec2& a = pt(t, 0);
    const Vec2& b = pt(t, 1);
    const Vec2& c = pt(t, 2);
    const double la = (b - c).norm(), lb = (c - a).norm(), lc = (a - b).norm();
    const double area2 = orient2d(a, b, c);
    if (area2 <= 0) return std::numeric_limits<double>::infinity();
    return la * lb * lc / (2.0 * area2);
  }

  bool tri_bad(int t) const {
    if (!T[t].alive || !T[t].inside) return false;
    const double R = circumradius(t);
    if (!std::isfinite(R)) return true;
    if (R > size_r) return true;
    const Vec2& a = pt(t, 0);
    const Vec2& b = pt(t, 1);
    const Vec2& c = pt(t, 2);
    const double shortest = std::min({(b - c).norm(), (c - a).norm(), (a - b).norm()});
    if (shortest <= 0) return true;
    return R / shortest > max_ratio;
  }

  bool seg_encroached(int si) const {
    const SubSeg& s = segs[static_cast<std::size_t>(si)];
    if (!s.alive) return false;
    const auto e = find_edge(s.a, s.b);
    if (e[0] < 0) return false;
    const Vec2& A = P[static_cast<std::size_t>(s.a)];
    const Vec2& B = P[static_cast<std::size_t>(s.b)];
    const int t = e[0], i = e[1];
    for (int side = 0; side < 2; ++side) {
      const int tt = side == 0 ? t : T[t].nb[i];
      if (tt < 0 || !T[tt].alive || !T[tt].inside) continue;
      const int ai = side == 0 ? i : mate(t, i);
      if (ai < 0) continue;
      const Vec2& X = pt(tt, ai);
      if ((A - X).dot(B - X) < 0) return true;
    }
    return false;
  }

  bool in_diametral(int si, const Vec2& p) const {
    const SubSeg& s = segs[static_cast<std::size_t>(si)];
    if (!s.alive) return false;
    return (P[static_cast<std::size_t>(s.a)] - p).dot(P[static_cast<std::size_t>(s.b)] - p) < 0;
  }

  // ---- diametral-disk grid of the sub-segments (Ruppert queries) ----------
  // The refinement asks, after every insertion, which sub-segments have the
  // new point inside their diametral circle: a scan of every sub-segment made
  // it O(insertions x sub-segments) = O(N^1.5). Each sub-segment is entered
  // (once, lazily: they are append-only, endpoints never move, a split one is
  // only marked dead) in the cells of a uniform grid its disk's bounding box
  // overlaps, inflated so that rounding can never put a point the exact test
  // accepts outside it; a query reads the point's cell, sorts the candidates
  // and applies the very same `in_diametral` — the same hits in the same
  // (index) order as the scan. Off-grid boxes or points use the scan.
  double gcell = 0;                                     // 0 = no grid
  std::unordered_map<long long, std::vector<int>> grid;
  std::vector<int> grid_big;                            // boxes too wide for cells
  std::size_t grid_done = 0;
  static constexpr double kGridSpan = 1e8;              // |coordinate / cell| bound
  static long long cell_key(long long ix, long long iy) { return ix * 4000000007LL + iy; }
  void grid_init(double cell) {
    gcell = cell > 0 && std::isfinite(cell) ? cell : 0;
    grid.clear();
    grid_big.clear();
    grid_done = 0;
  }
  void grid_sync() {
    for (; grid_done < segs.size(); ++grid_done) {
      const SubSeg& s = segs[grid_done];
      if (!s.alive) continue;  // a dead sub-segment never revives (the pre-split parents)
      const Vec2 A = P[static_cast<std::size_t>(s.a)], B = P[static_cast<std::size_t>(s.b)];
      const Vec2 c = 0.5 * (A + B);
      const double r = 0.5 * (B - A).norm() * (1 + 1e-6) + 1e-9 * (c.cwiseAbs().maxCoeff() + gcell);
      const double x0 = (c.x() - r) / gcell, x1 = (c.x() + r) / gcell;
      const double y0 = (c.y() - r) / gcell, y1 = (c.y() + r) / gcell;
      if (!(std::fabs(x0) < kGridSpan && std::fabs(x1) < kGridSpan && std::fabs(y0) < kGridSpan &&
            std::fabs(y1) < kGridSpan) ||
          (x1 - x0 + 1) * (y1 - y0 + 1) > 4096) {
        grid_big.push_back(static_cast<int>(grid_done));
        continue;
      }
      for (long long ix = static_cast<long long>(std::floor(x0)); ix <= static_cast<long long>(std::floor(x1)); ++ix)
        for (long long iy = static_cast<long long>(std::floor(y0)); iy <= static_cast<long long>(std::floor(y1)); ++iy)
          grid[cell_key(ix, iy)].push_back(static_cast<int>(grid_done));
    }
  }
  // the sub-segments i < below whose diametral circle holds p, ascending
  void diametral_hits(const Vec2& p, std::size_t below, std::vector<int>& out) {
    out.clear();
    const double px = p.x() / gcell, py = p.y() / gcell;
    if (gcell == 0 || !(std::fabs(px) < kGridSpan && std::fabs(py) < kGridSpan)) {
      for (std::size_t i = 0; i < below; ++i)
        if (in_diametral(static_cast<int>(i), p)) out.push_back(static_cast<int>(i));
      return;
    }
    grid_sync();
    auto it = grid.find(cell_key(static_cast<long long>(std::floor(px)), static_cast<long long>(std::floor(py))));
    if (it != grid.end())
      for (int i : it->second)
        if (static_cast<std::size_t>(i) < below && in_diametral(i, p)) out.push_back(i);
    for (int i : grid_big)
      if (static_cast<std::size_t>(i) < below && in_diametral(i, p)) out.push_back(i);
    std::sort(out.begin(), out.end());
  }
};

void check_finite(const std::vector<Vec2>& v, const char* what) {
  for (const auto& p : v)
    require(std::isfinite(p.x()) && std::isfinite(p.y()),
            std::string("mesh2d: non-finite coordinate in ") + what);
}

void check_simple(const std::vector<Vec2>& poly, const char* what) {
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = poly[i];
    const Vec2& b = poly[(i + 1) % n];
    require((b - a).norm() > 0, std::string("mesh2d: zero-length edge in ") + what);
    const Vec2& c = poly[(i + 2) % n];
    if (orient2d(a, b, c) == 0 && (a - b).dot(c - b) > 0)
      fail(std::string("mesh2d: self-intersecting ") + what + " (fold-back spike)");
    for (std::size_t j = i + 2; j < n; ++j) {
      if (i == 0 && j == n - 1) continue;  // adjacent through the wrap
      if (segments_conflict(a, b, poly[j], poly[(j + 1) % n]))
        fail(std::string("mesh2d: self-intersecting ") + what);
    }
  }
}

}  // namespace

double TriMesh::area() const {
  double s = 0;
  for (const auto& t : tris) {
    const Vec2& a = nodes[static_cast<std::size_t>(t[0])];
    const Vec2& b = nodes[static_cast<std::size_t>(t[1])];
    const Vec2& c = nodes[static_cast<std::size_t>(t[2])];
    s += 0.5 * ((b.x() - a.x()) * (c.y() - a.y()) - (b.y() - a.y()) * (c.x() - a.x()));
  }
  return s;
}

TriMesh triangulate(const MeshInput& in) {
  // ---- 0. validation ------------------------------------------------------
  require(in.outline.size() >= 3, "mesh2d: outline needs at least 3 vertices");
  check_finite(in.outline, "outline");
  for (const auto& hl : in.holes) {
    require(hl.size() >= 3, "mesh2d: a hole needs at least 3 vertices");
    check_finite(hl, "hole");
  }
  check_finite(in.points, "points");
  for (const auto& s : in.segments) {
    require(std::isfinite(s.first.x()) && std::isfinite(s.first.y()) &&
                std::isfinite(s.second.x()) && std::isfinite(s.second.y()),
            "mesh2d: non-finite coordinate in segments");
    require((s.second - s.first).norm() > 0, "mesh2d: zero-length constraint segment");
  }
  require(std::isfinite(in.target_h) && in.target_h >= 0,
          "mesh2d: target_h must be finite and >= 0");
  require(in.min_angle_deg >= 20.0 && in.min_angle_deg <= 34.0,
          "mesh2d: min_angle_deg must lie in [20, 34]");

  Vec2 lo = in.outline[0], hi = in.outline[0];
  auto grow = [&](const Vec2& p) {
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  };
  for (const auto& p : in.outline) grow(p);
  for (const auto& hl : in.holes)
    for (const auto& p : hl) grow(p);
  for (const auto& p : in.points) grow(p);
  for (const auto& s : in.segments) {
    grow(s.first);
    grow(s.second);
  }
  const double diag = (hi - lo).norm();
  require(diag > 0, "mesh2d: degenerate (zero-size) input");
  const double eps = kDupRel * diag;

  std::vector<Vec2> outline = in.outline;
  double a_out = poly_signed_area(outline);
  require(std::fabs(a_out) > 1e-24 * diag * diag, "mesh2d: zero-area outline");
  if (a_out < 0) {
    std::reverse(outline.begin(), outline.end());
    a_out = -a_out;
  }
  check_simple(outline, "outline");

  std::vector<std::vector<Vec2>> holes;
  double a_holes = 0;
  for (const auto& hin : in.holes) {
    std::vector<Vec2> hl = hin;
    double ah = poly_signed_area(hl);
    require(std::fabs(ah) > 1e-24 * diag * diag, "mesh2d: zero-area hole");
    if (ah > 0) {
      std::reverse(hl.begin(), hl.end());  // holes stored CW
      ah = -ah;
    }
    check_simple(hl, "hole");
    for (const auto& p : hl)
      require(point_vs_poly(outline, p, eps) >= 0, "mesh2d: hole vertex outside the outline");
    a_holes += -ah;
    holes.push_back(std::move(hl));
  }

  std::vector<std::array<Vec2, 2>> bedges;
  for (std::size_t i = 0; i < outline.size(); ++i)
    bedges.push_back({{outline[i], outline[(i + 1) % outline.size()]}});
  const std::size_t n_out = outline.size();
  for (const auto& hl : holes)
    for (std::size_t i = 0; i < hl.size(); ++i)
      bedges.push_back({{hl[i], hl[(i + 1) % hl.size()]}});

  for (std::size_t i = 0; i < bedges.size(); ++i)
    for (std::size_t j = i + 1; j < bedges.size(); ++j) {
      if (i < n_out && j < n_out) continue;  // covered by check_simple
      if (segments_conflict(bedges[i][0], bedges[i][1], bedges[j][0], bedges[j][1]))
        fail("mesh2d: boundary polygons intersect");
    }

  auto in_domain = [&](const Vec2& p) {
    if (point_vs_poly(outline, p, eps) < 0) return false;
    for (const auto& hl : holes)
      if (point_vs_poly(hl, p, eps) > 0) return false;
    return true;
  };
  for (const auto& p : in.points) require(in_domain(p), "mesh2d: extra point outside the domain");
  for (const auto& s : in.segments) {
    require(in_domain(s.first) && in_domain(s.second),
            "mesh2d: constraint segment endpoint outside the domain");
    require(in_domain(0.5 * (s.first + s.second)),
            "mesh2d: constraint segment leaves the domain");
    for (const auto& be : bedges)
      if (segments_properly_cross(s.first, s.second, be[0], be[1]))
        fail("mesh2d: constraint segment crosses a domain boundary");
  }

  // ---- 1. vertex list (deterministic order, duplicates rejected) ----------
  std::vector<Vec2> verts;
  auto push_input = [&](const Vec2& p, bool strict) -> int {
    for (std::size_t i = 0; i < verts.size(); ++i)
      if ((verts[i] - p).norm() <= eps) {
        require(!strict, "mesh2d: duplicate input points");
        return static_cast<int>(i);
      }
    verts.push_back(p);
    return static_cast<int>(verts.size()) - 1;
  };
  std::vector<int> out_id;
  for (const auto& p : outline) out_id.push_back(push_input(p, true));
  std::vector<std::vector<int>> hole_id;
  for (const auto& hl : holes) {
    std::vector<int> ids;
    for (const auto& p : hl) ids.push_back(push_input(p, true));
    hole_id.push_back(std::move(ids));
  }
  std::vector<std::array<int, 2>> seg_id;
  for (const auto& s : in.segments)
    seg_id.push_back({{push_input(s.first, false), push_input(s.second, false)}});
  for (const auto& p : in.points) push_input(p, true);

  // ---- 2. Delaunay triangulation of the vertex list -----------------------
  Mesher M;
  M.eps = eps;
  M.P.assign(3, Vec2::Zero());
  M.vtri.assign(3, -1);
  const Vec2 ctr = 0.5 * (lo + hi);
  const double Rs = 5.0 * diag;
  for (int k = 0; k < 3; ++k) {
    const double th = kPi * (0.5 + 2.0 * k / 3.0);
    M.P[static_cast<std::size_t>(k)] =
        Vec2(ctr.x() + Rs * std::cos(th), ctr.y() + Rs * std::sin(th));
  }
  M.add_tri(0, 1, 2);

  int hint = 0;
  for (std::size_t k = 0; k < verts.size(); ++k) {
    const int t = M.walk(hint, verts[k], false, nullptr);
    require(t >= 0, "mesh2d: point location failed during triangulation");
    const int c = M.classify(t, verts[k]);
    require(c >= 0 && c < 10, "mesh2d: degenerate or duplicate input point");
    const int nv = M.insert_vertex(verts[k], t, c == 0 ? -1 : c - 1);
    require(nv == static_cast<int>(k) + 3, "mesh2d: internal error, vertex numbering drifted");
    hint = M.vtri[static_cast<std::size_t>(nv)];
  }

  // ---- 3. constraint recovery --------------------------------------------
  auto add_chain = [&](const std::vector<int>& ids, bool boundary) {
    for (std::size_t i = 0; i < ids.size(); ++i) {
      const int a = ids[i] + 3;
      const int b = ids[(i + 1) % ids.size()] + 3;
      M.orig.push_back({M.P[static_cast<std::size_t>(a)], M.P[static_cast<std::size_t>(b)],
                        boundary});
      M.insert_constraint(a, b, static_cast<int>(M.orig.size()) - 1, 0.0, 1.0, 0);
    }
  };
  add_chain(out_id, true);
  for (const auto& ids : hole_id) add_chain(ids, true);
  for (const auto& s : seg_id) {
    const int a = s[0] + 3, b = s[1] + 3;
    if (a == b) continue;
    M.orig.push_back({M.P[static_cast<std::size_t>(a)], M.P[static_cast<std::size_t>(b)], false});
    M.insert_constraint(a, b, static_cast<int>(M.orig.size()) - 1, 0.0, 1.0, 0);
  }

  // ---- 4. domain classification ------------------------------------------
  M.flood_domain();

  // ---- 5. refinement parameters ------------------------------------------
  M.h = in.target_h > 0 ? in.target_h : std::max(hi.x() - lo.x(), hi.y() - lo.y()) / 10.0;
  require(M.h > 0, "mesh2d: could not derive a target edge length");
  M.size_r = M.h / std::sqrt(3.0);
  M.max_ratio = 1.0 / (2.0 * std::sin(in.min_angle_deg * kPi / 180.0));
  const double dom_area = a_out - a_holes;
  require(dom_area > 0, "mesh2d: holes consume the whole outline");
  const double exp_tris = dom_area / (0.25 * std::sqrt(3.0) * M.h * M.h);
  const double exp_nodes = std::max(16.0, std::ceil(0.5 * exp_tris));
  M.vert_budget = 3 + static_cast<int>(verts.size()) +
                  static_cast<int>(std::min(2.0e8, 200.0 * exp_nodes));

  // ---- 6. pre-split every constraint to at most h -------------------------
  if (in.refine) {
    const std::size_t n0 = M.segs.size();
    for (std::size_t si = 0; si < n0; ++si) {
      if (!M.segs[si].alive) continue;
      const double L = (M.P[static_cast<std::size_t>(M.segs[si].b)] -
                        M.P[static_cast<std::size_t>(M.segs[si].a)])
                           .norm();
      const int n = static_cast<int>(std::ceil(L / M.h - 1e-9));
      if (n < 2) continue;
      require(static_cast<int>(M.P.size()) + n <= M.vert_budget,
              "mesh2d: target_h too small for the size budget");
      const double t0 = M.segs[si].t0, t1 = M.segs[si].t1;
      int cur = static_cast<int>(si);
      for (int qi = 1; qi < n; ++qi) {
        require(M.split_subseg_at(cur, t0 + (t1 - t0) * static_cast<double>(qi) / n) >= 0,
                "mesh2d: target_h is below the floating-point resolution of the geometry");
        cur = static_cast<int>(M.segs.size()) - 1;
      }
    }
  }

  // ---- 7. Ruppert refinement ---------------------------------------------
  bool complete = true;
  if (in.refine) {
    std::deque<int> segq, triq;
    M.grid_init(M.h);
    std::vector<int> hits;
    auto after_insert = [&](const Vec2& p, std::size_t t_before, std::size_t s_before) {
      for (std::size_t i = t_before; i < M.T.size(); ++i)
        if (M.tri_bad(static_cast<int>(i))) triq.push_back(static_cast<int>(i));
      // the older sub-segments the point encroaches upon, then every new one
      M.diametral_hits(p, s_before, hits);
      for (int i : hits) segq.push_back(i);
      for (std::size_t i = s_before; i < M.segs.size(); ++i)
        if (M.segs[i].alive) segq.push_back(static_cast<int>(i));
    };
    auto split_seg = [&](int si) {
      const std::size_t tb = M.T.size(), sb = M.segs.size();
      const int nv = M.split_subseg_at(si, 0.5 * (M.segs[static_cast<std::size_t>(si)].t0 +
                                                  M.segs[static_cast<std::size_t>(si)].t1));
      if (nv < 0) return false;  // hit the splitting floor: leave it, stay honest
      after_insert(M.P[static_cast<std::size_t>(nv)], tb, sb);
      return true;
    };

    bool progress = true;
    for (int round = 0; progress && round < 1000; ++round) {
      progress = false;
      segq.clear();
      triq.clear();
      for (std::size_t i = 0; i < M.segs.size(); ++i)
        if (M.seg_encroached(static_cast<int>(i))) segq.push_back(static_cast<int>(i));
      for (std::size_t i = 0; i < M.T.size(); ++i)
        if (M.tri_bad(static_cast<int>(i))) triq.push_back(static_cast<int>(i));
      if (segq.empty() && triq.empty()) break;

      while (!segq.empty() || !triq.empty()) {
        if (static_cast<int>(M.P.size()) > M.vert_budget) {
          complete = false;
          break;
        }
        if (!segq.empty()) {
          const int si = segq.front();
          segq.pop_front();
          if (!M.segs[static_cast<std::size_t>(si)].alive || !M.seg_encroached(si)) continue;
          if (split_seg(si)) progress = true;
          continue;
        }
        const int t = triq.front();
        triq.pop_front();
        if (!M.tri_bad(t)) continue;
        const Vec2 cc = circumcenter(M.pt(t, 0), M.pt(t, 1), M.pt(t, 2));
        if (!cc.allFinite()) continue;

        // Ruppert: a circumcenter that encroaches upon sub-segments is
        // discarded and the sub-segments are split instead.
        std::vector<int> enc;
        M.diametral_hits(cc, M.segs.size(), enc);
        if (!enc.empty()) {
          for (const int si : enc)
            if (M.segs[static_cast<std::size_t>(si)].alive && split_seg(si)) progress = true;
          continue;
        }
        // Visibility: a sub-segment blocking the way to the circumcenter is
        // itself encroached by it, so split that one.
        std::array<int, 2> blocked{{-1, -1}};
        const int loc = M.walk(t, cc, true, &blocked);
        if (loc < 0) {
          if (blocked[0] >= 0) {
            const int bt = blocked[0], bi = blocked[1];
            const int si = M.seg_index(M.T[bt].v[(bi + 1) % 3], M.T[bt].v[(bi + 2) % 3]);
            if (si >= 0 && split_seg(si)) progress = true;
          }
          continue;
        }
        if (!M.T[loc].inside) continue;
        const int cl = M.classify(loc, cc);
        if (cl < 0 || cl >= 10) continue;
        bool too_close = false;
        for (int k = 0; k < 3; ++k)
          if ((M.pt(loc, k) - cc).norm() <= M.eps) too_close = true;
        if (too_close) continue;
        const std::size_t tb = M.T.size(), sb = M.segs.size();
        M.insert_vertex(cc, loc, cl == 0 ? -1 : cl - 1);
        after_insert(cc, tb, sb);
        progress = true;
      }
      if (!complete) break;
    }
    if (complete) {
      for (std::size_t i = 0; i < M.T.size() && complete; ++i)
        if (M.tri_bad(static_cast<int>(i))) complete = false;
      for (std::size_t i = 0; i < M.segs.size() && complete; ++i)
        if (M.seg_encroached(static_cast<int>(i))) complete = false;
    }
  }

  // ---- 8. output ----------------------------------------------------------
  TriMesh m;
  m.nodes.assign(M.P.begin() + 3, M.P.end());
  for (std::size_t i = 0; i < M.T.size(); ++i) {
    if (!M.T[i].alive || !M.T[i].inside) continue;
    require(M.T[i].v[0] >= 3 && M.T[i].v[1] >= 3 && M.T[i].v[2] >= 3,
            "mesh2d: internal error, a domain triangle touches the super-triangle");
    m.tris.push_back({{M.T[i].v[0] - 3, M.T[i].v[1] - 3, M.T[i].v[2] - 3}});
  }
  std::vector<int> order;
  for (std::size_t i = 0; i < M.segs.size(); ++i)
    if (M.segs[i].alive) order.push_back(static_cast<int>(i));
  std::stable_sort(order.begin(), order.end(), [&](int x, int y) {
    const SubSeg& sx = M.segs[static_cast<std::size_t>(x)];
    const SubSeg& sy = M.segs[static_cast<std::size_t>(y)];
    if (sx.orig != sy.orig) return sx.orig < sy.orig;
    return sx.t0 < sy.t0;
  });
  for (const int i : order) {
    const SubSeg& s = M.segs[static_cast<std::size_t>(i)];
    const std::array<int, 2> e{{s.a - 3, s.b - 3}};
    if (M.orig[static_cast<std::size_t>(s.orig)].boundary)
      m.boundary_edges.push_back(e);
    else
      m.constraint_edges.push_back(e);
  }
  m.refinement_complete = complete;
  require(complete || in.allow_partial,
          "mesh2d: refinement could not reach the angle/size policy within the size budget");
  return m;
}

int locate(const TriMesh& m, const Vec2& p) {
  const int nT = static_cast<int>(m.tris.size());
  const int nN = static_cast<int>(m.nodes.size());
  if (nT == 0 || nN == 0) return -1;

  // Vertex -> incident triangles, by counting sort (deterministic, no hashing).
  std::vector<int> start(static_cast<std::size_t>(nN) + 1, 0);
  for (const auto& t : m.tris)
    for (int k = 0; k < 3; ++k) ++start[static_cast<std::size_t>(t[k]) + 1];
  for (int i = 0; i < nN; ++i) start[static_cast<std::size_t>(i) + 1] += start[static_cast<std::size_t>(i)];
  std::vector<int> pos = start;
  std::vector<int> inc(static_cast<std::size_t>(3 * nT));
  for (int i = 0; i < nT; ++i)
    for (int k = 0; k < 3; ++k)
      inc[static_cast<std::size_t>(pos[static_cast<std::size_t>(m.tris[static_cast<std::size_t>(i)][k])]++)] = i;

  auto node = [&](int t, int k) -> const Vec2& {
    return m.nodes[static_cast<std::size_t>(m.tris[static_cast<std::size_t>(t)][k])];
  };
  auto contains = [&](int t) {
    for (int i = 0; i < 3; ++i)
      if (orient2d(node(t, (i + 1) % 3), node(t, (i + 2) % 3), p) < 0) return false;
    return true;
  };
  auto neighbor = [&](int t, int i) -> int {
    const int u = m.tris[static_cast<std::size_t>(t)][(i + 1) % 3];
    const int w = m.tris[static_cast<std::size_t>(t)][(i + 2) % 3];
    for (int q = start[static_cast<std::size_t>(u)]; q < start[static_cast<std::size_t>(u) + 1]; ++q) {
      const int s = inc[static_cast<std::size_t>(q)];
      if (s == t) continue;
      int hit = 0;
      for (int k = 0; k < 3; ++k) {
        const int vv = m.tris[static_cast<std::size_t>(s)][k];
        if (vv == u || vv == w) ++hit;
      }
      if (hit == 2) return s;
    }
    return -1;
  };

  int t = start[0] < start[1] ? inc[static_cast<std::size_t>(start[0])] : 0;
  int prev = -1;
  const int limit = 8 * nT + 64;
  for (int it = 0; it < limit; ++it) {
    bool moved = false;
    for (int k = 0; k < 3; ++k) {
      const int i = (it + k) % 3;
      if (orient2d(node(t, (i + 1) % 3), node(t, (i + 2) % 3), p) >= 0) continue;
      const int nx = neighbor(t, i);
      if (nx < 0 || nx == prev) continue;
      prev = t;
      t = nx;
      moved = true;
      break;
    }
    if (!moved) {
      if (contains(t)) return t;
      break;
    }
  }
  // The walk left the triangulated region (non-convex domain, hole, or the
  // point is simply outside): finish with a deterministic linear scan.
  for (int i = 0; i < nT; ++i)
    if (contains(i)) return i;
  return -1;
}

}  // namespace keta::mesh2

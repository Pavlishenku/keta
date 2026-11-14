// The exact section nodes:
//   Section (Exact)     an arbitrary holed polygon on wired pins;
//   Section (Editor)    the same polygon carried in the params (GUI-edited);
//   Section (Composite) a modular-ratio multi-material section (see below).
//
// Points are 3D; the section plane takes each point's (x, y) as (y, z) and
// ignores z. Coordinates in metres. A, Iy, Iz and Wel come from keta's closed
// form polygon integrals (exact for a polygon); the Tri6 warping FE replaces
// J, Avy and Avz and provides Iw and the shear centre, which have no home in
// the section properties and leave on their own pins.
//
// Mesh density: the target is ~3000 Tri6 elements at fineness 1 (clamped to
// 200..200000, proportional to the fineness). The mesher takes a target edge,
// so the count is converted with the equilateral-triangle area sqrt(3)/4 h^2
// over the net polygon area.
//
// Principal axes: keta's frame element requires principal section axes, so a
// polygon whose |Iyz| exceeds 5 % of max(Iy, Iz) is rejected (a typed error).
#include <algorithm>
#include <array>
#include <cmath>
#include <set>

#include "keta/param/components.h"
#include "keta/section_exact.h"

namespace keta::param {

namespace {

std::vector<Vec2> ring_of(const Branch<Value>& items) {
  std::vector<Vec2> ring;
  for (const auto& o : items)
    if (o)
      if (auto p = o->as_point()) ring.emplace_back(p->x(), p->y());
  return ring;
}

std::vector<Vec3> to_shape(const std::vector<Vec2>& ring) {
  std::vector<Vec3> out;
  out.reserve(ring.size());
  for (const Vec2& p : ring) out.emplace_back(p.x(), p.y(), 0.0);
  return out;
}

bool finite_ring(const std::vector<Vec2>& r) {
  return std::all_of(r.begin(), r.end(), [](const Vec2& p) { return std::isfinite(p.x()) && std::isfinite(p.y()); });
}

// Element budget of the mesh: 3000 * f elements, clamped to 200..200000.
double target_elements(double f) {
  if (std::isfinite(f) && f > 0) return std::clamp(std::round(3000.0 * f), 200.0, 200000.0);
  return 3000.0;
}

// Edge of the equilateral triangle that spends `n` elements on `area`.
double target_edge(double area, double n) { return std::sqrt(4.0 * area / (std::sqrt(3.0) * n)); }

// The exact-section core shared by Section (Exact) and Section (Editor): the
// same outer, holes, nu and fineness give identical outputs S, J, Iw, SC.
// `name` labels the produced section.
Outputs<Value> solve_exact_section(const SolveCtx& ctx, const std::vector<Vec2>& outer,
                                   std::vector<std::vector<Vec2>> holes, double nu, double fineness,
                                   const std::string& name) {
  if (outer.size() < 3) throw ComponentError("outer ring needs at least 3 points");
  if (!finite_ring(outer) || !std::all_of(holes.begin(), holes.end(), finite_ring))
    throw ComponentError("section geometry has a non-finite coordinate");
  const double n_target = target_elements(fineness);

  if (ctx.cancelled()) throw ComponentError::cancelled();
  Section s;
  s.name = name;
  ShapePolygon shape;
  shape.outer = to_shape(outer);
  for (const auto& h : holes) shape.holes.push_back(to_shape(h));
  s.shape = shape;
  std::shared_ptr<const ExactSection> ex;
  try {
    s.finalize();  // closed-form A, Iy, Iz, Wel (rejects non-principal axes)
    ExactInput xi;
    xi.outline = outer;
    xi.holes = std::move(holes);
    xi.materials = {ExactMaterial{"", 1.0, nu}};
    xi.target_h = target_edge(s.props.A, n_target);
    ex = std::make_shared<const ExactSection>(analyze_exact(xi));
    attach_exact(s, ex);
  } catch (const std::exception& e) {
    throw ComponentError(std::string("exact section: ") + e.what());
  }
  return {single(Value::of(SectionSpec{s})), single(Value::number(ex->J)), single(Value::number(ex->Iw)),
          single(Value::of(Point{Vec3(ex->y_sc, ex->z_sc, 0.0)}))};
}

}  // namespace

ComponentInfo SectionExactNode::info() const {
  return {"keta.section_exact", "Section (Exact)", "Model",
          "Exact J / Iw / shear areas of an arbitrary polygon (warping FE)"};
}
std::vector<PinSpec> SectionExactNode::inputs() const {
  return {PinSpec::list("P", "point", "Outer ring points (section y = x, z = y)"),
          PinSpec::tree("H", "point", "Hole rings, one branch each (optional)").opt(),
          PinSpec::item("N", "number", "Mesh fineness factor (default 1.0)").opt()};
}
std::vector<PinSpec> SectionExactNode::outputs() const {
  return {PinSpec::item("S", "section", "The exact section"),
          PinSpec::item("J", "number", "St-Venant torsion constant J (m^4)"),
          PinSpec::item("Iw", "number", "Warping constant Cw (m^6)"),
          PinSpec::item("SC", "point", "Shear-centre offset from centroid (y, z)")};
}
Outputs<Value> SectionExactNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  std::vector<Vec2> outer;
  for (const Value& v : input_or_empty(in, 0).non_nulls())
    if (auto p = v.as_point()) outer.emplace_back(p->x(), p->y());
  if (outer.size() < 3) throw ComponentError("P: outer ring needs at least 3 points");
  std::vector<std::vector<Vec2>> holes;
  for (const auto& kv : input_or_empty(in, 1).branches())
    if (auto ring = ring_of(kv.second); !ring.empty()) holes.push_back(std::move(ring));
  return solve_exact_section(ctx, outer, std::move(holes), nu, first_number(input_or_empty(in, 2)).value_or(fineness),
                             "Exact");
}
json::Value SectionExactNode::params() const {
  json::Value p = json::Value::object();
  p.set("fineness", fineness);
  p.set("nu", nu);
  return p;
}
void SectionExactNode::set_params(const json::Value& p) {
  if (const auto* x = p.get("fineness"); x && x->is_number()) fineness = x->as_number();
  if (const auto* x = p.get("nu"); x && x->is_number()) nu = x->as_number();
}

// ===========================================================================
// Section (Editor)
// ===========================================================================

namespace {

// A ring [[y, z], ...] from JSON, skipping malformed entries.
std::vector<Vec2> ring_from_json(const json::Value& v) {
  std::vector<Vec2> ring;
  if (!v.is_array()) return ring;
  for (const json::Value& p : v.as_array()) {
    if (!p.is_array() || p.as_array().size() < 2) continue;
    const json::Value& y = p.as_array()[0];
    const json::Value& z = p.as_array()[1];
    if (y.is_number() && z.is_number()) ring.emplace_back(y.as_number(), z.as_number());
  }
  return ring;
}

json::Value ring_to_json(const std::vector<Vec2>& ring) {
  json::Value a = json::Value::array();
  for (const Vec2& p : ring) {
    json::Value q = json::Value::array();
    q.push_back(p.x());
    q.push_back(p.y());
    a.push_back(std::move(q));
  }
  return a;
}

}  // namespace

ComponentInfo SectionEditorNode::info() const {
  return {"keta.section_editor", "Section (Editor)", "Model",
          "Interactively-edited exact polygon section (J / Iw / shear, warping FE)"};
}
std::vector<PinSpec> SectionEditorNode::inputs() const { return {}; }
std::vector<PinSpec> SectionEditorNode::outputs() const {
  return {PinSpec::item("S", "section", "The exact section"),
          PinSpec::item("J", "number", "St-Venant torsion constant J (m^4)"),
          PinSpec::item("Iw", "number", "Warping constant Cw (m^6)"),
          PinSpec::item("SC", "point", "Shear-centre offset from centroid (y, z)")};
}
Outputs<Value> SectionEditorNode::solve(const SolveCtx& ctx, const std::vector<Tree>&) const {
  return solve_exact_section(ctx, outer, holes, nu, fineness, "Editor");
}
json::Value SectionEditorNode::params() const {
  json::Value p = json::Value::object();
  p.set("outer", ring_to_json(outer));
  json::Value h = json::Value::array();
  for (const auto& r : holes) h.push_back(ring_to_json(r));
  p.set("holes", std::move(h));
  p.set("nu", nu);
  p.set("fineness", fineness);
  p.set("grid_step", grid_step);
  p.set("snap", snap);
  return p;
}
void SectionEditorNode::set_params(const json::Value& p) {
  // A truncated outer (< 3 points) keeps the current polygon rather than
  // blanking the geometry.
  if (const auto* o = p.get("outer"))
    if (auto r = ring_from_json(*o); r.size() >= 3) outer = std::move(r);
  if (const auto* h = p.get("holes"); h && h->is_array()) {
    holes.clear();
    for (const json::Value& x : h->as_array())
      if (auto r = ring_from_json(x); !r.empty()) holes.push_back(std::move(r));
  }
  if (const auto* x = p.get("nu"); x && x->is_number()) nu = x->as_number();
  if (const auto* x = p.get("fineness"); x && x->is_number()) fineness = x->as_number();
  if (const auto* x = p.get("grid_step"); x && x->is_number()) grid_step = x->as_number();
  if (const auto* x = p.get("snap"); x && x->is_bool()) snap = x->as_bool();
}

// ===========================================================================
// Section (Composite)
// ===========================================================================
//
// Overlapping region rings resolve last-wins: the exact engine takes one
// domain (outline + holes) plus material regions assigned last-wins by
// element centroid, with every region boundary forced onto element edges. The
// domain is the UNION of the rings, computed here by a small arrangement:
// every ring edge is split at the vertices of the other rings lying on it and
// at proper crossings; a sub-edge is on the union boundary iff exactly one of
// its two sides is covered by a ring; the oriented boundary sub-edges chain
// into loops (CCW outer, CW holes). The split rings are handed over as the
// material regions, so every interface vertex is shared verbatim.
//
// J_eff = (GJ)_eff / G_ref: the torsion problem is G-weighted.
// With a common Poisson ratio the G and E weights are proportional and the
// E-weighted J is already that; otherwise a second solve on the same mesh
// with E_e replaced by G_e recovers it (J does not depend on the reference
// point, so the E-centroid shift of the mesh is harmless).
//
// Conventions: the shear areas use the E-weighted mean Poisson ratio as
// normaliser (equal to EA / 2GA - 1 for a common ratio); like Section (Exact),
// a section whose |Iyz| exceeds 5 % of max(Iy, Iz) is rejected.

namespace {

struct UnionResult {
  std::vector<Vec2> outer;
  std::vector<std::vector<Vec2>> holes;
  std::vector<std::vector<Vec2>> split_rings;  // input rings with every arrangement vertex inserted
  double area = 0;
};

double cross2(const Vec2& a, const Vec2& b) { return a.x() * b.y() - a.y() * b.x(); }

UnionResult ring_union(const std::vector<std::vector<Vec2>>& rings) {
  Vec2 lo = rings[0][0], hi = rings[0][0];
  for (const auto& r : rings)
    for (const Vec2& p : r) {
      lo = lo.cwiseMin(p);
      hi = hi.cwiseMax(p);
    }
  const double diag = (hi - lo).norm();
  const double eps = 1e-9 * diag;

  // Split parameters per ring edge.
  struct Edge {
    std::size_t r, i;
    Vec2 a, b;
    std::vector<double> ts;
  };
  std::vector<Edge> edges;
  for (std::size_t r = 0; r < rings.size(); ++r)
    for (std::size_t i = 0; i < rings[r].size(); ++i)
      edges.push_back({r, i, rings[r][i], rings[r][(i + 1) % rings[r].size()], {0.0, 1.0}});
  for (std::size_t e = 0; e < edges.size(); ++e) {
    Edge& E = edges[e];
    const Vec2 d = E.b - E.a;
    const double L2 = d.squaredNorm();
    if (!(L2 > 0)) continue;
    for (std::size_t f = 0; f < edges.size(); ++f) {
      if (f == e) continue;
      const Edge& F = edges[f];
      // Vertices of F lying on E.
      for (const Vec2& q : {F.a, F.b}) {
        const double t = (q - E.a).dot(d) / L2;
        if (t > 0 && t < 1 && (E.a + t * d - q).norm() <= eps) E.ts.push_back(t);
      }
      // Proper crossing.
      const Vec2 g = F.b - F.a;
      const double den = cross2(d, g);
      if (std::fabs(den) <= 1e-14 * std::sqrt(L2) * g.norm()) continue;
      const double t = cross2(F.a - E.a, g) / den;
      const double u = cross2(F.a - E.a, d) / den;
      if (t > 0 && t < 1 && u > 0 && u < 1) E.ts.push_back(t);
    }
    std::sort(E.ts.begin(), E.ts.end());
  }

  // Vertex pool with coincidence merging (the first spelling wins).
  std::vector<Vec2> verts;
  auto vid = [&](const Vec2& p) {
    for (std::size_t k = 0; k < verts.size(); ++k)
      if ((verts[k] - p).norm() <= eps) return static_cast<int>(k);
    verts.push_back(p);
    return static_cast<int>(verts.size() - 1);
  };
  UnionResult out;
  out.split_rings.assign(rings.size(), {});
  std::vector<std::array<int, 2>> subs;
  std::set<std::pair<int, int>> seen;
  std::vector<std::vector<int>> ring_ids(rings.size());
  for (const Edge& E : edges) {
    int prev = -1;
    for (double t : E.ts) {
      const int v = vid(E.a + t * (E.b - E.a));
      if (v == prev) continue;
      if (t < 1) ring_ids[E.r].push_back(v);
      if (prev >= 0 && seen.insert({std::min(prev, v), std::max(prev, v)}).second) subs.push_back({prev, v});
      prev = v;
    }
  }
  for (std::size_t r = 0; r < rings.size(); ++r) {
    auto& ids = ring_ids[r];
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
    if (ids.size() > 1 && ids.front() == ids.back()) ids.pop_back();
    for (int v : ids) out.split_rings[r].push_back(verts[static_cast<std::size_t>(v)]);
  }

  // Boundary sub-edges, oriented with the covered side on the left.
  auto covered = [&](const Vec2& p) {
    return std::any_of(rings.begin(), rings.end(), [&](const auto& r) { return point_in_polygon(r, p); });
  };
  std::vector<std::array<int, 2>> bnd;
  for (const auto& s : subs) {
    const Vec2 a = verts[static_cast<std::size_t>(s[0])], b = verts[static_cast<std::size_t>(s[1])];
    const Vec2 d = b - a;
    const double L = d.norm();
    const Vec2 n(-d.y() / L, d.x() / L);
    const double delta = std::min(1e-6 * diag, 1e-3 * L);
    const Vec2 m = 0.5 * (a + b);
    const bool left = covered(m + delta * n), right = covered(m - delta * n);
    if (left == right) continue;
    bnd.push_back(left ? s : std::array<int, 2>{s[1], s[0]});
  }

  // Chain into loops; at a pinch vertex take the outgoing edge with the
  // smallest clockwise turn from the reversed incoming edge.
  std::vector<char> used(bnd.size(), 0);
  std::vector<std::vector<Vec2>> loops;
  for (std::size_t s0 = 0; s0 < bnd.size(); ++s0) {
    if (used[s0]) continue;
    std::vector<Vec2> loop;
    std::size_t cur = s0;
    while (!used[cur]) {
      used[cur] = 1;
      const int a = bnd[cur][0], b = bnd[cur][1];
      loop.push_back(verts[static_cast<std::size_t>(a)]);
      const Vec2 back = verts[static_cast<std::size_t>(a)] - verts[static_cast<std::size_t>(b)];
      const double ang_back = std::atan2(back.y(), back.x());
      std::size_t best = bnd.size();
      double best_turn = 1e300;
      for (std::size_t k = 0; k < bnd.size(); ++k) {
        if (used[k] || bnd[k][0] != b) continue;
        const Vec2 d = verts[static_cast<std::size_t>(bnd[k][1])] - verts[static_cast<std::size_t>(b)];
        double turn = ang_back - std::atan2(d.y(), d.x());
        while (turn <= 0) turn += 2 * M_PI;
        if (turn < best_turn) {
          best_turn = turn;
          best = k;
        }
      }
      if (best == bnd.size()) break;
      cur = best;
    }
    if (loop.size() >= 3) loops.push_back(std::move(loop));
  }

  std::vector<std::vector<Vec2>> outers;
  for (auto& l : loops) {
    const double a = polygon_area(l);
    out.area += a;
    (a > 0 ? outers : out.holes).push_back(std::move(l));
  }
  if (outers.empty()) throw ComponentError("composite section: the region rings enclose no area");
  if (outers.size() > 1)
    throw ComponentError("composite section: the region rings form " + std::to_string(outers.size()) +
                         " disjoint parts — Saint-Venant torsion is not defined component by component; "
                         "analyse each part separately");
  out.outer = std::move(outers[0]);
  return out;
}

}  // namespace

ComponentInfo SectionCompositeNode::info() const {
  return {"keta.section_composite", "Section (Composite)", "Model",
          "Modular-ratio transformed multi-material section (warping FE)"};
}
std::vector<PinSpec> SectionCompositeNode::inputs() const {
  return {PinSpec::tree("R", "point", "Region rings, one branch each (section y=x, z=y)"),
          PinSpec::list("M", "material", "One material per region (first is reference)"),
          PinSpec::item("N", "number", "Mesh fineness factor (default 1.0)").opt()};
}
std::vector<PinSpec> SectionCompositeNode::outputs() const {
  return {PinSpec::item("S", "section", "The transformed section"),
          PinSpec::item("M", "material", "Reference material (first input)"),
          PinSpec::item("J", "number", "Effective torsion constant J_eff (m^4)"),
          PinSpec::item("Iw", "number", "Warping constant Cw (m^6)")};
}
Outputs<Value> SectionCompositeNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  std::vector<std::vector<Vec2>> rings;
  for (const auto& kv : input_or_empty(in, 0).branches())
    if (auto ring = ring_of(kv.second); !ring.empty()) rings.push_back(std::move(ring));
  if (rings.empty()) throw ComponentError("R: at least one region ring is required");
  std::vector<Material> mats;
  for (const Value& v : input_or_empty(in, 1).non_nulls())
    if (const auto* m = v.get<MaterialSpec>()) mats.push_back(m->m);
  if (mats.size() != rings.size())
    throw ComponentError("M: " + std::to_string(mats.size()) + " materials for " + std::to_string(rings.size()) +
                         " region rings (counts must match)");
  // Input checks of the composite backend.
  auto fail_c = [](const std::string& msg) { throw ComponentError("composite section: " + msg); };
  const double e_ref = mats[0].E;
  if (!(std::isfinite(e_ref) && e_ref > 0))
    fail_c("reference modulus E_ref must be > 0 and finite (got " + Value::number(e_ref).describe() + ")");
  for (std::size_t i = 0; i < rings.size(); ++i) {
    if (rings[i].size() < 3) fail_c("degenerate polygon: region ring needs at least 3 points");
    if (!finite_ring(rings[i])) fail_c("degenerate polygon: non-finite coordinate");
    if (!(std::isfinite(mats[i].E) && mats[i].E > 0))
      fail_c("region modulus E must be > 0 and finite (got " + Value::number(mats[i].E).describe() + ")");
  }
  const double n_target = target_elements(first_number(input_or_empty(in, 2)).value_or(fineness));
  if (ctx.cancelled()) throw ComponentError::cancelled();

  UnionResult u = ring_union(rings);
  std::shared_ptr<ExactSection> ex;
  try {
    ExactInput xi;
    xi.outline = u.outer;
    xi.holes = u.holes;
    for (std::size_t i = 0; i < mats.size(); ++i) {
      xi.materials.push_back(ExactMaterial{mats[i].name, mats[i].E, mats[i].nu});
      xi.regions.push_back(ExactRegion{u.split_rings[i], static_cast<int>(i)});
    }
    xi.target_h = target_edge(u.area, n_target);
    ex = std::make_shared<ExactSection>(analyze_exact(xi));
    const bool uniform_nu = std::all_of(mats.begin(), mats.end(),
                                        [&](const Material& m) { return std::fabs(m.nu - mats[0].nu) < 1e-12; });
    if (!uniform_nu) {
      std::vector<ExactMaterial> g_mats;
      for (const Material& m : mats) g_mats.push_back(ExactMaterial{m.name, m.G(), m.nu});
      ex->J = analyze_exact(ex->mesh, g_mats).J;
    }
  } catch (const ComponentError&) {
    throw;
  } catch (const std::exception& e) {
    throw ComponentError(std::string("composite section: ") + e.what());
  }
  if (std::fabs(ex->Iyz) > 0.05 * std::max(ex->Iy, ex->Iz))
    throw ComponentError(
        "composite section: section axes are not principal (|Iyz| > 5% of max(Iy,Iz)); rotate the input");

  // Display outline: the largest-area ring as the outer, the rest as holes
  // (a best-effort picture of a nested composite).
  std::vector<std::size_t> order(rings.size());
  for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
    return std::fabs(polygon_area(rings[a])) > std::fabs(polygon_area(rings[b]));
  });
  Section s;
  s.name = "Composite";
  ShapePolygon shape;
  shape.outer = to_shape(rings[order[0]]);
  for (std::size_t k = 1; k < order.size(); ++k) shape.holes.push_back(to_shape(rings[order[k]]));
  s.shape = shape;
  s.props = ex->props();  // transformed to the reference material
  for (std::size_t i = 0; i < rings.size(); ++i)
    for (std::size_t k = 0; k < rings[i].size(); ++k)
      s.stress_points.push_back({"r" + std::to_string(i) + "v" + std::to_string(k), rings[i][k].x() - ex->yc,
                                 rings[i][k].y() - ex->zc, 0, 0});
  s.exact = ex;
  return {single(Value::of(SectionSpec{s})), single(Value::of(MaterialSpec{mats[0]})), single(Value::number(ex->J)),
          single(Value::number(ex->Iw))};
}
json::Value SectionCompositeNode::params() const {
  json::Value p = json::Value::object();
  p.set("fineness", fineness);
  return p;
}
void SectionCompositeNode::set_params(const json::Value& p) {
  if (const auto* x = p.get("fineness"); x && x->is_number()) fineness = x->as_number();
}

}  // namespace keta::param

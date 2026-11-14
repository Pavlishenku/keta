#include "keta/slabmesh.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "keta/shell.h"

namespace keta {

namespace {

using mesh2::Vec2;

// ---------------------------------------------------------------------------
// 2-D predicates (tolerance semantics live HERE; the mesher has its own).
// ---------------------------------------------------------------------------

double dist_seg2(const Vec2& p, const Vec2& a, const Vec2& b) {
  Vec2 d = b - a;
  double l2 = d.squaredNorm();
  double t = l2 > 0 ? std::clamp((p - a).dot(d) / l2, 0.0, 1.0) : 0.0;
  return (p - (a + t * d)).norm();
}

// Crossing-parity interiority (open set of the ring).
bool in_ring(const std::vector<Vec2>& ring, const Vec2& q) {
  bool in = false;
  const std::size_t n = ring.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Vec2& a = ring[i];
    const Vec2& b = ring[(i + 1) % n];
    if ((a.y() > q.y()) != (b.y() > q.y())) {
      double x = a.x() + (q.y() - a.y()) / (b.y() - a.y()) * (b.x() - a.x());
      if (q.x() < x) in = !in;
    }
  }
  return in;
}

double dist_ring(const std::vector<Vec2>& ring, const Vec2& q) {
  double best = 1e300;
  const std::size_t n = ring.size();
  for (std::size_t i = 0; i < n; ++i)
    best = std::min(best, dist_seg2(q, ring[i], ring[(i + 1) % n]));
  return best;
}

// The meshed domain: outline minus holes, with an `eps` boundary band.
struct Domain2 {
  const std::vector<Vec2>& outline;
  const std::vector<std::vector<Vec2>>& holes;
  double eps;

  bool on_boundary(const Vec2& q) const {
    if (dist_ring(outline, q) <= eps) return true;
    for (const auto& h : holes)
      if (dist_ring(h, q) <= eps) return true;
    return false;
  }
  bool inside_or_boundary(const Vec2& q) const {
    if (on_boundary(q)) return true;
    if (!in_ring(outline, q)) return false;
    for (const auto& h : holes)
      if (in_ring(h, q)) return false;
    return true;
  }
  bool interior_strict(const Vec2& q) const {
    if (on_boundary(q)) return false;
    if (!in_ring(outline, q)) return false;
    for (const auto& h : holes)
      if (in_ring(h, q)) return false;
    return true;
  }
};

// True when segments [a,b] and [c,d] touch anywhere but at a common endpoint
// pair (endpoint of BOTH, within tol): a proper crossing, a T-junction onto an
// interior, or a collinear overlap longer than tol.
bool segments_conflict_tol(const Vec2& a, const Vec2& b, const Vec2& c, const Vec2& d,
                           double tol) {
  Vec2 d1 = b - a, d2 = d - c;
  double l1 = d1.norm(), l2 = d2.norm();
  double denom = d1.x() * d2.y() - d1.y() * d2.x();
  Vec2 ca = c - a;
  if (std::fabs(denom) <= 1e-12 * l1 * l2) {
    // Parallel: a conflict only when collinear AND overlapping by more than tol.
    double off_c = std::fabs(ca.x() * d1.y() - ca.y() * d1.x()) / l1;
    double off_d = std::fabs((d - a).x() * d1.y() - (d - a).y() * d1.x()) / l1;
    if (off_c > tol || off_d > tol) return false;
    double t3 = (c - a).dot(d1) / (l1 * l1);
    double t4 = (d - a).dot(d1) / (l1 * l1);
    double lo = std::max(0.0, std::min(t3, t4));
    double hi = std::min(1.0, std::max(t3, t4));
    return (hi - lo) * l1 > tol;
  }
  double t = (ca.x() * d2.y() - ca.y() * d2.x()) / denom;
  double u = (ca.x() * d1.y() - ca.y() * d1.x()) / denom;
  double e1 = tol / l1, e2 = tol / l2;
  if (t < -e1 || t > 1 + e1 || u < -e2 || u > 1 + e2) return false;  // no contact
  bool end1 = t <= e1 || t >= 1 - e1;
  bool end2 = u <= e2 || u >= 1 - e2;
  return !(end1 && end2);
}

// Mesh vertex indices lying on the axis [A, B], ordered by abscissa. The
// mesher's contract makes every split point collinear with the input segment
// (mesh2d.h), so a tight relative band suffices; candidates are read off the
// constraint edges only.
std::vector<int> axis_chain(const mesh2::TriMesh& m, const Vec2& A, const Vec2& B,
                            const std::string& at, Id bar) {
  Vec2 d = B - A;
  double L2 = d.squaredNorm(), L = std::sqrt(L2);
  std::vector<std::pair<double, int>> on;
  auto consider = [&](int idx) {
    const Vec2& p = m.nodes[static_cast<std::size_t>(idx)];
    double cross = (p - A).x() * d.y() - (p - A).y() * d.x();
    if (std::fabs(cross) / L > 1e-9 * L) return;
    double t = (p - A).dot(d) / L2;
    if (t < -1e-9 || t > 1 + 1e-9) return;
    on.emplace_back(t, idx);
  };
  for (const auto& e : m.constraint_edges) {
    consider(e[0]);
    consider(e[1]);
  }
  std::sort(on.begin(), on.end());
  on.erase(std::unique(on.begin(), on.end()), on.end());
  require(on.size() >= 2 && std::fabs(on.front().first) <= 1e-9 &&
              std::fabs(on.back().first - 1.0) <= 1e-9,
          at + ": internal error — the mesh does not cover the axis of bar " +
              std::to_string(bar));
  std::vector<int> chain;
  chain.reserve(on.size());
  for (const auto& [t, idx] : on) chain.push_back(idx);
  return chain;
}

// Re-slice one bar load onto the sibling segments of a split bar. UDL / linear
// loads become per-segment tranches (the tranche boundaries partition the
// original range, so the resultant is preserved to round-off); a concentrated
// load goes to the containing segment; thermal / prestrain fields are uniform
// per length and are copied per segment.
void slice_bar_load(const Load& l, const SplitBarRecord& sp, const std::string& case_name,
                    std::vector<Load>& out) {
  const double L = sp.s.back();
  const std::string ctx =
      "case '" + case_name + "': load on auto-split bar " + std::to_string(sp.original);
  switch (l.kind) {
    case Load::Kind::BarThermal:
    case Load::Kind::Prestrain: {
      for (Id seg : sp.segments) {
        Load c = l;
        c.bar = seg;
        out.push_back(c);
      }
      return;
    }
    case Load::Kind::BarPoint: {
      double a = l.relative ? l.a * L : l.a;
      require(finite(a) && a >= -1e-9 * L && a <= L * (1 + 1e-9),
              ctx + ": point-load abscissa outside the bar");
      a = std::clamp(a, 0.0, L);
      for (std::size_t k = 0; k < sp.segments.size(); ++k) {
        double s1 = sp.s[k + 1];
        if (a <= s1 + 1e-12 * L || k + 1 == sp.segments.size()) {
          Load c = l;
          c.bar = sp.segments[k];
          c.relative = false;
          c.a = std::clamp(a - sp.s[k], 0.0, s1 - sp.s[k]);
          out.push_back(c);
          return;
        }
      }
      return;  // unreachable
    }
    case Load::Kind::BarUdl:
    case Load::Kind::BarLinear:
    case Load::Kind::BarParabolic: {
      double a0 = l.relative ? l.a * L : l.a;
      double b0 = l.b < 0 ? L : (l.relative ? l.b * L : l.b);
      require(finite(a0) && finite(b0) && a0 >= -1e-9 * L && b0 <= L * (1 + 1e-9) && b0 >= a0,
              ctx + ": abscissa range outside the bar");
      a0 = std::clamp(a0, 0.0, L);
      b0 = std::clamp(b0, 0.0, L);
      double span = b0 - a0;
      for (std::size_t k = 0; k < sp.segments.size(); ++k) {
        double c = std::max(a0, sp.s[k]);
        double d = std::min(b0, sp.s[k + 1]);
        if (d <= c) continue;
        Load nl = l;
        nl.bar = sp.segments[k];
        nl.relative = false;
        nl.a = c - sp.s[k];
        nl.b = d - sp.s[k];
        if (l.kind == Load::Kind::BarLinear && span > 0) {
          double fa = (c - a0) / span, fb = (d - a0) / span;
          nl.F = l.F + (l.F2 - l.F) * fa;
          nl.F2 = l.F + (l.F2 - l.F) * fb;
        }
        if (l.kind == Load::Kind::BarParabolic && span > 0) {
          // a parabola restricted to [c, d] is the parabola through its
          // values at c, (c + d) / 2 and d
          const Vec3 qc = l.Fmid - 0.5 * (l.F + l.F2);
          auto q = [&](double x) {
            const double t = (x - a0) / span;
            return Vec3(l.F + (l.F2 - l.F) * t + 4.0 * t * (1.0 - t) * qc);
          };
          nl.F = q(c);
          nl.Fmid = q(0.5 * (c + d));
          nl.F2 = q(d);
        }
        out.push_back(nl);
      }
      return;
    }
    case Load::Kind::PrestrainLinear:
    case Load::Kind::PrestrainParabolic: {
      // an imposed field over the whole beam: each segment [c, d] takes the
      // field restricted to it — a linear field through its values at c and
      // d, a parabola through its values at c, (c + d) / 2 and d (exact)
      const bool para = l.kind == Load::Kind::PrestrainParabolic;
      require(L > 0, ctx + ": zero-length bar");
      auto field = [&](double fi, double fm, double fj, double x) {
        const double t = x / L, c = para ? fm - 0.5 * (fi + fj) : 0.0;
        return fi + (fj - fi) * t + 4.0 * c * t * (1.0 - t);
      };
      for (std::size_t k = 0; k < sp.segments.size(); ++k) {
        const double c = sp.s[k], d = sp.s[k + 1], m = 0.5 * (c + d);
        Load nl = l;
        nl.bar = sp.segments[k];
        nl.strain = field(l.strain, l.strain_m, l.strain2, c);
        nl.kappa_y = field(l.kappa_y, l.kappa_y_m, l.kappa_y2, c);
        nl.kappa_z = field(l.kappa_z, l.kappa_z_m, l.kappa_z2, c);
        nl.kappa_x = field(l.kappa_x, l.kappa_x_m, l.kappa_x2, c);
        nl.strain2 = field(l.strain, l.strain_m, l.strain2, d);
        nl.kappa_y2 = field(l.kappa_y, l.kappa_y_m, l.kappa_y2, d);
        nl.kappa_z2 = field(l.kappa_z, l.kappa_z_m, l.kappa_z2, d);
        nl.kappa_x2 = field(l.kappa_x, l.kappa_x_m, l.kappa_x2, d);
        if (para) {
          nl.strain_m = field(l.strain, l.strain_m, l.strain2, m);
          nl.kappa_y_m = field(l.kappa_y, l.kappa_y_m, l.kappa_y2, m);
          nl.kappa_z_m = field(l.kappa_z, l.kappa_z_m, l.kappa_z2, m);
          nl.kappa_x_m = field(l.kappa_x, l.kappa_x_m, l.kappa_x2, m);
        }
        out.push_back(nl);
      }
      return;
    }
    case Load::Kind::PrestrainPoint:
      // Same policy as a tendon on a split beam: a concentrated dislocation
      // authored on the WHOLE beam is refused by type, never silently
      // re-sliced.
      fail(ctx + ": DEFOR CON on an auto-split beam is not supported — author the "
                 "dislocation on the mesh-conforming segment instead");
    default:
      fail(ctx + ": unsupported load kind for re-slicing");  // defensive
  }
}

}  // namespace

SlabMeshResult mesh_slab(Model& model, Id slab_id, const SlabMeshOptions& opt) {
  require(finite(opt.weld_tol) && opt.weld_tol > 0,
          "mesh_slab: weld_tol must be finite and positive");
  const Slab& s = model.slab(slab_id);  // typed lookup
  const std::string at = "slab " + std::to_string(slab_id);
  require(finite(s.t) && s.t > 0, at + ": thickness must be finite and positive");
  require(finite(s.mesh_size) && s.mesh_size > 0, at + ": mesh size must be finite and positive");
  require(opt.weld_tol < 0.5 * s.mesh_size,
          at + ": weld_tol must be well below the mesh size");
  (void)model.material(s.material);  // typed dangling-reference check
  require(s.outline.size() >= 3, at + ": outline needs at least 3 points");
  for (const Vec3& p : s.outline) require(finite(p), at + ": non-finite outline point");
  for (const auto& h : s.holes) {
    require(h.size() >= 3, at + ": a hole ring needs at least 3 points");
    for (const Vec3& p : h) require(finite(p), at + ": non-finite hole point");
  }
  for (const ShellTri& sh : model.shells)
    require(sh.slab != slab_id,
            at + ": already meshed (shell " + std::to_string(sh.id) + " carries its id)");

  // ---- plane fit (Newell normal through outline[0]) ----
  const Vec3 origin = s.outline.front();
  Vec3 nrm = Vec3::Zero();
  Vec3 lo = origin, hi = origin;
  for (std::size_t i = 0; i < s.outline.size(); ++i) {
    const Vec3& a = s.outline[i];
    const Vec3& b = s.outline[(i + 1) % s.outline.size()];
    nrm += a.cross(b);
    lo = lo.cwiseMin(a);
    hi = hi.cwiseMax(a);
  }
  double diag = (hi - lo).norm();
  require(nrm.norm() > 1e-12 * std::max(diag * diag, 1e-12),
          at + ": degenerate (zero-area) outline");
  const Vec3 unit_n = nrm.normalized();
  auto plane_dev = [&](const Vec3& p) { return std::fabs((p - origin).dot(unit_n)); };
  // Planarity gate = the WELDING tolerance (stricter than Model::validate's
  // 1e-6*diag record gate): the raise step projects every ring vertex exactly
  // into the plane, and a corner farther off-plane than weld_tol could no
  // longer weld back onto its own pre-existing node.
  {
    double dev = 0;
    for (const Vec3& p : s.outline) dev = std::max(dev, plane_dev(p));
    for (const auto& h : s.holes)
      for (const Vec3& p : h) dev = std::max(dev, plane_dev(p));
    require(dev <= opt.weld_tol,
            at + ": outline/hole ring is not planar within the welding tolerance (deviation " +
                std::to_string(dev) + " m > weld_tol " + std::to_string(opt.weld_tol) + " m)");
  }
  // The ONE shared frame definition (shell.h). The mesh frame deliberately
  // ignores Slab::local_x: that override steers RESULT axes, not geometry.
  const Mat3 fr = plane_frame(unit_n);
  const Vec3 xh = fr.row(0), yh = fr.row(1);
  auto to2 = [&](const Vec3& p) { return Vec2((p - origin).dot(xh), (p - origin).dot(yh)); };
  auto to3 = [&](const Vec2& q) { return Vec3(origin + q.x() * xh + q.y() * yh); };

  // Project the rings straight into the mesher input; `dom` references them.
  mesh2::MeshInput min2;
  min2.outline.reserve(s.outline.size());
  for (const Vec3& p : s.outline) min2.outline.push_back(to2(p));
  for (const auto& h : s.holes) {
    std::vector<Vec2> h2;
    h2.reserve(h.size());
    for (const Vec3& p : h) h2.push_back(to2(p));
    min2.holes.push_back(std::move(h2));
  }
  const Domain2 dom{min2.outline, min2.holes, opt.weld_tol};

  // ---- in-plane bar detection (model.bars order) ----
  struct Det {
    Id id = 0;
    Id ni = 0, nj = 0;
    Vec2 a, b;
    std::vector<int> chain;  // mesh vertex indices along the axis
  };
  std::vector<Det> det;
  for (const Bar& b : model.bars) {
    if (b.kind != BarKind::Frame) continue;
    const Vec3& pa = model.node(b.ni).x;
    const Vec3& pb = model.node(b.nj).x;
    if (plane_dev(pa) > opt.weld_tol || plane_dev(pb) > opt.weld_tol) continue;
    Vec2 a2 = to2(pa), b2 = to2(pb);
    if ((b2 - a2).norm() <= opt.weld_tol) continue;  // degenerate in-plane projection
    if (!dom.inside_or_boundary(a2) || !dom.inside_or_boundary(b2)) continue;
    // A boundary-running bar (midpoint on a ring) is excluded from mixing.
    if (!dom.interior_strict(0.5 * (a2 + b2))) continue;
    det.push_back(Det{b.id, b.ni, b.nj, a2, b2, {}});
  }

  // ---- conflict rejection (typed, actionable, BEFORE any mutation) ----
  for (std::size_t i = 0; i < det.size(); ++i)
    for (std::size_t j = i + 1; j < det.size(); ++j)
      if (segments_conflict_tol(det[i].a, det[i].b, det[j].a, det[j].b, opt.weld_tol))
        fail(at + ": bars " + std::to_string(det[i].id) + " and " + std::to_string(det[j].id) +
             " cross inside the slab without a shared node — split them at a common node or "
             "move one out of the slab plane");

  // ---- imposed points: pre-existing in-plane interior nodes (nodes order) ----
  std::vector<Vec2> verts;  // every already-collected input vertex, for dedup
  for (const Vec2& p : min2.outline) verts.push_back(p);
  for (const auto& h : min2.holes)
    for (const Vec2& p : h) verts.push_back(p);
  for (const Det& d : det) {
    verts.push_back(d.a);
    verts.push_back(d.b);
  }
  std::vector<Vec2> extra;
  for (const Node& nd : model.nodes) {
    if (plane_dev(nd.x) > opt.weld_tol) continue;
    Vec2 q = to2(nd.x);
    if (!dom.interior_strict(q)) continue;
    bool dup = false;
    for (const Vec2& v : verts)
      if ((v - q).norm() <= opt.weld_tol) {
        dup = true;
        break;
      }
    if (dup) continue;
    verts.push_back(q);
    extra.push_back(q);
  }

  // ---- triangulate ----
  for (const Det& d : det) min2.segments.emplace_back(d.a, d.b);
  min2.points = std::move(extra);
  min2.target_h = s.mesh_size;
  min2.min_angle_deg = opt.min_angle_deg;
  min2.allow_partial = opt.allow_partial;

  SlabMeshResult res;
  res.slab = slab_id;
  res.origin = origin;
  res.frame = fr;
  try {
    res.mesh = mesh2::triangulate(min2);
  } catch (const Error& e) {
    fail(at + ": " + e.what());
  }
  res.refinement_complete = res.mesh.refinement_complete;
  if (!res.mesh.refinement_complete)
    res.warnings.push_back(
        {SlabMeshWarning::Kind::RefinementIncomplete,
         at + ": the mesher could not reach the " + std::to_string(opt.min_angle_deg) +
             " deg quality policy within its size budget — the mesh is valid but partial "
             "(allow_partial was set)"});

  // ---- 2-D axis chains + split blockers (still BEFORE any mutation) ----
  for (Det& d : det) {
    d.chain = axis_chain(res.mesh, d.a, d.b, at, d.id);
    if (d.chain.size() < 3) continue;  // the bar will stay whole
    const std::string who = at + ": bar " + std::to_string(d.id) +
                            " must be auto-split by the slab mesh but ";
    for (const Tendon& t : model.tendons)
      for (Id bid : t.bars)
        require(bid != d.id, who + "tendon '" + t.name +
                                 "' runs through it — keep the bar out of the slab plane or "
                                 "route the tendon through the split segments afterwards");
    require(model.tapered_of(d.id) == nullptr,
            who + "it carries an exact non-prismatic (CARA VAR) kernel, which cannot be sliced");
    for (const BarMass& bm : model.bar_masses)
      require(bm.bar != d.id,
              who + "it carries a modal-mass entry, which cannot be sliced");
  }

  // ---- raise + weld (mesher node order; candidates = pre-existing nodes) ----
  const std::size_t n_pre = model.nodes.size();
  std::vector<char> claimed(n_pre, 0);
  res.nodes.reserve(res.mesh.nodes.size());
  for (const Vec2& q : res.mesh.nodes) {
    Vec3 p3 = to3(q);
    std::size_t best = n_pre;
    double bestd = 1e300;
    for (std::size_t k = 0; k < n_pre; ++k) {
      double dd = (model.nodes[k].x - p3).norm();
      if (dd < bestd) {  // strict: ties keep the lower index
        bestd = dd;
        best = k;
      }
    }
    if (best < n_pre && bestd <= opt.weld_tol) {
      require(!claimed[best],
              at + ": two mesh nodes weld onto model node " +
                  std::to_string(model.nodes[best].id) +
                  " — the welding tolerance exceeds the local mesh spacing");
      claimed[best] = 1;
      res.nodes.push_back(model.nodes[best].id);
      ++res.welded_nodes;
    } else {
      res.nodes.push_back(model.add_node(p3.x(), p3.y(), p3.z()));
    }
  }

  // ---- shells (mesher triangle order) ----
  res.shells.reserve(res.mesh.tris.size());
  for (const auto& tr : res.mesh.tris)
    res.shells.push_back(model.add_shell({res.nodes[static_cast<std::size_t>(tr[0])],
                                          res.nodes[static_cast<std::size_t>(tr[1])],
                                          res.nodes[static_cast<std::size_t>(tr[2])]},
                                         s.t, s.material, slab_id));

  // ---- auto-split (in-place replacement) + load re-slicing ----
  Id next_bar_id = 0;
  for (const Bar& b : model.bars) next_bar_id = std::max(next_bar_id, b.id);
  ++next_bar_id;

  for (const Det& d : det) {
    SplitBarRecord rec;
    rec.original = d.id;
    rec.chain_nodes.reserve(d.chain.size());
    for (int ix : d.chain) rec.chain_nodes.push_back(res.nodes[static_cast<std::size_t>(ix)]);
    require(rec.chain_nodes.front() == d.ni && rec.chain_nodes.back() == d.nj,
            at + ": the axis of bar " + std::to_string(d.id) +
                " did not weld back onto its own end nodes — near-duplicate nodes within "
                "weld_tol of the bar ends?");
    rec.s.reserve(rec.chain_nodes.size());
    double run = 0;
    rec.s.push_back(0.0);
    for (std::size_t k = 1; k < rec.chain_nodes.size(); ++k) {
      run += (model.node(rec.chain_nodes[k]).x - model.node(rec.chain_nodes[k - 1]).x).norm();
      rec.s.push_back(run);
    }

    if (rec.chain_nodes.size() == 2) {
      rec.segments = {d.id};  // untouched — the mesh put no vertex inside the axis
      res.split_bars.push_back(std::move(rec));
      continue;
    }

    auto pos = std::find_if(model.bars.begin(), model.bars.end(),
                            [&](const Bar& b) { return b.id == d.id; });
    require(pos != model.bars.end(), at + ": internal error — split bar vanished");
    const Bar orig = *pos;
    const double L_tot = rec.s.back();
    std::vector<Bar> segs;
    segs.reserve(rec.chain_nodes.size() - 1);
    for (std::size_t k = 0; k + 1 < rec.chain_nodes.size(); ++k) {
      Bar sb = orig;  // material, section, kind, roll, soil copied
      sb.id = next_bar_id++;
      sb.ni = rec.chain_nodes[k];
      sb.nj = rec.chain_nodes[k + 1];
      // End releases belong to the TRUE ends only.
      sb.releases = Releases{};
      if (k == 0)
        for (int c = 0; c < 6; ++c)
          sb.releases[static_cast<std::size_t>(c)] = orig.releases[static_cast<std::size_t>(c)];
      if (k + 2 == rec.chain_nodes.size())
        for (int c = 6; c < 12; ++c)
          sb.releases[static_cast<std::size_t>(c)] = orig.releases[static_cast<std::size_t>(c)];
      // Rigid-arm eccentricity: linear interpolation along the axis (exact for
      // the constant offsets of a composite girder).
      double t0 = rec.s[k] / L_tot, t1 = rec.s[k + 1] / L_tot;
      sb.ecc[0] = orig.ecc[0] * (1.0 - t0) + orig.ecc[1] * t0;
      sb.ecc[1] = orig.ecc[0] * (1.0 - t1) + orig.ecc[1] * t1;
      segs.push_back(sb);
      rec.segments.push_back(sb.id);
    }
    pos = model.bars.erase(pos);
    model.bars.insert(pos, segs.begin(), segs.end());
    model.invalidate_index();
    res.split_bars.push_back(std::move(rec));
  }

  // Re-slice the existing bar loads of every split bar (ordered map: lookups
  // only, but determinism costs nothing).
  std::map<Id, const SplitBarRecord*> split_of;
  for (const SplitBarRecord& r : res.split_bars)
    if (r.segments.size() >= 2) split_of[r.original] = &r;
  if (!split_of.empty()) {
    for (LoadCase& lc : model.cases) {
      std::vector<Load> out;
      out.reserve(lc.loads.size());
      for (const Load& l : lc.loads) {
        bool bar_load = l.kind == Load::Kind::BarUdl || l.kind == Load::Kind::BarLinear ||
                        l.kind == Load::Kind::BarPoint || l.kind == Load::Kind::BarThermal ||
                        l.kind == Load::Kind::Prestrain ||
                        l.kind == Load::Kind::PrestrainPoint ||
                        l.kind == Load::Kind::PrestrainLinear || l.kind == Load::Kind::BarParabolic ||
                        l.kind == Load::Kind::PrestrainParabolic;
        auto it = bar_load ? split_of.find(l.bar) : split_of.end();
        if (it == split_of.end()) {
          out.push_back(l);
          continue;
        }
        slice_bar_load(l, *it->second, lc.name, out);
      }
      lc.loads = std::move(out);
    }
  }

  return res;
}

}  // namespace keta

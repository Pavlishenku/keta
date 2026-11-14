// Analyze components (Assemble, Analyze, Combination) and the registry,
// lowered onto keta::Model and keta::solve_static_cached.
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <map>

#include "keta/param/components.h"
#include "keta/prestress.h"
#include "keta/slabmesh.h"

namespace keta::param {

namespace {

std::string fmt3(const Vec3& p) {
  char buf[96];
  std::snprintf(buf, sizeof buf, "(%.3f, %.3f, %.3f)", p.x(), p.y(), p.z());
  return buf;
}

double dist_point_segment(const Vec3& p, const Vec3& a, const Vec3& b) {
  Vec3 ab = b - a;
  double l2 = ab.squaredNorm();
  if (l2 == 0) return (p - a).norm();
  double t = std::clamp((p - a).dot(ab) / l2, 0.0, 1.0);
  return (p - (a + t * ab)).norm();
}

template <class T>
std::vector<const T*> collect(const Tree& t) {
  std::vector<const T*> out;
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (const T* x = o->get<T>()) out.push_back(x);
  return out;
}

bool same_support(const Support& a, const Support& b) {
  return a.node == b.node && a.fixed == b.fixed && a.spring == b.spring && a.oneside == b.oneside;
}

}  // namespace

// ============================ Assemble ============================

ComponentInfo AssembleNode::info() const {
  return {"keta.assemble", "Assemble", "Analyze", "Weld beams, supports and loads into a model"};
}
std::vector<PinSpec> AssembleNode::inputs() const {
  return {PinSpec::list("B", "beam", "Elements").opt(),
          PinSpec::list("S", "support", "Supports").opt(),
          PinSpec::list("L", "load", "Loads").opt(),
          PinSpec::list("K", "cable", "Cable elements").opt(),
          PinSpec::list("P", "pointmass", "Point masses").opt(),
          PinSpec::list("C", "combo", "Load combinations").opt(),
          PinSpec::list("T", "tendon", "Prestressing tendons").opt(),
          PinSpec::list("G", "foundation", "Winkler foundations (soil beds)").opt(),
          PinSpec::list("PL", "slab", "Slabs (meshed into shells)").opt()};
}
std::vector<PinSpec> AssembleNode::outputs() const {
  return {PinSpec::item("M", "model", "The assembled model"),
          PinSpec::list("W", "text", "One warning per meshing issue (empty when clean)")};
}

Outputs<Value> AssembleNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  const auto beams = collect<BeamSpec>(input_or_empty(in, 0));
  const auto supports = collect<SupportSpec>(input_or_empty(in, 1));
  const auto loads = collect<LoadSpec>(input_or_empty(in, 2));
  const auto masses = collect<PointMassSpec>(input_or_empty(in, 4));
  const auto combos = collect<ComboSpec>(input_or_empty(in, 5));
  const auto foundations = collect<FoundationSpec>(input_or_empty(in, 7));
  const auto tendons = collect<TendonSpec>(input_or_empty(in, 6));
  const auto slab_specs = collect<SlabSpec>(input_or_empty(in, 8));
  if (beams.empty() && slab_specs.empty())
    throw ComponentError("a structure needs beams or slabs — wire elements into B or slabs into PL");
  // Two slabs with EQUAL outlines would mesh onto identical welded nodes
  // (doubled stiffness) while a surface load resolves to the first only.
  for (std::size_t i = 0; i < slab_specs.size(); ++i)
    for (std::size_t j = i + 1; j < slab_specs.size(); ++j)
      if (slab_specs[i]->outline == slab_specs[j]->outline)
        throw ComponentError("slabs " + std::to_string(i + 1) + " and " + std::to_string(j + 1) +
                             " share an outline — surface-load targeting would be ambiguous and stiffness would "
                             "double; offset or merge them");
  // With slabs, everything positioned on NODES (supports, point masses, nodal
  // loads, settlements) and the surface loads resolve AFTER meshing, so they
  // see the mesh nodes; bar loads and tendons resolve before it and the slab
  // mesher re-slices them onto the auto-split segments.
  const bool has_slabs = !slab_specs.empty();
  bool meshed = false;
  if (ctx.cancelled()) throw ComponentError::cancelled();

  const double tol = tolerance;
  NodeMerger merger(tol);
  Model m;

  // Bars, with sections and materials deduplicated by value (ids renumbered).
  for (std::size_t bi = 0; bi < beams.size(); ++bi) {
    const BeamSpec& b = *beams[bi];
    Id n1 = merger.id_for(b.line.a), n2 = merger.id_for(b.line.b);
    if (n1 == n2)
      throw ComponentError("beam " + std::to_string(bi + 1) + " has zero length within tolerance " +
                           std::to_string(tol));
    auto sec = std::find_if(m.sections.begin(), m.sections.end(), [&](const Section& s) {
      Section probe = b.section.s;
      probe.id = s.id;
      // a section whose shear areas were resolved for its material (energy,
      // Model::resolve_shear) is the same catalogue entry
      if (!std::isnan(s.shear_nu) && probe.props.Avy == s.avy_shape && probe.props.Avz == s.avz_shape) {
        probe.props.Avy = s.props.Avy;
        probe.props.Avz = s.props.Avz;
      }
      return sections_equal(s, probe);
    });
    Id sec_id;
    if (sec != m.sections.end()) {
      sec_id = sec->id;
    } else {
      Section s = b.section.s;
      s.id = static_cast<Id>(m.sections.size() + 1);
      m.sections.push_back(s);
      sec_id = s.id;
    }
    auto mat = std::find_if(m.materials.begin(), m.materials.end(),
                            [&](const Material& x) { return materials_equal(x, b.material.m); });
    Id mat_id;
    if (mat != m.materials.end()) {
      mat_id = mat->id;
    } else {
      Material x = b.material.m;
      x.id = static_cast<Id>(m.materials.size() + 1);
      m.materials.push_back(x);
      mat_id = x.id;
    }
    Bar bar;
    bar.id = static_cast<Id>(bi + 1);
    bar.ni = n1;
    bar.nj = n2;
    bar.section = sec_id;
    bar.material = mat_id;
    bar.kind = b.truss ? BarKind::Truss : BarKind::Frame;
    bar.roll = b.roll;
    bar.releases = b.releases;
    if (b.ecc) bar.ecc = {*b.ecc, *b.ecc};
    m.bars.push_back(bar);
    m.resolve_shear(m.bars.back());
  }

  // Cables (pin K): Cable-kind bars whose anchors weld like beam ends, ids
  // after the beams (error messages name a cable by its 1-based index among the
  // cables). EA, w and the L0 route ride on CableProps; keta's element
  // constructor still wants a material and a section, so all cables share one
  // inert placeholder pair (E = 1 Pa, 1 m x 1 m, rho = 0) that the cable driver
  // never reads (the CableProps overrides win).
  const auto cables = collect<CableSpec>(input_or_empty(in, 3));
  if (!cables.empty()) {
    Material pm;
    pm.id = static_cast<Id>(m.materials.size() + 1);
    pm.name = "cable (placeholder)";
    pm.E = 1.0;
    pm.nu = 0.3;
    pm.rho = 0.0;
    m.materials.push_back(pm);
    Section ps;
    ps.id = static_cast<Id>(m.sections.size() + 1);
    ps.name = "cable (placeholder)";
    ps.shape = ShapeRect{1.0, 1.0};
    ps.finalize();
    m.sections.push_back(ps);
    for (std::size_t ci = 0; ci < cables.size(); ++ci) {
      const CableSpec& c = *cables[ci];
      const std::string who = "invalid model: cable " + std::to_string(ci + 1);
      auto positive = [&](double v, const char* what) {
        if (!(v > 0)) throw ComponentError(who + ": " + what + " must be > 0 (got " + Value::number(v).describe() + ")");
      };
      positive(c.ea, "axial rigidity EA");
      if (!(c.w >= 0))
        throw ComponentError(who + ": weight per length must be >= 0 (got " + Value::number(c.w).describe() + ")");
      positive(c.value, c.route == CableRoute::Length               ? "unstrained length L0"
                        : c.route == CableRoute::HorizontalTension ? "horizontal tension H"
                                                                   : "target end tension");
      Id n1 = merger.id_for(c.a), n2 = merger.id_for(c.b);
      if (n1 == n2)
        throw ComponentError("cable " + std::to_string(ci + 1) + " has zero length within tolerance " +
                             Value::number(tol).describe());
      Bar bar;
      bar.id = static_cast<Id>(m.bars.size() + 1);
      bar.ni = n1;
      bar.nj = n2;
      bar.material = pm.id;
      bar.section = ps.id;
      bar.kind = BarKind::Cable;
      CableProps cp;
      cp.ea = c.ea;
      cp.w = c.w;
      switch (c.route) {
        case CableRoute::Length: cp.unstressed_length = c.value; break;
        case CableRoute::HorizontalTension: cp.pretension = c.value; break;
        case CableRoute::EndTension: cp.end_tension = c.value; break;
      }
      bar.cable = cp;
      m.bars.push_back(bar);
      m.resolve_shear(m.bars.back());
    }
  }

  // Nodes (every node is created by a bar end: supports and loads never spawn one).
  const auto& pos = merger.positions();
  for (std::size_t i = 0; i < pos.size(); ++i) m.nodes.push_back({static_cast<Id>(i + 1), pos[i]});

  auto bar_for = [&](const Line& l, const std::string& what) {
    auto hit = find_bar(m, l, tol);
    if (!hit) throw ComponentError(what + ": no bar matches the loaded line");
    return *hit;
  };
  auto node_at = [&](const Vec3& p, const std::string& what) {
    std::optional<Id> id;
    if (!meshed) {
      id = merger.find(p);
    } else {  // the final node set (mesh nodes included): nearest within tol
      double best = tol;
      for (const keta::Node& n : m.nodes)
        if (const double dd = (n.x - p).norm(); dd <= best) {
          best = dd;
          id = n.id;
        }
    }
    if (!id) throw ComponentError(what + " at " + fmt3(p) + " lands on no node (tolerance " + std::to_string(tol) + ")");
    return *id;
  };

  // Winkler beds: a property of the matched bars; conflicting moduli are typed.
  for (std::size_t fi = 0; fi < foundations.size(); ++fi) {
    const FoundationSpec& f = *foundations[fi];
    SoilProps soil{f.kfx, f.kfy, f.kfz, f.kmx};
    for (const Line& l : f.lines) {
      auto hit = find_bar(m, l, tol);
      if (!hit) throw ComponentError("foundation " + std::to_string(fi + 1) + ": no bar matches the founded line");
      SoilProps& cur = m.bars[hit->index].soil;
      if (cur.any() && (cur.kfx != soil.kfx || cur.kfy != soil.kfy || cur.kfz != soil.kfz || cur.kmx != soil.kmx))
        throw ComponentError("foundation " + std::to_string(fi + 1) + ": bar " + std::to_string(hit->index + 1) +
                             " already carries different soil moduli (conflicting foundations on one bar)");
      cur = soil;
    }
  }

  auto resolve_supports_and_masses = [&] {
  // Supports: point supports land on the nearest node; a line support applies
  // its mask to EVERY node within tol of the segment. Exact duplicates collapse.
  // The supports of each node (positions in m.supports): the duplicate test
  // only meets the supports of its node, not every support (O(S^2) on a long
  // line support).
  std::unordered_map<Id, std::vector<std::size_t>> sup_at;
  for (std::size_t q = 0; q < m.supports.size(); ++q) sup_at[m.supports[q].node].push_back(q);
  for (std::size_t si = 0; si < supports.size(); ++si) {
    const SupportSpec& s = *supports[si];
    Support base;
    base.fixed = s.blocked;
    for (std::size_t k = 0; k < 6; ++k)
      if (s.springs[k] && !s.blocked[k]) base.spring[k] = *s.springs[k];
    // A support that blocks nothing and carries no spring is mechanically
    // inert and keta's model validation rejects it, so it is located (errors
    // still surface) but omitted.
    const bool inert = std::none_of(base.fixed.begin(), base.fixed.end(), [](bool b) { return b; }) &&
                       std::all_of(base.spring.begin(), base.spring.end(), [](double k) { return k == 0.0; });
    auto add = [&](Id node) {
      if (inert) return;
      Support x = base;
      x.node = node;
      std::vector<std::size_t>& here = sup_at[x.node];
      for (std::size_t q : here)
        if (same_support(m.supports[q], x)) return;
      here.push_back(m.supports.size());
      m.supports.push_back(x);
    };
    if (!s.line) {
      add(node_at(s.point, "support " + std::to_string(si + 1)));
      continue;
    }
    const Line& l = *s.line;
    if (!(std::isfinite(l.length()) && l.length() > tol))
      throw ComponentError("support " + std::to_string(si + 1) +
                           " (line): zero-length line segment — use a point support instead");
    std::size_t matched = 0;
    for (const keta::Node& n : m.nodes)
      if (dist_point_segment(n.x, l.a, l.b) <= tol) {
        add(n.id);
        ++matched;
      }
    if (matched == 0)
      throw ComponentError("support " + std::to_string(si + 1) + " (line): no model node within tolerance of the segment");
  }

  for (std::size_t pi = 0; pi < masses.size(); ++pi) {
    PointMass pm;
    pm.node = node_at(masses[pi]->point, "point mass " + std::to_string(pi + 1));
    pm.mass = masses[pi]->mass;
    pm.inertia = masses[pi]->inertia;
    m.point_masses.push_back(pm);
  }

  };
  if (!has_slabs) resolve_supports_and_masses();

  // Load cases: grouped by case id (sorted), an empty default case 1 if none.
  std::vector<Id> case_ids;
  for (const LoadSpec* l : loads) case_ids.push_back(l->case_id);
  for (const TendonSpec* t : tendons) case_ids.push_back(t->case_id);  // a tendon-only case exists too
  std::sort(case_ids.begin(), case_ids.end());
  case_ids.erase(std::unique(case_ids.begin(), case_ids.end()), case_ids.end());
  if (case_ids.empty()) case_ids.push_back(1);
  std::map<Id, std::size_t> case_pos;
  for (Id id : case_ids) {
    case_pos[id] = m.cases.size();
    m.cases.push_back({id, "LC" + std::to_string(id), {}});
  }
  std::map<Id, std::string> labels;

  // after_mesh: the loads positioned on nodes or slabs, resolved once the
  // mesh exists (all of them in one pass when there is no slab).
  auto after_mesh = [&](const LoadSpec& spec) {
    return std::holds_alternative<SurfacePressureLoad>(spec.kind) ||
           std::holds_alternative<SurfaceThermalLoad>(spec.kind) ||
           (has_slabs && (std::holds_alternative<NodalLoad>(spec.kind) || std::holds_alternative<SettlementLoad>(spec.kind)));
  };
  auto slab_of = [&](const std::vector<Vec3>& outline) -> std::optional<Id> {
    for (std::size_t k = 0; k < slab_specs.size(); ++k)
      if (slab_specs[k]->outline == outline) return static_cast<Id>(k + 1);
    return std::nullopt;
  };
  auto resolve_loads = [&](bool late) {
  for (std::size_t li = 0; li < loads.size(); ++li) {
    const LoadSpec& spec = *loads[li];
    if (after_mesh(spec) != late) continue;
    LoadCase& lc = m.cases[case_pos.at(spec.case_id)];
    const std::string tag = "load " + std::to_string(li + 1);
    std::visit(
        [&](const auto& l) {
          using T = std::decay_t<decltype(l)>;
          if constexpr (std::is_same_v<T, NodalLoad>) {
            lc.loads.push_back(Load::nodal(node_at(l.point, "point " + tag), l.F, l.M));
          } else if constexpr (std::is_same_v<T, UdlLoad>) {
            for (const Line& line : l.lines) {
              BarMatch hit = bar_for(line, "line " + tag);
              lc.loads.push_back(Load::udl(m.bars[hit.index].id, l.w, l.frame));
            }
          } else if constexpr (std::is_same_v<T, GravityLoad>) {
            lc.loads.push_back(Load::self_weight(l.g));
          } else if constexpr (std::is_same_v<T, LinearLoad>) {
            for (const Line& line : l.lines) {
              BarMatch hit = bar_for(line, "linear " + tag);
              // Positions are authored along the LINE; mirror them onto a bar
              // that runs the other way.
              Load x = hit.reversed ? Load::linear(m.bars[hit.index].id, l.w2, l.w1, l.frame, 1 - l.x2, 1 - l.x1)
                                    : Load::linear(m.bars[hit.index].id, l.w1, l.w2, l.frame, l.x1, l.x2);
              x.relative = true;
              lc.loads.push_back(x);
            }
          } else if constexpr (std::is_same_v<T, BarPointLoad>) {
            for (const Line& line : l.lines) {
              BarMatch hit = bar_for(line, "bar point " + tag);
              Load x = Load::point(m.bars[hit.index].id, hit.reversed ? 1 - l.at : l.at, l.F, Vec3::Zero(), l.frame);
              x.relative = true;
              lc.loads.push_back(x);
            }
          } else if constexpr (std::is_same_v<T, ThermalLoad>) {
            for (const Line& line : l.lines) {
              BarMatch hit = bar_for(line, "thermal " + tag);
              const Bar& bar = m.bars[hit.index];
              if (m.materials[static_cast<std::size_t>(bar.material - 1)].alpha == 0.0)
                throw ComponentError("thermal " + tag + ": bar material has no thermal dilation coefficient (alpha)");
              lc.loads.push_back(Load::thermal(bar.id, l.dt, l.gy, l.gz));
            }
          } else if constexpr (std::is_same_v<T, SettlementLoad>) {
            lc.loads.push_back(Load::support_disp(node_at(l.point, "settlement " + tag), l.d));
          } else if constexpr (std::is_same_v<T, CaseLabel>) {
            auto it = labels.find(spec.case_id);
            if (it != labels.end() && it->second != l.name)
              throw ComponentError("case " + std::to_string(spec.case_id) + ": conflicting labels \"" + it->second +
                                   "\" and \"" + l.name + "\"");
            labels[spec.case_id] = l.name;
            lc.name = l.name;
          } else if constexpr (std::is_same_v<T, SurfacePressureLoad>) {
            const auto sid = slab_of(l.outline);
            if (!sid)
              throw ComponentError("surface load " + std::to_string(li + 1) +
                                   ": target slab is not on Assemble's PL pin (wire the SAME slab into both the "
                                   "Surface Load and Assemble)");
            lc.loads.push_back(Load::surface_pressure(*sid, l.w, l.frame));
          } else if constexpr (std::is_same_v<T, SurfaceThermalLoad>) {
            const auto sid = slab_of(l.outline);
            if (!sid)
              throw ComponentError("surface thermal load " + std::to_string(li + 1) +
                                   ": target slab is not on Assemble's PL pin (wire the SAME slab into both the "
                                   "Surface Thermal and Assemble)");
            if (m.material(m.slab(*sid).material).alpha == 0.0)
              throw ComponentError("surface thermal load " + std::to_string(li + 1) +
                                   ": slab material has no thermal dilation coefficient (alpha_temp)");
            lc.loads.push_back(Load::surface_thermal(*sid, l.dt, l.grad));
          }
        },
        spec.kind);
  }
  };
  resolve_loads(false);

  // Tendons: resolve each chain to ordered (bar, flip) pairs, record the
  // tendon and lower its equivalent loads into its case.
  for (std::size_t ti = 0; ti < tendons.size(); ++ti) {
    const TendonSpec& spec = *tendons[ti];
    const std::string tn = "tendon " + spec.name.value_or(std::to_string(ti + 1));
    const std::vector<Line>& lines = spec.lines;
    if (lines.empty()) throw ComponentError(tn + ": chain has no lines");
    for (std::size_t i = 0; i < lines.size(); ++i)
      if (!(std::isfinite(lines[i].length()) && lines[i].length() > 0))
        throw ComponentError(tn + ": chain line " + std::to_string(i + 1) + " has zero or non-finite length");
    Tendon t = tendon_record(spec);
    t.id = static_cast<Id>(ti + 1);
    std::vector<std::size_t> bar_ix;
    for (std::size_t i = 0; i < lines.size(); ++i) {
      auto hit = find_bar(m, lines[i], tol);
      if (!hit)
        throw ComponentError(tn + ": chain line " + std::to_string(i + 1) + " matches no bar in the model (tolerance " +
                             Value::number(tol).describe() + ")");
      bar_ix.push_back(hit->index);
    }
    // Travel start of each line by a contiguity walk.
    auto close = [tol](const Vec3& a, const Vec3& b) { return (a - b).norm() <= tol; };
    std::vector<Vec3> starts;
    if (lines.size() == 1) {
      starts.push_back(lines[0].a);
    } else {
      const Line &l0 = lines[0], &l1 = lines[1];
      Vec3 shared;
      if (close(l0.a, l1.a) || close(l0.a, l1.b)) shared = l0.a;
      else if (close(l0.b, l1.a) || close(l0.b, l1.b)) shared = l0.b;
      else
        throw ComponentError(tn + ": chain lines 1 and 2 do not share an endpoint (within tolerance " +
                             Value::number(tol).describe() + ")");
      starts.push_back(close(shared, l0.a) ? l0.b : l0.a);
      Vec3 prev = shared;
      for (std::size_t i = 1; i < lines.size(); ++i) {
        if (close(lines[i].a, prev)) {
          starts.push_back(lines[i].a);
          prev = lines[i].b;
        } else if (close(lines[i].b, prev)) {
          starts.push_back(lines[i].b);
          prev = lines[i].a;
        } else {
          throw ComponentError(tn + ": chain lines " + std::to_string(i) + " and " + std::to_string(i + 1) +
                               " do not share an endpoint (within tolerance " + Value::number(tol).describe() + ")");
        }
      }
    }
    for (std::size_t i = 0; i < bar_ix.size(); ++i) {
      const Bar& b = m.bars[bar_ix[i]];
      const Vec3& p1 = m.nodes[static_cast<std::size_t>(b.ni - 1)].x;
      const Vec3& p2 = m.nodes[static_cast<std::size_t>(b.nj - 1)].x;
      if (!close(starts[i], p1) && !close(starts[i], p2))
        throw ComponentError(tn + ": chain line " + std::to_string(i + 1) + " does not align with bar " +
                             std::to_string(b.id));
      t.bars.push_back(b.id);
      t.flips.push_back(!close(starts[i], p1));
    }
    if (std::none_of(t.flips.begin(), t.flips.end(), [](bool f) { return f; })) t.flips.clear();
    m.tendons.push_back(t);
    m.invalidate_index();
    std::vector<Load> eq;
    try {
      eq = tendon_loads(m, m.tendons.back());
    } catch (const Error& e) {
      throw ComponentError(tn + ": lowering onto its bars failed: " + e.what());
    }
    LoadCase& lc = m.cases[case_pos.at(spec.case_id)];
    lc.loads.insert(lc.loads.end(), eq.begin(), eq.end());
  }

  // Slabs (pin PL): each meshed by keta's slab mesher (conforming to the
  // in-plane bars, which it auto-splits and whose loads it re-slices), welded
  // onto the existing nodes within the Assemble tolerance. An incomplete
  // refinement is a W warning, never silent.
  Tree warnings;
  for (std::size_t si = 0; si < slab_specs.size(); ++si) {
    const SlabSpec& spec = *slab_specs[si];
    auto mat = std::find_if(m.materials.begin(), m.materials.end(),
                            [&](const Material& x) { return materials_equal(x, spec.material.m); });
    Id mat_id;
    if (mat != m.materials.end()) {
      mat_id = mat->id;
    } else {
      Material x = spec.material.m;
      x.id = static_cast<Id>(m.materials.size() + 1);
      m.materials.push_back(x);
      mat_id = x.id;
    }
    Slab sl;
    sl.id = static_cast<Id>(si + 1);
    sl.outline = spec.outline;
    sl.holes = spec.holes;
    sl.t = spec.thickness;
    sl.material = mat_id;
    sl.mesh_size = spec.mesh_size;
    sl.local_x = spec.local_x;
    m.slabs.push_back(sl);
    m.invalidate_index();
    if (ctx.cancelled()) throw ComponentError::cancelled();
    SlabMeshOptions so;
    so.weld_tol = tol;
    so.allow_partial = true;
    try {
      const SlabMeshResult r = mesh_slab(m, sl.id, so);
      for (const SlabMeshWarning& w : r.warnings)
        warnings.append(Path::single(0), Value::text("slab " + std::to_string(si + 1) + ": " + w.message));
    } catch (const Error& e) {
      throw ComponentError("slab " + std::to_string(si + 1) + ": meshing failed: " + e.what());
    }
  }
  if (has_slabs) {
    meshed = true;
    m.invalidate_index();
    resolve_supports_and_masses();
  }
  resolve_loads(true);

  for (std::size_t ci = 0; ci < combos.size(); ++ci) {
    const ComboSpec& c = *combos[ci];
    std::vector<std::pair<Id, double>> terms;
    for (const auto& [factor, cid] : c.terms) {
      if (!case_pos.count(cid))
        throw ComponentError("combination " + std::to_string(ci + 1) + " references unknown load case " +
                             std::to_string(cid));
      terms.emplace_back(cid, factor);
    }
    m.add_combo(c.name.value_or("C" + std::to_string(ci + 1)), std::move(terms));
  }

  m.invalidate_index();
  try {
    m.validate();
  } catch (const Error& e) {
    throw ComponentError(std::string("invalid model: ") + e.what());
  }
  return {single(Value::of(ModelRef(std::make_shared<const Model>(std::move(m))))), warnings};
}

json::Value AssembleNode::params() const {
  json::Value p = json::Value::object();
  p.set("tolerance", tolerance);
  return p;
}
void AssembleNode::set_params(const json::Value& p) {
  if (const auto* t = p.get("tolerance"); t && t->is_number()) tolerance = t->as_number();
}

// ============================ Analyze ============================

ComponentInfo AnalyzeNode::info() const { return {"keta.analyze", "Analyze", "Analyze", "Linear elastic static analysis"}; }
std::vector<PinSpec> AnalyzeNode::inputs() const { return {PinSpec::item("M", "model", "The model")}; }
std::vector<PinSpec> AnalyzeNode::outputs() const { return {PinSpec::item("R", "results", "Model + results")}; }

FactorCache::Stats AnalyzeNode::cache_stats() const { return cache->stats(); }

Outputs<Value> AnalyzeNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model;
  for (const auto& kv : in[0].branches())
    for (const auto& o : kv.second)
      if (o && !model)
        if (const auto* r = o->get<ModelRef>()) model = *r;
  if (!model) throw ComponentError::missing("M");
  if (ctx.cancelled()) throw ComponentError::cancelled();
  StaticOptions opt;
  opt.stations_per_bar = stations;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  auto out = std::make_shared<AnalyzedModel>();
  out->model = model;
  {
    std::lock_guard<std::mutex> lk(cache->mu);
    try {
      out->results = solve_static_cached(*model, opt, cache->cache);
    } catch (const Cancelled&) {
      throw ComponentError::cancelled();  // polls sit outside the factorization: the cache stays coherent
    } catch (const std::exception& e) {
      cache->cache.invalidate();
      throw ComponentError(std::string("analysis failed: ") + e.what());
    }
  }
  // Section-point stresses, on demand. They never fail: a section whose
  // exact solve fails degrades to its stress points, listed with the reason.
  out->stress_source = std::make_shared<const LazySectionStress>(stress_cache);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(ResultsRef(std::move(out))))};
}

std::shared_ptr<const SectionStressBundle> section_stress_bundle(const Model& model,
                                                                 const StaticResults& results,
                                                                 SharedSectionCache& cache,
                                                                 const SolveHooks* hooks) {
  auto b = std::make_shared<SectionStressBundle>();
  {
    std::lock_guard<std::mutex> lk(cache.mu);
    b->engine = std::make_shared<const SectionStressEngine>(
        SectionStressEngine::build(model, MeshResolution::Standard, cache.cache, hooks));
  }
  b->results = recover_stress_results(model, results, *b->engine);
  return b;
}

json::Value AnalyzeNode::params() const {
  json::Value p = json::Value::object();
  p.set("stations", stations);
  return p;
}
void AnalyzeNode::set_params(const json::Value& p) {
  if (const auto* s = p.get("stations"); s && s->is_number() && s->as_int() >= 2)
    stations = static_cast<int>(s->as_int());
}

// ============================ Combination ============================

ComponentInfo CombinationNode::info() const {
  return {"keta.combination", "Combination", "Analyze", "Linear load combination sum factor*case"};
}
std::vector<PinSpec> CombinationNode::inputs() const {
  return {PinSpec::list("C", "int", "Load-case ids"), PinSpec::list("F", "number", "Factors (paired with cases; last repeats)"),
          PinSpec::item("N", "text", "Combination name").opt()};
}
std::vector<PinSpec> CombinationNode::outputs() const { return {PinSpec::item("CB", "combo", "The combination")}; }
Outputs<Value> CombinationNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::vector<Id> cases;
  std::vector<double> factors;
  for (const Value& v : in[0].non_nulls())
    if (auto i = v.as_int()) cases.push_back(*i);
  for (const Value& v : in[1].non_nulls())
    if (auto x = v.as_number()) factors.push_back(*x);
  if (cases.empty()) throw ComponentError::missing("C");
  if (factors.empty()) throw ComponentError::missing("F");
  ComboSpec c;
  for (std::size_t i = 0; i < cases.size(); ++i)
    c.terms.emplace_back(i < factors.size() ? factors[i] : factors.back(), cases[i]);
  c.name = first_text(input_or_empty(in, 2));
  return {single(Value::of(c))};
}

// ============================ Registry ============================

namespace {
// Retired type id -> current id (one hop; targets must be live ids).
const std::pair<const char*, const char*> kAliases[] = {
    {"keta.selectbysection", "keta.select_by_section"},
    {"keta.selectbymaterial", "keta.select_by_material"},
};

template <class T>
ComponentFactory factory() {
  return {T().info(), [] { return std::unique_ptr<Comp>(std::make_unique<T>()); }};
}
// Op-parametrized node families (one component type id per op).
template <class T, class Op>
ComponentFactory factory(Op op) {
  return {T(op).info(), [op] { return std::unique_ptr<Comp>(std::make_unique<T>(op)); }};
}
}  // namespace

const ComponentFactory* Registry::by_type_id(const std::string& id) const {
  for (const auto& f : factories_)
    if (f.info.type_id == id) return &f;
  return nullptr;
}
const ComponentFactory* Registry::resolve_type_id(const std::string& id) const {
  if (const auto* f = by_type_id(id)) return f;
  for (const auto& [old_id, new_id] : kAliases)
    if (id == old_id) return by_type_id(new_id);
  return nullptr;
}
std::unique_ptr<Comp> Registry::make(const std::string& id) const {
  const auto* f = resolve_type_id(id);
  return f ? f->make() : nullptr;
}
std::vector<std::string> Registry::categories() const {
  std::vector<std::string> cats;
  for (const auto& f : factories_)
    if (std::find(cats.begin(), cats.end(), f.info.category) == cats.end()) cats.push_back(f.info.category);
  return cats;
}

Registry standard_registry() {
  Registry r;
  r.push(factory<NumberSlider>());
  r.push(factory<IntSlider>());
  r.push(factory<Panel>());
  r.push(factory<BakedDataNode>());
  r.push(factory<ClusterNode>());
  r.push(factory<ClusterInNode>());
  r.push(factory<ClusterOutNode>());
  // Tree: Grasshopper "Sets -> Tree" restructuring ops (value-agnostic).
  r.push(factory<GraftNode>());
  r.push(factory<FlattenNode>());
  r.push(factory<SimplifyNode>());
  r.push(factory<TrimTreeNode>());
  r.push(factory<FlipMatrixNode>());
  r.push(factory<PruneTreeNode>());
  r.push(factory<ShiftPathsNode>());
  r.push(factory<MergeTreesNode>());
  r.push(factory<EntwineNode>());
  r.push(factory<CleanTreeNode>());
  r.push(factory<ParamViewerNode>());
  r.push(factory<DataDamNode>());
  for (BinaryOp op : {BinaryOp::Add, BinaryOp::Subtract, BinaryOp::Multiply, BinaryOp::Divide, BinaryOp::Power,
                      BinaryOp::Modulo})
    r.push(factory<BinaryMathNode>(op));
  for (UnaryOp op : {UnaryOp::Negate, UnaryOp::Abs, UnaryOp::Sqrt, UnaryOp::Sin, UnaryOp::Cos, UnaryOp::Tan,
                     UnaryOp::Floor, UnaryOp::Ceil, UnaryOp::Round})
    r.push(factory<UnaryMathNode>(op));
  for (BinaryOp op : {BinaryOp::Larger, BinaryOp::Smaller, BinaryOp::And, BinaryOp::Or}) r.push(factory<BinaryMathNode>(op));
  r.push(factory<NotNode>());
  r.push(factory<EqualsNode>());
  r.push(factory<RangeNode>());
  r.push(factory<RemapNode>());
  r.push(factory<MassAdditionNode>());
  r.push(factory<ExpressionNode>());
  // Lists: the Grasshopper "Sets -> List" family (per-branch).
  r.push(factory<ListItemNode>());
  r.push(factory<ListLengthNode>());
  r.push(factory<SubListNode>());
  r.push(factory<ShiftListNode>());
  r.push(factory<ListReverseNode>());
  r.push(factory<ListInsertNode>());
  r.push(factory<SortListNode>());
  r.push(factory<CullPatternNode>());
  r.push(factory<DispatchNode>());
  // Sets: selection, grouping, set algebra (per-branch; flat-index sets).
  r.push(factory<CullIndexNode>());
  r.push(factory<CullDuplicatesNode>());
  r.push(factory<SortPointsNode>());
  r.push(factory<SelectBoxNode>());
  r.push(factory<SelectLevelNode>());
  r.push(factory<SelectLengthNode>());
  r.push(factory<SelectOrientationNode>());
  r.push(factory<SelectProximityNode>());
  r.push(factory<SelectBySectionNode>());
  r.push(factory<SelectByMaterialNode>());
  r.push(factory<GroupByNode>());
  r.push(factory<CreateSetNode>());
  r.push(factory<SetOpNode>());
  r.push(factory<SetMembersNode>());
  r.push(factory<FilterBySetNode>());
  r.push(factory<SelectionNode>());
  r.push(factory<PointXyz>());
  r.push(factory<VectorXyz>());
  r.push(factory<LineAB>());
  r.push(factory<Series>());
  r.push(factory<FrameGrid>());
  r.push(factory<GridNode>());
  r.push(factory<PolarGridNode>());
  r.push(factory<ConnectNode>());
  r.push(factory<WeldNode>());
  r.push(factory<IntersectNode>());
  r.push(factory<DivideLineNode>());
  r.push(factory<ShatterNode>());
  r.push(factory<CullShortNode>());
  r.push(factory<EndPointsNode>());
  r.push(factory<CsvPointsNode>());
  r.push(factory<CsvBarsNode>());
  r.push(factory<MaterialNode>());
  r.push(factory<SectionCatalogNode>());
  r.push(factory<SectionRect>());
  r.push(factory<SectionExactNode>());
  r.push(factory<SectionEditorNode>());
  r.push(factory<SectionCompositeNode>());
  r.push(factory<BeamNode>());
  r.push(factory<CableNode>());
  r.push(factory<SlabNode>());
  r.push(factory<SurfaceLoadNode>());
  r.push(factory<SurfaceThermalNode>());
  r.push(factory<SupportNode>());
  r.push(factory<PointMassNode>());
  r.push(factory<PointLoadNode>());
  r.push(factory<LineLoadNode>());
  r.push(factory<GravityNode>());
  r.push(factory<LinearLoadNode>());
  r.push(factory<BarPointLoadNode>());
  r.push(factory<ThermalLoadNode>());
  r.push(factory<SettlementNode>());
  r.push(factory<FoundationNode>());
  r.push(factory<CaseLabelNode>());
  r.push(factory<TendonNode>());
  r.push(factory<PhaseNode>());
  r.push(factory<AssembleNode>());
  r.push(factory<AnalyzeNode>());
  r.push(factory<PhasedAnalyzeNode>());
  r.push(factory<CombinationNode>());
  r.push(factory<ComboFamilyNode>());
  r.push(factory<EnvelopeNode>());
  r.push(factory<LoadPathNode>());
  r.push(factory<InfluenceLineNode>());
  r.push(factory<ConvoyNode>());
  r.push(factory<MovingLoadNode>());
  r.push(factory<IlExtremesNode>());
  r.push(factory<InfluenceSurfaceNode>());
  r.push(factory<BucklingNode>());
  r.push(factory<SecondOrderNode>());
  r.push(factory<ModalNode>());
  r.push(factory<ArcLengthNode>());
  r.push(factory<CableStaticNode>());
  r.push(factory<HarmonicNode>());
  r.push(factory<ImperfectionNode>());
  r.push(factory<NonlinearNode>());
  r.push(factory<PsdNode>());
  r.push(factory<SpectrumNode>());
  return r;
}

}  // namespace keta::param

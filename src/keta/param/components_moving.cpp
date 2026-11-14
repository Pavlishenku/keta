// Moving-load family over keta's influence engine (influence.h,
// influence_surface.h).
//
// Physics / units: SI throughout. An ordinate is the response per unit force
// along D (|D| scales it); the canonical direction is gravity (0, 0, -1).
// Geometric tolerance: max(ctx.tolerance, 1 mm), so line -> bar matching agrees
// with the Assemble welder. The Moving Load abscissa is cumulative along the
// path in line order; a line authored against the travel direction is still
// emitted in its own orientation with its local position mirrored, so the
// convoy lands on the same physical station the influence line uses.
//
// Engine route: the sweep is by the adjoint (one solve K w = g, then
// eta_k = w.r_k per station, see keta::unit_load_sweep) instead of one
// back-substitution per station. The factorization persists across recomputes
// in a SharedFactorCache; a convoy sweep takes S >= 1 steps (S + 1 positions).
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "keta/influence_surface.h"
#include "keta/param/components.h"

namespace keta::param {

namespace {

const Tree& in_at(const std::vector<Tree>& in, std::size_t i) { return input_or_empty(in, i); }

bool finite3(const Vec3& p) { return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()); }
bool close(const Vec3& a, const Vec3& b, double tol) { return (a - b).norm() <= tol; }
double geo_tol(const SolveCtx& ctx) { return std::max(ctx.tolerance, 1e-3); }
std::string num_str(double x) { return Value::number(x).describe(); }

template <class T>
std::optional<T> first_spec(const Tree& t) {
  if (const T* x = first_of<T>(t)) return *x;
  return std::nullopt;
}

ModelRef model_of(const Tree& t) {
  if (const ModelRef* m = first_of<ModelRef>(t); m && *m) return *m;
  throw ComponentError::missing("M");
}

// ---- response quantities ----

enum class Cat { BarForce, Reaction, Displacement };
struct Quantity {
  const char* name;
  Cat cat;
  int comp;
};
constexpr Quantity kQuantities[18] = {
    {"N", Cat::BarForce, 0},      {"Vy", Cat::BarForce, 1},     {"Vz", Cat::BarForce, 2},
    {"Mt", Cat::BarForce, 3},     {"My", Cat::BarForce, 4},     {"Mz", Cat::BarForce, 5},
    {"Rx", Cat::Reaction, 0},     {"Ry", Cat::Reaction, 1},     {"Rz", Cat::Reaction, 2},
    {"Rmx", Cat::Reaction, 3},    {"Rmy", Cat::Reaction, 4},    {"Rmz", Cat::Reaction, 5},
    {"dx", Cat::Displacement, 0}, {"dy", Cat::Displacement, 1}, {"dz", Cat::Displacement, 2},
    {"rx", Cat::Displacement, 3}, {"ry", Cat::Displacement, 4}, {"rz", Cat::Displacement, 5}};

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

// Exact case first ("Rz" reaction vs "rz" rotation), then case-insensitive.
const Quantity& parse_quantity(const std::string& raw) {
  std::string t = raw;
  t.erase(0, t.find_first_not_of(" \t\r\n"));
  t.erase(t.find_last_not_of(" \t\r\n") + 1);
  for (const Quantity& q : kQuantities)
    if (t == q.name) return q;
  for (const Quantity& q : kQuantities)
    if (lower(t) == lower(q.name)) return q;
  throw ComponentError("unknown influence quantity \"" + raw +
                       "\" (expected one of: N, Vy, Vz, Mt, My, Mz, Rx, Ry, Rz, Rmx, Rmy, Rmz, dx, dy, dz, rx, ry, rz)");
}

IlDim dim_of(Cat c, int comp) {
  if (c == Cat::Displacement) return comp < 3 ? IlDim::DispTranslation : IlDim::DispRotation;
  return comp < 3 ? IlDim::Force : IlDim::Moment;
}

Id nearest_node(const Model& m, const Vec3& p) {
  Id best = 0;
  double bd = HUGE_VAL;
  for (const keta::Node& n : m.nodes) {
    const double d = (n.x - p).norm();
    if (d < bd) {
      bd = d;
      best = n.id;
    }
  }
  return best;
}

// The bar whose segment passes nearest p, with the clamped projection xi.
std::pair<Id, double> nearest_bar_station(const Model& m, const Vec3& p) {
  Id best = 0;
  double bd = HUGE_VAL, bt = 0;
  for (const Bar& b : m.bars) {
    const Vec3 a = m.node(b.ni).x, c = m.node(b.nj).x, ab = c - a;
    const double len2 = ab.squaredNorm();
    const double t = len2 <= 0 ? 0.0 : std::clamp((p - a).dot(ab) / len2, 0.0, 1.0);
    const double d = (p - (a + t * ab)).norm();
    if (d < bd) {
      bd = d;
      best = b.id;
      bt = t;
    }
  }
  return {best, bt};
}

struct Target {
  ResponseKind response;
  std::string label;
  IlDim dim;
};

Target resolve_target(const Model& m, const Tree& q_in, const Tree& t_in, const char* who) {
  const auto t = first_point(t_in);
  if (!t) throw ComponentError::missing("T");
  if (!finite3(*t)) throw ComponentError(std::string(who) + ": target point is not finite");
  const std::string q_text = first_text(q_in).value_or("My");
  const Quantity& q = parse_quantity(q_text);
  Target out;
  out.dim = dim_of(q.cat, q.comp);
  if (q.cat == Cat::BarForce) {
    const auto [bar, xi] = nearest_bar_station(m, *t);
    out.response = ResponseKind::bar_force(bar, xi * m.bar_chord_length(m.bar(bar)), q.comp);
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.2f", xi);
    out.label = std::string(q.name) + " @ ξ=" + buf + " bar " + std::to_string(bar);
  } else {
    const Id node = nearest_node(m, *t);
    out.response = q.cat == Cat::Reaction ? ResponseKind::reaction(node, q.comp)
                                          : ResponseKind::displacement(node, q.comp);
    out.label = std::string(q.name) + " @ node " + std::to_string(node);
  }
  return out;
}

Vec3 direction(const Tree& d_in, const char* who) {
  const Vec3 d = first_vector(d_in).value_or(Vec3(0, 0, -1));
  if (!finite3(d)) throw ComponentError(std::string(who) + ": load direction is not finite");
  return d;
}

// Travel flags of a contiguous chain of lines: flipped[i] = the path runs from
// the line's b to its a.
std::vector<bool> path_flips(const std::vector<Line>& lines, double tol) {
  std::vector<bool> flipped(lines.size(), false);
  if (lines.size() < 2) return flipped;
  const Line &l0 = lines[0], &l1 = lines[1];
  Vec3 shared;
  if (close(l0.a, l1.a, tol) || close(l0.a, l1.b, tol)) shared = l0.a;
  else if (close(l0.b, l1.a, tol) || close(l0.b, l1.b, tol)) shared = l0.b;
  else throw ComponentError("moving load: path lines 1 and 2 do not share an endpoint");
  flipped[0] = close(shared, l0.a, tol);
  Vec3 prev = shared;
  for (std::size_t i = 1; i < lines.size(); ++i) {
    if (close(lines[i].a, prev, tol)) {
      prev = lines[i].b;
    } else if (close(lines[i].b, prev, tol)) {
      flipped[i] = true;
      prev = lines[i].a;
    } else {
      throw ComponentError("moving load: path lines " + std::to_string(i) + " and " + std::to_string(i + 1) +
                           " do not share an endpoint");
    }
  }
  return flipped;
}

void check_convoy(const ConvoySpec& c, const char* who) {
  for (const auto& [o, f] : c.axles)
    if (!(std::isfinite(o) && std::isfinite(f))) throw ComponentError(std::string(who) + ": a convoy axle is not finite");
  for (const auto& s : c.segments)
    if (!(std::isfinite(s[0]) && std::isfinite(s[1]) && std::isfinite(s[2])))
      throw ComponentError(std::string(who) + ": a convoy segment is not finite");
}

}  // namespace

// ============================ Load Path ============================

ComponentInfo LoadPathNode::info() const {
  return {"keta.loadpath", "Load Path", "Analyze", "Ordered contiguous beams forming a moving-load travel line"};
}
std::vector<PinSpec> LoadPathNode::inputs() const {
  return {PinSpec::list("B", "beam", "Ordered beams forming the travel line"),
          PinSpec::item("N", "text", "Path name").opt()};
}
std::vector<PinSpec> LoadPathNode::outputs() const { return {PinSpec::item("P", "path", "The travel path")}; }
Outputs<Value> LoadPathNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  PathSpec p;
  for (const Value& v : in_at(in, 0).non_nulls())
    if (const auto* b = v.get<BeamSpec>()) p.lines.push_back(b->line);
  if (p.lines.empty()) return {Tree()};  // empty in, empty out (not an error)
  for (std::size_t i = 0; i < p.lines.size(); ++i)
    if (!finite3(p.lines[i].a) || !finite3(p.lines[i].b))
      throw ComponentError("load path: beam " + std::to_string(i + 1) + " has a non-finite endpoint");
  const double tol = geo_tol(ctx);
  for (std::size_t i = 0; i + 1 < p.lines.size(); ++i) {
    const Line &a = p.lines[i], &b = p.lines[i + 1];
    if (!(close(a.a, b.a, tol) || close(a.a, b.b, tol) || close(a.b, b.a, tol) || close(a.b, b.b, tol)))
      throw ComponentError("load path: beams " + std::to_string(i + 1) + " and " + std::to_string(i + 2) +
                           " do not share an endpoint (within tolerance " + num_str(tol) + ")");
  }
  p.name = first_text(in_at(in, 1)).value_or("");
  return {single(Value::of(std::move(p)))};
}

// ============================ Influence Line ============================

ComponentInfo InfluenceLineNode::info() const {
  return {"keta.influenceline", "Influence Line", "Analyze", "Influence line of a response as a unit load travels a path"};
}
std::vector<PinSpec> InfluenceLineNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"), PinSpec::item("P", "path", "The travel path"),
          PinSpec::item("Q", "text", "Response quantity (N,Vy,Vz,Mt,My,Mz; Rx..Rmz; dx..rz)").opt(),
          PinSpec::item("T", "point", "Target point (bar section or node)"),
          PinSpec::item("D", "vector", "Unit-load direction (default (0,0,-1))").opt(),
          PinSpec::item("S", "int", "Samples per bar (>= 2, default 5)").opt()};
}
std::vector<PinSpec> InfluenceLineNode::outputs() const {
  return {PinSpec::item("IL", "influence", "The influence line"),
          PinSpec::list("X", "number", "Abscissae along the path (m)"),
          PinSpec::list("V", "number", "Ordinates (response per unit load)")};
}
Outputs<Value> InfluenceLineNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  const ModelRef model = model_of(in_at(in, 0));
  const auto path = first_spec<PathSpec>(in_at(in, 1));
  if (!path) throw ComponentError::missing("P");
  if (path->lines.empty()) throw ComponentError("influence line: the path has no lines");
  const double tol = geo_tol(ctx);
  std::vector<Id> bars;
  for (std::size_t i = 0; i < path->lines.size(); ++i) {
    auto hit = find_bar(*model, path->lines[i], tol);
    if (!hit)
      throw ComponentError("influence line: path line " + std::to_string(i + 1) +
                           " matches no bar in the model (tolerance " + num_str(tol) + ")");
    bars.push_back(model->bars[hit->index].id);
  }
  const Target tgt = resolve_target(*model, in_at(in, 2), in_at(in, 3), "influence line");
  const Vec3 dir = direction(in_at(in, 4), "influence line");
  if (dir.squaredNorm() == 0) throw ComponentError("influence line: load direction is the zero vector — no load to sweep");
  const std::int64_t s = first_int(in_at(in, 5)).value_or(5);
  if (s < 2)
    throw ComponentError("influence line: samples per bar S=" + std::to_string(s) + " is too small (need >= 2)");
  if (ctx.cancelled()) throw ComponentError::cancelled();

  auto out = std::make_shared<InfluenceModel>();
  out->model = model;
  out->label = tgt.label;
  out->dir = dir;
  out->dim = tgt.dim;
  LoadPath lp;
  try {
    lp = build_path(*model, bars, static_cast<int>(s));
  } catch (const std::exception& e) {
    throw ComponentError(std::string("influence line: load path: ") + e.what());
  }
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  {
    std::lock_guard<std::mutex> lk(cache->mu);
    try {
      out->line = influence_line(*model, lp, dir, tgt.response, cache->cache, &hooks);
    } catch (const Cancelled&) {
      throw ComponentError::cancelled();  // polls sit outside the factorization: the cache stays coherent
    } catch (const std::exception& e) {
      cache->cache.invalidate();
      throw ComponentError(std::string("influence line: ") + e.what());
    }
  }
  Tree xs, vs;
  xs.ensure_path(Path::single(0));
  vs.ensure_path(Path::single(0));
  for (std::size_t k = 0; k < out->line.s.size(); ++k) {
    xs.append(Path::single(0), Value::number(out->line.s[k]));
    vs.append(Path::single(0), Value::number(out->line.eta[k]));
  }
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(InfluenceRef(std::move(out)))), xs, vs};
}

// ============================ Convoy ============================

ComponentInfo ConvoyNode::info() const {
  return {"keta.convoy", "Convoy", "Analyze", "Axle + uniform-segment convoy load pattern"};
}
std::vector<PinSpec> ConvoyNode::inputs() const {
  return {PinSpec::list("O", "number", "Axle offsets (m)"), PinSpec::list("F", "number", "Axle loads (N)"),
          PinSpec::item("SS", "number", "Segment start offset (m, default 0)").opt(),
          PinSpec::item("SL", "number", "Segment length (m)").opt(),
          PinSpec::item("SW", "number", "Segment total load (N)").opt(),
          PinSpec::item("N", "text", "Convoy name").opt()};
}
std::vector<PinSpec> ConvoyNode::outputs() const { return {PinSpec::item("C", "convoy", "The convoy pattern")}; }
Outputs<Value> ConvoyNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::vector<double> off, load;
  for (const Value& v : in_at(in, 0).non_nulls())
    if (auto x = v.as_number()) off.push_back(*x);
  for (const Value& v : in_at(in, 1).non_nulls())
    if (auto x = v.as_number()) load.push_back(*x);
  if (off.size() != load.size())
    throw ComponentError("convoy: " + std::to_string(off.size()) + " axle offset(s) but " + std::to_string(load.size()) +
                         " axle load(s) (O and F must be equal length)");
  for (std::size_t i = 0; i < off.size(); ++i)
    if (!std::isfinite(off[i])) throw ComponentError("convoy: axle offset " + std::to_string(i + 1) + " is not finite");
  for (std::size_t i = 0; i < load.size(); ++i)
    if (!std::isfinite(load[i])) throw ComponentError("convoy: axle load " + std::to_string(i + 1) + " is not finite");
  ConvoySpec c;
  for (std::size_t i = 0; i < off.size(); ++i) c.axles.emplace_back(off[i], load[i]);
  const auto sl = first_number(in_at(in, 3)), sw = first_number(in_at(in, 4));
  const double ss = first_number(in_at(in, 2)).value_or(0.0);
  if (sl && sw) {
    if (!std::isfinite(ss) || !std::isfinite(*sl) || !std::isfinite(*sw))
      throw ComponentError("convoy: segment start / length / total load must all be finite");
    if (*sl > 0) c.segments.push_back({ss, *sl, *sw});  // a non-positive length is inert
  }
  if (c.axles.empty() && c.segments.empty())
    throw ComponentError("convoy is empty: wire axle offsets/loads (O, F) or a segment (SL and SW)");
  c.name = first_text(in_at(in, 5)).value_or("");
  return {single(Value::of(std::move(c)))};
}

// ============================ Moving Load ============================

ComponentInfo MovingLoadNode::info() const {
  return {"keta.movingload", "Moving Load", "Analyze", "Sweep a convoy along a path into per-position load cases"};
}
std::vector<PinSpec> MovingLoadNode::inputs() const {
  return {PinSpec::item("P", "path", "The travel path"), PinSpec::item("C", "convoy", "The convoy"),
          PinSpec::item("S", "int", "Sweep steps (positions = S+1, default 20)").opt(),
          PinSpec::item("B", "int", "Base load-case id (default 1000)").opt(),
          PinSpec::item("D", "vector", "Unit-load direction (default (0,0,-1))").opt(),
          PinSpec::list("G", "int", "Permanent case ids for generated combos").opt(),
          PinSpec::item("GF", "number", "Permanent factor (default 1.35)").opt(),
          PinSpec::item("QF", "number", "Variable factor (default 1.5)").opt()};
}
std::vector<PinSpec> MovingLoadNode::outputs() const {
  return {PinSpec::list("L", "load", "Generated loads (flat list; case id = B+k)"),
          PinSpec::list("CB", "combo", "One ULS combo per loaded position (empty without G; skipped where no load lands)"),
          PinSpec::list("X", "number", "Reference abscissae of the sweep (m)")};
}
Outputs<Value> MovingLoadNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  const auto path = first_spec<PathSpec>(in_at(in, 0));
  if (!path) throw ComponentError::missing("P");
  const auto convoy = first_spec<ConvoySpec>(in_at(in, 1));
  if (!convoy) throw ComponentError::missing("C");
  if (path->lines.empty()) throw ComponentError("moving load: the path has no lines");

  std::vector<double> len, cum;
  double ltot = 0;
  for (std::size_t i = 0; i < path->lines.size(); ++i) {
    const double l = path->lines[i].length();
    if (!(std::isfinite(l) && l > 0))
      throw ComponentError("moving load: path line " + std::to_string(i + 1) + " has zero or non-finite length");
    cum.push_back(ltot);
    len.push_back(l);
    ltot += l;
  }
  const std::vector<bool> flipped = path_flips(path->lines, geo_tol(ctx));
  check_convoy(*convoy, "moving load");
  const Vec3 dir = direction(in_at(in, 4), "moving load");
  const std::int64_t steps = first_int(in_at(in, 2)).value_or(20);
  if (steps < 1)
    throw ComponentError("moving load: sweep steps S=" + std::to_string(steps) + " is too small (need >= 1)");
  const std::int64_t base = first_int(in_at(in, 3)).value_or(1000);
  const double gf = first_number(in_at(in, 6)).value_or(1.35), qf = first_number(in_at(in, 7)).value_or(1.5);
  if (!(std::isfinite(gf) && std::isfinite(qf))) throw ComponentError("moving load: combo factors GF / QF must be finite");
  std::vector<Id> perm;
  for (const Value& v : in_at(in, 5).non_nulls())
    if (auto i = v.as_int()) perm.push_back(*i);

  // Reference sweep: every axle enters and leaves.
  double min_off = 0, max_off = 0;
  for (const auto& [o, f] : convoy->axles) {
    min_off = std::min(min_off, o);
    max_off = std::max(max_off, o);
  }
  for (const auto& s : convoy->segments)
    if (s[1] > 0) {
      min_off = std::min(min_off, s[0]);
      max_off = std::max(max_off, s[0] + s[1]);
    }
  const double start = -max_off, end = ltot - min_off, dr = (end - start) / static_cast<double>(steps);
  // Earliest line whose span contains p (a junction lands at the end of the earlier line).
  auto locate = [&](double p) {
    for (std::size_t i = 0; i < len.size(); ++i)
      if (p <= cum[i] + len[i] || i + 1 == len.size()) return std::make_pair(i, std::clamp((p - cum[i]) / len[i], 0.0, 1.0));
    return std::make_pair(len.size() - 1, 1.0);
  };

  const Path b0 = Path::single(0);
  Tree loads, combos, xs;
  for (Tree* t : {&loads, &combos, &xs}) t->ensure_path(b0);
  for (std::int64_t k = 0; k <= steps; ++k) {
    const double r = start + dr * static_cast<double>(k);
    xs.append(b0, Value::number(r));
    const Id case_id = base + k;
    bool emitted = false;
    for (const auto& [o, f] : convoy->axles) {
      const double p = r + o;
      if (p < 0 || p > ltot) continue;  // off-path axles are dropped
      const auto [li, t] = locate(p);
      LoadSpec l;
      l.case_id = case_id;
      l.kind = BarPointLoad{{path->lines[li]}, flipped[li] ? 1 - t : t, dir * f, LoadFrame::Global};
      loads.append(b0, Value::of(std::move(l)));
      emitted = true;
    }
    for (const auto& s : convoy->segments) {
      if (s[1] <= 0) continue;
      const double w = s[2] / s[1], lo = std::max(r + s[0], 0.0), hi = std::min(r + s[0] + s[1], ltot);
      if (hi <= lo) continue;
      for (std::size_t li = 0; li < len.size(); ++li) {
        const double a = std::max(lo, cum[li]), b = std::min(hi, cum[li] + len[li]);
        if (b <= a) continue;
        double x1 = (a - cum[li]) / len[li], x2 = (b - cum[li]) / len[li];
        if (flipped[li]) std::tie(x1, x2) = std::make_pair(1 - x2, 1 - x1);
        LoadSpec l;
        l.case_id = case_id;
        l.kind = LinearLoad{{path->lines[li]}, dir * w, dir * w, LoadFrame::Global, x1, x2};
        loads.append(b0, Value::of(std::move(l)));
        emitted = true;
      }
    }
    // A position without load makes no combination (its case would not exist).
    if (!perm.empty() && emitted) {
      ComboSpec c;
      char buf[96];
      std::snprintf(buf, sizeof buf, "pos %lld @ %.3f m", static_cast<long long>(k), r);
      c.name = buf;
      for (Id g : perm) c.terms.emplace_back(gf, g);
      c.terms.emplace_back(qf, case_id);
      combos.append(b0, Value::of(std::move(c)));
    }
  }
  return {loads, combos, xs};
}

// ============================ IL Extremes ============================

ComponentInfo IlExtremesNode::info() const {
  return {"keta.ilextremes", "IL Extremes", "Analyze", "Worst-case response of an influence line under a convoy / UDL"};
}
std::vector<PinSpec> IlExtremesNode::inputs() const {
  return {PinSpec::item("IL", "influence", "The influence line"),
          PinSpec::item("C", "convoy", "Convoy load pattern").opt(),
          PinSpec::item("Q", "number", "UDL intensity (N/m)").opt(),
          PinSpec::item("S", "int", "Convoy sweep steps (default 200)").opt()};
}
std::vector<PinSpec> IlExtremesNode::outputs() const {
  return {PinSpec::item("MAX", "number", "Most positive response"),
          PinSpec::item("MIN", "number", "Most negative response")};
}
Outputs<Value> IlExtremesNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const InfluenceRef* il = first_of<InfluenceRef>(in_at(in, 0));
  if (!il || !*il) throw ComponentError::missing("IL");
  const auto convoy = first_spec<ConvoySpec>(in_at(in, 1));
  const auto q = first_number(in_at(in, 2));
  if (!convoy && !q) throw ComponentError("IL Extremes: wire a convoy (C) or a UDL intensity (Q)");
  const int steps = static_cast<int>(std::clamp<std::int64_t>(first_int(in_at(in, 3)).value_or(200), 1, 1'000'000));
  double mx = 0, mn = 0;
  if (convoy) {
    check_convoy(*convoy, "IL Extremes");
    ConvoyPattern pat;
    pat.axles = convoy->axles;
    pat.segments = convoy->segments;
    const Extremes e = convoy_extreme((*il)->line, pat, steps);
    // An optional variable action: leaving it off is always possible.
    mx += std::max(e.max, 0.0);
    mn += std::min(e.min, 0.0);
  }
  if (q) {
    if (!std::isfinite(*q)) throw ComponentError("IL Extremes: UDL intensity Q is not finite");
    const Extremes e = udl_extreme((*il)->line, *q);
    mx += e.max;
    mn += e.min;
  }
  return {single(Value::number(mx)), single(Value::number(mn))};
}

// ============================ Influence Surface ============================

ComponentInfo InfluenceSurfaceNode::info() const {
  return {"keta.influence_surface", "Influence Surface", "Analyze", "Influence surface of a response over a deck node grid"};
}
std::vector<PinSpec> InfluenceSurfaceNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::tree("G", "point", "Deck grid points (one branch per row)"),
          PinSpec::item("Q", "text", "Response quantity (N,Vy,Vz,Mt,My,Mz; Rx..Rmz; dx..rz)").opt(),
          PinSpec::item("T", "point", "Target point (bar section or node)"),
          PinSpec::item("D", "vector", "Unit-load direction (default (0,0,-1))").opt()};
}
std::vector<PinSpec> InfluenceSurfaceNode::outputs() const {
  return {PinSpec::tree("V", "number", "Ordinates (one branch per grid row, one item per column)"),
          PinSpec::tree("P", "point", "Resolved grid node positions (same shape)")};
}
Outputs<Value> InfluenceSurfaceNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  const ModelRef model = model_of(in_at(in, 0));
  // One branch per grid row; each point snapped to the nearest model node.
  std::vector<std::vector<Id>> rows;
  for (const auto& kv : in_at(in, 1).branches()) {
    std::vector<Id> row;
    for (const auto& o : kv.second)
      if (o)
        if (auto p = o->as_point()) row.push_back(nearest_node(*model, *p));
    if (!row.empty()) rows.push_back(std::move(row));
  }
  if (rows.empty()) throw ComponentError::missing("G");
  const Target tgt = resolve_target(*model, in_at(in, 2), in_at(in, 3), "influence surface");
  const Vec3 dir = direction(in_at(in, 4), "influence surface");
  for (std::size_t r = 1; r < rows.size(); ++r)
    if (rows[r].size() != rows[0].size())
      throw ComponentError("influence surface: grid rows must all have the same length (row 1 has " +
                           std::to_string(rows[0].size()) + ", row " + std::to_string(r + 1) + " has " +
                           std::to_string(rows[r].size()) + ")");
  if (ctx.cancelled()) throw ComponentError::cancelled();
  std::vector<Id> flat;
  for (const auto& row : rows) flat.insert(flat.end(), row.begin(), row.end());
  InfluenceSurface surf;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  try {
    const NodeGrid g = build_node_grid(*model, flat, static_cast<int>(rows[0].size()), static_cast<int>(rows.size()));
    std::lock_guard<std::mutex> lk(cache->mu);
    try {
      surf = influence_surface(*model, g, dir, tgt.response, cache->cache, &hooks);
    } catch (const Cancelled&) {
      throw;  // polls sit outside the factorization: the cache stays coherent
    } catch (...) {
      cache->cache.invalidate();
      throw;
    }
  } catch (const Cancelled&) {
    throw ComponentError::cancelled();
  } catch (const std::exception& e) {
    throw ComponentError(std::string("influence surface: ") + e.what());
  }
  Tree vs, ps;
  for (std::size_t r = 0; r < rows.size(); ++r) {
    const Path p = Path::single(static_cast<int>(r));
    for (std::size_t c = 0; c < rows[r].size(); ++c) {
      vs.append(p, Value::number(surf.eta[r * rows[r].size() + c]));
      ps.append(p, Value::of(Point{model->node(rows[r][c]).x}));
    }
  }
  return {vs, ps};
}

}  // namespace keta::param

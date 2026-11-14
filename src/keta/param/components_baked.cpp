// Data (Internalized), the bake operation and the Data Dam.
//
// Baking freezes a node's SOLVED output tree into an independent source node
// that carries the data in its params (Grasshopper's "Internalise data"); the
// graph that produced it can then be edited or deleted. The params are a
// tagged-JSON tree whose item tags are append-only; an item this build does
// not understand is preserved verbatim so a load -> save keeps it, and only a
// solve of it fails (typed error). Non-finite floats are refused at bake time
// (JSON has no NaN / inf) and degrade to a preserved unknown item on read.
//
// The Data Dam never inspects what an item is: it compares and copies whole
// trees, so one dam works on numbers, beams or results alike.
#include <cmath>

#include "keta/param/components.h"

namespace keta::param {

// ---------------------------------------------------------------------------
// Pin-tag universe
// ---------------------------------------------------------------------------

const std::vector<std::string>& all_type_tags() {
  static const std::vector<std::string> tags = {
      "number", "int",     "bool",     "text",       "point",    "vector",   "line",      "material", "section",
      "beam",   "cable",   "support",  "pointmass",  "load",     "combo",    "model",     "results",  "buckling",
      "modal",  "nonlinear", "arctrace", "cablestatic", "harmonic", "psd",    "spectrum",  "envelope", "set",
      "path",   "convoy",  "influence", "tendon",    "phase",    "phased",   "foundation", "check",   "slab"};
  return tags;
}

std::string known_type_tag(const std::string& tag) {
  const auto& all = all_type_tags();
  return std::find(all.begin(), all.end(), tag) != all.end() ? tag : "any";
}

// ---------------------------------------------------------------------------
// Data (Internalized)
// ---------------------------------------------------------------------------

namespace {

using Item = BakedDataNode::Item;
using BranchRec = BakedDataNode::BranchRec;

bool finite3(const Vec3& v) { return std::isfinite(v.x()) && std::isfinite(v.y()) && std::isfinite(v.z()); }

// The Value-universe tag of a real value item, empty for a null / unknown.
std::string item_tag(const Item& it) { return it.value ? it.value->type_tag() : std::string(); }

std::string dominant_tag(const std::vector<BranchRec>& branches) {
  std::string seen;
  for (const BranchRec& b : branches)
    for (const Item& it : b.items) {
      const std::string t = item_tag(it);
      if (t.empty()) continue;
      if (seen.empty())
        seen = t;
      else if (seen != t)
        return "any";  // mixed
    }
  return seen.empty() ? "any" : seen;
}

// Keys in alphabetical order, the document's key layout.
json::Value item_to_json(const Item& it) {
  if (it.raw) return *it.raw;
  if (!it.value) return json::Value();
  const Value& v = *it.value;
  json::Value o = json::Value::object();
  auto xyz = [&](const char* t, const Vec3& p) {
    o.set("t", t);
    o.set("x", p.x());
    o.set("y", p.y());
    o.set("z", p.z());
  };
  if (const auto* x = v.get<double>()) {
    o.set("t", "num");
    o.set("v", *x);
  } else if (const auto* i = v.get<std::int64_t>()) {
    o.set("t", "int");
    o.set("v", *i);
  } else if (const auto* b = v.get<bool>()) {
    o.set("t", "bool");
    o.set("v", *b);
  } else if (const auto* s = v.get<std::string>()) {
    o.set("t", "text");
    o.set("v", *s);
  } else if (const auto* p = v.get<Point>()) {
    xyz("point", p->p);
  } else if (const auto* q = v.get<Vector>()) {
    xyz("vector", q->v);
  } else if (const auto* l = v.get<Line>()) {
    o.set("ax", l->a.x());
    o.set("ay", l->a.y());
    o.set("az", l->a.z());
    o.set("bx", l->b.x());
    o.set("by", l->b.y());
    o.set("bz", l->b.z());
    o.set("t", "line");
  }
  return o;
}

// A JSON number that is an exact integer (keta's JSON numbers are doubles).
std::optional<std::int64_t> json_int(const json::Value* x) {
  if (!x || !x->is_number()) return std::nullopt;
  const double d = x->as_number();
  if (!std::isfinite(d) || std::trunc(d) != d || std::fabs(d) >= 9.2e18) return std::nullopt;
  return static_cast<std::int64_t>(d);
}

// Defensive parse: null -> null slot; a known well-formed tag -> its value;
// anything else (unknown tag, missing / wrong-kind / non-finite field) -> an
// unknown item kept verbatim.
Item item_from_json(const json::Value& v) {
  Item it;
  if (v.is_null()) return it;
  const Item unknown{std::nullopt, v};
  if (!v.is_object()) return unknown;
  auto f = [&](const char* k) -> std::optional<double> {
    const json::Value* x = v.get(k);
    if (!x || !x->is_number() || !std::isfinite(x->as_number())) return std::nullopt;
    return x->as_number();
  };
  const json::Value* t = v.get("t");
  const std::string tag = t && t->is_string() ? t->as_string() : std::string();
  if (tag == "num") {
    if (auto x = f("v")) return {Value::number(*x), std::nullopt};
  } else if (tag == "int") {
    if (auto i = json_int(v.get("v"))) return {Value::integer(*i), std::nullopt};
  } else if (tag == "bool") {
    if (const auto* b = v.get("v"); b && b->is_bool()) return {Value::boolean(b->as_bool()), std::nullopt};
  } else if (tag == "text") {
    if (const auto* s = v.get("v"); s && s->is_string()) return {Value::text(s->as_string()), std::nullopt};
  } else if (tag == "point" || tag == "vector") {
    auto x = f("x"), y = f("y"), z = f("z");
    if (x && y && z) {
      const Vec3 p(*x, *y, *z);
      return {tag == "point" ? Value::of(Point{p}) : Value::of(Vector{p}), std::nullopt};
    }
  } else if (tag == "line") {
    auto ax = f("ax"), ay = f("ay"), az = f("az"), bx = f("bx"), by = f("by"), bz = f("bz");
    if (ax && ay && az && bx && by && bz)
      return {Value::of(Line{Vec3(*ax, *ay, *az), Vec3(*bx, *by, *bz)}), std::nullopt};
  }
  return unknown;
}

}  // namespace

ComponentInfo BakedDataNode::info() const {
  return {"keta.baked", "Data (Internalized)", "Params",
          "Data frozen (internalized) into this node — independent of any upstream graph"};
}
std::vector<PinSpec> BakedDataNode::inputs() const { return {}; }
std::vector<PinSpec> BakedDataNode::outputs() const { return {PinSpec::tree("D", dominant, "The internalized data")}; }

Outputs<Value> BakedDataNode::solve(const SolveCtx&, const std::vector<Tree>&) const {
  Tree tree;
  for (const BranchRec& b : branches) {
    std::vector<int> idx;
    for (std::int64_t i : b.path) idx.push_back(static_cast<int>(i));  // narrows to the Path index type
    const Path path(idx);
    tree.ensure_path(path);  // an empty branch stays present
    for (const Item& it : b.items) {
      if (it.raw) {
        const json::Value* t = it.raw->is_object() ? it.raw->get("t") : nullptr;
        throw ComponentError("unknown baked item tag '" + (t && t->is_string() ? t->as_string() : "<none>") + "'");
      }
      if (it.value)
        tree.append(path, *it.value);
      else
        tree.append_null(path);
    }
  }
  return {tree};
}

json::Value BakedDataNode::params() const {
  json::Value arr = json::Value::array();
  for (const BranchRec& b : branches) {
    json::Value items = json::Value::array();
    for (const Item& it : b.items) items.push_back(item_to_json(it));
    json::Value path = json::Value::array();
    for (std::int64_t i : b.path) path.push_back(i);
    json::Value o = json::Value::object();
    o.set("items", std::move(items));
    o.set("path", std::move(path));
    arr.push_back(std::move(o));
  }
  json::Value p = json::Value::object();
  p.set("branches", std::move(arr));
  return p;
}

void BakedDataNode::set_params(const json::Value& p) {
  std::vector<BranchRec> out;
  if (const auto* arr = p.get("branches"); arr && arr->is_array())
    for (const json::Value& b : arr->as_array()) {
      BranchRec r;
      if (const auto* path = b.is_object() ? b.get("path") : nullptr; path && path->is_array())
        for (const json::Value& x : path->as_array())
          if (auto i = json_int(&x)) r.path.push_back(*i);
      if (const auto* items = b.is_object() ? b.get("items") : nullptr; items && items->is_array())
        for (const json::Value& x : items->as_array()) r.items.push_back(item_from_json(x));
      out.push_back(std::move(r));
    }
  dominant = dominant_tag(out);
  branches = std::move(out);
}

BakedDataNode BakedDataNode::from_tree(const Tree& t) {
  std::vector<std::string> kinds;  // offending, deduplicated, first-seen order
  bool non_finite = false;
  BakedDataNode n;
  for (const auto& [path, branch] : t.branches()) {
    BranchRec r;
    for (int i : path.indices()) r.path.push_back(i);
    for (const auto& slot : branch) {
      Item it;
      if (slot) {
        const Value& v = *slot;
        if (const auto* x = v.get<double>()) {
          non_finite |= !std::isfinite(*x);
        } else if (const auto* p = v.get<Point>()) {
          non_finite |= !finite3(p->p);
        } else if (const auto* q = v.get<Vector>()) {
          non_finite |= !finite3(q->v);
        } else if (const auto* l = v.get<Line>()) {
          non_finite |= !finite3(l->a) || !finite3(l->b);
        } else if (!v.get<std::int64_t>() && !v.get<bool>() && !v.get<std::string>()) {
          const std::string tag = v.type_tag();
          if (std::find(kinds.begin(), kinds.end(), tag) == kinds.end()) kinds.push_back(tag);
          r.items.push_back(it);
          continue;
        }
        it.value = v;
      }
      r.items.push_back(std::move(it));
    }
    n.branches.push_back(std::move(r));
  }
  if (!kinds.empty()) {
    std::string list;
    for (const std::string& k : kinds) list += (list.empty() ? "" : ", ") + k;
    throw BakeError("cannot bake: these value kinds are not bakeable: " + list);
  }
  if (non_finite) throw BakeError("cannot bake: a value is not finite (NaN or infinity)");
  n.dominant = dominant_tag(n.branches);
  return n;
}

NodeId bake_node(Graph<Value>& g, NodeId id) {
  const Node<Value>* node = g.node_ptr(id);
  if (!node) throw BakeError("cannot bake: the node no longer exists");
  const std::size_t n_out = node->component->outputs().size();
  if (n_out != 1)
    throw BakeError("cannot bake: only single-output nodes can be baked (this node has " + std::to_string(n_out) +
                    " outputs)");
  const auto pos = node->pos;
  const auto* outs = g.outputs(id);
  if (!outs || outs->empty()) throw BakeError("cannot bake: the node has no solved output yet");
  BakedDataNode baked = BakedDataNode::from_tree(outs->front());  // validates; graph untouched on error
  std::vector<std::pair<NodeId, std::size_t>> consumers;
  for (const Wire& w : g.wires())
    if (w.from.first == id && w.from.second == 0) consumers.push_back(w.to);
  const NodeId fresh = g.add_node(std::make_unique<BakedDataNode>(std::move(baked)), pos);
  g.remove_node(id);
  for (const auto& to : consumers) g.connect({fresh, 0}, to);  // fresh source: no cycle, no duplicate
  return fresh;
}

// ---------------------------------------------------------------------------
// Data Dam
// ---------------------------------------------------------------------------

std::optional<std::chrono::steady_clock::time_point> DataDamNode::stale_since() const {
  std::lock_guard<std::mutex> lock(state->m);
  return state->stale_since;
}

ComponentInfo DataDamNode::info() const {
  return {"keta.data_dam", "Data Dam", "Tree", "Hold upstream data until released manually or after a delay"};
}
std::vector<PinSpec> DataDamNode::inputs() const { return {PinSpec::tree("D", "any", "Data to dam").opt()}; }
std::vector<PinSpec> DataDamNode::outputs() const { return {PinSpec::tree("D", "any", "Held (or released) data")}; }

Outputs<Value> DataDamNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::vector<Tree> fresh{input_or_empty(in, 0)};
  State& s = *state;
  std::lock_guard<std::mutex> lock(s.m);
  auto clear_stale = [&] {
    s.stale = false;
    s.stale_since.reset();
  };
  // First solve ever (or after a document load): capture, hold and pass
  // through once. Passing through IS a release, so a release queued before
  // this first solve is consumed here (a leftover flag would turn the NEXT
  // upstream change into a release instead of a hold).
  if (!s.held) {
    s.release_pending = false;
    s.held = fresh;
    clear_stale();
    return fresh;
  }
  if (s.release_pending.exchange(false)) {
    s.held = fresh;
    clear_stale();
    return fresh;
  }
  // Holding: return the dammed data unchanged and record whether fresh input
  // is waiting (stamping the moment it first diverged).
  if (*s.held == fresh)
    clear_stale();
  else if (!s.stale.exchange(true))
    s.stale_since = std::chrono::steady_clock::now();
  return *s.held;
}

json::Value DataDamNode::params() const {
  json::Value p = json::Value::object();
  if (after_secs) {
    p.set("mode", "after_secs");
    p.set("secs", static_cast<double>(*after_secs));
  } else {
    p.set("mode", "manual");
  }
  return p;
}

void DataDamNode::set_params(const json::Value& p) {
  const auto* mode = p.get("mode");
  if (mode && mode->is_string() && mode->as_string() == "after_secs") {
    const auto* secs = p.get("secs");
    after_secs = static_cast<float>(secs && secs->is_number() ? secs->as_number() : 2.0);
  } else {
    after_secs.reset();
  }
}

}  // namespace keta::param

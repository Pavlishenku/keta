// Clusters. Hydration and the fold operation live in the document layer
// (document.cpp), which holds the registry.
#include "keta/param/components.h"

namespace keta::param {

namespace {

json::Value name_tag(const std::string& name, const std::string& tag) {
  json::Value p = json::Value::object();
  p.set("name", name);
  p.set("tag", tag);
  return p;
}

void read_name_tag(const json::Value& p, std::string& name, std::string& tag) {
  if (const auto* x = p.get("name"); x && x->is_string()) name = x->as_string();
  if (const auto* x = p.get("tag"); x && x->is_string()) tag = x->as_string();
}

}  // namespace

// ---------------------------------------------------------------------------
// Markers
// ---------------------------------------------------------------------------

ComponentInfo ClusterInNode::info() const {
  return {kClusterInTypeId, "Cluster Input", "Params", "Cluster plumbing: an input boundary of a cluster's inner graph"};
}
std::vector<PinSpec> ClusterInNode::inputs() const { return {}; }
std::vector<PinSpec> ClusterInNode::outputs() const {
  return {PinSpec::tree("V", known_type_tag(tag), "The injected value")};
}
Outputs<Value> ClusterInNode::solve(const SolveCtx&, const std::vector<Tree>&) const {
  return {injected.value_or(Tree{})};
}
json::Value ClusterInNode::params() const { return name_tag(name, tag); }
void ClusterInNode::set_params(const json::Value& p) { read_name_tag(p, name, tag); }

ComponentInfo ClusterOutNode::info() const {
  return {kClusterOutTypeId, "Cluster Output", "Params",
          "Cluster plumbing: an output boundary of a cluster's inner graph"};
}
std::vector<PinSpec> ClusterOutNode::inputs() const {
  return {PinSpec::tree("V", known_type_tag(tag), "The value to expose").opt()};
}
std::vector<PinSpec> ClusterOutNode::outputs() const {
  return {PinSpec::tree("V", known_type_tag(tag), "The exposed value")};
}
Outputs<Value> ClusterOutNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {input_or_empty(in, 0)};
}
json::Value ClusterOutNode::params() const { return name_tag(name, tag); }
void ClusterOutNode::set_params(const json::Value& p) { read_name_tag(p, name, tag); }

// ---------------------------------------------------------------------------
// Cluster
// ---------------------------------------------------------------------------

ComponentInfo ClusterNode::info() const {
  return {kClusterTypeId, "Cluster", "Params", "A subgraph packaged as one node (fold a selection, or load a preset)"};
}
std::vector<PinSpec> ClusterNode::inputs() const {
  std::vector<PinSpec> v;
  for (const ClusterPin& p : pins)
    if (!p.out) v.push_back(PinSpec::tree(p.name, known_type_tag(p.tag), "Cluster input"));
  return v;
}
std::vector<PinSpec> ClusterNode::outputs() const {
  std::vector<PinSpec> v;
  for (const ClusterPin& p : pins)
    if (p.out) v.push_back(PinSpec::tree(p.name, known_type_tag(p.tag), "Cluster output"));
  return v;
}

Outputs<Value> ClusterNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  if (!inner) {
    if (pins.empty() && doc.empty()) return {};  // a blank cluster is harmless
    throw ComponentError("cluster not hydrated");
  }
  // Work on a copy: the template stays clean and the engine memoizes THIS
  // node (one inner evaluation per outer-input change).
  Graph<Value> g = *inner;
  std::size_t k = 0;
  for (NodeId id : g.node_ids())
    if (g.node_ref(id).component->info().type_id == kClusterInTypeId) {
      Tree t = k < in.size() ? in[k] : Tree{};
      ++k;
      g.with_component(id, [&](Comp& c) {
        if (auto* ci = dynamic_cast<ClusterInNode*>(&c)) ci->injected = std::move(t);
      });
    }
  const EvalReport report = g.evaluate(ctx);
  if (report.cancelled) throw ComponentError::cancelled();
  for (NodeId id : g.node_ids())
    if (const ComponentError* e = g.error(id))
      throw ComponentError("cluster inner node '" + g.node_ref(id).component->info().name + "': " + e->what());
  std::vector<Tree> outs;
  for (NodeId id : g.node_ids())
    if (g.node_ref(id).component->info().type_id == kClusterOutTypeId) {
      const auto* o = g.outputs(id);
      outs.push_back(o && !o->empty() ? o->front() : Tree{});
    }
  return outs;
}

// Keys in alphabetical order (doc, name, pins; dir, name, tag), the document's
// key layout.
json::Value ClusterNode::params() const {
  json::Value arr = json::Value::array();
  for (const ClusterPin& p : pins) {
    json::Value o = json::Value::object();
    o.set("dir", p.out ? "out" : "in");
    o.set("name", p.name);
    o.set("tag", p.tag);
    arr.push_back(std::move(o));
  }
  json::Value v = json::Value::object();
  v.set("doc", doc);
  v.set("name", name);
  v.set("pins", std::move(arr));
  return v;
}

void ClusterNode::set_params(const json::Value& p) {
  if (const auto* x = p.get("name"); x && x->is_string()) name = x->as_string();
  if (const auto* x = p.get("doc"); x && x->is_string()) doc = x->as_string();
  if (const auto* arr = p.get("pins"); arr && arr->is_array()) {
    pins.clear();
    for (const json::Value& e : arr->as_array()) {
      ClusterPin cp;
      cp.tag = "any";
      if (e.is_object()) {
        if (const auto* x = e.get("name"); x && x->is_string()) cp.name = x->as_string();
        if (const auto* x = e.get("tag"); x && x->is_string()) cp.tag = x->as_string();
        if (const auto* x = e.get("dir"); x && x->is_string()) cp.out = x->as_string() == "out";
      }
      pins.push_back(std::move(cp));
    }
  }
  inner.reset();  // the document layer re-hydrates after the load
}

}  // namespace keta::param

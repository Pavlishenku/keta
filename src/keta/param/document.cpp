#include "keta/param/document.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <string_view>

#ifndef KETA_VERSION
#define KETA_VERSION "dev"
#endif

namespace keta::param {

namespace {

// ---------------- migration table ----------------

struct Step {
  int from_version;
  const char* label;
  void (*apply)(json::Value&);
};

// Schemas 1 and 2 name the component types `karamaru.<name>`; schema 3 names
// them `keta.<name>`. A cluster's inner document is a string that migrates on
// its own when the cluster is hydrated.
void rename_type_ids(json::Value& doc) {
  constexpr std::string_view kFormerPrefix = "karamaru.", kPrefix = "keta.";
  if (!doc.is_object() || !doc.has("nodes") || !doc.at("nodes").is_array()) return;
  json::Array nodes = doc.at("nodes").as_array();
  for (json::Value& node : nodes) {
    if (!node.is_object() || !node.has("type_id") || !node.at("type_id").is_string()) continue;
    const std::string id = node.at("type_id").as_string();
    if (id.starts_with(kFormerPrefix)) node.set("type_id", std::string(kPrefix) + id.substr(kFormerPrefix.size()));
  }
  doc.set("nodes", std::move(nodes));
}

// One row per schema hop, sorted by from_version.
const Step kSteps[] = {
    // v2 only appended defaulted blocks (meta, view): no structural rewrite.
    {1, "1->2: versioned document formalized (meta/view appended)", [](json::Value&) {}},
    {2, "2->3: component type ids renamed to `keta.<name>`", rename_type_ids},
};

int read_version(const json::Value& doc) {
  const json::Value* v = doc.get("version");
  if (!v || v->is_null()) return 1;
  if (!v->is_number())
    throw DocError(DocError::Kind::Parse, "document `version` must be a number, found " + v->dump());
  double n = v->as_number();
  if (!std::isfinite(n) || std::trunc(n) != n)
    throw DocError(DocError::Kind::Parse, "document `version` is not a valid schema number: " + v->dump());
  if (n < 1) throw DocError(DocError::Kind::Parse, "invalid document schema " + v->dump());
  if (n > 1e9) return static_cast<int>(1e9);
  return static_cast<int>(n);
}

DocError parse_error(const std::string& what) { return DocError(DocError::Kind::Parse, what); }

// ---------------- small json readers (typed failures, never a crash) ----------------

const json::Value& required(const json::Value& o, const char* key, const char* where) {
  const json::Value* v = o.get(key);
  if (!v) throw parse_error(std::string(where) + ": missing field `" + key + "`");
  return *v;
}

std::int64_t to_index(const json::Value& v, const char* what) {
  if (!v.is_number()) throw parse_error(std::string(what) + " must be a number");
  double n = v.as_number();
  if (!(n >= 0) || std::trunc(n) != n || n > 1e12) throw parse_error(std::string(what) + " must be a non-negative integer");
  return static_cast<std::int64_t>(n);
}

std::pair<std::int64_t, std::int64_t> index_pair(const json::Value& v, const char* what) {
  if (!v.is_array() || v.as_array().size() != 2) throw parse_error(std::string(what) + " must be a [key, pin] pair");
  return {to_index(v.as_array()[0], what), to_index(v.as_array()[1], what)};
}

std::pair<float, float> read_pos(const json::Value& v, const char* what) {
  if (!v.is_array() || v.as_array().size() != 2 || !v.as_array()[0].is_number() || !v.as_array()[1].is_number())
    throw parse_error(std::string(what) + " must be an [x, y] pair");
  return {static_cast<float>(v.as_array()[0].as_number()), static_cast<float>(v.as_array()[1].as_number())};
}

json::Value write_pos(std::pair<float, float> p) {
  // A non-finite position would poison the whole reload: pin it to the origin.
  float x = std::isfinite(p.first) ? p.first : 0.f, y = std::isfinite(p.second) ? p.second : 0.f;
  return json::Array{json::Value(static_cast<double>(x)), json::Value(static_cast<double>(y))};
}

bool read_bool(const json::Value& o, const char* key, bool dflt) {
  const json::Value* v = o.get(key);
  return v && v->is_bool() ? v->as_bool() : dflt;
}

PinFlags read_pin_flags(const json::Value& o) {
  return {read_bool(o, "flatten", false), read_bool(o, "graft", false), read_bool(o, "simplify", false),
          read_bool(o, "reverse", false)};
}

json::Value write_pin_flags(const PinFlags& f) {
  json::Value o = json::Value::object();
  o.set("flatten", f.flatten);
  o.set("graft", f.graft);
  o.set("simplify", f.simplify);
  o.set("reverse", f.reverse);
  return o;
}

const char* const kKnownKeys[] = {"version", "nodes", "wires", "units", "meta", "view", "groups", "annotations"};

// ---------------- graph rebuild ----------------

// Rebuild nodes + wires; returns key -> NodeId.
std::map<std::int64_t, NodeId> build_nodes_and_wires(const json::Value& doc, const Registry& reg, Graph<Value>& g,
                                                     std::pair<float, float> offset = {0.f, 0.f}) {
  const json::Value& nodes = required(doc, "nodes", "document");
  const json::Value& wires = required(doc, "wires", "document");
  if (!nodes.is_array() || !wires.is_array()) throw parse_error("document: `nodes` and `wires` must be arrays");
  std::map<std::int64_t, NodeId> keymap;
  for (const json::Value& dn : nodes.as_array()) {
    if (!dn.is_object()) throw parse_error("document: a node must be an object");
    std::int64_t key = to_index(required(dn, "key", "node"), "node key");
    const json::Value& tid = required(dn, "type_id", "node");
    if (!tid.is_string()) throw parse_error("node type_id must be a string");
    auto comp = reg.make(tid.as_string());
    if (!comp) throw DocError::unknown_type(tid.as_string());
    if (const json::Value* p = dn.get("params"); p && !p->is_null()) comp->set_params(*p);
    auto pos = read_pos(required(dn, "pos", "node"), "node pos");
    NodeId id = g.add_node(std::move(comp), {pos.first + offset.first, pos.second + offset.second});
    if (const json::Value* fl = dn.get("input_flags"); fl && fl->is_array())
      for (const json::Value& e : fl->as_array()) {
        if (!e.is_object()) continue;
        const json::Value* pin = e.get("pin");
        const json::Value* f = e.get("flags");
        if (pin && f && f->is_object()) g.set_input_flags(id, static_cast<std::size_t>(to_index(*pin, "pin")), read_pin_flags(*f));
      }
    if (const json::Value* nf = dn.get("flags"); nf && nf->is_object()) {
      g.set_node_preview(id, read_bool(*nf, "preview", true));
      g.set_node_enabled(id, read_bool(*nf, "enabled", true));
    }
    if (!keymap.emplace(key, id).second) throw parse_error("document: duplicate node key " + std::to_string(key));
  }
  for (const json::Value& w : wires.as_array()) {
    if (!w.is_object()) throw parse_error("document: a wire must be an object");
    auto from = index_pair(required(w, "from", "wire"), "wire from");
    auto to = index_pair(required(w, "to", "wire"), "wire to");
    auto fi = keymap.find(from.first), ti = keymap.find(to.first);
    if (fi == keymap.end()) throw DocError(DocError::Kind::Wire, "dangling wire source");
    if (ti == keymap.end()) throw DocError(DocError::Kind::Wire, "dangling wire target");
    try {
      g.connect({fi->second, static_cast<std::size_t>(from.second)}, {ti->second, static_cast<std::size_t>(to.second)});
    } catch (const GraphError& e) {
      throw DocError(DocError::Kind::Wire, "wire [" + std::to_string(from.first) + ", " + std::to_string(from.second) +
                                               "] -> [" + std::to_string(to.first) + ", " + std::to_string(to.second) +
                                               "]: " + e.what());
    }
  }
  return keymap;
}

json::Value parse_and_migrate(const std::string& text, MigrationReport* report) {
  json::Value root;
  try {
    root = json::Value::parse(text);
  } catch (const std::exception& e) {
    throw parse_error(e.what());
  }
  if (!root.is_object()) throw parse_error("document: the top level must be an object");
  MigrationReport r = migrate(root);
  if (report) *report = r;
  return root;
}

// ---------------- serialization ----------------

std::string serialize(const Graph<Value>& g, const std::vector<NodeId>& ids, const DocSettings* s) {
  std::map<NodeId, std::int64_t> index;
  for (std::size_t i = 0; i < ids.size(); ++i) index[ids[i]] = static_cast<std::int64_t>(i);

  json::Value root = json::Value::object();
  root.set("version", kCurrentSchema);
  json::Value nodes = json::Value::array();
  for (std::size_t i = 0; i < ids.size(); ++i) {
    const auto& n = g.node_ref(ids[i]);
    json::Value o = json::Value::object();
    o.set("key", static_cast<std::int64_t>(i));
    o.set("type_id", n.component->info().type_id);
    o.set("pos", write_pos(n.pos));
    json::Value p = n.component->params();
    if (!p.is_null()) o.set("params", std::move(p));
    json::Value flags = json::Value::array();
    const std::size_t n_in = n.component->inputs().size();
    for (std::size_t pin = 0; pin < n_in; ++pin) {
      PinFlags f = g.input_flags(ids[i], pin);
      if (f.is_default()) continue;
      json::Value e = json::Value::object();
      e.set("pin", static_cast<std::int64_t>(pin));
      e.set("flags", write_pin_flags(f));
      flags.push_back(std::move(e));
    }
    if (!flags.as_array().empty()) o.set("input_flags", std::move(flags));
    NodeFlags nf = g.node_flags(ids[i]);
    if (!nf.is_default()) {
      json::Value f = json::Value::object();
      f.set("preview", nf.preview);
      f.set("enabled", nf.enabled);
      o.set("flags", std::move(f));
    }
    nodes.push_back(std::move(o));
  }
  root.set("nodes", std::move(nodes));

  json::Value wires = json::Value::array();
  for (const Wire& w : g.wires()) {
    auto fi = index.find(w.from.first), ti = index.find(w.to.first);
    if (fi == index.end() || ti == index.end()) continue;  // a subset keeps internal wires only
    json::Value o = json::Value::object();
    o.set("from", json::Array{json::Value(fi->second), json::Value(static_cast<std::int64_t>(w.from.second))});
    o.set("to", json::Array{json::Value(ti->second), json::Value(static_cast<std::int64_t>(w.to.second))});
    wires.push_back(std::move(o));
  }
  root.set("wires", std::move(wires));
  root.set("units", (s ? s->units : DisplayUnits{}).to_json());

  if (s) {
    json::Value meta = json::Value::object();
    meta.set("app_version", std::string("keta ") + KETA_VERSION);
    root.set("meta", std::move(meta));
    if (s->view) root.set("view", *s->view);
    json::Value groups = json::Value::array();
    for (const CanvasGroup& gr : s->groups) {
      json::Value members = json::Value::array();
      for (NodeId id : gr.nodes)
        if (auto it = index.find(id); it != index.end()) members.push_back(json::Value(it->second));
      if (members.as_array().empty()) continue;  // stale membership is pruned
      json::Value o = json::Value::object();
      o.set("name", gr.name);
      o.set("color", json::Array{json::Value(gr.color[0]), json::Value(gr.color[1]), json::Value(gr.color[2])});
      o.set("nodes", std::move(members));
      groups.push_back(std::move(o));
    }
    if (!groups.as_array().empty()) root.set("groups", std::move(groups));
    json::Value notes = json::Value::array();
    for (const CanvasNote& n : s->notes) {
      json::Value o = json::Value::object();
      o.set("text", n.text);
      o.set("pos", write_pos(n.pos));
      notes.push_back(std::move(o));
    }
    if (!notes.as_array().empty()) root.set("annotations", std::move(notes));
    for (const auto& [k, v] : s->extra) root.set(k, v);
  }
  return root.dump(2);
}

}  // namespace

DocError DocError::newer_schema(int found, int supported) {
  return {Kind::NewerSchema, "this document was written by a newer keta (schema " + std::to_string(found) +
                                 "; this build reads up to " + std::to_string(supported) +
                                 ") — update keta to open it"};
}
DocError DocError::unknown_type(const std::string& type_id) {
  return {Kind::UnknownType, "unknown component type '" + type_id + "'"};
}


MigrationReport migrate(json::Value& doc) {
  MigrationReport r;
  r.from_schema = read_version(doc);
  if (r.from_schema > kCurrentSchema) throw DocError::newer_schema(r.from_schema, kCurrentSchema);
  for (const Step& s : kSteps)
    if (s.from_version >= r.from_schema) {
      s.apply(doc);
      r.steps.emplace_back(s.label);
    }
  if (doc.is_object()) doc.set("version", kCurrentSchema);
  return r;
}

LoadedDocument load_document(const std::string& text, const Registry& reg) {
  LoadedDocument out;
  json::Value doc = parse_and_migrate(text, &out.report);
  auto keymap = build_nodes_and_wires(doc, reg, out.graph);
  hydrate_clusters(out.graph, reg);

  DocSettings& s = out.settings;
  if (const json::Value* u = doc.get("units")) {
    // Missing families take the frozen defaults; an unknown spelling is a
    // typed parse error.
    try {
      s.units = DisplayUnits::from_json(*u);
    } catch (const keta::Error& e) {
      throw parse_error(e.what());
    }
  }
  if (const json::Value* v = doc.get("view"); v && !v->is_null()) s.view = *v;
  if (const json::Value* gs = doc.get("groups"); gs && gs->is_array())
    for (const json::Value& gr : gs->as_array()) {
      if (!gr.is_object()) continue;
      CanvasGroup cg;
      if (const json::Value* n = gr.get("name"); n && n->is_string()) cg.name = n->as_string();
      if (const json::Value* c = gr.get("color"); c && c->is_array() && c->as_array().size() == 3)
        for (int k = 0; k < 3; ++k)
          if (c->as_array()[static_cast<std::size_t>(k)].is_number())
            cg.color[static_cast<std::size_t>(k)] = static_cast<int>(c->as_array()[static_cast<std::size_t>(k)].as_number());
      if (const json::Value* ns = gr.get("nodes"); ns && ns->is_array())
        for (const json::Value& k : ns->as_array())
          if (k.is_number())
            if (auto it = keymap.find(static_cast<std::int64_t>(k.as_number())); it != keymap.end())
              cg.nodes.push_back(it->second);
      if (!cg.nodes.empty()) s.groups.push_back(std::move(cg));
    }
  if (const json::Value* ns = doc.get("annotations"); ns && ns->is_array())
    for (const json::Value& n : ns->as_array()) {
      if (!n.is_object()) continue;
      CanvasNote cn;
      if (const json::Value* t = n.get("text"); t && t->is_string()) cn.text = t->as_string();
      if (const json::Value* p = n.get("pos")) cn.pos = read_pos(*p, "annotation pos");
      s.notes.push_back(std::move(cn));
    }
  for (const auto& [k, v] : doc.members()) {
    bool known = false;
    for (const char* kk : kKnownKeys) known = known || k == kk;
    if (!known) s.extra.emplace_back(k, v);
  }
  return out;
}

std::string to_json(const Graph<Value>& g) { return serialize(g, g.node_ids(), nullptr); }

std::string to_json_with_settings(const Graph<Value>& g, const DocSettings& s) {
  return serialize(g, g.node_ids(), &s);
}

std::string subset_to_json(const Graph<Value>& g, const std::set<NodeId>& ids) {
  std::vector<NodeId> kept;
  for (NodeId id : g.node_ids())
    if (ids.count(id)) kept.push_back(id);
  return serialize(g, kept, nullptr);
}

std::vector<NodeId> paste_from_json(Graph<Value>& g, const Registry& reg, const std::string& text,
                                    std::pair<float, float> at) {
  json::Value doc = parse_and_migrate(text, nullptr);
  // Centre the pasted nodes on `at` (bounding-box centre of the stored positions).
  const json::Value& nodes = required(doc, "nodes", "document");
  if (!nodes.is_array()) throw parse_error("document: `nodes` must be an array");
  float min_x = 0, min_y = 0, max_x = 0, max_y = 0;
  bool first = true;
  for (const json::Value& n : nodes.as_array()) {
    if (!n.is_object()) continue;
    const json::Value* pv = n.get("pos");
    if (!pv) continue;
    auto p = read_pos(*pv, "node pos");
    if (first) {
      min_x = max_x = p.first;
      min_y = max_y = p.second;
      first = false;
    }
    min_x = std::min(min_x, p.first);
    max_x = std::max(max_x, p.first);
    min_y = std::min(min_y, p.second);
    max_y = std::max(max_y, p.second);
  }
  std::pair<float, float> offset{at.first - 0.5f * (min_x + max_x), at.second - 0.5f * (min_y + max_y)};
  // Build into a scratch graph first so a failing paste leaves `g` untouched,
  // then move the nodes over.
  Graph<Value> scratch;
  auto keymap = build_nodes_and_wires(doc, reg, scratch, offset);
  hydrate_clusters(scratch, reg);  // the clones below share the hydrated templates
  std::map<NodeId, NodeId> moved;
  std::vector<NodeId> out;
  for (NodeId sid : scratch.node_ids()) {
    auto& n = *scratch.node_ptr(sid);
    NodeId nid = g.add_node(n.component->clone(), n.pos);
    for (std::size_t pin = 0; pin < n.component->inputs().size(); ++pin)
      if (!scratch.input_flags(sid, pin).is_default()) g.set_input_flags(nid, pin, scratch.input_flags(sid, pin));
    g.set_node_preview(nid, n.flags.preview);
    g.set_node_enabled(nid, n.flags.enabled);
    moved[sid] = nid;
    out.push_back(nid);
  }
  for (const Wire& w : scratch.wires())
    g.connect({moved.at(w.from.first), w.from.second}, {moved.at(w.to.first), w.to.second});
  (void)keymap;
  return out;
}

// ---------------- clusters ----------------

namespace {

void hydrate_depth(Graph<Value>& g, const Registry& reg, int depth) {
  if (depth > kMaxClusterDepth)
    throw DocError(DocError::Kind::Cluster, "cluster nesting exceeds the maximum depth of " +
                                                std::to_string(kMaxClusterDepth) + " (a self-referencing cluster?)");
  for (NodeId id : g.node_ids()) {
    auto* c = dynamic_cast<ClusterNode*>(g.node_ptr(id)->component.get());
    if (!c) continue;
    json::Value inner_doc;
    try {
      inner_doc = parse_and_migrate(c->doc, nullptr);
    } catch (const DocError& e) {
      throw DocError(DocError::Kind::Cluster, std::string("cluster inner document failed to load: ") + e.what());
    }
    Graph<Value> inner;
    build_nodes_and_wires(inner_doc, reg, inner);
    hydrate_depth(inner, reg, depth + 1);
    auto shared = std::make_shared<const Graph<Value>>(std::move(inner));
    g.with_component(id, [&](Comp& comp) { dynamic_cast<ClusterNode&>(comp).inner = shared; });
  }
}

}  // namespace

void hydrate_clusters(Graph<Value>& g, const Registry& reg) { hydrate_depth(g, reg, 0); }

NodeId fold_selection(Graph<Value>& g, const std::set<NodeId>& selection, const std::string& name) {
  if (selection.empty()) throw FoldError("cannot fold: nothing is selected");
  for (NodeId id : selection)
    if (const auto* n = g.node_ptr(id); n && n->component->info().type_id == kClusterTypeId)
      throw FoldError("cannot fold: the selection contains a cluster (v1 does not nest via fold)");

  // Selection in deterministic (node) order.
  std::vector<NodeId> sel;
  std::map<NodeId, std::size_t> sel_index;
  for (NodeId id : g.node_ids())
    if (selection.count(id)) {
      sel_index[id] = sel.size();
      sel.push_back(id);
    }
  using End = std::pair<NodeId, std::size_t>;
  std::vector<End> in_ends, out_ends;
  for (const Wire& w : g.wires()) {
    const bool from_in = selection.count(w.from.first) > 0, to_in = selection.count(w.to.first) > 0;
    if (!from_in && to_in && std::find(in_ends.begin(), in_ends.end(), w.to) == in_ends.end())
      in_ends.push_back(w.to);
    else if (from_in && !to_in && std::find(out_ends.begin(), out_ends.end(), w.from) == out_ends.end())
      out_ends.push_back(w.from);
  }
  auto by_sel = [&](const End& a, const End& b) {
    const auto ka = std::make_pair(sel_index.count(a.first) ? sel_index[a.first] : SIZE_MAX, a.second);
    const auto kb = std::make_pair(sel_index.count(b.first) ? sel_index[b.first] : SIZE_MAX, b.second);
    return ka < kb;
  };
  std::stable_sort(in_ends.begin(), in_ends.end(), by_sel);
  std::stable_sort(out_ends.begin(), out_ends.end(), by_sel);

  // Inner graph: the selection (params via clone, node and pin flags — keep
  // aligned with `serialize`) and its internal wires.
  Graph<Value> inner;
  std::map<NodeId, NodeId> inner_map;
  for (NodeId old : sel) {
    const Node<Value>& n = g.node_ref(old);
    const NodeId nid = inner.add_node(n.component->clone(), n.pos);
    inner.set_node_preview(nid, n.flags.preview);
    inner.set_node_enabled(nid, n.flags.enabled);
    for (std::size_t pin = 0; pin < n.component->inputs().size(); ++pin)
      if (!g.input_flags(old, pin).is_default()) inner.set_input_flags(nid, pin, g.input_flags(old, pin));
    inner_map[old] = nid;
  }
  for (const Wire& w : g.wires())
    if (selection.count(w.from.first) && selection.count(w.to.first))
      inner.connect({inner_map[w.from.first], w.from.second}, {inner_map[w.to.first], w.to.second});

  // Boundary pins from the REAL PinSpecs, one marker each.
  std::vector<ClusterPin> pins;
  for (std::size_t i = 0; i < in_ends.size(); ++i) {
    const auto [target, pin] = in_ends[i];
    const auto spec = g.node_ref(target).component->inputs();
    ClusterPin cp{"in" + std::to_string(i), "any", false};
    if (pin < spec.size()) cp = {spec[pin].name, spec[pin].type_tag, false};
    auto marker = std::make_unique<ClusterInNode>();
    marker->name = cp.name;
    marker->tag = cp.tag;
    const NodeId m = inner.add_node(std::move(marker), {-160.f, static_cast<float>(i) * 60.f});
    inner.connect({m, 0}, {inner_map[target], pin});
    pins.push_back(cp);
  }
  for (std::size_t j = 0; j < out_ends.size(); ++j) {
    const auto [source, pin] = out_ends[j];
    const auto spec = g.node_ref(source).component->outputs();
    ClusterPin cp{"out" + std::to_string(j), "any", true};
    if (pin < spec.size()) cp = {spec[pin].name, spec[pin].type_tag, true};
    auto marker = std::make_unique<ClusterOutNode>();
    marker->name = cp.name;
    marker->tag = cp.tag;
    const NodeId m = inner.add_node(std::move(marker), {400.f, static_cast<float>(j) * 60.f});
    inner.connect({inner_map[source], pin}, {m, 0});
    pins.push_back(cp);
  }

  // External reconnections, collected before the outer graph changes.
  std::vector<std::pair<End, std::size_t>> in_reconnect;
  std::vector<std::pair<std::size_t, End>> out_reconnect;
  for (const Wire& w : g.wires()) {
    const bool from_in = selection.count(w.from.first) > 0, to_in = selection.count(w.to.first) > 0;
    if (!from_in && to_in) {
      const auto it = std::find(in_ends.begin(), in_ends.end(), w.to);
      if (it != in_ends.end()) in_reconnect.emplace_back(w.from, static_cast<std::size_t>(it - in_ends.begin()));
    } else if (from_in && !to_in) {
      const auto it = std::find(out_ends.begin(), out_ends.end(), w.from);
      if (it != out_ends.end()) out_reconnect.emplace_back(static_cast<std::size_t>(it - out_ends.begin()), w.to);
    }
  }
  float sx = 0, sy = 0;
  for (NodeId id : sel) {
    sx += g.node_ref(id).pos.first;
    sy += g.node_ref(id).pos.second;
  }
  const float cnt = static_cast<float>(sel.size());

  auto cluster = std::make_unique<ClusterNode>();
  cluster->name = name;
  cluster->doc = to_json(inner);
  cluster->pins = std::move(pins);
  cluster->inner = std::make_shared<const Graph<Value>>(std::move(inner));
  const NodeId cid = g.add_node(std::move(cluster), {sx / cnt, sy / cnt});
  // A refused reconnection is ignored (it would close a cycle through the
  // selection, which is removed next).
  for (const auto& [src, k] : in_reconnect) try {
      g.connect(src, {cid, k});
    } catch (const GraphError&) {
    }
  for (const auto& [k, dst] : out_reconnect) try {
      g.connect({cid, k}, dst);
    } catch (const GraphError&) {
    }
  for (NodeId id : sel) g.remove_node(id);
  return cid;
}

}  // namespace keta::param

// keta::param — the versioned parametric document: graph <-> JSON, the schema
// version gate and the JSON-level migration engine. UI-free: the headless load
// path needs only a Registry.
//
// ON-DISK SHAPE:
//
//   { "version": 3,
//     "nodes": [ { "key": 0, "type_id": "keta.slider", "pos": [x, y],
//                  "params": {...},                       // omitted when null
//                  "input_flags": [ {"pin": i, "flags": {flatten, graft, simplify, reverse}} ],
//                  "flags": { "preview": b, "enabled": b } } ],  // omitted when default
//     "wires": [ { "from": [key, pin], "to": [key, pin] } ],
//     "units": { "force": "Kn", ... },                    // display units (units.h)
//     "meta":  { "app_version": "..." },                  // informational, Save form only
//     "view":  { ... },                                   // viewport state (opaque here)
//     "groups": [ { "name", "color": [r,g,b], "nodes": [key...] } ],
//     "annotations": [ { "text", "pos": [x, y] } ],
//     ...unknown top-level keys preserved verbatim... }
//
// THE GRASSHOPPER STRATEGY. Migrations rewrite the RAW json BEFORE any typed
// reading, so old shapes never survive as C++ types: text -> json -> migrate()
// -> graph. An absent / null version means schema 1; a version above
// kCurrentSchema is a typed NewerSchema error; 0, negatives, non-integers and
// non-numbers are Parse errors. Adding a schema bump = raise kCurrentSchema and
// append one row to the step table (document.cpp).
//
// Wires bind (node key, pin INDEX): pins are append-only in every component.
// Node types resolve through Registry::resolve_type_id (one alias hop), so a
// renamed component keeps loading.
#pragma once
#include <array>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "keta/param/components.h"
#include "keta/param/units.h"

namespace keta::param {

inline constexpr int kCurrentSchema = 3;

struct MigrationReport {
  int from_schema = kCurrentSchema;  // as written (absent version => 1)
  std::vector<std::string> steps;    // labels of the applied steps (empty => current)
};

// Typed document failure. The what() wording is part of the contract.
class DocError : public std::runtime_error {
 public:
  enum class Kind { NewerSchema, Parse, UnknownType, Wire, Cluster };
  DocError(Kind k, std::string msg) : std::runtime_error(std::move(msg)), kind_(k) {}
  static DocError newer_schema(int found, int supported);
  static DocError unknown_type(const std::string& type_id);
  Kind kind() const { return kind_; }

 private:
  Kind kind_;
};

struct CanvasGroup {
  std::string name;
  std::array<int, 3> color{128, 128, 128};
  std::vector<NodeId> nodes;
};

struct CanvasNote {
  std::string text;
  std::pair<float, float> pos{0.f, 0.f};
};

// Document-level settings riding alongside the graph. `units` and `view` are
// kept as raw json (their typed forms belong to the display layer); `extra`
// holds the unknown top-level keys of a same-schema file written by a newer
// build, written back verbatim on save. `units` is typed (units.h); `view`
// stays raw json (the display layer's).
struct DocSettings {
  DisplayUnits units;  // the frozen kN / kN.m / mm / mrad / MPa / kN/m / kN.m/m set by default
  std::optional<json::Value> view;
  std::vector<CanvasGroup> groups;
  std::vector<CanvasNote> notes;
  json::Members extra;
};

struct LoadedDocument {
  Graph<Value> graph;
  DocSettings settings;
  MigrationReport report;
};

// Migrate a raw document up to kCurrentSchema in place (stamps "version").
MigrationReport migrate(json::Value& doc);

// text -> migrate -> graph (+ settings + report). Throws DocError.
LoadedDocument load_document(const std::string& text, const Registry& reg);

// Graph-only form (default units, no meta / view / groups): clipboard, undo.
std::string to_json(const Graph<Value>& g);
// Save form: units, meta.app_version, view, groups, notes and extra keys.
std::string to_json_with_settings(const Graph<Value>& g, const DocSettings& s);
// A subset of nodes and the wires internal to it (copy).
std::string subset_to_json(const Graph<Value>& g, const std::set<NodeId>& ids);
// Paste a (subset) document centred on `at`; returns the new node ids.
std::vector<NodeId> paste_from_json(Graph<Value>& g, const Registry& reg, const std::string& text,
                                    std::pair<float, float> at);

// Clusters.
//
// hydrate_clusters: parse every cluster node's stored inner document
// (through the migration gate), build its graph, hydrate ITS clusters
// recursively and install it. Run at the end of load_document and
// paste_from_json. Nesting deeper than kMaxClusterDepth (8 real levels; a
// self-referencing preset) and an unreadable inner document are typed
// DocError::Kind::Cluster errors, never a hang or a silent dead node.
inline constexpr int kMaxClusterDepth = 8;
void hydrate_clusters(Graph<Value>& g, const Registry& reg);

// Why a fold was refused; the wording of what() is contractual.
class FoldError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Fold the selected nodes into ONE hydrated cluster: each distinct incoming
// (target, pin) becomes a cluster input and each distinct outgoing (source,
// pin) an output, named and tagged from the real PinSpecs and ordered by
// (selection node order, pin); the inner graph is the selection (params,
// node and pin flags carried) plus one marker per boundary pin, stored as
// graph-only JSON and hydrated live. The selection is replaced by the
// cluster, placed at its centroid and wired to the old partners. Throws
// FoldError (graph untouched) on an empty selection or one containing a
// cluster.
NodeId fold_selection(Graph<Value>& g, const std::set<NodeId>& selection, const std::string& name);

}  // namespace keta::param

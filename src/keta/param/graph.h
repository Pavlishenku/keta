// keta::param — the dataflow node graph and its incremental evaluator.
//
// The graph is the artifact. Components declare typed pins once and implement
// a pure `solve(ctx, inputs) -> outputs` over DataTrees; per-node state (a
// slider's value) lives in the component instance each node owns and
// round-trips through `params()` / `set_params()` for persistence.
//
// EVALUATION (demand-driven, red/green with early cutoff). A node re-solves when
// (a) it has never solved, (b) it was expired, or (c) an upstream node's outputs
// actually CHANGED since this node last verified. After a re-solve the outputs
// are compared with the cache: if equal, `changed_at` is not bumped, so the
// downstream cone stays green — slider scrubbing on wide graphs stays cheap.
// `expire()` touches one node only; staleness propagates through the revision
// bookkeeping, never by recursive expiration.
//
// Order: topological (Kahn), ties broken by node insertion sequence; multi-wire
// inputs merge in wire insertion order — evaluation is deterministic.
//
// Per-pin flags (flatten -> graft -> simplify -> reverse) are applied by the
// engine to the merged input BEFORE the component sees it. Per-node flags:
// `preview` is display-only (bumps the revision, never expires); `enabled =
// false` skips the node AND its whole downstream cone (grey: no error, outputs
// cleared) — disabling is not a failure.
//
// Failures: a solve throwing ComponentError (or any std::exception, e.g. a typed
// keta::Error from an engine) marks the node failed (orange); downstream nodes
// are skipped with "upstream error". Cooperative cancellation (the token fires,
// or a component throws ComponentError::cancelled()) stops the walk, leaves the
// node's last outputs and error untouched and keeps it stale so a resumed pass
// re-solves it — cancellation never creates nor destroys error state.
//
// Node ids are slot-map keys {index, version}: stable across deletions, and a
// copied Graph keeps them identical (the async worker evaluates a COPY and merges
// results back by id: `adopt_results` + `touched_since`).
#pragma once
#include <algorithm>
#include <cstdint>
#include <memory>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "keta/json.h"
#include "keta/param/datatree.h"

namespace keta::param {

// ---------------------------------------------------------------------------
// Component contract
// ---------------------------------------------------------------------------

struct ComponentInfo {
  std::string type_id;      // stable across releases (documents bind to it)
  std::string name;         // canvas display name
  std::string category;     // palette category
  std::string description;  // one-line tooltip
};

struct PinSpec {
  std::string name;
  std::string type_tag;  // "number", "point", "model", "results", "any", ...
  std::string description;
  AccessMode access = AccessMode::Item;
  bool optional = false;

  static PinSpec item(std::string n, std::string tag, std::string d) {
    return {std::move(n), std::move(tag), std::move(d), AccessMode::Item, false};
  }
  static PinSpec list(std::string n, std::string tag, std::string d) {
    return {std::move(n), std::move(tag), std::move(d), AccessMode::List, false};
  }
  static PinSpec tree(std::string n, std::string tag, std::string d) {
    return {std::move(n), std::move(tag), std::move(d), AccessMode::Tree, false};
  }
  PinSpec opt() && {
    optional = true;
    return std::move(*this);
  }
};

// A solve failure. `Cancelled` is special-cased by the engine (not a failure).
class ComponentError : public std::runtime_error {
 public:
  enum class Kind { Message, MissingInput, Cancelled };
  ComponentError(Kind k, std::string msg) : std::runtime_error(std::move(msg)), kind_(k) {}
  explicit ComponentError(std::string msg) : ComponentError(Kind::Message, std::move(msg)) {}
  static ComponentError missing(const std::string& pin) {
    return {Kind::MissingInput, "input '" + pin + "' is required but empty"};
  }
  static ComponentError cancelled() { return {Kind::Cancelled, "solve cancelled"}; }
  Kind kind() const { return kind_; }
  friend bool operator==(const ComponentError& a, const ComponentError& b) {
    return a.kind_ == b.kind_ && std::string(a.what()) == b.what();
  }

 private:
  Kind kind_;
};

// Sink for long-phase progress reports from a heavy solve (frac in [0, 1]).
class ProgressSink {
 public:
  virtual ~ProgressSink() = default;
  virtual void report(const std::string& node_label, const std::string& phase, double frac) = 0;
};

struct SolveCtx {
  double tolerance = 1e-6;                  // geometric merge tolerance [m]
  std::optional<CancelToken> cancel;        // none = never cancelled
  std::shared_ptr<ProgressSink> progress;   // none = no reporting
  bool cancelled() const { return cancel && cancel->is_cancelled(); }
};

template <class V>
class Component {
 public:
  virtual ~Component() = default;
  virtual ComponentInfo info() const = 0;
  virtual std::vector<PinSpec> inputs() const = 0;
  virtual std::vector<PinSpec> outputs() const = 0;
  // One input tree per input pin (unwired optional pins arrive empty); returns
  // one tree per output pin. Throws ComponentError (or keta::Error) on failure.
  virtual Outputs<V> solve(const SolveCtx& ctx, const std::vector<DataTree<V>>& inputs) const = 0;
  // Persistent per-node state; null = stateless.
  virtual json::Value params() const { return json::Value(); }
  virtual void set_params(const json::Value&) {}
  // Every component is clonable (the whole graph copies per async dispatch).
  virtual std::unique_ptr<Component> clone() const = 0;
};

// CRTP helper giving `clone()` for free to any copyable component.
template <class Derived, class V>
class ComponentBase : public Component<V> {
 public:
  std::unique_ptr<Component<V>> clone() const override {
    return std::make_unique<Derived>(static_cast<const Derived&>(*this));
  }
};

// ---------------------------------------------------------------------------
// Graph
// ---------------------------------------------------------------------------

struct NodeId {
  std::uint32_t index = UINT32_MAX;
  std::uint32_t version = 0;
  bool valid() const { return index != UINT32_MAX; }
  friend bool operator==(const NodeId&, const NodeId&) = default;
  friend auto operator<=>(const NodeId&, const NodeId&) = default;
};

struct Wire {
  std::pair<NodeId, std::size_t> from;  // output pin
  std::pair<NodeId, std::size_t> to;    // input pin
  friend bool operator==(const Wire&, const Wire&) = default;
};

struct PinFlags {
  bool flatten = false, graft = false, simplify = false, reverse = false;
  bool is_default() const { return !flatten && !graft && !simplify && !reverse; }
  friend bool operator==(const PinFlags&, const PinFlags&) = default;
};

struct NodeFlags {
  bool preview = true;  // display only
  bool enabled = true;  // false = the node and its downstream cone are skipped
  bool is_default() const { return preview && enabled; }
  friend bool operator==(const NodeFlags&, const NodeFlags&) = default;
};

enum class GraphErrorKind { NodeNotFound, PinOutOfRange, Cycle, DuplicateWire };

class GraphError : public std::runtime_error {
 public:
  explicit GraphError(GraphErrorKind k)
      : std::runtime_error(k == GraphErrorKind::NodeNotFound    ? "node not found"
                           : k == GraphErrorKind::PinOutOfRange ? "pin index out of range"
                           : k == GraphErrorKind::Cycle         ? "connecting would create a cycle"
                                                                : "wire already exists"),
        kind_(k) {}
  GraphErrorKind kind() const { return kind_; }

 private:
  GraphErrorKind kind_;
};

enum class NodeStatus {
  Cached,                // cache valid, not recomputed
  Recomputed,            // recomputed, outputs changed
  Unchanged,             // recomputed, outputs identical (cutoff)
  Failed,                // solve failed (orange)
  SkippedUpstreamError,  // an upstream node failed
  Disabled,              // the node is disabled (grey, no error)
  SkippedDisabled,       // a disabled node lies upstream (grey, no error)
  Cancelled              // cooperative cancellation; node stays stale
};

struct EvalReport {
  std::vector<std::pair<NodeId, NodeStatus>> statuses;
  std::size_t recomputed = 0, cutoffs = 0, failures = 0;
  bool cancelled = false;
  NodeStatus status_of(NodeId id) const {
    for (const auto& [n, s] : statuses)
      if (n == id) return s;
    throw std::out_of_range("node not in report");
  }
};

template <class V>
struct Node {
  std::unique_ptr<Component<V>> component;
  std::pair<float, float> pos{0.f, 0.f};  // canvas position (ignored by the engine)
  std::uint64_t seq = 0;                  // insertion sequence (deterministic tie-break)
  std::vector<PinFlags> input_flags;      // sized lazily; missing = default
  NodeFlags flags;
  bool skipped_disabled = false;          // runtime grey cell, recomputed each pass
  std::optional<std::vector<DataTree<V>>> outputs;
  std::uint64_t changed_at = 0, verified_at = 0, expired_at = 0;
  std::optional<ComponentError> error;

  Node() = default;
  Node(const Node& o)
      : component(o.component ? o.component->clone() : nullptr), pos(o.pos), seq(o.seq),
        input_flags(o.input_flags), flags(o.flags), skipped_disabled(o.skipped_disabled),
        outputs(o.outputs), changed_at(o.changed_at), verified_at(o.verified_at),
        expired_at(o.expired_at), error(o.error) {}
  Node& operator=(const Node& o) {
    if (this != &o) {
      Node tmp(o);
      *this = std::move(tmp);
    }
    return *this;
  }
  Node(Node&&) noexcept = default;
  Node& operator=(Node&&) noexcept = default;
};

template <class V>
class Graph {
 public:
  std::uint64_t revision() const { return revision_; }

  NodeId add_node(std::unique_ptr<Component<V>> c, std::pair<float, float> pos = {0.f, 0.f}) {
    Node<V> n;
    n.component = std::move(c);
    n.pos = pos;
    n.seq = next_seq_++;
    ++revision_;
    std::uint32_t idx;
    if (!free_.empty()) {
      idx = free_.back();
      free_.pop_back();
    } else {
      idx = static_cast<std::uint32_t>(slots_.size());
      slots_.push_back({});
    }
    Slot& s = slots_[idx];
    s.node = std::move(n);
    return {idx, s.version};
  }

  // Remove a node and every wire touching it; its former consumers expire.
  void remove_node(NodeId id) {
    if (!contains(id)) return;
    Slot& s = slots_[id.index];
    s.node.reset();
    ++s.version;
    free_.push_back(id.index);
    ++revision_;
    std::vector<NodeId> downstream;
    for (const Wire& w : wires_)
      if (w.from.first == id) downstream.push_back(w.to.first);
    std::erase_if(wires_, [&](const Wire& w) { return w.from.first == id || w.to.first == id; });
    for (NodeId n : downstream) expire(n);
  }

  // Connect output pin -> input pin. Multiple wires into one input are legal.
  void connect(std::pair<NodeId, std::size_t> from, std::pair<NodeId, std::size_t> to) {
    if (!contains(from.first) || !contains(to.first)) throw GraphError(GraphErrorKind::NodeNotFound);
    if (from.second >= node_ref(from.first).component->outputs().size() ||
        to.second >= node_ref(to.first).component->inputs().size())
      throw GraphError(GraphErrorKind::PinOutOfRange);
    Wire w{from, to};
    if (std::find(wires_.begin(), wires_.end(), w) != wires_.end())
      throw GraphError(GraphErrorKind::DuplicateWire);
    if (from.first == to.first || reaches(to.first, from.first)) throw GraphError(GraphErrorKind::Cycle);
    wires_.push_back(w);
    ++revision_;
    expire(to.first);
  }

  void disconnect(const Wire& w) {
    auto before = wires_.size();
    std::erase(wires_, w);
    if (wires_.size() != before) {
      ++revision_;
      expire(w.to.first);
    }
  }

  // Mark one node stale (its next evaluation re-solves unconditionally).
  void expire(NodeId id) {
    ++revision_;
    if (Node<V>* n = node_ptr(id)) {
      n->verified_at = 0;
      n->expired_at = revision_;
    }
  }

  void bump_revision() { ++revision_; }

  // Nodes expired strictly after `dispatch_revision` (edited after a copy was taken).
  std::set<NodeId> touched_since(std::uint64_t dispatch_revision) const {
    std::set<NodeId> out;
    for (NodeId id : node_ids())
      if (node_ref(id).expired_at > dispatch_revision) out.insert(id);
    return out;
  }

  // Mutate a node's component (slider drag, text edit, ...) and expire it.
  template <class F>
  void with_component(NodeId id, F&& f) {
    Node<V>* n = node_ptr(id);
    if (!n) throw GraphError(GraphErrorKind::NodeNotFound);
    f(*n->component);
    expire(id);
  }

  // Merge the results of an asynchronously evaluated COPY back into this graph,
  // then re-expire the nodes edited since the copy was taken.
  void adopt_results(const Graph& evaluated, const std::set<NodeId>& touched) {
    for (NodeId id : node_ids()) {
      const Node<V>* src = evaluated.node_ptr(id);
      if (!src) continue;
      Node<V>& n = *node_ptr(id);
      n.outputs = src->outputs;
      n.changed_at = src->changed_at;
      n.verified_at = src->verified_at;
      n.error = src->error;
    }
    for (NodeId id : touched) expire(id);
    bump_revision();
  }

  PinFlags input_flags(NodeId id, std::size_t pin) const {
    const Node<V>* n = node_ptr(id);
    return n && pin < n->input_flags.size() ? n->input_flags[pin] : PinFlags{};
  }
  void set_input_flags(NodeId id, std::size_t pin, PinFlags f) {
    if (Node<V>* n = node_ptr(id)) {
      if (n->input_flags.size() <= pin) n->input_flags.resize(pin + 1);
      n->input_flags[pin] = f;
    }
    expire(id);
  }

  NodeFlags node_flags(NodeId id) const {
    const Node<V>* n = node_ptr(id);
    return n ? n->flags : NodeFlags{};
  }
  void set_node_preview(NodeId id, bool preview) {
    Node<V>* n = node_ptr(id);
    if (n && n->flags.preview != preview) {
      n->flags.preview = preview;
      ++revision_;
    }
  }
  void set_node_enabled(NodeId id, bool enabled) {
    Node<V>* n = node_ptr(id);
    if (n && n->flags.enabled != enabled) {
      n->flags.enabled = enabled;
      expire(id);
    }
  }
  bool is_skipped_disabled(NodeId id) const {
    const Node<V>* n = node_ptr(id);
    return n && n->skipped_disabled;
  }

  bool contains(NodeId id) const {
    return id.index < slots_.size() && slots_[id.index].version == id.version &&
           slots_[id.index].node.has_value();
  }
  const Node<V>* node_ptr(NodeId id) const { return contains(id) ? &*slots_[id.index].node : nullptr; }
  Node<V>* node_ptr(NodeId id) { return contains(id) ? &*slots_[id.index].node : nullptr; }
  const Node<V>& node_ref(NodeId id) const { return *node_ptr(id); }

  // Live node ids in insertion order.
  std::vector<NodeId> node_ids() const {
    std::vector<NodeId> ids;
    for (std::uint32_t i = 0; i < slots_.size(); ++i)
      if (slots_[i].node) ids.push_back({i, slots_[i].version});
    std::sort(ids.begin(), ids.end(),
              [&](NodeId a, NodeId b) { return node_ref(a).seq < node_ref(b).seq; });
    return ids;
  }
  std::size_t node_count() const { return slots_.size() - free_.size(); }
  const std::vector<Wire>& wires() const { return wires_; }

  // Cached outputs (nullptr until the first successful solve).
  const std::vector<DataTree<V>>* outputs(NodeId id) const {
    const Node<V>* n = node_ptr(id);
    return n && n->outputs ? &*n->outputs : nullptr;
  }
  const ComponentError* error(NodeId id) const {
    const Node<V>* n = node_ptr(id);
    return n && n->error ? &*n->error : nullptr;
  }

  // Evaluate the whole graph; only stale nodes recompute.
  EvalReport evaluate(const SolveCtx& ctx = {}) {
    EvalReport report;
    const std::uint64_t rev = revision_;
    for (NodeId id : topo_order()) {
      if (ctx.cancelled()) {
        report.cancelled = true;
        break;
      }
      NodeStatus st = eval_node(id, ctx, rev);
      report.statuses.emplace_back(id, st);
      if (st == NodeStatus::Recomputed) {
        ++report.recomputed;
      } else if (st == NodeStatus::Unchanged) {
        ++report.recomputed;
        ++report.cutoffs;
      } else if (st == NodeStatus::Failed) {
        ++report.failures;
      } else if (st == NodeStatus::Cancelled) {
        report.cancelled = true;
        break;
      }
    }
    return report;
  }

 private:
  struct Slot {
    std::optional<Node<V>> node;
    std::uint32_t version = 0;
  };

  bool reaches(NodeId start, NodeId target) const {
    std::vector<NodeId> stack{start};
    std::set<NodeId> seen;
    while (!stack.empty()) {
      NodeId n = stack.back();
      stack.pop_back();
      if (n == target) return true;
      if (!seen.insert(n).second) continue;
      for (const Wire& w : wires_)
        if (w.from.first == n) stack.push_back(w.to.first);
    }
    return false;
  }

  // Kahn, ready set ordered by insertion sequence.
  std::vector<NodeId> topo_order() const {
    std::vector<NodeId> ids = node_ids();
    std::vector<std::size_t> indeg(slots_.size(), 0);
    for (const Wire& w : wires_) ++indeg[w.to.first.index];
    auto by_seq = [&](NodeId a, NodeId b) { return node_ref(a).seq > node_ref(b).seq; };
    std::vector<NodeId> ready;
    for (NodeId id : ids)
      if (indeg[id.index] == 0) ready.push_back(id);
    std::sort(ready.begin(), ready.end(), by_seq);
    std::vector<NodeId> order;
    order.reserve(ids.size());
    while (!ready.empty()) {
      NodeId n = ready.back();
      ready.pop_back();
      order.push_back(n);
      for (const Wire& w : wires_)
        if (w.from.first == n && --indeg[w.to.first.index] == 0) ready.push_back(w.to.first);
      std::sort(ready.begin(), ready.end(), by_seq);
    }
    return order;
  }

  NodeStatus eval_node(NodeId id, const SolveCtx& ctx, std::uint64_t rev) {
    Node<V>& node = *node_ptr(id);
    // 1. Already verified this revision.
    if (node.verified_at == rev && node.outputs && !node.error) return NodeStatus::Cached;

    // 1a. Disabled: never solve; clear outputs and error (not a failure).
    if (!node.flags.enabled) {
      node.skipped_disabled = false;
      node.error.reset();
      if (node.outputs) {
        node.outputs.reset();
        node.changed_at = rev;
      }
      node.verified_at = rev;
      return NodeStatus::Disabled;
    }

    // 1b. Transitive grey state: any disabled / skipped source skips this node.
    for (const Wire& w : wires_) {
      if (w.to.first != id) continue;
      const Node<V>& src = node_ref(w.from.first);
      if (!src.flags.enabled || src.skipped_disabled) {
        node.skipped_disabled = true;
        node.error.reset();
        if (node.outputs) {
          node.outputs.reset();
          node.changed_at = rev;
        }
        node.verified_at = rev;
        return NodeStatus::SkippedDisabled;
      }
    }
    node.skipped_disabled = false;

    // 2. Gather inputs: merge source trees per pin in wire order, then apply the
    //    pin flags (flatten -> graft -> simplify -> reverse).
    const std::size_t n_in = node.component->inputs().size();
    std::vector<DataTree<V>> inputs;
    inputs.reserve(n_in);
    std::uint64_t deps_changed = 0;
    for (std::size_t pin = 0; pin < n_in; ++pin) {
      DataTree<V> merged;
      for (const Wire& w : wires_) {
        if (w.to.first != id || w.to.second != pin) continue;
        const Node<V>& src = node_ref(w.from.first);
        if (src.error) {
          node.error = ComponentError("upstream error");
          node.outputs.reset();
          node.verified_at = rev;
          return NodeStatus::SkippedUpstreamError;
        }
        deps_changed = std::max(deps_changed, src.changed_at);
        if (src.outputs && w.from.second < src.outputs->size()) {
          const DataTree<V>& t = (*src.outputs)[w.from.second];
          merged = merged.empty() ? t : merged.merge(t);
        }
      }
      PinFlags f = pin < node.input_flags.size() ? node.input_flags[pin] : PinFlags{};
      if (f.flatten) merged = merged.flatten();
      if (f.graft) merged = merged.graft();
      if (f.simplify) merged = merged.simplify();
      if (f.reverse) merged = merged.reverse_items();
      inputs.push_back(std::move(merged));
    }

    // 3. Up to date (verified after every dependency last changed, not expired).
    if (node.outputs && !node.error && node.verified_at != 0 && node.verified_at >= deps_changed) {
      node.verified_at = rev;
      return NodeStatus::Cached;
    }

    // 4. Re-solve.
    try {
      std::vector<DataTree<V>> outs = node.component->solve(ctx, inputs).release();
      node.verified_at = rev;
      node.error.reset();
      if (node.outputs && *node.outputs == outs) return NodeStatus::Unchanged;
      node.outputs = std::move(outs);
      node.changed_at = rev;
      return NodeStatus::Recomputed;
    } catch (const ComponentError& e) {
      if (e.kind() == ComponentError::Kind::Cancelled) {
        node.verified_at = 0;  // stays stale; last outputs / error untouched
        return NodeStatus::Cancelled;
      }
      fail_node(node, e, rev);
    } catch (const std::exception& e) {
      // A hooked engine driver aborted by the shared token (keta::Cancelled)
      // that a component did not map itself: still a cancellation, never a
      // failure frozen into the node.
      if (ctx.cancelled()) {
        node.verified_at = 0;
        return NodeStatus::Cancelled;
      }
      fail_node(node, ComponentError(e.what()), rev);
    }
    return NodeStatus::Failed;
  }

  static void fail_node(Node<V>& node, const ComponentError& e, std::uint64_t rev) {
    node.verified_at = rev;
    node.error = e;
    node.outputs.reset();
    node.changed_at = rev;
  }

  std::vector<Slot> slots_;
  std::vector<std::uint32_t> free_;
  std::vector<Wire> wires_;
  std::uint64_t revision_ = 0;
  std::uint64_t next_seq_ = 0;
};

}  // namespace keta::param

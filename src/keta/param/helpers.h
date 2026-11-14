// keta::param — data-matching helpers shared by the first-party components
// and the node welder.
#pragma once
#include <functional>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "keta/hooks.h"
#include "keta/param/graph.h"
#include "keta/param/value.h"

namespace keta::param {

using Tree = DataTree<Value>;
using Comp = Component<Value>;

// N-way item-wise map with Grasshopper longest-list semantics: positional
// branch alignment with last-branch repeat, then per-branch longest list with
// last-item repeat. Output paths come from the principal input (most branches;
// ties -> the earliest input). An entirely empty input is MissingInput (named
// by its pin); an item slot where ANY input is null is skipped.
Tree map_items(const std::vector<std::pair<std::string, const Tree*>>& inputs,
               const std::function<Value(const std::vector<const Value*>&)>& f);

// Engine hooks fed by a node's SolveCtx: the shared cancel token and the
// progress sink, reports labelled with the node's display name. The caller
// keeps the returned object alive across the hooked solve (options hold a
// non-owning pointer).
inline SolveHooks ctx_hooks(const SolveCtx& ctx, const std::string& node_label) {
  SolveHooks h;
  h.cancel = ctx.cancel;
  if (ctx.progress)
    h.progress = [sink = ctx.progress, node_label](double frac, const std::string& phase) {
      sink->report(node_label, phase, frac);
    };
  return h;
}

// A one-item tree {0}[v].
inline Tree single(Value v) { return Tree::item(std::move(v)); }

// The tree, or {0}[default] when it carries no value.
Tree number_or(const Tree& t, double dflt);

// First non-null value of a whole tree read through an accessor, else default.
std::optional<double> first_number(const Tree& t);
std::optional<std::int64_t> first_int(const Tree& t);
std::optional<Vec3> first_vector(const Tree& t);
std::optional<std::string> first_text(const Tree& t);
std::optional<bool> first_bool(const Tree& t);
std::optional<Vec3> first_point(const Tree& t);

// Positional branch access of a secondary input, built once per solve: at(bi)
// is the bi-th branch repeating the last one, nullptr for an empty tree — the
// per-branch alignment of the list family.
class BranchIndex {
 public:
  explicit BranchIndex(const Tree& t) {
    for (const auto& kv : t.branches()) b_.push_back(&kv.second);
  }
  const Branch<Value>* at(std::size_t bi) const { return b_.empty() ? nullptr : b_[std::min(bi, b_.size() - 1)]; }

 private:
  std::vector<const Branch<Value>*> b_;
};
// A boolean pattern read from a branch (nulls -> false), never empty: falls
// back to the classic alternating {true, false}.
std::vector<bool> pattern_of(const Branch<Value>* b);

// Optional input pin (absent in an older document = empty).
inline const Tree& input_or_empty(const std::vector<Tree>& in, std::size_t i) {
  static const Tree empty;
  return i < in.size() ? in[i] : empty;
}
inline bool has_values(const Tree& t) { return t.first() != nullptr; }

// The first value of alternative T in the tree, BORROWED (nullptr when none).
// Lvalue trees only, like DataTree::first(): the replacement for the
// `for (v : t.non_nulls()) p = v.get<T>();` pattern, whose pointer outlives
// the owning copy non_nulls() returns.
template <class T>
const T* first_of(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (const T* x = o->template get<T>()) return x;
  return nullptr;
}
template <class T>
const T* first_of(const Tree&&) = delete;

// Lines of beams (or bare lines) flowing into a bar-load "B" pin.
std::vector<Line> beam_lines(const Tree& t, bool accept_lines = false);

// Typed accessor failure (names the pin and the received type).
[[noreturn]] void wrong_type(const std::string& pin, const std::string& expected, const Value& got);

// ---------------------------------------------------------------------------
// NodeMerger — the canonical coincident-point clusterer.
//
// Merges points within `tol` and hands out stable 1-based node ids. Welding is
// ~O(n) through a spatial hash of cubic cells of side `tol`: any match lies in
// the 3x3x3 block of cells around the query, so the result is identical to a
// linear scan: `id_for` returns the smallest matching index, `find` the nearest
// match with ties to the smallest index. A non-positive / NaN tol falls back to
// the linear scan. Non-finite coordinates map to cell 0 and never deduplicate;
// Model::validate() rejects them downstream.
// ---------------------------------------------------------------------------
class NodeMerger {
 public:
  explicit NodeMerger(double tol) : tol_(tol) {}
  Id id_for(const Vec3& p);                    // existing id within tol, else a new node
  std::optional<Id> find(const Vec3& p) const;  // nearest existing node within tol
  const std::vector<Vec3>& positions() const { return pos_; }
  double tol() const { return tol_; }

 private:
  using Cell = std::tuple<std::int64_t, std::int64_t, std::int64_t>;
  Cell cell_of(const Vec3& p) const;
  std::vector<std::uint32_t> candidates(const Vec3& p) const;

  double tol_;
  std::vector<Vec3> pos_;
  std::map<Cell, std::vector<std::uint32_t>> grid_;
};

// The bar whose welded end nodes coincide with the line ends (either
// orientation) within tol, with `reversed` = the line runs from nj to ni.
// Assumes node ids are 1-based dense (Assemble's numbering).
struct BarMatch {
  std::size_t index = 0;
  bool reversed = false;
};
std::optional<BarMatch> find_bar(const Model& m, const Line& l, double tol);

}  // namespace keta::param

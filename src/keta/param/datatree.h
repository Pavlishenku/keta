// keta::param — the DataTree carried by every wire of the parametric layer
// (Grasshopper `GH_Structure<T>` semantics).
//
// One container, three shapes: a single item is a tree with one branch of one
// item, a list is one branch of N items, a tree is many branches. A branch is a
// dense list of NULLABLE items (std::nullopt is a first-class Grasshopper null).
//
// Invariants: paths are kept sorted (std::map, lexicographic order = GH_Path
// comparison); paths may be sparse ({0} and {5} with nothing between is legal);
// depths may be ragged (constant depth is conventional, not enforced).
//
// Every operation is deterministic: emission order is path-sorted, then item
// order. Restructuring operations are immutable (const -> new tree).
//
// Semantics notes (verified against the McNeel documentation in the source):
//   * simplify never strips so much that a path becomes empty, and a tree with
//     <= 1 branch is returned unchanged (nothing is "shared by all branches");
//   * shift_paths with a positive offset removes from the outer (left) end — the
//     exact GH reordering is under-documented, this is the documented
//     approximation; a negative offset is exactly trim(|offset|);
//   * matching: longest list repeats the LAST item (and the last branch), shortest
//     list drops the excess, cross reference is the cartesian product with
//     concatenated source paths. Output paths come from the principal input (the
//     tree with more branches).
#pragma once
#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "keta/cancel.h"

namespace keta::param {

// ---------------------------------------------------------------------------
// Path
// ---------------------------------------------------------------------------

// A branch address in structure space, e.g. {0;1;2}.
class Path {
 public:
  Path() = default;
  Path(std::initializer_list<int> idx) : idx_(idx) {}
  explicit Path(std::vector<int> idx) : idx_(std::move(idx)) {}

  static Path single(int i) { return Path({i}); }

  std::size_t size() const { return idx_.size(); }
  bool empty() const { return idx_.empty(); }
  int operator[](std::size_t i) const { return idx_[i]; }
  const std::vector<int>& indices() const { return idx_; }

  // GH `Valid`: false if empty or any index negative.
  bool is_valid() const {
    return !idx_.empty() && std::all_of(idx_.begin(), idx_.end(), [](int i) { return i >= 0; });
  }

  // GH AppendElement / PrependElement / CullElement (non-mutating).
  Path appended(int i) const {
    Path p = *this;
    p.idx_.push_back(i);
    return p;
  }
  Path prepended(int i) const {
    Path p;
    p.idx_.reserve(idx_.size() + 1);
    p.idx_.push_back(i);
    p.idx_.insert(p.idx_.end(), idx_.begin(), idx_.end());
    return p;
  }
  Path culled() const {
    Path p = *this;
    if (!p.idx_.empty()) p.idx_.pop_back();
    return p;
  }
  // Sub-path [from, to).
  Path slice(std::size_t from, std::size_t to) const {
    return Path(std::vector<int>(idx_.begin() + static_cast<std::ptrdiff_t>(from),
                                 idx_.begin() + static_cast<std::ptrdiff_t>(to)));
  }
  Path concat(const Path& other) const {
    Path p = *this;
    p.idx_.insert(p.idx_.end(), other.idx_.begin(), other.idx_.end());
    return p;
  }

  // GH IsAncestor: strict prefix.
  bool is_ancestor(const Path& other) const {
    return idx_.size() < other.idx_.size() &&
           std::equal(idx_.begin(), idx_.end(), other.idx_.begin());
  }

  // {0;1;2}; the empty path renders as {}.
  std::string to_string() const {
    std::string s = "{";
    for (std::size_t k = 0; k < idx_.size(); ++k) {
      if (k) s += ';';
      s += std::to_string(idx_[k]);
    }
    return s + "}";
  }

  friend bool operator==(const Path& a, const Path& b) { return a.idx_ == b.idx_; }
  friend auto operator<=>(const Path& a, const Path& b) { return a.idx_ <=> b.idx_; }

 private:
  std::vector<int> idx_;
};

template <class T>
using Branch = std::vector<std::optional<T>>;

// ---------------------------------------------------------------------------
// DataTree
// ---------------------------------------------------------------------------

template <class T>
class DataTree {
 public:
  using Map = std::map<Path, Branch<T>>;

  DataTree() = default;

  // Convenience: a one-item tree {0}[v] (the shape of a single value).
  static DataTree item(T v) {
    DataTree t;
    t.append(Path::single(0), std::move(v));
    return t;
  }
  // Convenience: a one-branch list {0}[v...].
  static DataTree list(std::vector<T> vs) {
    DataTree t;
    auto& b = t.ensure_path(Path::single(0));
    for (auto& v : vs) b.emplace_back(std::move(v));
    return t;
  }

  bool empty() const { return branches_.empty(); }
  std::size_t path_count() const { return branches_.size(); }  // GH PathCount
  std::size_t item_count() const {                              // GH DataCount (nulls included)
    std::size_t n = 0;
    for (const auto& kv : branches_) n += kv.second.size();
    return n;
  }

  // GH EnsurePath: create an empty branch if absent.
  Branch<T>& ensure_path(const Path& p) { return branches_[p]; }
  // Append one non-null item, creating the branch if needed.
  void append(const Path& p, T item) { branches_[p].emplace_back(std::move(item)); }
  void append_null(const Path& p) { branches_[p].emplace_back(std::nullopt); }

  // Borrowing accessors (branch, branches, first): they hand out pointers /
  // references INTO the tree and are therefore callable on lvalues only — a
  // call on a temporary tree does not compile.
  // nullptr when the path is absent (an existing empty branch is non-null).
  const Branch<T>* branch(const Path& p) const& {
    auto it = branches_.find(p);
    return it == branches_.end() ? nullptr : &it->second;
  }
  const Branch<T>* branch(const Path& p) const&& = delete;

  const Map& branches() const& { return branches_; }
  const Map& branches() const&& = delete;
  std::vector<Path> paths() const {
    std::vector<Path> v;
    v.reserve(branches_.size());
    for (const auto& kv : branches_) v.push_back(kv.first);
    return v;
  }

  // GH NonNulls: every present value in (path, item) order, as an OWNING copy.
  // Deliberately not a view: under C++20 a range-for over
  // `solve(...)[0].non_nulls()` destroys the temporary container before the
  // loop body runs, so a view would dangle there while the copy is safe
  // (P2718, C++23, fixes the rule). Never keep a pointer or reference into the
  // returned vector beyond the statement that produced it — copy what must
  // outlive it.
  [[nodiscard]] std::vector<T> non_nulls() const {
    std::vector<T> v;
    for (const auto& kv : branches_)
      for (const auto& o : kv.second)
        if (o) v.push_back(*o);
    return v;
  }
  // First present value: a BORROW into this tree (nullptr when none), valid
  // while the tree is alive and unmodified — no copy is made. Returning a copy
  // made `const X& r = t.first()->member;` dangle (the copy dies at the end of
  // the full-expression). Calling it on a
  // temporary tree is a compile error, the C++ stand-in for the borrow check.
  [[nodiscard]] const T* first() const& {
    for (const auto& kv : branches_)
      for (const auto& o : kv.second)
        if (o) return &*o;
    return nullptr;
  }
  const T* first() const&& = delete;

  // Map every item, preserving the branch structure and the null positions.
  template <class F>
  auto map(F&& f) const -> DataTree<std::decay_t<decltype(f(std::declval<const T&>()))>> {
    using U = std::decay_t<decltype(f(std::declval<const T&>()))>;
    DataTree<U> out;
    for (const auto& [p, b] : branches_) {
      auto& ob = out.ensure_path(p);
      ob.reserve(b.size());
      for (const auto& o : b) {
        if (o)
          ob.emplace_back(f(*o));
        else
          ob.emplace_back(std::nullopt);
      }
    }
    return out;
  }

  // --- restructuring -------------------------------------------------------

  // Graft: each item into its own child branch {path;i}. Empty branches vanish.
  DataTree graft() const {
    DataTree out;
    for (const auto& [p, b] : branches_)
      for (std::size_t i = 0; i < b.size(); ++i)
        out.branches_[p.appended(static_cast<int>(i))] = Branch<T>{b[i]};
    return out;
  }

  // Flatten: every item, path-sorted then item order, into the single branch {0}.
  // An empty tree (no items) stays empty.
  DataTree flatten() const {
    Branch<T> items;
    items.reserve(item_count());
    for (const auto& kv : branches_) items.insert(items.end(), kv.second.begin(), kv.second.end());
    DataTree out;
    if (!items.empty()) out.branches_[Path::single(0)] = std::move(items);
    return out;
  }

  // Trim: drop the last `depth` indices of every path, merging collisions.
  DataTree trim(std::size_t depth) const {
    return path_map([depth](const Path& p) {
      std::size_t keep = p.size() > depth ? p.size() - depth : 0;
      return p.slice(0, keep);
    });
  }

  // Simplify: strip the leading (and, unless `front`, trailing) index runs
  // shared by ALL branches; at least one index is always kept.
  DataTree simplify(bool front = false) const {
    if (branches_.size() <= 1) return *this;
    std::vector<const Path*> ps;
    std::size_t min_len = SIZE_MAX;
    for (const auto& kv : branches_) {
      ps.push_back(&kv.first);
      min_len = std::min(min_len, kv.first.size());
    }
    std::size_t lead = common_leading(ps);
    std::size_t trail = front ? 0 : common_trailing(ps);
    std::size_t allowed = min_len > 0 ? min_len - 1 : 0;
    if (lead > allowed) lead = allowed;
    if (lead + trail > allowed) trail = allowed - lead;
    return path_map([lead, trail](const Path& p) { return p.slice(lead, p.size() - trail); });
  }

  // Flip Matrix: transpose branches x items; ragged branches pad with nulls;
  // output paths renumbered {0}..{L-1}.
  DataTree flip_matrix() const {
    std::size_t l = 0;
    for (const auto& kv : branches_) l = std::max(l, kv.second.size());
    DataTree out;
    for (std::size_t j = 0; j < l; ++j) {
      Branch<T> col;
      col.reserve(branches_.size());
      for (const auto& kv : branches_)
        col.push_back(j < kv.second.size() ? kv.second[j] : std::nullopt);
      out.branches_[Path::single(static_cast<int>(j))] = std::move(col);
    }
    return out;
  }

  // Clean Tree: remove nulls (compacting), then empty branches. Paths unchanged.
  DataTree clean_tree(bool remove_nulls = true, bool remove_empty = true) const {
    DataTree out;
    for (const auto& [p, b] : branches_) {
      Branch<T> c;
      if (remove_nulls) {
        for (const auto& o : b)
          if (o) c.push_back(o);
      } else {
        c = b;
      }
      if (remove_empty && c.empty()) continue;
      out.branches_[p] = std::move(c);
    }
    return out;
  }

  // Prune Tree: drop a branch of count c if c < min or (max > 0 and c > max).
  DataTree prune_tree(std::size_t min, std::size_t max) const {
    DataTree out;
    for (const auto& [p, b] : branches_) {
      std::size_t c = b.size();
      if (c < min || (max > 0 && c > max)) continue;
      out.branches_[p] = b;
    }
    return out;
  }

  // Path Mapper (path level): rewrite every path, appending colliding branches
  // in path-sorted then item order (never overwriting).
  template <class F>
  DataTree path_map(F&& f) const {
    DataTree out;
    for (const auto& [p, b] : branches_) {
      auto& ob = out.branches_[f(p)];
      ob.insert(ob.end(), b.begin(), b.end());
    }
    return out;
  }

  // Shift Paths: positive offset removes from the outer end, negative = trim.
  DataTree shift_paths(int offset) const {
    std::size_t n = static_cast<std::size_t>(offset < 0 ? -static_cast<std::int64_t>(offset) : offset);
    return path_map([offset, n](const Path& p) {
      if (offset >= 0) return p.slice(std::min(n, p.size()), p.size());
      return p.slice(0, p.size() > n ? p.size() - n : 0);
    });
  }

  // Merge: union of paths; shared paths concatenate this ++ other.
  DataTree merge(const DataTree& other) const {
    DataTree out = *this;
    for (const auto& [p, b] : other.branches_) {
      auto& ob = out.branches_[p];
      ob.insert(ob.end(), b.begin(), b.end());
    }
    return out;
  }

  // Entwine: each input flattened to {0;i} (an empty input -> empty branch).
  static DataTree entwine(const std::vector<DataTree>& inputs) {
    DataTree out;
    for (std::size_t i = 0; i < inputs.size(); ++i) {
      auto& ob = out.branches_[Path({0, static_cast<int>(i)})];
      for (const auto& kv : inputs[i].branches_) ob.insert(ob.end(), kv.second.begin(), kv.second.end());
    }
    return out;
  }

  // Reverse the item order within each branch (nulls travel).
  DataTree reverse_items() const {
    DataTree out = *this;
    for (auto& kv : out.branches_) std::reverse(kv.second.begin(), kv.second.end());
    return out;
  }

  // Unflatten: pour this tree's items into the branch shape of `guide`.
  template <class U>
  DataTree unflatten(const DataTree<U>& guide) const {
    Branch<T> data;
    for (const auto& kv : branches_) data.insert(data.end(), kv.second.begin(), kv.second.end());
    std::size_t k = 0;
    DataTree out;
    for (const auto& [p, gb] : guide.branches()) {
      auto& ob = out.branches_[p];
      for (std::size_t i = 0; i < gb.size() && k < data.size(); ++i) ob.push_back(data[k++]);
    }
    return out;
  }

  friend bool operator==(const DataTree& a, const DataTree& b) { return a.branches_ == b.branches_; }

 private:
  static std::size_t common_leading(const std::vector<const Path*>& ps) {
    const Path& f = *ps[0];
    std::size_t l = 0;
    for (; l < f.size(); ++l)
      for (const Path* p : ps)
        if (l >= p->size() || (*p)[l] != f[l]) return l;
    return l;
  }
  static std::size_t common_trailing(const std::vector<const Path*>& ps) {
    const Path& f = *ps[0];
    std::size_t t = 0;
    for (; t < f.size(); ++t)
      for (const Path* p : ps)
        if (t >= p->size() || (*p)[p->size() - 1 - t] != f[f.size() - 1 - t]) return t;
    return t;
  }

  Map branches_;
};

// ---------------------------------------------------------------------------
// Data matching
// ---------------------------------------------------------------------------

// How a component pin consumes its incoming tree.
enum class AccessMode { Item, List, Tree };
enum class MatchMode { Shortest, Longest, CrossReference };

template <class A, class B>
using Pair = std::pair<std::optional<A>, std::optional<B>>;

// Longest list on two branches: length = max, saturating (repeat-last) indexing;
// an absent / empty side is null throughout.
template <class A, class B>
std::vector<Pair<A, B>> zip_longest_branch(const Branch<A>* a, const Branch<B>* b) {
  std::size_t al = a ? a->size() : 0, bl = b ? b->size() : 0;
  std::size_t n = std::max(al, bl);
  std::vector<Pair<A, B>> out;
  out.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    std::optional<A> av = al ? (*a)[std::min(i, al - 1)] : std::nullopt;
    std::optional<B> bv = bl ? (*b)[std::min(i, bl - 1)] : std::nullopt;
    out.emplace_back(std::move(av), std::move(bv));
  }
  return out;
}

// Shortest list on two branches: length = min.
template <class A, class B>
std::vector<Pair<A, B>> zip_shortest_branch(const Branch<A>& a, const Branch<B>& b) {
  std::size_t n = std::min(a.size(), b.size());
  std::vector<Pair<A, B>> out;
  out.reserve(n);
  for (std::size_t i = 0; i < n; ++i) out.emplace_back(a[i], b[i]);
  return out;
}

// Tree-to-tree longest list (the Grasshopper default): positional branch
// alignment with last-branch repeat, then per-branch longest list.
template <class A, class B>
DataTree<Pair<A, B>> zip_longest(const DataTree<A>& a, const DataTree<B>& b) {
  std::vector<std::pair<const Path*, const Branch<A>*>> av;
  std::vector<std::pair<const Path*, const Branch<B>*>> bv;
  for (const auto& kv : a.branches()) av.emplace_back(&kv.first, &kv.second);
  for (const auto& kv : b.branches()) bv.emplace_back(&kv.first, &kv.second);
  std::size_t n = std::max(av.size(), bv.size());
  DataTree<Pair<A, B>> out;
  if (n == 0) return out;
  bool use_b = bv.size() > av.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Branch<A>* ab = av.empty() ? nullptr : av[std::min(i, av.size() - 1)].second;
    const Branch<B>* bb = bv.empty() ? nullptr : bv[std::min(i, bv.size() - 1)].second;
    const Path& p = use_b ? *bv[std::min(i, bv.size() - 1)].first : *av[std::min(i, av.size() - 1)].first;
    auto& ob = out.ensure_path(p);
    for (auto& pr : zip_longest_branch(ab, bb)) ob.emplace_back(std::move(pr));
  }
  return out;
}

// Tree-to-tree shortest list: branch alignment stops at the smaller count.
template <class A, class B>
DataTree<Pair<A, B>> zip_shortest(const DataTree<A>& a, const DataTree<B>& b) {
  std::vector<std::pair<const Path*, const Branch<A>*>> av;
  std::vector<std::pair<const Path*, const Branch<B>*>> bv;
  for (const auto& kv : a.branches()) av.emplace_back(&kv.first, &kv.second);
  for (const auto& kv : b.branches()) bv.emplace_back(&kv.first, &kv.second);
  std::size_t n = std::min(av.size(), bv.size());
  bool use_b = bv.size() > av.size();
  DataTree<Pair<A, B>> out;
  for (std::size_t i = 0; i < n; ++i) {
    auto& ob = out.ensure_path(use_b ? *bv[i].first : *av[i].first);
    for (auto& pr : zip_shortest_branch(*av[i].second, *bv[i].second)) ob.emplace_back(std::move(pr));
  }
  return out;
}

// Cross reference: every item of a against every item of b; one output branch
// per pair of non-empty input branches at the concatenated path; a outer, b inner.
template <class A, class B>
DataTree<Pair<A, B>> cross_reference(const DataTree<A>& a, const DataTree<B>& b) {
  DataTree<Pair<A, B>> out;
  for (const auto& [pa, ba] : a.branches())
    for (const auto& [pb, bb] : b.branches()) {
      if (ba.empty() || bb.empty()) continue;
      auto& ob = out.ensure_path(pa.concat(pb));
      for (const auto& ia : ba)
        for (const auto& ib : bb) ob.emplace_back(Pair<A, B>(ia, ib));
    }
  return out;
}

// ---------------------------------------------------------------------------
// Cooperative cancellation
// ---------------------------------------------------------------------------

// A shared cancel flag polled cooperatively. Copies share the same flag.
// The trees a component's solve() returns, one per output pin. A thin
// wrapper over std::vector<DataTree<T>> whose element access is ref-qualified:
// on a NAMED Outputs `out[i]` is a reference into it, while on a temporary one
// (`node.solve(ctx, in)[i]`) it MOVES the tree out and returns it by value.
// Combined with the lvalue-only borrowing accessors of DataTree, a pointer
// into the result of an unnamed solve() cannot be written — it was the one
// borrow the compiler could not see.
template <class T>
class Outputs {
 public:
  using Tree = DataTree<T>;
  Outputs() = default;
  Outputs(std::initializer_list<Tree> trees) : v_(trees) {}
  Outputs(std::vector<Tree> trees) : v_(std::move(trees)) {}  // NOLINT: implicit by design

  std::size_t size() const { return v_.size(); }
  bool empty() const { return v_.empty(); }
  const Tree& operator[](std::size_t i) const& { return v_[i]; }
  Tree& operator[](std::size_t i) & { return v_[i]; }
  Tree operator[](std::size_t i) && { return std::move(v_[i]); }
  const Tree& at(std::size_t i) const& { return v_.at(i); }
  Tree at(std::size_t i) && { return std::move(v_.at(i)); }
  const Tree& front() const& { return v_.front(); }
  Tree front() && { return std::move(v_.front()); }
  const Tree& back() const& { return v_.back(); }
  Tree back() && { return std::move(v_.back()); }
  // Iteration stays available on temporaries: a range-for keeps its range
  // object alive for the whole loop.
  auto begin() const { return v_.begin(); }
  auto end() const { return v_.end(); }
  const std::vector<Tree>& trees() const& { return v_; }
  std::vector<Tree> release() && { return std::move(v_); }
  friend bool operator==(const Outputs&, const Outputs&) = default;

 private:
  std::vector<Tree> v_;
};

// The engine's token (keta/cancel.h): one flag cancels the graph walk and
// the hooked driver inside a node.
using CancelToken = keta::CancelToken;

}  // namespace keta::param

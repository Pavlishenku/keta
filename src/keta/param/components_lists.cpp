// Lists and Trees component families.
#include <algorithm>
#include <cmath>
#include <map>

#include "keta/param/components.h"

namespace keta::param {

namespace {

const char* const kLists = "Lists";
const char* const kTree = "Tree";

double num(const Value& v) {
  if (auto x = v.as_number()) return *x;
  throw ComponentError("expected a number, got " + v.type_tag());
}

std::int64_t opt_int(const Tree& t, std::int64_t dflt) { return first_int(t).value_or(dflt); }
bool opt_bool(const Tree& t, bool dflt) { return first_bool(t).value_or(dflt); }
std::size_t opt_usize(const Tree& t, std::int64_t dflt) {
  return static_cast<std::size_t>(std::max<std::int64_t>(0, opt_int(t, dflt)));
}

const Tree& in_at(const std::vector<Tree>& in, std::size_t i) { return input_or_empty(in, i); }

}  // namespace

// ============================ Lists ============================

ComponentInfo ListItemNode::info() const {
  return {"keta.list_item", "List Item", kLists, "The item at index i of each branch (out-of-range -> null)"};
}
std::vector<PinSpec> ListItemNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::item("i", "int", "Index (default 0)").opt()};
}
std::vector<PinSpec> ListItemNode::outputs() const { return {PinSpec::item("V", "any", "The item at index i")}; }
Outputs<Value> ListItemNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::int64_t idx = opt_int(in_at(in, 1), 0);
  Tree out;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    if (idx >= 0 && static_cast<std::size_t>(idx) < branch.size())
      out.ensure_path(path).push_back(branch[static_cast<std::size_t>(idx)]);
    else
      out.append_null(path);
  }
  return {out};
}

ComponentInfo ListLengthNode::info() const {
  return {"keta.list_length", "List Length", kLists, "The number of items in each branch"};
}
std::vector<PinSpec> ListLengthNode::inputs() const { return {PinSpec::list("L", "any", "List").opt()}; }
std::vector<PinSpec> ListLengthNode::outputs() const { return {PinSpec::item("N", "int", "Item count of each branch")}; }
Outputs<Value> ListLengthNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree out;
  for (const auto& [path, branch] : in_at(in, 0).branches())
    out.append(path, Value::integer(static_cast<std::int64_t>(branch.size())));
  return {out};
}

ComponentInfo SubListNode::info() const {
  return {"keta.sub_list", "Sub List", kLists, "The inclusive, clamped slice [S0, S1] of each branch"};
}
std::vector<PinSpec> SubListNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::item("S0", "int", "Start index (default 0)").opt(),
          PinSpec::item("S1", "int", "End index, inclusive (default: last)").opt()};
}
std::vector<PinSpec> SubListNode::outputs() const { return {PinSpec::list("L", "any", "The sub-list of each branch")}; }
Outputs<Value> SubListNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::int64_t start = opt_int(in_at(in, 1), 0);
  const bool end_wired = has_values(in_at(in, 2));
  const std::int64_t end_opt = opt_int(in_at(in, 2), 0);
  Tree out;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    const auto len = static_cast<std::int64_t>(branch.size());
    std::int64_t s = std::clamp<std::int64_t>(start, 0, len);
    std::int64_t e = std::clamp<std::int64_t>(end_wired ? end_opt : len - 1, -1, len - 1);
    auto& dst = out.ensure_path(path);
    for (std::int64_t k = s; k <= e; ++k) dst.push_back(branch[static_cast<std::size_t>(k)]);
  }
  return {out};
}

ComponentInfo ShiftListNode::info() const {
  return {"keta.list_shift", "Shift List", kLists, "Shift each branch by S (W wraps, default true)"};
}
std::vector<PinSpec> ShiftListNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::item("S", "int", "Shift offset (default 1)").opt(),
          PinSpec::item("W", "bool", "Wrap around (default true)").opt()};
}
std::vector<PinSpec> ShiftListNode::outputs() const { return {PinSpec::list("L", "any", "The shifted branch")}; }
Outputs<Value> ShiftListNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::int64_t shift = opt_int(in_at(in, 1), 1);
  const bool wrap = opt_bool(in_at(in, 2), true);
  Tree out;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    auto& dst = out.ensure_path(path);
    const std::size_t len = branch.size();
    if (len == 0) continue;
    if (wrap) {
      const auto n = static_cast<std::int64_t>(len);
      for (std::size_t i = 0; i < len; ++i) {
        std::int64_t src = (static_cast<std::int64_t>(i) + shift % n + n) % n;  // rem_euclid
        dst.push_back(branch[static_cast<std::size_t>(src)]);
      }
    } else if (shift >= 0) {
      for (std::size_t i = std::min<std::size_t>(static_cast<std::size_t>(shift), len); i < len; ++i)
        dst.push_back(branch[i]);
    } else {
      // Negating INT64_MIN would overflow: anything that large drops everything.
      std::size_t k = shift < -static_cast<std::int64_t>(len) ? len : static_cast<std::size_t>(-shift);
      for (std::size_t i = 0; i < len - k; ++i) dst.push_back(branch[i]);
    }
  }
  return {out};
}

ComponentInfo ListReverseNode::info() const {
  return {"keta.list_reverse", "Reverse List", kLists, "Reverse the items of each branch (nulls travel)"};
}
std::vector<PinSpec> ListReverseNode::inputs() const { return {PinSpec::list("L", "any", "List").opt()}; }
std::vector<PinSpec> ListReverseNode::outputs() const { return {PinSpec::list("L", "any", "The reversed branch")}; }
Outputs<Value> ListReverseNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).reverse_items()};
}

ComponentInfo ListInsertNode::info() const {
  return {"keta.list_insert", "Insert Items", kLists, "Insert items I into each branch at the (clamped) indices i"};
}
std::vector<PinSpec> ListInsertNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::list("I", "any", "Items to insert").opt(),
          PinSpec::list("i", "int", "Indices to insert at").opt()};
}
std::vector<PinSpec> ListInsertNode::outputs() const {
  return {PinSpec::list("L", "any", "The branch with items inserted")};
}
Outputs<Value> ListInsertNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const BranchIndex items_ix(in_at(in, 1)), idx_ix(in_at(in, 2));
  Tree out;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    const Branch<Value>* items = items_ix.at(bi);
    const Branch<Value>* idxs = idx_ix.at(bi);
    ++bi;
    auto& dst = out.ensure_path(path);
    if (!items || !idxs || items->empty() || idxs->empty()) {
      dst.insert(dst.end(), branch.begin(), branch.end());
      continue;
    }
    const std::size_t n = std::max(items->size(), idxs->size());
    const auto len = static_cast<std::int64_t>(branch.size());
    std::vector<Branch<Value>> buckets(branch.size() + 1);  // items inserted BEFORE each position
    for (std::size_t k = 0; k < n; ++k) {
      const auto& item = (*items)[std::min(k, items->size() - 1)];
      const auto& ix = (*idxs)[std::min(k, idxs->size() - 1)];
      std::optional<std::int64_t> i = ix ? ix->as_int() : std::nullopt;
      if (!i) throw ComponentError("insert index must be an integer");
      buckets[static_cast<std::size_t>(std::clamp<std::int64_t>(*i, 0, len))].push_back(item);
    }
    for (std::size_t j = 0; j < branch.size(); ++j) {
      dst.insert(dst.end(), buckets[j].begin(), buckets[j].end());
      dst.push_back(branch[j]);
    }
    dst.insert(dst.end(), buckets[branch.size()].begin(), buckets[branch.size()].end());
  }
  return {out};
}

ComponentInfo SortListNode::info() const {
  return {"keta.sort_list", "Sort List", kLists, "Sort each branch ascending (by keys K when wired; nulls last)"};
}
std::vector<PinSpec> SortListNode::inputs() const {
  return {PinSpec::list("L", "number", "List to sort").opt(), PinSpec::list("K", "number", "Optional sort keys").opt()};
}
std::vector<PinSpec> SortListNode::outputs() const {
  return {PinSpec::list("L", "any", "The sorted branch"),
          PinSpec::tree("I", "int", "Source index of each sorted item, per branch (nulls carry their source position)")};
}
Outputs<Value> SortListNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const bool has_keys = !in_at(in, 1).empty();
  const BranchIndex keys_ix(in_at(in, 1));
  Tree out, idx_out;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    const Branch<Value>* keys = has_keys ? keys_ix.at(bi) : nullptr;
    ++bi;
    struct Entry {
      double key;
      std::size_t j;
    };
    std::vector<Entry> e;
    e.reserve(branch.size());
    for (std::size_t j = 0; j < branch.size(); ++j) {
      const std::optional<Value>& slot = keys && !keys->empty() ? (*keys)[std::min(j, keys->size() - 1)] : branch[j];
      double k = slot ? num(*slot) : HUGE_VAL;  // nulls sort last
      if (std::isnan(k)) k = HUGE_VAL;           // keeps a strict weak order
      e.push_back({k, j});
    }
    // Stable: equal keys keep their source order (the I map stays deterministic).
    std::stable_sort(e.begin(), e.end(), [](const Entry& a, const Entry& b) { return a.key < b.key; });
    auto& dst = out.ensure_path(path);
    auto& idst = idx_out.ensure_path(path);
    for (const Entry& x : e) {
      dst.push_back(branch[x.j]);
      idst.emplace_back(Value::integer(static_cast<std::int64_t>(x.j)));
    }
  }
  return {out, idx_out};
}

ComponentInfo CullPatternNode::info() const {
  return {"keta.cull_pattern", "Cull Pattern", kLists, "Keep items where the repeated boolean pattern P is true"};
}
std::vector<PinSpec> CullPatternNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::list("P", "bool", "Cull pattern (default true,false)").opt()};
}
std::vector<PinSpec> CullPatternNode::outputs() const { return {PinSpec::list("L", "any", "The culled branch")}; }
Outputs<Value> CullPatternNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const BranchIndex pat_ix(in_at(in, 1));
  Tree out;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    auto pattern = pattern_of(pat_ix.at(bi++));
    auto& dst = out.ensure_path(path);
    for (std::size_t j = 0; j < branch.size(); ++j)
      if (pattern[j % pattern.size()]) dst.push_back(branch[j]);
  }
  return {out};
}

ComponentInfo DispatchNode::info() const {
  return {"keta.dispatch", "Dispatch", kLists, "Split each branch by the pattern P into A (true) and B (false)"};
}
std::vector<PinSpec> DispatchNode::inputs() const {
  return {PinSpec::list("L", "any", "List").opt(), PinSpec::list("P", "bool", "Dispatch pattern (default true,false)").opt()};
}
std::vector<PinSpec> DispatchNode::outputs() const {
  return {PinSpec::list("A", "any", "Items at a true slot"), PinSpec::list("B", "any", "Items at a false slot")};
}
Outputs<Value> DispatchNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const BranchIndex pat_ix(in_at(in, 1));
  Tree a, b;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    auto pattern = pattern_of(pat_ix.at(bi++));
    auto& da = a.ensure_path(path);
    auto& db = b.ensure_path(path);
    for (std::size_t j = 0; j < branch.size(); ++j) (pattern[j % pattern.size()] ? da : db).push_back(branch[j]);
  }
  return {a, b};
}

// ============================ Trees ============================

std::string tree_synopsis(const Tree& t) {
  const std::size_t branches = t.path_count();
  if (branches == 0) return "empty tree";
  const std::size_t items = t.item_count();
  std::size_t nulls = 0;
  std::map<std::string, std::size_t> counts;
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second) {
      if (o)
        ++counts[o->type_tag()];
      else
        ++nulls;
    }
  std::string out = std::to_string(branches) + (branches == 1 ? " branch, " : " branches, ") + std::to_string(items) +
                    (items == 1 ? " item" : " items");
  if (!counts.empty()) {
    std::vector<std::pair<std::string, std::size_t>> kinds(counts.begin(), counts.end());
    std::stable_sort(kinds.begin(), kinds.end(), [](const auto& x, const auto& y) { return x.second > y.second; });
    std::string dist;
    if (kinds.size() == 1) {
      dist = kinds[0].first;
    } else {
      for (std::size_t k = 0; k < kinds.size(); ++k)
        dist += (k ? ", " : "") + std::to_string(kinds[k].second) + " " + kinds[k].first;
    }
    out += " (" + dist + ")";
  }
  if (nulls > 0) out += ", " + std::to_string(nulls) + " null";
  constexpr std::size_t kCap = 30;
  std::size_t shown = 0;
  for (const auto& [path, branch] : t.branches()) {
    if (shown++ == kCap) break;
    out += "\n" + path.to_string() + " (N=" + std::to_string(branch.size()) + ")";
  }
  if (branches > kCap) out += "\n… and " + std::to_string(branches - kCap) + " more";
  return out;
}

ComponentInfo GraftNode::info() const {
  return {"keta.graft", "Graft", kTree, "Graft a tree: place every item in its own new child branch"};
}
std::vector<PinSpec> GraftNode::inputs() const { return {PinSpec::tree("T", "any", "Tree to graft").opt()}; }
std::vector<PinSpec> GraftNode::outputs() const { return {PinSpec::tree("T", "any", "Grafted tree")}; }
Outputs<Value> GraftNode::solve(const SolveCtx&, const std::vector<Tree>& in) const { return {in_at(in, 0).graft()}; }

ComponentInfo FlattenNode::info() const {
  return {"keta.flatten", "Flatten", kTree, "Flatten a tree: merge every branch into one list at {0}"};
}
std::vector<PinSpec> FlattenNode::inputs() const { return {PinSpec::tree("T", "any", "Tree to flatten").opt()}; }
std::vector<PinSpec> FlattenNode::outputs() const { return {PinSpec::tree("T", "any", "Flattened tree")}; }
Outputs<Value> FlattenNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).flatten()};
}

ComponentInfo SimplifyNode::info() const {
  return {"keta.simplify", "Simplify", kTree, "Simplify a tree: remove path segments shared by all branches"};
}
std::vector<PinSpec> SimplifyNode::inputs() const {
  return {PinSpec::tree("T", "any", "Tree to simplify").opt(),
          PinSpec::item("F", "bool", "Front-only: strip only the shared leading run (default false = both ends)").opt()};
}
std::vector<PinSpec> SimplifyNode::outputs() const { return {PinSpec::tree("T", "any", "Simplified tree")}; }
Outputs<Value> SimplifyNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).simplify(opt_bool(in_at(in, 1), false))};
}

ComponentInfo TrimTreeNode::info() const {
  return {"keta.trim_tree", "Trim Tree", kTree, "Trim a tree: drop the innermost D path levels, merging branches"};
}
std::vector<PinSpec> TrimTreeNode::inputs() const {
  return {PinSpec::tree("T", "any", "Tree to trim").opt(),
          PinSpec::item("D", "int", "Number of innermost path levels to remove (default 1)").opt()};
}
std::vector<PinSpec> TrimTreeNode::outputs() const { return {PinSpec::tree("T", "any", "Trimmed tree")}; }
Outputs<Value> TrimTreeNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).trim(opt_usize(in_at(in, 1), 1))};
}

ComponentInfo FlipMatrixNode::info() const {
  return {"keta.flip_matrix", "Flip Matrix", kTree, "Flip a matrix tree: swap branches with items (ragged rows padded null)"};
}
std::vector<PinSpec> FlipMatrixNode::inputs() const { return {PinSpec::tree("T", "any", "Tree to transpose").opt()}; }
std::vector<PinSpec> FlipMatrixNode::outputs() const { return {PinSpec::tree("T", "any", "Transposed tree")}; }
Outputs<Value> FlipMatrixNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).flip_matrix()};
}

ComponentInfo PruneTreeNode::info() const {
  return {"keta.prune_tree", "Prune Tree", kTree, "Prune a tree: remove branches whose item count is below N0 or above N1"};
}
std::vector<PinSpec> PruneTreeNode::inputs() const {
  return {PinSpec::tree("T", "any", "Tree to prune").opt(),
          PinSpec::item("N0", "int", "Minimum items to keep a branch (default 1)").opt(),
          PinSpec::item("N1", "int", "Maximum items to keep a branch (default 0 = no upper bound)").opt()};
}
std::vector<PinSpec> PruneTreeNode::outputs() const { return {PinSpec::tree("T", "any", "Pruned tree")}; }
Outputs<Value> PruneTreeNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).prune_tree(opt_usize(in_at(in, 1), 1), opt_usize(in_at(in, 2), 0))};
}

ComponentInfo ShiftPathsNode::info() const {
  return {"keta.shift_paths", "Shift Paths", kTree,
          "Shift paths: drop path levels and merge (default O=1 drops the outermost level)"};
}
std::vector<PinSpec> ShiftPathsNode::inputs() const {
  return {PinSpec::tree("T", "any", "Tree to shift").opt(),
          PinSpec::item("O", "int",
                        "Path levels to shift: +N drops the N outermost (left) levels, -N the innermost (default 1)")
              .opt()};
}
std::vector<PinSpec> ShiftPathsNode::outputs() const { return {PinSpec::tree("T", "any", "Shifted tree")}; }
Outputs<Value> ShiftPathsNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::int64_t o = std::clamp<std::int64_t>(opt_int(in_at(in, 1), 1), -1'000'000, 1'000'000);
  return {in_at(in, 0).shift_paths(static_cast<int>(o))};
}

ComponentInfo MergeTreesNode::info() const {
  return {"keta.merge_trees", "Merge", kTree, "Merge two trees (A then B); several wires into A merge into A as well"};
}
std::vector<PinSpec> MergeTreesNode::inputs() const {
  return {PinSpec::tree("A", "any", "First tree (wiring several sources here merges them all)").opt(),
          PinSpec::tree("B", "any", "Second tree, appended after A on shared paths").opt()};
}
std::vector<PinSpec> MergeTreesNode::outputs() const { return {PinSpec::tree("T", "any", "Merged tree")}; }
Outputs<Value> MergeTreesNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).merge(in_at(in, 1))};
}

ComponentInfo EntwineNode::info() const {
  return {"keta.entwine", "Entwine", kTree, "Entwine inputs: flatten each wired input into its own branch {0;i}"};
}
std::vector<PinSpec> EntwineNode::inputs() const {
  return {PinSpec::tree("A", "any", "First input (branch {0;0})").opt(),
          PinSpec::tree("B", "any", "Second input (branch {0;1})").opt(),
          PinSpec::tree("C", "any", "Third input (branch {0;2})").opt()};
}
std::vector<PinSpec> EntwineNode::outputs() const { return {PinSpec::tree("T", "any", "Entwined tree")}; }
Outputs<Value> EntwineNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::vector<Tree> wired;  // unwired (empty) inputs are skipped so the stamps pack
  for (const Tree& t : in)
    if (!t.empty()) wired.push_back(t);
  return {Tree::entwine(wired)};
}

ComponentInfo CleanTreeNode::info() const {
  return {"keta.clean_tree", "Clean Tree", kTree, "Clean a tree: drop null items (N) then empty branches (E)"};
}
std::vector<PinSpec> CleanTreeNode::inputs() const {
  return {PinSpec::tree("T", "any", "Tree to clean").opt(), PinSpec::item("N", "bool", "Remove null items (default true)").opt(),
          PinSpec::item("E", "bool", "Remove branches left empty (default true)").opt()};
}
std::vector<PinSpec> CleanTreeNode::outputs() const { return {PinSpec::tree("T", "any", "Cleaned tree")}; }
Outputs<Value> CleanTreeNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {in_at(in, 0).clean_tree(opt_bool(in_at(in, 1), true), opt_bool(in_at(in, 2), true))};
}

ComponentInfo ParamViewerNode::info() const {
  return {"keta.param_viewer", "Param Viewer", kTree,
          "Inspect a tree: branch/item counts, per-type distribution and per-branch sizes"};
}
std::vector<PinSpec> ParamViewerNode::inputs() const { return {PinSpec::tree("T", "any", "Tree to inspect").opt()}; }
std::vector<PinSpec> ParamViewerNode::outputs() const {
  return {PinSpec::item("S", "text", "Textual synopsis of the input tree")};
}
Outputs<Value> ParamViewerNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {single(Value::text(tree_synopsis(in_at(in, 0))))};
}

}  // namespace keta::param

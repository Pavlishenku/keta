// Sets family and the Catalog Section component.
//
// FLAT-INDEX SEMANTICS (shared by every set-aware node): a SetSpec addresses
// items by FLAT position — walk the branches in path order and count every slot,
// nulls included (DataTree::flatten order). The selectors' I output, Create
// Set, Filter by Set and Selection all count the same way, so a selection built
// against a list filters that same list exactly.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <set>

#include "keta/catalog.h"
#include "keta/param/components.h"

namespace keta::param {

namespace {

const char* const kSets = "Sets";
constexpr double kDefaultWeldTol = 1e-3;  // = Assemble's default tolerance

const Tree& in_at(const std::vector<Tree>& in, std::size_t i) { return input_or_empty(in, i); }
std::int64_t opt_int(const Tree& t, std::int64_t d) { return first_int(t).value_or(d); }
bool opt_bool(const Tree& t, bool d) { return first_bool(t).value_or(d); }
double opt_number(const Tree& t, double d) { return first_number(t).value_or(d); }
std::string opt_text(const Tree& t, const std::string& d) { return first_text(t).value_or(d); }

bool point_in_box(const Vec3& p, const Vec3& lo, const Vec3& hi) {
  return (p.array() >= lo.array()).all() && (p.array() <= hi.array()).all();
}

// Does the segment a-b touch the closed AABB [lo, hi]? (Liang-Barsky slab clip;
// a fully-inside segment counts.)
bool segment_intersects_box(const Vec3& a, const Vec3& b, const Vec3& lo, const Vec3& hi) {
  Vec3 d = b - a;
  double t0 = 0, t1 = 1;
  for (int k = 0; k < 3; ++k) {
    if (std::fabs(d[k]) < 1e-12) {
      if (a[k] < lo[k] || a[k] > hi[k]) return false;
    } else {
      double tn = (lo[k] - a[k]) / d[k], tf = (hi[k] - a[k]) / d[k];
      if (tn > tf) std::swap(tn, tf);
      t0 = std::max(t0, tn);
      t1 = std::min(t1, tf);
      if (t0 > t1) return false;
    }
  }
  return true;
}

double point_seg_dist(const Vec3& p, const Vec3& a, const Vec3& b) {
  Vec3 ab = b - a;
  double ab2 = ab.squaredNorm();
  if (ab2 <= 0) return (p - a).norm();
  double t = std::clamp((p - a).dot(ab) / ab2, 0.0, 1.0);
  return (p - (a + t * ab)).norm();
}

// The selector engine: S = matching items (paths kept, misses dropped, empty
// branches kept), I = flat indices of the matches at {0}, M = bool mask aligned
// 1:1 with the input (a null slot stays null). The flat counter advances over
// every slot, nulls included.
template <class Pred>
std::vector<Tree> select_by(const Tree& t, Pred pred) {
  Tree s, idx, mask;
  std::int64_t flat = 0;
  for (const auto& [path, branch] : t.branches()) {
    s.ensure_path(path);
    mask.ensure_path(path);
    for (const auto& slot : branch) {
      if (slot) {
        const bool hit = pred(*slot);
        mask.append(path, Value::boolean(hit));
        if (hit) {
          s.append(path, *slot);
          idx.append(Path::single(0), Value::integer(flat));
        }
      } else {
        mask.append_null(path);
      }
      ++flat;
    }
  }
  return {s, idx, mask};
}

std::vector<PinSpec> selector_outputs() {
  return {PinSpec::tree("S", "any", "Selected items (paths preserved, non-selected dropped)"),
          PinSpec::list("I", "int", "Flat indices of the selected items (branch {0})"),
          PinSpec::tree("M", "bool", "Selection mask aligned 1:1 with G (nulls stay null)")};
}

std::vector<Tree> empty_select() { return {Tree(), Tree(), Tree()}; }

const Vec3* point_like(const Value& v) {
  if (const auto* p = v.get<Point>()) return &p->p;
  if (const auto* q = v.get<Vector>()) return &q->v;
  return nullptr;
}

bool lines_match(const Line& a, const Line& b, double tol) {
  return ((a.a - b.a).norm() <= tol && (a.b - b.b).norm() <= tol) ||
         ((a.a - b.b).norm() <= tol && (a.b - b.a).norm() <= tol);
}

// Duplicate test: same type only; numbers / points / vectors / lines within
// tol (lines direction-agnostic), text and bool exact.
bool values_match(const Value& a, const Value& b, double tol) {
  if (const auto* x = a.get<double>())
    if (const auto* y = b.get<double>()) return std::fabs(*x - *y) <= tol;
  if (const auto* x = a.get<std::int64_t>())
    if (const auto* y = b.get<std::int64_t>()) return std::fabs(static_cast<double>(*x) - static_cast<double>(*y)) <= tol;
  if (const auto* x = a.get<Point>())
    if (const auto* y = b.get<Point>()) return (x->p - y->p).norm() <= tol;
  if (const auto* x = a.get<Vector>())
    if (const auto* y = b.get<Vector>()) return (x->v - y->v).norm() <= tol;
  if (const auto* x = a.get<Line>())
    if (const auto* y = b.get<Line>()) return lines_match(*x, *y, tol);
  if (const auto* x = a.get<std::string>())
    if (const auto* y = b.get<std::string>()) return *x == *y;
  if (const auto* x = a.get<bool>())
    if (const auto* y = b.get<bool>()) return *x == *y;
  return false;
}

bool iequals(const std::string& a, const std::string& b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin(), [](char x, char y) {
           return std::tolower(static_cast<unsigned char>(x)) == std::tolower(static_cast<unsigned char>(y));
         });
}

// Group keys: numbers and ints collapse (Int(2) and 2.0 share a group, keyed
// 2.0); total order numbers (total_cmp-like) < bools < texts.
struct GroupKey {
  int rank = 0;  // 0 number, 1 bool, 2 text
  double num = 0;
  bool b = false;
  std::string text;
};
GroupKey key_of(const Value& v) {
  if (const auto* x = v.get<double>()) return {0, *x, false, {}};
  if (const auto* i = v.get<std::int64_t>()) return {0, static_cast<double>(*i), false, {}};
  if (const auto* b = v.get<bool>()) return {1, 0, *b, {}};
  if (const auto* s = v.get<std::string>()) return {2, 0, false, *s};
  throw ComponentError("Group By: cannot group by a key of type " + v.type_tag());
}
// IEEE-754 totalOrder on doubles: -0 < +0, NaNs ordered.
std::int64_t total_key(double x) {
  std::int64_t bits;
  std::memcpy(&bits, &x, sizeof bits);
  return bits ^ static_cast<std::int64_t>(static_cast<std::uint64_t>(bits >> 63) >> 1);
}
int key_cmp(const GroupKey& a, const GroupKey& b) {
  if (a.rank != b.rank) return a.rank < b.rank ? -1 : 1;
  if (a.rank == 0) {
    auto x = total_key(a.num), y = total_key(b.num);
    return x < y ? -1 : (x > y ? 1 : 0);
  }
  if (a.rank == 1) return a.b == b.b ? 0 : (a.b ? 1 : -1);
  return a.text.compare(b.text) < 0 ? -1 : (a.text == b.text ? 0 : 1);
}
Value key_value(const GroupKey& k) {
  if (k.rank == 0) return Value::number(k.num);
  if (k.rank == 1) return Value::boolean(k.b);
  return Value::text(k.text);
}

const SetSpec& require_set(const Value& v, const char* what) {
  if (const auto* s = v.as_set()) return *s;
  throw ComponentError(std::string(what) + " must be a set");
}

}  // namespace

// ============================ cull / sort ============================

ComponentInfo CullIndexNode::info() const {
  return {"keta.cull_index", "Cull Index", kSets, "Drop the items at the given per-branch indices (K keeps them instead)"};
}
std::vector<PinSpec> CullIndexNode::inputs() const {
  return {PinSpec::tree("L", "any", "Items").opt(),
          PinSpec::list("I", "int", "Indices to cull, per branch (out-of-range ignored)").opt(),
          PinSpec::item("K", "bool", "Keep only the indices instead of culling (default false)").opt()};
}
std::vector<PinSpec> CullIndexNode::outputs() const { return {PinSpec::tree("L", "any", "The culled branch")}; }
Outputs<Value> CullIndexNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const bool keep = opt_bool(in_at(in, 2), false);
  const BranchIndex idx_ix(in_at(in, 1));
  Tree out;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    std::set<std::int64_t> cull;
    if (const Branch<Value>* idxs = idx_ix.at(bi))
      for (const auto& slot : *idxs) {
        if (!slot) continue;
        auto i = slot->as_int();
        if (!i) throw ComponentError("Cull Index: indices must be integers");
        if (*i < 0) throw ComponentError("Cull Index: negative index (indices are 0-based, use 0..len)");
        cull.insert(*i);
      }
    ++bi;
    auto& dst = out.ensure_path(path);
    for (std::size_t j = 0; j < branch.size(); ++j)
      if (keep == (cull.count(static_cast<std::int64_t>(j)) > 0)) dst.push_back(branch[j]);
  }
  return {out};
}

ComponentInfo CullDuplicatesNode::info() const {
  return {"keta.cull_duplicates", "Cull Duplicates", kSets, "Cull coincident items per branch (points/lines/numbers within tol)"};
}
std::vector<PinSpec> CullDuplicatesNode::inputs() const {
  return {PinSpec::tree("L", "any", "Items").opt(),
          PinSpec::item("T", "number", "Tolerance (m) for points/lines/numbers (default 0.001)").opt()};
}
std::vector<PinSpec> CullDuplicatesNode::outputs() const {
  return {PinSpec::tree("U", "any", "Unique items (paths preserved)"),
          PinSpec::tree("I", "int", "Representative index in U for each input item"),
          PinSpec::tree("C", "int", "Cluster count per unique item")};
}
Outputs<Value> CullDuplicatesNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const double tol = opt_number(in_at(in, 1), kDefaultWeldTol);
  Tree unique, remap, counts;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    std::vector<const Value*> reps;
    std::vector<std::int64_t> cluster;
    auto& rdst = remap.ensure_path(path);
    for (const auto& slot : branch) {
      if (!slot) {
        rdst.push_back(std::nullopt);
        continue;
      }
      std::size_t ri = 0;
      while (ri < reps.size() && !values_match(*reps[ri], *slot, tol)) ++ri;
      if (ri == reps.size()) {
        reps.push_back(&*slot);
        cluster.push_back(0);
      }
      ++cluster[ri];
      rdst.emplace_back(Value::integer(static_cast<std::int64_t>(ri)));
    }
    auto& udst = unique.ensure_path(path);
    for (const Value* r : reps) udst.emplace_back(*r);
    auto& cdst = counts.ensure_path(path);
    for (std::int64_t c : cluster) cdst.emplace_back(Value::integer(c));
  }
  return {unique, remap, counts};
}

ComponentInfo SortPointsNode::info() const {
  return {"keta.sort_points", "Sort Points", kSets, "Sort each branch of points (lexicographic or single-axis, stable)"};
}
std::vector<PinSpec> SortPointsNode::inputs() const {
  return {PinSpec::tree("P", "point", "Points").opt(),
          PinSpec::item("A", "int", "-1 = lexicographic (x,y,z); 0/1/2 = that axis only (default -1)").opt()};
}
std::vector<PinSpec> SortPointsNode::outputs() const {
  return {PinSpec::tree("P", "point", "The sorted branch"), PinSpec::tree("I", "int", "Source index of each sorted item, per branch")};
}
Outputs<Value> SortPointsNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::int64_t axis = opt_int(in_at(in, 1), -1);
  Tree sorted, idxmap;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    struct Item {
      std::size_t j;
      Vec3 p;
      const Value* v;
    };
    std::vector<Item> items;
    for (std::size_t j = 0; j < branch.size(); ++j) {
      if (!branch[j]) continue;
      auto p = branch[j]->as_point();
      if (!p) throw ComponentError("Sort Points: expected a point, got " + branch[j]->type_tag());
      items.push_back({j, *p, &*branch[j]});
    }
    auto less = [axis](const Item& a, const Item& b) {
      if (axis >= 0 && axis <= 2) return a.p[axis] < b.p[axis];
      for (int k = 0; k < 3; ++k)
        if (a.p[k] != b.p[k]) return a.p[k] < b.p[k];
      return false;
    };
    std::stable_sort(items.begin(), items.end(), less);
    auto& pd = sorted.ensure_path(path);
    auto& id = idxmap.ensure_path(path);
    for (const Item& it : items) {
      pd.emplace_back(*it.v);
      id.emplace_back(Value::integer(static_cast<std::int64_t>(it.j)));
    }
  }
  return {sorted, idxmap};
}

// ============================ selectors ============================

ComponentInfo SelectBoxNode::info() const {
  return {"keta.select_box", "Select by Box", kSets, "Select points/lines inside the axis-aligned box A-B (C = crossing mode)"};
}
std::vector<PinSpec> SelectBoxNode::inputs() const {
  return {PinSpec::tree("G", "any", "Points or lines").opt(), PinSpec::item("A", "point", "Box corner A"),
          PinSpec::item("B", "point", "Box corner B"),
          PinSpec::item("C", "bool", "Crossing mode: keep items touching the box, not only inside (default false)").opt()};
}
std::vector<PinSpec> SelectBoxNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectBoxNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  auto a = first_point(in_at(in, 1));
  if (!a) throw ComponentError("Select by Box: corner A is required");
  auto b = first_point(in_at(in, 2));
  if (!b) throw ComponentError("Select by Box: corner B is required");
  const bool crossing = opt_bool(in_at(in, 3), false);
  const Vec3 lo = a->cwiseMin(*b), hi = a->cwiseMax(*b);
  return select_by(in_at(in, 0), [&](const Value& v) {
    if (const Vec3* p = point_like(v)) return point_in_box(*p, lo, hi);
    if (const Line* l = v.as_line())
      return crossing ? segment_intersects_box(l->a, l->b, lo, hi) : point_in_box(l->a, lo, hi) && point_in_box(l->b, lo, hi);
    return false;
  });
}

ComponentInfo SelectLevelNode::info() const {
  return {"keta.select_level", "Select by Level", kSets, "Select points/lines within a coordinate band on one axis"};
}
std::vector<PinSpec> SelectLevelNode::inputs() const {
  return {PinSpec::tree("G", "any", "Points or lines").opt(), PinSpec::item("A", "int", "Axis: 0=X 1=Y 2=Z (default 2)").opt(),
          PinSpec::item("V", "number", "Level coordinate in m (default 0)").opt(),
          PinSpec::item("T", "number", "Band half-width in m (default 0.001)").opt()};
}
std::vector<PinSpec> SelectLevelNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectLevelNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  const std::int64_t axis = opt_int(in_at(in, 1), 2);
  if (axis < 0 || axis > 2) throw ComponentError("Select by Level: axis must be 0, 1 or 2");
  const double level = opt_number(in_at(in, 2), 0.0), band = opt_number(in_at(in, 3), kDefaultWeldTol);
  const auto k = static_cast<int>(axis);
  return select_by(in_at(in, 0), [&](const Value& v) {
    if (const Vec3* p = point_like(v)) return std::fabs((*p)[k] - level) <= band;
    if (const Line* l = v.as_line()) return std::fabs(l->a[k] - level) <= band && std::fabs(l->b[k] - level) <= band;
    return false;
  });
}

ComponentInfo SelectLengthNode::info() const {
  return {"keta.select_length", "Select by Length", kSets, "Select lines (or |numbers|) with length in [MIN, MAX]"};
}
std::vector<PinSpec> SelectLengthNode::inputs() const {
  return {PinSpec::tree("G", "any", "Lines or numbers").opt(), PinSpec::item("MIN", "number", "Minimum length in m (default 0)").opt(),
          PinSpec::item("MAX", "number", "Maximum length in m (unwired = no cap)").opt()};
}
std::vector<PinSpec> SelectLengthNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectLengthNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  const double mn = opt_number(in_at(in, 1), 0.0), mx = opt_number(in_at(in, 2), HUGE_VAL);
  return select_by(in_at(in, 0), [&](const Value& v) {
    if (const Line* l = v.as_line()) {
      double len = l->length();
      return len >= mn && len <= mx;
    }
    if (auto x = v.as_number()) return std::fabs(*x) >= mn && std::fabs(*x) <= mx;
    return false;
  });
}

ComponentInfo SelectOrientationNode::info() const {
  return {"keta.select_orientation", "Select by Orientation", kSets, "Select lines whose axis is within A degrees of +-D (undirected)"};
}
std::vector<PinSpec> SelectOrientationNode::inputs() const {
  return {PinSpec::tree("G", "any", "Lines").opt(), PinSpec::item("D", "vector", "Reference direction (default 0,0,1)").opt(),
          PinSpec::item("A", "number", "Angular tolerance in deg (default 1)").opt()};
}
std::vector<PinSpec> SelectOrientationNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectOrientationNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  const Vec3 d = first_point(in_at(in, 1)).value_or(Vec3(0, 0, 1));
  const double cos_tol = std::cos(opt_number(in_at(in, 2), 1.0) * kPi / 180.0);
  const double dn = d.norm();
  return select_by(in_at(in, 0), [&](const Value& v) {
    const Line* l = v.as_line();
    if (!l) return false;
    Vec3 axis = l->b - l->a;
    double an = axis.norm();
    if (an <= 0 || dn <= 0) return false;  // never a NaN comparison
    return std::fabs(axis.dot(d) / (an * dn)) >= cos_tol;
  });
}

ComponentInfo SelectProximityNode::info() const {
  return {"keta.select_proximity", "Select by Proximity", kSets, "Select points/lines within radius R of any target point"};
}
std::vector<PinSpec> SelectProximityNode::inputs() const {
  return {PinSpec::tree("G", "any", "Points or lines").opt(), PinSpec::list("P", "point", "Target points"),
          PinSpec::item("R", "number", "Radius in m (default 0.001)").opt()};
}
std::vector<PinSpec> SelectProximityNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectProximityNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  std::vector<Vec3> targets;
  for (const Value& v : in_at(in, 1).non_nulls())
    if (auto p = v.as_point()) targets.push_back(*p);
  constexpr unsigned long long kMaxPairs = 5'000'000;
  const unsigned long long pairs =
      static_cast<unsigned long long>(in_at(in, 0).item_count()) * static_cast<unsigned long long>(targets.size());
  if (pairs > kMaxPairs)
    throw ComponentError("Select by Proximity: " + std::to_string(pairs) + " item x target pairs exceed the cap of " +
                         std::to_string(kMaxPairs));
  const double r = opt_number(in_at(in, 2), kDefaultWeldTol);
  return select_by(in_at(in, 0), [&](const Value& v) {
    if (const Vec3* p = point_like(v))
      return std::any_of(targets.begin(), targets.end(), [&](const Vec3& t) { return (*p - t).norm() <= r; });
    if (const Line* l = v.as_line())
      return std::any_of(targets.begin(), targets.end(), [&](const Vec3& t) { return point_seg_dist(t, l->a, l->b) <= r; });
    return false;
  });
}

ComponentInfo SelectBySectionNode::info() const {
  return {"keta.select_by_section", "Select by Section", kSets, "Select beams whose section designation matches S (case-insensitive)"};
}
std::vector<PinSpec> SelectBySectionNode::inputs() const {
  return {PinSpec::tree("E", "beam", "Beams").opt(), PinSpec::item("S", "text", "Section designation to match (e.g. \"IPE300\")")};
}
std::vector<PinSpec> SelectBySectionNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectBySectionNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  const std::string want = opt_text(in_at(in, 1), "");
  if (want.empty()) throw ComponentError("Select by Section: S (section designation) is required");
  // Same key as the catalog resolver: "IPE 300" == "IPE300" == "ipe-300".
  const std::string norm = normalize_designation(want);
  return select_by(in_at(in, 0), [&](const Value& v) {
    const auto* b = v.get<BeamSpec>();
    return b && normalize_designation(b->section.s.name) == norm;
  });
}

ComponentInfo SelectByMaterialNode::info() const {
  return {"keta.select_by_material", "Select by Material", kSets, "Select beams whose material name matches S (case-insensitive)"};
}
std::vector<PinSpec> SelectByMaterialNode::inputs() const {
  return {PinSpec::tree("E", "beam", "Beams").opt(), PinSpec::item("S", "text", "Material name to match (e.g. \"Steel\")")};
}
std::vector<PinSpec> SelectByMaterialNode::outputs() const { return selector_outputs(); }
Outputs<Value> SelectByMaterialNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return empty_select();
  const std::string want = opt_text(in_at(in, 1), "");
  if (want.empty()) throw ComponentError("Select by Material: S (material name) is required");
  return select_by(in_at(in, 0), [&](const Value& v) {
    const auto* b = v.get<BeamSpec>();
    return b && !b->material.m.name.empty() && iequals(b->material.m.name, want);
  });
}

// ============================ group / sets ============================

ComponentInfo GroupByNode::info() const {
  return {"keta.group_by", "Group By", kSets, "Group each branch by key into child branches (keys sorted ascending)"};
}
std::vector<PinSpec> GroupByNode::inputs() const {
  return {PinSpec::tree("L", "any", "Items").opt(), PinSpec::tree("K", "any", "Keys (matched per branch, repeat-last)")};
}
std::vector<PinSpec> GroupByNode::outputs() const {
  return {PinSpec::tree("G", "any", "Grouped items, one child branch {p;g} per key"),
          PinSpec::tree("K", "any", "The group key of each child branch {p;g}")};
}
Outputs<Value> GroupByNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const BranchIndex keys_ix(in_at(in, 1));
  Tree grouped, keys_out;
  std::size_t bi = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    const Branch<Value>* kb = keys_ix.at(bi++);
    std::vector<std::pair<GroupKey, Branch<Value>>> groups;
    for (std::size_t j = 0; j < branch.size(); ++j) {
      if (!branch[j]) continue;
      if (!kb || kb->empty()) throw ComponentError("Group By: keys K are required to group items");
      const auto& ks = (*kb)[std::min(j, kb->size() - 1)];
      if (!ks) continue;
      GroupKey key = key_of(*ks);
      auto it = std::find_if(groups.begin(), groups.end(), [&](const auto& g) { return key_cmp(g.first, key) == 0; });
      if (it == groups.end())
        groups.push_back({key, {branch[j]}});
      else
        it->second.push_back(branch[j]);
    }
    std::stable_sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) { return key_cmp(a.first, b.first) < 0; });
    for (std::size_t g = 0; g < groups.size(); ++g) {
      Path gp = path.appended(static_cast<int>(g));
      grouped.ensure_path(gp) = groups[g].second;
      keys_out.append(gp, key_value(groups[g].first));
    }
  }
  return {grouped, keys_out};
}

ComponentInfo CreateSetNode::info() const {
  return {"keta.create_set", "Create Set", kSets, "Build a named set of flat indices (sorted + deduped)"};
}
std::vector<PinSpec> CreateSetNode::inputs() const {
  return {PinSpec::tree("I", "int", "Flat indices (branch-sorted positions)").opt(),
          PinSpec::item("N", "text", "Set name (default \"set\")").opt()};
}
std::vector<PinSpec> CreateSetNode::outputs() const { return {PinSpec::item("S", "set", "The named set")}; }
Outputs<Value> CreateSetNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  std::vector<std::int64_t> idx;
  for (const Value& v : in_at(in, 0).non_nulls()) {
    auto i = v.as_int();
    if (!i) throw ComponentError("Create Set: indices must be integers");
    if (*i < 0) throw ComponentError("Create Set: negative index");
    idx.push_back(*i);
  }
  return {single(Value::of(SetSpec::make(opt_text(in_at(in, 1), "set"), std::move(idx))))};
}

ComponentInfo SetOpNode::info() const {
  return {"keta.set_op", "Set Algebra", kSets, "Union / intersection / difference / symmetric difference of two sets"};
}
std::vector<PinSpec> SetOpNode::inputs() const { return {PinSpec::item("A", "set", "First set"), PinSpec::item("B", "set", "Second set")}; }
std::vector<PinSpec> SetOpNode::outputs() const { return {PinSpec::item("R", "set", "The combined set")}; }
Outputs<Value> SetOpNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {map_items({{"A", &in_at(in, 0)}, {"B", &in_at(in, 1)}}, [&](const std::vector<const Value*>& it) {
    const SetSpec& a = require_set(*it[0], "Set Algebra: A");
    const SetSpec& b = require_set(*it[1], "Set Algebra: B");
    std::vector<std::int64_t> r;
    const char* sym;
    if (op == "union") {
      std::set_union(a.indices.begin(), a.indices.end(), b.indices.begin(), b.indices.end(), std::back_inserter(r));
      sym = "∪";
    } else if (op == "intersection") {
      std::set_intersection(a.indices.begin(), a.indices.end(), b.indices.begin(), b.indices.end(), std::back_inserter(r));
      sym = "∩";
    } else if (op == "difference") {
      std::set_difference(a.indices.begin(), a.indices.end(), b.indices.begin(), b.indices.end(), std::back_inserter(r));
      sym = "∖";
    } else if (op == "symmetric difference") {
      std::set_symmetric_difference(a.indices.begin(), a.indices.end(), b.indices.begin(), b.indices.end(),
                                    std::back_inserter(r));
      sym = "△";
    } else {
      throw ComponentError("Set Algebra: unknown op '" + op + "'");
    }
    return Value::of(SetSpec::make(a.name + " " + sym + " " + b.name, std::move(r)));
  })};
}
json::Value SetOpNode::params() const {
  json::Value p = json::Value::object();
  p.set("op", op);
  return p;
}
void SetOpNode::set_params(const json::Value& p) {
  if (const auto* o = p.get("op"); o && o->is_string()) op = o->as_string();
}

ComponentInfo SetMembersNode::info() const {
  return {"keta.set_members", "Set Members", kSets, "Unpack a set into its indices, name and count"};
}
std::vector<PinSpec> SetMembersNode::inputs() const { return {PinSpec::item("S", "set", "The set").opt()}; }
std::vector<PinSpec> SetMembersNode::outputs() const {
  return {PinSpec::list("I", "int", "The flat indices (sorted)"), PinSpec::item("N", "text", "The set name"),
          PinSpec::item("C", "int", "The member count")};
}
Outputs<Value> SetMembersNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree idx, name, count;
  for (const auto& [path, branch] : in_at(in, 0).branches())
    for (const auto& slot : branch) {
      if (!slot) continue;
      const SetSpec& s = require_set(*slot, "Set Members: S");
      for (std::int64_t i : s.indices) idx.append(path, Value::integer(i));
      name.append(path, Value::text(s.name));
      count.append(path, Value::integer(static_cast<std::int64_t>(s.indices.size())));
    }
  return {idx, name, count};
}

ComponentInfo FilterBySetNode::info() const {
  return {"keta.filter_by_set", "Filter by Set", kSets, "Split a list into in-set (A) and not-in-set (B) by flat index"};
}
std::vector<PinSpec> FilterBySetNode::inputs() const {
  return {PinSpec::tree("L", "any", "Items (flat indexing over branch-sorted order)").opt(),
          PinSpec::item("S", "set", "The set of flat indices")};
}
std::vector<PinSpec> FilterBySetNode::outputs() const {
  return {PinSpec::tree("A", "any", "Items whose flat index is in the set"),
          PinSpec::tree("B", "any", "Items whose flat index is not in the set")};
}
Outputs<Value> FilterBySetNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (in_at(in, 0).empty()) return {Tree(), Tree()};
  const SetSpec* set = first_of<SetSpec>(in_at(in, 1));  // borrowed from the input tree
  if (!set) throw ComponentError("Filter by Set: S is required");
  Tree a, b;
  std::int64_t flat = 0;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    a.ensure_path(path);
    b.ensure_path(path);
    for (const auto& slot : branch) {
      const bool in_set = std::binary_search(set->indices.begin(), set->indices.end(), flat);
      if (slot) (in_set ? a : b).append(path, *slot);
      ++flat;
    }
  }
  return {a, b};
}

// ============================ Selection ============================

ComponentInfo SelectionNode::info() const {
  return {"keta.selection", "Selection", kSets,
          "A named set from a 3D viewport pick (nodes/bars -> flat-index Set). Wire E (the beam list that fed "
          "Assemble) so bar ids map to the right flat slots when the list has nulls; without E, index = id-1 "
          "(exact only for a null-free list). Node ids always map index = id-1."};
}
std::vector<PinSpec> SelectionNode::inputs() const {
  return {PinSpec::tree("E", "any", "Elements that fed Assemble (beams for bars)").opt()};
}
std::vector<PinSpec> SelectionNode::outputs() const {
  return {PinSpec::item("S", "set", "The named set (flat indices)"), PinSpec::list("I", "int", "The picked model ids"),
          PinSpec::item("N", "text", "The set name"), PinSpec::list("W", "text", "Warnings (empty when the pick maps cleanly)")};
}
Outputs<Value> SelectionNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  if (kind != "nodes" && kind != "bars")
    throw ComponentError("Selection: kind must be \"nodes\" or \"bars\" (got \"" + kind + "\")");
  const Tree& elements = in_at(in, 0);
  std::vector<std::int64_t> idx;
  Tree warn;
  if (kind == "bars" && !elements.empty()) {
    std::vector<std::int64_t> nonnull_flat;  // flat slot of the k-th non-null element
    std::int64_t flat = 0;
    for (const auto& kv : elements.branches())
      for (const auto& slot : kv.second) {
        if (slot) nonnull_flat.push_back(flat);
        ++flat;
      }
    std::size_t dropped = 0;
    for (std::int64_t id : ids) {
      if (id < 1) continue;
      if (static_cast<std::size_t>(id) <= nonnull_flat.size())
        idx.push_back(nonnull_flat[static_cast<std::size_t>(id - 1)]);
      else
        ++dropped;
    }
    if (nonnull_flat.empty() && !ids.empty())
      warn.append(Path::single(0), Value::text("E has no non-null elements — all " + std::to_string(ids.size()) +
                                                " picked ids dropped"));
    else if (dropped > 0)
      warn.append(Path::single(0), Value::text(std::to_string(dropped) + " of " + std::to_string(ids.size()) +
                                                " picked ids exceed the element list — dropped"));
  } else {
    for (std::int64_t id : ids)
      if (id >= 1) idx.push_back(id - 1);
  }
  std::vector<std::int64_t> sorted = ids;
  std::sort(sorted.begin(), sorted.end());
  sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
  Tree id_out;
  for (std::int64_t i : sorted) id_out.append(Path::single(0), Value::integer(i));
  return {single(Value::of(SetSpec::make(name, std::move(idx)))), id_out, single(Value::text(name)), warn};
}
json::Value SelectionNode::params() const {
  json::Value arr = json::Value::array();
  for (std::int64_t i : ids) arr.push_back(json::Value(i));
  json::Value p = json::Value::object();
  p.set("kind", kind);
  p.set("ids", std::move(arr));
  p.set("name", name);
  return p;
}
void SelectionNode::set_params(const json::Value& p) {
  if (const auto* k = p.get("kind"); k && k->is_string()) kind = k->as_string();
  if (const auto* a = p.get("ids"); a && a->is_array()) {
    ids.clear();
    for (const auto& v : a->as_array())
      if (v.is_number() && std::trunc(v.as_number()) == v.as_number()) ids.push_back(static_cast<std::int64_t>(v.as_number()));
  }
  if (const auto* n = p.get("name"); n && n->is_string()) name = n->as_string();
}

// ============================ Catalog Section ============================

namespace {
Value catalog_value(const std::string& designation) {
  auto s = catalog_section(designation);
  if (!s) throw ComponentError("unknown catalog section '" + designation + "'");
  return Value::of(SectionSpec{*s});
}
}  // namespace

ComponentInfo SectionCatalogNode::info() const {
  return {"keta.section_catalog", "Catalog Section", "Model", "European steel profile by designation"};
}
std::vector<PinSpec> SectionCatalogNode::inputs() const {
  return {PinSpec::item("D", "text",
                        "Per-item designations (override the param when wired; a null item falls back to the param "
                        "designation)")
              .opt()};
}
std::vector<PinSpec> SectionCatalogNode::outputs() const { return {PinSpec::item("S", "section", "The section")}; }
Outputs<Value> SectionCatalogNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const Tree& d = in_at(in, 0);
  if (!has_values(d)) return {single(catalog_value(designation))};
  // Wired: one section per D SLOT, nulls included (never dropped — that would
  // misalign every following section against the lines Beam pairs them with).
  Tree out;
  for (const auto& [path, slots] : d.branches()) {
    auto& dst = out.ensure_path(path);
    for (const auto& slot : slots) {
      if (!slot) {
        dst.emplace_back(catalog_value(designation));
        continue;
      }
      const std::string* t = slot->as_text();
      if (!t) throw ComponentError("D must be designation text, got " + slot->type_tag());
      dst.emplace_back(catalog_value(*t));
    }
  }
  return {out};
}
json::Value SectionCatalogNode::params() const {
  json::Value p = json::Value::object();
  p.set("designation", designation);
  return p;
}
void SectionCatalogNode::set_params(const json::Value& p) {
  if (const auto* d = p.get("designation"); d && d->is_string()) designation = d->as_string();
}

}  // namespace keta::param

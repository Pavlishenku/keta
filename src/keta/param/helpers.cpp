#include "keta/param/helpers.h"

#include <algorithm>
#include <cmath>

namespace keta::param {

Tree map_items(const std::vector<std::pair<std::string, const Tree*>>& inputs,
               const std::function<Value(const std::vector<const Value*>&)>& f) {
  for (const auto& [pin, t] : inputs)
    if (!has_values(*t)) throw ComponentError::missing(pin);

  using View = std::pair<const Path*, const Branch<Value>*>;
  std::vector<std::vector<View>> views;
  views.reserve(inputs.size());
  for (const auto& in : inputs) {
    std::vector<View> v;
    for (const auto& kv : in.second->branches()) v.emplace_back(&kv.first, &kv.second);
    views.push_back(std::move(v));
  }
  std::size_t n_branches = 0, principal = 0;
  for (std::size_t i = 0; i < views.size(); ++i) {
    if (views[i].size() > n_branches) {  // strict: ties keep the earliest input
      n_branches = views[i].size();
      principal = i;
    }
  }

  Tree out;
  std::vector<const Value*> items(inputs.size());
  for (std::size_t bi = 0; bi < n_branches; ++bi) {
    std::vector<const Branch<Value>*> brs;
    brs.reserve(views.size());
    std::size_t n_items = 0;
    for (const auto& v : views) {
      const Branch<Value>* b = v[std::min(bi, v.size() - 1)].second;
      brs.push_back(b);
      n_items = std::max(n_items, b->size());
    }
    const Path& path = *views[principal][std::min(bi, views[principal].size() - 1)].first;
    for (std::size_t ii = 0; ii < n_items; ++ii) {
      bool any_null = false;
      for (std::size_t k = 0; k < brs.size(); ++k) {
        const Branch<Value>& b = *brs[k];
        if (b.empty() || !b[std::min(ii, b.size() - 1)]) {
          any_null = true;
          break;
        }
        items[k] = &*b[std::min(ii, b.size() - 1)];
      }
      if (any_null) continue;
      out.append(path, f(items));
    }
  }
  return out;
}

Tree number_or(const Tree& t, double dflt) { return has_values(t) ? t : single(Value::number(dflt)); }

std::optional<double> first_number(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (auto x = o->as_number()) return x;
  return std::nullopt;
}
std::optional<std::int64_t> first_int(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (auto x = o->as_int()) return x;
  return std::nullopt;
}
std::optional<Vec3> first_vector(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (auto x = o->as_vector()) return x;
  return std::nullopt;
}
std::optional<std::string> first_text(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (const auto* s = o->as_text()) return *s;
  return std::nullopt;
}

std::optional<bool> first_bool(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (auto x = o->as_bool()) return x;
  return std::nullopt;
}
std::optional<Vec3> first_point(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (auto x = o->as_point()) return x;
  return std::nullopt;
}

std::vector<bool> pattern_of(const Branch<Value>* b) {
  std::vector<bool> p;
  if (b)
    for (const auto& o : *b) p.push_back(o && o->as_bool().value_or(false));
  if (p.empty()) p = {true, false};
  return p;
}

std::vector<Line> beam_lines(const Tree& t, bool accept_lines) {
  std::vector<Line> out;
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second) {
      if (!o) continue;
      if (const auto* b = o->get<BeamSpec>())
        out.push_back(b->line);
      else if (accept_lines)
        if (const auto* l = o->as_line()) out.push_back(*l);
    }
  return out;
}

void wrong_type(const std::string& pin, const std::string& expected, const Value& got) {
  throw ComponentError(pin + " must be a " + expected + ", got " + got.type_tag());
}

// ---------------------------------------------------------------------------
// NodeMerger
// ---------------------------------------------------------------------------

NodeMerger::Cell NodeMerger::cell_of(const Vec3& p) const {
  auto cell = [&](double v) -> std::int64_t {
    double c = std::floor(v / tol_);
    return std::isfinite(c) && std::fabs(c) < 9.0e18 ? static_cast<std::int64_t>(c) : 0;
  };
  return {cell(p.x()), cell(p.y()), cell(p.z())};
}

std::vector<std::uint32_t> NodeMerger::candidates(const Vec3& p) const {
  auto [cx, cy, cz] = cell_of(p);
  std::vector<std::uint32_t> c;
  for (int dx = -1; dx <= 1; ++dx)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dz = -1; dz <= 1; ++dz) {
        auto it = grid_.find({cx + dx, cy + dy, cz + dz});
        if (it != grid_.end()) c.insert(c.end(), it->second.begin(), it->second.end());
      }
  std::sort(c.begin(), c.end());
  return c;
}

Id NodeMerger::id_for(const Vec3& p) {
  if (!(tol_ > 0.0)) {  // exact linear scan fallback
    for (std::size_t i = 0; i < pos_.size(); ++i)
      if ((pos_[i] - p).norm() <= tol_) return static_cast<Id>(i + 1);
    pos_.push_back(p);
    return static_cast<Id>(pos_.size());
  }
  for (std::uint32_t idx : candidates(p))  // ascending: first match = smallest index
    if ((pos_[idx] - p).norm() <= tol_) return static_cast<Id>(idx) + 1;
  auto idx = static_cast<std::uint32_t>(pos_.size());
  pos_.push_back(p);
  grid_[cell_of(p)].push_back(idx);
  return static_cast<Id>(idx) + 1;
}

std::optional<Id> NodeMerger::find(const Vec3& p) const {
  std::optional<std::pair<std::uint32_t, double>> best;
  auto consider = [&](std::uint32_t idx) {
    double d = (pos_[idx] - p).norm();
    if (d <= tol_ && (!best || d < best->second)) best = std::make_pair(idx, d);
  };
  if (!(tol_ > 0.0)) {
    for (std::uint32_t i = 0; i < pos_.size(); ++i) consider(i);
  } else {
    for (std::uint32_t idx : candidates(p)) consider(idx);
  }
  if (!best) return std::nullopt;
  return static_cast<Id>(best->first) + 1;
}

std::optional<BarMatch> find_bar(const Model& m, const Line& l, double tol) {
  for (std::size_t i = 0; i < m.bars.size(); ++i) {
    const Vec3& p1 = m.nodes[static_cast<std::size_t>(m.bars[i].ni - 1)].x;
    const Vec3& p2 = m.nodes[static_cast<std::size_t>(m.bars[i].nj - 1)].x;
    if ((p1 - l.a).norm() <= tol && (p2 - l.b).norm() <= tol) return BarMatch{i, false};
    if ((p1 - l.b).norm() <= tol && (p2 - l.a).norm() <= tol) return BarMatch{i, true};
  }
  return std::nullopt;
}

}  // namespace keta::param

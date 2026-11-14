// Snapshot change describer (see diff.h).
#include "keta/param/diff.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <map>
#include <optional>
#include <set>
#include <tuple>

#include "keta/json.h"

namespace keta::param {

namespace {

constexpr std::size_t kCap = 6;
const std::string kMinus = "−";

struct NodeInfo {
  std::uint64_t key = 0;
  std::string type_id;
  json::Value params, pos;
};

// Structural equality; object members compare as maps.
bool same(const json::Value& a, const json::Value& b) {
  if (a.type() != b.type()) return false;
  switch (a.type()) {
    case json::Value::Type::Null: return true;
    case json::Value::Type::Bool: return a.as_bool() == b.as_bool();
    case json::Value::Type::Number: return a.as_number() == b.as_number();
    case json::Value::Type::String: return a.as_string() == b.as_string();
    case json::Value::Type::Array: {
      const auto &x = a.as_array(), &y = b.as_array();
      if (x.size() != y.size()) return false;
      for (std::size_t i = 0; i < x.size(); ++i)
        if (!same(x[i], y[i])) return false;
      return true;
    }
    case json::Value::Type::Object: {
      if (a.members().size() != b.members().size()) return false;
      for (const auto& [k, v] : a.members()) {
        const json::Value* w = b.get(k);
        if (!w || !same(v, *w)) return false;
      }
      return true;
    }
  }
  return false;
}

std::optional<std::uint64_t> as_u64(const json::Value* v) {
  if (!v || !v->is_number()) return std::nullopt;
  const double d = v->as_number();
  if (!(d >= 0) || std::trunc(d) != d || d >= 1.8e19) return std::nullopt;
  return static_cast<std::uint64_t>(d);
}

std::optional<std::vector<NodeInfo>> extract_nodes(const json::Value& doc) {
  const json::Value* arr = doc.is_object() ? doc.get("nodes") : nullptr;
  if (!arr || !arr->is_array()) return std::nullopt;
  std::vector<NodeInfo> out;
  std::size_t i = 0;
  for (const json::Value& n : arr->as_array()) {
    NodeInfo info;
    const bool obj = n.is_object();
    info.key = as_u64(obj ? n.get("key") : nullptr).value_or(i);
    if (const json::Value* t = obj ? n.get("type_id") : nullptr; t && t->is_string()) info.type_id = t->as_string();
    if (const json::Value* p = obj ? n.get("params") : nullptr) info.params = *p;
    if (const json::Value* p = obj ? n.get("pos") : nullptr) info.pos = *p;
    out.push_back(std::move(info));
    ++i;
  }
  return out;
}

using End = std::pair<std::uint64_t, std::uint64_t>;
std::optional<End> endpoint(const json::Value* v) {
  if (!v || !v->is_array() || v->as_array().size() < 2) return std::nullopt;
  auto a = as_u64(&v->as_array()[0]), b = as_u64(&v->as_array()[1]);
  if (!a || !b) return std::nullopt;
  return End{*a, *b};
}

std::set<std::pair<End, End>> extract_wires(const json::Value& doc) {
  std::set<std::pair<End, End>> out;
  const json::Value* arr = doc.get("wires");
  if (!arr || !arr->is_array()) return out;
  for (const json::Value& w : arr->as_array()) {
    if (!w.is_object()) continue;
    auto f = endpoint(w.get("from")), t = endpoint(w.get("to"));
    if (f && t) out.insert({*f, *t});
  }
  return out;
}

const json::Array& top_array(const json::Value& doc, const char* key) {
  static const json::Array empty;
  const json::Value* v = doc.get(key);
  return v && v->is_array() ? v->as_array() : empty;
}

std::string plural(std::size_t n, const std::string& w) { return n == 1 ? w : w + "s"; }

std::string short_name(const std::string& type_id) {
  const auto dot = type_id.rfind('.');
  std::string seg = dot == std::string::npos ? type_id : type_id.substr(dot + 1);
  std::replace(seg.begin(), seg.end(), '_', ' ');
  if (!seg.empty() && seg[0] >= 'a' && seg[0] <= 'z') seg[0] = static_cast<char>(seg[0] - 'a' + 'A');
  return seg;
}

bool arrays_equal(const json::Array& a, const json::Array& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (!same(a[i], b[i])) return false;
  return true;
}

std::optional<std::string> str_field(const json::Value& v, const char* k) {
  const json::Value* x = v.is_object() ? v.get(k) : nullptr;
  if (x && x->is_string()) return x->as_string();
  return std::nullopt;
}

const json::Value& field(const json::Value& v, const char* k) {
  static const json::Value null;
  const json::Value* x = v.is_object() ? v.get(k) : nullptr;
  return x ? *x : null;
}

std::vector<std::string> group_entries(const json::Value& od, const json::Value& nd) {
  const json::Array &o = top_array(od, "groups"), &n = top_array(nd, "groups");
  std::vector<std::string> out;
  if (n.size() > o.size()) {
    const std::size_t added = n.size() - o.size();
    if (added == 1) {
      const auto name = str_field(n.back(), "name");
      out.push_back(name && !name->empty() ? "+1 group '" + *name + "'" : "+1 group");
    } else {
      out.push_back("+" + std::to_string(added) + " " + plural(added, "group"));
    }
  } else if (n.size() < o.size()) {
    const std::size_t removed = o.size() - n.size();
    out.push_back(kMinus + std::to_string(removed) + " " + plural(removed, "group"));
  } else if (!arrays_equal(o, n)) {
    bool renamed = false;
    for (std::size_t i = 0; i < o.size(); ++i) renamed |= str_field(o[i], "name") != str_field(n[i], "name");
    out.push_back(renamed ? "group renamed" : "group edited");
  }
  return out;
}

std::vector<std::string> note_entries(const json::Value& od, const json::Value& nd) {
  const json::Array &o = top_array(od, "annotations"), &n = top_array(nd, "annotations");
  std::vector<std::string> out;
  if (n.size() > o.size()) {
    const std::size_t added = n.size() - o.size();
    out.push_back("+" + std::to_string(added) + " " + plural(added, "note"));
  } else if (n.size() < o.size()) {
    const std::size_t removed = o.size() - n.size();
    out.push_back(kMinus + std::to_string(removed) + " " + plural(removed, "note"));
  } else if (!arrays_equal(o, n)) {
    bool text_changed = false;
    std::size_t moved = 0;
    for (std::size_t i = 0; i < o.size(); ++i) {
      text_changed |= !same(field(o[i], "text"), field(n[i], "text"));
      moved += !same(field(o[i], "pos"), field(n[i], "pos"));
    }
    if (text_changed)
      out.push_back("note edited");
    else if (moved > 0)
      out.push_back("moved " + std::to_string(moved) + " " + plural(moved, "note"));
  }
  return out;
}

bool is_primitive(const json::Value& v) {
  return v.is_null() || v.is_bool() || v.is_number() || v.is_string();
}

// Shortest round-trip decimal form, never an exponent.
std::string render_number(double f) {
  if (std::isfinite(f) && std::trunc(f) == f && std::fabs(f) < 9.2e18) return std::to_string(static_cast<long long>(f));
  char buf[512];
  auto r = std::to_chars(buf, buf + sizeof buf, f, std::chars_format::fixed);
  return std::string(buf, r.ptr);
}

std::string render_scalar(const json::Value& v) {
  if (v.is_null()) return "null";
  if (v.is_bool()) return v.as_bool() ? "true" : "false";
  if (v.is_string()) return "\"" + v.as_string() + "\"";
  if (v.is_number()) return render_number(v.as_number());
  return v.dump();
}

std::string param_entry(const std::string& name, const json::Value& o, const json::Value& n) {
  if (o.is_object() && n.is_object()) {
    std::set<std::string> keys;
    for (const auto& [k, v] : o.members()) keys.insert(k);
    for (const auto& [k, v] : n.members()) keys.insert(k);
    std::vector<std::pair<const json::Value*, const json::Value*>> diffs;
    for (const std::string& k : keys) {
      const json::Value& ov = field(o, k.c_str());
      const json::Value& nv = field(n, k.c_str());
      if (!same(ov, nv)) diffs.emplace_back(&ov, &nv);
    }
    if (diffs.size() == 1 && is_primitive(*diffs[0].first) && is_primitive(*diffs[0].second))
      return name + ": " + render_scalar(*diffs[0].first) + " → " + render_scalar(*diffs[0].second);
  }
  return name + ": params changed";
}

}  // namespace

std::vector<std::string> describe_change(const std::string& old_json, const std::string& new_json) {
  if (old_json == new_json) return {};
  json::Value od, nd;
  try {
    od = json::Value::parse(old_json);
    nd = json::Value::parse(new_json);
  } catch (...) {
    return {"edit"};
  }
  const auto on = extract_nodes(od), nn = extract_nodes(nd);
  if (!on || !nn) return {"edit"};

  std::vector<std::string> entries;
  std::map<std::string, std::size_t> oc, nc;
  for (const NodeInfo& n : *on) ++oc[n.type_id];
  for (const NodeInfo& n : *nn) ++nc[n.type_id];
  const bool types_identical = oc == nc;
  {
    std::set<std::string> all;
    for (const auto& [t, c] : oc) all.insert(t);
    for (const auto& [t, c] : nc) all.insert(t);
    std::vector<std::pair<std::string, std::size_t>> adds, removes;
    for (const std::string& t : all) {
      const std::size_t o = oc.count(t) ? oc[t] : 0, n = nc.count(t) ? nc[t] : 0;
      if (n > o) adds.emplace_back(short_name(t), n - o);
      if (n < o) removes.emplace_back(short_name(t), o - n);
    }
    std::sort(adds.begin(), adds.end());
    std::sort(removes.begin(), removes.end());
    for (const auto& [name, n] : adds) entries.push_back("+" + std::to_string(n) + " " + name);
    for (const auto& [name, n] : removes) entries.push_back(kMinus + std::to_string(n) + " " + name);
  }
  for (auto& e : group_entries(od, nd)) entries.push_back(std::move(e));
  for (auto& e : note_entries(od, nd)) entries.push_back(std::move(e));

  if (types_identical) {
    std::map<std::uint64_t, const NodeInfo*> old_by_key;
    for (const NodeInfo& n : *on) old_by_key[n.key] = &n;  // last wins, like collect() into a HashMap
    std::size_t moved = 0;
    for (const NodeInfo& n : *nn) {
      auto it = old_by_key.find(n.key);
      if (it == old_by_key.end() || it->second->type_id != n.type_id) continue;
      if (!same(it->second->params, n.params))
        entries.push_back(param_entry(short_name(n.type_id), it->second->params, n.params));
      else if (!same(it->second->pos, n.pos))
        ++moved;
    }
    auto key_map = [](const std::vector<NodeInfo>& v) {
      std::map<std::uint64_t, std::string> m;
      for (const NodeInfo& n : v) m[n.key] = n.type_id;  // last wins, like collect() into a map
      return m;
    };
    if (key_map(*on) == key_map(*nn)) {
      const auto ow = extract_wires(od), nw = extract_wires(nd);
      std::size_t added = 0, removed = 0;
      for (const auto& w : nw) added += !ow.count(w);
      for (const auto& w : ow) removed += !nw.count(w);
      if (added > 0) entries.push_back("+" + std::to_string(added) + " " + plural(added, "wire"));
      if (removed > 0) entries.push_back(kMinus + std::to_string(removed) + " " + plural(removed, "wire"));
    }
    if (moved > 0) entries.push_back("moved " + std::to_string(moved) + " " + plural(moved, "node"));
  }
  if (entries.empty()) return {"edit"};
  if (entries.size() > kCap) {
    const std::size_t extra = entries.size() - (kCap - 1);
    entries.resize(kCap - 1);
    entries.push_back("… and " + std::to_string(extra) + " more");
  }
  return entries;
}

}  // namespace keta::param

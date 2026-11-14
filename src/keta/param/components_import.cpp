// Weld and the CSV Points / CSV Bars importers.
//
// CSV reading (dependency-free, RFC 4180): records end at LF, CRLF or a lone
// CR; a field that starts with `"` is quoted (`""` is a literal quote;
// delimiters and line breaks inside are data); text after a closing quote is
// appended to the field; completely blank lines are skipped (so the warning
// line numbers match the spreadsheet only without blank interior lines); a
// UTF-8 BOM on the first record is stripped; rows may have any width.
//
// Units: the core is SI; CSV Points applies the one `unit` scale at the parse
// edge. Flat-index contract: one item per data row, in row order; a row that
// cannot be parsed is a visible null mirrored by a `W` warning — never
// dropped, never zeroed.
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <array>
#include <limits>
#include <map>
#include <set>
#include <sstream>

#include "keta/param/components.h"

namespace keta::param {

namespace {

const Tree& in_at(const std::vector<Tree>& in, std::size_t i) { return input_or_empty(in, i); }

using Record = std::vector<std::string>;

// Read a whole file, or a typed error naming it.
std::string read_file(const std::string& path, const char* who) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw ComponentError(std::string(who) + ": cannot read '" + path + "': " + std::strerror(errno));
  std::ostringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

std::vector<Record> parse_csv(const std::string& text, char delim) {
  std::vector<Record> rows;
  Record rec;
  std::string field;
  bool in_quotes = false, quoted = false, any = false;
  auto end_field = [&] {
    rec.push_back(std::move(field));
    field.clear();
    quoted = false;
  };
  auto end_record = [&] {
    if (any) {
      end_field();
      rows.push_back(std::move(rec));
    }
    rec.clear();
    field.clear();
    quoted = false;
    any = false;
  };
  for (std::size_t i = 0; i < text.size(); ++i) {
    const char c = text[i];
    if (in_quotes) {
      if (c == '"') {
        if (i + 1 < text.size() && text[i + 1] == '"') {
          field.push_back('"');
          ++i;
        } else {
          in_quotes = false;
        }
      } else {
        field.push_back(c);
      }
      continue;
    }
    if (c == '"' && field.empty() && !quoted) {
      in_quotes = quoted = any = true;
    } else if (c == delim) {
      any = true;
      end_field();
    } else if (c == '\n' || c == '\r') {
      if (c == '\r' && i + 1 < text.size() && text[i + 1] == '\n') ++i;
      end_record();
    } else {
      any = true;
      field.push_back(c);
    }
  }
  end_record();
  if (!rows.empty() && !rows[0].empty() && rows[0][0].rfind("\xEF\xBB\xBF", 0) == 0) rows[0][0].erase(0, 3);
  return rows;
}

// "auto" sniffs the first non-empty line: most occurrences outside quotes,
// ties toward comma; anything unreadable falls back to comma.
char resolve_delimiter(const std::string& text, const std::string& param) {
  if (param == ",") return ',';
  if (param == ";") return ';';
  if (param == "\t") return '\t';
  std::size_t start = 0;
  std::string line;
  while (start < text.size()) {
    std::size_t e = text.find_first_of("\r\n", start);
    line = text.substr(start, e == std::string::npos ? std::string::npos : e - start);
    if (line.find_first_not_of(" \t") != std::string::npos) break;
    line.clear();
    if (e == std::string::npos) break;
    start = e + 1;
  }
  std::size_t counts[3] = {0, 0, 0};
  bool q = false;
  for (char c : line) {
    if (c == '"') q = !q;
    else if (!q && c == ',') ++counts[0];
    else if (!q && c == ';') ++counts[1];
    else if (!q && c == '\t') ++counts[2];
  }
  const char cand[3] = {',', ';', '\t'};
  char best = ',';
  std::size_t best_n = 0;
  for (int k = 0; k < 3; ++k)
    if (counts[k] > best_n) {
      best_n = counts[k];
      best = cand[k];
    }
  return best;
}

double scale_for_unit(const std::string& unit) {
  if (unit == "cm") return 0.01;
  if (unit == "mm") return 0.001;
  if (unit == "in") return 0.0254;
  if (unit == "ft") return 0.3048;
  return 1.0;  // "m" and anything unrecognised
}

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n\v\f");
  if (a == std::string::npos) return {};
  const auto b = s.find_last_not_of(" \t\r\n\v\f");
  return s.substr(a, b - a + 1);
}

// A trimmed cell as a finite number: empty / non-numeric / nan / inf -> none.
std::optional<double> parse_num(const std::string& cell) {
  std::string t = trim(cell);
  if (t.empty()) return std::nullopt;
  if (t[0] == '+' && t.size() > 1 && t[1] != '-' && t[1] != '+') t.erase(0, 1);
  double v = 0;
  const auto r = std::from_chars(t.data(), t.data() + t.size(), v);
  if (r.ec != std::errc() || r.ptr != t.data() + t.size() || !std::isfinite(v)) return std::nullopt;
  return v;
}

const std::string* cell(const Record& r, std::size_t col) { return col < r.size() ? &r[col] : nullptr; }

// Ids print without a decimal when integral; otherwise the shortest form.
std::string fmt_id(double x) {
  if (std::trunc(x) == x && std::fabs(x) < 9.007e15) return std::to_string(static_cast<long long>(x));
  char buf[32];
  auto r = std::to_chars(buf, buf + sizeof buf, x);
  return std::string(buf, r.ptr);
}

std::uint64_t id_key(double x) {
  if (x == 0.0) x = 0.0;  // fold -0
  std::uint64_t k;
  std::memcpy(&k, &x, sizeof k);
  return k;
}

json::Value col_json(const std::optional<std::size_t>& c) {
  return c ? json::Value(static_cast<double>(*c)) : json::Value();
}
std::optional<std::size_t> col_of(const json::Value& map, const char* key) {
  const auto* x = map.get(key);
  if (!x || !x->is_number()) return std::nullopt;
  const double v = x->as_number();
  if (v < 0 || v != std::trunc(v)) return std::nullopt;
  return static_cast<std::size_t>(v);
}

void read_common(const json::Value& p, std::string& path, std::string& delimiter, bool& has_header) {
  if (const auto* x = p.get("path"); x && x->is_string()) path = x->as_string();
  if (const auto* x = p.get("delimiter"); x && x->is_string()) delimiter = x->as_string();
  if (const auto* x = p.get("has_header"); x && x->is_bool()) has_header = x->as_bool();
}

}  // namespace

// ============================ Weld ============================

ComponentInfo WeldNode::info() const {
  return {"keta.weld", "Weld", "Topology", "Weld coincident line endpoints into a shared node cloud (remap out)"};
}
std::vector<PinSpec> WeldNode::inputs() const {
  return {PinSpec::tree("L", "line", "Lines to weld").opt(), PinSpec::item("T", "number", "Tolerance (m) (default 0.001)").opt()};
}
std::vector<PinSpec> WeldNode::outputs() const {
  return {PinSpec::tree("L", "line",
                        "Welded lines: same paths, endpoints snapped to their node; a line whose two ends weld together "
                        "becomes a null slot (collapsed lines become nulls, nothing silently dropped)"),
          PinSpec::list("N", "point",
                        "Unique node cloud, first-appearance order (0-based here; Assemble assigns the same order "
                        "1-based at the same tolerance)"),
          PinSpec::tree("A", "int", "Node index (into N) of each line's start, aligned with input L"),
          PinSpec::tree("B", "int", "Node index (into N) of each line's end, aligned with input L")};
}
Outputs<Value> WeldNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  // The SAME welder Assemble uses, fed in the same deterministic order
  // (branches sorted, items in order, a then b): what Weld merges, Assemble merges.
  NodeMerger merger(first_number(in_at(in, 1)).value_or(1e-3));
  Tree l_out, a_out, b_out;
  for (const auto& [path, branch] : in_at(in, 0).branches()) {
    l_out.ensure_path(path);
    a_out.ensure_path(path);
    b_out.ensure_path(path);
    for (const auto& slot : branch) {
      if (!slot) {
        l_out.append_null(path);
        a_out.append_null(path);
        b_out.append_null(path);
        continue;
      }
      const Line* l = slot->as_line();
      if (!l) throw ComponentError("Weld: expected a line, got " + slot->type_tag());
      const Id ia = merger.id_for(l->a) - 1, ib = merger.id_for(l->b) - 1;
      a_out.append(path, Value::integer(ia));
      b_out.append(path, Value::integer(ib));
      if (ia == ib) {
        l_out.append_null(path);  // both ends on one node: the line collapsed
      } else {
        const auto& pos = merger.positions();
        l_out.append(path, Value::of(Line{pos[static_cast<std::size_t>(ia)], pos[static_cast<std::size_t>(ib)]}));
      }
    }
  }
  Tree n_out;
  n_out.ensure_path(Path::single(0));
  for (const Vec3& p : merger.positions()) n_out.append(Path::single(0), Value::of(Point{p}));
  return {l_out, n_out, a_out, b_out};
}

// ============================ CSV Points ============================

ComponentInfo CsvPointsNode::info() const {
  return {"keta.csv_points", "CSV Points", "Model", "Import points from a CSV: column mapping, units, row order preserved"};
}
std::vector<PinSpec> CsvPointsNode::inputs() const { return {}; }
std::vector<PinSpec> CsvPointsNode::outputs() const {
  return {PinSpec::list("P", "point", "Points, in row order (null on a bad row)"),
          PinSpec::list("I", "number", "Ids (mapped id column, else 1-based row number)"),
          PinSpec::list("W", "text", "One warning per import issue")};
}
Outputs<Value> CsvPointsNode::solve(const SolveCtx&, const std::vector<Tree>&) const {
  if (trim(path).empty()) throw ComponentError("CSV Points: no file selected");
  if (!map_x && !map_y && !map_z)
    throw ComponentError("CSV Points: no coordinate column mapped (map at least one of X/Y/Z)");
  const std::string text = read_file(path, "CSV Points");
  const double scale = scale_for_unit(unit);
  const std::vector<Record> rows = parse_csv(text, resolve_delimiter(text, delimiter));
  const std::size_t hdr = has_header ? 1 : 0;
  if (rows.size() <= hdr) throw ComponentError("CSV Points: no data rows in '" + path + "'");

  const Path b0 = Path::single(0);
  Tree points, ids;
  points.ensure_path(b0);
  ids.ensure_path(b0);
  std::vector<std::string> warnings;
  std::set<std::uint64_t> seen;
  std::vector<double> dups;
  auto note_id = [&](double v) {
    if (!seen.insert(id_key(v)).second &&
        std::none_of(dups.begin(), dups.end(), [&](double d) { return id_key(d) == id_key(v); }))
      dups.push_back(v);
  };
  for (std::size_t i = hdr; i < rows.size(); ++i) {
    const Record& rec = rows[i];
    const std::size_t rown = i - hdr + 1, line = rown + hdr;
    double coord[3] = {0, 0, 0};
    bool bad = false;
    const std::pair<const char*, std::optional<std::size_t>> axes[3] = {{"x", map_x}, {"y", map_y}, {"z", map_z}};
    for (int a = 0; a < 3; ++a) {
      if (!axes[a].second) continue;
      const std::size_t col = *axes[a].second;
      if (const std::string* c = cell(rec, col)) {
        if (auto v = parse_num(*c)) {
          coord[a] = *v * scale;
        } else {
          warnings.push_back("row " + std::to_string(line) + ": " + axes[a].first + " '" + trim(*c) + "' is not a number");
          bad = true;
        }
      } else {
        warnings.push_back("row " + std::to_string(line) + ": " + axes[a].first + " column " + std::to_string(col) + " missing");
        bad = true;
      }
    }
    if (bad)
      points.append_null(b0);
    else
      points.append(b0, Value::of(Point{Vec3(coord[0], coord[1], coord[2])}));

    double id = static_cast<double>(rown);
    if (map_id) {
      const std::string* c = cell(rec, *map_id);
      if (auto v = c ? parse_num(*c) : std::nullopt) {
        id = *v;
      } else {
        warnings.push_back("row " + std::to_string(line) + ": id '" + (c ? trim(*c) : std::string()) + "' is not a number");
      }
      note_id(id);
    }
    ids.append(b0, Value::number(id));
  }
  if (!dups.empty()) {
    std::sort(dups.begin(), dups.end());
    std::string list;
    for (std::size_t k = 0; k < dups.size(); ++k) list += (k ? ", " : "") + fmt_id(dups[k]);
    warnings.push_back("duplicate ids: " + list);
  }
  Tree w;
  for (auto& s : warnings) w.append(b0, Value::text(std::move(s)));
  return {points, ids, w};
}
json::Value CsvPointsNode::params() const {
  json::Value map = json::Value::object();
  map.set("id", col_json(map_id));
  map.set("x", col_json(map_x));
  map.set("y", col_json(map_y));
  map.set("z", col_json(map_z));
  json::Value p = json::Value::object();
  p.set("path", path);
  p.set("delimiter", delimiter);
  p.set("has_header", has_header);
  p.set("unit", unit);
  p.set("map", std::move(map));
  return p;
}
void CsvPointsNode::set_params(const json::Value& p) {
  read_common(p, path, delimiter, has_header);
  if (const auto* x = p.get("unit"); x && x->is_string()) unit = x->as_string();
  if (const auto* m = p.get("map")) {
    map_id = col_of(*m, "id");
    map_x = col_of(*m, "x");
    map_y = col_of(*m, "y");
    map_z = col_of(*m, "z");
  }
}

// ============================ CSV Bars ============================

ComponentInfo CsvBarsNode::info() const {
  return {"keta.csv_bars", "CSV Bars", "Model", "Import bars from a CSV: id→point lookup, section/material tags"};
}
std::vector<PinSpec> CsvBarsNode::inputs() const {
  return {PinSpec::list("P", "point", "Points to reference (from CSV Points)"),
          PinSpec::list("I", "number", "Ids parallel to P (from CSV Points)")};
}
std::vector<PinSpec> CsvBarsNode::outputs() const {
  return {PinSpec::list("L", "line", "Bar axis lines, row order (null on a bad ref)"),
          PinSpec::list("I", "number", "Bar ids (mapped column, else row number)"),
          PinSpec::list("S", "text", "Section tags (null where unmapped/blank)"),
          PinSpec::list("M", "text", "Material tags (null where unmapped/blank)"),
          PinSpec::list("W", "text", "One warning per import issue")};
}
Outputs<Value> CsvBarsNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  // P and I flattened with their nulls: positional, flat-index aligned.
  std::vector<std::optional<Vec3>> pts;
  std::vector<std::optional<double>> ids;
  for (const auto& kv : in_at(in, 0).branches())
    for (const auto& o : kv.second) pts.push_back(o ? o->as_point() : std::nullopt);
  for (const auto& kv : in_at(in, 1).branches())
    for (const auto& o : kv.second) ids.push_back(o ? o->as_number() : std::nullopt);
  if (pts.empty()) throw ComponentError::missing("P");
  if (ids.empty()) throw ComponentError::missing("I");
  if (pts.size() != ids.size())
    throw ComponentError("CSV Bars: P and I must be equal-length lists (" + std::to_string(pts.size()) + " points, " +
                         std::to_string(ids.size()) + " ids)");
  if (trim(path).empty()) throw ComponentError("CSV Bars: no file selected");

  // id -> point (last wins on a duplicate id, which is warned, never silent).
  std::vector<std::string> warnings;
  std::map<std::uint64_t, Vec3> lookup;
  std::vector<double> dup_nodes;
  for (std::size_t k = 0; k < pts.size(); ++k) {
    if (!pts[k] || !ids[k]) continue;
    const bool existed = lookup.count(id_key(*ids[k])) > 0;
    lookup[id_key(*ids[k])] = *pts[k];
    if (existed && std::none_of(dup_nodes.begin(), dup_nodes.end(), [&](double d) { return id_key(d) == id_key(*ids[k]); }))
      dup_nodes.push_back(*ids[k]);
  }
  for (double d : dup_nodes) warnings.push_back("duplicate node id " + fmt_id(d) + ": last occurrence wins");

  const std::string text = read_file(path, "CSV Bars");
  const std::vector<Record> rows = parse_csv(text, resolve_delimiter(text, delimiter));
  const std::size_t hdr = has_header ? 1 : 0;
  if (rows.size() <= hdr) throw ComponentError("CSV Bars: no data rows in '" + path + "'");

  const Path b0 = Path::single(0);
  Tree lines, bar_ids, sections, materials;
  for (Tree* t : {&lines, &bar_ids, &sections, &materials}) t->ensure_path(b0);
  auto resolve_end = [&](const Record& rec, std::size_t col, const char* end, std::size_t line,
                         std::string& warn) -> std::optional<Vec3> {
    const std::string* c = cell(rec, col);
    if (auto v = c ? parse_num(*c) : std::nullopt) {
      auto it = lookup.find(id_key(*v));
      if (it != lookup.end()) return it->second;
      warn = "row " + std::to_string(line) + ": node id " + fmt_id(*v) + " not found";
      return std::nullopt;
    }
    warn = "row " + std::to_string(line) + ": " + end + " '" + (c ? trim(*c) : std::string()) + "' is not a valid node id";
    return std::nullopt;
  };
  for (std::size_t i = hdr; i < rows.size(); ++i) {
    const Record& rec = rows[i];
    const std::size_t rown = i - hdr + 1, line = rown + hdr;
    std::string wa, wb;
    const auto a = resolve_end(rec, map_a, "a", line, wa), b = resolve_end(rec, map_b, "b", line, wb);
    if (a && b) {
      lines.append(b0, Value::of(Line{*a, *b}));
    } else {
      if (!wa.empty()) warnings.push_back(wa);
      if (!wb.empty()) warnings.push_back(wb);
      lines.append_null(b0);
    }
    double id = static_cast<double>(rown);
    if (map_id) {
      const std::string* c = cell(rec, *map_id);
      if (auto v = c ? parse_num(*c) : std::nullopt)
        id = *v;
      else
        warnings.push_back("row " + std::to_string(line) + ": id '" + (c ? trim(*c) : std::string()) + "' is not a number");
    }
    bar_ids.append(b0, Value::number(id));
    // A mapped but blank section cell means the catalog default applies
    // downstream (Catalog Section's null-slot fallback): never silent.
    if (map_section) {
      const std::string* c = cell(rec, *map_section);
      const std::string tag = c ? trim(*c) : std::string();
      if (tag.empty()) {
        warnings.push_back("row " + std::to_string(line) + ": no section tag — catalog default applies");
        sections.append_null(b0);
      } else {
        sections.append(b0, Value::text(tag));
      }
    } else {
      sections.append_null(b0);
    }
    const std::string* mc = map_material ? cell(rec, *map_material) : nullptr;
    const std::string mtag = mc ? trim(*mc) : std::string();
    if (mtag.empty())
      materials.append_null(b0);
    else
      materials.append(b0, Value::text(mtag));
  }
  Tree w;
  for (auto& s : warnings) w.append(b0, Value::text(std::move(s)));
  return {lines, bar_ids, sections, materials, w};
}
json::Value CsvBarsNode::params() const {
  json::Value map = json::Value::object();
  map.set("id", col_json(map_id));
  map.set("a", json::Value(static_cast<double>(map_a)));
  map.set("b", json::Value(static_cast<double>(map_b)));
  map.set("section", col_json(map_section));
  map.set("material", col_json(map_material));
  json::Value p = json::Value::object();
  p.set("path", path);
  p.set("delimiter", delimiter);
  p.set("has_header", has_header);
  p.set("map", std::move(map));
  return p;
}
void CsvBarsNode::set_params(const json::Value& p) {
  read_common(p, path, delimiter, has_header);
  if (const auto* m = p.get("map")) {
    map_id = col_of(*m, "id");
    if (auto a = col_of(*m, "a")) map_a = *a;
    if (auto b = col_of(*m, "b")) map_b = *b;
    map_section = col_of(*m, "section");
    map_material = col_of(*m, "material");
  }
}

// ============================ Topology: the rest of the family ============================

namespace {

constexpr double kWeldTol = 1e-3;  // default weld tolerance [m] (= Assemble's)

const Line& line_of(const Value& v, const char* who) {
  const Line* l = v.as_line();
  if (!l) throw ComponentError(std::string(who) + ": expected a line, got " + v.type_tag());
  return *l;
}

// Closest points between segments p0-p1 and q0-q1 (Ericson, Real-Time
// Collision Detection 5.1.9): (s, t, cp, cq), robust to zero-length segments.
struct Closest {
  double s, t;
  Vec3 cp, cq;
};
Closest closest_seg_seg(const Vec3& p0, const Vec3& p1, const Vec3& q0, const Vec3& q1) {
  const Vec3 d1 = p1 - p0, d2 = q1 - q0, r = p0 - q0;
  const double a = d1.dot(d1), e = d2.dot(d2), f = d2.dot(r);
  constexpr double eps = 1e-18;
  double s, t;
  if (a <= eps && e <= eps) {
    s = t = 0;
  } else if (a <= eps) {
    s = 0;
    t = std::clamp(f / e, 0.0, 1.0);
  } else {
    const double c = d1.dot(r);
    if (e <= eps) {
      t = 0;
      s = std::clamp(-c / a, 0.0, 1.0);
    } else {
      const double b = d1.dot(d2), denom = a * e - b * b;
      s = denom > eps ? std::clamp((b * f - c * e) / denom, 0.0, 1.0) : 0.0;
      t = (b * s + f) / e;
      if (t < 0) {
        t = 0;
        s = std::clamp(-c / a, 0.0, 1.0);
      } else if (t > 1) {
        t = 1;
        s = std::clamp((b - c) / a, 0.0, 1.0);
      }
    }
  }
  return {s, t, p0 + d1 * s, q0 + d2 * t};
}

// Sort split records by parameter and drop any within ptol of the last kept.
std::vector<Vec3> dedup_splits(std::vector<std::pair<double, Vec3>> sp, double ptol) {
  std::stable_sort(sp.begin(), sp.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
  std::vector<Vec3> out;
  double last = -std::numeric_limits<double>::infinity();
  for (const auto& [p, pt] : sp)
    if (p - last > ptol) {
      out.push_back(pt);
      last = p;
    }
  return out;
}

void emit_segments(Tree& l_out, Tree& s_out, const Path& path, std::int64_t src, const Vec3& a, const Vec3& b,
                   const std::vector<Vec3>& pts) {
  Vec3 prev = a;
  for (const Vec3& pt : pts) {
    l_out.append(path, Value::of(Line{prev, pt}));
    s_out.append(path, Value::integer(src));
    prev = pt;
  }
  l_out.append(path, Value::of(Line{prev, b}));
  s_out.append(path, Value::integer(src));
}

std::int64_t cell(double x, double h) {  // saturating float -> int cast (NaN -> 0)
  const double v = std::floor(x / h);
  if (std::isnan(v)) return 0;
  if (v >= 9.2e18) return std::numeric_limits<std::int64_t>::max();
  if (v <= -9.2e18) return std::numeric_limits<std::int64_t>::min();
  return static_cast<std::int64_t>(v);
}

}  // namespace

ComponentInfo IntersectNode::info() const {
  return {"keta.intersect_node", "Intersect & Node", "Topology",
          "Split crossing / near-missing lines at shared points so they connect"};
}
std::vector<PinSpec> IntersectNode::inputs() const {
  return {PinSpec::tree("L", "line", "Lines").opt(), PinSpec::item("T", "number", "Tolerance (m) (default 0.001)").opt()};
}
std::vector<PinSpec> IntersectNode::outputs() const {
  return {PinSpec::tree("L", "line", "Split segments, in the source branch, ordered by (source line, then t)"),
          PinSpec::tree("S", "int", "Per-branch index of each segment's source line (the remap for re-attaching data)"),
          PinSpec::list("X", "point", "Intersection points, deduped, deterministic order")};
}
Outputs<Value> IntersectNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::size_t kMaxPairs = 2'000'000;
  constexpr std::uint64_t kMaxInsertions = 16'000'000;
  constexpr double kParallelEps = 1e-9;
  const double tol = first_number(input_or_empty(in, 1)).value_or(kWeldTol);
  struct Flat {
    Vec3 a, b;
    Path path;
    std::int64_t local;
    bool degenerate;
  };
  std::vector<Flat> segs;
  for (const auto& [path, branch] : input_or_empty(in, 0).branches())
    for (std::size_t j = 0; j < branch.size(); ++j) {
      if (!branch[j]) continue;
      const Line& l = line_of(*branch[j], "Intersect & Node");
      segs.push_back({l.a, l.b, path, static_cast<std::int64_t>(j), !((l.b - l.a).norm() > 0)});
    }
  if (segs.empty()) return {Tree{}, Tree{}, Tree{}};

  // Broad phase: tol-padded AABBs bucketed on a uniform grid (cell = mean
  // segment length), i < j co-cell pairs sorted + deduped.
  std::vector<std::pair<std::size_t, std::size_t>> pairs;
  std::vector<std::size_t> active;
  for (std::size_t i = 0; i < segs.size(); ++i)
    if (!segs[i].degenerate) active.push_back(i);
  if (active.size() >= 2) {
    double sum_len = 0;
    for (std::size_t i : active) sum_len += (segs[i].b - segs[i].a).norm();
    const double h = std::max({sum_len / static_cast<double>(active.size()), tol, 1e-9});
    std::map<std::array<std::int64_t, 3>, std::vector<std::size_t>> grid;
    std::uint64_t insertions = 0;
    for (std::size_t i : active) {
      const Flat& s = segs[i];
      const Vec3 lo = s.a.cwiseMin(s.b) - Vec3::Constant(tol), hi = s.a.cwiseMax(s.b) + Vec3::Constant(tol);
      std::array<std::int64_t, 3> c0{cell(lo.x(), h), cell(lo.y(), h), cell(lo.z(), h)};
      std::array<std::int64_t, 3> c1{cell(hi.x(), h), cell(hi.y(), h), cell(hi.z(), h)};
      std::uint64_t cells = 1;
      for (int k = 0; k < 3; ++k) {
        const double span = std::max(0.0, static_cast<double>(c1[static_cast<std::size_t>(k)]) -
                                              static_cast<double>(c0[static_cast<std::size_t>(k)]) + 1.0);
        cells = static_cast<std::uint64_t>(std::min(1.8e19, static_cast<double>(cells) * span));
      }
      insertions = insertions + cells < insertions ? std::numeric_limits<std::uint64_t>::max() : insertions + cells;
      if (insertions > kMaxInsertions)
        throw ComponentError(
            "Intersect & Node: line set too large or spread for the spatial grid; raise the tolerance or reduce the "
            "line count");
      for (std::int64_t cx = c0[0]; cx <= c1[0]; ++cx)
        for (std::int64_t cy = c0[1]; cy <= c1[1]; ++cy)
          for (std::int64_t cz = c0[2]; cz <= c1[2]; ++cz) grid[{cx, cy, cz}].push_back(i);
    }
    for (const auto& [key, bucket] : grid) {
      const std::size_t n = bucket.size();
      if (n < 2) continue;
      const std::uint64_t add = static_cast<std::uint64_t>(n) * (n - 1) / 2;
      if (pairs.size() + add > 4 * kMaxPairs)
        throw ComponentError(
            "Intersect & Node: too many candidate intersection pairs; raise the tolerance or reduce the line count");
      for (std::size_t x = 0; x < n; ++x)
        for (std::size_t y = x + 1; y < n; ++y) pairs.emplace_back(std::min(bucket[x], bucket[y]), std::max(bucket[x], bucket[y]));
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    if (pairs.size() > kMaxPairs)
      throw ComponentError("Intersect & Node: " + std::to_string(pairs.size()) +
                           " candidate pairs exceed the cap of " + std::to_string(kMaxPairs));
  }

  std::vector<std::vector<std::pair<double, Vec3>>> splits(segs.size());
  NodeMerger xmerger(tol);
  for (const auto& [i, j] : pairs) {
    const Flat &si = segs[i], &sj = segs[j];
    const Vec3 d1 = si.b - si.a, d2 = sj.b - sj.a;
    if (d1.cross(d2).norm() <= kParallelEps * d1.norm() * d2.norm()) continue;  // parallel: Weld's job
    const Closest c = closest_seg_seg(si.a, si.b, sj.a, sj.b);
    if ((c.cp - c.cq).norm() > tol) continue;
    const double mi = tol / d1.norm(), mj = tol / d2.norm();
    const bool ii = c.s > mi && c.s < 1 - mi, ji = c.t > mj && c.t < 1 - mj;
    if (!ii && !ji) continue;  // endpoint-to-endpoint touch
    // One shared point: the closest-approach midpoint for a crossing, the
    // toucher's own vertex for a T-junction (so Weld can connect the stem).
    const Vec3 shared = ii && ji ? Vec3((c.cp + c.cq) / 2.0)
                        : ii     ? (c.t <= 0.5 ? sj.a : sj.b)
                                 : (c.s <= 0.5 ? si.a : si.b);
    xmerger.id_for(shared);
    if (ii) splits[i].emplace_back(c.s, shared);
    if (ji) splits[j].emplace_back(c.t, shared);
  }
  Tree l_out, s_out, x_out;
  for (std::size_t k = 0; k < segs.size(); ++k) {
    const double len = (segs[k].b - segs[k].a).norm();
    const double ptol = len > 0 ? tol / len : std::numeric_limits<double>::infinity();
    emit_segments(l_out, s_out, segs[k].path, segs[k].local, segs[k].a, segs[k].b, dedup_splits(std::move(splits[k]), ptol));
  }
  for (const Vec3& p : xmerger.positions()) x_out.append(Path::single(0), Value::of(Point{p}));
  return {l_out, s_out, x_out};
}

ComponentInfo DivideLineNode::info() const {
  return {"keta.divide_line", "Divide", "Topology", "Divide each line into N equal segments (one child branch per line)"};
}
std::vector<PinSpec> DivideLineNode::inputs() const {
  return {PinSpec::tree("L", "line", "Lines").opt(), PinSpec::item("N", "int", "Segments per line (default 2)").opt()};
}
std::vector<PinSpec> DivideLineNode::outputs() const {
  return {PinSpec::tree("P", "point", "Division points: branch {p;n} for line n (its item index) of input branch {p}, N+1 points"),
          PinSpec::tree("S", "line", "Sub-segments: branch {p;n}, N segments")};
}
Outputs<Value> DivideLineNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::int64_t kMaxN = 10'000;
  constexpr std::uint64_t kMaxPoints = 1'000'000;
  const std::int64_t n = first_int(input_or_empty(in, 1)).value_or(2);
  if (n < 1 || n > kMaxN)
    throw ComponentError("Divide: N must be in 1..=" + std::to_string(kMaxN) + ", got " + std::to_string(n));
  const Tree& L = input_or_empty(in, 0);
  std::uint64_t lines = 0;
  for (const Value& v : L.non_nulls()) lines += v.as_line() != nullptr;
  const std::uint64_t total = lines * static_cast<std::uint64_t>(n + 1);
  if (total > kMaxPoints)
    throw ComponentError("Divide: " + std::to_string(total) + " emitted points exceed the cap of " +
                         std::to_string(kMaxPoints));
  Tree p_out, s_out;
  for (const auto& [path, branch] : L.branches())
    for (std::size_t j = 0; j < branch.size(); ++j) {
      if (!branch[j]) continue;
      const Line& l = line_of(*branch[j], "Divide");
      const Path child = path.appended(static_cast<int>(j));
      const Vec3 ab = l.b - l.a;
      std::vector<Vec3> pts;
      for (std::int64_t k = 0; k <= n; ++k)
        pts.push_back(k == 0 ? l.a : k == n ? l.b : Vec3(l.a + ab * (static_cast<double>(k) / static_cast<double>(n))));
      for (const Vec3& p : pts) p_out.append(child, Value::of(Point{p}));
      for (std::size_t k = 0; k + 1 < pts.size(); ++k) s_out.append(child, Value::of(Line{pts[k], pts[k + 1]}));
    }
  return {p_out, s_out};
}

ComponentInfo ShatterNode::info() const {
  return {"keta.shatter", "Shatter", "Topology", "Split lines at the projection of nearby points onto them"};
}
std::vector<PinSpec> ShatterNode::inputs() const {
  return {PinSpec::tree("L", "line", "Lines").opt(), PinSpec::list("P", "point", "Split points (projected onto each line)"),
          PinSpec::item("T", "number", "Tolerance (m) (default 0.001)").opt()};
}
std::vector<PinSpec> ShatterNode::outputs() const {
  return {PinSpec::tree("L", "line", "Sub-segments, in the source branch (as Intersect & Node)"),
          PinSpec::tree("S", "int", "Per-branch index of each segment's source line")};
}
Outputs<Value> ShatterNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::uint64_t kMaxPairs = 5'000'000;
  std::vector<Vec3> targets;
  for (const Value& v : input_or_empty(in, 1).non_nulls())
    if (auto p = v.as_point()) targets.push_back(*p);
  const double tol = first_number(input_or_empty(in, 2)).value_or(kWeldTol);
  const Tree& L = input_or_empty(in, 0);
  std::uint64_t lines = 0;
  for (const Value& v : L.non_nulls()) lines += v.as_line() != nullptr;
  const std::uint64_t pairs = lines * targets.size();
  if (pairs > kMaxPairs)
    throw ComponentError("Shatter: " + std::to_string(pairs) + " line×point pairs exceed the cap of " +
                         std::to_string(kMaxPairs));
  Tree l_out, s_out;
  for (const auto& [path, branch] : L.branches())
    for (std::size_t j = 0; j < branch.size(); ++j) {
      if (!branch[j]) continue;
      const Line& l = line_of(*branch[j], "Shatter");
      const Vec3 ab = l.b - l.a;
      const double len2 = ab.dot(ab);
      std::vector<std::pair<double, Vec3>> sp;
      if (len2 > 0) {
        const double margin = tol / std::sqrt(len2);
        for (const Vec3& tgt : targets) {
          const double tp = (tgt - l.a).dot(ab) / len2;
          if (tp <= margin || tp >= 1 - margin) continue;  // no slivers at the ends
          const Vec3 proj = l.a + ab * tp;
          if ((proj - tgt).norm() <= tol) sp.emplace_back(tp, proj);
        }
      }
      const double ptol = len2 > 0 ? tol / std::sqrt(len2) : std::numeric_limits<double>::infinity();
      emit_segments(l_out, s_out, path, static_cast<std::int64_t>(j), l.a, l.b, dedup_splits(std::move(sp), ptol));
    }
  return {l_out, s_out};
}

ComponentInfo CullShortNode::info() const {
  return {"keta.cull_short", "Cull Short Lines", "Topology", "Drop lines shorter than T (e.g. slivers left by welding)"};
}
std::vector<PinSpec> CullShortNode::inputs() const {
  return {PinSpec::tree("L", "line", "Lines").opt(),
          PinSpec::item("T", "number", "Minimum length (m) (default 0.001)").opt()};
}
std::vector<PinSpec> CullShortNode::outputs() const {
  return {PinSpec::tree("L", "line", "Kept lines (paths preserved, dropped items removed)"),
          PinSpec::item("C", "int", "Number of lines culled"),
          PinSpec::tree("M", "bool", "Mask aligned with input (true = kept, nulls stay null)")};
}
Outputs<Value> CullShortNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const double tol = first_number(input_or_empty(in, 1)).value_or(kWeldTol);
  Tree l_out, m_out;
  std::int64_t culled = 0;
  for (const auto& [path, branch] : input_or_empty(in, 0).branches()) {
    l_out.ensure_path(path);
    m_out.ensure_path(path);
    for (const auto& slot : branch) {
      if (!slot) {
        m_out.append_null(path);
        continue;
      }
      const Line& l = line_of(*slot, "Cull Short Lines");
      const bool keep = l.length() >= tol;
      m_out.append(path, Value::boolean(keep));
      if (keep)
        l_out.append(path, Value::of(l));
      else
        ++culled;
    }
  }
  return {l_out, single(Value::integer(culled)), m_out};
}

ComponentInfo EndPointsNode::info() const {
  return {"keta.end_points", "End Points", "Topology", "Extract each line's start and end points (aligned with input)"};
}
std::vector<PinSpec> EndPointsNode::inputs() const { return {PinSpec::tree("L", "line", "Lines").opt()}; }
std::vector<PinSpec> EndPointsNode::outputs() const {
  return {PinSpec::tree("A", "point", "Start points (aligned with input, nulls stay null)"),
          PinSpec::tree("B", "point", "End points (aligned with input, nulls stay null)")};
}
Outputs<Value> EndPointsNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree a_out, b_out;
  for (const auto& [path, branch] : input_or_empty(in, 0).branches()) {
    a_out.ensure_path(path);
    b_out.ensure_path(path);
    for (const auto& slot : branch) {
      if (!slot) {
        a_out.append_null(path);
        b_out.append_null(path);
        continue;
      }
      const Line& l = line_of(*slot, "End Points");
      a_out.append(path, Value::of(Point{l.a}));
      b_out.append(path, Value::of(Point{l.b}));
    }
  }
  return {a_out, b_out};
}

}  // namespace keta::param

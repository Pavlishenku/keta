// Frame Grid: a whole building frame's geometry from six numbers. Downstream
// Beam / Support / Load nodes pair items by index, so the generation ORDER is
// part of the contract. Grid, Polar Grid and Connect: the composable lattice +
// wiring-rule family.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numbers>

#include "keta/param/components.h"

namespace keta::param {

ComponentInfo FrameGrid::info() const {
  return {"keta.frame_grid", "Frame Grid", "Geometry", "Multi-storey multi-bay frame lines"};
}
std::vector<PinSpec> FrameGrid::inputs() const {
  return {PinSpec::item("NX", "int", "Bays along X (default 4)").opt(),
          PinSpec::item("NY", "int", "Bays along Y (default 3)").opt(),
          PinSpec::item("NS", "int", "Storeys (default 6)").opt(),
          PinSpec::item("LX", "number", "Bay width X in m (default 6)").opt(),
          PinSpec::item("LY", "number", "Bay width Y in m (default 6)").opt(),
          PinSpec::item("H", "number", "Storey height in m (default 3.5)").opt()};
}
std::vector<PinSpec> FrameGrid::outputs() const {
  return {PinSpec::list("C", "line", "Column lines"), PinSpec::list("X", "line", "Beam lines along X"),
          PinSpec::list("Y", "line", "Beam lines along Y"), PinSpec::list("B", "point", "Base points (z = 0)"),
          PinSpec::list("W", "point", "Windward-face points (x = 0, above ground)")};
}
Outputs<Value> FrameGrid::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::int64_t nx = first_int(input_or_empty(in, 0)).value_or(4);
  const std::int64_t ny = first_int(input_or_empty(in, 1)).value_or(3);
  const std::int64_t ns = first_int(input_or_empty(in, 2)).value_or(6);
  const double lx = first_number(input_or_empty(in, 3)).value_or(6.0);
  const double ly = first_number(input_or_empty(in, 4)).value_or(6.0);
  const double h = first_number(input_or_empty(in, 5)).value_or(3.5);

  // Guards: a bad parameter is a typed component error, never a crash.
  for (auto [name, v] : {std::pair<const char*, std::int64_t>{"NX", nx}, {"NY", ny}, {"NS", ns}})
    if (v < 1 || v > kMaxDim)
      throw ComponentError(std::string(name) + " must be in 1..=" + std::to_string(kMaxDim) + ", got " +
                           std::to_string(v));
  for (auto [name, v] : {std::pair<const char*, double>{"LX", lx}, {"LY", ly}, {"H", h}})
    if (!(std::isfinite(v) && v > 0.0))
      throw ComponentError(std::string(name) + " must be finite and > 0, got " + Value::number(v).describe());

  const auto NX = static_cast<std::size_t>(nx), NY = static_cast<std::size_t>(ny), NS = static_cast<std::size_t>(ns);
  // Bar-count guard, checked BEFORE materializing anything.
  const std::size_t total = (NX + 1) * (NY + 1) * NS + NS * (NY + 1) * NX + NS * (NX + 1) * NY;
  if (total > kMaxBars)
    throw ComponentError("frame grid would generate " + std::to_string(total) + " bars, over the cap of " +
                         std::to_string(kMaxBars));

  const Path p0 = Path::single(0);
  auto line = [](Vec3 a, Vec3 b) { return Value::of(Line{a, b}); };
  auto at = [](std::size_t k, double step) { return static_cast<double>(k) * step; };

  // Columns: for j, i, s -> (i*LX, j*LY, s*H) up to (.., (s+1)*H).
  Tree columns;
  columns.ensure_path(p0);
  for (std::size_t j = 0; j <= NY; ++j)
    for (std::size_t i = 0; i <= NX; ++i)
      for (std::size_t s = 0; s < NS; ++s)
        columns.append(p0, line(Vec3(at(i, lx), at(j, ly), at(s, h)), Vec3(at(i, lx), at(j, ly), at(s + 1, h))));

  // X-beams: at each floor level, one bay along +X.
  Tree xb;
  xb.ensure_path(p0);
  for (std::size_t s = 1; s <= NS; ++s)
    for (std::size_t j = 0; j <= NY; ++j)
      for (std::size_t i = 0; i < NX; ++i)
        xb.append(p0, line(Vec3(at(i, lx), at(j, ly), at(s, h)), Vec3(at(i + 1, lx), at(j, ly), at(s, h))));

  // Y-beams: at each floor level, one bay along +Y.
  Tree yb;
  yb.ensure_path(p0);
  for (std::size_t s = 1; s <= NS; ++s)
    for (std::size_t i = 0; i <= NX; ++i)
      for (std::size_t j = 0; j < NY; ++j)
        yb.append(p0, line(Vec3(at(i, lx), at(j, ly), at(s, h)), Vec3(at(i, lx), at(j + 1, ly), at(s, h))));

  // Base points: the whole ground grid (support seats).
  Tree base;
  base.ensure_path(p0);
  for (std::size_t j = 0; j <= NY; ++j)
    for (std::size_t i = 0; i <= NX; ++i) base.append(p0, Value::of(Point{Vec3(at(i, lx), at(j, ly), 0.0)}));

  // Windward-face points: the x = 0 face above ground (wind load seats).
  Tree windward;
  windward.ensure_path(p0);
  for (std::size_t s = 1; s <= NS; ++s)
    for (std::size_t j = 0; j <= NY; ++j) windward.append(p0, Value::of(Point{Vec3(0.0, at(j, ly), at(s, h))}));

  return {columns, xb, yb, base, windward};
}

// ============================ Grid / Polar Grid / Connect ============================

namespace {

// Defensive pin readers: the FIRST non-null item of the whole tree,
// converted (a first item of the wrong kind falls back to the default, it is
// not skipped) — except bool and point, which take the first convertible one.
const Value* head(const std::vector<Tree>& in, std::size_t i) {
  if (i >= in.size()) return nullptr;
  for (const auto& kv : in[i].branches())
    for (const auto& o : kv.second)
      if (o) return &*o;
  return nullptr;
}
std::int64_t int_at(const std::vector<Tree>& in, std::size_t i, std::int64_t dflt) {
  const Value* v = head(in, i);
  return v ? v->as_int().value_or(dflt) : dflt;
}
double number_at(const std::vector<Tree>& in, std::size_t i, double dflt) {
  const Value* v = head(in, i);
  return v ? v->as_number().value_or(dflt) : dflt;
}
std::vector<double> numbers_at(const std::vector<Tree>& in, std::size_t i) {
  std::vector<double> out;
  if (i < in.size())
    for (const Value& v : in[i].non_nulls())
      if (auto x = v.as_number()) out.push_back(*x);
  return out;
}

// A double as the error messages print it: "NaN", "inf", shortest form.
std::string rs(double x) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
  char buf[40];
  for (int prec = 1; prec <= 17; ++prec) {
    std::snprintf(buf, sizeof buf, "%.*g", prec, x);
    if (std::strtod(buf, nullptr) == x) break;
  }
  return buf;
}

}  // namespace

ComponentInfo GridNode::info() const {
  return {"keta.grid", "Grid", "Geometry", "Rectilinear point lattice (line / rectangle / box)"};
}
std::vector<PinSpec> GridNode::inputs() const {
  return {PinSpec::item("NX", "int", "Points along X (default 2)").opt(),
          PinSpec::item("NY", "int", "Points along Y (default 1)").opt(),
          PinSpec::item("NZ", "int", "Points along Z / levels (default 1)").opt(),
          PinSpec::item("DX", "number", "Spacing along X in m (default 6)").opt(),
          PinSpec::item("DY", "number", "Spacing along Y in m (default 6)").opt(),
          PinSpec::item("DZ", "number", "Spacing along Z in m (default 3.5)").opt(),
          PinSpec::item("O", "point", "Origin (default world origin)").opt()};
}
std::vector<PinSpec> GridNode::outputs() const {
  return {PinSpec::tree("P", "point",
                        "Point lattice: branch {i;j} is the vertical fiber at column (i,j), holding NZ points; item "
                        "index k is the level, at O + (i·DX, j·DY, k·DZ). j varies fastest in branch order ({0;0} < "
                        "{0;1} < {1;0}). Graft for {i;j;k}; Simplify drops degenerate dimensions (NY=1 gives {i;0}). "
                        "Spacings may be zero or negative (mirrored / coincident grids); only NaN/inf is rejected.")};
}
Outputs<Value> GridNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::int64_t kMaxN = 1000, kMaxPoints = 1'000'000;
  const std::int64_t nx = int_at(in, 0, 2), ny = int_at(in, 1, 1), nz = int_at(in, 2, 1);
  const double dx = number_at(in, 3, 6.0), dy = number_at(in, 4, 6.0), dz = number_at(in, 5, 3.5);
  const Vec3 o = first_point(input_or_empty(in, 6)).value_or(Vec3::Zero());
  for (auto [name, v] : {std::pair<const char*, std::int64_t>{"NX", nx}, {"NY", ny}, {"NZ", nz}})
    if (v < 1 || v > kMaxN)
      throw ComponentError(std::string("Grid: ") + name + " must be in 1..=" + std::to_string(kMaxN) + ", got " +
                           std::to_string(v));
  for (auto [name, v] : {std::pair<const char*, double>{"DX", dx}, {"DY", dy}, {"DZ", dz}})
    if (!std::isfinite(v)) throw ComponentError(std::string("Grid: ") + name + " must be finite, got " + rs(v));
  const std::int64_t total = nx * ny * nz;
  if (total > kMaxPoints)
    throw ComponentError("Grid: " + std::to_string(total) + " points exceed the cap of " + std::to_string(kMaxPoints));
  Tree p;
  for (std::int64_t i = 0; i < nx; ++i)
    for (std::int64_t j = 0; j < ny; ++j) {
      const Path path{static_cast<int>(i), static_cast<int>(j)};
      for (std::int64_t k = 0; k < nz; ++k)
        p.append(path, Value::of(Point{Vec3(o.x() + static_cast<double>(i) * dx, o.y() + static_cast<double>(j) * dy,
                                            o.z() + static_cast<double>(k) * dz)}));
    }
  return {p};
}

ComponentInfo PolarGridNode::info() const {
  return {"keta.polar_grid", "Polar Grid", "Geometry", "Radial point lattice (ring / arc / tapered / twisted tower)"};
}
std::vector<PinSpec> PolarGridNode::inputs() const {
  return {PinSpec::item("N", "int", "Points per ring (default 8)").opt(),
          PinSpec::item("NZ", "int", "Levels (default 1)").opt(),
          PinSpec::list("R", "number", "Ring radius per level in m (repeat-last; default 5)").opt(),
          PinSpec::item("H", "number", "Level height in m (default 3.5)").opt(),
          PinSpec::item("A", "number", "Sweep in deg (default 360)").opt(),
          PinSpec::item("TW", "number", "Twist per level in deg (default 0)").opt(),
          PinSpec::item("O", "point", "Center (default world origin)").opt()};
}
std::vector<PinSpec> PolarGridNode::outputs() const {
  return {PinSpec::tree("P", "point",
                        "Point lattice: branch {i} is angular fiber i, holding NZ level points. Full sweep (|A| within "
                        "1e-9 of 360) steps θ = A·i/N with no duplicate seam; a partial sweep steps θ = A·i/(N−1) with "
                        "both ends inclusive (needs N ≥ 2). Level k has radius R[min(k, len−1)], height k·H and angle "
                        "offset k·TW: point = O + (r·cosθ, r·sinθ, k·H).")};
}
Outputs<Value> PolarGridNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::int64_t kMaxN = 10'000, kMaxNz = 1000, kMaxPoints = 1'000'000;
  constexpr double kFullTol = 1e-9;
  const std::int64_t n = int_at(in, 0, 8), nz = int_at(in, 1, 1);
  std::vector<double> radii = numbers_at(in, 2);
  if (radii.empty()) radii.push_back(5.0);
  const double h = number_at(in, 3, 3.5), a = number_at(in, 4, 360.0), tw = number_at(in, 5, 0.0);
  const Vec3 o = first_point(input_or_empty(in, 6)).value_or(Vec3::Zero());
  if (n < 1 || n > kMaxN)
    throw ComponentError("Polar Grid: N must be in 1..=" + std::to_string(kMaxN) + ", got " + std::to_string(n));
  if (nz < 1 || nz > kMaxNz)
    throw ComponentError("Polar Grid: NZ must be in 1..=" + std::to_string(kMaxNz) + ", got " + std::to_string(nz));
  const std::int64_t total = n * nz;
  if (total > kMaxPoints)
    throw ComponentError("Polar Grid: " + std::to_string(total) + " points exceed the cap of " +
                         std::to_string(kMaxPoints));
  for (auto [name, v] : {std::pair<const char*, double>{"H", h}, {"A", a}, {"TW", tw}})
    if (!std::isfinite(v)) throw ComponentError(std::string("Polar Grid: ") + name + " must be finite, got " + rs(v));
  for (std::size_t k = 0; k < radii.size(); ++k)
    if (!(std::isfinite(radii[k]) && radii[k] >= 0))
      throw ComponentError("Polar Grid: R[" + std::to_string(k) + "] must be finite and ≥ 0, got " + rs(radii[k]));
  const bool full = std::fabs(a - 360.0) < kFullTol || std::fabs(a + 360.0) < kFullTol;
  if (!full && n < 2)
    throw ComponentError("Polar Grid: a partial sweep (A ≠ ±360°) needs N ≥ 2 to include both ends");
  const double denom = full ? static_cast<double>(n) : static_cast<double>(n - 1);
  Tree p;
  for (std::int64_t i = 0; i < n; ++i) {
    const Path path = Path::single(static_cast<int>(i));
    const double theta_base = a * static_cast<double>(i) / denom;
    for (std::int64_t k = 0; k < nz; ++k) {
      const double r = radii[std::min<std::size_t>(static_cast<std::size_t>(k), radii.size() - 1)];
      const double theta = (theta_base + static_cast<double>(k) * tw) * std::numbers::pi / 180.0;
      p.append(path, Value::of(Point{Vec3(o.x() + r * std::cos(theta), o.y() + r * std::sin(theta),
                                          o.z() + static_cast<double>(k) * h)}));
    }
  }
  return {p};
}

ComponentInfo ConnectNode::info() const {
  return {"keta.connect", "Connect", "Geometry", "Connect lattice points into lines by offsetting the tree address"};
}
std::vector<PinSpec> ConnectNode::inputs() const {
  return {PinSpec::tree("P", "point", "Grid points (branch = fiber, item = level)").opt(),
          PinSpec::item("DI", "int", "Offset on the 2nd-to-last path index (default 0)").opt(),
          PinSpec::item("DJ", "int", "Offset on the last path index (default 0)").opt(),
          PinSpec::item("DK", "int", "Offset on the item / level index (default 1)").opt(),
          PinSpec::item("W", "bool", "Wrap the last path index (polar seam) (default false)").opt()};
}
std::vector<PinSpec> ConnectNode::outputs() const {
  return {PinSpec::tree("L", "line",
                        "Connecting lines, in the source branch, in (branch, item) order. Zero-length results "
                        "(coincident points) pass through as-is.")};
}
Outputs<Value> ConnectNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  constexpr std::size_t kMaxLines = 1'000'000;
  const Tree& p = input_or_empty(in, 0);
  const std::int64_t di = int_at(in, 1, 0), dj = int_at(in, 2, 0), dk = int_at(in, 3, 1);
  const bool w = first_bool(input_or_empty(in, 4)).value_or(false);
  if (di == 0 && dj == 0 && dk == 0)
    throw ComponentError("Connect: DI, DJ and DK are all zero — nothing to connect (set at least one offset)");
  std::size_t src_count = 0;
  for (const auto& kv : p.branches())
    for (const auto& o : kv.second) src_count += o.has_value();
  if (src_count > kMaxLines)
    throw ComponentError("Connect: " + std::to_string(src_count) + " source points exceed the emitted-line cap of " +
                         std::to_string(kMaxLines));

  // Positional addressing: the ascending distinct last index per sibling
  // group, and the ascending distinct second-to-last index per parent group.
  using Key = std::vector<int>;
  std::map<Key, std::vector<int>> last_by_prefix, second_by_prefix2;
  for (const auto& kv : p.branches()) {
    const std::vector<int>& ix = kv.first.indices();
    const std::size_t l = ix.size();
    if (l >= 1) last_by_prefix[Key(ix.begin(), ix.end() - 1)].push_back(ix[l - 1]);
    if (l >= 2) {
      std::vector<int>& g = second_by_prefix2[Key(ix.begin(), ix.end() - 2)];
      if (g.empty() || g.back() != ix[l - 2]) g.push_back(ix[l - 2]);
    }
  }
  auto as_pt = [](const Value& v) {
    auto pt = v.as_point();
    if (!pt) throw ComponentError("Connect: expected a point, got " + v.type_tag());
    return *pt;
  };
  Tree out;
  for (const auto& [path, branch] : p.branches()) {
    const std::vector<int>& ix = path.indices();
    const std::size_t l = ix.size();
    if (di != 0 && l < 2)
      throw ComponentError("Connect: tree paths have no second dimension (DI requires a {i;j}+ path)");
    if (l == 0) continue;
    std::optional<int> target_second;
    if (l >= 2) {
      const std::vector<int>& g = second_by_prefix2.at(Key(ix.begin(), ix.end() - 2));
      const auto q = static_cast<std::int64_t>(std::find(g.begin(), g.end(), ix[l - 2]) - g.begin());
      const std::int64_t ti = q + di;
      if (ti < 0 || ti >= static_cast<std::int64_t>(g.size())) continue;
      target_second = g[static_cast<std::size_t>(ti)];
    }
    const std::vector<int>& lasts = last_by_prefix.at(Key(ix.begin(), ix.end() - 1));
    const auto len = static_cast<std::int64_t>(lasts.size());
    const auto q = static_cast<std::int64_t>(std::find(lasts.begin(), lasts.end(), ix[l - 1]) - lasts.begin());
    std::int64_t tj = q + dj;
    if (w) tj = ((tj % len) + len) % len;  // rem_euclid
    if (tj < 0 || tj >= len) continue;
    std::vector<int> target = ix;
    if (target_second) target[l - 2] = *target_second;
    target[l - 1] = lasts[static_cast<std::size_t>(tj)];
    const Branch<Value>* tb = p.branch(Path(std::move(target)));
    if (!tb) continue;
    for (std::size_t k = 0; k < branch.size(); ++k) {
      if (!branch[k]) continue;
      const Vec3 a = as_pt(*branch[k]);
      const std::int64_t tk = static_cast<std::int64_t>(k) + dk;
      if (tk < 0 || tk >= static_cast<std::int64_t>(tb->size())) continue;
      const auto& dst = (*tb)[static_cast<std::size_t>(tk)];
      if (!dst) continue;
      out.append(path, Value::of(Line{a, as_pt(*dst)}));
    }
  }
  return {out};
}

}  // namespace keta::param

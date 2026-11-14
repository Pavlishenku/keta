// Params, Geometry and Model components.
#include <algorithm>
#include <cctype>
#include <cmath>

#include "keta/param/components.h"

namespace keta::param {

namespace {

json::Value obj() { return json::Value::object(); }

double num_param(const json::Value& p, const char* key, double cur) {
  const auto* v = p.get(key);
  return v && v->is_number() ? v->as_number() : cur;
}
std::int64_t int_param(const json::Value& p, const char* key, std::int64_t cur) {
  const auto* v = p.get(key);
  return v && v->is_number() ? v->as_int() : cur;
}
std::string str_param(const json::Value& p, const char* key, const std::string& cur) {
  const auto* v = p.get(key);
  return v && v->is_string() ? v->as_string() : cur;
}

Vec3 point_of(const Value& v, const char* pin) {
  if (auto p = v.as_point()) return *p;
  wrong_type(pin, "point", v);
}
Vec3 vector_of(const Value& v, const char* pin) {
  if (auto p = v.as_vector()) return *p;
  wrong_type(pin, "vector", v);
}
double number_of(const Value& v, const char* pin) {
  if (auto x = v.as_number()) return *x;
  wrong_type(pin, "number", v);
}

// Three optional number pins (default 0) -> one Vec3 per matched item.
template <class Make>
Tree xyz(const std::vector<Tree>& in, Make make) {
  Tree x = number_or(in[0], 0.0), y = number_or(in[1], 0.0), z = number_or(in[2], 0.0);
  return map_items({{"X", &x}, {"Y", &y}, {"Z", &z}}, [&](const std::vector<const Value*>& it) {
    return make(Vec3(number_of(*it[0], "X"), number_of(*it[1], "Y"), number_of(*it[2], "Z")));
  });
}

std::vector<Line> require_lines(const Tree& t, const char* pin, bool accept_lines = false) {
  auto lines = beam_lines(t, accept_lines);
  if (lines.empty()) throw ComponentError::missing(pin);
  return lines;
}

}  // namespace

LoadFrame parse_frame(const std::string& s) {
  if (s == "local") return LoadFrame::Local;
  if (s == "projected") return LoadFrame::GlobalProjected;
  return LoadFrame::Global;
}

// ============================ Params ============================

ComponentInfo NumberSlider::info() const { return {"keta.slider", "Slider", "Params", "A draggable number"}; }
std::vector<PinSpec> NumberSlider::inputs() const { return {}; }
std::vector<PinSpec> NumberSlider::outputs() const { return {PinSpec::item("N", "number", "The current value")}; }
Outputs<Value> NumberSlider::solve(const SolveCtx&, const std::vector<Tree>&) const {
  return {single(Value::number(value))};
}
json::Value NumberSlider::params() const {
  json::Value p = obj();
  p.set("value", value);
  p.set("min", min);
  p.set("max", max);
  p.set("label", label);
  return p;
}
void NumberSlider::set_params(const json::Value& p) {
  value = num_param(p, "value", value);
  min = num_param(p, "min", min);
  max = num_param(p, "max", max);
  label = str_param(p, "label", label);
}

ComponentInfo IntSlider::info() const { return {"keta.intslider", "Int", "Params", "A draggable integer"}; }
std::vector<PinSpec> IntSlider::inputs() const { return {}; }
std::vector<PinSpec> IntSlider::outputs() const { return {PinSpec::item("N", "int", "The current value")}; }
Outputs<Value> IntSlider::solve(const SolveCtx&, const std::vector<Tree>&) const {
  return {single(Value::integer(value))};
}
json::Value IntSlider::params() const {
  json::Value p = obj();
  p.set("value", value);
  p.set("min", min);
  p.set("max", max);
  p.set("label", label);
  return p;
}
void IntSlider::set_params(const json::Value& p) {
  value = int_param(p, "value", value);
  min = int_param(p, "min", min);
  max = int_param(p, "max", max);
  label = str_param(p, "label", label);
}

ComponentInfo Panel::info() const { return {"keta.panel", "Panel", "Params", "Display whatever flows in"}; }
std::vector<PinSpec> Panel::inputs() const { return {PinSpec::tree("In", "any", "Anything").opt()}; }
std::vector<PinSpec> Panel::outputs() const { return {PinSpec::tree("Out", "any", "The input, unchanged")}; }
Outputs<Value> Panel::solve(const SolveCtx&, const std::vector<Tree>& in) const { return {in[0]}; }

// ============================ Geometry ============================

ComponentInfo PointXyz::info() const { return {"keta.point", "Point", "Geometry", "Point from x, y, z"}; }
std::vector<PinSpec> PointXyz::inputs() const {
  return {PinSpec::item("X", "number", "X coordinate (default 0)").opt(),
          PinSpec::item("Y", "number", "Y coordinate (default 0)").opt(),
          PinSpec::item("Z", "number", "Z coordinate (default 0)").opt()};
}
std::vector<PinSpec> PointXyz::outputs() const { return {PinSpec::item("P", "point", "The point")}; }
Outputs<Value> PointXyz::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {xyz(in, [](const Vec3& v) { return Value::of(Point{v}); })};
}

ComponentInfo VectorXyz::info() const { return {"keta.vector", "Vector", "Geometry", "Vector from x, y, z"}; }
std::vector<PinSpec> VectorXyz::inputs() const {
  return {PinSpec::item("X", "number", "X component (default 0)").opt(),
          PinSpec::item("Y", "number", "Y component (default 0)").opt(),
          PinSpec::item("Z", "number", "Z component (default 0)").opt()};
}
std::vector<PinSpec> VectorXyz::outputs() const { return {PinSpec::item("V", "vector", "The vector")}; }
Outputs<Value> VectorXyz::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {xyz(in, [](const Vec3& v) { return Value::of(Vector{v}); })};
}

ComponentInfo LineAB::info() const { return {"keta.line", "Line", "Geometry", "Line between two points"}; }
std::vector<PinSpec> LineAB::inputs() const {
  return {PinSpec::item("A", "point", "Start point"), PinSpec::item("B", "point", "End point")};
}
std::vector<PinSpec> LineAB::outputs() const { return {PinSpec::item("L", "line", "The line")}; }
Outputs<Value> LineAB::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {map_items({{"A", &in[0]}, {"B", &in[1]}}, [](const std::vector<const Value*>& it) {
    return Value::of(Line{point_of(*it[0], "A"), point_of(*it[1], "B")});
  })};
}

ComponentInfo Series::info() const {
  return {"keta.series", "Series", "Geometry", "start + step*i for i in 0..count"};
}
std::vector<PinSpec> Series::inputs() const {
  return {PinSpec::item("S", "number", "Start (default 0)").opt(),
          PinSpec::item("D", "number", "Step (default 1)").opt(),
          PinSpec::item("N", "int", "Count (default 10)").opt()};
}
std::vector<PinSpec> Series::outputs() const { return {PinSpec::list("R", "number", "The series")}; }
Outputs<Value> Series::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  double start = first_number(in[0]).value_or(0.0), step = first_number(in[1]).value_or(1.0);
  std::int64_t count = first_int(in[2]).value_or(10);
  if (count < 0 || count > 1'000'000)
    throw ComponentError("count " + std::to_string(count) + " out of range 0..=1e6");
  Tree t;
  t.ensure_path(Path::single(0));
  for (std::int64_t i = 0; i < count; ++i)
    t.append(Path::single(0), Value::number(start + step * static_cast<double>(i)));
  return {t};
}

// ============================ Model ============================

Material MaterialNode::material() const {
  Material m;
  if (preset == "steel") {
    m.name = "Steel";
    m.E = 210e9;
    m.nu = 0.3;
    m.rho = 7850.0;
    m.alpha = 1.2e-5;
  } else if (preset == "concrete") {
    m.name = "Concrete";
    m.E = 33e9;
    m.nu = 0.2;
    m.rho = 2500.0;
    m.alpha = 1.0e-5;
  } else {
    m.name = "Custom";
    m.E = e;
    m.nu = nu;
    m.rho = rho;
    m.alpha = 0.0;  // no dilation coefficient: a thermal load on it is a typed error
  }
  return m;
}
ComponentInfo MaterialNode::info() const {
  return {"keta.material", "Material", "Model", "Steel / concrete / custom material"};
}
std::vector<PinSpec> MaterialNode::inputs() const { return {}; }
std::vector<PinSpec> MaterialNode::outputs() const { return {PinSpec::item("M", "material", "The material")}; }
Outputs<Value> MaterialNode::solve(const SolveCtx&, const std::vector<Tree>&) const {
  return {single(Value::of(MaterialSpec{material()}))};
}
json::Value MaterialNode::params() const {
  json::Value p = obj();
  p.set("preset", preset);
  p.set("e", e);
  p.set("nu", nu);
  p.set("rho", rho);
  return p;
}
void MaterialNode::set_params(const json::Value& p) {
  preset = str_param(p, "preset", preset);
  e = num_param(p, "e", e);
  nu = num_param(p, "nu", nu);
  rho = num_param(p, "rho", rho);
}

ComponentInfo SectionRect::info() const { return {"keta.section_rect", "Rect Section", "Model", "Solid rectangle b x h"}; }
std::vector<PinSpec> SectionRect::inputs() const {
  return {PinSpec::item("B", "number", "Width along y (m)"), PinSpec::item("H", "number", "Height along z (m)")};
}
std::vector<PinSpec> SectionRect::outputs() const { return {PinSpec::item("S", "section", "The section")}; }
Outputs<Value> SectionRect::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {map_items({{"B", &in[0]}, {"H", &in[1]}}, [](const std::vector<const Value*>& it) {
    double b = number_of(*it[0], "B"), h = number_of(*it[1], "H");
    Section s;
    char name[64];
    std::snprintf(name, sizeof name, "Rect %.3fx%.3f", b, h);
    s.name = name;
    s.shape = ShapeRect{b, h};
    try {
      s.finalize();
    } catch (const std::exception& e) {
      throw ComponentError(std::string("section error: ") + e.what());
    }
    return Value::of(SectionSpec{s});
  })};
}

ComponentInfo BeamNode::info() const { return {"keta.beam", "Beam", "Model", "Line -> beam (or truss bar)"}; }
std::vector<PinSpec> BeamNode::inputs() const {
  return {PinSpec::item("L", "line", "Axis line"), PinSpec::item("S", "section", "Cross-section"),
          PinSpec::item("M", "material", "Material"),
          PinSpec::item("E", "vector", "Section eccentricity [ex, ey, ez] (m; global)").opt()};
}
std::vector<PinSpec> BeamNode::outputs() const { return {PinSpec::item("B", "beam", "The element")}; }
Outputs<Value> BeamNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const Tree& e = input_or_empty(in, 3);
  std::vector<std::pair<std::string, const Tree*>> pins = {{"L", &in[0]}, {"S", &in[1]}, {"M", &in[2]}};
  const bool has_e = has_values(e);
  if (has_e) pins.emplace_back("E", &e);
  const bool tr = truss;
  return {map_items(pins, [&](const std::vector<const Value*>& it) {
    BeamSpec b;
    const Line* l = it[0]->as_line();
    if (!l) wrong_type("L", "line", *it[0]);
    const auto* s = it[1]->get<SectionSpec>();
    if (!s) wrong_type("S", "section", *it[1]);
    const auto* m = it[2]->get<MaterialSpec>();
    if (!m) wrong_type("M", "material", *it[2]);
    b.line = *l;
    b.section = *s;
    b.material = *m;
    b.truss = tr;
    if (has_e) b.ecc = vector_of(*it[3], "E");
    return Value::of(b);
  })};
}
json::Value BeamNode::params() const {
  json::Value p = obj();
  p.set("truss", truss);
  return p;
}
void BeamNode::set_params(const json::Value& p) {
  if (const auto* t = p.get("truss"); t && t->is_bool()) truss = t->as_bool();
}

ComponentInfo CableNode::info() const {
  return {"keta.cable", "Cable", "Model", "Line -> elastic-catenary cable (EA, w, L0 route)"};
}
std::vector<PinSpec> CableNode::inputs() const {
  return {PinSpec::item("L", "line", "Anchor line (the two cable ends)"),
          PinSpec::item("EA", "number", "Axial rigidity EA (N)").opt(),
          PinSpec::item("w", "number", "Weight per unstrained length (N/m)").opt(),
          PinSpec::item("V", "number", "Route value: L0 (m) or H / T (N)").opt()};
}
std::vector<PinSpec> CableNode::outputs() const { return {PinSpec::item("C", "cable", "The cable element")}; }
Outputs<Value> CableNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  CableSpec proto;
  proto.ea = first_number(input_or_empty(in, 1)).value_or(ea);
  proto.w = first_number(input_or_empty(in, 2)).value_or(w);
  proto.value = first_number(input_or_empty(in, 3)).value_or(value);
  proto.route = route == "length"        ? CableRoute::Length
                : route == "end_tension" ? CableRoute::EndTension
                                         : CableRoute::HorizontalTension;
  return {map_items({{"L", &in[0]}}, [&](const std::vector<const Value*>& it) {
    const Line* l = it[0]->as_line();
    if (!l) throw ComponentError("L must be a line");
    CableSpec c = proto;
    c.a = l->a;
    c.b = l->b;
    return Value::of(c);
  })};
}
json::Value CableNode::params() const {
  json::Value p = obj();
  p.set("route", route);
  p.set("ea", ea);
  p.set("w", w);
  p.set("value", value);
  return p;
}
void CableNode::set_params(const json::Value& p) {
  if (const auto* x = p.get("route"); x && x->is_string()) route = x->as_string();
  if (const auto* x = p.get("ea"); x && x->is_number()) ea = x->as_number();
  if (const auto* x = p.get("w"); x && x->is_number()) w = x->as_number();
  if (const auto* x = p.get("value"); x && x->is_number()) value = x->as_number();
}

// ============================ Slab / Surface Load / Surface Thermal ============================

namespace {
bool finite3(const Vec3& p) { return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z()); }
std::string rs_num(double x) {
  if (std::isnan(x)) return "NaN";
  if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
  return Value::number(x).describe();
}
std::vector<Vec3> slab_outline_of(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (const auto* s = o->get<SlabSpec>()) return s->outline;
  throw ComponentError::missing("S");
}
SurfaceFrame parse_surface_frame(std::string f) {
  f.erase(0, f.find_first_not_of(" \t\r\n"));
  f.erase(f.find_last_not_of(" \t\r\n") + 1);
  std::transform(f.begin(), f.end(), f.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (f == "localnormal" || f == "local" || f == "normal") return SurfaceFrame::LocalNormal;
  if (f == "projected" || f == "proj") return SurfaceFrame::Projected;
  return SurfaceFrame::Global;
}
const char* surface_frame_str(SurfaceFrame f) {
  switch (f) {
    case SurfaceFrame::LocalNormal: return "localnormal";
    case SurfaceFrame::Projected: return "projected";
    default: return "global";
  }
}
}  // namespace

ComponentInfo SlabNode::info() const {
  return {"keta.slab", "Slab", "Model", "A meshable slab: outline + holes + thickness + material"};
}
std::vector<PinSpec> SlabNode::inputs() const {
  return {PinSpec::list("P", "point", "Outline ring points (>= 3, not closed)"),
          PinSpec::tree("H", "point", "Hole rings (one branch per hole)").opt(),
          PinSpec::item("T", "number", "Thickness t (m, > 0)"),
          PinSpec::item("M", "material", "Slab material"),
          PinSpec::item("S", "number", "Target mesh edge length (m, > 0; default 0.5)").opt(),
          PinSpec::item("LX", "vector", "Local-x axis override (global frame, direction; optional)").opt()};
}
std::vector<PinSpec> SlabNode::outputs() const { return {PinSpec::item("S", "slab", "The slab")}; }
Outputs<Value> SlabNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  SlabSpec spec;
  for (const Value& v : input_or_empty(in, 0).non_nulls())
    if (auto p = v.as_point()) spec.outline.push_back(*p);
  if (spec.outline.size() < 3)
    throw ComponentError("slab: outline needs at least 3 points, got " + std::to_string(spec.outline.size()));
  if (!std::all_of(spec.outline.begin(), spec.outline.end(), finite3))
    throw ComponentError("slab: outline has a non-finite (NaN/inf) coordinate");
  for (const auto& kv : input_or_empty(in, 1).branches()) {
    std::vector<Vec3> ring;
    for (const auto& o : kv.second)
      if (o)
        if (auto p = o->as_point()) ring.push_back(*p);
    if (ring.empty()) continue;
    if (ring.size() < 3)
      throw ComponentError("slab: a hole ring needs at least 3 points, got " + std::to_string(ring.size()));
    if (!std::all_of(ring.begin(), ring.end(), finite3))
      throw ComponentError("slab: a hole ring has a non-finite (NaN/inf) coordinate");
    spec.holes.push_back(std::move(ring));
  }
  const auto t = first_number(input_or_empty(in, 2));
  if (!t) throw ComponentError::missing("T");
  if (!(std::isfinite(*t) && *t > 0)) throw ComponentError("slab: thickness must be finite and > 0, got " + rs_num(*t));
  spec.thickness = *t;
  const MaterialSpec* mat = first_of<MaterialSpec>(input_or_empty(in, 3));
  if (!mat) throw ComponentError::missing("M");
  spec.material = *mat;
  spec.mesh_size = first_number(input_or_empty(in, 4)).value_or(0.5);
  if (!(std::isfinite(spec.mesh_size) && spec.mesh_size > 0))
    throw ComponentError("slab: mesh size must be finite and > 0, got " + rs_num(spec.mesh_size));
  if (auto lx = first_vector(input_or_empty(in, 5))) {
    if (!finite3(*lx)) throw ComponentError("slab: local-x axis override has a non-finite (NaN/inf) component");
    if (!(lx->norm() > 0)) throw ComponentError("slab: local-x axis override is zero-length (it carries no direction)");
    spec.local_x = *lx;
  }
  return {single(Value::of(std::move(spec)))};
}

ComponentInfo SurfaceLoadNode::info() const {
  return {"keta.surfaceload", "Surface Load", "Model", "Uniform pressure on a slab (N/m²)"};
}
std::vector<PinSpec> SurfaceLoadNode::inputs() const {
  return {PinSpec::item("S", "slab", "Target slab"),
          PinSpec::item("W", "vector", "Pressure vector (N/m²; LocalNormal uses w.z)"),
          PinSpec::item("C", "int", "Load-case id (default 1)").opt()};
}
std::vector<PinSpec> SurfaceLoadNode::outputs() const { return {PinSpec::item("L", "load", "The surface load")}; }
Outputs<Value> SurfaceLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  SurfacePressureLoad l;
  l.frame = parse_surface_frame(frame);
  l.outline = slab_outline_of(input_or_empty(in, 0));
  const auto w = first_vector(input_or_empty(in, 1));
  if (!w) throw ComponentError::missing("W");
  l.w = *w;
  LoadSpec spec;
  spec.kind = l;
  spec.case_id = first_int(input_or_empty(in, 2)).value_or(1);
  return {single(Value::of(std::move(spec)))};
}
json::Value SurfaceLoadNode::params() const {
  json::Value p = obj();
  p.set("frame", frame);
  return p;
}
void SurfaceLoadNode::set_params(const json::Value& p) {
  if (const auto* x = p.get("frame"); x && x->is_string()) frame = surface_frame_str(parse_surface_frame(x->as_string()));
}

ComponentInfo SurfaceThermalNode::info() const {
  return {"keta.surfacethermal", "Surface Thermal", "Model", "Thermal (ΔT + gradient) load on a slab"};
}
std::vector<PinSpec> SurfaceThermalNode::inputs() const {
  return {PinSpec::item("S", "slab", "Target slab"), PinSpec::item("DT", "number", "Temperature change ΔT (K)"),
          PinSpec::item("G", "number", "Gradient dT/dz (K/m; default 0)").opt(),
          PinSpec::item("C", "int", "Load-case id (default 1)").opt()};
}
std::vector<PinSpec> SurfaceThermalNode::outputs() const {
  return {PinSpec::item("L", "load", "The surface thermal load")};
}
Outputs<Value> SurfaceThermalNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  SurfaceThermalLoad l;
  l.outline = slab_outline_of(input_or_empty(in, 0));
  const auto dt = first_number(input_or_empty(in, 1));
  if (!dt) throw ComponentError::missing("DT");
  if (!std::isfinite(*dt)) throw ComponentError("surface thermal: ΔT must be finite, got " + rs_num(*dt));
  l.dt = *dt;
  l.grad = first_number(input_or_empty(in, 2)).value_or(0.0);
  if (!std::isfinite(l.grad))
    throw ComponentError("surface thermal: gradient dT/dz must be finite, got " + rs_num(l.grad));
  LoadSpec spec;
  spec.kind = l;
  spec.case_id = first_int(input_or_empty(in, 3)).value_or(1);
  return {single(Value::of(std::move(spec)))};
}

ComponentInfo SupportNode::info() const { return {"keta.support", "Support", "Model", "Block DOFs at a point"}; }
std::vector<PinSpec> SupportNode::inputs() const {
  return {PinSpec::item("P", "point", "Support location (point support)").opt(),
          PinSpec::list("K", "number",
                        "Elastic spring stiffness per DOF [kx, ky, kz, krx, kry, krz] (N/m, N.m/rad); "
                        "missing entries = free")
              .opt(),
          PinSpec::list("L", "line", "Support lines (line support along each segment)").opt()};
}
std::vector<PinSpec> SupportNode::outputs() const { return {PinSpec::item("S", "support", "The support")}; }
Outputs<Value> SupportNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  // Springs: up to six numbers read POSITIONALLY in DOF order; a missing or
  // null slot leaves that DOF free.
  std::array<std::optional<double>, 6> springs{};
  std::size_t i = 0;
  for (const auto& kv : input_or_empty(in, 1).branches())
    for (const auto& o : kv.second) {
      if (i >= 6) break;
      springs[i++] = o ? o->as_number() : std::nullopt;
    }
  const auto bl = blocked;
  Tree out;
  if (has_values(in[0]))
    out = map_items({{"P", &in[0]}}, [&](const std::vector<const Value*>& it) {
      SupportSpec s;
      s.point = point_of(*it[0], "P");
      s.blocked = bl;
      s.springs = springs;
      return Value::of(s);
    });
  for (const Line& l : beam_lines(input_or_empty(in, 2), /*accept_lines=*/true)) {
    SupportSpec s;
    s.point = l.a;
    s.blocked = bl;
    s.springs = springs;
    s.line = l;
    out.append(Path::single(0), Value::of(s));
  }
  return {out};
}
json::Value SupportNode::params() const {
  json::Value arr = json::Value::array();
  for (bool b : blocked) arr.push_back(json::Value(b));
  json::Value p = obj();
  p.set("blocked", std::move(arr));
  return p;
}
void SupportNode::set_params(const json::Value& p) {
  const auto* a = p.get("blocked");
  if (!a || !a->is_array()) return;
  const auto& arr = a->as_array();
  for (std::size_t i = 0; i < arr.size() && i < 6; ++i)
    if (arr[i].is_bool()) blocked[i] = arr[i].as_bool();
}

ComponentInfo PointMassNode::info() const {
  return {"keta.point_mass", "Point Mass", "Model", "Lumped nodal mass (modal + self-weight)"};
}
std::vector<PinSpec> PointMassNode::inputs() const {
  return {PinSpec::item("P", "point", "Mass location"), PinSpec::item("M", "number", "Translational mass (kg)").opt(),
          PinSpec::item("I", "vector", "Rotary inertia [Ixx, Iyy, Izz] (kg m^2)").opt()};
}
std::vector<PinSpec> PointMassNode::outputs() const { return {PinSpec::item("P", "pointmass", "The point mass")}; }
Outputs<Value> PointMassNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  double m = first_number(in[1]).value_or(mass);
  Vec3 inertia = first_vector(in[2]).value_or(Vec3::Zero());
  return {map_items({{"P", &in[0]}}, [&](const std::vector<const Value*>& it) {
    return Value::of(PointMassSpec{point_of(*it[0], "P"), m, inertia});
  })};
}
json::Value PointMassNode::params() const {
  json::Value p = obj();
  p.set("mass", mass);
  return p;
}
void PointMassNode::set_params(const json::Value& p) { mass = num_param(p, "mass", mass); }

ComponentInfo PointLoadNode::info() const { return {"keta.pointload", "Point Load", "Model", "Force at a point (N)"}; }
std::vector<PinSpec> PointLoadNode::inputs() const {
  return {PinSpec::item("P", "point", "Application point"), PinSpec::item("F", "vector", "Force vector (N)"),
          PinSpec::item("M", "vector", "Moment vector (N.m)").opt()};
}
std::vector<PinSpec> PointLoadNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> PointLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Vec3 M = first_vector(input_or_empty(in, 2)).value_or(Vec3::Zero());
  const Id c = case_id;
  return {map_items({{"P", &in[0]}, {"F", &in[1]}}, [&](const std::vector<const Value*>& it) {
    return Value::of(LoadSpec{NodalLoad{point_of(*it[0], "P"), vector_of(*it[1], "F"), M}, c});
  })};
}

ComponentInfo LineLoadNode::info() const { return {"keta.lineload", "Line Load", "Model", "Distributed load on beams (N/m)"}; }
std::vector<PinSpec> LineLoadNode::inputs() const {
  return {PinSpec::list("B", "beam", "Beams to load"), PinSpec::item("W", "vector", "Load per unit length (N/m)")};
}
std::vector<PinSpec> LineLoadNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> LineLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  auto lines = require_lines(in[0], "B", /*accept_lines=*/true);
  auto w = first_vector(in[1]);
  if (!w) throw ComponentError::missing("W");
  return {single(Value::of(LoadSpec{UdlLoad{std::move(lines), *w, parse_frame(frame)}, case_id}))};
}

ComponentInfo GravityNode::info() const { return {"keta.gravity", "Gravity", "Model", "Self-weight (rho*A*g on every bar)"}; }
std::vector<PinSpec> GravityNode::inputs() const { return {}; }
std::vector<PinSpec> GravityNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> GravityNode::solve(const SolveCtx&, const std::vector<Tree>&) const {
  return {single(Value::of(LoadSpec{GravityLoad{Vec3(0, 0, -kGravity * factor)}, case_id}))};
}
json::Value GravityNode::params() const {
  json::Value p = CaseNode::params();
  p.set("factor", factor);
  return p;
}
void GravityNode::set_params(const json::Value& p) {
  CaseNode::set_params(p);
  factor = num_param(p, "factor", factor);
}

ComponentInfo LinearLoadNode::info() const {
  return {"keta.linearload", "Linear Load", "Model", "Linearly-varying distributed load on beams (N/m)"};
}
std::vector<PinSpec> LinearLoadNode::inputs() const {
  return {PinSpec::list("B", "beam", "Beams to load"), PinSpec::item("W1", "vector", "Load per unit length at x1 (N/m)"),
          PinSpec::item("W2", "vector", "Load per unit length at x2 (N/m)"),
          PinSpec::item("x1", "number", "Range start (0..1 along the bar; default 0)").opt(),
          PinSpec::item("x2", "number", "Range end (0..1 along the bar; default 1)").opt()};
}
std::vector<PinSpec> LinearLoadNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> LinearLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  auto lines = require_lines(in[0], "B");
  auto w1 = first_vector(in[1]);
  if (!w1) throw ComponentError::missing("W1");
  auto w2 = first_vector(in[2]);
  if (!w2) throw ComponentError::missing("W2");
  double x1 = first_number(input_or_empty(in, 3)).value_or(0.0);
  double x2 = first_number(input_or_empty(in, 4)).value_or(1.0);
  if (!(x1 >= 0 && x2 <= 1 && x1 < x2))
    throw ComponentError("Linear Load: the range must satisfy 0 <= x1 < x2 <= 1");
  return {single(Value::of(LoadSpec{LinearLoad{std::move(lines), *w1, *w2, parse_frame(frame), x1, x2}, case_id}))};
}

ComponentInfo BarPointLoadNode::info() const {
  return {"keta.barpointload", "Point Load (Bar)", "Model", "Concentrated force on a bar interior (N)"};
}
std::vector<PinSpec> BarPointLoadNode::inputs() const {
  return {PinSpec::list("B", "beam", "Beams to load"), PinSpec::item("t", "number", "Position along the bar (0..1)"),
          PinSpec::item("F", "vector", "Force vector (N)")};
}
std::vector<PinSpec> BarPointLoadNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> BarPointLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  auto lines = require_lines(in[0], "B");
  auto t = first_number(in[1]);
  if (!t) throw ComponentError::missing("t");
  if (!(*t >= 0 && *t <= 1)) throw ComponentError("Point Load (Bar): t must lie in [0, 1]");
  auto f = first_vector(in[2]);
  if (!f) throw ComponentError::missing("F");
  return {single(Value::of(LoadSpec{BarPointLoad{std::move(lines), *t, *f, parse_frame(frame)}, case_id}))};
}

ComponentInfo ThermalLoadNode::info() const {
  return {"keta.thermalload", "Thermal Load", "Model", "Thermal load on beams (uniform dT + gradients)"};
}
std::vector<PinSpec> ThermalLoadNode::inputs() const {
  return {PinSpec::list("B", "beam", "Beams to load"), PinSpec::item("dT", "number", "Uniform temperature change (K)"),
          PinSpec::item("gy", "number", "Gradient dT/dy (K/m; default 0)").opt(),
          PinSpec::item("gz", "number", "Gradient dT/dz (K/m; default 0)").opt()};
}
std::vector<PinSpec> ThermalLoadNode::outputs() const { return {PinSpec::item("L", "load", "The load")}; }
Outputs<Value> ThermalLoadNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  auto lines = require_lines(in[0], "B");
  auto dt = first_number(in[1]);
  if (!dt) throw ComponentError::missing("dT");
  double gy = first_number(input_or_empty(in, 2)).value_or(0.0);
  double gz = first_number(input_or_empty(in, 3)).value_or(0.0);
  return {single(Value::of(LoadSpec{ThermalLoad{std::move(lines), *dt, gy, gz}, case_id}))};
}

ComponentInfo SettlementNode::info() const {
  return {"keta.settlement", "Settlement", "Model", "Imposed support displacement / rotation at a point (m, rad)"};
}
std::vector<PinSpec> SettlementNode::inputs() const {
  return {PinSpec::item("P", "point", "Support location (matched to the nearest node)"),
          PinSpec::item("D", "vector", "Imposed translation [dx, dy, dz] (m); nonzero components prescribe").opt(),
          PinSpec::item("R", "vector", "Imposed rotation [rx, ry, rz] (rad); nonzero components prescribe").opt()};
}
std::vector<PinSpec> SettlementNode::outputs() const {
  return {PinSpec::item("L", "load", "The settlement (imposed-displacement load)")};
}
Outputs<Value> SettlementNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const Tree& d = input_or_empty(in, 1);
  const Tree& r = input_or_empty(in, 2);
  if (!has_values(d) && !has_values(r))
    throw ComponentError("Settlement prescribes nothing: wire D (translation) and/or R (rotation)");
  std::vector<std::pair<std::string, const Tree*>> pins = {{"P", &in[0]}};
  int di = -1, ri = -1;
  if (has_values(d)) {
    pins.emplace_back("D", &d);
    di = static_cast<int>(pins.size()) - 1;
  }
  if (has_values(r)) {
    pins.emplace_back("R", &r);
    ri = static_cast<int>(pins.size()) - 1;
  }
  const Id c = case_id;
  return {map_items(pins, [&](const std::vector<const Value*>& it) {
    SettlementLoad s;
    s.point = point_of(*it[0], "P");
    auto fold = [&](int idx, int base, const char* pin) {
      if (idx < 0) return;
      Vec3 v = vector_of(*it[static_cast<std::size_t>(idx)], pin);
      for (int k = 0; k < 3; ++k) {
        if (!std::isfinite(v[k]))
          throw ComponentError(std::string("Settlement: ") + pin + " component " + std::to_string(k) + " is not finite");
        if (v[k] != 0.0) s.d[static_cast<std::size_t>(base + k)] = v[k];
      }
    };
    fold(di, 0, "D");
    fold(ri, 3, "R");
    bool any = false;
    for (const auto& x : s.d) any = any || x.has_value();
    if (!any) throw ComponentError("Settlement prescribes nothing: D and R are both all-zero");
    return Value::of(LoadSpec{s, c});
  })};
}

ComponentInfo FoundationNode::info() const {
  return {"keta.foundation", "Foundation", "Model", "Winkler elastic foundation (soil bed) on beams"};
}
std::vector<PinSpec> FoundationNode::inputs() const {
  return {PinSpec::list("B", "beam", "Beams to found"), PinSpec::item("KY", "number", "Bed modulus along local y (N/m per m)"),
          PinSpec::item("KZ", "number", "Bed modulus along local z (N/m per m; default 0)").opt(),
          PinSpec::item("KX", "number", "Bed modulus along local x (N/m per m; default 0)").opt(),
          PinSpec::item("KM", "number", "Torsional bed modulus about local x (N.m/m per rad; default 0)").opt()};
}
std::vector<PinSpec> FoundationNode::outputs() const { return {PinSpec::item("F", "foundation", "The foundation")}; }
Outputs<Value> FoundationNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  FoundationSpec f;
  f.lines = require_lines(in[0], "B");
  auto ky = first_number(in[1]);
  if (!ky) throw ComponentError::missing("KY");
  f.kfy = *ky;
  f.kfz = first_number(input_or_empty(in, 2)).value_or(0.0);
  f.kfx = first_number(input_or_empty(in, 3)).value_or(0.0);
  f.kmx = first_number(input_or_empty(in, 4)).value_or(0.0);
  for (auto [name, v] : {std::pair{"KX", f.kfx}, {"KY", f.kfy}, {"KZ", f.kfz}, {"KM", f.kmx}}) {
    if (!std::isfinite(v)) throw ComponentError(std::string("Foundation: ") + name + " modulus is not finite");
    if (v < 0) throw ComponentError(std::string("Foundation: ") + name + " modulus is negative; bed moduli must be >= 0");
  }
  return {single(Value::of(f))};
}

ComponentInfo CaseLabelNode::info() const {
  return {"keta.caselabel", "Case Label", "Model", "Name a load case (rides the L pin into Assemble)"};
}
std::vector<PinSpec> CaseLabelNode::inputs() const { return {PinSpec::item("N", "text", "Case name")}; }
std::vector<PinSpec> CaseLabelNode::outputs() const { return {PinSpec::item("L", "load", "The case label")}; }
Outputs<Value> CaseLabelNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  auto name = first_text(in[0]);
  if (!name) throw ComponentError::missing("N");
  return {single(Value::of(LoadSpec{CaseLabel{*name}, case_id}))};
}

}  // namespace keta::param

// Display units (see units.h).
#include "keta/param/units.h"

#include <array>
#include <cmath>
#include <utility>

#include "keta/core.h"

namespace keta::param {

namespace {

template <class E, std::size_t N>
using Names = std::array<std::pair<E, const char*>, N>;

const Names<ForceUnit, 3> kForce{{{ForceUnit::N, "N"}, {ForceUnit::Kn, "Kn"}, {ForceUnit::Mn, "Mn"}}};
const Names<MomentUnit, 3> kMoment{{{MomentUnit::Nm, "Nm"}, {MomentUnit::Knm, "Knm"}, {MomentUnit::Mnm, "Mnm"}}};
const Names<LengthUnit, 3> kLength{{{LengthUnit::M, "M"}, {LengthUnit::Mm, "Mm"}, {LengthUnit::Cm, "Cm"}}};
const Names<RotationUnit, 3> kRotation{
    {{RotationUnit::Mrad, "Mrad"}, {RotationUnit::Rad, "Rad"}, {RotationUnit::Deg, "Deg"}}};
const Names<StressUnit, 4> kStress{
    {{StressUnit::Pa, "Pa"}, {StressUnit::Kpa, "Kpa"}, {StressUnit::Mpa, "Mpa"}, {StressUnit::Gpa, "Gpa"}}};
const Names<LineForceUnit, 3> kLineForce{
    {{LineForceUnit::Npm, "Npm"}, {LineForceUnit::Knpm, "Knpm"}, {LineForceUnit::Mnpm, "Mnpm"}}};
const Names<LineMomentUnit, 3> kLineMoment{
    {{LineMomentUnit::Nmpm, "Nmpm"}, {LineMomentUnit::Knmpm, "Knmpm"}, {LineMomentUnit::Mnmpm, "Mnmpm"}}};

template <class E, std::size_t N>
const char* name_of(const Names<E, N>& t, E e) {
  for (const auto& [v, n] : t)
    if (v == e) return n;
  return "";
}

// Error wording: unknown variant `X`, expected one of `A`, `B`, `C`.
template <class E, std::size_t N>
void read(const json::Value& v, const char* key, const Names<E, N>& t, E& out) {
  const json::Value* x = v.get(key);
  if (!x) return;
  std::string list;
  for (const auto& [e, n] : t) list += (list.empty() ? "`" : ", `") + std::string(n) + "`";
  if (!x->is_string()) fail(std::string("units.") + key + ": expected a unit name (one of " + list + ")");
  for (const auto& [e, n] : t)
    if (x->as_string() == n) {
      out = e;
      return;
    }
  fail(std::string("units.") + key + ": unknown variant `" + x->as_string() + "`, expected one of " + list);
}

// Scale tables indexed by the enum value, in declaration order.
constexpr double k3[] = {1.0, 1e-3, 1e-6};

}  // namespace

double factor(ForceUnit u) { return k3[static_cast<int>(u)]; }
double factor(MomentUnit u) { return k3[static_cast<int>(u)]; }
double factor(LengthUnit u) { return u == LengthUnit::M ? 1.0 : u == LengthUnit::Mm ? 1e3 : 1e2; }
double factor(RotationUnit u) {
  return u == RotationUnit::Mrad ? 1e3 : u == RotationUnit::Rad ? 1.0 : 180.0 / M_PI;
}
double factor(StressUnit u) {
  constexpr double s[] = {1.0, 1e-3, 1e-6, 1e-9};
  return s[static_cast<int>(u)];
}
double factor(LineForceUnit u) { return k3[static_cast<int>(u)]; }
double factor(LineMomentUnit u) { return k3[static_cast<int>(u)]; }

const char* suffix(ForceUnit u) {
  constexpr const char* s[] = {"N", "kN", "MN"};
  return s[static_cast<int>(u)];
}
const char* suffix(MomentUnit u) {
  constexpr const char* s[] = {"N·m", "kN·m", "MN·m"};
  return s[static_cast<int>(u)];
}
const char* suffix(LengthUnit u) {
  constexpr const char* s[] = {"m", "mm", "cm"};
  return s[static_cast<int>(u)];
}
const char* suffix(RotationUnit u) {
  constexpr const char* s[] = {"mrad", "rad", "deg"};
  return s[static_cast<int>(u)];
}
const char* suffix(StressUnit u) {
  constexpr const char* s[] = {"Pa", "kPa", "MPa", "GPa"};
  return s[static_cast<int>(u)];
}
const char* suffix(LineForceUnit u) {
  constexpr const char* s[] = {"N/m", "kN/m", "MN/m"};
  return s[static_cast<int>(u)];
}
const char* suffix(LineMomentUnit u) {
  constexpr const char* s[] = {"N·m/m", "kN·m/m", "MN·m/m"};
  return s[static_cast<int>(u)];
}

int decimals(ForceUnit u) { return std::array{0, 2, 4}[static_cast<std::size_t>(u)]; }
int decimals(MomentUnit u) { return std::array{0, 2, 4}[static_cast<std::size_t>(u)]; }
int decimals(LengthUnit u) { return std::array{4, 2, 3}[static_cast<std::size_t>(u)]; }
int decimals(RotationUnit u) { return std::array{3, 5, 3}[static_cast<std::size_t>(u)]; }
int decimals(StressUnit u) { return std::array{0, 1, 1, 3}[static_cast<std::size_t>(u)]; }
int decimals(LineForceUnit u) { return std::array{0, 2, 4}[static_cast<std::size_t>(u)]; }
int decimals(LineMomentUnit u) { return std::array{0, 2, 4}[static_cast<std::size_t>(u)]; }
int glyph_decimals(ForceUnit u) { return std::array{0, 0, 3}[static_cast<std::size_t>(u)]; }
int glyph_decimals(MomentUnit u) { return std::array{0, 1, 3}[static_cast<std::size_t>(u)]; }
int glyph_decimals(LengthUnit u) { return std::array{3, 0, 1}[static_cast<std::size_t>(u)]; }
int glyph_decimals(StressUnit u) { return std::array{0, 1, 1, 3}[static_cast<std::size_t>(u)]; }

json::Value DisplayUnits::to_json() const {
  json::Value u = json::Value::object();
  u.set("force", name_of(kForce, force));
  u.set("moment", name_of(kMoment, moment));
  u.set("length", name_of(kLength, length));
  u.set("rotation", name_of(kRotation, rotation));
  u.set("stress", name_of(kStress, stress));
  u.set("line_force", name_of(kLineForce, line_force));
  u.set("line_moment", name_of(kLineMoment, line_moment));
  return u;
}

DisplayUnits DisplayUnits::from_json(const json::Value& v) {
  DisplayUnits d;
  if (!v.is_object()) fail("units: expected an object");
  read(v, "force", kForce, d.force);
  read(v, "moment", kMoment, d.moment);
  read(v, "length", kLength, d.length);
  read(v, "rotation", kRotation, d.rotation);
  read(v, "stress", kStress, d.stress);
  read(v, "line_force", kLineForce, d.line_force);
  read(v, "line_moment", kLineMoment, d.line_moment);
  return d;
}

}  // namespace keta::param

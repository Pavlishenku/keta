// keta::param — project display units.
//
// The core is SI everywhere (m, N, Pa, rad); a document chooses how values are
// SHOWN. Each family maps an SI value to its display value by multiplication
// (`display = si * factor()`), with a suffix and the decimal counts used for
// tables and for viewport glyphs.
//
// Persisted form (the `units` block of a document): the variant NAMES,
// spelled as the enumerators below —
//   {"force": "Kn", "moment": "Knm", "length": "Mm", "rotation": "Mrad",
//    "stress": "Mpa", "line_force": "Knpm", "line_moment": "Knmpm"}
// which is also the frozen default set. A missing family takes its default;
// an unknown spelling is a typed parse error.
#pragma once
#include <string>

#include "keta/json.h"

namespace keta::param {

enum class ForceUnit { N, Kn, Mn };
enum class MomentUnit { Nm, Knm, Mnm };
enum class LengthUnit { M, Mm, Cm };
enum class RotationUnit { Mrad, Rad, Deg };
enum class StressUnit { Pa, Kpa, Mpa, Gpa };
enum class LineForceUnit { Npm, Knpm, Mnpm };       // shell resultant forces
enum class LineMomentUnit { Nmpm, Knmpm, Mnmpm };   // shell resultant moments

// SI -> display multiplier, display suffix, table decimals, glyph decimals
// (rotation and the line families have no glyph form).
double factor(ForceUnit);
double factor(MomentUnit);
double factor(LengthUnit);
double factor(RotationUnit);
double factor(StressUnit);
double factor(LineForceUnit);
double factor(LineMomentUnit);
const char* suffix(ForceUnit);
const char* suffix(MomentUnit);
const char* suffix(LengthUnit);
const char* suffix(RotationUnit);
const char* suffix(StressUnit);
const char* suffix(LineForceUnit);
const char* suffix(LineMomentUnit);
int decimals(ForceUnit);
int decimals(MomentUnit);
int decimals(LengthUnit);
int decimals(RotationUnit);
int decimals(StressUnit);
int decimals(LineForceUnit);
int decimals(LineMomentUnit);
int glyph_decimals(ForceUnit);
int glyph_decimals(MomentUnit);
int glyph_decimals(LengthUnit);
int glyph_decimals(StressUnit);

struct DisplayUnits {
  ForceUnit force = ForceUnit::Kn;
  MomentUnit moment = MomentUnit::Knm;
  LengthUnit length = LengthUnit::Mm;
  RotationUnit rotation = RotationUnit::Mrad;
  StressUnit stress = StressUnit::Mpa;
  LineForceUnit line_force = LineForceUnit::Knpm;
  LineMomentUnit line_moment = LineMomentUnit::Knmpm;
  friend bool operator==(const DisplayUnits&, const DisplayUnits&) = default;

  json::Value to_json() const;
  // Missing families take their defaults; throws keta::Error naming the
  // family and the accepted spellings on an unknown one.
  static DisplayUnits from_json(const json::Value& v);
};

}  // namespace keta::param

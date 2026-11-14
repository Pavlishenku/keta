// Export to script: a solved keta Model as an executable Python program
// driving `import keta` (the emitted calls are the builder of src/keta/py, so
// the script also documents the model in keta's scripting dialect).
//
// The export is MODEL-LEVEL: a frozen, diffable snapshot of the model as
// currently parameterized (a parametric document's sliders and rules baked
// into concrete nodes, bars and loads). With results supplied, an assertion
// block pins — per case and per combination — the reaction sums, the peak
// nodal translation and the two spot member forces (global max |My| and
// |N|), each to 1e-6 of that case's own scale, so the script PROVES it
// re-solves to matching numbers by running.
//
// Deterministic: a pure function of (model, results, options), byte-identical
// for identical inputs (model vectors in order, no timestamp, shortest
// round-trip floats). Pure ASCII (non-ASCII text is escaped).
//
// Blocked features (scan the whole model, fail once listing every blocker,
// never a silent drop): what the Python builder cannot express faithfully —
// cables (linear analyze() rejects them), tendons, slabs / shells / surface
// loads, tapered bars, prestrains, nodal / bar mass records, rotated /
// unilateral / coupled supports, releases other than ry / rz, unequal end
// eccentricities, span moments and non-finite numbers. Partial-span uniform
// loads are emitted as constant linear loads over their relative range.
#pragma once
#include <stdexcept>
#include <string>
#include <vector>

#include "keta/statics.h"

namespace keta {

struct ScriptOptions {
  // Stations per bar passed to m.analyze(stations=...): must match the
  // results' recovery so the asserted station indices line up.
  int stations_per_bar = 11;
};

class ScriptError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
  std::vector<std::string> features;  // the blockers (empty for a station mismatch)
};

// Throws ScriptError; its messages read:
//   "cannot export to Python script: N unsupported model feature(s)\n  - ...\n"
//   "cannot export to Python script: results station-count mismatch — ..."
std::string export_script(const Model& model, const StaticResults* results = nullptr,
                          const ScriptOptions& opt = {});

}  // namespace keta

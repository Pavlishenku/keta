// Model + results JSON serialization (own dependency-free JSON DOM).
#pragma once
#include "keta/json.h"
#include "keta/statics.h"

namespace keta {

json::Value model_to_json(const Model& model);
Model model_from_json(const json::Value& v);

// ---------------------------------------------------------------------------
// Result store: a self-describing JSON document of the static
// results that reads back bit for bit (numbers are written with 17
// significant digits, which round-trip a double exactly).
// Layout: schema_version, provenance {generator, note}, sign_convention,
// cases[], combinations[], auto_fixed[]; a case carries id, name,
// non_superposable, node_disp (6 per node), reactions, bars (end_forces: 12,
// stations: rows x N Vy Vz T My Mz u v w rx ry rz, soil), shells, links,
// constraint_forces, skew_reactions.
// ---------------------------------------------------------------------------
inline constexpr int kResultsSchemaVersion = 1;
// The convention of the stored forces, stamped in every document.
extern const char* const kResultsSignConvention;

struct ResultsProvenance {
  std::string generator;  // producing tool and version; empty = "keta <version>"
  std::string note;       // the caller's note: model, date, command line
};
struct StoredResults {
  int schema_version = kResultsSchemaVersion;
  ResultsProvenance provenance;
  std::string sign_convention;
  StaticResults results;
};
json::Value results_to_json(const StaticResults& results, const ResultsProvenance& provenance = {});
// Typed errors: a document without schema_version (not a result store), a
// schema newer than kResultsSchemaVersion, a missing or malformed record.
StoredResults results_from_json(const json::Value& v);

// File helpers.
Model load_model_file(const std::string& path);
StoredResults load_results_file(const std::string& path);
void save_json_file(const std::string& path, const json::Value& v, int indent = 1);

}  // namespace keta

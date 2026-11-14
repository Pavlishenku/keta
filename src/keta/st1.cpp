// import_st1 is a thin facade over the full ST1 session (interpreter +
// runtime, st1_lang.h / st1_run.h): the deck is executed — pseudo-programming
// resolved, blocks dispatched, EXEC statements actually solved — and the
// resulting model snapshot plus the run metadata are returned. Every command
// goes through the one dispatch used by the CLI session path.
//
// LIRE resolution: include paths are tried verbatim first, then relative to
// the including deck's directory (multi-file projects à la OA108).
#include "keta/st1.h"

#include <fstream>
#include <sstream>

#include "keta/st1_run.h"

namespace keta {

namespace {

std::string read_file(const std::string& path) {
  std::ifstream f(path);
  require(f.good(), "cannot open '" + path + "'");
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

St1Import run_deck(const std::string& text, const std::string& base_dir,
                   const St1Options& opt) {
  st1::St1Session session(opt);
  if (!base_dir.empty()) {
    session.interp->set_file_loader([base_dir](const std::string& p) {
      std::ifstream direct(p);
      if (direct.good()) {
        std::stringstream ss;
        ss << direct.rdbuf();
        return ss.str();
      }
      return read_file(base_dir + "/" + p);
    });
  }
  session.run(text);
  st1::St1Runtime& rt = *session.runtime;
  St1Import out;
  out.model = rt.snapshot();
  // Import compatibility: the returned model carries every defined loading
  // (the runtime's snapshot only holds the last EXEC set) and the COMB
  // definitions as model combinations.
  out.model.cases.clear();
  for (Id cid : rt.charg_order_) {
    const st1::ChargDef& c = rt.chargs_.at(cid);
    if (c.is_comb) continue;
    LoadCase lc;
    lc.id = c.id;
    lc.name = c.name;
    lc.loads = rt.lower_loads(c.loads);
    out.model.cases.push_back(lc);
  }
  for (const auto& [id, cb] : rt.combs_) {
    Combination k;
    k.id = id;
    k.name = cb.name;
    for (const auto& t : cb.terms)
      if (t.kind == "CHARG") k.terms.emplace_back(t.id, t.coef);
    out.model.combos.push_back(k);
  }
  out.model.invalidate_index();
  out.title = rt.title();
  out.warnings = rt.warnings();
  out.exec_requested = rt.exec_charg_requested;
  out.exec_modal = rt.exec_modal_requested;
  out.exec_buckling = rt.exec_buckling_requested;
  out.etude_pas = rt.etude_pas_hint;
  return out;
}

}  // namespace

St1Import import_st1(const std::string& text, const St1Options& opt) {
  return run_deck(text, "", opt);
}

St1Import import_st1_file(const std::string& path, const St1Options& opt) {
  std::string base;
  auto slash = path.find_last_of("/\\");
  if (slash != std::string::npos) base = path.substr(0, slash);
  return run_deck(read_file(path), base, opt);
}

}  // namespace keta

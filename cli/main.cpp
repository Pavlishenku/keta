// keta CLI: solve a JSON model, import an ST1 deck, run modal/buckling or
// diagnostics. Results go to stdout (summary) and optionally to a JSON file.
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>

#include "keta/diagnostics.h"
#include "keta/eigsolve.h"
#include "keta/io.h"
#include "keta/param/document.h"
#include "keta/script.h"
#include "keta/st1.h"
#include "keta/st1_run.h"

using namespace keta;

namespace {

int usage() {
  std::printf(
      "keta 1.0.0 — bar-structure analysis kernel\n"
      "usage:\n"
      "  keta run <model.json> [--out results.json] [--stations N]\n"
      "  keta st1 <deck.st1> [--gravity G] [--out model.json] [--run]\n"
      "  keta session <deck.st1> [--gravity G] [--edition st1|keta] [--pas-conv]\n"
      "                           [--fluage-st1] [--pertes-st1] [--threads N]\n"
      "                           [--erreur FILE] [--progress]\n"
      "                           (full ST1 session: the deck drives\n"
      "                            EXEC/RESU/SORTIE/GET/graphics)\n"
      "      --gravity G: POIDS PROPRE weighs G x RO x A (ST1: RO is a weight\n"
      "                   density, G = 1, the default)\n"
      "      --edition: RESU tables as ST1 prints them (st1, the default) or\n"
      "                 keta's own tables (keta)\n"
      "      --pas-conv: SURCH convoys moved on the PAS_CONV grid only, as ST1\n"
      "                  (default: grid + exact influence-line breakpoints)\n"
      "      --fluage-st1: PHASAGE creep by ST1's equivalent-time method\n"
      "                  (default: the exact linear superposition)\n"
      "      --pertes-st1: PHASAGE deferred tendon losses by EN 1992-1-1\n"
      "                  eq. 5.46 on each interval, as ST1 (default: solved\n"
      "                  by the structure)\n"
      "      --threads N: workers of the PHASAGE station recovery (default 0 =\n"
      "                  all cores; results identical for any N)\n"
      "      --erreur FILE: the run's log and messages, located by file and\n"
      "                  line, as ST1's erreur.txt (default erreur.txt in the\n"
      "                  working directory; \"\" writes none)\n"
      "      --progress: report each EXEC/RESU/SORTIE on stderr as it starts\n"
      "                  (and each EXEC as it ends): `progress: ligne N : ...`\n"
      "  keta modal <model.json> [--modes N]\n"
      "  keta buckling <model.json> [--modes N]\n"
      "  keta diagnose <model.json>\n"
      "  keta script <model.json | document.json> [--out script.py] [--stations N]\n"
      "            (export a self-checking Python script of a model or of a parametric document)\n");
  return 2;
}

const char* opt_value(int argc, char** argv, const char* name) {
  for (int i = 0; i + 1 < argc; ++i)
    if (std::strcmp(argv[i], name) == 0) return argv[i + 1];
  return nullptr;
}
bool has_flag(int argc, char** argv, const char* name) {
  for (int i = 0; i < argc; ++i)
    if (std::strcmp(argv[i], name) == 0) return true;
  return false;
}

void print_static_summary(const Model& m, const StaticResults& r) {
  for (const CaseResult& c : r.cases) {
    double wmax = 0;
    for (const Vec6& d : c.node_disp) wmax = std::max(wmax, d.segment<3>(0).norm());
    Vec3 rsum = Vec3::Zero();
    for (const Reaction& re : c.reactions) rsum += re.R.segment<3>(0);
    std::printf("case %lld '%s': |u|max = %.6g m, sum R = (%.6g, %.6g, %.6g)\n",
                static_cast<long long>(c.id), c.name.c_str(), wmax, rsum.x(), rsum.y(), rsum.z());
  }
  (void)m;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) return usage();
  std::string cmd = argv[1];
  std::string path = argv[2];
  try {
    if (cmd == "run") {
      Model m = load_model_file(path);
      StaticOptions so;
      if (const char* s = opt_value(argc, argv, "--stations")) so.stations_per_bar = std::atoi(s);
      StaticResults r = solve_static(m, so);
      print_static_summary(m, r);
      if (const char* o = opt_value(argc, argv, "--out")) {
        save_json_file(o, results_to_json(r, {"", "keta run " + path}));
        std::printf("results written to %s\n", o);
      }
    } else if (cmd == "st1") {
      St1Options so;
      if (const char* g = opt_value(argc, argv, "--gravity")) so.gravity = std::atof(g);
      St1Import imp = import_st1_file(path, so);
      std::printf("ST1 '%s': %zu nodes, %zu bars, %zu cases\n", imp.title.c_str(),
                  imp.model.nodes.size(), imp.model.bars.size(), imp.model.cases.size());
      for (const std::string& w : imp.warnings) std::printf("  warning: %s\n", w.c_str());
      if (const char* o = opt_value(argc, argv, "--out")) {
        save_json_file(o, model_to_json(imp.model));
        std::printf("model written to %s\n", o);
      }
      if (has_flag(argc, argv, "--run") || imp.exec_requested) {
        StaticResults r = solve_static(imp.model);
        print_static_summary(imp.model, r);
      }
    } else if (cmd == "session") {
      St1Options so;
      if (const char* g = opt_value(argc, argv, "--gravity")) so.gravity = std::atof(g);
      so.st1_edition = true;
      if (const char* e = opt_value(argc, argv, "--edition")) {
        require(std::strcmp(e, "st1") == 0 || std::strcmp(e, "keta") == 0,
                "--edition: st1 or keta");
        so.st1_edition = std::strcmp(e, "st1") == 0;
      }
      so.pas_conv_grid_only = has_flag(argc, argv, "--pas-conv");
      so.creep_superposition = !has_flag(argc, argv, "--fluage-st1");
      so.exact_losses = !has_flag(argc, argv, "--pertes-st1");
      if (const char* t = opt_value(argc, argv, "--threads")) so.threads = std::atoi(t);
      std::ifstream f(path);
      require(f.good(), "cannot open '" + path + "'");
      std::stringstream ss;
      ss << f.rdbuf();
      st1::St1Session session(so);
      {  // LIRE resolves relative to the deck's directory too
        std::string base;
        auto slash = path.find_last_of("/\\");
        if (slash != std::string::npos) base = path.substr(0, slash);
        if (!base.empty())
          session.interp->set_file_loader([base](const std::string& p) {
            std::ifstream direct(p);
            std::stringstream buf;
            if (direct.good()) { buf << direct.rdbuf(); return buf.str(); }
            std::ifstream rel(base + "/" + p);
            require(rel.good(), "cannot open '" + p + "'");
            buf << rel.rdbuf();
            return buf.str();
          });
      }
      // Console SORTIE lines are printed as they are written (a long run
      // shows its progress); file SORTIE goes to its files.
      session.runtime->output().on_console = [](const std::string& l) {
        std::fputs(l.c_str(), stdout);
        std::fputc('\n', stdout);
        std::fflush(stdout);
      };
      // erreur.txt, as ST1 writes it (default: the working directory).
      const char* erreur_opt = opt_value(argc, argv, "--erreur");
      const std::string erreur_path = erreur_opt ? erreur_opt : "erreur.txt";
      std::ofstream erreur;
      if (!erreur_path.empty()) {
        erreur.open(erreur_path, std::ios::binary | std::ios::trunc);
        if (!erreur) std::fprintf(stderr, "  warning: cannot write '%s'\n", erreur_path.c_str());
      }
      st1::ErreurLog log([&erreur](const std::string& text) {
        if (erreur.is_open()) erreur << text << std::flush;
      });
      log.header(path);
      const bool progress = has_flag(argc, argv, "--progress");
      std::size_t warnings_logged = 0;
      auto log_warnings = [&]() {
        const auto& w = session.runtime->warnings();
        for (; warnings_logged < w.size(); ++warnings_logged)
          log.message(warnings_logged < session.warning_where.size()
                          ? session.warning_where[warnings_logged]
                          : st1::Interp::Where{path, 0},
                      false, w[warnings_logged]);
      };
      using clock = std::chrono::steady_clock;
      const auto t0 = clock::now();
      auto t_step = t0;
      session.on_statement = [&](const st1::Statement& st, const st1::Interp::Where& at,
                                 bool done) {
        std::string head = st.empty() || st[0].str ? "" : st[0].s;
        for (char& c : head) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        const bool step = head == "EXEC" || head == "RESU" || head == "SORTIE";
        if (!step) {
          if (done) log_warnings();
          return;
        }
        // The command as the deck reads it, keywords upper-cased, `1,2,3`.
        std::string what;
        for (const st1::Tok& t : st) {
          if (t.s == ",") {
            what += ',';
            continue;
          }
          if (!what.empty() && what.back() != ',') what += ' ';
          if (t.str) {
            what += "'" + t.s + "'";
          } else {
            for (char c : t.s) what += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
          }
        }
        if (what.size() > 80) what = what.substr(0, 77) + "...";
        const std::string label = head == "EXEC"   ? "Execution " + what
                                  : head == "RESU" ? "Ecriture des resultats (RESU)"
                                                   : what;
        if (!done) {
          t_step = clock::now();
          log.step(label);
          if (progress) {
            std::fprintf(stderr, "progress: ligne %d : %s\n", at.line, label.c_str());
            std::fflush(stderr);
          }
          return;
        }
        log_warnings();
        if (head != "EXEC") return;
        const double dt = std::chrono::duration<double>(clock::now() - t_step).count();
        log.elapsed(dt);
        if (progress) {
          std::fprintf(stderr, "progress: ligne %d : %s terminee en %.2f s\n", at.line,
                       label.c_str(), dt);
          std::fflush(stderr);
        }
      };
      try {
        session.run(ss.str(), path);
      } catch (const std::exception& e) {
        log_warnings();
        const auto& w = session.interp->error_where();
        const st1::Interp::Where at = w.empty() ? st1::Interp::Where{path, 0} : w[0];
        log.message(at, true, e.what());
        for (const std::string& m : session.runtime->warnings())
          std::printf("  warning: %s\n", m.c_str());
        std::fflush(stdout);
        // The error names its line: `error: ligne 12 : ST1: ...`.
        if (at.line > 0)
          throw std::runtime_error("ligne " + std::to_string(at.line) +
                                   (at.file != path ? " de " + at.file : std::string()) +
                                   " : " + e.what());
        throw;
      }
      log_warnings();
      log.elapsed(std::chrono::duration<double>(clock::now() - t0).count());
      for (const std::string& w : session.runtime->warnings())
        std::printf("  warning: %s\n", w.c_str());
      if (progress) std::fprintf(stderr, "progress: fin\n");
    } else if (cmd == "modal") {
      Model m = load_model_file(path);
      ModalOptions mo;
      if (const char* n = opt_value(argc, argv, "--modes")) mo.n_modes = std::atoi(n);
      ModalResults r = solve_modal(m, mo);
      for (std::size_t i = 0; i < r.frequencies_hz.size(); ++i)
        std::printf("mode %zu: f = %.6g Hz (Meff z = %.4g)\n", i + 1, r.frequencies_hz[i],
                    r.effective_mass[i][2]);
    } else if (cmd == "buckling") {
      Model m = load_model_file(path);
      BucklingOptions bo;
      if (const char* n = opt_value(argc, argv, "--modes")) bo.n_modes = std::atoi(n);
      BucklingResults r = solve_buckling(m, bo);
      if (r.load_factors.empty())
        std::printf("no positive buckling factor (all-tension or zero axial state)\n");
      for (std::size_t i = 0; i < r.load_factors.size(); ++i)
        std::printf("mode %zu: lambda = %.6g\n", i + 1, r.load_factors[i]);
    } else if (cmd == "script") {
      // A keta model JSON is solved here; a parametric document is evaluated and
      // its first Analyze node's model + results are exported (the coherent
      // pair the document solved to), else its first Assemble node's model
      // with linear statics solved here.
      std::ifstream in(path, std::ios::binary);
      if (!in) throw std::runtime_error("cannot open " + path);
      std::stringstream ss;
      ss << in.rdbuf();
      const json::Value root = json::Value::parse(ss.str());
      ScriptOptions so;
      std::string text;
      if (root.has("nodes") && root.has("wires")) {
        param::Registry reg = param::standard_registry();
        param::LoadedDocument d = param::load_document(ss.str(), reg);
        d.graph.evaluate(param::SolveCtx{});
        const param::AnalyzedModel* am = nullptr;
        for (param::NodeId id : d.graph.node_ids()) {
          if (d.graph.node_ref(id).component->info().type_id != "keta.analyze") continue;
          if (const auto* outs = d.graph.outputs(id); outs && !outs->empty() && !(*outs)[0].empty())
            if (const auto* r = (*outs)[0].first()->get<param::ResultsRef>()) am = r->get();
          if (am) break;
        }
        if (am) {
          for (const CaseResult& c : am->results.cases)
            if (!c.bars.empty()) {
              so.stations_per_bar = static_cast<int>(c.bars[0].stations.size());
              break;
            }
          text = export_script(*am->model, &am->results, so);
        } else {
          // No Analyze node (a second-order / nonlinear / dynamic document):
          // export the first assembled model, linear statics solved here.
          const Model* model = nullptr;
          for (param::NodeId id : d.graph.node_ids()) {
            if (d.graph.node_ref(id).component->info().type_id != "keta.assemble") continue;
            if (const auto* outs = d.graph.outputs(id); outs && !outs->empty() && !(*outs)[0].empty())
              if (const auto* mr = (*outs)[0].first()->get<param::ModelRef>()) model = mr->get();
            if (model) break;
          }
          if (!model) throw std::runtime_error("the document has no solved Analyze or Assemble node");
          if (const char* st = opt_value(argc, argv, "--stations")) so.stations_per_bar = std::atoi(st);
          (void)export_script(*model, nullptr, so);  // an unsupported feature refuses before the solve
          StaticOptions o;
          o.stations_per_bar = so.stations_per_bar;
          const StaticResults r = solve_static(*model, o);
          text = export_script(*model, &r, so);
        }
      } else {
        Model m = model_from_json(root);
        if (const char* st = opt_value(argc, argv, "--stations")) so.stations_per_bar = std::atoi(st);
        (void)export_script(m, nullptr, so);  // an unsupported feature refuses before the solve
        StaticOptions o;
        o.stations_per_bar = so.stations_per_bar;
        const StaticResults r = solve_static(m, o);
        text = export_script(m, &r, so);
      }
      if (const char* o = opt_value(argc, argv, "--out")) {
        std::ofstream out(o, std::ios::binary);
        out << text;
      } else {
        std::fwrite(text.data(), 1, text.size(), stdout);
      }
    } else if (cmd == "diagnose") {
      Model m = load_model_file(path);
      auto diags = diagnose(m);
      if (diags.empty()) std::printf("no findings\n");
      for (const auto& d : diags)
        std::printf("[%s] %s\n        -> %s\n",
                    d.severity == Severity::Error ? "ERROR" : "warn", d.message.c_str(),
                    d.suggestion.c_str());
      return diags.empty() ? 0 : 1;
    } else {
      return usage();
    }
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error: %s\n", e.what());
    return 1;
  }
  return 0;
}

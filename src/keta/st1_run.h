// ST1 runtime: executes the full ST1 v24 command language over the keta
// engine. The interpreter (st1_lang.h) resolves the pseudo-programming layer
// (variables, expressions, POUR/SI, LIRE) and hands engine statements to
// St1Runtime::execute() one at a time, in program order — so an EXEC solves
// immediately and a following GET/SI sees its results.
//
// State model: ST1 descriptions are re-definable at any time (the LAST
// description wins, manual p.42); the runtime therefore keeps raw builder
// maps and materializes a keta Model snapshot lazily whenever an EXEC needs
// one. Results are stored per kind (CHARG, COMB, ETAT, PHASAGE, SPECTRE,
// INCENDIE, ENV, SURCH, modal, buckling, historique) and are the source for
// GET, RESU, COMB, ENV and the drawing commands.
//
// Conventions: ST1 units are user-consistent (gravity option converts weight
// densities); ST1 lists efforts with N POSITIVE IN COMPRESSION (manual p.60)
// while keta is tension-positive — the flip happens at the GET/RESU/listing
// boundary only, never inside the engine.
#pragma once
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>

#include "keta/eig_ext.h"
#include "keta/eigsolve.h"
#include "keta/etude.h"
#include "keta/model.h"
#include "keta/soilzones.h"
#include "keta/springmass.h"
#include "keta/st1.h"
#include "keta/st1_lang.h"
#include "keta/statics.h"
#include "keta/surch.h"
#include "keta/tapered.h"

namespace keta::st1 {

// One ST1 bar-load line after its bar list (manual p.110-115):
// <UNI, CON XL x, LIN XL x1 x2> with the forces / couples FX..MZ (one value,
// or two for LIN) and the modifiers GLO LOC PRO REL ABS anywhere on the line
// (the manual writes the frame AFTER the values). Appends one load per bar;
// false when the kind is not a bar load. A word it cannot read stops the
// line; it is returned in `ignored` (when given) for the caller to warn.
bool st1_bar_load_line(const Statement& st, std::size_t i, const std::vector<Id>& ids,
                       std::vector<Load>& out, std::string* ignored = nullptr);

// ---------------------------------------------------------------- builder --
struct CaraProps {
  double sx = 0, sy = 0, sz = 0, ix = 0, iy = 0, iz = 0, ipol = 0;
  double vy = 0, wy = 0, vz = 0, wz = 0;
  double exty = 0, extz = 0;  // torsion-centre offsets (stored; stresses only)
  bool operator<(const CaraProps& o) const;
  bool operator==(const CaraProps& o) const;
};
struct ConsProps {
  double E = 0, nu = 0.3, ro = 0, alpha = 0, g_dyn = 0;
  int amortissement = 0;  // damping-set id (0 = none)
  int mat = 0;            // MAT reference (0 = direct constants)
  int mat_origin = 0;     // MAT materialized into these constants (keeps two MATs
                          // of equal E apart: their time-dependent laws differ)
  bool operator<(const ConsProps& o) const;
};

// ---------------------------------------------------- CARA VAR / PSE / SPECIAL --

// `CARA VAR <LIN Y, LIN Z, PARA> liste` + its OR / MI / EX control lines
// (manual p.59-64). The characteristics are kept in ST1 form and lowered to a
// `keta::TaperedVar` (SectionProps: A<-SX, Avy<-SY, Avz<-SZ, J<-IX, Iy<-IY,
// Iz<-IZ) when the model snapshot is built.
struct CaraVarDef {
  TaperLaw law = TaperLaw::LinY;
  CaraProps orig, mid, ext;                 // OR / MI (PARA only) / EX
  bool has_or = false, has_mi = false, has_ex = false;
};

// `CARA SPECIAL j` (manual p.210): a classic R.D.M. bar whose force components
// are capped. `cara` is the RIGIDITE CARA line, `lim` the EFFORT line.
struct CaraSpecialDef {
  CaraProps cara;
  SpecialLimits lim;
  bool has_cara = false, has_effort = false;
};

// ------------------------------------------------------------ MASSE PROPRE --
// `MASSE PROPRE ... FIN` (manual p.69-71) lowered onto the model's
// `NodalMass` / `BarMass` tables. `lumped` records ST1's `MASSE DI` on the
// AUTOMATIC bar masses; keta's lumped/consistent switch is whole-model
// (`ModalOptions::lumped_mass`), so a deck mixing MASSE and MASSE DI takes the
// diagonal form for every bar and is warned about — documented deviation.
struct MasseDef {
  std::vector<NodalMass> nodal;
  std::vector<BarMass> bars;
  bool lumped = false, consistent = false;
  // `BARRE liste` without a MASSE form: ST1's default automatic masses are
  // diagonal, translations only (p.69: "diagonales et definies en deplacement
  // seulement par defaut"); an explicit `MASSE` (consistent) elsewhere wins.
  bool lumped_default = false;
  // a MASSE PROPRE block was given: only its listed bars get automatic masses
  bool selective = false;
  bool diagonal() const { return lumped || (lumped_default && !consistent); }
  bool any() const {
    return !nodal.empty() || !bars.empty() || lumped || consistent || lumped_default;
  }
};

// ------------------------------------------------------------- ENV / ENV COMB --

// One `ENV i (COMB) ... FIN` block (manual p.231 / p.234).
struct EnvDef {
  Id id = 0;
  std::string name;
  bool comb = false;  // ENV i COMB: sources carry the [cfmin, cfmax] box
  // Privileged component (absent = plain field-by-field envelope, p.232 §13.3.2).
  bool has_priv = false;
  EnvWhat what = EnvWhat::Effort;
  int comp = 0;
  ContrFiber fiber = ContrFiber::VY;
  Id contrg = 0;  // ENV ... CONTRG ic
  // Optional REAC concomitance (p.231).
  bool has_reac = false;
  Id reac_appui = 0;
  int reac_comp = 0;
  struct Source {
    std::string kind;  // CHARG / SURCH / COMB / ETAT / PHASAGE / SPECTRE / HISTORIQUE / ENV
    Id id = 0;
    double cmin = 1.0, cmax = 1.0;
  };
  std::vector<Source> sources;
};

// The stored result of one ENV: the privileged run section by section, plus the
// plain field-by-field bounds that feed ENV-of-ENV, RESU and GET on the
// non-privileged components.
struct EnvResult {
  Id id = 0;
  std::string name;
  bool has_priv = false;
  EnvWhat what = EnvWhat::Effort;
  Id contrg = 0;  // the privileged CONTRG (what == Contrg)
  ContrFiber fiber = ContrFiber::VY;  // the privileged fibre (what == Contr)
  // reaction bounds per APPUI (lo, hi) over the sources, each source with
  // the APPUI it was computed with (absent = 0); the union is listed
  std::map<Id, std::pair<Vec6, Vec6>> appui_bounds;
  int comp = 0;
  std::vector<EnvPoint> points;  // one per (bar, study section)
  Envelope bounds;               // lo / hi pseudo-cases
  EnvRun run;                    // privileged run at the node / reaction, when declared
  bool has_run = false;
  // REAC comp (APPUI i) (p.231): per APPUI, the reaction states concomitant
  // with the min / max of its component (every APPUI when none is named)
  std::map<Id, std::pair<Vec6, Vec6>> reac_concom;
  int reac_comp = -1;
  // CONTR bounds (ST1 edition): the fibre stresses of every source, each with
  // the sections it was computed with, enveloped fibre by fibre — the stress
  // of the bounded efforts is not the bound of the stresses.
  struct ContrBound {
    double x = 0;
    std::vector<double> lo, hi;
  };
  std::map<Id, std::vector<ContrBound>> contr;
  // CONTRG bounds likewise (ST1 edition): CONTRG id -> bar -> the bounds of
  // the generalized stress of every source (lo / hi of one value).
  std::map<Id, std::map<Id, std::vector<ContrBound>>> contrg_b;
};

struct SupportDef {
  std::array<bool, 6> fix{};
  std::array<double, 6> spring{};   // EL DI diagonal terms
  std::vector<double> spring_lt;    // EL / EL DI raw terms (full matrix path)
  bool spring_di = false;           // the terms above are the DIAGONAL form
  double rot = 0;                   // ROT alpha about Z
  int decol = 0;                    // 0 none, +1 POS, -1 NEG
  int decol_axis = -1;              // axis of the DECOL criterion (0/1/2)
  Id node = 0;                      // supported node (APPUI liste NOEUD liste2)
};

// One CHARG block, kept as lowered keta loads plus the pieces that need
// deferred processing (CABLE actions, COMB references).
struct ChargDef {
  Id id = 1;
  std::string name;
  std::vector<Load> loads;
  bool is_comb = false;                      // CHARG i COMB
  std::vector<std::pair<Id, double>> parts;  // combined-load references
  // CABLE sub-command requests (lowered at EXEC when prestress is wired):
  struct CableAction {
    std::vector<Id> cables;
    // TENSION: the tension without loss; PERTE: the loss effects only (p.122).
    // A bare `CABLE liste` means TENSION.
    bool tension = false, perte = false, perte_inst = false, relax = false;
    double retrait = 0, uni = 0, relax_cf = 5.0 / 6.0;
  };
  std::vector<CableAction> cable_actions;
};

// ETUDE study sections (manual p.97-98): the blocks are kept in declaration
// order and replayed onto an `EtudePlan` (etude.h) whenever the study grid is
// needed, so ST1's accumulation rules (default / SUP / SEUL) are applied by the
// module that owns them.
struct EtudeDef {
  std::vector<EtudeSpec> blocks;
  bool any() const { return !blocks.empty(); }
};

// COMB combination definition (block form, manual p.230).
struct CombDef {
  Id id = 0;
  std::string name;
  // (kind, id, coef) with kind in {CHARG, COMB, ETAT, PHASAGE, SPECTRE, INCENDIE}
  struct Term {
    std::string kind;
    Id id;
    double coef;
  };
  std::vector<Term> terms;
};

// ---------------------------------------------------------------- results --
struct ResultsStore {
  std::map<Id, CaseResult> charg;
  std::map<Id, CaseResult> comb;
  std::map<Id, CaseResult> etat;
  std::map<Id, CaseResult> phasage;
  std::map<Id, CaseResult> spectre;
  std::map<Id, CaseResult> incendie;
  // essential steps of the phasages (ETAT, DATE): the sources of their ENV
  std::map<Id, CaseResult> phase_step;
  ModalResults modal;
  bool has_modal = false;
  BucklingResults buckling;
  bool has_buckling = false;
  const CaseResult* find(const std::string& kind, Id id) const;
};

// ----------------------------------------------------------------- output --
// SORTIE routing: everything RESU/LISTER/ECRIRE produces goes through here.
struct Output {
  std::string console;                       // accumulated console text
  std::unique_ptr<std::ofstream> file;       // active SORTIE file
  bool to_file = false;
  std::string last_file_path;
  int message_level = 2;                     // 0 AUCUN, 1 REDUIT, 2 TOUT
  void line(const std::string& s);
  // Receives each console line as it is written (live CLI output). Optional.
  std::function<void(const std::string&)> on_console;
};

// ---------------------------------------------------------------- runtime --
class St1Runtime {
 public:
  St1Runtime(Env& env, const St1Options& opt);
  ~St1Runtime();

  // Interpreter sink: execute one engine statement.
  void execute(const Statement& st);

  // Lazily build the keta Model from the current description (active
  // elements only). Throws on incomplete descriptions (missing CARA/CONS).
  Model& snapshot();

  // ST1's RO is a WEIGHT density (poids volumique): POIDS PROPRE weighs
  // RO x A per unit length (times St1Options::gravity, a pure units
  // multiplier, 1 by default), while the built model carries the MASS
  // density RO / G_DYN (G_DYN = 9.81 unless the deck sets it) so the modal
  // and dynamic masses are right. Every CHARG load list goes through this on
  // its way into a Model: the self-weight vectors are scaled by G_DYN when
  // all bar materials share one, else expanded into per-bar uniform loads
  // RO x A (a CARA VAR bar in such a model is a typed error).
  // `target` (default: the current model) is the model the loads go into.
  std::vector<Load> lower_loads(const std::vector<Load>& loads, const Model* target = nullptr);

  // POIDS PROPRE (liste) in direction `dir`: the model-wide self-weight when
  // the list is empty or names every bar, else one span load RO x SX per
  // listed bar (linear on a CARA VAR bar). Shared by CHARG and PHASAGE CHARG.
  std::vector<Load> self_weight_loads(const std::vector<Id>& list, const Vec3& dir);

  // the model with EXACTLY these bars active (a PHASAGE starts empty, p.180)
  Model snapshot_exact_bars(const std::set<Id>& bars);

  // Introspection for import_st1 compatibility.
  const std::string& title() const { return title_; }
  const std::vector<std::string>& warnings() const { return warnings_; }
  ResultsStore& results() { return store_; }
  Output& output() { return out_; }
  bool exec_charg_requested = false;
  bool exec_modal_requested = false;
  bool exec_buckling_requested = false;
  double etude_pas_hint = 0;
  bool finished = false;  // QUITTER

  // --- domain-extension seam (st1_domains.cpp) ---
  // Consulted for statement heads the core dispatch does not own (MAT, PREC,
  // CABLE, TABLIER/SURCH/ALG/CONV/EUROCODE, PHASAGE, the dynamics blocks,
  // the fire chain, DESS graphics) and for EXEC variants the core does not
  // implement. Returns true when the statement was consumed — including
  // while one of the domain's own FIN-blocks is open. Installed by
  // install_st1_domains(); domain_state keeps the extension alive.
  std::function<bool(const Statement&)> domain_handler;
  std::shared_ptr<void> domain_state;
  // Cross-domain result stores (SURCH envelopes feed ENV; ENV results feed
  // RESU/GET MAX/MIN).
  std::map<Id, SurchResults> surch_results;
  std::map<Id, std::string> surch_titles;  // SURCH i 'titre' (the ST1 edition prints it)
  // The CARA in force when a result (kind, id) was computed: ST1 edits the
  // CONTR stresses of a result with the sections it was computed with (a deck
  // may redefine a CARA between two EXEC). Every store point stamps.
  void stamp_cara(const std::string& kind, Id id);
  std::map<std::pair<std::string, Id>, std::shared_ptr<const std::map<Id, CaraProps>>> result_cara_;
  std::shared_ptr<const std::map<Id, CaraProps>> cara_snap_;
  // The APPUI active when a result (kind, id) was computed (stamped with the
  // CARA; a phasing state gets its phase's): ST1 keeps reactions PER DECLARED
  // APPUI — several may sit on one node (p.48), a support of a phasing is
  // another APPUI than the one of the loadings (Fontenoy node 200: APPUI 1
  // for CHARG / SURCH, APPUI 100 for the states).
  std::map<std::pair<std::string, Id>, std::set<Id>> result_appuis_;
  // APPUI redefined by a phasage (MODIFIER APPUI) as they stood in a result
  std::map<std::pair<std::string, Id>, std::map<Id, SupportDef>> result_appui_defs_;
  // The APPUI descriptions of a result: appui_ with its redefinitions.
  std::map<Id, SupportDef> appui_defs_of(const std::string& kind, Id id) const;
  std::set<Id> appuis_of(const std::string& kind, Id id) const;
  // The reaction of each active APPUI: every restrained component of a
  // node's reaction goes to the first active APPUI (declaration order)
  // blocking it.
  // `zero_at_held_node`: an inactive APPUI on a node another active APPUI
  // holds reads 0 instead of being absent (ENV sources, as ST1).
  std::map<Id, Vec6> appui_reactions(const CaseResult& r, const std::set<Id>& active,
                                     bool zero_at_held_node = false,
                                     const std::map<Id, SupportDef>* defs = nullptr) const;
  const std::map<Id, CaraProps>& cara_of(const std::string& kind, Id id) const;
  // CONTR fibre stresses (compression-positive, OPTION order) of station `s`
  // of bar `bid` with the sections `cara` (CARA VAR interpolated).
  std::vector<double> contr_fibres(const std::map<Id, CaraProps>& cara, Id bid,
                                   const Station& s) const;
  // ST1 edition of each SURCH (p.137): per study section / support, the two
  // extremes of every component (or the concomitant states of the declared
  // privileged component). When present it is what RESU prints and what ENV
  // folds (surch_bounds).
  std::map<Id, SurchEnvelope> surch_envelopes;
  // CONTR bounds of each SURCH (ST1 edition): the fibre extremes of its study
  // sections (SurchEnvSpec::effort_fibres), per bar.
  std::map<Id, std::map<Id, std::vector<EnvResult::ContrBound>>> surch_contr;
  // SURCH id -> CONTRG id -> bar -> bounds, each section on the influence
  // line of the generalized stress itself (without its constant, p.99)
  std::map<Id, std::map<Id, std::map<Id, std::vector<EnvResult::ContrBound>>>> surch_contrg;
  std::map<Id, EnvRun> env_results;
  // EXEC HISTORIQUE envelopes {min, max} filled by the dynamics domain so an
  // ENV can fold a HISTORIQUE source (manual p.231: "l'enveloppe d'un
  // historique = enveloppe sur tous les pas de temps").
  std::map<Id, std::pair<CaseResult, CaseResult>> historique_env;
  // ST1 `GET TENSION CABLE i SEC n ...` (manual p.276). GET is a CORE command
  // but the tendon results belong to the prestress domain: the domain
  // extension fills this store, one row per cable section, in the manual's own
  // buffer order — {$nsec, $xl, $sigp, $dsigp, $x, $y, $z}. An empty store
  // makes `GET TENSION CABLE` a typed `$err = 1` warning, never a silent zero.
  std::map<Id, std::vector<std::array<double, 7>>> cable_tension_store;
  // ST1 cable tables (RESU ... CABLE, ST1 edition): per result (kind, id) the
  // cables it carries, one row every 4 x PAS_CABLE of trace length from the
  // origin anchor plus the end — position (global), stress and, for the
  // states, its change from the jacking stress. Filled by the domain.
  struct CableRow {
    double x = 0, y = 0, z = 0, sigma = 0, dsigma = 0;
  };
  struct CableTable {
    Id cable = 0;
    std::vector<CableRow> rows;
  };
  std::map<std::pair<std::string, Id>, std::vector<CableTable>> cable_tables;
  // The time series of each EXEC HISTORIQUE (unless TEMPS NON): per step its
  // time and the state recovered from the displacement field at the study
  // sections (bar efforts, reactions) — RESU HISTORIQUE and the
  // concomitant states an ENV takes from a HISTORIQUE source. Filled by the
  // domain.
  struct HistorySeries {
    std::string title;
    std::vector<double> times;
    std::vector<CaseResult> states;
  };
  std::map<Id, HistorySeries> historique_series;
  // The spectral responses of the last EXEC SPECTRE REPONSE: what a RESU
  // SPECTRE without a list edits (ST1 keeps only them, Exemple6).
  std::vector<Id> spectre_last;

  // --- ETUDE lowering, shared with the domain extension ---
  // Per bar id, the sorted study abscissas [m] of the current ETUDE
  // declaration (manual p.97-98), i.e. exactly what `StaticOptions::stations_abs`
  // expects. EMPTY when the deck declared no ETUDE block: the uniform
  // `stations_per_bar` grid is then kept (it already carries both bar ends,
  // ST1's default study sections), so a deck without ETUDE is bit-identical to
  // the pre-ETUDE behaviour. st1_domains.cpp must call this (or
  // `etude_stations_ref()`) and pass the result to every solve it drives, so
  // the domain results carry the same sections as EXEC CHARG.
  std::map<Id, std::vector<double>> etude_stations(StudyFamily which = StudyFamily::Any);
  // Same study sections replayed on another model (a PHASAGE base, whose bars
  // are not those of the current activation).
  std::map<Id, std::vector<double>> etude_stations_for(const Model& m,
                                                       StudyFamily which = StudyFamily::Any) const;
  // The recovery options every EXEC CHARG solve uses (ETUDE stations, PAS).
  StaticOptions recovery_options();
  // Same map, cached inside the runtime so its address is stable: this is what
  // `StaticOptions::stations_abs` must point at (the option does not own it).
  // Returns nullptr when no ETUDE block was declared.
  const std::map<Id, std::vector<double>>* etude_stations_ref();

  // --- generalized stresses and envelopes (public: RESU / GET / domains) ---
  std::map<Id, ContrgDef> contrgs;    // CONTRG i ... FIN (manual p.99-101)
  std::map<Id, EnvDef> envs;          // ENV definitions, in declaration order
  std::map<Id, EnvResult> env_full;   // their evaluated results
  // `EXEC CHARG PSE` (manual p.121): per CHARG id, the converged soil state of
  // every bar (the displacement/force part is stored in `results().charg`).
  std::map<Id, std::vector<PseBarSoil>> pse_soil;

  // The whole builder/runtime state is deliberately public: the runtime is
  // an internal tool layer and the domain extension is a co-owner of it.
 public:
  // --- statement helpers ---
  // `comma_only`: items continue only through ',' (or an 'a' range) — a
  // number after a space then belongs to the caller (ENV `charg 4 0 1`: case
  // 4, coefficients 0 and 1).
  std::vector<Id> id_list(const Statement& st, std::size_t& i, const std::vector<Id>& all,
                          bool comma_only = false) const;
  std::vector<Id> all_bar_ids() const;
  // Every DEFINED bar, whatever the current activation: a PHASAGE manages
  // its own activations, so its `TOUT` means all of them (p.186).
  std::vector<Id> defined_bar_ids() const { return bar_order_; }
  std::vector<Id> all_node_ids() const;
  std::vector<Id> all_appui_ids() const;
  // the chord length of a declared bar (0 if unknown): CONTRG abscissae x/L
  double bar_length(Id bar) const;

  // --- per-family handlers (Top level) ---
  void cmd_noeud(const Statement&);
  void cmd_gener(const Statement&);
  void cmd_barre(const Statement&);
  void cmd_appui(const Statement&);
  void cmd_exc(const Statement&);
  void cmd_art(const Statement&);
  void cmd_beta(const Statement&);
  void cmd_cara(const Statement&);
  void cmd_cons(const Statement&);
  void cmd_charg_open(const Statement&);
  void cmd_comb(const Statement&);
  void cmd_etude(const Statement&);
  void cmd_activer(const Statement&, bool activate);
  void cmd_exec(const Statement&);
  void cmd_get(const Statement&);
  void cmd_sortie(const Statement&);
  void cmd_lister(const Statement&);
  void cmd_resu_line(const Statement&);
  void cmd_projet(const Statement&, bool save);

  void cmd_masse(const Statement&);
  void cmd_contrg(const Statement&);
  void cmd_env(const Statement&);

  // --- block content handlers ---
  void charg_stmt(const Statement&);
  void comb_stmt(const Statement&);
  void etude_stmt(const Statement&);
  void cara_stmt(const Statement&);
  void cons_stmt(const Statement&);
  void appui_stmt(const Statement&);
  void mat_stmt(const Statement&);
  // Commit one parsed APPUI description onto the builder map (shared by the
  // simple form, the block form and the multi-statement `EL` continuation).
  void apply_appui(const std::vector<Id>& ids, const std::vector<Id>& node_ids,
                   const std::array<bool, 6>& fix, int decol, int decol_axis, double rot,
                   const std::vector<double>& el, bool el_di, bool has_el);
  void cara_var_stmt(const Statement&);
  void cara_pse_stmt(const Statement&);
  void cara_special_stmt(const Statement&);
  void masse_stmt(const Statement&);
  void contrg_stmt(const Statement&);
  void env_stmt(const Statement&);
  void compute_comb(const CombDef&);
  void compute_env(const EnvDef&);
  void render_resu();
  void save_model_json(const std::string& name);
  std::size_t model_node_index(Id nid) const;

  // Field-by-field lo/hi bounds of a stored SURCH (its per-section governing
  // concomitant states enveloped as two pseudo-cases). Cached; typed error
  // when the SURCH stored no state.
  const Envelope& surch_bounds(Id id);

  // --- EXEC implementations ---
  void exec_charg(const Statement&);
  void exec_mode(const Statement&);
  void exec_flambement(const Statement&);
  // CARA SPECIAL secant driver (manual p.210): iterate solve_static with the
  // per-bar secant scaling of tapered.h until every capped component is within
  // its threshold.
  void exec_charg_special(const std::vector<Id>& cases);

  void build_model();  // materialize model_ from builder maps
  void invalidate() {
    built_ = false;
    etude_cache_valid_ = false;
    ++description_gen_;
  }
  // Generation of the description: every change of a builder map (node, bar,
  // CARA, CONS, APPUI, EXC, ART, BETA, activation, OPTION, MASSE, the
  // domain's MAT / PREC / CABLE blocks) advances it; a snapshot does not. A
  // result computed under one generation can be CONTINUED (a PHASAGE
  // resumed by SUITE PHASAGE from its saved state) only under the same one.
  unsigned long description_gen() const { return description_gen_; }
  void touch_description() { ++description_gen_; }
  unsigned long description_gen_ = 0;
  void warn(const std::string& s);

  // --- state ---
  Env& env_;
  St1Options opt_;
  // 0 PLANE, 1 GRILL, 2 SPATIALE. ST1's own default is PLANE: the real decks
  // omit `OPTION` for 2-D frames and rely on the PLANE auto-blocked DOFs.
  int option_ = 0;
  std::string title_;
  std::vector<std::string> warnings_;
  Output out_;

  enum class Mode {
    Top, Nodes, Bars, Charg, ChargBars, ChargNodes, ChargDefor, ChargTemp, ChargAppui,
    Comb, Etude, Pagination, Resu,
    Cara, Cons, Appui, Exc, CaraVar, CaraPse, CaraSpecial, Masse, Contrg, Env, Mat
  };
  Mode mode_ = Mode::Top;
  // One load line of a CHARG body read into `loads` by the grammar of
  // charg_stmt — the phasing's CHARG blocks share it. `sub` carries the open
  // sub-list (BARRE / NOEUD / TEMP / DEFOR / APPUI alone on a line) from one
  // line to the next; lists cover every defined bar (no selective activation).
  void charg_line(std::vector<Load>& loads, Mode& sub, const Statement& st);

  // Builder maps (last description wins; ids are ST1 user ids).
  struct RawBar {
    Id id, n1, n2;
  };
  std::map<Id, Vec3> nodes_;
  std::vector<Id> node_order_;
  std::map<Id, RawBar> bars_;
  std::vector<Id> bar_order_;
  std::map<Id, CaraProps> cara_;
  std::map<Id, ConsProps> cons_;
  std::map<Id, SupportDef> appui_;
  std::vector<Id> appui_order_;
  std::map<Id, std::array<Vec3, 2>> exc_;
  std::map<Id, Releases> art_;
  std::map<Id, double> beta_;
  std::set<Id> active_bars_, active_appuis_;
  bool selective_activation_ = false;

  // CARA VAR / PSE / SPECIAL and MASSE PROPRE descriptions.
  std::map<Id, CaraVarDef> cara_var_;
  SoilZoneMap soil_zones_;                     // CARA PSE (soilzones.h)
  // Bars of `soil_zones_` carried by the exact bedded block (bef.h) in model_;
  // the others (CARA VAR, releases) keep the Hermite zoned bed
  // of the PSE driver.
  std::set<Id> bef_bars_;
  std::map<Id, CaraSpecialDef> cara_special_;  // CARA SPECIAL j (manual p.210)
  MasseDef masse_;
  // Open-block scratch (CARA VAR/PSE bar list, CARA SPECIAL bar + section).
  std::vector<Id> cara_block_bars_;
  TaperLaw cara_block_law_ = TaperLaw::LinY;
  Id cara_special_bar_ = 0;
  int cara_special_section_ = 0;  // 0 = none, 1 = RIGIDITE CARA, 2 = EFFORT
  std::vector<double> cara_special_nums_;
  std::vector<Id> masse_pending_nodes_, masse_pending_bars_;
  Vec3 masse_exc_or_ = Vec3::Zero(), masse_exc_ex_ = Vec3::Zero();
  // Bare `CARA` / `CONS` block instruction: the bar list of the current line,
  // so a continuation line carrying only keywords keeps targeting it
  // (manual p.59 block form, used by the real decks).
  std::vector<Id> cara_pending_bars_;
  // `APPUI … EL m1 … m21` spread over several statements (the real decks lay
  // the lower triangle out row by row): the description is held here until the
  // term count is complete (manual p.49-50).
  bool appui_el_open_ = false;
  std::vector<Id> appui_el_ids_, appui_el_nodes_;
  std::array<bool, 6> appui_el_fix_{};
  int appui_el_decol_ = 0, appui_el_axis_ = -1;
  double appui_el_rot_ = 0;
  std::vector<double> appui_el_vals_;
  bool appui_el_di_ = false;
  std::size_t appui_el_need_ = 0;
  // Minimal fallback material store — see `cmd_mat` in st1_run.cpp.
  std::map<int, ConsProps> mat_props_;
  int open_mat_ = 0;

  std::map<Id, ChargDef> chargs_;
  std::vector<Id> charg_order_;
  ChargDef* open_charg_ = nullptr;
  CombDef open_comb_;
  bool comb_open_ = false;
  std::map<Id, CombDef> combs_;
  EtudeDef etude_;
  EtudeSpec open_etude_;  // the ETUDE block being read
  bool etude_open_ = false;
  std::map<Id, std::vector<double>> etude_cache_;
  bool etude_cache_valid_ = false;
  ContrgDef open_contrg_;
  bool contrg_open_ = false;
  ContrgKind contrg_kind_ = ContrgKind::Cte;
  std::vector<Id> contrg_pending_bars_;
  ContrgCoeffs contrg_or_, contrg_mi_, contrg_ex_;
  EnvDef open_env_;
  bool env_open_ = false;
  std::string env_cur_kind_;  // `CHARG` alone then `id (coef)` lines
  std::map<Id, Envelope> surch_bounds_;
  std::string project_;
  std::string pending_block_;  // family recognized but engine-wired later
  std::vector<std::vector<std::string>> resu_req_;

  // Model snapshot
  Model model_;
  bool built_ = false;
  std::map<Id, Id> bar_index_;  // ST1 bar id -> model bar id (identity today)

  ResultsStore store_;
};

// Session: interpreter + runtime glued together.
struct St1Session {
  Env env;
  std::unique_ptr<St1Runtime> runtime;
  std::unique_ptr<Interp> interp;
  explicit St1Session(const St1Options& opt = {});
  // `file` names the deck in the locations below (erreur.txt, progress).
  void run(const std::string& text, const std::string& file = "");

  // Called around every engine statement, before (done = false) and after
  // (done = true) the runtime executes it, with the statement's location:
  // the CLI's live progress report. Optional.
  std::function<void(const Statement&, const Interp::Where&, bool done)> on_statement;
  // Location of each runtime warning, parallel to runtime->warnings(): the
  // statement that was executing when it was raised.
  std::vector<Interp::Where> warning_where;
};

// erreur.txt as ST1 writes it (manual p.34; layout of the erreur.txt of
// ST1 v2.24): a header, then the
// run's log — a ` - TEXT - ` line per step (EXEC, RESU, SORTIE) and its
// computing time — and every message as
//     --- Fichier : >FILE          (when the file changes)
//     --- Ligne :   207 ---        (when the line changes)
//     attention : TEXT             (or `erreur : TEXT`)
// keta's "ST1: " message prefix is dropped. Each piece goes to `sink` as
// soon as it is known, so the file can be followed during a long run.
class ErreurLog {
 public:
  explicit ErreurLog(std::function<void(const std::string&)> sink) : sink_(std::move(sink)) {}
  void header(const std::string& deck);
  void step(const std::string& text);
  void elapsed(double seconds);
  void message(const Interp::Where& at, bool error, const std::string& text);

 private:
  std::function<void(const std::string&)> sink_;
  std::string file_;
  bool file_known_ = false;
  int line_ = -1;
};

}  // namespace keta::st1

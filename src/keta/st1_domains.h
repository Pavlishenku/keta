// ST1 v24 domain extension: the command families the core runtime dispatch
// (st1_run.cpp) deliberately leaves to a plug-in, wired onto the keta engines.
//
// The core exposes ONE seam — `St1Runtime::domain_handler` — consulted for
// every statement its own block modes do not claim, and again from `cmd_exec`
// for the EXEC variants it does not implement.  `install_st1_domains()` fills
// that seam with the handler built here and parks its state in
// `St1Runtime::domain_state`, so the extension lives exactly as long as the
// runtime.  Nothing in st1_run.* / st1_lang.* / st1.* is touched.
//
// Families wired (manual v24 page numbers are cited on every parser):
//
//   MAT / MATERIAU  (p.74-88)   -> keta/matlaws.h  MatLaw  (+ the fire branch
//                                  `E FEU EC2`, p.214-215 -> FireConcrete)
//   MASSE PROPRE    (p.69-71)   -> Model::nodal_masses / bar_masses
//   PREC (p.89-90) + CABLE (p.91-94) + CHARG CABLE (p.118-120)
//                               -> keta/prestress.h Tendon / tendon_loads
//   TABLIER (p.126-135), SURCH (p.136-167), EUROCODE (p.165), ALG (p.156-159),
//   CONV (p.160-162), EXEC SURCH (p.174)
//                               -> keta/surch.h run_surch
//   PHASAGE (p.175-192) + EXEC PHASAGE
//                               -> keta/phasing.h solve_phased
//   AMORTISSEMENT (p.196-197), SPECTRE (p.198-200), SPECTRE REPONSE
//   (p.201-202), ACCELEROGRAMME (p.203-204), DYNAMIQUE CONVOI (p.205-206),
//   DYNAMIQUE CHARGE (p.207-209)
//                               -> keta/dynamics.h  (+ keta/influence.h for
//                                  the moving convoy)
//   fire chain (p.213-228): MAT E FEU, EXPOSITION FEU, SECTION, ACIER,
//   ARMATURE, CARA ib SECTION i, FEU, INCENDIE, ANALYSE FEU, EXEC INCENDIE
//                               -> keta/heat2d.h + keta/firechain.h solve_fire
//   graphics DESS ... RETOUR (p.243-258)
//                               -> keta/draw.h DrawState
//   LISTER of the families above (p.236-239) -> St1Runtime::output()
//
// CONVENTIONS
//   * SI units in the engine; ST1 decks are user-consistent, so nothing is
//     rescaled here except where the manual names a unit (MAT `MPA/TM2/KNM2`
//     strengths).  The French SURCH catalogue is left in TONNES
//     (`SurchDef::tonne_force = 1`) exactly as ST1 does — the deck's own
//     `POND 9.81` is what converts to kN (manual p.123).
//   * ST1 lists N positive in COMPRESSION (p.60); that flip stays where the
//     core put it — at the GET/RESU boundary — never inside these wirings.
//   * Every ST1 action the underlying engine cannot express is a TYPED
//     WARNING naming the exact gap (`St1Runtime::warn`), never a silent drop.
//   * Keywords are matched case-insensitively and the truncations attested in
//     the reference decks are accepted as aliases: `MATERIAU`=`MAT`,
//     `ACCELEROGRAM`=`ACCELEROGRAMME`, `CABLES`=`CABLE`, `CONVOI`=`CONV`,
//     `ARMA`=`ARMATURE`, `BARRES`=`BARRE`.
//
// STATEMENT-ORDER CONTRACT.  `St1Runtime::execute` consults the domain BEFORE
// its own top-level dispatch, so:
//   * a statement head this extension owns is consumed here and never reaches
//     the core (that is how `TABLIER`/`MAT`/... stop being "pending blocks");
//   * `EXEC CHARG` is INSPECTED (pending CHARG CABLE actions are lowered into
//     the `ChargDef` load list, MAT constants are pushed into the CONS records)
//     and then RELEASED with `false`, so the core's own `exec_charg` solves the
//     enriched cases;
//   * `CARA` is claimed only in its fire form (`CARA liste SECTION i ...`, or
//     the bare `CARA ib` header once a fire SECTION exists), so the R.D.M.
//     `CARA` of p.59-68 stays the core's.
//
// WHAT THE CORE CANNOT READ BACK (documented, never silent)
//   * `RESU HISTORIQUE` / `ENV ... HISTORIQUE j`: the core's `ResultsStore` has
//     no time-history slot, so the per-field envelope of every EXEC HISTORIQUE
//     run is exposed here as `St1Domains::historique_env[j] = {max, min}`
//     (two pseudo-CaseResults).  A warning says so at every run.
//   * `RESU INCENDIE i BARRE b CARA CHAUD`: the hot section state has no slot
//     either; `St1Domains::hot_cara_lines` carries the rendered table (it is
//     also printed to the active SORTIE at EXEC INCENDIE time).
#pragma once
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <vector>

#include "keta/draw.h"
#include "keta/dynamics.h"
#include "keta/firechain.h"
#include "keta/matlaws.h"
#include "keta/phasing.h"
#include "keta/st1_run.h"
#include "keta/surch.h"

namespace keta::st1 {

// ===========================================================================
// MAT — material with time-dependent laws (p.74-88) and its fire branch (p.214)
// ===========================================================================
struct MatDef {
  Id id = 0;
  std::string title;
  MatLaw law;                 // law.code == Direct when `E value` was given
  bool has_retrait = false;   // RETRAIT declared -> eps_sh() is meaningful
  bool has_fluage = false;    // FLUAGE declared  -> phi() is meaningful
  bool creep_traction = true;   // FLUAGE ... TRACTION NON => compression only
  bool creep_linear = true;     // FLUAGE ... LINEAIRE OUI (ST1 default)
  // `RETRAIT <loi>` / `FLUAGE <loi>` naming a law other than E's (p.75-90)
  std::optional<MatLaw::Code> code_retrait, code_fluage;
  // Physical constants that pass straight through to the keta Material.
  double nu = 0.2;
  double ro = 0;       // weight density (ST1 unit) — the core converts with G_DYN
  double g_dyn = 0;
  double alpha = 0;    // TEMP: thermal expansion coefficient
  int amortissement = 0;
  bool has_nu = false, has_ro = false, has_alpha = false;
  // Fire branch (`E FEU EC2`): the material is a FireConcrete, not a MatLaw.
  bool fire = false;
  FireConcrete fire_mat;
  // The deck's stress-unit system (MPA/TM2/KNM2 tag, p.75): fck is stored in
  // SI Pa for the laws; every modulus/strength EMITTED back into the model
  // must be divided by `stress_scale()` so it lands in the USER units the
  // rest of the deck (loads, sections) is written in.
  StrengthUnit unit = StrengthUnit::MPa;
  double stress_scale() const { return strength_to_pa(1.0, unit); }
  // T_THERM — heat cure of a pretensioned element (p.80, annex B4.6): the
  // mean strength at its end, the share lambda of the thermal loss the strand
  // does not take, the cure steps (days) and their temperatures (deg C).
  bool therm = false;
  double therm_fcmp = 0, therm_lambda = 0.5;
  std::vector<double> therm_dt, therm_theta;
  double therm_end() const;        // age at the end of the cure [days]
  double therm_dtheta_max() const; // theta_max - 20 [deg C]
  double therm_teb() const;        // relaxation equivalent time t_eb [days] (EC2 eq. 10.2)
};

// ===========================================================================
// PREC — prestress family (p.89-90)
// ===========================================================================
struct PrecDef {
  Id id = 0;
  std::string title;
  double ap = 0;        // SECTION  [m^2]
  double sigma = 0;     // TENSION  [Pa]
  double ep = 195e9;    // E        [Pa]
  double fpk = 0;       // RG / FPK [Pa]
  double rho1000 = 0.0; // R1000    [%] (absent: no relaxation, manual p.93)
  double recul = 0;     // RECUL    [m]
  std::optional<double> seuil;  // SEUIL [h]
  double f = 0;         // PERTE F  / PERTE MU : friction coefficient [1/rad]
  double phi = 0;       // PERTE PHI: wobble    [rad/m]  (BPEL form)
  double kphi = 0;      // PERTE KPHI (Eurocode form; kphi = phi / f)
  bool eurocode_form = false;
  double perte_inst = 0;  // PERTE INST: uniform instantaneous loss [Pa]
  // LOI <BPEL, EC2, EC2_2027> (KP kp) (p.89-90): the relaxation law of the
  // steel; EC2_2027 is the power law of annex B4.10.2 (eq. B.43) in a
  // phasing, the BPEL rule with mu0 = 0.35 outside (annex B2.2)
  bool loi_ec2_2027 = false;
  int loi = 0;  // 0 = the law of the concrete (p.93), 1 = BPEL, 2 = EC2, 3 = EC2_2027
  double kp = 0.16;
  bool pretension = false;
  // PRETENSION transmission law (p.92-93): the force grows from 0 at the
  // bonded end over the transmission length. 0 = none given (full force at
  // the bonded end), 1 LINEAIRE, 2 PARABOL, 3 BILINEAIRE (lg1, lg2, alpha).
  int pret_law = 0;
  double lg_scel1 = 0, lg_scel2 = 0, pret_alpha = 1;
  // Fraction of the force at distance d from the bonded end.
  double transmission(double d) const;
  double jack_force() const { return sigma * ap; }
};

// ===========================================================================
// CABLE — tendon trace (p.91-94)
// ===========================================================================
struct TracePoint {
  Vec3 p = Vec3::Zero();
  bool aligne = false;             // straight segment from the previous point
  std::optional<double> pente;     // PENTE  [rad] -> de_z/ds
  std::optional<double> gis;       // GIS    [rad] -> de_y/ds
  std::optional<double> rayon;     // RAYON  [m]
};

struct CableDef {
  Id id = 0;
  std::string title;
  Id prec = 0;
  int simul = 1;
  bool non_simul = false;
  int anchor = 0;         // 0 = TENSION OR, 1 = TENSION EX, 2 = OR EX, 3 = EX OR
  std::vector<Id> bars;   // BARRE liste (empty = every active bar)
  double pas_cable = 0.25;
  double lg_gaine_or = 0, lg_gaine_ex = 0;
  bool exterieur = false;
  std::vector<TracePoint> trace;
};

// Curvilinear tension table of a cable, the source of `GET TENSION CABLE`
// (`$nsec $xl $sigp $dsigp $x $y $z`, manual p.276).
struct CableTension {
  std::vector<double> s;      // chain abscissa [m]
  std::vector<double> P;      // force after losses [N]
  std::vector<double> sigma;  // P / Ap [Pa]
  std::vector<Vec3> pos;      // global position of the tendon point
  double sigma0 = 0;          // jacking stress (PREC TENSION)
};

// ===========================================================================
// Dynamics blocks
// ===========================================================================
struct DampingDef {
  Id id = 0;
  DampingSpec spec;
  double pond = 1.0;
  bool energie = false;  // ENERGIE KSI: weighted by each bar's modal strain energy
};

// One EN 1998-1 branch of a SPECTRE block. ST1 lets ONE `SPECTRE isp` carry a
// HORIZONTAL and a VERTICAL line (p.199), so the parameters are per axis.
struct Ec8Branch {
  bool set = false;
  double ag = 0, s_soil = 1, st_topo = 1, alpha = 2.5, eta = 1;
  double q = 1, beta0 = 0.2, gamma = 2.0 / 3.0;
  double tb = 0.15, tc = 0.4, td = 2.0, te = 0;
};

struct SpectrumDef {
  Id id = 0;
  double pond = 1.0;
  double beta = 0;  // horizontal-plane angle of the local Ox [rad] (SPATIALE)
  bool eurocode = false, reduit = false;
  // Tabulated form: one (period, Sa) table per global direction.
  std::vector<std::pair<double, double>> ax, ay, az;
  std::array<Ec8Branch, 3> ec8;  // X, Y, Z
  ResponseSpectrum build(int dir) const;
};

struct SpectreResponseDef {
  Id id = 0;
  Id spectrum = 0;
  ComboRule rule = ComboRule::Srss;
  bool linear = false;  // LIN: algebraic sum of the signed modal peaks
  std::array<double, 3> pond{{0, 0, 0}};
};

struct AccelDef {
  Id id = 0;
  std::string title;
  int methode = 1;       // 1 Newmark moyenne, 2 Newmark lineaire, 3 HHT
  double hht_alpha = 0;  // METHODE 3 ALPHA
  double dt = 0;         // PAS_TEMPS [s]
  double pond = 1.0;
  double beta = 0;
  std::vector<std::array<double, 3>> accel;  // AX AY (AZ) per step
};

struct DynConvoiDef {
  Id id = 0;
  int methode = 1;
  double hht_alpha = 0;
  Id convoi = 0;
  double dt = 0, speed = 0, duration = 0, pond = 1.0;
};

struct DynChargeDef {
  Id id = 0;
  int methode = 1;
  double hht_alpha = 0;
  double dt = 0, duration = 0;
  int mode = 1;
  double f0 = 0;  // 0 = use the mode's own frequency
  std::vector<double> phi, coef;
  struct BarLine {
    std::vector<Id> bars;
    Vec3 F = Vec3::Zero(), M = Vec3::Zero();
    LoadFrame frame = LoadFrame::Global;
  };
  std::vector<BarLine> lines;
};

// Envelope of one time-history run (EXEC HISTORIQUE / DYNAMIQUE ...).
struct HistoryResult {
  Id id = 0;
  std::string name;
  CaseResult max, min;   // per-field extremes over the whole history
  int steps = 0;
  double dt = 0;
  bool keep_time = true;
  std::vector<double> peak_disp;  // |u|_inf per step (kept even for TEMPS NON)
};

// ===========================================================================
// PHASAGE (p.175-192)
// ===========================================================================
struct PhasageDef {
  Id id = 0;
  std::string title;
  PhasePlan plan;
  std::vector<LoadCase> cases;                    // CHARG blocks nested in the phasage
  std::vector<std::pair<std::size_t, Id>> etats;  // (phase index, ETAT id)
  std::vector<std::string> etat_names;
  std::vector<Id> suite;  // SUITE PHASAGE liste (p.182): phasages continued by this one
  // per phase (parallel to plan.phases): the APPUI it activates / deactivates
  std::vector<std::pair<std::vector<Id>, std::vector<Id>>> appui_ops;
  // per phase (parallel to plan.phases): the APPUI it redefines (MODIFIER
  // APPUI, p.195) with their new description
  std::vector<std::map<Id, SupportDef>> appui_mods;
  // ENV i of the phasage (p.183): their definitions (sources filled at EXEC)
  // and, per phase, those active when the phase ends on an essential step
  // (ETAT, DATE, end of the phasage) — empty for an intermediate one.
  std::vector<EnvDef> envs;
  std::vector<std::set<Id>> phase_envs;
};

// ===========================================================================
// Fire chain (p.213-228)
// ===========================================================================
struct FireSectionDef {
  Id id = 0;
  std::string title;
  FireSection section;
  Id mat = 0;                     // MAT m (a MatDef with `fire` set)
  std::array<Id, 4> expo_ids{{0, 0, 0, 0}};  // ST1 exposure ids: VY WY VZ WZ
};

struct ArmatureDef {
  Id id = 0;
  std::vector<RebarLayer> layers;
  std::vector<Id> steels;  // ACIER id per layer
};

struct FireBarDef {
  Id bar = 0;
  Id section = 0;
  Id armature = 0;
};

struct FeuDef {
  Id id = 0;
  std::string title;
  double duree = 60;
  double temp = 20;
  double maille = 0;  // 0 = keep each section's own MAILLE
};

struct IncendieDef {
  Id id = 0;
  std::string title;
  Id feu = 0;
  bool local = false;
  Id charg = 0;
};

// MASSE PROPRE is owned by the CORE runtime (st1_run parses the block into
// the Model's NodalMass/BarMass tables directly); the domain keeps no state
// for it. `DomMasseDef` remains only as the retired container type.
struct DomMasseDef {
  std::vector<NodalMass> nodes;
  std::vector<BarMass> bars;
  bool declared = false;
};

// ===========================================================================
// The extension state
// ===========================================================================
class St1Domains {
 public:
  explicit St1Domains(St1Runtime& rt) : rt_(&rt) {}

  // The seam. Returns true when the statement was consumed.
  bool handle(const Statement& st);

  // ---- descriptions (public: the runtime layer is an internal tool layer) ----
  std::map<Id, MatDef> mats;
  // ST1 bar -> MAT id, captured before sync_materials() materializes the
  // constants into ConsProps (the phasing creep path needs the law identity).
  std::map<Id, int> bar_mat;
  std::map<Id, PrecDef> precs;
  std::map<Id, CableDef> cables;
  std::map<Id, CableTension> cable_tension;
  Tablier tablier;
  bool has_tablier = false;
  std::map<Id, SurchDef> surchs;
  // SURCH privileged fibre (p.137): `CONTR <VY,WY,...>` -> {"CONTR", name},
  // `CONTRG j` -> {"CONTRG", "j"}; PRESS stays unsupported (warned)
  std::map<Id, std::pair<std::string, std::string>> surch_priv_fibre;
  std::map<Id, AlgDef> algs;
  std::map<Id, ConvoyDef> convs;
  // ST1 states the predefined Eurocode loads in kN (p.127): the deck's units.
  EurocodeParams eurocode = EurocodeParams{}.scaled_forces(1e-3);
  // EUROCODE `CLASSE n (AJUSTEMENT ...)` (p.170): the adjustment factors of
  // each traffic class; TABLIER `CLASSE TRAFIC icltra` (default 2) picks one
  // at EXEC SURCH. Classes not redefined keep the manual's table (p.164).
  std::map<int, EurocodeParams> ec_class_;
  int ec_cur_class_ = 0;
  EurocodeParams eurocode_for_class(int cls) const;
  std::map<Id, PhasageDef> phasages;
  // One executed PHASAGE kept alive: the model its chain ran on, the options
  // (and what they point at), the chain itself and the engine session at its
  // end state. A later `SUITE PHASAGE i` continues from it — the common
  // prefix is never replayed from the first date (keta/phasing.h
  // PhasingSession::fork) — when the description is at the same generation,
  // its chain is the stored one and the engine accepts the new model.
  struct PhasageRun {
    Model model;
    std::set<Id> st1_creep_bars, st1_tension_creep_bars;
    std::map<Id, std::vector<double>> sections;
    SoilZoneMap pse_soil;
    PhasingOptions po;
    bool pse = false;
    bool resumed = false;  // continued from a kept run (else replayed from its first date)
    unsigned long description_gen = 0;
    // the chain run: per phase, the phase, its APPUI (de)activations, its
    // MODIFIER APPUI; and the lowered CHARG blocks (ids renumbered)
    std::vector<std::tuple<Phase, std::pair<std::vector<Id>, std::vector<Id>>,
                           std::map<Id, SupportDef>>> seq;
    std::vector<LoadCase> cases;
    std::unique_ptr<PhasingSession> session;
  };
  std::map<Id, std::shared_ptr<PhasageRun>> phasage_runs;
  std::map<Id, DampingDef> dampings;
  std::map<Id, SpectrumDef> spectra;
  std::map<Id, SpectreResponseDef> spectre_reps;
  std::map<Id, AccelDef> accels;
  std::map<Id, DynConvoiDef> dyn_convois;
  std::map<Id, DynChargeDef> dyn_charges;
  DomMasseDef masse;  // retired: never filled (core owns MASSE PROPRE)

  // Results the core's ResultsStore has no slot for (see the header note).
  std::map<Id, HistoryResult> historique;
  std::map<Id, std::pair<CaseResult, CaseResult>> historique_env;  // {max, min}
  std::map<Id, FireOutcome> fire_outcomes;
  std::vector<std::string> hot_cara_lines;  // RESU INCENDIE ... CARA CHAUD TEMP

  // fire chain descriptions
  std::map<Id, Exposure> exposures;
  std::map<Id, FireSectionDef> fire_sections;
  std::map<Id, RebarSteel> aciers;
  std::map<Id, ArmatureDef> armatures;
  std::map<Id, FireBarDef> fire_bars;
  std::map<Id, FeuDef> feux;
  std::map<Id, IncendieDef> incendies;
  FireAnalysisOptions fire_opt;

  // graphics
  DrawState draw;
  std::string last_drawing_path;

  // Tendons built at the last EXEC (id = CABLE id).
  std::map<Id, Tendon> tendons;

  // The keta tendon of a CABLE on model `m` (its trace projected on the bar
  // chain, its PREC losses when `with_losses`).
  Tendon build_tendon(const Model& m, const CableDef& c, const PrecDef& p, bool with_losses);
  // The rows of the ST1 cable table of tendon `t` (St1Runtime::cable_tables):
  // every 4 x `pas` of 3-D trace length between its anchors plus the end,
  // `sigma_at(s)` the stress at chain abscissa s, `sigma0` the jacking stress.
  std::vector<St1Runtime::CableRow> cable_rows(const Model& m, const Tendon& t, double pas,
                                               const std::function<double(double)>& sigma_at,
                                               double sigma0,
                                               const std::function<double(double)>& ref = nullptr) const;

  // Push the MAT-derived constants into the runtime's CONS records (and
  // synthesise CARA/CONS for the fire bars, which carry neither, p.224).
  void sync_materials();

 private:
  enum class Block {
    None, Mat, Prec, Cable, Tablier, Surch, Alg, Conv, Eurocode, Phasage,
    Amortissement, Spectre, SpectreRep, Accel, DynConvoi, DynCharge, Masse,
    Exposition, FireSection, Acier, Armature, CaraFeu, Feu, Incendie, AnalyseFeu
  };

  // --- block plumbing ---
  bool open_block(const Statement& st, const std::string& head);
  void block_stmt(const Statement& st, const std::string& head);
  void close_block();

  // --- per-family parsers ---
  void mat_stmt(const Statement&, const std::string& head);
  void prec_stmt(const Statement&, const std::string& head);
  void cable_stmt(const Statement&, const std::string& head);
  void tablier_stmt(const Statement&, const std::string& head);
  void surch_stmt(const Statement&, const std::string& head);
  void alg_stmt(const Statement&, const std::string& head);
  void conv_stmt(const Statement&, const std::string& head);
  void eurocode_stmt(const Statement&, const std::string& head);
  void phasage_stmt(const Statement&, const std::string& head);
  void amortissement_stmt(const Statement&, const std::string& head);
  void spectre_stmt(const Statement&, const std::string& head);
  void spectre_rep_stmt(const Statement&, const std::string& head);
  void accel_stmt(const Statement&, const std::string& head);
  void dyn_convoi_stmt(const Statement&, const std::string& head);
  void dyn_charge_stmt(const Statement&, const std::string& head);
  void masse_stmt(const Statement&, const std::string& head);
  void exposition_stmt(const Statement&, const std::string& head);
  void fire_section_stmt(const Statement&, const std::string& head);
  void acier_stmt(const Statement&, const std::string& head);
  void armature_stmt(const Statement&, const std::string& head);
  void cara_feu_stmt(const Statement&, const std::string& head);
  void feu_stmt(const Statement&, const std::string& head);
  void incendie_stmt(const Statement&, const std::string& head);
  void analyse_feu_stmt(const Statement&, const std::string& head);

  // CHARG nested in a PHASAGE (p.186) — the chapter-7 load syntax, re-parsed.
  void phas_charg_stmt(const Statement&, const std::string& head);

  // --- EXEC ---
  bool exec_stmt(const Statement& st);
  void exec_surch(const Statement&);
  void exec_phasage(const Statement&);
  void exec_spectre_reponse(const Statement&);
  void exec_historique(const Statement&);
  void exec_incendie(const Statement&);
  void lower_cable_actions();

  // --- graphics ---
  void dess_stmt(const Statement&, const std::string& head);
  std::vector<const CaseResult*> draw_sources(const Statement& st, std::size_t i) const;

  // --- listings ---
  bool lister(const Statement& st, const std::string& what);
  void get_tension(const Statement& st);

  // --- helpers ---
  Model& model();
  double deferred_fraction(const PrecDef& p, const ChargDef::CableAction& a) const;
  // The deferred loss (RETRAIT, UNI, RELAX) at a point whose stress after the
  // instantaneous losses is `sigma` [Pa].
  double deferred_loss(const PrecDef& p, const ChargDef::CableAction& a, double sigma) const;
  DampingSpec damping_for() const;
  Support support_of(Id appui_id) const;
  void flush_phase(bool essential = false);
  void warn(const std::string& s) const { rt_->warn(s); }

  St1Runtime* rt_;
  Block block_ = Block::None;
  Id block_id_ = 0;
  // Sub-mode inside the open block, block-specific:
  //   Tablier 1 = REP_TRANS body, 2 = ZONE_TRANS body
  //   Cable   1 = TRACE point lines
  //   Alg     1 = DENS pairs        Conv 1 = ESSIEU lines
  //   Spectre 1 = ACCELERATION table, 2 = EUROCODE axis lines
  //   Accel   1 = ACCELERATION rows  FireSection 1 = PERIMETRE vertices
  int sub_ = 0;
  bool in_dess_ = false;
  bool phas_charg_ = false;      // a CHARG block is open inside the PHASAGE
  St1Runtime::Mode phas_charg_mode_ = St1Runtime::Mode::Top;  // open CHARG sub-list
  LoadCase phas_case_;
  Phase phas_cur_;               // phase under construction
  std::vector<Id> phas_appui_on_, phas_appui_off_;  // its APPUI (de)activations
  std::map<Id, SupportDef> phas_appui_mods_;  // its MODIFIER APPUI
  std::map<Id, SupportDef> phas_appui_defs_;  // APPUI redefined so far in the PHASAGE
  std::set<Id> phas_env_on_;                  // its ENV currently taking the extremes
  // ST1 runs a phasage's commands in sequence; keta's engine orders a phase
  // as activations, deactivations, support changes, loads, tensioning. A
  // command of an EARLIER kind than one already in the phase opens a new
  // phase at the same date (no creep in between): an APPUI activated after
  // a load does not carry it.
  int phas_kind_ = 0;
  // a DATE was just given: its creep step is computed on its own before the
  // next command (an ENV then sees, e.g., the hogging creep under prestress
  // from 10 to 33 days before the superstructures are applied)
  bool phas_date_step_ = false;
  int phas_tension_group_ = 0;  // TENDRE / DETENDRE lines of the phasage
  bool phas_et_ = false;        // the command being read is chained by ET (p.180)
  void phase_step(int kind);
  double phas_time_ = 0;
  bool phas_dirty_ = false;      // the current phase carries at least one action
  Id phas_id_ = 0;
  std::vector<Id> cara_ids_;     // bars of the open fire `CARA liste` block
  std::set<Id> cable_lowered_;   // CHARG ids whose cable actions were lowered
};

// Install the extension on a runtime (idempotent: replaces any previous one).
void install_st1_domains(St1Runtime& rt);
// The installed extension, or nullptr.
St1Domains* st1_domains(St1Runtime& rt);

}  // namespace keta::st1

// ST1 "domaine d'étude" + result-algebra glue, three families in one module:
//
//  * ETUDE  (manual v24 p.97-98) — study sections: per bar, the sorted set of
//    abscissas where results are evaluated, each flagged EFFORT (internal
//    forces / stresses) and/or DEPLA (displacements / soil pressures).
//    `EtudePlan` implements ST1's exact accumulation rules (defaults, plain
//    ETUDE, ETUDE SUP, ETUDE SEUL) and lowers to the per-bar abscissa map that
//    `StaticOptions::stations_abs` feeds to the statics recovery.
//
//  * CONTRG (p.99-101) — generalized stress: a linear combination of the six
//    internal forces plus a constant, with per-bar coefficients that are
//    constant (CTE), linearly interpolated (VAR LIN, OR->EX) or parabolically
//    interpolated (VAR PARA, through OR/MI/EX) along the bar.
//
//  * ENV / ENV COMB (p.230-235) — the general envelope evaluator over stored
//    results: a privileged quantity (EFFORT / DEPLA / CONTR fiber stress /
//    CONTRG / PRESS, plus an optional REAC component), sources carrying a
//    [cfmin, cfmax] coefficient box, extremes at box corners with the ENTIRE
//    companion state assembled from the SAME corner (a consistent state, never
//    a mix of independent extrema).
//
// SI units throughout (m, N, Pa). Every unsupported composition is a typed
// `Error`, never a silent approximation.
#pragma once
#include <map>

#include "keta/envelope.h"

namespace keta {

// ===================== 6.1 — study sections (ETUDE) =====================

// Which result family a study section carries (ST1 EFFORT / DEPLA keywords,
// p.98: EFFORT = R.D.M. forces and stresses, DEPLA = displacements and, for
// beams on elastic soil, the soil pressures).
enum class StudyFamily { Any, Effort, Depla };

struct StudySection {
  double x = 0;         // abscissa along the bar [m], measured from end i
  bool effort = false;  // ST1 EFFORT
  bool depla = false;   // ST1 DEPLA
};

// Accumulation mode of one ETUDE block (p.97, the three documented headers):
//   Default — `ETUDE …`      : addition to the sections defined BY DEFAULT
//                              (previously defined extra sections are dropped);
//   Sup     — `ETUDE SUP …`  : addition to the sections ALREADY defined;
//   Seul    — `ETUDE SEUL …` : definition WITHOUT the defaults and without any
//                              previously defined section.
enum class EtudeMode { Default, Sup, Seul };

// One `liste1 SE liste2 (PAS lpas) (ABS,REL)` line of an ETUDE block.
struct EtudeLine {
  std::vector<Id> bars;           // liste1
  bool all = false;               // liste1 replaced by TOUT (p.98)
  std::vector<double> abscissas;  // liste2; with `step`, the range [lo, hi]
  double step = 0;                // PAS lpas (0 = no step)
  bool absolute = false;          // ABS; the default is REL (p.98)
};

// One ETUDE block: mode + result-family flags in the header, then its lines.
// ST1 syntax puts (EFFORT)(DEPLA) on the header and (ABS,REL) on each line.
// Neither flag given = both families (keta reading; the manual leaves the
// empty-flag case unspecified and a section carrying no family is pointless).
struct EtudeSpec {
  EtudeMode mode = EtudeMode::Default;
  bool effort = false, depla = false;
  std::vector<EtudeLine> lines;
};

// Legality of a RELATIVE step (p.97): "la valeur de lpas doit être une valeur
// fractionnelle entière. Les valeurs 1/2, 1/5, 1/10 sont correctes tandis que
// 1/3, 1/7 sont incorrectes."  Extracted rule: 1/lpas must be an integer n
// (the step divides the bar evenly) AND n must be a terminating-decimal
// denominator, i.e. n = 2^a·5^b — which is exactly what separates the legal
// 1/2, 1/4, 1/5, 1/10, 1/20 from the illegal 1/3, 1/6, 1/7, 1/9. An absolute
// step (ABS) carries no such restriction.
bool etude_step_is_legal(double lpas);

// The set of study sections of a model, built by ST1's accumulation rules.
//
// Abscissas are stored in METERS along the element axis, measured from end i,
// on the SAME length the statics recovery uses (node distance including the
// rigid-arm eccentricities, see `element.cpp`), so a relative abscissa lands
// exactly where the recovery expects it.
class EtudePlan {
 public:
  // ST1 default (p.97): "toutes les barres contiennent une section d'étude des
  // efforts à chaque extrémité" — EFFORT at x = 0 and x = L of every bar.
  void add_default(const Model& model);

  // Apply one ETUDE block with its mode semantics.
  void add(const Model& model, const EtudeSpec& spec);

  void clear() { bars_.clear(); }
  bool empty() const { return bars_.empty(); }

  // Sections of one bar (nullptr when the bar carries none).
  const std::vector<StudySection>* sections(Id bar) const;
  const std::map<Id, std::vector<StudySection>>& all() const { return bars_; }

  // Lowering to the statics hook: per-bar sorted abscissas [m].
  //
  // keta's `Station` carries the forces AND the deflection, so one recovery
  // pass serves both families: the EFFORT/DEPLA flags are kept as reporting
  // metadata and `StudyFamily::Any` (the default) returns their union. Asking
  // for a single family narrows the grid — a documented deviation from ST1,
  // which evaluates the two families independently.
  std::map<Id, std::vector<double>> stations_abs(StudyFamily which = StudyFamily::Any) const;

 private:
  void insert(Id bar, double L, double x, bool effort, bool depla);
  std::map<Id, std::vector<StudySection>> bars_;
};

// Element-axis length of a bar (chord including the rigid-arm eccentricities) —
// the length the station abscissas are measured along.
double bar_axis_length(const Model& model, const Bar& bar);

// ================ 6.2 — generalized stress (CONTRG) ================

// One coefficient record: the weights of the six R.D.M. forces plus the
// constant (p.99). `cst` is documented as "applicable uniquement pour les
// charges, non-applicable aux surcharges" — see `contrg_at(..., include_cst)`.
struct ContrgCoeffs {
  double cfn = 0;   // N  — axial
  double cfty = 0;  // ty — shear along y
  double cftz = 0;  // tz — shear along z
  double cfmx = 0;  // mx — torsion
  double cfmy = 0;  // my — bending about y
  double cfmz = 0;  // mz — bending about z
  double cst = 0;   // CST
};

// Variation of the coefficients along the bar (p.99): constant, linear between
// the origin (OR) and the far end (EX), or parabolic through OR / MI / EX.
enum class ContrgKind { Cte, VarLin, VarPara };

struct ContrgRule {
  ContrgKind kind = ContrgKind::Cte;
  ContrgCoeffs or_;  // xx = 0   (CTE uses this record alone)
  ContrgCoeffs mi;   // xx = 1/2 (VAR PARA only)
  ContrgCoeffs ex;   // xx = 1   (VAR LIN / VAR PARA)
};

// Sign convention the CONTRG coefficients were authored in.
//
// keta stations are TENSION-POSITIVE (N > 0 in tension, +My tensions the +z
// fiber; see model.h). ST1 lists the same physical state with the opposite
// normal-stress convention: p.66 "l'effort normal N est positif en compression
// […] la contrainte normale est positive en compression", and p.68 gives the
// stress of a section as sigma = N/Sx + My/Iy·z − Mz/Iz·y with that
// compression-positive sigma. keta's tension-positive sigma obeys the same
// formula shape, so term by term (the fiber functions 1, z, y are independent):
//
//     n_ST1 = −N_keta,  my_ST1 = −My_keta,  mz_ST1 = −Mz_keta
//
// The manual does NOT fix the sense of ty / tz / mx anywhere in the ETUDE /
// CONTRG / ENV chapters, so `St1` passes those three through unchanged (the
// classic French mixed convention). A deck whose shear/torsion listings also
// run the other way must simply negate cfty / cftz / cfmx — never guessed here.
enum class ForceConvention { Keta, St1 };

// A generalized stress: one id, one coefficient rule per bar. The rule may
// cover only part of the model (p.101: "Les contraintes généralisées peuvent
// ne pas être définies sur toutes les barres").
struct ContrgDef {
  Id id = 0;
  std::string title;
  ForceConvention convention = ForceConvention::Keta;
  std::map<Id, ContrgRule> bars;

  void set_cte(const std::vector<Id>& bar_ids, const ContrgCoeffs& c);
  void set_var_lin(const std::vector<Id>& bar_ids, const ContrgCoeffs& or_,
                   const ContrgCoeffs& ex);
  void set_var_para(const std::vector<Id>& bar_ids, const ContrgCoeffs& or_,
                    const ContrgCoeffs& mi, const ContrgCoeffs& ex);
  bool defined_for(Id bar) const { return bars.find(bar) != bars.end(); }
};

// Station forces re-signed into ST1's listing convention (see ForceConvention).
Station st1_signs(const Station& st);

// Coefficients of `def` on `bar` at the fractional abscissa xx = x/L in [0, 1].
//   VAR LIN  (p.100): c = c_or + (c_ex − c_or)·xx
//   VAR PARA (p.101): c = a·xx² + b·xx + c0 with
//                     a = 2(c_or − 2c_mi + c_ex), b = −3c_or + 4c_mi − c_ex,
//                     c0 = c_or.
//   NOTE — the manual prints "c = cfex" for the constant term; that is a typo:
//   the parabola through (0, c_or), (1/2, c_mi), (1, c_ex) requires c0 = c_or
//   (with c0 = c_ex the curve misses OR unless c_or = c_ex). keta implements
//   the interpolating parabola, i.e. the manual's own stated intent.
ContrgCoeffs contrg_coeffs(const ContrgDef& def, Id bar, double xx);

// Generalized stress at one station. `xx` is the fractional abscissa (x/L);
// the station forces are re-signed per `def.convention` before weighting.
// `include_cst` mirrors p.99: the constant applies to load cases ("charges")
// and NOT to influence-line surcharges ("surcharges") — pass false there.
double contrg_at(const ContrgDef& def, Id bar, double xx, const Station& st,
                 bool include_cst = true);

// ================ 13.2 / 13.4 — ENV and ENV COMB ================

// The eight CONTR fiber components (p.231, option SPATIALE). ST1's extreme
// fibers come from CARA (p.60): vy / vz are the ordinates of the UPPER fibers
// on the local y / z axes, wy / wz the ABSOLUTE ordinates of the lower ones —
// so the sampled point is (y, z) = (±vy or −wy, +vz or −wz). The four crossed
// components are the section corners.
enum class ContrFiber { VY, WY, VZ, WZ, VYVZ, VYWZ, WYVZ, WYWZ };

struct SectionFibers {
  double vy = 0, wy = 0, vz = 0, wz = 0;  // all >= 0, ST1 CARA meaning
};

// Extreme-fiber ordinates of a section: taken from the named stress points when
// the section has them, else derived from the elastic moduli (symmetric
// fallback z_max = Iy/Wely, y_max = Iz/Welz). Typed error when neither exists —
// ST1 p.66 likewise makes VY/WY/VZ/WZ mandatory for normal-stress output.
SectionFibers section_fibers(const Section& sec);

// Normal stress at a fiber: sigma = N/A + My·z/Iy − Mz·y/Iz (p.68), returned
// TENSION-POSITIVE like every keta stress. `st1_sign` returns it
// compression-positive, the way ST1 lists it.
double fiber_sigma(const Model& model, const Bar& bar, const Station& st, ContrFiber fiber,
                   bool st1_sign = false);

// The privileged quantity of an envelope (p.231 / p.234).
//
// ST1 §13.3.2 "cas où aucune composante n'est déclarée" — an envelope with NO
// privileged component is the plain field-by-field envelope, which lives in
// envelope.h (`envelope` / `envelope_factored`); `env_run` always envelopes on
// one declared component. ST1's coherence warning (mixing sources enveloped on
// different components, p.233) has no equivalent here for the same reason: one
// run, one component.
enum class EnvWhat { Effort, Depla, Contr, Contrg, Press };

struct EnvSpec {
  EnvWhat what = EnvWhat::Effort;

  // Location. Station quantities (Effort / Contr / Contrg / Press) use
  // (bar, station); Depla uses (node). `env_run_sections` ignores bar/station
  // and sweeps every section of every bar.
  Id bar = 0;
  std::size_t station = 0;
  Id node = 0;

  // Effort: 0..5 = N TY TZ MX MY MZ (keta station order N Vy Vz T My Mz).
  // Depla : 0..5 = DX DY DZ RX RY RZ (global node displacement).
  // Press : 0..3 = PX PY PZ MX (soil pressure of a Winkler bed).
  int comp = 0;

  ContrFiber fiber = ContrFiber::VY;     // Contr
  bool sigma_st1_sign = false;           // report CONTR compression-positive
  const ContrgDef* contrg = nullptr;     // Contrg (must outlive the call)
  bool contrg_include_cst = true;        // p.99: false for "surcharges"

  // Optional REAC component (p.231): "Une concomitance peut être demandée à la
  // fois sur les réactions d'appui et sur les sections d'étude" — the reaction
  // is enveloped on its own component, with its own corner choice.
  bool has_reac = false;
  Id reac_node = 0;
  int reac_comp = 0;  // 0..5 = FX FY FZ MX MY MZ
};

// One source of an envelope. Plain ENV gives one coefficient (cmax = cmin);
// ENV COMB gives the [cfmin, cfmax] box (p.234: "cfmax=cfmin" by default).
//
// ENV of ENV: an envelope is a SET of possible states, not an additive action.
// When `alts` is non-empty it replaces `r`, and the term contributes the single
// most favourable alternative (never their sum) — so wrapping a previous
// envelope's lo/hi bounds as alternatives reproduces ST1's "enveloppe
// d'enveloppes" (p.232-233), whose concomitance type is inherited from the
// inner envelope: keta expresses that by re-running the outer envelope on the
// SAME privileged component.
struct EnvSource {
  const CaseResult* r = nullptr;
  double cmin = 1.0, cmax = 1.0;
  std::vector<const CaseResult*> alts;
  // The alternatives of a plain ENV: at a SECTION, one whose RAW value (its
  // result before the source coefficient, `raw[a]`) is null there does not
  // compete unless they all are — as in ST1 (a horizontal
  // force listed with span loads is not the minimum of Mz;
  // `CHARG 2 0.` competes with 0 where CHARG 2 is not null, not where it
  // is); the reactions keep them.
  bool skip_null = false;
  std::vector<const CaseResult*> raw;
};

inline EnvSource env_source(const CaseResult& r, double coef = 1.0) {
  return EnvSource{&r, coef, coef, {}, false, {}};
}
inline EnvSource env_source_range(const CaseResult& r, double cmin, double cmax) {
  return EnvSource{&r, cmin, cmax, {}, false, {}};
}
// Wrap a stored envelope (its two bound pseudo-cases) as ONE alternative-set
// source. The Envelope must outlive the call.
EnvSource env_source_of_env(const Envelope& e, double cmin = 1.0, double cmax = 1.0);

// One envelope point: the two extremes of the privileged quantity, the box
// corner behind each (coefficient and alternative index per source) and the
// concomitant station state assembled from that same corner.
struct EnvPoint {
  Id bar = 0;
  std::size_t station = 0;
  double x = 0;  // abscissa [m]
  double max = 0, min = 0;
  std::vector<double> coef_max, coef_min;
  std::vector<std::size_t> alt_max, alt_min;
  Station st_max, st_min;
};

// Full-state run at ONE privileged location.
struct EnvRun {
  double max = 0, min = 0;
  std::vector<double> coef_max, coef_min;
  std::vector<std::size_t> alt_max, alt_min;
  CaseResult state_max, state_min;  // whole consistent companion states
  bool has_reac = false;
  double reac_max = 0, reac_min = 0;
  std::vector<double> reac_coef_max, reac_coef_min;
  Vec6 reac_state_max = Vec6::Zero(), reac_state_min = Vec6::Zero();
};

// The privileged scalar of one stored result (no coefficients applied).
double env_value(const Model& model, const EnvSpec& spec, const CaseResult& r);

// Envelope at the spec's location: extremes of the privileged quantity over the
// coefficient box, plus the ENTIRE companion state on each side. The state is
// `combine_results` of every source taken with ITS chosen corner coefficient
// and alternative — one consistent load state, not a mix of per-field extrema.
EnvRun env_run(const Model& model, const EnvSpec& spec, const std::vector<EnvSource>& sources);

// Same envelope run over EVERY study section of every bar (ST1 evaluates the
// privileged component section by section). Node quantities (Depla) are a typed
// reject here — they have no section — and bars where the quantity does not
// exist are skipped: no Winkler bed for PRESS, no rule for CONTRG (p.101).
std::vector<EnvPoint> env_run_sections(const Model& model, const EnvSpec& spec,
                                       const std::vector<EnvSource>& sources);

}  // namespace keta

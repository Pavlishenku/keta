// Staged construction by incremental superposition: a PhasePlan is walked
// phase by phase; every increment is solved LINEARLY on the currently active
// substructure and accumulated component-wise into a full-layout master state.
// There is no global re-solve — the locked-in stresses and redistribution come
// from what structure each increment sees.
//   * facet shells: ACTIVE FROM THE FIRST PHASE by design — a plan names BARS
//     only; the slabs belong to the initial layout (cast at t = 0), every
//     working substructure carries the full shell set, and the per-phase
//     shell resultants are accumulated tensor-wise with the principals
//     RECOMPUTED (staged casting of shells would need its own
//     stress-free-cast bookkeeping and is out of scope). Facets of aging
//     concrete CREEP and SHRINK between the dates: Boltzmann superposition of
//     their stress events with the law itself, each event's creep field the
//     whole elastic strain field of the event (membrane, curvature,
//     transverse shear) imposed exactly as an element field
//     (Load::ShellImposed), the free shrinkage as the in-plane field
//     eps_cs (x, y); the implicit step of the section parts (E_eff, the
//     step's stress increment an event at its start). Orthotropic facets and
//     facets on a bed are refused when they creep;
//   * shell tendons (Model::shell_tendons, TensionAction::shell): jacked by
//     fractions of their P0 profile (K = 0 within a group), bonded — every
//     recorded increment strains them by the mean concrete strain along
//     their length (Ep Ap of it lost or gained), the loss loads lose again,
//     iterated; between the dates they lose by creep and shrinkage THROUGH
//     the structure and by EN 1992-1-1 relaxation (equivalent time on the
//     current stress). A loss scales the whole profile (PhaseResult::
//     shell_tendons). The shell tendons' losses do not feed back into the
//     bar tendons' K = 1 passes (a mixed model of bar and shell tendons);
//   * activation: cast stress-free in the current geometry. RE-activating a
//     previously deactivated bar is legal and returns it stress-free at its
//     new date: its accumulated state was released at deactivation
//     (spring-back) and its creep memory wiped with it;
//   * deactivation: spring-back — the removed bar's accumulated nodal forces
//     (and a removed support's accumulated reaction) are released onto the
//     remaining structure as one increment;
//   * case loads: applied on the active substructure of their phase;
//   * tensioning: tendons lowered at their fraction increments (simultaneous
//     within a phase = the jack compensates their own shortening, K = 0);
//     after every later increment, previously-tensioned tendons with live
//     Ep/Ap lose Ep Ap d_eps_mean (EN 1992-1-1 §5.10.5.1 K = 1, one pass);
//   * detensioning: a TensionAction whose target fraction is BELOW the
//     applied one (negative = full release, the ST1 `DETENDRE` marker)
//     releases the corresponding share of the CURRENT EFFECTIVE force — the
//     spring-back symmetric to tensioning (an inert tension -> detension
//     cycle cancels exactly);
//   * creep — decided by the MATERIALS and the live tendons, never by a flag
//     of the plan: between phase
//     dates each active concrete bar (Material::concrete / concrete_law)
//     receives the imposed-strain increment (creep of the accumulated stress
//     + free shrinkage); a tendon on bars whose material does not age (steel)
//     still relaxes; every LIVE tendon with Ap > 0 additionally takes the
//     COUPLED DEFERRED LOSS of the interval — EN 1992-1-1 eq. 5.46 applied
//     incrementally (shrinkage + creep at tendon level from the same engine +
//     relaxation pro-rata of the time since tensioning, full interaction
//     denominator), applied to the cable force BEFORE the K = 1
//     redistribution of the concrete increment; then K = 1 losses run again.
//     `PhasingOptions::creep_engine` selects HOW the creep increment is
//     computed — the single AAEM step, or the rate-type Kelvin chains of
//     `keta/hereditary.h` (see `CreepEngine` below).
#pragma once
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>

#include "keta/linsolve.h"
#include "keta/prestress.h"
#include "keta/soilzones.h"

namespace keta {

struct TendonPhaseForce {
  Id tendon = 0;
  double jacked_to = 0;   // fraction * P0 applied so far [N]
  double loss = 0;        // cumulative elastic (K=1) + deferred (eq. 5.46) loss [N]
  double effective = 0;   // jacked_to - loss
  bool bed = false;       // pretension strand still held by its bed (no force on the structure)
  // PhasingOptions::local_losses: the local loss field (chain abscissa s,
  // force lost there) to take off `effective / jack_force` times the
  // tendon's force profile — empty otherwise.
  std::vector<std::pair<double, double>> local_loss;
};

struct PhaseResult {
  std::string name;
  double time = 0;
  CaseResult cumulated;    // full layout (inactive bars carry zeros)
  CaseResult incremental;
  std::vector<TendonPhaseForce> tendons;
  // the shell tendons acted on (keta/shell_tendon.h, TensionAction::shell):
  // jacked_to = applied fraction x P0, effective = the fraction left x P0
  // (the losses scale the whole profile: the force at a point is
  // effective / P0 times the P0 profile there)
  std::vector<TendonPhaseForce> shell_tendons;
  // Composite bars (keta/composite.h): the cumulated resultants of each cast
  // part per station of the bar (N, My, Mz about the part's own centroid;
  // fibre stress N/A + My z'/Iy - Mz y'/Iz in the part's axes). The bar's
  // own record in `cumulated` is about its NODE LINE.
  struct Part {
    Id bar = 0;
    int part = 0;
    std::vector<Vec3> stations;
  };
  std::vector<Part> parts;
  // --- appended: PhasingOptions::origin_ledgers ---
  // The cumulated state split by the ORIGIN of its increments (sum =
  // `cumulated`): [0] S(G) — case loads, the structural creep and shrinkage
  // (the redistribution practice carries with the permanent effects), and
  // the releases; [1] S(Po) — tensioning, detensioning and the instantaneous
  // (K = 1) losses the non-deferred increments cause; [2] S(dP) — the
  // deferred tendon losses (eq. 5.46: creep, shrinkage, relaxation) and the
  // elastic re-losses they, and the shortening of the tendons by the
  // structural creep, cause (with their sign: a loss is a negative prestress
  // effect). Empty unless the option is set.
  std::vector<CaseResult> origins;
  // The primary (isostatic) prestress of the tendons acting on the structure
  // at the end of the phase, at the stations: N, Vy, Vz, My, Mz (ETAT
  // PRECONTRAINT); end forces, torsion, displacements and reactions zero.
  CaseResult isostatic;
  // The first event whose increment mixed the origins (a bar, link or
  // support release, a support modification or jacking): the per-origin
  // attribution is no longer exact after it, and the weighted states refuse.
  std::string origins_invalid;
  // --- appended: catenary cables (Phase::cable_actions) ---
  struct CableState {
    Id bar = 0;
    bool active = false;
    double L0 = 0, H = 0, T1 = 0, T2 = 0;  // unstrained length, horizontal and end tensions
  };
  std::vector<CableState> cables;  // every cable bar the phasing has seen
  // a rupture of this phase with a dynamic factor phi: (phi - 1) x the static
  // release increment — the accidental peak is cumulated + this
  CaseResult cable_rupture_peak;
  // --- appended: PhasingOptions::keep_structures ---
  // The active substructure at the end of the phase, as the increments were
  // solved on it (without its catenary cables, which ride their own overlay).
  std::shared_ptr<const Model> structure;
};

struct PhasedResults {
  std::vector<PhaseResult> phases;
  // Number of distinct NNLS chain fits performed (Kelvin engine only). Fits
  // are shared through a (material, quantized loading age) cache, so two bars
  // loaded at the same age cost ONE fit — exposed so a driver can
  // audit the cost of a phasing run.
  int kelvin_fits = 0;
  // Factor-reuse tallies of the PhasedFactorCache after this solve, summed
  // over every substructure slot (cumulative over the cache's lifetime).
  FactorCache::Stats cache_stats;
};

// One factorization slot per distinct ACTIVE SUBSTRUCTURE. A plan alternating
// active sets A, B, A reuses A's factor on its second visit instead of
// thrashing one shared slot, and a caller that recomputes the same plan (a
// load edit, a creep option) holds the cache across solves so every unchanged
// substructure skips its factorization.
// The key — sorted active bar ids, sorted active support node ids — is a pure
// performance hint: FactorCache exact-compares the assembled K before any
// reuse, so a stale or colliding slot only costs a re-factorization, never a
// wrong number. Results are bit-identical to the uncached path.
// Eviction: none. The map grows with the number of distinct active sets the
// owner has ever seen (a handful of phases, edited occasionally); the owner
// bounds its lifetime (a node's cache dies with its document).
class PhasedFactorCache {
 public:
  using Key = std::pair<std::vector<Id>, std::vector<Id>>;  // (bars, support nodes)
  FactorCache& slot(const Key& key) { return slots_[key]; }
  FactorCache::Stats stats() const;  // summed over the slots
  std::size_t size() const { return slots_.size(); }
  void clear() { slots_.clear(); }

 private:
  std::map<Key, FactorCache> slots_;
};

// Which time-dependent engine produces the imposed-strain increment applied
// between two phase dates.
//
//  * Aaem   — one age-adjusted effective-modulus step per interval:
//             d_eps = phi(a2, a1) * sigma / E, with a1/a2 the concrete ages at
//             the interval ends. Cheap and stateless, but NOT telescoping: the
//             creep coefficient of the second interval is the one of a load
//             applied at a1, so splitting one interval in two OVERSTATES the
//             total creep (EN 1992-1-1 Annex B knows only phi(t, t0)).
//  * Kelvin — the rate-type Kelvin chains of `keta/hereditary.h`, fitted on
//             J(t, t0) = 1/E + phi(t, t0)/E (the compliance the model's own
//             constant `Material::E` implies). Each chain carries one internal
//             strain per unit and per generalized stress component
//             (N/A, My/Iy, Mz/Iz), so the delayed strain telescopes EXACTLY
//             over any subdivision of the interval.
//             ONE CHAIN PER LOADING EVENT (Boltzmann superposition — McHenry
//             1943; Bazant & Wu, "Rate-type creep law of aging concrete based
//             on Maxwell chain", 1974): at every creep interval where a bar's
//             generalized stress has changed by more than 1e-8 (relative to
//             the larger of the new and carried stress, per component) since
//             the last event, the INCREMENT is handed to a NEW chain fitted at
//             the age of that event; sub-threshold drift is folded into the
//             latest chain. A stress applied late therefore creeps with the
//             compliance of the age it was applied at — the frozen-age error
//             (~13 % on a two-step history) is gone,
//             what remains is the chain fit error. Fits are shared through a
//             (material, quantized loading age) cache — 4096 log-steps per
//             decade, a 0.056 % age resolution far below the fit error — so
//             all bars loaded at one phase date cost ONE NNLS fit
//             (`PhasedResults::kelvin_fits` reports the count).
//             The stress is taken PIECEWISE CONSTANT over each interval, at its
//             value at the interval start, exactly like the AAEM path (the
//             redistribution of the interval itself is explicit, one pass).
// Creep engines of a phasing:
//   Exact  — Boltzmann superposition with the law itself: every stress
//            increment ds_j applied at the age t_j adds, over [a1, a2],
//            ds_j (phi(a2, t_j) - phi(a1, t_j)) / E_ref(t_j) — no fit, no
//            threshold (the exact linear-creep response to the piecewise
//            constant stress history);
//   Aaem   — one age-adjusted step per interval (historical);
//   Kelvin — Boltzmann superposition of the stress increments (EN 1992-1-1),
//            one Kelvin chain per loading event;
//   St1    — ST1 2.24's incremental method (manual annex B3.1), identified on
//            ST1's own results and run FIBRE BY FIBRE: 2 fibres
//            in a plane section, 4 in space, at the radii of gyration. With
//            phi(t, t0) = A(t0) f(t - t0) (A the final value, f the time
//            function), a fibre keeps the amplitude S = sum dsigma_j A(t_j) /
//            E_ref(t_j) of its stress steps and its creep strain e; at a step,
//            with r = e / S:
//              loading   (e = 0 or 0 < r < 1): e = S f(te), the equivalent
//                         duration te set so that e is continuous;
//              recovery  (r >= 1, r <= 0, or a compression creep going
//                         into tension): e = e1 + (S - e1) f(tr), from the
//                         strain e1 of the step, tr likewise continuous across
//                         later steps; a recovery towards S goes on while e
//                         lies between e1 and S;
//              reversal  (BPEL, tension creep e > 0 under compression, S < e):
//                         a loading of S from the equivalent duration
//                         f(te) = e / (2e - S) (exact); if S stays
//                         positive, a recovery towards S from the next step;
//              tension   (EC2 family, no tension creep): a virgin fibre keeps
//                         no amplitude; a crept one recovers towards S while
//                         S keeps the sign of its creep, towards 0 once it
//                         changes; the BPEL family creeps in tension
//                         (st1_tension_creep_bars).
//            With 4 fibres (SPATIALE) the part of their creep strains a plane
//            section cannot follow is restrained: a self-equilibrated stress
//            -E (e - plane fit) is added to each fibre.
//            The final values are the superposition's; right after a step the
//            creep is much slower (e.g. 0.03 mm in a day, 0.24 by
//            superposition). For a law whose time function depends on the
//            loading age (EC2 2027) f is taken at the fibre's first loading
//            age (an effect of the order of 5 %).
enum class CreepEngine { Aaem, Kelvin, St1, Exact };

struct PhasingOptions {
  int stations_per_bar = 11;
  CreepEngine creep_engine = CreepEngine::Exact;
  // Kelvin chain size (Kelvin engine only). The retardation times always span
  // 1e-2 .. 1e5 days (decades = 7 / (n - 1)), so more units = finer spacing:
  // 8 units fit J to ~1.2 %; 16 to < 0.1 % beyond one day of loading, ~1 %
  // in the first hours.
  int kelvin_units = 8;
  // Cooperative cancel + progress (keta/hooks.h, not owned): polled at entry
  // and before every phase ("phase i/n"). Null: unhooked.
  const SolveHooks* hooks = nullptr;
  // --- appended ---
  // Creep along the bar (ST1's per-section fibres, manual annex B3.1): the
  // creep strain is computed at EVERY station from the stress there and
  // lowered as the imposed-strain field it is (concentrated dislocations at
  // the segment midpoints; a least-squares linear field on a bar with its own
  // kernel), and the tendon terms (K = 1 strain, eq. 5.46 creep / shrinkage /
  // section factors) are length-weighted over the stations of the tendon's
  // ANCHORED length, its relaxation base being the mean stress of its
  // friction / draw-in profile. false (default): the mid-bar scheme — one
  // stress per bar, a uniform imposed strain, mid-bar tendon terms over the
  // whole chain, the jack stress as relaxation base — exact only for bars
  // short against the stress variation.
  bool creep_along_bar = false;
  // Sub-steps per decade of the elapsed time of a creep interval (exact
  // engine; 0 = one step per interval, the stress held at its start value).
  int creep_substeps_per_decade = 4;
  // Bonded-tendon losses solved through the structure (the loss loads lose
  // again in every tendon, iterated; creep under the loss by the exact
  // engine's sub-steps; EC2 relaxation by equivalent time on the current
  // stress) instead of the EN 1992-1-1 eq. 5.46 section formula per interval
  // (factor 1 + rho (1 + 0.8 dphi), relaxation x 0.8) and a single K = 1 pass.
  bool exact_losses = true;
  // Station abscissae per bar [m] (keta/statics.h StaticOptions::stations_abs,
  // not owned); bars absent from it use `stations_per_bar`.
  const std::map<Id, std::vector<double>>* stations_abs = nullptr;
  // A bar activated between nodes that have already moved is born straight
  // on their displaced positions (stations start from the chord of the end
  // translations), so its station displacements are totals, continuous with
  // the nodes in translation (ST1). false: the stations count only the
  // increments since activation.
  bool born_deformed = false;
  // Instantaneous modulus of an aging concrete bar at its age [days] on the
  // date of each increment (ST1: the modulus E_j of the age at loading —
  // Exemple5, a release at the end of a heat cure deflects as E(fcmp)), in
  // the model's units. Null (default): every increment uses Material::E.
  std::function<double(Id bar, double age_days)> modulus_at;
  // Tendon losses POINT BY POINT (ST1): the K = 1 elastic losses, the
  // deferred ones (eq. 5.46) and the relaxation are evaluated at every station
  // of the tendon's bars from the concrete strain and stress THERE, and kept
  // as a loss field along the tendon (TendonPhaseForce::local_loss) — a
  // pretension strand loses its plateau shortening, not the mean of its
  // transmission zones (Exemple5: 394.5, the mean gave 330). Needs
  // `creep_along_bar`. false: one loss fraction of the whole profile.
  bool local_losses = false;
  // EXEC PHASAGE PSE (ST1 manual p.192): the bars of this map (not owned)
  // whose CARA PSE carries a PLASTIQUE law rest on a pseudo-plastic soil. They
  // must carry the exact bed (keta/bef.h) on Model::tapered; the others are
  // ignored. Each increment iterates on the state of every force component of
  // their bed, decided as ST1 does on the MEAN pressure k.u along the bar
  // (manual p.60): elastic, capped at the threshold (a constant line force) or
  // lifted off (POS / NEG). The soil keeps no memory (r = clamp(k.u) of the
  // total displacement since the bar's activation): a change of state releases
  // the difference of soil force as a line load on the bar. The torsion
  // component (KMX) stays elastic. Null: every bed is linear.
  const SoilZoneMap* pse = nullptr;
  // CreepEngine::St1: the bars that take it (not owned); the others creep by
  // superposition (CreepEngine::Kelvin). Null: every bar.
  const std::set<Id>* st1_creep_bars = nullptr;
  // CreepEngine::St1: the bars whose law creeps in tension (the BPEL family;
  // not owned). The others (EC2 family) take no tension creep. Null: none.
  const std::set<Id>* st1_tension_creep_bars = nullptr;
  // CreepEngine::St1: the fibre layout of a section — 0: two fibres along y
  // (PLANE), 1: two along z (GRILL), 2: four (SPATIALE).
  int st1_fibre_layout = 0;
  // Workers for the station recovery of each increment
  // (StaticOptions::bar_threads, keta/parallel.h): 1 = sequential (default),
  // 0 = hardware concurrency. Bit-identical results for any count.
  int threads = 1;
  // --- appended ---
  // Keep the cumulated state split by origin and the isostatic prestress
  // (PhaseResult::origins, isostatic; weighted_state, weighted_envelope,
  // hyperstatic_prestress below). Costs one accumulation per origin and per
  // increment and three more states per recorded phase.
  bool origin_ledgers = false;
  // --- appended ---
  // Record in PhaseResult::structure the active substructure of the end of
  // each phase (bars, supports, links, sections and moduli of its date): the
  // structure a saved tangent linearizes (keta/tangent.h). One model copy per
  // recorded phase.
  bool keep_structures = false;
};

// ---- weighted (characteristic) states from the origin ledgers ----
// ETAT PONDERE: cumulated + (w_g - 1) S(G) + (w_po - 1) S(Po) + (w_dp - 1)
// S(dP) — weights (1, 1, 1) return `cumulated` exactly; BPEL 91 rev. 99 art.
// 4.1.2: P1 = 1.02 P0 - 0.80 dP is (1, 1.02, 0.80), P2 = 0.98 P0 - 1.20 dP is
// (1, 0.98, 1.20) (S(dP) carries the minus sign of the losses); EN 1992-1-1
// §5.10.9: (1, r, r). Over node displacements, reactions, bar end forces and
// stations, link forces (not the facet shells). Throws when the phase was
// run without the ledgers or after an event that mixed the origins.
CaseResult weighted_state(const PhaseResult& phase, double w_g, double w_po, double w_dp);
// The favourable / unfavourable coefficient pairs of ETAT PONDERE (KGFAV
// KGDEF KPOFAV KPODEF KPEFAV KPEDEF | KPTFAV KPTDEF): the prestress either
// decomposed (kpo on S(Po), kpe on S(dP), BPEL) or total (kpt on S(Po) +
// S(dP), EN 1992-1-1 §5.10.9).
struct WeightingCoefficients {
  std::pair<double, double> kg{1, 1}, kpo{1, 1}, kpe{1, 1}, kpt{1, 1};
  bool total = false;
  static WeightingCoefficients bpel() { return {{1, 1}, {0.98, 1.02}, {0.80, 1.20}, {1, 1}, false}; }
  static WeightingCoefficients ec2(double r_inf, double r_sup) { return {{1, 1}, {1, 1}, {1, 1}, {r_inf, r_sup}, true}; }
};
// The envelope (min, max) of the weighted states: every scalar component
// takes, independently and per origin, the member of each pair that
// minimises / maximises it. Equal pairs give weighted_state exactly.
std::pair<CaseResult, CaseResult> weighted_envelope(const PhaseResult& phase, const WeightingCoefficients& c);
// The hyperstatic (secondary) prestress effect at the stations: S(Po) +
// S(dP) - isostatic (EN 1992-1-1 §5.10.7 (2)): zero on a statically
// determinate structure, linear between the supports of a continuous beam.
CaseResult hyperstatic_prestress(const PhaseResult& phase);

// Two supports with the same node and restraints, field by field.
bool equal_support(const Support& a, const Support& b);

PhasedResults solve_phased(const Model& model, const PhasePlan& plan,
                           const PhasingOptions& opt = {});
// Same, against a caller-owned cache persisting across solves.
PhasedResults solve_phased_cached(const Model& model, const PhasePlan& plan,
                                  const PhasingOptions& opt, PhasedFactorCache& cache);

// The state of a phasing between two phases, run phase by phase and carried
// EXACTLY from one plan to the next: `solve_phased` is one session run on one
// plan; a driver that continues a phasing (ST1 `SUITE PHASAGE`: a chain
// resumed from the end state of another) runs the common prefix ONCE and
// forks the session for every continuation instead of replaying the prefix
// from the first date.
//
// State carried: the master (cumulated) state, the active substructure and
// its supports, every bar's accumulated nodal forces (spring-back), the
// supports' accumulated reactions, the tendon states (fractions, losses,
// local loss fields, relaxation clocks), the creep memories of every station
// (superposition events, Kelvin chains, ST1 fibres), the casting dates and
// the pseudo-plastic soil states — all of it keyed by node / bar / tendon
// ids. The factorisation cache, the tendon profiles and the loss bases are
// derived (pure functions of the model) and are rebuilt or inherited.
//
// fork(model2, opt2): a new session holding this one's state on `model2`,
// a model that CONTAINS this one's footprint — the same nodes (same ids,
// coordinates, relative order), the same bars (same ends, section, material
// constants, kernels; relative order kept), the same supports, the same
// tendons for those already acted on — and may add bars, tendons and load
// cases the continuation uses. Every increment of the prefix was solved on
// the active substructure only, whose free-DOF numbering, assembled system
// and recovery are the same on both models, so the forked session continues
// BIT-IDENTICALLY to a session that would have run the whole chain on
// `model2` from the start. What the engine
// cannot compare — the time-dependent laws (functions), `modulus_at`, the
// exact kernels of CARA VAR bars — is the caller's contract: they must be
// the same laws. A model that does not contain the footprint, a changed
// station layout, a changed support set or a changed option the increments
// depend on is a TYPED REJECT (`fork_reject` gives the reason without
// throwing), never a silent approximation.
//
// Lifetime: the session copies `model` and `opt`; the objects `opt` points at
// (stations_abs, pse, st1_creep_bars, st1_tension_creep_bars, hooks) are not
// owned and must outlive the session.
class PhasingSession {
 public:
  // A pristine session on `model`: `active_from_start` are the bars active
  // before the first phase (solve_phased: the bars no phase activates; an
  // ST1 phasage starts empty), cast at t = 0. `cache`: a caller-owned
  // factorisation cache (null: the session owns one).
  PhasingSession(const Model& model, const PhasingOptions& opt,
                 const std::set<Id>& active_from_start, PhasedFactorCache* cache = nullptr);
  ~PhasingSession();
  PhasingSession(PhasingSession&&) noexcept;
  PhasingSession& operator=(PhasingSession&&) noexcept;

  // Run the phases of `plan` in order from the current state. The results of
  // the phases run are appended to results().
  void run(const PhasePlan& plan);

  // The results of the phases this session ran itself (a forked session
  // starts empty: the prefix results belong to the session it came from).
  const std::vector<PhaseResult>& results() const;
  std::vector<PhaseResult> take_results();
  // Phases run so far along the chain, the inherited ones included.
  std::size_t phases_run() const;
  // The date of the last phase run.
  double time() const;

  // Links (elastic / rigid ids) inactive at the start of the session: a
  // phase activates them later (solve_phased: the links some phase of the
  // plan activates). Before the first run only.
  void defer_links(const std::set<Id>& elastic, const std::set<Id>& rigid);
  // Section parts (bar id, part index) cast later by a phase (solve_phased:
  // those some phase of the plan casts). Before the first run only.
  void defer_parts(const std::set<std::pair<Id, int>>& later);

  // See the class comment. Throws Error with the reason of a reject.
  std::unique_ptr<PhasingSession> fork(const Model& model, const PhasingOptions& opt) const;
  // The reason fork(model, opt) would be refused; empty when it is exact.
  std::string fork_reject(const Model& model, const PhasingOptions& opt) const;

  // The friction / draw-in force profile of a tendon of the model, computed
  // once per session (inherited by a fork): what the drivers read back for
  // the cable tables.
  const ForceProfile& force_profile(Id tendon);

  int kelvin_fits() const;                // see PhasedResults::kelvin_fits
  FactorCache::Stats cache_stats() const;
  const Model& model() const;

  struct Impl;

 private:
  explicit PhasingSession(std::unique_ptr<Impl> impl);
  std::unique_ptr<Impl> impl_;
};

}  // namespace keta

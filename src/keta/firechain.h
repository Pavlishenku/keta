// Fire chain: hot mechanical state of a heated cross-section, and the ST1
// `INCENDIE` driver (a preloaded frame re-solved with degraded, thermally
// pre-strained members).
//
// The chain, following ST1 chapter 12 and its annex B5 (manual p213-229, p380-391):
//
//   MAT ... E FEU EC2   ->  FireConcrete   (fck, TRACTION, thermal properties)
//   ACIER / ARMATURE    ->  RebarSteel / RebarLayer
//   EXPOSITION FEU      ->  Exposure       (heat2d.h)
//   SECTION             ->  FireSection    (heat2d.h)
//   FEU  DUREE / TEMP   ->  FireProblem::duration / t_init
//   INCENDIE  FEU j (LOCAL) CHARG k  ->  FireProblem::local / preload_case
//   ANALYSE FEU         ->  FireAnalysisOptions
//   EXEC INCENDIE       ->  solve_fire()
//
// Mechanics (annex B5 p382-384)
// ----------------------------
// Every fibre is treated as thermo-elastic with a SECANT modulus taken on its
// own EN 1992-1-2 law at its own temperature:
//
//   eps_e,i = eps_i - eps_th,i          sigma_i = f_i(eps_e,i)      E_i = sigma_i / eps_e,i
//
// so the member rigidities and the thermal pre-stress resultants are
//
//   EA    = int E dS                        N_th   = int E eps_th dS
//   EI_y  = int E (z - z_c)^2 dS            M_th,y = int E eps_th (z - z_c) dS
//   EI_z  = int E (y - y_c)^2 dS            M_th,z = - int E eps_th (y - y_c) dS
//
// about the STIFFNESS centroid (y_c, z_c) = (int E y dS, int E z dS) / EA, which
// fire moves away from the cold centroid (ST1 `Y VAR`, the default; `Y CTE`
// pins it at the cold position). The equivalent imposed strain handed to the
// frame element is then eps_0 = N_th/EA and kappa = M_th/EI, which is exactly
// keta's `Load::Prestrain` — a restrained bar therefore develops
// N = -EA_hot * eps_0, the classical restrained-expansion force.
//
// Frame iteration (`solve_fire`)
// ------------------------------
// 1. the section heat problem is solved once per fire bar (heat2d.h);
// 2. every fire bar's Section is replaced, in a working copy of the model, by a
//    generic section whose A / I reproduce the hot EA / EI against the bar's
//    material modulus (and whose cloned material density is rescaled so that
//    self-weight is UNCHANGED by the substitution);
// 3. the preload case (ST1 `CHARG ... COMB`) plus the thermal `Prestrain` loads
//    are solved; the mean member forces over the bar give a new strain plane;
// 4. repeat until the nodal displacement field stops moving (`tol`) or
//    `max_iter` is reached.
//
// NON-CONVERGENCE IS A TYPED FAILURE, AND THAT FAILURE IS THE COLLAPSE SIGNAL.
// Annex B5 p388: "L'issue du calcul peut etre une raideur axiale nulle signalant
// la rupture totale de la section. La raideur peut aussi devenir si faible que la
// methode de Newton ne parvient pas a converger". keta never returns a degraded
// "best effort" state for a structure that failed to find equilibrium.
//
// Known limitations (documented, never silent)
//   * SECANT approach: the frame is solved with secant member rigidities, so the
//     result is an equilibrium state, not a load-deflection history; no unloading
//     path, no plastic memory, no creep at temperature (EN 1992-1-2 transient
//     creep is only implicitly inside the §3.2.2.1 law).
//   * no SPALLING model — the user removes the spalled thickness by hand
//     (manual p219).
//   * mean-of-the-bar deformation: as ST1 (manual p213), one strain plane per
//     BAR, taken from the mean member forces, so a span must be discretised into
//     several bars.
//   * the thermal problem is uncoupled from the mechanical one (no cracking
//     feedback on conduction) and prestressing tendons are not re-tensioned by
//     the fire in this driver.
#pragma once
#include "keta/heat2d.h"
#include "keta/statics.h"

namespace keta {

// ---------------------------------------------------------------------------
// EN 1992-1-2 §3.2.2.1 concrete stress-strain relationship at temperature
// ---------------------------------------------------------------------------
// sigma(theta, eps) = 3 eps f_c,theta / ( eps_c1,theta (2 + (eps/eps_c1,theta)^3) )
// up to eps_c1,theta, then a linear descending branch down to zero at eps_cu1.
// Compression is NEGATIVE here (keta's fibre sign convention). The tangent at the
// origin is E(theta) = 3 f_c,theta / (2 eps_c1,theta), which is the value ST1
// assigns through `E FEU EC2` (manual p215).
struct ConcreteFireLaw {
  double fc = 30e6;     // f_c,theta [Pa]
  double ec1 = 0.0025;  // eps_c1,theta
  double ecu1 = 0.02;   // eps_cu1,theta
  // ST1 `TRACTION OUI` re-activates the origin tangent on the tension side. The
  // EN default is `TRACTION NON` — no tensile strength, i.e. a cracked analysis
  // (manual p215).
  bool tension = false;

  double e_origin() const { return 3.0 * fc / (2.0 * ec1); }
  double sigma(double eps) const;
  double tangent(double eps) const;
  double secant(double eps) const;  // sigma/eps, e_origin() in the limit eps -> 0
};

// ---------------------------------------------------------------------------
// Reinforcement — ST1 `ACIER` / `ARMATURE` (manual p223-224)
// ---------------------------------------------------------------------------
struct RebarSteel {
  double diameter = 0.0;  // DIAMETRE phi [m]
  double area = 0.0;      // SECTION as [m^2]; 0 => pi phi^2 / 4
  SteelClass cls = SteelClass::N;                   // CLASSE <N, X>
  SteelForming forming = SteelForming::ColdWorked;  // LAMINE <CHAUD, FROID>; ST1 default: cold
  double E = 200e9;    // E   [Pa]
  double fyk = 500e6;  // FYK [Pa]

  double bar_area() const;  // as, or pi phi^2/4
};

struct RebarLayer {
  // ENROBAGE SUP/INF measures the cover to the rebar CENTROID (manual p224:
  // "les enrobages ne tiennent pas compte des diametres des aciers").
  enum class Placement { CoverSup, CoverInf, Local };
  RebarSteel steel;
  int count = 1;  // NOMBRE
  Placement placement = Placement::CoverSup;
  double cover = 0.0;   // ENROBAGE [m]
  double y = 0, z = 0;  // LOCAL Y / Z [m]
};

// A rebar resolved into the section frame.
struct Rebar {
  double y = 0, z = 0, area = 0;
  RebarSteel steel;
};

// SUP/INF layers are spread evenly across the section width; LOCAL layers are
// stacked at their given point. Throws when a cover puts a bar outside the
// section or when SUP/INF is used on a PERIMETRE section (manual p223 reserves
// ENROBAGE for DALLE / RECTANGLE).
std::vector<Rebar> resolve_rebars(const FireSection& section,
                                  const std::vector<RebarLayer>& layers);

// ---------------------------------------------------------------------------
// Fire material + fire bar — ST1 `MAT ... E FEU EC2` / `CARA ib SECTION i ARMATURE j`
// ---------------------------------------------------------------------------
struct FireConcrete {
  double fck = 30e6;      // FCK [Pa] (characteristic cylinder strength at 20 C)
  bool tension = false;   // TRACTION OUI / NON
  ThermalMaterial thermal;
};

struct FireBarSpec {
  Id bar = 0;
  FireSection section;
  FireConcrete material;
  std::vector<RebarLayer> rebars;
};

// ---------------------------------------------------------------------------
// Fibre discretisation of a heated section
// ---------------------------------------------------------------------------
struct HotFibers {
  struct Fibre {
    double y = 0, z = 0, area = 0;
    double theta = 20.0;   // temperature [degC]
    double eps_th = 0.0;   // EN 1992-1-2 §3.3.1 free thermal elongation
    int law = 0;           // index into `concrete` or `steel`
    bool is_concrete = true;
  };
  std::vector<Fibre> fibres;
  std::vector<ConcreteFireLaw> concrete;  // one per 1 C temperature bin used
  std::vector<UniaxialLaw> steel;

  double sigma(const Fibre& f, double eps_e) const;
  double tangent(const Fibre& f, double eps_e) const;
  double secant(const Fibre& f, double eps_e) const;
  double max_rebar_temp() const;
};

// Build the fibre set of a fire bar at one heat snapshot. Concrete fibres are the
// heat cells; rebars are added at their resolved position and take the
// temperature of the hottest concrete cell within phi/2 of their centre
// (manual p223: the bar diameter selects the concrete cell that sets its
// temperature).
HotFibers build_hot_fibers(const FireBarSpec& spec, const HeatResult& heat, int snap = -1);

// ---------------------------------------------------------------------------
// Hot section state
// ---------------------------------------------------------------------------
struct HotSectionState {
  Id bar = 0;
  double EA = 0;             // int E dS [N]
  double EIy = 0, EIz = 0;   // about (yc, zc) [N m^2]
  double yc = 0, zc = 0;     // stiffness centroid in the section frame [m]
  double N_th = 0;           // int E eps_th dS [N]
  double M_th_y = 0, M_th_z = 0;  // [N m]
  double eps_th = 0;         // N_th / EA — equivalent imposed axial strain
  double kappa_th_y = 0, kappa_th_z = 0;  // M_th / EI [1/m]
  double theta_max = 0, theta_mean = 0, theta_rebar_max = 0;  // [degC]
  bool secant = true;        // rigidities from the secant (true) or tangent moduli
};

// Integrate the hot rigidities for a given strain plane. `strain_plane` is
// (eps0, kappa_y, kappa_z) referred to (yref, zref):
//   eps(y, z) = eps0 + kappa_y (z - zref) - kappa_z (y - yref)
// `ignore_thermal` (ST1 `INCENDIE ... LOCAL`) drops eps_th everywhere: the
// materials stay degraded but the structure feels no thermal deformation.
HotSectionState hot_section_state(const HotFibers& fib, const Vec3& strain_plane, double yref,
                                  double zref, bool secant, bool centroid_var,
                                  bool ignore_thermal);

// Section-level Newton on the strain plane: find (eps0, kappa_y, kappa_z) about
// (yref, zref) whose fibre resultants match `target` = (N, My, Mz), using the
// TANGENT 3x3. Returns the (possibly partially converged) plane — this is the
// tangent half of ST1's `NEWTON_RAPHSON ... PERIODE` alternation.
Vec3 section_newton(const HotFibers& fib, const Vec3& target, const Vec3& start, double yref,
                    double zref, bool ignore_thermal, int max_iter = 25, double tol = 1e-9);

// ---------------------------------------------------------------------------
// INCENDIE driver
// ---------------------------------------------------------------------------
// ST1 `OPTION`. Fire analysis is available in PLANE and SPATIALE only
// (manual p213: "il n'est pas disponible en option grille").
enum class ModelOption { Plane, Grill, Spatial };

struct FireProblem {
  std::vector<Exposure> exposures;  // EXPOSITION FEU i — indexed by FireSection
  std::vector<FireBarSpec> bars;    // the fire-exposed bars
  double duration = 60.0;           // FEU DUREE [minutes]
  double t_init = 20.0;             // FEU TEMP  [degC]
  double dt = 30.0;                 // conduction time step [s]
  Id preload_case = 0;              // INCENDIE ... CHARG k
  bool local = false;               // INCENDIE ... LOCAL
  ModelOption option = ModelOption::Spatial;
};

struct FireAnalysisOptions {
  int max_iter = 500;         // NEWTON_RAPHSON MAX_ITER (ST1 default 500)
  // PERIODE (ST1 default 5): the frame is re-solved with SECANT rigidities every
  // iteration; between two frame solves, `period - 1` section-level TANGENT
  // Newton corrections refine the strain plane at frozen member forces. This is
  // the "les autres iterations se faisant avec des rigidites tangentes" of
  // manual p226, applied where it converges: on the section state, not on the
  // fixed-point operator itself (whose fixed point must be the secant one, cf.
  // annex B5 p382).
  int period = 5;
  bool centroid_var = true;   // Y VAR (default) / Y CTE
  bool keep_hot = false;      // CARA GARDER (keep the degraded section in `hot_model`)
  // Fold the hot neutral-axis shift into the bar's rigid-arm eccentricity
  // (manual p229: "les excentrements a chaud des barres sont pris en compte").
  // Off by default: the shift is always REPORTED in HotSectionState.
  bool fold_centroid_ecc = false;
  // Relative infinity-norm change of the nodal displacements. The plain
  // secant iteration on a cracking reinforced section can lock into a
  // 2-cycle; a stalled residual triggers an under-relaxation of the strain
  // plane (halved at every further stall, floor 1/16) that breaks it, so a
  // non-convergence that survives it is the collapse signal. 1e-3 is the
  // ST1-faithful default, tighten via ANALYSE FEU when a study demands it.
  double tol = 1e-3;
  int stations_per_bar = 11;
};

struct FireOutcome {
  CaseResult result;                    // hot equilibrium state under the preload
  std::vector<HotSectionState> states;  // parallel to FireProblem::bars
  std::vector<HeatResult> fields;       // parallel to FireProblem::bars
  Model hot_model;                      // CARA GARDER: model with the degraded sections
  bool hot_kept = false;
  int iterations = 0;
  double residual = 0;
};

FireOutcome solve_fire(const Model& model, const FireProblem& problem,
                       const FireAnalysisOptions& opt = {});

}  // namespace keta

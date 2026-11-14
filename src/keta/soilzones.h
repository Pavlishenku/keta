// keta soil zones — ST1 `CARA PSE`: zoned Winkler beds and pseudo-plastic /
// one-sided soil.
//
// ---------------------------------------------------------------------------
// What this implements
// ---------------------------------------------------------------------------
// The base kernel (element.cpp) carries ONE uniform Winkler bed per bar
// (`Bar::soil`), consistent over the whole bar length. ST1's `CARA PSE`
// command is richer (manual §2.10, pp. 59-66):
//
//   CARA PSE liste SX sx ... (IZ iz) ...
//      ZONE 1 (KFX kfx) (KFY kfy) (KFZ kfz) (KMX kmx) XL xl
//      ...
//      ZONE n (KFX kfx) (KFY kfy) (KFZ kfz) (KMX kmx)            [p. 59]
//
//   CARA PSE liste SX sx ...
//      PLASTIQUE (POS)(NEG) (KFX ffx)(KFY ffy)(KFZ ffz)(KMX fmx)
//      ZONE 1 (KFX kfx) (KFY kfy) (KFZ kfz) (KMX kmx)            [pp. 59-60]
//
// with (p. 60):
//   kfx/kfy/kfz  "réaction d'un sol élastique suivant l'axe ox/oy/oz local"
//   kmx          "réaction d'un sol élastique à la torsion"
//   ffx/ffy/ffz  "effort seuil plastique du sol élastique pour la barre
//                 suivant l'axe ox/oy/oz local"
//   fmx          "effort seuil plastique ... à la torsion"
//   xl           "abscisse de la fin d'une zone de sol élastique"
//
// Zone geometry (p. 61, p. 66): `XL` is the END abscissa of a zone; the value
// given for the LAST zone is ignored — "les propriétés de la dernière zone
// s'appliquent jusqu'à l'extrémité de la poutre, y compris au-delà de
// l'abscisse XL spécifiée" (p. 61) and "la fin de la dernière zone comprend
// obligatoirement l'extrémité de la barre" (p. 66). `resolve_zones()`
// implements exactly that.
//
// ---------------------------------------------------------------------------
// Units and sign conventions (SI, keta base units)
// ---------------------------------------------------------------------------
//   kfx, kfy, kfz  distributed translational bed [N/m per m = Pa]
//   kmx            distributed rotational bed    [N·m/rad per m]
//   ffx, ffy, ffz  plastic threshold, a LINE FORCE intensity   [N/m]
//   fmx            plastic threshold, a LINE MOMENT intensity  [N·m/m]
//   xl_end         abscissa along the bar, from end i          [m]
//
// The thresholds are FORCES PER UNIT LENGTH, not pressures: they are directly
// comparable to `k·u` (a modulus [Pa] times a displacement [m] is [N/m]).
// This is what the manual's own example does (p. 65): a limit pressure
// `ffy0` [t/m²] is multiplied by the pile diameter `diam1` [m] and the pile
// count before being handed to `KFY` on the `PLASTIQUE` line —
// `PLASTIQUE NEG KFY ffy0*diam1*nbp` — exactly as the elastic modulus
// `KI1(i) = kfil(i)*diam1` is (its own comment: "rigidité sol/pieu ramenée à
// l'unité de longueur de pieu, (t/m)/m1 = t/m2"). Same reduction, same
// per-unit-length family.
//
// The reported reaction intensity `soil_p` / `soil_m` follows the existing
// keta convention (statics.cpp): `p = k · u_local`, i.e. the sign of the local
// displacement component. The ST1 "pression" used by the POS/NEG test carries
// that same sign:
//
//     p_st1(x) = k · u(x) = soil_p(x)
//
// so `POS` keeps the bed only where the bar moved towards +local and `NEG`
// only where it moved towards -local. Both flags, or neither, means
// bilateral. The manual only says "la commande POS et NEG indique que la
// réaction du sol est prise en compte que si la pression est positive et/ou
// négative respectivement" (p. 60); the sign follows from the manual's
// Exemple11: a beam loaded along -y on a `PLASTIQUE NEG KFY 300` bed keeps
// its bed where it moves down, and keta reproduces ST1's deflections there to
// 0.05 % with this sign (the opposite one dropped the whole bed).
//
// ---------------------------------------------------------------------------
// Zoned bed matrix
// ---------------------------------------------------------------------------
// The partial-span bed matrix over a zone [x1, x2] is the consistent
//
//     k_bed = sum_zones  k_zone * INTEGRAL_{x1}^{x2} N(x)^T N(x) dx           (12x12)
//
// with EXACTLY the shape functions the full-bar bed uses in element.cpp:
//   * kfx (DOFs 0, 6) and kmx (DOFs 3, 9): linear   N = [1-xi, xi]
//   * kfy (DOFs 1, 5, 7, 11, l = +L) and kfz (DOFs 2, 4, 8, 10, l = -L):
//     Hermite cubics N = [1-3xi²+2xi³, l(xi-2xi²+xi³), 3xi²-2xi³, l(-xi²+xi³)]
//     (the l = -L on the x-z plane is the theta_y = -dw/dx convention).
// The integral is evaluated by 8-point Gauss-Legendre per zone. The bending
// integrand N_a N_b is a degree-6 polynomial, so a 4-point rule would already
// be exact (an n-point Gauss rule is exact through degree 2n-1); 8 points
// (exact through degree 15) is used as one fixed rule for every component and
// every sub-segment — it costs nothing measurable, it keeps machine-precision
// exactness when a zone is split into many sub-segments for the active set, and
// it leaves room for a modulus varying polynomially inside a zone. A single
// zone spanning [0, L] therefore reproduces the closed-form full-bar bed of
// element.cpp to ~1e-16.
//
// KETA INVARIANT, preserved here: the bed enters `ke_global` ONLY. It is added
// AFTER release condensation (a hinge frees the member moment but the soil
// still reacts to the physical displacement field) and it NEVER enters
// member-force recovery — the soil reaction is recovered separately as
// p = k·u. Consequently the station diagrams INSIDE a bedded bar do not see the
// bed reaction: this is a MESHING-BOUNDED APPROXIMATION, and the member END
// forces are exact. ST1 says the same about its own implementation (p. 60):
// "Lorsque des valeurs fortes de rigidité du sol sont utilisées, si l'on
// définit la poutre par l'intermédiaire d'une seule barre, ST1 peut fournir des
// valeurs erronées liées à la théorie mise en œuvre. Pour éliminer ce problème,
// il convient de discrétiser la barre en plusieurs tronçons."
//
// Also from p. 60: `CARA PSE` forbids reduced shear areas ("il ne faut pas
// donner de valeurs aux sections réduites d'effort tranchant"). keta does not
// enforce that — the bed shape functions stay Hermite (shear-free) while the
// member block may be Timoshenko, which is the same documented mismatch the
// existing full-bar bed already carries.
//
// ---------------------------------------------------------------------------
// Pseudo-plastic driver
// ---------------------------------------------------------------------------
// Elastic-perfectly-plastic bed, solved by an ACTIVE SET on sub-segments:
//
//   r(x) = clamp(k·u(x), -ff, +ff),  and r(x) = 0 where POS/NEG excludes it
//
// Each zone is split into `PseOptions::sub_per_zone` equal sub-segments (10 by
// default). For every sub-segment and every one of the four bed components the
// MEAN reaction intensity over the sub-segment decides the state:
//
//   Inactive : the POS/NEG side test fails       -> no stiffness, no force
//   Capped   : |mean| > ff (ff > 0)              -> no stiffness, a CONSTANT
//              line force of intensity ff opposing the motion
//   Elastic  : otherwise                         -> the consistent bed matrix
//
// then K is reassembled and the case re-solved, until the state vector is
// stable (at most `max_iterations` = 50 sweeps, else a typed Error).
//
// DISCRETIZATION BOUND. The yield front is resolved to ONE SUB-SEGMENT: with a
// zone of length Lz split into n sub-segments the front position carries an
// error of up to Lz/n, and because the criterion uses the sub-segment MEAN, a
// sub-segment straddling the front is taken wholly elastic or wholly capped.
// This is a refinement of ST1, which evaluates plastification on the mean
// pressure over the WHOLE BAR (p. 60: "la plastification du sol n'est calculée
// que sur la pression moyenne le long de la barre. Seule une discrétisation
// importante permet d'avoir des résultats fiables ... La discrétisation doit
// être plus fine, si ces barres sont connectées directement à des appuis
// rigides"). The same advice applies here: the displacement field inside a bar
// is still a single cubic, so a plastified bar must be meshed — sub-segments
// refine the FRONT, not the kinematics.
//
// DEVIATION FROM ST1, deliberate: ST1 allows a single soil layer per bar in the
// plastic case ("on ne peut definir qu'une seule couche de sol par barre",
// p. 60). keta allows a plastic spec over any number of zones — strictly more
// general, and it degrades to ST1's case when one zone is given.
//
// Determinism: zones are resolved in order, bars are visited in model order,
// the per-bar soil map is an ordered std::map, and the quadrature loop is
// serial — the assembled matrix is bit-reproducible run to run.
#pragma once
#include <array>
#include <map>

#include "keta/statics.h"

namespace keta {

// ---------------------------------------------------------------------------
// Input data
// ---------------------------------------------------------------------------

// One ST1 `ZONE` of elastic soil along a bar. `xl_end` is the END abscissa of
// the zone measured from end i [m]; the value carried by the LAST zone of a bar
// is ignored (that zone always runs to the bar end — p. 61, p. 66).
struct SoilZone {
  double kfx = 0, kfy = 0, kfz = 0, kmx = 0;  // [Pa], [Pa], [Pa], [N·m/rad per m]
  double xl_end = 0;                          // ST1 XL [m]
};

// ST1 `PLASTIQUE (POS)(NEG) ...` line. Thresholds are line intensities:
// ffx/ffy/ffz [N/m], fmx [N·m/m]. A threshold <= 0 means "keyword absent" —
// that component is NOT capped (but POS/NEG still applies to it).
struct SoilPlastic {
  bool pos = false, neg = false;
  double ffx = 0, ffy = 0, ffz = 0, fmx = 0;

  // Neither flag given, or both given, means the soil reacts on both sides.
  bool bilateral() const { return pos == neg; }
  bool has_cap() const { return ffx > 0 || ffy > 0 || ffz > 0 || fmx > 0; }
  bool nonlinear() const { return !bilateral() || has_cap(); }
};

// Zoned bed of ONE bar (+ its optional pseudo-plastic law).
struct BarSoil {
  std::vector<SoilZone> zones;
  std::optional<SoilPlastic> plastic;
};

// Per-bar soil, keyed by bar id. Ordered: iteration order is deterministic.
using SoilZoneMap = std::map<Id, BarSoil>;

// ---------------------------------------------------------------------------
// Bed components
// ---------------------------------------------------------------------------
// Fixed component indexing shared by every array in this module:
//   0 = local x force  (kfx / ffx)   linear shapes on DOFs 0, 6
//   1 = local x moment (kmx / fmx)   linear shapes on DOFs 3, 9
//   2 = local y force  (kfy / ffy)   Hermite shapes on DOFs 1, 5, 7, 11
//   3 = local z force  (kfz / ffz)   Hermite shapes on DOFs 2, 4, 8, 10
inline constexpr int kBedComps = 4;
inline constexpr int kBedFx = 0, kBedMx = 1, kBedFy = 2, kBedFz = 3;

// A resolved piece of bed: constant moduli over [x0, x1] along the bar.
struct BedSegment {
  double x0 = 0, x1 = 0;                    // [m]
  std::array<double, kBedComps> k{};        // moduli, component-indexed
  int zone = 0;                             // originating zone index
  double length() const { return x1 - x0; }
};

// Turn the ST1 zone list into contiguous segments covering [0, L].
// Throws on an empty list, a non-monotonic / out-of-range XL, or a bad modulus.
std::vector<BedSegment> resolve_zones(const std::vector<SoilZone>& zones, double L);

// Consistent bed matrix in the ELEMENT LOCAL frame (12x12) — see the header
// comment for the shape functions and the quadrature.
Mat12 bed_matrix(const std::vector<BedSegment>& segs, double L);
Mat12 zoned_bed_matrix(const std::vector<SoilZone>& zones, double L);

// Same, in the nodal (global, rigid-arm folded) frame: T^t k_bed T. This is the
// contribution added to `ke_global`, and nothing else.
Mat12 zoned_bed_matrix_nodal(const BarElement& el, const std::vector<SoilZone>& zones);

// Whole-map validation: bars exist, are frames, carry no `Bar::soil` of their
// own (double bed), zones are well ordered, every number finite and >= 0.
void validate_soil_zones(const Model& model, const SoilZoneMap& soil);

// ---------------------------------------------------------------------------
// Pseudo-plastic driver
// ---------------------------------------------------------------------------

enum class BedState { Elastic, Capped, Inactive };

// Converged active-set state of one sub-segment.
struct BedSegState {
  double x0 = 0, x1 = 0;
  std::array<BedState, kBedComps> state{};  // default Elastic
  std::array<double, kBedComps> p_mean{};   // mean k·u over the sub-segment [N/m], [N·m/m]
  std::array<double, kBedComps> r_mean{};   // mean ACTUAL reaction intensity (after capping)
};

// Soil results of one bar, mirroring BarResult::soil_p / soil_m.
struct PseBarSoil {
  Id bar = 0;
  bool bedded = false;
  std::vector<Vec3> soil_p;    // per station, [N/m], keta sign (p = k·u)
  std::vector<double> soil_m;  // per station, [N·m/m]
  // Integral of the reaction intensity over the bar, keta sign. The force the
  // soil applies to the structure is the NEGATIVE of this.
  Vec3 total_force_local = Vec3::Zero();   // [N]
  double total_torque_local = 0;           // [N·m]
  Vec3 total_force_global = Vec3::Zero();  // [N]
  std::vector<BedSegState> segments;       // active-set states (empty when purely elastic)
};

struct PseCaseResult {
  CaseResult result;             // displacements / reactions / member forces
  std::vector<PseBarSoil> soil;  // parallel to model.bars
  int iterations = 1;
  bool plastic = false;          // true when any bed segment left the elastic state
};

struct PseResults {
  std::vector<PseCaseResult> cases;
  // Only produced when NO bar carries a nonlinear (POS/NEG or capped) spec —
  // pseudo-plastic results are not superposable. Soil fields are not combined.
  std::vector<CaseResult> combos;
};

struct PseOptions {
  int sub_per_zone = 10;      // active-set sub-segments per zone (>= 1)
  int max_iterations = 50;    // active-set sweeps before a typed failure
  double tol = 1e-8;          // relative hysteresis band on the state tests
  // --- appended ---
  // The exact beam on elastic foundation (bef.h) for every bedded bar with
  // neither releases nor a CARA VAR kernel (all of them, else the Hermite bed
  // for the whole model): exact deflection, pressures and member forces.
  bool exact = true;
};

// Zoned / pseudo-plastic statics over every load case of the model.
PseResults solve_static_pse(const Model& model, const SoilZoneMap& soil,
                            const StaticOptions& opt = {}, const PseOptions& pse = {});

}  // namespace keta

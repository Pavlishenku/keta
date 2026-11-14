// keta domain model: nodes, materials, sections, bars, supports, loads, cases,
// combinations, tendons, phases. Pure typed data + validation — no analysis here.
//
// Conventions (shared with the whole kernel):
//   * SI base units: m, N, Pa, kg, rad. Z up, right-handed global XYZ.
//   * Element local frame: x from node i to node j; y/z principal section axes.
//     The local frame is built from the "roll" angle (rotation of y/z about x),
//     with the reference orientation: local z in the global X-Z plane tilted up
//     for non-vertical bars, and local y = global Y for vertical bars.
//   * Axial force N is positive in TENSION (ST1's compression-positive N is
//     converted at the ST1 import boundary).
#pragma once
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>

#include "keta/core.h"
#include "keta/sections.h"

namespace keta {

struct BarSpanLoads;  // keta/assembly.h
struct BarResult;     // keta/statics.h


// ---------------- materials ----------------

// Time-dependent parameters for EC2 creep/shrinkage (Annex B) — optional.
struct ConcreteTime {
  double fck = 35e6;      // characteristic cylinder strength [Pa]
  double rh = 70.0;       // ambient relative humidity [%]
  double h0 = 0.3;        // notional size 2*Ac/u [m]
  double cement_s = 0.25; // cement class coefficient s: 0.20 (R), 0.25 (N), 0.38 (S)
  double t_shrink0 = 3.0; // age at start of drying shrinkage [days]
  // --- appended ---
  // false: no autogenous shrinkage (heat-cured concrete, ST1 annex B4.6:
  // eps_ca = 0, the drying one from the end of the cure)
  bool autogenous = true;
};

// Any other aging-concrete law (keta/matlaws.h `MatLaw`: BPEL, BHP, CEB, fib,
// B3 ...), as the three functions the phasing creep loop needs. Ages in days;
// moduli in the model's own units.
struct ConcreteLaw {
  std::function<double(double, double)> phi;  // creep coefficient phi(t, t0)
  std::function<double(double)> E_ref;        // modulus phi is referenced to, loaded at t0
  std::function<double(double)> shrinkage;    // shrinkage strain at age t (< 0 = shortening)
  // --- appended ---
  // false => a fibre in tension takes no creep increment (ST1 FLUAGE ... TRACTION NON)
  bool creep_tension = true;
  // set => non-linear creep (EN 1992-1-1 §3.1.4(4), ST1 LINEAIRE NON): each
  // loading event's creep is multiplied by exp(1.5 (k - 0.45)), k = max(0.45,
  // |sigma_c| / fcm(t0)); this returns fcm(t0) in the model's stress unit.
  std::function<double(double)> nonlinear_fcm;
  // the age shift [days] the laws above run on (BPEL heat cure: t_eb); the
  // concrete on a pretension bed shrinks on its real age (as in ST1)
  double age_shift = 0.0;
};

struct Material {
  Id id = 0;
  std::string name;
  double E = 210e9;    // Young's modulus [Pa]
  double nu = 0.3;     // Poisson ratio
  double rho = 7850;   // density [kg/m^3]
  double alpha = 1.2e-5;  // thermal expansion [1/K]; 0 = "no alpha" (thermal load on it is a typed error)
  double fy = 0.0;     // yield strength [Pa] (0 = not defined; needed by fiber/fire)
  double damping = 0.0;  // material damping ratio (energy-weighted modal damping)
  std::optional<ConcreteTime> concrete;  // set => EC2 creep/shrinkage laws apply
  // --- appended ---
  // Set => this law replaces the EC2 one of `concrete` in the phasing creep
  // loop (keta/phasing.h).
  std::shared_ptr<const ConcreteLaw> concrete_law;
  // Shear modulus given directly [Pa] (0 = isotropic, E / 2(1 + nu)): the
  // bars of an orthotropic material such as timber (E/G ~ 16, which no
  // isotropic nu gives). Bars only — shells and continua take (E, nu).
  double G_given = 0.0;

  double G() const { return G_given > 0 ? G_given : E / (2.0 * (1.0 + nu)); }
};

// ---------------- bars ----------------

enum class BarKind { Frame, Truss, Cable };

// End releases: 6 per end, order ux uy uz rx ry rz in the LOCAL frame.
// true = released (hinge / slider). Condensed out of the element matrices.
using Releases = std::array<bool, 12>;  // [0..5] end i, [6..11] end j

// Winkler elastic foundation moduli, bar-local frame.
// kfx/kfy/kfz: distributed translational bed [N/m per m = Pa], kmx: rotational [N·m/rad per m].
struct SoilProps {
  double kfx = 0, kfy = 0, kfz = 0, kmx = 0;
  bool any() const { return kfx != 0 || kfy != 0 || kfz != 0 || kmx != 0; }
};

// Cable-specific data (exact elastic catenary; used by cable statics).
struct CableProps {
  double pretension = 0.0;    // initial horizontal tension estimate [N]
  double unstressed_length = 0.0;  // L0 [m]; 0 = derive from chord + pretension
  // Appended: direct properties overriding the bar's material / section, and
  // the end-tension L0 route.
  double ea = 0.0;              // axial rigidity [N]; 0 = E*A of the bar
  std::optional<double> w;      // weight per unstrained length [N/m]; none = rho*A*g
  double end_tension = 0.0;     // > 0: L0 from this tension magnitude at end j [N]
  static CableProps with_length(double l0) {
    CableProps c;
    c.unstressed_length = l0;
    return c;
  }
};

struct Bar {
  Id id = 0;
  Id ni = 0, nj = 0;
  Id material = 0;
  Id section = 0;
  BarKind kind = BarKind::Frame;
  double roll = 0.0;                 // rotation of local y/z about local x [rad] (ST1 BETA)
  Releases releases{};               // all false = rigidly connected
  std::array<Vec3, 2> ecc{Vec3::Zero(), Vec3::Zero()};  // rigid-arm node-line offsets (global frame), per end
  SoilProps soil;                    // Winkler bed (zero = none)
  std::optional<CableProps> cable;   // required when kind == Cable

  bool has_ecc() const { return ecc[0].squaredNorm() > 0 || ecc[1].squaredNorm() > 0; }
};

// ---------------- nodes / supports ----------------

struct Node {
  Id id = 0;
  Vec3 x = Vec3::Zero();
};

// ---------------- tendons in slabs and shells ----------------
// A tendon carried by facet shells (keta/shell_tendon.h): its polyline at
// its true position (inside the thickness), the jacking force, the friction
// mu [1/rad] on the cumulated deviation and the wobble k [1/m], the jacked
// end(s), the carrying shells (slab id, 0: every shell).
struct ShellTendon {
  Id id = 0;
  std::vector<Vec3> pts;
  double P0 = 0;
  double mu = 0, k = 0;
  enum class Jack { Start, End, Both } jack = Jack::Start;
  Id slab = 0;
  // --- appended: the steel, for the losses of a phasing (keta/phasing.h):
  // area (0 = a force-only tendon, no loss), modulus, and the EN 1992-1-1
  // §3.3.2 relaxation (fpk = 0: none)
  double Ap = 0.0, Ep = 195e9;
  double fpk = 0.0, rho1000 = 2.5;
  int relaxation_class = 2;
};

// ---------------- composite bars (section parts) ----------------
// A part of a composite bar (keta/composite.h): its material, its section
// about its own centroid (principal axes parallel to the bar's), and the
// local z of that centroid measured from the bar's node line (the parts lie
// on the local z axis).
struct SectionPart {
  std::string name;
  Id material = 0;
  SectionProps props;
  double z = 0.0;
};
struct CompositeBar {
  std::vector<SectionPart> parts;
};

// ---------------- node-to-node links ----------------
//
// Rigid link (master-slave kinematic coupling, small displacements): for each
// slave node s and each coupled global component c,
//     translations  u_s[c] = (u_m + theta_m x (x_s - x_m))[c],
//     rotations     theta_s[c] = theta_m[c],
// eliminated EXACTLY (u = T u_r, K_r = T^t K T; no penalty, no multiplier).
// An all-true mask is a rigid body; a partial mask couples only those global
// components (e.g. {ux, uy, rz}: an in-plane floor diaphragm), the
// rigid-body formula still giving each coupled translation.
// A slave may be the master of another link (chains resolve transitively); a
// component constrained twice, a cycle, or a supported slave component is
// rejected. Supported by the linear static, modal, response-spectrum,
// buckling and phasing analyses
// (keta::assemble_elements rejects them elsewhere).
struct RigidLink {
  Id id = 0;
  Id master = 0;
  std::vector<Id> slaves;
  std::array<bool, 6> dofs{true, true, true, true, true, true};  // global ux uy uz rx ry rz
};
// General linear constraint (multi-point constraint, constraint
// equations): sum_k coef_k u(node_k, comp_k) = 0 over global DOF
// components (0..5 = ux uy uz rx ry rz). The FIRST term's DOF is the one
// eliminated (it must be free, carry no other constraint, and have a
// non-zero coefficient): u_1 = -sum_{k>1} (coef_k / coef_1) u_k, joined to
// the exact master-slave reduction of the rigid links (keta/assembly.h).
struct ConstraintTerm {
  Id node = 0;
  int comp = 0;
  double coef = 0;
};
struct Constraint {
  Id id = 0;
  std::vector<ConstraintTerm> terms;
};
// A support in any frame: `axes` rows = local x, y, z (global components,
// orthonormal); `fixed` blocks the local translations / rotations. Each
// blocked local component is the constraint a_k . u = 0 on the node's
// translations (or rotations); the global component with the largest
// |a_k| still free is eliminated. Reactions are reported in the local frame
// (CaseResult::skew_reactions).
struct SkewSupport {
  Id id = 0;
  Id node = 0;
  Mat3 axes = Mat3::Identity();
  std::array<bool, 6> fixed{};
};
// Elastic link: a 6x6 spring between two nodes, acting at the point
// P = x_i + s (x_j - x_i) reached from both nodes by rigid arms, so a rigid
// motion of the pair strains nothing. Its relative displacement there,
//     d = (u_j + theta_j x (P - x_j)) - (u_i + theta_i x (P - x_i)),
//     dtheta = theta_j - theta_i,
// expressed in the link frame, meets the stiffness `k` (6x6 symmetric,
// [N/m, N, N·m/rad] blocks): the link force on node j is -k (d, dtheta).
// Link frame: `axes` (rows = local x, y, z in global) when given; otherwise
// the bar frame of i -> j with `roll` for distinct nodes, the global frame
// for coincident ones.
// Non-linear law of one link component (force against the link-frame
// relative displacement of that component): piecewise linear through the
// points (d_i, F_i), d strictly increasing, F non-decreasing (a convex
// potential: the static solve is then a convex minimisation, keta/link_nl.h);
// beyond the end points the end segments continue. F(0) may differ from 0:
// a force present at zero deformation (an at-rest earth pressure, a
// preloaded buffer), which the structure carries from the start; the
// potential int_0^d F stays convex. Positive F = tension (the force on node j
// is -F along the component).
struct LinkLaw {
  std::vector<std::pair<double, double>> pts;
  double force(double d) const;
  double tangent(double d) const;   // slope of the segment holding d (the right one at a kink)
  double energy(double d) const;    // int_0^d F
  // tension only (slack in compression), compression only, a gap of g
  // either side before the stiffness k acts, elastic-perfectly plastic
  // (k up to |F| = Fy, then flat; kp > 0 a hardening slope)
  static LinkLaw tension_only(double k);
  static LinkLaw compression_only(double k);
  static LinkLaw gap(double g, double k);
  static LinkLaw bilinear(double k, double Fy, double kp = 0.0);
};

struct ElasticLink {
  Id id = 0;
  Id ni = 0, nj = 0;
  Mat6 k = Mat6::Zero();
  // per link-frame component, a non-linear law replacing k on it (k must
  // then have no coupling term with that component). Solved by
  // solve_static_link_laws (keta/link_nl.h); the linear analyses reject them.
  std::array<std::optional<LinkLaw>, 6> law{};
  bool has_law() const {
    for (const auto& l : law)
      if (l) return true;
    return false;
  }
  std::optional<Mat3> axes;
  double roll = 0.0;
  double s = 0.5;
  static Mat6 diagonal(double kx, double ky, double kz, double krx, double kry, double krz) {
    Mat6 m = Mat6::Zero();
    m.diagonal() << kx, ky, kz, krx, kry, krz;
    return m;
  }
};

struct PointMass {
  Id node = 0;
  double mass = 0;               // [kg]
  Vec3 inertia = Vec3::Zero();   // rotary inertia about global axes [kg·m^2]
};

struct Support {
  Id node = 0;
  std::array<bool, 6> fixed{};       // rigidly blocked DOFs (global frame)
  std::array<double, 6> spring{};    // elastic support stiffness per DOF (0 = none)
  // A modelling restraint rather than a bearing (ST1 PLANE / GRILL: the
  // out-of-plane DOFs of a node carrying no active APPUI): it blocks its
  // DOFs but reports no reaction, so an envelope over construction phases
  // never sees a "zero reaction" of a deactivated bearing.
  bool implicit = false;
  // One-sided (uplift) per DOF: 0 = bilateral, +1 = reaction must be
  // >= 0 (releases on uplift), -1 = reaction must be <= 0. Active-set solved.
  std::array<int, 6> oneside{};
  // ST1 `APPUI EL` full elastic support: the symmetric 6x6 stiffness stored
  // VERBATIM in the support frame (see `springmass.h`, manual p. 49-51). The
  // per-DOF `spring` diagonal above stays independent — both are summed into K.
  std::optional<Mat6> spring6;
  double rot = 0.0;  // ST1 `ROT alpha`: support-frame rotation about global Z [rad]
  // ST1 `DECOL` (manual p. 47, p. 49): all-or-nothing uplift. The criterion
  // bears on ONE global axis; when it is met EVERY DOF of the support is
  // released. Statics only.
  int decol_axis = -1;  // -1 = none, else 0/1/2 = global X/Y/Z
  int decol_sign = 1;   // +1 = POS (glued while the reaction is >= 0), -1 = NEG
  // --- appended ---
  // Force cap of a blocked DOF (bearing yield) [N or N.m], 0 = none: past
  // |R| = cap the support yields and carries exactly +-cap (the sign of the
  // reaction it yielded under) while the node keeps moving into it
  // (rigid-plastic support, active set); on a one-sided DOF the cap acts on
  // the admissible side only.
  std::array<double, 6> cap{};
  bool any() const {
    for (bool b : fixed)
      if (b) return true;
    for (double k : spring)
      if (k != 0) return true;
    if (spring6 && spring6->squaredNorm() > 0) return true;
    return false;
  }
};

// ---------------- modal masses (ST1 MASSE PROPRE, manual p. 69-71) ----------------

// Nodal mass matrix, global frame. `add` mirrors ST1's `SUP`: false = this
// entry REPLACES the running nodal matrix built by the earlier entries at the
// same node, true = it is added to it. It never touches the automatic bar
// masses (ST1 p. 69: `MASSE NUL` does that, and only for bars).
struct NodalMass {
  Id node = 0;
  Mat6 m = Mat6::Zero();
  bool add = false;
};

// Per-bar modal-mass override / addition.
//   Auto    — explicit no-op: keep the automatic rho*A consistent mass.
//   Nul     — ST1 `MASSE NUL`: this bar contributes no automatic mass.
//   Uni     — mu(x) = m1                                   [kg/m]
//   Lin     — mu(x) = m1 + (m2 - m1) * x/L
//   Parabol — quadratic through m1 (x=0), m2 (x=L/2), m3 (x=L)
// Uni/Lin/Parabol entries ADD a translational line mass (ST1 p. 69: "l'option
// UNI permet d'ajouter des masses lineiques"); several entries may target the
// same bar, each with its own offsets. ST1's non-`SUP` "replace" reading is a
// parser-level concern: emit a `Nul` entry alongside to drop the automatic mass.
// `ecc_or` / `ecc_ex` are the ST1 `EXC OR` / `EXC EX` offsets of the mass line
// from the element ends, GLOBAL frame [m], composed with the bar's own `ecc`.
enum class BarMassKind { Auto, Nul, Uni, Lin, Parabol };

struct BarMass {
  Id bar = 0;
  BarMassKind kind = BarMassKind::Auto;
  double m1 = 0, m2 = 0, m3 = 0;  // line mass [kg/m]
  Vec3 ecc_or = Vec3::Zero(), ecc_ex = Vec3::Zero();
  // Auto only: the form of this bar's automatic mass, overriding the
  // model-wide choice (ST1 `BARRE liste MASSE DI` / `MASSE`, set bar by bar).
  enum class Form { ModelWide, Diagonal, Consistent } form = Form::ModelWide;
};

// ---------------- flat shells (plates) ----------------

// Frame in which a `Load::SurfacePressure` vector is expressed. A dedicated
// enum (NOT `LoadFrame`, whose bar-axis semantics do not map to a surface): a
// surface pressure is either a global force-per-area vector, a scalar pressure
// along each element's own normal, or a plan-area (snow) vector.
enum class SurfaceFrame {
  Global,       // `w` is a global force/area vector [N/m^2]
  LocalNormal,  // the scalar `w.z` is the pressure along each element's +normal
  Projected,    // `w` is per m^2 of HORIZONTAL projection; scaled by |n.z|
};

// A descriptive slab record: a planar surface with an outer boundary, optional
// holes, a constant thickness, a material and a target mesh size, all SI. This
// is the USER-FACING entity; the mesher (mesh2d.h) turns it into the `ShellTri`
// units the engine actually assembles, exactly as a `Tendon` record is lowered
// into per-bar equivalent loads. Thickness is its own scalar — deliberately NOT
// a `SectionProps` set (a plate's defining property is a scalar t).
struct Slab {
  Id id = 0;
  std::vector<Vec3> outline;             // outer ring, global frame [m], not closed
  std::vector<std::vector<Vec3>> holes;  // interior rings, same convention
  double t = 0;                          // constant thickness [m], > 0
  Id material = 0;
  double mesh_size = 0;  // target mesh edge length [m], > 0
  // Optional per-slab RESULT-AXES override (a direction, magnitude ignored):
  // every facet builds its local x by projecting this vector onto the facet
  // plane instead of the default global-X rule, so the recovered mx/my/mxy read
  // along a user-chosen in-plane direction (e.g. a reinforcement grid). The
  // principals m1/m2 are frame-invariant. Rejected by validate() when
  // non-finite, zero, or (near-)parallel to the slab normal.
  std::optional<Vec3> local_x;
};

// One lowered facet-shell triangle — the per-element unit the engine assembles.
// It lives IN the model exactly like a `Bar` does. The winding order fixes the
// element normal ((p2-p1) x (p3-p1)).
//
// SENTINEL POLICY for `slab`: it is a GROUPING id, never a 0/-1 "none" magic
// value. A ShellTri whose `slab` matches an existing `Slab::id` is owned by that
// slab; one whose `slab` matches no Slab is a STANDALONE shell (legal —
// hand-built meshes carry no Slab record). validate() therefore enforces
// the shell->slab reference only in a slab-carrying model.
// Enriched facet (keta/shell.h): corner thicknesses (linear over the facet:
// haunches), an offset of the mid-surface from the nodes along the facet
// normal (a rigid offset, exact), an orthotropic material (axes at `angle`
// from the facet's local x: Ex, Ey, nu_xy, Gxy in plane, Gxz, Gyz transverse;
// the bending stiffness that of the homogeneous plate t^3/12 Q), a Winkler
// bed normal to the facet [N/m^3]. All zero / unset: the plain facet.
struct ShellOrtho {
  double Ex = 0, Ey = 0, nu_xy = 0, Gxy = 0, Gxz = 0, Gyz = 0;
  double angle = 0;  // [rad] from the facet local x to the material axis 1
};
// A laminate (appended, G11): orthotropic plies listed bottom to top, each
// with its thickness and its material axis 1 at `angle` [rad] from the
// facet local x; the bottom face at `z_bottom` above the facet mid-plane
// (`centred`: -sum t / 2). Classical laminate theory — N = A e + B k,
// M = B e + D k, A, B, D = sum Qbar_k (z), (z^2)/2, (z^3)/3 over each ply —
// plus first-order transverse shear H = sum t_k G_k (rotated), the shear
// correction 5/6. The plies' thicknesses must add up to the facet's.
struct ShellPly {
  double t = 0, E1 = 0, E2 = 0, nu12 = 0, G12 = 0, G13 = 0, G23 = 0, angle = 0;
};
struct ShellLaminate {
  std::vector<ShellPly> plies;
  bool centred = true;
  double z_bottom = 0;
};
struct ShellEnrich {
  std::array<double, 3> tn{0.0, 0.0, 0.0};
  double ecc = 0.0;
  std::optional<ShellOrtho> ortho;
  double kbed = 0.0;
  std::optional<ShellLaminate> laminate;  // (appended, G11)
  bool plain() const {
    return tn[0] == 0 && tn[1] == 0 && tn[2] == 0 && ecc == 0 && !ortho && kbed == 0 && !laminate;
  }
};

// Continuum element (appended, G15; keta/continuum.h): plane stress / plane
// strain (thickness t) and axisymmetric elements in the global XY plane
// (T3, T6, Q4, Q8), solids (T4, H8, H20); node order in continuum.h.
struct Continuum {
  enum class Kind { PlaneStress, PlaneStrain, Axisymmetric, Solid };
  Id id = 0;
  Kind kind = Kind::PlaneStress;
  std::vector<Id> n;
  Id material = 0;
  double t = 1.0;
};

struct ShellTri {
  Id id = 0;
  std::array<Id, 3> n{};  // corner node ids, winding order
  double t = 0;           // thickness [m], > 0
  Id material = 0;
  Id slab = 0;  // owning slab id (grouping key)
  ShellEnrich enrich;  // (appended) haunch, offset, orthotropy, bed
};

// ---------------- loads ----------------

enum class LoadFrame { Global, Local, GlobalProjected };  // ST1 GLO / LOC / PRO

struct Load {
  enum class Kind {
    NodalForce,     // node, F (3), M (3), global frame
    BarUdl,         // bar, q force/moment per length over [a,b] (abscissa along bar, m)
    BarLinear,      // bar, q varies linearly q_a -> q_b over [a,b]
    BarPoint,       // bar, concentrated F/M at abscissa a
    BarThermal,     // bar, uniform dt [K] + gradients gy, gz [K/m] (dT/dy, dT/dz)
    SelfWeight,     // whole model, gravity vector g [m/s^2] (default (0,0,-9.80665))
    SupportDisp,    // node, prescribed displacement per DOF (global frame); requires blocked DOF
    Prestrain,      // bar, imposed generalized strain (ST1 DEFOR): axial strain
                    // `strain` + curvatures kappa_y / kappa_z (creep, lack of fit)
    // --- flat shells (appended at the end of the enum) ---
    SurfacePressure,  // slab, uniform pressure `F` in `sframe` on every ShellTri
                      // carrying the slab id
    SurfaceThermal,   // slab, uniform dt [K] + through-thickness gradient
                      // `grad` = dT/dz [K/m] about the element normal
    // --- imposed-deformation variants (ST1 DEFOR CON / LIN, appended) ---
    PrestrainPoint,   // bar, CONCENTRATED dislocation at abscissa `a`: axial
                      // gap `strain` [m] + rotation kinks kappa_y/kappa_z [rad]
                      // (a Dirac of the Prestrain fields — a hinge-like "coude")
    PrestrainLinear,  // bar, strain/curvature fields varying LINEARLY over the
                      // full bar: strain/kappa_y/kappa_z at end i ->
                      // strain2/kappa_y2/kappa_z2 at end j
    // --- free polygonal area load (appended) ---
    PolygonPressure,  // `F` [N/m^2 of PROJECTED area, global vector] on the
                      // polygon `poly` minus its `holes`, projected along
                      // `proj` onto the facets of slab `slab` (0: every shell)
    ShellPrestress,   // shell tendon `node` (an id into Model::shell_tendons)
                      // at `strain` x its force (keta/shell_tendon.h)
    ShellImposed,     // shell `node`: an imposed strain field given as the
                      // element displacement d* (18 global DOFs in `poly`,
                      // 6 x (u, theta)) that generates it — f = Ke d*, the
                      // resultants recovered from d - d* (creep / shrinkage
                      // of the facets in the phasing)
    // --- variable hydrostatic pressure (appended, G11) ---
    Hydrostatic,      // slab (0: every shell), p = hydro_gamma x depth below
                      // the free surface through `hydro_point` along the
                      // unit `g` (zero above it), along hydro_side x the
                      // facet normal; exact consistent loads (facets cut
                      // by the free surface integrated over their wet part)
    // --- parabolic bar load (appended) ---
    BarParabolic,     // `bar`: intensities F at a, Fmid at (a + b) / 2, F2 at b
                      // [N/m]; the parabola through them (no distributed couple)
    // --- parabolic imposed strain (appended) ---
    PrestrainParabolic,  // bar, strain / curvature / twist-rate fields varying
                         // as the parabola through their values at end i
                         // (strain, kappa_*), mid-bar (strain_m, kappa_*_m)
                         // and end j (strain2, kappa_*2) — the creep curvature
                         // of a bar under a distributed load
  };
  Kind kind = Kind::NodalForce;

  Id node = 0;                // NodalForce / SupportDisp
  Id bar = 0;                 // Bar* loads
  Vec3 F = Vec3::Zero();      // force  (NodalForce, BarPoint) or start intensity (BarUdl/BarLinear)
  Vec3 M = Vec3::Zero();      // moment (NodalForce, BarPoint) or start moment intensity
  Vec3 F2 = Vec3::Zero();     // end intensity (BarLinear)
  LoadFrame frame = LoadFrame::Global;
  double a = 0.0, b = -1.0;   // abscissa range along the bar [m]; b < 0 = full length
  bool relative = false;      // a/b given as fractions of L (ST1 REL)
  double dt = 0.0, gy = 0.0, gz = 0.0;  // BarThermal
  Vec3 g = Vec3(0, 0, -kGravity);       // SelfWeight gravity vector
  std::array<std::optional<double>, 6> disp{};  // SupportDisp: prescribed value per DOF
  double strain = 0.0;                // Prestrain: imposed axial strain
  double kappa_y = 0.0, kappa_z = 0.0;  // Prestrain: imposed curvatures [1/m]
  // --- flat shells (appended fields) ---
  Id slab = 0;                                     // SurfacePressure / SurfaceThermal target
  SurfaceFrame sframe = SurfaceFrame::Global;      // SurfacePressure frame
  double grad = 0.0;                               // SurfaceThermal dT/dz [K/m]
  // --- PrestrainLinear end-j values (appended fields; PrestrainPoint reuses
  // `a` for its abscissa and strain/kappa_* as gap [m] / kinks [rad]) ---
  double strain2 = 0.0;
  double kappa_y2 = 0.0, kappa_z2 = 0.0;
  // BarUdl / BarLinear DISTRIBUTED MOMENT (appended): `M` its start intensity
  // [N·m/m], `M2` its end one (BarLinear; BarUdl reads `M` throughout). A
  // couple field, not a transverse load: it bends the bar without shear (the
  // bond of a prestressing strand at an eccentricity, keta/prestress.h).
  Vec3 M2 = Vec3::Zero();
  // PrestrainPoint TRANSVERSE dislocations (appended): jumps of the local
  // deflections v / w [m] and of the twist [rad] at `a` (ST1 DEFOR CON DY DZ
  // RX), beside the axial gap `strain` and the kinks kappa_y / kappa_z.
  double dvy = 0.0, dvz = 0.0, dthx = 0.0;
  // PolygonPressure (appended): the load polygon and its openings (3-D
  // points, projected along `proj`), and that direction.
  std::vector<Vec3> poly;
  std::vector<std::vector<Vec3>> holes;
  Vec3 proj = Vec3(0, 0, 1);
  // Hydrostatic (appended): fluid unit weight [N/m^3], a point of the free
  // surface, the side (+1: pushing along the facet normal, -1 against it)
  double hydro_gamma = 0.0, hydro_side = 1.0;
  Vec3 hydro_point = Vec3::Zero();
  // BarParabolic (appended): the intensity at the middle of the loaded range.
  Vec3 Fmid = Vec3::Zero();
  // Prestrain / PrestrainLinear (appended): imposed twist rate kappa_x =
  // d(theta_x)/dx [rad/m] — the free twist of the bar (positive about +x);
  // kappa_x2 its end-j value (PrestrainLinear).
  double kappa_x = 0.0, kappa_x2 = 0.0;
  // PrestrainParabolic (appended): the values at mid-bar; the field is
  // f(u) = f_i + (f_j - f_i) u/L + 4 (f_m - (f_i + f_j)/2) u (L - u)/L^2.
  double strain_m = 0.0, kappa_y_m = 0.0, kappa_z_m = 0.0, kappa_x_m = 0.0;

  // Convenience constructors.
  static Load nodal(Id node, const Vec3& F, const Vec3& M = Vec3::Zero());
  static Load udl(Id bar, const Vec3& q, LoadFrame fr = LoadFrame::Global, double a = 0,
                  double b = -1);
  static Load linear(Id bar, const Vec3& qa, const Vec3& qb, LoadFrame fr = LoadFrame::Global,
                     double a = 0, double b = -1);
  static Load point(Id bar, double a, const Vec3& F, const Vec3& M = Vec3::Zero(),
                    LoadFrame fr = LoadFrame::Global);
  // Parabolic distributed load through qa (at a), qmid (at (a + b) / 2) and
  // qb (at b) [N/m]; b < 0 = to the end of the bar.
  static Load parabolic(Id bar, const Vec3& qa, const Vec3& qmid, const Vec3& qb,
                        LoadFrame fr = LoadFrame::Global, double a = 0, double b = -1);
  static Load thermal(Id bar, double dt, double gy = 0, double gz = 0);
  static Load self_weight(const Vec3& g = Vec3(0, 0, -kGravity));
  static Load support_disp(Id node, const std::array<std::optional<double>, 6>& d);
  // Uniform pressure on every ShellTri carrying `slab`. In `Global` /
  // `Projected` frames `w` is a force/area vector; in `LocalNormal` the scalar
  // `w.z` is the pressure along each element's own +normal.
  static Load surface_pressure(Id slab, const Vec3& w,
                               SurfaceFrame fr = SurfaceFrame::Global);
  // Uniform temperature change `dt` [K] + through-thickness gradient
  // `grad` = dT/dz [K/m] (about the element normal) on every ShellTri carrying
  // `slab`. A material with alpha == 0 is a typed error, never a silent zero.
  static Load surface_thermal(Id slab, double dt, double grad = 0.0);
  // Free polygonal area load: `q` [N/m^2 of the area projected along
  // `proj`, global vector] on `poly` minus `holes`, carried by the facets of
  // `slab` (0: every shell) it covers — each facet takes q times the
  // projected area of its overlap, split between its nodes by the linear
  // shape functions at the overlap's centroid (exact for the facets' own
  // consistent pressure load: full coverage = surface_pressure, Projected).
  // The equivalent loads of shell tendon `tendon` times `factor` (1: the
  // jacked state with its friction losses).
  static Load shell_prestress(Id tendon, double factor = 1.0);
  static Load shell_imposed(Id shell, const Eigen::Matrix<double, 18, 1>& d);
  // Hydrostatic pressure gamma x depth below the free surface through
  // `surface_point`, depth along `down`, on the shells of `slab` (0: all),
  // pushing along side x the facet normal.
  static Load hydrostatic(double gamma, const Vec3& surface_point, Id slab = 0, double side = 1.0,
                          const Vec3& down = Vec3(0, 0, -1));
  static Load polygon_pressure(const std::vector<Vec3>& poly, const Vec3& q,
                               const std::vector<std::vector<Vec3>>& holes = {}, Id slab = 0,
                               const Vec3& proj = Vec3(0, 0, 1));
  // Uniform imposed generalized strain over the whole bar (ST1 DEFOR UNI).
  static Load prestrain(Id bar, double strain, double kappa_y = 0.0, double kappa_z = 0.0);
  // Concentrated dislocation at abscissa `a` [m] (ST1 DEFOR CON): axial gap
  // `dux` [m] + rotation kinks `theta_y`/`theta_z` [rad] about local y/z, the
  // Dirac counterparts of the uniform strain/curvature fields.
  static Load prestrain_point(Id bar, double a, double dux, double theta_y = 0.0,
                              double theta_z = 0.0);
  // Linearly varying imposed strain/curvatures over the full bar length
  // (ST1 DEFOR LIN): end-i values -> end-j values.
  static Load prestrain_linear(Id bar, double strain_i, double kappa_y_i, double kappa_z_i,
                               double strain_j, double kappa_y_j, double kappa_z_j);
  // Parabolic imposed strain / curvatures over the full bar: the values at
  // end i, at mid-bar and at end j (kappa_x through the appended fields).
  static Load prestrain_parabolic(Id bar, const Vec3& at_i, const Vec3& at_mid, const Vec3& at_j);
};

struct LoadCase {
  Id id = 0;
  std::string name;
  std::vector<Load> loads;
};

// ---------------- non-prismatic bars (ST1 `CARA VAR`) ----------------

// EXACT member kernel of a bar whose section varies along its length.
//
// The tapered kernel (`tapered.h`, ST1 manual v24 p.59-64) is layered ABOVE
// this file — it needs `Material` and `SectionProps` together — so its
// descriptor (`TaperedVar`) cannot be a member of `Bar`. A bar instead carries
// the two RAW LOCAL quantities that the prismatic `BarElement` would otherwise
// build itself, installed by the layer that owns the descriptor (the ST1
// runtime, `st1_run.cpp`):
//
//   ke — raw local 12x12 member stiffness, BEFORE release condensation, before
//        the rigid-arm eccentricity transform and before the Winkler bed; the
//        consumers fold those in exactly as they do for the prismatic block;
//   f0 — raw local work-equivalent nodal vector of ONE span load, with
//        `BarElement::f0_local`'s convention (member end forces = k u - f0).
//        Unsupported load kinds are a typed error raised by the installer.
//        EMPTY means "lower the span loads with the prismatic shape functions"
//        — what the ST1 `CARA SPECIAL` secant driver (manual p.210) wants: it
//        only rescales the member stiffness, the loads stay elastic.
//
// Three sites consult the table — the assembled stiffness (`assemble_stiffness`
// / `fold_prescribed`), the load lowering (`lower_case`) and the member-force
// recovery (`recover_case`) — through `Model::tapered_of()`. `Model::tapered`
// is EMPTY by default and every prismatic path is then bit-identical.
//
// Documented limitations (checked in `Model::validate`, never silent):
//   * a tapered bar must be a `BarKind::Frame` with NO end release (the release
//     condensation operator is built from the prismatic block);
//   * the consistent MASS and the geometric stiffness Kg still use the constant
//     `SectionProps` of the bar's `Section` — modal / buckling on a tapered bar
//     are section-averaged, not exact;
//   * station DEFLECTIONS inside the bar come from the prismatic double
//     integration of `statics.cpp`; the station internal FORCES are exact
//     (pure equilibrium from the exact end forces) and so are the end
//     displacements.
struct TaperedBar {
  Mat12 ke = Mat12::Zero();
  std::function<Vec12(const Load&)> f0;
  // Optional, for custom member blocks beyond CARA VAR (the exact beam on
  // elastic foundation, keta/bef.h): when set, the equivalent loads are
  // computed from the bar's span record of the load (`f0_span`, used instead
  // of `f0`) and the stations are recovered by `recover` (u_local = local end
  // displacements, xs = ascending abscissas) instead of the prismatic
  // equilibrium / curvature integration.
  std::function<Vec12(const BarSpanLoads&)> f0_span;
  std::function<void(const BarSpanLoads&, const Vec12& u_local, const std::vector<double>& xs,
                     BarResult&)>
      recover;
  // Optional section law along the bar, t = x / L in [0, 1] (CARA VAR): the
  // prismatic station recovery then integrates the curvature M / EI(x), the
  // shear angle V / GA(x) and the strain N / EA(x) instead of the constant
  // properties of the bar's Section.
  std::function<SectionProps(double t)> props;
};

// Linear combination of case results: sum(factor * case).
struct Combination {
  Id id = 0;
  std::string name;
  std::vector<std::pair<Id, double>> terms;  // (case id, factor)
};

// ---------------- prestress tendons ----------------

// Tendon eccentricity profile along a chain of bars: piecewise polynomial e_y(s), e_z(s)
// (local frame, s = curvilinear abscissa from the anchor at the start).
struct PolySegment {
  double s0 = 0, s1 = 0;              // segment range [m]
  std::array<double, 4> cy{};         // e_y(s) = cy0 + cy1*ds + cy2*ds^2 + cy3*ds^3, ds = s - s0
  std::array<double, 4> cz{};
  // --- appended ---
  // The EXACT curve of the segment when it is not a cubic in s (cy / cz then
  // hold its Hermite cubic in s, same ends and slopes, for readers of the
  // coefficients): e(s) is read on the curve by inverting s(u).
  //  1 = parametric cubic (s, e_y, e_z)(u) = sum p[k] u^k, u in [0, 1]
  //      (ST1's trace between two points);
  //  2 = circular arc X(phi) = C - R cos(phi) U + R sin(phi) D, phi in
  //      [0, th] (a RAYON fillet of a trace by poles).
  int exact = 0;
  std::array<Vec3, 4> p{Vec3::Zero(), Vec3::Zero(), Vec3::Zero(), Vec3::Zero()};  // kind 1: u^k
  Vec3 C = Vec3::Zero(), U = Vec3::Zero(), D = Vec3::Zero();  // kind 2
  double R = 0, th = 0;
};

struct TendonLosses {
  double mu = 0.19;      // friction coefficient [1/rad]
  double k = 0.01;       // wobble [rad/m]
  double slip = 0.006;   // anchorage draw-in [m]
  double Ep = 195e9;     // tendon modulus [Pa]
  double Ap = 0.0;       // tendon area [m^2] (0 = force-only tendon, no elastic-shortening coupling)
  double rho1000 = 2.5;  // relaxation class loss at 1000h [%]
  // Deferred-loss coupling in the phasing creep loop (EN 1992-1-1 eq. 5.46
  // incremental — see keta/phasing.h). fpk = 0 (the default) disables the
  // RELAXATION term only; the creep + shrinkage terms need just Ap > 0.
  double fpk = 0.0;          // tendon characteristic strength [Pa] (0 = no relaxation)
  int relaxation_class = 2;  // EC2 §3.3.2 relaxation class 1 / 2 / 3
  // --- appended ---
  // Relaxation by the BPEL 91 annex 2 law (fpk = fprg) under varying strain by
  // the equivalent-time method — what ST1 runs in a phasing whose concrete is
  // BPEL (manual annex B4.3.2). false = the EC2 §3.3.2 curve.
  bool relaxation_bpel = false;
  // --- appended ---
  // Heat-cured pretension (EN 1992-1-1 §10.3.2 / §10.5.2, ST1 annex B4.6):
  // the equivalent time t_eb added to the relaxation clock [days], and the
  // thermal loss of the strand applied when the bed is released [Pa].
  double relax_teb_days = 0.0;
  double thermal_loss = 0.0;
  // --- appended ---
  // Line-friction exponent per unit length independent of mu [1/m]: ST1's
  // BPEL `PERTE F 0 PHI phi` (sigma = sigma_a exp(-(f alpha + phi s)), f = 0).
  double phi_abs = 0.0;
  // --- appended ---
  // Elementary units tensioned ONE AFTER THE OTHER (ST1 `PREC j NON SIMUL n`,
  // annex B1.1.3): in a phasing the tendon takes (n - 1) / (2 n) of the
  // elastic loss its own jacking causes. 1 = simultaneous (no such loss).
  int non_simul = 1;
};

struct Tendon {
  Id id = 0;
  std::string name;
  std::vector<Id> bars;             // chain of bars the tendon runs through, in order
  std::vector<PolySegment> profile; // e(s) profile over the chain
  double jack_force = 0.0;          // P0 at the active anchor [N]
  bool jack_both_ends = false;
  TendonLosses losses;
  // --- appended (parametric layer tendons); defaults keep ST1 ---
  // AUTHORED force profile: (s, P) breakpoints over the chain, piecewise
  // linear, constant beyond the ends (one point = a constant force). When
  // non-empty it REPLACES the jack-force friction / draw-in law (the losses
  // are part of what was authored); `jack_force` then carries its chain mean,
  // the reference the phasing fractions and loss bookkeeping scale.
  std::vector<std::pair<double, double>> force_profile;
  // flips[i] = the chain runs through bars[i] from its nj to its ni (empty =
  // every bar along the chain). The eccentricity stays in each bar's OWN
  // section frame; only the arc parameter is reversed.
  std::vector<bool> flips;
  // Load case carrying this tendon's equivalent loads in a static analysis
  // (0 = none). Phased analysis tensions tendons through TensionActions, so a
  // phase applying this case would double-count them. Informational: not
  // validated against `cases` (the phasing sub-models carry a subset of them).
  Id case_id = 0;
  // Anchor abscissae along the chain [m]: a tendon may start and stop INSIDE
  // its bars (ST1 CABLE: the trace extent, e.g. a cable anchored 4 m into its
  // first bar). The force is zero outside [s_begin, s_end]; friction and
  // draw-in run from these anchors. s_end < 0 = the chain end.
  double s_begin = 0.0, s_end = -1.0;
  // The bar takes the force's AXIAL component P cos(alpha) (and P cos(alpha)
  // e' across), alpha the cable's angle to the bar — ST1's external method;
  // false keeps the small-slope theory (N = P, P e' across).
  bool exact_angle = false;
  // --- appended ---
  // One active anchor at the chain END (ST1 `TENSION EX`) instead of its origin.
  bool jack_from_end = false;
  // Double-end jacking released at the end first (ST1 `TENSION EX OR`); false
  // = origin first (`OR EX`). Matters only when the two draw-ins overlap.
  bool release_end_first = false;
  // Chain intervals [s0, s1] where the tendon runs FREE of the concrete (the
  // straight runs of an external tendon, ST1 TRACE EXTERIEUR): the wobble
  // term of the friction law does not act there, only the angular one.
  std::vector<std::pair<double, double>> free_runs;
};

// ---------------- staged construction ----------------

struct TensionAction {
  Id tendon = 0;
  // TARGET cumulated fraction of the jack force after this phase. Above the
  // already-applied fraction = tensioning by the difference; below it =
  // DETENSION down to the target (a negative value = full release, the ST1
  // `DETENDRE` marker) — see keta/phasing.h for the release semantics.
  double fraction = 1.0;
  // --- appended ---
  // Pretension on a bed (ST1 TENDRE_BANC / RELACHER_BANC, p.193-194): `Bed`
  // tensions the strand against the bed — no force on the structure, the
  // strand relaxes in the bed; `Release` transfers it by bond (the tendon's
  // force profile = its transmission law) after the thermal loss, the strand
  // then bonded: its own elastic shortening is a loss. `fraction` is unused.
  enum class Mode { Jack, Bed, Release };
  Mode mode = Mode::Jack;
  // ST1 `TENDRE CABLE ... NON INJECTE` (p.191): the tendon slides in its
  // duct until `Phase::inject` — the strain changes along it reach it as one
  // force change over its whole length (no local loss), not point by point.
  bool unbonded = false;
  // Command group within the phase: actions of one group are simultaneous
  // (K = 0 among them); the groups follow one another, a later group's
  // increment giving the K = 1 elastic loss to the tendons of the earlier
  // ones (ST1 p.181: without ET, each TENDRE / DETENDRE line is a step).
  int group = 0;
  // `tendon` names a SHELL tendon (Model::shell_tendons): jacked to
  // `fraction` of its P0 profile (Mode::Jack only, bonded)
  bool shell = false;
};

struct Phase {
  std::string name;
  double time = 0.0;                // phase date [days] (drives creep/shrinkage between phases)
  std::vector<Id> activate_bars;    // cast stress-free in this phase
  std::vector<Id> deactivate_bars;  // removed, accumulated state springs back
  std::vector<Support> add_supports;
  std::vector<Id> remove_support_nodes;
  std::vector<Id> cases;            // load cases applied (and kept) from this phase on
  std::vector<TensionAction> tension;
  // --- appended ---
  // Concrete age [days] of a bar AT its activation (ST1 `ACTIVER BARRE AGE a`,
  // p.187: the initial age is independent of the activation date). A bar not
  // listed is cast at activation (age 0).
  std::map<Id, double> activate_age;
  // Bars activated IN THE PROLONGATION of the deformed structure (ST1
  // `ACTIVER BARRE ... MODE REL`, its default, p.187): when exactly one end
  // node already belongs to the structure (an active bar or support), the new
  // node starts with the rigid continuation of that end (translation +
  // rotation) instead of the undeformed geometry. Displacements only — the
  // bar is cast stress-free either way.
  std::set<Id> activate_rel;
  // ST1 `MODIFIER APPUI` (p.195): the support at `node` takes the restraints
  // of the given Support, after the deactivations. A component no longer
  // blocked releases its accumulated reaction (spring-back); a newly blocked
  // one is set on the deformed structure (no reaction yet).
  std::vector<Support> modify_supports;
  // ST1 `VERINER APPUI` (p.186): the TOTAL displacement of every blocked
  // component of the support at these nodes is brought back to zero (an
  // imposed displacement increment equal to minus the accumulated one).
  std::vector<Id> jack_support_nodes;
  // ST1 `INJECTER CABLE` (p.191): these unbonded tendons are grouted — from
  // this phase on their losses are local again.
  std::vector<Id> inject;
  // Node-to-node links (model.h ElasticLink / RigidLink): activated
  // stress-free in this phase (they constrain the LATER increments only); an
  // elastic link is deactivated by spring-back of its accumulated force. A
  // link never named in an activate list of the plan is active from the
  // start. Rigid links cannot be deactivated.
  std::vector<Id> activate_links, deactivate_links;
  std::vector<Id> activate_rigid_links;
  // Parts of composite bars (bar id, part index) cast stress-free in this
  // phase (keta/composite.h). A part no phase of the plan casts is active
  // with its bar.
  std::vector<std::pair<Id, int>> activate_parts;
  // --- appended ---
  // Section changes of ACTIVE bars (bar, new section), after the
  // deactivations and before the support modifications: the later increments see the new stiffness, the
  // accumulated forces stay (a swap moves no force). Refused for a creeping
  // material (the locked-in stress of the old section would creep as the new
  // one — cast the added concrete as a composite part, `activate_parts`), a
  // CARA VAR (tapered) bar or a composite bar.
  std::vector<std::pair<Id, Id>> modify_sections;
  // Phase hinges at one end of an active bar (end 0 = i, 1 = j; dofs = the
  // LOCAL components 0..5, e.g. dofs[5] = rz). `place_hinges` releases them
  // for the later increments — only where the accumulated end force on those
  // components is zero (releasing a loaded component would move the
  // equilibrium). `close_hinges` (clavage) makes them continuous again: the
  // relative rotation accumulated while open stays locked in, at zero force;
  // `jump` (axial gap, twist, rotation about local y, about local z) then
  // forces an additional relative jump into the closed joint (a concentrated
  // dislocation at that end: a free offset on a determinate structure,
  // self-stress on an indeterminate one).
  struct HingeAction {
    Id bar = 0;
    int end = 1;
    std::array<bool, 6> dofs{};
    std::array<double, 4> jump{};
  };
  std::vector<HingeAction> place_hinges, close_hinges;
  // ST1 `INJECTER CABLE ... PARTIELLEMENT`: an unbonded tensioned
  // tendon grouted on the chain interval [s0, s1] only — bonded there (its
  // losses local again), still sliding with Coulomb friction elsewhere: the
  // grouted zone anchors the free stretches next to it (the force passes to
  // the bond over a short transfer length), so each free stretch slides
  // between its own anchors, its force set by the compatibility with its duct.
  // Needs PhasingOptions::local_losses. Cumulative with later intervals.
  struct PartialGrout {
    Id tendon = 0;
    double s0 = 0, s1 = 0;
  };
  std::vector<PartialGrout> inject_partial;
  // Catenary cables (BarKind::Cable bars). A cable enters the
  // phasing when an `activate_bars` list names it (a cable never named stays
  // OUT of the phasing); from then on EVERY increment is the coupled cable +
  // frame one: K du + sum_c [q_c(x + du; L0') - q_c(x; L0)] = dF, q_c the
  // exact elastic catenary (keta/cable.h) at the true anchor positions x (the
  // accumulated displacements), L0 its unstrained length — the only state of
  // a cable. Its installation is one coupled increment (after the hinges),
  // L0 from its CableProps (unstressed_length, else end_tension, else
  // pretension) resolved at the anchors AS THEY STAND then; the anchors must
  // belong to the active structure. Not supported with nonlinear soil
  // springs, unilateral supports or rigid links (typed refusals).
  // `cable_actions`, after the installations, one coupled increment each:
  //   Length      L0 := value
  //   Horizontal  L0 such that H = value at the EQUILIBRIUM of the increment
  //               (the anchors move with the re-tension: a secant on L0)
  //   EndTension  the same for the end-j tension magnitude (taut branch)
  //   Thermal     L0 := L0 (1 + value), value = alpha dT
  //   Rupture     the cable leaves: its force released statically; value
  //               = the dynamic factor phi >= 1 of the accidental state
  //               (PhaseResult::cable_rupture_peak = (phi - 1) x the release)
  struct CableAction {
    enum class Kind { Length, Horizontal, EndTension, Thermal, Rupture };
    Id bar = 0;
    Kind kind = Kind::Length;
    double value = 0;
  };
  std::vector<CableAction> cable_actions;
  // Analysis regime of the increments of this phase: Linear, or
  // SecondOrder — every increment (loads, creep, losses, cables excepted)
  // solves the equilibrium of the TOTAL state with the geometric stiffness of
  // the TOTAL axial force,
  //   [K + Kg(N_new)] du = dF - sum_e [kg_e(N_new) - kg_e(N_old)] d_e,
  // N_new = N_acc + N(du) by fixed point, d_e the displacement of bar e since
  // its stress-free activation plus the injected imperfection — so a later
  // axial load amplifies an earlier deflection (path-consistent: phases
  // that are all second order end on the one-shot solve_second_order).
  // Refused beyond the critical load (LDLT inertia), with prescribed
  // displacements, live catenary cables, nonlinear soil, unilateral supports,
  // rigid links and facet shells.
  enum class Regime { Linear, SecondOrder };
  Regime regime = Regime::Linear;
  // Initial geometric imperfection injected at the start of the phase (node,
  // translations + rotations), stress-free, cumulative: it acts only through
  // Kg in the second-order increments (a linear increment ignores geometry).
  std::vector<std::pair<Id, Vec6>> imperfection;
};

// The time-dependent effects between two dates are decided by the materials
// (a Material with `concrete` or `concrete_law` creeps and shrinks) and the
// live tendons (relaxation data), never by the plan — keta/phasing.h.
struct PhasePlan {
  std::vector<Phase> phases;
};

// ---------------- model ----------------

struct Model {
  std::vector<Node> nodes;
  std::vector<Material> materials;
  std::vector<Section> sections;
  std::vector<Bar> bars;
  std::vector<Support> supports;
  std::vector<PointMass> point_masses;
  std::vector<NodalMass> nodal_masses;  // ST1 MASSE PROPRE NOEUD
  std::vector<BarMass> bar_masses;      // ST1 MASSE PROPRE BARRE
  // true => only the bars carrying a BarMassKind::Auto entry get their
  // automatic rho*A element mass (ST1 MASSE PROPRE BARRE listeb, p.69:
  // "les barres pour lesquelles il souhaite affecter des masses")
  bool bar_mass_selective = false;
  std::vector<LoadCase> cases;
  std::vector<Combination> combos;
  std::vector<Tendon> tendons;
  // Flat shells: the descriptive `Slab` records and the lowered `ShellTri`
  // facets the engine assembles (empty in a bar-only model — every shell path
  // is then bit-identical to the pre-shell kernel).
  std::vector<Slab> slabs;
  std::vector<ShellTri> shells;
  std::vector<Continuum> continua;  // (appended, G15)
  // ST1 CARA VAR: exact member kernels of the non-prismatic bars, keyed by bar
  // id (see `TaperedBar`). Empty = every bar prismatic.
  std::map<Id, TaperedBar> tapered;
  // Node-to-node links (see RigidLink / ElasticLink above).
  std::vector<RigidLink> rigid_links;
  std::vector<ElasticLink> links;
  std::vector<ShellTendon> shell_tendons;
  // Composite bars by bar id (keta/composite.h): handled by the phasing;
  // the other analyses take `composite_lowered(model)`.
  std::map<Id, CompositeBar> composite;
  // General linear constraints and skew supports (eliminated with the rigid
  // links' master-slave reduction).
  std::vector<Constraint> constraints;
  std::vector<SkewSupport> skew_supports;
  bool has_links() const {
    return !rigid_links.empty() || !links.empty() || !constraints.empty() || !skew_supports.empty();
  }

  // --- builder helpers (return the created id) ---
  Id add_node(double x, double y, double z);
  Id add_material(const Material& m);
  Id add_section(const Section& s);
  Id add_bar(Id ni, Id nj, Id mat, Id sec, BarKind kind = BarKind::Frame);
  // Energy-consistent shear areas of the bar's section for its material's
  // Poisson ratio (Section::shear_nu); called by add_bar.
  void resolve_shear(Bar& b);
  void add_support(Id node, const std::array<bool, 6>& fixed,
                   const std::array<double, 6>& spring = {});
  Id add_case(const std::string& name, std::vector<Load> loads);
  // a continuum element (G15): validated node count for its kind
  Id add_continuum(Continuum::Kind kind, const std::vector<Id>& nodes, Id material, double t = 1.0);
  Id add_combo(const std::string& name, std::vector<std::pair<Id, double>> terms);
  Id add_slab(const Slab& s);
  Id add_shell(const std::array<Id, 3>& n, double t, Id mat, Id slab);
  Id add_rigid_link(Id master, std::vector<Id> slaves,
                    const std::array<bool, 6>& dofs = {true, true, true, true, true, true});
  Id add_elastic_link(Id ni, Id nj, const Mat6& k, double s = 0.5);
  Id add_constraint(std::vector<ConstraintTerm> terms);
  Id add_skew_support(Id node, const Mat3& axes, const std::array<bool, 6>& fixed);

  // --- lookups (validated maps built lazily by index()) ---
  const Node& node(Id id) const;
  const Material& material(Id id) const;
  const Section& section(Id id) const;
  const Bar& bar(Id id) const;
  const LoadCase& load_case(Id id) const;
  const Combination& combo(Id id) const;
  const Tendon& tendon(Id id) const;
  const Slab& slab(Id id) const;
  int node_index(Id id) const;  // dense index into nodes[]
  int bar_index(Id id) const;
  // True when `id` names a loadable surface: a `Slab` record OR a `ShellTri`
  // group carrying it (a standalone hand-built mesh has no Slab record).
  bool has_surface(Id id) const;
  // Exact tapered kernel of a bar, or nullptr when the bar is prismatic.
  const TaperedBar* tapered_of(Id bar) const {
    auto it = tapered.find(bar);
    return it == tapered.end() ? nullptr : &it->second;
  }

  double bar_chord_length(const Bar& b) const;  // node-to-node distance (before ecc)

  // Whole-model validation: finiteness sweep, dangling references, degenerate bars,
  // thermal-on-alpha-less material, cable/truss constraints... Throws Error.
  void validate() const;

  // Rebuild the id->index maps (called by validate(); analyses call it once up front).
  void index() const;
  // Direct mutation of the entity vectors (bypassing add_*) leaves the lazy
  // maps stale — call this afterwards so the next lookup rebuilds them.
  void invalidate_index() const { indexed_ = false; }

 private:
  mutable std::unordered_map<Id, int> node_ix_, mat_ix_, sec_ix_, bar_ix_, case_ix_, combo_ix_,
      tendon_ix_, slab_ix_, shell_ix_;
  mutable bool indexed_ = false;
};

}  // namespace keta

// Flat facet-triangle SHELL element — keta's plates-and-shells layer. Since
// "coques acte II" the element carries a TRUE drilling-DOF membrane (Allman),
// a consistent mass matrix and a membrane geometric stiffness, so shells
// participate in statics, modal, buckling, second-order, dynamics, moving
// loads, cables, fire, zoned soils and phasing.
// The corotational large-displacement, arc-length and ruin drivers carry them
// as EICR facets — no driver rejects a shell-bearing model any more (see the
// gate notes at the bottom of this header).
//
// ===========================================================================
// PHYSICS
// ===========================================================================
//
// A flat, constant-thickness triangular facet carrying membrane + bending +
// transverse-shear action with full 3-D orientation (plates, walls, ramps,
// folded plates, and FACETED CURVED surfaces — see the membrane below). Six
// engineering DOFs per node (ux, uy, uz, rx, ry, rz) in the GLOBAL frame => an
// 18x18 symmetric stiffness.
//
// --- Local frame ---
// From the node winding n1 -> n2 -> n3 (right-handed):
//   * local z = unit element normal (p2-p1) x (p3-p1) / ||.||;
//   * local x = global X projected onto the element plane and normalised; when
//     that projection collapses (the normal is within `kShellNearParallelTol`
//     of +-global X) global Y is projected instead;
//   * local y = z x x.
// An OPTIONAL per-slab `local_x` override replaces the global-X seed with a
// user-chosen in-plane direction (projected onto the facet plane), so the
// recovered mx/my/mxy resolve along an authored axis (e.g. a reinforcement
// grid). The principals m1/m2 are frame-invariant; only the axes move.
// The 3x3 rotation R (rows = local x,y,z in global coords) maps a global vector
// to local (v_local = R*v_global); the 18x18 transform T is block-diagonal (R
// repeated over the six translation/rotation triples) and K_global = T'*K_l*T.
//
// --- Membrane: ALLMAN triangle (true drilling DOFs) ---
// Allman's triangle with vertex drilling rotations (Allman, Comput. Struct. 19
// (1984) 1-8; the 1988 quadratic-edge form), implemented as the CONSTRAINED
// LST identified by Cook (Comput. Struct. 22 (1986) 1065-1067): the quadratic
// 6-node triangle whose midside displacements are slaved to the corner DOFs by
//   u_m = (u_i + u_j)/2 + (y_j - y_i)(theta_j - theta_i)/8
//   v_m = (v_i + v_j)/2 - (x_j - x_i)(theta_j - theta_i)/8
// on each edge i->j. The membrane block is therefore 9x9 on [u, v, rz] per
// node: B_Allman(r,s) = B_LST(r,s) * T_A (linear in (r,s)), integrated with the
// 3-point mid-edge rule (exact for the quadratic B'DB). Properties:
//   * the CONSTANT-STRAIN PATCH TEST passes exactly (equal corner rotations
//     make every hierarchical midside term vanish => pure linear field);
//   * the drilling stiffness is PHYSICAL (edge-normal quadratic displacement),
//     which removes the assembly-level membrane LOCKING that a CST facet shows
//     on faceted curved surfaces — Scordelis-Lo / hemisphere / pinched
//     cylinder converge to their reference values, so curved shells are
//     within the domain;
//   * ONE spurious zero-energy mode remains: the equal-rotation mode
//     (theta_1 = theta_2 = theta_3, translations zero). It is stabilised by
//     the classical penalty on the difference between the MEAN drilling
//     rotation and the CST (linear-field) continuum spin (Cook, Malkus,
//     Plesha & Witt, "Concepts and Applications of FEA", 4th ed., sec. 8.2),
//     with the Hughes-Brezzi modulus gamma = G (kShellDrillStab = 1):
//         E_stab = 1/2 * kShellDrillStab * G t A * (theta_bar - omega_CST)^2,
//     omega_CST = (dv/dx - du/dy)/2 of the corner (linear) field. The penalty
//     vanishes identically on rigid rotations and on every constant-strain
//     state whose nodal rotations equal the continuum spin, so the patch test
//     stays exact and NO recovered resultant ever reads from it. This is the
//     ONLY fictitious stiffness left in the element (the old alpha = 1e-5
//     diagonal drilling spring is gone — drilling stiffness is now real), and
//     it is never a design-moment path.
//
// --- Bending + shear: Mindlin-Reissner, MITC3 ---
// Mindlin-Reissner plate with MITC3-style ASSUMED transverse shear (Lee & Bathe
// 2004): the covariant transverse shear is held constant along each edge and
// tied at the edge midpoints, so the element does not shear-lock in the thin
// limit. The bending curvature Bb is constant (linear rotations) => bending is
// integrated 1-point (exact); the assumed shear field is (r,s)-linear => shear
// is integrated with the 3-point mid-edge rule (exact for the quadratic Bs'Bs).
// Shear correction kappa = 5/6.
//
// Rotational-DOF convention (documented, self-consistent — the eternal plate
// sign trap): rx, ry are the components of the mid-surface-normal rotation
// VECTOR about local x, y (right-hand rule). The through-thickness
// displacement of a fibre at height z is theta x (z*zhat) = (z*ry, -z*rx, 0), so
//   kx = d(ry)/dx,  ky = -d(rx)/dy,  kxy = d(ry)/dy - d(rx)/dx;
//   gxz = dw/dx + ry,  gyz = dw/dy - rx.
// Moments per unit width m = Db*kappa with Db = E t^3/(12(1-nu^2))*[...];
// transverse shears per unit width q = kappa_s*G*t*gamma.
//
// --- Consistent mass ---
// 18x18, built in the local frame and congruent-transformed:
//   * translational: rho*t*A with the exact linear-triangle consistent matrix
//     (A/12) [[2,1,1],[1,2,1],[1,1,2]] per direction (the 1-point-per-monomial
//     closed form — exact for the product of linear shape functions);
//   * rotary inertia rho*t^3/12 on the BENDING rotations rx, ry (same
//     consistent shape), the classical Mindlin plate rotary term;
//   * NO mass on the drilling rz. Harmless by construction: rz now carries
//     REAL (Allman) stiffness, so a massless rz is a static DOF condensed by
//     the eigen pencil (mu -> 0 pairs are skipped by `modal_from_pairs`), never
//     a spurious rigid mode.
// The lumped variant mirrors the bar convention: tributary A/3 per node on the
// translations, A/3 * rho t^3/12 on rx, ry.
//
// --- Membrane geometric stiffness ---
// The classical plate initial-stress (von Karman) term
//     Kg = integral [dw/dx dw/dy] [[nx, nxy],[nxy, ny]] [dw/dx dw/dy]' dA
// with n the membrane resultants [N/m] of the prestress state (tension
// stiffens, compression softens — same sign convention as the bar Kg(N)).
// w is linear on the facet, so grad w is CONSTANT and the 1-point rule is
// EXACT; the block lives on the three local w DOFs. This is the term that
// yields the k = 4 uniaxial SS plate buckling coefficient and the
// 1/(1 - P/Pcr) second-order amplification. In-plane and rotational initial-stress couplings
// are not modelled (documented limitation of the flat-facet formulation).
//
// --- Resultant recovery ---
// At the element CENTROID, stress-based convention: mx is the
// moment PRODUCING sigma_x stresses, positive = tension on the +z face.
// Membrane nx, ny, nxy [N/m] (Allman strain at the centroid), bending mx, my,
// mxy [N.m/m], transverse shear qx, qy [N/m], plus principals m1 >= m2 /
// n1 >= n2 and their angles.
//
// --- Consistent loads ---
// A uniform pressure (global force/area vector, local-normal scalar, or
// plan-area "projected" vector) and self-weight lower to tributary A/3 nodal
// translational forces — the consistent vector for a uniform field on a linear
// triangle; rotational terms are dropped (documented standard practice).
// `Projected` scales by |n.z| (the cosine between the element normal and
// vertical): a flat panel gets the full q, a panel tilted by theta gets
// q*cos(theta) (its plan area is A*cos theta), a VERTICAL panel gets exactly 0.
//
// --- Thermal: initial-strain method ---
// A uniform temperature change and a through-thickness gradient enter as
// INITIAL STRAINS eps0 (membrane) and initial curvatures kappa0 (bending),
// already multiplied by the material's alpha (the element is alpha-oblivious —
// the caller owns the coefficient, mirroring the bar thermal path). The
// equivalent nodal load is the standard initial-strain vector
// f = int B' D eps0 dV, split membrane + bending:
//   f_m = (t*A) * Bm' * (Dm*eps0)   on the local [u, v, rz] DOFs (Bm at the
//                                   centroid — exact, since B_Allman is linear
//                                   in (r,s) and the mean of a linear field is
//                                   its centroid value);
//   f_b = A     * Bb' * (Db*kappa0) on the local [w,rx,ry] DOFs.
// Curvature sign: the in-plane strain of a fibre at height z is eps(z) = z*kappa
// in this convention, so a gradient grad = dT/dz gives eps0(z) = alpha*grad*z
// => kappa0 = alpha*grad*[1,1,0] with NO sign flip.
// RECOVERY MUST SUBTRACT the initial strains (`recover` takes eps0/kappa0):
//   n = t*Dm*(eps - eps0),  m = Db*(kappa - kappa0),
// otherwise a FREE thermal expansion reports a spurious resultant (the classic
// "deformations right, stresses wrong" bug): a free expansion is
// stress-free, a restrained one carries the full thermal resultant.
//
// ===========================================================================
// UNITS / LIMITS
// ===========================================================================
// SI throughout: coordinates m, t m, E/nu Pa/-, rho kg/m^3, gravity m/s^2;
// outputs N, N/m, N.m/m, rad.
//
// FORMULATION DOMAIN (not "deferred features" — the honest boundary of a flat
// facet, small-displacement, small-strain element, the way an Euler beam does
// not do shear):
//   * geometry is FLAT per facet; curved surfaces are represented by faceting,
//     validated on the MacNeal-Harder curved obstacle course (Scordelis-Lo,
//     pinched hemisphere, pinched cylinder);
//   * kinematics are LINEAR in the element; LARGE displacements are carried by
//     the corotational drivers (solve_corotational, arc_length_trace,
//     ruin_load), which wrap this element in an element-independent
//     corotational frame (EICR, nonlinear.cpp: CorotShell);
//   * the geometric stiffness carries the membrane-force/out-of-plane-gradient
//     term only (linearised buckling and P-Delta about a membrane state).
//
// References:
//   * Allman, "A compatible triangular element including vertex rotations for
//     plane elasticity analysis", Comput. Struct. 19 (1984) 1-8.
//   * Cook, "On the Allman triangle and a related quadrilateral element",
//     Comput. Struct. 22 (1986) 1065-1067 (the constrained-LST identity).
//   * Cook, Malkus, Plesha & Witt, "Concepts and Applications of Finite
//     Element Analysis", 4th ed., sec. 8.2 (drilling DOFs, spurious-mode
//     penalty).
//   * Lee & Bathe, "Development of MITC isotropic triangular shell finite
//     elements", Comput. Struct. 82 (2004) 945-962 (assumed transverse shear).
//   * Timoshenko & Woinowsky-Krieger, "Theory of Plates and Shells"
//     (the plate-bending closed forms).
//   * MacNeal & Harder, "A proposed standard set of problems to test finite
//     element accuracy", Finite Elem. Anal. Des. 1 (1985) 3-20 (patch tests,
//     the curved-shell problems and their reference values).
#pragma once
#include <array>
#include <optional>

#include "keta/model.h"

namespace keta {

using Vec18 = Eigen::Matrix<double, 18, 1>;
using Mat18 = Eigen::Matrix<double, 18, 18>;

// Penalty coefficient of the Allman spurious-mode stabilisation:
// K_stab = kShellDrillStab * G t A * g'g with g d = theta_bar - omega_CST.
// It vanishes on rigid rotations and constant-strain states and no recovered
// resultant reads from it. The value is the Hughes-Brezzi gamma = G (Hughes &
// Brezzi, CMAME 72 (1989) 105-121; Ibrahimbegovic, Taylor & Wilson, IJNME 30
// (1990) 445-457: the drilling-rotation/skew-strain constraint with modulus
// G). It was 1e-3: the equal-rotation mode of the Allman membrane, nearly
// free, coupled to the bending of neighbouring NON-coplanar facets and the
// MacNeal-Harder twisted beam diverged under refinement (5.43 -> 6.85 mm from
// 12x2 to 192x32, reference 5.424); with gamma = G it converges monotonically
// (5.4181 mm at 384x64). Flat stiffened plates converge more slowly from
// above.
inline constexpr double kShellDrillStab = 1.0;
// Shear correction factor for the Mindlin transverse-shear stiffness.
inline constexpr double kShellShearCorrection = 5.0 / 6.0;
// The in-plane projection of the frame seed axis must be longer than this (on
// unit vectors) or the fallback seed is used.
inline constexpr double kShellNearParallelTol = 1e-6;
// Degeneracy: reject the triangle when 2A <= kShellAreaDegenTol * (longest
// edge)^2. Relative (area scales as length^2) so it tracks the mesh scale.
inline constexpr double kShellAreaDegenTol = 1e-12;

// Orthonormal frame of a PLANE from its (not necessarily unit) normal: rows =
// local (x, y, z) expressed in global coordinates, local z = the normalised
// normal. Local x seeds from `local_x` (projected onto the plane) when given
// and usable, else from global X, else from global Y — THE shared seed rule.
// This is the single geometric definition used by the facet element AND the
// slab mesher (slabmesh.h), so a meshed slab and its shells agree on axes by
// construction. Typed error on a non-finite / zero normal.
Mat3 plane_frame(const Vec3& normal, const std::optional<Vec3>& local_x = std::nullopt);

// Local frame of a flat facet triangle (winding p0 -> p1 -> p2): rows = local
// (x, y, z) expressed in global coordinates, so v_local = R * v_global.
// Delegates to `plane_frame` with the winding normal (p1-p0) x (p2-p0).
// `local_x`, when set, seeds local x with that direction projected onto the
// facet plane; a non-finite / zero / normal-parallel override falls back to the
// standard rule (defensive only — Model::validate rejects such an authoring).
// Typed error on non-finite coordinates or a degenerate (collinear) triangle.
Mat3 facet_frame(const std::array<Vec3, 3>& p,
                 const std::optional<Vec3>& local_x = std::nullopt);

// Principal values of the symmetric 2x2 tensor [[a,c],[c,b]]: p1 >= p2, and
// `theta` (rad) the angle from local x to the p1 axis. ONE canonical closed
// form, shared by recovery and by every consumer that re-derives principals of
// a combined / averaged tensor.
void principal2(double a, double b, double c, double& p1, double& p2, double& theta);

// Recovered per-unit-width stress resultants at the element centroid.
// Stress-based convention: mx produces sigma_x, positive = tension on +z.
struct ShellResultants {
  double nx = 0, ny = 0, nxy = 0;      // membrane forces [N/m]
  double mx = 0, my = 0, mxy = 0;      // bending moments [N.m/m]
  double qx = 0, qy = 0;               // transverse shears [N/m]
  double m1 = 0, m2 = 0, m_theta = 0;  // principal moments + angle [rad]
  double n1 = 0, n2 = 0, n_theta = 0;  // principal membrane forces + angle
  // Recompute m1/m2/n1/n2 and the angles from the eight tensor components.
  // Principals are NEVER combined/enveloped linearly — they are recomputed.
  void recompute_principals();
};

// The element kernel: geometry + material + the 18x18 global stiffness, mass
// and geometric stiffness, and the data recovery needs. DOF order per node
// [ux, uy, uz, rx, ry, rz], nodes in winding order.
// Strain-displacement operators of a facet for MATERIAL-NONLINEAR
// integration (keta/shell_plastic.h), local frame: membrane (Allman) B at the
// three mid-edge points of the membrane rule (weights A/3), the constant
// bending B (curvature (kx, ky, kxy)), and the parts that stay elastic — the
// MITC3 transverse shear block and the drilling stabilisation penalty, both
// on their 9-DOF blocks. mem_idx / bend_idx place the 9-DOF blocks in the
// local 18-DOF vector; T maps global element DOFs to local ones. For an
// elastic material, sum_p w_p Bm' (t_p Dm) Bm + sum_p w_p t_p^3/12 Bb' Dm Bb
// + shear + drill (+ the offset and the bed) is ke() exactly.
struct ShellNlOps {
  // the in-plane integration points: weight w (area share), thickness t,
  // Allman membrane B (membrane DOFs) and MITC3 assumed transverse shear B
  // (gxz, gyz; bending DOFs) at each. A plain facet: the three mid-edge
  // points, w = A/3 (sum_p w Bs' (k G t) Bs is k_shear exactly); an enriched
  // facet whose thickness varies: the 7-point Dunavant rule (degree 5: t
  // linear times the quadratic B'DB, t^3 of the bending)
  std::vector<Eigen::Matrix<double, 3, 9>> bm;
  std::vector<Eigen::Matrix<double, 2, 9>> bs;
  std::vector<double> w, tp;
  Eigen::Matrix<double, 3, 9> bb;
  Eigen::Matrix<double, 9, 9> k_drill, k_shear;
  std::array<int, 9> mem_idx{}, bend_idx{};
  Mat18 T;
  // enriched facets: the offset (node DOFs -> mid-surface DOFs, local) and
  // the elastic Winkler bed on the node DOFs (local)
  Mat18 te = Mat18::Identity(), k_bed = Mat18::Zero();
  bool offset = false;
  double area = 0, t = 0;
  // Appended: an orthotropic facet or a laminate — the elastic matrix of
  // each ply (sx, sy, txy, txz, tyz from ex, ey, gxy, gxz, gyz; local axes)
  // between the unit heights z0 < z1 (z / t from the nodes' plane: a
  // laminate off the mid-plane spans its own range); empty: isotropic E, nu
  struct Ply {
    double z0 = -0.5, z1 = 0.5;
    Eigen::Matrix<double, 5, 5> C = Eigen::Matrix<double, 5, 5>::Zero();
  };
  std::vector<Ply> plies;
  bool laminate = false;
};

class ShellElement {
 public:
  ShellElement(const std::array<Vec3, 3>& coords, double t, double E, double nu,
               const std::optional<Vec3>& local_x = std::nullopt, const ShellEnrich& enrich = {});
  const ShellEnrich& enrichment() const { return en_; }

  const Mat18& ke() const { return ke_; }
  double area() const { return area_; }
  double thickness() const { return t_; }
  const Mat3& R() const { return r_; }
  Mat18 local_transform() const;  // global element DOFs -> local (18x18 block of R)
  // The in-plane strain (ex, ey, gxy; local axes) of the displacement d
  // (global element DOFs) at the point whose barycentric coordinates are
  // l, at the height z above the mid-surface: Allman membrane + z kappa.
  Vec3 strain_at(const Vec18& d_global, const std::array<double, 3>& l, double z) const;
  Vec3 normal() const { return r_.row(2); }
  const std::array<Vec3, 3>& coords() const { return coords_; }

  // Consistent (or lumped) 18x18 mass for density rho [kg/m^3]: rho*t*A
  // translational + rho*t^3/12 rotary on the bending rotations, NO drilling
  // mass (see the header physics block). GLOBAL element DOF order.
  Mat18 me(double rho, bool lumped = false) const;

  // Membrane geometric stiffness for the prestress membrane resultants
  // (nx, ny, nxy) [N/m] of this element (tension-positive, same convention as
  // the recovered resultants): the exact 1-point von Karman
  // grad(w)' [n] grad(w) block on the transverse DOFs. GLOBAL DOF order.
  Mat18 kg(double nx, double ny, double nxy) const;

  // --- consistent load vectors (GLOBAL element DOF order) ---
  // Uniform pressure given as a global force-per-area vector [N/m^2].
  Vec18 pressure_load_global(const Vec3& q) const;
  // Uniform pressure `p` [N/m^2] along the element's +local normal.
  Vec18 pressure_load_normal(double p) const;
  // Hydrostatic pressure gamma x depth below the free surface through x0
  // (depth along `down`, zero above), along side x the normal: the exact
  // consistent nodal forces int p N_i dA, the facet clipped to its wet part
  // (depths at the node plane).
  Vec18 hydrostatic_load(double gamma, const Vec3& x0, const Vec3& down, double side) const;
  // Plan-area ("projected"/snow) pressure: the global vector `q` is per m^2 of
  // HORIZONTAL projection, so the element receives q * |n.z|. A vertical panel
  // receives exactly zero.
  Vec18 pressure_load_projected(const Vec3& q) const;
  // Self-weight: force/area = rho*t*g, tributary A/3 per node.
  Vec18 self_weight_load(double rho, const Vec3& g) const;
  // Thermal initial-strain load; eps0 / kappa0 already multiplied by alpha.
  // Exact (both B integrate to their centroid values). Self-equilibrated.
  Vec18 thermal_load(const Vec3& eps0, const Vec3& kappa0) const;

  // Centroid resultants from the 18 GLOBAL nodal displacements, SUBTRACTING the
  // thermal initial strains (pass zeros for a non-thermal case).
  ShellResultants recover(const Vec18& d_global, const Vec3& eps0 = Vec3::Zero(),
                          const Vec3& kappa0 = Vec3::Zero()) const;

  // operators for material-nonlinear integration (see ShellNlOps)
  ShellNlOps nl_ops() const;

 private:
  std::array<Vec3, 3> coords_;
  double t_ = 0, area_ = 0;
  Mat3 r_ = Mat3::Identity();
  Mat18 ke_ = Mat18::Zero();
  // recovery data, local frame
  Eigen::Matrix<double, 3, 9> bm_;   // Allman membrane B at the centroid
  Eigen::Matrix<double, 3, 9> bb_;   // bending B
  Eigen::Matrix<double, 2, 9> bsc_;  // shear B at the centroid
  Eigen::Matrix<double, 2, 3> gw_;   // grad(w) operator (constant) — feeds kg
  Mat3 dm_ = Mat3::Zero(), db_ = Mat3::Zero();
  double ks_coef_ = 0;  // kappa_s * G * t
  ShellEnrich en_;
  double tc_ = 0;                                       // thickness at the centroid
  Eigen::Matrix2d gs_ = Eigen::Matrix2d::Identity();    // transverse shear moduli (local axes)
  Mat18 te_ = Mat18::Identity();                        // node DOFs -> mid-surface DOFs (local)
  std::array<double, 3> tcorner_{};
  std::array<std::array<double, 2>, 3> xy_{};  // local nodal coordinates
  double G_ = 0;
  double gd_ = 0;  // in-plane shear modulus of the drilling stabilisation
  // laminate (G11): A, B, D about the mid-plane
  bool lam_ = false;
  Mat3 lam_A_ = Mat3::Zero(), lam_B_ = Mat3::Zero(), lam_D_ = Mat3::Zero();
};

// ---------------- model-side plumbing ----------------

// One lowered facet shell: the kernel plus the dense node indices and the
// per-element bookkeeping the assembly / load / recovery paths need. Built in
// `model.shells` order (assembly determinism is load-bearing).
struct ShellElem {
  std::size_t index = 0;  // dense index into model.shells
  Id id = 0;
  Id slab = 0;
  std::array<int, 3> nodes{};  // dense node indices
  double rho = 0;              // mass density [kg/m^3] (self-weight + mass)
  double alpha = 0;            // thermal dilation [1/K] (SurfaceThermal)
  ShellElement el;
};

// Lower every ShellTri of the model. Typed errors on dangling references and
// degenerate / non-positive geometry, naming the shell id.
std::vector<ShellElem> shell_elems(const Model& model);

// Per-shell thermal initial strain / curvature of one lowered case, in the
// element local frame. Zeroed in `model.shells` order and ACCUMULATED as
// SurfaceThermal loads are lowered, so two thermal loads on one slab compose.
struct ShellRecovery {
  Vec3 eps0 = Vec3::Zero();
  Vec3 kappa0 = Vec3::Zero();
  Vec18 dstar = Vec18::Zero();  // Load::ShellImposed: recovered from d - d*
  bool any() const { return eps0.squaredNorm() > 0 || kappa0.squaredNorm() > 0 || dstar.squaredNorm() > 0; }
};

// One element's recovered resultants for one case / combination / envelope
// bound (parallel to BarResult).
struct ShellResult {
  Id shell = 0;
  Id slab = 0;
  ShellResultants r;
};

// ---------------- line supports ----------------

// Ids of every model node lying within `tol` of the 3-D segment [a, b]
// (distance to the segment, endpoints included). Deterministic: model node
// order. A zero-length segment is a typed error — a "line support" degenerating
// to a point silently masking one node is exactly the trap this rejects.
std::vector<Id> nodes_on_segment(const Model& model, const Vec3& a, const Vec3& b,
                                 double tol);

// Convenience: add one Support with the given fixed mask on every node within
// `tol` of the segment [a, b] (a line support = a mask applied per node).
// Returns the node ids that received it.
std::vector<Id> add_line_support(Model& model, const Vec3& a, const Vec3& b,
                                 const std::array<bool, 6>& fixed, double tol);

}  // namespace keta

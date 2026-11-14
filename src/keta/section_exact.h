// Exact cross-section engine — Saint-Venant warping/torsion finite elements on a
// quadratic (Tri6) triangulation of an arbitrary holed, multi-material polygon,
// plus the ASCII-DXF outline importer.
//
// This module EXTENDS `keta/sections.h`: the closed-form `SectionProps` layer stays
// the default (it is instantaneous and covers the catalogue shapes), and an exact
// result is ATTACHED to a `Section` on demand (`attach_exact`) when the heuristic
// J / Avy / Avz of a generic polygon are not good enough, or when the warping
// constant Iw, the shear centre or a real stress field are needed.
//
// ---------------------------------------------------------------------------
// Formulation (Saint-Venant torsion + Timoshenko shear), references
// ---------------------------------------------------------------------------
// Notation: the section lies in the local (y, z) plane, x is the bar axis. All
// integrals below are over the section area, in CENTROIDAL coordinates, and are
// MODULUS-WEIGHTED for a composite (every element integral carries its own
// n_e = E_e / E_ref, i.e. the classic transformed-section weighting).
//
// 1. Warping function.  Saint-Venant torsion assumes u_x = theta' * omega(y, z),
//    u_y = -theta' x z, u_z = theta' x y.  Equilibrium reduces to the Laplace
//    problem with a Neumann condition on the traction-free boundary
//
//        div(n grad omega) = 0 in A,        n grad omega . nu = n (z n_y - y n_z) on dA
//
//    whose weak form is  integral n grad(w) . grad(omega) dA
//                      = integral n grad(w) . (z, -y) dA   for all test w.
//    The torsion constant then follows from the strain energy
//
//        J = Iy + Iz - omega^T f_tors                                        (1)
//
//    with f_tors the discrete right-hand side above.  `omega` is defined up to a
//    constant; one degree of freedom is pinned (a rigid-body warping shift, which
//    (1) is invariant to because sum(f_tors) = 0).
//    References: Timoshenko & Goodier, *Theory of Elasticity*, 3rd ed. ch. 10;
//    Pilkey, *Analysis and Design of Elastic Beams* (2002) ch. 5;
//    Gruttmann & Wagner, "Shear correction factors in Timoshenko's beam theory
//    for arbitrary shaped cross-sections", Comput. Mech. 27 (2001) 199-207.
//
// 2. Transverse shear.  Two further Laplace-type problems (right-hand sides
//    depending on nu and on the second moments) give the Prandtl-type potentials
//    psi (for Vy) and psi_z (for Vz); the shear stresses are
//
//        tau = Vy/Delta * (grad psi_y - nu/2 * d)  +  Vz/Delta * (grad psi_z - nu/2 * h)
//        Delta = 2 (1 + nu) (Iy Iz - Iyz^2)                                   (2)
//
//    For a composite the nu/2 terms take each element's own nu (compatibility)
//    while the source term keeps the global 2 (1 + nu_ref) of Delta
//    (equilibrium: div tau = n (Iy y - Iyz z) V / ip2), which keeps the Neumann
//    load self-equilibrated, so the result does not depend on the pinned DOF.
//    The ENERGY-CONSISTENT (Timoshenko) shear areas come from the shear
//    strain energy of that field,  Avy = Delta^2 / kappa_y  (Gruttmann & Wagner
//    eq. 22; the same form `sectionproperties` implements after Pilkey).  For a
//    rectangle this yields ~0.83 A, i.e. the classical 5/6 to better than 1 %.
//
// 3. Shear centre.  Both the elasticity (CTM) definition, from the coupling of
//    the shear potentials with the torsion load, and the geometric (Trefftz)
//    definition, from the sectorial products of omega, are recovered; `y_sc`,
//    `z_sc` below are the elasticity ones (the pair that makes the shear field
//    torsion-free).  Warping constant about the shear centre:
//
//        Iw = integral omega^2 dA - (integral omega dA)^2 / A
//             - z_sc integral y omega dA + y_sc integral z omega dA            (3)
//
// 4. Discretisation.  Six-node quadratic triangles with STRAIGHT sides (the
//    Jacobian is then constant over the element), 4-point symmetric Gauss rule
//    (exact to cubic; the centroid weight -27/48 is negative, which is correct
//    for this rule).  The assembled matrix is symmetric positive definite once
//    one DOF is pinned and is factored ONCE with Eigen's `SimplicialLDLT`; the
//    three right-hand sides (torsion, Vy, Vz) are solved against that single
//    factorisation.
//
// 5. Stresses.  Straight-sided Tri6 fields are exactly-representable quadratics,
//    so the per-action stress fields are evaluated directly (no extrapolation).
//    sigma_xx = N/A + (My (Iz z - Iyz y) - Mz (Iy y - Iyz z)) / (Iy Iz - Iyz^2)
//               + B omega_sc / Iw  (bimoment, Vlasov)
//    tau      = T/J (grad omega + (-z, y)) + the shear field of (2).
//
// ---------------------------------------------------------------------------
// Units, conventions, determinism
// ---------------------------------------------------------------------------
//   * SI throughout: coordinates and lengths in m, areas m^2, second moments
//     m^4, warping constant m^6, moduli Pa, actions N and N.m, stresses Pa.
//   * Section frame: y horizontal, z vertical, right-handed with the bar axis x;
//     positive My bends fibres at +z into tension (sigma = ... + My z / Iy),
//     positive Mz puts +y into compression (sigma = ... - Mz y / Iz), matching
//     `keta/sections.h`.
//   * Every public evaluator takes coordinates RELATIVE TO THE ELASTIC CENTROID,
//     like `StressPoint`.  `ExactSection::yc, zc` give the centroid in the input
//     frame so a caller can convert.
//   * Deterministic: the mesher and the solver contain no randomness, no
//     address-derived ordering and no parallelism; the same input always yields
//     the same node order and bit-identical properties.
//   * Failures are typed (`keta::fail`), never a silent NaN.
//
// ---------------------------------------------------------------------------
// Known limitations
// ---------------------------------------------------------------------------
//   * The triangulator of this engine (`mesh_section`) requires the target
//     element size to resolve the thinnest wall and it VERIFIES that it did
//     (typed error otherwise). The stress evaluators (section_stress.cpp)
//     mesh with `keta::mesh2` first and fall back on it.
//   * Composite (modular-ratio) sections are weighted by E only: the recovered
//     J and the shear areas are the transformed-section values relative to
//     E_ref and are the true rigidities divided by G_ref only when every
//     material shares the same Poisson ratio.
//   * A section made of several DISJOINT regions is rejected: Saint-Venant
//     torsion is not defined component-by-component about a common centroid.
//   * Curved boundaries must be supplied polygonised; the Tri6 elements are
//     straight-sided (sub-parametric), so a circle's J carries the polygon's
//     own geometric error (~pi^2/(3 n^2) on the area for an inscribed n-gon).
#pragma once
#include <array>
#include <memory>
#include <string>
#include <vector>

#include "keta/mesh2d.h"
#include "keta/sections.h"

namespace keta {

using Vec2 = Eigen::Vector2d;  // section-plane point (y, z) [m]

// ---------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------

// One material of a composite section. Only the ratio E/E_ref matters for the
// geometry-like properties; `nu` enters the shear problems per element.
struct ExactMaterial {
  std::string name;
  double E = 1.0;   // Young's modulus [Pa]
  double nu = 0.3;  // Poisson's ratio [-]
};

// A material patch: every element whose centroid falls inside `outline` takes
// `material`. Later regions win over earlier ones; elements matching no region
// keep material 0.
struct ExactRegion {
  std::vector<Vec2> outline;  // (y, z) [m], either orientation
  int material = 0;
};

struct ExactInput {
  std::vector<Vec2> outline;              // outer boundary, either orientation
  std::vector<std::vector<Vec2>> holes;   // inner boundaries
  std::vector<ExactMaterial> materials;   // empty => one default material
  std::vector<ExactRegion> regions;       // empty => uniform material 0
  double target_h = 0;                    // target element edge [m]; 0 => bbox/20
  int max_nodes = 400000;                 // typed reject above this
  bool smooth = true;                     // one smart-Laplacian pass on interior nodes
};

// ---------------------------------------------------------------------------
// Mesh
// ---------------------------------------------------------------------------

// Quadratic (6-node) triangle mesh. Node order per element is
// [c0, c1, c2, m01, m12, m20] with the corner triangle counter-clockwise, which
// matches the shape functions N1..N6 = eta(2eta-1), xi(2xi-1), zeta(2zeta-1),
// 4 eta xi, 4 xi zeta, 4 eta zeta (zeta = 1 - eta - xi).
struct Tri6Mesh {
  std::vector<Vec2> nodes;                    // (y, z) [m]
  std::vector<std::array<int, 6>> tris;
  std::vector<int> mat;                       // material index per triangle
  std::vector<int> boundary_nodes;            // corner nodes on the domain boundary
  double area() const;                        // sum of element areas [m^2]
};

// The triangulator of analyze_exact: a Bowyer-Watson Delaunay triangulation
// of (boundary samples at ~target_h) + (a triangular interior lattice at
// target_h kept clear of the boundary), filtered to the domain by a centroid
// point-in-polygon test. Under the stated sampling condition every boundary
// sub-edge is a GABRIEL edge and therefore belongs to the Delaunay
// triangulation, so the filtered mesh conforms to the boundary exactly; the
// mesher VERIFIES that (edge recovery + mesh area against the exact Green
// integral) and fails typed when the target size does not resolve a wall.
Tri6Mesh mesh_section(const ExactInput& in);

// Forward-compatible seam: promote a linear mesh2 triangulation to Tri6 and tag
// its elements with `in`'s material regions.
Tri6Mesh tri6_upgrade(const mesh2::TriMesh& linear, const ExactInput& in);

// ---------------------------------------------------------------------------
// Actions and stresses
// ---------------------------------------------------------------------------

// Section actions in the bar local frame. `B` is the Vlasov bimoment [N.m^2].
struct SectionActions {
  double N = 0, Vy = 0, Vz = 0, T = 0, My = 0, Mz = 0, B = 0;
};

struct StressState {
  double sig = 0;                // sigma_xx [Pa]
  double tau_y = 0, tau_z = 0;   // transverse shear components [Pa]
  double tau = 0;                // hypot(tau_y, tau_z)
  double vm = 0;                 // von Mises sqrt(sig^2 + 3 tau^2)
};

struct SigmaExtreme {
  double value = 0;      // sigma_xx [Pa]
  double y = 0, z = 0;   // location, centroid-relative [m]
};
struct SigmaRange {
  SigmaExtreme min, max;
};

// ---------------------------------------------------------------------------
// Result
// ---------------------------------------------------------------------------

// Everything the warping solve produced. Node coordinates in `mesh` are stored
// SHIFTED to the elastic centroid, so every evaluator works in centroid-relative
// coordinates.
struct ExactSection {
  Tri6Mesh mesh;
  std::vector<ExactMaterial> materials;
  double E_ref = 1;   // materials[0].E — the modulus every property is divided by
  double nu_ref = 0;  // modulus-weighted mean Poisson ratio (shear normaliser)

  double A = 0;             // geometric area [m^2] (unweighted)
  double EA_ref = 0;        // transformed area = sum n_e dA [m^2]
  double yc = 0, zc = 0;    // elastic centroid, in the INPUT frame [m]
  double Iy = 0, Iz = 0, Iyz = 0;  // centroidal transformed second moments [m^4]
  double I1 = 0, I2 = 0;    // principal second moments [m^4]
  double phi = 0;           // principal angle [rad], atan2(Iy - I1, Iyz)
  double J = 0;             // Saint-Venant torsion constant [m^4]
  double Iw = 0;            // warping constant about the shear centre [m^6]
  double y_sc = 0, z_sc = 0;    // elasticity shear centre, centroid-relative [m]
  double y_st = 0, z_st = 0;    // Trefftz (geometric) shear centre [m]
  double Avy = 0, Avz = 0;      // energy-consistent shear areas [m^2]
  double Avyz = 0;              // shear coupling area [m^2]; 0 = uncoupled
  double zmax = 0, ymax = 0;    // extreme centroid-relative fibre distances [m]

  // Solved nodal fields (length = mesh.nodes.size()).
  std::vector<double> omega;   // warping function [m^2]
  std::vector<double> psi_y;   // shear potential driving Vy
  std::vector<double> psi_z;   // shear potential driving Vz
  double omega_mean = 0;       // integral n omega dA / EA_ref (zero-mean shift)
  double delta_s = 0;          // 2 (1 + nu_ref) (Iy Iz - Iyz^2)

  // Vlasov sectorial coordinate about the shear centre, zero-mean, at a node.
  // Sign convention: d(omega) = (y - y_sc) dz - (z - z_sc) dy, i.e. the classical
  // one the thin-wall module uses; the raw finite-element field `omega` is its
  // NEGATIVE (its gradient is (z, -y) + correction).
  double omega_sc(int node) const;

  // --- evaluators, all in CENTROID-RELATIVE coordinates ---

  // Index of the element containing (y, z), or -1 when the point is outside.
  int locate(double y, double z) const;
  // Full stress state at (y, z). Throws when the point is outside the section.
  StressState stress_at(double y, double z, const SectionActions& a) const;
  // sigma_xx only. Same element lookup and range check as `stress_at`.
  double sigma_at(double y, double z, const SectionActions& a) const;
  // Extremes of the LINEAR normal stress (N, My, Mz only) over the section.
  // sigma is affine in (y, z), so its extremes sit on the boundary and scanning
  // the boundary nodes — which contain every input vertex verbatim — is EXACT.
  SigmaRange sigma_extremes(const SectionActions& a) const;
  // Largest von Mises stress over the mesh nodes (approximate: nodal sampling).
  double von_mises_max(const SectionActions& a) const;
  // Nodal stress field, one state per mesh node: the mean of the incident
  // elements' values at that node (the transverse-shear field is
  // element-discontinuous; sigma is continuous within one material).
  std::vector<StressState> nodal_stresses(const SectionActions& a) const;
  // The six unit-action nodal fields [N, Vy, Vz, My, Mz, T] in one sweep (the
  // element geometry and interpolated potentials computed once per node, not
  // once per action); each equals nodal_stresses() of its unit action.
  std::array<std::vector<StressState>, 6> nodal_unit_stresses() const;

  // Fill a closed-form-compatible property block (transformed values; identical
  // to the geometric ones for a single material).
  SectionProps props() const;
};

// Mesh and solve in one call.
ExactSection analyze_exact(const ExactInput& in);
// Solve on a mesh built elsewhere (hand-made or from mesh2).
ExactSection analyze_exact(const Tri6Mesh& mesh, const std::vector<ExactMaterial>& materials);

// ---------------------------------------------------------------------------
// Wire-in to `Section` (smallest possible edit: one shared_ptr field)
// ---------------------------------------------------------------------------

// Solve `in` and attach the result to `s`, overwriting the heuristic torsion and
// shear areas of `s.props` with the exact ones (A / Iy / Iz stay as the shape's
// own closed form: they already are exact for every keta shape). Statics is
// untouched — a caller that wants exact stresses queries the evaluators below.
void attach_exact(Section& s, const ExactInput& in);
// Attach an already solved result (shared between bars using the same section).
void attach_exact(Section& s, std::shared_ptr<const ExactSection> ex);

// Energy-consistent shear areas (Saint-Venant flexure, Gruttmann & Wagner) of
// a parametric shape for Poisson's ratio `nu`: the closed form for a solid
// circle, the Tri6 warping solve on the exact outline otherwise (root fillets
// of an I included). False for ShapeGeneric (no geometry).
bool energy_shear_areas(const Shape& shape, double nu, double& Avy, double& Avz);

// Null when no exact backend was attached.
std::shared_ptr<const ExactSection> exact_of(const Section& s);

// Free stress-recovery hooks (throw when `s` carries no exact backend).
StressState exact_stress_at(const Section& s, double y, double z, const SectionActions& a);
SigmaRange exact_sigma_extremes(const Section& s, const SectionActions& a);

// ---------------------------------------------------------------------------
// DXF import (ASCII, R12 and later)
// ---------------------------------------------------------------------------
// Reads LWPOLYLINE, POLYLINE/VERTEX, LINE and ARC entities out of the ENTITIES
// section, chains loose LINEs/ARCs into closed loops, then classifies the loops
// by area: the largest is the outline, the rest are holes.  Splines, blocks and
// binary DXF are NOT supported and are typed errors, never silently
// straightened.  `scale` converts the file's length unit to metres (e.g. 1e-3
// for a drawing in millimetres).
//
// Curved edges — ARC entities (groups 10/20 centre, 40 radius, 50/51 start/end
// angles in degrees, always counter-clockwise) and LWPOLYLINE / POLYLINE-VERTEX
// bulges (group 42, bulge = tan(theta/4), positive counter-clockwise) — are
// POLYGONISED under a sagitta (mid-ordinate) tolerance: the angular step dphi
// satisfies  R (1 - cos(dphi/2)) <= tol,  i.e.  dphi = 2 acos(1 - tol/R),
// n = ceil(|sweep| / dphi) equal steps.  Junction nodes are exact (the arc end
// points are the original vertices / the exact angle-50/51 points, never a
// re-sampled duplicate).  The resulting inscribed polygon UNDER-estimates each
// arc's area by exactly  n * (R^2/2) * (dphi - sin dphi)  per arc (circular
// segment areas), a bound known in closed form.
struct DxfOptions {
  double scale = 1.0;  // drawing length unit -> metres
  // Sagitta tolerance for arc polygonisation, in DRAWING units (pre-`scale`):
  //   > 0  : absolute cap on the mid-ordinate flèche of every chord;
  //   <= 0 : DEFAULT, relative — tol = 1e-3 * R per arc.  Chosen relative so
  //          the import is invariant to the drawing unit and the polygon area
  //          error stays ~dphi^2/6 ≈ 1.3e-3 of each arc's sector area
  //          (dphi = 2 acos(1 - 1e-3) ≈ 0.0894 rad, ~71 chords per full turn)
  //          — well below the Tri6 meshing error of the downstream engine.
  double arc_sagitta = 0.0;
};
struct DxfSection {
  std::vector<Vec2> outline;             // CCW [m]
  std::vector<std::vector<Vec2>> holes;  // CCW [m]
};
DxfSection import_dxf(const std::string& text, const DxfOptions& opt);
DxfSection load_dxf(const std::string& path, const DxfOptions& opt);
// Historical entry points: identical to {scale, default sagitta} (arc-free
// files import byte-identically to the pre-arc importer).
DxfSection import_dxf(const std::string& text, double scale = 1.0);
DxfSection load_dxf(const std::string& path, double scale = 1.0);

// ---------------------------------------------------------------------------
// Small geometry helpers shared with the thin-wall module
// ---------------------------------------------------------------------------
double polygon_area(const std::vector<Vec2>& ring);         // signed, CCW positive
bool point_in_polygon(const std::vector<Vec2>& ring, const Vec2& p);

}  // namespace keta

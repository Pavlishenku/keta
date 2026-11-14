// Thin-walled cross-sections: piecewise-straight midline with a thickness per
// wall — open (Vlasov) and closed (Bredt) single- or multi-cell.
//
// This is the analytical companion of `keta/section_exact.h`: the exact Tri6
// engine is the reference, this module is the closed-form theory used for fast
// design work and as an independent cross-check of it.
//
// ---------------------------------------------------------------------------
// Model
// ---------------------------------------------------------------------------
// The section is a planar GRAPH: nodes are midline points (y, z), segments are
// straight walls of constant thickness t.  Every wall is idealised as a
// rectangle of length L and thickness t centred on the midline, so
//
//   A_s = t L,   and the exact strip second moments about the global axes are
//   Iy_s = A_s (z_m^2 + (L^2 e_z^2 + t^2 e_y^2)/12),  e = wall unit direction
//   Iz_s = A_s (y_m^2 + (L^2 e_y^2 + t^2 e_z^2)/12)
//   Iyz_s = A_s (y_m z_m + (L^2 - t^2) e_y e_z / 12)
//
// Overlapping material at a junction is counted twice — the classical thin-wall
// idealisation; it is a ~(t/L) effect and it is the reason a thin-wall J of a
// rolled I-section runs BELOW the finite-element value (no junction fillet, no
// junction stiffening).  Documented, never hidden.
//
// ---------------------------------------------------------------------------
// Torsion — Bredt, any topology
// ---------------------------------------------------------------------------
// The independent cells are the CYCLE SPACE of the graph: with V nodes and E
// segments in one connected component there are n = E - V + 1 of them, and the
// fundamental cycles of any spanning tree are a basis.  Both sides of the
// compatibility relation are linear in the circuit, so working in that basis is
// exactly equivalent to working with the planar faces (and needs no planar
// embedding).  With a unit twist rate G theta' = 1 the circulating flows x solve
//
//     M x = 2 A,     M_kl = sum_s c_k[s] c_l[s] L_s / t_s,   A_k = cell area   (1)
//
// (c_k the +-1 incidence of cycle k), the wall flow is q = sum_k x_k c_k and
//
//     J_closed = T(x) = 2 sum_k x_k A_k,        J = J_closed + sum_open b t^3/3 (2)
//
// which reduces to Bredt's 4 A_m^2 / (contour integral ds/t) for a single cell.
// Walls that are BRIDGES of the graph (in no cycle) carry the open Saint-Venant
// term b t^3/3; walls inside cells carry the closed term only.
// References: Bredt (1896); Kollbrunner & Basler, *Torsion in Structures*;
// Murray, *Introduction to the Theory of Thin-Walled Structures*.
//
// ---------------------------------------------------------------------------
// Transverse shear — one square system for open and closed alike
// ---------------------------------------------------------------------------
// Wall equilibrium gives dq/ds = -(Kz z + Ky y) t with
//
//     Kz = (Vz Iz - Vy Iyz)/D,  Ky = (Vy Iy - Vz Iyz)/D,  D = Iy Iz - Iyz^2  (3)
//
// i.e. the familiar q = V S(s) / I once the axes are principal.  Taking the E
// wall-start flows as unknowns, the V-1 independent node balances plus the n
// zero-twist circuit conditions (a shear through the shear centre twists
// nothing) make a SQUARE system — which degenerates to the classical
// leaf-to-root accumulation for an open section and to the classical "cut the
// cells and restore compatibility" for a closed one, with no special-casing.
// The energy-consistent shear areas follow from the flow energy
//
//     Av = V^2 / (contour integral q^2/t ds).                                (4)
//
// ---------------------------------------------------------------------------
// Warping — Vlasov
// ---------------------------------------------------------------------------
// The sectorial coordinate about a pole P integrates
//
//     d(omega) = (y - y_P) dz - (z - z_P) dy - (q_1/t) ds                    (5)
//
// where q_1 is the UNIT-TWIST torsion flow of (1); the correction term vanishes
// for an open section and makes omega single-valued around every cell of a
// closed one (its circuit integral is 2A_k - contour(q_1/t) = 0 by (1)), so (5)
// covers single- AND multi-cell sections.  omega is linear along each straight
// wall.  The shear centre nulls the sectorial products,
//
//     [Iyz  -Iz ; Iy  -Iyz] (y_P - y_S, z_P - z_S)^T = -(I_omega_y, I_omega_z) (6)
//
// and, once omega is re-referred to that pole and given zero mean,
//
//     Iw = integral omega_S^2 dA,   sigma_w = B omega_S / Iw.                (7)
//
// References: Vlasov, *Thin-Walled Elastic Beams*, 2nd ed. (1961);
// Timoshenko & Gere, *Theory of Elastic Stability*, ch. 5.
//
// ---------------------------------------------------------------------------
// Units, conventions, limits
// ---------------------------------------------------------------------------
//   * SI: m, m^2, m^4, m^6, N, N.m, Pa. Section frame (y, z) as in sections.h.
//   * Station coordinates are CENTROID-RELATIVE, like the exact engine's.
//   * Deterministic, allocation-light, no dependency beyond core.h.
//   * The graph must be connected; a disconnected midline is a typed error.
//   * Second-order (t/L) effects, junction fillets and shear lag are outside the
//     model by construction: cross-check against `analyze_exact` when they matter.
//
// ---------------------------------------------------------------------------
// Linearly tapered walls
// ---------------------------------------------------------------------------
// A wall may carry a thickness that varies linearly from node a to node b
// (`WallSegment::t_b`), symmetric about the midline. Every formula above then
// holds with t = t(s) under its integral: the area and the second moments are
// polynomial integrals (the strip term becomes integral t^3 / 12), Bredt's
// contour integral of ds/t is L ln(t_b / t_a) / (t_b - t_a), the open
// torsion term is integral t^3 ds / 3 = L (t_a + t_b)(t_a^2 + t_b^2) / 12, the
// shear flow is cubic along the wall, and inside a cell the sectorial
// coordinate follows integral ds/t (a logarithm) instead of a straight line.
// The rational and logarithmic integrals (flow energy, zero-twist conditions,
// sectorial moments) are taken by Gauss-Legendre quadrature on pieces of the
// wall whose thickness ratio stays under 1.5 — exact to rounding.
// TorsionFlow::tau_seg reports the stress at the thin end, tau_open at the
// thick end (their largest values on the wall).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "keta/core.h"

namespace keta {

using Vec2 = Eigen::Vector2d;  // midline point (y, z) [m]

// One straight wall between two midline nodes.
struct WallSegment {
  int a = 0, b = 0;  // node indices
  double t = 0;      // wall thickness [m]; at node a when the wall tapers
  // --- appended ---
  double t_b = 0;    // thickness at node b of a linearly tapered wall; 0 = constant t
  double thickness_at_b() const { return t_b > 0 ? t_b : t; }
  // thickness at the fraction u of the length from node a
  double thickness(double u) const { return t + u * (thickness_at_b() - t); }
};

struct ThinWall {
  std::vector<Vec2> nodes;         // midline nodes (y, z) [m]
  std::vector<WallSegment> segs;
  void validate() const;           // typed rejects: finiteness, t > 0, connectivity
  double length(int seg) const;
};

// One independent cell: the fundamental cycle of a spanning tree, as a signed
// walk over segments (`dir` = +1 when traversed a -> b).
struct WallCycle {
  std::vector<int> seg;
  std::vector<int> dir;
  double area = 0;  // signed enclosed area [m^2]
};

struct ThinWallProps {
  double A = 0;                    // [m^2]
  double yc = 0, zc = 0;           // centroid in the input frame [m]
  double Iy = 0, Iz = 0, Iyz = 0;  // centroidal [m^4]
  double I1 = 0, I2 = 0, phi = 0;  // principal moments [m^4] and angle [rad]
  double J = 0, J_open = 0, J_closed = 0;  // [m^4]
  int n_cells = 0;
  double y_sc = 0, z_sc = 0;       // shear centre, centroid-relative [m]
  double Iw = 0;                   // warping constant about the shear centre [m^6]
  double Avy = 0, Avz = 0;         // shear areas from the flow energy [m^2]
  double Wely = 0, Welz = 0;       // elastic moduli over the midline extremes [m^3]
};

// A sampling point on the midline. `s` is the arc length from the segment's
// first node; y, z are CENTROID-RELATIVE.
struct WallStation {
  int seg = 0;
  double s = 0;
  double y = 0, z = 0;
  double q = 0;        // shear flow [N/m], positive along a -> b
  double tau = 0;      // q / t [Pa]
  double omega = 0;    // sectorial coordinate about the shear centre [m^2]
  double sigma_w = 0;  // warping normal stress B omega / Iw [Pa]
};

// Full property set (mesh-free, closed form).
ThinWallProps thinwall_props(const ThinWall& w);

// The independent cells (empty for an open section).
std::vector<WallCycle> thinwall_cycles(const ThinWall& w);

struct TorsionFlow {
  double J = 0, J_open = 0, J_closed = 0;
  std::vector<double> q_cell;  // circulating flow per cell [N/m]
  std::vector<double> q_seg;   // resulting wall flow [N/m], positive along a -> b
  std::vector<double> tau_seg; // closed-cell shear stress q/t [Pa]
  std::vector<double> tau_open;// open (Saint-Venant) peak G theta' t per wall [Pa]
};
// Bredt cell flows and stresses under a torque T. Walls that belong to no cell
// carry the open response only; walls in a cell carry the closed one.
TorsionFlow thinwall_torsion(const ThinWall& w, double T);

// Shear flow q(s) under Vy / Vz applied at the shear centre. Open and closed
// alike (see the header). `per_segment` >= 2 stations per wall.
std::vector<WallStation> thinwall_shear_flow(const ThinWall& w, double Vy, double Vz,
                                             int per_segment = 5);

// Sectorial coordinate and the warping normal stress under a bimoment B.
std::vector<WallStation> thinwall_warping(const ThinWall& w, double B, int per_segment = 2);

// ---------------------------------------------------------------------------
// Morphology: the wall graph from the contour of a section
// ---------------------------------------------------------------------------
// The midline of a contour is not found unaided (pairing the two faces of a
// wall is ambiguous at junctions): the caller says how the contour points
// group into walls, in one of two ways, and gets the midline nodes and the
// thickness of each wall end — a tapered WallSegment when the ends differ.
// A wall without a determinate midline is a typed error: degenerate length or
// thickness, or (by quads) two faces whose lengths differ by more than half
// the longer. The thicknesses of the two ends may differ by any ratio.
//
// By quads: each wall is four contour points [p, q, r, s] in order around
// it, p -> q and s -> r its two faces, p-s and q-r its end caps. The midline
// runs between the cap midpoints; the thickness of an end is its cap projected
// on the normal of the midline (a mitred corner gives the true thickness, not
// the diagonal). Cap midpoints closer than `snap` merge into one node
// (snap = 0: 1e-6 of the extent of the contour).
ThinWall thinwall_from_quads(const std::vector<Vec2>& contour, const std::vector<std::array<int, 4>>& walls,
                             double snap = 0.0);
// By nodes and branches: a node lists the contour points gathered at one
// junction (at least two), a branch joins two nodes. Three steps: (1) the
// centroid of the points of each node gives the axis of each branch; (2) the
// points of a node projected on the normal of that axis give the thickness of
// the branch at that end (their spread) and the point of its midline there
// (the middle of the spread); (3) each node moves to the least-squares
// intersection of the midlines of its branches — the mitre point of a
// corner; a node with one branch, or with parallel branches, keeps the mean
// of its midline points.
ThinWall thinwall_from_nodes(const std::vector<Vec2>& contour, const std::vector<std::vector<int>>& nodes,
                             const std::vector<std::array<int, 2>>& branches);

}  // namespace keta

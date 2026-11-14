// Concentrated plastic hinges and the EXACT event-to-event pushover of a
// frame carrying them (a bilinear moment-rotation work law at a bar end;
// pushover under a fixed load pattern scaled by lambda).
//
// Hinge: a zero-length rotational spring in series between the node and the
// bar end, about one local bending axis, bilinear with kinematic hardening
// (initial stiffness k1, yield moment My, post-yield stiffness k2 < k1):
//   M = k1 (theta - theta_p),  |M - beta| <= My,  beta_dot = H theta_p_dot,
// H = k1 k2 / (k1 - k2) (post-yield slope k1 H / (k1 + H) = k2) — the exact
// bilinear kinematic law, symmetric, unloading with k1.
//
// Element: the bar's own (release-condensed, local) stiffness with the hinged
// end rotations replaced by internal DOFs t_h, plus the springs
// (u_q - t_h)^T K_h (u_q - t_h) / 2 (tangent of the hinge), then the
// internal DOFs condensed statically —
// the exact semi-rigid member, whatever the bar kernel (Timoshenko,
// rigid arms). Hinge rotation recovered as theta_h = u_q - t_h.
//
// Solver: with first-order kinematics and piecewise-linear hinges the
// response to P(lambda) = lambda P is piecewise LINEAR in lambda. Between
// two events the tangent is constant: one factorization, one solve gives the
// rates, and the next event (an elastic hinge reaching its yield bound, i.e.
// a linear equation in lambda) is located in closed form. Plastic hinges
// whose rotation rate reverses unload elastically (the active set is
// re-solved until consistent: plastic <=> loading). No Newton iteration, no
// step-size error: every point of the curve is exact to round-off. Cost per
// event: one sparse factorization + O(#hinges).
//
// AXIAL FORCE - MOMENT INTERACTION (NMHinge): the hinge then has two
// deformations at the bar end, the axial slip and the rotation (springs kN,
// kM), and a CONVEX POLYGONAL yield surface in the (N, M) plane (facets
// aN N + aM M <= r), associated flow (the plastic slip/rotation rate is
// normal to the active facets) and linear kinematic hardening H (2 x 2,
// SPD or zero). Multi-surface plasticity (Koiter): with k active facets
// (rows n of N), lambda_dot = G^-1 N Ke e_dot, G = N (Ke + H) N^T, tangent
// Ke - Ke N^T G^-1 N Ke. The model is piecewise LINEAR — between events
// the tangent is constant and an inactive facet is reached when a linear
// function of lambda vanishes — so the event-to-event solution stays EXACT;
// a polygon reaching a vertex (two facets) is a corner of the surface, not
// an approximation of it. A curved interaction curve is represented by the
// polygon the caller chooses (InteractionSurface::from_vertices), never
// silently. The bilinear moment hinge is the special case of one deformation
// and the two facets M = +-My (H = k1 k2 / (k1 - k2)).
//
// PERMANENT LOAD: an optional case applied first (its own event-to-event
// path, lambda_perm 0 -> 1), then held while the controlled case grows —
// gravity under a lateral pushover, the axial force of an N-M hinge.
//
// SPAN LOADS: on a hinged bar the fixed-end forces split between the nodes
// and the internal hinge DOFs — equivalent nodal load f_b - K_bi K_ii^-1 f_i
// and hinge-deformation offset -K_ii^-1 f_i, with the current hinge tangents —
// exact in rate form between events (the response stays piecewise linear).
// Unhinged bars take their span loads as linear statics does.
//
// Typed rejections: prescribed displacements, imposed strains / dislocations
// on a hinged bar, tapered or cable bars carrying a hinge, a hinge on a
// released DOF.
#pragma once
#include <optional>
#include <utility>
#include <vector>

#include "keta/statics.h"

namespace keta {

struct BilinearHinge {
  double k1 = 0;  // initial rotational stiffness [N m/rad]
  double My = 0;  // yield moment [N m]
  double k2 = 0;  // post-yield stiffness [N m/rad], 0 <= k2 < k1
};

// Convex polygon in the (N, M) plane containing the origin (N tension
// positive at the hinged bar end, M the hinge moment): facets
// (aN, aM, r) meaning aN N + aM M <= r, r > 0.
struct InteractionSurface {
  std::vector<Vec3> facets;
  // from the vertices (N, M) of a convex polygon around the origin, in order
  static InteractionSurface from_vertices(const std::vector<std::pair<double, double>>& v);
  // |N| / Npl + |M| / Mpl <= 1 (linear interaction, distinct tension and
  // compression capacities allowed)
  static InteractionSurface diamond(double Npl_tension, double Npl_compression, double Mpl);
};

struct NMHinge {
  double kN = 0, kM = 0;              // elastic axial [N/m] and rotational [N m/rad] stiffness
  InteractionSurface surface;
  Eigen::Matrix2d H = Eigen::Matrix2d::Zero();  // kinematic hardening in (N, M) [N/m, N; N, N m]
};

struct HingeSpec {
  Id bar = 0;
  int end = 0;   // 0: at node i, 1: at node j
  int axis = 1;  // local bending axis: 0 = x (torsion), 1 = y, 2 = z
  BilinearHinge law;          // moment-only hinge (ignored when nm is set)
  std::optional<NMHinge> nm;  // axial force - moment interaction hinge
};

struct HingeState {
  double theta = 0;    // hinge rotation (node minus bar end), local axis [rad]
  double theta_p = 0;  // plastic rotation
  double M = 0;        // hinge moment [N m]
  double delta = 0, delta_p = 0;  // hinge elongation (conjugate to N) and its plastic part [m] (N-M hinges)
  double N = 0;        // axial force at the hinge, tension positive [N] (N-M hinges)
  bool plastic = false;  // a facet is active
};

struct HingePushoverPoint {
  double lambda = 0;
  std::vector<Vec6> node_disp;      // per model node, global
  std::vector<HingeState> hinges;   // per HingeSpec
};

struct HingePushoverResult {
  std::vector<HingePushoverPoint> points;  // one per requested lambda
  std::vector<double> events;              // lambda of every state change, in order
  std::vector<double> permanent_events;    // lambda_perm of the changes under the permanent case
  HingePushoverPoint permanent;            // state under the full permanent load (lambda = 0)
};

struct HingePushoverOptions {
  std::optional<Id> permanent_case;  // applied fully before the controlled case
  // SECOND ORDER (Theory II, EXACT beam-column): every frame bar with the
  // stability-function stiffness k_exact(N) and fixed-end forces f0(N) of its
  // axial force N (beam_column_f0), N from the current state (fixed point).
  // Inside a segment of constant active set the hinge law is linear in total
  // form, s = s0 + Kt (e - e0), so each load level is a nonlinear equation in
  // u alone; the events — a facet reached, an active multiplier at its
  // maximum (unloading) — are located by bisection / golden section on the
  // load factor to round-off, the path being probed every 1/substeps of each
  // requested interval. A limit point (no Theory-II equilibrium, an
  // indefinite K(N)) or a mechanism is a typed error.
  bool second_order = false;
  int substeps = 20;
};

// Exact pushover of `model` under the nodal loads of case `case_id` scaled by
// each lambda of `lambdas` (non-decreasing, >= 0), monotonic loading from
// the unstressed state (or from the permanent-load state).
HingePushoverResult hinge_pushover(const Model& model, Id case_id, const std::vector<HingeSpec>& hinges,
                                   const std::vector<double>& lambdas, const HingePushoverOptions& opt = {});

}  // namespace keta

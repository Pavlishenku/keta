// Design moments of reinforced-concrete slabs by Wood-Armer, their assessment
// against given capacities, and the Clark-Nielsen steel forces of a shell.
//
// --- The one kernel: two reinforcement directions resisting a 2x2 tensor ---
// Bars in direction x (unit e1) and in direction y' at the angle alpha from x
// (unit a = (cos alpha, sin alpha), counter-clockwise, 0 < alpha < pi; alpha =
// pi/2 is the orthogonal mesh). A layer of bars of capacity X per unit width
// (measured perpendicular to the bars) resists, across a line of normal n, the
// normal action X (e1.n)^2 (Johansen's stepped yield criterion). A demand
// tensor S (a moment field, or the membrane forces of a sandwich layer) is
// safe on EVERY line when
//      X e1 e1' + Y a a' - S  is positive semi-definite,  X, Y >= 0.
// Wood (1968) chooses the (X, Y) that minimises the steel X + Y. With
// T = [e1 a] (invertible for 0 < alpha < pi) the condition is congruent to
//      diag(X, Y) - S'  >= 0,   S' = T^-1 S T^-T,
//      s'11 = Sxx - 2 Sxy cot(alpha) + Syy cot(alpha)^2,
//      s'12 = (Sxy - Syy cot(alpha)) / sin(alpha),
//      s'22 = Syy / sin(alpha)^2
// (Armer's skew formulas), i.e. the orthogonal problem on S'. Its EXACT
// minimiser is closed-form: with A = X - s'11, B = Y - s'22, min A + B over
// A B >= s'12^2, A, B >= 0 is A = B = |s'12|; when that makes X < 0 the bound
// X = 0 is active and Y = s'22 + s'12^2 / |s'11| (A + s'12^2 / A grows past
// |s'12|), symmetrically for Y; when both stay negative S' is negative
// semi-definite and no steel is needed. These are Wood's rules, now proven
// optimal for any angle.
//
// --- Wood-Armer moments (keta's shell convention, keta/shell.h) ---
// mx, my, mxy per unit width, m > 0 stretches the +z ("top") face. The top
// steel resists S = M, the bottom steel S = -M. Reported signed:
// mx_top, my_top >= 0, mx_bot, my_bot <= 0; my is the moment
// of the y' bars. A result is resolved along the reinforcement first: the
// facet's local x, or the angle `beta` from it (rotate2).
// Combination order: Wood-Armer is NONLINEAR — combine
// (mx, my, mxy) first, then design; the envelope over combinations is then
// the max of the top values and the min of the bottom ones
// (wood_armer_envelope).
//
// --- Assessment (Denton & Burgoyne 1996; Clark, BS 5400 Appendix A) ---
// Given the capacities of the four layers, the utilisation of a face is the
// smallest lambda >= 0 with lambda C - S >= 0 (C = X e1 e1' + Y a a'), i.e.
// the largest ratio demand / capacity over every yield-line direction — the
// 1-D search on the angle done EXACTLY: the feasible set of lambda is
// [lambda*, inf) (C >= 0), its end is 0 or a root of one of the three
// polynomial conditions of a 2x2 PSD matrix, and those candidates are tested
// in increasing order. infinity when a direction loaded in tension has no
// capacity. theta = the critical yield-line normal (angle from x).
//
// --- Clark-Nielsen (sandwich of Morley + Nielsen membrane design) ---
// The six resultants are carried by two layers at the steel levels z_top > 0
// and z_bot < 0 from the mid-surface (equilibrium of force and moment:
// n_top = (m - n z_bot) / (z_top - z_bot), n_bot = (n z_top - m) /
// (z_top - z_bot), componentwise), each layer designed by the kernel on its
// membrane tensor (tension positive): steel forces nx, ny' >= 0 per unit
// width, and the concrete force fc = smallest principal value of the
// remaining concrete tensor S - X e1 e1' - Y a a' (<= 0; uniaxial at the
// optimum). Example: 122.3 / -122.3 kN/m, then 74.74 /
// -169.86 kN/m under 95.11 kN/m of compression.
#pragma once
#include <vector>

#include "keta/statics.h"

namespace keta {

// Rotate a symmetric 2x2 tensor [[a, c], [c, b]] to axes turned by beta
// (counter-clockwise): the components along x' = (cos beta, sin beta).
void rotate2(double& a, double& b, double& c, double beta);

// The kernel: minimal steel (x, y) >= 0 of bars along x and along alpha
// resisting the tensor [[sxx, sxy], [sxy, syy]] on every line.
struct TwoWaySteel {
  double x = 0, y = 0;
};
TwoWaySteel two_way_design(double sxx, double syy, double sxy, double alpha);

struct WoodArmerMoments {
  double mx_top = 0, my_top = 0;  // >= 0, top (+z) bars
  double mx_bot = 0, my_bot = 0;  // <= 0, bottom bars
};
WoodArmerMoments wood_armer(double mx, double my, double mxy, double alpha = kPi / 2);

// Field over a result's facets (model.shells order); beta = angle of the x
// bars from each facet's local x.
std::vector<WoodArmerMoments> wood_armer_field(const CaseResult& r, double alpha = kPi / 2, double beta = 0);
// Envelope of the design moments over several COMBINED results.
std::vector<WoodArmerMoments> wood_armer_envelope(const std::vector<const CaseResult*>& results,
                                                  double alpha = kPi / 2, double beta = 0);

// Capacities of the four layers, all >= 0 (moments per unit width).
struct SlabCapacity {
  double mx_top = 0, my_top = 0, mx_bot = 0, my_bot = 0;
};
struct WoodArmerAssessment {
  double u_top = 0, u_bot = 0, u = 0;  // utilisation (inf when unbounded)
  double theta_top = 0, theta_bot = 0; // critical yield-line normal [rad] (0 when u = 0)
};
WoodArmerAssessment wood_armer_assess(double mx, double my, double mxy, const SlabCapacity& cap,
                                      double alpha = kPi / 2);
// The utilisation of one face: min lambda >= 0, lambda C - S >= 0; theta the
// critical normal.
double face_utilisation(double sxx, double syy, double sxy, double cx, double cy, double alpha,
                        double* theta = nullptr);

struct ClarkNielsenForces {
  double nx_top = 0, ny_top = 0, nx_bot = 0, ny_bot = 0;  // steel forces >= 0 [N/m]
  double fc_top = 0, fc_bot = 0;                          // concrete principal forces <= 0 [N/m]
};
ClarkNielsenForces clark_nielsen(double nx, double ny, double nxy, double mx, double my, double mxy,
                                 double z_top, double z_bot, double alpha = kPi / 2);

}  // namespace keta

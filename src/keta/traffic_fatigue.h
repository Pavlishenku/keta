// Road fatigue load models and effect histories.
//
//  * Rainflow counting by the ASTM E1049-85 §5.4.4 THREE-POINT method on the
//    turning points of the series: X the newest range, Y the previous one;
//    while X >= Y, Y is a cycle — a half cycle when Y holds the start of the
//    history, else a full cycle whose two points leave the stack; the ranges
//    left at the end are half cycles. O(n): the start point leaves the stack
//    by an index shift, the two interior points by a top rewrite (erasing
//    from the front of a vector would be O(n^2) in the worst case).
//  * Effect history of a convoy crossing an influence line: E(r) = sum_j P_j
//    eta(r + a_j) + segment integrals, r the reference position. A uniform
//    sampling can step over a peak; `convoy_history_exact`
//    samples at every BREAKPOINT r = s_i - a_j (an axle on an influence-line
//    node, a segment end on a node): for a piecewise-linear line and point
//    axles E is linear between breakpoints, so the turning points — hence the
//    rainflow ranges — are exact; a cubic line or uniform segments add
//    `per_interval` interior samples per breakpoint interval.
//  * EN 1991-2 §4.6 fatigue load models (loads in N, offsets in m, the
//    keta ConvoyPattern of influence.h): FLM1 = 0.7 Q_ik tandem + 0.3 q_ik
//    (§4.6.2), FLM2 the five frequent lorries (Table 4.6), FLM3 the 4 x 120 kN
//    vehicle and its 36 kN follower 40 m behind (§4.6.4), FLM4 the five
//    equivalent lorries and their share per traffic type (Table 4.7).
#pragma once
#include <vector>

#include "keta/influence.h"

namespace keta {

// ---------------- rainflow ----------------

// Local peaks and valleys: plateaus collapsed, the first and last points kept.
std::vector<double> turning_points(const std::vector<double>& series);

struct RainflowCycle {
  double range = 0;  // |peak - valley|
  double mean = 0;   // (peak + valley) / 2
  double count = 0;  // 1 (closed cycle) or 0.5 (half cycle)
};
// ASTM E1049-85 §5.4.4, in the order the cycles are found.
std::vector<RainflowCycle> rainflow(const std::vector<double>& series);

// ---------------- effect histories ----------------

struct EffectHistory {
  std::vector<double> r;       // reference positions of the convoy, increasing
  std::vector<double> effect;  // E(r)
  double max() const;
  double min() const;
};
// E at reference position r (axles off the line contribute nothing).
double convoy_effect_at(const InfluenceLine& il, const IlPrefix& prefix, const ConvoyPattern& convoy, double r);
// Uniform sweep over every position where the convoy touches the line
// (steps + 1 samples, the convoy_extreme bounds).
EffectHistory convoy_history(const InfluenceLine& il, const ConvoyPattern& convoy, int steps = 500);
// The breakpoint sweep (see the header note); per_interval < 0: 0 for a
// linear line with point axles only, 8 otherwise. O(N_s N_a log(N_s N_a)
// + samples (N_a log N_s)).
EffectHistory convoy_history_exact(const InfluenceLine& il, const ConvoyPattern& convoy, int per_interval = -1);

// ---------------- EN 1991-2 §4.6 fatigue load models ----------------

// FLM1 (§4.6.2): tandem of lane `lane` (1..3, 0 beyond) at 0.7 Q_ik, axles 1.2 m apart.
ConvoyPattern flm1_tandem(int lane);
// FLM1 UDL of lane `lane` over its width: 0.3 q_ik [N/m] (q_1k = 9, q_ik = 2.5 kN/m^2).
double flm1_udl(int lane, double lane_width);
// FLM2 frequent lorry i = 1..5 (Table 4.6).
ConvoyPattern flm2_lorry(int i);
// FLM3 vehicle (§4.6.4): 4 x 120 kN at 0, 1.2, 7.2, 8.4 m; the double
// variant adds the 4 x 36 kN follower whose first axle is `gap` (>= 40 m)
// behind the first axle of the leader.
ConvoyPattern flm3_vehicle();
ConvoyPattern flm3_double(double gap = 40.0);
// FLM4 equivalent lorry i = 1..5 and its share [%] in the traffic type (Table 4.7).
enum class FatigueTraffic { LongDistance, MediumDistance, Local };
ConvoyPattern flm4_lorry(int i);
double flm4_share(int i, FatigueTraffic t);

}  // namespace keta

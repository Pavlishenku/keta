// Road traffic load models on influence lines.
//  * EN 1991-2: carriageway lane division (§4.2.3), LM1 (TS + UDL with alpha
//    adjustment factors, classes 1/2 + French tonnage classes), LM2 single
//    axle (beta_Q), LM4 crowd.
//  * Fascicule 61 titre II: system A(l) with class coefficients a1/a2,
//    Bc / Bt / Mc120 convoys, dynamic amplification delta.
// The transverse distribution is the caller's (Courbon coefficient or
// transverse IL ordinate as a per-lane weight) — v1 longitudinal engine.
#pragma once
#include "keta/influence.h"

namespace keta {

struct LaneDivision {
  int n_lanes = 0;
  double lane_width = 0;     // w_l [m]
  double remaining = 0;      // residual width [m]
};
// EN 1991-2 §4.2.3: w < 5.4 -> 1 lane of 3 m; 5.4 <= w < 6 -> 2 lanes of w/2;
// else n = floor(w/3) lanes of 3 m.
LaneDivision en_lane_division(double carriageway_width);

// LM1 adjustment factors per lane rank (1-based; rank >= 4 -> 0 for TS).
struct Lm1Factors {
  double aQ[3] = {1, 1, 1};  // alpha_Q per lane 1..3
  double aq1 = 1, aq_other = 1;
  static Lm1Factors class1();  // EN recommended
  static Lm1Factors class2();  // French NA 2nd class
};
// Tandem system of lane `rank` (1..3): two axles alpha_Q*Q_k at 1.2 m spacing,
// as a convoy pattern in load units [N].
ConvoyPattern lm1_tandem(int rank, const Lm1Factors& f = {});
// UDL intensity of a lane [N/m]: alpha_q * q_k [N/m^2] * lane width.
double lm1_udl(int rank, double lane_width, const Lm1Factors& f = {});
// LM2 single axle beta_Q * 400 kN (French NA: 320 kN).
ConvoyPattern lm2_axle(double beta_q = 1.0, bool french_na = false);
// LM4 crowd 5 kN/m^2 over a lane width -> line load [N/m].
double lm4_crowd(double lane_width);

// Worst LM1 effect of one lane on an influence line: TS extreme + UDL on the
// sign-matched zones. Returns (max, min).
Extremes lm1_lane_extreme(const InfluenceLine& il, int rank, double lane_width,
                          const Lm1Factors& f = {}, int steps = 500);

// ---------------- Fascicule 61 titre II ----------------

// System A: A(l) = 0.23 + 36/(l+12) t/m^2 (SI: returns N/m^2), l in m.
double fasc61_a_of_l(double loaded_length);
// Class coefficient a1 (class 1..3, loaded lanes 1..4+) and lane width factor
// a2 = v0 / v (v0 = 3.5 / 3.0 / 2.75 m).
double fasc61_a1(int bridge_class, int n_loaded_lanes);
double fasc61_a2(int bridge_class, double lane_width);
// Design lane load [N/m]: max(a1*A(l), 0.4 - 0.0002 l [MPa-scale floor]) * a2 * width.
double fasc61_a_lane_load(double loaded_length, int bridge_class, double lane_width);

// Convoys (axle loads in N, offsets in m, one file).
ConvoyPattern fasc61_bc();     // two 30 t trucks in file
ConvoyPattern fasc61_bt();     // two 16 t tandems
ConvoyPattern fasc61_mc120();  // 110 t tracked vehicle (2 tracks x 6.1 m)
// Dynamic amplification (art. 5.5): delta = 1 + 0.4/(1+0.2L) + 0.6/(1+4G/S).
double fasc61_delta(double element_length, double permanent_load, double max_axle_load);

}  // namespace keta

namespace keta {
// ---------------- EN 1991-2 rail ----------------
// Classification factor alpha scales LM71 and SW/0 only (not SW/2, not the
// unloaded train). Valid values: 0.75..1.46 per §6.3.2.
ConvoyPattern rail_lm71(double alpha = 1.0);      // 4 x 250 kN + 80 kN/m tails (20 m each side)
ConvoyPattern rail_sw0(double alpha = 1.0);       // 133 kN/m over 15 m, gap 5.3 m, twice
ConvoyPattern rail_sw2();                          // 150 kN/m over 25 m, gap 7 m (no alpha)
double rail_unloaded_train();                      // 10 kN/m (sign-matched-zone UDL)
ConvoyPattern rail_hslm_b(int n_axles, double spacing);  // N x 170 kN at spacing d
}  // namespace keta

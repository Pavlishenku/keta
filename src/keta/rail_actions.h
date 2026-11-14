// Railway actions other than the vertical load models. Forces in N.
//
//  * Traction and braking (EN 1991-2 §6.5.3): Q_lak = 33 kN/m x L_ab <= 1000
//    kN; Q_lbk = 20 kN/m x L_ab <= 6000 kN (LM71, SW/0, HSLM) or 35 kN/m x
//    L_ab (SW/2, no cap); times alpha for LM71 and SW/0 (the caller's alpha).
//  * Nosing (§6.5.2): 100 kN at rail level, times alpha only when alpha > 1.
//  * Centrifugal force (§6.5.1): Q_tk = V^2 / (127 r) f Q_vk (V in km/h, r in
//    m), f = 1 - (V - 120)/1000 (814/V + 1.75)(1 - sqrt(2.88 / L_f)) for V >
//    120 km/h and L_f > 2.88 m, else 1; V above 300 km/h counts as 300.
//  * Eccentric vertical loads (§6.3.5): the wheel loads of the two rails in
//    the ratio q2/q1 <= 1.25 — the resultant then sits at e = (q2 - q1)/(q2
//    + q1) r/2, which is r/18 for 1.25.
//  * Distribution through the ballast (§6.3.6): an axle load shared over
//    three sleepers Q/4, Q/2, Q/4 (§6.3.6.1, Figure 6.4, spacing a), each
//    spread 4:1 through the ballast depth h (§6.3.6.2): a footprint of
//    (sleeper width + h/2) along the track by (sleeper length + h/2) across
//    for a monobloc sleeper, or two blocks (block width + h/2) under the
//    rails for a twin-block sleeper.
#pragma once
#include <vector>

#include "keta/influence_surface.h"

namespace keta {

// Dynamic factor of the load models (§6.4.5.2): careful maintenance
// Phi2 = 1.44 / (sqrt(L_phi) - 0.2) + 0.82 in [1, 1.67], standard maintenance
// Phi3 = 2.16 / (sqrt(L_phi) - 0.2) + 0.73 in [1, 2].
double rail_phi2(double L_phi);
double rail_phi3(double L_phi);

enum class RailModel { LM71, SW0, SW2, HSLM };
double rail_traction(double loaded_length, double alpha = 1.0);
double rail_braking(double loaded_length, RailModel model, double alpha = 1.0);
double rail_nosing(double alpha = 1.0);
// f(V, L_f), V in km/h, L_f [m]
double rail_centrifugal_factor(double v_kmh, double L_f);
// Q_tk [N] of a vertical characteristic load q_vertical [N] (or N/m)
double rail_centrifugal_force(double v_kmh, double radius, double q_vertical, double L_f);

struct EccentricRail {
  double q1 = 0, q2 = 0;  // lighter rail at -r/2, heavier at +r/2
  double gauge = 1.5;     // r: distance between the two wheel loads [m]
  double e() const { return (q2 - q1) / (q2 + q1) * gauge / 2; }
};
// a track load split in the ratio q2/q1 = ratio (1 <= ratio <= 1.25)
EccentricRail eccentric_rail(double total, double ratio = 1.25, double gauge = 1.5);

// A uniformly loaded rectangle: centre offset (dx along, dy across) from the
// load's reference, total load, sides ax (along) by ay (across).
struct LoadPatch {
  double dx = 0, dy = 0, load = 0, ax = 0, ay = 0;
};
enum class Sleeper { Monobloc, TwinBlock };
struct BallastOptions {
  double depth = 0.3;            // ballast under the sleeper h [m] (0: slab track)
  Sleeper sleeper = Sleeper::Monobloc;
  double sleeper_width = 0.25;   // along the track [m]
  double sleeper_length = 2.6;   // monobloc, across [m]
  double block_width = 0.6;      // twin block, across, per rail [m]
  double gauge = 1.5;            // between the rails [m]
  double spacing = 0.6;          // sleeper spacing a [m]
  bool three_sleepers = true;    // §6.3.6.1 sharing Q/4, Q/2, Q/4
};
std::vector<LoadPatch> ballast_patches(double axle_load, const BallastOptions& o = {});
// the effect of patches at reference (rx, ry): sum load x mean of eta over each patch
double patches_value(const InfluenceSurface& s, const std::vector<LoadPatch>& p, double rx, double ry);

}  // namespace keta

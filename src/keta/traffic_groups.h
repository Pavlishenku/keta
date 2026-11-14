// Vehicles and EN 1991-2 load groups on influence SURFACES
// (keta/influence_surface.h, keta/traffic_surface.h).
//
// --- A vehicle placed exactly ---
// A vehicle is a set of wheels (dx, dy, P) about its reference point. Its
// effect f(x, y) = sum P_i eta(x + dx_i, y + dy_i) is bilinear between the
// positions where a wheel crosses a gridline: on every rectangle of those
// candidate abscissas x = x_g - dx_i and ordinates y = y_g - dy_i (and the
// bounds) each wheel stays in one cell, so f is bilinear there and its
// extreme sits at a corner. All corners are evaluated: the extreme over the
// band is EXACT (wheels off the grid count zero; both directions of travel
// optionally).
//
// --- Load models ---
// LM2: one axle beta_Q x 400 kN on two wheels 2.00 m apart (EN 1991-2
// 4.3.3; the contact areas taken as points). LM3 (Annex A): the standard
// special vehicles as axle lines — 150 / 200 / 240 kN lines at 1.50 m, the
// double-bogie vehicles with a 12 m gap, the make-weight lines of 1500/200
// and 3000/240 trailing — and the TRANSVERSE layout of an axle line supplied
// by the caller as (dy, fraction) wheels (Annex A figure A.1 and the
// national annexes; not built in: it could not be checked against the text
// here). The special vehicles of national annexes (UK SV, ...) are vehicles
// of the same kind, given by their wheels.
//
// --- Groups (EN 1991-2 table 4.4a), vertical part ---
// gr1a: LM1 characteristic + footways at the combination value 3 kN/m2 on
// their unfavourable parts; gr4: crowd 5 kN/m2 on the unfavourable parts of
// the carriageway and footways; gr5: a special vehicle + LM1 at its frequent
// values (psi1: tandem 0.75, UDL 0.40, EN 1990 A2) kept out of a box around
// the vehicle — its band across, `clear` beyond its first and last axles
// along (the national annex sets it). gr5 is solved EXACTLY by branch and
// bound: vehicle placements in decreasing order of their own effect, each
// completed by the optimal LM1 outside its box, stopped once the vehicle
// effect plus the unrestricted LM1 cannot beat the best total.
// Factors of other national annexes: Lm1Factors (alpha_Q, alpha_q) given by
// the caller.
#pragma once
#include "keta/traffic_surface.h"

namespace keta {

using Wheels = std::vector<std::array<double, 3>>;  // (dx, dy, P [N])

struct VehiclePlacement {
  double value = 0;
  double x = 0, y = 0;   // reference point
  bool reversed = false; // travelling the other way (dx -> -dx)
};
struct VehicleBand {
  double y0 = -1e300, y1 = 1e300;  // every wheel within [y0, y1] across
  bool both_directions = true;
  bool minimise = false;
};
VehiclePlacement vehicle_surface_extreme(const InfluenceSurface& s, const Wheels& w, const VehicleBand& b);
// every corner placement with its value (the search space of the exact
// extreme), in decreasing order of the (signed) value
std::vector<VehiclePlacement> vehicle_surface_placements(const InfluenceSurface& s, const Wheels& w, const VehicleBand& b);

// LM2: beta_Q x 400 kN, two wheels 2.00 m apart (reference: the axle centre)
Wheels lm2_wheels(double beta_q = 1.0);
// Annex A special vehicle "600/150" ... "3600/240/240": axle lines (dx, load
// [N]) from the first line
std::vector<std::pair<double, double>> annex_a_axle_lines(const std::string& name);
// axle lines spread over wheels: layout (dy, fraction of the line), fractions summing to 1
Wheels axle_lines_to_wheels(const std::vector<std::pair<double, double>>& lines,
                            const std::vector<std::pair<double, double>>& layout);

struct Footway {
  double y0 = 0, y1 = 0;  // across
};
struct GroupResult {
  double value = 0;
  Lm1SurfaceResult lm1;            // gr1a, gr5 (frequent)
  double footways = 0;             // gr1a
  VehiclePlacement vehicle;        // gr5
  int vehicle_positions_tried = 0; // gr5
};
GroupResult group_1a(const InfluenceSurface& s, const Lm1SurfaceOptions& lm1, const std::vector<Footway>& footways,
                     double q_footway = 3e3);
GroupResult group_4(const InfluenceSurface& s, double y0, double y1, double q = 5e3, bool minimise = false);
GroupResult group_5(const InfluenceSurface& s, const Wheels& vehicle, const VehicleBand& band, const Lm1SurfaceOptions& lm1,
                    double clear);

// --- Concomitant effects ---
// The effect on `target` (another response's surface on the same grid) of
// the arrangement found for `driver`: the same vehicle, lanes and tandem
// positions, the UDLs on the zones where the DRIVER is unfavourable (outside
// the exclusion boxes) — the companion values of a design extreme (the
// moment with the extreme shear, the reactions with the extreme moment).
double vehicle_value(const InfluenceSurface& target, const Wheels& w, const VehiclePlacement& p);
enum class LoadGroup { Gr1a, Gr4, Gr5 };
struct GroupSpec {
  LoadGroup group = LoadGroup::Gr1a;
  Lm1SurfaceOptions lm1;            // gr1a, gr5 (characteristic factors: gr5 applies psi1)
  std::vector<Footway> footways;    // gr1a
  double q_footway = 3e3;
  double y0 = 0, y1 = 0, q_crowd = 5e3;  // gr4
  Wheels vehicle;                   // gr5
  VehicleBand band;                 // gr5
  double clear = 0;                 // gr5
};
GroupResult group_extreme(const InfluenceSurface& driver, const GroupSpec& g);
double group_concomitant(const InfluenceSurface& target, const InfluenceSurface& driver, const GroupResult& r,
                         const GroupSpec& g);

}  // namespace keta

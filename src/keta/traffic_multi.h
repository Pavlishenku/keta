// Several vehicles on one influence line at once.
//
// Each vehicle (point axles) is placed at a reference position r of a grid
// over the line; its OCCUPANCY interval [r + occ_lo, r + occ_hi] may not
// overlap another vehicle's. A concomitant UDL q loads the sign-matched zones
// of the line except each vehicle's REMOVAL interval [r + rem_lo, r + rem_hi]
// (inside its occupancy, so the removals are disjoint): the objective is then
//   q A_sign + sum_v [E_v(r_v) - q A_sign(removal of v at r_v)],
// separable in the vehicles. A greedy search in the order of entry is only a
// heuristic. Here the optimum over the grid is EXACT, by dynamic programming over (set
// of vehicles placed, rightmost occupancy end): placements sorted by start,
// each one's best predecessor value per subset read from a running maximum
// of the placements already ended — O(P k 2^k) for k vehicles on P
// positions (k <= 12), which makes pruning by influence ratio unnecessary.
// Every vehicle is placed (on the grid of the deck); a deck too short for them
// all is a typed error, not a silent overlap.
#pragma once
#include <optional>
#include <vector>

#include "keta/influence.h"

namespace keta {

struct MobileVehicle {
  std::vector<std::pair<double, double>> axles;  // (offset from the reference [m], load)
  double occ_lo = 0, occ_hi = 0;                 // occupancy, relative to the reference
  double rem_lo = 0, rem_hi = 0;                 // UDL removal, relative (inside the occupancy)
  // occupancy / removal as the axle span extended back (below the least
  // offset) and front (beyond the largest)
  static MobileVehicle from_margins(std::vector<std::pair<double, double>> axles, double occ_back, double occ_front,
                                    double rem_back, double rem_front);
};

struct MultiVehicleResult {
  Extremes extremes;               // max / min of the joint effect, UDL included
  std::vector<double> max_r, min_r;  // reference of each vehicle (input order) at the max / min
};
// positions: steps + 1 references over [s.front(), s.back()].
MultiVehicleResult multi_vehicle_extreme(const InfluenceLine& il, const std::vector<MobileVehicle>& vehicles,
                                         double udl = 0, int steps = 200);

// Sign-split window integrals of a line: the integrals of max(eta, 0) and
// min(eta, 0) over any [a, b] in O(log n), the zero crossings found once
// (exact for a linear line; bisection on the cubic Hermite pieces).
class IlSignPrefix {
 public:
  explicit IlSignPrefix(const InfluenceLine& il);
  explicit IlSignPrefix(const InfluenceLine&&) = delete;
  std::pair<double, double> integrals(double a, double b) const;  // (positive part, negative part)
 private:
  const InfluenceLine& il_;
  std::vector<double> x_;              // the sample abscissae with the zero crossings inserted
  std::vector<double> cpos_, cneg_;    // cumulative parts at x_
  double part(std::size_t k, double a, double b, bool pos) const;  // inside piece k
};

}  // namespace keta

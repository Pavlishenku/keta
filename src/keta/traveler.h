// Form travelers (équipages mobiles).
//
// A traveler (or a launching gantry) is a RIGID body of weight W acting at its
// centre of gravity, bearing on the deck at N >= 3 points. Its reactions are
// those of a rigid body:
//   * N = 3 non-collinear points: the unique statics solution (force and two
//     moment balances);
//   * N > 3: the rigid body on IDENTICAL elastic supports — its displacement
//     is a plane, so the forces are linear in the in-plane position, F_i = a
//     + b x_i + c y_i, which with the three balances is the minimum-norm
//     distribution (centred normal equations, O(N)). Supports of different
//     stiffness or a flexible traveler need the structural model instead.
// A bearing that would pull (F_i < 0) means the traveler tips: reported by
// `traveler_stability`; an anchor (tie-down) may pull. The deck loads are the
// nodal forces -F_i along gravity: placing, moving and removing a traveler
// are ordinary load cases of the phasing (cases stay applied from their
// phase on, so a move is one case: removal at the old position + placement
// at the new one).
#pragma once
#include <vector>

#include "keta/model.h"

namespace keta {

struct TravelerSupport {
  Id node = 0;       // deck node the support bears on
  double x = 0, y = 0;  // position in the plane normal to gravity (the two other global axes, in order)
  bool anchor = false;  // tie-down: may carry uplift
};
struct FormTraveler {
  double weight = 0;          // W > 0 [N]
  double cg_x = 0, cg_y = 0;  // centre of gravity in the same plane
  std::vector<TravelerSupport> supports;
  int gravity = 2;  // global axis the weight acts along, downwards (0 x, 1 y, 2 z)
};
struct TravelerStability {
  std::vector<double> reactions;  // per support, > 0 = the deck pushes up on the traveler
  std::vector<Id> uplifting;      // bearings with a negative reaction: the traveler tips
  bool stable() const { return uplifting.empty(); }
};

// the reactions (throws: fewer than 3 supports, non-finite data, collinear points)
std::vector<double> traveler_reactions(const FormTraveler& t);
TravelerStability traveler_stability(const FormTraveler& t);
// the deck loads: -F_i along gravity at each support node
std::vector<Load> traveler_loads(const FormTraveler& t);
// moving: the removal at `from` and the placement at `to`, one case
std::vector<Load> traveler_move_loads(const FormTraveler& from, const FormTraveler& to);

}  // namespace keta

// Track-structure interaction (UIC 774-3, EN 1991-2 §6.5.4): the model of a
// continuous welded track on embankments and a row of decks, and the
// code's sequenced analysis — temperature with the unloaded ballast, then
// every train position from that state with the loaded resistance under the
// train (keta/link_history.h, a "material change" of the link laws).
//
// --- Model (2-D, in the x-z plane, x along the track, z up) ---
// * track: one bar line (both rails) at z = 0 from x = 0 to the end of the
//   right embankment, elements of `element` length, its far ends held
//   longitudinally (the welded rail continues);
// * decks one after the other, from x = emb_left; each deck: a bar line at
//   its neutral axis z = -e_track, rigid arms to the track level (z = 0) and
//   to its bearings (z = -e_track - e_bearing) at its two ends; each bearing
//   held vertically and longitudinally by a spring (k_left, k_right; 0 =
//   sliding);
// * interaction: at every track node, elasto-plastic longitudinal springs
//   to the ground (embankments) or to the deck top (decks), stiffness and
//   resistance per unit length times the node's tributary length (a node at
//   a deck end shares its two halves between the two supports); vertically
//   rigid (the train's vertical load reaches the deck, whose end rotations
//   move its top and its bearings).
// Out-of-plane freedoms held everywhere.
//
// --- Procedure ---
// thermal stage (rail and deck temperature changes, `steps` increments),
// then for each train position (front from front_first to front_last):
// copy the thermal state, give the springs under the train the loaded law
// (force continuous), add the train (braking or acceleration per unit
// length along +x or -x, vertical load) in `steps` increments. Results:
// the thermal state and every position, and their envelopes — the rail
// stresses (N / A), the bearing forces, the track/deck (and track/ground)
// relative displacement of the interaction springs.
#pragma once
#include <vector>

#include "keta/link_history.h"

namespace keta {

struct TrackDeck {
  double length = 0;
  double A = 0, I = 0, E = 210e9, alpha = 1e-5;
  double e_track = 0;    // neutral axis below the track level
  double e_bearing = 0;  // bearings below the neutral axis
  double k_left = 0, k_right = 0;  // longitudinal bearing stiffness [N/m] (0 = sliding)
  double dT = 0;                   // uniform temperature change [K]
};

struct TrackBridgeInput {
  double element = 2.0;
  double emb_left = 300, emb_right = 300;
  double rail_A = 0, rail_I = 0, rail_E = 210e9, rail_alpha = 1.2e-5, rail_dT = 0;
  // ballast per unit length of track: elastic stiffness [N/m/m], resistance [N/m]
  double k_unloaded = 0, r_unloaded = 0, k_loaded = 0, r_loaded = 0;
  std::vector<TrackDeck> decks;
};

struct TrainPass {
  double length = 300;
  double horizontal = 20e3;  // braking / traction per unit length [N/m], signed along +x
  double vertical = 80e3;    // downward [N/m]
  // the train occupies [front - length, front] when moving toward +x
  // (forward = true), [front, front + length] when moving toward -x
  bool forward = true;
  double front_first = 0, front_last = 0;
  int positions = 16;
};

struct TrackState {
  double rail_min = 0, rail_max = 0;   // extreme rail stresses [Pa] (compression < 0)
  std::vector<double> bearing;         // longitudinal bearing forces, 2 per deck (left, right)
  double rel_deck = 0;                 // max |track - deck| relative displacement
  double rel_ground = 0;               // max |track - embankment| relative displacement
};

struct TrackBridgeResult {
  TrackState thermal;
  std::vector<TrackState> positions;
  TrackState envelope;  // extremes over the positions (bearing: max |force|)
};

TrackBridgeResult track_bridge_uic774(const TrackBridgeInput& in, const TrainPass& train, int steps = 10);

}  // namespace keta

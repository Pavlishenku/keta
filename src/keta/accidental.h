// Accidental rupture during a time history: members (a stay, a
// suspender, a strut) vanish at a given instant of a linear transient.
//
// Two-segment integration: the intact structure from t = 0 to the rupture
// step, then the structure without the severed bars resumed from the state
// reached there. The nodal displacements and velocities are continuous
// through the rupture; the acceleration jumps — the severed members' forces
// leave the internal force at once, so the resumed state is out of
// equilibrium and M a = F - C v - K' u starts the transient (the released
// force applied as a step). The severed members leave with their mass and
// with their own span loads (a cut stay falls with its weight).
//
// By default the history starts at the static equilibrium of the intact
// structure under the loads at t = 0 (the usual cable-loss scenario: the
// bridge at rest under its permanent loads, a stay cut at t_r). Undamped and
// for one degree of freedom the peak of the displacement increment is twice
// the static one; with a damping ratio xi it is 1 + exp(-pi xi / sqrt(1 -
// xi^2)) times it.
#pragma once
#include <functional>

#include "keta/dynamics.h"
#include "keta/nonlinear.h"

namespace keta {

// One load case times a function of time.
struct TimedCase {
  Id case_id = 0;
  std::function<double(double)> factor;  // null = 1
};

struct RuptureHistoryOptions {
  double dt = 0.01;
  int steps = 100;
  int rupture_step = 0;           // the bars vanish at t_r = rupture_step dt (0 <= rupture_step < steps)
  std::vector<Id> bars;           // the severed bars (non-empty)
  DampingSpec damping = rayleigh_none();  // Rayleigh (direct integration)
  IntegrationScheme scheme;
  bool lumped_mass = false;
  bool start_static = true;       // t = 0 at the intact static equilibrium under the loads at t = 0 (else at rest)

  static DampingSpec rayleigh_none() {
    DampingSpec d;
    d.kind = DampingSpec::Kind::Rayleigh;
    return d;
  }
};

struct RuptureHistoryResults {
  double dt = 0;
  int rupture_step = 0;
  // per step (0..steps), per model node, global frame: the fields of the
  // structure standing at that step (nodes only carried by severed bars
  // take zero after the rupture)
  std::vector<std::vector<Vec6>> u, v, a;
  std::vector<Vec6> static_intact;   // static equilibrium of the intact structure under the loads at t = 0
  std::vector<Vec6> static_severed;  // static equilibrium of the severed structure under the loads at t_end
};

RuptureHistoryResults rupture_time_history(const Model& model, const std::vector<TimedCase>& loads,
                                           const RuptureHistoryOptions& opt);

// The same in large displacements (corot_newmark): catenary cables
// may be severed (a stay, a suspender). Segment 1 integrates the intact
// structure from opt.initial_displacement (typically its converged
// solve_corotational state) to t_r = rupture_step dt; segment 2 integrates
// the structure without `bars` from the displacements and velocities reached
// there, under lam(t) continued. The observed histories (opt.observe) and
// the final state are returned as by corot_newmark (steps + 1 rows).
CorotDynResult corot_rupture_history(const Model& model, Id case_id, const std::function<double(double)>& lam,
                                     const CorotDynOptions& opt, int rupture_step, const std::vector<Id>& bars);

}  // namespace keta

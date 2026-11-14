// History-dependent static solve of links with ELASTO-PLASTIC components
// (ballast and track resistance, friction bearings, yielding fuses), loaded
// in stages from the committed state of the previous one, with law changes
// between stages — the "material change" of track-structure interaction
// (UIC 774-3: thermal state with the unloaded ballast, then each train
// position with the loaded resistance under the train).
//
// --- Law ---
// Per link component: elastic stiffness k, yield force fy (symmetric about
// a back force beta), post-yield slope kp (linear kinematic hardening,
// 0 <= kp < k; H = k kp / (k - kp)). State: plastic slip up, back force
// beta; F = k (d - up), elastic while |F - beta| <= fy.
//
// --- Steps (backward Euler, exact) ---
// A stage adds the loads of one case in n equal increments. Over one
// increment, the return mapping of every component is, as a function of
// its end deformation d, the monotone piecewise-linear law through
//   (d-, beta - fy), (d+, beta + fy),  d+- = up + (beta +- fy) / k,
// slopes kp outside — so the increment is the convex problem of
// keta/link_nl.h with these laws, solved exactly (finite Newton), then
// committed (up, beta updated by the closed-form return mapping). The
// model's own non-linear ELASTIC laws (LinkLaw) take part unchanged.
//
// --- Law change (set_law) ---
// The force is kept and the new law continues from it — the end of the
// old curve is the origin of the new one (a "material change" of the
// track-bridge interaction: the thermal slip is locked in when the train
// arrives): up = d - F / k', beta = 0; a force beyond the new yield
// (|F| > fy', e.g. the train leaving) is clamped to it, the unbalance
// being resolved by the next apply().
//
// --- Branching ---
// The object is a value: copy it to run several continuations from one
// committed state (each train position from the same thermal state).
// Loads accumulate stage after stage; result() recovers the cumulative
// state (bar stations, reactions, link forces).
#pragma once
#include <memory>

#include "keta/link_nl.h"

namespace keta {

struct PlasticLink {
  double k = 0;   // elastic stiffness
  double fy = 0;  // yield force (> 0)
  double kp = 0;  // post-yield slope, 0 <= kp < k
};

class LinkHistory {
 public:
  struct Component {
    Id link = 0;
    int comp = 0;
    PlasticLink law;
  };
  LinkHistory(const Model& model, const std::vector<Component>& plastic);
  // change the law of one plastic component (force continuous, see above)
  void set_law(Id link, int comp, const PlasticLink& law);
  // add the loads of `case_id` in `steps` equal increments
  void apply(Id case_id, int steps = 10);
  CaseResult result(int stations = 11) const;
  double force(Id link, int comp) const;
  double slip(Id link, int comp) const;  // plastic slip up
  double deformation(Id link, int comp) const;

 private:
  struct State {
    PlasticLink law;
    double up = 0, beta = 0;
  };
  std::size_t index(Id link, int comp) const;
  LinkLaw step_law(const State& s) const;
  void step(const VecX& rhs, const std::vector<double>& d0);

  std::shared_ptr<const Model> model_;
  std::shared_ptr<const LinkLawCore> core_;
  std::size_t n_elastic_ = 0;  // the model's own laws come first in core order
  std::vector<State> st_;      // plastic components (core order after n_elastic_)
  std::vector<Id> applied_;
  VecX rhs_, v_;
  std::vector<double> d0_;
};

}  // namespace keta

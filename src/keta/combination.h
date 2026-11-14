// EN 1990 combinations solved as EXTREMES, field by field, with the national
// annex's coefficients as data, exclusive actions and the provenance of every
// extreme (a "smart combination" with the contribution of every action).
//
// --- Rules ---
// Each action has a kind (permanent, settlement, variable), its partial
// factors (sup / inf: unfavourable / favourable) and psi0, psi1, psi2. An
// EXCLUSIVE tag groups actions of which at most one acts (settlement
// variants, traffic groups gr1a / gr5...; for permanent actions exactly
// one). The national annex (CombinationRules) gives the permanent and
// settlement factors, xi and the 6.10 / 6.10a+b choice; the actions carry
// their own gamma_Q and psi (no annex is built in: they are data).
// Limit states (EN 1990 6.4.3.2, 6.5.3):
//   ULS 6.10   sum gG G + gP P + gQ1 Q1 + sum gQi psi0i Qi
//   ULS 6.10ab the more adverse of 6.10a (every Q at psi0, no leading one)
//              and 6.10b (xi on the unfavourable permanents)
//   SLS characteristic G + Q1 + sum psi0 Qi; frequent G + psi1 Q1 + sum psi2
//   Qi; quasi-permanent G + sum psi2 Qi  (settlements 1 / 0)
//
// --- Exact extreme of a scalar ---
// For one response, effect e_a of each action at unit factor: for every
// choice of the leading action (and none), each exclusive group contributes
// independently its most adverse admissible member and factor — the
// favourable factor of a permanent action, nothing for a favourable
// variable one; the best leading choice wins. EXACT (the optimum of a sum
// of independent group choices), and the winning terms are returned: their
// sum reproduces the value.
#pragma once
#include "keta/envelope.h"

namespace keta {

enum class ActionKind { Permanent, Settlement, Variable };
struct CombAction {
  Id case_id = 0;
  std::string name;
  ActionKind kind = ActionKind::Variable;
  double gamma_sup = 1.5, gamma_inf = 0.0;  // permanent / settlement: <0 = take the annex's
  double psi0 = 0.7, psi1 = 0.5, psi2 = 0.3;
  std::string exclusive;                   // non-empty: at most one of the tag acts
  bool may_lead = true;
};
struct CombinationRules {
  double gamma_G_sup = 1.35, gamma_G_inf = 1.0;
  double gamma_set_sup = 1.2, gamma_set_inf = 0.0;
  double xi = 0.85;
  bool eq610ab = false;
  // EN 1990 Annex A2, table A2.4(B), recommended values
  static CombinationRules en_recommended() { return {}; }
};
enum class LimitState { ULS, SLS_characteristic, SLS_frequent, SLS_quasi_permanent };
struct CombTerm {
  Id case_id = 0;
  double factor = 0;
};
struct CombinationExtreme {
  double value = 0;
  std::vector<CombTerm> terms;  // the acting actions and their factors (provenance)
  std::string label;            // e.g. "ULS 6.10b, leading: traffic"
};
CombinationExtreme combination_extreme(const std::vector<double>& effects, const std::vector<CombAction>& actions,
                                       const CombinationRules& rules, LimitState ls, bool maximise);

// The same field by field over CaseResults (one per action, same order, same
// shape): the bounds as an EnvelopeBound whose source indices point into
// `combinations` (each distinct winning combination once, with its label).
struct CombinationEnvelope {
  EnvelopeBound lo, hi;
  std::vector<CombinationExtreme> combinations;  // value unused: the terms and label
};
CombinationEnvelope combination_envelope(const std::vector<const CaseResult*>& results, const std::vector<CombAction>& actions,
                                         const CombinationRules& rules, LimitState ls);

}  // namespace keta

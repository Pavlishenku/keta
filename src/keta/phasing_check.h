// Dry-run check of a construction plan (the VERIFIER of staged-construction
// decks): every phase is
// walked on a lightweight simulated state — active bars, supports, links,
// cast parts, tendon states, the date — in the engine's own order
// (keta/phasing.cpp step: supports, bars, links, parts, deactivations,
// support removals, section changes, hinges, modifications, jacking,
// grouting, loads, tensioning),
// WITHOUT solving, and EVERY problem is reported (not the first one).
//
// `Erreur`: the engine would refuse the plan at that point (the same
// condition as its require) or the operation is meaningless (a live tendon
// on a removed bar, a date going backwards). `Avertissement`: legal but
// silently without effect (a span load on an inactive bar, an age given to a
// bar this phase does not cast, grouting a tendon never tensioned unbonded).
// A clean report proves the references resolve and the order is consistent;
// it does not prove the convergence of the solve, nor that a hinge is placed
// where no force is locked in (a run-time state). O(items log items).
#pragma once
#include <string>
#include <vector>

#include "keta/model.h"

namespace keta {

struct PlanIssue {
  enum class Severity { Erreur, Avertissement };
  Severity severity = Severity::Erreur;
  std::size_t phase = 0;  // 0-based index in the plan
  std::string kind;       // short stable code, e.g. "barre_inconnue"
  std::string message;    // French, with the ids involved
};

std::vector<PlanIssue> check_plan(const Model& model, const PhasePlan& plan);

}  // namespace keta

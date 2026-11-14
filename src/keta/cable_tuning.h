// Cable (stay) tuning: the actions — prestress of stays, factors of load
// cases — that bring chosen results to targets (deck profile, pylon
// verticality, forces, reactions): a linear tuning in three modes (exact,
// optimised, best fit) and its nonlinear follow-up.
//
// --- Linear tuning (tune_linear) ---
// Unknowns: the factors x_j of load cases (typically one unit-prestress
// case per stay: add_unit_prestress_case, an imposed shortening P / EA that
// tensions the stay by P where it is free to shorten — a "unit
// tensile load case"). Each target is a result r_i (displacement, reaction,
// bar force at a station, link force): r = b + A x, b from the base cases
// (summed), column j of the influence matrix A from case j — every case of
// one static solve, so the influence is exact (linear superposition).
// Targets are equalities, bounds (at most / at least) or least-squares
// terms (fit, weight w); unknowns may carry bounds and a quadratic cost
// q (x - x_ref)^2. The problem
//     min sum_fit w (a_i x + b_i - t_i)^2 + sum_j q_j (x_j - ref_j)^2
//     s.t. equal targets, bound targets, bounds on x
// is a convex QP solved exactly (keta/qp.h, Goldfarb-Idnani). With as many
// equal targets as unknowns and nothing else it is the square solve
// ("exact"); more unknowns than targets need costs ("optimised");
// more targets than unknowns are fitted ("best fit"). The objective must be
// strictly convex on the unknowns the targets leave free: a typed error
// otherwise.
//
// --- Nonlinear tuning of catenary cables (tune_cables) ---
// Unknowns: the unstrained lengths L0 of Cable bars and the lack of fit
// (shortening) of frame or truss bars (TuningUnknown::Kind: the
// general initial-state problem — a target geometry, forces or reactions ->
// the unstrained lengths that realise it); responses from the exact coupled
// cable statics (keta/cable.h), plus the cable end tensions.
// Sequential QP: linearise the targets by finite differences of the true
// solve, solve the QP for the step, re-solve, until every equal target is
// met to `tol` (relative to its scale: the target, the value, the starting
// value and, for a translation, the bounding-box diagonal of the model —
// a zero target is met to tol x the size of the structure) and, when fitted or bounded targets
// are present, the step itself has shrunk below `tol` (a least-squares fit
// of a target geometry converges like Gauss-Newton). The tuned state is the
// converged nonlinear one; the finite-difference Jacobian only drives the
// iteration. Rank: when the equal targets are at least as many as the
// unknowns and none carries a cost, a singular value ratio of their
// Jacobian below 1e-10 is a typed error that names the dependent targets
// and the undetermined unknowns (never a silent minimum-norm choice).
// The frame is linear (solve_cable_static) or, with large_displacements,
// in large displacements (solve_corotational).
#pragma once
#include <limits>
#include <vector>

#include "keta/cable.h"

namespace keta {

struct TuningTarget {
  enum class Kind { Displacement, Reaction, BarForce, LinkForce, CableTension };
  enum class Mode { Equal, AtMost, AtLeast, Fit };
  Kind kind = Kind::Displacement;
  Id id = 0;     // node, bar or link
  int comp = 0;  // 0..5 (displacement / reaction / link force); bar: N, Vy, Vz, T, My, Mz; cable: end 0 / 1
  int station = 0;  // bar force: station index (-1 = last)
  double value = 0;
  Mode mode = Mode::Equal;
  double weight = 1;  // Fit
};

struct TuningUnknown {
  Id case_id = 0;  // tune_linear: the case whose factor is sought; tune_cables: the Cable bar
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  double cost = 0;  // q in q (x - ref)^2
  double ref = 0;
  // --- appended ---
  // tune_cables: CableLength = the unstrained length L0 of the Cable bar
  // `case_id`; BarShortening = the lack of fit delta of the frame or truss bar
  // `case_id` (fabricated short by delta > 0: the imposed axial strain
  // -delta / L of the bar, exactly as a Prestrain load), starting from 0.
  enum class Kind { CableLength, BarShortening };
  Kind kind = Kind::CableLength;
};

struct TuningResult {
  std::vector<double> x;         // factors (linear) or unstrained lengths (cables)
  std::vector<double> achieved;  // target results in the tuned state
  MatX influence;                // d results / d x at the solution
  std::vector<double> base;      // linear: the results of the base cases
  CaseResult state;              // the tuned state
  int iterations = 0;            // cables: outer iterations
};

// A case tensioning `bar` by P where it can shorten freely (imposed axial
// strain -P / EA). Returns the case id.
Id add_unit_prestress_case(Model& m, Id bar, double P = 1.0, const std::string& name = "");

TuningResult tune_linear(const Model& model, const std::vector<Id>& base_cases,
                         const std::vector<TuningUnknown>& unknowns, const std::vector<TuningTarget>& targets);

struct CableTuningOptions {
  double tol = 1e-9;   // relative tolerance on the equal targets
  int max_iter = 50;
  CableStaticOptions statics;
  // true: the responses come from the corotational solve of the whole
  // structure (keta/nonlinear.h: frame in large displacements, exact
  // catenaries) — displacement, reaction and cable-tension targets
  bool large_displacements = false;
  int substeps = 1;  // of the corotational solve
};
TuningResult tune_cables(const Model& model, Id case_id, const std::vector<TuningUnknown>& cables,
                         const std::vector<TuningTarget>& targets, const CableTuningOptions& opt = {});

// one target's value in a result (cable tensions need the model and the
// unstrained lengths: tune_cables evaluates them)
double target_value(const Model& model, const CaseResult& r, const TuningTarget& t);

}  // namespace keta

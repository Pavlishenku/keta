// Saved tangent state: the structure linearized about one equilibrium state,
// and the eigen analyses run on it.
//
// A TangentState holds
//   * the structure that was linearized: the model itself, or the active
//     substructure of a construction phase (its bars, supports, links and
//     moduli at the date of the phase);
//   * the state: nodal displacements and the axial force of every bar
//     (tension +), plus the shell membrane resultants of a linear state;
//   * K, the tangent stiffness on the free DOFs of
//     assemble_elements(*structure, true) — K + Kg(N, n_shell) about a linear
//     state, the consistent corotational tangent about a large-displacement
//     state (exact catenary pair tangents included), K_frame + Kg(N_frame) +
//     the catenary pair tangents about a cable-static state;
//   * provenance: the source, the case or phase, a fingerprint of the model
//     the state came from, whether it is converged, and the asymmetry of the
//     raw tangent (follower loads give a non-symmetric load stiffness, which
//     the symmetric eigen analyses refuse).
// The consumers:
//   * solve_modal_on_state: (K_T, M) — the modes of the loaded structure
//     (tension stiffens, compression softens, omega_1 -> 0 at the state's own
//     buckling load);
//   * solve_buckling_on_state: K_T + lambda Kg(N_add) singular — the critical
//     multiplier of an ADDITIONAL load on the structure in its locked-in state
//     (a construction-stage check); N_add of the additional load solved
//     linearly on the elastic operator of the same structure.
// The mass is the consistent (or lumped) mass of the reference geometry.
#pragma once
#include <cstdint>
#include <memory>
#include <string>

#include "keta/cable.h"
#include "keta/eigsolve.h"
#include "keta/nonlinear.h"
#include "keta/phasing.h"

namespace keta {

enum class TangentSource { StaticCase, CableState, Corotational, PhasingStage };

struct TangentState {
  TangentSource source = TangentSource::StaticCase;
  Id case_id = 0;                          // the load case of the state (0 for a phase)
  std::string label;                       // the phase name (PhasingStage)
  std::uint64_t model_fingerprint = 0;     // model_fingerprint() of the model the state came from
  std::shared_ptr<const Model> structure;  // what was linearized (the model, or the phase's substructure)
  std::vector<Vec6> u;                     // nodal displacements of the state, per node of `structure`
  std::vector<double> axial;               // N per bar of `structure` (tension +); 0 on cable bars
  std::vector<Vec3> shell_membrane;        // (nx, ny, nxy) per shell of a linear state (empty otherwise)
  SpMat K;                                 // tangent, free DOFs of assemble_elements(*structure, true)
  bool converged = true;
  // |K - K^T|_max / |K|_max of the raw tangent; K itself is stored symmetrized.
  double asymmetry = 0;
};

// FNV-1a over the canonical JSON of the model's structure (load cases and
// combinations excluded): equal for the same structure in one build, used to
// refuse a state handed with another model.
std::uint64_t model_fingerprint(const Model& model);

// About the linear static solution of `case_id`: K + Kg(N, n_shell).
TangentState tangent_of_static_case(const Model& model, Id case_id);
// About a converged solve_cable_static state: K_frame + Kg(N_frame) + the
// catenary pair tangents at the deformed anchor positions. The frame's
// geometric stiffness is included (without it a mast compressed by its stays
// would vibrate as if unloaded).
TangentState tangent_of_cable_state(const Model& model, Id case_id, const CableStaticResults& state);
// About a converged solve_corotational state (load factor 1): the consistent
// tangent the Newton driver factorizes, cables included.
TangentState tangent_of_corotational(const Model& model, Id case_id, const CorotResults& state,
                                     CorotRotationMeasure measure = CorotRotationMeasure::Log);
// About the end of a construction phase: its active substructure (recorded
// with PhasingOptions::keep_structures) at the moduli of its date, K + Kg of
// the cumulated axial field. Catenary cables active in the phase are refused.
TangentState tangent_of_phase(const Model& model, const PhaseResult& phase);

// Modal analysis on the saved tangent. The results' shapes are per node of
// tangent.structure. opt.axial_state must be null (the state is the tangent's).
ModalResults solve_modal_on_state(const Model& model, const TangentState& tangent, ModalOptions opt = {});

// Buckling of an additional load about the saved state: K_T + lambda Kg(N_add)
// singular, smallest positive lambda first. An additional load without any
// axial force is a typed reject, so is a K_T that is not positive definite
// (the state is already unstable on its own).
BucklingResults solve_buckling_on_state(const Model& model, const TangentState& tangent,
                                        const std::vector<Load>& additional, int n_modes = 3,
                                        int dense_threshold = 600);

}  // namespace keta

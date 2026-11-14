// Substructuring (superelements): static condensation of Guyan and the
// fixed-interface modal synthesis of Craig-Bampton, for modal analysis.
//
// A substructure is a set of bars and shells. Its INTERIOR nodes are the
// nodes of its elements that no other element (bar, shell, link, rigid link,
// constraint, skew support) touches and that are not declared masters; the
// free DOFs of every other node form the boundary (the masters and the
// residual structure). Per substructure s, interior I and boundary B:
//   u_I = Psi u_B + Phi q,  Psi = -K_II^-1 K_IB  (constraint modes, exact
//   static condensation), Phi = the first n_modes fixed-interface modes of
//   (K_II, M_II), mass-normalized (Craig-Bampton; n_modes = 0 is Guyan);
// the reduced pencil T^t K T, T^t M T over (u_B, q_1, q_2, ...) is solved
// densely and its modes expanded by T to the whole model. EXACT limits:
// every fixed-interface mode kept = the full model's spectrum; massless
// interiors = the full spectrum by Guyan alone (static condensation is exact
// for statics). Otherwise the reduced frequencies are upper bounds
// (Rayleigh-Ritz on the subspace of T).
#pragma once
#include "keta/eigsolve.h"

namespace keta {

struct Substructure {
  std::vector<Id> bars, shells;  // its elements
  int n_modes = 0;               // fixed-interface modes kept (0 = Guyan)
};

struct SubstructureInfo {
  int n_interior = 0;                       // interior free DOFs condensed
  std::vector<double> fixed_interface_hz;   // the kept fixed-interface frequencies
};

struct SubstructuredModal {
  ModalResults modal;            // modes expanded to the whole model
  int n_boundary = 0;            // boundary free DOFs kept
  int n_reduced = 0;             // n_boundary + the modal coordinates
  std::vector<SubstructureInfo> subs;
  // the superelement matrices: reduced K, M over (u_B in free-DOF order,
  // then each substructure's modal coordinates), and the expansion T
  MatX K, M;
  SpMat T;
  std::vector<int> boundary_free;  // free index of each boundary coordinate
};

// `masters`: nodes kept even when interior to a substructure. Rejected:
// rigid-link / constraint reductions (master-slave models), a bar or shell in
// two substructures, a substructure without interior DOF.
SubstructuredModal solve_modal_substructured(const Model& model, const std::vector<Substructure>& subs,
                                             const std::vector<Id>& masters = {},
                                             const ModalOptions& opt = {});

}  // namespace keta

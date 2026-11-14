// Modal analysis and linear (eigenvalue) buckling.
//
// Both drivers put the guaranteed-SPD reduced stiffness K on the factorized
// side (the flipped pencil):
//   modal:    M phi = mu K phi   (mu largest  -> lowest frequencies, om^2 = 1/mu)
//   buckling: Kg x  = mu K x     (mu most negative -> smallest positive load
//                                 factor lambda = -1/mu, since K + lambda Kg
//                                 goes singular; Kg here is tension-stiffening)
// Solved by a K-inner-product Lanczos with full reorthogonalization (own
// implementation, deterministic start vectors) above a size threshold, dense
// Cholesky-reduced eigensolve below. Every returned pair is certified by its
// pencil residual.
//
// Modal completeness (solve_modal). A single-vector Krylov space holds one
// direction per eigenspace: the second copy of a repeated frequency (square
// or circular sections, rings, symmetric towers, the rigid modes of a free
// body) only appears if rounding happens to feed it, and the Lanczos stopping
// test, which bounds the error of the pairs it HAS, cannot see a pair it
// lacks. solve_modal therefore follows the shifted block Lanczos of Grimes,
// Lewis & Simon (SIAM J. Matrix Anal. Appl. 15(1), 1994):
//   * spectral shift: the pencil is (K - sigma M, M), mu = 1 / (omega^2 - sigma);
//     sigma = 0 when K is positive definite, a small negative shift when it is
//     singular (free bodies: their rigid modes come out at 0 Hz);
//   * block Krylov space (ModalOptions::block_size columns), which captures a
//     multiplicity up to the block size in one pass;
//   * Sturm certificate: the inertia of K - c M (Sylvester) counts the
//     eigenvalues below a cut c just above the highest mode returned; while
//     some are missing, the Lanczos restarts on a fresh block deflated against
//     the modes already found. A count still unmatched after the restarts is
//     reported in ModalResults::warning, never returned silently.
#pragma once
#include <cstdint>
#include <limits>
#include <string>

#include "keta/linsolve.h"
#include "keta/statics.h"

namespace keta {

struct ModalOptions {
  int n_modes = 10;
  bool lumped_mass = false;
  int dense_threshold = 600;  // free-DOF count above which Lanczos is used
  // Optional prestress state: axial force per bar folds Kg into the pencil
  // (vibration about a loaded state; tension raises frequencies). Per-bar by
  // signature: on a shell-bearing model the shells enter K and M but carry no
  // prestress Kg here (mirrors the eig_ext `_fields` convention).
  const std::vector<double>* axial_state = nullptr;
  const SolveHooks* hooks = nullptr;  // cancellation + progress (keta/hooks.h), not owned
  // --- appended ---
  // Participation factors with the consistent-mass coupling to the moving
  // supports, Gamma = phi^T (M_ff r_f + M_fs r_s) (ground_coupling): the exact
  // earthquake load of a consistent mass (Chopra; without it the error is
  // -1.5 % at 10 elements, O(h^2)). ON by default (the exact solution, the
  // rule of keta); false = phi^T M_ff r_f, what ST1 takes.
  // total_mass is then c^T M_ff^-1 c, c = M_ff r_f + M_fs r_s — the sum of
  // Gamma^2 over the whole spectrum, so the effective-mass ratios sum to 1
  // (r_f^T M_ff r_f without the coupling). No effect on a lumped mass.
  bool support_coupling = true;
  // Spectral shift sigma [rad^2/s^2] of the pencil (K - sigma M, M); NaN =
  // automatic: 0 when K is positive definite, else sigma_auto (free bodies).
  // Must keep K - sigma M positive definite (sigma below the lowest omega^2).
  double shift = std::numeric_limits<double>::quiet_NaN();
  // Block size of the Lanczos (columns per Krylov step): repeated frequencies
  // up to this multiplicity are found in one pass, the Sturm certificate and
  // the deflated restarts cover larger ones. 6 = the rigid modes of a free body.
  int block_size = 6;
  // --- appended ---
  // A stiffness that REPLACES K (+ Kg of axial_state) on the free DOFs of
  // assemble_elements(model, true, axial_state): a saved tangent
  // (keta/tangent.h). axial_state then only shapes the zero-stiffness
  // auto-fix (the DOFs a prestress braces). Not owned.
  const SpMat* stiffness = nullptr;
};

struct ModalResults {
  std::vector<double> frequencies_hz;
  std::vector<double> omegas;                 // rad/s
  std::vector<VecX> modes_free;               // free-DOF mode shapes, M-normalized
  std::vector<std::vector<Vec6>> modes;       // per node, unit-crest scaled
  std::vector<std::array<double, 6>> participation;   // Gamma per direction
  std::vector<std::array<double, 6>> effective_mass;  // Gamma^2
  std::array<double, 6> total_mass{};
  std::vector<double> residuals;              // pencil residual per mode
  // --- appended ---
  double shift = 0;        // sigma of the pencil solved (0 = unshifted)
  int sturm_count = -1;    // eigenvalues below the certification cut (-1: not certified, dense path)
  std::string warning;     // non-empty iff fewer modes than asked could be certified (never silent)
};

ModalResults solve_modal(const Model& model, const ModalOptions& opt = {});

// Automatic spectral shift for a singular stiffness: sigma = -1e-8 max_i
// K_ii / M_ii over the dofs carrying mass. K_ii / M_ii is the Rayleigh
// quotient of a unit vector: below the highest eigenvalue omega_max^2 of the
// pencil and within a small factor of it. Two bounds fix the factor:
//   * conditioning: K - sigma M has its smallest eigenvalues ~ |sigma| m along
//     the rigid modes, a condition number ~ omega_max^2 / |sigma| = 1e8; the
//     solves then keep ~1e-8 of accuracy and every pair passes its residual
//     test (1e-12 left 1e-4: the pairs were rejected);
//   * convergence: the shift-invert gap ratios are those of the unshifted
//     pencil while |sigma| < omega_1^2, i.e. omega_max / omega_1 < 1e4, which
//     covers the usual meshes; beyond, the low modes are compressed and
//     converge more slowly (never wrongly: the Sturm certificate stands).
double automatic_shift(const SpMat& K, const SpMat& M);

// Number of eigenvalues omega^2 < c of (K, M), by the inertia of K - c M
// (Sylvester's law); c exactly at an eigenvalue is nudged up by 1e-9 relative.
int eigenvalue_count_below(const SpMat& K, const SpMat& M, double c);

struct BucklingOptions {
  Id case_id = 0;  // 0 = first case
  int n_modes = 3;
  int dense_threshold = 600;
  const SolveHooks* hooks = nullptr;  // cancellation + progress (keta/hooks.h), not owned
};

struct BucklingResults {
  std::vector<double> load_factors;           // ascending positive
  std::vector<std::vector<Vec6>> modes;       // unit-crest
  std::vector<double> residuals;
  // Appended: ST1's modal quantities (see BucklingCteVar in eig_ext.h):
  // generalized geometric stiffness x^T KG x of each unit-crest shape,
  // participation phi^T KG R on the unit vector of every free DOF, total
  // R^T KG R — KG = -Kg(N_ref), positive under compression.
  std::vector<double> gen_stiffness;
  std::vector<double> participation;
  double total_kg = 0;
};

// Reference-case buckling: solve the case linearly, freeze its axial field,
// find the smallest positive multipliers. All-tension -> empty result.
BucklingResults solve_buckling(const Model& model, const BucklingOptions& opt = {});

// Critical multiplier by LDLT inertia bisection on K + lambda*Kg(N_ref):
// bracket-double then bisect on the first negative-pivot flip. Exact by
// Sylvester; cross-checks the eigen route.
double buckling_factor_bisection(const Model& model, Id case_id = 0, double lambda_max = 1e12,
                                 double tol_rel = 1e-9, bool exact_stability = false);

// Number of (K, M) eigenvalues with omega in [f_lo, f_hi] Hz — Sturm count via
// the signed LDLT diagonal of K - (2 pi f)^2 M.
int mode_count_in_band(const Model& model, double f_lo_hz, double f_hi_hz,
                       bool lumped_mass = false);

// Geometric imperfection: clone the model with node coordinates offset by
// scale * the mode's translational components.
Model apply_imperfection(const Model& model, const std::vector<Vec6>& shape, double scale);

// ---- generic kernel (exposed for direct use) ----
// Solve B x = mu K x for the extremal mu (both ends), K SPD prefactored.
// Returns pairs sorted by descending mu. n is the problem size.
struct PencilPair {
  double mu;
  VecX x;
  double residual;
};
// Early stop (appended; null = the full 8 n_wanted + 60 iterations, every
// Ritz pair returned): every `check_every` iterations the Ritz estimates
// beta_j |y_j,i| of the n_wanted LARGEST mu (the modal end) are tested
// against tol |mu_max|; once all pass, the iteration stops (the converged
// pairs agree with the full run to that relative accuracy). With
// `lazy_vectors`, Ritz vectors are formed from the top down only until
// n_wanted pairs have a pencil residual <= accept (and mu above the
// massless floor 1e-12 max(1, mu_max)) — what modal_from_pairs consumes.
struct LanczosStop {
  double tol = 1e-12;
  int check_every = 5;
  bool lazy_vectors = true;
  double accept = 1e-6;
};
std::vector<PencilPair> k_inner_lanczos(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                        int n_wanted, int max_iters = 0, const LanczosStop* stop = nullptr);
// Block version (k_inner_lanczos = block_size 1, seed 0, no locked vectors,
// with the same arithmetic). `locked` (K-orthonormal, with their K-images
// `k_locked`) are deflated: the Krylov space is kept K-orthogonal to them, so
// the pairs returned are new ones. Start block: column 0 of seed 0 is the
// historical 1, 1/2, 1/3, ...; the others are reproducible pseudo-random
// columns (splitmix64 of seed and column). A column that collapses during the
// block orthonormalization (the invariant subspace is exhausted along it)
// leaves the block; the iteration stops when the block is empty. `max_columns`
// bounds the basis size (0 = min(n, 8 n_wanted + 60 + block_size)).
struct BlockLanczosOptions {
  int n_wanted = 1;
  int block_size = 1;
  int max_columns = 0;
  std::uint64_t seed = 0;
  const LanczosStop* stop = nullptr;
};
std::vector<PencilPair> k_inner_block_lanczos(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                              const BlockLanczosOptions& opt,
                                              const std::vector<VecX>& locked = {},
                                              const std::vector<VecX>& k_locked = {});
std::vector<PencilPair> dense_pencil(const SpMat& K, const SpMat& B);
// Scale-free residual of one pencil pair: |B x - mu K x| / (|B x| + |mu| |K x|).
double pencil_residual(const SpMat& K, const SpMat& B, double mu, const VecX& x);

// Rigid-body influence vector {R} for direction d over the free DOFs
// (rotations are taken about the global origin).
VecX rigid_influence_vector(const Model& model, const DofMap& dofs, int d);
// Support coupling of the earthquake load for a rigid ground motion in
// direction d: M_fs r_s, the consistent mass between the free dofs and the
// support-blocked dofs that move with the ground (Chopra: the effective load
// is p_eff = -(M_ff r_f + M_fs r_s) u_g). Zero for a lumped mass (M_fs = 0)
// and for a model without blocked dofs.
VecX ground_coupling(const Model& model, const Assembled& asm_, bool lumped, int d);
// Expand a free-DOF vector to per-node Vec6, scaled so its largest-magnitude
// component is exactly +1 (signed).
std::vector<Vec6> expand_unit_crest(const DofMap& dofs, const VecX& x);
// Modal post-processing shared by every modal driver: take the pencil pairs of
// M phi = mu K phi (descending mu), keep the first n_modes converged ones,
// M-normalize, and fill frequencies / shapes / participation / effective mass.
// `lumped`: the mass flavour of M, for the support coupling of the
// participation factors (ground_inertia).
// `shift`: the sigma of the pencil the pairs solve (omega^2 = 1/mu + sigma; an
// omega^2 within 1e-6 |sigma| of zero is a rigid mode, 0 Hz). A shortfall
// (fewer converged pairs than n_modes while the pencil has more) is written
// to ModalResults::warning.
ModalResults modal_from_pairs(const Model& model, const Assembled& asm_, const SpMat& M,
                              const std::vector<PencilPair>& pairs, int n_modes,
                              bool lumped = false, double shift = 0.0);

}  // namespace keta

namespace keta::ltb {
// Lateral-torsional buckling, closed forms (EN 1993-1-1 background):
// doubly-symmetric section under uniform moment, fork supports:
//   M_cr = (pi/L) sqrt(E Iz G J) sqrt(1 + pi^2 E Iw / (L^2 G J))
// and the monosymmetric root with the Wagner asymmetry beta_x. The
// moment-gradient factor C1 is applied multiplicatively by the caller.
double critical_uniform_moment(double L, double E, double Iz, double G, double J, double Iw);
double critical_uniform_moment_monosym(double L, double E, double Iz, double G, double J,
                                       double Iw, double beta_x);
}  // namespace keta::ltb

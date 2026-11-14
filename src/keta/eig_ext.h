// ST1-flavoured eigen drivers on top of eigsolve.h:
//   * CTE/VAR linear buckling      (manual p. 211-212, "EXEC FLAMBEMENT"),
//   * METHODE 1 = inverse iteration with deflation (p. 194, p. 211),
//   * modal stopping criteria FREQUENCE / RATIO ... AVANCE (p. 194-195).
//
// ---------------- CTE/VAR buckling ----------------
// ST1 splits the axial state in two parts: a FIXED one (CTE <case>, typically
// the permanent loads) kept as it is, and a VARIABLE one (VAR <case>, typically
// a live-load case) amplified until the tangent goes singular. The manual
// (p. 212) writes
//     [K - KG(Ni)] x = c [KG(Nj)] x
// with ST1's compression-positive KG. keta's convention is tangent = K + Kg(N)
// with N tension-positive (assembly.h), so the same statement reads
//     ( K + Kg(N_cte) ) + lambda Kg(N_var)   singular,
// and on the flipped pencil that keeps the guaranteed-SPD side factorized
// (eigsolve.h doctrine)
//     Kg(N_var) x = mu K' x,   K' = K + Kg(N_cte),   lambda = -1/mu,
// so the smallest positive lambda is the most negative mu — exactly the mapping
// solve_buckling() uses, with K' in place of K.
// K' must be positive definite: a CTE state that already buckles the structure
// is a typed failure, verified on the LDLT inertia (Sylvester's law), never a
// silently meaningless eigenvalue. CTE case id 0 means "no fixed axial state"
// (the manual's « il suffit d'appliquer une charge 0 pour la composante fixe »),
// i.e. K' = K, which reproduces solve_buckling() on the VAR case alone.
//
// Tension handling — the manual, p. 212, verbatim:
//   « Il est possible de prendre en compte l'effet du raidissage par traction
//     des cas de charges appliqués (Ni, Nj) (attention car il peut rigidifier
//     excessivement la structure dans les parties tendues et induire des
//     problèmes de convergence numérique). Il suffit de modifier le code
//     correspondant dans le fichier paramètre st1.par. Par défaut l'effort de
//     traction du aux charges appliquées n'est pas pris en compte. »
// keta mirrors that DEFAULT: `tension_stiffening = false` clamps BOTH axial
// fields to N <= 0, so bars in tension contribute nothing to Kg (neither to the
// fixed part nor to the amplified one). Set it to true for the full pencil —
// what solve_buckling() does — knowing it can stiffen the tensioned parts a lot.
// The other documented ST1 fallback (« Si on met une charge 0 sur la composante
// variable, ST1 considère un effort normal uniforme de compression des barres »)
// is opt-in here through `uniform_var_if_zero`; without it a zero variable
// axial field is a typed error, since lambda would otherwise be meaningless.
//
// ---------------- METHODE 1 / METHODE 2 ----------------
// « La méthode 1 est la méthode par itération inverse. Cette méthode calcule
//   mode par mode. La méthode 2 est la méthode par itération du sous espace des
//   vecteurs propres [...] Elle est définie par défaut. » (p. 194)
// METHODE 2 here is keta's existing subspace machinery (K-inner-product Lanczos
// above the dense threshold, Cholesky-reduced dense eigensolve below).
// METHODE 1 is a deflated inverse iteration x <- K^-1 B x, mode by mode, with
// Rayleigh-quotient convergence — an independent cross-check of METHODE 2.
//
// ---------------- modal stopping criteria ----------------
// « EXEC MODE (METHODE <1,2>) (<FREQUENCE fq, RATIO rt>)(AVANCE nm) » (p. 194):
// fq = cutoff frequency (all modes below it), rt = cumulative participation
// target in %, nm = mode-count increment (default 10). Participation follows
// the manual's definitions (p. 195): gamma_j = phi_j^t M R with M-normalized
// modes, modal mass M_j = gamma_j^2, total M_tot = R^t M R, ratio alpha_j =
// M_j / M_tot, cumulated over the modes computed so far.
#pragma once
#include "keta/eigsolve.h"

namespace keta {

// ST1 METHODE selector (values are the manual's).
enum class EigMethod { InverseIteration = 1, Subspace = 2 };

// ---------------- generic kernel ----------------

// Deflated inverse iteration on the flipped pencil B x = mu K x (K SPD, already
// factorized in `Kfac`): power iteration on K^-1 B, so it converges to the
// eigenvalues of LARGEST |mu| first — the lowest frequencies for B = M, the
// smallest |lambda| for B = Kg. Deflation is a two-pass Gram-Schmidt against
// the converged vectors in the K inner product (eigenvectors of a symmetric
// pencil with SPD K are K-orthogonal *and* B-orthogonal; the K product is the
// numerically safe one — M can be singular with lumped/massless DOFs and Kg is
// indefinite by nature).
// Convergence: relative change of the Rayleigh quotient mu = x^t B x / x^t K x
// below `tol`. Returns the converged pairs by descending |mu|; it stops early
// (fewer pairs than asked) when an iteration stalls, which happens for pairs of
// eigenvalues with equal |mu| and opposite signs — every returned pair is
// certified by its pencil residual, none is returned unconverged.
std::vector<PencilPair> inverse_iteration(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                          int n_wanted, int max_iters = 2000,
                                          double tol = 1e-12);

// ---------------- CTE/VAR buckling ----------------

struct BucklingCteVarOptions {
  int dense_threshold = 600;                    // free-DOF count above which Lanczos is used
  EigMethod method = EigMethod::Subspace;       // ST1 METHODE (2 = default)
  bool tension_stiffening = false;              // st1.par default: tension NOT accounted
  bool uniform_var_if_zero = false;             // zero VAR field -> uniform unit compression
  int max_iters = 2000;                         // inverse-iteration cap (METHODE 1)
};

struct BucklingCteVar {
  std::vector<double> load_factors;       // ascending positive lambda (VAR multipliers)
  std::vector<std::vector<Vec6>> modes;   // unit-crest mode shapes
  std::vector<double> residuals;          // pencil residual per mode
  std::vector<double> N_cte, N_var;       // frozen axial fields actually used [N], tension +
  Inertia cte_inertia;                    // inertia of K' = K + Kg(N_cte) (n_neg == 0 enforced)
  // Appended: ST1's modal quantities of a critical mode (manual p.200 with
  // the mass replaced by the geometric stiffness KG = -Kg(N_var), p.215):
  // generalized geometric stiffness x^T KG x of each unit-crest shape (ST1
  // prints the modes normalized to 1: divide by the square root), the
  // participation factor gamma = phi^T KG R of the normalized mode on the
  // unit vector R of every free DOF, and the total R^T KG R
  // (rigidite modale = gamma^2, taux = gamma^2 / total).
  std::vector<double> gen_stiffness;
  std::vector<double> participation;
  double total_kg = 0;
};

// Solve both cases, freeze their axial fields and amplify the variable one.
// `case_cte_id` == 0 means "no fixed axial state" (K' = K); `case_var_id` must
// name a real load case. Throws when the CTE state is already beyond critical.
// Facet shells participate fully: each resolved case also freezes the per-shell
// membrane resultants, which enter Kg through the plate initial-stress term
// (the st1.par tension policy stays a BAR-axial concept — a shell's membrane
// TENSOR enters whole, its principal directions not being separable by a
// scalar sign).
BucklingCteVar solve_buckling_cte_var(const Model& model, Id case_cte_id, Id case_var_id,
                                      int n_modes = 2, const BucklingCteVarOptions& opt = {});

// Same, on axial fields the caller froze itself (combination / phase / envelope
// states — ST1's CTE|VAR <COMB,ETAT,PHASAGE>). Fields are per bar, tension-positive,
// in model bar order; an empty field means "all zero". The per-bar signature
// cannot carry a shell membrane state: on a shell-bearing model the shells
// enter with stiffness + mass only (no shell destabilising term) — use the
// case-id entry points when the shells themselves are prestressed.
BucklingCteVar solve_buckling_cte_var_fields(const Model& model, const std::vector<double>& N_cte,
                                             const std::vector<double>& N_var, int n_modes = 2,
                                             const BucklingCteVarOptions& opt = {});

struct BucklingCteVarBisectOptions {
  double lambda_max = 1e12;
  double tol_rel = 1e-9;
  // Exact stability-function member tangents at N = N_cte + lambda N_var (one
  // element per member is then enough for the exact critical factor).
  bool exact_stability = false;
  bool tension_stiffening = false;
  bool uniform_var_if_zero = false;
};

// Critical VAR multiplier by LDLT inertia bisection on K + Kg(N_cte) + lambda Kg(N_var)
// (or on the exact stability-function tangent assembled at N_cte + lambda N_var).
// Bracket-double then bisect on the first negative-pivot flip — exact by Sylvester,
// an independent cross-check of the eigen route.
double buckling_cte_var_bisection(const Model& model, Id case_cte_id, Id case_var_id,
                                  const BucklingCteVarBisectOptions& opt = {});
double buckling_cte_var_bisection_fields(const Model& model, const std::vector<double>& N_cte,
                                         const std::vector<double>& N_var,
                                         const BucklingCteVarBisectOptions& opt = {});

// ---------------- modal drivers ----------------

struct ModalExtOptions : ModalOptions {
  EigMethod method = EigMethod::Subspace;
  int max_iters = 2000;  // inverse-iteration cap (METHODE 1)
};

// EXEC MODE (nb) (METHODE <1,2>): same result struct as solve_modal(), which is
// exactly what METHODE 2 delegates to.
ModalResults solve_modal_method(const Model& model, const ModalExtOptions& opt = {});

struct ModalFrequencyResults {
  ModalResults modal;
  double cutoff_hz = 0;
  int expected_count = 0;  // Sturm count of eigenvalues in [0, cutoff]
  bool complete = true;    // expected_count modes were actually converged
  std::string warning;     // non-empty iff !complete (never silent)
};

// EXEC MODE FREQUENCE fq: size the problem with the Sturm count in [0, fq],
// solve exactly that many modes and verify the count. Frequencies above the
// cutoff are dropped (they can only come from a count/solver disagreement) and
// reported in `warning`.
ModalFrequencyResults solve_modal_frequency(const Model& model, double f_cutoff_hz,
                                            const ModalExtOptions& opt = {});

struct ModalRatioResults {
  ModalResults modal;
  std::array<double, 6> reached_pct{};              // cumulated ratio per direction [%]
  std::vector<std::array<double, 6>> cumulative_pct;  // running cumul after each mode [%]
  std::array<bool, 6> has_mass{};                   // direction carries mass (target applies)
  int n_modes = 0;
  bool target_reached = false;
  std::string warning;  // non-empty iff !target_reached (never silent)
};

// EXEC MODE RATIO rt AVANCE nm: grow the mode count by `nm_increment` until every
// TRANSLATIONAL direction that carries mass has cumulated `target_pct` % of its
// total mass, or the free-DOF count is exhausted. Like ST1's AVANCE, each
// increment restarts the eigensolve on the bigger mode count (the mode set is
// nested, so `cumulative_pct` is monotone in the mode index).
ModalRatioResults solve_modal_ratio(const Model& model, double target_pct, int nm_increment = 10,
                                    const ModalExtOptions& opt = {});

}  // namespace keta

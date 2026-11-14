// Sparse symmetric solver with a three-rung factorization reuse ladder:
//   1. exact same pattern AND values  -> reuse the numeric factor (pure back-substitution)
//   2. same pattern, new values       -> reuse the symbolic analysis, refactor numerics
//   3. new pattern                    -> full symbolic + numeric factorization
// The compare is exact (element-wise on the CSC arrays), never a hash.
//
// Backend: Eigen SimplicialLDLT (LDL^t without pivoting), fill-reducing
// ordering KetaOrdering (keta/ordering.h: AMD, or a nested dissection on
// large 3-D graphs when it predicts less work); from kSupernodalMin
// unknowns on, the supernodal multifrontal LDL^t of keta/supernodal.h (the
// same factorization with dense kernels). The signed diagonal D
// gives the matrix inertia by Sylvester's law — used by P-Delta beyond-critical
// detection, buckling bisection and arc-length crossing counts.
#pragma once
#include <memory>

#include "keta/core.h"
#include "keta/ordering.h"
#include "keta/supernodal.h"

namespace keta {

struct Inertia {
  int n_neg = 0, n_zero = 0, n_pos = 0;
};

class FactorCache {
 public:
  // Factor (or reuse) for the given matrix. Throws Error on a singular or
  // non-finite factorization. Returns which rung was used (1/2/3).
  int factorize(const SpMat& K);

  // Back-substitute one or many RHS against the current factor.
  VecX solve(const VecX& rhs) const;
  MatX solve(const MatX& rhs) const;
  void require_spd() const;  // throws if any negative pivot (mechanism/beyond-critical)
  Inertia inertia() const;
  // Signed-pivot extremes (local elimination conditioning; a ratio below
  // ~1e-10 marks an effectively singular system).
  double min_abs_pivot() const;
  double max_abs_pivot() const;
  // Smallest ratio |d_k| / |a_kk| of a pivot to the diagonal term of the
  // (reduced) matrix it eliminates. A mechanism whose pivot cancels only to
  // round-off (aligned hinges) escapes the zero-pivot test; its ratio is of
  // the order of the machine precision, where a sound model, however stiff
  // its parts, keeps ratios many orders of magnitude above it.
  double min_pivot_ratio() const;

  // Exact master-slave reduction (rigid links, keta/assembly.h
  // LinkReduction): with a T (n x n_red) set, factorize(K) factors T^t K T
  // and solve(rhs) returns T (T^t K T)^-1 T^t rhs — the extended-space
  // displacement of the constrained problem. Inertia and pivots are those of
  // the reduced matrix. An empty T (0 x 0) clears it. A changed T invalidates the factor (the same
  // T, compared exactly, keeps it).
  void set_reduction(const SpMat& T);
  bool reduced() const { return has_T_; }
  void invalidate() { valid_ = false; }

  struct Stats {
    int full = 0, symbolic_reused = 0, factor_reused = 0;
  };
  const Stats& stats() const { return stats_; }

 public:
  FactorCache() : ldlt_ptr_(std::make_unique<Eigen::SimplicialLDLT<SpMat, Eigen::Lower, KetaOrdering>>()) {}
  FactorCache(FactorCache&&) = default;
  FactorCache& operator=(FactorCache&&) = default;

 private:
  Eigen::SimplicialLDLT<SpMat, Eigen::Lower, KetaOrdering>& ldlt_() const { return *ldlt_ptr_; }
  std::unique_ptr<Eigen::SimplicialLDLT<SpMat, Eigen::Lower, KetaOrdering>> ldlt_ptr_;
  std::unique_ptr<SupernodalLDLT> sn_ = std::make_unique<SupernodalLDLT>();
  bool use_sn_ = false;  // the current factor is the supernodal one
  VecX D_() const { return use_sn_ ? VecX(sn_->vectorD()) : VecX(ldlt_().vectorD()); }  // (Eigen returns D by value)
  SpMat stored_;  // exact copy of the last factored matrix (pattern + values)
  SpMat T_;
  bool has_T_ = false;
  bool valid_ = false;
  bool analyzed_ = false;
  Stats stats_;
};

}  // namespace keta

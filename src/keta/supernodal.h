// Supernodal multifrontal LDL^t (no pivoting) for large sparse symmetric
// matrices — the backend of FactorCache (keta/linsolve.h) from
// kSupernodalMin unknowns on.
//
// Same factorization as Eigen's SimplicialLDLT (unit lower L, diagonal D,
// no pivoting: the signed D gives the inertia; a zero pivot is a
// mechanism), organised for dense kernels: the fill-reducing order of
// KetaOrdering (AMD or nested dissection, keta/ordering.h), the elimination
// tree and the column counts (up-looking symbolic elimination), the
// FUNDAMENTAL supernodes (a column and its only child with the child's
// count = its count + 1: identical structure below), and a multifrontal
// numeric phase in postorder — each supernode's dense frontal matrix gathers
// its original entries and its children's update matrices (extend-add), its
// pivot block is factored densely (unpivoted LDL^t), the panel below by a
// triangular solve and the update matrix by GEMMs, both on fixed blocks.
// The solve runs supernode by supernode (dense triangular solve +
// GEMV/GEMM). Deterministic (fixed operation order, also in parallel:
// set_solver_threads); results equal the simplicial factorization's to
// round-off.
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta {

constexpr int kSupernodalMin = 10000;

// Threads of the supernodal factorization (process-wide; 1 = sequential,
// the default; 0 = every hardware thread). Independent subtrees of the
// supernodal tree run concurrently, and the panel and update of a large
// front are cut into FIXED blocks of 256 rows / columns that run
// concurrently: the result is bit-identical whatever the count.
void set_solver_threads(int n);
int solver_threads();

class SupernodalLDLT {
 public:
  // symbolic analysis of the pattern (both triangles or the lower one; the
  // lower triangle's values are used)
  void analyze(const SpMat& A);
  // numeric factorization (the pattern of the last analyze); false on a
  // zero pivot or a non-finite value
  bool factorize(const SpMat& A);
  VecX solve(const VecX& b) const;
  MatX solve(const MatX& B) const;
  const VecX& vectorD() const { return Dp_; }  // in the permuted order (inertia, pivot extremes)
  const Eigen::VectorXi& permutation() const { return perm_; }  // perm(k) = original index eliminated k-th
  Eigen::Index rows() const { return n_; }

 private:
  Eigen::Index n_ = 0;
  Eigen::VectorXi perm_, iperm_;  // perm_(k) = original index eliminated k-th; iperm_ its inverse
  std::vector<int> sn_start_;     // supernode s: columns [sn_start_[s], sn_start_[s + 1])
  std::vector<std::vector<int>> rows_;  // per supernode: its rows (permuted), the pivot columns first
  std::vector<int> sn_parent_;    // supernodal elimination tree (-1: root)
  std::vector<MatX> L_;           // per supernode: |rows| x width (unit lower on the pivot block)
  VecX Dp_;
  // the permuted lower pattern of A, column by column (original entry
  // positions, for the numeric gather)
  std::vector<std::vector<std::pair<int, int>>> colmap_;  // per permuted column: (permuted row >= col, original (row, col) position index)
  std::vector<std::pair<int, int>> origpos_;              // original (row, col) per position
};

}  // namespace keta

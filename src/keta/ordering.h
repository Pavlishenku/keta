// Fill-reducing ordering of the sparse LDL^t (keta/linsolve.h).
//
// Eigen's AMD is near-optimal on the plane-like graphs of bridge decks and
// frames, but on genuinely 3-D graphs (space frames, solid meshes) its
// greedy minimum degree leaves far more fill than a nested dissection: a
// 14 x 14 x 14 space frame (18 900 DOFs) factors 4x faster with the
// dissection below. KetaOrdering therefore computes AMD and, from
// kNestedDissectionMin DOFs on, also a nested dissection by level-set
// separators (BFS from a pseudo-peripheral vertex, the smallest level
// within the 35-65 % balance band as separator, recursion down to 64
// vertices, separators eliminated last), predicts the factorization work
// of both exactly (symbolic up-looking elimination: column counts c_j,
// work sum c_j^2) and keeps the cheaper. Deterministic (no randomness, ties
// to AMD); below the threshold the ordering is AMD unchanged, so small
// models are bit-identical to the plain AMD factorization.
#pragma once
#include <Eigen/OrderingMethods>

#include "keta/core.h"

namespace keta {

constexpr int kNestedDissectionMin = 10000;

// elimination order (old indices, k-th eliminated first) of a symmetric
// pattern given in full (both triangles)
Eigen::VectorXi nested_dissection_order(const SpMat& full_pattern);
// the exact symbolic factorization of P A P^t for an elimination order:
// strictly-lower nonzeros of L and the work sum c_j^2
struct SymbolicFill {
  double nnz = 0, work = 0;
};
SymbolicFill symbolic_fill(const SpMat& full_pattern, const Eigen::VectorXi& order);

struct KetaOrdering {
  using PermutationType = Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int>;
  template <typename MatrixType>
  void operator()(const MatrixType& mat, PermutationType& perm) {
    Eigen::AMDOrdering<int>()(mat, perm);
    if (mat.rows() < kNestedDissectionMin) return;
    const SpMat A = mat;
    const Eigen::VectorXi nd = nested_dissection_order(A);
    const Eigen::VectorXi amd = perm.indices();
    if (symbolic_fill(A, nd).work < 0.9 * symbolic_fill(A, amd).work) perm.indices() = nd;
  }
};

}  // namespace keta

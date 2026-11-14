// Dense strictly convex quadratic programme, solved EXACTLY by the dual
// active-set method of Goldfarb & Idnani (Math. Programming 27, 1983):
//
//     min 1/2 x' G x + a' x   s.t.   E x = e,   I x >= d,
//
// G symmetric positive definite. The method starts from the unconstrained
// minimum, makes the equalities active, then adds the most violated
// inequality, stepping in the primal and dual spaces at once and dropping
// an active inequality whose multiplier would turn negative: every
// iterate is dual feasible, the objective increases strictly, and it ends
// after finitely many steps at the KKT point (or proves the constraints
// infeasible: a typed error). The small dense systems of each step are
// recomputed from the active set (Cholesky of G, then M = N' G^-1 N) — sizes
// of tuning problems (tens of unknowns) make that cheap and keep the method
// free of factor updates.
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta {

struct QpResult {
  VecX x;
  VecX mult_eq;        // Lagrange multipliers of the equalities
  VecX mult_in;        // of the inequalities (>= 0; 0 when inactive)
  std::vector<int> active;  // active inequality indices
  double objective = 0;
  int iterations = 0;
};

// E (me x n), I (mi x n): rows are the constraints. Empty matrices allowed.
QpResult solve_qp(const MatX& G, const VecX& a, const MatX& E, const VecX& e, const MatX& I, const VecX& d);

}  // namespace keta

#include "keta/qp.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace keta {

QpResult solve_qp(const MatX& G, const VecX& a, const MatX& E, const VecX& e, const MatX& I, const VecX& d) {
  const Eigen::Index n = G.rows();
  require(n > 0 && G.cols() == n && a.size() == n, "qp: G must be square and match a");
  require((E.rows() == 0 || E.cols() == n) && E.rows() == e.size(), "qp: equality rows / right-hand side mismatch");
  require((I.rows() == 0 || I.cols() == n) && I.rows() == d.size(), "qp: inequality rows / right-hand side mismatch");
  require(G.allFinite() && a.allFinite() && E.allFinite() && e.allFinite() && I.allFinite() && d.allFinite(),
          "qp: non-finite data");
  const Eigen::LLT<MatX> L(G);
  require(L.info() == Eigen::Success, "qp: G must be symmetric positive definite");
  const Eigen::Index me = E.rows(), mi = I.rows();
  VecX x = -L.solve(a);
  std::vector<int> act;  // >= 0: inequality; < 0: equality -(k + 1)
  std::vector<double> u;
  auto normal = [&](int c) -> VecX { return c < 0 ? VecX(E.row(-c - 1).transpose()) : VecX(I.row(c).transpose()); };
  // primal step z = H n and dual step r = N* n for the active set
  auto directions = [&](const VecX& np, VecX& z, VecX& r, VecX& Gn) {
    Gn = L.solve(np);
    const Eigen::Index q = static_cast<Eigen::Index>(act.size());
    if (q == 0) {
      z = Gn;
      r = VecX(0);
      return;
    }
    MatX N(n, q);
    for (Eigen::Index k = 0; k < q; ++k) N.col(k) = normal(act[static_cast<std::size_t>(k)]);
    const MatX GN = L.solve(N);
    const MatX M = N.transpose() * GN;
    r = M.ldlt().solve(N.transpose() * Gn);
    z = Gn - GN * r;
  };
  auto drop = [&](std::size_t k) {
    act.erase(act.begin() + static_cast<std::ptrdiff_t>(k));
    u.erase(u.begin() + static_cast<std::ptrdiff_t>(k));
  };
  auto step_dual = [&](const VecX& r, double t) {
    for (std::size_t k = 0; k < u.size(); ++k) u[k] -= t * r(static_cast<Eigen::Index>(k));
  };
  QpResult out;
  // the equalities first (their multipliers are free)
  for (Eigen::Index k = 0; k < me; ++k) {
    const VecX np = E.row(k).transpose();
    VecX z, r, Gn;
    directions(np, z, r, Gn);
    const double zn = z.dot(np);
    require(z.norm() > 1e-12 * Gn.norm() && zn > 0, "qp: dependent equality constraints");
    const double t = -(np.dot(x) - e(k)) / zn;
    x += t * z;
    step_dual(r, t);
    act.push_back(-static_cast<int>(k) - 1);
    u.push_back(t);
    ++out.iterations;
  }
  const int max_it = 50 * static_cast<int>(me + mi + n) + 100;
  for (int it = 0;; ++it) {
    require(it < max_it, "qp: no termination (degenerate constraints)");
    // the most violated inequality, measured along its normal
    int p = -1;
    double worst = 0;
    for (Eigen::Index j = 0; j < mi; ++j) {
      if (std::find(act.begin(), act.end(), static_cast<int>(j)) != act.end()) continue;
      const double nn = I.row(j).norm();
      if (nn == 0) {
        require(d(j) <= 0, "qp: infeasible constraint 0 >= d");
        continue;
      }
      const double s = (I.row(j).dot(x) - d(j)) / nn;
      const double tol = 1e-11 * (1 + std::fabs(d(j)) / nn + x.cwiseAbs().maxCoeff());
      if (s < -tol && s < worst) {
        worst = s;
        p = static_cast<int>(j);
      }
    }
    if (p < 0) break;
    const VecX np = I.row(p).transpose();
    double up = 0;
    for (int inner = 0;; ++inner) {
      require(inner < max_it, "qp: no termination (degenerate constraints)");
      VecX z, r, Gn;
      directions(np, z, r, Gn);
      // the dual blocking: an active inequality whose multiplier reaches 0
      double t1 = std::numeric_limits<double>::infinity();
      std::size_t kdrop = 0;
      for (std::size_t k = 0; k < act.size(); ++k)
        if (act[k] >= 0 && r(static_cast<Eigen::Index>(k)) > 0) {
          const double ratio = u[k] / r(static_cast<Eigen::Index>(k));
          if (ratio < t1) {
            t1 = ratio;
            kdrop = k;
          }
        }
      ++out.iterations;
      if (z.norm() <= 1e-12 * Gn.norm()) {
        // n_p depends on the active normals: a pure dual step
        require(std::isfinite(t1), "qp: infeasible constraints");
        step_dual(r, t1);
        up += t1;
        drop(kdrop);
        continue;
      }
      const double t2 = -(np.dot(x) - d(p)) / z.dot(np);
      const double t = std::min(t1, t2);
      x += t * z;
      step_dual(r, t);
      up += t;
      if (t2 <= t1) {
        act.push_back(p);
        u.push_back(up);
        break;
      }
      drop(kdrop);
    }
  }
  out.x = x;
  out.mult_eq = VecX::Zero(me);
  out.mult_in = VecX::Zero(mi);
  for (std::size_t k = 0; k < act.size(); ++k) {
    if (act[k] < 0) {
      out.mult_eq(-act[k] - 1) = u[k];
    } else {
      out.mult_in(act[k]) = u[k];
      out.active.push_back(act[k]);
    }
  }
  out.objective = 0.5 * x.dot(G * x) + a.dot(x);
  return out;
}

}  // namespace keta

#include "keta/linsolve.h"

#include <limits>

namespace keta {

namespace {

bool same_pattern(const SpMat& a, const SpMat& b) {
  if (a.rows() != b.rows() || a.cols() != b.cols() || a.nonZeros() != b.nonZeros()) return false;
  const auto n = a.cols();
  for (Eigen::Index c = 0; c <= n; ++c)
    if (a.outerIndexPtr()[c] != b.outerIndexPtr()[c]) return false;
  for (Eigen::Index i = 0; i < a.nonZeros(); ++i)
    if (a.innerIndexPtr()[i] != b.innerIndexPtr()[i]) return false;
  return true;
}

bool same_values(const SpMat& a, const SpMat& b) {
  for (Eigen::Index i = 0; i < a.nonZeros(); ++i)
    if (a.valuePtr()[i] != b.valuePtr()[i]) return false;  // bitwise
  return true;
}

}  // namespace

void FactorCache::set_reduction(const SpMat& T) {
  const bool want = T.rows() > 0 || T.cols() > 0;
  if (want == has_T_) {
    if (!want) return;  // no reduction before, none now: the factor stays
    if (T.rows() == T_.rows() && T.cols() == T_.cols() && T.nonZeros() == T_.nonZeros()) {
      SpMat a = T, b = T_;
      a.makeCompressed();
      b.makeCompressed();
      if (same_pattern(a, b) && same_values(a, b)) return;  // the same T: the factor stays
    }
  }
  has_T_ = want;
  T_ = T;
  valid_ = false;
  analyzed_ = false;
  stored_ = SpMat();
}

int FactorCache::factorize(const SpMat& K) {
  SpMat Kc = has_T_ ? SpMat((T_.transpose() * K * T_).pruned()) : K;
  Kc.makeCompressed();
  if (Kc.rows() == 0) {  // fully clamped model: nothing to factor
    stored_ = std::move(Kc);
    valid_ = true;
    return 1;
  }
  int rung;
  const bool big = Kc.rows() >= kSupernodalMin;
  if (valid_ && same_pattern(stored_, Kc)) {
    if (same_values(stored_, Kc)) {
      ++stats_.factor_reused;
      return 1;
    }
    stored_ = std::move(Kc);
    valid_ = false;
    if (use_sn_) {
      if (!sn_->factorize(stored_)) fail("solver: singular stiffness (zero pivot) — the model contains a mechanism");
    } else {
      ldlt_().factorize(stored_);
    }
    ++stats_.symbolic_reused;
    rung = 2;
  } else {
    stored_ = std::move(Kc);
    valid_ = false;  // drop before the fallible factorization (failure hygiene)
    use_sn_ = big;
    if (use_sn_) {
      sn_->analyze(stored_);
      if (!sn_->factorize(stored_)) fail("solver: singular stiffness (zero pivot) — the model contains a mechanism");
    } else {
      ldlt_().analyzePattern(stored_);
      ldlt_().factorize(stored_);
    }
    analyzed_ = true;
    ++stats_.full;
    rung = 3;
  }
  if (!use_sn_ && ldlt_().info() != Eigen::Success) {
    valid_ = false;
    fail("solver: factorization failed (singular or structurally deficient stiffness)");
  }
  // A zero or non-finite pivot means an exact mechanism.
  const VecX Dv = D_();
  for (Eigen::Index i = 0; i < Dv.size(); ++i)
    if (Dv(i) == 0.0 || !std::isfinite(Dv(i))) {
      valid_ = false;
      fail("solver: singular stiffness (zero pivot) — the model contains a mechanism");
    }
  valid_ = true;
  return rung;
}

VecX FactorCache::solve(const VecX& rhs) const {
  require(valid_, "solver: no valid factorization");
  if (rhs.size() == 0) return rhs;
  if (has_T_) {
    if (T_.cols() == 0) return VecX::Zero(rhs.size());
    VecX xr = use_sn_ ? sn_->solve(VecX(T_.transpose() * rhs)) : VecX(ldlt_().solve(VecX(T_.transpose() * rhs)));
    require(xr.allFinite(), "solver: non-finite solution");
    return T_ * xr;
  }
  VecX x = use_sn_ ? sn_->solve(rhs) : VecX(ldlt_().solve(rhs));
  require(x.allFinite(), "solver: non-finite solution");
  return x;
}

MatX FactorCache::solve(const MatX& rhs) const {
  require(valid_, "solver: no valid factorization");
  if (rhs.rows() == 0) return rhs;
  if (has_T_) {
    if (T_.cols() == 0) return MatX::Zero(rhs.rows(), rhs.cols());
    MatX xr = use_sn_ ? sn_->solve(MatX(T_.transpose() * rhs)) : MatX(ldlt_().solve(MatX(T_.transpose() * rhs)));
    require(xr.allFinite(), "solver: non-finite solution");
    return T_ * xr;
  }
  MatX x = use_sn_ ? sn_->solve(rhs) : MatX(ldlt_().solve(rhs));
  require(x.allFinite(), "solver: non-finite solution");
  return x;
}

Inertia FactorCache::inertia() const {
  require(valid_, "solver: no valid factorization");
  Inertia in;
  const VecX Dv = D_();
  for (Eigen::Index i = 0; i < Dv.size(); ++i) {
    if (Dv(i) < 0)
      ++in.n_neg;
    else if (Dv(i) > 0)
      ++in.n_pos;
    else
      ++in.n_zero;
  }
  return in;
}

double FactorCache::min_abs_pivot() const {
  require(valid_, "solver: no valid factorization");
  const VecX Dv = D_();
  double m = Dv.size() ? std::fabs(Dv(0)) : 0.0;
  for (Eigen::Index i = 1; i < Dv.size(); ++i) m = std::min(m, std::fabs(Dv(i)));
  return m;
}
double FactorCache::min_pivot_ratio() const {
  require(valid_, "solver: no valid factorization");
  const VecX Dv = D_();
  const VecX diag = stored_.diagonal();
  VecX dp(diag.size());
  if (use_sn_) {
    const Eigen::VectorXi& p = sn_->permutation();
    for (Eigen::Index k = 0; k < dp.size(); ++k) dp(k) = diag(p(k));
  } else {
    dp = ldlt_().permutationP() * diag;  // P A P^-1 = L D L^t
  }
  double r = std::numeric_limits<double>::infinity();
  for (Eigen::Index k = 0; k < Dv.size(); ++k)
    if (dp(k) != 0.0) r = std::min(r, std::fabs(Dv(k)) / std::fabs(dp(k)));
  return r;
}
double FactorCache::max_abs_pivot() const {
  require(valid_, "solver: no valid factorization");
  const VecX Dv = D_();
  double m = 0.0;
  for (Eigen::Index i = 0; i < Dv.size(); ++i) m = std::max(m, std::fabs(Dv(i)));
  return m;
}

void FactorCache::require_spd() const {
  Inertia in = inertia();
  require(in.n_neg == 0 && in.n_zero == 0,
          "solver: stiffness is not positive definite (" + std::to_string(in.n_neg) +
              " negative pivots) — mechanism or beyond critical load");
}

}  // namespace keta

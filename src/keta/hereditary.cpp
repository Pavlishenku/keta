#include "keta/hereditary.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {

// (1 - e^-y)/y with the cancellation-free small-y branch.
double lambda_of(double y) {
  if (y < 1e-6) return 1.0 - y / 2.0 + y * y / 6.0;
  return -std::expm1(-y) / y;
}
double one_minus_exp(double y) { return -std::expm1(-y); }

}  // namespace

// ---------------- compliance functions ----------------

ComplianceFn ec2_compliance(const ConcreteTime& c, double E) {
  require(finite(E) && E > 0, "hereditary: ec2_compliance needs a finite positive E");
  ConcreteTime cc = c;
  return [cc, E](double t, double t0) {
    require(t0 > 0 && t >= t0, "hereditary: J(t,t0) needs t >= t0 > 0");
    return 1.0 / E + ec2_phi(t, t0, cc) / E;
  };
}

ComplianceFn matlaw_compliance(const MatLaw& law) {
  law.validate();
  MatLaw m = law;
  return [m](double t, double t0) {
    require(t0 > 0 && t >= t0, "hereditary: J(t,t0) needs t >= t0 > 0");
    double e0 = m.E(t0);
    require(e0 > 0, "hereditary: matlaw modulus is not positive");
    return 1.0 / e0 + m.phi(t, t0) / m.E_ref(t0);
  };
}

// ---------------- NNLS (Lawson & Hanson 1974, algorithm 23.10) ----------------

VecX nnls(const MatX& A, const VecX& b, int max_iter) {
  const int n = static_cast<int>(A.cols());
  require(A.rows() == b.size(), "nnls: shape mismatch");
  if (max_iter <= 0) max_iter = 3 * n + 10;
  VecX x = VecX::Zero(n);
  std::vector<bool> in_p(static_cast<std::size_t>(n), false);
  double scale = A.cwiseAbs().maxCoeff() * b.cwiseAbs().maxCoeff();
  double tol = 1e-12 * std::max(1.0, scale);

  auto ls_on_p = [&](VecX& s) {
    std::vector<int> cols;
    for (int j = 0; j < n; ++j)
      if (in_p[static_cast<std::size_t>(j)]) cols.push_back(j);
    s = VecX::Zero(n);
    if (cols.empty()) return;
    MatX Ap(A.rows(), static_cast<Eigen::Index>(cols.size()));
    for (std::size_t k = 0; k < cols.size(); ++k)
      Ap.col(static_cast<Eigen::Index>(k)) = A.col(cols[k]);
    VecX sp = Ap.colPivHouseholderQr().solve(b);
    for (std::size_t k = 0; k < cols.size(); ++k) s(cols[k]) = sp(static_cast<Eigen::Index>(k));
  };

  for (int outer = 0; outer < max_iter; ++outer) {
    VecX w = A.transpose() * (b - A * x);
    int jmax = -1;
    double wmax = tol;
    for (int j = 0; j < n; ++j)
      if (!in_p[static_cast<std::size_t>(j)] && w(j) > wmax) {
        wmax = w(j);
        jmax = j;
      }
    if (jmax < 0) break;
    in_p[static_cast<std::size_t>(jmax)] = true;
    for (int inner = 0; inner < 3 * n + 10; ++inner) {
      VecX s;
      ls_on_p(s);
      bool any_neg = false;
      for (int j = 0; j < n; ++j)
        if (in_p[static_cast<std::size_t>(j)] && s(j) <= 0) any_neg = true;
      if (!any_neg) {
        x = s;
        break;
      }
      double alpha = 1.0;
      for (int j = 0; j < n; ++j)
        if (in_p[static_cast<std::size_t>(j)] && s(j) <= 0) {
          double den = x(j) - s(j);
          if (den > 0) alpha = std::min(alpha, x(j) / den);
        }
      x += alpha * (s - x);
      for (int j = 0; j < n; ++j)
        if (in_p[static_cast<std::size_t>(j)] && x(j) <= 1e-14 * std::max(1.0, x.maxCoeff())) {
          in_p[static_cast<std::size_t>(j)] = false;
          x(j) = 0;
        }
    }
  }
  for (int j = 0; j < n; ++j) x(j) = std::max(0.0, x(j));
  return x;
}

// ---------------- Kelvin chain ----------------

double KelvinChain::J_delayed(double dt) const {
  require(dt >= 0, "KelvinChain: negative elapsed time");
  double s = 0;
  for (std::size_t m = 0; m < tau.size(); ++m) s += a[m] * one_minus_exp(dt / tau[m]);
  return s;
}
double KelvinChain::J(double dt) const { return a0 + J_delayed(dt); }

KelvinChain fit_kelvin_chain(const ComplianceFn& J, double t0, const KelvinFitOptions& o) {
  require(t0 > 0, "fit_kelvin_chain: loading age must be > 0");
  require(o.n_units >= 1 && o.n_units <= 64, "fit_kelvin_chain: n_units out of range");
  require(o.tau_min > 0 && o.decades > 0, "fit_kelvin_chain: bad retardation grid");
  require(o.fit_min > 0 && o.fit_max > o.fit_min && o.n_samples >= o.n_units,
          "fit_kelvin_chain: bad fit window");

  KelvinChain c;
  c.t0 = t0;
  c.a0 = J(t0, t0);
  require(finite(c.a0) && c.a0 > 0, "fit_kelvin_chain: non-positive instantaneous compliance");
  for (int m = 0; m < o.n_units; ++m)
    c.tau.push_back(o.tau_min * std::pow(10.0, o.decades * m));

  const int ns = o.n_samples;
  MatX A(ns, o.n_units);
  VecX b(ns);
  std::vector<double> dt(static_cast<std::size_t>(ns));
  double lmin = std::log10(o.fit_min), lmax = std::log10(o.fit_max);
  for (int i = 0; i < ns; ++i) {
    double d = std::pow(10.0, lmin + (lmax - lmin) * i / (ns - 1));
    dt[static_cast<std::size_t>(i)] = d;
    double jt = J(t0 + d, t0);
    require(finite(jt) && jt > 0, "fit_kelvin_chain: non-finite compliance sample");
    // relative weighting: divide the row by J so the residual is |Jfit/J - 1|
    double w = 1.0 / jt;
    b(i) = (jt - c.a0) * w;
    for (int m = 0; m < o.n_units; ++m)
      A(i, m) = one_minus_exp(d / c.tau[static_cast<std::size_t>(m)]) * w;
  }
  VecX x = nnls(A, b);
  c.a.assign(static_cast<std::size_t>(o.n_units), 0.0);
  for (int m = 0; m < o.n_units; ++m) c.a[static_cast<std::size_t>(m)] = x(m);

  double err = 0;
  for (int i = 0; i < ns; ++i) {
    double d = dt[static_cast<std::size_t>(i)];
    double jt = J(t0 + d, t0);
    err = std::max(err, std::fabs(c.J(d) / jt - 1.0));
  }
  c.fit_max_rel_error = err;
  return c;
}

// ---------------- exponential algorithm ----------------

KelvinState kelvin_state(const KelvinChain& c, double t_start) {
  KelvinState s;
  s.gamma.assign(c.tau.size(), 0.0);
  s.t = t_start;
  s.sigma = 0;
  return s;
}

KelvinStep kelvin_step(KelvinState& st, const KelvinChain& c, double t_end, double sigma_end) {
  require(st.gamma.size() == c.tau.size(), "kelvin_step: state/chain mismatch");
  require(finite(t_end) && finite(sigma_end), "kelvin_step: non-finite input");
  double dt = t_end - st.t;
  require(dt >= 0, "kelvin_step: time must not go backwards");
  double dsig = sigma_end - st.sigma;
  KelvinStep out;
  double flex_delayed = 0;
  for (std::size_t m = 0; m < c.tau.size(); ++m) {
    double y = dt / c.tau[m];
    double e1 = one_minus_exp(y);      // 1 - e^-y
    double lam = lambda_of(y);         // (1 - e^-y)/y
    double dg_hist = e1 * (st.sigma * c.a[m] - st.gamma[m]);
    double dg_lin = dsig * c.a[m] * (1.0 - lam);
    st.gamma[m] += dg_hist + dg_lin;
    out.d_eps_delayed += dg_hist + dg_lin;
    flex_delayed += c.a[m] * (1.0 - lam);
  }
  out.d_eps_total = dsig * c.a0 + out.d_eps_delayed;
  double flex = c.a0 + flex_delayed;
  out.e_algorithmic = flex > 0 ? 1.0 / flex : 0.0;
  st.sigma = sigma_end;
  st.t = t_end;
  return out;
}

double kelvin_creep(KelvinState& st, const KelvinChain& c, double t_end) {
  return kelvin_step(st, c, t_end, st.sigma).d_eps_delayed;
}

// ---------------- Volterra reference ----------------

std::vector<double> log_time_grid(double t0, double dt_min, double dt_max, int n) {
  require(t0 > 0 && dt_min > 0 && dt_max > dt_min && n >= 2, "log_time_grid: bad input");
  std::vector<double> t;
  t.push_back(t0);
  double lmin = std::log10(dt_min), lmax = std::log10(dt_max);
  for (int i = 0; i < n; ++i)
    t.push_back(t0 + std::pow(10.0, lmin + (lmax - lmin) * i / (n - 1)));
  return t;
}

namespace {
// Trapezoidal kernel: the average compliance seen by the stress increment
// applied over the step (t_{j-1}, t_j], observed at t_i.
double jbar(const ComplianceFn& J, double ti, double tjm1, double tj) {
  return 0.5 * (J(ti, tjm1) + J(ti, tj));
}
void check_grid(const std::vector<double>& t) {
  require(t.size() >= 1 && t[0] > 0, "volterra: grid must start at a positive age");
  for (std::size_t i = 1; i < t.size(); ++i)
    require(t[i] > t[i - 1], "volterra: grid must be strictly increasing");
}
}  // namespace

std::vector<double> volterra_strain(const ComplianceFn& J, const std::vector<double>& t,
                                    const std::vector<double>& sigma) {
  check_grid(t);
  require(t.size() == sigma.size(), "volterra_strain: size mismatch");
  const std::size_t n = t.size();
  std::vector<double> eps(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    double e = J(t[i], t[0]) * sigma[0];
    for (std::size_t j = 1; j <= i; ++j)
      e += jbar(J, t[i], t[j - 1], t[j]) * (sigma[j] - sigma[j - 1]);
    eps[i] = e;
  }
  return eps;
}

std::vector<double> volterra_stress(const ComplianceFn& J, const std::vector<double>& t,
                                    const std::vector<double>& eps) {
  check_grid(t);
  require(t.size() == eps.size(), "volterra_stress: size mismatch");
  const std::size_t n = t.size();
  std::vector<double> sig(n, 0.0);
  sig[0] = eps[0] / J(t[0], t[0]);
  for (std::size_t i = 1; i < n; ++i) {
    double rhs = eps[i] - J(t[i], t[0]) * sig[0];
    for (std::size_t j = 1; j < i; ++j)
      rhs -= jbar(J, t[i], t[j - 1], t[j]) * (sig[j] - sig[j - 1]);
    double dii = jbar(J, t[i], t[i - 1], t[i]);
    require(dii > 0, "volterra_stress: non-positive diagonal compliance");
    sig[i] = sig[i - 1] + rhs / dii;
  }
  return sig;
}

std::vector<double> relaxation_function(const ComplianceFn& J, double t0,
                                        const std::vector<double>& t) {
  check_grid(t);
  require(std::fabs(t[0] - t0) <= 1e-12 * t0, "relaxation_function: grid must start at t0");
  std::vector<double> unit(t.size(), 1.0);
  return volterra_stress(J, t, unit);
}

}  // namespace keta

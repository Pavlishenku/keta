// Hereditary (Volterra) linear viscoelasticity of aging concrete — the
// rate-type Kelvin-chain engine and its integral-form reference.
//
// PHYSICS. Under a uniaxial stress history sigma(t') the delayed strain of a
// linear-viscoelastic aging material is the Boltzmann/Volterra integral
//
//     eps(t) = INT_{0}^{t} J(t, t') dsigma(t')          (+ shrinkage/thermal)
//
// with J(t, t') the COMPLIANCE FUNCTION [1/Pa]: the strain at age t under a
// unit stress applied at age t' and held constant. Every creep code of
// `keta/creep.h` / `keta/matlaws.h` writes J in the same shape,
//
//     J(t, t') = 1/E(t') + phi(t, t') / E_ref(t')
//
// (EN 1992-1-1 eq. 3.6 for the EC2 law, with E_ref = Ec = 1.05 Ecm; see the
// "WHAT phi() IS REFERENCED TO" table of `matlaws.h` for the other codes).
//
// TWO SOLVERS, ONE REFERENCE. This module provides
//
//   * the INTEGRAL form (`volterra_strain` / `volterra_stress` /
//     `relaxation_function`): a direct trapezoidal quadrature of the
//     hereditary integral on a user-supplied (log) time grid. It stores the
//     whole stress history and costs O(n^2) — it is the REFERENCE against
//     which the rate-type engine below is verified, not a production solver;
//
//   * the RATE-TYPE form (`KelvinChain` + `kelvin_step`): the compliance is
//     approximated by a Kelvin (retardation) chain
//
//         J(t, t') ~ a0 + SUM_mu a_mu * (1 - exp(-(t - t')/tau_mu))
//
//     with FIXED, log-spaced retardation times tau_mu and non-negative unit
//     compliances a_mu = 1/E_mu. A chain has a finite internal state (one
//     strain gamma_mu per unit), so a step-by-step analysis needs no stress
//     history at all: the memory lives in gamma. The step update is EXACT for
//     a stress varying linearly over the step (the "exponential algorithm",
//     Bazant & Wu 1974; Jirasek & Bazant, *Inelastic Analysis of Structures*,
//     ch. 27):
//
//         gamma_mu^{n+1} = e^{-y} gamma_mu^n + (sigma_n/E_mu)(1 - e^{-y})
//                          + (dsigma/E_mu)(1 - lambda),
//         y = dt/tau_mu,  lambda = (1 - e^{-y})/y,
//
//     giving the incremental stress-strain law
//
//         deps = dsigma / E'' + deps''
//         1/E'' = a0 + SUM_mu a_mu (1 - lambda_mu)
//         deps'' = SUM_mu (1 - e^{-y_mu}) (sigma_n a_mu - gamma_mu^n).
//
// FITTING THE CHAIN. `fit_kelvin_chain` fits the a_mu by NON-NEGATIVE least
// squares (Lawson & Hanson, *Solving Least Squares Problems*, 1974, ch. 23) on
// a log-spaced sample of (t - t'), with a RELATIVE weighting 1/J^2 so the fit
// error is uniform in relative terms over the whole window. Non-negativity is
// not cosmetic: a negative unit compliance is a thermodynamically inadmissible
// chain (negative delayed stiffness) and produces non-monotone creep.
// The alternative classical route is the CONTINUOUS RETARDATION SPECTRUM of
// Bazant & Xi, "Continuous retardation spectrum for solidification theory of
// concrete creep", J. Eng. Mech. 121(2), 1995, 281-288 — the Widder inversion
//     L(tau) = lim_{k->inf} (-k tau)^k / (k-1)! * J^{(k)}(k tau),
//     a_mu = L(tau_mu) ln(10) Dlog(tau),
// evaluated in practice at k = 3. It needs the k-th derivative of J in closed
// form (available for the power-type B3 / MC90 terms, NOT for the EC2
// beta_c((t-t0)/(betaH + t-t0))^0.3 product), so keta fits by least squares
// instead and reports the achieved error in `KelvinChain::fit_max_rel_error`.
// Both routes place the SAME log-spaced tau grid; only the a_mu differ.
//
// AGING. The unit compliances of a real concrete depend on the age at loading.
// A chain fitted at one loading age t0 is EXACT (to the fit error) for any
// stress history applied at that age, and is what `fit_kelvin_chain` returns.
// A history with SEVERAL loading events superposes, per Boltzmann (McHenry
// 1943), one chain per event, each fitted at its own event age — that is what
// `phasing.cpp` does (one chain per loading event, fits shared through a
// (material, age) cache; see the CreepEngine::Kelvin note there). Against the
// Volterra integral on a stress history, forcing a whole history through ONE
// frozen-age chain errs by ~13 %, an error the per-event superposition
// removes. Stepping ONE state set
// through a chain whose moduli are re-evaluated at the CURRENT age is the
// classical aging Kelvin chain and is NOT done here: it does not reproduce
// J(t, t') for t' later than the fit age, and the consistent treatment needs
// solidification theory (Bazant & Prasannan, J. Eng. Mech. 115(8), 1989).
//
// UNITS. Ages and times in DAYS (the creep.h convention), stresses in Pa,
// compliances in 1/Pa, strains dimensionless. Every entry point rejects a
// non-finite / out-of-order input with a typed `Error`.
#pragma once
#include <functional>

#include "keta/matlaws.h"

namespace keta {

// ---------------- compliance functions ----------------

// J(t, t0) [1/Pa], ages in days, requires t >= t0 > 0.
using ComplianceFn = std::function<double(double, double)>;

// EC2 compliance anchored on ONE elastic modulus E [Pa]:
//     J(t, t0) = 1/E + phi(t, t0)/E        (EN 1992-1-1 eq. 3.6 with Ec := E)
// This is the form a keta FE model implies: the bar stiffness uses the single
// constant `Material::E`, so both the instantaneous and the creep term must be
// referenced to it for the rate-type engine to agree with the AAEM increment
// it replaces. Use `matlaw_compliance` for the age-dependent modulus.
ComplianceFn ec2_compliance(const ConcreteTime& c, double E);

// Compliance of any `MatLaw` (keta/matlaws.h):
//     J(t, t0) = 1/E(t0) + phi(t, t0)/E_ref(t0)
// with each law's own reference modulus (`MatLaw::E_ref`).
ComplianceFn matlaw_compliance(const MatLaw& law);

// ---------------- Kelvin chain ----------------

struct KelvinFitOptions {
  int n_units = 8;          // number of Kelvin units (Bazant: >= 1 per decade)
  double tau_min = 1e-2;    // shortest retardation time [days]
  double decades = 1.0;     // log10 spacing between retardation times
  double fit_min = 1e-2;    // fit window on the elapsed time (t - t0) [days]
  double fit_max = 1e5;
  int n_samples = 120;      // log-spaced samples in the fit window
};

// An aging-free Kelvin chain fitted at ONE loading age.
struct KelvinChain {
  double t0 = 0;                 // loading age the fit was made at [days]
  double a0 = 0;                 // instantaneous compliance 1/E0 [1/Pa]
  std::vector<double> tau;       // retardation times [days]
  std::vector<double> a;         // unit compliances 1/E_mu [1/Pa], >= 0
  double fit_max_rel_error = 0;  // max |J_chain/J - 1| over the fit window

  // Chain compliance at elapsed time dt = t - t0 [days].
  double J(double dt) const;
  // Delayed part only (without a0).
  double J_delayed(double dt) const;
};

// Fit `n_units` non-negative unit compliances to J(t0 + dt, t0). a0 is set
// EXACTLY to J(t0, t0) (the instantaneous compliance is never fitted).
KelvinChain fit_kelvin_chain(const ComplianceFn& J, double t0, const KelvinFitOptions& o = {});

// ---------------- exponential algorithm ----------------

// Internal state of one chain driven by ONE stress component.
struct KelvinState {
  std::vector<double> gamma;  // unit strains [-]
  double sigma = 0;           // stress at the end of the last step [Pa]
  double t = 0;               // time at the end of the last step [days]
};
KelvinState kelvin_state(const KelvinChain& c, double t_start = 0);

struct KelvinStep {
  double d_eps_total = 0;    // instantaneous + delayed strain increment [-]
  double d_eps_delayed = 0;  // delayed (Kelvin units) part only [-]
  double e_algorithmic = 0;  // E'' = dsigma / (deps - deps'') [Pa]
};

// Exact one-step update under a stress varying LINEARLY from `st.sigma` to
// `sigma_end` over [st.t, t_end]. Updates `st`.
KelvinStep kelvin_step(KelvinState& st, const KelvinChain& c, double t_end, double sigma_end);

// Delayed strain increment over [st.t, t_end] under a stress held CONSTANT at
// `st.sigma` (the piecewise-constant history a staged-construction driver
// sees between two phase dates). Updates `st`; returns deps''.
double kelvin_creep(KelvinState& st, const KelvinChain& c, double t_end);

// ---------------- Volterra reference ----------------

// Log-spaced grid t0, t0 + dt_min, ... , t0 + dt_max with `n` points after t0.
std::vector<double> log_time_grid(double t0, double dt_min, double dt_max, int n);

// Strain history for a given stress history (trapezoidal hereditary integral).
// t must be strictly increasing with t[0] > 0; sigma[0] is applied AT t[0].
std::vector<double> volterra_strain(const ComplianceFn& J, const std::vector<double>& t,
                                    const std::vector<double>& sigma);

// Inverse problem: stress history for a given strain history (same
// quadrature, solved by forward substitution — Bazant's step-by-step
// relaxation algorithm).
std::vector<double> volterra_stress(const ComplianceFn& J, const std::vector<double>& t,
                                    const std::vector<double>& eps);

// Relaxation function R(t, t0): stress under a UNIT strain imposed at t0 and
// held. R[0] = 1/J(t0, t0). `t` must start at t0.
std::vector<double> relaxation_function(const ComplianceFn& J, double t0,
                                        const std::vector<double>& t);

// Non-negative least squares min ||A x - b||, x >= 0 (Lawson & Hanson 1974,
// algorithm NNLS). Exposed for direct use; deterministic (no random restarts).
VecX nnls(const MatX& A, const VecX& b, int max_iter = 0);

}  // namespace keta

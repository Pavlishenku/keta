// Material-nonlinear reinforced-concrete beam, solved EXACTLY where the
// statics allow it (e.g. the VDI 6201 example beams):
//
//  * TabulatedLaw — a uniaxial stress-strain table: points flagged POL
//    (polygon) or SPL (spline); each maximal
//    run of SPL points, extended by its POL neighbour on each side, is
//    interpolated by ONE natural cubic spline through that extended run,
//    every other segment linearly. Beyond the table: typed
//    rejection. Evaluation O(log n) (binary search).
//  * concrete_sargin — EN 1992-1-1 §3.1.5 eq. (3.14), compression only (no
//    tensile strength), valid up to eps_cu1; section integrals in CLOSED FORM
//    (polynomial division + logarithm), no integration-point error.
//  * ts_steel_heft525 — the tension-stiffened ("mean strain") steel law of
//    DIN 1045-1 / DAfStb Heft 525 (DIN EN 1992-1-1/NA):
//      sigma_sr = fct,eff (1 / rho_eff + alpha_e - 1),  eps_sr1 = fct,eff / Ec,
//      eps_sr2 = sigma_sr / Es, alpha_e = Es / Ec;
//      uncracked       0 <= sigma <= sigma_sr : eps = eps_sr1 sigma / sigma_sr
//      crack formation sigma_sr..1.3 sigma_sr : eps = eps_s2 - [beta_t (sigma -
//                       sigma_sr) + (1.3 sigma_sr - sigma)] / (0.3 sigma_sr) dEps
//      stabilised      1.3 sigma_sr..fy      : eps = eps_s2 - beta_t dEps
//      yielding        fy..ft                : eps = eps_sy - beta_t dEps +
//                       delta_d (1 - sigma_sr / fy)(eps_s2 - eps_sy)
//    (dEps = eps_sr2 - eps_sr1, eps_s2 the bare steel strain). Returned as a
//    polygonal TabulatedLaw (each stage is linear in sigma for a bilinear bare
//    steel). Compression: the bare steel.
//  * RcSection — rectangle b x h of a concrete law + discrete bars, N = 0
//    state for a given curvature (exact root, the concrete in closed form),
//    M(kappa), and its inverse on the ascending branch.
//  * SimplySupportedRcBeam — statically determinate: M(x) = M_perm(x) +
//    P m_P(x) is known exactly for a FIXED permanent load (uniform w) under
//    the controlled point load P, so the deflection under the load is the
//    virtual-work integral  delta(P) = int kappa(M(x)) mbar(x) dx — integrated
//    by adaptive Gauss-Legendre on the sub-intervals between the abscissae
//    where M(x) crosses a kink of kappa(M) (a breakpoint of a law reached by
//    a fibre): exact to round-off, no element discretisation. Euler-Bernoulli
//    (N = 0, so second-order effects vanish identically).
#pragma once
#include <functional>
#include <vector>

#include "keta/core.h"

namespace keta::rcnl {

class TabulatedLaw {
 public:
  struct Point {
    double eps = 0, sig = 0;
    bool spline = false;  // "SPL" (else "POL")
  };
  // Points in increasing strain over the whole range (tension and
  // compression), the origin included or bridged linearly by its neighbours.
  explicit TabulatedLaw(std::vector<Point> pts);
  // Symmetric table from the tension branch (eps > 0), mirrored.
  static TabulatedLaw symmetric(const std::vector<Point>& tension);
  double sigma(double eps) const;
  double tangent(double eps) const;
  double eps_min() const { return x_.front(); }
  double eps_max() const { return x_.back(); }
  const std::vector<double>& breakpoints() const { return x_; }

 private:
  std::vector<double> x_, y_;
  std::vector<double> m2_;  // spline second derivatives at the points (0 on linear segments)
  std::vector<char> seg_spline_;  // segment k (x_k, x_k+1) is a spline piece
};

// EN 1992-1-1 eq. (3.14), compression negative: sigma(eps) for eps in
// [-eps_cu1, 0], 0 in tension; k = 1.05 Ecm eps_c1 / fcm.
struct Sargin {
  double fcm = 0, Ecm = 0, eps_c1 = 0, eps_cu1 = 0;
  double k() const { return 1.05 * Ecm * eps_c1 / fcm; }
  double sigma(double eps) const;
  // closed-form antiderivatives in eta = -eps / eps_c1 (compressive):
  // F0(eta) = int_0^eta s(t) dt, F1(eta) = int_0^eta t s(t) dt, s = sigma / -fcm
  void integrals(double eta, double& F0, double& F1) const;
};
// EN 1992-1-1 Table 3.1 from fcm (Ecm in Pa, strains positive): Ecm = 22
// (fcm / 10)^0.3 GPa, eps_c1 = 0.7 fcm^0.31 per mille (<= 2.8), eps_cu1 = 3.5
// per mille up to fck = 50 MPa (fck = fcm - 8).
Sargin sargin_from_fcm(double fcm);

struct Heft525 {
  double Es = 200e9, fy = 500e6, ft = 525e6, eps_su = 0.025;  // bare bilinear steel
  double fct_eff = 0, rho_eff = 0, Ec = 0;                   // cracking
  double beta_t = 0.25, delta_d = 0.8;
};
TabulatedLaw ts_steel_heft525(const Heft525& p);

struct Bar {
  double z = 0;   // from the centroid, upward [m]
  double As = 0;  // [m2]
  const TabulatedLaw* law = nullptr;
};

class RcSection {
 public:
  RcSection(double b, double h, Sargin concrete, std::vector<Bar> bars);
  // N = 0 state for a sagging curvature kappa >= 0 (top compressed):
  // strain eps(z) = eps0 - kappa z. Returns M (sagging positive) and eps0.
  double moment(double kappa, double* eps0 = nullptr) const;
  // curvature of the moment peak on [0, kappa_limit], kappa_limit the largest
  // admissible curvature (top fibre at -eps_cu1 or a bar at the end of its
  // table, whichever first); = kappa_limit when M still grows there: the
  // load-controlled limit of the section
  double kappa_peak() const { return kpeak_; }
  // inverse on [0, kappa_peak], where M(kappa) is verified increasing
  double kappa_of(double M) const;
  // curvatures at which some fibre reaches a law breakpoint (kinks of M)
  std::vector<double> kink_curvatures() const;

 private:
  double axial(double eps0, double kappa) const;
  double b_, h_;
  Sargin c_;
  std::vector<Bar> bars_;
  double klim_ = 0, kpeak_ = 0;
};

class SimplySupportedRcBeam {
 public:
  // span L, permanent uniform load w [N/m] (fixed), point load P at a
  SimplySupportedRcBeam(const RcSection& s, double L, double w, double a);
  // deflection under the point load for P, total (permanent included)
  double deflection(double P) const;
  // load-controlled limit: the section under the load at its moment peak
  double P_limit() const;
  double moment(double P, double x) const;

 private:
  double integrate(double P, double x0, double x1, double tol) const;
  const RcSection& s_;
  double L_, w_, a_;
  std::vector<double> Mk_;  // moments at the kinks of kappa(M)
};

}  // namespace keta::rcnl

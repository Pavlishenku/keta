// Reinforced-concrete section analysis to EN 1992-1-1 — RAW VALUES ONLY.
//
// This module computes numbers (neutral axes, stresses, crack widths,
// curvatures, resistances). It NEVER returns a verdict: no utilisation ratio,
// no pass/fail, no partial-factor policy beyond the explicit `gamma_c` /
// `gamma_s` / `alpha_cc` the caller hands in. Design checking belongs to the
// layer above.
//
// WHAT IS IMPLEMENTED
//   * SLS state I  (uncracked transformed section)  — §7.1, §7.2
//   * SLS state II (fully cracked, concrete in tension ignored) — §7.1(2)
//   * cracking moment M_cr from f_ct,eff on the state-I section
//   * crack width w_k — §7.3.4 eq. (7.8), (7.9), (7.11)/(7.14) with the
//     rho_p,eff / h_c,ef machinery of §7.3.2(3) eq. (7.10)
//   * tension-stiffening curvature interpolation — §7.4.3 eq. (7.18)/(7.19)
//     and a simply-supported deflection integrator built on it
//   * ULS bending resistance M_Rd with the rectangular stress block — §3.1.7(3)
//
// GEOMETRY AND SIGNS (fixed once, used by every entry point)
//   * z is measured DOWNWARD from the top (compression) fibre, in metres;
//   * the section is a rectangle b x h, optionally with a FLANGE bf x hf at the
//     top (T beam, flange in compression). bf <= b or hf <= 0 means rectangular;
//   * tension steel `As` at depth `d`, optional compression steel `As2` at `d2`;
//   * actions: N tension-POSITIVE, applied at the GROSS CONCRETE centroid;
//     M positive = SAGGING (compression at the top fibre);
//   * stresses tension-positive: a compressed top fibre reports sigma_c < 0,
//     the tension bars report sigma_s > 0;
//   * curvature kappa > 0 for sagging (strain increasing downward).
//
// MODULAR RATIO. Every SLS entry point takes a creep coefficient `phi`;
// alpha_e = Es / Ec,eff with Ec,eff = Ecm/(1 + phi) (EN 1992-1-1 eq. 7.20).
// phi = 0 is the short-term (instantaneous) analysis. In the cracked
// transformed section the tension steel is counted with alpha_e and the
// compression steel with (alpha_e - 1), the concrete of the block being
// counted once — the classical convention.
//
// DOCUMENTED LIMITS
//   * uniaxial bending about the horizontal axis, symmetric reinforcement
//     layers reduced to two (tension / compression);
//   * state II assumes a linear concrete law in compression (SLS) — the
//     nonlinear §3.1.5 branch is not used;
//   * the state-II neutral axis is found by bisection on a sign-changing
//     bracket scanned over (0, 4h]; a section that cannot equilibrate the
//     given (N, M) in a cracked state (e.g. pure tension across the whole
//     depth) is a typed `Error`, never a silent number;
//   * the flange is assumed to be on the COMPRESSION side. A T beam in hogging
//     is modelled by passing its web as the rectangle.
#pragma once
#include <functional>

#include "keta/creep.h"

namespace keta {

// ---------------- inputs ----------------

struct RcSection {
  double b = 0;    // web width [m]
  double h = 0;    // total depth [m]
  double bf = 0;   // effective flange width [m] (<= b or hf <= 0 -> rectangular)
  double hf = 0;   // flange thickness [m]
  double d = 0;    // effective depth of the tension steel [m]
  double As = 0;   // tension steel area [m^2]
  double d2 = 0;   // depth of the compression steel [m]
  double As2 = 0;  // compression steel area [m^2]

  void validate() const;
  double width(double z) const;          // section width at depth z [m]
  double conc_area(double t) const;      // INT_0^t width dz [m^2]
  double conc_first_moment(double t) const;  // INT_0^t width z dz [m^3]
  double conc_second_moment(double t) const; // INT_0^t width z^2 dz [m^4]
  double area() const { return conc_area(h); }
  double centroid() const;               // gross concrete centroid depth [m]
};

struct RcConcrete {
  double fck = 0;        // characteristic cylinder strength [Pa]
  double fctm = 0;       // 0 -> EN 1992-1-1 Table 3.1
  double Ecm = 0;        // 0 -> ec2_ecm(fck)
  double gamma_c = 1.5;  // ULS partial factor
  double alpha_cc = 1.0; // EN 1992-1-1 §3.1.6(1)
};

struct RcSteel {
  double fyk = 500e6;
  double Es = 200e9;
  double gamma_s = 1.15;
};

struct RcMaterials {
  RcConcrete c;
  RcSteel s;
};

// EN 1992-1-1 Table 3.1: fctm = 0.30 fck^(2/3) (<= C50/60),
// else 2.12 ln(1 + fcm/10), fcm = fck + 8 MPa. fck in Pa, result in Pa.
double ec2_fctm(double fck);
// Table 3.1 ultimate strain of the rectangular-block design law eps_cu3.
double ec2_eps_cu3(double fck);
// §3.1.7(3): block depth factor lambda and strength factor eta.
double ec2_block_lambda(double fck);
double ec2_block_eta(double fck);

double rc_fctm(const RcConcrete& c);           // c.fctm or the Table 3.1 value
double rc_ecm(const RcConcrete& c);            // c.Ecm or ec2_ecm(fck)
double rc_ec_eff(const RcConcrete& c, double phi);   // Ecm / (1 + phi)  (eq. 7.20)
double rc_alpha_e(const RcMaterials& m, double phi); // Es / Ec,eff

// ---------------- SLS state I (uncracked) ----------------

struct RcStateI {
  double alpha_e = 0;
  double A = 0;            // transformed area [m^2]
  double x = 0;            // transformed centroid depth [m]
  double I = 0;            // transformed second moment about it [m^4]
  double sigma_c_top = 0;  // [Pa] (negative = compression)
  double sigma_c_bot = 0;
  double sigma_s = 0;      // tension steel [Pa]
  double sigma_s2 = 0;     // compression steel [Pa]
  double kappa = 0;        // curvature [1/m]
};
RcStateI rc_uncracked(const RcSection& s, const RcMaterials& m, double N, double M, double phi);

// Cracking moment: the sagging moment that brings the bottom fibre of the
// state-I section to f_ct,eff (f_ct,eff = fctm unless overridden), at the given
// axial force.
double rc_cracking_moment(const RcSection& s, const RcMaterials& m, double N, double phi,
                          double fct_eff = 0);

// ---------------- SLS state II (fully cracked) ----------------

struct RcStateII {
  double alpha_e = 0;
  double x = 0;            // neutral axis depth from the top fibre [m]
  double I = 0;            // cracked transformed second moment about x [m^4]
  double sigma_c_top = 0;  // [Pa] (negative = compression)
  double sigma_s = 0;      // tension steel [Pa]
  double sigma_s2 = 0;     // compression steel [Pa]
  double kappa = 0;        // curvature [1/m]
  double lever_arm = 0;    // z = M / (As sigma_s) [m] (0 when As = 0)
};
RcStateII rc_cracked(const RcSection& s, const RcMaterials& m, double N, double M, double phi);

// ---------------- crack width (§7.3) ----------------

struct RcCrackParams {
  double cover = 0.030;    // c, cover to the longitudinal reinforcement [m]
  double bar_dia = 0.016;  // equivalent bar diameter phi_eq [m]
  double spacing = 0;      // bar spacing [m]; > 5(c + phi/2) switches to eq. 7.14
  double k1 = 0.8;         // bond: 0.8 high bond, 1.6 plain bars
  double k2 = 0.5;         // strain distribution: 0.5 bending, 1.0 pure tension
  double k3 = 3.4;         // eq. 7.11 (recommended value)
  double k4 = 0.425;
  double kt = 0.6;         // 0.6 short-term, 0.4 long-term loading
  double fct_eff = 0;      // 0 -> fctm
};

struct RcCrackResult {
  double x = 0;            // cracked neutral axis [m]
  double sigma_s = 0;      // steel stress in the cracked state [Pa]
  double hc_ef = 0;        // effective tension height [m]  (§7.3.2(3))
  double Ac_eff = 0;       // effective concrete area in tension [m^2]
  double rho_p_eff = 0;    // eq. 7.10
  double d_eps = 0;        // eps_sm - eps_cm, eq. 7.9
  double s_r_max = 0;      // eq. 7.11 or 7.14 [m]
  bool wide_spacing = 0;   // true when eq. 7.14 was used
  double wk = 0;           // crack width [m]
};
RcCrackResult rc_crack_width(const RcSection& s, const RcMaterials& m, const RcCrackParams& p,
                             double N, double M, double phi);

// ---------------- tension stiffening / deflection (§7.4.3) ----------------

struct RcCurvature {
  double M_cr = 0;
  double kappa_I = 0;   // uncracked curvature [1/m]
  double kappa_II = 0;  // fully cracked curvature [1/m]
  double zeta = 0;      // eq. 7.19 distribution coefficient
  double kappa = 0;     // eq. 7.18 interpolated curvature [1/m]
};
// beta: 1.0 single short-term load, 0.5 sustained / many cycles (§7.4.3(3)).
RcCurvature rc_curvature(const RcSection& s, const RcMaterials& m, double N, double M, double phi,
                         double beta = 0.5);

// Midspan deflection of a simply-supported beam of span L under the moment
// distribution M_of_x(x), by the unit-load integral
//     delta = INT_0^L kappa(x) m(x) dx,   m = x/2 | (L-x)/2,
// with the eq. 7.18 curvature at every station (Simpson on each half span,
// `n_per_half` intervals, forced even). N is held constant along the span.
double rc_deflection_ss(const RcSection& s, const RcMaterials& m, double L,
                        const std::function<double(double)>& M_of_x, double N, double phi,
                        double beta = 0.5, int n_per_half = 40);

// Shrinkage curvature, EN 1992-1-1 eq. 7.21: 1/r_cs = eps_cs alpha_e S / I,
// alpha_e = Es / Ec,eff(phi), S the first moment of the reinforcement about the
// centroid of the transformed section (state I) or about the neutral axis
// (state II), I its second moment. eps_cs < 0 (shortening) gives a sagging
// curvature for bottom steel; returned with the sign of a sagging curvature
// (positive when it adds to a sagging load curvature).
double rc_shrinkage_curvature(const RcSection& s, const RcMaterials& m, double eps_cs, double phi,
                              bool cracked);
// Mid-span deflection of a simply supported beam with the load AND shrinkage
// curvatures interpolated by the same zeta (eq. 7.18 applied to the total
// curvature, the tension-stiffening of shrinkage): delta = INT m(x) [zeta
// (k_II + k_cs,II) + (1 - zeta)(k_I + k_cs,I)] dx, integrated on the pieces
// between the cracking points |M| = M_cr (located on a 256-interval scan of
// M_of_x then bisected to round-off) by adaptive 5-point Gauss-Legendre to
// 1e-12 relative — exact to that tolerance for the smooth integrand there.
// The cracking moment M_cr of zeta (eq. 7.19, sigma_sr / sigma_s = M_cr / M)
// is taken on the transformed section at creep coefficient phi_cr (< 0: phi;
// 0: the short-term section on which the first cracking occurs).
double rc_deflection_ss_cs(const RcSection& s, const RcMaterials& m, double L,
                           const std::function<double(double)>& M_of_x, double phi, double eps_cs,
                           double beta = 0.5, double phi_cr = -1.0);

// ---------------- ULS bending (§3.1.7(3) + §6.1) ----------------

struct RcUls {
  double fcd = 0, fyd = 0;
  double lambda = 0, eta = 0, eps_cu3 = 0;
  double x = 0;         // neutral axis depth [m]
  double block = 0;     // block depth lambda * x [m]
  double eps_s = 0;     // tension steel strain [-]
  double sigma_s = 0;   // tension steel stress [Pa]
  double sigma_s2 = 0;  // compression steel stress [Pa] (negative = compression)
  double Fc = 0;        // concrete block force [N] (negative = compression)
  double M_Rd = 0;      // resisting moment about the gross concrete centroid [N.m]
  bool steel_yields = false;
};
// Bending resistance at a given axial force N (tension +). Throws when the
// section cannot equilibrate N at all (N beyond the squash/tension capacity).
RcUls rc_mrd(const RcSection& s, const RcMaterials& m, double N = 0.0);

}  // namespace keta

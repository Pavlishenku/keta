// EN 1993-1-1 / EN 1993-1-5 steel member checks — a VERIFICATION layer on
// top of the solver: raw resistances, reduction factors, interaction factors
// and utilisation ratios, every formula cited by its clause. No verdict is
// taken here (a ratio is returned, the caller compares it with 1). SI units
// throughout (m, N, Pa); fy is the caller's (EN 1993-1-1 Table 3.1), so are
// the partial factors.
//
// Exactness: these are the code's own formulas (definitions of the norm, not
// keta simplifications). Inputs the code leaves to "an appropriate method"
// (Mcr of a real moment diagram, Ncr,T, delta_x of Table A.2) come from keta's
// solvers (warping.h, statics.h) — the caller passes them in.
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta::ec3 {

// ---------------------------------------------------------------- materials
// Table 5.2: epsilon = sqrt(235 MPa / fy).
double epsilon(double fy);

// ------------------------------------------------- classification (§5.5)
// Table 5.2 sheet 1, internal compression part of slenderness c/t:
//   alpha = compressed fraction of c under the PLASTIC distribution (classes
//   1, 2), psi = sigma_2 / sigma_1 of the ELASTIC one (class 3), sigma_1 the
//   largest compression. Returns 1..4 (4 = beyond the class-3 limit, 5.5.2(8)).
int class_internal(double c_t, double eps, double alpha, double psi);
// Table 5.2 sheet 2, outstand flange (rolled or welded `rolled`): uniform
// compression limits 9 / 10 / 14 eps; a stress gradient through alpha (tip in
// compression: 9 eps/alpha, 10 eps/alpha; tip in tension: 9 eps/(alpha
// sqrt(alpha)), ...) and class 3 by 21 eps sqrt(k_sigma) (EN 1993-1-5 Table 4.2).
int class_outstand(double c_t, double eps, double alpha, double psi, bool tip_in_compression);

// Rolled doubly-symmetric I/H section (root radius r).
struct ISection {
  double h = 0, b = 0, tw = 0, tf = 0, r = 0;
};
// Table 5.2: c = b/2 - tw/2 - r (outstand flange), c = h - 2 tf - 2 r (web).
double c_flange(const ISection& s);
double c_web(const ISection& s);

// -------------------------------------------- EN 1993-1-5 plate buckling
// Table 4.1 (internal) and Table 4.2 (outstand) buckling factors k_sigma.
double ksigma_internal(double psi);
double ksigma_outstand(double psi, bool max_compression_at_free_edge);
// §4.4(2): lambda_p = (b/t) / (28.4 eps sqrt(k_sigma)).
double lambda_p(double b_t, double eps, double ksigma);
// Eq. (4.2) internal, (4.3) outstand: reduction factor rho <= 1.
double rho_internal(double lambda_p, double psi);
double rho_outstand(double lambda_p);

// Annex A.2.2(1): elastic critical stress of a plate stiffener as a column on
// the elastic foundation of the plate, a <= a_c = 4.33 (I b1^2 b2^2 / (t^3 b))^(1/4):
//   sigma_cr,sl = pi^2 E I_sl / (A_sl a^2) + E t^3 b a^2 / (4 pi^2 (1 - nu^2) A_sl b1^2 b2^2)
// b1, b2: distances from the stiffener to the longitudinal plate edges (b =
// b1 + b2). a > a_c is a typed rejection (A.2.2(1) second formula not here).
double sigma_cr_sl(double E, double nu, double t, double a, double b1, double b2, double A_sl,
                   double I_sl);
double a_c(double t, double b1, double b2, double I_sl);
// A flat stiffener hs x ts welded on one face of a plate of thickness t, with
// the plate strip `bp` (A.2.1(3) Figure A.1: (3 - psi)/(5 - psi) b1 +
// 2/(5 - psi) b2 plus ts): gross A_sl and I_sl about the combined centroid.
void flat_stiffener(double t, double bp, double hs, double ts, double& A_sl, double& I_sl);

// Effective section of plates under UNIFORM compression (EN 1993-1-5 §4.3,
// §4.4, psi = 1): the section is a set of rectangles (y0 < y1, z0 < z1) in the
// section axes; each declared internal plate (Table 4.1) loses its central
// strip (1 - rho) b of the clear width b centred at (yc, zc), oriented along y
// (`along_y`) or z. Exact rectangle integrals.
struct Rect {
  double y0 = 0, y1 = 0, z0 = 0, z1 = 0;
};
struct InternalPlate {
  double b = 0, t = 0;   // clear width and thickness (Table 4.1 b-bar)
  double yc = 0, zc = 0; // centre of the clear width
  bool along_y = true;   // the plate's width runs along y (a flange) or z (a web)
};
struct EffectiveSection {
  double A = 0, Aeff = 0;
  double yg = 0, zg = 0;          // gross centroid
  double yg_eff = 0, zg_eff = 0;  // effective centroid
  double Iy_eff = 0, Iz_eff = 0;  // about the effective centroid
  std::vector<double> rho, lambda_p, removed;  // per plate
};
EffectiveSection effective_uniform_compression(const std::vector<Rect>& gross,
                                               const std::vector<InternalPlate>& plates,
                                               double fy);
// Second moments of a set of rectangles minus removed strips about (y, z).
struct RectProps {
  double A = 0, yg = 0, zg = 0, Iy = 0, Iz = 0;
};
RectProps rect_props(const std::vector<Rect>& add, const std::vector<Rect>& remove = {});

// --------------------------------------------- cross-section resistance
// §6.2.6(2) eq. 6.18: Vpl,Rd = Av fy / (sqrt 3 gamma_M0).
double vpl_rd(double Av, double fy, double gM0);
// §6.2.6(3)a: rolled I/H, load parallel to the web:
// Av = A - 2 b tf + (tw + 2 r) tf, not less than eta hw tw.
double av_z_rolled_i(const ISection& s, double A, double eta = 1.0);
// Load parallel to the flanges: the two flanges, Av = 2 b tf (EN 1993-1-1
// §6.2.6(3) gives no rolled-I rule for this direction; the flange area is the
// shear area of the classical tables).
double av_y_flanges(const ISection& s);
// §6.2.8(3) / §6.2.10(3): rho = (2 VEd / Vpl,Rd - 1)^2 when VEd > 0.5 Vpl,Rd, else 0.
double rho_shear(double VEd, double VplRd);
// §6.2.8(5) eq. 6.30: My,V,Rd = (Wpl,y - rho Aw^2 / (4 tw)) fy / gamma_M0, Aw = hw tw.
double my_v_rd_i(double Wply, double rho, double hw, double tw, double fy, double gM0);
// §6.2.8(3): reduced yield (1 - rho) fy on the shear area of Vy, the two
// flanges: Mz,V,Rd = (Wpl,z - rho 2 tf b^2 / 4) fy / gamma_M0.
double mz_v_rd_i(double Wplz, double rho, const ISection& s, double fy, double gM0);
// §6.2.10(3): N with (1 - rho) fy on the shear area: NV,Rd = (A - rho Av) fy / gamma_M0.
double nv_rd(double A, double Av, double rho, double fy, double gM0);
// §6.2.9.1(5) eq. 6.36 (I/H, major axis): MN,y,Rd = M (1 - n) / (1 - 0.5 a) <= M,
// a = min((A - 2 b tf) / A, 0.5).
double mn_y_rd_i(double M, double n, double a);
// §6.2.9.1(4) eqs 6.33 / 6.34 with (1 - rho) fy on the web (§6.2.10(3)):
// the axial force may be neglected about y-y when NEd <= 0.25 N and
// NEd <= 0.5 hw tw (1 - rho) fy / gamma_M0.
bool axial_negligible_y(double NEd, double N, double hw, double tw, double rho, double fy,
                        double gM0);

// ------------------------------------------------- member buckling (§6.3)
enum class Curve { a0, a, b, c, d };
double imperfection(Curve c);  // Table 6.1 / 6.3: 0.13, 0.21, 0.34, 0.49, 0.76
// Table 6.2, rolled I sections (fy < 460 MPa column); welded I sections.
Curve curve_flexural_i(bool rolled, double h, double b, double tf, bool about_y);
// Table 6.4 (general case) and Table 6.5 (rolled or equivalent welded, §6.3.2.3).
Curve curve_lt_general_i(bool rolled, double h, double b);
Curve curve_lt_rolled_i(bool rolled, double h, double b);
// Eq. 6.49 / 6.56: Phi = 0.5 (1 + alpha (lambda - 0.2) + lambda^2),
// chi = 1 / (Phi + sqrt(Phi^2 - lambda^2)) <= 1.
double chi(double lambda, double alpha);
// Eq. 6.57 (§6.3.2.3): Phi = 0.5 (1 + alpha (lambda - lambda0) + beta lambda^2),
// chi = 1 / (Phi + sqrt(Phi^2 - beta lambda^2)) <= 1 and <= 1 / lambda^2.
double chi_lt_rolled(double lambda, double alpha, double lambda0 = 0.4, double beta = 0.75);
// Eq. 6.58: f = 1 - 0.5 (1 - kc) (1 - 2 (lambda - 0.8)^2) <= 1.
double f_modification(double kc, double lambda_lt);
// Ncr = pi^2 E I / Lcr^2; torsional, doubly symmetric section, fork supports:
// Ncr,T = A / (Iy + Iz) (G IT + pi^2 E Iw / LT^2).
double ncr(double EI, double Lcr);
double ncr_torsional(double A, double Iy, double Iz, double G, double IT, double E, double Iw,
                     double LT);
// DIN 18800-2 eq. (20) (Schneider, Bautabellen): Mcr = zeta Ncr,z (sqrt(c^2 +
// 0.25 zp^2) + 0.5 zp), c^2 = (Iw + 0.039 L^2 IT) / Iz, zp the load point
// from the shear centre (negative above it). A code's own closed form.
double mcr_din18800(double zeta, double Ncrz, double Iw, double IT, double Iz, double L, double zp);

// ------------------------------- Annex A (method 1), Tables A.1 and A.2
// Table A.2: Cm,0 for a linear moment (psi = M_small / M_large) and for
// transverse loading with the deflection delta_x and moment M_x at the
// critical section.
double cm0_linear(double psi, double NEd, double Ncr);
double cm0_transverse(double EI, double L, double delta, double M, double NEd, double Ncr);

struct AnnexAInput {
  // section (class 1 or 2: plastic moduli)
  double A = 0, Iy = 0, IT = 0, Wply = 0, Wplz = 0, Wely = 0, Welz = 0;
  double fy = 0, gM1 = 1.0;
  // actions (magnitudes; NEd compression > 0)
  double NEd = 0, MyEd = 0, MzEd = 0;
  // critical values and reduction factors
  double Ncry = 0, Ncrz = 0, NcrT = 0, Mcr0 = 0, C1 = 1;
  double chiy = 1, chiz = 1, chiLT = 1;  // chiLT as used (chi_LT,mod when modified)
  double lambday = 0, lambdaz = 0;
  double Cmy0 = 1, Cmz0 = 1;
};
struct AnnexAResult {
  double mu_y = 0, mu_z = 0, wy = 0, wz = 0, npl = 0, aLT = 0;
  double lambda0 = 0, lambda0_lim = 0, eps_y = 0;
  double Cmy = 0, Cmz = 0, CmLT = 0;
  double bLT = 0, cLT = 0, dLT = 0, eLT = 0;
  double Cyy = 0, Cyz = 0, Czy = 0, Czz = 0;
  double kyy = 0, kyz = 0, kzy = 0, kzz = 0;
  double eq661 = 0, eq662 = 0;  // utilisation of eqs. 6.61 / 6.62
};
AnnexAResult annex_a(const AnnexAInput& in);

// --------------------------------------------- Annex B (method 2)
// Table B.3 Cm for a linear moment: 0.6 + 0.4 psi >= 0.4.
double cm_b3_linear(double psi);
// Table B.3, uniform load between supports with end moments psi Mh, Ms:
// alpha_s = Ms / Mh < 0 ... here the class used by the sheets: alpha_h =
// Mh / Ms with |Mh| <= |Ms|: Cm = 0.95 + 0.05 alpha_h (psi >= 0 ... uniform load)
// and for the LT factor 0.2 + 0.8 alpha_s >= 0.4 (Table B.3 row 2, alpha_s >= 0).
double cm_b3_uniform_load_ah(double alpha_h);
double cm_b3_uniform_load_as(double alpha_s);
struct AnnexBInput {
  double NEd = 0, NRk = 0, gM1 = 1.0;
  double chiy = 1, chiz = 1, lambday = 0, lambdaz = 0;
  double Cmy = 1, CmLT = 1;
  bool class34 = false;  // Table B.1/B.2 class 3, 4 rows (elastic)
};
// Table B.1 kyy (I sections, class 1/2: Cmy (1 + (lambda_y - 0.2) n) <= Cmy (1 + 0.8 n);
// class 3/4: Cmy (1 + 0.6 lambda_y n) <= Cmy (1 + 0.6 n)), n = NEd / (chi NRk / gM1).
double kyy_b(const AnnexBInput& in);
// Table B.2 (members susceptible to torsional deformation), class 1/2:
// lambda_z >= 0.4: kzy = 1 - 0.1 lambda_z n / (CmLT - 0.25) >= 1 - 0.1 n / (CmLT - 0.25);
// lambda_z < 0.4: kzy = 0.6 + lambda_z <= 1 - 0.1 lambda_z n / (CmLT - 0.25).
double kzy_b2(const AnnexBInput& in);

}  // namespace keta::ec3

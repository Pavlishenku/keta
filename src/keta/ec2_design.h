// EN 1992-1-1 (and its German National Annex, DIN EN 1992-1-1/NA) ULS design
// formulas beyond bending: shear §6.2, torsion §6.3, web-flange shear §6.2.4,
// interface shear §6.2.5, compression-reinforcement bending design. Raw
// values, no verdict — like rc.h. SI units (m, N, Pa); fck etc. in Pa. Every
// function cites its clause; the DIN NA variants are separate functions
// named *_din (the NA's own equations: 6.7aDE, 6.7bDE, NA.6.28.1, ...).
#pragma once
#include "keta/core.h"

namespace keta::ec2 {

// ------------------------------------------------------------- shear §6.2
// §6.2.2(1) eq. 6.2a/6.2b: VRd,c = max(CRd,c k (100 rho_l fck)^(1/3) + k1 sigma_cp,
// v_min + k1 sigma_cp) bw d. k = 1 + sqrt(200 mm / d) <= 2, rho_l <= 0.02,
// sigma_cp = NEd / Ac (compression > 0) < 0.2 fcd. EN: CRd,c = 0.18 / gamma_c,
// k1 = 0.15, v_min = 0.035 k^1.5 fck^0.5 (6.3N). DIN NA: CRd,c = 0.15 / gamma_c,
// k1 = 0.12, v_min = 0.0525 / gamma_c k^1.5 fck^0.5 for d <= 600 mm (6.3aDE).
struct VrdcInput {
  double bw = 0, d = 0, Asl = 0, fck = 0, gamma_c = 1.5;
  double sigma_cp = 0, fcd = 0;
};
double shear_k(double d);
double vrd_c_en(const VrdcInput& in);
double vrd_c_din(const VrdcInput& in);

// §6.2.2(6) eq. 6.6N: nu = 0.6 (1 - fck / 250 MPa). DIN NA: nu1 = 0.75 nu2,
// nu2 = 1.1 - fck / 500 MPa <= 1.
double nu_en(double fck);
double nu1_din(double fck);
// Lever arm. EN §6.2.3(1): z = 0.9 d. DIN NA §6.2.3(1): z = 0.9 d but not
// more than max(d - c_v,l - 30 mm, d - 2 c_v,l).
double z_din(double d, double c_vl);
// DIN NA eq. 6.7bDE: VRd,cc = c 0.48 fck^(1/3) (1 - 1.2 sigma_cd / fcd) bw z,
// c = 0.5, fck in MPa (sigma_cd > 0 compression); eq. 6.7aDE:
// cot(theta) = (1.2 + 1.4 sigma_cd / fcd) / (1 - VRd,cc / VEd), limited to [1, 3].
double vrd_cc_din(double fck, double sigma_cd, double fcd, double bw, double z);
// (upper bound 3.0 of DIN EN 1992-1-1/NA; DIN EN 1992-2/NA eq. NA.6.107aDE
// caps it at 1.75: pass cot_max)
double cot_theta_din(double sigma_cd, double fcd, double VRdcc, double VEd, double cot_max = 3.0);
// §6.2.4(4) eq. 6.21: transverse flange steel Asf / sf = vEd hf / (cot theta_f fyd);
// §6.3.2(1) eq. 6.26: shear flow of torsion tau_t t_ef = TEd / (2 A_k).
double asf_flange(double vEd, double hf, double cot_theta_f, double fyd);
double torsion_shear_flow(double TEd, double Ak);
// §6.2.3(3) eq. 6.8 (vertical links): Asw / s = VEd / (z fywd cot theta);
// eq. 6.9: VRd,max = alpha_cw bw z nu1 fcd / (cot theta + tan theta).
double asw_s(double VEd, double z, double fywd, double cot_theta);
double vrd_max(double bw, double z, double nu1, double fcd, double cot_theta, double alpha_cw = 1.0);
// §9.2.2(5) eq. 9.5N: rho_w,min = 0.08 sqrt(fck) / fyk (MPa); Asw,min / s =
// rho_w,min bw sin(alpha).
double asw_min_s_en(double fck, double fyk, double bw, double alpha = kPi / 2);

// ------------------------------------------------------------ torsion §6.3
// §6.3.2(1) equivalent thin-walled rectangle: t_ef given (DIN NA: twice the
// distance from the surface to the longitudinal bar axis), A_k = (b - t_ef)
// (h - t_ef), u_k = 2 (b + h - 2 t_ef).
struct TorsionBox {
  double tef = 0, Ak = 0, uk = 0;
};
TorsionBox torsion_box(double b, double h, double tef);
// DIN NA eq. NA.6.28.1: Asw / s = TEd tan(theta) / (fyd 2 A_k);
// eq. 6.28: Asl = TEd u_k cot(theta) / (fyd 2 A_k);
// eq. 6.30: TRd,max = 2 nu alpha_cw fcd A_k t_ef sin(theta) cos(theta).
double asw_s_torsion(double TEd, const TorsionBox& t, double fyd, double cot_theta);
double asl_torsion(double TEd, const TorsionBox& t, double fyd, double cot_theta);
double trd_max(double nu, double fcd, const TorsionBox& t, double cot_theta, double alpha_cw = 1.0);


// ------------------------------------- bending design of a rectangle §6.1
// Parabola-rectangle concrete (§3.1.7(1), eq. 3.17; n, eps_c2, eps_cu2 of
// Table 3.1 for fck), steel of Figure 3.8 with the inclined top branch up to
// eps_ud (fyd at eps_yd = fyd / Es, ftd at eps_ud). Exact stress-block
// integrals (closed forms of the parabola); the strain plane is the root of
// mu(state) = MEds / (b d^2 fcd) in domain 2 (eps_s1 = eps_ud, eps_c <
// eps_cu2) or domain 3 (eps_c = eps_cu2, quadratic in xi). Beyond xi_lim
// (DIN NA 5.4 (NA.5): 0.45 up to C50/60 for a linear-elastic analysis) the
// compression zone is frozen at xi_lim and the remainder dM is carried by a
// steel couple: As1 = (M_lim / (zeta d) + dM / (d - d2) + NEd) / sigma_s1,
// As2 = dM / ((d - d2) |sigma_s2|). MEds = MEd - NEd zs1 (caller).
// eps_lim (> eps_ud, optional) extends the law by a horizontal branch at
// ftd up to the design strain limit eps_lim, which then bounds domain 2
// (e.g. EN 1992 B500B: inclined to 25 per mille, flat to 45).
struct SteelDesign {
  double fyd = 500e6 / 1.15, Es = 200e9, ftd = 525e6 / 1.15, eps_ud = 0.025;
  double eps_lim = 0;  // 0: the law ends at eps_ud
  double limit() const { return eps_lim > eps_ud ? eps_lim : eps_ud; }
};
double steel_stress(const SteelDesign& s, double eps);
// alpha_R and k_a of the parabola-rectangle block for a top strain eps_c:
// resultant = alpha_R b x fcd at k_a x from the compressed edge.
struct Block {
  double alpha = 0, ka = 0;
};
Block parabola_rectangle_block(double fck, double eps_c);
struct RectDesign {
  double mu = 0, mu_lim = 0, xi = 0, zeta = 0, eps_c = 0, eps_s1 = 0, eps_s2 = 0;
  double sigma_s1 = 0, sigma_s2 = 0, As1 = 0, As2 = 0;
};
RectDesign design_rect_bending(double b, double d, double d2, double fck, double fcd,
                               const SteelDesign& st, double MEds, double NEd, double xi_lim);


// ---------------------------------------- interface shear §6.2.5
// Eq. 6.25 (DIN NA: 1.2 mu in the steel term): vRdi = c fctd + mu sigma_n +
// rho fyd (1.2 mu sin(alpha) + cos(alpha)) <= 0.5 nu fcd (the cap is the
// caller's: vrdi_max). `din` selects the NA's 1.2 mu (EN: mu). Required ratio
// rho = (vEdi - c fctd - mu sigma_n) / (fyd (k mu sin(alpha) + cos(alpha))), >= 0.
double rho_interface(double vEdi, double c, double fctd, double mu, double sigma_n, double fyd,
                     double alpha, bool din);
double vrdi_max(double nu, double fcd);


// ----------------------------------- crack control, minimum reinforcement
// §7.3.2(2) eq. 7.1: As,min sigma_s = kc k fct,eff Act. Eq. 7.2 (rectangles,
// webs): kc = 0.4 (1 - sigma_c / (k1 (h / h*) fct,eff)) in [0, 1], sigma_c =
// NEd / (b h) compression > 0, h* = min(h, 1 m), k1 = 1.5 (compression) or
// 2 h* / (3 h) (tension).
double kc_rect(double sigma_c, double h, double fct_eff);
// k for non-uniform self-equilibrating stresses, h = min(b, h) of the part,
// DIN NA: 0.8 (h <= 300 mm) .. 0.5 (h >= 800 mm), linear between.
double k_selfeq_din(double h);
// DIN NA Table 7.2DE (a): sigma_s = sqrt(wk 3.48e6 / phi_s*) [MPa, wk and
// phi_s* in mm, the formula of the table]; phi_s* from the bar diameter:
// for As,min (NA.7.5.2) phi_s = phi_s* fct,eff / 2.9 MPa; under load
// (7.7.1DE) phi_s = phi_s* sigma_s As / (4 (h - d) b 2.9 MPa) >= phi_s* fct,eff / 2.9.
double sigma_s_tab72de(double wk, double phi_star);
double phi_star_min_reinf(double phi_s, double fct_eff);
double phi_star_load(double phi_s, double sigma_s, double As, double h, double d, double b, double fct_eff);
double as_min(double kc, double k, double fct_eff, double Act, double sigma_s);
// §7.3.4 crack width with the DIN NA crack spacing: w_k = s_r,max (eps_sm -
// eps_cm), eq. 7.9: eps_sm - eps_cm = (sigma_s - kt fct,eff / rho_eff (1 +
// alpha_e rho_eff)) / Es >= 0.6 sigma_s / Es, rho_eff = As / Ac,eff (eq. 7.10,
// no prestressing steel); NA eq. 7.11: s_r,max = phi / (3.6 rho_eff) <=
// sigma_s phi / (3.6 fct,eff). Ac,eff = h_c,ef b with h_c,ef of the caller
// (DIN NA Fig. 7.1DE).
double crack_width_din(double sigma_s, double As, double Ac_eff, double phi, double fct_eff, double kt,
                       double alpha_e, double Es);


// ------------------------------------------------------------ fatigue §6.8
// Eq. 6.76: fcd,fat = k1 beta_cc(t0) fcd (1 - fck / 250 MPa) (k1 = 0.85 EN
// recommended, 1.0 DIN NA; beta_cc(t0) of eq. 3.2).
double fcd_fat(double k1, double beta_cc_t0, double fcd, double fck);
// Table 6.3N (EN) / 6.3DE (DIN NA), straight and bent bars: Delta sigma_Rsk
// at N* = 1e6 cycles: 162.5 MPa (EN), 175 MPa (DIN).
double dsigma_rsk_bars(bool din);
// §6.8.5 / NA: equivalent damage stress range of the simplified method,
// Delta sigma_s,equ(N*) = lambda_s (sigma_s,max - sigma_s,min).
double dsigma_s_equ(double lambda_s, double sigma_max, double sigma_min);


// ------------------------------------ slender columns, §5.2 and §5.8.8
// Eq. 5.1: theta_i = theta_0 alpha_h alpha_m, theta_0 = 1/200 (recommended),
// alpha_h = 2 / sqrt(l) in [2/3, 1] (l in m), alpha_m = sqrt(0.5 (1 + 1/m));
// eq. 5.2 (isolated member): e_i = theta_i l0 / 2.
double imperfection_theta(double l, double m, double theta0 = 1.0 / 200);
double imperfection_ei(double theta_i, double l0);
// §5.8.8.3 eq. 5.34: 1/r = Kr K_phi eps_yd / (0.45 d); eq. 5.33: e2 = (1/r)
// l0^2 / c with c = 10 (pi^2 for a sinusoidal curvature, K1 = 1 above
// lambda = 35 — the caller's K1 for 15 < lambda < 35).
double nominal_curvature(double Kr, double Kphi, double eps_yd, double d);
double second_order_e2(double K1, double inv_r, double l0, double c = 10.0);


// --------------------------------------------------------- punching §6.4
// §6.4.2 basic control perimeter u1 at 2d of a rectangular loaded area a x b
// (interior column): 2 (a + b) + 4 pi d; u0 = 2 (a + b).
double punching_u1_rect(double a, double b, double d);
// §6.4.4(1) eq. 6.47: vRd,c = CRd,c k (100 rho_l fck)^(1/3) + k1 sigma_cp >=
// v_min + k1 sigma_cp; rho_l = sqrt(rho_ly rho_lz) <= 0.02 (DIN NA also <=
// 0.5 fcd / fyd); k = 1 + sqrt(200 mm / d) <= 2; k1 = 0.1. EN: CRd,c = 0.18 /
// gamma_c, v_min = 0.035 k^1.5 fck^0.5; DIN NA: CRd,c = 0.18 / gamma_c for
// u0 / d >= 4, else 0.18 / gamma_c (0.1 u0 / d + 0.6), v_min of 6.3aDE/6.3bDE
// (0.0525 / gamma_c for d <= 600 mm, 0.0375 / gamma_c for d > 800 mm).
struct PunchingInput {
  double d = 0, rho_ly = 0, rho_lz = 0, fck = 0, fcd = 0, fyd = 0, gamma_c = 1.5;
  double sigma_cp = 0, u0 = 0;
};
double vrd_c_punching_en(const PunchingInput& in);
double vrd_c_punching_din(const PunchingInput& in);
// DIN NA 6.4.3(6) (braced systems, adjacent spans within 25 %): beta = 1.10
// interior, 1.40 edge, 1.50 corner column, 1.35 wall end, 1.20 wall corner.
enum class PunchingPosition { Interior, Edge, Corner, WallEnd, WallCorner };
double beta_punching_din(PunchingPosition p);
// eq. 6.38: vEd = beta VEd / (u_i d).
double ved_punching(double beta, double VEd, double u, double d);
// DIN NA eq. NA.6.53.1: vRd,max = 1.4 vRd,c (first perimeter u1).
double vrd_max_punching_din(double vrdc);

// --------------------------------- shell reinforcement, Baumann method
// Sandwich panel forces of a layer (the plate resultants carried by the two
// outer layers at the lever arm z): N_i = n_i / 2 + m_i / z (bottom / tension
// layer of a sagging m; the caller flips the sign of m for the top layer).
double sandwich_panel_force(double n, double m, double z);
// Orthogonal reinforcement x/y for the membrane forces (Nx, Ny, Nxy) of a
// panel, Baumann (1972) with the minimum-reinforcement strut (theta = 45 deg,
// Nielsen; EN 1992-1-1 Annex F, eq. F.2/F.3 and F.6/F.7 for tension):
//   Nx >= -|Nxy| and Ny >= -|Nxy|: Zx = Nx + |Nxy|, Zy = Ny + |Nxy|;
//   Nx < -|Nxy|: Zx = 0, Zy = Ny + Nxy^2 / |Nx| (symmetrically for y).
// Z <= 0 means no reinforcement in that direction (returned as 0).
struct MembraneSteel {
  double Zx = 0, Zy = 0, Fc = 0;  // tie forces and concrete strut force per unit width
};
MembraneSteel membrane_reinforcement(double Nx, double Ny, double Nxy);
}  // namespace keta::ec2

// SIA 262:2013 punching (§4.3.6), the Swiss code's own formulas.
namespace keta::sia262 {
// eq. (3): tau_cd = 0.3 eta_t sqrt(fck) / gamma_c (MPa).
double tau_cd(double fck, double gamma_c = 1.5, double eta_t = 1.0);
// §4.3.6.2: control perimeter at d/2 of a rectangular column a x b:
// u = 2 a + 2 b + pi d; eq. (55): k_e = 1 / (1 + e_u / b_u), b_u = sqrt(4 A_c / pi)
// of the area A_c inside the perimeter (rounded corners), = 0.9 for interior
// columns (4.3.6.2.5) when the caller takes that simplification.
double perimeter_u(double a, double b, double d);
double area_inside_u(double a, double b, double d);
double ke(double e_u, double b_u);
// eq. (59), level of approximation III: psi = 1.2 (r_s / d) (f_sd / E_s)
// (m_sd / m_Rd)^1.5.
double psi_loa3(double rs, double d, double fsd, double Es, double msd, double mRd);
// eq. (58): k_r = 1 / (0.45 + 0.18 psi d k_g) (d in mm), k_g = 48 / (16 + D_max) (mm).
double kr(double psi, double d, double Dmax = 0.032);
// eq. (57): V_Rd,c = k_r tau_cd d u_red.
double vrd_c(double kr, double tau_cd, double d, double u_red);
}  // namespace keta::sia262

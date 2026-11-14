// EN 1994-1-1 composite beams — a VERIFICATION layer on top of the solver
// (like ec3.h): raw resistances and the intermediate quantities of each
// clause, every formula cited; the caller compares the utilisations with 1.
// SI units throughout (m, N, Pa); partial factors are the caller's.
//
//   * headed studs: §6.6.3.1 (eqs 6.18, 6.19, alpha 6.20/6.21) and the
//     reduction of profiled sheeting with ribs transverse to the beam,
//     §6.6.4.2 (eq 6.23, kt,max of Table 6.2);
//   * effective width of the concrete flange, §5.4.1.2 (eq 5.3);
//   * plastic resistance moment in sagging bending, §6.2.1.2 (rectangular
//     stress blocks, plastic neutral axis in the slab / the steel flange /
//     the web) and partial shear connection by the linear interpolation of
//     §6.2.1.3(5) (eq 6.1, ductile connectors);
//   * plastic resistance moment in HOGGING bending, §6.2.1.2(1)(a)-(c)
//     (slab reinforcement at fsd, concrete in tension neglected) with the
//     classes of the compressed web and bottom flange (§5.5.2, EN 1993-1-1
//     Table 5.2) that the plastic resistance requires (class 1 or 2);
//   * elastic properties of the transformed section (uncracked in sagging,
//     cracked in hogging), long-term modular ratios §5.4.2.2(2) eq (5.6), and
//     the SLS deflection of a simply supported beam, §7.3.1: construction
//     stage (steel alone under the wet concrete when unpropped), permanent
//     actions on the composite section (nL, psiL = 1.1), variable actions
//     (n0), primary effects of shrinkage (nL, psiL = 0.55, §5.4.2.2(1)/(6));
//   * minimum degree of shear connection, §6.6.1.2 (eqs 6.12-6.17);
//   * composite beams with large WEB OPENINGS by the method of SCI P355
//     (R. M. Lawson, S. J. Hicks, "Design of composite beams with large web
//     openings", SCI 2011) — the design guide of EN 1993-1-1 / EN 1994-1-1
//     practice: Tee properties of the perforated section (circular openings
//     by the equivalent rectangle 0.9 ho x 0.45 ho, §3.1.4), global bending
//     at the opening (eqs 5-11), shear of the Tees and of the slab (§3.3,
//     eqs 12-13), Vierendeel bending (§3.4, eqs 14-25, 61) and web-posts
//     between adjacent openings (§3.5, eqs 29-53).
//
// The rolled section enters with its real dimensions (h, b, tw, tf, r) and
// its tabulated A and Wpl (root fillets included); hw = h - 2 tf.
#pragma once
#include <utility>
#include <vector>

#include "keta/ec3.h"

namespace keta::ec4 {

// ------------------------------------------------------- shear connectors
// §6.6.3.1: PRd = min(0.8 fu pi d^2 / 4, 0.29 alpha d^2 sqrt(fck Ecm)) / gamma_v
// alpha = 0.2 (hsc/d + 1) for 3 <= hsc/d <= 4, 1 for hsc/d > 4 (eqs 6.20, 6.21).
double stud_prd_solid(double d, double fu, double hsc, double fck, double Ecm, double gamma_v);
// §6.6.4.2, sheeting with ribs transverse to the beam: kt = 0.7 / sqrt(nr)
// b0/hp (hsc/hp - 1) (eq 6.23) <= kt,max (Table 6.2: nr = 1: 0.85 / 1.0
// (t <= 1 / > 1 mm) through-deck welded d < 20, 0.75 / 0.75 with holes;
// nr = 2: 0.70 / 0.80, 0.60 / 0.60).
double kt_transverse(int nr, double b0, double hp, double hsc, double t_sheet, bool through_welded);

// ------------------------------------------------------- effective width
// §5.4.1.2 eq (5.3): beff = b0 + sum_i bei, bei = min(Le / 8, bi).
double beff(double b0, double Le, double b1, double b2);

// ------------------------------------------------------- sections
struct SteelI {
  double h = 0, b = 0, tw = 0, tf = 0, r = 0;  // real dimensions
  double A = 0, Wpl = 0;                        // tabulated (root fillets included)
  double I = 0;                                 // tabulated major-axis second moment (elastic checks)
  double fy = 0, gM0 = 1.0, gM1 = 1.0, E = 210e9;
  double hw() const { return h - 2 * tf; }
};
struct Slab {
  double hs = 0;  // total depth of the slab (sheeting included)
  double hp = 0;  // depth of the profiled sheeting (ribs transverse: no concrete in compression there)
  double fck = 0, gamma_c = 1.5;
  double hc() const { return hs - hp; }
  double fcd() const { return fck / gamma_c; }
};

// ------------------------------------------------- sagging bending §6.2.1
enum class PnaIn { slab, flange, web };
struct SaggingResult {
  double Ncf = 0;     // compressive resistance of the slab 0.85 fcd beff hc (§6.2.1.2(1))
  double Npla = 0;    // tensile resistance of the steel section A fy / gM0
  double Mpla = 0;    // plastic moment of the steel section Wpl fy / gM0
  PnaIn pna = PnaIn::slab;
  double x = 0;       // PNA depth: in the slab from its top, else in the steel from its top
  double Mpl = 0;     // full interaction
  double Nc = 0;      // compressive force the connectors transfer (<= Ncf)
  double eta = 1;     // degree of shear connection Nc / Ncf
  double MRd = 0;     // eq (6.1): Mpla + (Mpl - Mpla) eta
};
// Nc: the connectors' force between the critical sections (n PRd); Nc >= Ncf
// is full interaction.
SaggingResult sagging(const SteelI& s, const Slab& sl, double beff, double Nc);

// ------------------------------------------------- hogging bending §6.2.1
// Longitudinal reinforcement of the effective slab width: total area As,
// characteristic strength fsk, its centroid zs ABOVE the top face of the steel
// section. Es is taken equal to Ea (§3.2(2)).
struct Rebar {
  double As = 0, fsk = 0, gamma_s = 1.15, zs = 0;
  double fsd() const { return fsk / gamma_s; }
};
struct HoggingResult {
  double Ns = 0;      // As fsd, the reinforcement in tension (§6.2.1.2(1)(b))
  double Npla = 0, Mpla = 0;
  PnaIn pna = PnaIn::web;  // slab: the whole steel section in compression (Ns >= Npla)
  double x = 0;       // PNA depth from the top of the steel section (0 for `slab`)
  double Mpl = 0;     // plastic resistance moment (positive, top in tension)
  // classification of the compressed parts under the plastic distribution
  double alpha = 0;   // compressed fraction of the web's c = hw - 2 r
  double psi = 0;     // sigma_top / sigma_bottom of c in the cracked ELASTIC section (needs I)
  int web_class = 0, flange_class = 0;
  bool plastic_valid() const { return web_class <= 2 && flange_class <= 2; }
};
// The web case keeps the tabulated Wpl: the PNA rises from mid-depth by
// Ns / (2 tw fyd) (root fillets counted with the web, as in `sagging`). The
// web class needs s.I > 0 for psi (class 3); without it, a web beyond the
// class-2 limit reports 3. §5.5.2(3) ("hole in the web") is not applied.
HoggingResult hogging(const SteelI& s, const Rebar& r);

// ------------------------------------------- elastic section (SLS, §5.4.2)
// Transformed to steel: A, I about the elastic neutral axis, z of that axis
// ABOVE the centroid of the steel section (symmetric I: its mid-depth).
struct ElasticSection {
  double A = 0, I = 0, z = 0;
};
// Uncracked (sagging): the concrete of the slab above the ribs (hc x beff)
// with modular ratio n, its centroid at hp + hc / 2 above the steel; concrete
// in tension is NOT discounted (EN 1994 "uncracked" analysis); optional
// reinforcement (Es = Ea).
ElasticSection elastic_uncracked(const SteelI& s, const Slab& sl, double beff, double n,
                                 const Rebar* r = nullptr);
// Cracked (hogging): steel section + reinforcement.
ElasticSection elastic_cracked(const SteelI& s, const Rebar& r);
// §5.4.2.2(2) eq (5.6): nL = n0 (1 + psiL phi_t); psiL = 1.1 permanent,
// 0.55 shrinkage (primary and secondary effects), 1.5 imposed deformations.
double n_long_term(double n0, double psiL, double phi_t);

// ------------------------------------- SLS deflection, simply supported span
struct DeflectionInput {
  double L = 0;          // span
  double beff = 0;       // effective width (sagging, mid-span)
  double g_a = 0;        // steel self-weight [N/m]
  double g_c = 0;        // wet concrete + sheeting (construction stage) [N/m]
  double g_s = 0;        // superimposed permanent on the composite section [N/m]
  double q = 0;          // variable action (SLS combination) [N/m]
  bool propped = false;  // propped during casting: g_a + g_c go to the composite section
  double n0 = 0;         // Ea / Ecm (§5.4.2.2(1))
  double phi_t = 0;      // creep coefficient phi(t, t0) for permanent loads (0 = short term only)
  double phi_sh = 0;     // creep coefficient for shrinkage (age of 1 day at its start, §5.4.2.2(4))
  double eps_cs = 0;     // free shrinkage strain, magnitude (>= 0)
};
struct DeflectionResult {
  double nL = 0, n_sh = 0;
  ElasticSection perm, var, sh;   // transformed sections (nL, n0, n_sh)
  double d_construction = 0;      // steel alone (unpropped), 5 (g_a + g_c) L^4 / 384 Ea Ia
  double d_perm = 0;              // composite, nL (+ g_a + g_c when propped)
  double d_var = 0;               // composite, n0
  double Nsh = 0, e_sh = 0, Msh = 0;  // primary shrinkage force eps Ea / n_sh Ac at the slab centroid
  double d_sh = 0;                // Msh L^2 / (8 Ea I_sh)
  double d_total = 0;
};
DeflectionResult deflection_simply_supported(const SteelI& s, const Slab& sl, const DeflectionInput& in);

// ------------------------------------ minimum degree of shear connection
// §6.6.1.2(1): Le the length in sagging between zero-moment points [m], fy
// the nominal yield strength. Equal flanges (Abf / Atf = 1): eqs (6.12),
// (6.13); bottom flange 3 x top flange: eqs (6.14), (6.15); between the two,
// linear interpolation (§6.6.1.2(1)); Abf / Atf > 3 is outside the clause
// (rejected). eta >= 0.4 in every case.
double eta_min(double fy, double Le, double Abf_over_Atf = 1.0);
// §6.6.1.2(3), eqs (6.16), (6.17): headed studs d = 19 mm, hsc >= 76 mm after
// welding; rolled or welded I/H section with equal flanges; profiled sheeting
// continuous across the beam, ribs transverse, b0 / hp >= 2, hp <= 60 mm;
// one stud per rib; Nc from the linear interpolation of §6.2.1.3(5). The
// caller asserts the conditions; the geometric ones it passes are checked.
double eta_min_sheeting(double fy, double Le, double d, double hsc, double b0, double hp, int studs_per_rib);

// ------------------------------- continuous beams, global analysis §5.4.2.3
// (3): cracked regions of 15 % of the span each side of an internal support
// for spans with Lmin / Lmax >= 0.6 (the simplified method); rejected
// beyond, where the cracked analysis must be run on the stresses (5.4.2.3(2)).
double cracked_length(double L_left, double L_right);

// --------------------------------- lateral-torsional buckling §6.4.2/6.4.3
// Inverted U-frame stiffness per unit length (6.8): ks = k1 k2 / (k1 + k2),
// k1 = alpha Ea I2 / a (6.9: alpha 2 edge beam, 3 inner beam, 4 inner beam of
// four or more similar ones; I2 the cracked flexural stiffness of the slab per
// unit width, a the beam spacing), k2 = Ea tw^3 / (4 (1 - nu^2) hs) (6.10,
// unencased web, hs between the flange shear centres).
double uframe_ks(double alpha, double Ea, double I2, double a, double tw, double hs, double nu = 0.3);
// Table 6.1: the largest depth h [m] of a rolled IPE (or similar) / HE
// section, unencased or partially encased, that needs no LTB check by
// §6.4.3 (S235, S275, S355, S420 / S460). fy outside these grades: rejected.
double table61_max_depth(double fy, bool he, bool encased);

// --------------------------------------------- composite columns §6.7
// Rectangle-based cross-section about its y axis (z up): steel plates, a
// concrete core, reinforcement bars, each a rectangle (b x h centred at z)
// of a design strength (steel fyd, bars fsd, concrete 0.85 fcd or fcd for a
// filled hollow section, §6.7.3.2(1)); concrete takes no tension.
struct ColumnRect {
  double b = 0, h = 0, z = 0;
  enum Kind { Steel, Concrete, Bar } kind = Steel;
};
struct CompositeColumn {
  std::vector<ColumnRect> parts;
  double fyd = 0, fsd = 0, fcd_eff = 0;  // fcd_eff: 0.85 fcd (encased) or fcd (filled)
  double Ea = 210e9, Es = 210e9, Ecm = 33e9;
};
// (6.30) Npl,Rd = Aa fyd + Ac fcd_eff + As fsd.
double npl_rd(const CompositeColumn& c);
// the rigid-plastic resistance with the neutral axis at z0 (compression
// above): (N, M) about the section's y axis (N < 0 compression; M > 0
// compressing the top); exact rectangle integrals
std::pair<double, double> plastic_nm(const CompositeColumn& c, double z0);
// the plastic M at axial force N (bisection on the neutral axis, exact)
double mpl_n(const CompositeColumn& c, double N);
// (6.40) / (6.42): (EI)eff = Ke,II (Ea Ia + Es Is + Ke Ecm,eff Ic); Ke = 0.6,
// Ke,II = 1 (buckling) or 0.9 / Ke,II = 0.5 (second-order analysis);
// Ecm,eff = Ecm / (1 + (NG,Ed / NEd) phi_t) (6.41)
double ei_eff(const CompositeColumn& c, double Ke, double Ke2, double phi_t = 0.0, double NG_over_N = 0.0);
struct ColumnCheck {
  double Npl_Rk = 0, Npl_Rd = 0, Ncr = 0, lambda = 0, chi = 0;  // (6.39), Table 6.5 curve
  double MEd2 = 0;       // first-order moment amplified by k = beta / (1 - NEd / Ncr,eff) (6.44) and e0 NEd
  double mu_d = 0, MRd = 0, util_N = 0, util_M = 0;  // 6.7.3.6: MEd <= alpha_M mu_d Mpl,Rd
};
// Member check by the simplified method of §6.7.3 (doubly symmetric section,
// end moments ratio r for beta = max(0.66 + 0.44 r, 0.44) (Table 6.4; beta =
// 1 for a transverse load: pass r = 1), member imperfection e0 = L / d0 of
// Table 6.5 (d0 = 300 curve a, 200 b, 150 c), alpha_M = 0.9 (S235-S355) /
// 0.8 (S420 / S460).
ColumnCheck column_check(const CompositeColumn& c, double L_cr, double NEd, double MEd1, double r_end,
                         ec3::Curve curve, double d0, double alpha_M, double phi_t = 0.0, double NG_over_N = 0.0);

// ----------------------------------------- SCI P355 web openings (§3)
struct Opening {
  bool circular = false;
  double ho = 0;   // depth (diameter)
  double lo = 0;   // length (rectangular)
  double eo = 0;   // eccentricity of the centre above the mid-depth of the web
  // horizontal stiffeners (one pair, both Tees), §3.4.5
  bool stiffened = false;
  double tr = 0, br = 0, er = 0, fyr = 0;  // thickness, width, distance from the opening edge
  // ko of the local composite action: eq (61) (35) for a stiffened opening,
  // or eq (25) (25) as the worked example of P355 Appendix A applies it
  bool ko_eq61 = true;
  double x = 0;    // abscissa of the opening centre from the nearest support
};
struct Connection {
  double PRd = 0;     // design resistance of one stud
  int nr = 1;         // studs per rib
  double s = 0;       // rib (stud) spacing along the beam
  bool floor_studs = false;  // count studs over a length by whole ribs or pro rata (P355)
};
struct BeamData {
  SteelI steel;
  Slab slab;
  double Le = 0;            // span (sagging: the effective length)
  double bl = 0, br = 0;    // half distances to the adjacent beams (effective width bounds)
  double dl = 0, sl = 0;    // slab reinforcement: bar diameter and spacing (slab shear, §3.3.2)
  Connection conn;
};
struct TeeProps {
  double hw = 0;  // web depth of the Tee
  double A = 0;   // area (stiffener excluded)
  double z = 0;   // depth of its elastic neutral axis from the outer face of the flange
  double Ar = 0;  // stiffener area (effective, 10 eps tr), 0 without
  double NRd = 0; // A fy / gM0 + Ar fyr / gM0
};
struct OpeningResult {
  double ho_eq = 0, lo_eq = 0;   // equivalent rectangle (§3.1.4)
  TeeProps top, bottom;
  double heff = 0;               // distance between the Tees' centroids
  // global bending at the opening centre
  double beff_o = 0, NcRd1 = 0, NcRd2 = 0, NcRd = 0, nsc = 0;
  bool pna_in_top_tee = true;
  double NbTEd = 0, AtTcap = 0, NtTRd = 0, MoRd = 0;
  // shear (§3.3)
  double AvT = 0, VtRd = 0, VbRd = 0;  // Tee shear areas (the real opening depth)
  double sigma_cp = 0, VcRd_a = 0, VcRd_min = 0, VcRd = 0, VRd = 0;
  // Vierendeel bending (§3.4)
  double lo_eff = 0, mu = 0;
  double hwt_v = 0, hwb_v = 0;        // Tee webs after the Table 3.2 class limit
  double MtTVRd = 0, MbTVRd = 0;      // eq (19)
  double NtEd = 0, MtTNVRd = 0, MbTNVRd = 0;  // eq (20) (linear for a stiffened Tee)
  double nsco = 0, dNcRd = 0, ko = 1, MvcRd = 0;  // eqs (24), (25) / (61)
  double MVbRd = 0;                   // eq (14): 2 MtT + 2 MbT + Mvc
  double MVbEd = 0;                   // demand |V| lo_eq at the opening edge
};
// MEd, VEd at the opening centre; VEd_edge the shear at its support-side
// edge (x - ho/2 circular, x - lo/2 rectangular).
OpeningResult opening(const BeamData& b, const Opening& o, double MEd, double VEd, double VEd_edge);

// web-post between two adjacent openings (§3.5): s0 the post width, s the
// centre distance, ho their depth; heff, zt and MbTNVRd those of the
// adjacent Tees, nscs the studs over s.
struct WebPostInput {
  double s0 = 0, s = 0, ho = 0, eo = 0;
  double heff = 0, zt = 0, MbTNVRd = 0, lo = 0;
  double VEd = 0;   // vertical shear at the web-post
  double nscs = 0;  // studs over s
};
struct WebPostResult {
  double VwpEd = 0;    // eq (29) longitudinal shear
  double dNcsRd = 0;   // eq (30)
  double VwpRd = 0;    // eq (36)
  double MwpRd = 0;    // eq (44)
  double lambda = 0, chi = 0, NwpRd = 0;  // eqs (41), (42), buckling curve c
  double VRdM = 0, VRdN = 0;  // eqs (51), (53): the vertical shear the post allows
};
WebPostResult web_post(const BeamData& b, const WebPostInput& in);

}  // namespace keta::ec4

// Unified aging concrete material laws — the compute side of the ST1 `MAT`
// command (ST1 v24 manual §3.3, pp. 74-88; formulas from annexes B3.3 to B3.9,
// pp. 337-371). One struct, `MatLaw`, answers the three time-dependent
// questions a phasing/prestress driver asks a concrete:
//
//     E(t)          Young's modulus at age t          [Pa]
//     phi(t, t0)    creep coefficient                 [-]
//     eps_sh(t, ts) shrinkage strain, NEGATIVE = shortening   [-]
//
// The laws (`MatLaw::Code`) and their authoritative sources:
//
//  * Direct    — no law: a constant modulus (ST1 `CONS E ...` / `MAT E value`).
//  * BPEL      — BPEL 91 chap. 2.1 (modulus) + annexe 1 §3 (shrinkage) and
//                §4.2 (creep). Validity fc28 <= 60 MPa (manual p. 74/75).
//  * BHP1/BHP2 — AFREM/BHP proposal of 16 May 1995 ("BHP99" in the annexes),
//                without (BHP1) / with (BHP2) silica fume. 40 <= fc28 <= 80 MPa.
//  * BPEL_BL   — BPEL 99 annexe 5 lightweight concrete (modulus, shrinkage and
//                creep variants of the BPEL law).
//  * CEB       — CEB-FIP Model Code 1990 (bulletin 213/214), eq. 2.1-15/54/57
//                (modulus), 2.1-65..71 (creep), 2.1-74..79 (shrinkage).
//  * EC2       — EN 1992-1-1 Table 3.1 + §3.1.4 + Annex B. This law DELEGATES
//                to the existing keta `creep` module (`ec2_ecm`, `ec2_phi`,
//                `ec2_shrinkage`) — it is not re-derived here.
//  * EC2_BHP1/2— EN 1992-2 Annex B (B.103 shrinkage, B.104-B.116 HP creep),
//                without / with silica fume. 40 <= fck <= 90 MPa.
//  * FIB2010   — fib Model Code 2010 eq. 5.1-20..24, 5.1-56/57 (modulus),
//                5.1-63..73 (creep, basic + drying split), 5.1-75..83
//                (shrinkage, autogenous + drying split).
//  * ACI_B3    — Bazant-Baweja model B3 as published in ACI 209.2R-21:
//                compliance q1..q5 (eq. A-37..A-50), shrinkage eq. A-31..A-36.
//  * IP2, MC78, COURBON — the French legacy laws:
//      IP2     — Instruction Provisoire n°2 (1973): E28 = 76650 sqrt(sigma28)
//                [T/m2], Ej = E28 sqrt(0.685 log10(j + 1)) below 28 d;
//                shrinkage kb eps_c(H) kc1(hm) kD r(t'), t' = t sqrt(1/(10 hm));
//                creep J_def = kb kc(H) kd(t'0) ker(2 hm) r(t - t0) / Ei(t0),
//                t'0 = t0 (theta + 10)/30 (phi = J_def E28).
//      MC78    — CEB-FIP Model Code 1978/80: E28 = 204671 (fck + 800)^(1/3)
//                [T/m2], aging (fcm(j')/fcm(28))^(1/3) to 400 fictitious
//                days; shrinkage k1 eps_s1(H) eps_s2(h0) beta_s(t'); creep
//                0.4 beta_d(t' - t'0) + phi_f (beta_f(t') - beta_f(t'0)),
//                phi_f = k1 phi_f1(H) phi_f2, t' = alpha t (theta + 10)/30,
//                h0 = 2 lambda(H) hm (MC78's lambda 2 Ac/u). phi_f2 is read on
//                MC78 fig. e.2, which no accessible source
//                tabulates: it MUST be given (`phi_f2`, no default). The fit of
//                beta_s (MC78 fig. e.4) is valid for 0.05 <= h0 < 0.99 m: its
//                kinetics are evaluated with h0 clamped to that range (the
//                magnitude eps_s2 takes the true h0). No reinforcement factor
//                (the source's MC78 has none).
//      COURBON — J. Courbon's non-aging law: phi = 2 (1 -
//                exp(-0.006 (t - t0))), constant modulus `e_direct`, no
//                shrinkage.
//              The source supports neither a default phi_f2 table nor a
//              steel factor on the MC78 shrinkage: keta has none.
//
// UNITS. keta is SI: `fc` and every modulus are in Pa, ages in DAYS, RH in %.
// The mean radius keeps the two vocabularies the laws themselves use:
// `rm_cm` = Ac/u in CENTIMETRES (BPEL family) and `h_mm` = 2*Ac/u in
// MILLIMETRES (CEB/EC2/fib/B3) — they describe the same section and satisfy
// h_mm = 20 * rm_cm; each law reads only its own. An ST1 file gives the
// characteristic strength in MPa, t/m^2 or kN/m^2 (`MPA`/`TM2`/`KNM2`
// selectors): `strength_to_pa()` converts, with 1 t/m^2 = 9806.65 Pa
// (g = 9.80665, keta's kGravity) and 1 kN/m^2 = 1000 Pa. NOTE: ST1 itself uses
// g = 9.81 for that conversion (k_unite = 1000/9.81); keta uses the exact
// standard gravity — a 5.5e-5 relative difference on a TM2/KNM2 input.
//
// WHAT phi() IS REFERENCED TO. Every code writes the creep strain as
// eps_cc = phi * sigma / E_ref, but E_ref differs from code to code, so phi
// values are NOT interchangeable between laws. `E_ref(t0)` returns it:
//     BPEL / BHP / BPEL_BL : Ei28, the 28-day instantaneous modulus
//     CEB                  : Eci, the 28-day modulus (MC90 2.1-15)
//     EC2 / EC2_BHP        : Ec = 1.05 Ecm (EN 1992-1-1 §3.1.4(2), eq. 3.6)
//     FIB2010              : Eci, the 28-day tangent modulus
//     ACI_B3               : 1/(q1 + Delta_0.01), the conventional instantaneous
//                            compliance — the only t0-dependent one.
//
// AVANCE MULTIPLIERS (manual §3.3.7 p. 88) are applied per component. Where a
// law has no endogenous/drying split, the manual's own formula shows which
// single coefficient it carries, and the others are documented no-ops:
//
//   law        E    shrinkage              creep                    kinetics
//   BPEL/BL    km   kre (whole, annex      kfe (whole, annex        kfc unused
//                   B3.3 "eps_r = kre ks   B3.3 "eps_fl = kfe ...")
//                   eps_0"); krd unused
//   BHP1/2     km   krd*eps_rd + kre*eps_re  kfd*eps_fd + kfe*eps_fp  kfc active
//   CEB        km   kre (whole, B3.5)      kfe (whole, B3.5)        kfc unused
//   EC2        km   krd*eps_cd + kre*eps_ca  kfe (whole, eq. 3.6)   kfc unused
//   EC2_BHP1/2 km   krd*eps_cd + kre*eps_ca  kfd*phi_d + kfe*phi_b  kfc active
//   FIB2010    km   krd*eps_cd + kre*eps_ca  kfd*phi_dc + kfe*phi_bc kfc active
//   ACI_B3     km   krd (drying only; B3    kfd*Cd + kfe*C0          kfc unused
//                   has no autogenous term
//                   in ST1 — kre unused)
//
// BONDED-STEEL REDUCTION. BHP99, EC2-2 HP, fib 2010 and B3 divide the total
// deferred strain by (1 + n*rho_s) (rho_s = As/B, ST1 `ROS`). n is fixed by the
// BHP99 annex (15 for 40 <= fc28 < 60, 9 for 60 <= fc28 < 80) and is Es/Ec(28)
// for the others (Es = `es`, default 200 GPa). BPEL instead carries the same
// idea inside ks = 1/(1 + 20 rho_s); CEB and EC2 have no such term in ST1.
//
// KNOWN DEVIATIONS / LIMITS (never silent):
//  * Thermal treatment (`T_THERM`, the t_eb equivalent-age shift of BPEL
//    annexe 6 / EC2 §10.3) is NOT implemented — it only ever applies inside a
//    construction phasing, and ST1 itself asks the user to enter the extra
//    losses by hand outside phasing. No parameter of it is accepted here.
//  * BPEL/BHP fcj is used as written by the annex, WITHOUT the "fcj = fc28 for
//    j > 28" cap: at 28 days the ordinary branch is exactly fc28
//    (4.76 + 0.83*28 = 28.0) while the high-strength branch overshoots by 3.7 %,
//    and both keep growing beyond — this is what ST1 prints (the annex states
//    the j > 28 cap only for the thermal-treatment variant).
//  * BPEL_BL modulus: the annex line reads "Eij = km k_unite 1650 cbrt(fcj(rho_bs))"
//    in an OCR-damaged form. keta reads it as Eij = 1650 * rho_bs^2 * fcj^(1/3)
//    with rho_bs in t/m^3 — the only reading of that line that degenerates to
//    the normal-weight 11000 * fcj^(1/3) family at rho_bs = 2.4 t/m^3 (it gives
//    9504 there, i.e. 86 % of 11000, the usual lightweight abatement). Flagged
//    as the one law here whose constant is reconstructed rather than read.
//  * MC90's loading-age adjustment (eq. 2.1-72) is applied for CEB although
//    ST1's annex B3.5 does not restate it; at t0 = 28 d with an N/R cement it
//    is the identity.
//  * ACI_B3 q5 (eq. A-49) needs |eps_sh_inf| in units of 1e-6, not as a bare
//    strain — with a bare strain the drying compliance comes out ~4000x too
//    large. keta scales explicitly; see `b3_q5` in matlaws.cpp.
//  * Relaxation is not re-derived: `relaxation_ratio_seuil` reuses the EC2
//    `ec2_relaxation_ratio` of keta/creep.h and only adds ST1's `SEUIL`
//    threshold semantics (manual p. 89).
#pragma once
#include <memory>
#include <optional>

#include "keta/creep.h"

namespace keta {

// ---------------- units ----------------

// Unit selector of an ST1 strength (`MPA` / `TM2` / `KNM2` after `E <law>`).
enum class StrengthUnit { MPa, TM2, KNM2 };

// Convert a characteristic strength expressed in one of ST1's units to Pa.
// 1 t/m^2 = 9806.65 Pa (kGravity), 1 kN/m^2 = 1000 Pa.
double strength_to_pa(double value, StrengthUnit unit);
// Same, to MPa (the unit every law formula below is written in).
double strength_to_mpa(double value, StrengthUnit unit);

// ---------------- cement ----------------

// Cement class. The three ST1 vocabularies map onto one enum:
//   EC2 / fib / EC2_BHP : S, N, R          (CEM 32.5N / 32.5R-42.5N / 42.5R+)
//   CEB90               : SL == S, N, R, RS
// The class drives the beta_cc exponent s, the loading-age adjustment exponent
// alpha, and the shrinkage coefficients (alpha_ds1/2, beta_sc, alpha_as).
enum class CementClass { S, N, R, RS };

// ---------------- advanced multipliers (MAT ... AVANCE ...) ----------------

struct AvanceFactors {
  double k_module = 1.0;               // K_MODULE
  double k_retrait_endogene = 1.0;     // K_RETRAIT_ENDOGENE
  double k_retrait_dessication = 1.0;  // K_RETRAIT_DESSICATION
  double k_fluage_endogene = 1.0;      // K_FLUAGE_ENDOGENE (basic / own creep)
  double k_fluage_dessication = 1.0;   // K_FLUAGE_DESSICATION
  double k_fluage_cinetique = 1.0;     // K_FLUAGE_CINETIQUE (creep kinetics)
  // EN 1992-1-1 (2027) kinetics (annex B4.10.1): they divide the age in the
  // time functions of their term.
  double k_retrait_endogene_cin = 1.0;     // K_RETRAIT_ENDOGENE_CINETIQUE (krec)
  double k_retrait_dessication_cin = 1.0;  // K_RETRAIT_DESSIC(ATION)_CINETIQUE (krdc)
  double k_fluage_endogene_cin = 1.0;      // K_FLUAGE_ENDOGENE_CINETIQUE (kfec)
  double k_fluage_dessication_cin = 1.0;   // K_FLUAGE_DESSIC(ATION)_CINETIQUE (kfdc)
};

// ---------------- the law ----------------

struct MatLaw {
  enum class Code {
    Direct,     // constant modulus, no creep, no shrinkage
    BPEL,       // BPEL 91, ordinary concrete (fc28 <= 60 MPa)
    BHP1,       // AFREM/BHP 1995 without silica fume (40..80 MPa)
    BHP2,       // AFREM/BHP 1995 with silica fume    (40..80 MPa)
    BPEL_BL,    // BPEL 99 annexe 5, lightweight
    CEB,        // CEB-FIP MC90
    EC2,        // EN 1992-1-1 (delegates to keta/creep.h)
    EC2_BHP1,   // EN 1992-2 Annex B without silica fume (40..90 MPa)
    EC2_BHP2,   // EN 1992-2 Annex B with silica fume    (40..90 MPa)
    FIB2010,    // fib Model Code 2010
    ACI_B3,     // Bazant B3 per ACI 209.2R-21
    EC2_2027,   // EN 1992-1-1:2023 annex B (ST1 annex B4.10.1)
    IP2,        // Instruction Provisoire n°2 (1973) (appended)
    MC78,       // CEB-FIP Model Code 1978/80 (appended)
    COURBON     // J. Courbon's non-aging creep (appended)
  };

  Code code = Code::EC2;

  // Strength: fc28 for the BPEL/BHP family, fck for CEB/EC2/fib/B3 [Pa].
  double fc = 35e6;
  // Modulus for Code::Direct [Pa].
  double e_direct = 35e9;

  // Environment and geometry.
  double rh = 70.0;      // ambient relative humidity [%]        (ST1 ROH / RH)
  double rm_cm = 12.5;   // BPEL mean radius Ac/u [cm]           (ST1 RM)
  double h_mm = 250.0;   // notional size 2Ac/u [mm]             (ST1 H / RM_EC)
  double ros = 0.0;      // bonded steel ratio As/B [-]          (ST1 ROS)
  double ts_dry = 3.0;   // age at start of drying [days]; for ACI_B3 this is
                         // the end-of-cure age jc (ST1 TEMPS_CURE)
  CementClass cement = CementClass::N;
  double es = 200e9;     // reinforcement modulus for n = Es/Ec [Pa]

  // Explicit final shrinkage (ST1 `RETRAIT <law> ER er`), a positive magnitude.
  // What it replaces follows the manual's own wording per family:
  //   BPEL / BPEL_BL : the TOTAL final shrinkage eps_r    (p. 75)
  //   BHP1/2, EC2, EC2_BHP1/2, FIB2010, ACI_B3 : the final DRYING value only —
  //     "le terme de retrait endogene est toujours calcule par ST1"
  //     (pp. 80, 82, 84, 86)
  //   CEB            : the final total shrinkage eps_cs0  (p. 78)
  std::optional<double> er;

  // Lightweight concrete (BPEL_BL).
  double robs = 1.8;  // dry density [t/m^3]                     (ST1 ROBS)
  double ea = 0.0;    // water absorbed by the aggregate [l/m^3] (ST1 EA)

  // fib 2010: apply the secant reduction alpha_i (eq. 5.1-24). ST1 applies the
  // tangent modulus (alpha_i = 1) — keep false to match it.
  bool fib_secant = false;

  // ACI_B3 composition and curing.
  int b3_cement_type = 1;    // 1 = type I, 2 = type II, 3 = type III (CIMENT)
  double k_forme = 1.0;      // K_FORME ks: 1.00 infinite slab, 1.15 infinite
                             // cylinder, 1.25 infinite prism, 1.30 sphere,
                             // 1.55 cube (ACI 209.2R-21 tab. A-9)
  double k_cure = 1.0;       // K_CURE: 0.75 evaporation, 1.20 protected,
                             // 1.00 immersed (tab. A-8)
  double water = 205.0;      // QUANTITE_EAU w [kg/m^3], 0.35 <= w/c <= 0.85
  double cement_qty = 409.0; // QUANTITE_CIMENT c [kg/m^3], 160 <= c <= 720
  double aggregate = 1730.0; // QUANTITE_AGREGAT a [kg/m^3], 2.5 <= a/c <= 13.5

  AvanceFactors avance;
  // EC2_2027: silica fume (ST1 SILICE OUI) lowers the creep exponents.
  bool silica = false;
  // IP2 / MC78 / COURBON (the mean radius hm = S/p is `rm_cm`):
  double theta = 20.0;          // ambient temperature [deg C] of the fictitious ages
  double kb = 1.0;              // IP2 concrete coefficient kb (0.6..2; ip2_kb_from_mix)
  double k1_consistency = 1.0;  // MC78 k1: 1.25 soft, 1 normal, 0.75 firm
  double phi_f2 = 0.0;          // MC78 phi_f2(h0) of fig. e.2 — REQUIRED for MC78
  // IP2 kb from the mix: (w/c)(7.9e-3 c + 0.766) - (9.4e-4 c + 0.397)
  static double ip2_kb_from_mix(double cement_kg_m3, double water_kg_m3);

  // Typed rejection of an out-of-range / non-finite parameter set. Called at
  // the head of every query below, so a bad law can never return a number.
  void validate() const;

  // Young's modulus at age t [days] -> Pa. Includes K_MODULE.
  double E(double t_days) const;
  // BPEL/BHP deferred ("differed") modulus Ev = 3700 fcj^(1/3) MPa
  // (BPEL 91 §2.1.3); typed error for any other law.
  double Ev(double t_days) const;
  // Modulus phi() is referenced to (see the header note). t0 [days] matters
  // only for ACI_B3.
  double E_ref(double t0 = 28.0) const;

  // Mean strength at age t [days] -> Pa: fcj for the BPEL/BHP family,
  // fcm(t) = beta_cc(t) * fcm for the others.
  double fc_mean(double t_days) const;

  // Creep coefficient phi(t, t0), ages in days, t >= t0 > 0.
  double phi(double t, double t0) const;

  // Shrinkage strain at age t with drying starting at ts (end of cure), ages in
  // days. NEGATIVE = shortening (the keta/creep.h convention). The DRYING part
  // is counted from ts, the AUTOGENOUS part is the law's own function of the
  // age since casting — exactly what EN 1992-1-1 eq. 3.8 does; laws without an
  // autogenous term (BPEL, CEB, ACI_B3) therefore return 0 for t <= ts.
  double eps_sh(double t, double ts) const;

  // Steel/concrete modular ratio n used by the (1 + n rho_s) reductions.
  double n_equiv() const;
};

// ---------------- composition ----------------

// The phasing law of a concrete composed of independently chosen components
// (the C++ face of ST1's MAT ... FLUAGE <loi> RETRAIT <loi>): creep phi and
// the modulus it is referenced to from `creep`, the shrinkage from
// `shrinkage` (drying counted from age `ts`); `stress_scale` = Pa per model
// stress unit (1 in SI). The instantaneous modulus of the increments is a
// third choice, PhasingOptions::modulus_at (e.g. a MatLaw's E). The laws are
// validated here, once.
struct ConcreteLaw;
// `nonlinear`: the EN 1992-1-1 §3.1.4(4) non-linear creep under high
// compression (ST1 LINEAIRE NON); off by default, as ST1's LINEAIRE OUI.
std::shared_ptr<ConcreteLaw> concrete_law_from(const MatLaw& creep, const MatLaw& shrinkage, double ts = 0.0,
                                               double stress_scale = 1.0, bool nonlinear = false);

// ---------------- non-linear creep (ST1 FLUAGE ... LINEAIRE NON) ----------------

// EN 1992-1-1 §3.1.4(4): under a high compressive stress at loading the creep
// coefficient is multiplied by exp(1.5 (k_sigma - 0.45)) with
// k_sigma = max(0.45, |sigma_c| / fcm(t0)) — the manual states exactly this
// correction for EC2 and EC2_BHP (pp. 81 and 83) and it is OFF by default
// there ("LINEAIRE OUI"). `sigma_c` and `fcm_t0` in Pa; the factor is >= 1.
//
// ST1's other creep switch, `TRACTION OUI/NON`, is a driver-level policy (may
// creep develop under tension at all?), not a parameter of any law here — it
// belongs to whoever applies phi() to a stress state.
double nonlinear_creep_factor(double sigma_c, double fcm_t0);

// ---------------- relaxation threshold (ST1 PREC ... SEUIL) ----------------

// ST1 `SEUIL seuil_relax` (manual p. 89): relaxation threshold date in HOURS.
// No value  -> no threshold, relaxation keeps growing;
// value 0   -> ST1's 5e5 h (57 years) threshold;
// any other -> that value.
// Returns the time to feed the relaxation law, i.e. min(t, threshold).
inline constexpr double kRelaxSeuilDefaultHours = 5e5;
double relax_seuil_hours(double t_hours, std::optional<double> seuil);

// Relaxation loss ratio with the SEUIL cap applied: the EC2 §3.3.2 law of
// keta/creep.h (`ec2_relaxation_ratio`) evaluated at `relax_seuil_hours`.
double relaxation_ratio_seuil(double t_hours, double sigma_pi, double fpk,
                              int relaxation_class, double rho1000,
                              std::optional<double> seuil);

}  // namespace keta

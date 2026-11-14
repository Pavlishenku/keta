// Serviceability values of a concrete section, on top of its service states
// (keta/section_losses.h, keta/section_service.h) — RAW VALUES: limits,
// ratios of a computed stress to its limit, cracking moment and minimum
// reinforcement, never a verdict.
//
// STRESS LIMITS (EN 1992-1-1 §7.2). The clause names five stresses:
//   k1 fck  concrete compression under the characteristic combination, in
//           areas of exposure classes XD, XF and XS (§7.2(2));
//   k2 fck  concrete compression under the quasi-permanent loads below which
//           creep may be taken linear (§7.2(3));
//   k3 fyk  tension of the reinforcement under the characteristic combination,
//           k4 fyk where an imposed deformation causes the stress (§7.2(5));
//   k5 fpk  mean value of the stress in the prestressing tendons (§7.2(5):
//           "should not exceed k5 fpk" in EN 1992-1-1:2004 and in its French
//           and German versions; the consolidated BSI reprints of 2011 and
//           2014 print "k5 fyk" there, a misprint no corrigendum orders).
// The recommended values are k1 = 0,6, k2 = 0,45, k3 = 0,8, k4 = 1 and
// k5 = 0,75; a National Annex may give others, the caller hands them in.
// Which combination a state was computed under is the caller's knowledge:
// stress_ratios divides the stresses of ONE state by the five limits and the
// caller reads the ratios that belong to its combination.
#pragma once
#include "keta/section_service.h"

namespace keta::ec2 {

struct StressLimitFactors {
  double k1 = 0.6, k2 = 0.45, k3 = 0.8, k4 = 1.0, k5 = 0.75;
};
// The five limits [Pa], all positive.
struct StressLimits {
  double concrete_characteristic = 0;   // k1 fck
  double concrete_quasi_permanent = 0;  // k2 fck
  double bar = 0;                       // k3 fyk
  double bar_imposed = 0;               // k4 fyk
  double tendon = 0;                    // k5 fpk
};
// Typed errors: a strength or a factor that is not finite and positive.
StressLimits stress_limits(double fck, double fyk, double fpk, const StressLimitFactors& k = {});

// The stresses of a service state that §7.2 limits, and each over its limit.
// `concrete_compression` is the compression of the most compressed vertex of
// the rings (-fibre_min.stress, 0 when nothing is compressed) — the largest
// compression of the concrete wherever its law rises with the shortening,
// which a softening law past its peak does not. `bar_tension` and
// `tendon_tension` are the largest tensile stress among the bars and among
// the tendons of the state (0 without any, or when none is stretched).
struct StressRatios {
  double concrete_compression = 0, bar_tension = 0, tendon_tension = 0;  // [Pa], >= 0
  double concrete_characteristic = 0;   // concrete_compression / (k1 fck)
  double concrete_quasi_permanent = 0;  // concrete_compression / (k2 fck)
  double bar = 0;                       // bar_tension / (k3 fyk)
  double bar_imposed = 0;               // bar_tension / (k4 fyk)
  double tendon = 0;                    // tendon_tension / (k5 fpk)
};
// Typed error: a limit that is not finite and positive.
StressRatios stress_ratios(const PrestressServiceState& state, const StressLimits& limits);

// CHOICE OF THE STATE (EN 1992-1-1 §7.1(2)). Stresses are computed on the
// uncracked section provided that the flexural tensile stress does not exceed
// fct,eff — fctm or fctm,fl, the value the minimum tension reinforcement is
// based on. service_state computes the uncracked state
// (prestress_service_state), reads its largest concrete stress
// `flexural_tension` (fibre_max, negative while every fibre is compressed)
// and returns that state while flexural_tension <= fct_eff, the cracked state
// (cracked_service_state) beyond. The choice is made on THIS action alone: a
// section cracked by an earlier, larger action stays cracked, which the
// caller knows and this function does not — it then calls the cracked state
// itself. Typed errors: fct_eff negative or not finite, and those of the two
// states.
struct ServiceState {
  PrestressServiceState state;
  bool cracked = false;
  double flexural_tension = 0;  // [Pa], of the uncracked section
};
ServiceState service_state(const RcPolySection& s, double Ecm, double fct_eff,
                           const Vec3& sustained = Vec3::Zero(), const std::vector<Vec2>& points = {},
                           const PassiveSteel& bars = {});

// CRACKING MOMENT. The moment M_cr >= 0 of direction theta — (My, Mz) =
// M_cr (cos theta, sin theta), as moment_resistance reads it — that, added to
// `sustained`, brings the first vertex of the rings of the uncracked section
// (prestress_service_state) to the tensile strength fct. With sigma_cp the
// stress of that fibre under the prestress and `sustained` (negative where
// they compress it) and W the modulus of the homogenised section at that
// fibre for this direction, M_cr = (fct - sigma_cp) W: the prestress shifts
// the cracking moment of the reinforced section, fct W, by -sigma_cp W.
// `fibre` is that vertex with its stress sigma_cp, `modulus` is W.
// Typed errors: fct negative or not finite, a fibre already past fct under
// `sustained` alone, and those of prestress_service_state.
struct CrackingMoment {
  double moment = 0;   // M_cr [N m]
  FibreStress fibre;   // the vertex that cracks first; stress = sigma_cp [Pa]
  double modulus = 0;  // W [m^3]
};
CrackingMoment cracking_moment(const RcPolySection& s, double Ecm, double fct, double theta,
                               const Vec3& sustained = Vec3::Zero(), const PassiveSteel& bars = {});
// The section with its tendons slack: every initial strain eps_p0 at zero. A
// bonded tendon stays in the homogenised section as a steel, an unbonded one
// carries nothing. Its cracking moment ignores any effect of prestressing —
// M_rep of EN 1992-2 §6.1(109) b). The initial strains of the bars, which are
// no prestress, stay.
RcPolySection without_prestress(const RcPolySection& s);

// MINIMUM REINFORCEMENT AGAINST BRITTLE FAILURE (EN 1992-2 §6.1(109) b) and
// (110)). Expr. (6.101a): As,min = M_rep / (z_s fyk), z_s the lever arm at
// the ultimate limit state related to the reinforcing steel. With tendons
// counted as fully active (pretensioned, below 0,6 fpk under the
// characteristic combination), Expr. (6.101b): As,min fyk + Ap dsigma_p >=
// M_rep / z, dsigma_p the smaller of 0,4 fptk and 500 MPa — give their area
// `Ap` and strength `fptk`; the result is not below zero.
// Typed errors: M_rep or Ap negative, z or fyk not positive, Ap > 0 without
// fptk > 0.
double brittle_failure_steel(double M_rep, double z, double fyk, double Ap = 0.0, double fptk = 0.0);
// The lever arm z_s on the ultimate plane of moment_resistance(s, omega, NEd,
// theta): the distance between the point where the resultant of the concrete
// compression acts and the point where the resultant of the forces of the
// stretched bars acts — the arm of their couple, M / (As fyd) when those two
// forces are the only ones. Compressed bars and tendons enter the ultimate
// plane and neither resultant. Typed errors: no ultimate plane carrying NEd
// with a moment along theta (NEd outside the axial capacity, or a direction
// the moments about the origin of the axes never take), no stretched bar or
// no compressed concrete on that plane.
double ultimate_lever_arm(const RcPolySection& s, double omega, double NEd, double theta);

// CRACK WIDTH (EN 1992-1-1 §7.3.4, corrigendum AC:2010 included).
//   wk = sr,max (eps_sm - eps_cm)                                        (7.8)
//   eps_sm - eps_cm = [sigma_s - kt fct,eff / rho_p,eff (1 + alpha_e
//                      rho_p,eff)] / Es  >=  0,6 sigma_s / Es            (7.9)
//   rho_p,eff = (As + xi1 Ap') / Ac,eff                                  (7.10)
//   xi1 = sqrt(xi phi_s / phi_p), sqrt(xi) with tendons alone            (7.5)
//   sr,max = k3 c + k1 k2 k4 phi / rho_p,eff                             (7.11)
//   sr,max = 1,3 (h - x) where the bars are farther apart than
//            5 (c + phi / 2) or no bonded steel lies in Ac,eff           (7.14)
// Ac,eff is the concrete within hc,ef = min(2,5 (h - d), (h - x) / 3, h / 2)
// of the most stretched face (§7.3.2(3), Fig. 7.1), Ap' the bonded tendons
// within it; with a compressed zone (x > 0) h / 2 never rules, (h - x) / 3
// being below h / 3. Expr. (7.10) carries xi1 since the corrigendum
// EN 1992-1-1:2004/AC:2010 (modification 99; the official Czech version of
// the corrigendum and the Spanish Codigo Estructural of 2021 print it so).
// The text before it squared xi1, and the German version and most design
// software still do: `xi1_squared` gives that form.
struct CrackWidthParams {
  double fct_eff = 0;      // [Pa], > 0 (fctm for crack widths, §7.1(2))
  double cover = 0.030;    // c, cover to the longitudinal reinforcement [m]
  double bar_dia = 0.016;  // phi, (equivalent, Expr. (7.12)) bar diameter [m]; the largest, phi_s of (7.5)
  double spacing = 0;      // bar spacing [m]; > 5 (c + phi / 2) switches to Expr. (7.14)
  double k1 = 0.8;         // bond: 0,8 high bond bars, 1,6 plain surface
  double k2 = 0.5;         // strain distribution: 0,5 bending, 1,0 pure tension
  double k3 = 3.4, k4 = 0.425;
  double kt = 0.6;         // 0,6 short term, 0,4 long term loading
  double alpha_e = 0;      // Es / Ecm of Expr. (7.9); 0: Es over the modulus of the cracked state
  // bonded tendons: xi the ratio of bond strength of Table 6.2 (0: tendons
  // left out of rho_p,eff), tendon_dia the equivalent diameter phi_p of §6.8.2
  double xi = 0, tendon_dia = 0;
  bool xi1_squared = false;  // Expr. (7.10) as printed before AC:2010
};
struct CrackWidth {
  // geometry of the cracked state along its strain gradient, depths from the
  // most compressed fibre (zero when the chain is called on its own)
  double h = 0, d = 0, x = 0;
  double hc_ef = 0, Ac_eff = 0;  // [m], [m^2]
  double As = 0, Ap = 0;         // bars and bonded tendons within Ac,eff [m^2]
  double xi1 = 0, rho_p_eff = 0;
  double sigma_s = 0;         // [Pa]
  double d_eps = 0;           // eps_sm - eps_cm
  double s_r_max = 0;         // [m]
  bool wide_spacing = false;  // Expr. (7.14) was used
  double wk = 0;              // [m]
};
// Expr. (7.5) and (7.8) to (7.14) from their terms: the steel stress sigma_s
// of the cracked section, the bars As and bonded tendons Ap within Ac,eff,
// the stretched depth h - x (read by Expr. (7.14) alone) and the modulus Es
// of the bars. `alpha_e` of the parameters must be given here. Typed errors:
// sigma_s, Ac_eff, Es, fct_eff, alpha_e or bar_dia not positive, a negative
// area, cover or factor, tendons counted without their diameter.
CrackWidth crack_width_of(double sigma_s, double As, double Ap, double Ac_eff, double h_minus_x, double Es,
                          const CrackWidthParams& p);
// The crack width of the section under `action`, from its cracked state
// (cracked_service_state on Ecm — the effective modulus for a long-term
// state). Along the strain gradient of that state: h the depth of the rings,
// x the compressed depth, d the depth of the centroid of the stretched bars
// (of the stretched bonded tendons when no bar is stretched). Ac,eff is the
// area of the rings within hc,ef of the most stretched vertex — exact for a
// polygon with voids and an oblique neutral axis; As and Ap are the bars and
// bonded tendons whose centre lies in that band. sigma_s is the largest
// stress variation Es eps of a stretched bar from the zero strain of the
// concrete at its level (§7.3.4(1): the tension of the cracked state for a
// bar without initial strain) — without stretched bar, the largest variation
// Ep eps of a bonded tendon (§7.3.4(2)). Typed errors: a section the action leaves
// compressed everywhere or stretched everywhere (the effective area of a
// member in tension, Fig. 7.1 c), is the caller's: crack_width_of), no bonded
// steel in the stretched zone, and those of cracked_service_state and
// crack_width_of.
CrackWidth crack_width(const RcPolySection& s, double Ecm, const Vec3& action, const PassiveSteel& bars,
                       const CrackWidthParams& p);

}  // namespace keta::ec2

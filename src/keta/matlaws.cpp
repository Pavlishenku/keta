#include "keta/matlaws.h"
#include "keta/model.h"

#include <algorithm>
#include <cmath>

namespace keta {
namespace {

// ---------------- small shared helpers ----------------

constexpr double kPaPerTm2 = kGravity * 1000.0;  // 1 t/m^2 = 9806.65 Pa
constexpr double kPaPerKnm2 = 1000.0;            // 1 kN/m^2 = 1 kPa

using Code = MatLaw::Code;

bool bpel_family(Code c) {
  return c == Code::BPEL || c == Code::BHP1 || c == Code::BHP2 || c == Code::BPEL_BL;
}
bool pcp_legacy(Code c) { return c == Code::IP2 || c == Code::MC78 || c == Code::COURBON; }
// laws that read the mean radius rm_cm = Ac/u (cm) rather than h_mm
bool reads_rm(Code c) { return bpel_family(c) || c == Code::IP2 || c == Code::MC78; }
bool bhp99(Code c) { return c == Code::BHP1 || c == Code::BHP2; }
bool ec2_bhp(Code c) { return c == Code::EC2_BHP1 || c == Code::EC2_BHP2; }
bool silica_fume(Code c) { return c == Code::BHP2 || c == Code::EC2_BHP2; }
bool lightweight(Code c) { return c == Code::BPEL_BL; }
bool in_range(double v, double lo, double hi) { return v >= lo && v <= hi; }

// Nominal cement coefficient s of the class. Used BOTH as the beta_cc exponent
// and as the selector of `ec2_t0_adjusted`, whose alpha mapping is
// s <= 0.21 -> +1, s >= 0.30 -> -1, else 0 — i.e. exactly the alpha the codes
// prescribe (S/SL -1, N 0, R +1; MC90 lumps N and R at alpha = 0).
double class_s(Code code, CementClass c) {
  if (code == Code::CEB) {
    // MC90 2.1-58: SL 0.38, N and R 0.25, RS 0.20.
    switch (c) {
      case CementClass::S: return 0.38;
      case CementClass::N: return 0.25;
      case CementClass::R: return 0.25;
      case CementClass::RS: return 0.20;
    }
  }
  // EN 1992-1-1 §3.1.2(6) / fib 5.1-51: S 0.38, N 0.25, R 0.20.
  switch (c) {
    case CementClass::S: return 0.38;
    case CementClass::N: return 0.25;
    case CementClass::R: return 0.20;
    case CementClass::RS: return 0.20;
  }
  return 0.25;
}

// Mean strength fcm = fck + delta_f [MPa]; delta_f = 8 MPa (EC2/CEB/fib),
// 8.3 MPa for B3 (ACI 209.2R-21).
double fcm_mpa(const MatLaw& m) {
  double fck = m.fc / 1e6;
  return fck + (m.code == Code::ACI_B3 ? 8.3 : 8.0);
}

// ---- French legacy laws (IP2, MC78) ----
// The source writes in T/m2: constants kept verbatim, converted at the
// boundary (1 T/m2 = kPaPerTm2 Pa). hm = S/p in metres.
double pcp_hm(const MatLaw& m) { return m.rm_cm / 100.0; }
double ip2_r(double x) {  // IP2, shared shrinkage / creep evolution
  x = std::max(x, 0.0);
  return std::max(0.0, 1.0 - 0.195 * std::exp(-0.127 * x) - 0.51 * std::exp(-0.0144 * x) -
                           0.295 * std::exp(-0.00178 * x));
}
double ip2_e28(const MatLaw& m) {  // IP2 28-day modulus [Pa]
  return 76650.0 * std::sqrt(m.fc / kPaPerTm2) * kPaPerTm2;
}
double ip2_ei(const MatLaw& m, double t) {  // J = 2, capped at 28 days
  if (t >= 28.0) return ip2_e28(m);
  return ip2_e28(m) * std::sqrt(std::max(0.0, 0.685 * std::log10(t + 1.0)));
}
double mc78_e28(const MatLaw& m) {  // MC78 28-day modulus [Pa]
  return 204671.0 * std::cbrt(m.fc / kPaPerTm2 + 800.0) * kPaPerTm2;
}
double mc78_alpha(const MatLaw& m) {  // MC78: 1 normal / slow, 2 rapid, 3 rapid high strength
  switch (m.cement) {
    case CementClass::R: return 2.0;
    case CementClass::RS: return 3.0;
    default: return 1.0;
  }
}
double mc78_ratio(double j) {  // fcm(j')/fcm(28), fit of MC78 fig. e.1
  j = std::max(j, 0.0);
  return 1.428 * (1.0 - 1.0 / (2.692e-5 * j * j + 0.0756 * j + 1.235));
}
double mc78_fict(const MatLaw& m, double t) { return mc78_alpha(m) * t * (m.theta + 10.0) / 30.0; }
// linear interpolation over the MC78 anchor humidities 40 / 70 / 90 / 100 %, clamped
double mc78_h(double H, double v40, double v70, double v90, double v100) {
  const double h[4] = {40, 70, 90, 100}, v[4] = {v40, v70, v90, v100};
  if (H <= h[0]) return v[0];
  if (H >= h[3]) return v[3];
  int i = H <= 70 ? 0 : (H <= 90 ? 1 : 2);
  return v[i] + (v[i + 1] - v[i]) * (H - h[i]) / (h[i + 1] - h[i]);
}
double mc78_h0(const MatLaw& m) { return 2.0 * mc78_h(m.rh, 1.0, 1.5, 5.0, 30.0) * pcp_hm(m); }
double mc78_beta_d(double d) {  // MC78 (constants inside the exponent: continuous at 1 and 100 d)
  d = std::max(d, 0.0);
  if (d <= 1.0) return 1.0 - std::exp(-0.4415 * d);
  if (d <= 100.0) return 1.0 - std::exp(-0.00769 * d - 0.433);
  return 1.0 - std::exp(-0.00378 * d - 0.825);
}
double mc78_beta_f(double h0, double tf) {
  const double psi = h0 < 0.40 ? 2000.0 * h0 + 200.0 : (h0 < 1.6 ? 1500.0 * h0 + 400.0 : 2800.0);
  tf = std::max(tf, 0.0);
  return std::cbrt(tf / (tf + psi));
}

// beta_cc(t) = exp(s (1 - sqrt(28/t))) — EN 1992-1-1 eq. 3.2, MC90 2.1-54,
// fib 5.1-51. s follows the class, except fib which forces s = 0.20 when
// fcm > 60 MPa (annex B3.8).
// EN 1992-1-1 (2027) table B.2 (ST1 annex B4.10.1): s_C by cement class and
// strength — CR 0.3 / 0.2 / 0.1, CN 0.5 / 0.4 / 0.3, CS 0.6 / 0.5 / 0.4 for
// fck < 35, 35..60, > 60 MPa.
double s_c_2027(const MatLaw& m) {
  const double fck = m.fc / 1e6;
  const int band = fck < 35.0 ? 0 : (fck <= 60.0 ? 1 : 2);
  static const double tab[3][3] = {{0.6, 0.5, 0.4}, {0.5, 0.4, 0.3}, {0.3, 0.2, 0.1}};
  const int c = m.cement == CementClass::S ? 0 : (m.cement == CementClass::N ? 1 : 2);
  return tab[c][band];
}

double beta_cc(const MatLaw& m, double t) {
  // 2027 (eq. B.2): exp(s_C (1 - sqrt(28 / t))) below the reference age of
  // 28 days, 1 beyond
  if (m.code == Code::EC2_2027)
    return t < 28.0 ? std::exp(s_c_2027(m) * (1.0 - std::sqrt(28.0 / t))) : 1.0;
  double s = class_s(m.code, m.cement);
  if (m.code == Code::FIB2010 && fcm_mpa(m) > 60.0) s = 0.20;
  return std::exp(s * (1.0 - std::sqrt(28.0 / t)));
}

// ---------------- BPEL / BHP99 family ----------------

// fcj [MPa]. BPEL 91 §2.1 (annex B3.3): j/(4.76 + 0.83 j) fc28 for
// fc28 <= 40 MPa, j/(1.40 + 0.95 j) fc28 above (the BHP99 law always takes the
// second branch). Lightweight: BPEL 99 annexe 5, 0.680 log10(2j+1) fc28 below
// 7 days then [0.519 + 0.332 log10 j] fc28 (the two agree to 3e-4 at j = 7 and
// reach 1.000 at j = 28). No "fcj = fc28 for j > 28" cap — see matlaws.h.
double bpel_fcj(const MatLaw& m, double j) {
  double fc28 = m.fc / 1e6;
  if (lightweight(m.code)) {
    if (j >= 28.0) return fc28;
    if (j < 7.0) return 0.680 * std::log10(2.0 * j + 1.0) * fc28;
    return (0.519 + 0.332 * std::log10(j)) * fc28;
  }
  if (bhp99(m.code) || fc28 > 40.0) return j / (1.40 + 0.95 * j) * fc28;
  return j / (4.76 + 0.83 * j) * fc28;
}

// ks = 1/(1 + 20 rho_s) — BPEL annexe 1 §3 and §4.2,2 (the bonded-steel
// abatement shared by the BPEL shrinkage and creep coefficients).
double bpel_ks(const MatLaw& m) { return 1.0 / (1.0 + 20.0 * m.ros); }

// Final BPEL shrinkage magnitude eps_r (annexe 1 §3): kre ks eps_0, with
// eps_0 = (100 - rho_h)(6 + 80/(10 + 3 rm)) 1e-6 for ordinary concrete and
// (100 - rho_h)(9 + 120/(10 + 3 rm)) 1e-6 for lightweight. `er` (a TOTAL final
// shrinkage for this family, manual p. 75) overrides the computed value.
double bpel_eps_r(const MatLaw& m) {
  if (m.er) return m.avance.k_retrait_endogene * *m.er;
  double rm = m.rm_cm;
  double e0 = lightweight(m.code) ? (100.0 - m.rh) * (9.0 + 120.0 / (10.0 + 3.0 * rm)) * 1e-6
                                  : (100.0 - m.rh) * (6.0 + 80.0 / (10.0 + 3.0 * rm)) * 1e-6;
  return m.avance.k_retrait_endogene * bpel_ks(m) * e0;
}

// BPEL shrinkage kinetics r(t) (annexe 1 §3): t/(t + 9 rm); lightweight adds
// the absorbed-water term Ea((rho_h - 50)/25 + k), k = (rm - 4)/rm if rm > 4.
double bpel_r(const MatLaw& m, double t) {
  double rm = m.rm_cm, extra = 0.0;
  if (lightweight(m.code)) {
    double k = rm > 4.0 ? (rm - 4.0) / rm : 0.0;
    extra = m.ea * ((m.rh - 50.0) / 25.0 + k);
  }
  double den = t + 9.0 * rm + extra;
  return den > 0 ? t / den : 0.0;
}

// BHP99 drying shrinkage eps_rd(t) since ts, magnitude (annex B3.4). Same law
// as EN 1992-2 B.116 with rm in cm: 8.4 rm^2 == 0.021 h^2 (h = 20 rm mm), and
// 2.8 rm^2 == 0.007 h^2 with silica fume.
double bhp99_eps_rd(const MatLaw& m, double t, double ts) {
  if (t <= ts) return 0.0;
  double fc28 = m.fc / 1e6;
  double K = fc28 < 57.0 ? 18.0 : 30.0 - 0.21 * fc28;
  double amp = m.er ? *m.er : K * (72.0 * std::exp(-0.046 * fc28) + 75.0 - m.rh) * 1e-6;
  double beta = silica_fume(m.code) ? 2.8 : 8.4;
  double d = t - ts;
  return amp * d / (d + beta * m.rm_cm * m.rm_cm);
}

// BHP99 autogenous shrinkage eps_re(age), magnitude (annex B3.4).
double bhp99_eps_re(const MatLaw& m, double age) {
  double fc28 = m.fc / 1e6;
  if (age <= 0) return 0.0;
  if (age <= 28.0) {
    double ratio = bpel_fcj(m, age) / fc28;
    if (ratio < 0.1) return 0.0;
    return (fc28 - 20.0) * (2.2 * ratio - 0.2) * 1e-6;
  }
  return (fc28 - 20.0) * (2.8 - 1.1 * std::exp(-age / 96.0)) * 1e-6;
}

// ---------------- EN 1992-2 Annex B (HP concrete) ----------------

// Drying shrinkage eps_cd(t) since ts, magnitude — EN 1992-2 B.116 (annex B3.7).
double ec2bhp_eps_cd(const MatLaw& m, double t, double ts) {
  if (t <= ts) return 0.0;
  double fck = m.fc / 1e6;
  double K = fck < 55.0 ? 18.0 : 30.0 - 0.21 * fck;
  double amp = m.er ? *m.er : K * (72.0 * std::exp(-0.046 * fck) + 75.0 - m.rh) * 1e-6;
  double beta_cd = silica_fume(m.code) ? 0.007 : 0.021;
  double d = t - ts;
  return amp * d / (d + beta_cd * m.h_mm * m.h_mm);
}

// Autogenous shrinkage eps_ca(age), magnitude — EN 1992-2 B.103 (annex B3.7).
double ec2bhp_eps_ca(const MatLaw& m, double age) {
  double fck = m.fc / 1e6;
  if (age <= 0) return 0.0;
  // the early-age branch for t < 28 days only (B.103; manual annex B4.7):
  // at 28 days exactly ST1 takes the long-term one
  if (age < 28.0) {
    double ratio = fcm_mpa(m) * beta_cc(m, age) / fck;
    if (ratio < 0.1) return 0.0;
    return (fck - 20.0) * (2.2 * ratio - 0.2) * 1e-6;
  }
  return (fck - 20.0) * (2.8 - 1.1 * std::exp(-age / 96.0)) * 1e-6;
}

// ---------------- EN 1992-1-1 (delegation to keta/creep.h) ----------------

ConcreteTime ec2_concrete(const MatLaw& m, double ts) {
  ConcreteTime c;
  c.fck = m.fc;
  c.rh = m.rh;
  c.h0 = m.h_mm / 1000.0;
  c.cement_s = class_s(Code::EC2, m.cement);
  c.t_shrink0 = ts;
  return c;
}

// EN 1992-1-1 eq. 3.11-3.13 autogenous shrinkage, magnitude. Needed ONLY to
// split `ec2_shrinkage` into its two AVANCE-weighted parts (and to honour an
// explicit ER, which the manual defines as the final DRYING value); the
// default total stays the creep module's, bit for bit.
double ec2_eps_ca(const MatLaw& m, double t) {
  double fck = m.fc / 1e6;
  return (1.0 - std::exp(-0.2 * std::sqrt(t))) * 2.5 * std::max(0.0, fck - 10.0) * 1e-6;
}
// EN 1992-1-1 eq. 3.10 drying kinetics — only needed on the explicit-ER path.
double ec2_beta_ds(const MatLaw& m, double t, double ts) {
  if (t <= ts) return 0.0;
  double h = m.h_mm;
  return (t - ts) / ((t - ts) + 0.04 * std::sqrt(h * h * h));
}

// ---------------- fib Model Code 2010 ----------------

double fib_alpha_ds1(CementClass c) {
  return c == CementClass::S ? 3.0 : (c == CementClass::N ? 4.0 : 6.0);
}
double fib_alpha_ds2(CementClass c) { return c == CementClass::S ? 0.13 : 0.12; }
double fib_alpha_as(CementClass c) {
  return c == CementClass::S ? 800.0 : (c == CementClass::N ? 700.0 : 600.0);
}

// ---------------- ACI 209.2R-21 model B3 ----------------

double b3_ecm_mpa(const MatLaw& m) { return 4734.0 * std::sqrt(fcm_mpa(m)); }
// eq. A-34 modulus growth factor sqrt(t/(4 + 0.85 t)).
double b3_age(double t) { return std::sqrt(t / (4.0 + 0.85 * t)); }

// tau_sh [days], eq. A-36 with h in mm. (The manual annotates kt in "days/m^2";
// that is a slip — with h in mm and kt = 8.5e-2 jc^-0.08 fcm^-0.25 a 100 mm
// slab gets tau_sh ~ 300 d, the documented B3 order of magnitude, whereas the
// m^2 reading would put it at 8e4 years.)
double b3_tau_sh(const MatLaw& m, double jc) {
  double kt = 8.5e-2 * std::pow(std::max(jc, 1.0), -0.08) * std::pow(fcm_mpa(m), -0.25);
  double d = m.k_forme * m.h_mm;
  return kt * d * d;
}
double b3_S(const MatLaw& m, double t, double jc) {
  if (t <= jc) return 0.0;
  return std::tanh(std::sqrt((t - jc) / b3_tau_sh(m, jc)));
}
// eq. A-33 eps_s_inf and eq. A-32 eps_sh_inf (magnitudes).
double b3_eps_s_inf(const MatLaw& m) {
  double a1 = m.b3_cement_type == 2 ? 0.85 : (m.b3_cement_type == 3 ? 1.10 : 1.00);
  return a1 * m.k_cure *
         (1.9e-2 * std::pow(m.water, 2.1) * std::pow(fcm_mpa(m), -0.28) + 270.0) * 1e-6;
}
double b3_eps_sh_inf(const MatLaw& m, double jc) {
  if (m.er) return *m.er;
  double tau = b3_tau_sh(m, jc);
  return b3_eps_s_inf(m) * b3_age(607.0) / b3_age(jc + tau);
}
// eq. A-6 humidity factor (the two branches agree at RH = 98 %, -0.2 at 100 %).
double b3_kh(const MatLaw& m) {
  double h = m.rh / 100.0;
  if (h <= 0.98) return 1.0 - h * h * h;
  return -12.9404 * h + 12.7404;
}
// Compliance parameters, all in 1/MPa (eq. A-38, A-41, A-46, A-47, A-49).
double b3_q1(const MatLaw& m) { return 0.6 / (m.avance.k_module * b3_ecm_mpa(m)); }
double b3_q2(const MatLaw& m) {
  return 185.4e-6 * std::sqrt(m.cement_qty) * std::pow(fcm_mpa(m), -0.9);
}
double b3_q3(const MatLaw& m) {
  double wc = m.water / m.cement_qty;
  return 0.29 * wc * wc * wc * wc * b3_q2(m);
}
double b3_q4(const MatLaw& m) { return 20.3e-6 * std::pow(m.aggregate / m.cement_qty, -0.7); }
// eq. A-49. |eps_sh_inf| enters in units of 1e-6, NOT as a bare strain: with a
// bare strain (5e-4) the factor |.|^-0.6 is 4000x larger and the drying
// compliance dwarfs q1. See the note in matlaws.h.
double b3_q5(const MatLaw& m, double jc) {
  double e = std::fabs(b3_eps_sh_inf(m, jc)) * 1e6;
  return 0.757 / fcm_mpa(m) * std::pow(std::max(e, 1e-9), -0.6);
}
double b3_delta001(const MatLaw& m, double t0) {
  return (b3_q3(m) + b3_q2(m) / std::sqrt(t0)) * std::log1p(std::pow(0.01, 0.1));
}
// Basic creep compliance C0(t, t0) [1/MPa], eq. A-40..A-45.
double b3_C0(const MatLaw& m, double t, double t0) {
  double d = std::max(t - t0, 0.0);
  double lnz = std::log1p(std::pow(d, 0.1));
  double Qf = 1.0 / (0.086 * std::pow(t0, 2.0 / 9.0) + 1.21 * std::pow(t0, 4.0 / 9.0));
  double Z = std::pow(t0, -0.5) * lnz;
  double Q = Z > 0 ? Qf * std::pow(1.0 + std::pow(Qf / Z, 1.7 * std::pow(t0, 0.12) + 8.0),
                                   -1.0 / (1.7 * std::pow(t0, 0.12) + 8.0))
                   : 0.0;
  return b3_q2(m) * Q + b3_q3(m) * lnz + b3_q4(m) * std::log(t / t0) - b3_delta001(m, t0);
}
// Drying creep compliance Cd(t, t0, jc) [1/MPa], eq. A-48/A-50.
double b3_Cd(const MatLaw& m, double t, double t0, double jc) {
  double rh = m.rh / 100.0;
  double Ht = 1.0 - b3_S(m, t, jc) * (1.0 - rh);
  double H0 = 1.0 - b3_S(m, t0, jc) * (1.0 - rh);
  double br = std::exp(-8.0 * Ht) - std::exp(-8.0 * H0);
  if (br <= 0) return 0.0;
  return b3_q5(m, jc) * std::sqrt(br);
}

}  // namespace

// ---------------- units ----------------

double strength_to_pa(double value, StrengthUnit unit) {
  require(finite(value) && value > 0, "strength_to_pa: strength must be finite and positive");
  switch (unit) {
    case StrengthUnit::MPa: return value * 1e6;
    case StrengthUnit::TM2: return value * kPaPerTm2;
    case StrengthUnit::KNM2: return value * kPaPerKnm2;
  }
  fail("strength_to_pa: unknown unit");
}

double strength_to_mpa(double value, StrengthUnit unit) {
  return strength_to_pa(value, unit) / 1e6;
}

// ---------------- validation ----------------

void MatLaw::validate() const {
  const AvanceFactors& a = avance;
  require(finite(a.k_module) && finite(a.k_retrait_endogene) && finite(a.k_retrait_dessication) &&
              finite(a.k_fluage_endogene) && finite(a.k_fluage_dessication) &&
              finite(a.k_fluage_cinetique) && finite(a.k_retrait_endogene_cin) &&
              finite(a.k_retrait_dessication_cin) && finite(a.k_fluage_endogene_cin) &&
              finite(a.k_fluage_dessication_cin),
          "MatLaw: AVANCE multipliers must be finite");
  require(a.k_module > 0, "MatLaw: K_MODULE must be > 0");
  if (code == Code::Direct) {
    require(finite(e_direct) && e_direct > 0, "MatLaw(Direct): E must be finite and positive");
    return;
  }
  require(finite(fc) && fc > 0, "MatLaw: fc must be finite and positive");
  require(finite(rh) && rh > 0 && rh <= 100, "MatLaw: RH must lie in (0, 100] %");
  require(finite(ts_dry) && ts_dry >= 0, "MatLaw: drying start age must be finite and >= 0");
  require(finite(ros) && ros >= 0 && ros < 1, "MatLaw: ROS must lie in [0, 1)");
  require(finite(es) && es > 0, "MatLaw: Es must be finite and positive");
  if (er) require(finite(*er) && *er >= 0, "MatLaw: ER must be finite and >= 0");
  if (reads_rm(code))
    require(finite(rm_cm) && rm_cm > 0, "MatLaw: RM (mean radius, cm) must be finite and positive");
  else if (code != Code::COURBON)
    require(finite(h_mm) && h_mm > 0, "MatLaw: H (notional size, mm) must be finite and positive");
  if (pcp_legacy(code)) {
    require(finite(theta) && theta > -10, "MatLaw(IP2/MC78): ambient temperature must be finite and > -10 C");
    require(finite(kb) && kb > 0, "MatLaw(IP2): kb must be finite and positive");
    require(finite(k1_consistency) && k1_consistency > 0, "MatLaw(MC78): k1 must be finite and positive");
    if (code == Code::MC78)
      require(finite(phi_f2) && phi_f2 > 0,
              "MatLaw(MC78): phi_f2 (MC78 fig. e.2, not tabulated by the sources) must be given (> 0)");
    if (code == Code::COURBON)
      require(finite(e_direct) && e_direct > 0, "MatLaw(COURBON): the constant modulus e_direct must be > 0");
  }

  double f = fc / 1e6;
  switch (code) {
    case Code::BPEL:
    case Code::BPEL_BL:
      // Manual p. 74/75: "BPEL: beton courant (fc28 <= 60 MPa)".
      require(f <= 60.0, "MatLaw(BPEL): fc28 must be <= 60 MPa");
      if (lightweight(code))
        require(finite(robs) && robs > 0 && finite(ea) && ea >= 0,
                "MatLaw(BPEL_BL): ROBS must be positive and EA >= 0");
      break;
    case Code::BHP1:
    case Code::BHP2:
      // Manual p. 75: "40 MPa <= fc28 <= 80 MPa".
      require(in_range(f, 40.0, 80.0), "MatLaw(BHP): fc28 must lie in [40, 80] MPa");
      break;
    case Code::EC2_BHP1:
    case Code::EC2_BHP2:
      // Manual p. 82: "40 MPa <= fck <= 90 MPa".
      require(in_range(f, 40.0, 90.0), "MatLaw(EC2_BHP): fck must lie in [40, 90] MPa");
      break;
    case Code::ACI_B3: {
      require(finite(water) && finite(cement_qty) && finite(aggregate) && cement_qty > 0,
              "MatLaw(ACI_B3): composition must be finite with c > 0");
      require(in_range(cement_qty, 160.0, 720.0),
              "MatLaw(ACI_B3): cement content must lie in [160, 720] kg/m3");
      require(in_range(water / cement_qty, 0.35, 0.85),
              "MatLaw(ACI_B3): w/c must lie in [0.35, 0.85]");
      require(in_range(aggregate / cement_qty, 2.5, 13.5),
              "MatLaw(ACI_B3): a/c must lie in [2.5, 13.5]");
      require(finite(k_forme) && k_forme > 0 && finite(k_cure) && k_cure > 0,
              "MatLaw(ACI_B3): K_FORME and K_CURE must be finite and positive");
      require(b3_cement_type >= 1 && b3_cement_type <= 3,
              "MatLaw(ACI_B3): CIMENT must be 1, 2 or 3 (type I, II, III)");
      break;
    }
    default: break;
  }
}

// ---------------- strength and modulus ----------------

double MatLaw::ip2_kb_from_mix(double c, double w) {
  require(c > 0 && w > 0, "MatLaw::ip2_kb_from_mix: cement and water contents > 0");
  return (w / c) * (7.9e-3 * c + 0.766) - (9.4e-4 * c + 0.397);
}

double MatLaw::fc_mean(double t_days) const {
  validate();
  require(t_days > 0, "MatLaw::fc_mean: age must be > 0");
  switch (code) {
    case Code::Direct: return fc;
    case Code::BPEL:
    case Code::BHP1:
    case Code::BHP2:
    case Code::BPEL_BL: return bpel_fcj(*this, t_days) * 1e6;
    // ACI 209.2R-21 eq. A-34 is E ~ sqrt(t/(4 + 0.85 t)), i.e. strength grows
    // as t/(4 + 0.85 t).
    case Code::ACI_B3: return fcm_mpa(*this) * t_days / (4.0 + 0.85 * t_days) * 1e6;
    // IP2 / COURBON: no strength-growth law in the source
    case Code::IP2:
    case Code::COURBON: return fc;
    case Code::MC78: return mc78_ratio(std::min(mc78_fict(*this, t_days), 400.0)) * fcm_mpa(*this) * 1e6;
    default: return fcm_mpa(*this) * beta_cc(*this, t_days) * 1e6;
  }
}

double MatLaw::E(double t_days) const {
  validate();
  require(t_days > 0, "MatLaw::E: age must be > 0");
  double km = avance.k_module;
  switch (code) {
    case Code::Direct: return km * e_direct;
    // BPEL 91 §2.1: Eij = 11000 fcj^(1/3) [MPa].
    case Code::BPEL:
    case Code::BHP1:
    case Code::BHP2: return km * 11000.0 * std::cbrt(bpel_fcj(*this, t_days)) * 1e6;
    // BPEL 99 annexe 5 lightweight — see the reconstruction note in matlaws.h.
    case Code::BPEL_BL:
      return km * 1650.0 * robs * robs * std::cbrt(bpel_fcj(*this, t_days)) * 1e6;
    // MC90 2.1-15 + 2.1-57: Eci(t) = 21500 (fcm/10)^(1/3) sqrt(beta_cc(t)).
    case Code::CEB:
      return km * 21500.0 * std::cbrt(fcm_mpa(*this) / 10.0) * std::sqrt(beta_cc(*this, t_days)) *
             1e6;
    // EN 1992-1-1 Table 3.1 + eq. 3.5: Ecm(t) = Ecm beta_cc(t)^0.3.
    case Code::EC2:
    case Code::EC2_BHP1:
    case Code::EC2_BHP2: return km * ec2_ecm(fc) * std::pow(beta_cc(*this, t_days), 0.3);
    // fib 5.1-21/23 + 5.1-56/57: Eci(t) = 21500 (fcm/10)^(1/3) alpha_i beta_cc^0.5.
    case Code::FIB2010: {
      double ai = fib_secant ? std::min(0.8 + 0.2 * fcm_mpa(*this) / 88.0, 1.0) : 1.0;
      return km * 21500.0 * std::cbrt(fcm_mpa(*this) / 10.0) * ai *
             std::sqrt(beta_cc(*this, t_days)) * 1e6;
    }
    // ACI 209.2R-21 eq. A-39 + A-34.
    case Code::ACI_B3: return km * b3_ecm_mpa(*this) * b3_age(t_days) * 1e6;
    // 2027 (EN 1992-1-1:2023 eq. B.4, ST1 annex B4.10.1): Ecm = 9800 fcm^(1/3),
    // E(j) = beta_cc^(1/3) Ecm. (ST1 2.24 keeps the EC2 2004 exponent 0.3 —
    // E(7) = 31 588 MPa — an anomaly; keta follows the standard.)
    case Code::EC2_2027:
      return km * 9800.0 * std::cbrt(fcm_mpa(*this)) * std::cbrt(beta_cc(*this, t_days)) * 1e6;
    case Code::IP2: return km * ip2_ei(*this, t_days);
    case Code::MC78:
      return km * mc78_e28(*this) *
             std::cbrt(std::max(0.0, mc78_ratio(std::min(mc78_fict(*this, t_days), 400.0))));
    case Code::COURBON: return km * e_direct;
  }
  fail("MatLaw::E: unknown law");
}

double MatLaw::Ev(double t_days) const {
  validate();
  require(code == Code::BPEL || code == Code::BHP1 || code == Code::BHP2,
          "MatLaw::Ev: the deferred modulus 3700 fcj^(1/3) is a BPEL/BHP law only");
  require(t_days > 0, "MatLaw::Ev: age must be > 0");
  // BPEL 91 §2.1.3: Evj = 3700 fcj^(1/3) [MPa].
  return avance.k_module * 3700.0 * std::cbrt(bpel_fcj(*this, t_days)) * 1e6;
}

double MatLaw::E_ref(double t0) const {
  validate();
  switch (code) {
    case Code::Direct: return avance.k_module * e_direct;
    // EN 1992-1-1 eq. 3.6 references phi to the tangent modulus Ec = 1.05 Ecm.
    case Code::EC2:
    case Code::EC2_BHP1:
    case Code::EC2_BHP2: return 1.05 * avance.k_module * ec2_ecm(fc);
    // B3: the conventional instantaneous compliance q1 + Delta_0.01 (t0-dependent).
    case Code::ACI_B3:
      require(t0 > 0, "MatLaw::E_ref: t0 must be > 0");
      return 1e6 / (b3_q1(*this) + b3_delta001(*this, t0));
    // 2027 (eq. B.21/B.22): the tangent modulus Ec = Ecm / (0.8 + 0.2 fcm/88) >= Ecm.
    // (ST1 2.24 refers phi to Ecm — its phi equals fib's to 4 digits —
    // an anomaly; keta follows the standard.)
    case Code::EC2_2027: return E(28.0) / std::min(1.0, 0.8 + 0.2 * fcm_mpa(*this) / 88.0);
    // MC78: J = (...) / E28, the 28-day modulus of the code (not E(28) of the
    // fictitious-age fit, 0.15 % off it)
    case Code::MC78: return avance.k_module * mc78_e28(*this);
    // Everything else references phi to its own 28-day modulus.
    default: return E(28.0);
  }
}

double MatLaw::n_equiv() const {
  validate();
  // Annex B3.4: n = 15 for 40 <= fc28 < 60, n = 9 for 60 <= fc28 <= 80.
  if (bhp99(code)) return fc / 1e6 < 60.0 ? 15.0 : 9.0;
  // ACI_B3: ST1 takes n = 15 (delayed strain / analytic = 1 + 15 ROS =
  // 1.75000 at every state; Es/E28 would give 1.7419).
  if (code == Code::ACI_B3) return 15.0;
  if (ec2_bhp(code) || code == Code::FIB2010 || code == Code::EC2_2027 || code == Code::IP2)
    return es / E(28.0);
  return 0.0;  // BPEL folds it into ks; CEB and EC2 have no (1 + n rho_s) term
}

// ---------------- creep ----------------

double MatLaw::phi(double t, double t0) const {
  validate();
  if (code == Code::Direct) return 0.0;
  require(t0 > 0 && t >= t0, "MatLaw::phi: need t >= t0 > 0");
  const AvanceFactors& a = avance;
  double d = t - t0;

  switch (code) {
    // ---- BPEL 91 annexe 1 §4.2 (and its lightweight variant) ----
    // phi = Kfl(t0) f(t - t0), Kfl = ks [ke + kc k(t0)], referenced to Ei28.
    case Code::BPEL:
    case Code::BPEL_BL: {
      double rm = rm_cm;
      double ke = lightweight(code) ? 0.20 : 0.40;
      double kc = lightweight(code)
                      ? (120.0 - rh) / 50.0 + 0.4 * (100.0 - rh) / (20.0 + rm)
                      : (120.0 - rh) / 30.0 + (2.0 / 3.0) * (100.0 - rh) / (20.0 + rm);
      double kfl = bpel_ks(*this) * (ke + kc * 100.0 / (100.0 + t0));
      double f = std::sqrt(d) / (std::sqrt(d) + 5.0 * std::sqrt(rm));
      return a.k_fluage_endogene * kfl * f;
    }

    // ---- AFREM/BHP 1995, annex B3.4 ----
    // phi = [kfe Kfl f + kfd phi_d0 (eps_rd(t) - eps_rd(t0))] / (1 + n rho_s),
    // referenced to Ei28.
    case Code::BHP1:
    case Code::BHP2: {
      double fc28 = fc / 1e6, fcj0 = bpel_fcj(*this, t0);
      double kfl, beta;
      if (silica_fume(code)) {
        kfl = 3.6 / std::pow(fcj0, 0.37);
        beta = a.k_fluage_cinetique * 0.37 * std::exp(2.8 * fcj0 / fc28);
      } else {
        kfl = 1.4;
        beta = a.k_fluage_cinetique * 0.40 * std::exp(3.1 * fcj0 / fc28);
      }
      double sd = std::sqrt(d);
      double basic = sd + beta > 0 ? kfl * sd / (sd + beta) : 0.0;
      double phid0 = silica_fume(code) ? 1e3 : 3.2e3;
      double drying = phid0 * (bhp99_eps_rd(*this, t, ts_dry) - bhp99_eps_rd(*this, t0, ts_dry));
      return (a.k_fluage_endogene * basic + a.k_fluage_dessication * drying) /
             (1.0 + n_equiv() * ros);
    }

    // ---- CEB-FIP MC90, annex B3.5 (eq. 2.1-61..71) ----
    case Code::CEB: {
      double fcm = fcm_mpa(*this);
      double phi_rh = 1.0 + (1.0 - rh / 100.0) / (0.46 * std::cbrt(h_mm / 100.0));
      double beta_fcm = 5.3 / std::sqrt(fcm / 10.0);
      double t0a = ec2_t0_adjusted(t0, class_s(code, cement));  // MC90 2.1-72
      double beta_t0 = 1.0 / (0.1 + std::pow(t0a, 0.2));
      double bh = std::min(
          150.0 * (1.0 + std::pow(1.2 * rh / 100.0, 18)) * (h_mm / 100.0) + 250.0, 1500.0);
      double beta_c = std::pow(d / (bh + d), 0.3);
      return a.k_fluage_endogene * phi_rh * beta_fcm * beta_t0 * beta_c;
    }

    // ---- EN 1992-1-1 Annex B: delegated to keta/creep.h ----
    case Code::EC2: return a.k_fluage_endogene * ec2_phi(t, t0, ec2_concrete(*this, ts_dry));

    // ---- EN 1992-2 Annex B (HP), annex B3.7 ----
    case Code::EC2_BHP1:
    case Code::EC2_BHP2: {
      double fck = fc / 1e6, fcm0 = fcm_mpa(*this) * beta_cc(*this, t0);
      double phib0, beta_bc;
      if (silica_fume(code)) {
        phib0 = 3.6 / std::pow(fcm0, 0.37);
        beta_bc = a.k_fluage_cinetique * 0.37 * std::exp(2.8 * fcm0 / fck);
      } else {
        phib0 = 1.4;
        beta_bc = a.k_fluage_cinetique * 0.40 * std::exp(3.1 * fcm0 / fck);
      }
      double sd = std::sqrt(d);
      double basic = sd + beta_bc > 0 ? phib0 * sd / (sd + beta_bc) : 0.0;
      double phid0 = silica_fume(code) ? 1000.0 : 3200.0;
      double drying = phid0 * (ec2bhp_eps_cd(*this, t, ts_dry) - ec2bhp_eps_cd(*this, t0, ts_dry));
      return (a.k_fluage_endogene * basic + a.k_fluage_dessication * drying) /
             (1.0 + n_equiv() * ros);
    }

    // ---- fib Model Code 2010, annex B3.8 (eq. 5.1-63..73) ----
    case Code::FIB2010: {
      double fcm = fcm_mpa(*this);
      double t0a = ec2_t0_adjusted(t0, class_s(code, cement));  // eq. 5.1-73
      // basic creep, eq. 5.1-64..66
      double phib0 = 1.8 / std::pow(fcm, 0.7);
      double g = 30.0 / t0a + 0.035;
      double basic = phib0 * std::log1p(std::max(a.k_fluage_cinetique, 0.0) * g * g * d);
      // drying creep, eq. 5.1-67..71
      double beta_dc_fcm = 412.0 / std::pow(fcm, 1.4);
      double beta_rh = (1.0 - rh / 100.0) / std::cbrt(0.1 * h_mm / 100.0);
      double beta_dc_t0 = 1.0 / (0.1 + std::pow(t0a, 0.2));
      double alpha_fcm = std::sqrt(35.0 / fcm);
      double bh = std::min(1.5 * h_mm + 250.0 * alpha_fcm, 1500.0 * alpha_fcm);
      double gamma = 1.0 / (2.3 + 3.5 / std::sqrt(t0a));
      double drying = beta_dc_fcm * beta_rh * beta_dc_t0 * std::pow(d / (bh + d), gamma);
      return (a.k_fluage_endogene * basic + a.k_fluage_dessication * drying) /
             (1.0 + n_equiv() * ros);
    }

    // ---- EN 1992-1-1 (2027) annex B, ST1 annex B4.10.1 (eq. B.6..B.16) ----
    case Code::EC2_2027: {
      double fcm = fcm_mpa(*this);
      double t0a = ec2_t0_adjusted(t0, class_s(Code::EC2, cement));  // eq. B.17
      // basic creep, eq. B.6..B.8 (kfec divides the time under load)
      double alpha = silica ? std::min(0.7, 0.7 - 0.11 * (fcm - 30.0) / 90.0) : 0.7;
      double g = 0.035 + 30.0 / t0a;
      double basic = 1.8 / std::pow(fcm, alpha) *
                     std::log1p(g * g * d / std::max(a.k_fluage_endogene_cin, 1e-12));
      // drying creep, eq. B.9..B.16 (kfdc scales beta_H)
      double b = silica ? std::min(1.4, 1.4 - 0.25 * (fcm - 30.0) / 90.0) : 1.4;
      double beta_dc_fcm = 412.0 / std::pow(fcm, b);
      double beta_rh = (1.0 - rh / 100.0) / std::cbrt(h_mm / 1000.0);
      double beta_dc_t0 = 1.0 / (0.1 + std::pow(t0a, 0.2));
      double alpha_fcm = std::sqrt(35.0 / fcm);
      double bh = std::min(1.5 * h_mm + 250.0 * alpha_fcm, 1500.0 * alpha_fcm);
      double gamma = 1.0 / (2.3 + 3.5 / std::sqrt(t0a));
      double drying = beta_dc_fcm * beta_rh * beta_dc_t0 *
                      std::pow(d / (bh * a.k_fluage_dessication_cin + d), gamma);
      return (a.k_fluage_endogene * basic + a.k_fluage_dessication * drying) /
             (1.0 + n_equiv() * ros);
    }

    // ---- ACI 209.2R-21 model B3, annex B3.9 ----
    // phi = (kfe C0 + kfd Cd) / [(q1 + Delta_0.01)(1 + n rho_s)].
    case Code::ACI_B3: {
      double j0 = b3_q1(*this) + b3_delta001(*this, t0);
      double c0 = std::max(b3_C0(*this, t, t0), 0.0);
      double cd = b3_Cd(*this, t, t0, ts_dry);
      return (a.k_fluage_endogene * c0 + a.k_fluage_dessication * cd) / j0 /
             (1.0 + n_equiv() * ros);
    }

    // ---- IP2, phi = J_def E28, J_def = kb kc kd ker r / Ei(t0) ----
    case Code::IP2: {
      const double H = rh, hmean = 2.0 * pcp_hm(*this);  // ker reads the mean thickness 2 S/p
      const double kc = -2.73e-4 * H * H + 3.53e-3 * H + 3.38;
      const double t0f = t0 * (theta + 10.0) / 30.0;
      const double kd = (2.82 * t0f + 173.0) / (5.75 * t0f + 94.3);
      const double ker = (464.0 * hmean + 94.2) / (834.0 * hmean + 54.7);
      return a.k_fluage_endogene * kb * kc * kd * ker * ip2_r(d) * ip2_e28(*this) / ip2_ei(*this, t0);
    }
    // ---- MC78, J = (0.4 beta_d + phi_f dbeta_f) / E28 ----
    case Code::MC78: {
      const double tf = mc78_fict(*this, t), t0f = mc78_fict(*this, t0), h0 = mc78_h0(*this);
      const double phi_f = k1_consistency * mc78_h(rh, 3.0, 2.0, 1.0, 0.8) * phi_f2;
      return a.k_fluage_endogene *
             std::max(0.0, 0.4 * mc78_beta_d(tf - t0f) + phi_f * (mc78_beta_f(h0, tf) - mc78_beta_f(h0, t0f)));
    }
    // ---- Courbon, non-aging ----
    case Code::COURBON: return a.k_fluage_endogene * 2.0 * (1.0 - std::exp(-0.006 * d));

    default: break;
  }
  fail("MatLaw::phi: unknown law");
}

// ---------------- shrinkage ----------------
//
// Convention (documented in matlaws.h): the DRYING part is measured from ts,
// the AUTOGENOUS part is the law's own function of the age since casting —
// exactly what EN 1992-1-1 eq. 3.8 does. Negative = shortening.

double MatLaw::eps_sh(double t, double ts) const {
  validate();
  if (code == Code::Direct) return 0.0;
  require(t > 0 && ts >= 0, "MatLaw::eps_sh: need t > 0 and ts >= 0");
  const AvanceFactors& a = avance;

  switch (code) {
    // ---- BPEL 91 annexe 1 §3: eps(ts, t) = eps_r [r(t) - r(ts)] ----
    case Code::BPEL:
    case Code::BPEL_BL:
      if (t <= ts) return 0.0;
      return -bpel_eps_r(*this) * (bpel_r(*this, t) - bpel_r(*this, ts));

    // ---- AFREM/BHP 1995, annex B3.4 ----
    case Code::BHP1:
    case Code::BHP2:
      return -(a.k_retrait_dessication * bhp99_eps_rd(*this, t, ts) +
               a.k_retrait_endogene * bhp99_eps_re(*this, t)) /
             (1.0 + n_equiv() * ros);

    // ---- CEB-FIP MC90, annex B3.5 (eq. 2.1-74..79) ----
    case Code::CEB: {
      if (t <= ts) return 0.0;
      double fcm = fcm_mpa(*this);
      double beta_sc = cement == CementClass::S ? 4.0 : (cement == CementClass::RS ? 8.0 : 5.0);
      double beta_rh = rh >= 99.0 ? 0.25 : -1.55 * (1.0 - std::pow(rh / 100.0, 3));
      double eps_cs0 = er ? -*er : (160.0 + 10.0 * beta_sc * (9.0 - fcm / 10.0)) * beta_rh * 1e-6;
      double d = t - ts;
      double beta_s = std::sqrt(d / (350.0 * (h_mm / 100.0) * (h_mm / 100.0) + d));
      return a.k_retrait_endogene * eps_cs0 * beta_s;
    }

    // ---- EN 1992-1-1: delegated to keta/creep.h ----
    case Code::EC2: {
      double ca = ec2_eps_ca(*this, t);
      if (er) return -(a.k_retrait_dessication * *er * ec2_beta_ds(*this, t, ts) +
                       a.k_retrait_endogene * ca);
      ConcreteTime c = ec2_concrete(*this, ts);
      // Same multiplier on both parts => the creep module's total, bit for bit.
      if (a.k_retrait_dessication == a.k_retrait_endogene)
        return a.k_retrait_endogene * ec2_shrinkage(t, c);
      double total = -ec2_shrinkage(t, c);  // magnitude
      return -(a.k_retrait_dessication * (total - ca) + a.k_retrait_endogene * ca);
    }

    // ---- EN 1992-2 Annex B (HP), annex B3.7 ----
    case Code::EC2_BHP1:
    case Code::EC2_BHP2:
      return -(a.k_retrait_dessication * ec2bhp_eps_cd(*this, t, ts) +
               a.k_retrait_endogene * ec2bhp_eps_ca(*this, t)) /
             (1.0 + n_equiv() * ros);

    // ---- fib Model Code 2010, annex B3.8 (eq. 5.1-75..83) ----
    case Code::FIB2010: {
      double fcm = fcm_mpa(*this);
      // autogenous, eq. 5.1-76/78/79 (negative = shortening)
      double eps_cas0 =
          -fib_alpha_as(cement) * std::pow((fcm / 10.0) / (6.0 + fcm / 10.0), 2.5) * 1e-6;
      double ca = eps_cas0 * (1.0 - std::exp(-0.2 * std::sqrt(t)));
      // drying, eq. 5.1-77/80..83
      double beta_s1 = std::min(std::pow(35.0 / fcm, 0.1), 1.0);
      double beta_rh = rh >= 99.0 * beta_s1 ? 0.25 : -1.55 * (1.0 - std::pow(rh / 100.0, 3));
      double eps_cds0 = (220.0 + 110.0 * fib_alpha_ds1(cement)) *
                        std::exp(-fib_alpha_ds2(cement) * fcm / 10.0) * 1e-6;
      double d = std::max(t - ts, 0.0);
      double beta_ds = std::sqrt(d / (d + 0.035 * h_mm * h_mm));
      double cd = er ? -*er * beta_ds : eps_cds0 * beta_ds * beta_rh;
      return (a.k_retrait_dessication * cd + a.k_retrait_endogene * ca) / (1.0 + n_equiv() * ros);
    }

    // ---- EN 1992-1-1 (2027) annex B, ST1 annex B4.10.1 (eq. B.23..B.33) ----
    case Code::EC2_2027: {
      double fcm = fcm_mpa(*this);
      const int c = cement == CementClass::S ? 0 : (cement == CementClass::N ? 1 : 2);
      static const double alpha_bs[3] = {800.0, 700.0, 600.0};
      static const double alpha_ds[3] = {3.0, 4.0, 6.0};
      // autogenous, eq. B.24/B.26/B.27 (negative = shortening): beta_bs =
      // 1 - exp(-0.2 krec sqrt(t)), krec multiplying the exponent
      double ca = -alpha_bs[c] * std::pow(fcm / (60.0 + fcm), 2.5) * 1e-6 *
                  (1.0 - std::exp(-0.2 * a.k_retrait_endogene_cin * std::sqrt(t)));
      // drying, eq. B.25/B.28..B.33 (krdc scales the 0.035 h^2 term)
      double rh_eq = std::min(99.0, 90.0 * std::pow(35.0 / fcm, 0.1));
      double beta_rh = rh <= rh_eq ? 1.55 * (1.0 - std::pow(rh / rh_eq, 3))
                                   : 1.55 * (1.0 - std::pow(rh / rh_eq, 2)) -
                                         (rh >= 100.0 ? 0.25 : 0.0);
      double eps_cds = (220.0 + 110.0 * alpha_ds[c]) * std::exp(-0.012 * fcm) * 1e-6;
      double dd = std::max(t - ts, 0.0);
      double beta_ds = std::sqrt(dd / (dd + 0.035 * a.k_retrait_dessication_cin * h_mm * h_mm));
      double cd = er ? -*er * beta_ds : -eps_cds * beta_rh * beta_ds;
      return (a.k_retrait_dessication * cd + a.k_retrait_endogene * ca) / (1.0 + n_equiv() * ros);
    }

    // ---- ACI 209.2R-21 model B3, annex B3.9 (eq. A-31..A-36) ----
    // B3 as parameterized by ST1 carries a drying term only: K_RETRAIT_ENDOGENE
    // is a documented no-op here.
    case Code::ACI_B3: {
      double s = b3_S(*this, t, ts);
      double cd = er ? *er * s : b3_eps_sh_inf(*this, ts) * b3_kh(*this) * s;
      return -a.k_retrait_dessication * cd / (1.0 + n_equiv() * ros);
    }

    // ---- IP2, eps = kb eps_c kc1 kD r(t'), from casting ----
    case Code::IP2: {
      const double hm = pcp_hm(*this), H = rh;
      const double eps_c = (-8.68e-3 * H * H + 0.545 * H + 32.4) * 1e-5;
      const double kc1 = (47.6 * hm + 155.0) / (525.0 * hm + 101.0);
      const double fin = er ? *er : kb * eps_c * kc1 / (1.0 + n_equiv() * ros);
      return -a.k_retrait_endogene * fin * ip2_r(t * std::sqrt(1.0 / (10.0 * hm)));
    }
    // ---- MC78, eps = k1 eps_s1(H) eps_s2(h0) beta_s(t') ----
    case Code::MC78: {
      const double h0 = mc78_h0(*this), hk = std::clamp(h0, 0.05, 0.99);
      const double eps_s1 = mc78_h(rh, -5.2e-4, -3.2e-4, -1.3e-4, 1.0e-4);
      const double eps_s2 = 0.7166 + 0.0333 / h0;
      const double r1 = (-1.85e-3 * hk + 0.0879) / (258.0 * hk - 3.99);
      const double r2 = (-1.61 * hk + 1.6) / (334.0 * hk - 9.36);
      const double tp = t * (theta + 10.0) / 30.0;
      const double bs = 1.0 - 0.633 * std::exp(-r1 * tp) - 0.367 * std::exp(-r2 * tp);
      const double fin = er ? -*er : k1_consistency * eps_s1 * eps_s2;
      return a.k_retrait_endogene * fin * bs;
    }
    case Code::COURBON: return 0.0;

    default: break;
  }
  fail("MatLaw::eps_sh: unknown law");
}

// ---------------- composition ----------------

std::shared_ptr<ConcreteLaw> concrete_law_from(const MatLaw& creep, const MatLaw& shrinkage, double ts,
                                               double stress_scale, bool nonlinear) {
  creep.validate();
  shrinkage.validate();
  require(ts >= 0 && stress_scale > 0, "concrete_law_from: ts >= 0 and a positive stress scale");
  auto cl = std::make_shared<ConcreteLaw>();
  cl->phi = [creep](double t, double t0) { return creep.phi(t, t0); };
  cl->E_ref = [creep, stress_scale](double t0) { return creep.E_ref(t0) / stress_scale; };
  cl->shrinkage = [shrinkage, ts](double t) { return t > 0 ? shrinkage.eps_sh(t, ts) : 0.0; };
  if (nonlinear) cl->nonlinear_fcm = [creep, stress_scale](double t0) { return creep.fc_mean(t0) / stress_scale; };
  return cl;
}

// ---------------- non-linear creep ----------------

double nonlinear_creep_factor(double sigma_c, double fcm_t0) {
  require(finite(sigma_c) && finite(fcm_t0) && fcm_t0 > 0,
          "nonlinear_creep_factor: need a finite stress and fcm(t0) > 0");
  // EN 1992-1-1 eq. 3.7: k_sigma = sigma_c / fcm(t0), floored at 0.45 so the
  // factor is 1 below the threshold.
  double k = std::max(0.45, std::fabs(sigma_c) / fcm_t0);
  return std::exp(1.5 * (k - 0.45));
}

// ---------------- relaxation threshold ----------------

double relax_seuil_hours(double t_hours, std::optional<double> seuil) {
  require(finite(t_hours) && t_hours >= 0, "relax_seuil_hours: t must be finite and >= 0");
  if (!seuil) return t_hours;
  require(finite(*seuil) && *seuil >= 0, "relax_seuil_hours: SEUIL must be finite and >= 0");
  double cap = *seuil == 0.0 ? kRelaxSeuilDefaultHours : *seuil;
  return std::min(t_hours, cap);
}

double relaxation_ratio_seuil(double t_hours, double sigma_pi, double fpk, int relaxation_class,
                              double rho1000, std::optional<double> seuil) {
  return ec2_relaxation_ratio(relax_seuil_hours(t_hours, seuil), sigma_pi, fpk, relaxation_class,
                              rho1000);
}

}  // namespace keta

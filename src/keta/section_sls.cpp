// Serviceability values of a concrete section — see section_sls.h.
#include "keta/section_sls.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace keta::ec2 {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0; }

// The largest tensile stress of a list, 0 when it is empty or all compressed.
double largest_tension(const std::vector<double>& stresses) {
  double largest = 0;
  for (double s : stresses) largest = std::max(largest, s);
  return largest;
}

}  // namespace

StressLimits stress_limits(double fck, double fyk, double fpk, const StressLimitFactors& k) {
  require(positive(fck) && positive(fyk) && positive(fpk), "stress limits: fck, fyk and fpk > 0");
  require(positive(k.k1) && positive(k.k2) && positive(k.k3) && positive(k.k4) && positive(k.k5),
          "stress limits: factors k1 to k5 > 0");
  return {k.k1 * fck, k.k2 * fck, k.k3 * fyk, k.k4 * fyk, k.k5 * fpk};
}

StressRatios stress_ratios(const PrestressServiceState& state, const StressLimits& limits) {
  require(positive(limits.concrete_characteristic) && positive(limits.concrete_quasi_permanent) &&
              positive(limits.bar) && positive(limits.bar_imposed) && positive(limits.tendon),
          "stress ratios: limits > 0");
  StressRatios r;
  r.concrete_compression = std::max(0.0, -state.fibre_min.stress);
  r.bar_tension = largest_tension(state.bar_stress);
  r.tendon_tension = largest_tension(state.tendon_stress);
  r.concrete_characteristic = r.concrete_compression / limits.concrete_characteristic;
  r.concrete_quasi_permanent = r.concrete_compression / limits.concrete_quasi_permanent;
  r.bar = r.bar_tension / limits.bar;
  r.bar_imposed = r.bar_tension / limits.bar_imposed;
  r.tendon = r.tendon_tension / limits.tendon;
  return r;
}

ServiceState service_state(const RcPolySection& s, double Ecm, double fct_eff, const Vec3& sustained,
                           const std::vector<Vec2>& points, const PassiveSteel& bars) {
  require(std::isfinite(fct_eff) && fct_eff >= 0, "service state: fct_eff >= 0");
  ServiceState out;
  out.state = prestress_service_state(s, Ecm, sustained, points, bars);
  out.flexural_tension = out.state.fibre_max.stress;
  out.cracked = out.flexural_tension > fct_eff;
  if (out.cracked) out.state = cracked_service_state(s, Ecm, sustained, points, bars);
  return out;
}

RcPolySection without_prestress(const RcPolySection& s) {
  RcPolySection out = s;
  for (RcPolySection::Tendon& t : out.tendons) t.eps_p0 = 0.0;
  return out;
}

CrackingMoment cracking_moment(const RcPolySection& s, double Ecm, double fct, double theta, const Vec3& sustained,
                               const PassiveSteel& bars) {
  require(std::isfinite(fct) && fct >= 0, "cracking moment: fct >= 0");
  require(std::isfinite(theta), "cracking moment: finite direction theta");
  // The uncracked state is linear: the stress of a fibre under `sustained`
  // plus M along theta is its stress under `sustained` plus M times its
  // stress under the unit moment on the section without initial strains.
  const StrainPlane base = prestress_service_state(s, Ecm, sustained, {}, bars).plane;
  RcPolySection slack = without_prestress(s);
  for (RcPolySection::Bar& b : slack.bars) b.eps0 = 0.0;
  const StrainPlane unit =
      prestress_service_state(slack, Ecm, Vec3(0.0, std::cos(theta), std::sin(theta)), {}, bars).plane;
  CrackingMoment out;
  out.moment = std::numeric_limits<double>::infinity();
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      const double stress = Ecm * base.at(v.x(), v.y()), rate = Ecm * unit.at(v.x(), v.y());
      require(stress <= fct, "cracking moment: a fibre is already past fct under the sustained action");
      if (rate <= 0) continue;  // the moment compresses this fibre
      const double moment = (fct - stress) / rate;
      if (moment < out.moment) out = {moment, {v.x(), v.y(), stress}, 1.0 / rate};
    }
  return out;
}

double brittle_failure_steel(double M_rep, double z, double fyk, double Ap, double fptk) {
  require(std::isfinite(M_rep) && M_rep >= 0, "brittle failure steel: M_rep >= 0");
  require(positive(z) && positive(fyk), "brittle failure steel: z and fyk > 0");
  require(std::isfinite(Ap) && Ap >= 0 && (Ap == 0 || positive(fptk)),
          "brittle failure steel: Ap >= 0, with fptk > 0 when tendons are counted");
  constexpr double kShare = 0.4, kCap = 500e6;  // dsigma_p = min(0,4 fptk, 500 MPa)
  const double tendons = Ap > 0 ? Ap * std::min(kShare * fptk, kCap) : 0.0;
  return std::max(0.0, (M_rep / z - tendons) / fyk);
}

double ultimate_lever_arm(const RcPolySection& s, double omega, double NEd, double theta) {
  const SectionDesign r = moment_resistance(s, omega, NEd, theta);
  require(r.feasible, "ultimate lever arm: no ultimate plane carries NEd with a moment along theta");
  // the concrete alone under the ultimate plane
  RcPolySection concrete = s;
  concrete.bars.clear();
  concrete.tendons.clear();
  const Vec3 C = section_resultants(concrete, 0.0, r.plane);
  Vec3 T = Vec3::Zero();
  for (const RcPolySection::Bar& b : s.bars) {
    const double force = omega * b.w * steel_stress(s.steel, b.eps0 + r.plane.at(b.y, b.z));
    if (force > 0) T += force * Vec3(1.0, b.z, -b.y);
  }
  require(C(0) < 0 && T(0) > 0, "ultimate lever arm: no compressed concrete or no stretched bar on the ultimate plane");
  // points of application (y, z) = (-Mz / N, My / N)
  const Vec2 at_c(-C(2) / C(0), C(1) / C(0)), at_t(-T(2) / T(0), T(1) / T(0));
  return (at_t - at_c).norm();
}

CrackWidth crack_width_of(double sigma_s, double As, double Ap, double Ac_eff, double h_minus_x, double Es,
                          const CrackWidthParams& p) {
  require(positive(sigma_s) && positive(Ac_eff) && positive(Es), "crack width: sigma_s, Ac_eff and Es > 0");
  require(std::isfinite(As) && As >= 0 && std::isfinite(Ap) && Ap >= 0, "crack width: steel areas >= 0");
  require(positive(p.fct_eff) && positive(p.alpha_e) && positive(p.bar_dia),
          "crack width: fct_eff, alpha_e and bar_dia > 0");
  require(std::isfinite(p.cover) && p.cover >= 0 && std::isfinite(p.spacing) && p.spacing >= 0,
          "crack width: cover and spacing >= 0");
  require(positive(p.k1) && positive(p.k2) && positive(p.k3) && positive(p.k4) && positive(p.kt),
          "crack width: factors k1 to k4 and kt > 0");
  require(std::isfinite(p.xi) && p.xi >= 0, "crack width: bond ratio xi >= 0");
  CrackWidth r;
  r.sigma_s = sigma_s;
  r.Ac_eff = Ac_eff;
  r.As = As;
  r.Ap = Ap;
  if (Ap > 0 && p.xi > 0) {
    require(As == 0 || positive(p.tendon_dia), "crack width: tendons counted beside bars need their diameter");
    r.xi1 = std::sqrt(As > 0 ? p.xi * p.bar_dia / p.tendon_dia : p.xi);  // (7.5)
  }
  r.rho_p_eff = (As + (p.xi1_squared ? r.xi1 * r.xi1 : r.xi1) * Ap) / Ac_eff;  // (7.10)
  r.d_eps = 0.6 * sigma_s / Es;                                                // floor of (7.9)
  if (r.rho_p_eff > 0)
    r.d_eps = std::max(r.d_eps, (sigma_s - p.kt * p.fct_eff / r.rho_p_eff * (1 + p.alpha_e * r.rho_p_eff)) / Es);
  r.wide_spacing = r.rho_p_eff == 0 || p.spacing > 5 * (p.cover + p.bar_dia / 2);
  if (r.wide_spacing) {
    require(positive(h_minus_x), "crack width: stretched depth h - x > 0 for Expr. (7.14)");
    r.s_r_max = 1.3 * h_minus_x;  // (7.14)
  } else {
    r.s_r_max = p.k3 * p.cover + p.k1 * p.k2 * p.k4 * p.bar_dia / r.rho_p_eff;  // (7.11)
  }
  r.wk = r.s_r_max * r.d_eps;  // (7.8)
  return r;
}

namespace {

// Signed area of a ring within the half-plane p.u >= level (Sutherland-
// Hodgman against one line: exact for any simple ring, a void comes out
// negative).
double area_beyond(const std::vector<Vec2>& ring, const Vec2& u, double level) {
  std::vector<Vec2> kept;
  for (std::size_t i = 0; i < ring.size(); ++i) {
    const Vec2 &a = ring[i], &b = ring[(i + 1) % ring.size()];
    const double da = a.dot(u) - level, db = b.dot(u) - level;
    if (da >= 0) kept.push_back(a);
    if ((da >= 0) != (db >= 0)) kept.push_back(a + (b - a) * (da / (da - db)));
  }
  double twice = 0;
  for (std::size_t i = 0; i < kept.size(); ++i) {
    const Vec2 &a = kept[i], &b = kept[(i + 1) % kept.size()];
    twice += a.x() * b.y() - b.x() * a.y();
  }
  return twice / 2;
}

}  // namespace

CrackWidth crack_width(const RcPolySection& s, double Ecm, const Vec3& action, const PassiveSteel& bars,
                       const CrackWidthParams& p) {
  const PrestressServiceState state = cracked_service_state(s, Ecm, action, {}, bars);
  constexpr double kInf = std::numeric_limits<double>::infinity();
  double eps_min = kInf, eps_max = -kInf;
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      eps_min = std::min(eps_min, state.plane.at(v.x(), v.y()));
      eps_max = std::max(eps_max, state.plane.at(v.x(), v.y()));
    }
  require(eps_max > 0, "crack width: the action leaves the section compressed everywhere");
  require(eps_min < 0, "crack width: the action stretches the section everywhere (use crack_width_of)");
  const Vec2 gradient(-state.plane.kz, state.plane.ky);
  const double slope = gradient.norm();  // strain per unit depth, > 0 between the two signs
  const Vec2 u = gradient / slope;       // towards the stretched face
  double low = kInf, high = -kInf;
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      low = std::min(low, v.dot(u));
      high = std::max(high, v.dot(u));
    }
  CrackWidth geometry;
  geometry.h = high - low;
  geometry.x = -state.plane.e0 / slope - low;  // the strain e0 + slope (p.u) vanishes there

  // the stretched bonded steels: bars first, the bonded tendons without any
  struct Stretched {
    double level = 0, area = 0, stress = 0;
    bool tendon = false;
  };
  std::vector<Stretched> steels;
  // (a steel is read by the strain of the plane at its level: the stress
  // variation from the zero strain of the concrete there, §7.3.4(1))
  for (std::size_t i = 0; i < state.bar_stress.size(); ++i) {
    const double strain = state.plane.at(s.bars[i].y, s.bars[i].z);
    if (strain > 0)
      steels.push_back({Vec2(s.bars[i].y, s.bars[i].z).dot(u), bars.omega * s.bars[i].w, s.steel.Es * strain, false});
  }
  const bool bars_stretched = !steels.empty();
  for (const RcPolySection::Tendon& t : s.tendons) {
    const double strain = state.plane.at(t.y, t.z);
    if (t.bonded && strain > 0) steels.push_back({Vec2(t.y, t.z).dot(u), t.Ap, s.strand.E * strain, true});
  }
  double area = 0, moment = 0, sigma_s = 0;
  for (const Stretched& b : steels) {
    if (b.tendon == bars_stretched) continue;  // the bars rule when some are stretched, else the tendons
    area += b.area;
    moment += b.area * b.level;
    sigma_s = std::max(sigma_s, b.stress);
  }
  require(area > 0, "crack width: no bonded steel in the stretched zone");
  geometry.d = moment / area - low;
  // h / 2 never rules here: with a compressed zone, (h - x) / 3 < h / 3
  geometry.hc_ef = std::min(2.5 * (geometry.h - geometry.d), (geometry.h - geometry.x) / 3);
  const double band = high - geometry.hc_ef;
  double Ac_eff = 0, As = 0, Ap = 0;
  for (const auto& ring : s.rings) Ac_eff += area_beyond(ring, u, band);
  for (const Stretched& b : steels)
    if (b.level >= band) (b.tendon ? Ap : As) += b.area;
  CrackWidthParams q = p;
  if (q.alpha_e == 0) q.alpha_e = s.steel.Es / Ecm;
  CrackWidth r = crack_width_of(sigma_s, As, Ap, Ac_eff, geometry.h - geometry.x, s.steel.Es, q);
  r.h = geometry.h;
  r.d = geometry.d;
  r.x = geometry.x;
  r.hc_ef = geometry.hc_ef;
  return r;
}

}  // namespace keta::ec2

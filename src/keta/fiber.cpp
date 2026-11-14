#include "keta/fiber.h"

#include <limits>
#include <optional>

namespace keta {

UniaxialLaw UniaxialLaw::elastic(double E) {
  UniaxialLaw l;
  l.kind = Kind::Elastic;
  l.E = E;
  return l;
}
UniaxialLaw UniaxialLaw::steel(double E, double fy, double eps_u) {
  UniaxialLaw l;
  l.kind = Kind::SteelBilinear;
  l.E = E;
  l.fy = fy;
  l.eps_u = eps_u;
  return l;
}
UniaxialLaw UniaxialLaw::concrete_pr(double fcd, double eps_c2, double eps_cu2, double n) {
  UniaxialLaw l;
  l.kind = Kind::ConcretePR;
  l.fcd = fcd;
  l.eps_c2 = eps_c2;
  l.eps_cu2 = eps_cu2;
  l.n_exp = n;
  return l;
}
UniaxialLaw UniaxialLaw::steel_elliptic(double E, double fp, double fy) {
  UniaxialLaw l;
  l.kind = Kind::SteelElliptic;
  l.E = E;
  l.fp = fp;
  l.fy = fy;
  return l;
}

double UniaxialLaw::sigma(double eps) const {
  switch (kind) {
    case Kind::Elastic: return E * eps;
    case Kind::SteelBilinear: {
      double ey = fy / E;
      if (eps_u > 0 && std::fabs(eps) > eps_u) return 0.0;  // ruptured
      if (eps > ey) return fy;
      if (eps < -ey) return -fy;
      return E * eps;
    }
    case Kind::ConcretePR: {
      double x = -eps;  // compression positive internally
      if (x <= 0) return 0.0;         // no tension
      if (x >= eps_cu2) return 0.0;   // crushed
      if (x >= eps_c2) return -fcd;   // plateau (sign: compression stress < 0)
      return -fcd * (1.0 - std::pow(1.0 - x / eps_c2, n_exp));
    }
    case Kind::SteelElliptic: {
      // EN 1993-1-2 §3.2: linear to eps_p, ellipse to eps_y = 0.02, plateau.
      double s = eps < 0 ? -1.0 : 1.0;
      double x = std::fabs(eps);
      double eps_p = fp / E, eps_y = 0.02;
      if (x <= eps_p) return s * E * x;
      if (x >= eps_y) return s * fy;
      double c = (fy - fp) * (fy - fp) / ((eps_y - eps_p) * E - 2 * (fy - fp));
      double a2 = (eps_y - eps_p) * (eps_y - eps_p + c / E);
      double b2 = c * (eps_y - eps_p) * E + c * c;
      double val = fp - c + (std::sqrt(b2) / std::sqrt(a2)) *
                                std::sqrt(std::max(0.0, a2 - (eps_y - x) * (eps_y - x)));
      return s * val;
    }
    case Kind::ConcreteTensionSoftening: {
      if (eps <= 0) {
        const double x = -eps;
        if (x >= eps_cu2) return 0.0;
        if (x >= eps_c2) return -fcd;
        return -fcd * (1.0 - std::pow(1.0 - x / eps_c2, n_exp));
      }
      const double eps_cr = ft / E_t;
      if (eps <= eps_cr) return E_t * eps;
      if (eps < eps_tu) return ft * (eps_tu - eps) / (eps_tu - eps_cr);
      return 0.0;
    }
    case Kind::SteelHardening:
    case Kind::Strand: {
      const double s = eps < 0 ? -1.0 : 1.0, x = std::fabs(eps);
      const double y = kind == Kind::Strand ? fpd : fy, ey = y / E;
      if (x <= ey) return E * eps;
      const double xb = eps_ud > 0 ? std::min(x, eps_ud) : x;  // clamp past the bound (finite probes)
      return s * (y + hardening * (xb - ey));
    }
    case Kind::Sargin: {
      if (eps >= 0) return 0.0;
      const double x = -eps;
      if (x > eps_cu1) return 0.0;
      const double eta = x / eps_c1;
      return -fcm * (k_sargin * eta - eta * eta) / (1.0 + (k_sargin - 2.0) * eta);
    }
    case Kind::Bfup: {
      if (eps <= 0) {
        const double x = -eps;
        if (x > eps_cud) return 0.0;
        return -std::min(E * x, fcd);
      }
      const double e_el = ft_el / E;
      if (eps <= e_el) return E * eps;
      if (eps <= eps_lim) return ft_el + (ft_f - ft_el) * (eps - e_el) / (eps_lim - e_el);
      return 0.0;
    }
  }
  return 0;
}

double UniaxialLaw::tangent(double eps) const {
  switch (kind) {
    case Kind::Elastic: return E;
    case Kind::SteelBilinear: {
      double ey = fy / E;
      if (eps_u > 0 && std::fabs(eps) > eps_u) return 0.0;
      return std::fabs(eps) <= ey ? E : 0.0;
    }
    case Kind::ConcretePR: {
      double x = -eps;
      if (x <= 0 || x >= eps_c2) return 0.0;
      return fcd * n_exp / eps_c2 * std::pow(1.0 - x / eps_c2, n_exp - 1.0);
    }
    case Kind::SteelElliptic: {
      double h = 1e-7;
      return (sigma(eps + h) - sigma(eps - h)) / (2 * h);
    }
    case Kind::ConcreteTensionSoftening: {
      if (eps <= 0) {
        const double x = -eps;
        if (x >= eps_c2) return 0.0;
        return fcd * n_exp / eps_c2 * std::pow(1.0 - x / eps_c2, n_exp - 1.0);
      }
      const double eps_cr = ft / E_t;
      if (eps <= eps_cr) return E_t;
      if (eps < eps_tu) return -ft / (eps_tu - eps_cr);
      return 0.0;
    }
    case Kind::SteelHardening:
    case Kind::Strand: {
      const double x = std::fabs(eps), ey = (kind == Kind::Strand ? fpd : fy) / E;
      if (x <= ey) return E;
      return eps_ud > 0 && x >= eps_ud ? 0.0 : hardening;
    }
    case Kind::Sargin: {
      if (eps >= 0) return 0.0;
      const double x = -eps;
      if (x > eps_cu1) return 0.0;
      const double eta = x / eps_c1, k = k_sargin, d = 1.0 + (k - 2.0) * eta;
      // d sigma / d eps = fcm / eps_c1 (k - 2 eta - (k - 2) eta^2) / d^2: zero at the
      // peak eta = 1, positive before, negative on the descending branch
      return fcm / eps_c1 * (k - 2.0 * eta - (k - 2.0) * eta * eta) / (d * d);
    }
    case Kind::Bfup: {
      if (eps <= 0) return -eps < fcd / E ? E : 0.0;
      const double e_el = ft_el / E;
      if (eps <= e_el) return E;
      if (eps <= eps_lim) return (ft_f - ft_el) / (eps_lim - e_el);
      return 0.0;
    }
  }
  return 0;
}

double UniaxialLaw::sup_tangent() const {
  switch (kind) {
    case Kind::Elastic:
    case Kind::SteelBilinear:
    case Kind::SteelElliptic: return E;
    case Kind::ConcretePR: return fcd * n_exp / eps_c2;
    case Kind::ConcreteTensionSoftening: return std::max(fcd * n_exp / eps_c2, E_t);
    case Kind::SteelHardening:
    case Kind::Strand:
    case Kind::Bfup: return E;
    case Kind::Sargin: return fcm * k_sargin / eps_c1;  // the initial slope 1.05 Ecm
  }
  return E;
}

double UniaxialLaw::eps_min() const {
  switch (kind) {
    case Kind::ConcretePR:
    case Kind::ConcreteTensionSoftening: return -eps_cu2;
    case Kind::SteelBilinear: return eps_u > 0 ? -eps_u : -std::numeric_limits<double>::infinity();
    case Kind::SteelHardening:
    case Kind::Strand: return eps_ud > 0 ? -eps_ud : -std::numeric_limits<double>::infinity();
    case Kind::Sargin: return -eps_cu1;
    case Kind::Bfup: return -eps_cud;
    default: return -std::numeric_limits<double>::infinity();
  }
}

double UniaxialLaw::eps_max() const {
  switch (kind) {
    case Kind::SteelBilinear: return eps_u > 0 ? eps_u : std::numeric_limits<double>::infinity();
    case Kind::SteelHardening:
    case Kind::Strand: return eps_ud > 0 ? eps_ud : std::numeric_limits<double>::infinity();
    case Kind::Bfup: return eps_lim;
    default: return std::numeric_limits<double>::infinity();
  }
}

UniaxialLaw UniaxialLaw::concrete_tension_softening(double fcd, double eps_c2, double eps_cu2, double n, double E_t,
                                                    double ft, double Gf, double lc) {
  for (double v : {fcd, eps_c2, eps_cu2, n, E_t, ft, Gf, lc})
    require(std::isfinite(v) && v > 0, "concrete_tension_softening: every parameter finite and > 0");
  require(eps_cu2 >= eps_c2, "concrete_tension_softening: eps_cu2 >= eps_c2");
  UniaxialLaw l = concrete_pr(fcd, eps_c2, eps_cu2, n);
  l.kind = Kind::ConcreteTensionSoftening;
  l.E_t = E_t;
  l.ft = ft;
  l.eps_tu = ft / E_t + 2.0 * Gf / (ft * lc);
  return l;
}

UniaxialLaw UniaxialLaw::steel_hardening(double E, double fy, double k, double eps_uk, double eps_ud) {
  for (double v : {E, fy, eps_uk, eps_ud})
    require(std::isfinite(v) && v > 0, "steel_hardening: E, fy, eps_uk, eps_ud finite and > 0");
  require(std::isfinite(k) && k > 1, "steel_hardening: k = ft / fy > 1");
  require(eps_uk > fy / E && eps_ud > fy / E, "steel_hardening: eps_uk and eps_ud beyond the yield strain");
  require(eps_ud <= eps_uk, "steel_hardening: eps_ud <= eps_uk");
  UniaxialLaw l;
  l.kind = Kind::SteelHardening;
  l.E = E;
  l.fy = fy;
  l.hardening = (k - 1.0) * fy / (eps_uk - fy / E);
  l.eps_ud = eps_ud;
  return l;
}

UniaxialLaw UniaxialLaw::sargin(double fcm, double Ecm, double eps_c1, double eps_cu1) {
  for (double v : {fcm, Ecm, eps_c1, eps_cu1})
    require(std::isfinite(v) && v > 0, "sargin: fcm, Ecm, eps_c1, eps_cu1 finite and > 0");
  require(eps_cu1 >= eps_c1, "sargin: eps_cu1 >= eps_c1");
  UniaxialLaw l;
  l.kind = Kind::Sargin;
  l.fcm = fcm;
  l.E = Ecm;
  l.eps_c1 = eps_c1;
  l.eps_cu1 = eps_cu1;
  l.k_sargin = 1.05 * Ecm * eps_c1 / fcm;
  // the denominator 1 + (k - 2) eta is linear, 1 at eta = 0: positive over
  // the range iff positive at eps_cu1 (else the stress blows up mid-curve)
  require(1.0 + (l.k_sargin - 2.0) * eps_cu1 / eps_c1 > 0, "sargin: the pole of Expr. (3.14) falls before eps_cu1");
  return l;
}

UniaxialLaw UniaxialLaw::sargin_ec2(double fck) {
  require(std::isfinite(fck) && fck >= 12e6 && fck <= 90e6, "sargin_ec2: fck in [12, 90] MPa (Table 3.1)");
  const double fcm = fck + 8e6, fm = fcm / 1e6;
  const double Ecm = 22e9 * std::pow(fm / 10.0, 0.3);
  const double eps_c1 = std::min(0.7 * std::pow(fm, 0.31), 2.8) / 1000.0;
  const double eps_cu1 = fck < 50e6 ? 3.5e-3 : (2.8 + 27.0 * std::pow((98.0 - fm) / 100.0, 4)) / 1000.0;
  return sargin(fcm, Ecm, eps_c1, eps_cu1);
}

UniaxialLaw UniaxialLaw::strand(double Ep, double fp01k, double fpk, double eps_uk, double eps_ud, bool inclined,
                                double gamma_s) {
  for (double v : {Ep, fp01k, fpk, eps_uk, gamma_s})
    require(std::isfinite(v) && v > 0, "strand: Ep, fp01k, fpk, eps_uk, gamma_s finite and > 0");
  require(fpk > fp01k, "strand: fpk > fp0,1k");
  UniaxialLaw l;
  l.kind = Kind::Strand;
  l.E = Ep;
  l.fpd = fp01k / gamma_s;
  require(eps_uk > l.fpd / Ep, "strand: eps_uk beyond the elastic limit");
  if (inclined) {
    require(std::isfinite(eps_ud) && eps_ud > l.fpd / Ep && eps_ud <= eps_uk,
            "strand: eps_ud between the elastic limit and eps_uk");
    l.hardening = (fpk - fp01k) / gamma_s / (eps_uk - l.fpd / Ep);
    l.eps_ud = eps_ud;
  }
  return l;
}

UniaxialLaw UniaxialLaw::bfup_t3(double fck, double fcm, double Ecm, double fctm, double K, double fctk_el,
                                 double fctfk, double eps_u_lim, double gamma_c, double gamma_cf, double alpha_cc) {
  for (double v : {fck, fcm, Ecm, fctm, K, fctk_el, fctfk, eps_u_lim, gamma_c, gamma_cf, alpha_cc})
    require(std::isfinite(v) && v > 0, "bfup_t3: every parameter finite and > 0");
  UniaxialLaw l;
  l.kind = Kind::Bfup;
  l.E = Ecm;
  l.fcd = alpha_cc * fck / gamma_c;
  l.eps_cud = (1.0 + 14.0 * fctm / (K * fcm)) * l.fcd / Ecm;
  l.ft_el = fctk_el / gamma_cf;
  l.ft_f = fctfk / (gamma_cf * K);
  l.eps_lim = eps_u_lim;
  require(eps_u_lim > l.ft_el / Ecm, "bfup_t3: eps_u,lim beyond the elastic tensile strain");
  return l;
}

void LawRegistry::add(const std::string& name, const UniaxialLaw& law) {
  require(!name.empty(), "law registry: an empty name");
  require(!laws_.count(name), "law registry: '" + name + "' is already registered");
  laws_[name] = law;
}

const UniaxialLaw& LawRegistry::get(const std::string& name) const {
  const auto it = laws_.find(name);
  require(it != laws_.end(), "law registry: no law '" + name + "'");
  return it->second;
}

std::vector<std::string> LawRegistry::names() const {
  std::vector<std::string> v;
  for (const auto& [k, l] : laws_) v.push_back(k);
  return v;
}

double FiberSection::area() const {
  double a = 0;
  for (const Fiber& f : fibers) a += f.A;
  return a;
}

Vec3 FiberSection::resultants(const Vec3& e) const {
  Vec3 r = Vec3::Zero();
  for (const Fiber& f : fibers) {
    double eps = f.eps0 + e(0) + e(1) * f.z - e(2) * f.y;
    double s = laws[static_cast<std::size_t>(f.law)].sigma(eps);
    r(0) += s * f.A;
    r(1) += s * f.A * f.z;
    r(2) -= s * f.A * f.y;
  }
  for (std::size_t i = 0; i < tendons.size(); ++i) {
    const Tendon& t = tendons[i];
    const double s = tendon_stress(i, e);
    r(0) += s * t.A;
    r(1) += s * t.A * t.z;
    r(2) -= s * t.A * t.y;
  }
  return r;
}

double FiberSection::tendon_stress(std::size_t i, const Vec3& e) const {
  const Tendon& t = tendons.at(i);
  const UniaxialLaw& law = laws[static_cast<std::size_t>(t.law)];
  return law.sigma(t.bonded ? t.eps_p0 + e(0) + e(1) * t.z - e(2) * t.y : t.eps_p0);
}

void FiberSection::add_tendon(double y, double z, double A, const UniaxialLaw& law, double eps_p0, bool bonded) {
  using K = UniaxialLaw::Kind;
  require(law.kind == K::Elastic || law.kind == K::SteelBilinear || law.kind == K::SteelHardening ||
              law.kind == K::Strand,
          "fiber tendon: the law must be elastic, steel, steel_hardening or strand");
  require(std::isfinite(y) && std::isfinite(z) && std::isfinite(A) && A > 0 && std::isfinite(eps_p0),
          "fiber tendon: finite position, area > 0 and finite initial strain");
  require(eps_p0 < law.eps_max(), "fiber tendon: initial strain at or past the law's bound");
  laws.push_back(law);
  tendons.push_back({y, z, A, static_cast<int>(laws.size()) - 1, eps_p0, bonded});
}

Mat3 FiberSection::tangent(const Vec3& e) const {
  Mat3 k = Mat3::Zero();
  for (const Fiber& f : fibers) {
    double eps = f.eps0 + e(0) + e(1) * f.z - e(2) * f.y;
    double Et = laws[static_cast<std::size_t>(f.law)].tangent(eps);
    Vec3 w(1.0, f.z, -f.y);
    k += Et * f.A * (w * w.transpose());
  }
  for (const Tendon& t : tendons) {
    if (!t.bonded) continue;  // a frozen force
    const double eps = t.eps_p0 + e(0) + e(1) * t.z - e(2) * t.y;
    const double Et = laws[static_cast<std::size_t>(t.law)].tangent(eps);
    Vec3 w(1.0, t.z, -t.y);
    k += Et * t.A * (w * w.transpose());
  }
  return k;
}

namespace {

struct StrainRange {
  double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
};

// Narrows `range` to the axial strains e0 that keep eps = e0 + offset
// strictly within [eps_min, eps_max]: a few ulps inside, a law may already
// have dropped its stress at the bound itself.
void keep_within_bounds(StrainRange& range, double eps_min, double eps_max, double offset) {
  constexpr double kUlps = 8 * std::numeric_limits<double>::epsilon();
  if (std::isfinite(eps_min))
    range.lo = std::max(range.lo, eps_min - offset + kUlps * (std::fabs(eps_min) + std::fabs(offset)));
  if (std::isfinite(eps_max))
    range.hi = std::min(range.hi, eps_max - offset - kUlps * (std::fabs(eps_max) + std::fabs(offset)));
}

// The axial strains e0 that keep every fibre and bonded tendon of `s` within
// the strain bounds of its law under the curvatures (ky, kz). `eps_steel_max`
// also bounds the bilinear and elliptic steel fibres, whose plateau has no
// end of its own. The bounds are those of the strain the law reads: the
// initial strain of a fibre or of a tendon plus the plane.
StrainRange strains_within_bounds(const FiberSection& s, double ky, double kz, double eps_steel_max) {
  StrainRange range;
  for (const Fiber& f : s.fibers) {
    const UniaxialLaw& law = s.laws[static_cast<std::size_t>(f.law)];
    const bool plateau_steel =
        law.kind == UniaxialLaw::Kind::SteelBilinear || law.kind == UniaxialLaw::Kind::SteelElliptic;
    const double cap = plateau_steel ? eps_steel_max : std::numeric_limits<double>::infinity();
    keep_within_bounds(range, std::max(law.eps_min(), -cap), std::min(law.eps_max(), cap),
                       f.eps0 + ky * f.z - kz * f.y);
  }
  for (const FiberSection::Tendon& t : s.tendons) {
    if (!t.bonded) continue;  // a frozen force, whatever the plane
    const UniaxialLaw& law = s.laws[static_cast<std::size_t>(t.law)];
    keep_within_bounds(range, law.eps_min(), law.eps_max(), t.eps_p0 + ky * t.z - kz * t.y);
  }
  return range;
}

// The largest e0 of `range` with N(e0) = n_target; none when the range does
// not bracket the target. An infinite end is searched geometrically.
// A softening branch (past the peak of Sargin, tension softening, a fibre
// past its bound) makes N non-monotone and the range may hold several roots:
// the physical one is the LARGEST e0 (the least compressed state, on the
// rising branch), so scan down from the top to the first point below the
// target and bisect there.
template <class AxialForce>
std::optional<double> largest_axial_root(const AxialForce& N, double n_target, const StrainRange& range) {
  constexpr double kFirstReach = 0.05, kGrowth = 1.5;
  constexpr int kExpansions = 60, kScan = 256, kBisections = 200;
  double lo = range.lo, hi = range.hi;
  if (!std::isfinite(lo)) {
    const double from = std::min(hi, 0.0);
    double reach = kFirstReach;
    lo = from - reach;
    for (int i = 0; i < kExpansions && N(lo) > n_target; ++i) lo = from - (reach *= kGrowth);
  }
  if (!std::isfinite(hi)) {
    const double from = std::max(range.lo, 0.0);
    double reach = kFirstReach;
    hi = from + reach;
    for (int i = 0; i < kExpansions && N(hi) < n_target; ++i) hi = from + (reach *= kGrowth);
  }
  if (!(N(lo) <= n_target && N(hi) >= n_target)) return std::nullopt;
  const double top = hi, step = (hi - lo) / kScan;
  for (int i = 1; i <= kScan; ++i) {
    const double e = i == kScan ? lo : top - step * i;
    if (N(e) <= n_target) {
      lo = e;
      break;
    }
    hi = e;
  }
  for (int i = 0; i < kBisections; ++i) {
    const double mid = 0.5 * (lo + hi);
    (N(mid) < n_target ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

// The axial equilibrium at the curvatures (ky, kz) that keeps the section
// within its strain bounds (see strains_within_bounds); none past the ultimate
// curvature of this N.
std::optional<double> equilibrium_within_bounds(const FiberSection& s, double n_target, double ky, double kz,
                                                double eps_steel_max) {
  const StrainRange range = strains_within_bounds(s, ky, kz, eps_steel_max);
  if (range.lo > range.hi) return std::nullopt;
  return largest_axial_root([&](double e0) { return s.resultants(Vec3(e0, ky, kz))(0); }, n_target, range);
}

}  // namespace

Vec3 FiberSection::solve_state(const Vec3& target, int max_iter, double tol) const {
  // scale from the sup-tangent stiffness
  Mat3 k0 = Mat3::Zero();
  for (const Fiber& f : fibers) {
    double Es = laws[static_cast<std::size_t>(f.law)].sup_tangent();
    Vec3 w(1.0, f.z, -f.y);
    k0 += Es * f.A * (w * w.transpose());
  }
  for (const Tendon& t : tendons) {
    if (!t.bonded) continue;
    const double Es = laws[static_cast<std::size_t>(t.law)].sup_tangent();
    Vec3 w(1.0, t.z, -t.y);
    k0 += Es * t.A * (w * w.transpose());
  }
  Vec3 scale(std::max(1.0, k0(0, 0) * 1e-4), std::max(1.0, k0(1, 1) * 1e-4),
             std::max(1.0, k0(2, 2) * 1e-4));
  Eigen::FullPivLU<Mat3> lu0(k0);
  require(lu0.isInvertible(), "fiber: degenerate section (all fibers on a line?)");

  // largest component of a residual, each relative to its target (or its floor)
  auto scaled = [&](const Vec3& rho) {
    double err = 0;
    for (int c = 0; c < 3; ++c) err = std::max(err, std::fabs(rho(c)) / std::max(std::fabs(target(c)), scale(c)));
    return err;
  };
  // A law drops its stress past its strain bound (crushing, rupture): a plane
  // out there, the steels left carrying the target alone, is a root of the
  // residual and no state of the section. The search stays where every fibre
  // and bonded tendon is within the bounds of its law.
  constexpr double kNoCap = std::numeric_limits<double>::infinity();
  auto within_bounds = [&](const Vec3& e) {
    const StrainRange range = strains_within_bounds(*this, e(1), e(2), kNoCap);
    return e(0) >= range.lo && e(0) <= range.hi;
  };
  // Start on the plane of the section at its initial moduli: at a zero strain
  // a concrete without tension has no tangent at all, and the first step would
  // be that of the steels alone.
  Vec3 e = lu0.solve(target - resultants(Vec3::Zero()));
  for (int it = 0; it < max_iter; ++it) {
    const Vec3 rho = target - resultants(e);
    const double r0 = scaled(rho);
    if (r0 < tol) return e;
    Mat3 k = tangent(e);
    Vec3 de;
    Eigen::FullPivLU<Mat3> lu(k);
    if (lu.isInvertible() && k.determinant() > 1e-12 * k0.determinant())
      de = lu.solve(rho);
    else
      de = lu0.solve(rho);  // sup-tangent gradient descent (monotone energy)
    // backtracking on the scaled residual, within the strain bounds
    double s = 1.0;
    bool ok = false;
    for (int back = 0; back < 40; ++back) {
      const Vec3 et = e + s * de;
      if (within_bounds(et) && scaled(target - resultants(et)) < r0) {
        e = et;
        ok = true;
        break;
      }
      s *= 0.5;
    }
    if (!ok) fail("fiber: state determination stalled (target outside capacity?)");
  }
  fail("fiber: state determination did not converge (target outside the capacity polytope)");
}

double FiberSection::solve_axial(double n_target, double ky, double kz) const {
  // A law drops its stress past its strain bound (crushing, rupture): N(e0)
  // steps back across the target right below the equilibrium that sits at
  // the bound, and a search over every e0 lands on a post-limit state.
  // Where no fibre is past a bound N has no such step, so an equilibrium
  // there is the physical one and is searched first.
  if (const auto e0 = equilibrium_within_bounds(*this, n_target, ky, kz, std::numeric_limits<double>::infinity()))
    return *e0;
  // past the ultimate curvature of this N: the post-limit equilibrium
  const auto e0 =
      largest_axial_root([&](double e) { return resultants(Vec3(e, ky, kz))(0); }, n_target, StrainRange{});
  require(e0.has_value(), "fiber: axial target unreachable at this curvature");
  return *e0;
}

std::vector<FiberSection::MKappaPoint> FiberSection::m_kappa(
    double n_target, const std::vector<double>& kappas) const {
  std::vector<MKappaPoint> out;
  for (double k : kappas) {
    MKappaPoint p;
    p.kappa = k;
    p.eps0 = solve_axial(n_target, k, 0.0);
    p.moment = resultants(Vec3(p.eps0, k, 0.0))(1);
    out.push_back(p);
  }
  return out;
}

FiberSection FiberSection::rect(double b, double h, int ny, int nz, const UniaxialLaw& law) {
  FiberSection s;
  s.laws.push_back(law);
  for (int i = 0; i < ny; ++i)
    for (int j = 0; j < nz; ++j) {
      Fiber f;
      f.y = -b / 2 + b * (i + 0.5) / ny;
      f.z = -h / 2 + h * (j + 0.5) / nz;
      f.A = b * h / (ny * nz);
      f.law = 0;
      s.fibers.push_back(f);
    }
  return s;
}

void FiberSection::add_rebar(double y, double z, double A, const UniaxialLaw& law) {
  laws.push_back(law);
  fibers.push_back({y, z, A, static_cast<int>(laws.size()) - 1});
}

void FiberSection::add_rebar(double y, double z, double A, const UniaxialLaw& law, double eps0) {
  // strictly inside the bounds of the law (which a NaN or an infinite value is not)
  require(eps0 > law.eps_min() && eps0 < law.eps_max(),
          "fiber rebar: an initial strain finite and within the strain bounds of the law");
  add_rebar(y, z, A, law);
  fibers.back().eps0 = eps0;
}

double FiberSection::fiber_stress(std::size_t i, const Vec3& e) const {
  const Fiber& f = fibers.at(i);
  return laws[static_cast<std::size_t>(f.law)].sigma(f.eps0 + e(0) + e(1) * f.z - e(2) * f.y);
}

}  // namespace keta

namespace keta {

UltimatePoint ultimate_point(const FiberSection& s, double n_target, double alpha, double eps_steel_max) {
  constexpr double kFirstCurvature = 1e-4, kLargestCurvature = 1e3;
  constexpr int kBisections = 60;
  const double cy = std::cos(alpha), cz = std::sin(alpha);
  auto equilibrium = [&](double k) { return equilibrium_within_bounds(s, n_target, k * cy, k * cz, eps_steel_max); };
  double lo = kFirstCurvature;
  require(equilibrium(lo).has_value(), "fiber: axial target out of reach even at a vanishing curvature");
  while (equilibrium(2 * lo)) {
    lo *= 2;
    require(lo < kLargestCurvature, "fiber: no strain bound is reached in this direction (no ultimate point)");
  }
  double hi = 2 * lo;
  for (int i = 0; i < kBisections; ++i) {
    const double mid = 0.5 * (lo + hi);
    (equilibrium(mid) ? lo : hi) = mid;
  }
  UltimatePoint p;
  p.kappa = lo;
  p.eps0 = *equilibrium(lo);
  const Vec3 r = s.resultants(Vec3(p.eps0, lo * cy, lo * cz));
  p.my = r(1);
  p.mz = r(2);
  return p;
}

std::vector<InteractionPoint> interaction_curve(const FiberSection& s, double n_target, int n_alpha,
                                                double eps_steel_max) {
  std::vector<InteractionPoint> out;
  for (int ia = 0; ia < n_alpha; ++ia) {
    const double alpha = 2 * kPi * ia / n_alpha;
    const UltimatePoint p = ultimate_point(s, n_target, alpha, eps_steel_max);
    out.push_back({alpha, p.my, p.mz});
  }
  return out;
}

}  // namespace keta

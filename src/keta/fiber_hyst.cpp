#include "keta/fiber_hyst.h"

#include <limits>

#include <cmath>
#include <optional>

namespace keta {

// ---------------- HystLaw ----------------

HystLaw HystLaw::wrap(const UniaxialLaw& l) {
  HystLaw h;
  h.kind = Kind::Memoryless;
  h.base = l;
  return h;
}

HystLaw HystLaw::steel_kinematic(double E, double fy, double b) {
  require(E > 0 && fy > 0, "hyst: steel_kinematic needs E > 0 and fy > 0");
  require(b >= 0 && b < 1, "hyst: steel_kinematic hardening ratio must be in [0, 1)");
  HystLaw h;
  h.kind = Kind::SteelKinematic;
  h.E = E;
  h.fy = fy;
  h.b = b;
  return h;
}

HystLaw HystLaw::concrete_kj(double fcd, double eps_c2, double eps_cu2, double n) {
  require(fcd > 0 && eps_c2 > 0 && eps_cu2 >= eps_c2 && n >= 1,
          "hyst: concrete_kj parameters out of range");
  HystLaw h;
  h.kind = Kind::ConcreteKJ;
  h.fcd = fcd;
  h.eps_c2 = eps_c2;
  h.eps_cu2 = eps_cu2;
  h.n_exp = n;
  return h;
}

HystLaw HystLaw::steel_mixed(double E, double fy, double H, double iso) {
  require(E > 0 && fy > 0 && H > 0, "hyst: steel_mixed needs E, fy, H > 0");
  require(iso >= 0 && iso <= 1, "hyst: steel_mixed isotropic share must be in [0, 1]");
  HystLaw h;
  h.kind = Kind::SteelMixed;
  h.E = E;
  h.fy = fy;
  h.H = H;
  h.iso = iso;
  return h;
}

HystLaw HystLaw::asym_von_mises(double E, double ft0, double Ht, double fc0, double Hc) {
  require(E > 0 && ft0 > 0 && fc0 > 0 && Ht > 0 && Hc > 0, "hyst: asym_von_mises needs positive E, yields, hardenings");
  HystLaw h;
  h.kind = Kind::AsymVonMises;
  h.E = E;
  h.ft0 = ft0;
  h.Ht = Ht;
  h.fc0 = fc0;
  h.Hc = Hc;
  return h;
}

HystLaw HystLaw::damage_simo_ju(double E, double r0, double A, double B) {
  require(E > 0 && r0 > 0 && A >= 0 && A <= 1 && B >= 0, "hyst: damage_simo_ju parameters out of range");
  HystLaw h;
  h.kind = Kind::DamageSimoJu;
  h.E = E;
  h.r0 = r0;
  h.A = A;
  h.B = B;
  return h;
}

HystLaw HystLaw::plastic_damage_oliver(double E, double fy0, double hh, double r0, double A, double n) {
  require(E > 0 && fy0 > 0 && hh > 0 && r0 > 0 && A >= 0 && n >= 1, "hyst: plastic_damage_oliver parameters out of range");
  HystLaw h;
  h.kind = Kind::PlasticDamageOliver;
  h.E = E;
  h.fy = fy0;
  h.H = hh;
  h.r0 = r0;
  h.A = A;
  h.n_ratio = n;
  return h;
}

namespace {
struct StressState {
  double eps = 0;
  HystState st;
};
}  // namespace

// one stress-driven step of the G10 laws: the strain and the new state
static StressState stress_step(const HystLaw& l, const HystState& s, double sig) {
  StressState o;
  o.st = s;
  const double sg = sig >= 0 ? 1.0 : -1.0;
  switch (l.kind) {
    case HystLaw::Kind::SteelMixed: {
      const double R = l.fy + l.iso * l.H * s.ebar, xi = sig - s.alpha;
      double dg = 0;
      if (std::fabs(xi) > R) dg = (std::fabs(xi) - R) / l.H;
      const double sx = xi >= 0 ? 1.0 : -1.0;
      o.st.eps_p = s.eps_p + dg * sx;
      o.st.alpha = s.alpha + (1 - l.iso) * l.H * dg * sx;
      o.st.ebar = s.ebar + dg;
      o.eps = sig / l.E + o.st.eps_p;
      return o;
    }
    case HystLaw::Kind::AsymVonMises: {
      const double e = std::max(s.ebar, sig > 0 ? (sig - l.ft0) / l.Ht : (-sig - l.fc0) / l.Hc);
      const double gam = (l.fc0 + l.Hc * e) - (l.ft0 + l.Ht * e);
      const double g = (sig + gam / 2) / std::sqrt(sig * sig + gam * gam / 2);
      o.st.ebar = e;
      o.st.eps_p = s.eps_p + (e - s.ebar) * (e > s.ebar ? g : 0.0);
      o.eps = sig / l.E + o.st.eps_p;
      return o;
    }
    case HystLaw::Kind::DamageSimoJu: {
      const double r = std::max({l.r0, s.r, std::fabs(sig) / std::sqrt(l.E)});
      const double d = 1 - l.r0 * (1 - l.A) / r - l.A * std::exp(l.B * (l.r0 - r));
      require(d < 1, "hyst: damage reached 1");
      o.st.r = r;
      o.eps = sig / ((1 - d) * l.E);
      return o;
    }
    case HystLaw::Kind::PlasticDamageOliver: {
      const double beta = sig > 0 ? 1.0 : 1.0 / l.n_ratio;
      const double r = std::max({l.r0, s.r, beta * std::fabs(sig) / std::sqrt(l.E)});
      const double d = 1 - l.r0 / r * std::exp(l.A * (1 - r / l.r0));
      require(d < 1, "hyst: damage reached 1");
      const double e = std::max(s.ebar, (std::fabs(sig) - l.fy) / l.H);
      o.st.r = r;
      o.st.ebar = e;
      o.st.eps_p = s.eps_p + sg * (e - s.ebar);
      o.eps = sig / ((1 - d) * l.E) + o.st.eps_p;
      return o;
    }
    default:
      throw Error("hyst: strain_at is defined for the stress-driven laws (steel_mixed, asym_von_mises, damage)");
  }
}

double HystLaw::strain_at(const HystState& st, double sig) const { return stress_step(*this, st, sig).eps; }
HystState HystLaw::update_stress(const HystState& st, double sig) const { return stress_step(*this, st, sig).st; }

double HystLaw::damage(const HystState& st) const {
  if (kind == Kind::DamageSimoJu) {
    const double r = std::max(r0, st.r);
    return 1 - r0 * (1 - A) / r - A * std::exp(B * (r0 - r));
  }
  if (kind == Kind::PlasticDamageOliver) {
    const double r = std::max(r0, st.r);
    return 1 - r0 / r * std::exp(A * (1 - r / r0));
  }
  return 0;
}

namespace {
// the stress at a trial strain: invert the strictly increasing eps(sigma)
double invert(const HystLaw& l, const HystState& st, double eps) {
  auto F = [&](double s) { return stress_step(l, st, s).eps - eps; };
  const double s0 = F(0);
  if (s0 == 0) return 0;
  double lo = 0, hi = 0, step = l.E * std::max(std::fabs(eps), 1e-12) + 1;
  if (s0 < 0) {
    hi = step;
    while (F(hi) < 0) hi *= 2, require(hi < 1e300, "hyst: no stress for this strain");
  } else {
    lo = -step;
    while (F(lo) > 0) lo *= 2, require(lo > -1e300, "hyst: no stress for this strain");
  }
  double flo = F(lo), fhi = F(hi);
  for (int it = 0; it < 300; ++it) {  // Illinois regula falsi, bracketed
    double s = (lo * fhi - hi * flo) / (fhi - flo);
    if (!(s > lo && s < hi)) s = 0.5 * (lo + hi);
    const double fs = F(s);
    if (fs == 0 || hi - lo <= 1e-15 * std::max(std::fabs(lo), std::fabs(hi))) return s;
    if ((fs < 0) == (flo < 0)) {
      lo = s;
      flo = fs;
      fhi *= 0.5;
    } else {
      hi = s;
      fhi = fs;
      flo *= 0.5;
    }
  }
  return 0.5 * (lo + hi);
}
}  // namespace

double HystLaw::kj_residual(double eps_min) const {
  double eta = std::max(0.0, -eps_min) / eps_c2;
  double ep = eta <= 2.0 ? 0.145 * eta * eta + 0.127 * eta : 0.707 * (eta - 2.0) + 0.834;
  return ep * eps_c2;
}

namespace {

// EC2 parabola-rectangle envelope (compression-negative), identical piecewise
// form to `UniaxialLaw::ConcretePR`.
void concrete_envelope(const HystLaw& l, double eps, double& sig, double& Et) {
  double x = -eps;  // compression positive internally
  if (x <= 0 || x >= l.eps_cu2) {
    sig = 0;
    Et = 0;
    return;
  }
  if (x >= l.eps_c2) {
    sig = -l.fcd;
    Et = 0;
    return;
  }
  sig = -l.fcd * (1.0 - std::pow(1.0 - x / l.eps_c2, l.n_exp));
  Et = l.fcd * l.n_exp / l.eps_c2 * std::pow(1.0 - x / l.eps_c2, l.n_exp - 1.0);
}

}  // namespace

void HystLaw::response_mech(const HystState& st, double eps, double& sig, double& Et) const {
  switch (kind) {
    case Kind::Memoryless:
      sig = base.sigma(eps);
      Et = base.tangent(eps);
      return;
    case Kind::SteelKinematic: {
      // 1-D return mapping, linear kinematic hardening (see header).
      double H = b * E / (1.0 - b);
      double sig_tr = E * (eps - st.eps_p);
      double alpha = H * st.eps_p;
      double xi = sig_tr - alpha;
      double f = std::fabs(xi) - fy;
      if (f <= 0) {
        sig = sig_tr;
        Et = E;
      } else {
        double dg = f / (E + H);
        double sgn = xi >= 0 ? 1.0 : -1.0;
        sig = sig_tr - dg * E * sgn;
        Et = E * H / (E + H);  // = b E
      }
      return;
    }
    case Kind::SteelMixed: {
      const double sig_tr = E * (eps - st.eps_p), xi = sig_tr - st.alpha;
      const double f = std::fabs(xi) - (fy + iso * H * st.ebar);
      if (f <= 0) {
        sig = sig_tr;
        Et = E;
      } else {
        sig = sig_tr - f / (E + H) * E * (xi >= 0 ? 1.0 : -1.0);
        Et = E * H / (E + H);
      }
      return;
    }
    case Kind::AsymVonMises:
    case Kind::DamageSimoJu:
    case Kind::PlasticDamageOliver: {
      sig = invert(*this, st, eps);
      // the tangent: the inverse slope of eps(sigma) on the loading side
      const double ds = 1e-7 * (std::fabs(sig) + 1e-6 * E);
      const double s2 = sig + (sig >= 0 ? ds : -ds);
      const double de = strain_at(st, s2) - strain_at(st, sig);
      Et = de != 0 ? (s2 - sig) / de : E;
      return;
    }
    case Kind::ConcreteKJ: {
      if (eps <= st.eps_min) {  // new compressive extreme: on the envelope
        concrete_envelope(*this, eps, sig, Et);
        return;
      }
      double ep = -kj_residual(st.eps_min);  // residual strain (<= 0)
      if (eps >= ep) {  // beyond the residual towards tension: open crack
        sig = 0;
        Et = 0;
        return;
      }
      // On the unload/reload line towards (eps_min, sigma_env(eps_min)).
      double sig_m, dummy;
      concrete_envelope(*this, st.eps_min, sig_m, dummy);
      double Eun = sig_m / (st.eps_min - ep);  // > 0 (both factors negative)
      sig = Eun * (eps - ep);
      Et = Eun;
      return;
    }
  }
  sig = 0;
  Et = 0;
}

HystState HystLaw::update_mech(const HystState& st, double eps) const {
  HystState out = st;
  switch (kind) {
    case Kind::Memoryless: return out;
    case Kind::SteelKinematic: {
      double H = b * E / (1.0 - b);
      double sig_tr = E * (eps - st.eps_p);
      double xi = sig_tr - H * st.eps_p;
      double f = std::fabs(xi) - fy;
      if (f > 0) out.eps_p += f / (E + H) * (xi >= 0 ? 1.0 : -1.0);
      return out;
    }
    case Kind::ConcreteKJ:
      out.eps_min = std::min(st.eps_min, eps);
      return out;
    case Kind::SteelMixed: {
      const double sig_tr = E * (eps - st.eps_p), xi = sig_tr - st.alpha;
      const double f = std::fabs(xi) - (fy + iso * H * st.ebar);
      if (f > 0) {
        const double dg = f / (E + H), sx = xi >= 0 ? 1.0 : -1.0;
        out.eps_p += dg * sx;
        out.alpha += (1 - iso) * H * dg * sx;
        out.ebar += dg;
      }
      return out;
    }
    case Kind::AsymVonMises:
    case Kind::DamageSimoJu:
    case Kind::PlasticDamageOliver:
      return update_stress(st, invert(*this, st, eps));
  }
  return out;
}


namespace {
// the creep step of a law: the creep strain at the end of the step is
// A + B s(sigma_end) (see fiber_hyst.h)
struct CreepStep {
  double A = 0, B = 0;
};
CreepStep creep_step(const HystLaw& l, const HystState& st) {
  CreepStep c;
  // the shrinkage at the end of the step (an imposed strain, in A)
  c.A += st.sh_mode == 1 ? st.sh_drv
                          : (st.cr_dt > 0 && !l.shrink.empty() ? l.shrink_at(st.age + st.cr_dt) - l.shrink_at(st.age0) : st.sh);
  for (std::size_t i = 0; i < l.creep_a.size(); ++i) {
    if (st.cr_dt <= 0) {
      c.A += st.cr_k[i];
      continue;
    }
    const double x = st.cr_dt / l.creep_tau[i], E = std::exp(-x);
    const double lin = x > 1e-8 ? 1 - (1 - E) / x : x / 2;  // 1 - (1 - e^-x) / x
    c.A += st.cr_k[i] * E + l.creep_a[i] * st.cr_s * (1 - E) - l.creep_a[i] * st.cr_s * lin;
    c.B += l.creep_a[i] * lin;
  }
  return c;
}
// the driving stress s = g(sigma) sigma and its derivative
void creep_drive(const HystLaw& l, double sig, double& s, double& ds) {
  double g = 1, dg = 0;
  if (l.creep_fck > 0 && -sig > 0.45 * l.creep_fck) {
    g = std::exp(1.5 * (-sig / l.creep_fck - 0.45));
    dg = -1.5 / l.creep_fck * g;
  }
  s = g * sig;
  ds = g + sig * dg;
}
// the end-of-step stress: h(sigma) = sigma - law(eps - A - B s(sigma)) = 0,
// increasing in sigma (Et >= 0, B >= 0, s' > 0): safeguarded Newton
double creep_solve(const HystLaw& l, const HystState& st, double eps, const CreepStep& c, double& Et, double& ds) {
  double sig, s;
  l.response_mech(st, eps - c.A - c.B * st.cr_s, sig, Et);  // start: the committed driving stress
  if (c.B == 0) {
    creep_drive(l, sig, s, ds);
    return sig;
  }
  double lo = -std::numeric_limits<double>::infinity(), hi = std::numeric_limits<double>::infinity();
  for (int it = 0; it < 200; ++it) {
    creep_drive(l, sig, s, ds);
    double sm;
    l.response_mech(st, eps - c.A - c.B * s, sm, Et);
    const double h = sig - sm;
    if (h == 0) break;
    (h > 0 ? hi : lo) = sig;
    double nx = sig - h / (1 + Et * c.B * ds);
    if (!(nx > lo && nx < hi)) {
      if (std::isfinite(lo) && std::isfinite(hi)) nx = 0.5 * (lo + hi);
      else nx = sig - h;  // (h' >= 1 away from the kinks: a damped step)
    }
    if (std::fabs(nx - sig) <= 1e-15 * std::max(std::fabs(sig), 1.0)) {
      sig = nx;
      break;
    }
    sig = nx;
  }
  creep_drive(l, sig, s, ds);
  double sm;
  l.response_mech(st, eps - c.A - c.B * s, sm, Et);
  return sig;
}
}  // namespace

void HystLaw::response(const HystState& st, double eps, double& sig, double& Et) const {
  if (!creeps() && st.sh_mode == 0) return response_mech(st, eps, sig, Et);
  const CreepStep c = creep_step(*this, st);
  double ds = 1;
  sig = creep_solve(*this, st, eps, c, Et, ds);
  Et = Et / (1 + Et * c.B * ds);
}

HystState HystLaw::update(const HystState& st, double eps) const {
  if (!creeps() && st.sh_mode == 0) return update_mech(st, eps);
  require(creep_a.size() == creep_tau.size() && creep_a.size() <= 6, "fiber creep: at most 6 Kelvin units, a and tau paired");
  for (double tau : creep_tau) require(tau > 0, "fiber creep: retardation times > 0");
  const CreepStep c = creep_step(*this, st);
  double Et = 0, ds = 1, s = 0;
  const double sig = creep_solve(*this, st, eps, c, Et, ds);
  creep_drive(*this, sig, s, ds);
  HystState out = update_mech(st, eps - c.A - c.B * s);
  if (st.cr_dt > 0)
    for (std::size_t i = 0; i < creep_a.size(); ++i) {
      const double x = st.cr_dt / creep_tau[i], E = std::exp(-x);
      const double lin = x > 1e-8 ? 1 - (1 - E) / x : x / 2;
      out.cr_k[i] = st.cr_k[i] * E + creep_a[i] * st.cr_s * (1 - E) + creep_a[i] * (s - st.cr_s) * lin;
    }
  out.cr_s = s;
  if (st.sh_mode == 1) out.sh = st.sh_drv;
  else if (st.cr_dt > 0 && !shrink.empty()) out.sh = shrink_at(st.age + st.cr_dt) - shrink_at(st.age0);
  return out;
}

double HystLaw::shrink_at(double age) const {
  require(!shrink.empty(), "fiber shrinkage: an empty table");
  if (age <= shrink.front().first) return shrink.front().second;
  for (std::size_t i = 1; i < shrink.size(); ++i)
    if (age <= shrink[i].first) {
      const auto& a = shrink[i - 1];
      const auto& b = shrink[i];
      require(b.first > a.first, "fiber shrinkage: ages must increase");
      return a.second + (b.second - a.second) * (age - a.first) / (b.first - a.first);
    }
  return shrink.back().second;
}

double HystLaw::sup_tangent() const {
  switch (kind) {
    case Kind::Memoryless: return base.sup_tangent();
    case Kind::SteelKinematic:
    case Kind::SteelMixed:
    case Kind::AsymVonMises:
    case Kind::DamageSimoJu:
    case Kind::PlasticDamageOliver: return E;
    case Kind::ConcreteKJ:
      // The unload/reload slope exceeds the initial envelope tangent
      // E0 = fcd n / eps_c2 by up to (2 - eta) / (2 (0.873 - 0.145 eta)) with
      // eta -> 0, i.e. ~1.15 E0; 1.25 E0 is a safe supremum.
      return 1.25 * fcd * n_exp / eps_c2;
  }
  return E;
}

// ---------------- HystFiberSection ----------------

HystFiberSection HystFiberSection::wrap(const FiberSection& s) {
  require(s.tendons.empty(), "HystFiberSection::wrap: prestressing tendons are not carried by the hysteretic section");
  for (const Fiber& f : s.fibers)
    require(f.eps0 == 0, "HystFiberSection::wrap: the initial strain of a fibre is not carried by the hysteretic section");
  HystFiberSection h;
  h.fibers = s.fibers;
  h.laws.reserve(s.laws.size());
  for (const UniaxialLaw& l : s.laws) h.laws.push_back(HystLaw::wrap(l));
  return h;
}

HystFiberSection HystFiberSection::rect(double b, double h, int ny, int nz, const HystLaw& law) {
  HystFiberSection s;
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

void HystFiberSection::add_rebar(double y, double z, double A, const HystLaw& law) {
  laws.push_back(law);
  fibers.push_back({y, z, A, static_cast<int>(laws.size()) - 1});
}

Vec3 HystFiberSection::resultants(const Vec3& e, const std::vector<HystState>& st) const {
  Vec3 r = Vec3::Zero();
  for (std::size_t i = 0; i < fibers.size(); ++i) {
    const Fiber& f = fibers[i];
    double eps = e(0) + e(1) * f.z - e(2) * f.y;
    double s, Et;
    laws[static_cast<std::size_t>(f.law)].response(st[i], eps, s, Et);
    r(0) += s * f.A;
    r(1) += s * f.A * f.z;
    r(2) -= s * f.A * f.y;
  }
  return r;
}

Mat3 HystFiberSection::tangent(const Vec3& e, const std::vector<HystState>& st) const {
  Mat3 k = Mat3::Zero();
  for (std::size_t i = 0; i < fibers.size(); ++i) {
    const Fiber& f = fibers[i];
    double eps = e(0) + e(1) * f.z - e(2) * f.y;
    double s, Et;
    laws[static_cast<std::size_t>(f.law)].response(st[i], eps, s, Et);
    Vec3 w(1.0, f.z, -f.y);
    k += Et * f.A * (w * w.transpose());
  }
  return k;
}

Mat3 HystFiberSection::sup_stiffness() const {
  Mat3 k0 = Mat3::Zero();
  for (const Fiber& f : fibers) {
    double Es = laws[static_cast<std::size_t>(f.law)].sup_tangent();
    Vec3 w(1.0, f.z, -f.y);
    k0 += Es * f.A * (w * w.transpose());
  }
  return k0;
}

Vec3 HystFiberSection::solve_state(const Vec3& target, const std::vector<HystState>& st,
                                   const Vec3& e_start, int max_iter, double tol) const {
  require(st.size() == fibers.size(), "hyst: state vector size mismatch");
  Mat3 k0 = sup_stiffness();
  Vec3 scale(std::max(1.0, k0(0, 0) * 1e-4), std::max(1.0, k0(1, 1) * 1e-4),
             std::max(1.0, k0(2, 2) * 1e-4));
  Eigen::FullPivLU<Mat3> lu0(k0);
  require(lu0.isInvertible(), "hyst: degenerate section (all fibers on a line?)");

  // Scaled 2-NORM of the residual — the line-search metric MUST be a norm
  // that decreases when any component improves while the others are exactly
  // constant: on a crack-open branch a whole residual component can be
  // exactly invariant along the search direction until another fiber branch
  // engages, and a max-component metric then stalls on exact ties.
  //
  // ABSOLUTE NOISE FLOOR per component: a component whose target AND residual
  // both sit below 1e-10 of the component scale counts as converged. Without
  // it, a roundoff-seeded nano-moment target (e.g. Mz ~ 1e-9 N.m leaking into
  // the out-of-plane component of a plane problem) would demand crossing a
  // MACROSCOPIC Karsan-Jirsa dead zone to be matched exactly — an engagement
  // problem no local iteration solves, for a quantity of no physical meaning.
  auto rerr = [&](const Vec3& e) {
    Vec3 rho = target - resultants(e, st);
    double s2 = 0;
    for (int c = 0; c < 3; ++c) {
      double floor_c = 1e-10 * scale(c);
      if (std::fabs(target(c)) < floor_c && std::fabs(rho(c)) < floor_c) continue;
      double x = rho(c) / std::max(std::fabs(target(c)), scale(c));
      s2 += x * x;
    }
    return std::sqrt(s2);
  };
  // Backtracking line search along `de`: walk DOWN the halving sequence to
  // the local optimum, not to the first decrease — an overshooting direction
  // (regularized null space of a plateaued/cracked section) shows a first,
  // marginal decrease at a wildly-too-long step; the residual keeps improving
  // as s shrinks towards the physical engagement scale, and we stop at the
  // first worsening past a found improvement.
  auto backtrack = [&](Vec3& e, double r0, const Vec3& de, int budget) {
    double s = 1.0;
    double best_err = r0;
    Vec3 best_e = e;
    bool found = false;
    for (int back = 0; back < budget; ++back) {
      Vec3 et = e + s * de;
      double r1 = rerr(et);
      if (r1 < best_err) {
        best_err = r1;
        best_e = et;
        found = true;
      } else if (found) {
        break;  // past the local optimum along the halving sequence
      }
      s *= 0.5;
    }
    if (found) e = best_e;
    return found;
  };
  // GRACEFUL FLOOR: a deeply plastified state chases the last decades of the
  // tolerance along a nearly flat response (plateau "valley crawl", linear
  // rate ~0.93) and can stagnate at a solver-precision floor far above `tol`
  // yet far below any structural relevance. When the iteration stops
  // improving there, return the best point; targets that are genuinely
  // outside the capacity stagnate at a relative error orders of magnitude
  // above this ceiling (~1e-3 and up) and still fail. The element-level
  // determination applies the same rule to its own residual.
  constexpr double kFloorCeiling = 1e-5;
  auto attempt = [&](Vec3 e) -> std::optional<Vec3> {
    double best_err = std::numeric_limits<double>::max();
    Vec3 best_e = e;
    int stall = 0;
    for (int it = 0; it < max_iter; ++it) {
      Vec3 rho = target - resultants(e, st);
      double err = rerr(e);
      if (err < tol) return e;
      if (err < 0.99 * best_err) {
        best_err = err;
        best_e = e;
        stall = 0;
      } else if (++stall >= 10) {
        break;
      }
      // Newton on the REGULARIZED tangent k + 1e-6 k0: identical to plain
      // Newton where k is well conditioned (the perturbation is 1e-6 of the
      // sup stiffness), and where k is rank-deficient (cracked RC: steel-only
      // rank 2; plateaued EPP: rank 0) the k0 part fills the null space with
      // a direction whose overshoot the walk-to-the-optimum backtracking
      // scales down to the engagement threshold — a pure-k0 direction there
      // would creep with a linear rate ~1 and never meet the tolerance. The
      // k0 direction stays as the monotone fallback on the next line.
      Mat3 kreg = tangent(e, st) + 1e-6 * k0;
      Eigen::FullPivLU<Mat3> lu(kreg);
      bool moved = lu.isInvertible() && backtrack(e, err, lu.solve(rho), 40);
      if (!moved && !backtrack(e, err, lu0.solve(rho), 40)) break;
    }
    if (best_err < kFloorCeiling) return best_e;
    return std::nullopt;
  };
  if (auto r = attempt(e_start)) return *r;
  if (e_start.squaredNorm() > 0)  // warm start dead-ended: retry cold
    if (auto r = attempt(Vec3::Zero())) return *r;
  fail("hyst: state determination did not converge (target outside the capacity polytope?)");
}

void HystFiberSection::commit(const Vec3& e, std::vector<HystState>& st) const {
  require(st.size() == fibers.size(), "hyst: state vector size mismatch");
  for (std::size_t i = 0; i < fibers.size(); ++i) {
    const Fiber& f = fibers[i];
    double eps = e(0) + e(1) * f.z - e(2) * f.y;
    st[i] = laws[static_cast<std::size_t>(f.law)].update(st[i], eps);
  }
}

}  // namespace keta

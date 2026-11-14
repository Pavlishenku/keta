#include "keta/soil_springs.h"

#include <algorithm>
#include <cmath>

namespace keta::soil {

namespace {
constexpr double kDeg45 = kPi / 4;

// 8-point Gauss-Legendre on [0, 1]
constexpr double kGx[8] = {0.019855071751231856, 0.10166676129318664, 0.2372337950418355, 0.4082826787521751,
                           0.5917173212478249,   0.7627662049581645,  0.8983332387068134, 0.9801449282487681};
constexpr double kGw[8] = {0.05061426814518813, 0.11119051722668724, 0.15685332293894363, 0.18134189168918100,
                           0.18134189168918100, 0.15685332293894363, 0.11119051722668724, 0.05061426814518813};

template <class F>
double gauss(double a, double b, F f) {
  double s = 0;
  for (int k = 0; k < 8; ++k) s += kGw[k] * f(a + (b - a) * kGx[k]);
  return s * (b - a);
}

// odd-symmetric law through (y_k, p_k), y ascending from 0, flat beyond the last
LinkLaw odd_law(const std::vector<double>& y, const std::vector<double>& p) {
  LinkLaw l;
  const std::size_t n = y.size();
  const double ext = y.back() > 0 ? y.back() : 1.0;
  l.pts.emplace_back(-y.back() - ext, -p.back());
  for (std::size_t k = n; k-- > 1;) l.pts.emplace_back(-y[k], -p[k]);
  for (std::size_t k = 0; k < n; ++k) l.pts.emplace_back(y[k], p[k]);
  l.pts.emplace_back(y.back() + ext, p.back());
  return l;
}
}  // namespace

double ka_rankine(double phi) { return std::pow(std::tan(kDeg45 - phi / 2), 2); }
double kp_rankine(double phi) { return std::pow(std::tan(kDeg45 + phi / 2), 2); }
double k0_jaky(double phi) { return 1 - std::sin(phi); }
double ka_coulomb(double phi, double delta) {
  const double r = std::sqrt(std::sin(phi + delta) * std::sin(phi) / std::cos(delta));
  return std::pow(std::cos(phi), 2) / (std::cos(delta) * (1 + r) * (1 + r));
}
double kp_coulomb(double phi, double delta) {
  const double s = std::sin(phi + delta) * std::sin(phi) / std::cos(delta);
  require(s < 1, "kp_coulomb: delta too large for the Coulomb passive wedge");
  const double r = std::sqrt(s);
  return std::pow(std::cos(phi), 2) / (std::cos(delta) * (1 - r) * (1 - r));
}

Id add_ground_spring(Model& m, Id node, const Vec3& dir, const LinkLaw& law) {
  require(dir.norm() > 0, "add_ground_spring: zero direction");
  const Vec3 x = dir.normalized();
  const Vec3 seed = std::fabs(x.z()) < 0.9 ? Vec3(0, 0, 1) : Vec3(1, 0, 0);
  const Vec3 y = seed.cross(x).normalized();
  const Vec3 z = x.cross(y);
  const Vec3 p = m.node(node).x;
  const Id g = m.add_node(p.x(), p.y(), p.z());
  m.add_support(g, {true, true, true, true, true, true});
  const Id l = m.add_elastic_link(g, node, Mat6::Zero());
  ElasticLink& e = m.links.back();
  Mat3 R;
  R.row(0) = x.transpose();
  R.row(1) = y.transpose();
  R.row(2) = z.transpose();
  e.axes = R;
  e.law[0] = law;
  return l;
}

LinkLaw law_sum(const LinkLaw& a, const LinkLaw& b) {
  std::vector<double> xs;
  for (const auto& q : a.pts) xs.push_back(q.first);
  for (const auto& q : b.pts) xs.push_back(q.first);
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  LinkLaw s;
  for (double x : xs) s.pts.emplace_back(x, a.force(x) + b.force(x));
  return s;
}

LinkLaw law_scaled(const LinkLaw& a, double f) {
  require(f >= 0, "law_scaled: negative factor");
  LinkLaw s = a;
  for (auto& q : s.pts) q.second *= f;
  return s;
}

// ---------------------------------------------------------- earth pressure
EarthPoint earth_pressure_at(const EarthLayer& l, double y) {
  const double depth = l.top - y;
  const double sv = l.sv_top + l.gamma * depth;
  EarthPoint e;
  e.p0 = l.K0 * sv;
  e.pa = std::max(0.0, l.Ka * sv - l.Kac * l.c);
  e.pp = l.Kp * sv + l.Kpc * l.c;
  e.kh = l.kh_A + l.kh_B * (depth > 0 ? std::pow(depth, l.kh_kappa) : 0.0);
  return e;
}

LinkLaw earth_pressure_law(double P0, double Pa, double Pp, double K) {
  require(K > 0, "earth pressure: the subgrade stiffness must be > 0");
  const double tol = 1e-12 * std::max({std::fabs(P0), std::fabs(Pa), std::fabs(Pp), 1e-300});
  require(Pa <= P0 + tol && P0 <= Pp + tol, "earth pressure: needs active <= at rest <= passive");
  const double dp = -(std::max(Pp, P0) - P0) / K, da = (P0 - std::min(Pa, P0)) / K;
  LinkLaw l;
  const double ext = std::max({-dp, da, 1e-3});
  if (dp < 0) {
    l.pts.emplace_back(dp - ext, -Pp);
    l.pts.emplace_back(dp, -Pp);
  } else {
    l.pts.emplace_back(-ext, -P0);
  }
  l.pts.emplace_back(0.0, -P0);
  if (da > 0) {
    l.pts.emplace_back(da, -Pa);
    l.pts.emplace_back(da + ext, -Pa);
  } else {
    l.pts.emplace_back(ext, -P0);
  }
  return l;
}

std::vector<Id> add_earth_pressure(Model& m, const std::vector<Id>& wall, const Vec3& away, const Vec3& up,
                                   const std::vector<EarthLayer>& layers, double width, EarthDistribution dist) {
  require(wall.size() >= 2, "add_earth_pressure: a chain of two nodes at least");
  require(width > 0, "add_earth_pressure: width must be > 0");
  require(up.norm() > 0 && away.norm() > 0, "add_earth_pressure: zero direction");
  const Vec3 u = up.normalized();
  for (const EarthLayer& l : layers) require(l.top > l.bottom, "add_earth_pressure: layer top must lie above its bottom");
  std::vector<Vec3> X;
  for (Id n : wall) X.push_back(m.node(n).x);
  std::vector<Id> out(wall.size(), 0);
  for (std::size_t i = 0; i < wall.size(); ++i) {
    LinkLaw law;
    bool any = false;
    for (const EarthLayer& l : layers) {
      double P0 = 0, Pa = 0, Pp = 0, K = 0;
      for (int side : {-1, 1}) {
        const long o = static_cast<long>(i) + side;
        if (o < 0 || o >= static_cast<long>(wall.size())) continue;
        const Vec3 &xi = X[i], &xo = X[static_cast<std::size_t>(o)];
        const double L = (xo - xi).norm(), yi = xi.dot(u), yo = xo.dot(u);
        // t in [0, 1] from node i; N_i = 1 - t; the part inside the layer
        double t0 = 0, t1 = 1;
        if (yo != yi) {
          double ta = (l.top - yi) / (yo - yi), tb = (l.bottom - yi) / (yo - yi);
          if (ta > tb) std::swap(ta, tb);
          t0 = std::max(0.0, ta);
          t1 = std::min(1.0, tb);
        } else if (yi > l.top || yi < l.bottom) {
          continue;
        }
        if (t1 <= t0) continue;
        if (dist == EarthDistribution::Lumped) {
          // the half-segment next to node i: t in [0, 1/2]
          const double h1 = std::min(t1, 0.5);
          if (h1 <= t0) continue;
          const EarthPoint e = earth_pressure_at(l, std::clamp(yi, l.bottom, l.top));
          const double len = L * (h1 - t0);
          P0 += e.p0 * len;
          Pa += e.pa * len;
          Pp += e.pp * len;
          K += e.kh * len;
          continue;
        }
        auto at = [&](double t) { return earth_pressure_at(l, yi + t * (yo - yi)); };
        // split at the active cut-off (Ka sv - Kac c = 0), linear in t
        std::vector<double> cuts{t0, t1};
        const double fa0 = l.Ka * (l.sv_top + l.gamma * (l.top - yi)) - l.Kac * l.c;
        const double fa1 = l.Ka * (l.sv_top + l.gamma * (l.top - yo)) - l.Kac * l.c;
        if ((fa0 < 0) != (fa1 < 0) && fa1 != fa0) {
          const double ts = fa0 / (fa0 - fa1);
          if (ts > t0 && ts < t1) cuts.insert(cuts.begin() + 1, ts);
        }
        for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
          const double a = cuts[k], b = cuts[k + 1];
          P0 += L * gauss(a, b, [&](double t) { return at(t).p0 * (1 - t); });
          Pa += L * gauss(a, b, [&](double t) { return at(t).pa * (1 - t); });
          Pp += L * gauss(a, b, [&](double t) { return at(t).pp * (1 - t); });
          K += L * gauss(a, b, [&](double t) { return at(t).kh * (1 - t); });
        }
      }
      if (K <= 0 && P0 == 0 && Pp == 0) continue;
      const LinkLaw li = earth_pressure_law(width * P0, width * Pa, width * Pp, width * K);
      law = any ? law_sum(law, li) : li;
      any = true;
    }
    if (any) out[i] = add_ground_spring(m, wall[i], away, law);
  }
  return out;
}

// ---------------------------------------------------------- piles
std::vector<double> tributary_lengths(const Model& m, const std::vector<Id>& chain) {
  std::vector<double> L(chain.size(), 0.0);
  for (std::size_t i = 0; i + 1 < chain.size(); ++i) {
    const double l = (m.node(chain[i + 1]).x - m.node(chain[i]).x).norm();
    L[i] += l / 2;
    L[i + 1] += l / 2;
  }
  return L;
}

double pu_clay(double c, double gamma_eff, double X, double D, double J) {
  require(c > 0 && D > 0 && X >= 0, "p-y clay: c, D > 0 and depth X >= 0");
  return std::min(3 * c + gamma_eff * X + J * c * X / D, 9 * c) * D;
}

SandCoefficients api_sand_coefficients(double phi) {
  require(phi > 0 && phi < kPi / 2, "p-y sand: friction angle in (0, 90) degrees");
  const double a = phi / 2, b = kDeg45 + phi / 2, K0 = 0.4, Ka = ka_rankine(phi);
  const double tb = std::tan(b), tbp = std::tan(b - phi), tp = std::tan(phi), ta = std::tan(a);
  SandCoefficients s;
  s.C1 = K0 * tp * std::sin(b) / (tbp * std::cos(a)) + tb * tb * ta / tbp + K0 * tb * (tp * std::sin(b) - ta);
  s.C2 = tb / tbp - Ka;
  s.C3 = Ka * (std::pow(tb, 8) - 1) + K0 * tp * std::pow(tb, 4);
  return s;
}

double pu_sand(double phi, double gamma_eff, double X, double D) {
  require(D > 0 && X >= 0 && gamma_eff >= 0, "p-y sand: D > 0, depth X >= 0");
  const SandCoefficients s = api_sand_coefficients(phi);
  return std::min((s.C1 * X + s.C2 * D) * gamma_eff * X, s.C3 * D * gamma_eff * X);
}

LinkLaw py_soft_clay(double c, double gamma_eff, double X, double D, double eps50, double J) {
  require(eps50 > 0, "p-y soft clay: eps50 must be > 0");
  const double pu = pu_clay(c, gamma_eff, X, D, J), yc = 2.5 * eps50 * D;
  const double r[6] = {0, 0.1, 0.3, 1.0, 3.0, 8.0}, q[6] = {0, 0.23, 0.33, 0.50, 0.72, 1.00};
  std::vector<double> y, p;
  for (int k = 0; k < 6; ++k) {
    y.push_back(r[k] * yc);
    p.push_back(q[k] * pu);
  }
  return odd_law(y, p);
}

LinkLaw py_stiff_clay(double c, double gamma_eff, double X, double D, double eps50, int n) {
  require(eps50 > 0 && n >= 2, "p-y stiff clay: eps50 > 0, n >= 2");
  const double pu = pu_clay(c, gamma_eff, X, D, 0.5), y50 = 2.5 * eps50 * D;
  std::vector<double> y, p;
  for (int k = 0; k <= n; ++k) {  // equal force steps on p / pu = 0.5 (y / y50)^(1/4)
    const double s = static_cast<double>(k) / n;
    y.push_back(16 * std::pow(s, 4) * y50);
    p.push_back(s * pu);
  }
  return odd_law(y, p);
}

LinkLaw py_sand(double phi, double gamma_eff, double X, double D, double k, bool cyclic, int n) {
  require(k > 0 && n >= 2, "p-y sand: k > 0, n >= 2");
  const double pu = pu_sand(phi, gamma_eff, X, D);
  if (pu <= 0) return LinkLaw{{{-1.0, 0.0}, {1.0, 0.0}}};
  const double A = cyclic ? 0.9 : std::max(0.9, 3.0 - 0.8 * X / D);
  const double s = k * X / (A * pu);  // p = A pu tanh(s y)
  std::vector<double> y, p;
  for (int i = 0; i <= n; ++i) {  // equal steps of s y up to 4 (tanh 4 = 0.99933), flat beyond
    const double a = 4.0 * i / n;
    y.push_back(a / s);
    p.push_back(A * pu * std::tanh(a));
  }
  return odd_law(y, p);
}

LinkLaw tz_clay(double tmax_line, double D) {
  require(tmax_line >= 0 && D > 0, "t-z clay: tmax >= 0, D > 0");
  const double r[6] = {0, 0.0016, 0.0031, 0.0057, 0.0080, 0.0100}, q[6] = {0, 0.30, 0.50, 0.75, 0.90, 1.00};
  std::vector<double> z, t;
  for (int k = 0; k < 6; ++k) {
    z.push_back(r[k] * D);
    t.push_back(q[k] * tmax_line);
  }
  return odd_law(z, t);
}

LinkLaw tz_sand(double tmax_line) {
  require(tmax_line >= 0, "t-z sand: tmax >= 0");
  return odd_law({0.0, 0.00254}, {0.0, tmax_line});
}

LinkLaw qz_tip(double Qp, double D) {
  require(Qp >= 0 && D > 0, "q-z: Qp >= 0, D > 0");
  const double r[6] = {0, 0.002, 0.013, 0.042, 0.073, 0.100}, q[6] = {0, 0.25, 0.50, 0.75, 0.90, 1.00};
  LinkLaw l;
  l.pts.emplace_back(-0.2 * D, -Qp);
  for (int k = 5; k >= 1; --k) l.pts.emplace_back(-r[k] * D, -q[k] * Qp);
  l.pts.emplace_back(0.0, 0.0);
  l.pts.emplace_back(0.1 * D, 0.0);
  return l;
}

double api_alpha(double c, double sv_eff) {
  require(c > 0 && sv_eff > 0, "api_alpha: c and sv' must be > 0");
  const double psi = c / sv_eff;
  return std::min(1.0, psi <= 1 ? 0.5 / std::sqrt(psi) : 0.5 * std::pow(psi, -0.25));
}

}  // namespace keta::soil

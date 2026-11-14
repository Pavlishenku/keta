// Contour integrals of a uniaxial law over a polygonal section — see
// contour_integrals.h.
#include "keta/contour_integrals.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace keta::contour {

namespace {

// Gauss-Legendre on [0, 1] (Newton on P_n): exact to degree 2 n - 1 = 31. The
// integrands of a polynomial law are of degree 4 at most along an edge, 14
// after the change of variable of add_piece.
constexpr int kGauss = 16;
struct GaussRule {
  std::array<double, kGauss> x{}, w{};
};
const GaussRule& gauss_rule() {
  static const GaussRule rule = [] {
    GaussRule r;
    for (int i = 0; i < kGauss; ++i) {
      double t = std::cos(kPi * (i + 0.75) / (kGauss + 0.5)), dp = 1;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1, p1 = t;
        for (int k = 2; k <= kGauss; ++k) {
          const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = kGauss * (t * p1 - p0) / (t * t - 1);
        const double dt = p1 / dp;
        t -= dt;
        if (std::fabs(dt) < 1e-16) break;
      }
      r.x[static_cast<std::size_t>(i)] = 0.5 * (1 - t);
      r.w[static_cast<std::size_t>(i)] = 1.0 / ((1 - t * t) * dp * dp);
    }
    return r;
  }();
  return rule;
}

// The contour integrals of a plane: int sigma {1, s, t} dA, then
// int Et {1, s, t, s^2, s t, t^2} dA.
using Moments = Eigen::Matrix<double, 9, 1>;

// One piece of boundary, from the parameter a to b: `curve` gives (s, t,
// ds / dparameter). `sign` is +1 on a ring, -1 round a displaced disc.
template <class Curve>
void add_piece(Moments& m, const Law& law, double e0, double gradient, const Curve& curve, double a, double b,
               double sign) {
  // parameter a + (b - a) xi^2 (3 - 2 xi): its derivative vanishes at both ends
  const GaussRule& rule = gauss_rule();
  for (std::size_t q = 0; q < kGauss; ++q) {
    const double xi = rule.x[q];
    const std::array<double, 3> c = curve(a + (b - a) * xi * xi * (3 - 2 * xi));
    const double s = c[0], t = c[1], eps = e0 + gradient * s;
    const double w = -sign * (b - a) * 6 * xi * (1 - xi) * rule.w[q] * c[2];
    const double sg = w * law.sigma(eps);
    m(0) += sg * t;
    m(1) += sg * s * t;
    m(2) += sg * t * t / 2;
    if (!law.tangent) continue;  // resultants alone
    const double et = w * law.tangent(eps);
    m(3) += et * t;
    m(4) += et * s * t;
    m(5) += et * t * t / 2;
    m(6) += et * s * s * t;
    m(7) += et * s * t * t / 2;
    m(8) += et * t * t * t / 3;
  }
}

// `cuts` sorted, each interval between two of them integrated on its own.
template <class Curve>
void add_curve(Moments& m, const Law& law, double e0, double gradient, const Curve& curve, std::vector<double>& cuts,
               double sign) {
  std::sort(cuts.begin(), cuts.end());
  for (std::size_t k = 0; k + 1 < cuts.size(); ++k)
    if (cuts[k + 1] > cuts[k]) add_piece(m, law, e0, gradient, curve, cuts[k], cuts[k + 1], sign);
}

}  // namespace

Response response(const std::vector<std::vector<Vec2>>& rings, const std::vector<Disc>& discs, const Law& law,
                  const Vec3& e) {
  // gradient of eps in (y, z): (-kz, ky); a component at the round-off of the
  // other is snapped to zero (an axis-parallel gradient stays axis-parallel)
  Vec2 grad(-e(2), e(1));
  const double raw = grad.norm();
  if (std::fabs(grad.x()) < 1e-14 * raw) grad.x() = 0;
  if (std::fabs(grad.y()) < 1e-14 * raw) grad.y() = 0;
  const double G = grad.norm(), e0 = e(0);
  const Vec2 u = G > 0 ? Vec2(grad / G) : Vec2(1, 0), v(-u.y(), u.x());

  Moments m = Moments::Zero();
  std::vector<double> cuts;
  for (const auto& ring : rings)
    for (std::size_t i = 0; i < ring.size(); ++i) {
      const Vec2 &p = ring[i], &q = ring[(i + 1) % ring.size()];
      const double sp = p.dot(u), tp = p.dot(v), sq = q.dot(u), tq = q.dot(v);
      if (sp == sq) continue;  // ds = 0 along the edge
      cuts = {0.0, 1.0};
      const double ep = e0 + G * sp, eq = e0 + G * sq;
      if (ep != eq)
        for (double kink : law.kinks) {
          const double at = (kink - ep) / (eq - ep);
          if (at > 0 && at < 1) cuts.push_back(at);
        }
      add_curve(
          m, law, e0, G,
          [&](double l) { return std::array<double, 3>{sp + l * (sq - sp), tp + l * (tq - tp), sq - sp}; }, cuts, 1.0);
    }
  for (const Disc& disc : discs) {
    // the circle counter-clockwise from the point of largest s; the strain is
    // monotonic on each half, a kink cuts each at most once
    const double r = std::sqrt(disc.area / kPi), sc = disc.centre.dot(u), tc = disc.centre.dot(v);
    cuts = {0.0, kPi / 2, kPi, 3 * kPi / 2, 2 * kPi};
    if (G > 0)
      for (double kink : law.kinks) {
        const double x = (kink - (e0 + G * sc)) / (G * r);
        if (std::fabs(x) >= 1) continue;
        cuts.push_back(std::acos(x));
        cuts.push_back(2 * kPi - std::acos(x));
      }
    add_curve(
        m, law, e0, G,
        [&](double a) { return std::array<double, 3>{sc + r * std::cos(a), tc + r * std::sin(a), -r * std::sin(a)}; },
        cuts, -1.0);
  }

  // back to (y, z): p = s u + t v
  Response out;
  const double Sy = m(1) * u.x() + m(2) * v.x(), Sz = m(1) * u.y() + m(2) * v.y();
  out.R = Vec3(m(0), Sz, -Sy);
  const double Ey = m(4) * u.x() + m(5) * v.x(), Ez = m(4) * u.y() + m(5) * v.y();
  const double Eyy = m(6) * u.x() * u.x() + 2 * m(7) * u.x() * v.x() + m(8) * v.x() * v.x();
  const double Ezz = m(6) * u.y() * u.y() + 2 * m(7) * u.y() * v.y() + m(8) * v.y() * v.y();
  const double Eyz = m(6) * u.x() * u.y() + m(7) * (u.x() * v.y() + u.y() * v.x()) + m(8) * v.x() * v.y();
  out.K << m(3), Ez, -Ey, Ez, Ezz, -Eyz, -Ey, -Eyz, Eyy;
  return out;
}

Law concrete_law(const UniaxialLaw& l, const std::string& who) {
  using K = UniaxialLaw::Kind;
  std::vector<double> kinks;
  switch (l.kind) {
    case K::Elastic: break;
    case K::ConcretePR: kinks = {0.0, -l.eps_c2, -l.eps_cu2}; break;
    case K::Sargin: kinks = {0.0, -l.eps_cu1}; break;
    case K::ConcreteTensionSoftening: kinks = {0.0, -l.eps_c2, -l.eps_cu2, l.ft / l.E_t, l.eps_tu}; break;
    default: fail(who + ": the concrete law must be elastic, concrete_pr, sargin or concrete_tension_softening");
  }
  return {[l](double eps) { return l.sigma(eps); }, [l](double eps) { return l.tangent(eps); }, std::move(kinks)};
}

}  // namespace keta::contour

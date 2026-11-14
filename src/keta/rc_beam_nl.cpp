// Material-nonlinear RC beam, exact determinate solution — see rc_beam_nl.h.
#include "keta/rc_beam_nl.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace keta::rcnl {

namespace {
// Root of a continuous f on [lo, hi] with f(lo) f(hi) <= 0: Illinois regula
// falsi with a bisection guard, to the floating-point resolution of x.
template <class F>
double root(F&& f, double lo, double hi) {
  double flo = f(lo), fhi = f(hi);
  if (flo == 0) return lo;
  if (fhi == 0) return hi;
  require((flo < 0) != (fhi < 0), "rc beam: no sign change (no N = 0 state in the admissible strain range)");
  int side = 0;
  for (int it = 0; it < 400; ++it) {
    double x = (lo * fhi - hi * flo) / (fhi - flo);
    if (!(x > lo && x < hi) || it % 4 == 3) x = 0.5 * (lo + hi);  // guard
    if (hi - lo <= 4 * std::numeric_limits<double>::epsilon() * std::max(std::fabs(lo), std::fabs(hi)) ||
        x == lo || x == hi)
      return x;
    const double fx = f(x);
    if (fx == 0) return x;
    if ((fx < 0) == (flo < 0)) {
      lo = x;
      flo = fx;
      if (side == -1) fhi *= 0.5;
      side = -1;
    } else {
      hi = x;
      fhi = fx;
      if (side == 1) flo *= 0.5;
      side = 1;
    }
  }
  return 0.5 * (lo + hi);
}

// natural cubic spline second derivatives through (x, y), tridiagonal solve
std::vector<double> natural_m2(const std::vector<double>& x, const std::vector<double>& y) {
  const std::size_t n = x.size();
  std::vector<double> m(n, 0.0);
  if (n < 3) return m;
  std::vector<double> a(n, 0), b(n, 0), c(n, 0), r(n, 0);
  for (std::size_t i = 1; i + 1 < n; ++i) {
    const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
    a[i] = h0;
    b[i] = 2 * (h0 + h1);
    c[i] = h1;
    r[i] = 6 * ((y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0);
  }
  for (std::size_t i = 2; i + 1 < n; ++i) {  // Thomas
    const double w = a[i] / b[i - 1];
    b[i] -= w * c[i - 1];
    r[i] -= w * r[i - 1];
  }
  for (std::size_t i = n - 2; i >= 1; --i) {
    m[i] = (r[i] - c[i] * m[i + 1]) / b[i];
    if (i == 1) break;
  }
  return m;
}

// 8-point Gauss-Legendre on [a, b]
template <class F>
double gl8(F&& f, double a, double b) {
  static const double xg[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363};
  static const double wg[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
  const double c = 0.5 * (a + b), h = 0.5 * (b - a);
  double s = 0;
  for (int i = 0; i < 4; ++i) s += wg[i] * (f(c - h * xg[i]) + f(c + h * xg[i]));
  return s * h;
}
}  // namespace

// ------------------------------------------------------------ TabulatedLaw

TabulatedLaw::TabulatedLaw(std::vector<Point> pts) {
  require(pts.size() >= 2, "tabulated law: at least two points");
  for (std::size_t i = 1; i < pts.size(); ++i)
    require(pts[i].eps > pts[i - 1].eps, "tabulated law: strains must increase");
  const std::size_t n = pts.size();
  for (const Point& p : pts) {
    x_.push_back(p.eps);
    y_.push_back(p.sig);
  }
  m2_.assign(n, 0.0);
  seg_spline_.assign(n - 1, 0);
  // maximal runs of SPL points, each extended by one neighbour on each side
  std::size_t i = 0;
  while (i < n) {
    if (!pts[i].spline) {
      ++i;
      continue;
    }
    std::size_t j = i;
    while (j + 1 < n && pts[j + 1].spline) ++j;
    const std::size_t a = i > 0 ? i - 1 : i, b = j + 1 < n ? j + 1 : j;
    if (b > a) {
      std::vector<double> xs(x_.begin() + static_cast<long>(a), x_.begin() + static_cast<long>(b) + 1);
      std::vector<double> ys(y_.begin() + static_cast<long>(a), y_.begin() + static_cast<long>(b) + 1);
      const std::vector<double> m = natural_m2(xs, ys);
      for (std::size_t k = a; k <= b; ++k) m2_[k] = m[k - a];
      for (std::size_t k = a; k < b; ++k) seg_spline_[k] = 1;
    }
    i = j + 1;
  }
}

TabulatedLaw TabulatedLaw::symmetric(const std::vector<Point>& t) {
  std::vector<Point> p;
  for (auto it = t.rbegin(); it != t.rend(); ++it)
    if (it->eps > 0) p.push_back({-it->eps, -it->sig, it->spline});
  if (t.empty() || t.front().eps != 0) p.push_back({0.0, 0.0, false});
  for (const Point& q : t) p.push_back(q);
  return TabulatedLaw(p);
}

double TabulatedLaw::sigma(double e) const {
  require(e >= x_.front() * (1 + 1e-12) - 1e-300 && e <= x_.back() * (1 + 1e-12) + 1e-300,
          "tabulated law: strain beyond the table (SSLA TYPE LIM)");
  e = std::clamp(e, x_.front(), x_.back());
  std::size_t k = static_cast<std::size_t>(std::upper_bound(x_.begin(), x_.end(), e) - x_.begin());
  k = std::clamp<std::size_t>(k, 1, x_.size() - 1) - 1;
  const double h = x_[k + 1] - x_[k], t = (e - x_[k]) / h;
  const double lin = y_[k] + t * (y_[k + 1] - y_[k]);
  if (!seg_spline_[k]) return lin;
  const double A = 1 - t, B = t;
  return lin + ((A * A * A - A) * m2_[k] + (B * B * B - B) * m2_[k + 1]) * h * h / 6;
}

double TabulatedLaw::tangent(double e) const {
  e = std::clamp(e, x_.front(), x_.back());
  std::size_t k = static_cast<std::size_t>(std::upper_bound(x_.begin(), x_.end(), e) - x_.begin());
  k = std::clamp<std::size_t>(k, 1, x_.size() - 1) - 1;
  const double h = x_[k + 1] - x_[k], t = (e - x_[k]) / h;
  const double s = (y_[k + 1] - y_[k]) / h;
  if (!seg_spline_[k]) return s;
  const double A = 1 - t, B = t;
  return s + (-(3 * A * A - 1) * m2_[k] + (3 * B * B - 1) * m2_[k + 1]) * h / 6;
}

// ------------------------------------------------------------------ Sargin

double Sargin::sigma(double eps) const {
  if (eps >= 0) return 0.0;
  require(-eps <= eps_cu1 * (1 + 1e-12), "Sargin law: compressive strain beyond eps_cu1");
  const double eta = -eps / eps_c1, kk = k();
  return -fcm * (kk * eta - eta * eta) / (1 + (kk - 2) * eta);
}

void Sargin::integrals(double eta, double& F0, double& F1) const {
  const double kk = k(), c = kk - 2;
  if (std::fabs(c * eta) < 0.25) {
    // series of (k t - t^2) sum (-c t)^n: no cancellation at small strains
    // (the closed form loses ~ eta^-2 relative digits there)
    F0 = F1 = 0;
    double p = 1;  // (-c)^n
    for (int n = 0; n < 80; ++n) {
      const double e2 = std::pow(eta, n + 2);
      const double t0 = p * (kk * e2 / (n + 2) - e2 * eta / (n + 3));
      const double t1 = p * (kk * e2 * eta / (n + 3) - e2 * eta * eta / (n + 4));
      F0 += t0;
      F1 += t1;
      if (std::fabs(t0) <= 1e-18 * std::fabs(F0) && std::fabs(t1) <= 1e-18 * std::fabs(F1)) break;
      p *= -c;
    }
    return;
  }
  if (std::fabs(c) < 1e-12) {
    F0 = kk * eta * eta / 2 - eta * eta * eta / 3;
    F1 = kk * eta * eta * eta / 3 - eta * eta * eta * eta / 4;
    return;
  }
  // (k eta - eta^2) / (1 + c eta) = A eta + B + R / (1 + c eta)
  const double A = -1 / c, B = (kk - A) / c, R = -B;
  const double l = std::log1p(c * eta);
  F0 = A * eta * eta / 2 + B * eta + R / c * l;
  F1 = A * eta * eta * eta / 3 + B * eta * eta / 2 + R * (eta / c - l / (c * c));
}

Sargin sargin_from_fcm(double fcm) {
  const double f = fcm / 1e6;
  Sargin s;
  s.fcm = fcm;
  s.Ecm = 22e9 * std::pow(f / 10, 0.3);
  s.eps_c1 = std::min(0.7 * std::pow(f, 0.31), 2.8) * 1e-3;
  const double fck = f - 8;
  s.eps_cu1 = (fck < 50 ? 3.5 : 2.8 + 27 * std::pow((98 - f) / 100, 4)) * 1e-3;
  return s;
}

// ---------------------------------------------------------- Heft 525 (TS I)

TabulatedLaw ts_steel_heft525(const Heft525& p) {
  require(p.fct_eff > 0 && p.rho_eff > 0 && p.Ec > 0 && p.Es > 0, "tension stiffening: invalid data");
  const double ae = p.Es / p.Ec, ssr = p.fct_eff * (1 / p.rho_eff + ae - 1);
  const double esr1 = p.fct_eff / p.Ec, esr2 = ssr / p.Es, dE = esr2 - esr1, esy = p.fy / p.Es;
  require(1.3 * ssr < p.fy, "tension stiffening: 1.3 sigma_sr must stay below fy");
  require(p.ft > p.fy && p.eps_su > esy, "tension stiffening: bare steel needs ft > fy, eps_su > eps_sy");
  const double s13 = 1.3 * ssr;
  const double e13 = s13 / p.Es - p.beta_t * dE;
  const double ey = esy - p.beta_t * dE;
  const double eu = ey + p.delta_d * (1 - ssr / p.fy) * (p.eps_su - esy);
  std::vector<TabulatedLaw::Point> pts = {{-p.eps_su, -p.ft, false}, {-esy, -p.fy, false}, {0.0, 0.0, false},
                                          {esr1, ssr, false},       {e13, s13, false},     {ey, p.fy, false},
                                          {eu, p.ft, false}};
  return TabulatedLaw(pts);
}

// --------------------------------------------------------------- RcSection

RcSection::RcSection(double b, double h, Sargin concrete, std::vector<Bar> bars)
    : b_(b), h_(h), c_(concrete), bars_(std::move(bars)) {
  require(b > 0 && h > 0 && c_.fcm > 0 && c_.eps_c1 > 0 && c_.eps_cu1 >= c_.eps_c1, "rc section: invalid data");
  for (const Bar& br : bars_) require(br.law != nullptr && br.As > 0, "rc section: bar without law");
  // limit: top fibre at -eps_cu1 with N = 0, or a bar at the end of its table
  // (a bar beyond its table counts as the sign of its side: the root found is
  // then the admissibility boundary, handled below)
  auto top_at_limit = [&](double k) {
    const double e0 = -c_.eps_cu1 + k * h_ / 2;
    for (const Bar& br : bars_) {
      const double e = e0 - k * br.z;
      if (e > br.law->eps_max()) return 1.0;
      if (e < br.law->eps_min()) return -1.0;
    }
    return axial(e0, k);
  };
  double khi = 1e-4;
  auto admissible = [&](double k, double e0) {
    for (const Bar& br : bars_) {
      const double e = e0 - k * br.z;
      if (e < br.law->eps_min() || e > br.law->eps_max()) return false;
    }
    return true;
  };
  // the concrete criterion: g(k) = N(top at -eps_cu1) grows from < 0
  while (top_at_limit(khi) < 0 && khi < 10) khi *= 2;
  require(top_at_limit(khi) >= 0, "rc section: no N = 0 state with the top at eps_cu1");
  const double kc = root(top_at_limit, 0.0, khi);
  // largest curvature in [0, kc] with an admissible N = 0 state (kc itself
  // when the concrete governs; below it when a bar reaches its table end)
  auto ok = [&](double k) {
    try {
      double e0 = 0;
      moment(k, &e0);
      return admissible(k, e0);
    } catch (const Error&) {
      return false;
    }
  };
  if (ok(kc)) {
    klim_ = kc;
  } else {
    double lo = 0, hi = kc;
    for (int it = 0; it < 200 && hi - lo > 4 * std::numeric_limits<double>::epsilon() * hi; ++it) {
      const double mid = 0.5 * (lo + hi);
      (ok(mid) ? lo : hi) = mid;
    }
    klim_ = lo;
  }
  // peak of M(kappa) on [0, klim]: dense scan, then golden-section refinement
  // in the bracket of the largest sample (M is unimodal there)
  const int ns = 2000;
  int imax = 0;
  double Mmax = 0;
  for (int i = 1; i <= ns; ++i) {
    const double M = moment(klim_ * i / ns);
    if (M > Mmax) {
      Mmax = M;
      imax = i;
    }
  }
  if (imax == ns) {
    kpeak_ = klim_;
  } else {
    double a = klim_ * (imax - 1) / ns, b = klim_ * (imax + 1) / ns;
    const double g = 0.5 * (std::sqrt(5.0) - 1);
    double c = b - g * (b - a), d = a + g * (b - a), Mc = moment(c), Md = moment(d);
    for (int it = 0; it < 200 && b - a > 4 * std::numeric_limits<double>::epsilon() * b; ++it) {
      if (Mc > Md) {
        b = d;
        d = c;
        Md = Mc;
        c = b - g * (b - a);
        Mc = moment(c);
      } else {
        a = c;
        c = d;
        Mc = Md;
        d = a + g * (b - a);
        Md = moment(d);
      }
    }
    kpeak_ = 0.5 * (a + b);
  }
  double Mp = 0;
  for (int i = 1; i <= ns; ++i) {
    const double M = moment(kpeak_ * i / ns);
    require(M > Mp, "rc section: M(kappa) is not increasing up to its peak (non-unimodal)");
    Mp = M;
  }
}

double RcSection::axial(double e0, double k) const {
  double N = 0;
  if (k > 0) {
    // compressed zone: z from max(-h/2, e0 / k) to h/2, eta = (k z - e0) / eps_c1
    const double eta_top = (k * h_ / 2 - e0) / c_.eps_c1;
    if (eta_top > 0) {
      const double eta_bot = std::max(0.0, (-k * h_ / 2 - e0) / c_.eps_c1);
      double F0t, F1t, F0b, F1b;
      c_.integrals(eta_top, F0t, F1t);
      c_.integrals(eta_bot, F0b, F1b);
      N -= b_ * c_.fcm * c_.eps_c1 / k * (F0t - F0b);
    }
  } else if (e0 < 0) {
    N += b_ * h_ * c_.sigma(e0);
  }
  for (const Bar& br : bars_) N += br.As * br.law->sigma(e0 - k * br.z);
  return N;
}

double RcSection::moment(double k, double* eps0) const {
  require(k >= 0, "rc section: sagging curvature only");
  // bracket: top at -eps_cu1 .. the tightest bar table end
  const double lo = -c_.eps_cu1 + k * h_ / 2;
  double hi = std::numeric_limits<double>::infinity();
  for (const Bar& br : bars_) hi = std::min(hi, br.law->eps_max() + k * br.z);
  require(hi > lo, "rc section: empty admissible strain range");
  const double e0 = root([&](double e) { return axial(e, k); }, lo, hi);
  if (eps0) *eps0 = e0;
  double M = 0;
  if (k > 0) {
    const double eta_top = (k * h_ / 2 - e0) / c_.eps_c1;
    if (eta_top > 0) {
      const double eta_bot = std::max(0.0, (-k * h_ / 2 - e0) / c_.eps_c1);
      double F0t, F1t, F0b, F1b;
      c_.integrals(eta_top, F0t, F1t);
      c_.integrals(eta_bot, F0b, F1b);
      // M = -int sigma z dA, z = (eps_c1 eta + e0) / k
      M += b_ * c_.fcm * c_.eps_c1 / (k * k) * (c_.eps_c1 * (F1t - F1b) + e0 * (F0t - F0b));
    }
  }
  for (const Bar& br : bars_) M -= br.As * br.law->sigma(e0 - k * br.z) * br.z;
  return M;
}

double RcSection::kappa_of(double M) const {
  if (M <= 0) {
    require(M == 0, "rc section: sagging moments only");
    return 0;
  }
  require(M <= moment(kpeak_) * (1 + 1e-12), "rc section: moment beyond the section's peak");
  if (M >= moment(kpeak_)) return kpeak_;
  return root([&](double k) { return moment(k) - M; }, 0.0, kpeak_);
}

std::vector<double> RcSection::kink_curvatures() const {
  // scan grid: geometric from 1e-9 kappa_peak (the uncracked range is tiny)
  // and uniform up to kappa_peak, sign changes refined to round-off
  std::vector<double> grid{0.0};
  for (int i = 0; i <= 400; ++i) grid.push_back(kpeak_ * std::pow(1e-9, 1.0 - i / 400.0));
  for (int i = 1; i <= 400; ++i) grid.push_back(kpeak_ * i / 400);
  std::sort(grid.begin(), grid.end());
  grid.erase(std::unique(grid.begin(), grid.end()), grid.end());
  std::vector<double> out;
  for (const Bar& br : bars_) {
    for (double ek : br.law->breakpoints()) {
      auto g = [&](double k) {
        double e0;
        moment(k, &e0);
        return e0 - k * br.z - ek;
      };
      double kp = grid[1], gp = g(kp);
      for (std::size_t i = 2; i < grid.size(); ++i) {
        const double kk = grid[i], gk = g(kk);
        if (gp != 0 && gk != 0 && (gp < 0) != (gk < 0)) out.push_back(root(g, kp, kk));
        kp = kk;
        gp = gk;
      }
    }
  }
  std::sort(out.begin(), out.end());
  return out;
}

// --------------------------------------------------- SimplySupportedRcBeam

SimplySupportedRcBeam::SimplySupportedRcBeam(const RcSection& s, double L, double w, double a)
    : s_(s), L_(L), w_(w), a_(a) {
  require(L > 0 && a > 0 && a < L && w >= 0, "rc beam: invalid span / load position");
  for (double k : s.kink_curvatures()) Mk_.push_back(s.moment(k));
}

double SimplySupportedRcBeam::moment(double P, double x) const {
  const double mp = x <= a_ ? x * (L_ - a_) / L_ : a_ * (L_ - x) / L_;
  return w_ * x * (L_ - x) / 2 + P * mp;
}

double SimplySupportedRcBeam::P_limit() const {
  const double Ml = s_.moment(s_.kappa_peak());
  const double mp = a_ * (L_ - a_) / L_;
  const double P = (Ml - w_ * a_ * (L_ - a_) / 2) / mp;
  // the maximum of M(x) must sit under the load: M'(a-) >= 0 >= M'(a+)
  require(w_ * (L_ - 2 * a_) / 2 + P * (L_ - a_) / L_ >= 0 && w_ * (L_ - 2 * a_) / 2 - P * a_ / L_ <= 0,
          "rc beam: the governing section is not under the load");
  return P;
}

double SimplySupportedRcBeam::integrate(double P, double x0, double x1, double tol) const {
  auto f = [&](double x) {
    const double mb = x <= a_ ? x * (L_ - a_) / L_ : a_ * (L_ - x) / L_;
    return s_.kappa_of(moment(P, x)) * mb;
  };
  // adaptive: accept when two halves agree with the whole
  struct Seg {
    double a, b, whole;
    int depth;
  };
  std::vector<Seg> stack{{x0, x1, gl8(f, x0, x1), 0}};
  double sum = 0;
  while (!stack.empty()) {
    const Seg s = stack.back();
    stack.pop_back();
    const double m = 0.5 * (s.a + s.b);
    const double l = gl8(f, s.a, m), r = gl8(f, m, s.b);
    if (std::fabs(l + r - s.whole) <= tol || s.depth > 30) {
      sum += l + r;
    } else {
      stack.push_back({s.a, m, l, s.depth + 1});
      stack.push_back({m, s.b, r, s.depth + 1});
    }
  }
  return sum;
}

double SimplySupportedRcBeam::deflection(double P) const {
  require(P >= 0 && P <= P_limit() * (1 + 1e-12), "rc beam: load beyond the limit load");
  // abscissae where M(x) crosses a kink moment (exact roots of the quadratics)
  std::vector<double> xs{0.0, a_, L_};
  for (double Mk : Mk_) {
    // left: -w/2 x^2 + (wL/2 + P (L-a)/L) x - Mk = 0 ; right: -w/2 x^2 + (wL/2 - P a/L) x + P a - Mk = 0
    const double cl[3] = {-w_ / 2, w_ * L_ / 2 + P * (L_ - a_) / L_, -Mk};
    const double cr[3] = {-w_ / 2, w_ * L_ / 2 - P * a_ / L_, P * a_ - Mk};
    for (int side = 0; side < 2; ++side) {
      const double* c = side == 0 ? cl : cr;
      std::vector<double> r;
      if (c[0] == 0) {
        if (c[1] != 0) r.push_back(-c[2] / c[1]);
      } else {
        const double disc = c[1] * c[1] - 4 * c[0] * c[2];
        if (disc >= 0) {
          const double q = -0.5 * (c[1] + std::copysign(std::sqrt(disc), c[1]));
          if (q != 0) r.push_back(c[2] / q);
          r.push_back(q / c[0]);
        }
      }
      for (double x : r)
        if (side == 0 ? (x > 0 && x < a_) : (x > a_ && x < L_)) xs.push_back(x);
    }
  }
  std::sort(xs.begin(), xs.end());
  // tolerance: round-off level of the whole integral, shared by length
  double scale = 0;
  for (std::size_t i = 0; i + 1 < xs.size(); ++i)
    if (xs[i + 1] > xs[i]) scale += std::fabs(integrate(P, xs[i], xs[i + 1], std::numeric_limits<double>::infinity()));
  double total = 0;
  for (std::size_t i = 0; i + 1 < xs.size(); ++i)
    if (xs[i + 1] > xs[i]) total += integrate(P, xs[i], xs[i + 1], 1e-13 * scale * (xs[i + 1] - xs[i]) / L_);
  return total;
}

}  // namespace keta::rcnl

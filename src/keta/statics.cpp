#include "keta/statics.h"

#include <map>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include "keta/bef.h"
#include "keta/parallel.h"
#include "keta/springmass.h"

namespace keta {

namespace {

// Exact load integrals along a bar for station recovery.
// Sum of distributed + point load components up to abscissa x (component c of
// the LOCAL force), and their first moment about the cut at x.
struct SpanIntegrals {
  const BarSpanLoads& sp;

  // Indexed form (many loads on the bar: tendon kinks and segments): the
  // load is piecewise linear between the sorted segment ends, so prefix
  // values at those breakpoints (force, moment about the breakpoint, couple)
  // and sorted point-load prefix sums answer every query in O(log n) instead
  // of a walk over all the loads — the same integrals, exactly.
  bool indexed = false;
  std::vector<double> bk;                      // breakpoints, ascending
  std::vector<std::array<double, 3>> qa_, qs_;  // q_c(s) = qa + qs (s - bk[i]) on [bk[i], bk[i+1]]
  std::vector<std::array<double, 3>> ma_, ms_;  // couple intensities, same form
  std::vector<std::array<double, 3>> F0, M0, C0;  // prefix at bk[i]: force, moment about bk[i], couple
  std::vector<double> pa;                      // point abscissae, ascending
  std::vector<std::size_t> pord;               // their indices in sp.points
  std::vector<std::array<double, 3>> pF, pFa, pM;  // prefix sums of F, F a, M (pa[0..k) )

  explicit SpanIntegrals(const BarSpanLoads& s) : sp(s) {
    // The prefix sums hold linear intensities: a parabolic segment takes the
    // exact walk below.
    if (sp.segs.size() + sp.points.size() < 16 || has_parabolic(sp.segs)) return;
    indexed = true;
    // points
    std::vector<std::size_t> ord(sp.points.size());
    for (std::size_t i = 0; i < ord.size(); ++i) ord[i] = i;
    std::stable_sort(ord.begin(), ord.end(),
                     [&](std::size_t i, std::size_t j) { return sp.points[i].a < sp.points[j].a; });
    pord = ord;
    pa.resize(ord.size());
    pF.assign(ord.size() + 1, {0, 0, 0});
    pFa.assign(ord.size() + 1, {0, 0, 0});
    pM.assign(ord.size() + 1, {0, 0, 0});
    for (std::size_t k = 0; k < ord.size(); ++k) {
      const PointAt& p = sp.points[ord[k]];
      pa[k] = p.a;
      for (int c = 0; c < 3; ++c) {
        pF[k + 1][c] = pF[k][c] + p.F(c);
        pFa[k + 1][c] = pFa[k][c] + p.F(c) * p.a;
        pM[k + 1][c] = pM[k][c] + p.M(c);
      }
    }
    // distributed: breakpoints and the linear load on each interval
    for (const DistSeg& g : sp.segs)
      if (g.b > g.a) {
        bk.push_back(g.a);
        bk.push_back(g.b);
      }
    std::sort(bk.begin(), bk.end());
    bk.erase(std::unique(bk.begin(), bk.end()), bk.end());
    const std::size_t m = bk.size();
    qa_.assign(m, {0, 0, 0});
    qs_.assign(m, {0, 0, 0});
    ma_.assign(m, {0, 0, 0});
    ms_.assign(m, {0, 0, 0});
    if (m >= 2) {
      // q on [bk_i, bk_i+1] = sum over active segs of qa + sl (s - a); written
      // as A + sl * s per seg and summed by difference arrays, then shifted
      std::vector<std::array<double, 3>> dA(m, {0, 0, 0}), dS(m, {0, 0, 0}), dmA(m, {0, 0, 0}),
          dmS(m, {0, 0, 0});
      for (const DistSeg& g : sp.segs) {
        if (!(g.b > g.a)) continue;
        const std::size_t i0 = static_cast<std::size_t>(std::lower_bound(bk.begin(), bk.end(), g.a) - bk.begin());
        const std::size_t i1 = static_cast<std::size_t>(std::lower_bound(bk.begin(), bk.end(), g.b) - bk.begin());
        const double len = g.b - g.a;
        for (int c = 0; c < 3; ++c) {
          const double sl = (g.qb(c) - g.qa(c)) / len, A = g.qa(c) - sl * g.a;
          const double msl = (g.mb(c) - g.ma(c)) / len, mA = g.ma(c) - msl * g.a;
          dA[i0][c] += A;
          dS[i0][c] += sl;
          dmA[i0][c] += mA;
          dmS[i0][c] += msl;
          dA[i1][c] -= A;
          dS[i1][c] -= sl;
          dmA[i1][c] -= mA;
          dmS[i1][c] -= msl;
        }
      }
      std::array<double, 3> A{0, 0, 0}, S{0, 0, 0}, mA{0, 0, 0}, mS{0, 0, 0};
      for (std::size_t i = 0; i < m; ++i)
        for (int c = 0; c < 3; ++c) {
          A[c] += dA[i][c];
          S[c] += dS[i][c];
          mA[c] += dmA[i][c];
          mS[c] += dmS[i][c];
          qa_[i][c] = A[c] + S[c] * bk[i];
          qs_[i][c] = S[c];
          ma_[i][c] = mA[c] + mS[c] * bk[i];
          ms_[i][c] = mS[c];
        }
    }
    F0.assign(m, {0, 0, 0});
    M0.assign(m, {0, 0, 0});
    C0.assign(m, {0, 0, 0});
    for (std::size_t i = 0; i + 1 < m; ++i) {
      const double h = bk[i + 1] - bk[i];
      for (int c = 0; c < 3; ++c) {
        const double q0 = qa_[i][c], sl = qs_[i][c];
        const double fi = q0 * h + sl * h * h / 2;              // int q
        const double mi = q0 * h * h / 2 + sl * h * h * h / 6;  // int q (b - s)
        F0[i + 1][c] = F0[i][c] + fi;
        M0[i + 1][c] = M0[i][c] + F0[i][c] * h + mi;
        C0[i + 1][c] = C0[i][c] + ma_[i][c] * h + ms_[i][c] * h * h / 2;
      }
    }
  }

  // distributed integrals up to x: force, moment about x, couple
  void dist(int c, double x, double* f, double* mo, double* cp) const {
    const std::size_t m = bk.size();
    if (m < 2 || x <= bk[0]) {
      if (f) *f = 0;
      if (mo) *mo = 0;
      if (cp) *cp = 0;
      return;
    }
    std::size_t i = static_cast<std::size_t>(std::upper_bound(bk.begin(), bk.end(), x) - bk.begin()) - 1;
    if (i >= m - 1) {  // past the last breakpoint: no load beyond
      i = m - 1;
      const double d = x - bk[i];  // beyond the loads: the moment arm only
      if (f) *f = F0[i][c];
      if (mo) *mo = M0[i][c] + F0[i][c] * d;
      if (cp) *cp = C0[i][c];
      return;
    }
    const double d = x - bk[i], q0 = qa_[i][c], sl = qs_[i][c];
    if (f) *f = F0[i][c] + q0 * d + sl * d * d / 2;
    if (mo) *mo = M0[i][c] + F0[i][c] * d + q0 * d * d / 2 + sl * d * d * d / 6;
    if (cp) *cp = C0[i][c] + ma_[i][c] * d + ms_[i][c] * d * d / 2;
  }
  // The cuts of (x0, x1): the abscissae of the point loads `keep` selects
  // and (with_segs) the ends of the distributed loads strictly inside,
  // ascending, between x0 and x1 (both included). Indexed: binary searches.
  template <class Keep>
  void cuts(double x0, double x1, Keep keep, bool with_segs, std::vector<double>& out) const {
    out.clear();
    out.push_back(x0);
    if (!indexed) {
      for (const PointAt& p : sp.points)
        if (keep(p) && p.a > x0 && p.a < x1) out.push_back(p.a);
      if (with_segs)
        for (const DistSeg& sg : sp.segs)
          for (double e : {sg.a, sg.b})
            if (e > x0 && e < x1) out.push_back(e);
      std::sort(out.begin() + 1, out.end());
    } else {
      auto p0 = std::upper_bound(pa.begin(), pa.end(), x0), p1 = std::lower_bound(pa.begin(), pa.end(), x1);
      auto b0 = std::upper_bound(bk.begin(), bk.end(), x0), b1 = std::lower_bound(bk.begin(), bk.end(), x1);
      if (!with_segs) b0 = b1;
      while (p0 < p1 || b0 < b1) {  // merge, ascending
        const bool take_p = b0 >= b1 || (p0 < p1 && *p0 <= *b0);
        if (take_p) {
          if (keep(sp.points[pord[static_cast<std::size_t>(p0 - pa.begin())]])) out.push_back(*p0);
          ++p0;
        } else {
          out.push_back(*b0);
          ++b0;
        }
      }
    }
    out.push_back(x1);
  }

  // number of sorted points counted at the cut x
  std::size_t npts(double x, bool left_side) const {
    const double tol = 1e-11 * std::max(1.0, std::abs(x));
    if (left_side)  // a < x - tol
      return static_cast<std::size_t>(std::lower_bound(pa.begin(), pa.end(), x - tol) - pa.begin());
    return static_cast<std::size_t>(std::upper_bound(pa.begin(), pa.end(), x + tol) - pa.begin());
  }

  // a point load at a counts at the cut x (left_side: only strictly before
  // it). The same abscissa reached by two roundings (the bar's chord length,
  // a load put on its end node) is ONE point: compared to a relative 1e-11.
  static bool at_or_before(double a, double x, bool left_side) {
    const double tol = 1e-11 * std::max(1.0, std::abs(x));
    return left_side ? a < x - tol : a <= x + tol;
  }

  double force(int c, double x, bool left_side = false) const {
    if (indexed) {
      double f;
      dist(c, x, &f, nullptr, nullptr);
      return f + pF[npts(x, left_side)][c];
    }
    // integral of q_c over [0, x] + point sums (left_side excludes a point
    // load exactly at x — the value just before the jump)
    double q = 0;
    for (const DistSeg& s : sp.segs) {
      double lo = s.a, hi = std::min(x, s.b);
      if (hi <= lo) continue;
      double len = s.b - s.a;
      auto qc = [&](double t) { return s.qa(c) + (s.qb(c) - s.qa(c)) * (len > 0 ? t / len : 0.0); };
      double t0 = lo - s.a, t1 = hi - s.a;
      q += 0.5 * (qc(t0) + qc(t1)) * (t1 - t0);  // exact for linear q
      if (s.qc(c) != 0.0 && len > 0) {  // parabolic part 4 qc u (l - u) / l^2
        auto P = [&](double u) { return u * u * (len / 2 - u / 3); };
        q += 4.0 * s.qc(c) / (len * len) * (P(t1) - P(t0));
      }
    }
    for (const PointAt& p : sp.points)
      if (at_or_before(p.a, x, left_side)) q += p.F(c);
    return q;
  }

  double moment(int c, double x) const {  // integral of q_c(s)*(x-s) ds + point sums
    if (indexed) {
      double mo;
      dist(c, x, nullptr, &mo, nullptr);
      const std::size_t k = static_cast<std::size_t>(std::upper_bound(pa.begin(), pa.end(), x) - pa.begin());
      return mo + pF[k][c] * x - pFa[k][c];
    }
    double m = 0;
    for (const DistSeg& s : sp.segs) {
      double lo = s.a, hi = std::min(x, s.b);
      if (hi <= lo) continue;
      double len = s.b - s.a;
      double slope = len > 0 ? (s.qb(c) - s.qa(c)) / len : 0.0;
      // q(t) = qa + slope*(t - a); integral over [lo, hi] of q(t)(x - t) dt, closed form.
      auto F = [&](double t) {
        double u = t - s.a;
        // antiderivative of (qa + slope*u)(x - t) dt with t = a + u:
        // = qa[(x-a)u - u^2/2] + slope[(x-a)u^2/2 - u^3/3]
        return s.qa(c) * ((x - s.a) * u - u * u / 2) + slope * ((x - s.a) * u * u / 2 - u * u * u / 3);
      };
      m += F(hi) - F(lo);
      if (s.qc(c) != 0.0 && len > 0) {
        // parabolic part: antiderivative of 4 qc u (l - u) (X - u) / l^2, X = x - a
        const double X = x - s.a;
        auto P = [&](double t) {
          const double u = t - s.a;
          return len * X * u * u / 2 - (len + X) * u * u * u / 3 + u * u * u * u / 4;
        };
        m += 4.0 * s.qc(c) / (len * len) * (P(hi) - P(lo));
      }
    }
    for (const PointAt& p : sp.points)
      if (p.a <= x) m += p.F(c) * (x - p.a);
    return m;
  }

  // applied point moments M_c with a <= x, plus the distributed moments
  // integrated up to x (a couple field: no shear)
  double point_moment(int c, double x, bool left_side = false) const {
    if (indexed) {
      double cp;
      dist(c, x, nullptr, nullptr, &cp);
      return cp + pM[npts(x, left_side)][c];
    }
    double m = 0;
    for (const PointAt& p : sp.points)
      if (at_or_before(p.a, x, left_side)) m += p.M(c);
    for (const DistSeg& s : sp.segs) {
      if (s.ma(c) == 0 && s.mb(c) == 0) continue;
      const double hi = std::min(x, s.b), len = s.b - s.a;
      if (hi <= s.a || len <= 0) continue;
      const double u = hi - s.a;
      m += s.ma(c) * u + (s.mb(c) - s.ma(c)) * u * u / (2 * len);
    }
    return m;
  }
};

// Adaptive 3-point Gauss of a vector integrand on [a, b]: a piece is kept
// when its two halves agree with it to round-off. A CARA VAR bar's
// integrands (M / EI(x), T / GJ(x), N / EA(x)) are smooth between the load
// cuts but not polynomial: no fixed sub-step count is exact.
template <class F>
Vec3 adaptive_gauss(const F& f, double a, double b) {
  auto g3 = [&](double u, double v) {
    static const double gx = 0.7745966692414834;
    const double m = 0.5 * (u + v), h = 0.5 * (v - u);
    return Vec3(h * (5.0 / 9.0 * f(m - h * gx) + 8.0 / 9.0 * f(m) + 5.0 / 9.0 * f(m + h * gx)));
  };
  if (!(b > a)) return Vec3::Zero();
  const Vec3 whole = g3(a, b);
  const double tol = 1e-14 * std::max(whole.cwiseAbs().maxCoeff(), 1e-300);
  struct Rec {
    const decltype(g3)& g;
    double tol;
    Vec3 run(double u, double v, const Vec3& w, int depth) const {
      const double m = 0.5 * (u + v);
      const Vec3 l = g(u, m), r = g(m, v);
      if (depth >= 24 || (l + r - w).cwiseAbs().maxCoeff() <= tol) return l + r;
      return run(u, m, l, depth + 1) + run(m, v, r, depth + 1);
    }
  };
  return Rec{g3, tol}.run(a, b, whole, 0);
}

void recover_stations(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                      const Vec12& f_end, const Vec12& u_local, int n_st, BarResult& out,
                      const std::vector<double>* custom_x,
                      const std::function<SectionProps(double)>* props_at = nullptr,
                      bool forces_only = false) {
  double L = el.L();
  const Bar& b = el.bar();
  const Material& mat = model.material(b.material);
  const SectionProps& sec = model.section(b.section).props;
  SpanIntegrals in{sp};

  // Study sections (ST1 ETUDE, see etude.h): when a custom abscissa list is
  // given for this bar the stations sit exactly there; otherwise the uniform
  // `stations_per_bar` grid is used, on the bit-identical historical path.
  int n = std::max(1, n_st);
  const bool uniform = custom_x == nullptr || custom_x->empty();
  if (!uniform) {
    n = static_cast<int>(custom_x->size());
    double prev = -1;
    for (double x : *custom_x) {
      require(finite(x) && x >= -1e-9 * L && x <= L * (1 + 1e-9),
              "study sections: abscissa outside bar " + std::to_string(b.id));
      require(x >= prev, "study sections: abscissas must be sorted ascending");
      prev = x;
    }
  }
  out.stations.resize(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k) {
    double x = uniform ? (n == 1 ? L / 2 : L * k / (n - 1))
                       : std::min(std::max((*custom_x)[static_cast<std::size_t>(k)], 0.0), L);
    Station& st = out.stations[static_cast<std::size_t>(k)];
    st.x = x;
    // Both END stations read the member's own end force: a point load sitting
    // exactly on the start is counted (a <= x), one exactly on the end is not
    // (it acts on the joint), so neither end reads the value outside the bar.
    const bool end = x >= L;
    st.N = -(f_end(0) + in.force(0, x, end));
    st.Vy = -(f_end(1) + in.force(1, x, end));
    st.Vz = -(f_end(2) + in.force(2, x, end));
    st.T = -(f_end(3) + in.point_moment(0, x, end));
    st.My = -f_end(4) - f_end(2) * x - in.moment(2, x) - in.point_moment(1, x, end);
    st.Mz = -f_end(5) + f_end(1) * x + in.moment(1, x) - in.point_moment(2, x, end);
  }

  // Deflection recovery by double integration of curvature (+ shear angle),
  // anchored on the END DISPLACEMENTS so released ends recover their true
  // deflected shape. Axial/twist single integration anchored at end i.
  if (forces_only) return;
  if (b.kind == BarKind::Frame && n >= 2) {
    double EIy = mat.E * sec.Iy, EIz = mat.E * sec.Iz;
    double EA = mat.E * sec.A;
    double GAy = sec.Avy > 0 ? mat.G() * sec.Avy : 0;
    double GAz = sec.Avz > 0 ? mat.G() * sec.Avz : 0;
    double h_uniform = L / (n - 1);
    // Variable section (CARA VAR): stiffnesses at x, and each station
    // interval integrated by a composite Simpson rule of `nsub` sub-steps
    // (the curvature M / EI(x) is no longer quadratic between stations).
    auto EI_at = [&](int plane, double x) {
      if (!props_at) return plane == 1 ? EIz : EIy;
      SectionProps p = (*props_at)(std::min(std::max(x / L, 0.0), 1.0));
      return mat.E * (plane == 1 ? p.Iz : p.Iy);
    };
    auto GA_at = [&](int plane, double x, double GA) {
      if (!props_at || GA <= 0) return GA;
      SectionProps p = (*props_at)(std::min(std::max(x / L, 0.0), 1.0));
      double Av = plane == 1 ? p.Avy : p.Avz;
      return Av > 0 ? mat.G() * Av : GA;
    };
    auto EA_at = [&](double x) {
      if (!props_at) return EA;
      return mat.E * (*props_at)(std::min(std::max(x / L, 0.0), 1.0)).A;
    };
    const int nsub = 1;  // a variable section integrates adaptively (adaptive_gauss)
    // Per-segment Simpson double integration: the bending moment is CONTINUOUS
    // (only shear jumps), so evaluating the curvature at segment midpoints is
    // side-unambiguous and the scheme is exact for our piecewise-linear load
    // vocabulary whenever point loads land on stations.
    // (`left`: the value just before an applied couple at x — the moment
    // itself jumps there, so the end of an interval takes its inner side)
    auto moment_at = [&](int plane, double x, bool left = false) {  // plane 1 = Mz, 2 = My
      if (plane == 1)
        return -f_end(5) + f_end(1) * x + in.moment(1, x) - in.point_moment(2, x, left);
      return -f_end(4) - f_end(2) * x - in.moment(2, x) - in.point_moment(1, x, left);
    };
    auto shear_at = [&](int c, double x, bool left) { return -(f_end(c) + in.force(c, x, left)); };
    // The imposed-curvature field is linear (kap0 + dkap*x, ST1 DEFOR LIN) and
    // a concentrated dislocation (DEFOR CON) is its Dirac: a slope jump
    // -sgn*theta at the kink abscissa, invisible AT its own station
    // (half-open, like a point load), superposed on the smooth Simpson
    // integration below.
    //
    // SIGN of the imposed-curvature term: the theta_y = -dw/dx flip applies to
    // the curvature exactly as it does to the moment, so kap sits INSIDE the
    // sgn factor: v'' = +(Mz/EIz - kz), w'' = -(My/EIy - ky). With -kap in
    // BOTH planes a single bar's INTERIOR deflections under an imposed
    // kappa_y would bow OPPOSITE to the nodal solution of the same beam
    // meshed in two — the f0 lowering, i.e. the nodal convention, is the
    // authority.
    auto integrate2 = [&](int plane, int c, double kap0, double dkap, double c2kap, double GA, double d0,
                          double dL, double r0, auto kink_of, auto jump_of, auto set,
                          auto set_rot) {
      std::vector<double> dp(static_cast<std::size_t>(n), 0.0), v(static_cast<std::size_t>(n), 0.0);
      double sgn = plane == 1 ? 1.0 : -1.0;
      double EI = plane == 1 ? EIz : EIy;
      for (int k = 1; k < n; ++k) {
        double x0 = out.stations[static_cast<std::size_t>(k - 1)].x;
        double x1 = out.stations[static_cast<std::size_t>(k)].x;
        double h = uniform ? h_uniform : (x1 - x0);  // per-segment on custom sections
        double vk = v[static_cast<std::size_t>(k - 1)], dk = dp[static_cast<std::size_t>(k - 1)];
        // Exact integration: the interval is split at the transverse point
        // loads and couples inside it and at the ends of the distributed
        // loads, so the curvature is a polynomial (cubic: M under a linear
        // load) on each piece; v(b) = v(a) + (b - a) v'(a) + int (b - s) c(s)
        // and v'(b) = v'(a) + int c(s) by 3-point Gauss, exact to degree 5
        // (Simpson was exact only when point loads landed on stations). A
        // variable section keeps its 16 sub-pieces per interval.
        (void)h;
        std::vector<double> cuts;
        in.cuts(x0, x1, [&](const PointAt& p) { return p.F(c) != 0.0 || p.M(plane == 1 ? 2 : 1) != 0.0; },
                true, cuts);
        static const double gx[3] = {-0.7745966692414834, 0.0, 0.7745966692414834};
        static const double gw[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
        auto curv = [&](double x) {
          const double EIx = props_at ? EI_at(plane, x) : EI;
          return sgn * (moment_at(plane, x) / EIx - (kap0 + (dkap + c2kap * x) * x));
        };
        for (std::size_t q = 1; q < cuts.size(); ++q) {
          const double pa = cuts[q - 1], pb = cuts[q];
          if (pb - pa <= 0) continue;
          if (props_at) {
            // v(b) = v(a) + (b - a) v'(a) + int (b - s) c(s) + int V / GA
            const Vec3 I = adaptive_gauss(
                [&](double x) {
                  const double cg = curv(x);
                  return Vec3(cg, (pb - x) * cg,
                              GA > 0 ? shear_at(c, x, false) / GA_at(plane, x, GA) : 0.0);
                },
                pa, pb);
            vk += (pb - pa) * dk + I(1) + I(2);
            dk += I(0);
            continue;
          }
          for (int j2 = 0; j2 < nsub; ++j2) {
            const double a = pa + (pb - pa) * j2 / nsub, b2 = pa + (pb - pa) * (j2 + 1) / nsub;
            const double hh = b2 - a, mid = 0.5 * (a + b2);
            double iv = 0, id = 0, ig = 0;
            for (int g = 0; g < 3; ++g) {
              const double x = mid + 0.5 * hh * gx[g], wg = 0.5 * hh * gw[g];
              const double cg = curv(x);
              id += wg * cg;
              iv += wg * (b2 - x) * cg;
              if (GA > 0) ig += wg * shear_at(c, x, false) / GA;
            }
            vk += hh * dk + iv + ig;
            dk += id;
          }
        }
        v[static_cast<std::size_t>(k)] = vk;
        dp[static_cast<std::size_t>(k)] = dk;
        for (const KinkAt& kk : sp.kinks) {  // dislocations inside (x0, x1]
          double th = kink_of(kk), jv = jump_of(kk);
          if ((th == 0.0 && jv == 0.0) || kk.a <= x0 || kk.a > x1) continue;
          v[static_cast<std::size_t>(k)] += -sgn * th * (x1 - kk.a) + jv;
          dp[static_cast<std::size_t>(k)] += -sgn * th;
        }
      }
      // Linear correction to anchor on the end displacements.
      double a0 = d0 - v[0];
      double a1 = (dL - (v[static_cast<std::size_t>(n - 1)] + a0)) / L;
      for (int k = 0; k < n; ++k) {
        double x = out.stations[static_cast<std::size_t>(k)].x;
        set(k, v[static_cast<std::size_t>(k)] + a0 + a1 * x);
        // Section rotation = end rotation + integrated curvature (the shear
        // angle is not a rotation): theta_z = +v', theta_y = -w' (bending).
        set_rot(k, r0 + sgn * dp[static_cast<std::size_t>(k)]);
      }
    };
    integrate2(1, 1, sp.kappa_z, sp.dkappa_z, sp.c2kappa_z, GAy, u_local(1), u_local(7), u_local(5),
               [](const KinkAt& kk) { return kk.theta_z; }, [](const KinkAt& kk) { return kk.dvy; },
               [&](int k, double val) { out.stations[static_cast<std::size_t>(k)].defl.y() = val; },
               [&](int k, double val) { out.stations[static_cast<std::size_t>(k)].rot.z() = val; });
    integrate2(2, 2, sp.kappa_y, sp.dkappa_y, sp.c2kappa_y, GAz, u_local(2), u_local(8), u_local(4),
               [](const KinkAt& kk) { return kk.theta_y; }, [](const KinkAt& kk) { return kk.dvz; },
               [&](int k, double val) { out.stations[static_cast<std::size_t>(k)].defl.z() = val; },
               [&](int k, double val) { out.stations[static_cast<std::size_t>(k)].rot.y() = val; });
    // Twist: theta_x = theta_i + int T / GJ, T = -(f_end(3) + applied
    // torques up to x) — quadratic under a distributed torque (Simpson per
    // piece, split at the point torques), anchored linearly on theta_j.
    const double GJ = mat.G() * sec.J;
    std::vector<double> tw(static_cast<std::size_t>(n), 0.0);
    auto GJ_at = [&](double x) {
      if (!props_at) return GJ;
      return mat.G() * (*props_at)(std::min(std::max(x / L, 0.0), 1.0)).J;
    };
    if (GJ > 0 && props_at) {  // CARA VAR: int T / GJ(x), adaptive between the point torques
      for (int k = 1; k < n; ++k) {
        const double xa = out.stations[static_cast<std::size_t>(k - 1)].x;
        const double xb = out.stations[static_cast<std::size_t>(k)].x;
        std::vector<double> cuts;
        in.cuts(xa, xb, [](const PointAt& p) { return p.M(0) != 0.0; }, true, cuts);
        double t = tw[static_cast<std::size_t>(k - 1)];
        for (std::size_t q = 1; q < cuts.size(); ++q)
          t += adaptive_gauss(
                   [&](double x) {
                     const double g = GJ_at(x);
                     return Vec3((g > 0 ? -(f_end(3) + in.point_moment(0, x, false)) / g : 0.0) + sp.kappa_x +
                                     (sp.dkappa_x + sp.c2kappa_x * x) * x,
                                 0, 0);
                   },
                   cuts[q - 1], cuts[q])(0);
        for (const KinkAt& kk : sp.kinks)
          if (kk.dthx != 0.0 && kk.a > xa && kk.a <= xb) t += kk.dthx;
        tw[static_cast<std::size_t>(k)] = t;
      }
    }
    if (GJ > 0 && !props_at) {
      auto curv_x = [&](double x, bool left) {  // theta_x' = T / GJ + imposed twist rate
        return -(f_end(3) + in.point_moment(0, x, left)) / GJ + sp.kappa_x + (sp.dkappa_x + sp.c2kappa_x * x) * x;
      };
      for (int k = 1; k < n; ++k) {
        const double xa = out.stations[static_cast<std::size_t>(k - 1)].x;
        const double xb = out.stations[static_cast<std::size_t>(k)].x;
        std::vector<double> cuts;
        in.cuts(xa, xb, [](const PointAt& p) { return p.M(0) != 0.0; }, false, cuts);
        double t = tw[static_cast<std::size_t>(k - 1)];
        for (std::size_t q = 1; q < cuts.size(); ++q) {
          const double a = cuts[q - 1], b2 = cuts[q];
          if (b2 - a > 0)
            t += (b2 - a) / 6.0 *
                 (curv_x(a, false) + 4.0 * curv_x(0.5 * (a + b2), false) + curv_x(b2, true));
        }
        for (const KinkAt& kk : sp.kinks)  // twist jumps inside (xa, xb]
          if (kk.dthx != 0.0 && kk.a > xa && kk.a <= xb) t += kk.dthx;
        tw[static_cast<std::size_t>(k)] = t;
      }
    }
    const double tw_err = u_local(9) - u_local(3) - tw[static_cast<std::size_t>(n - 1)];
    for (int k = 0; k < n; ++k) {
      double t = out.stations[static_cast<std::size_t>(k)].x / L;
      out.stations[static_cast<std::size_t>(k)].rot.x() =
          u_local(3) + tw[static_cast<std::size_t>(k)] + t * tw_err;
    }
    // Axial: u(x) = u_i + int (N/EA + eps0(x)) + axial gap jumps past their
    // abscissa (trapezoid; exact for the piecewise-linear N and linear eps0).
    // (accumulated interval by interval: each interval integrated once)
    double u_run = u_local(0);
    for (int k = 0; k < n; ++k) {
      double u = u_run;
      for (int j = k; j <= k && j >= 1; ++j) {
        double xa = out.stations[static_cast<std::size_t>(j - 1)].x;
        double xb = out.stations[static_cast<std::size_t>(j)].x;
        if (props_at) {  // N / EA(x): adaptive between the axial point loads
          std::vector<double> cuts;
          in.cuts(xa, xb, [](const PointAt& p) { return p.F(0) != 0.0; }, true, cuts);
          for (std::size_t q = 1; q < cuts.size(); ++q)
            u += adaptive_gauss(
                     [&](double x) {
                       return Vec3(-(f_end(0) + in.force(0, x, false)) / EA_at(x) + sp.eps0 +
                                       (sp.deps0 + sp.c2eps0 * x) * x,
                                   0, 0);
                     },
                     cuts[q - 1], cuts[q])(0);
          continue;
        }
        // N / EA piecewise: split at the axial point loads inside the
        // interval (and at the ends of the distributed loads) and take N
        // from the inner side of each piece (a point
        // load AT a station must not leak into the interval before it);
        // Simpson per piece is exact for the quadratic N of a linear load.
        std::vector<double> cuts;
        in.cuts(xa, xb, [](const PointAt& p) { return p.F(0) != 0.0; }, true, cuts);
        auto strain = [&](double x, bool left) {
          return -(f_end(0) + in.force(0, x, left)) / EA + sp.eps0 + (sp.deps0 + sp.c2eps0 * x) * x;
        };
        for (std::size_t q = 1; q < cuts.size(); ++q) {
          const double a = cuts[q - 1], b2 = cuts[q];
          if (b2 - a <= 0) continue;
          u += (b2 - a) / 6.0 *
               (strain(a, false) + 4.0 * strain(0.5 * (a + b2), false) + strain(b2, true));
        }
      }
      u_run = u;
      double xk = out.stations[static_cast<std::size_t>(k)].x;
      for (const KinkAt& kk : sp.kinks)
        if (kk.dux != 0.0 && kk.a < xk) u += kk.dux;
      out.stations[static_cast<std::size_t>(k)].defl.x() = u;
    }
  } else if (n >= 2) {
    // Truss: linear interpolation of the end displacements.
    for (int k = 0; k < n; ++k) {
      double t = out.stations[static_cast<std::size_t>(k)].x / L;
      out.stations[static_cast<std::size_t>(k)].defl =
          (1 - t) * u_local.segment<3>(0) + t * u_local.segment<3>(6);
      out.stations[static_cast<std::size_t>(k)].rot =
          (1 - t) * u_local.segment<3>(3) + t * u_local.segment<3>(9);
    }
  }

  // Soil pressure recovery p = k * u_local(x) (linear interpolation of the
  // displacement field between end values — meshing-bounded approximation).
  if (b.soil.any() && n >= 1) {
    out.soil_p.resize(static_cast<std::size_t>(n));
    out.soil_m.resize(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) {
      const Station& st = out.stations[static_cast<std::size_t>(k)];
      double t = L > 0 ? st.x / L : 0;
      out.soil_p[static_cast<std::size_t>(k)] =
          Vec3(b.soil.kfx * st.defl.x(), b.soil.kfy * st.defl.y(), b.soil.kfz * st.defl.z());
      double thx = (1 - t) * u_local(3) + t * u_local(9);
      out.soil_m[static_cast<std::size_t>(k)] = b.soil.kmx * thx;
    }
  }
}

}  // namespace

void recover_bar_stations(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                          const Vec12& f_end, const Vec12& u_local, const std::vector<double>& xs,
                          BarResult& out) {
  const TaperedBar* tb = model.tapered_of(el.bar().id);
  if (tb && tb->recover) {
    tb->recover(sp, u_local, xs, out);
    return;
  }
  recover_stations(model, el, sp, f_end, u_local, static_cast<int>(xs.size()), out, &xs,
                   tb && tb->props ? &tb->props : nullptr);
}

// ---------------- exact constant-N beam-column ----------------

namespace {

// F_m(s) = sum_{j>=0} alpha^j s^(2j+m) / (2j+m)!, m = 0..5: the transfer
// functions of y'' - alpha y = h — F_0 = C, F_1 = S the homogeneous pair
// (y(s) = y0 C + y0' S), F_(m+2) the response from rest to h = s^m / m!
// (F_m' = F_(m-1), F_0' = alpha F_1, F_m - alpha F_(m+2) = s^m / m!).
// |alpha s^2| <= 4: the series, 17 terms (truncation below 1e-25
// relative); beyond: cosh/sinh (tension) or cos/sin (compression) and the
// exact recurrence F_(m+2) = (F_m - s^m/m!) / alpha, whose cancellation is
// bounded by (m+1)(m+2) / |alpha s^2| < 5 per step.
std::array<double, 7> bc_transfer(double alpha, double s) {
  static const double inv_fact[7] = {1.0, 1.0, 0.5, 1.0 / 6, 1.0 / 24, 1.0 / 120, 1.0 / 720};
  std::array<double, 7> F{};
  double sm[7];
  sm[0] = 1.0;
  for (int m = 1; m < 7; ++m) sm[m] = sm[m - 1] * s;
  const double z = alpha * s * s;
  if (std::fabs(z) <= 4.0) {
    for (int m = 0; m < 7; ++m) {
      double term = sm[m] * inv_fact[m], sum = term;
      for (int j = 1; j <= 16; ++j) {
        term *= z / ((2.0 * j + m - 1) * (2.0 * j + m));
        sum += term;
      }
      F[static_cast<std::size_t>(m)] = sum;
    }
    return F;
  }
  const double k = std::sqrt(std::fabs(alpha));
  if (alpha > 0) {
    F[0] = std::cosh(k * s);
    F[1] = std::sinh(k * s) / k;
  } else {
    F[0] = std::cos(k * s);
    F[1] = std::sin(k * s) / k;
  }
  for (int m = 2; m < 7; ++m)
    F[static_cast<std::size_t>(m)] =
        (F[static_cast<std::size_t>(m - 2)] - sm[m - 2] * inv_fact[m - 2]) / alpha;
  return F;
}

// March y'' - alpha y = h from (y, y')(0) = (0, yp0) over the pieces
// [cuts[q-1], cuts[q]] (cuts ascending from 0 to L), h a cubic on each piece,
// recovered exactly from four interior Chebyshev samples (h is sampled
// strictly inside a piece, so a jump AT a cut is never straddled). Stations
// `xs` (ascending, may be null) receive y and y'. Returns (y, y')(L).
template <class H>
std::pair<double, double> bc_march(double alpha, double yp0, const std::vector<double>& cuts, H&& h,
                                   const std::vector<double>* xs, std::vector<double>* ys,
                                   std::vector<double>* yps) {
  static const std::array<double, 4> tau = [] {
    std::array<double, 4> t{};
    for (int i = 0; i < 4; ++i)
      t[static_cast<std::size_t>(i)] = 0.5 * (1.0 - std::cos(kPi * (2 * i + 1) / 8.0));
    return t;
  }();
  static const Eigen::Matrix4d Vinv = [] {
    Eigen::Matrix4d V;
    for (int i = 0; i < 4; ++i)
      for (int m = 0; m < 4; ++m) V(i, m) = std::pow(tau[static_cast<std::size_t>(i)], m);
    return Eigen::Matrix4d(V.inverse());
  }();
  static const double fact[4] = {1.0, 1.0, 2.0, 6.0};
  double y = 0.0, yp = yp0;
  std::size_t k = 0;
  const std::size_t ns = xs ? xs->size() : 0;
  if (xs) {
    ys->assign(ns, 0.0);
    yps->assign(ns, 0.0);
  }
  while (k < ns && (*xs)[k] <= cuts.front()) {
    (*yps)[k] = yp0;
    ++k;
  }
  for (std::size_t q = 1; q < cuts.size(); ++q) {
    const double a = cuts[q - 1], b = cuts[q], l = b - a;
    const bool last = q + 1 == cuts.size();
    if (l <= 0 && !last) continue;
    double c[4] = {0, 0, 0, 0};
    if (l > 0) {
      Eigen::Vector4d hv;
      for (int i = 0; i < 4; ++i) hv(i) = h(a + l * tau[static_cast<std::size_t>(i)]);
      const Eigen::Vector4d d = Vinv * hv;
      double lp = 1.0;
      for (int m = 0; m < 4; ++m) {
        c[m] = fact[m] * d(m) / lp;  // h^(m)(a+)
        lp *= l;
      }
    }
    const double ya = y, ypa = yp;
    auto at = [&](double t, double& yo, double& ypo) {
      const std::array<double, 7> F = bc_transfer(alpha, t);
      yo = ya * F[0] + ypa * F[1] + c[0] * F[2] + c[1] * F[3] + c[2] * F[4] + c[3] * F[5];
      ypo = ya * alpha * F[1] + ypa * F[0] + c[0] * F[1] + c[1] * F[2] + c[2] * F[3] + c[3] * F[4];
    };
    while (k < ns && ((*xs)[k] <= b || last)) {
      at(std::min(std::max((*xs)[k], a), b) - a, (*ys)[k], (*yps)[k]);
      ++k;
    }
    if (l > 0) at(l, y, yp);
  }
  return {y, yp};
}

struct BcBar {
  double L = 0, EIy = 0, EIz = 0;
  std::vector<double> cuts;  // 0 .. L, at every point load and segment end
};

BcBar bc_setup(const Model& model, const BarElement& el, const BarSpanLoads& sp, double N) {
  const Bar& b = el.bar();
  const std::string tag = "second order (exact stability): bar " + std::to_string(b.id);
  require(model.tapered_of(b.id) == nullptr,
          tag + " is a variable section (the stability functions are prismatic)");
  require(!b.soil.any(), tag + " has a Winkler bed (the beam-column ODE carries no foundation term)");
  require(sp.kinks.empty(),
          tag + " carries a concentrated dislocation (not in the exact beam-column member)");
  for (const DistSeg& g : sp.segs)
    require(g.qa(0) == 0.0 && g.qb(0) == 0.0,
            tag + " carries an axial span load: N varies along the member (mesh it, nodal loads)");
  require_linear_span(sp.segs, tag);
  require_no_imposed_twist(sp, tag);
  require_linear_strain(sp, tag);
  for (const PointAt& p : sp.points)
    require(p.F(0) == 0.0,
            tag + " carries an axial point load: N varies along the member (split the bar there)");
  const Material& mat = model.material(b.material);
  const SectionProps& sec = model.section(b.section).props;
  BcBar bc;
  bc.L = el.L();
  bc.EIy = mat.E * sec.Iy;
  bc.EIz = mat.E * sec.Iz;
  const double EImin = std::min(bc.EIy, bc.EIz);
  require(!(N > 0) || bc.L * std::sqrt(N / EImin) <= 18.0,
          tag + ": tension with L sqrt(N/EI) > 18 — the transfer from end i loses e^(kL) of "
                "precision; split the member");
  SpanIntegrals in{sp};
  in.cuts(0.0, bc.L, [](const PointAt&) { return true; }, true, bc.cuts);
  bc.cuts.push_back(bc.L);
  return bc;
}

// Curvature source h(x) of plane 1 (v, Mz) or 2 (w, My) for the end forces
// f: h = sgn (M1(x) / EI - kappa(x)), M1 the first-order moment of f and the
// span loads, sgn = +1 / -1 (theta_y = -w').
struct BcCurv {
  const SpanIntegrals& in;
  const BarSpanLoads& sp;
  double EIy, EIz;
  double operator()(int plane, const Vec12& f, double x) const {
    if (plane == 1)
      return (-f(5) + f(1) * x + in.moment(1, x) - in.point_moment(2, x)) / EIz -
             (sp.kappa_z + sp.dkappa_z * x);
    return -((-f(4) - f(2) * x - in.moment(2, x) - in.point_moment(1, x)) / EIy -
             (sp.kappa_y + sp.dkappa_y * x));
  }
};

}  // namespace

std::array<double, 7> beam_transfer(double alpha, double s) { return bc_transfer(alpha, s); }

Vec12 beam_column_f0(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                     const Vec12& f0_linear, double N) {
  if (el.bar().kind != BarKind::Frame) return f0_linear;
  bool loaded = sp.kappa_y != 0 || sp.kappa_z != 0 || sp.dkappa_y != 0 || sp.dkappa_z != 0 ||
                has_parabolic_strain(sp);
  for (const DistSeg& g : sp.segs)
    loaded |= !g.qc.isZero() || g.qa(1) != 0 || g.qb(1) != 0 || g.qa(2) != 0 || g.qb(2) != 0 || g.ma(1) != 0 ||
              g.mb(1) != 0 || g.ma(2) != 0 || g.mb(2) != 0;
  for (const PointAt& p : sp.points)
    loaded |= p.F(1) != 0 || p.F(2) != 0 || p.M(1) != 0 || p.M(2) != 0;
  Vec12 f0 = f0_linear;
  if (!loaded) {  // no bending span load: no bending fixed-end force
    for (int c : {1, 2, 4, 5, 7, 8, 10, 11}) f0(c) = 0.0;
    return f0;
  }
  BcBar bc = bc_setup(model, el, sp, N);
  SpanIntegrals in{sp};
  BcCurv curv{in, sp, bc.EIy, bc.EIz};
  const double L = bc.L;
  const Vec12 zero = Vec12::Zero();
  for (int plane : {1, 2}) {
    const double EI = plane == 1 ? bc.EIz : bc.EIy;
    const double sgn = plane == 1 ? 1.0 : -1.0;  // also the sign of the shear lever in M1
    const double alpha = N / EI;
    // load response from rest at a clamped end i
    const auto yl = bc_march(alpha, 0.0, bc.cuts, [&](double x) { return curv(plane, zero, x); },
                             nullptr, nullptr, nullptr);
    // end forces at i (moment f_r, shear f_v) closing y(L) = y'(L) = 0:
    // h += sgn (-f_r + sgn f_v x) / EI
    const std::array<double, 7> F = bc_transfer(alpha, L);
    Eigen::Matrix2d A;
    A << -F[2], sgn * F[3], -F[1], sgn * F[2];
    A *= sgn / EI;
    const Eigen::Vector2d fe = A.fullPivLu().solve(Eigen::Vector2d(-yl.first, -yl.second));
    const int iv = plane, ir = plane == 1 ? 5 : 4, jv = iv + 6, jr = ir + 6;
    // f0 = -(end forces of the clamped member); end j by equilibrium of the
    // increments over the linear f0 (same loads, the clamped chord is still)
    const double d_v = -fe(1) - f0_linear(iv), d_r = -fe(0) - f0_linear(ir);
    f0(iv) = f0_linear(iv) + d_v;
    f0(ir) = f0_linear(ir) + d_r;
    f0(jv) = f0_linear(jv) - d_v;
    f0(jr) = f0_linear(jr) - d_r + sgn * L * d_v;
  }
  return f0;
}

void recover_beam_column(const Model& model, const BarElement& el, const BarSpanLoads& sp,
                         const Vec12& f0_exact, double N, const Vec12& u_nodal,
                         const StaticOptions& opt, BarResult& br) {
  const Mat12 kx = el.k_exact_local(N);
  const Vec12 u_local = el.local_displacements_with(kx, u_nodal, f0_exact);
  br.end_forces_local = kx * u_local - f0_exact;
  const std::vector<double>* custom_x = nullptr;
  if (opt.stations_abs) {
    auto it = opt.stations_abs->find(el.bar().id);
    if (it != opt.stations_abs->end() && !it->second.empty()) custom_x = &it->second;
  }
  if (opt.stations_per_bar <= 0 && !custom_x) return;
  // first-order station fields of the exact end forces (N, V, T, axial, twist)
  recover_stations(model, el, sp, br.end_forces_local, u_local, opt.stations_per_bar, br, custom_x,
                   nullptr, opt.forces_only);
  if (el.bar().kind != BarKind::Frame) return;
  BcBar bc = bc_setup(model, el, sp, N);
  SpanIntegrals in{sp};
  BcCurv curv{in, sp, bc.EIy, bc.EIz};
  std::vector<double> xs;
  xs.reserve(br.stations.size());
  for (const Station& st : br.stations) xs.push_back(st.x);
  std::vector<double> ys, yps;
  for (int plane : {1, 2}) {
    const double EI = plane == 1 ? bc.EIz : bc.EIy;
    const double yp0 = plane == 1 ? u_local(5) : -u_local(4);
    bc_march(N / EI, yp0, bc.cuts, [&](double x) { return curv(plane, br.end_forces_local, x); },
             &xs, &ys, &yps);
    for (std::size_t k = 0; k < br.stations.size(); ++k) {
      Station& st = br.stations[k];
      if (plane == 1) {
        st.Mz += N * ys[k];  // M = M1 + N v (v counted from end i)
        if (!opt.forces_only) {
          st.defl.y() = u_local(1) + ys[k];
          st.rot.z() = yps[k];
        }
      } else {
        st.My -= N * ys[k];  // M = M1 - N w (theta_y = -w')
        if (!opt.forces_only) {
          st.defl.z() = u_local(2) + ys[k];
          st.rot.y() = -yps[k];
        }
      }
    }
  }
}

CaseResult recover_case(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                        const VecX& u_free, const StaticOptions& opt, Id id,
                        const std::string& name) {
  const DofMap& d = asm_.dofs;
  int nn = d.n_nodes;
  CaseResult r;
  r.id = id;
  r.name = name;

  // Full displacement vector: free + prescribed + zeros.
  VecX u_full = VecX::Zero(kDofPerNode * nn);
  const bool has_u0 = cl.link_u0.size() > 0;
  for (int f = 0; f < kDofPerNode * nn; ++f) {
    int fi = d.free_ix[static_cast<std::size_t>(f)];
    if (fi >= 0) u_full(f) = u_free(fi) + (has_u0 ? cl.link_u0(fi) : 0.0);
  }
  for (const auto& [pf, pv] : cl.prescribed) u_full(pf) = pv;

  r.node_disp.resize(static_cast<std::size_t>(nn));
  for (int ni = 0; ni < nn; ++ni)
    for (int c = 0; c < 6; ++c) r.node_disp[static_cast<std::size_t>(ni)](c) = u_full(d.full(ni, c));

  // Internal nodal force accumulation (for reactions).
  VecX fint = VecX::Zero(kDofPerNode * nn);

  r.bars.resize(asm_.elements.size());
  // per bar: K_e u_e for the reactions (summed below, in bar order)
  std::vector<Vec12> f_nodal_of(asm_.elements.size(), Vec12::Zero());
  auto recover_bar = [&](std::size_t e) {
    const BarElement& el = asm_.elements[e];
    const Bar& b = model.bars[e];
    BarResult& br = r.bars[e];
    br.bar = b.id;
    if (b.kind == BarKind::Cable) return;

    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    Vec12 u_nodal;
    for (int c = 0; c < 6; ++c) {
      u_nodal(c) = u_full(d.full(ni, c));
      u_nodal(6 + c) = u_full(d.full(nj, c));
    }
    // Reactions use R = K u - F: accumulate K u (below); the applied vector
    // (nodal + equivalent span loads) is subtracted once via cl.nodal_full.
    f_nodal_of[e] = bar_ke_assembled(model, el) * u_nodal;
    // Member forces: pure member matrices, releases re-injected, soil excluded.
    // ST1 CARA VAR: the tapered member block replaces the prismatic one (a
    // tapered bar carries no release, so u_local is simply T u_nodal).
    Vec12 u_local = el.local_displacements(u_nodal, cl.f0[e]);
    const TaperedBar* tb = model.tapered_of(b.id);
    br.end_forces_local = tb ? Vec12(tb->ke * u_local - cl.f0[e])
                             : el.end_forces_local(u_local, cl.f0[e]);
    const std::vector<double>* custom_x = nullptr;
    if (opt.stations_abs) {
      auto it = opt.stations_abs->find(b.id);
      if (it != opt.stations_abs->end() && !it->second.empty()) custom_x = &it->second;
    }
    if (tb && tb->recover && (opt.stations_per_bar > 0 || custom_x)) {
      // Custom member block (bef.h): its own exact station recovery.
      std::vector<double> xs;
      if (custom_x) {
        xs = *custom_x;
      } else {
        int n = std::max(1, opt.stations_per_bar);
        double L = el.L();
        for (int k = 0; k < n; ++k) xs.push_back(n == 1 ? L / 2 : L * k / (n - 1));
      }
      tb->recover(cl.span[e], u_local, xs, br);
    } else if (opt.stations_per_bar > 0 || custom_x) {
      recover_stations(model, el, cl.span[e], br.end_forces_local, u_local, opt.stations_per_bar,
                       br, custom_x, tb && tb->props ? &tb->props : nullptr, opt.forces_only);
    }
  };
  // (a few bars already pay: a staged deck's bars carry many stations and
  // hundreds of creep and tendon loads each, ~70 us of recovery a bar)
  if (opt.bar_threads == 1 || asm_.elements.size() < 4) {
    for (std::size_t e = 0; e < asm_.elements.size(); ++e) recover_bar(e);
  } else {
    if (!model.bars.empty()) (void)model.node_index(model.bars.front().ni);  // id index built before the workers read it
    parallel_for(asm_.elements.size(), opt.bar_threads, recover_bar);
  }
  for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) continue;
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    const Vec12& f_nodal = f_nodal_of[e];
    for (int c = 0; c < 6; ++c) {
      fint(d.full(ni, c)) += f_nodal(c);
      fint(d.full(nj, c)) += f_nodal(6 + c);
    }
  }

  // Facet shells: internal force accumulation (for reactions) + centroid
  // resultant recovery, SUBTRACTING the thermal initial strains so a free
  // thermal expansion reads exactly zero.
  r.shells.resize(asm_.shells.size());
  for (std::size_t s = 0; s < asm_.shells.size(); ++s) {
    const ShellElem& se = asm_.shells[s];
    Vec18 u18;
    for (int k = 0; k < 3; ++k)
      for (int c = 0; c < 6; ++c)
        u18(6 * k + c) = u_full(d.full(se.nodes[static_cast<std::size_t>(k)], c));
    Vec18 f18 = se.el.ke() * u18;
    for (int k = 0; k < 3; ++k)
      for (int c = 0; c < 6; ++c)
        fint(d.full(se.nodes[static_cast<std::size_t>(k)], c)) += f18(6 * k + c);
    ShellRecovery rec;
    if (!cl.shell_recovery.empty()) rec = cl.shell_recovery[se.index];
    ShellResult& sr = r.shells[s];
    sr.shell = se.id;
    sr.slab = se.slab;
    sr.r = se.el.recover(u18 - rec.dstar, rec.eps0, rec.kappa0);
  }

  // Continuum elements: internal forces for the reactions (stresses by
  // continuum_stress)
  for (const ContinuumKernel& k : asm_.continua) {
    const int n = static_cast<int>(k.nodes.size());
    VecX ue(k.dim * n);
    for (int i = 0; i < n; ++i)
      for (int a = 0; a < k.dim; ++a) ue(k.dim * i + a) = u_full(d.full(k.nodes[static_cast<std::size_t>(i)], a));
    const VecX fe = k.K * ue;
    for (int i = 0; i < n; ++i)
      for (int a = 0; a < k.dim; ++a) fint(d.full(k.nodes[static_cast<std::size_t>(i)], a)) += fe(k.dim * i + a);
  }

  // Elastic links: their nodal forces (reactions) and link-frame results.
  r.links.reserve(asm_.links.size());
  for (const LinkElem& e : asm_.links) {
    Vec12 u12;
    for (int c = 0; c < 6; ++c) {
      u12(c) = u_full(d.full(e.ni, c));
      u12(6 + c) = u_full(d.full(e.nj, c));
    }
    const Vec12 f12 = e.K * u12;
    for (int c = 0; c < 6; ++c) {
      fint(d.full(e.ni, c)) += f12(c);
      fint(d.full(e.nj, c)) += f12(6 + c);
    }
    LinkResult lr;
    lr.link = e.id;
    lr.deformation = e.B * u12;
    lr.force = e.k * lr.deformation;
    r.links.push_back(lr);
  }
  // General constraints and skew supports: their multipliers from the
  // unbalanced forces K u - F on the free DOFs they touch (least squares on
  // the constraint gradients: exact when those DOFs carry no rigid link).
  if (!model.constraints.empty() || !model.skew_supports.empty()) {
    std::vector<std::vector<std::pair<int, double>>> grads;  // per multiplier: (full dof, coef)
    for (const Constraint& c : model.constraints) {
      std::vector<std::pair<int, double>> g;
      for (const ConstraintTerm& t : c.terms) g.emplace_back(d.full(model.node_index(t.node), t.comp), t.coef);
      grads.push_back(g);
    }
    std::vector<std::pair<std::size_t, int>> skew_cols;  // (skew index, local comp)
    for (std::size_t k = 0; k < model.skew_supports.size(); ++k) {
      const SkewSupport& sk = model.skew_supports[k];
      const int ni = model.node_index(sk.node);
      for (int c = 0; c < 6; ++c) {
        if (!sk.fixed[static_cast<std::size_t>(c)]) continue;
        std::vector<std::pair<int, double>> g;
        for (int a = 0; a < 3; ++a) g.emplace_back(d.full(ni, (c < 3 ? 0 : 3) + a), sk.axes(c % 3, a));
        grads.push_back(g);
        skew_cols.emplace_back(k, c);
      }
    }
    std::map<int, int> row;
    for (const auto& g : grads)
      for (const auto& [f, c] : g)
        if (d.free_ix[static_cast<std::size_t>(f)] >= 0 && !row.count(f)) row.emplace(f, static_cast<int>(row.size()));
    MatX A = MatX::Zero(static_cast<Eigen::Index>(row.size()), static_cast<Eigen::Index>(grads.size()));
    VecX b(static_cast<Eigen::Index>(row.size()));
    for (std::size_t j = 0; j < grads.size(); ++j)
      for (const auto& [f, c] : grads[j])
        if (auto it = row.find(f); it != row.end()) A(it->second, static_cast<Eigen::Index>(j)) += c;
    for (const auto& [f, i] : row) b(i) = fint(f) - cl.nodal_full(f);
    const VecX mu = grads.empty() || row.empty() ? VecX::Zero(static_cast<Eigen::Index>(grads.size()))
                                                : VecX(A.colPivHouseholderQr().solve(b));
    std::size_t j = 0;
    for (; j < model.constraints.size(); ++j) r.constraint_forces.push_back(mu(static_cast<Eigen::Index>(j)));
    r.skew_reactions.resize(model.skew_supports.size());
    for (std::size_t k = 0; k < model.skew_supports.size(); ++k) r.skew_reactions[k].node = model.skew_supports[k].node;
    for (const auto& [k, c] : skew_cols) {
      const double m = mu(static_cast<Eigen::Index>(j++));
      CaseResult::SkewReaction& sr = r.skew_reactions[k];
      sr.local(c) = m;
      sr.global.segment<3>(c < 3 ? 0 : 3) += m * model.skew_supports[k].axes.row(c % 3).transpose();
    }
  }
  // Rigid links: the unbalanced force of each slave DOF (its elements'
  // K u minus its applied load) is the constraint force the link carries to
  // the DOFs of its expression — supported masters report it as reaction.
  for (const LinkReduction::Expr& e : asm_.red.slaves) {
    const double lam = fint(e.slave) - cl.nodal_full(e.slave);
    for (const auto& [k, c] : e.terms) fint(k) += c * lam;
    fint(e.slave) -= lam;
  }

  // Reactions R = K u - F at supported nodes (declaration order, deduped).
  std::vector<Id> seen;
  for (const Support& s : model.supports) {
    if (s.implicit) continue;  // modelling restraint: no reaction reported
    bool dup = false;
    for (Id x : seen)
      if (x == s.node) dup = true;
    if (dup) continue;
    seen.push_back(s.node);
    int ni = model.node_index(s.node);
    Reaction rea;
    rea.node = s.node;
    for (int c = 0; c < 6; ++c) rea.R(c) = fint(d.full(ni, c)) - cl.nodal_full(d.full(ni, c));
    r.reactions.push_back(rea);
  }
  return r;
}

namespace {
// A mechanism whose pivot cancels only to round-off (aligned hinges, a bar
// chain free to rotate) passes the zero-pivot test of the factorization and
// returns huge displacements of either sign. Linear statics refuses it by
// the ratio of each pivot to the diagonal term it eliminates
// (FactorCache::min_pivot_ratio): round-off leaves ratios near 1e-16, a sound
// model, however contrasted its stiffnesses, ratios far above the threshold.
constexpr double kMechanismPivotRatio = 1e-13;
void require_no_mechanism(const FactorCache& cache) {
  const double r = cache.min_pivot_ratio();
  if (r < kMechanismPivotRatio) {
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.1e", r);
    fail(std::string("solver: singular stiffness (pivot ") + buf +
         " of its diagonal term) — the model contains a mechanism, or stiffnesses more than 1e13 apart");
  }
}
}  // namespace

StaticContext StaticContext::build(const Model& model, bool links_supported) {
  StaticContext ctx;
  ctx.model = &model;
  ctx.asm_ = assemble_elements(model, links_supported);
  ctx.K = assemble_stiffness(model, ctx.asm_);
  if (ctx.asm_.red.active) ctx.cache.set_reduction(ctx.asm_.red.T);
  ctx.cache.factorize(ctx.K);
  require_no_mechanism(ctx.cache);
  return ctx;
}

VecX StaticContext::solve_displacements(const CaseLoads& cl) {
  VecX u = cache.solve(cl.rhs);
  if (cl.link_u0.size()) u += cl.link_u0;
  return u;
}

CaseResult StaticContext::solve_case(const LoadCase& lc, const StaticOptions& opt) {
  CaseLoads cl = lower_case(*model, asm_, lc);
  fold_prescribed(*model, asm_, cl);
  fold_link_offsets(*model, asm_, K, cl);
  VecX u = cache.solve(cl.rhs);
  return recover_case(*model, asm_, cl, u, opt, lc.id, lc.name);
}

namespace {

// Any unilateral support: per-DOF one-sided (uplift) OR an ST1 DECOL
// all-or-nothing switch. Both take the active-set path below.
bool has_unilateral(const Model& model) {
  for (const Support& s : model.supports)
    for (int c = 0; c < 6; ++c)
      if (s.oneside[static_cast<std::size_t>(c)] != 0 || s.cap[static_cast<std::size_t>(c)] > 0) return true;
  return has_decol(model);
}

// Active-set solve for unilateral supports until the Kuhn-Tucker conditions
// hold. Two independent flavours, per load case, non-superposable:
//   * per-DOF `oneside`: only the offending DOF is released / re-engaged;
//   * ST1 `DECOL` (manual p. 47, p. 49): the criterion bears on ONE global
//     axis (decol_axis) and, when it is met, EVERY DOF of that support is
//     released; the support re-engages when its displacement on that axis
//     would push back into the support.
CaseResult solve_case_uplift(const Model& model, const LoadCase& lc, const StaticOptions& opt) {
  // engaged[i][c] tracks whether one-sided DOF (support i, comp c) is active.
  std::vector<std::array<bool, 6>> engaged(model.supports.size());
  for (auto& a : engaged) a.fill(true);
  // engaged_sup[i] tracks whether the whole DECOL support i is still glued.
  std::vector<char> engaged_sup(model.supports.size(), 1);
  // capped[i][c]: the sign (+-1) of the reaction a capped DOF yielded under,
  // 0 while it holds: released, the support then pushes with
  // exactly capped x cap.
  std::vector<std::array<int, 6>> capped(model.supports.size());
  for (auto& a : capped) a.fill(0);
  int cap = 3 * static_cast<int>(model.supports.size()) * 6 + 10;

  CaseResult result;
  for (int it = 0; it < cap; ++it) {
    Model work = model;
    LoadCase wlc = lc;
    for (std::size_t i = 0; i < model.supports.size(); ++i)
      for (int c = 0; c < 6; ++c)
        if (capped[i][static_cast<std::size_t>(c)] != 0) {
          const Support& s = model.supports[i];
          const double F = capped[i][static_cast<std::size_t>(c)] * s.cap[static_cast<std::size_t>(c)];
          Vec3 f = Vec3::Zero(), mo = Vec3::Zero();
          (c < 3 ? f : mo)(c % 3) = F;
          wlc.loads.push_back(Load::nodal(s.node, f, mo));
        }
    for (std::size_t i = 0; i < work.supports.size(); ++i) {
      Support& s = work.supports[i];
      if (s.decol_axis >= 0 && !engaged_sup[i]) {
        // Separated: the whole support vanishes this iteration.
        s.fixed.fill(false);
        s.spring.fill(0.0);
        s.spring6.reset();
      }
      for (int c = 0; c < 6; ++c) {
        if (capped[i][static_cast<std::size_t>(c)] != 0) {
          s.fixed[static_cast<std::size_t>(c)] = false;
          s.spring[static_cast<std::size_t>(c)] = 0;
          continue;
        }
        if (s.oneside[static_cast<std::size_t>(c)] == 0) continue;
        if (!engaged[i][static_cast<std::size_t>(c)]) {
          s.fixed[static_cast<std::size_t>(c)] = false;
          s.spring[static_cast<std::size_t>(c)] = 0;
        }
      }
      s.oneside = {};      // the working model is bilateral
      s.cap = {};
      s.decol_axis = -1;
    }
    // Drop supports left with no active DOF (fully released this iteration).
    std::erase_if(work.supports, [](const Support& s) { return !s.any(); });
    StaticContext ctx = StaticContext::build(work);
    CaseResult r = ctx.solve_case(wlc, opt);
    // node -> its LAST entry in r.reactions (what a scan of the list keeping
    // the last match reads), O(1) per support instead of O(reactions)
    std::unordered_map<Id, std::size_t> reaction_of;
    for (std::size_t k = 0; k < r.reactions.size(); ++k) reaction_of[r.reactions[k].node] = k;
    auto reaction_at = [&](Id node, int c) {
      auto it = reaction_of.find(node);
      return it == reaction_of.end() ? 0.0 : r.reactions[it->second].R(c);
    };

    bool changed = false;
    for (std::size_t i = 0; i < model.supports.size(); ++i) {
      const Support& s = model.supports[i];
      int ni = model.node_index(s.node);
      if (s.decol_axis >= 0) {
        int ax = s.decol_axis, side = s.decol_sign;
        if (engaged_sup[i]) {
          // Glued while the reaction on the DECOL axis keeps the admissible
          // sign (POS: R >= 0). One offending component releases everything.
          const double R = reaction_at(s.node, ax);
          if (R * side < -1e-9 * (1.0 + std::fabs(R))) {
            engaged_sup[i] = 0;
            changed = true;
          }
        } else {
          // Separated: re-glue as soon as the node would penetrate again.
          double u = r.node_disp[static_cast<std::size_t>(ni)](ax);
          if (u * side < -1e-12) {
            engaged_sup[i] = 1;
            changed = true;
          }
        }
        continue;  // DECOL and per-DOF oneside are mutually exclusive (validate)
      }
      for (int c = 0; c < 6; ++c) {
        int side = s.oneside[static_cast<std::size_t>(c)];
        const double cp = s.cap[static_cast<std::size_t>(c)];
        if (side == 0 && cp == 0) continue;
        int& yc = capped[i][static_cast<std::size_t>(c)];
        if (yc != 0) {
          // Yielded: the node must keep moving INTO the support (the slip
          // has the sign of the penetration, opposite to the reaction's
          // push); else it holds again.
          double u = r.node_disp[static_cast<std::size_t>(ni)](c);
          if (u * yc > 1e-12) {
            yc = 0;
            changed = true;
          }
        } else if (engaged[i][static_cast<std::size_t>(c)]) {
          // Reaction must have the admissible sign, and stay below the cap.
          const double R = reaction_at(s.node, c);
          if (side != 0 && R * side < -1e-9 * (1.0 + std::fabs(R))) {
            engaged[i][static_cast<std::size_t>(c)] = false;
            changed = true;
          } else if (cp > 0 && std::fabs(R) > cp * (1 + 1e-12)) {
            yc = R > 0 ? 1 : -1;
            changed = true;
          }
        } else {
          // Released DOF must not penetrate: displacement moves away from support.
          double u = r.node_disp[static_cast<std::size_t>(ni)](c);
          if (u * side < -1e-12) {
            engaged[i][static_cast<std::size_t>(c)] = true;
            changed = true;
          }
        }
      }
    }
    if (!changed) {
      r.non_superposable = true;
      // A support whose DOFs are ALL released disappears from the working
      // model, hence from `recover_case`'s reaction list. Restore a zero entry
      // so the reported list always mirrors the user's support declaration
      // order — a separated support simply carries no reaction (every
      // component, moments included, is exactly zero).
      std::vector<Reaction> full;
      std::unordered_set<Id> seen;
      // the supports of each node, in declaration order
      std::unordered_map<Id, std::vector<std::size_t>> supports_of;
      for (std::size_t q = 0; q < model.supports.size(); ++q) supports_of[model.supports[q].node].push_back(q);
      for (const Support& s : model.supports) {
        if (!seen.insert(s.node).second) continue;
        Reaction rea;
        rea.node = s.node;
        if (auto it = reaction_of.find(s.node); it != reaction_of.end()) rea.R = r.reactions[it->second].R;
        // a yielded component carries exactly its cap
        for (std::size_t q : supports_of[s.node])
          for (int c = 0; c < 6; ++c)
            if (capped[q][static_cast<std::size_t>(c)] != 0)
              rea.R(c) = capped[q][static_cast<std::size_t>(c)] * model.supports[q].cap[static_cast<std::size_t>(c)];
        full.push_back(rea);
      }
      r.reactions = std::move(full);
      return r;
    }
    result = r;
  }
  fail("uplift: active-set iteration did not converge for case '" + lc.name + "'");
}

}  // namespace

StaticResults solve_static_cached(const Model& model, const StaticOptions& opt,
                                  FactorCache& cache) {
  hooks_check(opt.hooks);
  model.validate();
  // A uniform Winkler bed is solved exactly (bef.h) on a working copy: the
  // Hermite bed of the element is an approximation whose station forces
  // also leave out the soil reaction.
  if (opt.exact_soil_beds && std::any_of(model.bars.begin(), model.bars.end(),
                                         [&](const Bar& b) { return exact_bed_applies(model, b); })) {
    Model exact = model;
    install_exact_beds(exact);
    return solve_static_cached(exact, opt, cache);
  }
  StaticResults out;

  if (has_unilateral(model)) {
    for (const LoadCase& lc : model.cases) {
      hooks_check(opt.hooks);
      out.cases.push_back(solve_case_uplift(model, lc, opt));
    }
    // Combinations over non-superposable uplift results are rejected.
    require(model.combos.empty(),
            "uplift supports: results are non-superposable, combinations are not allowed");
    return out;
  }

  hooks_report(opt.hooks, 0.0, "assembly");
  Assembled asm_ = assemble_elements(model, true);
  out.auto_fixed = asm_.dofs.auto_fixed;
  SpMat K = assemble_stiffness(model, asm_);
  hooks_check(opt.hooks);
  hooks_report(opt.hooks, 0.5, "factorization");
  if (asm_.red.active)
    cache.set_reduction(asm_.red.T);
  else if (cache.reduced())
    cache.set_reduction(SpMat());
  cache.factorize(K);
  require_no_mechanism(cache);

  // Lower all cases, multi-RHS solve.
  std::vector<CaseLoads> lowered;
  lowered.reserve(model.cases.size());
  MatX rhs(asm_.dofs.n_free, static_cast<Eigen::Index>(model.cases.size()));
  for (std::size_t i = 0; i < model.cases.size(); ++i) {
    lowered.push_back(lower_case(model, asm_, model.cases[i]));
    fold_prescribed(model, asm_, lowered.back());
    fold_link_offsets(model, asm_, K, lowered.back());
    rhs.col(static_cast<Eigen::Index>(i)) = lowered.back().rhs;
  }
  MatX U;
  if (!model.cases.empty()) U = cache.solve(rhs);

  // Per-case recovery, optionally parallel (StaticOptions::threads): every
  // case writes its own pre-allocated slot, `recover_case` is a pure function
  // of (model, asm_, lowered[i], U.col(i)) — bit-determinism by construction,
  // the thread count and scheduling cannot change a single ULP. The model's
  // lazy id index is already built (validate + lower_case above), so the
  // workers only READ shared state. The Eigen multi-RHS solve above stays
  // sequential (already vectorized).
  hooks_check(opt.hooks);
  hooks_report(opt.hooks, 0.8, "case recovery");
  out.cases.resize(model.cases.size());
  int nthreads = opt.threads == 0
                     ? std::max(1, static_cast<int>(std::thread::hardware_concurrency()))
                     : std::max(1, opt.threads);
  auto recover_one = [&](std::size_t i) {
    out.cases[i] = recover_case(model, asm_, lowered[i], U.col(static_cast<Eigen::Index>(i)),
                                opt, model.cases[i].id, model.cases[i].name);
  };
  if (nthreads == 1 || model.cases.size() < 2) {
    for (std::size_t i = 0; i < model.cases.size(); ++i) recover_one(i);
  } else {
    std::vector<std::exception_ptr> errs(model.cases.size());
    std::vector<std::thread> pool;
    pool.reserve(static_cast<std::size_t>(nthreads));
    for (int t = 0; t < nthreads; ++t)
      pool.emplace_back([&, t]() {
        for (std::size_t i = static_cast<std::size_t>(t); i < model.cases.size();
             i += static_cast<std::size_t>(nthreads)) {
          try {
            recover_one(i);
          } catch (...) {
            errs[i] = std::current_exception();
          }
        }
      });
    for (std::thread& th : pool) th.join();
    for (const std::exception_ptr& e : errs)  // lowest failing case index wins
      if (e) std::rethrow_exception(e);
  }

  for (const Combination& cb : model.combos) {
    std::vector<const CaseResult*> srcs;
    std::vector<double> facs;
    for (const auto& [cid, f] : cb.terms) {
      for (std::size_t i = 0; i < model.cases.size(); ++i)
        if (model.cases[i].id == cid) srcs.push_back(&out.cases[i]);
      facs.push_back(f);
    }
    require(srcs.size() == cb.terms.size(), "combination '" + cb.name + "': missing case");
    out.combos.push_back(combine_results(srcs, facs, cb.id, cb.name));
  }
  hooks_report(opt.hooks, 1.0, "done");
  return out;
}

StaticResults solve_static(const Model& model, const StaticOptions& opt) {
  FactorCache cache;
  return solve_static_cached(model, opt, cache);
}

CaseResult combine_results(const std::vector<const CaseResult*>& sources,
                           const std::vector<double>& factors, Id id, const std::string& name) {
  require(!sources.empty() && sources.size() == factors.size(), "combine: bad inputs");
  for (const CaseResult* s : sources)
    require(!s->non_superposable, "combine: source '" + s->name + "' is non-superposable");
  CaseResult r = *sources[0];
  r.id = id;
  r.name = name;
  // zero all numeric fields, then accumulate
  for (auto& nd : r.node_disp) nd.setZero();
  for (auto& re : r.reactions) re.R.setZero();
  for (auto& br : r.bars) {
    br.end_forces_local.setZero();
    for (auto& st : br.stations) {
      st.N = st.Vy = st.Vz = st.T = st.My = st.Mz = 0;
      st.defl.setZero();
      st.rot.setZero();
    }
    for (auto& p : br.soil_p) p.setZero();
    for (auto& m : br.soil_m) m = 0;
  }
  for (auto& sh : r.shells) sh.r = ShellResultants{};
  for (auto& lk : r.links) {
    lk.deformation.setZero();
    lk.force.setZero();
  }
  for (double& m : r.constraint_forces) m = 0;
  for (auto& sr : r.skew_reactions) {
    sr.local.setZero();
    sr.global.setZero();
  }
  for (std::size_t s = 0; s < sources.size(); ++s) {
    const CaseResult& c = *sources[s];
    double f = factors[s];
    // The reactions are checked too: the loop below indexes `c.reactions` with
    // `r.reactions`'s range, so a shorter source would read out of bounds
    // instead of reporting a typed mismatch (a staged model whose active
    // support set changed between the two results).
    require(c.node_disp.size() == r.node_disp.size() && c.bars.size() == r.bars.size() &&
                c.reactions.size() == r.reactions.size() && c.shells.size() == r.shells.size() &&
                c.links.size() == r.links.size(),
            "combine: shape mismatch");
    for (std::size_t i = 0; i < r.node_disp.size(); ++i) r.node_disp[i] += f * c.node_disp[i];
    for (std::size_t i = 0; i < r.links.size(); ++i) {
      r.links[i].deformation += f * c.links[i].deformation;
      r.links[i].force += f * c.links[i].force;
    }
    for (std::size_t i = 0; i < r.reactions.size(); ++i) r.reactions[i].R += f * c.reactions[i].R;
    require(c.constraint_forces.size() == r.constraint_forces.size() && c.skew_reactions.size() == r.skew_reactions.size(),
            "combine: shape mismatch");
    for (std::size_t i = 0; i < r.constraint_forces.size(); ++i) r.constraint_forces[i] += f * c.constraint_forces[i];
    for (std::size_t i = 0; i < r.skew_reactions.size(); ++i) {
      r.skew_reactions[i].local += f * c.skew_reactions[i].local;
      r.skew_reactions[i].global += f * c.skew_reactions[i].global;
    }
    for (std::size_t i = 0; i < r.bars.size(); ++i) {
      r.bars[i].end_forces_local += f * c.bars[i].end_forces_local;
      require(c.bars[i].stations.size() == r.bars[i].stations.size(), "combine: station mismatch");
      for (std::size_t k = 0; k < r.bars[i].stations.size(); ++k) {
        Station& a = r.bars[i].stations[k];
        const Station& b = c.bars[i].stations[k];
        a.N += f * b.N;
        a.Vy += f * b.Vy;
        a.Vz += f * b.Vz;
        a.T += f * b.T;
        a.My += f * b.My;
        a.Mz += f * b.Mz;
        a.defl += f * b.defl;
        a.rot += f * b.rot;
      }
      for (std::size_t k = 0; k < r.bars[i].soil_p.size(); ++k) {
        r.bars[i].soil_p[k] += f * c.bars[i].soil_p[k];
        r.bars[i].soil_m[k] += f * c.bars[i].soil_m[k];
      }
    }
    // Facet shells: only the EIGHT tensor components superpose linearly.
    for (std::size_t i = 0; i < r.shells.size(); ++i) {
      ShellResultants& a = r.shells[i].r;
      const ShellResultants& b = c.shells[i].r;
      a.nx += f * b.nx;
      a.ny += f * b.ny;
      a.nxy += f * b.nxy;
      a.mx += f * b.mx;
      a.my += f * b.my;
      a.mxy += f * b.mxy;
      a.qx += f * b.qx;
      a.qy += f * b.qy;
    }
  }
  // Principal values are NOT linear in the sources: recompute them from the
  // combined tensor (combining m1/m2 directly is the classic wrong answer).
  for (auto& sh : r.shells) sh.r.recompute_principals();
  return r;
}

std::vector<double> axial_field(const CaseResult& r) {
  std::vector<double> N(r.bars.size(), 0.0);
  for (std::size_t i = 0; i < r.bars.size(); ++i)
    N[i] = -r.bars[i].end_forces_local(0);  // N at origin, tension-positive
  return N;
}

std::vector<Vec3> shell_membrane_field(const CaseResult& r) {
  std::vector<Vec3> n(r.shells.size(), Vec3::Zero());
  for (std::size_t i = 0; i < r.shells.size(); ++i)
    n[i] = Vec3(r.shells[i].r.nx, r.shells[i].r.ny, r.shells[i].r.nxy);
  return n;
}

std::vector<StressAtPoint> station_stresses(const Model& model, const Bar& bar,
                                            const Station& st) {
  const Section& sec = model.section(bar.section);
  const SectionProps& p = sec.props;
  std::vector<StressAtPoint> out;
  out.reserve(sec.stress_points.size());
  for (const StressPoint& sp : sec.stress_points) {
    StressAtPoint s;
    s.point = sp.name;
    s.sigma = st.N / p.A;
    if (p.Iy > 0) s.sigma += st.My * sp.z / p.Iy;
    if (p.Iz > 0) s.sigma -= st.Mz * sp.y / p.Iz;
    double tau = 0;
    if (p.Avy > 0) tau += sp.tau_factor_y * st.Vy / p.Avy;
    if (p.Avz > 0) tau += sp.tau_factor_z * st.Vz / p.Avz;
    if (p.J > 0) {
      // St-Venant torsion shear estimate at the outline: tau_t = T r / J
      double rr = std::hypot(sp.y, sp.z);
      tau += st.T * rr / p.J;
    }
    s.tau = tau;
    s.von_mises = std::sqrt(s.sigma * s.sigma + 3 * tau * tau);
    out.push_back(std::move(s));
  }
  return out;
}

}  // namespace keta

#include "keta/tapered.h"
#include "keta/element.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {

// ---------------------------------------------------------------------------
// Gauss-Legendre nodes/weights on [0, 1], generated deterministically by Newton
// iteration on the Legendre polynomials (no tables, any order).
// ---------------------------------------------------------------------------

struct GaussRule {
  std::vector<double> x, w;
};

GaussRule gauss_legendre_01(int n) {
  require(n >= 2 && n <= 256, "tapered: Gauss-Legendre order must be in [2, 256]");
  GaussRule g;
  g.x.resize(static_cast<std::size_t>(n));
  g.w.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    double z = std::cos(kPi * (i + 0.75) / (n + 0.5));  // classical asymptotic seed
    double p0 = 1, dp = 1;
    for (int it = 0; it < 100; ++it) {
      p0 = 1;
      double p1 = 0;
      for (int j = 0; j < n; ++j) {
        double p2 = p1;
        p1 = p0;
        p0 = ((2 * j + 1) * z * p1 - j * p2) / (j + 1);
      }
      dp = n * (z * p0 - p1) / (z * z - 1);
      double dz = p0 / dp;
      z -= dz;
      if (std::fabs(dz) <= 1e-16) break;
    }
    // Map [-1, 1] -> [0, 1]: x = (1 - z)/2, w_[0,1] = w_[-1,1]/2.
    g.x[static_cast<std::size_t>(i)] = 0.5 * (1.0 - z);
    g.w[static_cast<std::size_t>(i)] = 1.0 / ((1.0 - z * z) * dp * dp);
  }
  return g;
}

// ---------------------------------------------------------------------------
// Characteristic interpolation
// ---------------------------------------------------------------------------

// Power carried by the linearly varying rectangular dimension h(x) for each
// characteristic (manual p. 62-64; see the header doc for the citations).
struct Exponents {
  double A, Avy, Avz, Iy, Iz, J, Wely, Welz;
};

Exponents exponents_for(TaperLaw law) {
  // SX/SY/SZ are area-like (p = 1); IX is p = 1 by keta convention (the manual is
  // silent and a rectangle's J is not a power of h); Wel = I/v adds one power of h
  // in the denominator for the inertia that varies cubically.
  if (law == TaperLaw::LinY) return {1, 1, 1, 1, 3, 1, 1, 2};
  return {1, 1, 1, 3, 1, 1, 2, 1};  // LinZ
}

double lin_c(double c0, double c1, double t, double p, const char* what) {
  if (c0 == 0.0 && c1 == 0.0) return 0.0;  // ST1: characteristic not defined
  require(finite(c0) && finite(c1) && c0 > 0 && c1 > 0,
          std::string("CARA VAR LIN: ") + what +
              " must be > 0 at both OR and EX (or 0 at both to leave it undefined)");
  if (p == 1.0) return c0 + (c1 - c0) * t;
  if (p == 3.0) {
    double r = std::cbrt(c0) + (std::cbrt(c1) - std::cbrt(c0)) * t;
    return r * r * r;
  }
  double r0 = std::pow(c0, 1.0 / p), r1 = std::pow(c1, 1.0 / p);
  return std::pow(r0 + (r1 - r0) * t, p);
}

double para_c(double c0, double cm, double c1, double t, const char* what) {
  if (c0 == 0.0 && cm == 0.0 && c1 == 0.0) return 0.0;
  require(finite(c0) && finite(cm) && finite(c1) && c0 > 0 && cm > 0 && c1 > 0,
          std::string("CARA VAR PARA: ") + what +
              " must be > 0 at OR, MI and EX (or 0 at all three to leave it undefined)");
  double l0 = (1 - t) * (1 - 2 * t), l1 = 4 * t * (1 - t), l2 = t * (2 * t - 1);
  double inv = l0 / c0 + l1 / cm + l2 / c1;
  require(inv > 0, std::string("CARA VAR PARA: the interpolated 1/") + what +
                       " parabola becomes non-positive inside the bar (OR/MI/EX are "
                       "not interpolable by a parabola on their inverse)");
  return 1.0 / inv;
}

// Rigid equilibrium transfer G with r = (L, 0, 0): p_i = -G p_j.
Mat6 rigid_G(double L) {
  Mat6 g = Mat6::Identity();
  g(4, 2) = -L;  // M_y at i picks up -L * F_z at j
  g(5, 1) = L;   // M_z at i picks up +L * F_y at j
  return g;
}

// Unit-end-action internal-force fields at the cut x (rows n, vy, vz, t, my, mz;
// columns Fx, Fy, Fz, Mx, My, Mz applied at end j).
Mat6 unit_fields(double L, double x) {
  Mat6 u = Mat6::Identity();
  double xi = L - x;
  u(4, 2) = -xi;
  u(5, 1) = xi;
  return u;
}

// Sorted, de-duplicated quadrature breakpoints (load discontinuities + the ends).
std::vector<double> load_breakpoints(const std::vector<TaperedSpanLoad>& loads, double L) {
  std::vector<double> br{0.0, L};
  for (const TaperedSpanLoad& ld : loads) {
    double a = ld.a;
    double b = ld.kind == TaperedSpanLoad::Kind::Point ? ld.a : (ld.b < 0 ? L : ld.b);
    br.push_back(std::clamp(a, 0.0, L));
    br.push_back(std::clamp(b, 0.0, L));
  }
  std::sort(br.begin(), br.end());
  double tol = 1e-12 * std::max(L, 1.0);
  br.erase(std::unique(br.begin(), br.end(),
                       [tol](double p, double q) { return std::fabs(p - q) <= tol; }),
           br.end());
  return br;
}

}  // namespace

// ---------------------------------------------------------------------------
// Interpolation accessors
// ---------------------------------------------------------------------------

SectionProps tapered_props_at(const TaperedVar& var, double t) {
  require(finite(t) && t >= -1e-9 && t <= 1.0 + 1e-9,
          "tapered: normalized abscissa t = x/L outside [0, 1]");
  t = std::clamp(t, 0.0, 1.0);
  SectionProps p;
  if (var.law == TaperLaw::Para) {
    const SectionProps &o = var.orig, &m = var.mid, &e = var.ext;
    p.A = para_c(o.A, m.A, e.A, t, "SX");
    p.Avy = para_c(o.Avy, m.Avy, e.Avy, t, "SY");
    p.Avz = para_c(o.Avz, m.Avz, e.Avz, t, "SZ");
    p.Iy = para_c(o.Iy, m.Iy, e.Iy, t, "IY");
    p.Iz = para_c(o.Iz, m.Iz, e.Iz, t, "IZ");
    p.J = para_c(o.J, m.J, e.J, t, "IX");
    p.Wely = para_c(o.Wely, m.Wely, e.Wely, t, "Wely");
    p.Welz = para_c(o.Welz, m.Welz, e.Welz, t, "Welz");
    p.Wply = para_c(o.Wply, m.Wply, e.Wply, t, "Wply");
    p.Wplz = para_c(o.Wplz, m.Wplz, e.Wplz, t, "Wplz");
  } else {
    Exponents x = exponents_for(var.law);
    const SectionProps &o = var.orig, &e = var.ext;
    p.A = lin_c(o.A, e.A, t, x.A, "SX");
    p.Avy = lin_c(o.Avy, e.Avy, t, x.Avy, "SY");
    p.Avz = lin_c(o.Avz, e.Avz, t, x.Avz, "SZ");
    p.Iy = lin_c(o.Iy, e.Iy, t, x.Iy, "IY");
    p.Iz = lin_c(o.Iz, e.Iz, t, x.Iz, "IZ");
    p.J = lin_c(o.J, e.J, t, x.J, "IX");
    p.Wely = lin_c(o.Wely, e.Wely, t, x.Wely, "Wely");
    p.Welz = lin_c(o.Welz, e.Welz, t, x.Welz, "Welz");
    p.Wply = lin_c(o.Wply, e.Wply, t, x.Wely, "Wply");
    p.Wplz = lin_c(o.Wplz, e.Wplz, t, x.Welz, "Wplz");
  }
  return p;
}

Vec6 tapered_compliance(const TaperedVar& var, const Material& mat, double t) {
  require(finite(mat.E) && mat.E > 0, "tapered: material E must be > 0");
  double G = mat.G();
  require(finite(G) && G > 0, "tapered: material G = E/(2(1+nu)) must be > 0");
  SectionProps s = tapered_props_at(var, t);
  require(s.A > 0, "CARA VAR: SX must be > 0 everywhere along the bar");
  require(s.Iy > 0 && s.Iz > 0, "CARA VAR: IY and IZ must be > 0 everywhere along the bar");
  require(s.J > 0,
          "CARA VAR: IX (torsion constant) must be > 0 — the flexibility-integrated tapered "
          "stiffness inverts the end flexibility, and a zero torsion constant makes it singular");
  Vec6 d;
  d(0) = 1.0 / (mat.E * s.A);
  d(1) = s.Avy > 0 ? 1.0 / (G * s.Avy) : 0.0;
  d(2) = s.Avz > 0 ? 1.0 / (G * s.Avz) : 0.0;
  d(3) = 1.0 / (G * s.J);
  d(4) = 1.0 / (mat.E * s.Iy);
  d(5) = 1.0 / (mat.E * s.Iz);
  return d;
}

// ---------------------------------------------------------------------------
// Stiffness
// ---------------------------------------------------------------------------

Mat6 tapered_flexibility(const TaperedVar& var, const Material& mat, double L, int gauss) {
  require(finite(L) && L > 0, "tapered: bar length must be > 0");
  GaussRule g = gauss_legendre_01(gauss);
  Mat6 F = Mat6::Zero();
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    double t = g.x[k], w = g.w[k] * L, x = t * L;
    Vec6 d = tapered_compliance(var, mat, t);
    Mat6 u = unit_fields(L, x);
    Mat6 ud = u.transpose() * d.asDiagonal();
    F += w * (ud * u);
  }
  return 0.5 * (F + F.transpose());  // kill the round-off asymmetry
}

Mat6 tapered_end_stiffness(const TaperedVar& var, const Material& mat, double L, int gauss) {
  Mat6 F = tapered_flexibility(var, mat, L, gauss);
  Eigen::LLT<Mat6> llt(F);
  require(llt.info() == Eigen::Success,
          "tapered: the end flexibility is not positive definite — check that SX, SY/SZ, IX, IY "
          "and IZ stay physical along the bar");
  Mat6 K = llt.solve(Mat6::Identity());
  return 0.5 * (K + K.transpose());
}

Mat12 tapered_ke(const TaperedVar& var, const Material& mat, double L, int gauss) {
  Mat6 Kjj = tapered_end_stiffness(var, mat, L, gauss);
  Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
  B.block<6, 6>(0, 0) = -rigid_G(L).transpose();
  B.block<6, 6>(0, 6) = Mat6::Identity();
  Mat12 k = B.transpose() * Kjj * B;
  return 0.5 * (k + k.transpose());
}

// ---------------------------------------------------------------------------
// Span loads
// ---------------------------------------------------------------------------

TaperedSpanLoad TaperedSpanLoad::udl(const Vec3& q, double a, double b) {
  TaperedSpanLoad l;
  l.kind = Kind::Udl;
  l.qa = l.qb = q;
  l.a = a;
  l.b = b;
  return l;
}

TaperedSpanLoad TaperedSpanLoad::linear(const Vec3& qa, const Vec3& qb, double a, double b) {
  TaperedSpanLoad l;
  l.kind = Kind::Linear;
  l.qa = qa;
  l.qb = qb;
  l.a = a;
  l.b = b;
  return l;
}

TaperedSpanLoad TaperedSpanLoad::point(double a, const Vec3& F, const Vec3& M) {
  TaperedSpanLoad l;
  l.kind = Kind::Point;
  l.qa = F;
  l.ma = M;
  l.a = a;
  l.b = a;
  return l;
}

SpanResultant tapered_reduce_loads(const std::vector<TaperedSpanLoad>& loads, double L, double s0,
                                   double s1, double x, bool open_at_s1) {
  require(finite(L) && L > 0, "tapered: bar length must be > 0");
  double tol = 1e-12 * std::max(L, 1.0);
  SpanResultant out;
  for (const TaperedSpanLoad& ld : loads) {
    require(finite(ld.a) && ld.qa.allFinite() && ld.qb.allFinite() && ld.ma.allFinite(),
            "tapered span load: non-finite abscissa or intensity");
    if (ld.kind == TaperedSpanLoad::Kind::Point) {
      require(ld.a >= -tol && ld.a <= L + tol, "tapered span load: point abscissa outside the bar");
      bool inside = ld.a >= s0 - tol && (open_at_s1 ? ld.a < s1 - tol : ld.a <= s1 + tol);
      if (!inside) continue;
      out.F += ld.qa;
      out.M += Vec3(ld.a - x, 0, 0).cross(ld.qa) + ld.ma;
      continue;
    }
    double a = ld.a, b = ld.b < 0 ? L : ld.b;
    require(finite(b) && a >= -tol && b <= L + tol && b >= a - tol,
            "tapered span load: abscissa range outside the bar");
    a = std::clamp(a, 0.0, L);
    b = std::clamp(b, 0.0, L);
    double len = b - a;
    if (len <= tol) continue;  // zero-length distributed load carries nothing
    double lo = std::max(a, s0), hi = std::min(b, s1);
    if (hi <= lo) continue;
    Vec3 slope = (ld.kind == TaperedSpanLoad::Kind::Udl ? Vec3::Zero().eval()
                                                        : ((ld.qb - ld.qa) / len).eval());
    // q(s) = qa + slope * (s - a); integrate over u = s - a in [u0, u1].
    double u0 = lo - a, u1 = hi - a, x0 = x - a;
    double i0 = u1 - u0;
    double i1 = (u1 * u1 - u0 * u0) / 2.0;
    double i2 = (u1 * u1 * u1 - u0 * u0 * u0) / 3.0;
    Vec3 f = ld.qa * i0 + slope * i1;                                  // int q ds
    Vec3 lever = ld.qa * (i1 - x0 * i0) + slope * (i2 - x0 * i1);      // int (s - x) q ds
    out.F += f;
    // int (s - x) e_x x q ds = (0, -lever_z, +lever_y)
    out.M += Vec3(0, -lever.z(), lever.y());
  }
  return out;
}

Vec12 tapered_f0(const TaperedVar& var, const Material& mat, double L,
                 const std::vector<TaperedSpanLoad>& loads, int gauss) {
  require(finite(L) && L > 0, "tapered: bar length must be > 0");
  Mat6 Kjj = tapered_end_stiffness(var, mat, L, gauss);

  // Free-end displacements of the base system under the span load (unit-load
  // theorem), quadrature split at every load breakpoint.
  GaussRule g = gauss_legendre_01(gauss);
  std::vector<double> br = load_breakpoints(loads, L);
  Vec6 d0 = Vec6::Zero();
  for (std::size_t seg = 0; seg + 1 < br.size(); ++seg) {
    double s0 = br[seg], s1 = br[seg + 1], len = s1 - s0;
    if (len <= 0) continue;
    for (std::size_t k = 0; k < g.x.size(); ++k) {
      double x = s0 + g.x[k] * len, w = g.w[k] * len;
      Vec6 d = tapered_compliance(var, mat, x / L);
      Mat6 u = unit_fields(L, x);
      SpanResultant r = tapered_reduce_loads(loads, L, x, L, x);
      Vec6 s;
      s << r.F, r.M;
      d0 += w * (u.transpose() * (d.asDiagonal() * s));
    }
  }

  Vec6 pj = -(Kjj * d0);
  SpanResultant qres = tapered_reduce_loads(loads, L, 0.0, L, 0.0);
  Vec6 q;
  q << qres.F, qres.M;
  Vec6 pi = -(rigid_G(L) * pj) - q;  // member equilibrium: exact by construction
  Vec12 f0;
  f0 << -pi, -pj;
  return f0;
}

Vec12 tapered_f0_strain(const TaperedVar& var, const Material& mat, double L, const Vec3& ea,
                        const Vec3& eb, int gauss) {
  require(finite(L) && L > 0, "tapered: bar length must be > 0");
  require(ea.allFinite() && eb.allFinite(), "tapered: non-finite imposed strain");
  Mat6 Kjj = tapered_end_stiffness(var, mat, L, gauss);
  // Free-end displacement of the base system under the imposed generalized
  // strains (unit-load theorem, no stress): d0 = int U^T e dx, e = (eps, 0, 0,
  // 0, k_my, k_mz). The curvatures enter with the signs that make a constant
  // section reproduce BarElement's prismatic couples (+EI k at i, -EI k at j):
  // k_my = -kappa_y, k_mz = -kappa_z against the my / mz resultants of
  // unit_fields.
  GaussRule g = gauss_legendre_01(gauss);
  Vec6 d0 = Vec6::Zero();
  for (std::size_t k = 0; k < g.x.size(); ++k) {
    const double x = g.x[k] * L, w = g.w[k] * L, t = g.x[k];
    const Vec3 e = (1 - t) * ea + t * eb;
    Vec6 ev;
    ev << e(0), 0, 0, 0, -e(1), -e(2);
    d0 += w * (unit_fields(L, x).transpose() * ev);
  }
  Vec6 pj = -(Kjj * d0);
  Vec6 pi = -(rigid_G(L) * pj);  // self-equilibrated: no applied force
  Vec12 f0;
  f0 << -pi, -pj;
  return f0;
}

Vec12 tapered_f0_udl(const TaperedVar& var, const Material& mat, double L, const Vec3& q, double a,
                     double b, int gauss) {
  return tapered_f0(var, mat, L, {TaperedSpanLoad::udl(q, a, b)}, gauss);
}

Vec12 tapered_f0_linear(const TaperedVar& var, const Material& mat, double L, const Vec3& qa,
                        const Vec3& qb, double a, double b, int gauss) {
  return tapered_f0(var, mat, L, {TaperedSpanLoad::linear(qa, qb, a, b)}, gauss);
}

Vec12 tapered_f0_point(const TaperedVar& var, const Material& mat, double L, double a,
                       const Vec3& F, const Vec3& M, int gauss) {
  return tapered_f0(var, mat, L, {TaperedSpanLoad::point(a, F, M)}, gauss);
}

// ---------------------------------------------------------------------------
// Station recovery
// ---------------------------------------------------------------------------

std::vector<TaperedStation> tapered_stations(const TaperedVar& var, const Material& mat, double L,
                                             const std::vector<TaperedSpanLoad>& loads,
                                             const Vec12& f_end_local, const Vec12& u_local,
                                             int n_stations) {
  require(finite(L) && L > 0, "tapered: bar length must be > 0");
  require(f_end_local.allFinite() && u_local.allFinite(),
          "tapered stations: non-finite end forces or displacements");
  int n = std::max(1, n_stations);
  std::vector<TaperedStation> st(static_cast<std::size_t>(n));

  const Vec3 Fi = f_end_local.segment<3>(0);
  const Vec3 Mi = f_end_local.segment<3>(3);
  // Internal resultants at the cut x, reduced at (x, 0, 0): the NEGATIVE of the
  // reduction of everything acting on [0, x] (end-i member forces + span loads) —
  // exactly the convention of statics.cpp::recover_stations.
  auto cut = [&](double x, bool left) {
    SpanResultant l = tapered_reduce_loads(loads, L, 0.0, x, x, left);
    Vec3 F = Fi + l.F;
    Vec3 M = Mi + Vec3(-x, 0, 0).cross(Fi) + l.M;
    Vec6 s;
    s << -F, -M;
    return s;
  };

  for (int k = 0; k < n; ++k) {
    double x = n == 1 ? L / 2 : L * k / (n - 1);
    Vec6 s = cut(x, false);
    TaperedStation& t = st[static_cast<std::size_t>(k)];
    t.x = x;
    t.N = s(0);
    t.Vy = s(1);
    t.Vz = s(2);
    t.T = s(3);
    t.My = s(4);
    t.Mz = s(5);
  }

  if (n < 2) return st;

  double h = L / (n - 1);
  // Double integration of the curvature with a per-segment Simpson rule (the
  // bending moment is continuous, so the segment midpoint is side-unambiguous);
  // O(h^4) in the station spacing — with a variable EI there is no closed form.
  auto integrate2 = [&](int mom, int shr, double sgn, bool plane_xy, double dend0, double dendL,
                        auto set) {
    std::vector<double> dp(static_cast<std::size_t>(n), 0.0), v(static_cast<std::size_t>(n), 0.0);
    auto curv = [&](double x, bool left) {
      SectionProps s = tapered_props_at(var, x / L);
      double EI = mat.E * (plane_xy ? s.Iz : s.Iy);
      return sgn * cut(x, left)(mom) / EI;
    };
    auto gam = [&](double x, bool left) {
      SectionProps s = tapered_props_at(var, x / L);
      double Av = plane_xy ? s.Avy : s.Avz;
      if (Av <= 0) return 0.0;
      return cut(x, left)(shr) / (mat.G() * Av);
    };
    for (int k = 1; k < n; ++k) {
      std::size_t k0 = static_cast<std::size_t>(k - 1), k1 = static_cast<std::size_t>(k);
      double x0 = st[k0].x, x1 = st[k1].x, xm = 0.5 * (x0 + x1);
      double c0 = curv(x0, false), cm = curv(xm, false), c1 = curv(x1, false);
      v[k1] = v[k0] + h * dp[k0] + h * h / 6.0 * (c0 + 2.0 * cm);
      dp[k1] = dp[k0] + h / 6.0 * (c0 + 4.0 * cm + c1);
      // Shear deflection: V jumps at a concentrated load, so the right endpoint is
      // taken on its LEFT side (value just before the jump).
      v[k1] += h / 6.0 * (gam(x0, false) + 4.0 * gam(xm, false) + gam(x1, true));
    }
    double a0 = dend0 - v[0];
    double a1 = (dendL - (v[static_cast<std::size_t>(n - 1)] + a0)) / L;
    for (int k = 0; k < n; ++k)
      set(k, v[static_cast<std::size_t>(k)] + a0 + a1 * st[static_cast<std::size_t>(k)].x);
  };
  // v'' = +Mz/EIz (x-y plane), w'' = -My/EIy (x-z plane, theta_y = -dw/dx).
  integrate2(5, 1, +1.0, true, u_local(1), u_local(7), [&](int k, double val) {
    st[static_cast<std::size_t>(k)].defl.y() = val;
  });
  integrate2(4, 2, -1.0, false, u_local(2), u_local(8), [&](int k, double val) {
    st[static_cast<std::size_t>(k)].defl.z() = val;
  });

  // Axial: u(x) = u_i + int N/(E A(x)) dx (per-segment Simpson).
  auto eps = [&](double x, bool left) {
    SectionProps s = tapered_props_at(var, x / L);
    return cut(x, left)(0) / (mat.E * s.A);
  };
  double u = u_local(0);
  st[0].defl.x() = u;
  for (int k = 1; k < n; ++k) {
    double x0 = st[static_cast<std::size_t>(k - 1)].x, x1 = st[static_cast<std::size_t>(k)].x;
    double xm = 0.5 * (x0 + x1);
    u += h / 6.0 * (eps(x0, false) + 4.0 * eps(xm, false) + eps(x1, true));
    st[static_cast<std::size_t>(k)].defl.x() = u;
  }
  return st;
}

// ---------------------------------------------------------------------------
// CARA SPECIAL
// ---------------------------------------------------------------------------

bool SpecialLimits::any() const {
  for (double c : cap)
    if (c != 0.0) return true;
  return false;
}

void SpecialLimits::validate() const {
  for (std::size_t i = 0; i < cap.size(); ++i)
    require(finite(cap[i]) && cap[i] >= 0.0,
            "CARA SPECIAL: force threshold " + std::to_string(i) +
                " must be finite and >= 0 (0 = unlimited, the component stays elastic)");
}

std::array<double, 6> special_demands(const Vec12& f_local) {
  require(f_local.allFinite(), "CARA SPECIAL: non-finite member end forces");
  std::array<double, 6> d{};
  for (int c = 0; c < 6; ++c)
    d[static_cast<std::size_t>(c)] = std::max(std::fabs(f_local(c)), std::fabs(f_local(c + 6)));
  return d;
}

std::array<double, 6> special_ratios(const SpecialLimits& lim, const Vec12& f_local) {
  lim.validate();
  std::array<double, 6> d = special_demands(f_local);
  std::array<double, 6> r{};
  for (std::size_t c = 0; c < 6; ++c) {
    double cap = lim.cap[c];
    r[c] = (cap == 0.0 || d[c] <= cap) ? 1.0 : cap / d[c];
  }
  return r;
}

Vec12 special_clamp(const SpecialLimits& lim, const Vec12& f_local) {
  std::array<double, 6> r = special_ratios(lim, f_local);
  Vec12 out = f_local;
  for (int c = 0; c < 6; ++c) {
    double s = r[static_cast<std::size_t>(c)];
    out(c) *= s;
    out(c + 6) *= s;
  }
  return out;
}

std::array<double, 12> special_dof_scales(const SpecialLimits& lim, const Vec12& f_local) {
  std::array<double, 6> r = special_ratios(lim, f_local);
  std::array<double, 12> s{};
  for (int c = 0; c < 6; ++c) {
    double v = std::sqrt(r[static_cast<std::size_t>(c)]);
    s[static_cast<std::size_t>(c)] = v;
    s[static_cast<std::size_t>(c + 6)] = v;
  }
  return s;
}

Mat12 special_secant_stiffness(const Mat12& k_local, const std::array<double, 12>& scales) {
  for (double s : scales)
    require(finite(s) && s >= 0.0, "CARA SPECIAL: secant scale must be finite and >= 0");
  Mat12 out;
  for (int i = 0; i < 12; ++i)
    for (int j = 0; j < 12; ++j)
      out(i, j) = scales[static_cast<std::size_t>(i)] * scales[static_cast<std::size_t>(j)] *
                  k_local(i, j);
  return out;
}

double special_overshoot(const SpecialLimits& lim, const Vec12& f_local) {
  lim.validate();
  std::array<double, 6> d = special_demands(f_local);
  double worst = 0.0;
  for (std::size_t c = 0; c < 6; ++c) {
    if (lim.cap[c] == 0.0) continue;
    worst = std::max(worst, d[c] / lim.cap[c] - 1.0);
  }
  return worst;
}

void install_tapered(Model& m, Id bid, const TaperedVar& var) {
  // the bar through the model index (checked against its slot), else the
  // first bar of that id by a scan: O(1) per call instead of O(bars) — it is
  // called once per tapered bar
  const Bar* found = nullptr;
  try {
    const std::size_t ix = static_cast<std::size_t>(m.bar_index(bid));
    if (m.bars[ix].id == bid) found = &m.bars[ix];
  } catch (const Error&) {
  }
  if (!found) {
    auto it = std::find_if(m.bars.begin(), m.bars.end(), [&](const Bar& b) { return b.id == bid; });
    if (it != m.bars.end()) found = &*it;
  }
  require(found, "tapered bar: unknown bar " + std::to_string(bid));
  const Bar& bar = *found;
  const Material& mat = m.material(bar.material);
  double L = (m.node(bar.nj).x + bar.ecc[1] - m.node(bar.ni).x - bar.ecc[0]).norm();
  TaperedBar tb;
  tb.ke = tapered_ke(var, mat, L);
  tb.props = [var](double t) { return tapered_props_at(var, t); };
  Mat3 R = bar_rotation(m.node(bar.ni).x + bar.ecc[0], m.node(bar.nj).x + bar.ecc[1], bar.roll);
  double rho = mat.rho;
  Id id_for_msg = bid;
  tb.f0 = [var, mat, L, R, rho, id_for_msg, para = var.law == TaperLaw::Para,
           ke = tb.ke](const Load& l) -> Vec12 {
    auto to_local = [&](const Vec3& v) {
      return l.frame == LoadFrame::Local ? v : Vec3(R * v);
    };
    auto proj = [&](const Vec3& v) {
      if (l.frame != LoadFrame::GlobalProjected) return v;
      Vec3 ex = R.row(0);
      return Vec3(v.x() * std::sqrt(std::max(0.0, 1 - ex.x() * ex.x())),
                  v.y() * std::sqrt(std::max(0.0, 1 - ex.y() * ex.y())),
                  v.z() * std::sqrt(std::max(0.0, 1 - ex.z() * ex.z())));
    };
    double a = l.relative ? l.a * L : l.a;
    double b = l.b < 0 ? -1.0 : (l.relative ? l.b * L : l.b);
    switch (l.kind) {
      case Load::Kind::BarUdl:
        return tapered_f0_udl(var, mat, L, to_local(proj(l.F)), a, b);
      case Load::Kind::BarLinear:
        return tapered_f0_linear(var, mat, L, to_local(proj(l.F)), to_local(proj(l.F2)), a, b);
      case Load::Kind::BarPoint:
        return tapered_f0_point(var, mat, L, a, to_local(l.F), to_local(l.M));
      case Load::Kind::BarThermal: {
        const Vec3 e(mat.alpha * l.dt, -mat.alpha * l.gz, mat.alpha * l.gy);  // kappa_y = -alpha gz
        return tapered_f0_strain(var, mat, L, e, e);
      }
      case Load::Kind::Prestrain: {
        require(l.kappa_x == 0, "tapered bar: an imposed twist rate (kappa_x) is not supported on bar " +
                                    std::to_string(id_for_msg));
        const Vec3 e(l.strain, l.kappa_y, l.kappa_z);
        return tapered_f0_strain(var, mat, L, e, e);
      }
      case Load::Kind::PrestrainLinear:
        require(l.kappa_x == 0 && l.kappa_x2 == 0,
                "tapered bar: an imposed twist rate (kappa_x) is not supported on bar " + std::to_string(id_for_msg));
        return tapered_f0_strain(var, mat, L, Vec3(l.strain, l.kappa_y, l.kappa_z),
                                 Vec3(l.strain2, l.kappa_y2, l.kappa_z2));
      case Load::Kind::PrestrainParabolic:
        fail("tapered bar: a parabolic imposed strain is not supported on bar " + std::to_string(id_for_msg) +
             " (linear imposed fields only)");
      case Load::Kind::PrestrainPoint:
        if (l.frame == LoadFrame::Local) {
          // a concentrated dislocation: the stress-free displacement of
          // the free end j it causes (rigid motion beyond it), held by the
          // exact tapered member stiffness — as BarElement, any section law
          const double xa = std::clamp(l.relative ? l.a * L : l.a, 0.0, L);
          const double thz = -l.kappa_z, thy = -l.kappa_y;
          Vec6 tip;
          tip << l.strain, l.dvy + thz * (L - xa), l.dvz - thy * (L - xa), l.dthx, thy, thz;
          return Vec12(ke.block<12, 6>(0, 6) * tip);
        }
        fail("tapered bar: a concentrated dislocation in global axes is not carried on bar " +
             std::to_string(id_for_msg) + " non portee");
      case Load::Kind::SelfWeight: {
        // SX is LINEAR under both LIN laws, so self-weight lowers EXACTLY to
        // a linearly varying span load (tapered.h); under PARA, rho g A(x)
        // is integrated piecewise-linearly on 64 pieces (the interpolation
        // error of A, <= (1/64)^2 / 8 of its curvature, below any printed
        // digit).
        if (!para) {
          Vec3 qa = R * (rho * tapered_props_at(var, 0.0).A * l.g);
          Vec3 qb = R * (rho * tapered_props_at(var, 1.0).A * l.g);
          return tapered_f0_linear(var, mat, L, qa, qb, 0.0, -1.0);
        }
        Vec12 f0 = Vec12::Zero();
        constexpr int kPieces = 64;
        for (int k = 0; k < kPieces; ++k) {
          const double u0 = static_cast<double>(k) / kPieces, u1 = static_cast<double>(k + 1) / kPieces;
          Vec3 qa = R * (rho * tapered_props_at(var, u0).A * l.g);
          Vec3 qb = R * (rho * tapered_props_at(var, u1).A * l.g);
          f0 += tapered_f0_linear(var, mat, L, qa, qb, u0 * L, u1 * L);
        }
        return f0;
      }
      default:
        fail("tapered bar: load kind not carried on bar " + std::to_string(id_for_msg) +
             " (DEFOR CON: documented limitation, tapered.h)");
    }
  };
  m.tapered[bid] = std::move(tb);
}

}  // namespace keta

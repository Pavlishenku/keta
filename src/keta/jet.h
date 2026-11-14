// Forward-mode automatic differentiation to second order ("jets"), for exact
// gradients and Hessians of small closed-form functions of N variables (the
// corotational element energy of nonlinear.cpp: N = 12).
//
// Jet<N> carries a value, its gradient and its dense Hessian with respect to
// N seed variables; every operation applies the exact chain rule
//   unary  f(a):    g = f' ga,          H = f' Ha + f'' ga ga^T
//   product a b:    g = a gb + b ga,    H = a Hb + b Ha + ga gb^T + gb ga^T
// so ONE evaluation of a function written generically on the scalar type
// returns its value, gradient and Hessian to round-off — no truncation error,
// no step size, deterministic (fixed operation order). Grad<N> is the
// first-order variant (value + gradient) for the cheaper force-only
// evaluations. Both are plain aggregates of fixed-size arrays (no heap).
//
// Supported: + - * / (jet-jet, jet-double, double-jet), unary minus, sqrt,
// sin, cos, atan2, and value comparisons through `val`. Functions written for
// the jets also compile for plain double (val(double) is the identity), so one
// template serves the value, the force and the tangent.
#pragma once
#include <array>
#include <cmath>

namespace keta::ad {

template <int N>
struct Grad {
  double v = 0;
  std::array<double, N> g{};
  static Grad seed(double x, int i) {
    Grad r;
    r.v = x;
    r.g[static_cast<std::size_t>(i)] = 1;
    return r;
  }
  static Grad constant(double x) {
    Grad r;
    r.v = x;
    return r;
  }
  // r = f(a) with f' = d1 (the second derivative is not tracked)
  static Grad chain(const Grad& a, double f, double d1, double) {
    Grad r;
    r.v = f;
    for (int i = 0; i < N; ++i) r.g[static_cast<std::size_t>(i)] = d1 * a.g[static_cast<std::size_t>(i)];
    return r;
  }
  // r = f(a, b) with partials fa, fb
  static Grad chain2(const Grad& a, const Grad& b, double f, double fa, double fb, double, double, double) {
    Grad r;
    r.v = f;
    for (std::size_t i = 0; i < N; ++i) r.g[i] = fa * a.g[i] + fb * b.g[i];
    return r;
  }
  friend Grad operator*(const Grad& a, const Grad& b) {
    Grad r;
    r.v = a.v * b.v;
    for (std::size_t i = 0; i < N; ++i) r.g[i] = a.v * b.g[i] + b.v * a.g[i];
    return r;
  }
};

template <int N>
struct Jet {
  double v = 0;
  std::array<double, N> g{};
  std::array<double, N * N> h{};  // row-major, symmetric
  static Jet seed(double x, int i) {
    Jet r;
    r.v = x;
    r.g[static_cast<std::size_t>(i)] = 1;
    return r;
  }
  static Jet constant(double x) {
    Jet r;
    r.v = x;
    return r;
  }
  static Jet chain(const Jet& a, double f, double d1, double d2) {
    Jet r;
    r.v = f;
    for (std::size_t i = 0; i < N; ++i) r.g[i] = d1 * a.g[i];
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j) r.h[i * N + j] = d1 * a.h[i * N + j] + d2 * a.g[i] * a.g[j];
    return r;
  }
  // r = f(a, b): partials fa, fb and second partials faa, fbb, fab.
  static Jet chain2(const Jet& a, const Jet& b, double f, double fa, double fb, double faa, double fbb, double fab) {
    Jet r;
    r.v = f;
    for (std::size_t i = 0; i < N; ++i) r.g[i] = fa * a.g[i] + fb * b.g[i];
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        r.h[i * N + j] = fa * a.h[i * N + j] + fb * b.h[i * N + j] + faa * a.g[i] * a.g[j] + fbb * b.g[i] * b.g[j] +
                         fab * (a.g[i] * b.g[j] + b.g[i] * a.g[j]);
    return r;
  }
  friend Jet operator*(const Jet& a, const Jet& b) {
    Jet r;
    r.v = a.v * b.v;
    for (std::size_t i = 0; i < N; ++i) r.g[i] = a.v * b.g[i] + b.v * a.g[i];
    for (std::size_t i = 0; i < N; ++i)
      for (std::size_t j = 0; j < N; ++j)
        r.h[i * N + j] = a.v * b.h[i * N + j] + b.v * a.h[i * N + j] + a.g[i] * b.g[j] + b.g[i] * a.g[j];
    return r;
  }
};

// ---- generic operations over both jet kinds (T = Grad<N> or Jet<N>) ----

template <class T>
concept JetLike = requires(T t) {
  t.v;
  t.g;
  T::constant(0.0);
};

inline double val(double x) { return x; }
template <JetLike T>
double val(const T& a) {
  return a.v;
}

template <JetLike T>
T operator+(const T& a, const T& b) {
  T r = a;
  r.v += b.v;
  for (std::size_t i = 0; i < a.g.size(); ++i) r.g[i] += b.g[i];
  if constexpr (requires { r.h; })
    for (std::size_t i = 0; i < a.h.size(); ++i) r.h[i] += b.h[i];
  return r;
}
template <JetLike T>
T operator-(const T& a) {
  T r = a;
  r.v = -r.v;
  for (auto& x : r.g) x = -x;
  if constexpr (requires { r.h; })
    for (auto& x : r.h) x = -x;
  return r;
}
template <JetLike T>
T operator-(const T& a, const T& b) {
  return a + (-b);
}
template <JetLike T>
T operator*(const T& a, double s) {
  T r = a;
  r.v *= s;
  for (auto& x : r.g) x *= s;
  if constexpr (requires { r.h; })
    for (auto& x : r.h) x *= s;
  return r;
}
template <JetLike T>
T operator*(double s, const T& a) {
  return a * s;
}
template <JetLike T>
T operator+(const T& a, double s) {
  T r = a;
  r.v += s;
  return r;
}
template <JetLike T>
T operator+(double s, const T& a) {
  return a + s;
}
template <JetLike T>
T operator-(const T& a, double s) {
  return a + (-s);
}
template <JetLike T>
T operator-(double s, const T& a) {
  return (-a) + s;
}
template <JetLike T>
T operator/(const T& a, double s) {
  return a * (1.0 / s);
}
template <JetLike T>
T recip(const T& a) {
  const double x = a.v;
  return T::chain(a, 1.0 / x, -1.0 / (x * x), 2.0 / (x * x * x));
}
template <JetLike T>
T operator/(const T& a, const T& b) {
  return a * recip(b);
}
template <JetLike T>
T operator/(double s, const T& b) {
  return recip(b) * s;
}
template <JetLike T>
T sqrt(const T& a) {
  const double s = std::sqrt(a.v);
  return T::chain(a, s, 0.5 / s, -0.25 / (s * a.v));
}
template <JetLike T>
T sin(const T& a) {
  const double s = std::sin(a.v), c = std::cos(a.v);
  return T::chain(a, s, c, -s);
}
template <JetLike T>
T cos(const T& a) {
  const double s = std::sin(a.v), c = std::cos(a.v);
  return T::chain(a, c, -s, -c);
}
// atan2(y, x): partials y -> x/r2, x -> -y/r2; second partials
// yy = -2xy/r2^2, xx = 2xy/r2^2, xy = (y^2 - x^2)/r2^2.
template <JetLike T>
T atan2(const T& y, const T& x) {
  const double r2 = x.v * x.v + y.v * y.v, r4 = r2 * r2;
  return T::chain2(y, x, std::atan2(y.v, x.v), x.v / r2, -y.v / r2, -2 * x.v * y.v / r4, 2 * x.v * y.v / r4,
                   (y.v * y.v - x.v * x.v) / r4);
}

}  // namespace keta::ad

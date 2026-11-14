#include "keta/wood_armer.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace keta {

namespace {

void check_alpha(double alpha) {
  require(std::isfinite(alpha) && std::sin(alpha) > 1e-6 && alpha > 0 && alpha < kPi,
          "wood_armer: the angle between the two bar directions must lie in (0, pi)");
}

// cos / sin of the bar angle; cos(pi/2) rounds to 6e-17, snapped to the exact
// orthogonal mesh so that orthogonal designs carry no spurious coupling.
void trig(double alpha, double& ca, double& sa) {
  ca = std::cos(alpha);
  sa = std::sin(alpha);
  if (std::fabs(ca) < 1e-14) {
    ca = 0;
    sa = 1;
  }
}

// Capacity tensor X e1 e1' + Y a a' of the two layers.
void capacity_tensor(double X, double Y, double alpha, double& c11, double& c22, double& c12) {
  double ca, sa;
  trig(alpha, ca, sa);
  c11 = X + Y * ca * ca;
  c22 = Y * sa * sa;
  c12 = Y * ca * sa;
}

double min_eig(double a, double b, double c) {
  return 0.5 * (a + b) - std::sqrt(0.25 * (a - b) * (a - b) + c * c);
}

}  // namespace

void rotate2(double& a, double& b, double& c, double beta) {
  const double cs = std::cos(beta), sn = std::sin(beta);
  const double a2 = a * cs * cs + b * sn * sn + 2 * c * sn * cs;
  const double b2 = a * sn * sn + b * cs * cs - 2 * c * sn * cs;
  const double c2 = (b - a) * sn * cs + c * (cs * cs - sn * sn);
  a = a2;
  b = b2;
  c = c2;
}

TwoWaySteel two_way_design(double sxx, double syy, double sxy, double alpha) {
  check_alpha(alpha);
  double ca, sa;
  trig(alpha, ca, sa);
  const double cot = ca / sa;
  // S' = T^-1 S T^-T, T = [e1 a]
  const double p11 = sxx - 2 * sxy * cot + syy * cot * cot;
  const double p12 = (sxy - syy * cot) / sa;
  const double p22 = syy / (sa * sa);
  const double t = std::fabs(p12);
  TwoWaySteel r;
  r.x = p11 + t;
  r.y = p22 + t;
  if (r.x < 0) {
    r.x = 0;
    r.y = p22 + p12 * p12 / -p11;  // p11 < -t <= 0
  } else if (r.y < 0) {
    r.y = 0;
    r.x = p11 + p12 * p12 / -p22;
  }
  r.x = std::max(0.0, r.x);
  r.y = std::max(0.0, r.y);
  return r;
}

WoodArmerMoments wood_armer(double mx, double my, double mxy, double alpha) {
  const TwoWaySteel t = two_way_design(mx, my, mxy, alpha);
  const TwoWaySteel b = two_way_design(-mx, -my, -mxy, alpha);
  WoodArmerMoments r;
  r.mx_top = t.x;
  r.my_top = t.y;
  r.mx_bot = -b.x;
  r.my_bot = -b.y;
  return r;
}

std::vector<WoodArmerMoments> wood_armer_field(const CaseResult& res, double alpha, double beta) {
  std::vector<WoodArmerMoments> out;
  out.reserve(res.shells.size());
  for (const ShellResult& s : res.shells) {
    double a = s.r.mx, b = s.r.my, c = s.r.mxy;
    rotate2(a, b, c, beta);
    out.push_back(wood_armer(a, b, c, alpha));
  }
  return out;
}

std::vector<WoodArmerMoments> wood_armer_envelope(const std::vector<const CaseResult*>& results, double alpha,
                                                  double beta) {
  require(!results.empty(), "wood_armer_envelope: no result");
  std::vector<WoodArmerMoments> env;
  for (const CaseResult* r : results) {
    require(r != nullptr, "wood_armer_envelope: null result");
    const std::vector<WoodArmerMoments> f = wood_armer_field(*r, alpha, beta);
    if (env.empty()) {
      env = f;
      continue;
    }
    require(f.size() == env.size(), "wood_armer_envelope: results of different models");
    for (std::size_t i = 0; i < f.size(); ++i) {
      env[i].mx_top = std::max(env[i].mx_top, f[i].mx_top);
      env[i].my_top = std::max(env[i].my_top, f[i].my_top);
      env[i].mx_bot = std::min(env[i].mx_bot, f[i].mx_bot);
      env[i].my_bot = std::min(env[i].my_bot, f[i].my_bot);
    }
  }
  return env;
}

double face_utilisation(double sxx, double syy, double sxy, double cx, double cy, double alpha, double* theta) {
  check_alpha(alpha);
  require(cx >= 0 && cy >= 0, "wood_armer_assess: capacities must be >= 0");
  double c11, c22, c12;
  capacity_tensor(cx, cy, alpha, c11, c22, c12);
  const double smax = std::max({std::fabs(sxx), std::fabs(syy), std::fabs(sxy)});
  const double cmax = std::max({std::fabs(c11), std::fabs(c22), std::fabs(c12)});
  if (theta) *theta = 0;
  if (smax == 0) return 0;
  const double eps = 1e-9;
  auto feasible = [&](double l) {
    const double sc = smax + l * cmax;
    const double f1 = l * c11 - sxx, f2 = l * c22 - syy, h = l * c12 - sxy;
    return f1 >= -eps * sc && f2 >= -eps * sc && f1 * f2 - h * h >= -eps * sc * sc;
  };
  // candidates: 0, the diagonal roots, the roots of det(l C - S) = 0
  std::vector<double> cand{0.0};
  if (c11 > 0) cand.push_back(sxx / c11);
  if (c22 > 0) cand.push_back(syy / c22);
  const double qa = c11 * c22 - c12 * c12;
  const double qb = -(c11 * syy + c22 * sxx - 2 * c12 * sxy);
  const double qc = sxx * syy - sxy * sxy;
  if (std::fabs(qa) > 1e-14 * cmax * cmax) {
    const double disc = qb * qb - 4 * qa * qc;
    if (disc >= 0) {
      const double q = -0.5 * (qb + std::copysign(std::sqrt(disc), qb));
      if (q != 0) {
        cand.push_back(q / qa);
        cand.push_back(qc / q);
      } else {
        cand.push_back(0.0);
      }
    }
  } else if (qb != 0) {
    cand.push_back(-qc / qb);
  }
  std::sort(cand.begin(), cand.end());
  double lam = std::numeric_limits<double>::infinity();
  for (double l : cand)
    if (l >= 0 && std::isfinite(l) && feasible(l)) {
      lam = l;
      break;
    }
  if (theta && std::isfinite(lam) && lam > 0) {
    // the critical normal: the null direction of lam C - S
    const double a = lam * c11 - sxx, b = lam * c22 - syy, c = lam * c12 - sxy;
    const double e = min_eig(a, b, c);
    // eigenvector of [[a, c], [c, b]] for e: (c, e - a) or (e - b, c)
    double nx = c, ny = e - a;
    if (std::hypot(nx, ny) < std::hypot(e - b, c)) {
      nx = e - b;
      ny = c;
    }
    double t = std::atan2(ny, nx);
    if (t < 0) t += kPi;
    if (t >= kPi) t -= kPi;
    *theta = t;
  }
  return lam;
}

WoodArmerAssessment wood_armer_assess(double mx, double my, double mxy, const SlabCapacity& cap, double alpha) {
  WoodArmerAssessment r;
  r.u_top = face_utilisation(mx, my, mxy, cap.mx_top, cap.my_top, alpha, &r.theta_top);
  r.u_bot = face_utilisation(-mx, -my, -mxy, cap.mx_bot, cap.my_bot, alpha, &r.theta_bot);
  r.u = std::max(r.u_top, r.u_bot);
  return r;
}

ClarkNielsenForces clark_nielsen(double nx, double ny, double nxy, double mx, double my, double mxy, double z_top,
                                 double z_bot, double alpha) {
  require(z_top > z_bot, "clark_nielsen: the top steel level must lie above the bottom one");
  const double dz = z_top - z_bot;
  auto top = [&](double n, double m) { return (m - n * z_bot) / dz; };
  auto bot = [&](double n, double m) { return (n * z_top - m) / dz; };
  ClarkNielsenForces r;
  auto layer = [&](double sx, double sy, double sxy, double& fx, double& fy, double& fc) {
    const TwoWaySteel s = two_way_design(sx, sy, sxy, alpha);
    fx = s.x;
    fy = s.y;
    double c11, c22, c12;
    capacity_tensor(s.x, s.y, alpha, c11, c22, c12);
    fc = std::min(0.0, min_eig(sx - c11, sy - c22, sxy - c12));
  };
  layer(top(nx, mx), top(ny, my), top(nxy, mxy), r.nx_top, r.ny_top, r.fc_top);
  layer(bot(nx, mx), bot(ny, my), bot(nxy, mxy), r.nx_bot, r.ny_bot, r.fc_bot);
  return r;
}

}  // namespace keta

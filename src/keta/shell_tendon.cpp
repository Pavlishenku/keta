#include "keta/shell_tendon.h"

#include <algorithm>
#include <cmath>

namespace keta {

std::vector<Vec3> tendon_polyline(const std::vector<Eigen::Vector2d>& plan, const std::vector<ProfilePoint>& prof,
                                  int per_segment) {
  require(plan.size() >= 2 && prof.size() >= 2 && per_segment >= 1, "tendon profile: invalid data");
  std::vector<double> sp{0.0};
  for (std::size_t i = 1; i < plan.size(); ++i) sp.push_back(sp.back() + (plan[i] - plan[i - 1]).norm());
  for (std::size_t i = 1; i < prof.size(); ++i)
    require(prof[i].s > prof[i - 1].s, "tendon profile: abscissae must increase");
  require(prof.front().s >= -1e-12 && prof.back().s <= sp.back() * (1 + 1e-12),
          "tendon profile: control points outside the plan trace");
  auto z_at = [&](double s) {
    std::size_t i = 0;
    while (i + 2 < prof.size() && s > prof[i + 1].s) ++i;
    const ProfilePoint &a = prof[i], &b = prof[i + 1];
    if (a.flat && b.flat) {
      const double d1 = a.lambda * (b.s - a.s), d2 = b.s - a.s - d1;
      const double c1 = (b.z - a.z) / (d1 * (d1 + d2));
      if (s <= a.s + d1) return a.z + c1 * (s - a.s) * (s - a.s);
      const double c2 = -c1 * d1 / d2;
      return b.z + c2 * (s - b.s) * (s - b.s);
    }
    if (a.flat) return a.z + (b.z - a.z) / ((b.s - a.s) * (b.s - a.s)) * (s - a.s) * (s - a.s);
    if (b.flat) return b.z + (a.z - b.z) / ((a.s - b.s) * (a.s - b.s)) * (s - b.s) * (s - b.s);
    return a.z + (b.z - a.z) * (s - a.s) / (b.s - a.s);
  };
  auto xy_at = [&](double s) {
    std::size_t i = 0;
    while (i + 2 < plan.size() && s > sp[i + 1]) ++i;
    const double u = (s - sp[i]) / (sp[i + 1] - sp[i]);
    return Eigen::Vector2d(plan[i] + u * (plan[i + 1] - plan[i]));
  };
  std::vector<double> ss;
  for (std::size_t i = 0; i + 1 < prof.size(); ++i)
    for (int k = 0; k < per_segment; ++k)
      ss.push_back(prof[i].s + (prof[i + 1].s - prof[i].s) * k / per_segment);
  ss.push_back(prof.back().s);
  for (double s : sp)
    if (s > prof.front().s && s < prof.back().s) ss.push_back(s);
  std::sort(ss.begin(), ss.end());
  std::vector<Vec3> out;
  for (double s : ss) {
    const Eigen::Vector2d p = xy_at(s);
    const Vec3 v(p.x(), p.y(), z_at(s));
    if (out.empty() || (v - out.back()).norm() > 1e-12 * std::max(1.0, sp.back())) out.push_back(v);
  }
  return out;
}

namespace {
struct Profile {
  std::vector<double> s;         // vertex abscissae
  std::vector<Vec3> t;           // segment tangents
  std::vector<double> len;       // segment lengths
  std::vector<double> PSp, PSm;  // start-jacked force after / before each vertex
  std::vector<double> PEp, PEm;  // end-jacked
};
Profile profile_of(const ShellTendon& tn) {
  Profile p;
  const std::size_t n = tn.pts.size();
  p.s.assign(n, 0.0);
  for (std::size_t j = 0; j + 1 < n; ++j) {
    const Vec3 d = tn.pts[j + 1] - tn.pts[j];
    p.len.push_back(d.norm());
    p.t.push_back(d.normalized());
    p.s[j + 1] = p.s[j] + d.norm();
  }
  auto dev = [&](std::size_t i) {  // deviation angle at an interior vertex
    if (i == 0 || i + 1 >= n) return 0.0;
    return std::acos(std::clamp(p.t[i - 1].dot(p.t[i]), -1.0, 1.0));
  };
  const double inf = 0.0;
  p.PSp.assign(n, inf);
  p.PSm.assign(n, inf);
  p.PEp.assign(n, inf);
  p.PEm.assign(n, inf);
  const bool js = tn.jack != ShellTendon::Jack::End, je = tn.jack != ShellTendon::Jack::Start;
  if (js) {
    p.PSm[0] = p.PSp[0] = tn.P0;
    for (std::size_t i = 1; i < n; ++i) {
      p.PSm[i] = p.PSp[i - 1] * std::exp(-tn.k * p.len[i - 1]);
      p.PSp[i] = p.PSm[i] * std::exp(-tn.mu * dev(i));
    }
  }
  if (je) {
    p.PEp[n - 1] = p.PEm[n - 1] = tn.P0;
    for (std::size_t i = n - 1; i-- > 0;) {
      p.PEp[i] = p.PEm[i + 1] * std::exp(-tn.k * p.len[i]);
      p.PEm[i] = p.PEp[i] * std::exp(-tn.mu * dev(i));
    }
  }
  return p;
}
}  // namespace

std::vector<std::pair<double, double>> shell_tendon_forces(const ShellTendon& t) {
  const Profile p = profile_of(t);
  std::vector<std::pair<double, double>> out;
  for (std::size_t i = 0; i < t.pts.size(); ++i)
    out.emplace_back(std::max(p.PSm[i], p.PEm[i]), std::max(p.PSp[i], p.PEp[i]));
  return out;
}

std::vector<std::pair<std::size_t, Vec18>> shell_tendon_loads(const Model& model, const std::vector<ShellElem>& shells,
                                                              const ShellTendon& tn, double factor) {
  const Profile pr = profile_of(tn);
  const std::size_t n = tn.pts.size();
  struct Fac {
    std::size_t k;
    Vec3 x0, n;
    Vec3 x[3];
    double half_t, area2;
  };
  std::vector<Fac> fac;
  for (std::size_t k = 0; k < shells.size(); ++k) {
    const ShellElem& se = shells[k];
    if (tn.slab != 0 && se.slab != tn.slab) continue;
    Fac f;
    f.k = k;
    for (int c = 0; c < 3; ++c) f.x[c] = model.nodes[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(c)])].x;
    f.x0 = f.x[0];
    const Vec3 nn = (f.x[1] - f.x[0]).cross(f.x[2] - f.x[0]);
    f.area2 = nn.norm();
    if (f.area2 == 0) continue;
    f.n = nn / f.area2;
    f.half_t = 0.5 * model.shells[se.index].t;
    fac.push_back(f);
  }
  // barycentric coordinates of the projection of x on a facet
  auto bary = [](const Fac& f, const Vec3& x, double l[3]) {
    for (int c = 0; c < 3; ++c) {
      const Vec3 a = f.x[(c + 1) % 3], b = f.x[(c + 2) % 3];
      l[c] = (b - a).cross(x - a).dot(f.n) / f.area2;
    }
  };
  std::vector<Vec18> acc(shells.size(), Vec18::Zero());
  std::vector<char> touched(shells.size(), 0);
  double size = 0;
  for (const Vec3& x : tn.pts) size = std::max(size, (x - tn.pts[0]).norm());
  const double tol = 1e-10;
  auto add_point = [&](const Fac& f, const Vec3& x, const Vec3& F) {
    double l[3];
    bary(f, x, l);
    const Vec3 r = f.n * f.n.dot(x - f.x0);
    const Vec3 M = r.cross(F);
    for (int c = 0; c < 3; ++c) {
      acc[f.k].segment<3>(6 * c) += l[c] * F;
      acc[f.k].segment<3>(6 * c + 3) += l[c] * M;
    }
    touched[f.k] = 1;
  };
  auto locate = [&](const Vec3& x) -> const Fac* {
    const Fac* best = nullptr;
    double bd = 1e300;
    for (const Fac& f : fac) {
      const double d = std::fabs(f.n.dot(x - f.x0));
      if (d > f.half_t * (1 + 1e-9) + tol * size) continue;
      double l[3];
      bary(f, x, l);
      if (l[0] < -tol || l[1] < -tol || l[2] < -tol) continue;
      if (d < bd) {
        bd = d;
        best = &f;
      }
    }
    require(best != nullptr, "shell tendon " + std::to_string(tn.id) + ": a point lies outside the shells' thickness");
    return best;
  };
  const double s = factor;
  // anchors and deviations
  const std::vector<std::pair<double, double>> P = shell_tendon_forces(tn);
  add_point(*locate(tn.pts[0]), tn.pts[0], s * P[0].second * pr.t[0]);
  add_point(*locate(tn.pts[n - 1]), tn.pts[n - 1], -s * P[n - 1].first * pr.t[n - 2]);
  for (std::size_t i = 1; i + 1 < n; ++i) {
    const Vec3 F = s * (P[i].second * pr.t[i] - P[i].first * pr.t[i - 1]);
    if (F.squaredNorm() > 0) add_point(*locate(tn.pts[i]), tn.pts[i], F);
  }
  // wobble friction along the runs: P'(s) t, piece by piece over the facets
  if (tn.k > 0) {
    static constexpr double gx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                     -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                     0.7966664774136267,  0.9602898564975363};
    static constexpr double gw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                     0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                     0.2223810344533745, 0.1012285362903763};
    for (std::size_t j = 0; j + 1 < n; ++j) {
      const Vec3 a = tn.pts[j], b = tn.pts[j + 1];
      const double L = pr.len[j];
      // P(u) on the run: the larger of the start- and end-jacked curves
      auto PS = [&](double u) { return pr.PSp[j] * std::exp(-tn.k * L * u); };
      auto PE = [&](double u) { return pr.PEm[j + 1] * std::exp(-tn.k * L * (1 - u)); };
      auto dP = [&](double u) {  // dP/ds
        const double ps = PS(u), pe = PE(u);
        return ps >= pe ? -tn.k * ps : tn.k * pe;
      };
      std::vector<double> cuts{0.0, 1.0};
      if (pr.PSp[j] > 0 && pr.PEm[j + 1] > 0) {
        // PS = PE: exp(-k L u) PSp = exp(-k L (1 - u)) PEm
        const double u = 0.5 + std::log(pr.PSp[j] / pr.PEm[j + 1]) / (2 * tn.k * L);
        if (u > 0 && u < 1) cuts.push_back(u);
      }
      // the facet pieces: the parameter interval of each facet
      std::vector<std::tuple<double, double, const Fac*>> iv;
      for (const Fac& f : fac) {
        double la[3], lb[3];
        bary(f, a, la);
        bary(f, b, lb);
        double lo = 0, hi = 1;
        for (int c = 0; c < 3 && lo <= hi; ++c) {
          const double d = lb[c] - la[c];
          if (std::fabs(d) < 1e-15) {
            if (la[c] < -tol) hi = -1;
            continue;
          }
          const double u = -la[c] / d;
          if (d > 0)
            lo = std::max(lo, u);
          else
            hi = std::min(hi, u);
        }
        if (hi - lo <= 1e-14) continue;
        const Vec3 xm = a + 0.5 * (lo + hi) * (b - a);
        if (std::fabs(f.n.dot(xm - f.x0)) > f.half_t * (1 + 1e-9) + tol * size) continue;
        iv.emplace_back(lo, hi, &f);
      }
      std::sort(iv.begin(), iv.end(), [](const auto& x, const auto& y) { return std::get<0>(x) < std::get<0>(y); });
      // a disjoint cover of [0, 1]
      double cur = 0;
      std::vector<std::tuple<double, double, const Fac*>> pieces;
      while (cur < 1 - 1e-14) {
        const Fac* bf = nullptr;
        double reach = cur;
        for (const auto& [lo, hi, f] : iv)
          if (lo <= cur + 1e-12 && hi > reach) {
            reach = hi;
            bf = f;
          }
        require(bf != nullptr, "shell tendon " + std::to_string(tn.id) + ": a run leaves the shells");
        pieces.emplace_back(cur, reach, bf);
        cur = reach;
      }
      for (const auto& [lo, hi, f] : pieces) {
        std::vector<double> sub{lo, hi};
        for (double c : cuts)
          if (c > lo && c < hi) sub.push_back(c);
        std::sort(sub.begin(), sub.end());
        for (std::size_t q = 0; q + 1 < sub.size(); ++q) {
          const double u0 = sub[q], u1 = sub[q + 1], hm = 0.5 * (u1 - u0), um = 0.5 * (u0 + u1);
          for (int g = 0; g < 8; ++g) {
            const double u = um + hm * gx[g];
            const Vec3 x = a + u * (b - a);
            add_point(*f, x, s * dP(u) * pr.t[j] * (L * hm * gw[g]));
          }
        }
      }
    }
  }
  std::vector<std::pair<std::size_t, Vec18>> out;
  for (std::size_t k = 0; k < shells.size(); ++k)
    if (touched[k]) out.emplace_back(k, acc[k]);
  return out;
}

}  // namespace keta

namespace keta {

namespace {
constexpr double kGx[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
                           0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
constexpr double kGw[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
                           0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763};
}  // namespace

double shell_tendon_mean_force(const ShellTendon& tn) {
  const Profile pr = profile_of(tn);
  double sum = 0, len = 0;
  for (std::size_t j = 0; j + 1 < tn.pts.size(); ++j) {
    const double L = pr.len[j];
    auto P = [&](double u) {
      return std::max(pr.PSp[j] * std::exp(-tn.k * L * u), pr.PEm[j + 1] * std::exp(-tn.k * L * (1 - u)));
    };
    std::vector<double> cuts{0.0, 1.0};
    if (tn.k > 0 && pr.PSp[j] > 0 && pr.PEm[j + 1] > 0) {
      const double u = 0.5 + std::log(pr.PSp[j] / pr.PEm[j + 1]) / (2 * tn.k * L);
      if (u > 0 && u < 1) cuts.push_back(u);
    }
    std::sort(cuts.begin(), cuts.end());
    for (std::size_t q = 0; q + 1 < cuts.size(); ++q) {
      const double hm = 0.5 * (cuts[q + 1] - cuts[q]), um = 0.5 * (cuts[q] + cuts[q + 1]);
      for (int g = 0; g < 8; ++g) sum += P(um + hm * kGx[g]) * L * hm * kGw[g];
    }
    len += L;
  }
  return len > 0 ? sum / len : 0.0;
}

double shell_tendon_mean_strain(const Model& model, const std::vector<ShellElem>& shells, const ShellTendon& tn,
                                const std::vector<Vec6>& node_disp) {
  struct Fac {
    std::size_t k;
    Vec3 x0, n;
    Vec3 x[3];
    double half_t, area2, ecc;
  };
  std::vector<Fac> fac;
  for (std::size_t k = 0; k < shells.size(); ++k) {
    const ShellElem& se = shells[k];
    if (tn.slab != 0 && se.slab != tn.slab) continue;
    Fac f;
    f.k = k;
    for (int c = 0; c < 3; ++c) f.x[c] = model.nodes[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(c)])].x;
    f.x0 = f.x[0];
    const Vec3 nn = (f.x[1] - f.x[0]).cross(f.x[2] - f.x[0]);
    f.area2 = nn.norm();
    if (f.area2 == 0) continue;
    f.n = nn / f.area2;
    f.half_t = 0.5 * se.el.thickness();
    f.ecc = model.shells[se.index].enrich.ecc;
    fac.push_back(f);
  }
  auto bary = [](const Fac& f, const Vec3& x, double l[3]) {
    for (int c = 0; c < 3; ++c) {
      const Vec3 a = f.x[(c + 1) % 3], b = f.x[(c + 2) % 3];
      l[c] = (b - a).cross(x - a).dot(f.n) / f.area2;
    }
  };
  double size = 0;
  for (const Vec3& x : tn.pts) size = std::max(size, (x - tn.pts[0]).norm());
  const double tol = 1e-10;
  double sum = 0, len = 0;
  for (std::size_t j = 0; j + 1 < tn.pts.size(); ++j) {
    const Vec3 a = tn.pts[j], b = tn.pts[j + 1];
    const double L = (b - a).norm();
    if (L == 0) continue;
    const Vec3 t = (b - a) / L;
    std::vector<std::tuple<double, double, const Fac*>> iv;
    for (const Fac& f : fac) {
      double la[3], lb[3];
      bary(f, a, la);
      bary(f, b, lb);
      double lo = 0, hi = 1;
      for (int c = 0; c < 3 && lo <= hi; ++c) {
        const double d = lb[c] - la[c];
        if (std::fabs(d) < 1e-15) {
          if (la[c] < -tol) hi = -1;
          continue;
        }
        const double u = -la[c] / d;
        if (d > 0)
          lo = std::max(lo, u);
        else
          hi = std::min(hi, u);
      }
      if (hi - lo <= 1e-14) continue;
      const Vec3 xm = a + 0.5 * (lo + hi) * (b - a);
      if (std::fabs(f.n.dot(xm - f.x0) - f.ecc) > f.half_t * (1 + 1e-9) + tol * size) continue;
      iv.emplace_back(lo, hi, &f);
    }
    std::sort(iv.begin(), iv.end(), [](const auto& x, const auto& y) { return std::get<0>(x) < std::get<0>(y); });
    double cur = 0;
    while (cur < 1 - 1e-14) {
      const Fac* bf = nullptr;
      double reach = cur;
      for (const auto& [lo, hi, f] : iv)
        if (lo <= cur + 1e-12 && hi > reach) {
          reach = hi;
          bf = f;
        }
      require(bf != nullptr, "shell tendon " + std::to_string(tn.id) + ": a run leaves the shells");
      const ShellElem& se = shells[bf->k];
      Vec18 d;
      for (int c = 0; c < 3; ++c)
        d.segment<6>(6 * c) = node_disp[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(c)])];
      // the tendon direction in the facet axes (its in-plane components)
      const Vec3 tl = se.el.R() * t;
      const double hm = 0.5 * (reach - cur), um = 0.5 * (reach + cur);
      for (int g = 0; g < 8; ++g) {
        const double u = um + hm * kGx[g];
        const Vec3 x = a + u * (b - a);
        double l[3];
        bary(*bf, x, l);
        const Vec3 e = se.el.strain_at(d, {l[0], l[1], l[2]}, bf->n.dot(x - bf->x0) - bf->ecc);
        sum += (tl(0) * tl(0) * e(0) + tl(1) * tl(1) * e(1) + tl(0) * tl(1) * e(2)) * L * hm * kGw[g];
      }
      cur = reach;
    }
    len += L;
  }
  return len > 0 ? sum / len : 0.0;
}

}  // namespace keta

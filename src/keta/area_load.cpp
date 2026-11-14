#include "keta/area_load.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
using P2 = Eigen::Vector2d;
double cross2(const P2& a, const P2& b) { return a.x() * b.y() - a.y() * b.x(); }

// clip a convex CCW polygon by the half-plane a . p <= c
std::vector<P2> clip_half(const std::vector<P2>& poly, const P2& a, double c) {
  std::vector<P2> out;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const P2 &p = poly[i], &q = poly[(i + 1) % poly.size()];
    const double sp = c - a.dot(p), sq = c - a.dot(q);
    if (sp >= 0) out.push_back(p);
    if ((sp >= 0) != (sq >= 0)) out.push_back(p + (q - p) * (sp / (sp - sq)));
  }
  return out;
}
// the interval of the line {o + lambda dir} inside a convex polygon (empty: lo > hi)
void chord(const std::vector<P2>& poly, const P2& o, const P2& dir, double& lo, double& hi) {
  lo = -1e300;
  hi = 1e300;
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i) {
    const P2 &a = poly[i], &b = poly[(i + 1) % n];
    const P2 e = b - a;
    // inside: cross(e, p - a) >= 0 (CCW)
    const double c0 = cross2(e, o - a), c1 = cross2(e, dir);
    if (c1 == 0) {
      if (c0 < -1e-12 * e.norm() * (o - a).norm()) {
        lo = 1;
        hi = 0;
        return;
      }
      continue;
    }
    const double t = -c0 / c1;
    if (c1 > 0)
      lo = std::max(lo, t);
    else
      hi = std::min(hi, t);
  }
}
}  // namespace

std::vector<Load> area_load_to_bars(const Model& model, const AreaLoadZone& z) {
  model.index();
  const std::size_t n = z.poly.size();
  require(n >= 3, "area load: a panel needs 3 vertices at least");
  require(z.q.allFinite(), "area load: non-finite load");
  // plane: Newell normal, (e1, e2) in it
  Vec3 nrm = Vec3::Zero();
  for (std::size_t i = 0; i < n; ++i) nrm += z.poly[i].cross(z.poly[(i + 1) % n]);
  require(nrm.norm() > 0, "area load: degenerate panel");
  nrm.normalize();
  const Vec3 o = z.poly[0];
  const Vec3 e1 = (z.poly[1] - z.poly[0]).normalized(), e2 = nrm.cross(e1);
  double size = 0;
  for (const Vec3& p : z.poly) size = std::max(size, (p - o).norm());
  for (const Vec3& p : z.poly)
    require(std::fabs((p - o).dot(nrm)) <= 1e-9 * size, "area load: the panel vertices are not coplanar");
  std::vector<P2> P;
  for (const Vec3& p : z.poly) P.emplace_back((p - o).dot(e1), (p - o).dot(e2));  // CCW about nrm
  for (std::size_t i = 0; i < n; ++i)
    require(cross2(P[(i + 1) % n] - P[i], P[(i + 2) % n] - P[(i + 1) % n]) >= -1e-12 * size * size,
            "area load: the panel must be convex");
  P2 d2 = P2::Zero();
  if (z.mode == AreaLoadZone::Mode::OneWay) {
    d2 = P2(z.span.dot(e1), z.span.dot(e2));
    require(d2.norm() > 1e-12 * z.span.norm(), "area load: the span direction is normal to the panel");
    d2.normalize();
  }
  // candidate bars
  std::vector<const Bar*> cand;
  for (const Bar& b : model.bars)
    if (z.bars.empty() || std::find(z.bars.begin(), z.bars.end(), b.id) != z.bars.end()) cand.push_back(&b);
  std::vector<Load> out;
  const double tol = 1e-9 * size;
  for (std::size_t k = 0; k < n; ++k) {
    const P2 a = P[k], b = P[(k + 1) % n];
    const double len = (b - a).norm();
    if (len <= tol) continue;
    const P2 t = (b - a) / len, nn(-t.y(), t.x());  // inward normal (CCW)
    // breakpoints s and depths h(s)
    std::vector<double> sv{0.0, len};
    std::vector<P2> region;
    if (z.mode == AreaLoadZone::Mode::TwoWay) {
      region = P;
      for (std::size_t f = 0; f < n && !region.empty(); ++f) {
        if (f == k) continue;
        const P2 af = P[f], bf = P[(f + 1) % n];
        const double lf = (bf - af).norm();
        if (lf <= tol) continue;
        const P2 tf = (bf - af) / lf, nf(-tf.y(), tf.x());
        if ((nf - nn).norm() < 1e-14) continue;  // collinear edges: the same distance
        // nn . (p - a) <= nf . (p - af)
        region = clip_half(region, nn - nf, nn.dot(a) - nf.dot(af));
      }
      for (const P2& v : region) sv.push_back(t.dot(v - a));
    } else {
      if (std::fabs(cross2(d2, t)) < 1e-12) continue;  // an edge along the span takes nothing
      // the line through each panel vertex along d meets the edge at s
      for (const P2& v : P) {
        // a + s t = v + mu d  ->  s = cross(v - a, d) / cross(t, d)
        sv.push_back(cross2(v - a, d2) / cross2(t, d2));
      }
    }
    std::sort(sv.begin(), sv.end());
    std::vector<double> s2;
    for (double s : sv)
      if (s >= -tol && s <= len + tol && (s2.empty() || s - s2.back() > tol)) s2.push_back(std::clamp(s, 0.0, len));
    auto depth = [&](double s) {
      const P2 p = a + s * t;
      double lo, hi;
      if (z.mode == AreaLoadZone::Mode::TwoWay) {
        if (region.size() < 3) return 0.0;
        chord(region, p, nn, lo, hi);
        return hi > lo ? std::max(0.0, hi) - std::max(0.0, lo) : 0.0;
      }
      chord(P, p, d2, lo, hi);
      return hi > lo ? 0.5 * (hi - lo) * std::fabs(cross2(d2, t)) : 0.0;
    };
    std::vector<double> hv;
    for (double s : s2) hv.push_back(depth(s));
    // the bars on the edge, and their coverage
    std::vector<std::pair<double, double>> covered;
    for (const Bar* bar : cand) {
      const Vec3 xi = model.node(bar->ni).x, xj = model.node(bar->nj).x;
      auto loc = [&](const Vec3& x) { return P2((x - o).dot(e1), (x - o).dot(e2)); };
      const P2 pi = loc(xi), pj = loc(xj);
      if (std::fabs((xi - o).dot(nrm)) > tol || std::fabs((xj - o).dot(nrm)) > tol) continue;
      if (std::fabs(nn.dot(pi - a)) > tol || std::fabs(nn.dot(pj - a)) > tol) continue;
      const double sa = t.dot(pi - a), sb = t.dot(pj - a);
      const double lo = std::max(0.0, std::min(sa, sb)), hi = std::min(len, std::max(sa, sb));
      if (hi - lo <= tol) continue;
      covered.emplace_back(lo, hi);
      const double L = std::fabs(sb - sa);
      auto xbar = [&](double s) { return std::clamp(sa < sb ? s - sa : sa - s, 0.0, L); };
      // pieces between the breakpoints inside [lo, hi]
      std::vector<double> cuts{lo};
      for (double s : s2)
        if (s > lo + tol && s < hi - tol) cuts.push_back(s);
      cuts.push_back(hi);
      auto h_at = [&](double s) {
        // piecewise linear interpolation of (s2, hv)
        std::size_t i = 0;
        while (i + 2 < s2.size() && s > s2[i + 1]) ++i;
        const double u = s2[i + 1] > s2[i] ? (s - s2[i]) / (s2[i + 1] - s2[i]) : 0.0;
        return hv[i] + u * (hv[i + 1] - hv[i]);
      };
      for (std::size_t c = 0; c + 1 < cuts.size(); ++c) {
        const double s0 = cuts[c], s1 = cuts[c + 1];
        const double h0 = h_at(s0), h1 = h_at(s1);
        if (h0 == 0 && h1 == 0) continue;
        double x0 = xbar(s0), x1 = xbar(s1);
        Vec3 q0 = z.q * h0, q1 = z.q * h1;
        if (x0 > x1) {
          std::swap(x0, x1);
          std::swap(q0, q1);
        }
        out.push_back(Load::linear(bar->id, q0, q1, LoadFrame::Global, x0, x1));
      }
    }
    // every loaded stretch of the edge must be carried
    std::sort(covered.begin(), covered.end());
    double reach = 0;
    for (const auto& [lo, hi] : covered) {
      if (lo > reach + tol) break;
      reach = std::max(reach, hi);
    }
    double hmax = 0;
    for (double h : hv) hmax = std::max(hmax, h);
    require(reach >= len - tol || hmax == 0,
            "area load: edge " + std::to_string(k) + " of the panel is not carried by bars along its length");
  }
  return out;
}

}  // namespace keta

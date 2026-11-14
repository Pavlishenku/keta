#include "keta/prestress.h"

#include <algorithm>
#include <iterator>
#include <limits>

namespace keta {

namespace {
const PolySegment* find_seg(const Tendon& t, double s) {
  if (t.profile.empty()) return nullptr;
  // the segments run along s: the first one whose end reaches s (binary
  // search), the same segment the former linear scan returned
  auto it = std::lower_bound(t.profile.begin(), t.profile.end(), s,
                             [](const PolySegment& seg, double x) { return seg.s1 + 1e-12 < x; });
  if (it != t.profile.end() && s >= it->s0 - 1e-12) {
    // a joint shared with the previous segment belongs to the previous one
    if (it != t.profile.begin() && s <= std::prev(it)->s1 + 1e-12) return &*std::prev(it);
    return &*it;
  }
  for (const PolySegment& seg : t.profile)  // unsorted or gapped profile
    if (s >= seg.s0 - 1e-12 && s <= seg.s1 + 1e-12) return &seg;
  return s < t.profile.front().s0 ? &t.profile.front() : &t.profile.back();
}
double poly(const std::array<double, 4>& c, double x) {
  return c[0] + x * (c[1] + x * (c[2] + x * c[3]));
}
double dpoly(const std::array<double, 4>& c, double x) {
  return c[1] + x * (2 * c[2] + x * 3 * c[3]);
}
double d2poly(const std::array<double, 4>& c, double x) { return 2 * c[2] + 6 * c[3] * x; }

// Value, first and second derivative (d/ds) of e_y (c = 1) or e_z (c = 2) on
// a segment at s. An exact curve (PolySegment::exact) is read at the
// parameter where its abscissa is s (Newton on the monotone s(u), safeguarded
// by bisection): e' = Y'/X', e'' = (Y'' X' - Y' X'') / X'^3.
struct SegEval {
  double v[3], d1[3], d2[3];
};
SegEval seg_eval(const PolySegment& g, double s) {
  SegEval r{};
  if (g.exact == 0) {
    const double x = s - g.s0;
    r.v[1] = poly(g.cy, x);
    r.d1[1] = dpoly(g.cy, x);
    r.d2[1] = d2poly(g.cy, x);
    r.v[2] = poly(g.cz, x);
    r.d1[2] = dpoly(g.cz, x);
    r.d2[2] = d2poly(g.cz, x);
    return r;
  }
  // the curve X(u), X'(u), X''(u) on [lo, hi]
  double lo = 0, hi = g.exact == 1 ? 1.0 : g.th;
  auto curve = [&](double u, Vec3& x, Vec3& dx, Vec3& ddx) {
    if (g.exact == 1) {
      x = g.p[0] + u * (g.p[1] + u * (g.p[2] + u * g.p[3]));
      dx = g.p[1] + u * (2 * g.p[2] + u * 3 * g.p[3]);
      ddx = 2 * g.p[2] + 6 * u * g.p[3];
    } else {
      const double c = std::cos(u), sn = std::sin(u);
      x = g.C - g.R * c * g.U + g.R * sn * g.D;
      dx = g.R * sn * g.U + g.R * c * g.D;
      ddx = g.R * c * g.U - g.R * sn * g.D;
    }
  };
  Vec3 x, dx, ddx;
  const double span = hi - lo;
  double u = lo + span * std::clamp((s - g.s0) / (g.s1 - g.s0), 0.0, 1.0);
  for (int it = 0; it < 100; ++it) {
    curve(u, x, dx, ddx);
    const double f = x.x() - s;
    if (std::abs(f) <= 1e-14 * std::max(1.0, std::abs(s))) break;
    (f > 0 ? hi : lo) = u;
    double un = dx.x() > 0 ? u - f / dx.x() : 0.5 * (lo + hi);
    if (!(un > lo && un < hi)) un = 0.5 * (lo + hi);
    if (std::abs(un - u) <= 1e-15 * span) {
      u = un;
      curve(u, x, dx, ddx);
      break;
    }
    u = un;
  }
  const double xp = dx.x();
  for (int c = 1; c <= 2; ++c) {
    r.v[c] = x[c];
    r.d1[c] = dx[c] / xp;
    r.d2[c] = (ddx[c] * xp - dx[c] * ddx.x()) / (xp * xp * xp);
  }
  return r;
}

// Slope breaks of the trace strictly inside (a, b): the joints between
// profile segments where e_y' or e_z' jumps (ALIGNE / PENTE polygons). The
// cable is straight on each side, so the deviation there is CONCENTRATED: a
// point force P (t+ - t-) and the friction drop exp(-mu dtheta) at the joint.
std::vector<double> trace_kinks(const Tendon& t, double a, double b) {
  std::vector<double> out;
  const double tol = 1e-9 * std::max(1.0, b - a);
  for (std::size_t i = 0; i + 1 < t.profile.size(); ++i) {
    const PolySegment &l = t.profile[i], &r = t.profile[i + 1];
    if (std::abs(l.s1 - r.s0) > 1e-9 || l.s1 <= a + tol || l.s1 >= b - tol) continue;
    const SegEval el = seg_eval(l, l.s1), er = seg_eval(r, r.s0);
    const double jy = er.d1[1] - el.d1[1], jz = er.d1[2] - el.d1[2];
    if (std::abs(jy) + std::abs(jz) > 1e-9) out.push_back(l.s1);
  }
  std::sort(out.begin(), out.end());
  return out;
}
}  // namespace


double tendon_ey(const Tendon& t, double s) {
  const PolySegment* g = find_seg(t, s);
  require(g, "tendon: empty profile");
  return seg_eval(*g, s).v[1];
}
double tendon_ez(const Tendon& t, double s) {
  const PolySegment* g = find_seg(t, s);
  require(g, "tendon: empty profile");
  return seg_eval(*g, s).v[2];
}
double tendon_dey(const Tendon& t, double s) { return seg_eval(*find_seg(t, s), s).d1[1]; }
double tendon_dez(const Tendon& t, double s) { return seg_eval(*find_seg(t, s), s).d1[2]; }
double tendon_d2ey(const Tendon& t, double s) { return seg_eval(*find_seg(t, s), s).d2[1]; }
double tendon_d2ez(const Tendon& t, double s) { return seg_eval(*find_seg(t, s), s).d2[2]; }

double tendon_chain_length(const Model& model, const Tendon& t) {
  double L = 0;
  for (Id b : t.bars) L += model.bar_chord_length(model.bar(b));
  return L;
}

double ForceProfile::at(double x) const {
  if (s.empty()) return 0;
  if (x <= s.front()) return P.front();
  if (x >= s.back()) return P.back();
  // the first breakpoint at or beyond x (binary search)
  const std::size_t i = static_cast<std::size_t>(std::lower_bound(s.begin() + 1, s.end(), x) - s.begin());
  if (i < s.size()) {
    double t = (x - s[i - 1]) / (s[i] - s[i - 1]);
    return P[i - 1] + t * (P[i] - P[i - 1]);
  }
  return P.back();
}
double ForceProfile::slope(double x) const {
  if (s.size() < 2) return 0;
  for (std::size_t i = 1; i < s.size(); ++i)
    if (x <= s[i] || i + 1 == s.size()) return (P[i] - P[i - 1]) / (s[i] - s[i - 1]);
  return 0;
}

ForceProfile tendon_force_profile(const Model& model, const Tendon& t, int samples_per_meter) {
  ForceProfile fp;
  if (!t.force_profile.empty()) {
    // Authored profile, used as given (friction / draw-in are part of what
    // was authored). ForceProfile::at holds the end values beyond the ends.
    for (const auto& [ps, pp] : t.force_profile) {
      fp.s.push_back(ps);
      fp.P.push_back(pp);
    }
    return fp;
  }
  const double Lchain = tendon_chain_length(model, t);
  // The friction / draw-in law runs between the anchors [sb, se]; `L` below
  // is that anchored length and grid abscissae are offset by sb.
  const double sb = std::clamp(t.s_begin, 0.0, Lchain);
  const double se = t.s_end < 0 ? Lchain : std::clamp(t.s_end, sb, Lchain);
  const double L = se - sb;
  require(L > 0, "tendon: anchored length must be > 0");
  // The grid: samples_per_meter nodes per metre, plus every slope break of
  // the trace as a DOUBLE node (left and right limits at the same abscissa):
  // the friction drop of a concentrated deviation is a jump of the force,
  // not a ramp over the neighbouring step.
  std::vector<double> sg;
  {
    const int nu = std::max(4, static_cast<int>(L * samples_per_meter));
    const std::vector<double> kinks = trace_kinks(t, sb, se);
    const double snap = 1e-6 * L / nu;
    std::size_t q = 0;
    for (int i = 0; i <= nu; ++i) {
      const double su = sb + L * i / nu;
      while (q < kinks.size() && kinks[q] <= su + snap) {
        if (!sg.empty() && sg.back() >= kinks[q] - snap && sg.back() < kinks[q]) sg.pop_back();
        sg.push_back(kinks[q]);
        sg.push_back(kinks[q]);
        ++q;
      }
      if (!sg.empty() && std::abs(su - sg.back()) <= snap) continue;  // a kink took this node
      sg.push_back(su);
    }
  }
  const int n = static_cast<int>(sg.size()) - 1;
  fp.s.resize(static_cast<std::size_t>(n) + 1);
  fp.P.resize(static_cast<std::size_t>(n) + 1);
  const TendonLosses& lo = t.losses;
  auto node = [&](int i) { return sg[static_cast<std::size_t>(i)]; };

  // Cumulated angular deviation theta(s): the angle swept by the unit tangent
  // (1, e_y', e_z') / |.| of the 3-D trace, summed between closely spaced
  // tangents (16 per grid step). Exact for steep traces, where |e''| ds
  // overstates the angle by (1 + e'^2); a slope break is the zero-length
  // step between its two nodes (the left limit — a joint belongs to the
  // segment before it — then the right limit).
  std::vector<double> theta(static_cast<std::size_t>(n) + 1, 0.0);
  {
    // sampled strictly inside the anchors: the profile segment just outside
    // an anchor (a flat extension to the chain end) carries no force
    const double pad = 1e-9 * std::max(1.0, L);
    auto tangent = [&](double s) {
      s = std::clamp(s, sb + pad, se - pad);
      Vec3 v(1.0, tendon_dey(t, s), tendon_dez(t, s));
      return Vec3(v / v.norm());
    };
    constexpr int kSub = 16;
    Vec3 prev = tangent(sb);
    for (int i = 1; i <= n; ++i) {
      const double a = node(i - 1), b = node(i);
      double acc = 0;
      if (b - a <= 0) {  // a slope break: its right limit
        const Vec3 cur = tangent(b + pad);
        acc = std::atan2(prev.cross(cur).norm(), prev.dot(cur));
        prev = cur;
      } else {
        for (int k = 1; k <= kSub; ++k) {
          const Vec3 cur = tangent(a + (b - a) * k / kSub);
          acc += std::atan2(prev.cross(cur).norm(), prev.dot(cur));
          prev = cur;
        }
      }
      theta[static_cast<std::size_t>(i)] = theta[static_cast<std::size_t>(i - 1)] + acc;
    }
  }
  // Curvilinear length of the cable from the origin anchor, arc[i] = int
  // |(1, e_y', e_z')| dx — the cable is longer than its bars wherever it
  // slopes (no small-slope assumption: the wobble and the draw-in run on it).
  std::vector<double> arc(static_cast<std::size_t>(n) + 1, 0.0);
  {
    const double pad = 1e-9 * std::max(1.0, L);
    auto speed = [&](double s) {
      s = std::clamp(s, sb + pad, se - pad);
      const double dy = tendon_dey(t, s), dz = tendon_dez(t, s);
      return std::sqrt(1.0 + dy * dy + dz * dz);
    };
    // 3-point Gauss on 8 sub-steps: interior points only, so a step never
    // reads the far side of a slope break at its ends
    constexpr int kSub = 8;
    const double g = std::sqrt(0.6);
    for (int i = 1; i <= n; ++i) {
      const double a = node(i - 1), h = (node(i) - a) / kSub;
      double acc = 0;
      for (int k = 0; k < kSub && h > 0; ++k) {
        const double c = a + (k + 0.5) * h;
        acc += h / 18.0 * (5.0 * speed(c - g * h / 2) + 8.0 * speed(c) + 5.0 * speed(c + g * h / 2));
      }
      arc[static_cast<std::size_t>(i)] = arc[static_cast<std::size_t>(i - 1)] + acc;
    }
  }
  auto arc_at = [&](double sx) {  // the curvilinear length up to bar abscissa sx
    if (sx <= sg.front()) return arc.front();
    if (sx >= sg.back()) return arc.back();
    const std::size_t j = static_cast<std::size_t>(std::upper_bound(sg.begin(), sg.end(), sx) - sg.begin()) - 1;
    const double h = sg[j + 1] - sg[j];
    return h > 0 ? arc[j] + (sx - sg[j]) / h * (arc[j + 1] - arc[j]) : arc[j];
  };
  // Length rubbing against the concrete from the origin anchor: the anchored
  // length minus the free runs (external tendons), the wobble's lever.
  std::vector<double> wob(static_cast<std::size_t>(n) + 1, 0.0);
  for (int i = 0; i <= n; ++i) {
    const double si = node(i);
    double free_len = 0;
    for (const auto& [f0, f1] : t.free_runs) {
      const double a = std::max(f0, sb), b = std::min(f1, si);
      if (b > a) free_len += arc_at(b) - arc_at(a);
    }
    wob[static_cast<std::size_t>(i)] = arc[static_cast<std::size_t>(i)] - free_len;
  }
  auto friction_from = [&](double P0, bool from_end) {
    std::vector<double> P(static_cast<std::size_t>(n) + 1);
    for (int i = 0; i <= n; ++i) {
      double th = from_end ? theta[static_cast<std::size_t>(n)] - theta[static_cast<std::size_t>(i)]
                           : theta[static_cast<std::size_t>(i)];
      double x = from_end ? wob[static_cast<std::size_t>(n)] - wob[static_cast<std::size_t>(i)]
                          : wob[static_cast<std::size_t>(i)];
      P[static_cast<std::size_t>(i)] = P0 * std::exp(-lo.mu * (th + lo.k * x) - lo.phi_abs * x);
    }
    return P;
  };
  const std::vector<double> P1 = friction_from(t.jack_force, false);  // from the origin
  const std::vector<double> P2 = friction_from(t.jack_force, true);   // from the end
  std::vector<double> P = t.jack_both_ends ? P1 : (t.jack_from_end ? P2 : P1);
  if (t.jack_both_ends)
    for (int i = 0; i <= n; ++i)
      P[static_cast<std::size_t>(i)] = std::max(P1[static_cast<std::size_t>(i)], P2[static_cast<std::size_t>(i)]);

  // Anchorage draw-in (ST1 annex B1.1.2): over the zone [anchor, s1] the force
  // becomes sigma(s1) e^{-f alpha(s1 - s) - phi (s1 - s)}, the friction curve
  // reflected about s1 — with P_side the friction curve from the released
  // anchor, P'(s) = min(P(s), C / P_side(s)), C set so that exactly
  // slip * Ep * Ap is lost (a draw-in reaching the far end lowers the whole
  // reflected curve). Double-end jacking: the two anchors are released ONE
  // AFTER THE OTHER (TENSION OR EX: origin first; EX OR: end first), so the
  // overlapping draw-ins of a short tendon depend on the release order.
  if (lo.slip > 0 && lo.Ap > 0) {
    const double target = lo.slip * lo.Ep * lo.Ap;
    auto release = [&](const std::vector<double>& Ps) {
      // lost area for a level C: the positive part of P - C / Ps, its sign
      // change located linearly within an interval
      auto lost = [&](double C) {
        double sum = 0;
        auto d = [&](int i) {
          const std::size_t u = static_cast<std::size_t>(i);
          return P[u] - C / Ps[u];
        };
        for (int i = 1; i <= n; ++i) {
          // the curvilinear step of the interval
          const double h = arc[static_cast<std::size_t>(i)] - arc[static_cast<std::size_t>(i - 1)];
          const double a = d(i - 1), b = d(i);
          if (a >= 0 && b >= 0) sum += 0.5 * (a + b) * h;
          else if (a > 0) sum += 0.5 * a * h * a / (a - b);
          else if (b > 0) sum += 0.5 * b * h * b / (b - a);
        }
        return sum;
      };
      if (lost(0.0) <= target) {  // the whole tendon slips out: no force left
        std::fill(P.begin(), P.end(), 0.0);
        return;
      }
      double clo = 0, chi = 0;
      for (int i = 0; i <= n; ++i)
        chi = std::max(chi, P[static_cast<std::size_t>(i)] * Ps[static_cast<std::size_t>(i)]);
      for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (clo + chi);
        (lost(mid) > target ? clo : chi) = mid;
      }
      const double C = 0.5 * (clo + chi);
      for (int i = 0; i <= n; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        P[u] = std::min(P[u], C / Ps[u]);
      }
    };
    if (!t.jack_both_ends) {
      release(t.jack_from_end ? P2 : P1);
    } else if (t.release_end_first) {
      release(P2);
      release(P1);
    } else {
      release(P1);
      release(P2);
    }
  }

  for (int i = 0; i <= n; ++i) {
    fp.s[static_cast<std::size_t>(i)] = node(i);
    fp.P[static_cast<std::size_t>(i)] = P[static_cast<std::size_t>(i)];
  }
  // No force outside the anchors (ForceProfile::at holds the end values).
  const double eps = 1e-9 * std::max(1.0, Lchain);
  if (sb > eps) {
    fp.s.insert(fp.s.begin(), {0.0, sb - eps});
    fp.P.insert(fp.P.begin(), {0.0, 0.0});
  }
  if (se < Lchain - eps) {
    fp.s.insert(fp.s.end(), {se + eps, Lchain});
    fp.P.insert(fp.P.end(), {0.0, 0.0});
  }
  return fp;
}

std::vector<Load> tendon_loads(const Model& model, const Tendon& t, double fraction, int n_sub) {
  return tendon_loads(model, t, fraction, n_sub, -std::numeric_limits<double>::infinity(),
                      std::numeric_limits<double>::infinity());
}

std::vector<Load> tendon_loads(const Model& model, const Tendon& t, double fraction, int n_sub,
                               double s_lo, double s_hi) {
  require(fraction >= -1.0 && fraction <= 1.0 + 1e-12, "tendon: fraction out of range");
  ForceProfile fp = tendon_force_profile(model, t);
  std::vector<Load> out;
  double s0 = 0;
  double Lchain = tendon_chain_length(model, t);
  const double sb = std::clamp(t.s_begin, 0.0, Lchain);
  const double se = t.s_end < 0 ? Lchain : std::clamp(t.s_end, sb, Lchain);
  const std::vector<double> kinks = trace_kinks(t, sb, se);
  for (std::size_t bi = 0; bi < t.bars.size(); ++bi) {
    const Bar& b = model.bar(t.bars[bi]);
    double Lb = model.bar_chord_length(b);
    // A flipped bar is run from nj to ni: the chain abscissa at the bar-local
    // x is s0 + Lb - x. A cable's equivalent loads do not depend on the travel
    // direction, so every term below is written in the bar-local x with the
    // profile re-expressed in it (the eccentricity stays in the bar's frame).
    const bool flip = !t.flips.empty() && t.flips[bi];
    // The tendon force is P along the cable's tangent (1, e_y', e_z') / |…|:
    // with exact_angle the bar takes its AXIAL component P cos(alpha) (as
    // ST1: N = P cos alpha, 2.5 % below P at a 0.23 slope) and P cos(alpha)
    // e' across.
    // nil force all along this bar: zero loads, nothing evaluated
    const bool idle = s0 + Lb <= s_lo || s0 >= s_hi;
    auto Pax = [&](double s) {
      if (idle) return 0.0;
      if (!t.exact_angle) return fraction * fp.at(s);
      const double dy = tendon_dey(t, s), dz = tendon_dez(t, s);
      return fraction * fp.at(s) / std::sqrt(1.0 + dy * dy + dz * dz);
    };
    auto Pe = [&](double x, int comp) {  // (P_axial * e_comp) at bar-local abscissa x
      if (idle) return 0.0;
      double s = flip ? s0 + Lb - x : s0 + x;
      double P = Pax(s);
      return comp == 0 ? P : P * (comp == 1 ? tendon_ey(t, s) : tendon_ez(t, s));
    };
    // The occupied part [xa0, xb0] of the bar (set below); every quantity is
    // read inside it, so an interior anchor (where the force drops to zero)
    // never enters a derivative.
    double xa0 = 0, xb0 = Lb;
    // P e_comp' (d/dx): the transverse component of the tendon force — the
    // deviation the concrete takes as a force. The rest of d(P e)/dx, P' e,
    // is the bond force acting at the eccentricity: a distributed COUPLE, no
    // shear (ST1: Ty = 0 along a pretension transmission zone, Exemple5).
    auto PdE = [&](double x, int comp) {
      if (idle) return 0.0;
      const double sgn = flip ? -1.0 : 1.0;
      double s = flip ? s0 + Lb - x : s0 + x;
      double P = Pax(s);
      return sgn * P * (comp == 1 ? tendon_dey(t, s) : tendon_dez(t, s));
    };
    // Anchor closures ONLY at the true chain ends: at internal boundaries the
    // adjacent bars' terms cancel exactly in the load vector, and emitting
    // them would pollute the per-bar span records station recovery integrates
    // against (internal forces flow through the nodes there).
    // The chain ends sit at local x = 0 or x = Lb depending on the flip.
    // The anchors sit at chain abscissae sb / se (the chain ends by default):
    // in bar-local x, the part of this bar the tendon occupies is [xa0, xb0],
    // and an anchor closure is emitted where that part ends at an anchor.
    const double tol = 1e-9 * std::max(1.0, Lchain);
    auto local_x = [&](double sc) { return flip ? s0 + Lb - sc : sc - s0; };
    const bool has_sb = sb >= s0 - tol && sb <= s0 + Lb + tol;
    const bool has_se = se >= s0 - tol && se <= s0 + Lb + tol;
    if (se < s0 - tol || sb > s0 + Lb + tol) {  // the tendon does not reach this bar
      s0 += Lb;
      continue;
    }
    if (has_sb) (flip ? xb0 : xa0) = std::clamp(local_x(sb), 0.0, Lb);
    if (has_se) (flip ? xa0 : xb0) = std::clamp(local_x(se), 0.0, Lb);
    if (xb0 - xa0 <= tol) {  // an anchor on this bar's end node: the neighbour carries it
      s0 += Lb;
      continue;
    }
    // Anchor at the low local end of the occupied part (start anchor of an
    // unflipped bar, end anchor of a flipped one) and at the high end.
    const bool anchor_lo = flip ? has_se : has_sb;
    const bool anchor_hi = flip ? has_sb : has_se;
    // G = e_y F_z - e_z F_y: the torque about the bar axis of the tendon's
    // transverse force acting at the eccentricity (a cable offset sideways
    // twists the bar where it deviates vertically, as in ST1). Its density is
    // dG/dx exactly (e_y' F_z - e_z' F_y = P (e_y' e_z' - e_z' e_y') = 0).
    auto G = [&](double x) {
      if (idle) return 0.0;
      const double s = flip ? s0 + Lb - x : s0 + x;
      return tendon_ey(t, s) * PdE(x, 2) - tendon_ez(t, s) * PdE(x, 1);
    };
    if (anchor_lo) {
      const double x = xa0, xi = std::min(x + 1e-9 * Lb, xb0);
      Vec3 Fi(Pe(xi, 0), PdE(xi, 1), PdE(xi, 2));
      Vec3 Mi(G(xi), +Pe(xi, 2), -Pe(xi, 1));
      out.push_back(Load::point(b.id, x, Fi, Mi, LoadFrame::Local));
    }
    if (anchor_hi) {
      const double x = xb0, xj = std::max(x - 1e-9 * Lb, xa0);
      Vec3 Fj = -Vec3(Pe(xj, 0), PdE(xj, 1), PdE(xj, 2));
      Vec3 Mj = -Vec3(G(xj), +Pe(xj, 2), -Pe(xj, 1));
      out.push_back(Load::point(b.id, x, Fj, Mj, LoadFrame::Local));
    }
    // Distributed deviation + drag on n_sub spans of the occupied part. The
    // load density is the derivative of F_c = P (axial drag P') and
    // F_c = (P e_c)' (deviation (P e_c)''); on each span a LINEAR load is
    // fitted to the EXACT resultant F_c(xb) - F_c(xa) and first moment
    // h F_c(xb) - int F_c — so the tendon's loads are self-equilibrated to
    // round-off whatever the trace's kinks (sampling (P e)'' at the span ends
    // left MN-level vertical residuals on tight ALIGNE kinks, Fontenoy).
    // F_c is read at the same clamped points as the anchor closures above.
    const double Lo = xb0 - xa0;
    const double dclamp = 1e-9 * Lb;
    // Pulled inside only at an anchor (where the force jumps to zero); at an
    // internal bar end both bars read the SAME node value, so the per-bar
    // resultants telescope exactly along the chain.
    const double lo_pull = anchor_lo ? dclamp : 0.0, hi_pull = anchor_hi ? dclamp : 0.0;
    auto cl = [&](double x) { return std::clamp(x, xa0 + lo_pull, xb0 - hi_pull); };
    auto F = [&](double x, int c) { return c == 0 ? Pe(cl(x), 0) : PdE(cl(x), c); };
    // composite Simpson of g over [a, b] (piecewise-smooth integrands)
    auto simpson = [&](double a, double bb, auto g) {
      const int m = 8;
      double acc = 0;
      for (int j = 0; j < m; ++j) {
        double u = a + (bb - a) * j / m, v = a + (bb - a) * (j + 1) / m;
        acc += (v - u) / 6.0 * (g(u) + 4 * g(0.5 * (u + v)) + g(v));
      }
      return acc;
    };
    // linear density over [xa, xb] with resultant R and first moment M1 about xa
    auto fit = [](double h, double R, double M1, double& wa_, double& wb_) {
      wb_ = 6 * M1 / (h * h) - 2 * R / h;
      wa_ = 2 * R / h - wb_;
    };
    // Slope breaks of the trace on this bar (bar-local x). At a break the
    // deviation is concentrated: the point force F(x+) - F(x-) (P times the
    // turn of the tangent, the friction drop, the change of P cos alpha) and
    // the couple of the force jump at the eccentricity; the spans on either
    // side read the one-sided limits. A break on the bar's end node belongs
    // to the bar ending there (chain-wise); the next bar's first span reads
    // its right limit.
    const double dk = 1e-9 * std::max(1.0, Lb);
    std::vector<std::pair<double, bool>> cuts;  // span boundaries, is a break
    for (int k = 0; k <= n_sub; ++k) cuts.push_back({xa0 + Lo * k / n_sub, false});
    for (double sk : kinks) {
      if (sk < s0 - tol || sk > s0 + Lb + tol) continue;
      const double xk = std::clamp(local_x(sk), 0.0, Lb);
      if (xk < xa0 - tol || xk > xb0 + tol) continue;
      if (sk > s0 + tol) {  // emitted by this bar
        Vec3 Fk, Mk;
        for (int c = 0; c < 3; ++c)
          Fk(c) = c == 0 ? Pe(xk + dk, 0) - Pe(xk - dk, 0) : PdE(xk + dk, c) - PdE(xk - dk, c);
        Mk = Vec3(G(xk + dk) - G(xk - dk), Pe(xk + dk, 2) - Pe(xk - dk, 2),
                  -(Pe(xk + dk, 1) - Pe(xk - dk, 1)));
        out.push_back(Load::point(b.id, xk, Fk, Mk, LoadFrame::Local));
      }
      bool merged = false;
      for (auto& c : cuts)
        if (std::abs(c.first - xk) <= 1e-6 * Lo / n_sub) {
          c = {xk, true};
          merged = true;
        }
      if (!merged) cuts.push_back({xk, true});
    }
    std::sort(cuts.begin(), cuts.end());
    for (std::size_t k = 0; k + 1 < cuts.size() && Lo > tol; ++k) {
      const double xa = cuts[k].first, xb = cuts[k + 1].first;
      const double h = xb - xa;
      if (h <= tol) continue;
      // one-sided at a break: the span's own side of it
      const double ea = cuts[k].second ? dk : 0.0, eb = cuts[k + 1].second ? dk : 0.0;
      auto Fs = [&](double x, int c) {
        if (x <= xa) x = xa + ea;
        if (x >= xb) x = xb - eb;
        return F(x, c);
      };
      auto Pes = [&](double x, int c) {
        if (x <= xa) x = xa + ea;
        if (x >= xb) x = xb - eb;
        return Pe(cl(x), c);
      };
      auto Gs = [&](double x) {
        if (x <= xa) x = xa + ea;
        if (x >= xb) x = xb - eb;
        return G(cl(x));
      };
      Vec3 wa, wb, ma = Vec3::Zero(), mb = Vec3::Zero();
      for (int c = 0; c < 3; ++c) {
        // force density dF_c/dx: resultant F(xb) - F(xa), first moment
        // h F(xb) - int F
        double R = Fs(xb, c) - Fs(xa, c);
        double I = simpson(xa, xb, [&](double x) { return Fs(x, c); });
        fit(h, R, h * Fs(xb, c) - I, wa(c), wb(c));
        if (c == 0) continue;
        // couple density D = P' e (= d(P e)/dx - P e'), exact by parts:
        // int D = [P e] - int P e', int (x - xa) D = h (P e)(xb) - int P e
        // - int (x - xa) P e'
        const double ipe = simpson(xa, xb, [&](double x) { return Pes(x, c); });
        const double R_d = Pes(xb, c) - Pes(xa, c) - I;
        const double M_d = h * Pes(xb, c) - ipe -
                           simpson(xa, xb, [&](double x) { return (x - xa) * Fs(x, c); });
        double da, db;
        fit(h, R_d, M_d, da, db);
        // about local z the couple is -D_y, about local y +D_z (anchors: M =
        // (0, +P e_z, -P e_y))
        if (c == 1) {
          ma.z() = -da;
          mb.z() = -db;
        } else {
          ma.y() = da;
          mb.y() = db;
        }
      }
      {  // torque density dG/dx, same exact resultant / first-moment fit
        auto Gc = Gs;
        const double R = Gc(xb) - Gc(xa);
        const double I = simpson(xa, xb, Gc);
        fit(h, R, h * Gc(xb) - I, ma.x(), mb.x());
      }
      Load l = Load::linear(b.id, wa, wb, LoadFrame::Local, xa, xb);
      l.M = ma;
      l.M2 = mb;
      out.push_back(l);
    }
    s0 += Lb;
  }
  return out;
}

}  // namespace keta

#include "keta/influence.h"

#include <algorithm>
#include <memory>
#include <unordered_map>

namespace keta {

ResponseKind ResponseKind::displacement(Id node, int comp) {
  ResponseKind r;
  r.what = What::Displacement;
  r.node = node;
  r.comp = comp;
  return r;
}
ResponseKind ResponseKind::reaction(Id node, int comp) {
  ResponseKind r;
  r.what = What::Reaction;
  r.node = node;
  r.comp = comp;
  return r;
}
ResponseKind ResponseKind::bar_force(Id bar, double x, int comp) {
  ResponseKind r;
  r.what = What::BarForce;
  r.bar = bar;
  r.x = x;
  r.comp = comp;
  return r;
}

ResponseKind ResponseKind::bar_displacement(Id bar, double x, int comp) {
  ResponseKind r;
  r.what = What::BarDisplacement;
  r.bar = bar;
  r.x = x;
  r.comp = comp;
  return r;
}

LoadPath build_path(const Model& model, const std::vector<Id>& bars, int n_per_bar) {
  require(!bars.empty() && n_per_bar >= 2, "path: need bars and >= 2 stations per bar");
  LoadPath p;
  double s0 = 0;
  Id chain_end = 0;
  for (std::size_t bi = 0; bi < bars.size(); ++bi) {
    const Bar& b = model.bar(bars[bi]);
    double L = model.bar_chord_length(b);
    bool flip = false;
    if (bi == 0) {
      // orientation fixed by the second bar if any
      if (bars.size() > 1) {
        const Bar& b2 = model.bar(bars[1]);
        if (b.nj == b2.ni || b.nj == b2.nj) {
          flip = false;
          chain_end = b.nj;
        } else if (b.ni == b2.ni || b.ni == b2.nj) {
          flip = true;
          chain_end = b.ni;
        } else {
          fail("path: bars " + std::to_string(b.id) + " and " + std::to_string(b2.id) +
               " are not contiguous");
        }
      } else {
        chain_end = b.nj;
      }
    } else {
      if (b.ni == chain_end) {
        flip = false;
        chain_end = b.nj;
      } else if (b.nj == chain_end) {
        flip = true;
        chain_end = b.ni;
      } else {
        fail("path: bar " + std::to_string(b.id) + " does not continue the chain");
      }
    }
    for (int k = 0; k < n_per_bar; ++k) {
      if (bi > 0 && k == 0) continue;  // shared node: keep one station
      double t = static_cast<double>(k) / (n_per_bar - 1);
      double xl = flip ? L * (1 - t) : L * t;
      p.stations.push_back({b.id, xl, s0 + L * t});
    }
    s0 += L;
  }
  p.length = s0;
  return p;
}

double InfluenceLine::on(std::size_t i, double x) const {
  const double h = s[i + 1] - s[i];
  if (h <= 0) return eta[i];
  const double t = (x - s[i]) / h;
  if (!cubic()) return eta[i] + t * (eta[i + 1] - eta[i]);
  const double t2 = t * t, t3 = t2 * t;
  return (2 * t3 - 3 * t2 + 1) * eta[i] + (t3 - 2 * t2 + t) * h * dr[i] +
         (-2 * t3 + 3 * t2) * eta[i + 1] + (t3 - t2) * h * dl[i + 1];
}

namespace {
// exact integral of the line over [x0, x1] inside interval i (2-point Gauss:
// exact for the cubic)
double piece_integral(const InfluenceLine& il, std::size_t i, double x0, double x1) {
  if (x1 <= x0) return 0;
  if (!il.cubic()) return 0.5 * (il.on(i, x0) + il.on(i, x1)) * (x1 - x0);
  const double m = 0.5 * (x0 + x1), r = 0.5 * (x1 - x0) / std::sqrt(3.0);
  return 0.5 * (x1 - x0) * (il.on(i, m - r) + il.on(i, m + r));
}
// the real zeros of the line strictly inside interval i, ascending
std::vector<double> piece_zeros(const InfluenceLine& il, std::size_t i) {
  std::vector<double> z;
  const double a = il.s[i], b = il.s[i + 1];
  if (b <= a) return z;
  if (!il.cubic()) {
    const double ea = il.eta[i], eb = il.eta[i + 1];
    if ((ea > 0 && eb < 0) || (ea < 0 && eb > 0)) z.push_back(a + ea / (ea - eb) * (b - a));
    return z;
  }
  // split at the extrema of the cubic (roots of its derivative), then bisect
  // each monotone piece that changes sign
  const double h = b - a, e0 = il.eta[i], e1 = il.eta[i + 1], m0 = il.dr[i] * h, m1 = il.dl[i + 1] * h;
  // p(t) = c3 t^3 + c2 t^2 + c1 t + c0 ; p'(t) = 3 c3 t^2 + 2 c2 t + c1
  const double c3 = 2 * e0 + m0 - 2 * e1 + m1, c2 = -3 * e0 - 2 * m0 + 3 * e1 - m1, c1 = m0;
  std::vector<double> ts{0.0};
  const double A = 3 * c3, B = 2 * c2, C = c1;
  if (std::fabs(A) > 1e-300) {
    const double disc = B * B - 4 * A * C;
    if (disc > 0) {
      const double q = -0.5 * (B + std::copysign(std::sqrt(disc), B));
      for (double t : {q / A, q != 0 ? C / q : 2.0})
        if (t > 0 && t < 1) ts.push_back(t);
    }
  } else if (std::fabs(B) > 1e-300) {
    const double t = -C / B;
    if (t > 0 && t < 1) ts.push_back(t);
  }
  ts.push_back(1.0);
  std::sort(ts.begin(), ts.end());
  auto p = [&](double t) { return il.on(i, a + t * h); };
  for (std::size_t k = 1; k < ts.size(); ++k) {
    double lo = ts[k - 1], hi = ts[k], plo = p(lo), phi = p(hi);
    if (!((plo > 0 && phi < 0) || (plo < 0 && phi > 0))) continue;
    for (int it = 0; it < 100 && hi - lo > 1e-15; ++it) {
      const double mid = 0.5 * (lo + hi), pm = p(mid);
      if ((pm > 0) == (plo > 0)) {
        lo = mid;
        plo = pm;
      } else {
        hi = mid;
      }
    }
    z.push_back(a + 0.5 * (lo + hi) * h);
  }
  return z;
}
}  // namespace

double InfluenceLine::at(double x) const {
  if (s.empty()) return 0;
  if (x <= s.front()) return x < s.front() ? 0.0 : eta.front();
  if (x >= s.back()) return x > s.back() ? 0.0 : eta.back();
  // binary search, left-most interval wins on duplicates
  std::size_t lo = 0, hi = s.size() - 1;
  while (hi - lo > 1) {
    std::size_t mid = (lo + hi) / 2;
    if (s[mid] <= x)
      lo = mid;
    else
      hi = mid;
  }
  (void)hi;
  return on(lo, x);
}

double InfluenceLine::integral(double a, double b) const {
  if (s.size() < 2 || b <= a) return 0;
  a = std::max(a, s.front());
  b = std::min(b, s.back());
  if (b <= a) return 0;
  // Interpolate inside interval i itself (not through at(), which picks the
  // interval AFTER a duplicate abscissa and so read a shear jump's far side).
  double sum = 0;
  for (std::size_t i = 1; i < s.size(); ++i) {
    double x0 = std::max(a, s[i - 1]), x1 = std::min(b, s[i]);
    if (x1 <= x0) continue;
    sum += piece_integral(*this, i - 1, x0, x1);
  }
  return sum;
}

IlPrefix::IlPrefix(const InfluenceLine& il) : il_(il) {
  cum_.assign(il.s.size(), 0.0);
  for (std::size_t i = 1; i < il.s.size(); ++i)
    cum_[i] = cum_[i - 1] + piece_integral(il, i - 1, il.s[i - 1], il.s[i]);
}

double IlPrefix::prefix_at(double q) const {
  const auto& s = il_.s;
  if (s.size() < 2 || q <= s.front()) return 0;
  if (q >= s.back()) return cum_.back();
  // last station at or before q (the same interval InfluenceLine::at picks);
  // a zero-length interval at a duplicate abscissa adds nothing to cum_
  std::size_t lo = static_cast<std::size_t>(std::upper_bound(s.begin(), s.end(), q) - s.begin()) - 1;
  return cum_[lo] + piece_integral(il_, lo, s[lo], q);
}

double IlPrefix::prefix_in(std::size_t lo, double q) const {
  return cum_[lo] + piece_integral(il_, lo, il_.s[lo], q);
}

double IlPrefix::integral(double a, double b) const {
  if (b <= a) return 0;
  return prefix_at(b) - prefix_at(a);
}

namespace {
double signed_integral(const InfluenceLine& il, bool positive) {
  double sum = 0;
  for (std::size_t i = 1; i < il.s.size(); ++i) {
    if (il.s[i] - il.s[i - 1] <= 0) continue;
    std::vector<double> cuts{il.s[i - 1]};
    for (double z : piece_zeros(il, i - 1)) cuts.push_back(z);
    cuts.push_back(il.s[i]);
    for (std::size_t k = 1; k < cuts.size(); ++k) {
      const double v = piece_integral(il, i - 1, cuts[k - 1], cuts[k]);
      if (positive && v > 0) sum += v;
      if (!positive && v < 0) sum += v;
    }
  }
  return sum;
}
}  // namespace

double InfluenceLine::integral_pos() const { return signed_integral(*this, true); }
double InfluenceLine::integral_neg() const { return signed_integral(*this, false); }

std::pair<double, double> InfluenceLine::point_extreme() const {
  double mx = 0, mn = 0;
  for (double e : eta) {
    mx = std::max(mx, e);
    mn = std::min(mn, e);
  }
  if (cubic())  // the extrema of each cubic piece (roots of its derivative)
    for (std::size_t i = 0; i + 1 < s.size(); ++i) {
      const double h = s[i + 1] - s[i];
      if (h <= 0) continue;
      const double e0 = eta[i], e1 = eta[i + 1], m0 = dr[i] * h, m1 = dl[i + 1] * h;
      const double A = 3 * (2 * e0 + m0 - 2 * e1 + m1), B = 2 * (-3 * e0 - 2 * m0 + 3 * e1 - m1), C = m0;
      double ts[2];
      int nt = 0;
      if (std::fabs(A) > 1e-300) {
        const double disc = B * B - 4 * A * C;
        if (disc >= 0) {
          ts[nt++] = (-B + std::sqrt(disc)) / (2 * A);
          ts[nt++] = (-B - std::sqrt(disc)) / (2 * A);
        }
      } else if (std::fabs(B) > 1e-300) {
        ts[nt++] = -C / B;
      }
      for (int q = 0; q < nt; ++q)
        if (const double t = ts[q]; t > 0 && t < 1) {
          const double v = on(i, s[i] + t * h);
          mx = std::max(mx, v);
          mn = std::min(mn, v);
        }
    }
  return {mx, mn};
}

double extract_response(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                        const VecX& u, const ResponseKind& response) {
  double value = 0;
  switch (response.what) {
    case ResponseKind::What::Displacement: {
      int ni = model.node_index(response.node);
      int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, response.comp))];
      value = fi >= 0 ? u(fi) : 0.0;
      break;
    }
    case ResponseKind::What::Reaction: {
      // R = sum over incident elements of (K_e u)_dof - F_applied_dof.
      require(asm_.continua.empty(), "influence: a reaction response does not take continuum elements");
      int ni = model.node_index(response.node);
      int dof = asm_.dofs.full(ni, response.comp);
      double fint = 0;
      for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
        const Bar& b = model.bars[e];
        if (b.kind == BarKind::Cable) continue;
        int na = model.node_index(b.ni), nb = model.node_index(b.nj);
        if (na != ni && nb != ni) continue;
        Vec12 ue;
        for (int c = 0; c < 6; ++c) {
          int fa = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(na, c))];
          int fb = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nb, c))];
          ue(c) = fa >= 0 ? u(fa) : 0.0;
          ue(6 + c) = fb >= 0 ? u(fb) : 0.0;
        }
        // bar_ke_assembled: the CARA VAR tapered block when the bar has one.
        Vec12 fe = bar_ke_assembled(model, asm_.elements[e]) * ue;
        fint += (na == ni) ? fe(response.comp) : fe(6 + response.comp);
      }
      // Facet shells incident to the node contribute their internal forces
      // f_int = K_e u_e exactly like the bars — WITHOUT this term a reaction
      // influence on a shell-bearing model reads silently wrong ordinates
      // (the historical reason shells were gated out of the moving-load
      // engine).
      for (const ShellElem& se : asm_.shells) {
        int local = -1;
        for (int k = 0; k < 3; ++k)
          if (se.nodes[static_cast<std::size_t>(k)] == ni) local = k;
        if (local < 0) continue;
        Vec18 ue;
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 6; ++c) {
            int fi = asm_.dofs.free_ix[static_cast<std::size_t>(
                asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c))];
            ue(6 * k + c) = fi >= 0 ? u(fi) : 0.0;
          }
        Vec18 fe = se.el.ke() * ue;
        fint += fe(6 * local + response.comp);
      }
      // Elastic links incident to the node act on it like the elements:
      // without them, a support reached through a link read a zero reaction.
      for (const LinkElem& le : asm_.links) {
        if (le.ni != ni && le.nj != ni) continue;
        Vec12 ue;
        for (int c = 0; c < 6; ++c) {
          int fa = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(le.ni, c))];
          int fb = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(le.nj, c))];
          ue(c) = fa >= 0 ? u(fa) : 0.0;
          ue(6 + c) = fb >= 0 ? u(fb) : 0.0;
        }
        const Vec12 fe = le.K * ue;
        if (le.ni == ni) fint += fe(response.comp);
        if (le.nj == ni) fint += fe(6 + response.comp);
      }
      // R = K u - F restricted to the incident rows; an elastic spring at a
      // free DOF surfaces automatically as R = -k u through the equilibrium.
      value = fint - cl.nodal_full(dof);
      break;
    }
    case ResponseKind::What::BarForce: {
      std::size_t e = static_cast<std::size_t>(model.bar_index(response.bar));
      const BarElement& el = asm_.elements[e];
      const Bar& b = model.bars[e];
      int na = model.node_index(b.ni), nb = model.node_index(b.nj);
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        int fa = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(na, c))];
        int fb = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nb, c))];
        ue(c) = fa >= 0 ? u(fa) : 0.0;
        ue(6 + c) = fb >= 0 ? u(fb) : 0.0;
      }
      const Vec12& f0 = cl.f0[e];
      Vec12 ul = el.local_displacements(ue, f0);
      const TaperedBar* tb = model.tapered_of(b.id);
      Vec12 fe = tb ? Vec12(tb->ke * ul - f0) : el.end_forces_local(ul, f0);
      if (tb && tb->recover) {  // custom block (bef.h): its exact stations
        double xr = std::min(std::max(response.x, 0.0), el.L());
        BarResult br;
        recover_bar_stations(model, el, cl.span[e], fe, ul, {xr}, br);
        const Station& st = br.stations[0];
        const double comps[6] = {st.N, st.Vy, st.Vz, st.T, st.My, st.Mz};
        require(response.comp >= 0 && response.comp <= 5, "influence: bad component");
        value = comps[response.comp];
        break;
      }
      // single-station recovery at response.x with the unit load's span record
      double x = response.x;
      const BarSpanLoads& sp = cl.span[e];
      auto force_int = [&](int c) {
        double q = 0;
        for (const PointAt& pl : sp.points)
          if (pl.a <= x) q += pl.F(c);
        return q;
      };
      auto point_moment = [&](int c) {
        double m = 0;
        for (const PointAt& pl : sp.points)
          if (pl.a <= x) m += pl.M(c);
        return m;
      };
      auto moment_of = [&](int c) {
        double m = 0;
        for (const PointAt& pl : sp.points)
          if (pl.a <= x) m += pl.F(c) * (x - pl.a);
        return m;
      };
      switch (response.comp) {
        case 0: value = -(fe(0) + force_int(0)); break;
        case 1: value = -(fe(1) + force_int(1)); break;
        case 2: value = -(fe(2) + force_int(2)); break;
        case 3: value = -(fe(3) + point_moment(0)); break;
        case 4: value = -fe(4) - fe(2) * x - moment_of(2) - point_moment(1); break;
        case 5: value = -fe(5) + fe(1) * x + moment_of(1) - point_moment(2); break;
        default: fail("influence: bad component");
      }
      break;
    }
    case ResponseKind::What::BarDisplacement: {
      require(response.comp >= 0 && response.comp <= 5, "influence: bad component");
      std::size_t e = static_cast<std::size_t>(model.bar_index(response.bar));
      const BarElement& el = asm_.elements[e];
      const Bar& b = model.bars[e];
      int na = model.node_index(b.ni), nb = model.node_index(b.nj);
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        int fa = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(na, c))];
        int fb = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nb, c))];
        ue(c) = fa >= 0 ? u(fa) : 0.0;
        ue(6 + c) = fb >= 0 ? u(fb) : 0.0;
      }
      const Vec12& f0 = cl.f0[e];
      Vec12 ul = el.local_displacements(ue, f0);
      const TaperedBar* tb = model.tapered_of(b.id);
      Vec12 fe = tb ? Vec12(tb->ke * ul - f0) : el.end_forces_local(ul, f0);
      double L = el.L();
      double x = std::min(std::max(response.x, 0.0), L);
      // Recovery grid: the ends, the section, every concentrated load and
      // segment end of the bar (the Simpson curvature integration is exact
      // between kinks of the moment diagram) and a few interior points.
      const BarSpanLoads& sp = cl.span[e];
      std::vector<double> xs = {0.0, x, L};
      for (int k = 1; k < 8; ++k) xs.push_back(L * k / 8.0);
      for (const PointAt& pl : sp.points) xs.push_back(std::min(std::max(pl.a, 0.0), L));
      require_linear_span(sp.segs, "influence line");
      for (const DistSeg& sg : sp.segs) {
        xs.push_back(std::min(std::max(sg.a, 0.0), L));
        xs.push_back(std::min(std::max(sg.b, 0.0), L));
      }
      for (const KinkAt& kk : sp.kinks) xs.push_back(std::min(std::max(kk.a, 0.0), L));
      std::sort(xs.begin(), xs.end());
      xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
      BarResult br;
      recover_bar_stations(model, el, sp, fe, ul, xs, br);
      std::size_t ix = static_cast<std::size_t>(std::find(xs.begin(), xs.end(), x) - xs.begin());
      const Station& st = br.stations[ix];
      Vec3 g = el.R().transpose() * (response.comp < 3 ? st.defl : st.rot);
      value = g(response.comp % 3);
      break;
    }
  }
  require(std::isfinite(value), "influence: non-finite ordinate");
  return value;
}

void extract_bar_displacements(const Model& model, const Assembled& asm_, const CaseLoads& cl,
                               const VecX& u, Id bar, const std::vector<const ResponseKind*>& rs,
                               std::vector<double>& out) {
  out.assign(rs.size(), 0.0);
  if (rs.empty()) return;
  std::size_t e = static_cast<std::size_t>(model.bar_index(bar));
  const BarElement& el = asm_.elements[e];
  const Bar& b = model.bars[e];
  int na = model.node_index(b.ni), nb = model.node_index(b.nj);
  Vec12 ue;
  for (int c = 0; c < 6; ++c) {
    int fa = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(na, c))];
    int fb = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nb, c))];
    ue(c) = fa >= 0 ? u(fa) : 0.0;
    ue(6 + c) = fb >= 0 ? u(fb) : 0.0;
  }
  const Vec12& f0 = cl.f0[e];
  Vec12 ul = el.local_displacements(ue, f0);
  const TaperedBar* tb = model.tapered_of(b.id);
  Vec12 fe = tb ? Vec12(tb->ke * ul - f0) : el.end_forces_local(ul, f0);
  const double L = el.L();
  // the same recovery grid as extract_response, with every section asked
  const BarSpanLoads& sp = cl.span[e];
  std::vector<double> xs = {0.0, L};
  for (int k = 1; k < 8; ++k) xs.push_back(L * k / 8.0);
  for (const PointAt& pl : sp.points) xs.push_back(std::min(std::max(pl.a, 0.0), L));
  require_linear_span(sp.segs, "influence line");
  for (const DistSeg& sg : sp.segs) {
    xs.push_back(std::min(std::max(sg.a, 0.0), L));
    xs.push_back(std::min(std::max(sg.b, 0.0), L));
  }
  for (const KinkAt& kk : sp.kinks) xs.push_back(std::min(std::max(kk.a, 0.0), L));
  for (const ResponseKind* r : rs) xs.push_back(std::min(std::max(r->x, 0.0), L));
  std::sort(xs.begin(), xs.end());
  xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  BarResult br;
  recover_bar_stations(model, el, sp, fe, ul, xs, br);
  for (std::size_t q = 0; q < rs.size(); ++q) {
    const double x = std::min(std::max(rs[q]->x, 0.0), L);
    const std::size_t ix = static_cast<std::size_t>(std::find(xs.begin(), xs.end(), x) - xs.begin());
    const Station& st = br.stations[ix];
    Vec3 g = el.R().transpose() * (rs[q]->comp < 3 ? st.defl : st.rot);
    out[q] = g(rs[q]->comp % 3);
    require(std::isfinite(out[q]), "influence: non-finite ordinate");
  }
}

// Full DOF indices a response extractor reads (the support of its gradient).
std::vector<int> response_dofs(const Model& model, const Assembled& asm_, const ResponseKind& r) {
  const DofMap& d = asm_.dofs;
  std::vector<int> dofs;
  auto node_dofs = [&](int ni) {
    for (int c = 0; c < 6; ++c) dofs.push_back(d.full(ni, c));
  };
  switch (r.what) {
    case ResponseKind::What::Displacement:
      dofs.push_back(d.full(model.node_index(r.node), r.comp));
      break;
    case ResponseKind::What::BarForce:
    case ResponseKind::What::BarDisplacement: {
      const Bar& b = model.bars[static_cast<std::size_t>(model.bar_index(r.bar))];
      node_dofs(model.node_index(b.ni));
      node_dofs(model.node_index(b.nj));
      break;
    }
    case ResponseKind::What::Reaction: {
      int ni = model.node_index(r.node);
      for (const Bar& b : model.bars) {
        if (b.kind == BarKind::Cable) continue;
        int na = model.node_index(b.ni), nb = model.node_index(b.nj);
        if (na != ni && nb != ni) continue;
        node_dofs(na);
        node_dofs(nb);
      }
      for (const ShellElem& se : asm_.shells)
        if (std::find(se.nodes.begin(), se.nodes.end(), ni) != se.nodes.end())
          for (int k : se.nodes) node_dofs(k);
      for (const LinkElem& le : asm_.links)
        if (le.ni == ni || le.nj == ni) {
          node_dofs(le.ni);
          node_dofs(le.nj);
        }
      break;
    }
  }
  std::sort(dofs.begin(), dofs.end());
  dofs.erase(std::unique(dofs.begin(), dofs.end()), dofs.end());
  return dofs;
}

namespace {

// Full DOF indices one unit load writes (a bar load: its bar's 12; a nodal
// load: its node's 6).
void load_dofs(const Model& model, const Load& l, const DofMap& d, std::vector<int>& out) {
  out.clear();
  auto node_dofs = [&](Id node) {
    int ni = model.node_index(node);
    for (int c = 0; c < 6; ++c) out.push_back(d.full(ni, c));
  };
  if (l.kind == Load::Kind::NodalForce) {
    node_dofs(l.node);
  } else {
    const Bar& b = model.bars[static_cast<std::size_t>(model.bar_index(l.bar))];
    node_dofs(b.ni);
    node_dofs(b.nj);
  }
}

}  // namespace

std::vector<double> unit_load_sweep(const Model& model, const std::vector<Load>& unit,
                                    const ResponseKind& response, FactorCache& fac,
                                    const SolveHooks* hooks) {
  hooks_check(hooks);
  hooks_report(hooks, 0.0, "assembly");
  model.validate();
  Assembled asm_ = assemble_elements(model, true);
  // rigid links: the adjoint solve on the reduced space (w = T Kr^-1 T^t g);
  // a reaction ordinate at a master would miss its slaves' transferred
  // forces (extract_response reads the master's own elements): refused
  if (asm_.red.active && response.what == ResponseKind::What::Reaction)
    for (const RigidLink& l : model.rigid_links)
      require(l.master != response.node,
              "influence: the reaction of a rigid-link master is not available (its slaves' forces)");
  SpMat K = assemble_stiffness(model, asm_);
  hooks_check(hooks);
  hooks_report(hooks, 0.25, "factorization");
  fac.set_reduction(asm_.red.active ? asm_.red.T : SpMat());
  fac.factorize(K);
  const DofMap& d = asm_.dofs;

  CaseLoads cl = empty_case_loads(asm_);
  VecX probe = VecX::Zero(d.n_free);
  VecX g = VecX::Zero(d.n_free);
  for (int f : response_dofs(model, asm_, response)) {
    int fi = d.free_ix[static_cast<std::size_t>(f)];
    if (fi < 0) continue;
    probe(fi) = 1;
    g(fi) = extract_response(model, asm_, cl, probe, response);
    probe(fi) = 0;
  }
  hooks_check(hooks);
  hooks_report(hooks, 0.5, "adjoint solve");
  const VecX w = fac.solve(g);

  std::size_t target_bar = 0;
  if (response.what == ResponseKind::What::BarForce ||
      response.what == ResponseKind::What::BarDisplacement)
    target_bar = static_cast<std::size_t>(model.bar_index(response.bar));
  const int reaction_dof = response.what == ResponseKind::What::Reaction
                               ? d.full(model.node_index(response.node), response.comp)
                               : -1;
  std::vector<double> eta;
  eta.reserve(unit.size());
  std::vector<int> touched;
  for (std::size_t k = 0; k < unit.size(); ++k) {
    const Load& l = unit[k];
    if (hooks && k % 64 == 0) {
      hooks->check();
      hooks->report(0.5 + 0.5 * static_cast<double>(k) / static_cast<double>(unit.size()),
                    "ordinates " + std::to_string(k) + "/" + std::to_string(unit.size()));
    }
    require(l.kind == Load::Kind::NodalForce || l.kind == Load::Kind::BarPoint,
            "influence: unit loads are nodal or bar point loads");
    lower_load(model, asm_, "unit", l, cl);
    load_dofs(model, l, d, touched);
    double v = 0;
    for (int f : touched) {
      int fi = d.free_ix[static_cast<std::size_t>(f)];
      if (fi >= 0) v += w(fi) * cl.rhs(fi);
    }
    switch (response.what) {
      case ResponseKind::What::Displacement: break;
      case ResponseKind::What::Reaction: v -= cl.nodal_full(reaction_dof); break;
      case ResponseKind::What::BarForce:
      case ResponseKind::What::BarDisplacement:
        if (l.kind == Load::Kind::BarPoint &&
            static_cast<std::size_t>(model.bar_index(l.bar)) == target_bar)
          v += extract_response(model, asm_, cl, probe, response);  // probe is all-zero here
        break;
    }
    require(std::isfinite(v), "influence: non-finite ordinate");
    eta.push_back(v);
    // reset exactly what lower_load wrote
    for (int f : touched) {
      cl.nodal_full(f) = 0;
      int fi = d.free_ix[static_cast<std::size_t>(f)];
      if (fi >= 0) cl.rhs(fi) = 0;
    }
    if (l.kind == Load::Kind::BarPoint) {
      std::size_t e = static_cast<std::size_t>(model.bar_index(l.bar));
      cl.f0[e] = Vec12::Zero();
      cl.span[e] = BarSpanLoads{};
    }
  }
  hooks_report(hooks, 1.0, "done");
  return eta;
}

InfluenceLine influence_line(const Model& model, const LoadPath& path, const Vec3& dir,
                             const ResponseKind& response) {
  FactorCache fac;
  return influence_line(model, path, dir, response, fac);
}

InfluenceLine influence_line(const Model& model, const LoadPath& path, const Vec3& dir,
                             const ResponseKind& response, FactorCache& fac,
                             const SolveHooks* hooks) {
  std::vector<Load> unit;
  unit.reserve(path.stations.size());
  for (const PathStation& st : path.stations) unit.push_back(Load::point(st.bar, st.x_local, dir));
  // The slopes at both ends of each interval, one-sided second-order
  // differences of ordinates at e and 2e inside the interval, on the bar that
  // carries it (the line's kinks at nodes and its jump at the target section
  // stay at the samples), read through the same adjoint solution.
  const std::size_t n = path.stations.size();
  auto bar_end = [&](const PathStation& st) {
    const double Lb = model.bar_chord_length(model.bar(st.bar));
    return st.x_local <= 1e-9 * Lb || st.x_local >= Lb * (1 - 1e-9);
  };
  struct Iv {
    bool ok = false;
    double e = 0;
  };
  std::vector<Iv> ivs(n > 0 ? n - 1 : 0);
  // The stations of each bar, in path order: the direction search below
  // visits only the stations of the interval's bar (the first one, in path
  // order, at another abscissa), O(stations per bar) instead of a scan of the
  // whole path per interval (O(n^2) over the line).
  std::unordered_map<Id, std::vector<std::size_t>> on_bar;
  for (std::size_t q = 0; q < n; ++q) on_bar[path.stations[q].bar].push_back(q);
  for (std::size_t k = 0; k + 1 < n; ++k) {
    const PathStation &a = path.stations[k], &b = path.stations[k + 1];
    const double h = b.s - a.s;
    if (h <= 0) continue;
    // the bar of the interval and a reference station on it
    const PathStation& ref = a.bar == b.bar ? a : (bar_end(a) ? b : a);
    // x_local(s) on that bar: the sign of dx/ds from the reference station
    double dir_s = 0;
    for (std::size_t q : on_bar[ref.bar])
      if (path.stations[q].s != ref.s) {
        dir_s = (path.stations[q].x_local - ref.x_local) / (path.stations[q].s - ref.s) > 0 ? 1.0 : -1.0;
        break;
      }
    if (dir_s == 0) continue;
    const double e = 1e-3 * h;
    ivs[k].ok = true;
    ivs[k].e = e;
    for (double sp : {a.s + e, a.s + 2 * e, b.s - e, b.s - 2 * e})
      unit.push_back(Load::point(ref.bar, ref.x_local + dir_s * (sp - ref.s), dir));
  }
  const std::vector<double> all = unit_load_sweep(model, unit, response, fac, hooks);
  InfluenceLine il;
  il.eta.assign(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(n));
  il.s.reserve(n);
  for (const PathStation& st : path.stations) il.s.push_back(st.s);
  il.dl.assign(n, 0.0);
  il.dr.assign(n, 0.0);
  std::size_t base = n;
  for (std::size_t k = 0; k + 1 < n; ++k) {
    if (!ivs[k].ok) continue;
    const double e = ivs[k].e;
    const double f1 = all[base], f2 = all[base + 1], g1 = all[base + 2], g2 = all[base + 3];
    base += 4;
    // the ordinates of the interval's own ends, approached from inside
    // (second-order extrapolation: a jump at the target section stays out)
    const double fa = 2 * f1 - f2, gb = 2 * g1 - g2;
    il.dr[k] = (-3 * fa + 4 * f1 - f2) / (2 * e);
    il.dl[k + 1] = (3 * gb - 4 * g1 + g2) / (2 * e);
  }
  return il;
}

// ---------------- convoys ----------------

ConvoyPattern ConvoyPattern::mirrored() const {
  ConvoyPattern m;
  for (const auto& [a, p] : axles) m.axles.emplace_back(-a, p);
  for (const auto& seg : segments) m.segments.push_back({-(seg[0] + seg[1]), seg[1], seg[2]});
  return m;
}

Extremes convoy_extreme(const InfluenceLine& il, const ConvoyPattern& convoy, int steps) {
  require(steps >= 1, "convoy: need steps >= 1");
  const IlPrefix prefix(il);
  double a_min = 0, a_max = 0;
  for (const auto& [a, p] : convoy.axles) {
    a_min = std::min(a_min, a);
    a_max = std::max(a_max, a);
  }
  for (const auto& seg : convoy.segments) {
    a_min = std::min(a_min, seg[0]);
    a_max = std::max(a_max, seg[0] + seg[1]);
  }
  double r0 = il.s.front() - a_max, r1 = il.s.back() - a_min;
  Extremes ex;
  bool first = true;
  for (int k = 0; k <= steps; ++k) {
    double r = r0 + (r1 - r0) * k / steps;
    double E = 0;
    for (const auto& [a, p] : convoy.axles) E += p * il.at(r + a);
    for (const auto& seg : convoy.segments)
      if (seg[1] > 0) E += seg[2] / seg[1] * prefix.integral(r + seg[0], r + seg[0] + seg[1]);
    if (first || E > ex.max) {
      ex.max = std::max(first ? E : ex.max, E);
      if (first || E >= ex.max) ex.r_max = r;
    }
    if (first || E < ex.min) {
      ex.min = std::min(first ? E : ex.min, E);
      ex.r_min = r;
    }
    first = false;
  }
  return ex;
}

Extremes udl_extreme(const InfluenceLine& il, double q) {
  Extremes ex;
  ex.max = q * il.integral_pos();
  ex.min = q * il.integral_neg();
  return ex;
}

}  // namespace keta

namespace keta {

namespace {
// The equivalent nodal forces of a convoy at the reference position r(step):
// one CaseLoads per step, every axle lowered into it as a BarPoint (its
// station pair found by bisection), every segment as the UDL pieces of the
// bars it covers — O(n + (axles + segment pieces) log S) per step.
std::function<VecX(int)> convoy_force_factory(const Model& model, const LoadPath& path, const ConvoyPattern& convoy,
                                              const Vec3& dir, std::function<double(int)> position) {
  require(path.stations.size() >= 2, "moving convoy: the path needs 2 stations");
  auto asm_ptr = std::make_shared<Assembled>(assemble_elements(model, true));
  auto model_ptr = &model;
  LoadPath p = path;
  ConvoyPattern c = convoy;
  Vec3 d = dir;
  return [asm_ptr, model_ptr, p, c, d, position](int step) -> VecX {
    const Assembled& asm_ = *asm_ptr;
    CaseLoads cl = empty_case_loads(asm_);
    const double r = position(step);  // reference abscissa of the convoy
    auto place = [&](double s, double magnitude) {
      if (s < p.stations.front().s || s > p.stations.back().s) return;
      // the first station k >= 1 with s <= s_k (+ 1e-12)
      const auto it = std::lower_bound(p.stations.begin() + 1, p.stations.end(), s - 1e-12,
                                       [](const PathStation& st, double v) { return st.s < v; });
      if (it == p.stations.end()) return;
      const PathStation& a = *(it - 1);
      const PathStation& b = *it;
      const double t = (s - a.s) / std::max(1e-12, b.s - a.s);
      // stations within one bar share the bar; at bar joints pick the side
      const Id bar = (a.bar == b.bar || t > 0) ? b.bar : a.bar;
      const double xl = a.bar == b.bar ? a.x_local + t * (b.x_local - a.x_local) : (t == 0 ? a.x_local : b.x_local * t);
      lower_load(*model_ptr, asm_, "convoy", Load::point(bar, std::clamp(xl, 0.0, 1e9), magnitude * d), cl);
    };
    for (const auto& [off, P] : c.axles) place(r + off, P);
    // a segment: its uniform density over [s0, s1], cut along the path's
    // station intervals into the bars that carry it (consistent loads)
    for (const auto& seg : c.segments) {
      if (seg[1] <= 0) continue;
      const double s0 = r + seg[0], s1 = r + seg[0] + seg[1], q = seg[2] / seg[1];
      for (std::size_t k = 1; k < p.stations.size(); ++k) {
        const PathStation& a = p.stations[k - 1];
        const PathStation& b = p.stations[k];
        if (a.bar != b.bar || b.s <= a.s) continue;
        const double u0 = std::max(s0, a.s), u1 = std::min(s1, b.s);
        if (u1 <= u0) continue;
        auto xl = [&](double u) { return a.x_local + (u - a.s) / (b.s - a.s) * (b.x_local - a.x_local); };
        const double x0 = xl(u0), x1 = xl(u1);
        lower_load(*model_ptr, asm_, "convoy",
                   Load::udl(a.bar, q * d, LoadFrame::Global, std::min(x0, x1), std::max(x0, x1)), cl);
      }
    }
    return cl.rhs;
  };
}
}  // namespace

std::function<VecX(int)> moving_convoy_force(const Model& model, const LoadPath& path,
                                             const ConvoyPattern& convoy, const Vec3& dir,
                                             double speed, double dt) {
  require(speed > 0 && dt > 0, "moving convoy: bad speed or step");
  return convoy_force_factory(model, path, convoy, dir, [speed, dt](int step) { return speed * dt * step; });
}

double convoy_position_accel(double r0, double v0, double a, double t) {
  require(std::isfinite(r0) && std::isfinite(v0) && std::isfinite(a) && v0 >= 0 && t >= 0,
          "moving convoy: finite r0, v0 >= 0, a, t >= 0");
  // a braking convoy stops at t* = v0 / |a| and stays there (it does not reverse)
  if (a < 0) t = std::min(t, v0 / -a);
  return r0 + v0 * t + 0.5 * a * t * t;
}

std::function<VecX(int)> moving_convoy_force_accel(const Model& model, const LoadPath& path,
                                                   const ConvoyPattern& convoy, const Vec3& dir, double r0,
                                                   double v0, double a, double dt) {
  require(dt > 0 && v0 >= 0 && std::isfinite(a) && std::isfinite(r0), "moving convoy: dt > 0, v0 >= 0, finite a, r0");
  require(v0 > 0 || a > 0, "moving convoy: a convoy at rest needs a > 0 to move");
  return convoy_force_factory(model, path, convoy, dir,
                              [r0, v0, a, dt](int step) { return convoy_position_accel(r0, v0, a, dt * step); });
}

}  // namespace keta

#include "keta/thinwall.h"

#include <algorithm>
#include <cmath>
#include <deque>

namespace keta {

double ThinWall::length(int seg) const {
  const WallSegment& s = segs[static_cast<std::size_t>(seg)];
  return (nodes[static_cast<std::size_t>(s.b)] - nodes[static_cast<std::size_t>(s.a)]).norm();
}

void ThinWall::validate() const {
  require(nodes.size() >= 2, "thinwall: at least 2 midline nodes are needed");
  require(!segs.empty(), "thinwall: at least one wall segment is needed");
  for (const Vec2& p : nodes)
    require(std::isfinite(p.x()) && std::isfinite(p.y()), "thinwall: non-finite midline node");
  int n = static_cast<int>(nodes.size());
  for (std::size_t i = 0; i < segs.size(); ++i) {
    const WallSegment& s = segs[i];
    require(s.a >= 0 && s.a < n && s.b >= 0 && s.b < n,
            "thinwall: segment " + std::to_string(i) + " references a missing node");
    require(s.a != s.b, "thinwall: segment " + std::to_string(i) + " is a loop on one node");
    require(s.t > 0 && std::isfinite(s.t),
            "thinwall: segment " + std::to_string(i) + " needs a positive thickness");
    require(s.t_b >= 0 && std::isfinite(s.t_b),
            "thinwall: segment " + std::to_string(i) + " needs a thickness at node b that is positive, or 0 for none");
    require((nodes[static_cast<std::size_t>(s.b)] - nodes[static_cast<std::size_t>(s.a)]).norm() >
                0,
            "thinwall: segment " + std::to_string(i) + " has zero length");
  }
  // Connectivity: Saint-Venant torsion and the shear-flow balance couple the
  // whole midline; a disconnected graph is a modelling error, not a section.
  std::vector<std::vector<int>> adj(nodes.size());
  for (const WallSegment& s : segs) {
    adj[static_cast<std::size_t>(s.a)].push_back(s.b);
    adj[static_cast<std::size_t>(s.b)].push_back(s.a);
  }
  std::vector<char> seen(nodes.size(), 0);
  std::deque<int> queue{0};
  seen[0] = 1;
  std::size_t count = 1;
  while (!queue.empty()) {
    int u = queue.front();
    queue.pop_front();
    for (int v : adj[static_cast<std::size_t>(u)])
      if (!seen[static_cast<std::size_t>(v)]) {
        seen[static_cast<std::size_t>(v)] = 1;
        ++count;
        queue.push_back(v);
      }
  }
  require(count == nodes.size(), "thinwall: the midline graph is not connected");
}

namespace {

// Everything the four public entry points share: centroid-shifted geometry, the
// section moments, a spanning tree and the fundamental-cycle basis.
struct Frame {
  const ThinWall* w = nullptr;
  std::vector<Vec2> p;         // centroid-shifted nodes
  std::vector<double> L;       // segment lengths
  double A = 0, yc = 0, zc = 0;
  double Iy = 0, Iz = 0, Iyz = 0, det = 0;
  std::vector<WallCycle> cycles;
  std::vector<char> in_cycle;  // segment belongs to at least one cell
  std::vector<int> parent_node, parent_seg, depth;  // BFS spanning tree
  std::vector<int> order;      // BFS visit order

  int E() const { return static_cast<int>(w->segs.size()); }
  int V() const { return static_cast<int>(p.size()); }
  const WallSegment& seg(int e) const { return w->segs[static_cast<std::size_t>(e)]; }
  Vec2 pa(int e) const { return p[static_cast<std::size_t>(seg(e).a)]; }
  Vec2 pb(int e) const { return p[static_cast<std::size_t>(seg(e).b)]; }
  double length(int e) const { return L[static_cast<std::size_t>(e)]; }

  // Thickness integrals along wall e, in the fraction u of its length.
  // integral over [0, 1] of u^k t(u)
  double t_moment(int e, int k) const {
    const double ta = seg(e).t, tb = seg(e).thickness_at_b();
    return ta / (k + 1) + (tb - ta) / (k + 2);
  }
  // integral over [0, 1] of t(u)^3
  double t_cubed(int e) const {
    const double ta = seg(e).t, tb = seg(e).thickness_at_b();
    return (ta + tb) * (ta * ta + tb * tb) / 4.0;
  }
  // integral over [0, u] of dv / t(v): (u / t_a) ln(1 + x) / x with x the
  // relative thickness change up to u (the constant wall is its limit)
  double inverse_t(int e, double u) const {
    const double ta = seg(e).t, x = u * (seg(e).thickness_at_b() - ta) / ta;
    const double log_ratio = std::fabs(x) < 1e-6 ? 1.0 - x / 2.0 + x * x / 3.0 : std::log1p(x) / x;
    return u / ta * log_ratio;
  }
  // integral over [0, 1] of fn(u): 8-point Gauss-Legendre on pieces of the
  // wall whose thickness ratio stays under kPieceRatio (one piece for a
  // constant wall). Exact for polynomials up to degree 15; the integrands
  // divided by t(u) or carrying ln t(u) converge to rounding on such pieces.
  template <class Fn>
  double integrate(int e, const Fn& fn) const {
    static constexpr double kNode[4] = {0.1834346424956498, 0.5255324099163290, 0.7966664774136267,
                                        0.9602898564975363};
    static constexpr double kWeight[4] = {0.3626837833783620, 0.3137066458778873, 0.2223810344533745,
                                          0.1012285362903763};
    constexpr double kPieceRatio = 1.5;
    const double ta = seg(e).t, tb = seg(e).thickness_at_b();
    const double ratio = std::max(ta, tb) / std::min(ta, tb);
    const int pieces = std::max(1, static_cast<int>(std::ceil(std::log(ratio) / std::log(kPieceRatio))));
    // breakpoints in geometric progression of the thickness
    auto breakpoint = [&](int i) {
      if (i == 0) return 0.0;
      if (i == pieces) return 1.0;
      return (ta * std::pow(tb / ta, static_cast<double>(i) / pieces) - ta) / (tb - ta);
    };
    double sum = 0;
    for (int i = 0; i < pieces; ++i) {
      const double u0 = breakpoint(i), half = 0.5 * (breakpoint(i + 1) - u0), mid = u0 + half;
      for (int g = 0; g < 4; ++g)
        sum += half * kWeight[g] * (fn(mid - half * kNode[g]) + fn(mid + half * kNode[g]));
    }
    return sum;
  }
};

void build_tree(Frame& f) {
  int nv = f.V();
  f.parent_node.assign(static_cast<std::size_t>(nv), -1);
  f.parent_seg.assign(static_cast<std::size_t>(nv), -1);
  f.depth.assign(static_cast<std::size_t>(nv), 0);
  std::vector<std::vector<std::pair<int, int>>> adj(static_cast<std::size_t>(nv));
  for (int e = 0; e < f.E(); ++e) {
    adj[static_cast<std::size_t>(f.seg(e).a)].push_back({f.seg(e).b, e});
    adj[static_cast<std::size_t>(f.seg(e).b)].push_back({f.seg(e).a, e});
  }
  std::vector<char> seen(static_cast<std::size_t>(nv), 0);
  std::deque<int> queue{0};
  seen[0] = 1;
  f.order.clear();
  std::vector<char> used(static_cast<std::size_t>(f.E()), 0);
  while (!queue.empty()) {
    int u = queue.front();
    queue.pop_front();
    f.order.push_back(u);
    for (const auto& [v, e] : adj[static_cast<std::size_t>(u)]) {
      if (seen[static_cast<std::size_t>(v)]) continue;
      seen[static_cast<std::size_t>(v)] = 1;
      f.parent_node[static_cast<std::size_t>(v)] = u;
      f.parent_seg[static_cast<std::size_t>(v)] = e;
      f.depth[static_cast<std::size_t>(v)] = f.depth[static_cast<std::size_t>(u)] + 1;
      used[static_cast<std::size_t>(e)] = 1;
      queue.push_back(v);
    }
  }
  // Fundamental cycle of every non-tree segment.
  f.in_cycle.assign(static_cast<std::size_t>(f.E()), 0);
  auto step_up = [&](int n, int& out_seg) {
    int e = f.parent_seg[static_cast<std::size_t>(n)];
    out_seg = e;
    return f.seg(e).a == n ? +1 : -1;  // travelling n -> parent
  };
  for (int e = 0; e < f.E(); ++e) {
    if (used[static_cast<std::size_t>(e)]) continue;
    int a = f.seg(e).a, b = f.seg(e).b;
    std::vector<std::pair<int, int>> up_a, up_b;  // (segment, direction of travel)
    int na = a, nb = b;
    while (f.depth[static_cast<std::size_t>(na)] > f.depth[static_cast<std::size_t>(nb)]) {
      int s;
      int d = step_up(na, s);
      up_a.push_back({s, d});
      na = f.parent_node[static_cast<std::size_t>(na)];
    }
    while (f.depth[static_cast<std::size_t>(nb)] > f.depth[static_cast<std::size_t>(na)]) {
      int s;
      int d = step_up(nb, s);
      up_b.push_back({s, d});
      nb = f.parent_node[static_cast<std::size_t>(nb)];
    }
    while (na != nb) {
      int s;
      int d = step_up(na, s);
      up_a.push_back({s, d});
      na = f.parent_node[static_cast<std::size_t>(na)];
      d = step_up(nb, s);
      up_b.push_back({s, d});
      nb = f.parent_node[static_cast<std::size_t>(nb)];
    }
    // Walk: a -> b through the chord, b -> lca up the tree, lca -> a downwards.
    WallCycle c;
    c.seg.push_back(e);
    c.dir.push_back(+1);
    for (const auto& [s, d] : up_b) {
      c.seg.push_back(s);
      c.dir.push_back(d);
    }
    for (auto it = up_a.rbegin(); it != up_a.rend(); ++it) {
      c.seg.push_back(it->first);
      c.dir.push_back(-it->second);
    }
    double area = 0;
    for (std::size_t i = 0; i < c.seg.size(); ++i) {
      Vec2 from = c.dir[i] > 0 ? f.pa(c.seg[i]) : f.pb(c.seg[i]);
      Vec2 to = c.dir[i] > 0 ? f.pb(c.seg[i]) : f.pa(c.seg[i]);
      area += 0.5 * (from.x() * to.y() - to.x() * from.y());
    }
    c.area = area;
    for (int s : c.seg) f.in_cycle[static_cast<std::size_t>(s)] = 1;
    f.cycles.push_back(std::move(c));
  }
}

Frame build_frame(const ThinWall& w) {
  w.validate();
  Frame f;
  f.w = &w;
  f.p = w.nodes;
  f.L.resize(w.segs.size());
  double A = 0, Sy = 0, Sz = 0;
  for (int e = 0; e < f.E(); ++e) {
    f.L[static_cast<std::size_t>(e)] = w.length(e);
    const double a = f.length(e) * f.t_moment(e, 0);
    // centroid of the wall: at the first moment of its thickness along the midline
    const Vec2 m = f.pa(e) + f.t_moment(e, 1) / f.t_moment(e, 0) * (f.pb(e) - f.pa(e));
    A += a;
    Sz += a * m.x();
    Sy += a * m.y();
  }
  require(A > 0, "thinwall: non-positive wall area");
  f.A = A;
  f.yc = Sz / A;
  f.zc = Sy / A;
  for (Vec2& q : f.p) q -= Vec2(f.yc, f.zc);
  for (int e = 0; e < f.E(); ++e) {
    // midline point a + u d, thickness t(u) across the normal (-d.y, d.x) / L:
    // integral (coordinate)^2 t du along the midline plus the strip term t^3 / 12
    const double L = f.length(e);
    const Vec2 a = f.pa(e), d = f.pb(e) - a, unit = d / L;
    const double T0 = f.t_moment(e, 0), T1 = f.t_moment(e, 1), T2 = f.t_moment(e, 2);
    const double strip = f.t_cubed(e) / 12.0;
    f.Iy += L * (a.y() * a.y() * T0 + 2 * a.y() * d.y() * T1 + d.y() * d.y() * T2 + unit.x() * unit.x() * strip);
    f.Iz += L * (a.x() * a.x() * T0 + 2 * a.x() * d.x() * T1 + d.x() * d.x() * T2 + unit.y() * unit.y() * strip);
    f.Iyz += L * (a.x() * a.y() * T0 + (a.x() * d.y() + a.y() * d.x()) * T1 + d.x() * d.y() * T2 -
                  unit.x() * unit.y() * strip);
  }
  f.det = f.Iy * f.Iz - f.Iyz * f.Iyz;
  require(f.det > 0, "thinwall: degenerate second moments (collinear midline?)");
  build_tree(f);
  return f;
}

// Circulating cell flows for a UNIT twist rate (G theta' = 1); `J_closed` is the
// torque they carry, which is exactly the closed torsion constant.
std::vector<double> unit_cell_flow(const Frame& f, double& J_closed) {
  int nc = static_cast<int>(f.cycles.size());
  J_closed = 0;
  if (nc == 0) return {};
  MatX M = MatX::Zero(nc, nc);
  VecX rhs(nc);
  std::vector<std::vector<double>> c(static_cast<std::size_t>(nc),
                                     std::vector<double>(static_cast<std::size_t>(f.E()), 0.0));
  for (int k = 0; k < nc; ++k) {
    for (std::size_t i = 0; i < f.cycles[static_cast<std::size_t>(k)].seg.size(); ++i)
      c[static_cast<std::size_t>(k)][static_cast<std::size_t>(
          f.cycles[static_cast<std::size_t>(k)].seg[i])] +=
          f.cycles[static_cast<std::size_t>(k)].dir[i];
    rhs(k) = 2.0 * f.cycles[static_cast<std::size_t>(k)].area;
  }
  for (int k = 0; k < nc; ++k)
    for (int l = 0; l < nc; ++l) {
      double m = 0;
      for (int e = 0; e < f.E(); ++e)
        m += c[static_cast<std::size_t>(k)][static_cast<std::size_t>(e)] *
             c[static_cast<std::size_t>(l)][static_cast<std::size_t>(e)] *
             f.length(e) * f.inverse_t(e, 1.0);
      M(k, l) = m;
    }
  Eigen::FullPivLU<MatX> lu(M);
  require(lu.isInvertible(), "thinwall: singular cell compatibility system");
  VecX x = lu.solve(rhs);
  std::vector<double> q(static_cast<std::size_t>(f.E()), 0.0);
  for (int k = 0; k < nc; ++k) {
    J_closed += 2.0 * x(k) * f.cycles[static_cast<std::size_t>(k)].area;
    for (int e = 0; e < f.E(); ++e)
      q[static_cast<std::size_t>(e)] +=
          x(k) * c[static_cast<std::size_t>(k)][static_cast<std::size_t>(e)];
  }
  // The cell amplitudes are returned appended after the wall flows so the caller
  // can report both without a second solve.
  for (int k = 0; k < nc; ++k) q.push_back(x(k));
  return q;
}

double open_torsion(const Frame& f) {
  double j = 0;
  for (int e = 0; e < f.E(); ++e)
    if (!f.in_cycle[static_cast<std::size_t>(e)])
      j += f.length(e) * f.t_cubed(e) / 3.0;
  return j;
}

// Nodal sectorial coordinate about the CENTROID, eq. (5) of the header.
std::vector<double> sectorial(const Frame& f, const std::vector<double>& q1) {
  std::vector<double> om(static_cast<std::size_t>(f.V()), 0.0);
  // BFS order guarantees a node's parent already carries its value.
  for (std::size_t i = 1; i < f.order.size(); ++i) {
    int v = f.order[i];
    int e = f.parent_seg[static_cast<std::size_t>(v)];
    int u = f.parent_node[static_cast<std::size_t>(v)];
    // Travel u -> v: +1 when the segment is stored a -> b with a == u.
    int dir = f.seg(e).a == u ? +1 : -1;
    Vec2 from = dir > 0 ? f.pa(e) : f.pb(e);
    Vec2 to = dir > 0 ? f.pb(e) : f.pa(e);
    Vec2 m = 0.5 * (from + to);
    Vec2 d = to - from;
    double dom = m.x() * d.y() - m.y() * d.x();
    if (!q1.empty())
      dom -= dir * q1[static_cast<std::size_t>(e)] * f.length(e) * f.inverse_t(e, 1.0);
    om[static_cast<std::size_t>(v)] = om[static_cast<std::size_t>(u)] + dom;
  }
  return om;
}

// The sectorial coordinate at the fraction u of wall e, from its nodal values:
// a straight line, but for the torsion-flow term -q1 integral ds/t of eq. (5),
// which bends it inside a cell when the wall tapers.
double omega_on_wall(const Frame& f, const std::vector<double>& om, const std::vector<double>& q1, int e, double u) {
  const double oa = om[static_cast<std::size_t>(f.seg(e).a)], ob = om[static_cast<std::size_t>(f.seg(e).b)];
  const double flow = q1.empty() ? 0.0 : q1[static_cast<std::size_t>(e)] * f.length(e);
  return oa + u * (ob - oa + flow * f.inverse_t(e, 1.0)) - flow * f.inverse_t(e, u);
}

// Area integrals of the sectorial coordinate over the section.
struct SectorialInt {
  double one = 0;  // integral omega dA
  double y = 0;    // integral omega y dA
  double z = 0;    // integral omega z dA
  double sq = 0;   // integral omega^2 dA
};
SectorialInt sectorial_integrals(const Frame& f, const std::vector<double>& om, const std::vector<double>& q1) {
  SectorialInt out;
  for (int e = 0; e < f.E(); ++e) {
    const Vec2 a = f.pa(e), d = f.pb(e) - a;
    const double L = f.length(e);
    auto over_area = [&](auto&& weight) {
      return L * f.integrate(e, [&](double u) {
        const double omega = omega_on_wall(f, om, q1, e, u);
        return f.seg(e).thickness(u) * weight(omega, u);
      });
    };
    out.one += over_area([](double omega, double) { return omega; });
    out.sq += over_area([](double omega, double) { return omega * omega; });
    out.y += over_area([&](double omega, double u) { return omega * (a.x() + u * d.x()); });
    out.z += over_area([&](double omega, double u) { return omega * (a.y() + u * d.y()); });
  }
  return out;
}

// Shear-flow coefficients of eq. (3).
struct ShearK {
  double Kz = 0, Ky = 0;
};
ShearK shear_k(const Frame& f, double Vy, double Vz) {
  ShearK k;
  k.Kz = (Vz * f.Iz - Vy * f.Iyz) / f.det;
  k.Ky = (Vy * f.Iy - Vz * f.Iyz) / f.det;
  return k;
}

// Accumulated wall load between the a-end and the fraction u of segment e:
// -(Kz Qz + Ky Qy) with Q the first moment of the wall area up to u,
// L integral (a + v d) t(v) dv — cubic in u when the wall tapers.
double accum(const Frame& f, int e, double u, const ShearK& k) {
  const Vec2 a = f.pa(e), d = f.pb(e) - a;
  const double ta = f.seg(e).t, dt = f.seg(e).thickness_at_b() - ta;
  auto first_moment = [&](double c, double dc) {
    return f.length(e) * u * (ta * c + u * ((ta * dc + dt * c) / 2.0 + u * dt * dc / 3.0));
  };
  return -(k.Kz * first_moment(a.y(), d.y()) + k.Ky * first_moment(a.x(), d.x()));
}

// The E wall-start flows: V-1 node balances + one zero-twist condition per cell.
std::vector<double> shear_start_flows(const Frame& f, const ShearK& k) {
  int E = f.E(), V = f.V(), nc = static_cast<int>(f.cycles.size());
  require(E == V - 1 + nc, "thinwall: inconsistent midline topology");
  MatX M = MatX::Zero(E, E);
  VecX rhs = VecX::Zero(E);
  int row = 0;
  for (int n = 1; n < V; ++n, ++row) {
    for (int e = 0; e < E; ++e) {
      if (f.seg(e).b == n) {
        M(row, e) += 1.0;
        rhs(row) -= accum(f, e, 1.0, k);
      }
      if (f.seg(e).a == n) M(row, e) -= 1.0;
    }
  }
  for (int c = 0; c < nc; ++c, ++row) {
    const WallCycle& cy = f.cycles[static_cast<std::size_t>(c)];
    for (std::size_t i = 0; i < cy.seg.size(); ++i) {
      int e = cy.seg[i];
      // zero twist: the circuit integral of q / t, q = q0 + accum along the wall
      const double w = cy.dir[i] * f.length(e);
      M(row, e) += w * f.inverse_t(e, 1.0);
      rhs(row) -= w * f.integrate(e, [&](double u) { return accum(f, e, u, k) / f.seg(e).thickness(u); });
    }
  }
  Eigen::FullPivLU<MatX> lu(M);
  require(lu.isInvertible(), "thinwall: singular shear-flow system");
  VecX x = lu.solve(rhs);
  return std::vector<double>(x.data(), x.data() + E);
}

// Shear-flow energy, the contour integral of q^2 / t.
double flow_energy(const Frame& f, const std::vector<double>& q0, const ShearK& k) {
  double energy = 0;
  for (int e = 0; e < f.E(); ++e)
    energy += f.length(e) * f.integrate(e, [&](double u) {
      const double q = q0[static_cast<std::size_t>(e)] + accum(f, e, u, k);
      return q * q / f.seg(e).thickness(u);
    });
  return energy;
}

}  // namespace

std::vector<WallCycle> thinwall_cycles(const ThinWall& w) { return build_frame(w).cycles; }

TorsionFlow thinwall_torsion(const ThinWall& w, double T) {
  require(std::isfinite(T), "thinwall: non-finite torque");
  Frame f = build_frame(w);
  TorsionFlow out;
  double j_closed = 0;
  std::vector<double> unit = unit_cell_flow(f, j_closed);
  int nc = static_cast<int>(f.cycles.size());
  out.J_closed = j_closed;
  out.J_open = open_torsion(f);
  out.J = out.J_closed + out.J_open;
  require(out.J > 0, "thinwall: non-positive torsion constant");
  // The twist rate the applied torque produces: G theta' = T / J. The unit
  // fields above were built for G theta' = 1, so they scale by that factor.
  double rate = T / out.J;
  out.q_seg.assign(static_cast<std::size_t>(f.E()), 0.0);
  out.tau_seg.assign(static_cast<std::size_t>(f.E()), 0.0);
  out.tau_open.assign(static_cast<std::size_t>(f.E()), 0.0);
  for (int e = 0; e < f.E(); ++e) {
    if (nc > 0) {
      out.q_seg[static_cast<std::size_t>(e)] = rate * unit[static_cast<std::size_t>(e)];
      out.tau_seg[static_cast<std::size_t>(e)] =
          out.q_seg[static_cast<std::size_t>(e)] / std::min(f.seg(e).t, f.seg(e).thickness_at_b());
    }
    // Saint-Venant part: tau_max = G theta' t on the wall surface. Present on
    // every open branch; on a cell wall it is the (small) secondary term.
    if (!f.in_cycle[static_cast<std::size_t>(e)])
      out.tau_open[static_cast<std::size_t>(e)] = rate * std::max(f.seg(e).t, f.seg(e).thickness_at_b());
  }
  for (int k = 0; k < nc; ++k)
    out.q_cell.push_back(rate * unit[static_cast<std::size_t>(f.E() + k)]);
  return out;
}

std::vector<WallStation> thinwall_shear_flow(const ThinWall& w, double Vy, double Vz,
                                             int per_segment) {
  require(std::isfinite(Vy) && std::isfinite(Vz), "thinwall: non-finite shear force");
  require(per_segment >= 2, "thinwall: per_segment must be >= 2");
  Frame f = build_frame(w);
  ShearK k = shear_k(f, Vy, Vz);
  std::vector<double> q0 = shear_start_flows(f, k);
  std::vector<WallStation> out;
  out.reserve(static_cast<std::size_t>(f.E() * per_segment));
  for (int e = 0; e < f.E(); ++e) {
    for (int i = 0; i < per_segment; ++i) {
      double u = static_cast<double>(i) / (per_segment - 1);
      WallStation st;
      st.seg = e;
      st.s = u * f.L[static_cast<std::size_t>(e)];
      Vec2 p = f.pa(e) + u * (f.pb(e) - f.pa(e));
      st.y = p.x();
      st.z = p.y();
      st.q = q0[static_cast<std::size_t>(e)] + accum(f, e, u, k);
      st.tau = st.q / f.seg(e).thickness(u);
      out.push_back(st);
    }
  }
  return out;
}

std::vector<WallStation> thinwall_warping(const ThinWall& w, double B, int per_segment) {
  require(std::isfinite(B), "thinwall: non-finite bimoment");
  require(per_segment >= 2, "thinwall: per_segment must be >= 2");
  ThinWallProps pr = thinwall_props(w);
  require(pr.Iw > 0, "thinwall: the section has a vanishing warping constant");
  Frame f = build_frame(w);
  double j_closed = 0;
  std::vector<double> unit = unit_cell_flow(f, j_closed);
  unit.resize(static_cast<std::size_t>(f.E()));
  std::vector<double> om = sectorial(f, unit);
  // Re-refer to the shear centre and give zero mean (eqs. 6-7).
  for (int n = 0; n < f.V(); ++n)
    om[static_cast<std::size_t>(n)] += -pr.y_sc * f.p[static_cast<std::size_t>(n)].y() +
                                       pr.z_sc * f.p[static_cast<std::size_t>(n)].x();
  const double mean = sectorial_integrals(f, om, unit).one / f.A;
  for (double& v : om) v -= mean;

  std::vector<WallStation> out;
  for (int e = 0; e < f.E(); ++e) {
    for (int i = 0; i < per_segment; ++i) {
      double u = static_cast<double>(i) / (per_segment - 1);
      WallStation st;
      st.seg = e;
      st.s = u * f.L[static_cast<std::size_t>(e)];
      Vec2 p = f.pa(e) + u * (f.pb(e) - f.pa(e));
      st.y = p.x();
      st.z = p.y();
      st.omega = omega_on_wall(f, om, unit, e, u);
      st.sigma_w = B * st.omega / pr.Iw;
      out.push_back(st);
    }
  }
  return out;
}

ThinWallProps thinwall_props(const ThinWall& w) {
  Frame f = build_frame(w);
  ThinWallProps p;
  p.A = f.A;
  p.yc = f.yc;
  p.zc = f.zc;
  p.Iy = f.Iy;
  p.Iz = f.Iz;
  p.Iyz = f.Iyz;
  double avg = 0.5 * (f.Iy + f.Iz);
  double rad = std::sqrt(0.25 * (f.Iy - f.Iz) * (f.Iy - f.Iz) + f.Iyz * f.Iyz);
  p.I1 = avg + rad;
  p.I2 = avg - rad;
  p.phi = std::atan2(f.Iy - p.I1, f.Iyz);
  p.n_cells = static_cast<int>(f.cycles.size());

  double j_closed = 0;
  std::vector<double> unit = unit_cell_flow(f, j_closed);
  unit.resize(static_cast<std::size_t>(f.E()));
  p.J_closed = j_closed;
  p.J_open = open_torsion(f);
  p.J = p.J_closed + p.J_open;

  // Sectorial coordinate about the centroid -> shear centre -> Iw.
  std::vector<double> om = sectorial(f, unit);
  const SectorialInt li = sectorial_integrals(f, om, unit);
  double u = (f.Iyz * li.y - f.Iz * li.z) / f.det;   // y_P - y_S with y_P = 0
  double v = (f.Iy * li.y - f.Iyz * li.z) / f.det;   // z_P - z_S
  p.y_sc = -u;
  p.z_sc = -v;
  for (int n = 0; n < f.V(); ++n)
    om[static_cast<std::size_t>(n)] += u * f.p[static_cast<std::size_t>(n)].y() -
                                       v * f.p[static_cast<std::size_t>(n)].x();
  const SectorialInt ls = sectorial_integrals(f, om, unit);
  p.Iw = ls.sq - ls.one * ls.one / f.A;

  // Shear areas from the flow energy of a unit force, eq. (4).
  {
    ShearK ky = shear_k(f, 1.0, 0.0);
    std::vector<double> q0 = shear_start_flows(f, ky);
    double e_y = flow_energy(f, q0, ky);
    p.Avy = e_y > 0 ? 1.0 / e_y : 0.0;
    ShearK kz = shear_k(f, 0.0, 1.0);
    q0 = shear_start_flows(f, kz);
    double e_z = flow_energy(f, q0, kz);
    p.Avz = e_z > 0 ? 1.0 / e_z : 0.0;
  }

  double ymax = 0, zmax = 0;
  for (const Vec2& q : f.p) {
    ymax = std::max(ymax, std::fabs(q.x()));
    zmax = std::max(zmax, std::fabs(q.y()));
  }
  p.Wely = zmax > 0 ? p.Iy / zmax : 0;
  p.Welz = ymax > 0 ? p.Iz / ymax : 0;
  return p;
}

}  // namespace keta

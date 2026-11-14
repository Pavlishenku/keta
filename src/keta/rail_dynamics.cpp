#include "keta/rail_dynamics.h"

#include <algorithm>
#include <cmath>

namespace keta {

ConvoyPattern rail_hslm_a(int train) {
  // EN 1991-2 Table 6.3: N, D [m], d [m], P [kN]
  static const double tab[10][4] = {{18, 18, 2.0, 170}, {17, 19, 3.5, 200}, {16, 20, 2.0, 180}, {15, 21, 3.0, 190},
                                    {14, 22, 2.0, 170}, {13, 23, 2.0, 180}, {13, 24, 2.0, 190}, {12, 25, 2.5, 190},
                                    {11, 26, 2.0, 210}, {11, 27, 2.0, 210}};
  require(train >= 1 && train <= 10, "HSLM-A: train 1..10");
  const double* t = tab[train - 1];
  const int N = static_cast<int>(t[0]);
  const double D = t[1], d = t[2], P = t[3] * 1e3;
  std::vector<double> xs{0.0, 3.0, 14.0, 17.0};
  const double lead = 17.0 + 3.525;
  for (int j = 0; j <= N + 2; ++j) {
    xs.push_back(lead + j * D);
    xs.push_back(lead + j * D + d);
  }
  const double tail = lead + (N + 2) * D + d + 3.525;
  for (double x : {0.0, 3.0, 14.0, 17.0}) xs.push_back(tail + x);
  ConvoyPattern c;
  for (double x : xs) c.axles.emplace_back(-x, P);
  return c;
}


namespace {
// one bar of the path: s in [s0, s0 + L], the modal force of a unit load
// at the four sample abscissas (Lagrange-cubic in the path abscissa)
struct PathBar {
  double s0 = 0, L = 0;
  std::array<VecX, 4> g;  // per sample t = 0, 1/3, 2/3, 1: modal forces
};

struct SweepContext {
  int nm = 0;
  std::vector<PathBar> bars;
  double path_len = 0;
  MatX phi;  // observed DOFs x modes
  double wmax = 0;
  std::vector<double> zeta;

  SweepContext(const Model& model, const ModalResults& modal, const std::vector<Id>& path_bars, const Vec3& dir,
               const RailSweepOptions& opt) {
    require(!modal.modes_free.empty(), "rail sweep: no modes");
    require(!opt.speeds.empty() && !opt.observe.empty(), "rail sweep: speeds and observed DOFs needed");
    require(opt.steps_per_period >= 4 && opt.max_travel > 0 && opt.free_time >= 0, "rail sweep: bad step options");
    nm = static_cast<int>(modal.modes_free.size());
    const Assembled asm_ = assemble_elements(model, true);
    require(modal.modes_free.front().size() == asm_.dofs.n_free, "rail sweep: the modes belong to another model");
    const LoadPath lp = build_path(model, path_bars, 4);
    for (std::size_t bi = 0; bi < path_bars.size(); ++bi) {
      const Bar& b = model.bar(path_bars[bi]);
      const double L = model.bar_chord_length(b);
      std::vector<const PathStation*> st;
      for (const PathStation& p : lp.stations)
        if (p.bar == b.id) st.push_back(&p);
      require(st.size() >= 3, "rail sweep: path stations");
      const PathStation& last = *st.back();
      const bool flip = std::fabs(last.x_local) < 0.5 * L;
      PathBar pb;
      pb.L = L;
      pb.s0 = last.s - L;
      for (int k = 0; k < 4; ++k) {
        const double t = k / 3.0, x = flip ? L * (1 - t) : L * t;
        LoadCase lc;
        lc.id = 1;
        lc.loads.push_back(Load::point(b.id, std::clamp(x, 0.0, L), dir));
        const CaseLoads cl = lower_case(model, asm_, lc);
        pb.g[static_cast<std::size_t>(k)] = VecX(nm);
        for (int n = 0; n < nm; ++n)
          pb.g[static_cast<std::size_t>(k)](n) = modal.modes_free[static_cast<std::size_t>(n)].dot(cl.rhs);
      }
      bars.push_back(pb);
    }
    path_len = lp.length;
    const std::size_t no = opt.observe.size();
    phi = MatX::Zero(static_cast<Eigen::Index>(no), nm);
    for (std::size_t o = 0; o < no; ++o) {
      require(opt.observe[o].second >= 0 && opt.observe[o].second < 6, "rail sweep: dof 0..5");
      const int fi = asm_.dofs.free_of(model.node_index(opt.observe[o].first), opt.observe[o].second);
      if (fi < 0) continue;
      for (int n = 0; n < nm; ++n)
        phi(static_cast<Eigen::Index>(o), n) = modal.modes_free[static_cast<std::size_t>(n)](fi);
    }
    for (double w : modal.omegas) wmax = std::max(wmax, w);
    zeta.resize(static_cast<std::size_t>(nm));
    for (int n = 0; n < nm; ++n)
      zeta[static_cast<std::size_t>(n)] =
          opt.damping.zeta_for(modal.omegas[static_cast<std::size_t>(n)], static_cast<std::size_t>(n));
  }

  // the modal forces of a unit load at s and their first and second
  // derivatives along the path (zero off the path)
  bool shape(double s, VecX& g, VecX* g1 = nullptr, VecX* g2 = nullptr) const {
    g.setZero(nm);
    if (g1) g1->setZero(nm);
    if (g2) g2->setZero(nm);
    if (s < 0 || s > path_len) return false;
    std::size_t k = 0;
    while (k + 1 < bars.size() && s > bars[k].s0 + bars[k].L) ++k;
    const PathBar& b = bars[k];
    const double t = std::clamp((s - b.s0) / b.L, 0.0, 1.0);
    const double T[4] = {0, 1.0 / 3, 2.0 / 3, 1};
    for (int i = 0; i < 4; ++i) {
      double den = 1;
      // l_i(t) = prod_{j != i} (t - T_j) / (T_i - T_j), with its derivatives
      double a[3];
      int m = 0;
      for (int j = 0; j < 4; ++j)
        if (j != i) {
          den *= T[i] - T[j];
          a[m++] = t - T[j];
        }
      const double l = a[0] * a[1] * a[2];
      const double d1 = a[1] * a[2] + a[0] * a[2] + a[0] * a[1];
      const double d2 = 2 * (a[0] + a[1] + a[2]);
      g += (l / den) * b.g[static_cast<std::size_t>(i)];
      if (g1) *g1 += (d1 / den / b.L) * b.g[static_cast<std::size_t>(i)];
      if (g2) *g2 += (d2 / den / (b.L * b.L)) * b.g[static_cast<std::size_t>(i)];
    }
    return true;
  }
};

void finish(RailSweepResult& out, RailSweepPoint pt) {
  for (std::size_t o = 0; o < pt.u_max.size(); ++o) {
    const double st = std::max(std::fabs(out.static_max[o]), std::fabs(out.static_min[o]));
    const double dy = std::max(std::fabs(pt.u_max[o]), std::fabs(pt.u_min[o]));
    if (st > 0 && dy / st > out.dynamic_factor[o]) {
      out.dynamic_factor[o] = dy / st;
      out.critical_speed[o] = pt.speed;
    }
  }
  out.points.push_back(std::move(pt));
}

RailSweepResult static_envelope(const Model& model, const std::vector<Id>& path_bars, const ConvoyPattern& train,
                                const Vec3& dir, const RailSweepOptions& opt) {
  RailSweepResult out;
  const LoadPath ilp = build_path(model, path_bars, 21);
  for (const auto& [node, c] : opt.observe) {
    const InfluenceLine il = influence_line(model, ilp, dir, ResponseKind::displacement(node, c));
    const Extremes ex = convoy_extreme(il, train, 4000);
    out.static_max.push_back(ex.max);
    out.static_min.push_back(ex.min);
  }
  out.dynamic_factor.assign(opt.observe.size(), 0.0);
  out.critical_speed.assign(opt.observe.size(), 0.0);
  return out;
}

void record(RailSweepPoint& pt, const VecX& u, const VecX& a) {
  for (std::size_t o = 0; o < pt.u_max.size(); ++o) {
    pt.u_max[o] = std::max(pt.u_max[o], u(static_cast<Eigen::Index>(o)));
    pt.u_min[o] = std::min(pt.u_min[o], u(static_cast<Eigen::Index>(o)));
    pt.a_max[o] = std::max(pt.a_max[o], std::fabs(a(static_cast<Eigen::Index>(o))));
  }
}

RailSweepPoint new_point(double v, std::size_t no) {
  RailSweepPoint pt;
  pt.speed = v;
  pt.u_max.assign(no, -1e300);
  pt.u_min.assign(no, 1e300);
  pt.a_max.assign(no, 0.0);
  return pt;
}
}  // namespace

RailSweepResult rail_speed_sweep(const Model& model, const ModalResults& modal, const std::vector<Id>& path_bars,
                                 const ConvoyPattern& train, const Vec3& dir, const RailSweepOptions& opt) {
  require(!train.axles.empty() && train.segments.empty(), "rail sweep: a train of point axles");
  const SweepContext cx(model, modal, path_bars, dir, opt);
  const int nm = cx.nm;
  double off_max = -1e300, off_min = 1e300;
  for (const auto& [off, P] : train.axles) {
    off_max = std::max(off_max, off);
    off_min = std::min(off_min, off);
  }
  const double train_len = off_max - off_min;
  RailSweepResult out = static_envelope(model, path_bars, train, dir, opt);
  VecX g(nm);
  for (double v : opt.speeds) {
    require(v > 0, "rail sweep: speeds > 0");
    const double T = (cx.path_len + train_len) / v + opt.free_time;
    double dt = std::min(2 * kPi / cx.wmax / opt.steps_per_period, opt.max_travel / v);
    const int nsteps = static_cast<int>(std::ceil(T / dt));
    dt = T / nsteps;
    VecX q = VecX::Zero(nm), qd = VecX::Zero(nm), p0 = VecX::Zero(nm), p1(nm);
    RailSweepPoint pt = new_point(v, opt.observe.size());
    auto load_at = [&](double t, VecX& p) {
      p.setZero();
      for (const auto& [off, P] : train.axles)
        if (cx.shape(v * t + off - off_max, g)) p += P * g;
    };
    load_at(0, p0);
    for (int i = 1; i <= nsteps; ++i) {
      load_at(i * dt, p1);
      for (int n = 0; n < nm; ++n)
        oscillator_step(modal.omegas[static_cast<std::size_t>(n)], cx.zeta[static_cast<std::size_t>(n)], dt, p0(n),
                        p1(n), q(n), qd(n));
      VecX qdd(nm);
      for (int n = 0; n < nm; ++n) {
        const double w = modal.omegas[static_cast<std::size_t>(n)];
        qdd(n) = p1(n) - 2 * cx.zeta[static_cast<std::size_t>(n)] * w * qd(n) - w * w * q(n);
      }
      record(pt, cx.phi * q, cx.phi * qdd);
      std::swap(p0, p1);
    }
    finish(out, std::move(pt));
  }
  return out;
}

RailSweepResult rail_speed_sweep_vbi(const Model& model, const ModalResults& modal, const std::vector<Id>& path_bars,
                                     const VehicleTrain& train, const Vec3& dir, const RailSweepOptions& opt) {
  require(!train.axles.empty(), "rail sweep: a train of sprung axles");
  for (const SprungAxle& a : train.axles)
    require(a.m_wheel >= 0 && a.m_body > 0 && a.k > 0 && a.c >= 0, "rail sweep: m_wheel >= 0, m_body > 0, k > 0, c >= 0");
  const SweepContext cx(model, modal, path_bars, dir, opt);
  const int nm = cx.nm, na = static_cast<int>(train.axles.size()), n = nm + na;
  double off_max = -1e300, off_min = 1e300;
  ConvoyPattern weights;
  for (const SprungAxle& a : train.axles) {
    off_max = std::max(off_max, a.offset);
    off_min = std::min(off_min, a.offset);
    weights.axles.emplace_back(a.offset, (a.m_wheel + a.m_body) * train.g);
  }
  const double train_len = off_max - off_min;
  RailSweepResult out = static_envelope(model, path_bars, weights, dir, opt);
  double wveh = 0;
  for (const SprungAxle& a : train.axles) wveh = std::max(wveh, std::sqrt(a.k / a.m_body));
  for (double v : opt.speeds) {
    require(v > 0, "rail sweep: speeds > 0");
    const double T = (cx.path_len + train_len) / v + opt.free_time;
    double dt = std::min({2 * kPi / std::max(cx.wmax, wveh) / opt.steps_per_period, opt.max_travel / v});
    const int nsteps = static_cast<int>(std::ceil(T / dt));
    dt = T / nsteps;
    // state x = [q; y_s], Newmark average acceleration on M x'' + C x' + K x = f
    VecX x = VecX::Zero(n), xd = VecX::Zero(n), xdd = VecX::Zero(n);
    auto system = [&](double t, MatX& M, MatX& C, MatX& K, VecX& f) {
      M = MatX::Zero(n, n);
      C = MatX::Zero(n, n);
      K = MatX::Zero(n, n);
      f = VecX::Zero(n);
      for (int i = 0; i < nm; ++i) {
        const double w = modal.omegas[static_cast<std::size_t>(i)];
        M(i, i) = 1;
        C(i, i) = 2 * cx.zeta[static_cast<std::size_t>(i)] * w;
        K(i, i) = w * w;
      }
      VecX ph(nm), p1(nm), p2(nm);
      for (int k = 0; k < na; ++k) {
        const SprungAxle& a = train.axles[static_cast<std::size_t>(k)];
        const int s = nm + k;
        M(s, s) = a.m_body;
        C(s, s) = a.c;
        K(s, s) = a.k;
        if (!cx.shape(v * t + a.offset - off_max, ph, &p1, &p2)) continue;
        // deck rows: + phi [m_w y_w'' + c y_w' + k y_w - c y_s' - k y_s] = phi G
        M.topLeftCorner(nm, nm) += a.m_wheel * ph * ph.transpose();
        C.topLeftCorner(nm, nm) += a.c * ph * ph.transpose() + 2 * a.m_wheel * v * ph * p1.transpose();
        K.topLeftCorner(nm, nm) += a.k * ph * ph.transpose() + a.c * v * ph * p1.transpose() +
                                   a.m_wheel * v * v * ph * p2.transpose();
        C.block(0, s, nm, 1) -= a.c * ph;
        K.block(0, s, nm, 1) -= a.k * ph;
        f.head(nm) += (a.m_wheel + a.m_body) * train.g * ph;
        // sprung row: m_s y_s'' + c (y_s' - y_w') + k (y_s - y_w) = 0
        C.block(s, 0, 1, nm) -= a.c * ph.transpose();
        K.block(s, 0, 1, nm) -= a.k * ph.transpose() + a.c * v * p1.transpose();
      }
    };
    MatX M, C, K;
    VecX f;
    system(0, M, C, K, f);
    xdd = M.partialPivLu().solve(f - C * xd - K * x);
    RailSweepPoint pt = new_point(v, opt.observe.size());
    for (int i = 1; i <= nsteps; ++i) {
      system(i * dt, M, C, K, f);
      // x1 = x + dt xd + dt^2/4 (xdd + xdd1), xd1 = xd + dt/2 (xdd + xdd1)
      const VecX xp = x + dt * xd + dt * dt / 4 * xdd, vp = xd + dt / 2 * xdd;
      const MatX A = M + dt / 2 * C + dt * dt / 4 * K;
      const VecX xdd1 = A.partialPivLu().solve(f - C * vp - K * xp);
      x = xp + dt * dt / 4 * xdd1;
      xd = vp + dt / 2 * xdd1;
      xdd = xdd1;
      record(pt, cx.phi * x.head(nm), cx.phi * xdd.head(nm));
    }
    finish(out, std::move(pt));
  }
  return out;
}

}  // namespace keta

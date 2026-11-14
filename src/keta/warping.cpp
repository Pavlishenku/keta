// Non-uniform (Vlasov) torsion driver — see warping.h.
#include "keta/warping.h"

#include <algorithm>
#include <cmath>

#include "keta/eigsolve.h"
#include "keta/linsolve.h"

namespace keta {

namespace {

using Vec14 = Eigen::Matrix<double, 14, 1>;
using Mat14 = Eigen::Matrix<double, 14, 14>;
using Vec7 = Eigen::Matrix<double, 7, 1>;
constexpr int kW = 7;  // DOFs per node

// local torsion DOFs of the 14-vector: theta_i, psi_i, theta_j, psi_j
constexpr int kTor[4] = {3, 12, 9, 13};

struct Torsion {
  double G = 0;    // G_t: GJ (+ N r0^2 in second order)
  double EIw = 0;  // E Iw
  double L = 0;
  Eigen::Matrix4d H, Hinv, D;  // basis {1, x, F2, F3}: end values / end forces
  bool warping = false;        // EIw > 0
  // basis values at x: theta, theta', theta'', theta''' of {1, x, F2, F3}
  Eigen::Matrix4d at(double x) const {
    const double a = G / EIw;
    const std::array<double, 7> F = beam_transfer(a, x);
    Eigen::Matrix4d b;
    // rows: theta, theta', theta'', theta'''; cols: 1, x, F2, F3
    b << 1, x, F[2], F[3],  //
        0, 1, F[1], F[2],   //
        0, 0, F[0], F[1],   //
        0, 0, a * F[1], F[0];
    return b;
  }
  void build(double G_t, double EIw_, double L_, Id bar) {
    G = G_t;
    EIw = EIw_;
    L = L_;
    warping = EIw > 0;
    if (!warping) return;
    require(std::fabs(G) * L * L / EIw <= 400.0,
            "warping: bar " + std::to_string(bar) +
                ": |G_t| L^2 / (E Iw) > 400 — the element is too long for the exact transfer "
                "(cancellation e^(kL)); split it");
    const Eigen::Matrix4d b0 = at(0.0), bL = at(L);
    H.row(0) = b0.row(0);
    H.row(1) = b0.row(1);
    H.row(2) = bL.row(0);
    H.row(3) = bL.row(1);
    D.row(0) = -(G * b0.row(1) - EIw * b0.row(3));
    D.row(1) = -EIw * b0.row(2);
    D.row(2) = G * bL.row(1) - EIw * bL.row(3);
    D.row(3) = EIw * bL.row(2);
    Hinv = H.fullPivLu().inverse();
  }
  Eigen::Matrix4d stiffness() const {
    if (!warping) {
      Eigen::Matrix4d k = Eigen::Matrix4d::Zero();
      k(0, 0) = k(2, 2) = G / L;
      k(0, 2) = k(2, 0) = -G / L;
      return k;
    }
    Eigen::Matrix4d k = D * Hinv;
    return 0.5 * (k + k.transpose());
  }
};

// Bar data shared by the drivers.
struct WBar {
  std::size_t e = 0;
  int ni = 0, nj = 0;
  double L = 0, E = 0, G = 0, A = 0, Iy = 0, Iz = 0, J = 0, Iw = 0, r02 = 0;
  Eigen::Vector2d off = Eigen::Vector2d::Zero();
  Mat14 T = Mat14::Identity();  // local <- global (14)
};

struct WSystem {
  const Model& model;
  const WarpingSpec& spec;
  Assembled asm_;
  CaseLoads cl;
  std::vector<WBar> bars;
  int n_nodes = 0;
  std::vector<int> free_ix;  // full (7 n + c) -> free or -1
  int n_free = 0;
  VecX P_direct;  // loads applied at the nodes (full 7 n), span equivalents excluded

  WSystem(const Model& m, Id case_id, const WarpingSpec& s) : model(m), spec(s) {
    require(model.shells.empty(), "warping: facet shells are not supported");
    model.validate();
    asm_ = assemble_elements(model);
    cl = lower_case(model, asm_, model.load_case(case_id));
    require(cl.prescribed.empty(), "warping: settlements are not supported");
    n_nodes = asm_.dofs.n_nodes;
    // bars
    std::vector<std::vector<std::size_t>> at_node(static_cast<std::size_t>(n_nodes));
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      const Bar& b = model.bars[e];
      const std::string tag = "warping: bar " + std::to_string(b.id);
      require(b.kind == BarKind::Frame, tag + " is not a frame bar");
      bool rel = false;
      for (bool r : b.releases) rel |= r;
      require(!rel, tag + " has end releases");
      require(!b.has_ecc(), tag + " has a rigid-arm eccentricity");
      require(!b.soil.any(), tag + " has a Winkler bed");
      require(model.tapered_of(b.id) == nullptr, tag + " is a variable section");
      const BarElement& el = asm_.elements[e];
      const Material& mat = model.material(b.material);
      const SectionProps& sp = model.section(b.section).props;
      WBar w;
      w.e = e;
      w.ni = model.node_index(b.ni);
      w.nj = model.node_index(b.nj);
      w.L = el.L();
      w.E = mat.E;
      w.G = mat.G();
      w.A = sp.A;
      w.Iy = sp.Iy;
      w.Iz = sp.Iz;
      w.J = sp.J;
      auto it = spec.Iw.find(b.section);
      // the spec overrides; else the section's own warping constant
      w.Iw = it == spec.Iw.end() ? sp.Iw : it->second;
      require(w.Iw >= 0, tag + ": negative warping constant");
      require(sp.beta_x == 0.0, tag + ": monosymmetric section (beta_x != 0) — the Vlasov driver covers "
                                       "doubly symmetric sections only");
      w.r02 = (sp.Iy + sp.Iz) / sp.A;
      auto ot = spec.load_offset.find(b.id);
      if (ot != spec.load_offset.end()) w.off = ot->second;
      w.T.setIdentity();
      w.T.topLeftCorner<12, 12>() = el.T();
      at_node[static_cast<std::size_t>(w.ni)].push_back(bars.size());
      at_node[static_cast<std::size_t>(w.nj)].push_back(bars.size());
      bars.push_back(w);
    }
    // DOF map: supports on the 6 classic DOFs; warping fixed by the spec or
    // dropped when no incident bar has Iw > 0
    std::vector<char> fixed(static_cast<std::size_t>(kW * n_nodes), 0);
    for (const Support& s : model.supports) {
      require(!s.spring6, "warping: full elastic supports (APPUI EL) are not supported");
      for (int c = 0; c < 6; ++c) {
        require(s.oneside[static_cast<std::size_t>(c)] == 0, "warping: one-sided supports are not supported");
        if (s.fixed[static_cast<std::size_t>(c)])
          fixed[static_cast<std::size_t>(kW * model.node_index(s.node) + c)] = 1;
      }
    }
    for (Id nd : spec.warping_fixed) fixed[static_cast<std::size_t>(kW * model.node_index(nd) + 6)] = 1;
    for (int n = 0; n < n_nodes; ++n) {
      const auto& inc = at_node[static_cast<std::size_t>(n)];
      bool any_w = false;
      for (std::size_t k : inc) any_w |= bars[k].Iw > 0;
      if (!any_w) {
        fixed[static_cast<std::size_t>(kW * n + 6)] = 1;  // Saint-Venant only: no warping DOF
        continue;
      }
      if (fixed[static_cast<std::size_t>(kW * n + 6)]) continue;
      const Vec3 x0 = asm_.elements[bars[inc.front()].e].R().row(0);
      for (std::size_t k : inc) {
        const Vec3 xk = asm_.elements[bars[k].e].R().row(0);
        require(x0.cross(xk).norm() <= 1e-9,
                "warping: node " + std::to_string(model.nodes[static_cast<std::size_t>(n)].id) +
                    " joins non-collinear bars with free warping (continuity undefined): restrain it");
      }
    }
    // DOFs without any stiffness (rotations of an all-truss node do not occur:
    // frames only) are kept; a zero pivot is a typed singular error.
    free_ix.assign(fixed.size(), -1);
    for (std::size_t f = 0; f < fixed.size(); ++f)
      if (!fixed[f]) free_ix[f] = n_free++;
    // loads applied directly at the nodes
    P_direct = VecX::Zero(kW * n_nodes);
    for (int n = 0; n < n_nodes; ++n)
      for (int c = 0; c < 6; ++c) P_direct(kW * n + c) = cl.nodal_full(asm_.dofs.full(n, c));
    for (const WBar& w : bars) {
      const Vec12 eq = asm_.elements[w.e].T().transpose() * cl.f0[w.e];
      for (int c = 0; c < 6; ++c) {
        P_direct(kW * w.ni + c) -= eq(c);
        P_direct(kW * w.nj + c) -= eq(6 + c);
      }
    }
  }

  int dof(const WBar& w, int k) const {  // local 14 index -> full index
    if (k < 6) return kW * w.ni + k;
    if (k < 12) return kW * w.nj + (k - 6);
    return kW * (k == 12 ? w.ni : w.nj) + 6;
  }

  Vec14 gather(const WBar& w, const VecX& ufull) const {
    Vec14 g;
    for (int k = 0; k < 14; ++k) g(k) = ufull(dof(w, k));
    return w.T * g;  // local
  }
  VecX expand(const VecX& u) const {
    VecX full = VecX::Zero(kW * n_nodes);
    for (std::size_t f = 0; f < free_ix.size(); ++f)
      if (free_ix[f] >= 0) full(static_cast<Eigen::Index>(f)) = u(free_ix[f]);
    return full;
  }
};

// Torsion fixed-end forces (theta_i, psi_i, theta_j, psi_j) of the bar's
// span torques — distributed m_x and point torques plus the torque of the
// offset transverse loads — for the exact torsion block `t` (from rest:
// particular solutions m F4/EIw, P F3/EIw; then the homogeneous closure).
Eigen::Vector4d torsion_f0(const WBar& w, const Torsion& t, const BarSpanLoads& sp, const Vec12& f0_lin) {
  auto mx_of = [&](const Vec3& F, const Vec3& M) { return M(0) + w.off(0) * F(2) - w.off(1) * F(1); };
  Eigen::Vector4d f0 = Eigen::Vector4d::Zero();
  if (!t.warping) {
    // Saint-Venant: fixed-end torques of a clamped shaft, m(x) linear per
    // segment: T_i = int m (L - x)/L, T_j = int m x / L (theta'' = -m/GJ)
    double fi = 0, fj = 0;
    require_linear_span(sp.segs, "warping (Saint-Venant torques)");
    for (const DistSeg& g : sp.segs) {
      const double ma = mx_of(g.qa, g.ma), mb = mx_of(g.qb, g.mb), a = g.a, b = g.b;
      const double I0 = 0.5 * (ma + mb) * (b - a);
      const double I1 = (b - a) * (ma * (2 * a + b) + mb * (a + 2 * b)) / 6.0;  // int m x
      fi += I0 - I1 / t.L;
      fj += I1 / t.L;
    }
    for (const PointAt& p : sp.points) {
      const double P = mx_of(p.F, p.M);
      fi += P * (t.L - p.a) / t.L;
      fj += P * p.a / t.L;
    }
    // the eccentric part is new; the pure m_x part is the linear f0
    (void)f0_lin;
    f0 << fi, 0.0, fj, 0.0;
    return f0;
  }
  const double a_ = t.G / t.EIw, L = t.L;
  // particular solution values at L: theta, theta', theta'', theta'''
  Eigen::Vector4d p = Eigen::Vector4d::Zero();
  auto add_from = [&](double x0, double c0, double c1) {  // (c0 + c1 (x - x0)) / EIw from x0
    if (x0 >= L) return;
    const std::array<double, 7> F = beam_transfer(a_, L - x0);
    p(0) += (c0 * F[4] + c1 * F[5]) / t.EIw;
    p(1) += (c0 * F[3] + c1 * F[4]) / t.EIw;
    p(2) += (c0 * F[2] + c1 * F[3]) / t.EIw;
    p(3) += (c0 * F[1] + c1 * F[2]) / t.EIw;
  };
  require_linear_span(sp.segs, "warping");
  require_no_imposed_twist(sp, "warping");
  for (const DistSeg& g : sp.segs) {
    const double ma = mx_of(g.qa, g.ma), mb = mx_of(g.qb, g.mb);
    if (ma == 0 && mb == 0) continue;
    const double s = g.b > g.a ? (mb - ma) / (g.b - g.a) : 0.0;
    add_from(g.a, ma, s);
    add_from(g.b, -mb, -s);
  }
  for (const PointAt& pt : sp.points) {
    const double P = mx_of(pt.F, pt.M);
    if (P == 0 || pt.a >= L) continue;
    const std::array<double, 7> F = beam_transfer(a_, L - pt.a);
    p(0) += P * F[3] / t.EIw;
    p(1) += P * F[2] / t.EIw;
    p(2) += P * F[1] / t.EIw;
    p(3) += P * F[0] / t.EIw;
  }
  // homogeneous closure to zero end displacements
  const Eigen::Vector4d c = -t.Hinv * Eigen::Vector4d(0, 0, p(0), p(1));
  Eigen::Vector4d fp(0, 0, t.G * p(1) - t.EIw * p(3), t.EIw * p(2));
  const Eigen::Vector4d f_fixed = t.D * c + fp;
  return -f_fixed;
}

// Hermite cubic on [0, L]: N, N', N'' of (value_i, slope_i, value_j, slope_j)
struct Herm {
  Eigen::RowVector4d n, d1, d2;
};
Herm hermite(double x, double L) {
  const double s = x / L;
  Herm h;
  h.n << 1 - 3 * s * s + 2 * s * s * s, L * (s - 2 * s * s + s * s * s), 3 * s * s - 2 * s * s * s,
      L * (-s * s + s * s * s);
  h.d1 << (-6 * s + 6 * s * s) / L, 1 - 4 * s + 3 * s * s, (6 * s - 6 * s * s) / L, -2 * s + 3 * s * s;
  h.d2 << (-6 + 12 * s) / (L * L), (-4 + 6 * s) / L, (6 - 12 * s) / (L * L), (-2 + 6 * s) / L;
  return h;
}

// Second-order (geometric) matrix of one bar, local 14: the coupling
// My theta v'' + Mz theta w'', the load heights and, when `with_N`, the
// N (v'^2 + w'^2 + r0^2 theta'^2)/2 terms (consistent Hermite).
Mat14 geometric14(const Model& model, const BarElement& el, const WBar& w, const BarSpanLoads& sp,
                  const Vec12& f_end, bool with_N) {
  static const double gx[5] = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831,
                               0.9061798459386640};
  static const double gw[5] = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889,
                               0.4786286704993665, 0.2369268850561891};
  const double L = w.L;
  std::vector<double> xs(5);
  for (int g = 0; g < 5; ++g) xs[static_cast<std::size_t>(g)] = 0.5 * L * (1 + gx[g]);
  BarResult br;
  recover_bar_stations(model, el, sp, f_end, Vec12::Zero(), xs, br);
  const double N = -f_end(0);
  Mat14 kg = Mat14::Zero();
  static const int iv[4] = {1, 5, 7, 11}, iw[4] = {2, 4, 8, 10}, it[4] = {3, 12, 9, 13};
  static const double sw[4] = {1, -1, 1, -1};
  for (int g = 0; g < 5; ++g) {
    const double x = xs[static_cast<std::size_t>(g)], wt = 0.5 * L * gw[g];
    const Herm h = hermite(x, L);
    Vec14 v1 = Vec14::Zero(), v2 = Vec14::Zero(), w1 = Vec14::Zero(), w2 = Vec14::Zero(),
          t0 = Vec14::Zero(), t1 = Vec14::Zero();
    for (int k = 0; k < 4; ++k) {
      v1(iv[k]) = h.d1(k);
      v2(iv[k]) = h.d2(k);
      w1(iw[k]) = sw[k] * h.d1(k);
      w2(iw[k]) = sw[k] * h.d2(k);
      t0(it[k]) = h.n(k);
      t1(it[k]) = h.d1(k);
    }
    const double My = br.stations[static_cast<std::size_t>(g)].My;
    const double Mz = br.stations[static_cast<std::size_t>(g)].Mz;
    kg += wt * (My * (t0 * v2.transpose() + v2 * t0.transpose()) +
                Mz * (t0 * w2.transpose() + w2 * t0.transpose()));
    if (with_N) kg += wt * N * (v1 * v1.transpose() + w1 * w1.transpose() + w.r02 * t1 * t1.transpose());
  }
  // load heights: (q . e) theta^2 / 2
  if (w.off.squaredNorm() > 0) {
    require_linear_span(sp.segs, "warping (load heights)");
    for (const DistSeg& s : sp.segs) {
      const double ea = w.off(0) * s.qa(1) + w.off(1) * s.qa(2), eb = w.off(0) * s.qb(1) + w.off(1) * s.qb(2);
      if (ea == 0 && eb == 0) continue;
      for (int g = 0; g < 5; ++g) {
        const double x = s.a + 0.5 * (s.b - s.a) * (1 + gx[g]), wt = 0.5 * (s.b - s.a) * gw[g];
        const double q = ea + (eb - ea) * (x - s.a) / std::max(s.b - s.a, 1e-300);
        const Herm h = hermite(x, L);
        Vec14 t0 = Vec14::Zero();
        for (int k = 0; k < 4; ++k) t0(it[k]) = h.n(k);
        kg += wt * q * t0 * t0.transpose();
      }
    }
    for (const PointAt& p : sp.points) {
      const double q = w.off(0) * p.F(1) + w.off(1) * p.F(2);
      if (q == 0) continue;
      const Herm h = hermite(p.a, L);
      Vec14 t0 = Vec14::Zero();
      for (int k = 0; k < 4; ++k) t0(it[k]) = h.n(k);
      kg += q * t0 * t0.transpose();
    }
  }
  return kg;
}

// One bar's tangent and fixed-end vector for a state (N, end forces).
struct BarState {
  Mat14 k = Mat14::Zero();   // local
  Vec14 f0 = Vec14::Zero();  // local
};

BarState bar_state(const WSystem& sys, const WBar& w, bool second, double N, const Vec12* f_prev) {
  const Model& model = sys.model;
  const BarElement& el = sys.asm_.elements[w.e];
  const BarSpanLoads& sp = sys.cl.span[w.e];
  const Vec12& f0_lin = sys.cl.f0[w.e];
  BarState s;
  Torsion t;
  t.build(w.G * w.J + (second ? N * w.r02 : 0.0), w.E * w.Iw, w.L, model.bars[w.e].id);
  Vec12 f0b;
  if (second) {
    s.k.topLeftCorner<12, 12>() = el.k_exact_local(N);
    f0b = beam_column_f0(model, el, sp, f0_lin, N);
  } else {
    s.k.topLeftCorner<12, 12>() = el.k_local();
    f0b = f0_lin;
  }
  // torsion block replaces the Saint-Venant entries
  for (int a = 0; a < 14; ++a)
    for (int b : {3, 9}) {
      s.k(a, b) = 0;
      s.k(b, a) = 0;
    }
  const Eigen::Matrix4d kt = t.stiffness();
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) s.k(kTor[a], kTor[b]) += kt(a, b);
  s.f0.head<12>() = f0b;
  const Eigen::Vector4d ft = torsion_f0(w, t, sp, f0_lin);
  for (int a = 0; a < 4; ++a) s.f0(kTor[a]) = ft(a);
  if (second && f_prev) s.k += geometric14(model, el, w, sp, *f_prev, false);
  return s;
}

struct Assembly {
  SpMat K;
  VecX rhs;
};

Assembly assemble(const WSystem& sys, const std::vector<BarState>& st,
                  const std::vector<Mat14>* extra = nullptr) {
  std::vector<Triplet> trips;
  VecX rhs = VecX::Zero(sys.n_free);
  for (std::size_t f = 0; f < sys.free_ix.size(); ++f)
    if (sys.free_ix[f] >= 0) rhs(sys.free_ix[f]) += sys.P_direct(static_cast<Eigen::Index>(f));
  for (std::size_t k = 0; k < sys.bars.size(); ++k) {
    const WBar& w = sys.bars[k];
    Mat14 kl = extra ? (*extra)[k] : st[k].k;
    const Mat14 kg = w.T.transpose() * kl * w.T;
    const Vec14 fg = w.T.transpose() * st[k].f0;
    int g[14];
    for (int a = 0; a < 14; ++a) g[a] = sys.free_ix[static_cast<std::size_t>(sys.dof(w, a))];
    for (int a = 0; a < 14; ++a) {
      if (g[a] < 0) continue;
      if (!extra) rhs(g[a]) += fg(a);
      for (int b = 0; b < 14; ++b)
        if (g[b] >= 0 && kg(a, b) != 0) trips.emplace_back(g[a], g[b], kg(a, b));
    }
  }
  if (!extra)
    for (int n = 0; n < sys.n_nodes; ++n)
      for (int c = 0; c < 6; ++c) {
        const double ks = sys.asm_.spring[static_cast<std::size_t>(sys.asm_.dofs.full(n, c))];
        const int f = sys.free_ix[static_cast<std::size_t>(kW * n + c)];
        if (ks != 0 && f >= 0) trips.emplace_back(f, f, ks);
      }
  Assembly a;
  a.K = SpMat(sys.n_free, sys.n_free);
  a.K.setFromTriplets(trips.begin(), trips.end());
  a.rhs = rhs;
  return a;
}

}  // namespace

WarpingResults solve_warping(const Model& model, Id case_id, const WarpingSpec& spec,
                             const WarpingOptions& opt) {
  require(model.continua.empty(), "solve_warping: continuum elements are not supported (linear statics, modal and dynamics take them)");
  WSystem sys(model, case_id, spec);
  const bool second = opt.theory == WarpingTheory::Second;
  const std::size_t nb = sys.bars.size();
  std::vector<double> N(nb, 0.0);
  std::vector<Vec12> fprev(nb, Vec12::Zero());
  std::vector<BarState> st(nb);
  VecX u = VecX::Zero(sys.n_free);
  WarpingResults out;
  for (int it = 0;; ++it) {
    require(it < (second ? opt.max_iterations : 1) + 1, "warping: second-order iteration did not converge");
    for (std::size_t k = 0; k < nb; ++k)
      st[k] = bar_state(sys, sys.bars[k], second && it > 0, N[k], it > 0 ? &fprev[k] : nullptr);
    Assembly a = assemble(sys, st);
    FactorCache fac;
    fac.factorize(a.K);
    require(fac.inertia().n_neg == 0, second ? "warping: load beyond the critical state (indefinite tangent)"
                                             : "warping: stiffness not positive definite (mechanism)");
    VecX u1 = fac.solve(a.rhs);
    const double kappa = fac.max_abs_pivot() / fac.min_abs_pivot();
    const VecX uf = sys.expand(u1);
    for (std::size_t k = 0; k < nb; ++k) {
      const Vec14 ul = sys.gather(sys.bars[k], uf);
      const Vec14 f = st[k].k * ul - st[k].f0;
      fprev[k] = f.head<12>();
      N[k] = -f(0);
    }
    const double du = (u1 - u).lpNorm<Eigen::Infinity>(), un = std::max(u1.lpNorm<Eigen::Infinity>(), 1e-300);
    u = u1;
    ++out.iterations;
    if (!second) break;
    if (it == 0) continue;  // the first pass is the linear state
    // Round-off floor of a solve with Kt: ~ eps kappa, kappa = the LDLT
    // pivot ratio (near the critical load the tangent is ill-conditioned and
    // the fixed point stalls at that floor: 1e-8 for a member with C_M = 0
    // on 100 elements, kappa = 3e7). Capped at 1e-6: a worse conditioning is a
    // typed non-convergence, not a silently loose answer.
    const double floor = std::min(64.0 * 2.220446049250313e-16 * kappa, 1e-6);
    if (du <= std::max(opt.tol_rel, floor) * un) break;
  }
  // recovery
  const VecX uf = sys.expand(u);
  out.node_disp.assign(static_cast<std::size_t>(sys.n_nodes), Vec7::Zero());
  for (int n = 0; n < sys.n_nodes; ++n)
    for (int c = 0; c < kW; ++c) out.node_disp[static_cast<std::size_t>(n)](c) = uf(kW * n + c);
  VecX fint = VecX::Zero(kW * sys.n_nodes);
  for (std::size_t k = 0; k < nb; ++k) {
    const WBar& w = sys.bars[k];
    const Vec14 ul = sys.gather(w, uf);
    const Vec14 f = st[k].k * ul - st[k].f0;
    const Vec14 fg = w.T.transpose() * f;
    for (int a = 0; a < 14; ++a) fint(sys.dof(w, a)) += fg(a);
    WarpingBarResult r;
    r.bar = model.bars[w.e].id;
    r.end_forces = f;
    r.N = -f(0);
    const double r3 = second ? r.N * w.r02 : 0.0;
    r.MT = {-f(3), f(9)};
    r.Mw = {f(12), -f(13)};
    for (int s = 0; s < 2; ++s) {
      const double psi = ul(12 + s);
      r.MT1[static_cast<std::size_t>(s)] = w.G * w.J * psi;
      r.MT3[static_cast<std::size_t>(s)] = r3 * psi;
      r.MT2[static_cast<std::size_t>(s)] = r.MT[static_cast<std::size_t>(s)] - r.MT1[static_cast<std::size_t>(s)] -
                                           r.MT3[static_cast<std::size_t>(s)];
    }
    out.bars.push_back(r);
  }
  for (int n = 0; n < sys.n_nodes; ++n)
    for (int c = 0; c < 6; ++c)
      fint(kW * n + c) += sys.asm_.spring[static_cast<std::size_t>(sys.asm_.dofs.full(n, c))] * uf(kW * n + c);
  std::vector<Id> seen;
  for (const Support& s : model.supports) {
    if (s.implicit || std::find(seen.begin(), seen.end(), s.node) != seen.end()) continue;
    seen.push_back(s.node);
    const int n = model.node_index(s.node);
    Vec7 R;
    for (int c = 0; c < kW; ++c) R(c) = fint(kW * n + c) - sys.P_direct(kW * n + c);
    out.reactions.emplace_back(s.node, R);
  }
  return out;
}

WarpingBucklingResults warping_buckling(const Model& model, Id case_id, const WarpingSpec& spec,
                                        int n_modes) {
  require(model.continua.empty(), "warping_buckling: continuum elements are not supported (linear statics, modal and dynamics take them)");
  WarpingResults lin = solve_warping(model, case_id, spec);
  WSystem sys(model, case_id, spec);
  const std::size_t nb = sys.bars.size();
  std::vector<BarState> st(nb);
  std::vector<Mat14> kg(nb);
  for (std::size_t k = 0; k < nb; ++k) {
    st[k] = bar_state(sys, sys.bars[k], false, 0.0, nullptr);
    kg[k] = geometric14(model, sys.asm_.elements[sys.bars[k].e], sys.bars[k], sys.cl.span[sys.bars[k].e],
                        lin.bars[k].end_forces.head<12>(), true);
  }
  Assembly a = assemble(sys, st);
  Assembly g = assemble(sys, st, &kg);
  std::vector<PencilPair> pairs;
  if (sys.n_free <= 600) {
    pairs = dense_pencil(a.K, g.K);
  } else {
    FactorCache fac;
    fac.factorize(a.K);
    fac.require_spd();
    pairs = k_inner_lanczos(fac, a.K, g.K, std::max(2 * n_modes, n_modes + 8));
  }
  // Kg x = mu K x, lambda = -1/mu for mu < 0 (most negative first)
  std::vector<std::pair<double, VecX>> lam;
  for (const PencilPair& p : pairs)
    if (p.mu < 0) lam.emplace_back(-1.0 / p.mu, p.x);
  std::sort(lam.begin(), lam.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
  WarpingBucklingResults r;
  for (std::size_t k = 0; k < lam.size() && static_cast<int>(k) < n_modes; ++k) {
    r.load_factors.push_back(lam[k].first);
    const VecX full = sys.expand(lam[k].second);
    double crest = 0;
    for (int i = 0; i < full.size(); ++i)
      if (std::fabs(full(i)) > std::fabs(crest)) crest = full(i);
    std::vector<Vec7> mode(static_cast<std::size_t>(sys.n_nodes));
    for (int n = 0; n < sys.n_nodes; ++n)
      for (int c = 0; c < kW; ++c) mode[static_cast<std::size_t>(n)](c) = full(kW * n + c) / crest;
    r.modes.push_back(mode);
  }
  return r;
}

}  // namespace keta

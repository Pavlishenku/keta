// Concentrated plastic hinges (moment and N-M interaction), exact
// event-to-event pushover — see plastic_hinge.h.
#include "keta/plastic_hinge.h"

#include <algorithm>
#include <cmath>
#include <map>

#include "keta/assembly.h"
#include "keta/linsolve.h"

namespace keta {

// ------------------------------------------------------- interaction surface

InteractionSurface InteractionSurface::from_vertices(const std::vector<std::pair<double, double>>& v) {
  require(v.size() >= 3, "interaction surface: at least three vertices");
  InteractionSurface s;
  const std::size_t n = v.size();
  double orient = 0;
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    orient += a.first * b.second - b.first * a.second;
  }
  require(orient != 0, "interaction surface: degenerate polygon");
  for (std::size_t i = 0; i < n; ++i) {
    const auto& a = v[i];
    const auto& b = v[(i + 1) % n];
    // outward normal of edge a -> b for a counter-clockwise polygon
    double aN = b.second - a.second, aM = -(b.first - a.first);
    if (orient < 0) {
      aN = -aN;
      aM = -aM;
    }
    const double r = aN * a.first + aM * a.second;
    require(r > 0, "interaction surface: the origin must lie strictly inside");
    for (std::size_t k = 0; k < n; ++k)
      require(aN * v[k].first + aM * v[k].second <= r * (1 + 1e-12), "interaction surface: polygon not convex");
    s.facets.emplace_back(aN, aM, r);
  }
  return s;
}

InteractionSurface InteractionSurface::diamond(double Nt, double Nc, double Mpl) {
  require(Nt > 0 && Nc > 0 && Mpl > 0, "interaction surface: capacities must be positive");
  return from_vertices({{Nt, 0}, {0, Mpl}, {-Nc, 0}, {0, -Mpl}});
}

namespace {

// A facet counts as reached when f >= -kTie r: the hinge forces come out of a
// global solve whose relative accuracy is ~ cond(K) eps (a stiff elastic
// hinge, My / k1 ~ 1e-6 rad, makes cond(K) ~ 1e5 .. 1e7), so a reach within
// 1e-9 of the facet's capacity is a tie, not a missed event.
constexpr double kTie = 1e-9;

// One hinge in its spring coordinates: m deformations (rotation, or axial
// slip + rotation) at local DOFs q, elastic Ke, hardening H, facets n . s <= r
struct HingeModel {
  int m = 1;
  int q[2] = {0, 0};
  Eigen::Matrix2d Ke = Eigen::Matrix2d::Zero(), H = Eigen::Matrix2d::Zero();
  std::vector<Eigen::Vector2d> n;  // facet normals (first m components used)
  std::vector<double> r;
  double sgnN = 1;  // N = sgnN * s_axial (tension positive)
};

struct HingeVar {
  Eigen::Vector2d e = Eigen::Vector2d::Zero(), ep = Eigen::Vector2d::Zero(), beta = Eigen::Vector2d::Zero();
  std::vector<int> active;
};

Eigen::MatrixXd blk(const Eigen::Matrix2d& A, int m) { return A.topLeftCorner(m, m); }
Eigen::VectorXd vec(const Eigen::Vector2d& a, int m) { return a.head(m); }

// active-set matrices: N (k x m), G = N (Ke + H) N^T
struct ActiveOps {
  Eigen::MatrixXd N, Ginv;
  Eigen::MatrixXd Kt;
};
ActiveOps active_ops(const HingeModel& h, const std::vector<int>& act) {
  ActiveOps o;
  const Eigen::MatrixXd Ke = blk(h.Ke, h.m);
  if (act.empty()) {
    o.Kt = Ke;
    return o;
  }
  require(static_cast<int>(act.size()) <= h.m, "plastic hinge: more active facets than deformations");
  o.N.resize(static_cast<Eigen::Index>(act.size()), h.m);
  for (std::size_t k = 0; k < act.size(); ++k) o.N.row(static_cast<Eigen::Index>(k)) = vec(h.n[static_cast<std::size_t>(act[k])], h.m).transpose();
  const Eigen::MatrixXd G = o.N * (Ke + blk(h.H, h.m)) * o.N.transpose();
  const Eigen::FullPivLU<Eigen::MatrixXd> lu(G);
  require(lu.isInvertible(), "plastic hinge: dependent active facets");
  o.Ginv = lu.inverse();
  o.Kt = Ke - Ke * o.N.transpose() * o.Ginv * o.N * Ke;
  return o;
}

struct HingedBar {
  std::size_t e = 0;
  std::vector<int> hinges;   // hinge indices on this bar
  std::vector<int> qs;       // all hinged local DOFs, in hinge order
  int gdof[12];
  Mat12 kc_minus_k = Mat12::Zero();
  MatX G;                    // hinge deformations = G * u_nodal (+ load offset, rows: qs order)
  MatX Kii_inv, Kbi;         // internal block inverse, bar-end / internal coupling (span loads)
};

void condense(const BarElement& el, HingedBar& hb, const std::vector<HingeModel>& hm,
              const std::vector<Eigen::MatrixXd>& Kt) {
  const int m = static_cast<int>(hb.qs.size());
  const Mat12& k = el.k_local();
  MatX Ka = MatX::Zero(12 + m, 12 + m);
  std::vector<int> map(12);
  for (int a = 0; a < 12; ++a) map[static_cast<std::size_t>(a)] = a;
  for (int h = 0; h < m; ++h) map[static_cast<std::size_t>(hb.qs[static_cast<std::size_t>(h)])] = 12 + h;
  for (int a = 0; a < 12; ++a)
    for (int b = 0; b < 12; ++b) Ka(map[static_cast<std::size_t>(a)], map[static_cast<std::size_t>(b)]) += k(a, b);
  // springs: 1/2 (u_q - t)^T Kt (u_q - t), one block per hinge
  int off = 0;
  for (int hi : hb.hinges) {
    const HingeModel& h = hm[static_cast<std::size_t>(hi)];
    const Eigen::MatrixXd& K = Kt[static_cast<std::size_t>(hi)];
    for (int a = 0; a < h.m; ++a)
      for (int b = 0; b < h.m; ++b) {
        const int qa = h.q[a], qb = h.q[b], ta = 12 + off + a, tb = 12 + off + b;
        Ka(qa, qb) += K(a, b);
        Ka(ta, tb) += K(a, b);
        Ka(qa, tb) -= K(a, b);
        Ka(ta, qb) -= K(a, b);
      }
    off += h.m;
  }
  const MatX Kii = Ka.bottomRightCorner(m, m), Kib = Ka.bottomLeftCorner(m, 12);
  const Eigen::FullPivLU<MatX> f(Kii);
  require(f.isInvertible(), "plastic hinge: singular hinge condensation (a fully plastic hinge with no hardening "
                            "on an otherwise free bar end)");
  const MatX X = f.solve(Kib);  // t = -X u (+ Kii^-1 f_i under span loads)
  hb.Kii_inv = f.inverse();
  hb.Kbi = Ka.topRightCorner(12, m);
  const Mat12 Kc = Ka.topLeftCorner(12, 12) - Ka.topRightCorner(12, m) * X;
  hb.kc_minus_k = el.T().transpose() * (Kc - k) * el.T();
  MatX Gl = X;
  for (int h = 0; h < m; ++h) Gl(h, hb.qs[static_cast<std::size_t>(h)]) += 1.0;
  hb.G = Gl * el.T();
}

}  // namespace

HingePushoverResult hinge_pushover(const Model& model, Id case_id, const std::vector<HingeSpec>& hinges,
                                   const std::vector<double>& lambdas, const HingePushoverOptions& opt) {
  require(model.continua.empty(), "hinge_pushover: continuum elements are not supported (linear statics, modal and dynamics take them)");
  for (std::size_t k = 0; k < lambdas.size(); ++k)
    require(lambdas[k] >= 0 && (k == 0 || lambdas[k] >= lambdas[k - 1]),
            "hinge pushover: load factors must be non-negative and non-decreasing");
  const Assembled asm_ = assemble_elements(model);
  const DofMap& d = asm_.dofs;
  auto load_case = [&](Id cid) {
    CaseLoads cl = lower_case(model, asm_, model.load_case(cid));
    require(cl.prescribed.empty(), "hinge pushover: prescribed displacements are not supported");
    return cl;
  };
  const CaseLoads clP = load_case(case_id);
  const CaseLoads clPerm = opt.permanent_case ? load_case(*opt.permanent_case) : empty_case_loads(asm_);

  // hinge models
  std::vector<HingeModel> hm(hinges.size());
  std::map<std::size_t, std::size_t> bar_slot;
  std::vector<HingedBar> hbars;
  std::vector<std::size_t> hinge_bar(hinges.size());
  std::vector<int> hinge_row(hinges.size());  // first row of the hinge in its bar's G
  for (std::size_t i = 0; i < hinges.size(); ++i) {
    const std::size_t e = static_cast<std::size_t>(model.bar_index(hinges[i].bar));
    for (const CaseLoads* c : {&clP, &clPerm}) {
      const BarSpanLoads& sp = c->span[e];
      require(sp.kinks.empty() && sp.eps0 == 0 && sp.kappa_y == 0 && sp.kappa_z == 0 && sp.deps0 == 0 &&
                  sp.dkappa_y == 0 && sp.dkappa_z == 0 && sp.kappa_x == 0 && sp.dkappa_x == 0 && !has_parabolic_strain(sp),
              "hinge pushover: imposed strains / dislocations on a hinged bar are not supported");
    }
  }
  for (std::size_t i = 0; i < hinges.size(); ++i) {
    const HingeSpec& s = hinges[i];
    require(s.end == 0 || s.end == 1, "hinge pushover: end must be 0 or 1");
    require(s.axis >= 0 && s.axis <= 2, "hinge pushover: axis must be 0, 1 or 2");
    HingeModel& h = hm[i];
    const int qrot = 6 * s.end + 3 + s.axis;
    if (s.nm) {
      const NMHinge& nm = *s.nm;
      require(nm.kN > 0 && nm.kM > 0, "hinge pushover: N-M hinge stiffnesses must be positive");
      require(!nm.surface.facets.empty(), "hinge pushover: N-M hinge without an interaction surface");
      require((nm.H - nm.H.transpose()).norm() <= 1e-12 * nm.H.norm() &&
                  nm.H.selfadjointView<Eigen::Lower>().eigenvalues().minCoeff() >= 0,
              "hinge pushover: hardening must be symmetric positive semi-definite");
      h.m = 2;
      h.q[0] = 6 * s.end;  // local axial
      h.q[1] = qrot;
      h.sgnN = s.end == 1 ? 1.0 : -1.0;
      h.Ke << nm.kN, 0, 0, nm.kM;
      const Eigen::Matrix2d D = Eigen::Vector2d(h.sgnN, 1.0).asDiagonal();
      h.H = D * nm.H * D;
      for (const Vec3& f : nm.surface.facets) {
        h.n.push_back(D * Eigen::Vector2d(f(0), f(1)));
        h.r.push_back(f(2));
      }
    } else {
      const BilinearHinge& L = s.law;
      require(L.k1 > 0 && L.My > 0 && L.k2 >= 0 && L.k2 < L.k1,
              "hinge pushover: bilinear law needs k1 > 0, My > 0, 0 <= k2 < k1");
      h.m = 1;
      h.q[0] = qrot;
      h.Ke(0, 0) = L.k1;
      h.H(0, 0) = L.k2 == 0 ? 0.0 : L.k1 * L.k2 / (L.k1 - L.k2);  // plastic slope k1 H / (k1 + H) = k2
      h.n = {Eigen::Vector2d(1, 0), Eigen::Vector2d(-1, 0)};
      h.r = {L.My, L.My};
    }
    const std::size_t e = static_cast<std::size_t>(model.bar_index(s.bar));
    const Bar& b = model.bars[e];
    require(b.kind != BarKind::Cable, "hinge pushover: a hinge on a cable");
    require(model.tapered_of(b.id) == nullptr, "hinge pushover: a hinge on a tapered bar");
    for (int a = 0; a < h.m; ++a)
      require(!b.releases[static_cast<std::size_t>(h.q[a])], "hinge pushover: a hinge on a released DOF");
    auto it = bar_slot.find(e);
    if (it == bar_slot.end()) {
      it = bar_slot.emplace(e, hbars.size()).first;
      HingedBar hb;
      hb.e = e;
      const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      for (int c = 0; c < 6; ++c) {
        hb.gdof[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
        hb.gdof[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
      }
      hbars.push_back(hb);
    }
    HingedBar& hb = hbars[it->second];
    for (int a = 0; a < h.m; ++a)
      require(std::find(hb.qs.begin(), hb.qs.end(), h.q[a]) == hb.qs.end(), "hinge pushover: two hinges on one DOF");
    hinge_bar[i] = it->second;
    hinge_row[i] = static_cast<int>(hb.qs.size());
    hb.hinges.push_back(static_cast<int>(i));
    for (int a = 0; a < h.m; ++a) hb.qs.push_back(h.q[a]);
  }

  if (opt.second_order) {
    // ======================= second order (Theory II) =======================
    require(opt.substeps >= 1, "hinge pushover: substeps must be >= 1");
    std::vector<int> bar_hb(model.bars.size(), -1);
    for (std::size_t b = 0; b < hbars.size(); ++b) bar_hb[hbars[b].e] = static_cast<int>(b);
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      const Bar& b = model.bars[e];
      require(b.kind == BarKind::Frame, "hinge pushover (second order): truss / cable bars are not supported");
      require(model.tapered_of(b.id) == nullptr, "hinge pushover (second order): tapered bars are not supported");
      if (bar_hb[e] >= 0) require(!b.soil.any(), "hinge pushover (second order): a Winkler bed on a hinged bar");
    }
    // nodal-only parts of the phases (the bars' equivalent loads are
    // re-lowered with f0(N) at every state)
    auto nodal_only = [&](const CaseLoads& cl) {
      VecX F = cl.rhs;
      for (std::size_t e = 0; e < model.bars.size(); ++e) {
        const Vec12 g = asm_.elements[e].T().transpose() * asm_.elements[e].condense_f0(cl.f0[e]);
        const Bar& b = model.bars[e];
        const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
        for (int c = 0; c < 6; ++c) {
          const int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))], fj = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
          if (fi >= 0) F(fi) -= g(c);
          if (fj >= 0) F(fj) -= g(6 + c);
        }
      }
      return F;
    };
    const VecX FnP = nodal_only(clP), FnPerm = nodal_only(clPerm);
    std::vector<Triplet> spring_trips;
    for (std::size_t f = 0; f < asm_.spring.size(); ++f)
      if (asm_.spring[f] != 0.0 && d.free_ix[f] >= 0) spring_trips.emplace_back(d.free_ix[f], d.free_ix[f], asm_.spring[f]);
    for (const auto& [ni, K6] : asm_.spring6)
      for (int a = 0; a < 6; ++a)
        for (int c = 0; c < 6; ++c) {
          const int fa = d.free_ix[static_cast<std::size_t>(d.full(ni, a))], fc2 = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
          if (fa >= 0 && fc2 >= 0 && K6(a, c) != 0.0) spring_trips.emplace_back(fa, fc2, K6(a, c));
        }

    // segment base (constant active set): hinge vars at the segment start
    struct Base {
      std::vector<HingeVar> v;
      std::vector<Eigen::MatrixXd> Kt;
      std::vector<ActiveOps> ops;
    };
    Base base;
    base.v.assign(hinges.size(), HingeVar{});
    for (std::size_t i = 0; i < hinges.size(); ++i) {
      base.ops.push_back(active_ops(hm[i], {}));
      base.Kt.push_back(base.ops.back().Kt);
    }
    auto sbase = [&](std::size_t i) {
      const HingeModel& h = hm[i];
      Eigen::VectorXd s = blk(h.Ke, h.m) * (base.v[i].e - base.v[i].ep).head(h.m);
      return s;
    };

    struct Sol {
      VecX u;
      std::vector<double> N;
      std::vector<Eigen::Vector2d> e;  // hinge deformations
    };
    Sol warm;
    warm.u = VecX::Zero(d.n_free);
    warm.N.assign(model.bars.size(), 0.0);
    warm.e.assign(hinges.size(), Eigen::Vector2d::Zero());

    // Theory-II equilibrium at (mu, lam) for the current segment base
    auto solve_at = [&](double mu, double lam) {
      Sol S = warm;
      for (int it = 0;; ++it) {
        if (it >= 500) throw Error("hinge pushover (second order): no Theory-II equilibrium at load factor " +
                                   std::to_string(lam) + " (limit point)");
        std::vector<Triplet> trips = spring_trips;
        VecX F = mu * FnPerm + lam * FnP;
        struct BarEval {
          Mat12 kraw, kr;
          Vec12 f0c;
          MatX X, Kii_inv;
          VecX Fa;
          int m = 0;
        };
        std::vector<BarEval> be(model.bars.size());
        for (std::size_t e = 0; e < model.bars.size(); ++e) {
          const BarElement& el = asm_.elements[e];
          BarEval& B = be[e];
          const double Ne = S.N[e];
          B.kraw = el.k_exact_local(Ne);
          Vec12 f0 = Vec12::Zero();
          if (mu != 0 && clPerm.f0[e].cwiseAbs().maxCoeff() > 0)
            f0 += mu * beam_column_f0(model, el, clPerm.span[e], clPerm.f0[e], Ne);
          if (lam != 0 && clP.f0[e].cwiseAbs().maxCoeff() > 0)
            f0 += lam * beam_column_f0(model, el, clP.span[e], clP.f0[e], Ne);
          B.kr = el.condense_matrix(B.kraw);
          B.f0c = el.condense_f0_with(B.kraw, f0);
          const Bar& bar = model.bars[e];
          const int ni = model.node_index(bar.ni), nj = model.node_index(bar.nj);
          int g[12];
          for (int c = 0; c < 6; ++c) {
            g[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
            g[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
          }
          Mat12 Kg;
          Vec12 feq;
          if (bar_hb[e] < 0) {
            Kg = el.ke_exact(Ne);
            feq = el.T().transpose() * B.f0c;
          } else {
            const HingedBar& hb = hbars[static_cast<std::size_t>(bar_hb[e])];
            const int m = static_cast<int>(hb.qs.size());
            B.m = m;
            MatX Ka = MatX::Zero(12 + m, 12 + m);
            std::vector<int> map(12);
            for (int a = 0; a < 12; ++a) map[static_cast<std::size_t>(a)] = a;
            for (int h = 0; h < m; ++h) map[static_cast<std::size_t>(hb.qs[static_cast<std::size_t>(h)])] = 12 + h;
            for (int a = 0; a < 12; ++a)
              for (int b = 0; b < 12; ++b) Ka(map[static_cast<std::size_t>(a)], map[static_cast<std::size_t>(b)]) += B.kr(a, b);
            B.Fa = VecX::Zero(12 + m);
            for (int a = 0; a < 12; ++a) B.Fa(map[static_cast<std::size_t>(a)]) += B.f0c(a);
            int off = 0;
            for (int hi : hb.hinges) {
              const std::size_t i = static_cast<std::size_t>(hi);
              const HingeModel& h = hm[i];
              const Eigen::MatrixXd& K = base.Kt[i];
              // s = Kt e + c, c = s_base - Kt e_base
              const Eigen::VectorXd c = sbase(i) - K * base.v[i].e.head(h.m);
              for (int a = 0; a < h.m; ++a) {
                for (int b = 0; b < h.m; ++b) {
                  const int qa = h.q[a], qb = h.q[b], ta = 12 + off + a, tb = 12 + off + b;
                  Ka(qa, qb) += K(a, b);
                  Ka(ta, tb) += K(a, b);
                  Ka(qa, tb) -= K(a, b);
                  Ka(ta, qb) -= K(a, b);
                }
                B.Fa(h.q[a]) -= c(a);
                B.Fa(12 + off + a) += c(a);
              }
              off += h.m;
            }
            const MatX Kii = Ka.bottomRightCorner(m, m), Kib = Ka.bottomLeftCorner(m, 12), Kbi = Ka.topRightCorner(12, m);
            const Eigen::FullPivLU<MatX> lu(Kii);
            if (!lu.isInvertible()) throw Error("hinge pushover (second order): singular hinge condensation (mechanism)");
            B.Kii_inv = lu.inverse();
            B.X = B.Kii_inv * Kib;
            const Mat12 Kc = Ka.topLeftCorner(12, 12) - Kbi * B.X;
            const Vec12 fl = B.Fa.head(12) - Kbi * (B.Kii_inv * B.Fa.tail(m));
            Kg = el.T().transpose() * Kc * el.T();
            feq = el.T().transpose() * fl;
          }
          for (int a = 0; a < 12; ++a) {
            if (g[a] < 0) continue;
            F(g[a]) += feq(a);
            for (int b = 0; b < 12; ++b)
              if (g[b] >= 0 && Kg(a, b) != 0.0) trips.emplace_back(g[a], g[b], Kg(a, b));
          }
        }
        SpMat K(d.n_free, d.n_free);
        K.setFromTriplets(trips.begin(), trips.end());
        FactorCache f2;
        try {
          f2.factorize(K);
          f2.require_spd();
        } catch (const Error&) {
          throw Error("hinge pushover (second order): K(N) singular or indefinite at load factor " +
                      std::to_string(lam) + " (limit point or mechanism)");
        }
        S.u = f2.solve(F);
        // new axial forces and the hinge deformations
        double dN = 0, Nmax = 1.0;
        for (std::size_t e = 0; e < model.bars.size(); ++e) {
          const BarElement& el = asm_.elements[e];
          const BarEval& B = be[e];
          const Bar& bar = model.bars[e];
          const int ni = model.node_index(bar.ni), nj = model.node_index(bar.nj);
          Vec12 un;
          for (int c = 0; c < 6; ++c) {
            const int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))], fj = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
            un(c) = fi >= 0 ? S.u(fi) : 0.0;
            un(6 + c) = fj >= 0 ? S.u(fj) : 0.0;
          }
          double Nn;
          if (bar_hb[e] < 0) {
            const Vec12 ul = el.local_displacements_with(B.kraw, un, B.f0c);
            Nn = (B.kraw * ul - B.f0c)(6);
          } else {
            const HingedBar& hb = hbars[static_cast<std::size_t>(bar_hb[e])];
            const Vec12 ul = el.T() * un;
            const VecX t = B.Kii_inv * B.Fa.tail(B.m) - B.X * ul;
            Vec12 v = ul;
            for (int h = 0; h < B.m; ++h) v(hb.qs[static_cast<std::size_t>(h)]) = t(h);
            Nn = (B.kr * v - B.f0c)(6);
            for (int hi : hb.hinges) {
              const std::size_t i = static_cast<std::size_t>(hi);
              for (int a = 0; a < hm[i].m; ++a) {
                const int row = hinge_row[i] + a;
                S.e[i](a) = ul(hb.qs[static_cast<std::size_t>(row)]) - t(row);
              }
            }
          }
          dN = std::max(dN, std::fabs(Nn - S.N[e]));
          Nmax = std::max(Nmax, std::fabs(Nn));
          S.N[e] = Nn;
        }
        if (dN <= 1e-13 * Nmax) break;
      }
      return S;
    };
    // hinge quantities at a solution, from the segment base
    struct HQ {
      std::vector<Eigen::VectorXd> dlam;  // multiplier increments of the active facets
      std::vector<std::vector<double>> f; // all facets
    };
    auto hinge_q = [&](const Sol& S) {
      HQ q;
      q.dlam.resize(hinges.size());
      q.f.resize(hinges.size());
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        const HingeModel& h = hm[i];
        const Eigen::VectorXd de = (S.e[i] - base.v[i].e).head(h.m);
        Eigen::VectorXd dep = Eigen::VectorXd::Zero(h.m);
        if (!base.v[i].active.empty()) {
          q.dlam[i] = base.ops[i].Ginv * base.ops[i].N * blk(h.Ke, h.m) * de;
          dep = base.ops[i].N.transpose() * q.dlam[i];
        }
        const Eigen::VectorXd s = sbase(i) + base.Kt[i] * de;
        const Eigen::VectorXd beta = base.v[i].beta.head(h.m) + blk(h.H, h.m) * dep;
        for (std::size_t j = 0; j < h.n.size(); ++j) q.f[i].push_back(h.n[j].head(h.m).dot(s - beta) - h.r[j]);
      }
      return q;
    };
    auto is_active = [&](std::size_t i, int j) {
      return std::find(base.v[i].active.begin(), base.v[i].active.end(), j) != base.v[i].active.end();
    };
    // move the segment base to the solution S (commit the hinge states)
    auto rebase = [&](const Sol& S) {
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        const HingeModel& h = hm[i];
        const Eigen::VectorXd de = (S.e[i] - base.v[i].e).head(h.m);
        if (!base.v[i].active.empty()) {
          const Eigen::VectorXd dl = base.ops[i].Ginv * base.ops[i].N * blk(h.Ke, h.m) * de;
          const Eigen::VectorXd dep = base.ops[i].N.transpose() * dl;
          base.v[i].ep.head(h.m) += dep;
          base.v[i].beta.head(h.m) += blk(h.H, h.m) * dep;
        }
        base.v[i].e = S.e[i];
      }
      warm = S;
    };
    auto set_act = [&](std::size_t i, std::vector<int> act) {
      std::sort(act.begin(), act.end());
      base.v[i].active = act;
      base.ops[i] = active_ops(hm[i], act);
      base.Kt[i] = base.ops[i].Kt;
    };
    // consistent active set at the base, probed by a tiny step along the path
    auto make_consistent = [&](double mu, double lam, double dmu, double dlam) {
      for (std::size_t pass = 0;; ++pass) {
        require(pass < 4 * hinges.size() + 8, "hinge pushover (second order): no consistent active set");
        const Sol S = solve_at(mu + dmu, lam + dlam);
        const HQ q = hinge_q(S);
        bool changed = false;
        for (std::size_t i = 0; i < hinges.size() && !changed; ++i) {
          for (Eigen::Index k = 0; k < q.dlam[i].size() && !changed; ++k)
            if (q.dlam[i](k) < 0) {
              std::vector<int> act = base.v[i].active;
              act.erase(act.begin() + k);
              set_act(i, act);
              changed = true;
            }
          if (changed) break;
          const HingeModel& h = hm[i];
          for (int j = 0; j < static_cast<int>(h.n.size()) && !changed; ++j)
            if (!is_active(i, j) && q.f[i][static_cast<std::size_t>(j)] > 0 &&
                static_cast<int>(base.v[i].active.size()) < h.m) {
              std::vector<int> act = base.v[i].active;
              act.push_back(j);
              set_act(i, act);
              changed = true;
            }
        }
        if (!changed) return;
      }
    };
    // first violation in [a, b] of the path (phase direction), or b
    auto advance2 = [&](bool perm_phase, double& par, double target, std::vector<double>& events) {
      auto at = [&](double x) { return perm_phase ? solve_at(x, 0.0) : solve_at(1.0, x); };
      const double h0 = (target - par) / opt.substeps;
      while (par < target) {
        const double probe = 1e-7 * std::max(std::fabs(h0), 1e-300);
        if (perm_phase)
          make_consistent(par, 0.0, probe, 0.0);
        else
          make_consistent(1.0, par, 0.0, probe);
        const double a = par, b = std::min(target, par + h0);
        const Sol Sb = at(b);
        const HQ qb = hinge_q(Sb);
        // the multipliers' slope at b (a step back along the path)
        const HQ qbm = hinge_q(at(b - std::min(probe, 0.5 * (b - a))));
        // (i) an inactive facet reached; (ii) an active multiplier decreasing
        double xev = b;
        bool ev = false;
        for (std::size_t i = 0; i < hinges.size(); ++i) {
          const HingeModel& h = hm[i];
          for (int j = 0; j < static_cast<int>(h.n.size()); ++j) {
            if (is_active(i, j) || qb.f[i][static_cast<std::size_t>(j)] < -kTie * h.r[static_cast<std::size_t>(j)])
              continue;
            double lo = a, hi = b;
            for (int it = 0; it < 200 && hi - lo > 4e-16 * std::max(1.0, std::fabs(hi)); ++it) {
              const double mid = 0.5 * (lo + hi);
              (hinge_q(at(mid)).f[i][static_cast<std::size_t>(j)] >= 0 ? hi : lo) = mid;  // the root itself
            }
            if (hi < xev || !ev) xev = std::min(xev, hi);
            ev = true;
          }
          for (Eigen::Index k = 0; k < qb.dlam[i].size(); ++k) {
            if (qb.dlam[i](k) >= qbm.dlam[i](k)) continue;  // still loading at b
            // golden section for the maximum of the multiplier in [a, b]
            double l = a, r = b;
            const double g = 0.5 * (std::sqrt(5.0) - 1);
            double c = r - g * (r - l), dd = l + g * (r - l);
            double fc = hinge_q(at(c)).dlam[i](k), fd = hinge_q(at(dd)).dlam[i](k);
            for (int it = 0; it < 200 && r - l > 4e-16 * std::max(1.0, std::fabs(r)); ++it) {
              if (fc > fd) {
                r = dd;
                dd = c;
                fd = fc;
                c = r - g * (r - l);
                fc = hinge_q(at(c)).dlam[i](k);
              } else {
                l = c;
                c = dd;
                fc = fd;
                dd = l + g * (r - l);
                fd = hinge_q(at(dd)).dlam[i](k);
              }
            }
            xev = std::min(xev, 0.5 * (l + r));
            ev = true;
          }
        }
        if (!ev) {
          rebase(Sb);
          par = b;
          continue;
        }
        const Sol Se = at(xev);
        rebase(Se);
        par = xev;
        events.push_back(xev);
        // an event located (to bisection precision) at the requested load
        // factor: the state there is the target's
        if (target - par <= 1e-12 * std::max(1.0, std::fabs(target))) par = target;
        // facets reached join (the next consistency pass settles the rest)
        const HQ qe = hinge_q(Se);
        for (std::size_t i = 0; i < hinges.size(); ++i) {
          const HingeModel& h = hm[i];
          std::vector<int> act = base.v[i].active;
          for (int j = 0; j < static_cast<int>(h.n.size()); ++j)
            if (!is_active(i, j) && qe.f[i][static_cast<std::size_t>(j)] >= -kTie * h.r[static_cast<std::size_t>(j)] &&
                static_cast<int>(act.size()) < h.m)
              act.push_back(j);
          if (act.size() != base.v[i].active.size()) set_act(i, act);
        }
      }
    };
    auto snapshot2 = [&](double lam) {
      HingePushoverPoint pt;
      pt.lambda = lam;
      pt.node_disp.assign(model.nodes.size(), Vec6::Zero());
      for (std::size_t n = 0; n < model.nodes.size(); ++n)
        for (int c = 0; c < 6; ++c) {
          const int fi = d.free_ix[static_cast<std::size_t>(d.full(static_cast<int>(n), c))];
          if (fi >= 0) pt.node_disp[n](c) = warm.u(fi);
        }
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        const HingeModel& h = hm[i];
        const Eigen::VectorXd s = sbase(i);
        HingeState hs;
        hs.plastic = !base.v[i].active.empty();
        if (h.m == 1) {
          hs.theta = base.v[i].e(0);
          hs.theta_p = base.v[i].ep(0);
          hs.M = s(0);
        } else {
          hs.delta = h.sgnN * base.v[i].e(0);
          hs.delta_p = h.sgnN * base.v[i].ep(0);
          hs.N = h.sgnN * s(0);
          hs.theta = base.v[i].e(1);
          hs.theta_p = base.v[i].ep(1);
          hs.M = s(1);
        }
        pt.hinges.push_back(hs);
      }
      return pt;
    };
    HingePushoverResult out2;
    if (opt.permanent_case) {
      double mu = 0;
      advance2(true, mu, 1.0, out2.permanent_events);
    }
    out2.permanent = snapshot2(0.0);
    double lam2 = 0;
    for (double target : lambdas) {
      if (target > lam2) advance2(false, lam2, target, out2.events);
      out2.points.push_back(snapshot2(target));
    }
    return out2;
  }

  const SpMat K0 = assemble_stiffness(model, asm_);
  std::vector<HingeVar> st(hinges.size());
  std::vector<ActiveOps> ops(hinges.size());
  std::vector<Eigen::MatrixXd> Kt(hinges.size());
  std::vector<char> dirty(hbars.size(), 1);
  for (std::size_t i = 0; i < hinges.size(); ++i) {
    ops[i] = active_ops(hm[i], {});
    Kt[i] = ops[i].Kt;
  }
  FactorCache fc;
  bool need_factor = true;
  auto set_active = [&](std::size_t i, std::vector<int> act) {
    std::sort(act.begin(), act.end());
    st[i].active = act;
    ops[i] = active_ops(hm[i], act);
    Kt[i] = ops[i].Kt;
    dirty[hinge_bar[i]] = 1;
    need_factor = true;
  };
  auto factor = [&]() {
    for (std::size_t b = 0; b < hbars.size(); ++b)
      if (dirty[b]) {
        condense(asm_.elements[hbars[b].e], hbars[b], hm, Kt);
        dirty[b] = 0;
      }
    std::vector<Triplet> trips;
    trips.reserve(hbars.size() * 144);
    for (const HingedBar& hb : hbars)
      for (int a = 0; a < 12; ++a)
        for (int c = 0; c < 12; ++c)
          if (hb.gdof[a] >= 0 && hb.gdof[c] >= 0) trips.emplace_back(hb.gdof[a], hb.gdof[c], hb.kc_minus_k(a, c));
    SpMat dK(d.n_free, d.n_free);
    dK.setFromTriplets(trips.begin(), trips.end());
    const SpMat K = K0 + dK;
    try {
      fc.factorize(K);
      fc.require_spd();
    } catch (const Error&) {
      throw Error("hinge pushover: the tangent is singular or indefinite — a plastic mechanism has formed "
                  "(the limit load is reached)");
    }
    need_factor = false;
  };
  auto nodal_u = [&](const VecX& u, const HingedBar& hb) {
    Vec12 v;
    for (int a = 0; a < 12; ++a) v(a) = hb.gdof[a] >= 0 ? u(hb.gdof[a]) : 0.0;
    return v;
  };
  auto stress = [&](std::size_t i) -> Eigen::Vector2d {
    const HingeModel& h = hm[i];
    Eigen::Vector2d s = Eigen::Vector2d::Zero();
    s.head(h.m) = blk(h.Ke, h.m) * (st[i].e - st[i].ep).head(h.m);
    return s;
  };
  auto fval = [&](std::size_t i, int j) {
    const HingeModel& h = hm[i];
    return h.n[static_cast<std::size_t>(j)].head(h.m).dot((stress(i) - st[i].beta).head(h.m)) -
           h.r[static_cast<std::size_t>(j)];
  };

  // A load phase: its assembled vector and, per hinged bar, the span-load
  // fixed-end forces (release-condensed, local). With hinges the bar's span
  // load splits between the nodes and the internal hinge DOFs: the
  // equivalent nodal load is f_b - K_bi K_ii^-1 f_i and the hinge
  // deformations carry the offset -K_ii^-1 f_i — both with the CURRENT hinge
  // tangents, exact in rate form between events.
  struct Phase {
    const CaseLoads* cl = nullptr;
    std::vector<Vec12> cf0;  // per hinged bar
  };
  auto make_phase = [&](const CaseLoads& cl) {
    Phase ph;
    ph.cl = &cl;
    for (const HingedBar& hb : hbars) ph.cf0.push_back(asm_.elements[hb.e].condense_f0(cl.f0[hb.e]));
    return ph;
  };
  struct Rates {
    VecX du;
    std::vector<Eigen::Vector2d> de, dep;
  };
  auto rates = [&](const Phase& ph) {
    Rates R;
    R.de.assign(hinges.size(), Eigen::Vector2d::Zero());
    R.dep.assign(hinges.size(), Eigen::Vector2d::Zero());
    VecX D = ph.cl->rhs;
    std::vector<VecX> e0(hbars.size());
    for (std::size_t b = 0; b < hbars.size(); ++b) {
      const HingedBar& hb = hbars[b];
      const Vec12& f = ph.cf0[b];
      const int m = static_cast<int>(hb.qs.size());
      VecX fi(m);
      Vec12 fb = f;
      for (int h = 0; h < m; ++h) {
        fi(h) = f(hb.qs[static_cast<std::size_t>(h)]);
        fb(hb.qs[static_cast<std::size_t>(h)]) = 0;
      }
      e0[b] = -hb.Kii_inv * fi;
      if (f.cwiseAbs().maxCoeff() == 0.0) continue;
      const Vec12 feq = fb - hb.Kbi * (hb.Kii_inv * fi);
      const Vec12 corr = asm_.elements[hb.e].T().transpose() * (feq - f);
      for (int a = 0; a < 12; ++a)
        if (hb.gdof[a] >= 0) D(hb.gdof[a]) += corr(a);
    }
    R.du = fc.solve(D);
    for (std::size_t b = 0; b < hbars.size(); ++b) {
      const HingedBar& hb = hbars[b];
      const VecX g = hb.G * nodal_u(R.du, hb) + e0[b];
      for (int hi : hb.hinges) {
        const std::size_t i = static_cast<std::size_t>(hi);
        for (int a = 0; a < hm[i].m; ++a) R.de[i](a) = g(hinge_row[i] + a);
      }
    }
    for (std::size_t i = 0; i < hinges.size(); ++i) {
      if (st[i].active.empty()) continue;
      const HingeModel& h = hm[i];
      const Eigen::VectorXd lam = ops[i].Ginv * ops[i].N * blk(h.Ke, h.m) * R.de[i].head(h.m);
      R.dep[i].head(h.m) = ops[i].N.transpose() * lam;
    }
    return R;
  };

  VecX u = VecX::Zero(d.n_free);
  const std::size_t max_passes = 4 * hinges.size() + 8;
  // advance the load parameter from 0 to `target` along `dir` (base load held)
  auto advance = [&](const Phase& dir, double& par, double target, std::vector<double>& events) {
    while (par < target) {
      Rates R;
      for (std::size_t pass = 0;; ++pass) {
        require(pass < max_passes, "hinge pushover: no consistent active set (non-unique plastic response)");
        if (need_factor) factor();
        R = rates(dir);
        bool changed = false;
        for (std::size_t i = 0; i < hinges.size() && !changed; ++i) {
          const HingeModel& h = hm[i];
          // active facets must load (lambda_dot >= 0)
          if (!st[i].active.empty()) {
            const Eigen::VectorXd lam = ops[i].Ginv * ops[i].N * blk(h.Ke, h.m) * R.de[i].head(h.m);
            for (Eigen::Index k = 0; k < lam.size(); ++k)
              if (lam(k) < 0) {
                std::vector<int> act = st[i].active;
                act.erase(act.begin() + k);
                set_active(i, act);
                changed = true;
                break;
              }
            if (changed) break;
          }
          // inactive facets on the bound and loading join
          const Eigen::Vector2d ds = Eigen::Vector2d::Zero() + [&] {
            Eigen::Vector2d v = Eigen::Vector2d::Zero();
            v.head(h.m) = blk(h.Ke, h.m) * (R.de[i] - R.dep[i]).head(h.m);
            return v;
          }();
          const Eigen::Vector2d dbeta = h.H * R.dep[i];
          for (int j = 0; j < static_cast<int>(h.n.size()); ++j) {
            if (std::find(st[i].active.begin(), st[i].active.end(), j) != st[i].active.end()) continue;
            const double rj = h.r[static_cast<std::size_t>(j)];
            if (fval(i, j) >= -kTie * rj && h.n[static_cast<std::size_t>(j)].head(h.m).dot((ds - dbeta).head(h.m)) > 0) {
              std::vector<int> act = st[i].active;
              act.push_back(j);
              set_active(i, act);
              changed = true;
              break;
            }
          }
        }
        if (!changed) break;
      }
      // next event: an inactive facet reached
      double dl = target - par;
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        const HingeModel& h = hm[i];
        Eigen::Vector2d ds = Eigen::Vector2d::Zero();
        ds.head(h.m) = blk(h.Ke, h.m) * (R.de[i] - R.dep[i]).head(h.m);
        const Eigen::Vector2d dbeta = h.H * R.dep[i];
        for (int j = 0; j < static_cast<int>(h.n.size()); ++j) {
          if (std::find(st[i].active.begin(), st[i].active.end(), j) != st[i].active.end()) continue;
          const double rate = h.n[static_cast<std::size_t>(j)].head(h.m).dot((ds - dbeta).head(h.m));
          if (rate <= 0) continue;
          dl = std::min(dl, std::max(0.0, -fval(i, j) / rate));
        }
      }
      u += dl * R.du;
      par += dl;
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        st[i].e += dl * R.de[i];
        st[i].ep += dl * R.dep[i];
        st[i].beta += dl * (hm[i].H * R.dep[i]);
      }
      // facets reached (ties included) join at once; the next consistency
      // pass drops any whose multiplier rate turns out negative
      bool event = false;
      for (std::size_t i = 0; i < hinges.size(); ++i) {
        const HingeModel& h = hm[i];
        Eigen::Vector2d ds = Eigen::Vector2d::Zero();
        ds.head(h.m) = blk(h.Ke, h.m) * (R.de[i] - R.dep[i]).head(h.m);
        const Eigen::Vector2d dbeta = h.H * R.dep[i];
        std::vector<int> act = st[i].active;
        for (int j = 0; j < static_cast<int>(h.n.size()); ++j) {
          if (std::find(act.begin(), act.end(), j) != act.end()) continue;
          const double rate = h.n[static_cast<std::size_t>(j)].head(h.m).dot((ds - dbeta).head(h.m));
          if (rate > 0 && fval(i, j) >= -kTie * h.r[static_cast<std::size_t>(j)] &&
              static_cast<int>(act.size()) < h.m)
            act.push_back(j);
        }
        if (act.size() != st[i].active.size()) {
          set_active(i, act);
          event = true;
        }
      }
      if (event) events.push_back(par);
      if (target - par <= 1e-15 * std::max(1.0, target)) par = target;
    }
  };
  auto snapshot = [&](double lam) {
    HingePushoverPoint pt;
    pt.lambda = lam;
    pt.node_disp.assign(model.nodes.size(), Vec6::Zero());
    for (std::size_t n = 0; n < model.nodes.size(); ++n)
      for (int c = 0; c < 6; ++c) {
        const int fi = d.free_ix[static_cast<std::size_t>(d.full(static_cast<int>(n), c))];
        if (fi >= 0) pt.node_disp[n](c) = u(fi);
      }
    for (std::size_t i = 0; i < hinges.size(); ++i) {
      const HingeModel& h = hm[i];
      const Eigen::Vector2d s = stress(i);
      HingeState hs;
      hs.plastic = !st[i].active.empty();
      if (h.m == 1) {
        hs.theta = st[i].e(0);
        hs.theta_p = st[i].ep(0);
        hs.M = s(0);
      } else {
        hs.delta = h.sgnN * st[i].e(0);  // elongation, conjugate to N
        hs.delta_p = h.sgnN * st[i].ep(0);
        hs.N = h.sgnN * s(0);
        hs.theta = st[i].e(1);
        hs.theta_p = st[i].ep(1);
        hs.M = s(1);
      }
      pt.hinges.push_back(hs);
    }
    return pt;
  };

  HingePushoverResult out;
  if (opt.permanent_case) {
    double mu = 0;
    advance(make_phase(clPerm), mu, 1.0, out.permanent_events);
  }
  out.permanent = snapshot(0.0);
  double lam = 0;
  const Phase phP = make_phase(clP);
  for (double target : lambdas) {
    advance(phP, lam, target, out.events);
    out.points.push_back(snapshot(target));
  }
  return out;
}

}  // namespace keta

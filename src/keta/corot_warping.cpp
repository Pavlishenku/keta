#include "keta/corot_warping.h"

#include <cmath>
#include <limits>

#include "keta/corot_kin.h"
#include "keta/linsolve.h"

namespace keta {

namespace {
using corot::CorotElem;

struct WElem {
  CorotElem ce;
  Eigen::Matrix4d Ktw = Eigen::Matrix4d::Zero();  // (theta1, psi1, theta2, psi2)
  bool warp = false;
  double EIw = 0;
  // the element energy on the 14 element DOFs (q12, psi_i, psi_j)
  template <class T>
  T energy(const std::array<T, 14>& x) const {
    std::array<T, 12> q;
    for (int i = 0; i < 12; ++i) q[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)];
    const std::array<T, 7> d = ce.deform_t(q);
    T w = d[0] * (0.5 * ce.Kbar(0, 0)) * d[0];
    for (int i = 0; i < 7; ++i)
      for (int j = 0; j < 7; ++j) {
        if (i == 0 && j == 0) continue;
        if (warp && (i == 1 || i == 4) && (j == 1 || j == 4)) continue;  // the torsion: K_tw below
        if (ce.Kbar(i, j) != 0.0) w = w + d[static_cast<std::size_t>(i)] * (0.5 * ce.Kbar(i, j)) * d[static_cast<std::size_t>(j)];
      }
    if (warp) {
      const std::array<T, 4> t{d[1], x[12], d[4], x[13]};
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          if (Ktw(i, j) != 0.0) w = w + t[static_cast<std::size_t>(i)] * (0.5 * Ktw(i, j)) * t[static_cast<std::size_t>(j)];
    }
    return w;
  }
};
}  // namespace

CorotWarpingResults solve_corotational_warping(const Model& model, Id case_id, const WarpingSpec& spec,
                                               const CorotOptions& opt) {
  require(model.continua.empty(), "solve_corotational_warping: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  require(model.shells.empty() && !model.has_links() && model.composite.empty(),
          "corotational warping: bar models only (no shells, links or composite bars)");
  const Assembled asm_ = assemble_elements(model);
  const int nn = static_cast<int>(model.nodes.size());
  const int ndof = 7 * nn;
  // elements
  std::vector<WElem> el;
  std::vector<char> warp_node(static_cast<std::size_t>(nn), 0);
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    require(b.kind == BarKind::Frame, "corotational warping: frame bars only");
    for (bool r : b.releases) require(!r, "corotational warping: bars with releases are not supported");
    require(!b.soil.any() && !model.tapered_of(b.id), "corotational warping: no Winkler bed or CARA VAR");
    WElem w;
    w.ce = corot::make_corot_elem(model, asm_, e, opt.measure);
    const auto it = spec.Iw.find(b.section);
    // the spec overrides; else the section's own warping constant
    const SectionProps& sp0 = model.section(b.section).props;
    require(sp0.beta_x == 0.0, "corotational warping: monosymmetric section (beta_x != 0) not covered");
    const double Iw = it != spec.Iw.end() ? it->second : sp0.Iw;
    if (Iw > 0) {
      const Material& mat = model.material(b.material);
      const double L = w.ce.l0, EIw = mat.E * Iw, GJ = mat.G() * model.section(b.section).props.J;
      w.EIw = EIw;
      Eigen::Matrix4d kb, kg;
      kb << 12, 6 * L, -12, 6 * L, 6 * L, 4 * L * L, -6 * L, 2 * L * L, -12, -6 * L, 12, -6 * L, 6 * L, 2 * L * L, -6 * L,
          4 * L * L;
      kb *= EIw / (L * L * L);
      kg << 36, 3 * L, -36, 3 * L, 3 * L, 4 * L * L, -3 * L, -L * L, -36, -3 * L, 36, -3 * L, 3 * L, -L * L, -3 * L,
          4 * L * L;
      kg *= GJ / (30 * L);
      w.Ktw = kb + kg;
      w.warp = true;
      warp_node[static_cast<std::size_t>(w.ce.ni)] = 1;
      warp_node[static_cast<std::size_t>(w.ce.nj)] = 1;
    }
    el.push_back(w);
  }
  // DOFs: supports, restrained warping, warping where no warping bar reaches
  std::vector<char> blocked(static_cast<std::size_t>(ndof), 0);
  for (const Support& s : model.supports) {
    const int ni = model.node_index(s.node);
    for (int c = 0; c < 6; ++c) {
      require(s.spring[static_cast<std::size_t>(c)] == 0 && !s.spring6, "corotational warping: no elastic supports");
      if (s.fixed[static_cast<std::size_t>(c)]) blocked[static_cast<std::size_t>(7 * ni + c)] = 1;
    }
  }
  for (Id n : spec.warping_fixed) blocked[static_cast<std::size_t>(7 * model.node_index(n) + 6)] = 1;
  for (int ni = 0; ni < nn; ++ni)
    if (!warp_node[static_cast<std::size_t>(ni)]) blocked[static_cast<std::size_t>(7 * ni + 6)] = 1;
  std::vector<int> fr(static_cast<std::size_t>(ndof), -1);
  int nf = 0;
  for (int g = 0; g < ndof; ++g)
    if (!blocked[static_cast<std::size_t>(g)]) fr[static_cast<std::size_t>(g)] = nf++;
  // loads: dead nodal forces
  VecX F = VecX::Zero(ndof);
  for (const Load& l : model.load_case(case_id).loads) {
    require(l.kind == Load::Kind::NodalForce && l.M.squaredNorm() == 0,
            "corotational warping: nodal forces only (no moments, span loads or settlements)");
    const int ni = model.node_index(l.node);
    F.segment<3>(7 * ni) += l.F;
  }
  auto gather = [&](const WElem& w, const VecX& u, std::array<double, 14>& x) {
    for (int c = 0; c < 6; ++c) {
      x[static_cast<std::size_t>(c)] = u(7 * w.ce.ni + c);
      x[static_cast<std::size_t>(6 + c)] = u(7 * w.ce.nj + c);
    }
    x[12] = u(7 * w.ce.ni + 6);
    x[13] = u(7 * w.ce.nj + 6);
  };
  auto gdof = [&](const WElem& w, int a) {
    if (a < 6) return 7 * w.ce.ni + a;
    if (a < 12) return 7 * w.ce.nj + a - 6;
    return a == 12 ? 7 * w.ce.ni + 6 : 7 * w.ce.nj + 6;
  };
  // internal force (and tangent) on the full DOFs
  auto evaluate = [&](const VecX& u, VecX& fint, SpMat* K) {
    fint = VecX::Zero(ndof);
    std::vector<Triplet> t;
    for (const WElem& w : el) {
      std::array<double, 14> x;
      gather(w, u, x);
      if (K) {
        std::array<ad::Jet<14>, 14> xj;
        for (int i = 0; i < 14; ++i) xj[static_cast<std::size_t>(i)] = ad::Jet<14>::seed(x[static_cast<std::size_t>(i)], i);
        const ad::Jet<14> e = w.energy(xj);
        for (int a = 0; a < 14; ++a) {
          if (!w.warp && a >= 12) continue;
          fint(gdof(w, a)) += e.g[static_cast<std::size_t>(a)];
          const int fa = fr[static_cast<std::size_t>(gdof(w, a))];
          if (fa < 0) continue;
          for (int b = 0; b < 14; ++b) {
            if (!w.warp && b >= 12) continue;
            const int fb = fr[static_cast<std::size_t>(gdof(w, b))];
            if (fb >= 0) t.emplace_back(fa, fb, e.h[static_cast<std::size_t>(14 * a + b)]);
          }
        }
      } else {
        std::array<ad::Grad<14>, 14> xg;
        for (int i = 0; i < 14; ++i) xg[static_cast<std::size_t>(i)] = ad::Grad<14>::seed(x[static_cast<std::size_t>(i)], i);
        const ad::Grad<14> e = w.energy(xg);
        for (int a = 0; a < 14; ++a) {
          if (!w.warp && a >= 12) continue;
          fint(gdof(w, a)) += e.g[static_cast<std::size_t>(a)];
        }
      }
    }
    if (K) {
      K->resize(nf, nf);
      K->setFromTriplets(t.begin(), t.end());
    }
  };
  auto free_part = [&](const VecX& v) {
    VecX o(nf);
    for (int g = 0; g < ndof; ++g)
      if (fr[static_cast<std::size_t>(g)] >= 0) o(fr[static_cast<std::size_t>(g)]) = v(g);
    return o;
  };
  CorotWarpingResults out;
  VecX u = VecX::Zero(ndof);
  FactorCache fac;
  const int ns = std::max(1, opt.substeps);
  for (int st = 1; st <= ns; ++st) {
    const double lam = static_cast<double>(st) / ns;
    const VecX Fl = free_part(lam * F);
    const double scale = std::max(Fl.norm(), 1e-300);
    bool conv = false;
    for (int it = 0; it < opt.max_iterations; ++it) {
      VecX fint;
      SpMat K;
      evaluate(u, fint, &K);
      const VecX r = Fl - free_part(fint);
      // the round-off floor of K u (stiff axial terms: eps ||K|| |u| that no
      // iteration removes), as solve_corotational
      double kmax = 0;
      for (int k = 0; k < K.outerSize(); ++k)
        for (SpMat::InnerIterator iq(K, k); iq; ++iq) kmax = std::max(kmax, std::fabs(iq.value()));
      const double floor = 64 * std::numeric_limits<double>::epsilon() * kmax * u.lpNorm<Eigen::Infinity>() *
                           std::sqrt(static_cast<double>(nf));
      if (r.norm() <= std::max(opt.tol_rel * scale, floor)) {
        conv = true;
        break;
      }
      fac.factorize(K);
      const VecX du = fac.solve(r);
      // line search (Crisfield vol. 1 sec. 9.3, as solve_corotational):
      // accept once the directional derivative of the potential along du
      // has dropped below 0.8 of its start value, or the inf-norm residual
      // has dropped (a stiff load arm carries a large but harmless axial
      // residual on the linear predictor)
      const double g0 = std::fabs(du.dot(r)), r0 = r.lpNorm<Eigen::Infinity>();
      double s = 1.0;
      VecX un;
      for (int bt = 0; bt <= opt.max_backtracks; ++bt, s *= 0.5) {
        un = u;
        for (int g = 0; g < ndof; ++g)
          if (fr[static_cast<std::size_t>(g)] >= 0) un(g) += s * du(fr[static_cast<std::size_t>(g)]);
        VecX f2;
        evaluate(un, f2, nullptr);
        const VecX rt = Fl - free_part(f2);
        if (std::fabs(du.dot(rt)) <= 0.8 * g0 || rt.lpNorm<Eigen::Infinity>() < r0 || bt == opt.max_backtracks) break;
      }
      u = un;
      ++out.iterations;
    }
    require(conv, "corotational warping: no convergence at load factor " + std::to_string(lam));
  }
  out.node_disp.resize(static_cast<std::size_t>(nn));
  for (int ni = 0; ni < nn; ++ni) out.node_disp[static_cast<std::size_t>(ni)] = u.segment<7>(7 * ni);
  // end forces: Kbar d (the torsion of warping bars from K_tw: the end
  // torques and the bimoments)
  for (const WElem& w : el) {
    std::array<double, 14> x;
    gather(w, u, x);
    std::array<double, 12> q;
    for (int i = 0; i < 12; ++i) q[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)];
    const std::array<double, 7> d = w.ce.deform_t(q);
    Eigen::Matrix<double, 7, 1> dv;
    for (int i = 0; i < 7; ++i) dv(i) = d[static_cast<std::size_t>(i)];
    Eigen::Matrix<double, 7, 1> fl = w.ce.Kbar * dv;
    Eigen::Vector2d bm = Eigen::Vector2d::Zero();
    if (w.warp) {
      const Eigen::Vector4d tt(d[1], x[12], d[4], x[13]);
      const Eigen::Vector4d ft = w.Ktw * tt;
      fl(1) = ft(0);
      fl(4) = ft(2);
      bm = Eigen::Vector2d(ft(1), ft(3));
    }
    Vec12 f = Vec12::Zero();
    f(6) = fl(0);
    f(0) = -fl(0);
    f.segment<3>(3) = fl.segment<3>(1);
    f.segment<3>(9) = fl.segment<3>(4);
    const double ln = w.ce.l0;
    f(1) = -(fl(3) + fl(6)) / ln;
    f(7) = -f(1);
    f(2) = (fl(2) + fl(5)) / ln;
    f(8) = -f(2);
    out.bar_forces.push_back(f);
    out.bimoments.push_back(bm);
  }
  return out;
}

}  // namespace keta

#include "keta/link_nl.h"

#include <Eigen/SparseLU>
#include <functional>

#include <limits>

#include <cmath>
#include <map>

#include "keta/assembly.h"
#include "keta/linsolve.h"

namespace keta {

LinkLawCore::LinkLawCore(const Model& model, const std::vector<std::pair<std::size_t, int>>& extra, bool allow_softening)
    : w_(model) {
  model.validate();
  if (!allow_softening)
    for (const ElasticLink& l : model.links)
      for (const auto& law : l.law)
        if (law)
          for (std::size_t i = 1; i < law->pts.size(); ++i)
            require(law->pts[i].second >= law->pts[i - 1].second,
                    "link " + std::to_string(l.id) + ": a softening law (force decreasing) needs the path-following "
                    "driver solve_link_laws_path, not the convex solve");
  std::vector<std::pair<std::size_t, int>> list;
  for (std::size_t i = 0; i < model.links.size(); ++i)
    for (int c = 0; c < 6; ++c)
      if (model.links[i].law[static_cast<std::size_t>(c)]) list.emplace_back(i, c);
  for (const auto& e : extra) {
    require(e.first < model.links.size() && e.second >= 0 && e.second < 6, "link laws: bad extra component");
    for (const auto& q : list) require(q != e, "link laws: a component listed twice");
    for (int q = 0; q < 6; ++q)
      require(q == e.second || (model.links[e.first].k(e.second, q) == 0 && model.links[e.first].k(q, e.second) == 0),
              "link laws: a component with a law takes no coupling term in k");
    list.push_back(e);
  }
  // a unit marker while the DOF map is built: the zero-stiffness auto-fix
  // must keep the DOFs a law holds (a node held by a law alone was
  // eliminated, and an empty system crashed the solve)
  for (const auto& [i, c] : list) {
    comps_.push_back({i, c, {}, {}});
    w_.links[i].k(c, c) = 1.0;
  }
  for (ElasticLink& l : w_.links) l.law = {};
  asm_ = assemble_elements(w_, true);
  // the linear part proper: the laws' components out of the links
  for (const Comp& q : comps_) w_.links[q.link].k(q.c, q.c) = 0;
  for (std::size_t i = 0; i < w_.links.size(); ++i) asm_.links[i] = link_element(w_, w_.links[i]);
  K_ = assemble_stiffness(w_, asm_);
  const DofMap& d = asm_.dofs;
  for (Comp& q : comps_) {
    const LinkElem e = link_element(model, model.links[q.link]);
    for (int a = 0; a < 12; ++a) {
      const double coef = e.B(q.c, a);
      if (coef == 0) continue;
      const int full = d.full(a < 6 ? e.ni : e.nj, a % 6);
      q.bfull.emplace_back(full, coef);
      const int fi = d.free_ix[static_cast<std::size_t>(full)];
      if (fi >= 0) q.bfree.emplace_back(fi, coef);
    }
  }
}

CaseLoads LinkLawCore::lower(const LoadCase& lc) const {
  CaseLoads cl = lower_case(w_, asm_, lc);
  fold_prescribed(w_, asm_, cl);
  fold_link_offsets(w_, asm_, K_, cl);
  return cl;
}

std::vector<double> LinkLawCore::d0(const CaseLoads& cl) const {
  std::map<int, double> presc(cl.prescribed.begin(), cl.prescribed.end());
  std::vector<double> out(comps_.size(), 0.0);
  for (std::size_t j = 0; j < comps_.size(); ++j)
    for (const auto& [full, coef] : comps_[j].bfull) {
      const int fi = asm_.dofs.free_ix[static_cast<std::size_t>(full)];
      if (fi >= 0) {
        if (cl.link_u0.size()) out[j] += coef * cl.link_u0(fi);
      } else if (auto it = presc.find(full); it != presc.end()) {
        out[j] += coef * it->second;
      }
    }
  return out;
}

double LinkLawCore::deformation(std::size_t j, const VecX& v, const std::vector<double>& d0) const {
  double s = d0[j];
  for (const auto& [i, c] : comps_[j].bfree) s += c * v(i);
  return s;
}

int LinkLawCore::solve(const VecX& rhs, const std::vector<double>& d0, const std::vector<const LinkLaw*>& laws,
                       VecX& v, int max_iter) const {
  require(laws.size() == comps_.size() && d0.size() == comps_.size(), "link laws: one law per component");
  const int n = asm_.dofs.n_free;
  if (v.size() != n) v = VecX::Zero(n);
  FactorCache cache;
  if (asm_.red.active) cache.set_reduction(asm_.red.T);
  auto energy = [&](const VecX& x) {
    double e = 0.5 * x.dot(K_ * x) - rhs.dot(x);
    for (std::size_t j = 0; j < comps_.size(); ++j) e += laws[j]->energy(deformation(j, x, d0));
    return e;
  };
  auto residual = [&](const VecX& x) {
    VecX r = rhs - K_ * x;
    for (std::size_t j = 0; j < comps_.size(); ++j) {
      const double F = laws[j]->force(deformation(j, x, d0));
      for (const auto& [i, c] : comps_[j].bfree) r(i) -= c * F;
    }
    return r;
  };
  auto reduced_norm = [&](const VecX& r) { return asm_.red.active ? VecX(asm_.red.reduce(r)).norm() : r.norm(); };
  double fscale = n > 0 ? rhs.cwiseAbs().maxCoeff() : 0.0;
  for (std::size_t j = 0; j < comps_.size(); ++j) fscale = std::max(fscale, std::fabs(laws[j]->force(d0[j])));
  int it = 0;
  // every DOF prescribed: the state is known, nothing to solve
  for (; n > 0 && it < max_iter; ++it) {
    const VecX r = residual(v);
    // tangent: K + sum kt b b' (the slope of each law's current segment)
    std::vector<Triplet> trips;
    for (int k = 0; k < K_.outerSize(); ++k)
      for (SpMat::InnerIterator iter(K_, k); iter; ++iter) trips.emplace_back(iter.row(), iter.col(), iter.value());
    for (std::size_t j = 0; j < comps_.size(); ++j) {
      const double kt = laws[j]->tangent(deformation(j, v, d0));
      if (kt == 0) continue;
      for (const auto& [a, ca] : comps_[j].bfree)
        for (const auto& [b, cb] : comps_[j].bfree) trips.emplace_back(a, b, kt * ca * cb);
    }
    SpMat Kt(n, n);
    Kt.setFromTriplets(trips.begin(), trips.end());
    // the round-off floor of the residual: K v carries ~ eps max|K| |v| of
    // noise (a stiff contact law, 1e12 against a beam's 1e5, sits far above
    // a tolerance relative to the loads)
    double kmax = 0;
    for (int k = 0; k < Kt.outerSize(); ++k)
      for (SpMat::InnerIterator iter(Kt, k); iter; ++iter) kmax = std::max(kmax, std::fabs(iter.value()));
    // a slack component leaves its DOFs without stiffness only if nothing
    // else holds them: a mechanism is a typed error of the factorization
    cache.factorize(Kt);
    const VecX dv = cache.solve(r);
    // backtracking on the convex energy (the full step is exact once the
    // segments are right)
    const double e0 = energy(v);
    double alpha = 1.0;
    VecX vn = v + dv;
    for (int ls = 0; ls < 60 && energy(vn) > e0 + 1e-15 * std::max(1.0, std::fabs(e0)); ++ls) {
      alpha *= 0.5;
      vn = v + alpha * dv;
    }
    v = vn;
    const double rn = reduced_norm(residual(v));
    const double floor = 64 * std::numeric_limits<double>::epsilon() * kmax * v.cwiseAbs().maxCoeff();
    if (alpha == 1.0 && (rn <= 1e-11 * std::max(fscale, 1.0) || rn <= floor)) {
      ++it;
      break;
    }
    // the energy's own resolution: once the backtracking cannot decrease it
    // any more, the residual it can still certify is ~ sqrt(eps |E| kmax)
    // (a step of residual r lowers the energy by ~ r^2 / k): very stiff
    // components (1e12 against 1e7) put that floor above the round-off of
    // K v. The minimum is then reached to the precision of the energy.
    const double efloor =
        std::sqrt(64 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::fabs(e0)) * kmax);
    if (alpha < 1e-6 && rn <= efloor) {
      ++it;
      break;
    }
  }
  require(it < max_iter, "link laws: no convergence (non-monotone law, or a mechanism)");
  return it;
}

CaseResult LinkLawCore::recover(const CaseLoads& cl, const VecX& v, const std::vector<double>& F, int stations, Id id,
                                const std::string& name) const {
  // the laws' forces enter as nodal actions on the structure (the reactions
  // of the supported nodes and the rigid links' transfers see them)
  CaseLoads cr = cl;
  for (std::size_t j = 0; j < comps_.size(); ++j)
    for (const auto& [full, c] : comps_[j].bfull) cr.nodal_full(full) -= c * F[j];
  StaticOptions so;
  so.stations_per_bar = stations;
  CaseResult out = recover_case(w_, asm_, cr, v, so, id, name);
  out.non_superposable = true;
  for (std::size_t j = 0; j < comps_.size(); ++j) out.links[comps_[j].link].force(comps_[j].c) = F[j];
  return out;
}

LinkLawResult solve_static_link_laws(const Model& model, Id case_id, const LinkLawOptions& opt) {
  const LinkLawCore core(model);
  const LoadCase& lc = model.load_case(case_id);
  const CaseLoads cl = core.lower(lc);
  const std::vector<double> d0 = core.d0(cl);
  std::vector<const LinkLaw*> laws;
  for (const LinkLawCore::Comp& q : core.comps())
    laws.push_back(&*model.links[q.link].law[static_cast<std::size_t>(q.c)]);
  VecX v = VecX::Zero(core.n_free());
  LinkLawResult out;
  out.iterations = core.solve(cl.rhs, d0, laws, v, opt.max_iter);
  std::vector<double> F(laws.size());
  for (std::size_t j = 0; j < laws.size(); ++j) F[j] = laws[j]->force(core.deformation(j, v, d0));
  out.result = core.recover(cl, v, F, opt.stations_per_bar, lc.id, lc.name);
  return out;
}

LinkPathResult solve_link_laws_path(const Model& model, Id case_id, const LinkPathOptions& opt) {
  require((!opt.targets.empty() || opt.arc_steps > 0) && opt.steps_per_leg >= 1,
          "link path: targets (or arc_steps) and steps_per_leg >= 1 needed");
  require(opt.ctrl_comp >= 0 && opt.ctrl_comp < 6, "link path: ctrl_comp 0..5");
  // the family components (not carrying a model law) join as extra components
  std::vector<std::pair<std::size_t, int>> extra;
  auto link_ix = [&](Id id) {
    for (std::size_t i = 0; i < model.links.size(); ++i)
      if (model.links[i].id == id) return i;
    throw Error("link path: no link " + std::to_string(id));
  };
  for (const LinkHingeFamily& h : opt.families) {
    require(h.family != nullptr, "link path: a family is required");
    h.family->validate();
    require(h.comp >= 0 && h.comp < 6 && h.axial_comp >= 0 && h.axial_comp < 6 && h.comp != h.axial_comp,
            "link path: components 0..5, distinct");
    const std::size_t li = link_ix(h.link);
    require(!model.links[li].law[static_cast<std::size_t>(h.comp)], "link path: a family component carries no model law");
    require(static_cast<bool>(model.links[li].law[static_cast<std::size_t>(h.axial_comp)]),
            "link path: the axial component of a family needs a law (its force is N)");
    extra.emplace_back(li, h.comp);
  }
  LinkLawCore core(model, extra, true);
  const Assembled& as = core.assembled();
  require(!as.red.active, "link path: rigid links, constraints and skew supports are not supported");
  const int n = core.n_free();
  const auto& comps = core.comps();
  const std::size_t nm = comps.size() - extra.size();
  std::vector<const LinkLaw*> laws(comps.size(), nullptr);
  for (std::size_t j = 0; j < nm; ++j) laws[j] = &*model.links[comps[j].link].law[static_cast<std::size_t>(comps[j].c)];
  std::vector<std::size_t> axial_of(comps.size(), 0);
  for (std::size_t q = 0; q < opt.families.size(); ++q) {
    const std::size_t j = nm + q, li = comps[j].link;
    bool found = false;
    for (std::size_t a = 0; a < nm; ++a)
      if (comps[a].link == li && comps[a].c == opt.families[q].axial_comp) {
        axial_of[j] = a;
        found = true;
      }
    require(found, "link path: the axial component of a family must be a law component");
  }
  const CaseLoads cl = core.lower(model.load_case(case_id));
  require(cl.prescribed.empty(), "link path: prescribed displacements are not supported");
  VecX fdead = VecX::Zero(n);
  if (opt.dead_case != 0) {
    const CaseLoads cd = core.lower(model.load_case(opt.dead_case));
    require(cd.prescribed.empty(), "link path: prescribed displacements are not supported");
    fdead = cd.rhs;
  }
  const VecX f = cl.rhs;
  const std::vector<double> d0(comps.size(), 0.0);
  const int cf = as.dofs.free_of(model.node_index(opt.ctrl_node), opt.ctrl_comp);
  require(cf >= 0, "link path: the controlled DOF is not free");
  // forces and tangents of the components at u
  struct St {
    std::vector<double> d, F, k, dFdN;
  };
  auto state = [&](const VecX& u) {
    St s;
    const std::size_t m = comps.size();
    s.d.assign(m, 0);
    s.F.assign(m, 0);
    s.k.assign(m, 0);
    s.dFdN.assign(m, 0);
    for (std::size_t j = 0; j < m; ++j) s.d[j] = core.deformation(j, u, d0);
    for (std::size_t j = 0; j < nm; ++j) {
      s.F[j] = laws[j]->force(s.d[j]);
      s.k[j] = laws[j]->tangent(s.d[j]);
    }
    for (std::size_t q = 0; q < opt.families.size(); ++q) {
      const std::size_t j = nm + q;
      const HingeLawFamily& fam = *opt.families[q].family;
      const double N = s.F[axial_of[j]];
      const LinkLaw l = fam.at(N);
      s.F[j] = l.force(s.d[j]);
      s.k[j] = l.tangent(s.d[j]);
      if (fam.N.size() > 1) {
        const double h = 1e-6 * (fam.N.back() - fam.N.front());
        const double lo = std::max(fam.N.front(), N - h), hi = std::min(fam.N.back(), N + h);
        s.dFdN[j] = (fam.at(hi).force(s.d[j]) - fam.at(lo).force(s.d[j])) / (hi - lo);
      }
    }
    return s;
  };
  auto residual = [&](const VecX& u, double lam, const St& s) {
    VecX r = fdead + lam * f - core.K() * u;
    for (std::size_t j = 0; j < comps.size(); ++j)
      for (const auto& [i, c] : comps[j].bfree) r(i) -= c * s.F[j];
    return r;
  };
  auto tangent = [&](const St& s) {
    std::vector<Triplet> t;
    const SpMat& K = core.K();
    for (int k = 0; k < K.outerSize(); ++k)
      for (SpMat::InnerIterator it(K, k); it; ++it) t.emplace_back(it.row(), it.col(), it.value());
    for (std::size_t j = 0; j < comps.size(); ++j) {
      for (const auto& [a, ca] : comps[j].bfree)
        for (const auto& [b, cb] : comps[j].bfree) t.emplace_back(a, b, s.k[j] * ca * cb);
      if (s.dFdN[j] != 0) {
        const std::size_t ax = axial_of[j];
        for (const auto& [a, ca] : comps[j].bfree)
          for (const auto& [b, cb] : comps[ax].bfree) t.emplace_back(a, b, s.dFdN[j] * s.k[ax] * ca * cb);
      }
    }
    SpMat Kt(n, n);
    Kt.setFromTriplets(t.begin(), t.end());
    return Kt;
  };
  LinkPathResult out;
  VecX u = VecX::Zero(n);
  double lam = 0;
  const double fscale = std::max({f.cwiseAbs().maxCoeff(), fdead.cwiseAbs().maxCoeff(), 1e-30});
  // one bordered Newton to u_ctrl = tgt from (u, lam); false on failure.
  // The augmented system [Kt -f; e_c' 0] (du, dlam) = (r, tgt - u_c) is
  // solved whole: it stays regular on a plateau (Kt singular) and through
  // limit points. After the first (full) step the control holds exactly and
  // every later step keeps it, so it is damped by backtracking on |r| (a
  // Newton cycle between two segments of a law cannot persist).
  auto rnorm = [&](const VecX& uu, double ll) {
    const St s = state(uu);
    return residual(uu, ll, s).cwiseAbs().maxCoeff();
  };
  auto step_to = [&](double tgt, VecX& uu, double& ll) {
    for (int it = 0; it < opt.max_iter; ++it) {
      const St s = state(uu);
      const VecX r = residual(uu, ll, s);
      double fsc = fscale;
      for (double F : s.F) fsc = std::max(fsc, std::fabs(F));
      const double rn = r.cwiseAbs().maxCoeff();
      const bool on = std::fabs(uu(cf) - tgt) <= 1e-12 * std::max(1.0, std::fabs(tgt));
      if (rn <= opt.tol * fsc && on) return true;
      ++out.iterations;
      const SpMat Kt = tangent(s);
      std::vector<Triplet> t;
      for (int k = 0; k < Kt.outerSize(); ++k)
        for (SpMat::InnerIterator i2(Kt, k); i2; ++i2) t.emplace_back(i2.row(), i2.col(), i2.value());
      for (int i = 0; i < n; ++i)
        if (f(i) != 0) t.emplace_back(i, n, -f(i));
      t.emplace_back(n, cf, 1.0);
      SpMat A(n + 1, n + 1);
      A.setFromTriplets(t.begin(), t.end());
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(A);
      if (lu.info() != Eigen::Success) return false;
      VecX rhs(n + 1);
      rhs.head(n) = r;
      rhs(n) = tgt - uu(cf);
      const VecX x = lu.solve(rhs);
      if (!x.allFinite()) return false;
      double alpha = 1.0;
      if (on) {
        for (int ls = 0; ls < 30; ++ls) {
          if (rnorm(uu + alpha * x.head(n), ll + alpha * x(n)) < rn) break;
          alpha *= 0.5;
        }
      }
      uu += alpha * x.head(n);
      ll += alpha * x(n);
    }
    return false;
  };
  std::function<void(double, int)> advance = [&](double tgt, int depth) {
    VecX uu = u;
    double ll = lam;
    if (step_to(tgt, uu, ll)) {
      u = uu;
      lam = ll;
      return;
    }
    require(depth < opt.max_cuts, "link path: no equilibrium near u_ctrl = " + std::to_string(tgt) + " (step halved)");
    const double mid = 0.5 * (u(cf) + tgt);
    advance(mid, depth + 1);
    advance(tgt, depth + 1);
  };
  auto record = [&] {
    const St s = state(u);
    out.trace.push_back({lam, u(cf), s.d, s.F});
  };
  // the dead load first (lambda = 0, the controlled DOF free): plain Newton
  if (opt.dead_case != 0) {
    for (int it = 0;; ++it) {
      require(it < opt.max_iter, "link path: the dead load does not converge");
      const St s = state(u);
      const VecX r = residual(u, 0, s);
      if (r.cwiseAbs().maxCoeff() <= opt.tol * fscale) break;
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(tangent(s));
      require(lu.info() == Eigen::Success, "link path: singular tangent under the dead load");
      u += lu.solve(r);
      ++out.iterations;
    }
  }
  record();
  if (opt.arc_steps > 0) {
    require(opt.arc_length > 0 && opt.arc_psi >= 0 && (opt.arc_first_sign == 1 || opt.arc_first_sign == -1),
            "link path: arc_length > 0, arc_psi >= 0, arc_first_sign = +-1");
    const double w = opt.arc_psi * opt.arc_psi * f.squaredNorm();  // weight of dlam^2 in the arc norm
    const std::size_t m = comps.size();
    std::vector<const HingeLawFamily*> famof(m, nullptr);
    for (std::size_t q = 0; q < opt.families.size(); ++q) famof[nm + q] = opt.families[q].family.get();
    // ---- cells: per component its law segment (and, for a family, its N
    // interval); within a cell the forces are smooth (linear for a plain
    // law, rational in N for a family whose points move with N) ----
    struct Cell {
      std::vector<int> seg, iv;
    };
    auto law_pts = [&](std::size_t j, int iv, double t) {  // the points of the law (family: at N parameter t in interval iv)
      if (!famof[j]) return laws[j]->pts;
      const HingeLawFamily& fm = *famof[j];
      if (fm.N.size() == 1) return fm.laws[0].pts;
      std::vector<std::pair<double, double>> p;
      const auto& A = fm.laws[static_cast<std::size_t>(iv) - 1].pts;
      const auto& B = fm.laws[static_cast<std::size_t>(iv)].pts;
      for (std::size_t i = 0; i < A.size(); ++i)
        p.emplace_back(A[i].first + t * (B[i].first - A[i].first), A[i].second + t * (B[i].second - A[i].second));
      return p;
    };
    // the N interval parameter of a family component: t = (N - N_lo) / dN
    auto tN = [&](std::size_t j, int iv, double N, double& dN) {
      const HingeLawFamily& fm = *famof[j];
      if (fm.N.size() == 1) {
        dN = 1;
        return 0.0;
      }
      dN = fm.N[static_cast<std::size_t>(iv)] - fm.N[static_cast<std::size_t>(iv) - 1];
      return (N - fm.N[static_cast<std::size_t>(iv) - 1]) / dN;
    };
    // the state of the cell at u: deformations, forces, slopes, dF/dN (exact)
    auto eval = [&](const VecX& uu, const Cell& c) {
      St s;
      s.d.assign(m, 0);
      s.F.assign(m, 0);
      s.k.assign(m, 0);
      s.dFdN.assign(m, 0);
      for (std::size_t j = 0; j < m; ++j) s.d[j] = core.deformation(j, uu, d0);
      for (std::size_t j = 0; j < nm; ++j) {
        const auto& P = laws[j]->pts;
        const std::size_t i = static_cast<std::size_t>(c.seg[j]);
        s.k[j] = (P[i + 1].second - P[i].second) / (P[i + 1].first - P[i].first);
        s.F[j] = P[i].second + s.k[j] * (s.d[j] - P[i].first);
      }
      for (std::size_t j = nm; j < m; ++j) {
        const HingeLawFamily& fm = *famof[j];
        const double N = s.F[axial_of[j]];
        double dN = 1;
        const double t = tN(j, c.iv[j], N, dN);
        const std::size_t i = static_cast<std::size_t>(c.seg[j]);
        const auto P = law_pts(j, c.iv[j], t);
        const double th0 = P[i].first, th1 = P[i + 1].first, M0 = P[i].second, M1 = P[i + 1].second;
        const double sl = (M1 - M0) / (th1 - th0);
        s.k[j] = sl;
        s.F[j] = M0 + sl * (s.d[j] - th0);
        if (fm.N.size() > 1) {
          const auto& A = fm.laws[static_cast<std::size_t>(c.iv[j]) - 1].pts;
          const auto& B = fm.laws[static_cast<std::size_t>(c.iv[j])].pts;
          const double dth0 = B[i].first - A[i].first, dth1 = B[i + 1].first - A[i + 1].first;
          const double dM0 = B[i].second - A[i].second, dM1 = B[i + 1].second - A[i + 1].second;
          const double h = th1 - th0;
          const double dsl = ((dM1 - dM0) * h - (M1 - M0) * (dth1 - dth0)) / (h * h);
          s.dFdN[j] = (dM0 - dth0 * sl + (s.d[j] - th0) * dsl) / dN;
        }
      }
      return s;
    };
    // the boundaries of a cell: g >= 0 inside, with the gradient over u
    struct Bound {
      double g = 0;
      VecX grad;
      bool hard = false;  // the family's tabulated N range: no cell beyond (no continuation outward)
    };
    auto bvec = [&](std::size_t j) {
      VecX b = VecX::Zero(n);
      for (const auto& [i, cc] : comps[j].bfree) b(i) += cc;
      return b;
    };
    auto bounds = [&](const VecX& uu, const Cell& c) {
      const St s = eval(uu, c);
      std::vector<Bound> out_b;
      for (std::size_t j = 0; j < nm; ++j) {
        const auto& P = laws[j]->pts;
        const std::size_t i = static_cast<std::size_t>(c.seg[j]);
        if (i > 0) out_b.push_back({s.d[j] - P[i].first, bvec(j), false});
        if (i + 2 < P.size()) out_b.push_back({P[i + 1].first - s.d[j], -bvec(j), false});
      }
      for (std::size_t j = nm; j < m; ++j) {
        const HingeLawFamily& fm = *famof[j];
        const std::size_t ax = axial_of[j];
        const VecX gN = s.k[ax] * bvec(ax);  // dN/du
        const double N = s.F[ax];
        double dN = 1;
        const double t = tN(j, c.iv[j], N, dN);
        if (fm.N.size() > 1) {
          const std::size_t k = static_cast<std::size_t>(c.iv[j]);
          out_b.push_back({N - fm.N[k - 1], gN, k == 1});
          out_b.push_back({fm.N[k] - N, -gN, k + 1 == fm.N.size()});
        }
        const std::size_t i = static_cast<std::size_t>(c.seg[j]);
        const auto P = law_pts(j, c.iv[j], t);
        auto dth = [&](std::size_t p) {
          if (fm.N.size() == 1) return 0.0;
          return (fm.laws[static_cast<std::size_t>(c.iv[j])].pts[p].first - fm.laws[static_cast<std::size_t>(c.iv[j]) - 1].pts[p].first) / dN;
        };
        if (i > 0) out_b.push_back({s.d[j] - P[i].first, VecX(bvec(j) - dth(i) * gN), false});
        if (i + 2 < P.size()) out_b.push_back({P[i + 1].first - s.d[j], VecX(dth(i + 1) * gN - bvec(j)), false});
      }
      return out_b;
    };
    const double gtol = 1e-11;
    auto gscale = [&](const Bound& b) { return std::max(1.0, b.grad.cwiseAbs().maxCoeff()); };
    // the candidate cells of each component at u (those holding the point)
    auto candidates = [&](const VecX& uu, std::size_t j, const St& nat) {
      std::vector<std::pair<int, int>> c;  // (seg, iv)
      const double d = core.deformation(j, uu, d0);
      auto segs_holding = [&](const std::vector<std::pair<double, double>>& P) {
        std::vector<int> r;
        const double tol = gtol * std::max(1.0, std::fabs(P.back().first - P.front().first));
        for (std::size_t i = 0; i + 1 < P.size(); ++i) {
          const bool lo_ok = i == 0 || d >= P[i].first - tol;
          const bool hi_ok = i + 2 == P.size() || d <= P[i + 1].first + tol;
          if (lo_ok && hi_ok) r.push_back(static_cast<int>(i));
        }
        return r;
      };
      if (!famof[j]) {
        for (int i : segs_holding(laws[j]->pts)) c.emplace_back(i, -1);
        return c;
      }
      const HingeLawFamily& fm = *famof[j];
      const double N = nat.F[axial_of[j]];
      std::vector<int> ivs;
      if (fm.N.size() == 1) ivs.push_back(0);
      else {
        const double tol = gtol * std::max(1.0, std::fabs(fm.N.back() - fm.N.front()));
        for (std::size_t k = 1; k < fm.N.size(); ++k)
          if (N >= fm.N[k - 1] - tol && N <= fm.N[k] + tol) ivs.push_back(static_cast<int>(k));
      }
      for (int iv : ivs) {
        double dN = 1;
        const double t = fm.N.size() == 1 ? 0.0 : tN(j, iv, N, dN);
        for (int i : segs_holding(law_pts(j, iv, t))) c.emplace_back(i, iv);
      }
      return c;
    };
    // the plain-law forces at u (the axial forces of the families: plain laws)
    auto natural = [&](const VecX& uu) {
      St s;
      s.F.assign(m, 0);
      for (std::size_t j = 0; j < nm; ++j) s.F[j] = laws[j]->force(core.deformation(j, uu, d0));
      return s;
    };
    // bordered solve [Kt -f; a' b] x = rhs
    auto bsolve = [&](const St& s, const VecX& a_u, double a_l, const VecX& rhs, VecX& x) {
      const SpMat Kt = tangent(s);
      std::vector<Triplet> t;
      for (int k = 0; k < Kt.outerSize(); ++k)
        for (SpMat::InnerIterator i2(Kt, k); i2; ++i2) t.emplace_back(i2.row(), i2.col(), i2.value());
      for (int i = 0; i < n; ++i) {
        if (f(i) != 0) t.emplace_back(i, n, -f(i));
        if (a_u(i) != 0) t.emplace_back(n, i, a_u(i));
      }
      if (a_l != 0) t.emplace_back(n, n, a_l);
      SpMat A(n + 1, n + 1);
      A.setFromTriplets(t.begin(), t.end());
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(A);
      if (lu.info() != Eigen::Success) return false;
      x = lu.solve(rhs);
      return x.allFinite();
    };
    VecX pu = VecX::Zero(n);
    pu(cf) = opt.arc_first_sign;
    double pl = 0;
    Cell prev_cell;
    bool have_prev = false;
    double travelled = 0;
    int recorded = 0;
    for (int guard = 0; recorded < opt.arc_steps && guard < 100 * opt.arc_steps + 1000; ++guard) {
      // the cell to leave in, and the direction: every candidate combination,
      // the cells different from the previous one (continuing) first, both
      // signs of the tangent (the one along the previous direction first);
      // consistent = every boundary at zero opening into the cell
      const St nat = natural(u);
      std::vector<std::vector<std::pair<int, int>>> cand(m);
      std::size_t combos = 1;
      for (std::size_t j = 0; j < m; ++j) {
        cand[j] = candidates(u, j, nat);
        require(!cand[j].empty(), "link path: a family's axial force left its tabulated range");
        if (have_prev && cand[j].size() > 1)
          std::stable_sort(cand[j].begin(), cand[j].end(), [&](const auto& a, const auto& b) {
            const bool sa = a.first == prev_cell.seg[j] && a.second == prev_cell.iv[j];
            const bool sb = b.first == prev_cell.seg[j] && b.second == prev_cell.iv[j];
            return !sa && sb;
          });
        combos *= cand[j].size();
        require(combos <= 4096, "link path: too many components at a breakpoint at once");
      }
      Cell cell;
      VecX tu;
      double tl = 0;
      bool found = false;
      for (std::size_t mask = 0; mask < combos && !found; ++mask) {
        Cell c;
        c.seg.assign(m, 0);
        c.iv.assign(m, -1);
        std::size_t r = mask;
        for (std::size_t j = 0; j < m; ++j) {
          const auto& pr = cand[j][r % cand[j].size()];
          r /= cand[j].size();
          c.seg[j] = pr.first;
          c.iv[j] = pr.second;
        }
        const St s = eval(u, c);
        VecX rhs = VecX::Zero(n + 1), x;
        rhs(n) = 1;
        if (!bsolve(s, pu, w * pl, rhs, x)) continue;
        const double nt = std::sqrt(x.head(n).squaredNorm() + w * x(n) * x(n));
        if (!(nt > 0)) continue;
        const std::vector<Bound> B = bounds(u, c);
        for (double sg : {1.0, -1.0}) {
          const VecX t1 = sg * x.head(n) / nt;
          bool ok = true;
          for (const Bound& b : B)
            if (std::fabs(b.g) <= gtol * gscale(b) && b.grad.dot(t1) < -1e-12 * b.grad.norm() * t1.norm()) ok = false;
          if (ok) {
            cell = c;
            tu = t1;
            tl = sg * x(n) / nt;
            found = true;
            break;
          }
        }
      }
      require(found, "link path: no consistent continuation at a breakpoint");
      // the first boundary met (linear estimate), then the corrector: the
      // arc's normal plane, or the boundary itself (an exact event)
      const VecX u0 = u;
      const double l0 = lam;
      const double rest = opt.arc_length - travelled;
      const std::vector<Bound> B0 = bounds(u0, cell);
      int target = -1;
      double s_t = rest;
      for (std::size_t b = 0; b < B0.size(); ++b) {
        const double rate = B0[b].grad.dot(tu);
        if (rate < 0 && B0[b].g > gtol * gscale(B0[b])) {
          const double sb = -B0[b].g / rate;
          if (sb < s_t) {
            s_t = sb;
            target = static_cast<int>(b);
          }
        }
      }
      VecX un;
      double ln = 0;
      bool done = false;
      for (int attempt = 0; attempt < 8 && !done; ++attempt) {
        un = u0 + s_t * tu;
        ln = l0 + s_t * tl;
        bool conv = false;
        for (int it = 0; it < opt.max_iter; ++it) {
          const St s = eval(un, cell);
          const VecX r = residual(un, ln, s);
          double c;
          VecX a_u;
          double a_l = 0;
          if (target < 0) {
            c = tu.dot(un - u0) + w * tl * (ln - l0) - s_t;
            a_u = tu;
            a_l = w * tl;
          } else {
            const Bound bt = bounds(un, cell)[static_cast<std::size_t>(target)];
            c = bt.g;
            a_u = bt.grad;
          }
          double fsc = fscale;
          for (double F : s.F) fsc = std::max(fsc, std::fabs(F));
          if (r.cwiseAbs().maxCoeff() <= opt.tol * fsc && std::fabs(c) <= 1e-13 * std::max(1.0, std::fabs(s_t) + u0.norm())) {
            conv = true;
            break;
          }
          ++out.iterations;
          VecX rhs(n + 1), x;
          rhs.head(n) = r;
          rhs(n) = -c;
          if (!bsolve(s, a_u, a_l, rhs, x)) break;
          un += x.head(n);
          ln += x(n);
        }
        require(conv, "link path: the arc corrector did not converge");
        // a boundary crossed on the way (the path bends): retarget to the
        // one crossed first (linear interpolation of its opening)
        const std::vector<Bound> B1 = bounds(un, cell);
        int crossed = -1;
        double frac = 2;
        for (std::size_t b = 0; b < B1.size(); ++b) {
          if (static_cast<int>(b) == target) continue;
          if (B1[b].g < -gtol * gscale(B1[b])) {
            const double fr = B0[b].g / (B0[b].g - B1[b].g);
            if (fr < frac) {
              frac = fr;
              crossed = static_cast<int>(b);
            }
          }
        }
        if (crossed < 0) {
          done = true;
        } else {
          target = crossed;
          s_t = std::max(frac, 1e-6) * tu.dot(un - u0);
        }
      }
      require(done, "link path: the arc step could not land inside its cell");
      const double step = tu.dot(un - u0) + w * tl * (ln - l0);
      u = un;
      lam = ln;
      pu = tu;
      pl = tl;
      prev_cell = cell;
      have_prev = true;
      travelled += step;
      {
        const St s = state(u);
        const VecX r = residual(u, lam, s);
        double fsc = fscale;
        for (double F : s.F) fsc = std::max(fsc, std::fabs(F));
        require(r.cwiseAbs().maxCoeff() <= 1e-8 * fsc, "link path: the arc step left equilibrium");
      }
      if (travelled >= opt.arc_length * (1 - 1e-12) || target >= 0) {
        record();
        ++recorded;
        if (travelled >= opt.arc_length * (1 - 1e-12)) travelled = 0;
      }
      if (std::isfinite(opt.arc_stop_ctrl) && (u(cf) - opt.arc_stop_ctrl) * opt.arc_first_sign >= 0) break;
    }
  }
  double prev = u(cf);
  for (double target : opt.arc_steps > 0 ? std::vector<double>{} : opt.targets) {
    for (int k = 1; k <= opt.steps_per_leg; ++k) {
      advance(prev + (target - prev) * k / opt.steps_per_leg, 0);
      record();
    }
    prev = target;
  }
  out.node_disp.assign(model.nodes.size(), Vec6::Zero());
  for (std::size_t ni = 0; ni < model.nodes.size(); ++ni)
    for (int c = 0; c < 6; ++c) {
      const int fi = as.dofs.free_of(static_cast<int>(ni), c);
      if (fi >= 0) out.node_disp[ni](c) = u(fi);
    }
  return out;
}

}  // namespace keta

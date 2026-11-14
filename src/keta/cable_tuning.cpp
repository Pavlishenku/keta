#include "keta/cable_tuning.h"

#include <algorithm>
#include <cmath>

#include "keta/nonlinear.h"
#include "keta/qp.h"

namespace keta {

Id add_unit_prestress_case(Model& m, Id bar, double P, const std::string& name) {
  const Bar& b = m.bar(bar);
  const double EA = m.material(b.material).E * m.section(b.section).props.A;
  require(EA > 0, "add_unit_prestress_case: the bar needs E A > 0");
  return m.add_case(name.empty() ? "unit prestress " + std::to_string(bar) : name, {Load::prestrain(bar, -P / EA)});
}

double target_value(const Model& model, const CaseResult& r, const TuningTarget& t) {
  require(t.comp >= 0 && t.comp < 6, "tuning target: component 0..5");
  switch (t.kind) {
    case TuningTarget::Kind::Displacement:
      return r.node_disp.at(static_cast<std::size_t>(model.node_index(t.id)))(t.comp);
    case TuningTarget::Kind::Reaction:
      for (const Reaction& q : r.reactions)
        if (q.node == t.id) return q.R(t.comp);
      throw Error("tuning target: no reaction at node " + std::to_string(t.id));
    case TuningTarget::Kind::BarForce: {
      const BarResult& b = r.bars.at(static_cast<std::size_t>(model.bar_index(t.id)));
      require(!b.stations.empty(), "tuning target: the bar has no stations");
      const std::size_t s = t.station < 0 ? b.stations.size() - 1 : static_cast<std::size_t>(t.station);
      require(s < b.stations.size(), "tuning target: station out of range");
      const Station& st = b.stations[s];
      const double v[6] = {st.N, st.Vy, st.Vz, st.T, st.My, st.Mz};
      return v[t.comp];
    }
    case TuningTarget::Kind::LinkForce:
      for (std::size_t i = 0; i < model.links.size(); ++i)
        if (model.links[i].id == t.id) return r.links.at(i).force(t.comp);
      throw Error("tuning target: no link " + std::to_string(t.id));
    case TuningTarget::Kind::CableTension:
      throw Error("tuning target: cable tensions are evaluated by tune_cables");
  }
  return 0;
}

namespace {
// the QP of one tuning step: results r = b + A x
VecX solve_tuning_qp(const MatX& A, const VecX& b, const std::vector<TuningUnknown>& unk,
                     const std::vector<TuningTarget>& tg) {
  const Eigen::Index n = A.cols();
  MatX G = MatX::Zero(n, n);
  VecX a = VecX::Zero(n);
  std::vector<VecX> Erows, Irows;
  std::vector<double> ev, iv;
  bool other = false;
  for (std::size_t i = 0; i < tg.size(); ++i) {
    const VecX ai = A.row(static_cast<Eigen::Index>(i)).transpose();
    const double gap = tg[i].value - b(static_cast<Eigen::Index>(i));
    switch (tg[i].mode) {
      case TuningTarget::Mode::Equal:
        Erows.push_back(ai);
        ev.push_back(gap);
        break;
      case TuningTarget::Mode::AtLeast:
        Irows.push_back(ai);
        iv.push_back(gap);
        other = true;
        break;
      case TuningTarget::Mode::AtMost:
        Irows.push_back(-ai);
        iv.push_back(-gap);
        other = true;
        break;
      case TuningTarget::Mode::Fit:
        require(tg[i].weight > 0, "tuning: a fitted target needs a weight > 0");
        G += 2 * tg[i].weight * ai * ai.transpose();
        a -= 2 * tg[i].weight * gap * ai;
        other = true;
        break;
    }
  }
  for (Eigen::Index j = 0; j < n; ++j) {
    const TuningUnknown& u = unk[static_cast<std::size_t>(j)];
    require(u.lower <= u.upper && u.cost >= 0, "tuning: unknown bounds lower <= upper, cost >= 0");
    if (u.cost > 0) {
      G(j, j) += 2 * u.cost;
      a(j) -= 2 * u.cost * u.ref;
      other = true;
    }
    VecX e = VecX::Zero(n);
    e(j) = 1;
    if (std::isfinite(u.lower)) {
      Irows.push_back(e);
      iv.push_back(u.lower);
      other = true;
    }
    if (std::isfinite(u.upper)) {
      Irows.push_back(-e);
      iv.push_back(-u.upper);
      other = true;
    }
  }
  MatX E(static_cast<Eigen::Index>(Erows.size()), n), I(static_cast<Eigen::Index>(Irows.size()), n);
  for (std::size_t k = 0; k < Erows.size(); ++k) E.row(static_cast<Eigen::Index>(k)) = Erows[k].transpose();
  for (std::size_t k = 0; k < Irows.size(); ++k) I.row(static_cast<Eigen::Index>(k)) = Irows[k].transpose();
  const VecX e = Eigen::Map<const VecX>(ev.data(), static_cast<Eigen::Index>(ev.size()));
  const VecX d = Eigen::Map<const VecX>(iv.data(), static_cast<Eigen::Index>(iv.size()));
  if (!other) {
    // exact: as many independent equal targets as unknowns
    require(E.rows() == n, "tuning: " + std::to_string(E.rows()) + " equal targets for " + std::to_string(n) +
                               " unknowns — add costs (optimised) or fitted targets (best fit)");
    const Eigen::FullPivLU<MatX> lu(E);
    require(lu.isInvertible(), "tuning: the targets do not determine the unknowns (singular influence matrix)");
    return lu.solve(e);
  }
  // augment with the equalities' penalty (zero on the feasible set: the
  // optimum is unchanged) so that G is definite whenever the problem has a
  // unique solution
  if (E.rows() > 0) {
    const MatX EtE = E.transpose() * E;
    const double rho = std::max(G.diagonal().cwiseAbs().maxCoeff(), 1e-300) / EtE.diagonal().maxCoeff();
    G += rho * EtE;
    a -= rho * E.transpose() * e;
  }
  return solve_qp(G, a, E, e, I, d).x;
}
}  // namespace

TuningResult tune_linear(const Model& model, const std::vector<Id>& base_cases,
                         const std::vector<TuningUnknown>& unknowns, const std::vector<TuningTarget>& targets) {
  require(!unknowns.empty() && !targets.empty(), "tune_linear: unknowns and targets needed");
  const StaticResults r = solve_static(model);
  auto res = [&](Id c) -> const CaseResult& {
    for (const CaseResult& x : r.cases)
      if (x.id == c) return x;
    throw Error("tune_linear: no case " + std::to_string(c));
  };
  const Eigen::Index m = static_cast<Eigen::Index>(targets.size()), n = static_cast<Eigen::Index>(unknowns.size());
  VecX b = VecX::Zero(m);
  MatX A(m, n);
  for (Eigen::Index i = 0; i < m; ++i) {
    for (Id c : base_cases) b(i) += target_value(model, res(c), targets[static_cast<std::size_t>(i)]);
    for (Eigen::Index j = 0; j < n; ++j)
      A(i, j) = target_value(model, res(unknowns[static_cast<std::size_t>(j)].case_id), targets[static_cast<std::size_t>(i)]);
  }
  const VecX x = solve_tuning_qp(A, b, unknowns, targets);
  TuningResult out;
  out.x.assign(x.data(), x.data() + x.size());
  const VecX ach = b + A * x;
  out.achieved.assign(ach.data(), ach.data() + ach.size());
  out.base.assign(b.data(), b.data() + b.size());
  out.influence = A;
  std::vector<const CaseResult*> src;
  std::vector<double> f;
  for (Id c : base_cases) {
    src.push_back(&res(c));
    f.push_back(1.0);
  }
  for (Eigen::Index j = 0; j < n; ++j) {
    src.push_back(&res(unknowns[static_cast<std::size_t>(j)].case_id));
    f.push_back(x(j));
  }
  out.state = combine_results(src, f, 0, "tuned");
  return out;
}

TuningResult tune_cables(const Model& model, Id case_id, const std::vector<TuningUnknown>& cables,
                         const std::vector<TuningTarget>& targets, const CableTuningOptions& opt) {
  require(!cables.empty() && !targets.empty(), "tune_cables: cables and targets needed");
  Model w = model;
  std::vector<std::size_t> bi;
  std::vector<double> scale;  // per unknown: the length the finite-difference step is relative to
  const std::vector<double> l0 = resolve_cable_l0(w);
  VecX x(static_cast<Eigen::Index>(cables.size()));
  for (std::size_t j = 0; j < cables.size(); ++j) {
    const TuningUnknown& c = cables[j];
    require(w.bar_index(c.case_id) >= 0, "tune_cables: no bar " + std::to_string(c.case_id));
    const std::size_t k = static_cast<std::size_t>(w.bar_index(c.case_id));
    if (c.kind == TuningUnknown::Kind::CableLength) {
      require(w.bars[k].kind == BarKind::Cable && w.bars[k].cable,
              "tune_cables: unknown " + std::to_string(c.case_id) + " is not a cable bar");
      x(static_cast<Eigen::Index>(j)) = l0[k];
    } else {
      require(w.bars[k].kind != BarKind::Cable,
              "tune_cables: the lack of fit of bar " + std::to_string(c.case_id) +
                  " — a cable's length is its CableLength unknown");
      require(!opt.large_displacements,
              "tune_cables: a lack-of-fit unknown is an imposed strain, which the corotational solver does not take");
      x(static_cast<Eigen::Index>(j)) = 0.0;
    }
    bi.push_back(k);
    scale.push_back(w.bar_chord_length(w.bars[k]));
  }
  for (std::size_t j = 0; j < bi.size(); ++j)
    for (std::size_t i = j + 1; i < bi.size(); ++i)
      require(bi[i] != bi[j], "tune_cables: bar " + std::to_string(cables[j].case_id) + " is an unknown twice");
  bool has_cables = false;
  for (const Bar& b : w.bars) has_cables |= b.kind == BarKind::Cable;
  // the responses of a set of unknowns
  auto evaluate = [&](const VecX& L, CaseResult* state) {
    Model t = w;
    LoadCase& lc = t.cases[static_cast<std::size_t>(&t.load_case(case_id) - t.cases.data())];
    for (std::size_t j = 0; j < bi.size(); ++j) {
      const double xj = L(static_cast<Eigen::Index>(j));
      if (cables[j].kind == TuningUnknown::Kind::CableLength) {
        CableProps cp = *t.bars[bi[j]].cable;
        cp.unstressed_length = xj;
        cp.end_tension = 0;
        t.bars[bi[j]].cable = cp;
      } else if (xj != 0.0) {
        lc.loads.push_back(Load::prestrain(t.bars[bi[j]].id, -xj / scale[j]));
      }
    }
    VecX v(static_cast<Eigen::Index>(targets.size()));
    if (opt.large_displacements) {
      CorotOptions co;
      co.tol_rel = opt.statics.tol_rel;
      co.substeps = opt.substeps;
      const CorotResults cr = solve_corotational(t, case_id, co);
      CaseResult res;
      res.node_disp = cr.node_disp;
      res.reactions = cr.reactions;
      for (std::size_t i = 0; i < targets.size(); ++i) {
        const TuningTarget& tg = targets[i];
        if (tg.kind == TuningTarget::Kind::CableTension) {
          require(tg.comp == 0 || tg.comp == 1, "tuning target: cable end 0 or 1");
          v(static_cast<Eigen::Index>(i)) = cr.cable_T.at(static_cast<std::size_t>(t.bar_index(tg.id)))[static_cast<std::size_t>(tg.comp)];
        } else {
          require(tg.kind == TuningTarget::Kind::Displacement || tg.kind == TuningTarget::Kind::Reaction,
                  "tune_cables: large displacements take displacement, reaction and cable-tension targets");
          v(static_cast<Eigen::Index>(i)) = target_value(t, res, tg);
        }
      }
      if (state) *state = res;
      return v;
    }
    if (!has_cables) {  // lack-of-fit unknowns alone: the linear statics
      StaticOptions so;
      so.stations_per_bar = opt.statics.stations_per_bar;
      Model t1 = t;
      t1.combos.clear();
      std::erase_if(t1.cases, [&](const LoadCase& c) { return c.id != case_id; });
      const StaticResults r = solve_static(t1, so);
      for (std::size_t i = 0; i < targets.size(); ++i) {
        require(targets[i].kind != TuningTarget::Kind::CableTension, "tune_cables: a cable-tension target without cables");
        v(static_cast<Eigen::Index>(i)) = target_value(t1, r.cases[0], targets[i]);
      }
      if (state) *state = r.cases[0];
      return v;
    }
    const CableStaticResults r = solve_cable_static(t, case_id, opt.statics);
    const std::vector<double> l0t = resolve_cable_l0(t);
    for (std::size_t i = 0; i < targets.size(); ++i) {
      const TuningTarget& tg = targets[i];
      if (tg.kind == TuningTarget::Kind::CableTension) {
        const std::size_t k = static_cast<std::size_t>(t.bar_index(tg.id));
        const Bar& b = t.bars[k];
        require(b.kind == BarKind::Cable, "tuning target: bar " + std::to_string(tg.id) + " is not a cable");
        require(tg.comp == 0 || tg.comp == 1, "tuning target: cable end 0 or 1");
        const auto [EA, wgt] = cable_ea_w(t, b);
        const Vec3 pi = t.node(b.ni).x + r.result.node_disp[static_cast<std::size_t>(t.node_index(b.ni))].head<3>();
        const Vec3 pj = t.node(b.nj).x + r.result.node_disp[static_cast<std::size_t>(t.node_index(b.nj))].head<3>();
        const CableEval ev = eval_cable(pi, pj, l0t[k], EA, wgt);
        v(static_cast<Eigen::Index>(i)) = tg.comp == 0 ? ev.f1.norm() : ev.f2.norm();
      } else {
        v(static_cast<Eigen::Index>(i)) = target_value(t, r.result, tg);
      }
    }
    if (state) *state = r.result;
    return v;
  };
  TuningResult out;
  const Eigen::Index m = static_cast<Eigen::Index>(targets.size()), n = x.size();
  VecX v = evaluate(x, nullptr);
  // the scale of each target quantity (a zero target included): its start
  // value, and for a translation the size of the model (its bounding-box
  // diagonal) — a deviation is measured against the structure
  VecX v_start = v.cwiseAbs();
  {
    Vec3 lo = Vec3::Constant(std::numeric_limits<double>::infinity()), hi = -lo;
    for (const Node& nd : w.nodes) {
      lo = lo.cwiseMin(nd.x);
      hi = hi.cwiseMax(nd.x);
    }
    const double diag = w.nodes.empty() ? 0.0 : (hi - lo).norm();
    for (Eigen::Index i = 0; i < m; ++i)
      if (targets[static_cast<std::size_t>(i)].kind == TuningTarget::Kind::Displacement &&
          targets[static_cast<std::size_t>(i)].comp < 3)
        v_start(i) = std::max(v_start(i), diag);
  }
  MatX J(m, n);
  bool step_targets = false;  // fitted or bounded targets: convergence also on the step
  for (const TuningTarget& tg : targets) step_targets |= tg.mode != TuningTarget::Mode::Equal;
  bool any_cost = false;
  for (const TuningUnknown& c : cables) any_cost |= c.cost > 0;
  double last_step = std::numeric_limits<double>::infinity();
  auto describe = [&](const TuningTarget& tg) {
    static const char* kinds[] = {"displacement", "reaction", "bar force", "link force", "cable tension"};
    return std::string(kinds[static_cast<int>(tg.kind)]) + " " + std::to_string(tg.id) + "/" + std::to_string(tg.comp);
  };
  for (int it = 0;; ++it) {
    // converged: every equal target met (and the step settled)
    bool done = it > 0 && (!step_targets || last_step <= opt.tol);
    for (Eigen::Index i = 0; i < m && done; ++i)
      if (targets[static_cast<std::size_t>(i)].mode == TuningTarget::Mode::Equal) {
        const double sc = std::max({std::fabs(targets[static_cast<std::size_t>(i)].value), std::fabs(v(i)),
                                    std::fabs(v_start(i)), 1e-300});
        if (std::fabs(v(i) - targets[static_cast<std::size_t>(i)].value) > opt.tol * sc) done = false;
      }
    if (done) {
      out.iterations = it;
      break;
    }
    require(it < opt.max_iter, "tune_cables: no convergence");
    // central differences of the true coupled solve
    for (Eigen::Index j = 0; j < n; ++j) {
      const double h = cables[static_cast<std::size_t>(j)].kind == TuningUnknown::Kind::CableLength
                           ? 1e-6 * std::fabs(x(j))
                           : 1e-7 * scale[static_cast<std::size_t>(j)];
      VecX xp = x, xm = x;
      xp(j) += h;
      xm(j) -= h;
      J.col(j) = (evaluate(xp, nullptr) - evaluate(xm, nullptr)) / (2 * h);
    }
    // rank of the equal targets (square or over-determined, no cost): a
    // dependent set is named, never resolved by a silent choice
    std::vector<Eigen::Index> eq;
    for (Eigen::Index i = 0; i < m; ++i)
      if (targets[static_cast<std::size_t>(i)].mode == TuningTarget::Mode::Equal) eq.push_back(i);
    if (!any_cost && static_cast<Eigen::Index>(eq.size()) >= n) {
      MatX Je(static_cast<Eigen::Index>(eq.size()), n);
      for (std::size_t r = 0; r < eq.size(); ++r) Je.row(static_cast<Eigen::Index>(r)) = J.row(eq[r]);
      // columns scaled to unit norm: the units of the unknowns do not count
      VecX cn = Je.colwise().norm().transpose();
      for (Eigen::Index j = 0; j < n; ++j) Je.col(j) /= cn(j) > 0 ? cn(j) : 1.0;
      Eigen::JacobiSVD<MatX> svd(Je, Eigen::ComputeFullU | Eigen::ComputeFullV);
      const VecX sv = svd.singularValues();
      if (sv.size() > 0 && !(sv(sv.size() - 1) > 1e-10 * sv(0))) {
        const Eigen::Index k = sv.size() - 1;
        const VecX u = svd.matrixU().col(k), vv = svd.matrixV().col(k);
        std::string tg_names, un_names;
        for (std::size_t r = 0; r < eq.size(); ++r)
          if (std::fabs(u(static_cast<Eigen::Index>(r))) > 0.1)
            tg_names += (tg_names.empty() ? "" : ", ") + describe(targets[static_cast<std::size_t>(eq[r])]);
        for (Eigen::Index j = 0; j < n; ++j)
          if (std::fabs(vv(j)) > 0.1) un_names += (un_names.empty() ? "" : ", ") + std::to_string(cables[static_cast<std::size_t>(j)].case_id);
        fail("tune_cables: the equal targets are rank deficient (singular value ratio " +
             std::to_string(sv(k) / sv(0)) + "): the targets [" + tg_names +
             "] are dependent and the unknowns (bars) [" + un_names + "] are not determined by them");
      }
    }
    // the step: the tuning QP on the linearisation, unknowns shifted to x
    std::vector<TuningUnknown> du = cables;
    for (Eigen::Index j = 0; j < n; ++j) {
      TuningUnknown& u = du[static_cast<std::size_t>(j)];
      u.lower -= x(j);
      u.upper -= x(j);
      u.ref -= x(j);
    }
    const VecX dx = solve_tuning_qp(J, v, du, targets);
    // damped by the true response (halving while the equal-target misfit grows)
    auto misfit = [&](const VecX& vv) {
      double s = 0;
      for (Eigen::Index i = 0; i < m; ++i)
        if (targets[static_cast<std::size_t>(i)].mode == TuningTarget::Mode::Equal)
          s += std::pow(vv(i) - targets[static_cast<std::size_t>(i)].value, 2);
      return s;
    };
    double alpha = 1;
    VecX xn = x + dx, vn = evaluate(xn, nullptr);
    for (int ls = 0; ls < 30 && misfit(vn) > misfit(v) && misfit(v) > 0; ++ls) {
      alpha *= 0.5;
      xn = x + alpha * dx;
      vn = evaluate(xn, nullptr);
    }
    double xs = 0;
    for (Eigen::Index j = 0; j < n; ++j) xs = std::max(xs, scale[static_cast<std::size_t>(j)]);
    last_step = (xn - x).lpNorm<Eigen::Infinity>() / std::max(xs, 1e-300);
    x = xn;
    v = vn;
  }
  evaluate(x, &out.state);
  out.x.assign(x.data(), x.data() + n);
  out.achieved.assign(v.data(), v.data() + m);
  out.influence = J;
  return out;
}

}  // namespace keta

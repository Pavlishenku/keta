#include "keta/fiber_frame.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "keta/assembly.h"
#include "keta/fiber_elem.h"
#include "keta/nonlinear.h"  // rot_exp / rot_log

namespace keta {

using namespace keta::fiberel;

void gauss_lobatto01(int n, std::vector<double>& x, std::vector<double>& w) {
  switch (n) {
    case 3:
      x = {0.0, 0.5, 1.0};
      w = {1.0 / 6, 2.0 / 3, 1.0 / 6};
      break;
    case 4: {
      double a = 1.0 / std::sqrt(5.0);
      x = {0.0, 0.5 * (1 - a), 0.5 * (1 + a), 1.0};
      w = {1.0 / 12, 5.0 / 12, 5.0 / 12, 1.0 / 12};
      break;
    }
    case 5: {
      double a = std::sqrt(3.0 / 7.0);
      x = {0.0, 0.5 * (1 - a), 0.5, 0.5 * (1 + a), 1.0};
      w = {1.0 / 20, 49.0 / 180, 16.0 / 45, 49.0 / 180, 1.0 / 20};
      break;
    }
    case 6: {
      double s7 = std::sqrt(7.0);
      double a = std::sqrt((7.0 - 2.0 * s7) / 21.0), b = std::sqrt((7.0 + 2.0 * s7) / 21.0);
      x = {0.0, 0.5 * (1 - b), 0.5 * (1 - a), 0.5 * (1 + a), 0.5 * (1 + b), 1.0};
      w = {1.0 / 30, (14.0 - s7) / 60, (14.0 + s7) / 60,
           (14.0 + s7) / 60, (14.0 - s7) / 60, 1.0 / 30};
      break;
    }
    default: fail("fiber_frame: n_ip must be 3, 4, 5 or 6");
  }
}

namespace {

using Vec5 = Eigen::Matrix<double, 5, 1>;
struct BarInfo3 {
  bool fiber = false;
  double L0 = 0;
  Mat612 Blin = Mat612::Zero();  // a6 * blkdiag(R,R,R,R): linear basic map
  Mat12 ke = Mat12::Zero();      // geometrically linear elastic bars
  Mat6 kb6 = Mat6::Zero();       // corotational elastic bars: basic stiffness
  Mat3 E0 = Mat3::Identity();    // columns = initial local axes (global frame)
  Vec3 Xi = Vec3::Zero(), Xj = Vec3::Zero();
  int gi[12] = {0};
};

// Corotational basic deformation extraction — see the header block.
Vec6 corot_v(const BarInfo3& bi, const Vec12& ue) {
  Vec3 ui = ue.segment<3>(0), thi = ue.segment<3>(3);
  Vec3 uj = ue.segment<3>(6), thj = ue.segment<3>(9);
  Vec3 d0 = bi.Xj - bi.Xi;
  Vec3 dv = d0 + (uj - ui);
  double ln = dv.norm();
  require(ln > 1e-9 * bi.L0, "fiber_frame: corotational element collapsed to zero length");
  // Cancellation-safe elongation: (ln^2 - L0^2) / (ln + L0).
  double delta = ((uj - ui).dot(dv + d0)) / (ln + bi.L0);
  Mat3 Ri = rot_exp(thi), Rj = rot_exp(thj);
  Vec3 xhat = dv / ln;
  Vec3 ym = 0.5 * (Ri * bi.E0.col(1) + Rj * bi.E0.col(1));
  Vec3 zhat = xhat.cross(ym);
  double zn = zhat.norm();
  require(zn > 1e-8, "fiber_frame: degenerate corotational frame (chord parallel to mean y)");
  zhat /= zn;
  Vec3 yhat = zhat.cross(xhat);
  Mat3 E;
  E.col(0) = xhat;
  E.col(1) = yhat;
  E.col(2) = zhat;
  Vec3 Ti = rot_log(E.transpose() * Ri * bi.E0);
  Vec3 Tj = rot_log(E.transpose() * Rj * bi.E0);
  Vec6 v;
  v << delta, Ti(1), Tj(1), Ti(2), Tj(2), Tj(0) - Ti(0);
  return v;
}

// B = dv/du by 4th-order central differences (see header: FD is the house
// pattern for corotational kinematics; h = 1e-5 * scale keeps the roundoff on
// the equilibrium residual ~1e-11 relative).
Mat612 corot_B(const BarInfo3& bi, const Vec12& ue) {
  Mat612 B;
  for (int j = 0; j < 12; ++j) {
    double h = (j % 6) < 3 ? 1e-5 * std::max(1.0, bi.L0) : 1e-5;
    Vec12 up = ue, um = ue, up2 = ue, um2 = ue;
    up(j) += h;
    um(j) -= h;
    up2(j) += 2 * h;
    um2(j) -= 2 * h;
    B.col(j) =
        (8.0 * (corot_v(bi, up) - corot_v(bi, um)) - (corot_v(bi, up2) - corot_v(bi, um2))) /
        (12.0 * h);
  }
  return B;
}

// Geometric part of the tangent: Hessian of g(u) = q . v(u) at frozen q, by
// central second differences (tangent-only accuracy).
Mat12 corot_geom(const BarInfo3& bi, const Vec12& ue, const Vec6& q6) {
  auto g = [&](const Vec12& u) { return q6.dot(corot_v(bi, u)); };
  double g0 = g(ue);
  double h[12];
  for (int j = 0; j < 12; ++j) h[j] = (j % 6) < 3 ? 1e-4 * std::max(1.0, bi.L0) : 1e-4;
  Mat12 H;
  for (int j = 0; j < 12; ++j) {
    Vec12 up = ue, um = ue;
    up(j) += h[j];
    um(j) -= h[j];
    H(j, j) = (g(up) - 2 * g0 + g(um)) / (h[j] * h[j]);
    for (int k = j + 1; k < 12; ++k) {
      Vec12 upp = ue, upm = ue, ump = ue, umm = ue;
      upp(j) += h[j];
      upp(k) += h[k];
      upm(j) += h[j];
      upm(k) -= h[k];
      ump(j) -= h[j];
      ump(k) += h[k];
      umm(j) -= h[j];
      umm(k) -= h[k];
      double v = (g(upp) - g(upm) - g(ump) + g(umm)) / (4 * h[j] * h[k]);
      H(j, k) = v;
      H(k, j) = v;
    }
  }
  return H;
}

}  // namespace

// ---------------- spatial driver ----------------

FiberFrame3DResult solve_fiber_frame3d(const Model& model, Id case_id,
                                       const FiberFrame3DOptions& opt) {
  require(model.continua.empty(), "solve_fiber_frame3d: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  const bool staged = !opt.stages.empty();
  require(opt.drying.empty() || staged, "fiber_frame: drying runs in staged construction (stages with durations)");
  const int nn = static_cast<int>(model.nodes.size());
  const int ndof = 6 * nn;
  const bool corot = opt.corotational;
  require(model.shells.empty(), "fiber_frame: shells are not supported");

  Assembled asm_ = assemble_elements(model);

  // ---- supports ----
  std::vector<bool> blocked(static_cast<std::size_t>(ndof), false);
  std::vector<double> spring(static_cast<std::size_t>(ndof), 0.0);
  for (const Support& sp : model.supports) {
    int ix = model.node_index(sp.node);
    require(!sp.spring6 && sp.decol_axis < 0 && sp.rot == 0.0,
            "fiber_frame: only plain fixed / spring supports are supported");
    for (int c = 0; c < 6; ++c) {
      if (sp.fixed[static_cast<std::size_t>(c)]) blocked[static_cast<std::size_t>(6 * ix + c)] = true;
      spring[static_cast<std::size_t>(6 * ix + c)] += sp.spring[static_cast<std::size_t>(c)];
    }
  }

  // ---- per-bar setup ----
  std::vector<BarInfo3> info(model.bars.size());
  std::vector<FiberElement3D> elems(model.bars.size());
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    BarInfo3& bi = info[e];
    const BarElement& el = asm_.elements[e];
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (int c = 0; c < 6; ++c) {
      bi.gi[c] = 6 * ni + c;
      bi.gi[6 + c] = 6 * nj + c;
    }
    bi.L0 = el.L();
    bi.Xi = model.nodes[static_cast<std::size_t>(ni)].x;
    bi.Xj = model.nodes[static_cast<std::size_t>(nj)].x;
    bi.E0 = el.R().transpose();  // columns = local axes
    auto it = opt.sections.find(b.id);
    bi.fiber = (it != opt.sections.end());
    Mat12 T = Mat12::Zero();
    for (int blk = 0; blk < 4; ++blk) T.block<3, 3>(3 * blk, 3 * blk) = el.R();
    if (bi.fiber) {
      require(b.kind == BarKind::Frame, "fiber_frame: a fiber bar must be a Frame bar");
      require(!b.has_ecc(), "fiber_frame: eccentricity on a fiber bar is not supported");
      require(!b.soil.any(), "fiber_frame: a Winkler bed on a fiber bar is not supported");
      for (bool r : b.releases)
        require(!r, "fiber_frame: releases on a fiber bar are not supported");
      require(!model.tapered_of(b.id), "fiber_frame: a tapered fiber bar is not supported");
      double GJ = model.material(b.material).G() * model.section(b.section).props.J;
      bi.Blin = basic_transform6(bi.L0) * T;
      elems[e].init(it->second, bi.L0, GJ, opt.n_ip, opt.tangent_floor, opt.max_element_iter,
                    opt.tol_element);
    } else if (!corot) {
      bi.ke = el.ke();  // releases / ecc / soil / taper folded in
      if (model.tapered_of(b.id)) bi.ke = bar_ke_assembled(model, el);
    } else {
      require(b.kind != BarKind::Cable, "fiber_frame: cable bars in corotational mode");
      require(!b.has_ecc(), "fiber_frame: eccentric bars in corotational mode");
      require(!b.soil.any(), "fiber_frame: Winkler beds in corotational mode");
      require(!model.tapered_of(b.id), "fiber_frame: tapered bars in corotational mode");
      bi.Blin = basic_transform6(bi.L0) * T;
      Mat12 kl = el.k_local();
      // Rotate the LOCAL condensed matrix into the basic system.
      bi.kb6 = elastic_basic6(kl, bi.L0);
    }
  }

  // ---- external load pattern ----
  // `weigh`: the bars a SelfWeight load weighs (null: all)
  auto build_loads = [&](const LoadCase& lc, const std::vector<char>* weigh) {
  VecX fext = VecX::Zero(ndof);
  for (const Load& l : lc.loads) {
    switch (l.kind) {
      case Load::Kind::NodalForce: {
        int ix = model.node_index(l.node);
        for (int c = 0; c < 3; ++c) {
          fext(6 * ix + c) += l.F(c);
          fext(6 * ix + 3 + c) += l.M(c);
        }
        break;
      }
      case Load::Kind::SelfWeight: {
        require(!corot, "fiber_frame: corotational mode takes nodal loads only");
        for (std::size_t e = 0; e < model.bars.size(); ++e) {
          if (weigh && !(*weigh)[e]) continue;
          require(!info[e].fiber, "fiber_frame: span loads on a fiber bar are not supported");
          const Bar& b = model.bars[e];
          Vec12 f = asm_.elements[e].equivalent_nodal(l, model.material(b.material));
          for (int c = 0; c < 6; ++c) {
            fext(info[e].gi[c]) += f(c);
            fext(info[e].gi[6 + c]) += f(6 + c);
          }
        }
        break;
      }
      case Load::Kind::SupportDisp:
        fail("fiber_frame: prescribed support displacements are not supported");
        break;
      case Load::Kind::SurfacePressure:
      case Load::Kind::SurfaceThermal:
      case Load::Kind::Hydrostatic:
        fail("fiber_frame: surface loads are not supported");
        break;
      default: {
        require(!corot, "fiber_frame: corotational mode takes nodal loads only");
        std::size_t e = static_cast<std::size_t>(model.bar_index(l.bar));
        require(!info[e].fiber, "fiber_frame: span loads on a fiber bar are not supported");
        const Bar& b = model.bars[e];
        Vec12 f = asm_.elements[e].equivalent_nodal(l, model.material(b.material));
        for (int c = 0; c < 6; ++c) {
          fext(info[e].gi[c]) += f(c);
          fext(info[e].gi[6 + c]) += f(6 + c);
        }
        break;
      }
    }
  }
  for (int i = 0; i < ndof; ++i)
    if (blocked[static_cast<std::size_t>(i)]) fext(i) = 0;
  return fext;
  };
  VecX fext = staged ? VecX(VecX::Zero(ndof)) : build_loads(model.load_case(case_id), nullptr);
  VecX fbase = VecX::Zero(ndof);  // staged: the loads of the completed stages
  std::vector<char> active(model.bars.size(), staged ? 0 : 1);
  std::vector<Vec12> ue_act(model.bars.size(), Vec12::Zero());  // displacements at activation
  std::vector<Vec6> v_act(model.bars.size(), Vec6::Zero());     // corotational: basic deformations then

  // ---- residual / tangent at a trial displacement ----
  FiberFrame3DResult out;
  auto evaluate = [&](const VecX& u, MatX& K, VecX& fint) {
    K = MatX::Zero(ndof, ndof);
    fint = VecX::Zero(ndof);
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      if (!active[e]) continue;
      BarInfo3& bi = info[e];
      Vec12 ue;
      for (int c = 0; c < 12; ++c) ue(c) = u(bi.gi[c]);
      Vec12 fe;
      Mat12 kee;
      if (!bi.fiber && !corot) {
        kee = bi.ke;
        fe = bi.ke * (ue - ue_act[e]);
      } else {
        Vec6 v;
        Mat612 B;
        if (corot) {
          v = corot_v(bi, ue) - v_act[e];
          B = corot_B(bi, ue);
        } else {
          v = bi.Blin * (ue - ue_act[e]);
          B = bi.Blin;
        }
        Vec6 q6;
        Mat6 kb6;
        if (bi.fiber) {
          elems[e].determine(v);
          q6 = elems[e].q6();
          kb6 = elems[e].kb;
        } else {
          q6 = bi.kb6 * v;
          kb6 = bi.kb6;
        }
        fe = B.transpose() * q6;
        kee = B.transpose() * kb6 * B;
        if (corot) kee += corot_geom(bi, ue, q6);
      }
      for (int r = 0; r < 12; ++r) {
        fint(bi.gi[r]) += fe(r);
        for (int c = 0; c < 12; ++c) K(bi.gi[r], bi.gi[c]) += kee(r, c);
      }
    }
    for (int i = 0; i < ndof; ++i) {
      if (spring[static_cast<std::size_t>(i)] != 0) {
        K(i, i) += spring[static_cast<std::size_t>(i)];
        fint(i) += spring[static_cast<std::size_t>(i)] * u(i);
      }
    }
    for (int i = 0; i < ndof; ++i)
      if (blocked[static_cast<std::size_t>(i)]) {
        K.row(i).setZero();
        K.col(i).setZero();
        K(i, i) = 1.0;
        fint(i) = 0;
      }
    // a DOF that received no stiffness at all (isolated component) is fixed
    for (int i = 0; i < ndof; ++i)
      if (!blocked[static_cast<std::size_t>(i)] && K(i, i) == 0.0) {
        K(i, i) = 1.0;
        fint(i) = 0;
      }
  };

  std::vector<Vec5> q_saved(model.bars.size(), Vec5::Zero());
  auto save_q = [&] {
    for (std::size_t e = 0; e < model.bars.size(); ++e) q_saved[e] = elems[e].q;
  };
  auto restore_q = [&] {
    for (std::size_t e = 0; e < model.bars.size(); ++e) elems[e].q = q_saved[e];
  };
  auto commit_all = [&] {
    for (std::size_t e = 0; e < model.bars.size(); ++e)
      if (info[e].fiber) elems[e].commit_states();
  };

  VecX u = VecX::Zero(ndof);
  double lambda = 0;
  MatX K;
  VecX fint;
  double fref = fext.norm();
  auto zero_blocked = [&](VecX& r) {
    for (int i = 0; i < ndof; ++i)
      if (blocked[static_cast<std::size_t>(i)]) r(i) = 0;
  };
  auto record = [&](double u_ctrl) {
    FiberTracePoint p;
    p.lambda = lambda;
    p.u_ctrl = u_ctrl;
    p.node_disp.resize(model.nodes.size());
    for (int i = 0; i < nn; ++i) p.node_disp[static_cast<std::size_t>(i)] = u.segment<6>(6 * i);
    out.trace.push_back(std::move(p));
  };

  // One load-controlled Newton step to `lam_target` (line search + q restore
  // on backtracks; commits on success). Throws NonConvergence on failure.
  auto newton_to = [&](double lam_target, double res_scale) {
    VecX target = fbase + lam_target * fext;
    save_q();
    evaluate(u, K, fint);
    VecX res = target - fint;
    zero_blocked(res);
    double r0 = res.norm();
    // (scaled by the internal forces too: a step driven by imposed strains
    // alone — restrained shrinkage — has no external load to measure by)
    const double scale = std::max(res_scale, fint.norm());
    bool converged = r0 <= opt.tol_force * scale;
    for (int it = 0; it < opt.max_global_iter && !converged; ++it) {
      ++out.iterations;
      Eigen::FullPivLU<MatX> lu(K);
      if (!lu.isInvertible())
        throw NonConvergence("fiber_frame: singular tangent at load factor " +
                             std::to_string(lam_target));
      VecX du = lu.solve(res);
      if (!du.allFinite())
        throw NonConvergence("fiber_frame: non-finite Newton step at load factor " +
                             std::to_string(lam_target));
      double alpha = 1.0;
      bool accepted = false;
      for (int back = 0; back <= opt.max_backtracks; ++back) {
        VecX ut = u + alpha * du;
        restore_q();
        try {
          evaluate(ut, K, fint);
        } catch (const Error&) {
          alpha *= 0.5;
          continue;
        }
        VecX rt = target - fint;
        zero_blocked(rt);
        double rn = rt.norm();
        if (rn < r0 || rn <= opt.tol_force * scale) {
          u = ut;
          res = rt;
          r0 = rn;
          accepted = true;
          save_q();
          break;
        }
        alpha *= 0.5;
      }
      if (!accepted)
        throw NonConvergence("fiber_frame: line search failed at load factor " +
                             std::to_string(lam_target));
      converged = r0 <= opt.tol_force * scale;
    }
    if (!converged)
      throw NonConvergence("fiber_frame: no equilibrium at load factor " +
                           std::to_string(lam_target));
    lambda = lam_target;
    commit_all();
  };

  // ---- control dispatch ----
  if (staged) {
    require(opt.control == FiberFrame3DOptions::Control::Load, "fiber_frame: staged construction runs under load control");
    std::vector<char> in_place(static_cast<std::size_t>(nn), 0);
    // in place from the start: the nodes whose three translations are held
    for (const Support& sp : model.supports)
      if (sp.fixed[0] && sp.fixed[1] && sp.fixed[2]) in_place[static_cast<std::size_t>(model.node_index(sp.node))] = 1;
    record(0.0);
    double clock = 0;
    out.trace_time.push_back(clock);
    // the step duration seen by every committed fiber state (creep)
    std::vector<double> born(model.bars.size(), 0.0);  // the clock time of age 0 of each bar's concrete
    // drying bars: their humidity fields (from their activation)
    std::map<std::size_t, DryingField> fields;
    for (const auto& [bid, ds] : opt.drying) {
      const std::size_t e = static_cast<std::size_t>(model.bar_index(bid));
      require(info[e].fiber, "fiber_frame: drying of bar " + std::to_string(bid) + ", not a fiber bar");
      const HystFiberSection& hs = opt.sections.at(bid);
      (void)drying_init(ds);  // (validation)
      require(static_cast<int>(hs.fibers.size()) >= ds.ny * ds.nz, "fiber_frame: the drying grid exceeds the section's fibers");
      for (int i = 0; i < ds.ny; ++i)
        for (int j = 0; j < ds.nz; ++j) {
          const Fiber& fb = hs.fibers[static_cast<std::size_t>(i * ds.nz + j)];
          const double y = -ds.b / 2 + ds.b * (i + 0.5) / ds.ny, z = -ds.h / 2 + ds.h * (j + 0.5) / ds.nz;
          require(std::fabs(fb.y - y) <= 1e-9 * ds.b && std::fabs(fb.z - z) <= 1e-9 * ds.h,
                  "fiber_frame: the section's first fibers are not the drying grid (HystFiberSection::rect order)");
          require(hs.laws[static_cast<std::size_t>(fb.law)].shrink.empty(),
                  "fiber_frame: a drying fiber's law carries a shrinkage table too");
        }
    }
    // the drying step of every active drying bar, its fibers' shrinkage at the step's end
    auto dry = [&](double dt) {
      for (auto& [e, fld] : fields) {
        const DryingSection& ds = opt.drying.at(model.bars[e].id);
        drying_step(ds, fld, dt);
        for (auto& pt : elems[e].committed)
          for (int c = 0; c < ds.ny * ds.nz; ++c) {
            HystState& h = pt[static_cast<std::size_t>(c)];
            h.sh_mode = 1;
            h.sh_drv = ds.shrinkage(fld.hum[static_cast<std::size_t>(c)]);
          }
      }
    };
    auto set_dt = [&](double dt) {
      for (std::size_t e = 0; e < model.bars.size(); ++e)
        if (info[e].fiber)
          for (auto& pt : elems[e].committed)
            for (HystState& h : pt) {
              h.cr_dt = dt;
              h.age = clock - born[e];
            }
    };
    for (const FiberFrame3DOptions::Stage& st : opt.stages) {
      require(st.load_steps >= 1, "fiber_frame: a stage needs load_steps >= 1");
      std::vector<std::size_t> fresh;
      std::vector<char> weigh(model.bars.size(), 0);
      for (Id b : st.activate) {
        const std::size_t e = static_cast<std::size_t>(model.bar_index(b));
        require(!active[e], "fiber_frame: bar " + std::to_string(b) + " activated twice");
        fresh.push_back(e);
        weigh[e] = 1;
      }
      // place the new nodes by rigid continuation through the new bars
      std::vector<char> todo(fresh.size(), 1);
      for (bool grew = true; grew;) {
        grew = false;
        for (std::size_t q = 0; q < fresh.size(); ++q) {
          if (!todo[q]) continue;
          const Bar& b = model.bars[fresh[q]];
          const int i = model.node_index(b.ni), j = model.node_index(b.nj);
          const bool pi = in_place[static_cast<std::size_t>(i)], pj = in_place[static_cast<std::size_t>(j)];
          if (pi && pj) {
            todo[q] = 0;
            continue;
          }
          if (!pi && !pj) continue;
          const int from = pi ? i : j, to = pi ? j : i;
          const Vec3 th = u.segment<3>(6 * from + 3);
          const Vec3 X = model.nodes[static_cast<std::size_t>(to)].x - model.nodes[static_cast<std::size_t>(from)].x;
          u.segment<3>(6 * to) = u.segment<3>(6 * from) + th.cross(X);
          u.segment<3>(6 * to + 3) = th;
          for (int c = 0; c < 6; ++c)
            if (blocked[static_cast<std::size_t>(6 * to + c)]) u(6 * to + c) = 0;
          in_place[static_cast<std::size_t>(to)] = 1;
          todo[q] = 0;
          grew = true;
        }
      }
      for (char t : todo) require(!t, "fiber_frame: a stage activates a bar with neither end in place");
      // stress-free activation: the deformation counts from now
      for (std::size_t e : fresh) {
        born[e] = clock - st.activation_age;
        if (opt.drying.count(model.bars[e].id)) fields[e] = drying_init(opt.drying.at(model.bars[e].id));
        for (auto& pt : elems[e].committed)
          for (HystState& h : pt) h.age0 = st.activation_age;
        active[e] = 1;
        for (int c = 0; c < 12; ++c) ue_act[e](c) = u(info[e].gi[c]);
        if (corot) v_act[e] = corot_v(info[e], ue_act[e]);
      }
      fbase += fext;
      fext = st.case_id != 0 ? build_loads(model.load_case(st.case_id), &weigh) : VecX(VecX::Zero(ndof));
      if (st.case_id != 0)
        for (const Load& l : model.load_case(st.case_id).loads)
          if (l.kind == Load::Kind::NodalForce)
            require(in_place[static_cast<std::size_t>(model.node_index(l.node))], "fiber_frame: a nodal load on a node not in place");
      const double fscale = std::max({fbase.norm(), fext.norm(), 1.0});
      for (int step = 1; step <= st.load_steps; ++step) {
        newton_to(static_cast<double>(step) / st.load_steps, fscale);
        ++out.steps;
        record(0.0);
        out.trace_time.push_back(clock);
      }
      if (st.duration > 0) {
        require(st.time_steps >= 1, "fiber_frame: a stage's creep needs time_steps >= 1");
        const double dt = st.duration / st.time_steps;
        for (int k = 1; k <= st.time_steps; ++k) {
          set_dt(dt);  // (the ages advance with the clock)
          dry(dt);
          newton_to(1.0, fscale);
          ++out.steps;
          clock += dt;
          record(0.0);
          out.trace_time.push_back(clock);
        }
        set_dt(0);
      } else {
        require(st.duration == 0, "fiber_frame: a stage's duration must be >= 0");
      }
      out.stage_end.push_back(out.trace.size() - 1);
    }
    for (const auto& [e, fld] : fields) out.humidity[model.bars[e].id] = fld.hum;
  } else if (opt.control == FiberFrame3DOptions::Control::Load) {
    require(opt.load_steps >= 1, "fiber_frame: load_steps must be >= 1");
    double fscale = fref > 0 ? fref : 1.0;
    for (int step = 1; step <= opt.load_steps; ++step) {
      double lam = static_cast<double>(step) / opt.load_steps;
      newton_to(lam, lam * fscale);
      out.steps = step;
      record(0.0);
    }
  } else if (opt.control == FiberFrame3DOptions::Control::Displacement) {
    require(fref > 0, "fiber_frame: displacement control needs a non-zero load pattern");
    require(!opt.ctrl_targets.empty(), "fiber_frame: displacement control needs targets");
    require(opt.ctrl_comp >= 0 && opt.ctrl_comp < 6, "fiber_frame: ctrl_comp must be 0..5");
    require(opt.steps_per_leg >= 1, "fiber_frame: steps_per_leg must be >= 1");
    int cdof = 6 * model.node_index(opt.ctrl_node) + opt.ctrl_comp;
    require(!blocked[static_cast<std::size_t>(cdof)],
            "fiber_frame: the controlled DOF is blocked by a support");
    double lam_hist = 0;
    // One increment of the controlled DOF to `tgt` (bordered Newton).
    std::function<void(double, int)> disp_step = [&](double tgt, int depth) {
      VecX u_in = u;
      double lam_in = lambda;
      save_q();
      try {
        for (int it = 0; it < opt.max_global_iter; ++it) {
          evaluate(u, K, fint);
          VecX r = lambda * fext - fint;
          zero_blocked(r);
          double scale = fref * std::max({std::fabs(lambda), 0.05 * lam_hist, 1e-9});
          if (r.norm() <= opt.tol_force * scale &&
              std::fabs(u(cdof) - tgt) <= 1e-10 * std::max(1.0, std::fabs(tgt))) {
            commit_all();
            lam_hist = std::max(lam_hist, std::fabs(lambda));
            return;
          }
          ++out.iterations;
          // Near-singular tangents are legal here (the bordered system stays
          // regular through load limit points): solve and check finiteness.
          Eigen::FullPivLU<MatX> lu(K);
          VecX du1 = lu.solve(r);
          VecX du2 = lu.solve(fext);
          if (!du1.allFinite() || !du2.allFinite())
            throw NonConvergence("fiber_frame: singular tangent (displacement control)");
          require(std::fabs(du2(cdof)) > 1e-300,
                  "fiber_frame: the load pattern does not move the controlled DOF");
          double dlam = (tgt - u(cdof) - du1(cdof)) / du2(cdof);
          u += du1 + dlam * du2;
          lambda += dlam;
          if (!u.allFinite() || !std::isfinite(lambda))
            throw NonConvergence("fiber_frame: non-finite step (displacement control)");
        }
        throw NonConvergence("fiber_frame: no equilibrium (displacement control)");
      } catch (const Error&) {
        // cut the increment in half (restore, do the two halves)
        if (depth >= opt.max_cuts) throw;
        u = u_in;
        lambda = lam_in;
        restore_q();
        double mid = 0.5 * (u_in(cdof) + tgt);
        disp_step(mid, depth + 1);
        disp_step(tgt, depth + 1);
      }
    };
    record(0.0);
    double prev = 0;
    for (double target : opt.ctrl_targets) {
      for (int s = 1; s <= opt.steps_per_leg; ++s) {
        double tgt = prev + (target - prev) * s / opt.steps_per_leg;
        disp_step(tgt, 0);
        ++out.steps;
        record(u(cdof));
      }
      prev = target;
    }
  } else {  // Control::ArcLength
    require(fref > 0, "fiber_frame: arc length needs a non-zero load pattern");
    require(opt.dl0 > 0, "fiber_frame: dl0 must be > 0");
    double dl = opt.dl0;
    VecX du_prev = VecX::Zero(ndof);
    record(0.0);
    for (int step = 0; step < opt.max_steps; ++step) {
      bool converged = false;
      VecX u_new;
      double lam_new = 0;
      VecX Du;
      VecX u0 = u;
      double lam0 = lambda;
      save_q();
      for (int cut = 0; cut <= opt.max_cuts && !converged; ++cut) {
        restore_q();
        u = u0;
        lambda = lam0;
        try {
          evaluate(u, K, fint);
          // Near the limit point the tangent is near-singular BY DESIGN: solve
          // anyway and let the cylindrical constraint tame the large solution
          // (an exactly singular factorization yields non-finite values, which
          // are caught below). Gating on isInvertible() here would deadlock the
          // continuation exactly at the point it exists to cross.
          Eigen::FullPivLU<MatX> lu(K);
          VecX duF = lu.solve(fext);
          if (!duF.allFinite()) throw NonConvergence("fiber_frame: singular arc tangent");
          double dlam = dl / duF.norm();
          if (step > 0 && du_prev.dot(duF) < 0) dlam = -dlam;
          double lam = lambda + dlam;
          VecX uu = u + dlam * duF;
          Du = uu - u;
          bool ok = false;
          for (int it = 0; it < opt.max_arc_iterations; ++it) {
            ++out.iterations;
            evaluate(uu, K, fint);
            VecX r = lam * fext - fint;
            zero_blocked(r);
            double scale = std::max({std::fabs(lam) * fref, fint.norm(), 1e-30});
            if (r.norm() <= opt.tol_force * scale) {
              ok = true;
              break;
            }
            Eigen::FullPivLU<MatX> lu2(K);
            VecX dur = lu2.solve(r);
            VecX duF2 = lu2.solve(fext);
            if (!dur.allFinite() || !duF2.allFinite()) break;
            VecX t = Du + dur;
            double a = duF2.dot(duF2);
            double b = 2 * duF2.dot(t);
            double c = t.dot(t) - dl * dl;
            double disc = b * b - 4 * a * c;
            if (disc < 0) break;  // step cut
            double r1 = (-b + std::sqrt(disc)) / (2 * a);
            double r2 = (-b - std::sqrt(disc)) / (2 * a);
            // Crisfield root selection: maximize dot(Du_new, Du).
            VecX c1 = t + r1 * duF2, c2 = t + r2 * duF2;
            double pick = c1.dot(Du) >= c2.dot(Du) ? r1 : r2;
            uu += dur + pick * duF2;
            lam += pick;
            Du = uu - u;
          }
          if (ok) {
            u_new = uu;
            lam_new = lam;
            converged = true;
          } else {
            dl /= 2;
          }
        } catch (const Error&) {
          dl /= 2;
        }
      }
      if (!converged) throw NonConvergence("fiber_frame: arc-length step failed after cuts");
      du_prev = u_new - u0;
      u = u_new;
      lambda = lam_new;
      commit_all();
      ++out.steps;
      record(0.0);
      dl = std::clamp(dl * 1.2, opt.dl0 / 16, 4 * opt.dl0);
      if (opt.stop_lambda > 0 && std::fabs(lambda) >= opt.stop_lambda) break;
      if (opt.stop_disp_inf > 0 && u.lpNorm<Eigen::Infinity>() >= opt.stop_disp_inf) break;
    }
  }

  // ---- results ----
  out.lambda = lambda;
  out.node_disp.assign(model.nodes.size(), Vec6::Zero());
  for (int i = 0; i < nn; ++i) out.node_disp[static_cast<std::size_t>(i)] = u.segment<6>(6 * i);
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    FiberBarState3D st;
    st.bar = model.bars[e].id;
    st.fiber = info[e].fiber;
    if (st.fiber) {
      st.q = elems[e].q6();
      for (std::size_t k = 0; k < elems[e].xg.size(); ++k) {
        FiberStation3D fs;
        fs.x = elems[e].xg[k];
        fs.weight = elems[e].wg[k];
        fs.N = elems[e].s[k](0);
        fs.My = elems[e].s[k](1);
        fs.Mz = elems[e].s[k](2);
        fs.eps0 = elems[e].d[k](0);
        fs.kappa_y = elems[e].d[k](1);
        fs.kappa_z = elems[e].d[k](2);
        st.stations.push_back(fs);
      }
    }
    out.bars.push_back(std::move(st));
  }
  return out;
}

// ---------------- plane compatibility wrapper ----------------

FiberFrameResult solve_fiber_frame(const Model& model, Id case_id, const FiberFrameOptions& opt) {
  require(model.continua.empty(), "solve_fiber_frame: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  require(opt.load_steps >= 1, "fiber_frame: load_steps must be >= 1");
  const LoadCase& lc = model.load_case(case_id);

  // ---- plane-frame gate ----
  double extent = 0;
  for (const Node& n : model.nodes) extent = std::max(extent, n.x.norm());
  double ytol = 1e-9 * std::max(1.0, extent);
  for (const Node& n : model.nodes)
    require(std::fabs(n.x.y()) <= ytol,
            "fiber_frame: node " + std::to_string(n.id) + " is out of the global X-Z plane");

  // ---- plane load gate (out-of-plane components are typed errors) ----
  for (const Load& l : lc.loads) {
    switch (l.kind) {
      case Load::Kind::NodalForce:
        require(l.F.y() == 0 && l.M.x() == 0 && l.M.z() == 0,
                "fiber_frame: out-of-plane nodal load component");
        break;
      case Load::Kind::SelfWeight:
        require(l.g.y() == 0, "fiber_frame: out-of-plane gravity");
        break;
      case Load::Kind::SupportDisp:
        fail("fiber_frame: prescribed support displacements are not supported");
        break;
      default:
        require(l.F.y() == 0 && l.F2.y() == 0 && l.M.x() == 0 && l.M.z() == 0,
                "fiber_frame: out-of-plane span load component");
        break;
    }
  }

  // ---- plane fiber-bar gate ----
  for (const Bar& b : model.bars) {
    if (opt.sections.find(b.id) == opt.sections.end()) continue;
    Mat3 R = bar_rotation(model.node(b.ni).x, model.node(b.nj).x, b.roll);
    require(std::fabs(std::fabs(R(1, 1)) - 1.0) < 1e-9,
            "fiber_frame: bar " + std::to_string(b.id) +
                " has its local y out of the global Y axis (out-of-plane bar or roll != 0)");
  }

  FiberFrame3DOptions o3;
  for (const auto& [id, fs] : opt.sections) o3.sections.emplace(id, HystFiberSection::wrap(fs));
  o3.n_ip = opt.n_ip;
  o3.corotational = false;
  o3.control = FiberFrame3DOptions::Control::Load;
  o3.load_steps = opt.load_steps;
  o3.max_global_iter = opt.max_global_iter;
  o3.max_element_iter = opt.max_element_iter;
  o3.max_backtracks = opt.max_backtracks;
  o3.tol_force = opt.tol_force;
  o3.tol_element = opt.tol_element;
  o3.tangent_floor = opt.tangent_floor;

  FiberFrame3DResult r3 = solve_fiber_frame3d(model, case_id, o3);

  FiberFrameResult out;
  out.steps = r3.steps;
  out.iterations = r3.iterations;
  out.lambda = r3.lambda;
  out.node_disp.reserve(r3.node_disp.size());
  for (const Vec6& d : r3.node_disp) out.node_disp.emplace_back(d(0), d(2), d(4));
  out.bars.reserve(r3.bars.size());
  for (const FiberBarState3D& b3 : r3.bars) {
    FiberBarState st;
    st.bar = b3.bar;
    st.fiber = b3.fiber;
    st.q = Vec3(b3.q(0), b3.q(1), b3.q(2));
    for (const FiberStation3D& s3 : b3.stations) {
      FiberStation fs;
      fs.x = s3.x;
      fs.weight = s3.weight;
      fs.N = s3.N;
      fs.My = s3.My;
      fs.eps0 = s3.eps0;
      fs.kappa = s3.kappa_y;
      st.stations.push_back(fs);
    }
    out.bars.push_back(std::move(st));
  }
  return out;
}

}  // namespace keta

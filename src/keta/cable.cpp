#include "keta/cable.h"

#include <cmath>
#include <limits>

namespace keta {

void catenary_projections(const Catenary2D& c, double& lx, double& lz) {
  double H = c.H, V1 = c.V1, V2 = c.V2(), w = c.w, L0 = c.L0, EA = c.EA;
  lx = (H / w) * (std::asinh(V2 / H) - std::asinh(V1 / H)) + H * L0 / EA;
  lz = (c.T2() - c.T1()) / w + (L0 / EA) * (V1 + w * L0 / 2);
}

Mat3 catenary_flexibility22(const Catenary2D& c) {
  double H = c.H, V1 = c.V1, V2 = c.V2(), w = c.w, L0 = c.L0, EA = c.EA;
  double T1 = c.T1(), T2 = c.T2();
  Mat3 F = Mat3::Zero();
  F(0, 0) = (std::asinh(V2 / H) - std::asinh(V1 / H)) / w - (V2 / T2 - V1 / T1) / w + L0 / EA;
  F(0, 1) = (H / T2 - H / T1) / w;
  F(1, 0) = F(0, 1);
  F(1, 1) = (V2 / T2 - V1 / T1) / w + L0 / EA;
  return F;
}

Catenary2D solve_catenary(double dx, double dz, double L0, double EA, double w) {
  require(dx > 0 && L0 > 0 && EA > 0 && w >= 0, "catenary: invalid input");
  Catenary2D c;
  c.L0 = L0;
  c.EA = EA;
  c.w = std::max(w, 1e-9 * EA / (L0 * L0));  // strictly positive weight for the closed forms
  // Initial guess (Peyrot-Goulois style): lambda from slack.
  double chord = std::hypot(dx, dz);
  if (L0 <= chord) {
    // Taut: chord-aligned tension from the elastic stretch. BOTH components
    // must follow the chord — a steep, strongly stretched chord has |V1| of
    // the order of T0, and a weight-only V1 guess puts the Newton outside its
    // basin (the H > 0 damping then stalls it).
    double T0 = std::max({EA * (chord - L0) / L0, c.w * dx / 8, 1e-6 * EA});
    c.H = std::max(T0 * dx / chord, 1e-8 * c.w * L0);
    c.V1 = T0 * dz / chord - 0.5 * c.w * L0;
  } else {
    double sag = std::sqrt(std::max(1e-12, 3.0 * dx * (L0 - chord) / 8.0));
    c.H = std::max(c.w * dx * dx / (8 * sag), 1e-8 * c.w * L0);
    c.V1 = c.w * L0 * (dz / (2 * std::max(chord, 1e-12)) - 0.5);
  }

  double ref = std::max({dx, std::fabs(dz), L0});
  for (int it = 0; it < 200; ++it) {
    double lx, lz;
    catenary_projections(c, lx, lz);
    double rx = lx - dx, rz = lz - dz;
    // The projection integrals cancel terms of size ~T/w, so their FP noise
    // floor scales with the TENSION, not with the geometry: on a strongly
    // stretched chord (taut intermediate anchor states of a coupled solve)
    // the classic 1e-13 * geometry threshold sits BELOW that floor and the
    // Newton would spin on noise forever. Accepting the floor costs
    // delta_H ~ noise / F00 ~ noise * EA / L0, a relative error < 1e-9 here.
    double floor_ =
        8.0 * std::numeric_limits<double>::epsilon() * (c.T1() + c.T2()) / c.w;
    if (std::fabs(rx) + std::fabs(rz) < std::max(1e-13 * ref, floor_)) return c;
    Mat3 F = catenary_flexibility22(c);
    double det = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
    require(std::fabs(det) > 1e-300, "catenary: singular flexibility");
    double dH = -(F(1, 1) * rx - F(0, 1) * rz) / det;
    double dV = -(-F(1, 0) * rx + F(0, 0) * rz) / det;
    // Damped update keeping H > 0.
    double s = 1.0;
    while (c.H + s * dH <= 0) s *= 0.5;
    c.H += s * dH;
    c.V1 += s * dV;
  }
  fail("catenary: Newton did not converge");
}

double catenary_l0_for_tension(double dx, double dz, double EA, double w, double H_target) {
  require(H_target > 0, "catenary: H target must be positive");
  double chord = std::hypot(dx, dz);
  // Secant on f(L0) = H(L0) - H_target; H decreases with L0.
  double l_lo = chord * 0.5, l_hi = chord * 1.5;
  auto Hof = [&](double l0) { return solve_catenary(dx, dz, l0, EA, w).H; };
  // Expand bracket.
  for (int i = 0; i < 60 && Hof(l_lo) < H_target; ++i) l_lo *= 0.95;
  for (int i = 0; i < 60 && Hof(l_hi) > H_target; ++i) l_hi *= 1.05;
  double f_lo = Hof(l_lo) - H_target, f_hi = Hof(l_hi) - H_target;
  require(f_lo >= 0 && f_hi <= 0, "catenary: tension target outside reachable bracket");
  for (int i = 0; i < 200; ++i) {
    double mid = 0.5 * (l_lo + l_hi);
    double fm = Hof(mid) - H_target;
    if (std::fabs(fm) < 1e-10 * H_target) return mid;
    if (fm > 0)
      l_lo = mid;
    else
      l_hi = mid;
  }
  return 0.5 * (l_lo + l_hi);
}

double catenary_l0_for_end_tension(double dx, double dz, double EA, double w, double T_target) {
  return catenary_l0_for_end_tension(dx, dz, EA, w, T_target, 2);
}

double catenary_l0_for_end_tension(double dx, double dz, double EA, double w, double T_target, int end) {
  require(std::isfinite(T_target) && T_target > 0, "catenary: end tension target must be positive");
  require(end == 1 || end == 2, "catenary: the end tension is measured at end 1 or 2");
  const double chord = std::hypot(dx, dz);
  auto f = [&](double l0) {
    const Catenary2D c = solve_catenary(dx, dz, l0, EA, w);
    return (end == 1 ? c.T1() : c.T2()) - T_target;
  };
  const double lo = chord / (1.0 + 0.9 * 0.05);
  double a = lo, fa = f(lo);
  require(fa >= 0, "catenary: end tension target above the elastic-validity cap");
  double b = 0, fb = 0;
  bool bracketed = false;
  for (int k = 1; k <= 60; ++k) {
    double hi = lo * std::pow(1.05, k);
    double fhi = f(hi);
    if (fhi <= 0) {
      b = hi;
      fb = fhi;
      bracketed = true;
      break;
    }
    require(fhi <= fa, "catenary: end tension target below the tension minimum (no taut root)");
    a = hi;
    fa = fhi;
  }
  require(bracketed, "catenary: end tension target outside the reachable bracket");
  for (int i = 0; i < 200 && b - a > 1e-15 * b; ++i) {
    double mid = 0.5 * (a + b);
    double fm = f(mid);
    if (std::fabs(fm) <= 1e-12 * T_target) return mid;
    if (fm > 0)
      a = mid;
    else
      b = mid;
  }
  (void)fb;
  return 0.5 * (a + b);
}

CableEval eval_cable(const Vec3& p1, const Vec3& p2, double L0, double EA, double w) {
  Vec3 d = p2 - p1;
  double dxh = std::hypot(d.x(), d.y());
  require(dxh > 1e-9 * (std::fabs(d.z()) + 1.0), "cable: vertical chord is not supported");
  Vec3 uh(d.x() / dxh, d.y() / dxh, 0.0);
  Vec3 ez(0, 0, 1);
  Vec3 oop = uh.cross(ez);  // out-of-plane unit vector

  Catenary2D c = solve_catenary(dxh, d.z(), L0, EA, w);

  CableEval ev;
  ev.H = c.H;
  ev.T_max = std::max(c.T1(), c.T2());
  // Forces the cable exerts ON the nodes: node 1 is pulled toward the span
  // (+H along the chord) and carries V1 (negative for a hanging cable);
  // symmetric-span check: f1 = (H, -wL0/2), f2 = (-H, -wL0/2), and the
  // identity f1 + f2 = (0, 0, -w L0) holds by construction (V2 - V1 = w L0).
  ev.f1 = c.H * uh + c.V1 * ez;
  ev.f2 = -(c.H * uh + c.V2() * ez);

  Mat3 F = catenary_flexibility22(c);
  double det = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
  double Kxx = F(1, 1) / det, Kxz = -F(0, 1) / det, Kzz = F(0, 0) / det;
  ev.k3 = Kxx * (uh * uh.transpose()) + Kxz * (uh * ez.transpose() + ez * uh.transpose()) +
          Kzz * (ez * ez.transpose()) + (c.H / dxh) * (oop * oop.transpose());
  return ev;
}

std::pair<double, double> cable_ea_w(const Model& model, const Bar& b) {
  const CableProps& cp = *b.cable;
  const Material& mat = model.material(b.material);
  const SectionProps& sec = model.section(b.section).props;
  return {cp.ea > 0 ? cp.ea : mat.E * sec.A, cp.w ? *cp.w : mat.rho * sec.A * kGravity};
}

std::vector<double> resolve_cable_l0(const Model& model) {
  std::vector<double> l0(model.bars.size(), 0.0);
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    const Bar& b = model.bars[i];
    if (b.kind != BarKind::Cable) continue;
    const CableProps& cp = *b.cable;
    const auto [EA, w] = cable_ea_w(model, b);
    Vec3 d = model.node(b.nj).x - model.node(b.ni).x;
    double dxh = std::hypot(d.x(), d.y());
    if (cp.unstressed_length > 0) {
      l0[i] = cp.unstressed_length;
    } else if (cp.end_tension > 0) {
      l0[i] = catenary_l0_for_end_tension(dxh, d.z(), EA, w, cp.end_tension);
    } else {
      require(cp.pretension > 0, "cable bar " + std::to_string(b.id) +
                                     ": give unstressed_length or a positive pretension");
      l0[i] = catenary_l0_for_tension(dxh, d.z(), EA, w, cp.pretension);
    }
  }
  return l0;
}

CableStaticResults solve_cable_static(const Model& model, Id case_id,
                                      const CableStaticOptions& opt) {
  require(model.continua.empty(), "solve_cable_static: continuum elements are not supported (linear statics, modal and dynamics take them)");
  hooks_check(opt.hooks);
  model.validate();
  bool any_cable = false;
  for (const Bar& b : model.bars) any_cable |= (b.kind == BarKind::Cable);
  require(any_cable, "cable static: model has no cable bars");

  Assembled asm_ = assemble_elements(model);
  const DofMap& dofs = asm_.dofs;
  int n = dofs.n_free;
  SpMat Kf = assemble_stiffness(model, asm_, true);  // frame only (the cables are added below)
  std::vector<double> l0 = resolve_cable_l0(model);

  CaseLoads cl = lower_case(model, asm_, model.load_case(case_id));
  require(cl.prescribed.empty(), "cable static: settlements are not supported");

  struct CableRef {
    std::size_t bar_ix;
    int ni, nj;
    double EA, w;
  };
  std::vector<CableRef> cabs;
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    const Bar& b = model.bars[i];
    if (b.kind != BarKind::Cable) continue;
    const auto [EA, w] = cable_ea_w(model, b);
    cabs.push_back({i, model.node_index(b.ni), model.node_index(b.nj), EA, w});
  }

  auto cable_forces = [&](const VecX& u, std::vector<CableEval>& evs) {
    VecX q = VecX::Zero(n);
    evs.clear();
    for (const CableRef& cr : cabs) {
      const Bar& b = model.bars[cr.bar_ix];
      Vec3 p1 = model.node(b.ni).x, p2 = model.node(b.nj).x;
      for (int c = 0; c < 3; ++c) {
        int f1 = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.ni, c))];
        int f2 = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.nj, c))];
        if (f1 >= 0) p1(c) += u(f1);
        if (f2 >= 0) p2(c) += u(f2);
      }
      CableEval ev = eval_cable(p1, p2, l0[cr.bar_ix], cr.EA, cr.w);
      evs.push_back(ev);
      for (int c = 0; c < 3; ++c) {
        int f1 = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.ni, c))];
        int f2 = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.nj, c))];
        // q accumulates the forces cables exert on the nodes; the residual is
        // F_ext + q - K u = 0.
        if (f1 >= 0) q(f1) += ev.f1(c);
        if (f2 >= 0) q(f2) += ev.f2(c);
      }
    }
    return q;
  };

  VecX u = VecX::Zero(n);
  std::vector<CableEval> evs;
  int iterations = 0;
  // |Kf|, for the round-off floor of the residual (below)
  const SpMat Kabs = Kf.cwiseAbs();
  constexpr double kEps = std::numeric_limits<double>::epsilon();
  // One cache for every Newton step: the tangent keeps the pattern of Kf plus
  // the cable pairs, so each step after the first reuses the symbolic
  // analysis (fill-reducing ordering included) and refactors the numerics
  // only — the same deterministic ordering, the same factor as a fresh cache.
  FactorCache fac;
  for (int sub = 1; sub <= opt.substeps; ++sub) {
    double lam = static_cast<double>(sub) / opt.substeps;
    VecX Fext = lam * cl.rhs;
    hooks_report(opt.hooks, static_cast<double>(sub - 1) / opt.substeps,
                 "substep " + std::to_string(sub) + "/" + std::to_string(opt.substeps));
    for (int it = 0;; ++it) {
      hooks_check(opt.hooks);
      require(it < opt.max_iterations, "cable static: Newton did not converge");
      VecX q = cable_forces(u, evs);
      VecX r = Fext + q - Kf * u;
      // The scale must NOT collapse with the residual: include the cable
      // tension magnitude (an unloaded pure-cable net has Fext = Kf u = 0 and
      // q -> 0 at equilibrium, but the member forces stay finite).
      double t_scale = 0;
      for (const CableEval& ev : evs)
        t_scale = std::max({t_scale, ev.T_max, ev.f1.norm(), ev.f2.norm()});
      double scale = std::max({Fext.lpNorm<Eigen::Infinity>(), (Kf * u).lpNorm<Eigen::Infinity>(),
                               t_scale, 1e-30});
      // Converged when every component is within the relative tolerance OR
      // within the round-off floor of its own evaluation: r_i sums the terms
      // K_ij u_j, whose rounding error is ~ eps sum_j |K_ij u_j| — on a stiff
      // member (a pylon, EI / L^3 ~ 1e12 N/m, metre displacements) that floor
      // reaches 1e-3 N while the terms cancel to ~1e2 N, above 1e-10 |K u|,
      // and the Newton then spins on noise until the line search gives up on
      // an equilibrated state. The floor only relaxes such components; where
      // it is below tol_rel * scale the test is the plain relative one.
      const VecX floor_r = 8.0 * kEps * (Kabs * u.cwiseAbs());
      if (((r.cwiseAbs() - floor_r).array() <= opt.tol_rel * scale).all()) break;
      ++iterations;  // Newton steps (tangent solves)
      // Tangent: K_frame + cable pairs.
      std::vector<Triplet> trips;
      for (int kcol = 0; kcol < Kf.outerSize(); ++kcol)
        for (SpMat::InnerIterator itk(Kf, kcol); itk; ++itk)
          trips.emplace_back(static_cast<int>(itk.row()), static_cast<int>(itk.col()), itk.value());
      for (std::size_t ci = 0; ci < cabs.size(); ++ci) {
        const CableRef& cr = cabs[ci];
        const Mat3& k3 = evs[ci].k3;
        int g[6];
        for (int c = 0; c < 3; ++c) {
          g[c] = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.ni, c))];
          g[3 + c] = dofs.free_ix[static_cast<std::size_t>(dofs.full(cr.nj, c))];
        }
        for (int a = 0; a < 3; ++a)
          for (int bb = 0; bb < 3; ++bb) {
            double v = k3(a, bb);
            if (g[a] >= 0 && g[bb] >= 0) trips.emplace_back(g[a], g[bb], v);
            if (g[3 + a] >= 0 && g[3 + bb] >= 0) trips.emplace_back(g[3 + a], g[3 + bb], v);
            if (g[a] >= 0 && g[3 + bb] >= 0) trips.emplace_back(g[a], g[3 + bb], -v);
            if (g[3 + a] >= 0 && g[bb] >= 0) trips.emplace_back(g[3 + a], g[bb], -v);
          }
      }
      SpMat Kt(n, n);
      Kt.setFromTriplets(trips.begin(), trips.end());
      fac.factorize(Kt);
      VecX du = fac.solve(r);
      // Backtracking line search on ||r||_inf; invalid cable states backtrack too.
      double r0 = r.lpNorm<Eigen::Infinity>();
      double s = 1.0;
      for (int back = 0; back < 30; ++back) {
        VecX ut = u + s * du;
        try {
          std::vector<CableEval> tmp;
          VecX qt = cable_forces(ut, tmp);
          VecX rt = Fext + qt - Kf * ut;
          if (rt.lpNorm<Eigen::Infinity>() < r0) {
            u = ut;
            break;
          }
        } catch (const Error&) {
          // invalid trial state: halve
        }
        s *= 0.5;
        require(back < 29, "cable static: line search failed");
      }
    }
  }

  hooks_report(opt.hooks, 1.0, "done");
  // Recovery: frame part with the linear operator; cable columns reported apart.
  CableStaticResults out;
  StaticOptions sopt;
  sopt.stations_per_bar = opt.stations_per_bar;
  out.result = recover_case(model, asm_, cl, u, sopt, case_id, model.load_case(case_id).name);
  out.result.non_superposable = true;
  out.iterations = iterations;
  // Patch reactions with the cable anchor pull: R = K u + (-q) - F.
  VecX qfin = cable_forces(u, evs);
  for (Reaction& re : out.result.reactions) {
    int ni = model.node_index(re.node);
    for (const auto& [ci, cr] : [&] {
           std::vector<std::pair<std::size_t, CableRef>> v;
           for (std::size_t i = 0; i < cabs.size(); ++i) v.emplace_back(i, cabs[i]);
           return v;
         }()) {
      if (cr.ni == ni)
        for (int c = 0; c < 3; ++c) re.R(c) -= evs[ci].f1(c);
      if (cr.nj == ni)
        for (int c = 0; c < 3; ++c) re.R(c) -= evs[ci].f2(c);
    }
  }
  out.cable_H.assign(model.bars.size(), 0.0);
  out.cable_Tmax.assign(model.bars.size(), 0.0);
  for (std::size_t ci = 0; ci < cabs.size(); ++ci) {
    out.cable_H[cabs[ci].bar_ix] = evs[ci].H;
    out.cable_Tmax[cabs[ci].bar_ix] = evs[ci].T_max;
  }
  return out;
}

double ernst_modulus(double E, double gamma, double l, double sigma) {
  require(E > 0 && gamma >= 0 && l > 0 && sigma > 0 && finite(E) && finite(gamma) && finite(l) && finite(sigma),
          "ernst_modulus: E > 0, gamma >= 0, l > 0, sigma > 0");
  return E / (1.0 + std::pow(gamma * l, 2) * E / (12.0 * std::pow(sigma, 3)));
}

ErnstStiffness ernst_chord_stiffness(double EA, double L0, double dx, double dz, double w, double H) {
  require(EA > 0 && L0 > 0 && dx > 0 && w >= 0 && H > 0 && finite(EA) && finite(L0) && finite(dx) && finite(dz) &&
              finite(w) && finite(H),
          "ernst_chord_stiffness: EA, L0, dx, H > 0 and w >= 0");
  ErnstStiffness e;
  const double Lc = std::hypot(dx, dz), T = H * Lc / dx;
  e.k_elastic = EA / L0;
  e.chi = w * dx * dx / (Lc * H);
  if (w == 0) {
    e.k_sag = std::numeric_limits<double>::infinity();
    e.k = e.k_elastic;
  } else {
    e.k_sag = 12.0 * T * T * T / (w * w * dx * dx * Lc);
    e.k = 1.0 / (1.0 / e.k_elastic + 1.0 / e.k_sag);
  }
  e.note = e.chi > kErnstSagNote;
  return e;
}

ErnstSubstitute ernst_substitute(const Model& model, const CableStaticResults& state) {
  require(state.result.node_disp.size() == model.nodes.size() && state.cable_H.size() == model.bars.size(),
          "ernst_substitute: the cable-static state does not belong to this model");
  ErnstSubstitute out;
  out.model = model;
  out.axial = axial_field(state.result);
  out.cables.assign(model.bars.size(), ErnstStiffness{});
  const std::vector<double> l0 = resolve_cable_l0(model);
  Id next_mat = 0;
  for (const Material& m : model.materials) next_mat = std::max(next_mat, m.id);
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    const Bar& b = model.bars[i];
    if (b.kind != BarKind::Cable) continue;
    require(!b.has_ecc(), "ernst_substitute: a cable with end offsets is not supported");
    const Vec3 p1 = model.node(b.ni).x + state.result.node_disp[static_cast<std::size_t>(model.node_index(b.ni))].head<3>();
    const Vec3 p2 = model.node(b.nj).x + state.result.node_disp[static_cast<std::size_t>(model.node_index(b.nj))].head<3>();
    const Vec3 d = p2 - p1;
    const double dxh = std::hypot(d.x(), d.y());
    const auto [EA, w] = cable_ea_w(model, b);
    const double H = state.cable_H[i];
    require(H > 0, "ernst_substitute: cable bar " + std::to_string(b.id) + " is slack in the state");
    const ErnstStiffness e = ernst_chord_stiffness(EA, l0[i], dxh, d.z(), w, H);
    out.cables[i] = e;
    const double A = model.section(b.section).props.A;
    require(A > 0, "ernst_substitute: cable bar " + std::to_string(b.id) + " has no section area");
    const double Lc_model = (model.node(b.nj).x - model.node(b.ni).x).norm();
    Material m = model.material(b.material);
    m.id = ++next_mat;
    m.E = e.k * Lc_model / A;
    m.rho = w * l0[i] / kGravity / (A * Lc_model);
    out.model.materials.push_back(m);
    Bar& nb = out.model.bars[i];
    nb.kind = BarKind::Truss;
    nb.material = m.id;
    nb.cable.reset();
    out.axial[i] = H * d.norm() / dxh;
  }
  out.model.invalidate_index();
  return out;
}

}  // namespace keta

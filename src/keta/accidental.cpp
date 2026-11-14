#include "keta/accidental.h"

#include <set>
#include <string>

#include "keta/assembly.h"
#include "keta/statics.h"

namespace keta {

namespace {

std::vector<Vec6> to_nodes(const Model& m, const DofMap& d, const VecX& x) {
  std::vector<Vec6> out(m.nodes.size(), Vec6::Zero());
  for (std::size_t n = 0; n < m.nodes.size(); ++n)
    for (int c = 0; c < 6; ++c) {
      const int fi = d.free_of(static_cast<int>(n), c);
      if (fi >= 0) out[n](c) = x(fi);
    }
  return out;
}

VecX from_nodes(const DofMap& d, const std::vector<Vec6>& f) {
  VecX x = VecX::Zero(d.n_free);
  for (std::size_t n = 0; n < f.size(); ++n)
    for (int c = 0; c < 6; ++c) {
      const int fi = d.free_of(static_cast<int>(n), c);
      if (fi >= 0) x(fi) = f[n](c);
    }
  return x;
}

// The structure without the severed bars (their per-bar data and their own
// span loads go with them).
Model severed(const Model& model, const std::set<Id>& cut) {
  Model m = model;
  std::erase_if(m.bars, [&](const Bar& b) { return cut.count(b.id) != 0; });
  std::erase_if(m.bar_masses, [&](const BarMass& b) { return cut.count(b.bar) != 0; });
  for (Id b : cut) m.tapered.erase(b);
  for (LoadCase& c : m.cases)
    std::erase_if(c.loads, [&](const Load& l) { return l.bar != 0 && cut.count(l.bar) != 0; });
  m.invalidate_index();
  return m;
}

// The lowered load vectors of the timed cases on one structure.
struct Lowered {
  std::vector<VecX> rhs;
  std::vector<std::function<double(double)>> f;
  VecX at(double t, int n) const {
    VecX F = VecX::Zero(n);
    for (std::size_t k = 0; k < rhs.size(); ++k) F += (f[k] ? f[k](t) : 1.0) * rhs[k];
    return F;
  }
};

Lowered lower_all(const Model& m, const Assembled& a, const std::vector<TimedCase>& loads) {
  Lowered L;
  for (const TimedCase& tc : loads) {
    const CaseLoads cl = lower_case(m, a, m.load_case(tc.case_id));
    require(cl.prescribed.empty(), "rupture_time_history: case " + std::to_string(tc.case_id) +
                                       " prescribes support displacements (not supported)");
    L.rhs.push_back(cl.rhs);
    L.f.push_back(tc.factor);
  }
  return L;
}

std::vector<Vec6> static_state(const Model& m, const Lowered& L, double t) {
  StaticContext ctx = StaticContext::build(m, false);
  ctx.cache.require_spd();
  CaseLoads cl = empty_case_loads(ctx.asm_);
  cl.rhs = L.at(t, ctx.asm_.dofs.n_free);
  return to_nodes(m, ctx.asm_.dofs, ctx.solve_displacements(cl));
}

}  // namespace

RuptureHistoryResults rupture_time_history(const Model& model, const std::vector<TimedCase>& loads,
                                           const RuptureHistoryOptions& opt) {
  const std::string who = "rupture_time_history";
  model.validate();
  require(opt.dt > 0 && opt.steps > 0, who + ": dt > 0 and steps > 0");
  require(opt.rupture_step >= 0 && opt.rupture_step < opt.steps,
          who + ": rupture_step in 0.." + std::to_string(opt.steps - 1));
  require(!opt.bars.empty(), who + ": no severed bar (use time_integrate)");
  require(!model.has_links(), who + ": node-to-node links and constraints are not supported");
  require(model.continua.empty(), who + ": continuum elements are not supported");
  require(model.tendons.empty(), who + ": tendons are lowered by the phasing, not by this driver");
  std::set<Id> cut;
  for (Id b : opt.bars) {
    require(model.bar_index(b) >= 0, who + ": bar " + std::to_string(b) + " is not in the model");
    require(cut.insert(b).second, who + ": bar " + std::to_string(b) + " listed twice");
    require(model.bar(b).kind != BarKind::Cable,
            who + ": bar " + std::to_string(b) +
                " is a catenary cable — a nonlinear element, outside this linear integration");
    require(!model.composite.count(b), who + ": composite bar " + std::to_string(b) + " is not supported");
  }
  for (const Bar& b : model.bars)
    require(b.kind != BarKind::Cable, who + ": the model carries catenary cables (nonlinear elements)");

  const Model post = severed(model, cut);
  post.validate();
  const Assembled a1 = assemble_elements(model, true);
  const Assembled a2 = assemble_elements(post, true);
  const Lowered L1 = lower_all(model, a1, loads);
  const Lowered L2 = lower_all(post, a2, loads);
  const int n1 = a1.dofs.n_free, n2 = a2.dofs.n_free;

  RuptureHistoryResults out;
  out.dt = opt.dt;
  out.rupture_step = opt.rupture_step;
  out.static_intact = static_state(model, L1, 0.0);
  out.static_severed = static_state(post, L2, opt.steps * opt.dt);

  // segment 1: the intact structure, 0 -> rupture step
  VecX u0 = opt.start_static ? from_nodes(a1.dofs, out.static_intact) : VecX::Zero(n1);
  VecX v0 = VecX::Zero(n1);
  std::vector<Vec6> ur(model.nodes.size(), Vec6::Zero()), vr = ur;
  if (opt.rupture_step > 0) {
    const TimeHistoryResults s1 = time_integrate_from(
        model, opt.rupture_step, opt.dt, [&](int s) { return L1.at(s * opt.dt, n1); }, u0, v0, opt.damping,
        opt.scheme, opt.lumped_mass);
    for (int s = 0; s < opt.rupture_step; ++s) {
      out.u.push_back(to_nodes(model, a1.dofs, s1.u[static_cast<std::size_t>(s)]));
      out.v.push_back(to_nodes(model, a1.dofs, s1.v[static_cast<std::size_t>(s)]));
      out.a.push_back(to_nodes(model, a1.dofs, s1.a[static_cast<std::size_t>(s)]));
    }
    ur = to_nodes(model, a1.dofs, s1.u.back());
    vr = to_nodes(model, a1.dofs, s1.v.back());
  } else {
    ur = to_nodes(model, a1.dofs, u0);
  }

  // a DOF free before the rupture and fixed after it (a node left without
  // stiffness) must be at rest there: the state would be lost silently
  for (std::size_t n = 0; n < model.nodes.size(); ++n)
    for (int c = 0; c < 6; ++c)
      if (a1.dofs.free_of(static_cast<int>(n), c) >= 0 && a2.dofs.free_of(static_cast<int>(n), c) < 0 &&
          (ur[n](c) != 0 || vr[n](c) != 0))
        fail(who + ": node " + std::to_string(model.nodes[n].id) + " DOF " + std::to_string(c) +
             " loses its last stiffness at the rupture while moving (a member left hanging free)");

  // segment 2: the severed structure from the rupture state
  const int rest = opt.steps - opt.rupture_step;
  const VecX u2 = from_nodes(a2.dofs, ur), v2 = from_nodes(a2.dofs, vr);
  const TimeHistoryResults s2 = time_integrate_from(
      post, rest, opt.dt, [&](int s) { return L2.at((opt.rupture_step + s) * opt.dt, n2); }, u2, v2, opt.damping,
      opt.scheme, opt.lumped_mass);
  for (std::size_t s = 0; s < s2.u.size(); ++s) {
    out.u.push_back(to_nodes(post, a2.dofs, s2.u[s]));
    out.v.push_back(to_nodes(post, a2.dofs, s2.v[s]));
    out.a.push_back(to_nodes(post, a2.dofs, s2.a[s]));
  }
  return out;
}

CorotDynResult corot_rupture_history(const Model& model, Id case_id, const std::function<double(double)>& lam,
                                     const CorotDynOptions& opt, int rupture_step, const std::vector<Id>& bars) {
  const std::string who = "corot_rupture_history";
  model.validate();
  require(opt.steps > 0 && rupture_step >= 0 && rupture_step < opt.steps,
          who + ": rupture_step in 0.." + std::to_string(opt.steps - 1));
  require(!bars.empty(), who + ": no severed bar (use corot_newmark)");
  require(!model.has_links(), who + ": node-to-node links and constraints are not supported");
  std::set<Id> cut;
  for (Id b : bars) {
    require(model.bar_index(b) >= 0, who + ": bar " + std::to_string(b) + " is not in the model");
    require(cut.insert(b).second, who + ": bar " + std::to_string(b) + " listed twice");
  }
  const Model post = severed(model, cut);
  post.validate();

  // segment 1: the intact structure
  CorotDynResult out;
  std::vector<Vec6> ur = opt.initial_displacement, vr(model.nodes.size(), Vec6::Zero());
  if (ur.empty()) ur.assign(model.nodes.size(), Vec6::Zero());
  for (const auto& [node, vel] : opt.initial_velocity) vr[static_cast<std::size_t>(model.node_index(node))] = vel;
  if (rupture_step > 0) {
    CorotDynOptions o1 = opt;
    o1.steps = rupture_step;
    out = corot_newmark(model, case_id, lam, o1);
    ur = out.node_disp;
    vr = out.node_vel;
  }

  // a DOF free before the rupture and fixed after it must be at rest there
  const Assembled a1 = assemble_elements(model, true), a2 = assemble_elements(post, true);
  for (std::size_t n = 0; n < model.nodes.size(); ++n)
    for (int c = 0; c < 6; ++c)
      if (a1.dofs.free_of(static_cast<int>(n), c) >= 0 && a2.dofs.free_of(static_cast<int>(n), c) < 0 &&
          (ur[n](c) != 0 || vr[n](c) != 0))
        fail(who + ": node " + std::to_string(model.nodes[n].id) + " DOF " + std::to_string(c) +
             " loses its last stiffness at the rupture while moving (a member left hanging free)");

  // segment 2: the severed structure from the state reached
  CorotDynOptions o2 = opt;
  o2.steps = opt.steps - rupture_step;
  o2.initial_displacement = ur;
  for (std::size_t n = 0; n < model.nodes.size(); ++n)
    for (int c = 0; c < 6; ++c)
      if (a2.dofs.free_of(static_cast<int>(n), c) < 0) o2.initial_displacement[n](c) = 0.0;
  o2.initial_velocity.clear();
  for (std::size_t n = 0; n < model.nodes.size(); ++n) {
    Vec6 w = vr[n];
    for (int c = 0; c < 6; ++c)
      if (a2.dofs.free_of(static_cast<int>(n), c) < 0) w(c) = 0.0;
    if (w.squaredNorm() > 0) o2.initial_velocity[model.nodes[n].id] = w;
  }
  const double tr = rupture_step * opt.dt;
  const CorotDynResult s2 = corot_newmark(post, case_id, [&](double t) { return lam(t + tr); }, o2);
  if (rupture_step == 0) return s2;
  for (std::size_t i = 1; i < s2.t.size(); ++i) {
    out.t.push_back(s2.t[i] + tr);
    out.u.push_back(s2.u[i]);
  }
  out.node_disp = s2.node_disp;
  out.node_vel = s2.node_vel;
  out.iterations += s2.iterations;
  return out;
}

}  // namespace keta

#include "keta/tangent.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>

#include "keta/assembly.h"
#include "keta/statics.h"

namespace keta {

namespace {

struct Fnv {
  std::uint64_t h = 1469598103934665603ULL;
  void bytes(const void* p, std::size_t n) {
    const auto* c = static_cast<const unsigned char*>(p);
    for (std::size_t i = 0; i < n; ++i) {
      h ^= c[i];
      h *= 1099511628211ULL;
    }
  }
  void d(double v) { bytes(&v, sizeof v); }
  void i(long long v) { bytes(&v, sizeof v); }
  void v3(const Vec3& v) {
    for (int k = 0; k < 3; ++k) d(v(k));
  }
};

// The symmetric part of K, and |K - K^T|_max / |K|_max.
SpMat symmetrized(const SpMat& K, double& asymmetry) {
  const SpMat Kt = SpMat(K.transpose());
  const SpMat D = K - Kt;
  double dmax = 0, kmax = 0;
  for (int c = 0; c < D.outerSize(); ++c)
    for (SpMat::InnerIterator it(D, c); it; ++it) dmax = std::max(dmax, std::abs(it.value()));
  for (int c = 0; c < K.outerSize(); ++c)
    for (SpMat::InnerIterator it(K, c); it; ++it) kmax = std::max(kmax, std::abs(it.value()));
  asymmetry = kmax > 0 ? dmax / kmax : 0.0;
  SpMat S = (0.5 * (K + Kt)).eval();
  S.makeCompressed();
  return S;
}

void require_linear_structure(const Model& m, const std::string& who) {
  require(m.continua.empty(), who + ": continuum elements carry no geometric stiffness here");
  require(!links_have_arms(m), who + ": links with arms carry a second-order term that is not modelled");
}

// Linear elastic K + Kg(N, n_shell) of `m` on assemble_elements(m, true, &N).
SpMat linear_tangent(const Model& m, const std::vector<double>& N, const std::vector<Vec3>* sn,
                     bool cables_added_by_caller = false) {
  const Assembled a = assemble_elements(m, true, &N);
  SpMat K = assemble_stiffness(m, a, cables_added_by_caller);
  const SpMat Kg = assemble_geometric(m, a, N, sn);
  SpMat T = (K + Kg).eval();
  T.makeCompressed();
  return T;
}

void check_model(const Model& model, const TangentState& t, const std::string& who) {
  require(t.structure != nullptr, who + ": empty tangent state");
  require(model_fingerprint(model) == t.model_fingerprint,
          who + ": the tangent state was saved on another model (fingerprint mismatch)");
  require(t.converged, who + ": the tangent state is not a converged equilibrium");
  require(t.asymmetry <= 1e-8, who + ": the saved tangent is not symmetric (asymmetry " + std::to_string(t.asymmetry) +
                                   ", follower load stiffness): a symmetric eigen analysis does not apply");
}

}  // namespace

std::uint64_t model_fingerprint(const Model& model) {
  Fnv f;
  for (const Node& n : model.nodes) {
    f.i(n.id);
    f.v3(n.x);
  }
  for (const Material& m : model.materials) {
    f.i(m.id);
    f.d(m.E);
    f.d(m.nu);
    f.d(m.rho);
    f.d(m.G());
  }
  for (const Section& s : model.sections) {
    f.i(s.id);
    const SectionProps& p = s.props;
    for (double v : {p.A, p.Avy, p.Avz, p.Iy, p.Iz, p.J, p.Ip, p.Iw}) f.d(v);
  }
  for (const Bar& b : model.bars) {
    for (long long v : {static_cast<long long>(b.id), static_cast<long long>(b.ni), static_cast<long long>(b.nj),
                        static_cast<long long>(b.material), static_cast<long long>(b.section),
                        static_cast<long long>(b.kind)})
      f.i(v);
    f.d(b.roll);
    for (bool r : b.releases) f.i(r ? 1 : 0);
    f.v3(b.ecc[0]);
    f.v3(b.ecc[1]);
    if (b.cable) {
      f.d(b.cable->ea);
      f.d(b.cable->w.value_or(-1.0));
      f.d(b.cable->unstressed_length);
      f.d(b.cable->pretension);
      f.d(b.cable->end_tension);
    }
  }
  for (const Support& s : model.supports) {
    f.i(s.node);
    for (bool x : s.fixed) f.i(x ? 1 : 0);
    for (double k : s.spring) f.d(k);
  }
  f.i(static_cast<long long>(model.shells.size()));
  f.i(static_cast<long long>(model.links.size()));
  f.i(static_cast<long long>(model.rigid_links.size()));
  f.i(static_cast<long long>(model.point_masses.size()));
  return f.h;
}

TangentState tangent_of_static_case(const Model& model, Id case_id) {
  const std::string who = "tangent_of_static_case";
  model.validate();
  require_linear_structure(model, who);
  StaticContext ctx = StaticContext::build(model, true);
  StaticOptions sopt;
  sopt.stations_per_bar = 2;
  const CaseResult r = ctx.solve_case(model.load_case(case_id), sopt);
  TangentState t;
  t.source = TangentSource::StaticCase;
  t.case_id = case_id;
  t.model_fingerprint = model_fingerprint(model);
  t.structure = std::make_shared<const Model>(model);
  t.u = r.node_disp;
  t.axial = axial_field(r);
  t.shell_membrane = shell_membrane_field(r);
  t.K = linear_tangent(model, t.axial, t.shell_membrane.empty() ? nullptr : &t.shell_membrane);
  return t;
}

TangentState tangent_of_cable_state(const Model& model, Id case_id, const CableStaticResults& state) {
  const std::string who = "tangent_of_cable_state";
  model.validate();
  require_linear_structure(model, who);
  require(state.result.node_disp.size() == model.nodes.size(),
          who + ": the cable-static results do not belong to this model");
  TangentState t;
  t.source = TangentSource::CableState;
  t.case_id = case_id;
  t.model_fingerprint = model_fingerprint(model);
  t.structure = std::make_shared<const Model>(model);
  t.u = state.result.node_disp;
  t.axial = axial_field(state.result);
  for (std::size_t i = 0; i < model.bars.size(); ++i)
    if (model.bars[i].kind == BarKind::Cable) t.axial[i] = 0.0;
  // frame part: K + Kg(N_frame) (cables skipped by both assemblies)
  const Assembled a = assemble_elements(model, true, &t.axial);
  SpMat K = linear_tangent(model, t.axial, nullptr, true);
  // catenary pair tangents at the deformed anchors
  const std::vector<double> l0 = resolve_cable_l0(model);
  std::vector<Triplet> trip;
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    const Bar& b = model.bars[i];
    if (b.kind != BarKind::Cable) continue;
    require(!b.has_ecc(), who + ": a cable with end offsets is not supported");
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    const Vec3 p1 = model.node(b.ni).x + t.u[static_cast<std::size_t>(ni)].head<3>();
    const Vec3 p2 = model.node(b.nj).x + t.u[static_cast<std::size_t>(nj)].head<3>();
    const auto [EA, w] = cable_ea_w(model, b);
    const CableEval ev = eval_cable(p1, p2, l0[i], EA, w);
    const int nodes[2] = {ni, nj};
    for (int A = 0; A < 2; ++A)
      for (int B = 0; B < 2; ++B) {
        const double s = A == B ? 1.0 : -1.0;
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) {
            const int fr = a.dofs.free_of(nodes[A], r), fc = a.dofs.free_of(nodes[B], c);
            if (fr >= 0 && fc >= 0 && ev.k3(r, c) != 0) trip.emplace_back(fr, fc, s * ev.k3(r, c));
          }
      }
  }
  SpMat Kc(K.rows(), K.cols());
  Kc.setFromTriplets(trip.begin(), trip.end());
  t.K = symmetrized((K + Kc).eval(), t.asymmetry);
  return t;
}

TangentState tangent_of_corotational(const Model& model, Id case_id, const CorotResults& state,
                                     CorotRotationMeasure measure) {
  const std::string who = "tangent_of_corotational";
  model.validate();
  require(state.node_disp.size() == model.nodes.size(), who + ": the corotational results do not belong to this model");
  const CorotTangentProbe p = corot_tangent_probe(model, case_id, state.u_free, 1.0, measure, true);
  TangentState t;
  t.source = TangentSource::Corotational;
  t.case_id = case_id;
  t.model_fingerprint = model_fingerprint(model);
  t.structure = std::make_shared<const Model>(model);
  t.u = state.node_disp;
  t.axial.assign(model.bars.size(), 0.0);
  for (std::size_t i = 0; i < model.bars.size() && i < state.bar_forces.size(); ++i)
    if (model.bars[i].kind != BarKind::Cable) t.axial[i] = -state.bar_forces[i](0);
  // the probe's free DOFs are those of assemble_elements(model, true); the
  // modal pencil numbers them with the prestress auto-fix, which can only
  // keep MORE dofs free — require the two numberings to agree
  const Assembled a = assemble_elements(model, true, &t.axial);
  require(a.dofs.n_free == p.tangent.rows(),
          who + ": the corotational free DOFs differ from the prestressed numbering (a DOF braced only by the axial "
                "state)");
  const double scale = std::max(1.0, p.tangent.cwiseAbs().maxCoeff());
  t.converged = p.residual.size() == 0 || p.residual.cwiseAbs().maxCoeff() <= 1e-6 * scale;
  t.K = symmetrized(p.tangent.sparseView(), t.asymmetry);
  return t;
}

TangentState tangent_of_phase(const Model& model, const PhaseResult& phase) {
  const std::string who = "tangent_of_phase";
  require(phase.structure != nullptr,
          who + ": the phase did not record its structure (run the phasing with PhasingOptions::keep_structures)");
  for (const PhaseResult::CableState& c : phase.cables)
    require(!c.active, who + ": catenary cable " + std::to_string(c.bar) +
                           " is active in the phase — its coupled stiffness is not part of the recorded structure");
  require(phase.cumulated.bars.size() == model.bars.size(), who + ": the phase does not belong to this model");
  const Model& w = *phase.structure;
  require_linear_structure(w, who);
  TangentState t;
  t.source = TangentSource::PhasingStage;
  t.label = phase.name;
  t.model_fingerprint = model_fingerprint(model);
  t.structure = phase.structure;
  t.u.assign(w.nodes.size(), Vec6::Zero());
  for (std::size_t i = 0; i < w.nodes.size(); ++i)
    t.u[i] = phase.cumulated.node_disp[static_cast<std::size_t>(model.node_index(w.nodes[i].id))];
  t.axial.assign(w.bars.size(), 0.0);
  for (std::size_t i = 0; i < w.bars.size(); ++i)
    t.axial[i] = -phase.cumulated.bars[static_cast<std::size_t>(model.bar_index(w.bars[i].id))].end_forces_local(0);
  if (!w.shells.empty() && phase.cumulated.shells.size() == w.shells.size())
    t.shell_membrane = shell_membrane_field(phase.cumulated);
  t.K = linear_tangent(w, t.axial, t.shell_membrane.empty() ? nullptr : &t.shell_membrane);
  return t;
}

ModalResults solve_modal_on_state(const Model& model, const TangentState& tangent, ModalOptions opt) {
  check_model(model, tangent, "solve_modal_on_state");
  require(opt.axial_state == nullptr && opt.stiffness == nullptr,
          "solve_modal_on_state: the state is the tangent's (leave axial_state and stiffness unset)");
  opt.axial_state = &tangent.axial;
  opt.stiffness = &tangent.K;
  return solve_modal(*tangent.structure, opt);
}

BucklingResults solve_buckling_on_state(const Model& model, const TangentState& tangent,
                                        const std::vector<Load>& additional, int n_modes, int dense_threshold) {
  const std::string who = "solve_buckling_on_state";
  check_model(model, tangent, who);
  const Model& w = *tangent.structure;
  require_linear_structure(w, who);
  for (const Bar& b : w.bars)
    require(b.kind != BarKind::Cable, who + ": cable bars carry no linear geometric stiffness of an additional load");
  require(!additional.empty(), who + ": no additional load");
  // the axial field of the additional load, solved linearly on the elastic operator
  Model wa = w;
  wa.cases.clear();
  wa.combos.clear();
  const Id cid = wa.add_case("additional", additional);
  StaticContext ctx = StaticContext::build(wa, true);
  StaticOptions sopt;
  sopt.stations_per_bar = 2;
  const CaseResult r = ctx.solve_case(wa.load_case(cid), sopt);
  const std::vector<double> N = axial_field(r);
  const std::vector<Vec3> sn = shell_membrane_field(r);
  bool any = false;
  for (double v : N) any |= v != 0;
  for (const Vec3& v : sn) any |= v.squaredNorm() != 0;
  require(any, who + ": the additional load produces no axial or membrane force — no geometric stiffness to amplify");

  const Assembled a = assemble_elements(w, true, &tangent.axial);
  require(a.dofs.n_free == tangent.K.rows(), who + ": the saved tangent does not match its structure's DOFs");
  SpMat Kg = assemble_geometric(w, a, N, sn.empty() ? nullptr : &sn);
  SpMat KT = tangent.K;
  if (a.red.active) {
    KT = a.red.reduce(KT);
    Kg = a.red.reduce(Kg);
  }
  std::vector<PencilPair> pairs;
  if (KT.rows() <= dense_threshold) {
    require(Eigen::LLT<MatX>(MatX(KT)).info() == Eigen::Success,
            who + ": the tangent of the state is not positive definite — the structure is already unstable under "
                  "its locked-in state");
    pairs = dense_pencil(KT, Kg);
  } else {
    FactorCache fac;
    fac.factorize(KT);
    const Inertia in = fac.inertia();
    require(in.n_neg == 0 && in.n_zero == 0,
            who + ": the tangent of the state is not positive definite — the structure is already unstable under "
                  "its locked-in state");
    pairs = k_inner_lanczos(fac, KT, Kg, n_modes);
  }
  if (a.red.active)
    for (PencilPair& p : pairs) p.x = a.red.expand(p.x);
  // K_T + lambda Kg singular <=> Kg x = -(1/lambda) K_T x
  std::vector<std::pair<double, const PencilPair*>> lam;
  for (const PencilPair& p : pairs)
    if (p.mu < -1e-14 && p.residual < 1e-6) lam.emplace_back(-1.0 / p.mu, &p);
  std::sort(lam.begin(), lam.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
  BucklingResults out;
  for (const auto& [l, p] : lam) {
    if (static_cast<int>(out.load_factors.size()) >= n_modes) break;
    out.load_factors.push_back(l);
    out.modes.push_back(expand_unit_crest(a.dofs, p->x));
    out.residuals.push_back(p->residual);
  }
  return out;
}

}  // namespace keta

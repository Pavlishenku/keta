#include "keta/wind_equivalent.h"

#include <cmath>

#include "keta/nonlinear.h"

namespace keta {

Eigen::VectorXd lrc_loads(const Eigen::MatrixXd& C, const Eigen::VectorXd& I, double g) {
  require(C.rows() == C.cols() && C.rows() == I.size(), "LRC: covariance n x n and influence n");
  const double v = I.dot(C * I);
  require(v > 0, "LRC: the target response has no variance");
  return g * (C * I) / std::sqrt(v);
}

Eigen::VectorXd lrc_modal(const Eigen::MatrixXd& Cq, const Eigen::VectorXd& psi, double g) { return lrc_loads(Cq, psi, g); }

LoadCase modal_inertia_case(const Model& model, const ModalResults& modal, const Eigen::VectorXd& q, Id case_id,
                            const std::string& name) {
  require(q.size() == static_cast<Eigen::Index>(modal.modes_free.size()), "inertia loads: one amplitude per mode");
  const Assembled asm_ = assemble_elements(model);
  const SpMat M = assemble_mass(model, asm_);
  VecX u = VecX::Zero(asm_.dofs.n_free);
  for (std::size_t i = 0; i < modal.modes_free.size(); ++i)
    u += modal.omegas[i] * modal.omegas[i] * q(static_cast<Eigen::Index>(i)) * modal.modes_free[i];
  const VecX f = M * u;
  LoadCase lc;
  lc.id = case_id;
  lc.name = name;
  for (std::size_t n = 0; n < model.nodes.size(); ++n) {
    Vec3 F = Vec3::Zero(), Mo = Vec3::Zero();
    bool any = false;
    for (int c = 0; c < 6; ++c) {
      const int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(static_cast<int>(n), c))];
      if (fi < 0) continue;
      (c < 3 ? F(c) : Mo(c - 3)) = f(fi);
      any = any || f(fi) != 0;
    }
    if (any) lc.loads.push_back(Load::nodal(model.nodes[n].id, F, Mo));
  }
  return lc;
}

namespace {
// the 8 vertices of the tangent octagon in units of g, deviations (a, b)
std::vector<std::pair<double, double>> octagon(double sa, double sb, double rho) {
  const double p = std::sqrt(2 * (1 + rho)), q = std::sqrt(2 * (1 - rho));
  const std::pair<double, double> A{sa, sb * (p - 1)}, B{sa * (p - 1), sb}, C{sa * (1 - q), sb}, D{sa, sb * (1 - q)};
  auto neg = [](std::pair<double, double> v) { return std::make_pair(-v.first, -v.second); };
  return {A, B, C, neg(D), neg(A), neg(B), neg(C), D};  // around the polygon
}
}  // namespace

std::vector<Eigen::VectorXd> equiprobable_polyhedron(const Eigen::VectorXd& m, const Eigen::MatrixXd& S, double g) {
  const Eigen::Index k = m.size();
  require(k == 2 || k == 3, "equiprobable polyhedron: 2 or 3 effects");
  require(S.rows() == k && S.cols() == k && (S - S.transpose()).cwiseAbs().maxCoeff() <= 1e-12 * S.cwiseAbs().maxCoeff(),
          "equiprobable polyhedron: symmetric covariance k x k");
  require(g > 0, "equiprobable polyhedron: g > 0");
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(S);
  require(es.eigenvalues().minCoeff() >= -1e-12 * std::max(es.eigenvalues().maxCoeff(), 1e-300),
          "equiprobable polyhedron: the covariance is not positive semi-definite");
  const double lmax = std::max(es.eigenvalues().maxCoeff(), 0.0);
  std::vector<Eigen::Index> act;
  for (Eigen::Index i = 0; i < k; ++i)
    if (es.eigenvalues()(i) > 1e-12 * lmax) act.push_back(i);
  std::vector<Eigen::VectorXd> out;
  const Eigen::Index r = static_cast<Eigen::Index>(act.size());
  if (r == 0) return {m};
  if (r == 1) {
    const Eigen::VectorXd a = es.eigenvectors().col(act[0]) * std::sqrt(es.eigenvalues()(act[0]));
    return {m + g * a, m - g * a};
  }
  if (r < k) {  // rank 2 of 3: the tangent octagon of the planar ellipse, principal frame
    const Eigen::VectorXd a1 = es.eigenvectors().col(act[0]) * std::sqrt(es.eigenvalues()(act[0]));
    const Eigen::VectorXd a2 = es.eigenvectors().col(act[1]) * std::sqrt(es.eigenvalues()(act[1]));
    for (const auto& [x, y] : octagon(1, 1, 0)) out.push_back(m + g * (x * a1 + y * a2));
    return out;
  }
  // full rank: the octagon of every coordinate plane, the rest concomitant
  for (Eigen::Index a = 0; a < k; ++a)
    for (Eigen::Index b = a + 1; b < k; ++b) {
      const double sa = std::sqrt(S(a, a)), sb = std::sqrt(S(b, b)), rho = S(a, b) / (sa * sb);
      Eigen::Matrix2d Sab;
      Sab << S(a, a), S(a, b), S(a, b), S(b, b);
      const Eigen::Matrix2d inv = Sab.inverse();
      for (const auto& [x, y] : octagon(sa, sb, rho)) {
        Eigen::VectorXd v = m;
        const Eigen::Vector2d d(g * x, g * y);
        v(a) += d(0);
        v(b) += d(1);
        for (Eigen::Index c = 0; c < k; ++c)
          if (c != a && c != b) v(c) += Eigen::RowVector2d(S(c, a), S(c, b)) * inv * d;
        out.push_back(v);
      }
    }
  return out;
}

WindLoopResult wind_nonlinear_loop(const Model& model, Id permanent_case, const BuffetingSpec& spec,
                                   const BuffetingRequest& target, int n_modes, int max_iterations, double tol) {
  require(max_iterations >= 1 && tol > 0 && n_modes >= 1, "wind loop: iterations >= 1, tol > 0, modes >= 1");
  require(model.node_index(target.node) >= 0 && target.comp >= 0 && target.comp < 6, "wind loop: target node / component");
  std::vector<Load> base;
  if (permanent_case != 0) base = model.load_case(permanent_case).loads;
  const LoadCase mean = mean_buffeting_case(model, spec, 1, "mean");
  base.insert(base.end(), mean.loads.begin(), mean.loads.end());
  // a second-order static of `loads` on a copy of the model: (result, target value)
  auto so_static = [&](const std::vector<Load>& loads) {
    Model w = model;
    w.cases.clear();
    w.cases.push_back({1, "wind loop", loads});
    SecondOrderOptions so;
    so.tol_rel_du = 1e-13;
    const SecondOrderResults r = solve_second_order(w, 1, so);
    return std::make_pair(axial_field(r.result),
                          r.result.node_disp[static_cast<std::size_t>(model.node_index(target.node))](target.comp));
  };
  WindLoopResult out;
  std::vector<Load> loads = base;
  double prev_sigma = -1, prev_ext = 0;
  for (int it = 1; it <= max_iterations; ++it) {
    const auto [N, ext_prev_state] = so_static(loads);
    (void)ext_prev_state;
    ModalOptions mo;
    mo.n_modes = n_modes;
    mo.axial_state = &N;
    const ModalResults modal = solve_modal(model, mo);
    const MultiBuffetingResults b = spectral_buffeting_multi(model, modal, spec, {target});
    // the peak: modal LRC on the target, as inertia loads
    const Assembled asm_ = assemble_elements(model);
    const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(model.node_index(target.node), target.comp))];
    require(f >= 0, "wind loop: the target degree of freedom is restrained");
    Eigen::VectorXd psi(static_cast<Eigen::Index>(modal.omegas.size()));
    for (std::size_t i = 0; i < modal.omegas.size(); ++i) psi(static_cast<Eigen::Index>(i)) = modal.modes_free[i](f);
    const Eigen::VectorXd qe = lrc_modal(b.modal_cov, psi, b.stats[0].peak_factor);
    const LoadCase peak = modal_inertia_case(model, modal, qe, 1, "peak");
    loads = base;
    loads.insert(loads.end(), peak.loads.begin(), peak.loads.end());
    const auto [N2, ext] = so_static(loads);
    out.sigma = b.stats[0].sigma;
    out.peak_factor = b.stats[0].peak_factor;
    out.extreme = ext;
    out.iterations = it;
    out.axial = N2;
    out.frequencies_hz = modal.frequencies_hz;
    if (prev_sigma >= 0 && std::abs(out.sigma - prev_sigma) <= tol * out.sigma &&
        std::abs(ext - prev_ext) <= tol * std::max(std::abs(ext), 1e-300)) {
      out.converged = true;
      break;
    }
    prev_sigma = out.sigma;
    prev_ext = ext;
  }
  return out;
}

}  // namespace keta

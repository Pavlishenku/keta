#include "keta/rotating.h"

#include "keta/assembly.h"

namespace keta {

namespace {
// the bar's displacement interpolation in the local frame at xi: u = S q
// (3 x 12, local DOFs u v w rx ry rz per end) — linear axial; Hermite
// transverse for a frame (v' = rz, w' = -ry), linear for a truss
Eigen::Matrix<double, 3, 12> interp(double xi, double L, bool frame) {
  Eigen::Matrix<double, 3, 12> S = Eigen::Matrix<double, 3, 12>::Zero();
  S(0, 0) = 1 - xi;
  S(0, 6) = xi;
  if (!frame) {
    for (int c : {1, 2}) {
      S(c, c) = 1 - xi;
      S(c, 6 + c) = xi;
    }
    return S;
  }
  const double h1 = 1 - 3 * xi * xi + 2 * xi * xi * xi, h2 = L * (xi - 2 * xi * xi + xi * xi * xi);
  const double h3 = 3 * xi * xi - 2 * xi * xi * xi, h4 = L * (-xi * xi + xi * xi * xi);
  S(1, 1) = h1, S(1, 5) = h2, S(1, 7) = h3, S(1, 11) = h4;
  S(2, 2) = h1, S(2, 4) = -h2, S(2, 8) = h3, S(2, 10) = -h4;
  return S;
}
// d/dx of the transverse rows
Eigen::Matrix<double, 2, 12> slope(double xi, double L, bool frame) {
  Eigen::Matrix<double, 2, 12> G = Eigen::Matrix<double, 2, 12>::Zero();
  if (!frame) {
    for (int c : {0, 1}) {
      G(c, 1 + c) = -1 / L;
      G(c, 7 + c) = 1 / L;
    }
    return G;
  }
  const double d1 = (-6 * xi + 6 * xi * xi) / L, d2 = 1 - 4 * xi + 3 * xi * xi;
  const double d3 = (6 * xi - 6 * xi * xi) / L, d4 = -2 * xi + 3 * xi * xi;
  G(0, 1) = d1, G(0, 5) = d2, G(0, 7) = d3, G(0, 11) = d4;
  G(1, 2) = d1, G(1, 4) = -d2, G(1, 8) = d3, G(1, 10) = -d4;
  return G;
}
const double kGx[4] = {0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970263};
const double kGw[4] = {0.1739274225687269, 0.3260725774312731, 0.3260725774312731, 0.1739274225687269};
}  // namespace

RotatingModal solve_modal_rotating(const Model& model, const Spin& spin, const ModalOptions& opt) {
  require(model.continua.empty(), "solve_modal_rotating: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  require(spin.axis.norm() > 0, "rotating: the spin axis needs a direction");
  require(std::isfinite(spin.omega), "rotating: non-finite spin rate");
  require(model.shells.empty(), "rotating: shells are not taken (bars and point masses only)");
  require(opt.axial_state == nullptr, "rotating: the prestress is the centrifugal steady state");
  for (const Bar& b : model.bars) {
    require(!b.has_ecc(), "rotating: bar " + std::to_string(b.id) + " has rigid arms");
    require(model.tapered_of(b.id) == nullptr, "rotating: bar " + std::to_string(b.id) + " is tapered");
    for (bool r : b.releases) require(!r, "rotating: bar " + std::to_string(b.id) + " has end releases");
  }
  const Vec3 e = spin.axis.normalized();
  const Mat3 P = Mat3::Identity() - e * e.transpose();
  const double w2 = spin.omega * spin.omega;

  // the steady state under the exact centrifugal loads
  Model mc = model;
  mc.cases.clear();
  std::vector<Load> loads;
  for (const Bar& b : model.bars) {
    const double rA = model.material(b.material).rho * model.section(b.section).props.A;
    if (rA == 0) continue;
    const Vec3 qa = rA * w2 * P * (model.node(b.ni).x - spin.point), qb = rA * w2 * P * (model.node(b.nj).x - spin.point);
    loads.push_back(Load::linear(b.id, qa, qb));
  }
  for (const PointMass& pm : model.point_masses)
    if (pm.mass != 0) loads.push_back(Load::nodal(pm.node, pm.mass * w2 * P * (model.node(pm.node).x - spin.point)));
  mc.add_case("centrifugal", loads);
  RotatingModal out;
  StaticContext ctx = StaticContext::build(mc, false);
  StaticOptions so;
  so.stations_per_bar = 2;
  out.steady = ctx.solve_case(mc.cases.front(), so);

  // the steady axial forces hold the transverse DOFs of truss nodes (the
  // zero-stiffness auto-fix must see them, as in a prestressed modal)
  const std::vector<double> Nst = axial_field(out.steady);
  const Assembled as = assemble_elements(model, false, &Nst);
  const DofMap& d = as.dofs;
  SpMat K = assemble_stiffness(model, as);
  const SpMat M = assemble_mass(model, as, opt.lumped_mass);
  require(M.nonZeros() > 0, "rotating: model has no mass");
  std::vector<Triplet> t;
  auto scatter = [&](const Bar& b, const Mat12& g) {
    const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    int f[12];
    for (int c = 0; c < 6; ++c) {
      f[c] = d.free_of(ni, c);
      f[6 + c] = d.free_of(nj, c);
    }
    for (int r = 0; r < 12; ++r)
      for (int c = 0; c < 12; ++c)
        if (f[r] >= 0 && f[c] >= 0 && g(r, c) != 0) t.emplace_back(f[r], f[c], g(r, c));
  };
  for (std::size_t k = 0; k < model.bars.size(); ++k) {
    const Bar& b = model.bars[k];
    const BarElement& el = as.elements[k];
    const bool frame = b.kind == BarKind::Frame;
    const double L = el.L(), rA = model.material(b.material).rho * model.section(b.section).props.A;
    const Mat3& R = el.R();
    Mat12 g = Mat12::Zero();
    if (spin.stiffening) {
      // N(x) = N_i - int_0^x q_x, q_x the local axial centrifugal load (linear)
      const double Ni = -out.steady.bars[k].end_forces_local(0);
      const Vec3 xa = model.node(b.ni).x - spin.point, xb = model.node(b.nj).x - spin.point;
      const double qa = rA * w2 * (R * (P * xa))(0), qb = rA * w2 * (R * (P * xb))(0);
      for (int q = 0; q < 4; ++q) {
        const double xi = kGx[q], N = Ni - L * (qa * xi + (qb - qa) * xi * xi / 2);
        const Eigen::Matrix<double, 2, 12> G = slope(xi, L, frame);
        g += (kGw[q] * L * N) * (G.transpose() * G);
      }
    }
    if (spin.softening && rA != 0 && opt.lumped_mass) {
      // the lumped mass: half the bar's at each end, m P on its translations
      for (int end : {0, 6})
        for (int r = 0; r < 3; ++r)
          for (int c = 0; c < 3; ++c) g(end + r, end + c) -= rA * L / 2 * w2 * (R * P * R.transpose())(r, c);
    } else if (spin.softening && rA != 0) {
      const Mat3 Pl = R * P * R.transpose();
      for (int q = 0; q < 4; ++q) {
        const Eigen::Matrix<double, 3, 12> S = interp(kGx[q], L, frame);
        g -= (kGw[q] * L * rA * w2) * (S.transpose() * Pl * S);
      }
    }
    scatter(b, Mat12(el.T().transpose() * g * el.T()));
  }
  if (spin.softening)
    for (const PointMass& pm : model.point_masses) {
      const int ni = model.node_index(pm.node);
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c) {
          const int fr = d.free_of(ni, r), fc = d.free_of(ni, c);
          if (fr >= 0 && fc >= 0 && P(r, c) != 0) t.emplace_back(fr, fc, -pm.mass * w2 * P(r, c));
        }
    }
  SpMat G(d.n_free, d.n_free);
  G.setFromTriplets(t.begin(), t.end());
  K = (K + G).eval();

  std::vector<PencilPair> pairs;
  if (K.rows() <= opt.dense_threshold) {
    Eigen::LLT<MatX> llt{MatX(K)};
    require(llt.info() == Eigen::Success, "rotating: the spinning structure is beyond divergence (K_eff not positive definite)");
    pairs = dense_pencil(K, M);
  } else {
    FactorCache f;
    f.factorize(K);
    f.require_spd();
    pairs = k_inner_lanczos(f, K, M, opt.n_modes);
  }
  out.modal = modal_from_pairs(model, as, M, pairs, opt.n_modes, opt.support_coupling && !opt.lumped_mass);
  return out;
}

}  // namespace keta

#include "keta/nldyn_model.h"

#include <cmath>
#include <limits>

#include "keta/assembly.h"

namespace keta {

std::shared_ptr<Hysteresis> lrb_law(double k1, double Qd, double k2) {
  require(k1 > 0 && Qd > 0 && k2 >= 0 && k2 < k1, "lrb_law: k1 > 0, Qd > 0, 0 <= k2 < k1");
  const double dy = Qd / (k1 - k2);
  return std::make_shared<ElastoPlastic>(k1, Qd + k2 * dy, k2 / k1);
}

namespace {
// a path-independent elastic law of the deformation (keta/model.h LinkLaw)
class ElasticLawHyst : public Hysteresis {
 public:
  explicit ElasticLawHyst(LinkLaw l) : l_(std::move(l)) {}
  void trial(double d, double& f, double& kt) const override {
    f = l_.force(d);
    kt = l_.tangent(d);
  }
  void commit(double) override {}
  std::unique_ptr<Hysteresis> clone() const override { return std::make_unique<ElasticLawHyst>(*this); }

 private:
  LinkLaw l_;
};
}  // namespace

std::shared_ptr<Hysteresis> elastic_law(const LinkLaw& l) {
  require(l.pts.size() >= 2, "elastic_law: a law needs two points");
  return std::make_shared<ElasticLawHyst>(l);
}

std::shared_ptr<Hysteresis> fps_law(double W, double R, double mu, double dy) {
  require(W > 0 && R > 0 && mu > 0 && dy > 0, "fps_law: W, R, mu, dy > 0");
  const double k2 = W / R, k1 = mu * W / dy + k2;
  return std::make_shared<ElastoPlastic>(k1, k1 * dy, k2 / k1);
}

ModelNlDynResult newmark_model(const Model& model, const std::vector<LinkDevice>& devices_in, const NlExcitation& load,
                               const ModelNlDynOptions& opt) {
  model.validate();
  require(opt.dt > 0 && opt.steps > 0 && opt.beta > 0, "newmark_model: dt, steps, beta > 0");
  require(model.rigid_links.empty(), "newmark_model: rigid links are not supported (use stiff elastic links)");
  // the linear part: the devices' components out of their links (a unit
  // marker keeps their DOFs while the DOF map is built)
  // the model's own link laws (gaps, contact, compression / tension only,
  // multilinear elastic) join as conservative devices, after the caller's
  std::vector<LinkDevice> devices = devices_in;
  for (const ElasticLink& l : model.links)
    for (int c = 0; c < 6; ++c)
      if (l.law[static_cast<std::size_t>(c)]) {
        for (const LinkDevice& d : devices_in)
          require(!(d.link == l.id && d.comp == c), "newmark_model: a device on a component that carries a link law");
        LinkDevice d;
        d.link = l.id;
        d.comp = c;
        d.law = elastic_law(*l.law[static_cast<std::size_t>(c)]);
        devices.push_back(d);
      }
  Model w = model;
  for (ElasticLink& l : w.links)
    for (auto& q : l.law) q.reset();
  std::vector<std::size_t> li;
  for (const LinkDevice& d : devices) {
    std::size_t k = model.links.size();
    for (std::size_t i = 0; i < model.links.size(); ++i)
      if (model.links[i].id == d.link) k = i;
    require(k < model.links.size(), "newmark_model: no link " + std::to_string(d.link));
    require(d.comp >= 0 && d.comp < 6, "newmark_model: component 0..5");
    require(d.law || d.C > 0, "newmark_model: a device needs a law or a damper");
    require(d.C >= 0 && d.alpha > 0 && d.v_lin > 0, "newmark_model: C >= 0, alpha > 0, v_lin > 0");
    for (int q = 0; q < 6; ++q)
      require(q == d.comp || (model.links[k].k(d.comp, q) == 0 && model.links[k].k(q, d.comp) == 0),
              "newmark_model: a device component takes no coupling term in k");
    w.links[k].k(d.comp, d.comp) = 1.0;
    li.push_back(k);
  }
  Assembled asm_ = assemble_elements(w, true);
  require(!asm_.red.active, "newmark_model: rigid links are not supported");
  for (std::size_t j = 0; j < devices.size(); ++j) w.links[li[j]].k(devices[j].comp, devices[j].comp) = 0;
  for (std::size_t i = 0; i < w.links.size(); ++i) asm_.links[i] = link_element(w, w.links[i]);
  const SpMat K = assemble_stiffness(w, asm_);
  const SpMat M = assemble_mass(w, asm_, opt.lumped_mass);
  const SpMat C = (opt.rayleigh_a * M + opt.rayleigh_b * K).pruned();
  const int n = asm_.dofs.n_free;
  require(n > 0, "newmark_model: no free DOF");
  // device operators
  std::vector<std::vector<std::pair<int, double>>> B(devices.size());
  for (std::size_t j = 0; j < devices.size(); ++j) {
    const LinkElem e = link_element(model, model.links[li[j]]);
    for (int a = 0; a < 12; ++a) {
      const double coef = e.B(devices[j].comp, a);
      if (coef == 0) continue;
      const int fi = asm_.dofs.free_of(a < 6 ? e.ni : e.nj, a % 6);
      if (fi >= 0) B[j].emplace_back(fi, coef);
    }
  }
  std::vector<std::unique_ptr<Hysteresis>> law(devices.size());
  for (std::size_t j = 0; j < devices.size(); ++j)
    if (devices[j].law) law[j] = devices[j].law->clone();
  auto dot = [&](std::size_t j, const VecX& x) {
    double s = 0;
    for (const auto& [i, c] : B[j]) s += c * x(i);
    return s;
  };
  auto damper = [&](const LinkDevice& d, double v, double& f, double& ct) {
    f = 0;
    ct = 0;
    if (d.C <= 0) return;
    if (std::fabs(v) < d.v_lin) {
      ct = d.C * std::pow(d.v_lin, d.alpha - 1);
      f = ct * v;
    } else {
      f = d.C * std::copysign(std::pow(std::fabs(v), d.alpha), v);
      ct = d.alpha * d.C * std::pow(std::fabs(v), d.alpha - 1);
    }
  };
  // device forces and their stiffness / damping tangents at (u, v)
  auto devices_at = [&](const VecX& u, const VecX& v, std::vector<double>& F, std::vector<double>& kt,
                        std::vector<double>& ct) {
    F.assign(devices.size(), 0);
    kt.assign(devices.size(), 0);
    ct.assign(devices.size(), 0);
    for (std::size_t j = 0; j < devices.size(); ++j) {
      double f = 0, k = 0;
      if (law[j]) law[j]->trial(dot(j, u), f, k);
      double fd, cd;
      damper(devices[j], dot(j, v), fd, cd);
      F[j] = f + fd;
      kt[j] = k;
      ct[j] = cd;
    }
  };
  // loads
  std::vector<std::pair<VecX, std::function<double(double)>>> cases;
  for (const auto& [cid, f] : load.cases) {
    const CaseLoads cl = lower_case(w, asm_, w.load_case(cid));
    require(cl.prescribed.empty(), "newmark_model: prescribed displacements are not supported");
    cases.emplace_back(cl.rhs, f);
  }
  VecX Mr;
  if (load.ground_dir >= 0) {
    require(load.ground_dir < 3 && load.ground_accel, "newmark_model: ground direction 0..2 and its acceleration");
    VecX r = VecX::Zero(n);
    for (int ni = 0; ni < asm_.dofs.n_nodes; ++ni) {
      const int fi = asm_.dofs.free_of(ni, load.ground_dir);
      if (fi >= 0) r(fi) = 1;
    }
    Mr = M * r;
  }
  auto p_at = [&](double t) {
    VecX p = VecX::Zero(n);
    for (const auto& [v, f] : cases) p += f(t) * v;
    if (Mr.size()) p -= load.ground_accel(t) * Mr;
    return p;
  };
  // observed DOFs
  std::vector<int> obs;
  for (const auto& [node, c] : opt.observe) {
    require(c >= 0 && c < 6, "newmark_model: observed dof 0..5");
    obs.push_back(asm_.dofs.free_of(model.node_index(node), c));
  }
  ModelNlDynResult out;
  out.device_work.assign(devices.size(), 0);
  VecX u = VecX::Zero(n), v = VecX::Zero(n), a;
  std::vector<double> F, kt, ct;
  devices_at(u, v, F, kt, ct);
  VecX p = p_at(0);
  {
    VecX fd = VecX::Zero(n);
    for (std::size_t j = 0; j < devices.size(); ++j)
      for (const auto& [i, c] : B[j]) fd(i) += c * F[j];
    Eigen::SimplicialLDLT<SpMat> lm(M);
    require(lm.info() == Eigen::Success, "newmark_model: singular mass (massless free DOFs)");
    a = lm.solve(p - C * v - K * u - fd);
  }
  auto record = [&](double t) {
    out.t.push_back(t);
    std::vector<double> uo;
    for (int fi : obs) uo.push_back(fi >= 0 ? u(fi) : 0.0);
    out.u.push_back(uo);
    out.device_force.push_back(F);
    std::vector<double> dd;
    for (std::size_t j = 0; j < devices.size(); ++j) dd.push_back(dot(j, u));
    out.device_deform.push_back(dd);
  };
  record(0);
  const double dt = opt.dt, g = opt.gamma, b = opt.beta;
  double mmax = 0, kmax = 0;
  for (int k = 0; k < M.outerSize(); ++k)
    for (SpMat::InnerIterator itr(M, k); itr; ++itr) mmax = std::max(mmax, std::fabs(itr.value()));
  for (int k = 0; k < K.outerSize(); ++k)
    for (SpMat::InnerIterator itr(K, k); itr; ++itr) kmax = std::max(kmax, std::fabs(itr.value()));
  Eigen::SimplicialLDLT<SpMat> solver;
  bool analysed = false;
  for (int step = 1; step <= opt.steps; ++step) {
    const double t = step * dt;
    const VecX pn = p_at(t);
    VecX un = u + dt * v + dt * dt * (0.5 - b) * a;
    VecX an, vn;
    auto kin = [&](const VecX& x) {
      an = (x - u - dt * v - dt * dt * (0.5 - b) * a) / (b * dt * dt);
      vn = v + dt * ((1 - g) * a + g * an);
    };
    std::vector<double> Fn;
    bool done = false;
    // the residual of a trial displacement (devices' trial state, no commit)
    auto residual = [&](const VecX& x, double& scale) {
      kin(x);
      devices_at(x, vn, Fn, kt, ct);
      VecX fd = VecX::Zero(n);
      for (std::size_t j = 0; j < devices.size(); ++j)
        for (const auto& [i, c] : B[j]) fd(i) += c * Fn[j];
      const VecX Ma = M * an, Ku = K * x;
      scale = pn.lpNorm<Eigen::Infinity>() + Ma.lpNorm<Eigen::Infinity>() + (Ku + fd).lpNorm<Eigen::Infinity>() + 1e-300;
      return VecX(pn - Ma - C * vn - Ku - fd);
    };
    for (int it = 0; it < opt.max_iterations; ++it) {
      double scale;
      const VecX r = residual(un, scale);
      ++out.newton_iterations;
      // the round-off floor of the residual: a = (x - ...) / (beta dt^2)
      // carries eps |x| / (beta dt^2) of noise, times the mass
      double kd = 0;
      for (double k : kt) kd = std::max(kd, std::fabs(k));
      const double floor = 64 * std::numeric_limits<double>::epsilon() * (mmax / (b * dt * dt) + kmax + kd) *
                           std::max(un.lpNorm<Eigen::Infinity>(), u.lpNorm<Eigen::Infinity>());
      if (r.lpNorm<Eigen::Infinity>() <= std::max(opt.tol * scale, floor)) {
        done = true;
        break;
      }
      // J = M / (b dt^2) + (g / (b dt)) C + K + sum b (kt + g / (b dt) ct) b'
      std::vector<Triplet> trips;
      auto add = [&](const SpMat& A, double s) {
        for (int k = 0; k < A.outerSize(); ++k)
          for (SpMat::InnerIterator itr(A, k); itr; ++itr) trips.emplace_back(itr.row(), itr.col(), s * itr.value());
      };
      add(M, 1 / (b * dt * dt));
      add(C, g / (b * dt));
      add(K, 1);
      for (std::size_t j = 0; j < devices.size(); ++j) {
        const double kk = kt[j] + g / (b * dt) * ct[j];
        if (kk == 0) continue;
        for (const auto& [i1, c1] : B[j])
          for (const auto& [i2, c2] : B[j]) trips.emplace_back(i1, i2, kk * c1 * c2);
      }
      SpMat J(n, n);
      J.setFromTriplets(trips.begin(), trips.end());
      if (!analysed) {
        solver.analyzePattern(J);
        analysed = true;
      }
      solver.factorize(J);
      require(solver.info() == Eigen::Success, "newmark_model: singular tangent");
      // the full Newton step, halved while the residual grows (a bilinear
      // device flipping between its branches makes plain Newton cycle)
      const VecX du = solver.solve(r);
      const double r0 = r.norm();
      double s = 1;
      for (int ls = 0; ls < 30; ++ls, s *= 0.5) {
        double sc;
        if (residual(un + s * du, sc).norm() < r0) break;
      }
      un += s * du;
    }
    require(done, "newmark_model: Newton did not converge at t = " + std::to_string(t));
    kin(un);
    devices_at(un, vn, Fn, kt, ct);
    // energy terms of the step (average acceleration: exact discrete balance)
    const VecX du = un - u;
    out.input_work += 0.5 * (p + pn).dot(du);
    out.rayleigh_work += 0.5 * (v + vn).dot(C * du) ;
    for (std::size_t j = 0; j < devices.size(); ++j)
      out.device_work[j] += 0.5 * (F[j] + Fn[j]) * (dot(j, un) - dot(j, u));
    for (std::size_t j = 0; j < devices.size(); ++j)
      if (law[j]) law[j]->commit(dot(j, un));
    u = un;
    v = vn;
    a = an;
    p = pn;
    F = Fn;
    record(t);
  }
  out.kinetic = 0.5 * v.dot(M * v);
  out.strain = 0.5 * u.dot(K * u);
  return out;
}

}  // namespace keta

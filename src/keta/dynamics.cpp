#include "keta/dynamics.h"

#include <complex>

namespace keta {

// ---------------- damping ----------------

double DampingSpec::zeta_for(double omega, std::size_t mode_index) const {
  switch (kind) {
    case Kind::Modal: return zeta;
    case Kind::PerMode:
      require(mode_index < per_mode.size(), "damping: per-mode list shorter than mode count");
      return per_mode[mode_index];
    case Kind::Rayleigh: return alpha / (2 * omega) + beta * omega / 2;
  }
  return zeta;
}

void oscillator_step(double w, double zeta, double dt, double p0, double p1, double& u, double& v) {
  require(dt > 0 && zeta >= 0 && zeta < 1, "oscillator: dt > 0, 0 <= zeta < 1");
  if (w <= 0) {  // a free mass: u'' = p, p linear
    u += v * dt + dt * dt * (2 * p0 + p1) / 6;
    v += dt * (p0 + p1) / 2;
    return;
  }
  // particular solution of p0 + s t, then the homogeneous part (exact)
  const double s = (p1 - p0) / dt, w2 = w * w, wd = w * std::sqrt(1 - zeta * zeta);
  const double C1 = u - (p0 / w2 - 2 * zeta * s / (w2 * w));
  const double C2 = (v - s / w2 + zeta * w * C1) / wd;
  const double e = std::exp(-zeta * w * dt), c = std::cos(wd * dt), sn = std::sin(wd * dt);
  u = (p1 / w2 - 2 * zeta * s / (w2 * w)) + e * (C1 * c + C2 * sn);
  v = s / w2 + e * ((-zeta * w * C1 + wd * C2) * c + (-zeta * w * C2 - wd * C1) * sn);
}

DampingSpec DampingSpec::rayleigh_two_mode(double om1, double z1, double om2, double z2) {
  DampingSpec d;
  d.kind = Kind::Rayleigh;
  double den = om2 * om2 - om1 * om1;
  require(std::fabs(den) > 0, "rayleigh calibration: equal frequencies");
  d.alpha = 2 * om1 * om2 * (z1 * om2 - z2 * om1) / den;
  d.beta = 2 * (z2 * om2 - z1 * om1) / den;
  return d;
}

IntegrationScheme IntegrationScheme::hht(double alpha) {
  require(alpha >= -1.0 / 3.0 && alpha <= 0.0, "HHT alpha must be in [-1/3, 0]");
  IntegrationScheme s;
  s.hht_alpha = alpha;
  s.gamma = (1 - 2 * alpha) / 2;
  s.beta = (1 - alpha) * (1 - alpha) / 4;
  return s;
}

// ---------------- direct integration ----------------

namespace {
TimeHistoryResults integrate_newmark(const Model& model, int n_steps, double dt,
                                     const std::function<VecX(int)>& force, const VecX* u0, const VecX* v0,
                                     const DampingSpec& damping, const IntegrationScheme& scheme, bool lumped_mass);
}  // namespace

TimeHistoryResults time_integrate(const Model& model, int n_steps, double dt,
                                  const std::function<VecX(int)>& force,
                                  const DampingSpec& damping, const IntegrationScheme& scheme,
                                  bool lumped_mass) {
  return integrate_newmark(model, n_steps, dt, force, nullptr, nullptr, damping, scheme, lumped_mass);
}

TimeHistoryResults time_integrate_from(const Model& model, int n_steps, double dt,
                                       const std::function<VecX(int)>& force, const VecX& u0, const VecX& v0,
                                       const DampingSpec& damping, const IntegrationScheme& scheme,
                                       bool lumped_mass) {
  return integrate_newmark(model, n_steps, dt, force, &u0, &v0, damping, scheme, lumped_mass);
}

namespace {
TimeHistoryResults integrate_newmark(const Model& model, int n_steps, double dt,
                                     const std::function<VecX(int)>& force, const VecX* u0, const VecX* v0,
                                     const DampingSpec& damping, const IntegrationScheme& scheme, bool lumped_mass) {
  require(dt > 0 && n_steps > 0, "time_integrate: bad dt or step count");
  model.validate();
  Assembled asm_ = assemble_elements(model, true);
  require(asm_.dofs.n_free > 0, "time_integrate: no free DOFs");
  SpMat K = assemble_stiffness(model, asm_);
  SpMat M = assemble_mass(model, asm_, lumped_mass);
  // rigid links: the integration runs on the reduced coordinates u = T u_r
  // (exact elimination), the forces reduced by T^t, the states expanded
  const bool red = asm_.red.active;
  const int n_ext = asm_.dofs.n_free;
  if (red) {
    K = asm_.red.reduce(K);
    M = asm_.red.reduce(M);
  }
  const int n = static_cast<int>(K.rows());
  auto force_r = [&](int s) {
    const VecX f = force(s);
    require(f.size() == n_ext, "time_integrate: force vector size mismatch");
    return red ? VecX(asm_.red.reduce(f)) : f;
  };
  require(M.nonZeros() > 0, "time_integrate: model has no mass");

  require(damping.kind == DampingSpec::Kind::Rayleigh || damping.zeta == 0,
          "time_integrate: direct integration supports Rayleigh damping only "
          "(use time_history_modal for modal damping)");
  double ca = damping.kind == DampingSpec::Kind::Rayleigh ? damping.alpha : 0.0;
  double cb = damping.kind == DampingSpec::Kind::Rayleigh ? damping.beta : 0.0;
  SpMat C = (ca * M + cb * K).eval();

  double g = scheme.gamma, b = scheme.beta, ah = scheme.hht_alpha;
  double op = 1 + ah;  // (1 + alpha)

  SpMat Aeff = (M + (op * g * dt) * C + (op * b * dt * dt) * K).eval();
  FactorCache fac;
  fac.factorize(Aeff);

  TimeHistoryResults out;
  out.dt = dt;
  VecX u = VecX::Zero(n), v = VecX::Zero(n);
  if (u0 || v0) {
    require(!red, "time_integrate_from: rigid links are not supported with initial conditions");
    require((!u0 || u0->size() == n) && (!v0 || v0->size() == n),
            "time_integrate_from: initial conditions of " + std::to_string(n) + " free DOFs expected");
    if (u0) u = *u0;
    if (v0) v = *v0;
  }
  VecX F0 = force_r(0);
  // initial acceleration M a0 = F0 - C v0 - K u0; M is singular when some
  // DOFs carry no mass (rotations of a diagonal translational mass, ST1's
  // default): those take the effective matrix's (quasi-static) answer —
  // exact a0 = 0 from rest under a zero initial force
  VecX a;
  try {
    FactorCache mfac;
    mfac.factorize(M);
    a = mfac.solve(VecX(F0 - C * v - K * u));
  } catch (const Error&) {
    a = fac.solve(VecX(F0 - C * v - K * u));
  }

  auto push = [&](const VecX& uu, const VecX& vv, const VecX& aa) {
    out.u.push_back(red ? asm_.red.expand(uu) : uu);
    out.v.push_back(red ? asm_.red.expand(vv) : vv);
    out.a.push_back(red ? asm_.red.expand(aa) : aa);
    out.energy.push_back(0.5 * vv.dot(M * vv) + 0.5 * uu.dot(K * uu));
  };
  push(u, v, a);

  VecX Fn = F0;
  for (int s = 1; s <= n_steps; ++s) {
    VecX Fn1 = force_r(s);
    VecX ut = u + dt * v + dt * dt * (0.5 - b) * a;
    VecX vt = v + dt * (1 - g) * a;
    // HHT: M a1 + (1+ah)(C v1 + K u1) - ah (C v0 + K u0) = (1+ah)F1 - ah F0
    VecX rhs = op * Fn1 - ah * Fn - op * (C * vt + K * ut) + ah * (C * v + K * u);
    VecX a1 = fac.solve(rhs);
    u = ut + b * dt * dt * a1;
    v = vt + g * dt * a1;
    a = a1;
    Fn = Fn1;
    push(u, v, a);
  }
  return out;
}
}  // namespace

TimeHistoryResults time_history_modal(const Model& model, const ModalResults& modal, int n_steps,
                                      double dt, const std::function<VecX(int)>& force,
                                      const DampingSpec& damping,
                                      const IntegrationScheme& scheme) {
  (void)model;
  require(dt > 0 && n_steps > 0, "time_history_modal: bad dt or step count");
  std::size_t nm = modal.omegas.size();
  require(nm > 0, "time_history_modal: no modes");
  int n = static_cast<int>(modal.modes_free.front().size());

  double g = scheme.gamma, b = scheme.beta, ah = scheme.hht_alpha;
  double op = 1 + ah;

  std::vector<double> q(nm, 0), qd(nm, 0), qdd(nm, 0);
  VecX F = force(0);
  std::vector<double> pn(nm);
  for (std::size_t i = 0; i < nm; ++i) {
    pn[i] = modal.modes_free[i].dot(F);
    qdd[i] = pn[i];  // from rest: q'' = p(t0), unit modal mass
  }

  TimeHistoryResults out;
  out.dt = dt;
  auto recombine = [&]() {
    VecX u = VecX::Zero(n), v = VecX::Zero(n), a = VecX::Zero(n);
    for (std::size_t i = 0; i < nm; ++i) {
      u += q[i] * modal.modes_free[i];
      v += qd[i] * modal.modes_free[i];
      a += qdd[i] * modal.modes_free[i];
    }
    out.u.push_back(u);
    out.v.push_back(v);
    out.a.push_back(a);
    out.energy.push_back(0.0);
  };
  recombine();

  for (int s = 1; s <= n_steps; ++s) {
    VecX F1 = force(s);
    for (std::size_t i = 0; i < nm; ++i) {
      double om = modal.omegas[i];
      double z = damping.zeta_for(om, i);
      double c = 2 * z * om, k = om * om;
      double p1 = modal.modes_free[i].dot(F1);
      double aeff = 1 + op * g * dt * c + op * b * dt * dt * k;
      double ut = q[i] + dt * qd[i] + dt * dt * (0.5 - b) * qdd[i];
      double vt = qd[i] + dt * (1 - g) * qdd[i];
      double rhs = op * p1 - ah * pn[i] - op * (c * vt + k * ut) + ah * (c * qd[i] + k * q[i]);
      double a1 = rhs / aeff;
      q[i] = ut + b * dt * dt * a1;
      qd[i] = vt + g * dt * a1;
      qdd[i] = a1;
      pn[i] = p1;
    }
    recombine();
  }
  return out;
}

std::function<VecX(int)> seismic_force(const Model& model, int dir,
                                       const std::vector<double>& accel, bool lumped_mass,
                                       bool support_coupling) {
  Assembled asm_ = assemble_elements(model, true);
  SpMat M = assemble_mass(model, asm_, lumped_mass);
  VecX r = VecX::Zero(asm_.dofs.n_free);
  for (int ni = 0; ni < asm_.dofs.n_nodes; ++ni) {
    const Vec3& p = model.nodes[static_cast<std::size_t>(ni)].x;
    Vec6 rn = Vec6::Zero();
    if (dir < 3) {
      rn(dir) = 1;
    } else {
      Vec3 e = Vec3::Zero();
      e(dir - 3) = 1;
      rn.segment<3>(0) = e.cross(p);
      rn(dir) = 1;
    }
    for (int c = 0; c < 6; ++c) {
      int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
      if (fi >= 0) r(fi) += rn(c);
    }
  }
  VecX Mr = M * r;
  // the consistent mass coupling to the moving supports (eigsolve.h)
  if (support_coupling) {
    const VecX cpl = ground_coupling(model, asm_, lumped_mass, dir);
    if (cpl.squaredNorm() > 0) Mr += cpl;
  }
  return [Mr, accel](int step) -> VecX {
    double ag = step < static_cast<int>(accel.size()) ? accel[static_cast<std::size_t>(step)] : 0.0;
    return -ag * Mr;
  };
}

std::function<VecX(int)> harmonic_force(const Model& model, Id node, int comp, double amplitude,
                                        double omega, double dt) {
  Assembled asm_ = assemble_elements(model, true);
  int ni = model.node_index(node);
  int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, comp))];
  require(fi >= 0, "harmonic_force: target DOF is blocked");
  int n = asm_.dofs.n_free;
  return [=](int step) -> VecX {
    VecX f = VecX::Zero(n);
    f(fi) = amplitude * std::sin(omega * step * dt);
    return f;
  };
}

// ---------------- harmonic ----------------

HarmonicResults solve_harmonic(const Model& model, Id case_id, const std::vector<double>& freqs_hz,
                               const DampingSpec& damping, int n_modes, int stations,
                               const SolveHooks* hooks) {
  hooks_check(hooks);
  hooks_report(hooks, 0.0, "modal analysis");
  model.validate();
  Assembled asm_ = assemble_elements(model, true);  // links: the modes are solve_modal's (reduced pencil)
  ModalOptions mo;
  mo.n_modes = n_modes > 0 ? n_modes : (asm_.red.active ? asm_.red.n_red : asm_.dofs.n_free);
  ModalResults modal = solve_modal(model, mo);
  require(!modal.omegas.empty(), "harmonic: no modes");
  for (std::size_t i = 0; i < modal.omegas.size(); ++i)
    require(damping.zeta_for(modal.omegas[i], i) > 0,
            "harmonic: an undamped resonance would blow up the FRF - give every mode zeta > 0");

  CaseLoads cl = lower_case(model, asm_, model.load_case(case_id));
  require(cl.prescribed.empty(), "harmonic: settlements are not supported in the FRF driver");
  const VecX& P = cl.rhs;
  // Load-free lowering for the imaginary (quadrature) part: the recovery
  // operator is affine and the applied span loads are IN PHASE with the drive,
  // so they belong to the real part only.
  CaseLoads cl0;
  StaticOptions sopt;
  if (stations > 0) {
    cl0 = lower_case(model, asm_, LoadCase{});
    sopt.stations_per_bar = stations;
  }

  HarmonicResults out;
  out.freqs_hz = freqs_hz;
  out.omegas = modal.omegas;
  for (std::size_t i = 0; i < modal.omegas.size(); ++i) out.zetas.push_back(damping.zeta_for(modal.omegas[i], i));
  int nn = asm_.dofs.n_nodes;
  const std::size_t nf = freqs_hz.size();
  for (std::size_t k = 0; k < nf; ++k) {
    const double f = freqs_hz[k];
    hooks_check(hooks);
    hooks_report(hooks, static_cast<double>(k) / static_cast<double>(nf),
                 "frequency " + std::to_string(k + 1) + "/" + std::to_string(nf));
    double om = 2 * kPi * f;
    VecX ur = VecX::Zero(asm_.dofs.n_free), ui = VecX::Zero(asm_.dofs.n_free);
    for (std::size_t i = 0; i < modal.omegas.size(); ++i) {
      double omi = modal.omegas[i];
      double z = damping.zeta_for(omi, i);
      double p = modal.modes_free[i].dot(P);
      std::complex<double> H = p / std::complex<double>(omi * omi - om * om, 2 * z * omi * om);
      ur += H.real() * modal.modes_free[i];
      ui += H.imag() * modal.modes_free[i];
    }
    std::vector<Vec6> dr(static_cast<std::size_t>(nn), Vec6::Zero()),
        di(static_cast<std::size_t>(nn), Vec6::Zero()),
        da(static_cast<std::size_t>(nn), Vec6::Zero());
    for (int ni = 0; ni < nn; ++ni)
      for (int c = 0; c < 6; ++c) {
        int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        if (fi >= 0) {
          dr[static_cast<std::size_t>(ni)](c) = ur(fi);
          di[static_cast<std::size_t>(ni)](c) = ui(fi);
          da[static_cast<std::size_t>(ni)](c) = std::hypot(ur(fi), ui(fi));
        }
      }
    out.disp_re.push_back(std::move(dr));
    out.disp_im.push_back(std::move(di));
    out.disp_amp.push_back(std::move(da));

    if (stations > 0) {
      // Complex station diagrams: the STATIC recovery applied to each part
      // (span loads on the real part only; see the header for the distributed
      // inertia caveat).
      CaseResult cre = recover_case(model, asm_, cl, ur, sopt, case_id, "re");
      CaseResult cim = recover_case(model, asm_, cl0, ui, sopt, case_id, "im");
      std::vector<BarResult> amp = cre.bars;
      for (std::size_t e = 0; e < amp.size(); ++e)
        for (std::size_t k = 0; k < amp[e].stations.size(); ++k) {
          Station& a = amp[e].stations[k];
          const Station& sr = cre.bars[e].stations[k];
          const Station& si = cim.bars[e].stations[k];
          a.N = std::hypot(sr.N, si.N);
          a.Vy = std::hypot(sr.Vy, si.Vy);
          a.Vz = std::hypot(sr.Vz, si.Vz);
          a.T = std::hypot(sr.T, si.T);
          a.My = std::hypot(sr.My, si.My);
          a.Mz = std::hypot(sr.Mz, si.Mz);
          for (int c = 0; c < 3; ++c) a.defl(c) = std::hypot(sr.defl(c), si.defl(c));
          for (int c = 0; c < 3; ++c) a.rot(c) = std::hypot(sr.rot(c), si.rot(c));
        }
      out.bars_re.push_back(std::move(cre.bars));
      out.bars_im.push_back(std::move(cim.bars));
      out.bars_amp.push_back(std::move(amp));
    }
  }
  hooks_report(hooks, 1.0, "done");
  return out;
}

// ---------------- PSD ----------------

PsdResults solve_psd(const Model& model, Id case_id, const std::vector<double>& freqs_hz,
                     const std::vector<double>& s_alpha, const DampingSpec& damping, int n_modes,
                     int stations, const SolveHooks* hooks) {
  require(freqs_hz.size() == s_alpha.size() && freqs_hz.size() >= 2, "psd: bad input grids");
  HarmonicResults h = solve_harmonic(model, case_id, freqs_hz, damping, n_modes, stations, hooks);
  hooks_check(hooks);
  PsdResults out;
  out.freqs_hz = freqs_hz;
  std::size_t nn = h.disp_amp.front().size();
  out.disp_psd.assign(freqs_hz.size(), std::vector<Vec6>(nn, Vec6::Zero()));
  for (std::size_t k = 0; k < freqs_hz.size(); ++k)
    for (std::size_t ni = 0; ni < nn; ++ni)
      for (int c = 0; c < 6; ++c) {
        double amp = h.disp_amp[k][ni](c);
        out.disp_psd[k][ni](c) = amp * amp * s_alpha[k];
      }
  out.rms.assign(nn, Vec6::Zero());
  for (std::size_t ni = 0; ni < nn; ++ni)
    for (int c = 0; c < 6; ++c) {
      double var = 0;
      for (std::size_t k = 1; k < freqs_hz.size(); ++k)
        var += 0.5 * (out.disp_psd[k - 1][ni](c) + out.disp_psd[k][ni](c)) *
               (freqs_hz[k] - freqs_hz[k - 1]);
      out.rms[ni](c) = std::sqrt(var);
    }

  if (stations > 0 && !h.bars_amp.empty()) {
    // Station-force RMS by the same linearity: S_q(f) = |H_q(f)|^2 S_alpha,
    // trapezoid over the SAME frequency grid as the nodal RMS.
    out.rms_bars = h.bars_amp.front();  // shape template (abscissas, ids)
    std::size_t nb = out.rms_bars.size();
    for (std::size_t e = 0; e < nb; ++e)
      for (std::size_t s = 0; s < out.rms_bars[e].stations.size(); ++s) {
        auto var_of = [&](auto get) {
          double var = 0;
          for (std::size_t k = 1; k < freqs_hz.size(); ++k) {
            double a0 = get(h.bars_amp[k - 1][e].stations[s]);
            double a1 = get(h.bars_amp[k][e].stations[s]);
            var += 0.5 * (a0 * a0 * s_alpha[k - 1] + a1 * a1 * s_alpha[k]) *
                   (freqs_hz[k] - freqs_hz[k - 1]);
          }
          return std::sqrt(var);
        };
        Station& st = out.rms_bars[e].stations[s];
        st.N = var_of([](const Station& q) { return q.N; });
        st.Vy = var_of([](const Station& q) { return q.Vy; });
        st.Vz = var_of([](const Station& q) { return q.Vz; });
        st.T = var_of([](const Station& q) { return q.T; });
        st.My = var_of([](const Station& q) { return q.My; });
        st.Mz = var_of([](const Station& q) { return q.Mz; });
        for (int c = 0; c < 3; ++c)
          st.defl(c) = var_of([c](const Station& q) { return q.defl(c); });
        for (int c = 0; c < 3; ++c)
          st.rot(c) = var_of([c](const Station& q) { return q.rot(c); });
      }
  }
  return out;
}

// ---------------- response spectrum ----------------

double ResponseSpectrum::sa(double T, double zeta) const {
  switch (kind) {
    case Kind::Table: {
      require(!table.empty(), "spectrum: empty table");
      // damping correction of a table given at table_zeta: the EN 1998-1
      // eta(zeta) / eta(table_zeta), eta = max(0.55, sqrt(10 / (5 + 100 zeta)))
      double f = 1;
      if (table_zeta > 0) {
        auto eta = [](double z) { return std::max(0.55, std::sqrt(10.0 / (5.0 + 100 * z))); };
        f = T >= table_correction_from ? eta(zeta) / eta(table_zeta) : 1.0;
      }
      if (T <= table.front().first) return f * table.front().second;
      if (T >= table.back().first) return f * table.back().second;
      for (std::size_t i = 1; i < table.size(); ++i)
        if (T <= table[i].first) {
          double t0 = table[i - 1].first, t1 = table[i].first;
          double s0 = table[i - 1].second, s1 = table[i].second;
          return f * (s0 + (s1 - s0) * (T - t0) / (t1 - t0));
        }
      return f * table.back().second;
    }
    case Kind::Ec8Elastic: {
      double eta = eta_fixed > 0 ? eta_fixed : std::max(0.55, std::sqrt(10.0 / (5.0 + 100 * zeta)));
      if (T <= TB) return ag * S * (1 + T / TB * (eta * plateau - 1));
      if (T <= TC) return ag * S * eta * plateau;
      if (T <= TD) return ag * S * eta * plateau * TC / T;
      return ag * S * eta * plateau * TC * TD / (T * T);
    }
    case Kind::Ec8Design: {
      double v;
      if (T <= TB)
        v = ag * S * (2.0 / 3.0 + T / TB * (plateau / q - 2.0 / 3.0));
      else if (T <= TC)
        v = ag * S * plateau / q;
      else if (T <= TD)
        v = std::max(ag * S * plateau / q * TC / T, beta_floor * ag);
      else
        v = std::max(ag * S * plateau / q * TC * TD / (T * T), beta_floor * ag);
      return v;
    }
  }
  return 0;
}

ResponseSpectrum ResponseSpectrum::ec8_elastic(double ag, double S, double TB, double TC,
                                               double TD) {
  ResponseSpectrum s;
  s.kind = Kind::Ec8Elastic;
  s.ag = ag;
  s.S = S;
  s.TB = TB;
  s.TC = TC;
  s.TD = TD;
  return s;
}
ResponseSpectrum ResponseSpectrum::ec8_design(double ag, double S, double TB, double TC, double TD,
                                              double q) {
  ResponseSpectrum s = ec8_elastic(ag, S, TB, TC, TD);
  s.kind = Kind::Ec8Design;
  s.q = q;
  return s;
}

double cqc_rho(double om_i, double om_j, double zi, double zj) {
  // Rosenblueth-Elishakoff two-damping form (reduces to Der Kiureghian at zi == zj).
  double r = om_j / om_i;
  double num = 8 * std::sqrt(zi * zj) * (zi + r * zj) * std::pow(r, 1.5);
  double den = std::pow(1 - r * r, 2) + 4 * zi * zj * r * (1 + r * r) +
               4 * (zi * zi + zj * zj) * r * r;
  return num / den;
}

double double_sum_rho(double om_i, double om_j, double zi, double zj, double td) {
  double omdi = om_i * std::sqrt(std::max(0.0, 1 - zi * zi));
  double omdj = om_j * std::sqrt(std::max(0.0, 1 - zj * zj));
  double zpi = zi + 2.0 / (td * om_i);
  double zpj = zj + 2.0 / (td * om_j);
  double eps = (omdi - omdj) / (zpi * om_i + zpj * om_j);
  return 1.0 / (1.0 + eps * eps);
}

SpectrumResults solve_spectrum(const Model& model, const ResponseSpectrum& spectrum,
                               const DampingSpec& damping, const SpectrumOptions& opt) {
  ModalOptions mo;
  mo.n_modes = opt.n_modes;
  mo.support_coupling = opt.support_coupling;
  mo.lumped_mass = opt.lumped_mass;
  ModalResults modal = solve_modal(model, mo);
  std::size_t nm = modal.omegas.size();
  require(nm > 0, "spectrum: no modes");

  // node-to-node links: the modes come from solve_modal (reduced pencil);
  // the per-mode recovery below only reads displacement fields
  Assembled asm_ = assemble_elements(model, true);
  SpectrumResults out;
  out.frequencies_hz = modal.frequencies_hz;
  std::size_t nn = modal.modes.front().size();

  // Per-mode station recovery setup (load-free static recovery of the signed
  // modal peak field Gamma_i Sd_i phi_i).
  CaseLoads cl0;
  StaticOptions sopt;
  if (opt.stations > 0) {
    cl0 = lower_case(model, asm_, LoadCase{});
    sopt.stations_per_bar = opt.stations;
    sopt.stations_abs = opt.stations_abs;
  }

  // excitation along a vector of translations: Gamma_v = v . Gamma (the
  // participation is linear in the influence vector)
  const bool use_vec = opt.dir_vector.squaredNorm() > 0;
  Vec3 v = Vec3::Zero();
  if (use_vec) {
    require(opt.dir_vector.tail<3>().squaredNorm() == 0, "spectrum: the direction vector holds translations only");
    v = opt.dir_vector.head<3>().normalized();
  } else {
    require(opt.dir >= 0 && opt.dir < 6, "spectrum: direction 0..5");
  }
  // DEGENERATE modes (equal frequencies to 1e-9: the paired modes of a
  // symmetric structure) have an arbitrary basis — set by round-off — and
  // an SRSS of per-mode quantities is not invariant under its rotation:
  // the cluster is rotated (orthogonally: the M-orthonormality and the
  // eigenspace are kept) so that its first mode carries ALL the
  // participation in the excitation direction and the others none — the
  // result is then independent of round-off, and equal to CQC's for the
  // cluster (rho = 1 at equal frequencies)
  {
    auto gamma_of = [&](std::size_t i) {
      if (!use_vec) return modal.participation[i][static_cast<std::size_t>(opt.dir)];
      double g = 0;
      for (int c = 0; c < 3; ++c) g += v(c) * modal.participation[i][static_cast<std::size_t>(c)];
      return g;
    };
    for (std::size_t a = 0; a < nm;) {
      std::size_t b = a + 1;
      while (b < nm && std::fabs(modal.omegas[b] - modal.omegas[a]) <= 1e-9 * std::max(modal.omegas[a], 1e-300)) ++b;
      const std::size_t c = b - a;
      if (c >= 2) {
        VecX g(static_cast<Eigen::Index>(c));
        for (std::size_t k = 0; k < c; ++k) g(static_cast<Eigen::Index>(k)) = gamma_of(a + k);
        const double gn = g.norm();
        if (gn > 0) {
          // an orthonormal Q with first column g / |g| (Householder)
          VecX w = g / gn;
          w(0) += w(0) >= 0 ? 1.0 : -1.0;
          MatX Q = MatX::Identity(static_cast<Eigen::Index>(c), static_cast<Eigen::Index>(c)) - 2 * w * w.transpose() / w.squaredNorm();
          if (Q.col(0).dot(g) < 0) Q.col(0) *= -1;  // first column = +g / |g|
          std::vector<VecX> phi(c);
          std::vector<std::array<double, 6>> part(c);
          for (std::size_t k = 0; k < c; ++k) {
            phi[k] = VecX::Zero(modal.modes_free[a].size());
            part[k].fill(0);
          }
          for (std::size_t l = 0; l < c; ++l)
            for (std::size_t k = 0; k < c; ++k) {
              const double q = Q(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(l));
              phi[l] += q * modal.modes_free[a + k];
              for (int d = 0; d < 6; ++d) part[l][static_cast<std::size_t>(d)] += q * modal.participation[a + k][static_cast<std::size_t>(d)];
            }
          for (std::size_t l = 0; l < c; ++l) {
            modal.modes_free[a + l] = phi[l];
            modal.participation[a + l] = part[l];
            for (int d = 0; d < 6; ++d)
              modal.effective_mass[a + l][static_cast<std::size_t>(d)] = part[l][static_cast<std::size_t>(d)] * part[l][static_cast<std::size_t>(d)];
          }
        }
      }
      a = b;
    }
  }
  std::vector<double> vals(nm);  // Gamma * Sd per mode
  for (std::size_t i = 0; i < nm; ++i) {
    double om = modal.omegas[i];
    double z = damping.zeta_for(om, i);
    out.zetas.push_back(z);
    double gamma = 0;
    if (use_vec)
      for (int c = 0; c < 3; ++c) gamma += v(c) * modal.participation[i][static_cast<std::size_t>(c)];
    else
      gamma = modal.participation[i][static_cast<std::size_t>(opt.dir)];
    out.participation.push_back(gamma);
    double sd = 0, sa = 0;
    if (om > 1e-9) {  // rigid-mode gating on displacement AND base shear
      double T = 2 * kPi / om;
      sa = spectrum.sa(T, z);
      sd = sa / (om * om);
    }
    vals[i] = gamma * sd;
    out.modal_base_shear.push_back(gamma * gamma * sa);
    std::vector<Vec6> field(nn, Vec6::Zero());
    for (int ni = 0; ni < asm_.dofs.n_nodes; ++ni)
      for (int c = 0; c < 6; ++c) {
        int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        if (fi >= 0) field[static_cast<std::size_t>(ni)](c) = vals[i] * modal.modes_free[i](fi);
      }
    out.modal_disp.push_back(std::move(field));
    if (opt.stations > 0) {
      CaseResult cm = recover_case(model, asm_, cl0, VecX(vals[i] * modal.modes_free[i]), sopt,
                                   static_cast<Id>(i + 1), "mode");
      out.modal_bars.push_back(std::move(cm.bars));
      out.modal_reactions.push_back(std::move(cm.reactions));
    }
  }

  auto combine = [&](auto value_of) {
    if (opt.rule == ComboRule::Srss) {
      double s = 0;
      for (std::size_t i = 0; i < nm; ++i) s += value_of(i) * value_of(i);
      return std::sqrt(s);
    }
    if (opt.rule == ComboRule::Abs) {
      double s = 0;
      for (std::size_t i = 0; i < nm; ++i) s += std::fabs(value_of(i));
      return s;
    }
    if (opt.rule == ComboRule::Nrl) {
      std::size_t imax = 0;
      for (std::size_t i = 1; i < nm; ++i)
        if (std::fabs(value_of(i)) > std::fabs(value_of(imax))) imax = i;
      double s = 0;
      for (std::size_t i = 0; i < nm; ++i)
        if (i != imax) s += value_of(i) * value_of(i);
      return nm == 0 ? 0.0 : std::fabs(value_of(imax)) + std::sqrt(s);
    }
    if (opt.rule == ComboRule::TenPercent) {
      double s = 0;
      for (std::size_t i = 0; i < nm; ++i) s += value_of(i) * value_of(i);
      for (std::size_t i = 0; i < nm; ++i)
        for (std::size_t j = i + 1; j < nm; ++j)
          if ((modal.omegas[j] - modal.omegas[i]) / modal.omegas[i] <= 0.10)
            s += 2 * std::fabs(value_of(i) * value_of(j));
      return std::sqrt(s);
    }
    double s = 0;
    for (std::size_t i = 0; i < nm; ++i)
      for (std::size_t j = 0; j < nm; ++j) {
        double rho = opt.rule == ComboRule::Cqc
                         ? cqc_rho(modal.omegas[i], modal.omegas[j], out.zetas[i], out.zetas[j])
                         : double_sum_rho(modal.omegas[i], modal.omegas[j], out.zetas[i],
                                          out.zetas[j], opt.duration);
        s += rho * value_of(i) * value_of(j);
      }
    return std::sqrt(std::max(0.0, s));
  };

  out.combined_disp.assign(nn, Vec6::Zero());
  for (std::size_t ni = 0; ni < nn; ++ni)
    for (int c = 0; c < 6; ++c)
      out.combined_disp[ni](c) = combine([&](std::size_t i) { return out.modal_disp[i][ni](c); });
  out.base_shear = combine([&](std::size_t i) { return out.modal_base_shear[i]; });

  if (opt.stations > 0 && !out.modal_reactions.empty()) {
    out.combined_reactions = out.modal_reactions.front();  // node template
    for (std::size_t r = 0; r < out.combined_reactions.size(); ++r)
      for (int c = 0; c < 6; ++c)
        out.combined_reactions[r].R(c) = combine([&](std::size_t i) {
          return r < out.modal_reactions[i].size() ? out.modal_reactions[i][r].R(c) : 0.0;
        });
  }
  if (opt.stations > 0 && !out.modal_bars.empty()) {
    // Station-level combination: same rule, per station and per component,
    // componentwise >= 0 (the nodal convention; signs live in modal_bars).
    out.combined_bars = out.modal_bars.front();  // shape template
    for (std::size_t e = 0; e < out.combined_bars.size(); ++e)
      for (std::size_t s = 0; s < out.combined_bars[e].stations.size(); ++s) {
        auto comb = [&](auto get) {
          return combine([&](std::size_t i) { return get(out.modal_bars[i][e].stations[s]); });
        };
        Station& st = out.combined_bars[e].stations[s];
        st.N = comb([](const Station& q) { return q.N; });
        st.Vy = comb([](const Station& q) { return q.Vy; });
        st.Vz = comb([](const Station& q) { return q.Vz; });
        st.T = comb([](const Station& q) { return q.T; });
        st.My = comb([](const Station& q) { return q.My; });
        st.Mz = comb([](const Station& q) { return q.Mz; });
        for (int c = 0; c < 3; ++c)
          st.defl(c) = comb([c](const Station& q) { return q.defl(c); });
        for (int c = 0; c < 3; ++c)
          st.rot(c) = comb([c](const Station& q) { return q.rot(c); });
      }
    for (std::size_t e = 0; e < out.combined_bars.size(); ++e) {
      for (int c = 0; c < 12; ++c)
        out.combined_bars[e].end_forces_local(c) = combine(
            [&](std::size_t i) { return out.modal_bars[i][e].end_forces_local(c); });
      for (std::size_t s = 0; s < out.combined_bars[e].soil_p.size(); ++s) {
        for (int c = 0; c < 3; ++c)
          out.combined_bars[e].soil_p[s](c) = combine(
              [&](std::size_t i) { return out.modal_bars[i][e].soil_p[s](c); });
        out.combined_bars[e].soil_m[s] =
            combine([&](std::size_t i) { return out.modal_bars[i][e].soil_m[s]; });
      }
    }
  }

  if (opt.missing_mass) {
    // along a vector: the rigid translations' mass is m I, so the total is
    // sum v_c^2 m_c and the effective mass of a mode Gamma_v^2
    double m_tot = 0, m_eff = 0;
    if (use_vec) {
      for (int c = 0; c < 3; ++c) m_tot += v(c) * v(c) * modal.total_mass[static_cast<std::size_t>(c)];
      for (double g : out.participation) m_eff += g * g;
    } else {
      m_tot = modal.total_mass[static_cast<std::size_t>(opt.dir)];
      for (std::size_t i = 0; i < nm; ++i) m_eff += modal.effective_mass[i][static_cast<std::size_t>(opt.dir)];
    }
    double dm = std::max(0.0, m_tot - m_eff);
    double zpa = spectrum.sa(0.0, out.zetas.empty() ? 0.05 : out.zetas[0]);
    out.missing_mass_shear = dm * zpa;
    out.base_shear = std::hypot(out.base_shear, out.missing_mass_shear);
  }
  return out;
}

std::vector<Vec6> combine_directional_srss(const std::vector<std::vector<Vec6>>& dirs) {
  require(!dirs.empty(), "directional combination: no inputs");
  std::vector<Vec6> out(dirs[0].size(), Vec6::Zero());
  for (const auto& d : dirs)
    for (std::size_t i = 0; i < out.size(); ++i)
      for (int c = 0; c < 6; ++c) out[i](c) += d[i](c) * d[i](c);
  for (auto& v : out)
    for (int c = 0; c < 6; ++c) v(c) = std::sqrt(v(c));
  return out;
}

std::vector<Vec6> combine_directional_30pct(const std::vector<std::vector<Vec6>>& dirs) {
  require(!dirs.empty(), "directional combination: no inputs");
  std::vector<Vec6> out(dirs[0].size(), Vec6::Zero());
  for (std::size_t g = 0; g < dirs.size(); ++g) {
    for (std::size_t i = 0; i < out.size(); ++i)
      for (int c = 0; c < 6; ++c) {
        double v = 0;
        for (std::size_t d = 0; d < dirs.size(); ++d)
          v += (d == g ? 1.0 : 0.3) * dirs[d][i](c);
        out[i](c) = std::max(out[i](c), v);
      }
  }
  return out;
}

std::vector<double> response_spectrum_of(const std::vector<double>& accel, double dt, const std::vector<double>& periods,
                                         double zeta) {
  require(accel.size() >= 2 && dt > 0, "response spectrum: an acceleration history and dt > 0");
  require(zeta >= 0 && zeta < 1, "response spectrum: 0 <= zeta < 1");
  std::vector<double> sa;
  for (double T : periods) {
    require(T > 0, "response spectrum: periods > 0");
    const double w = 2 * kPi / T;
    double u = 0, v = 0, umax = 0;
    for (std::size_t i = 1; i < accel.size(); ++i) {
      oscillator_step(w, zeta, dt, -accel[i - 1], -accel[i], u, v);
      umax = std::max(umax, std::fabs(u));
    }
    sa.push_back(w * w * umax);
  }
  return sa;
}

FloorSpectrum floor_response_spectrum(const Model& model, const ModalResults& modal, Id node, int comp, int dir,
                                      const std::vector<double>& ground_accel, double dt, const DampingSpec& damping,
                                      const std::vector<double>& periods, double zeta_floor) {
  require(dir >= 0 && dir < 3 && comp >= 0 && comp < 6, "floor spectrum: direction 0..2, component 0..5");
  require(!modal.modes_free.empty(), "floor spectrum: no modes");
  const Assembled asm_ = assemble_elements(model, true);
  const int f = asm_.dofs.free_of(model.node_index(node), comp);
  const std::size_t nm = modal.omegas.size(), ns = ground_accel.size();
  FloorSpectrum out;
  out.accel.assign(ns, 0.0);
  // relative modal responses: q'' + 2 z w q' + w^2 q = -Gamma a_g (exact steps)
  for (std::size_t i = 0; i < nm; ++i) {
    const double w = modal.omegas[i], z = damping.zeta_for(w, i), g = modal.participation[i][static_cast<std::size_t>(dir)];
    const double phi = f >= 0 ? modal.modes_free[i](f) : 0.0;
    if (phi == 0 || g == 0) continue;
    double q = 0, qd = 0;
    for (std::size_t k = 1; k < ns; ++k) {
      oscillator_step(w, z, dt, -g * ground_accel[k - 1], -g * ground_accel[k], q, qd);
      const double qdd = -g * ground_accel[k] - 2 * z * w * qd - w * w * q;
      out.accel[k] += phi * qdd;
    }
  }
  // absolute: the ground's own acceleration along the observed translation
  if (comp == dir)
    for (std::size_t k = 0; k < ns; ++k) out.accel[k] += ground_accel[k];
  out.periods = periods;
  out.sa = response_spectrum_of(out.accel, dt, periods, zeta_floor);
  return out;
}

DampingSpec strain_energy_damping(const Model& model, const ModalResults& modal, const std::map<Id, double>& zeta_of_material,
                                  double zeta_links) {
  require(model.continua.empty(), "strain_energy_damping: continuum elements are not supported");
  const Assembled asm_ = assemble_elements(model, true);
  const DofMap& d = asm_.dofs;
  DampingSpec out;
  out.kind = DampingSpec::Kind::PerMode;
  auto zeta_mat = [&](Id m) {
    auto it = zeta_of_material.find(m);
    require(it != zeta_of_material.end(), "strain-energy damping: no damping ratio for material " + std::to_string(m));
    require(it->second >= 0 && it->second < 1, "strain-energy damping: 0 <= zeta < 1");
    return it->second;
  };
  for (std::size_t i = 0; i < modal.modes_free.size(); ++i) {
    const VecX& ph = modal.modes_free[i];
    auto val = [&](int ni, int c) {
      const int fi = d.free_of(ni, c);
      return fi >= 0 ? ph(fi) : 0.0;
    };
    double e = 0, ze = 0;
    for (std::size_t k = 0; k < asm_.elements.size(); ++k) {
      const Bar& b = model.bars[k];
      if (b.kind == BarKind::Cable) continue;
      const int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      Vec12 u;
      for (int c = 0; c < 6; ++c) {
        u(c) = val(ni, c);
        u(6 + c) = val(nj, c);
      }
      const double en = 0.5 * u.dot(bar_ke_assembled(model, asm_.elements[k]) * u);
      e += en;
      ze += zeta_mat(b.material) * en;
    }
    for (const ShellElem& se : asm_.shells) {
      Vec18 u;
      for (int k = 0; k < 3; ++k)
        for (int c = 0; c < 6; ++c) u(6 * k + c) = val(se.nodes[static_cast<std::size_t>(k)], c);
      const double en = 0.5 * u.dot(se.el.ke() * u);
      e += en;
      ze += zeta_mat(model.shells[se.index].material) * en;
    }
    for (const LinkElem& le : asm_.links) {
      Vec12 u;
      for (int c = 0; c < 6; ++c) {
        u(c) = val(le.ni, c);
        u(6 + c) = val(le.nj, c);
      }
      const double en = 0.5 * u.dot(le.K * u);
      e += en;
      ze += zeta_links * en;
    }
    require(e > 0, "strain-energy damping: a mode without strain energy");
    out.per_mode.push_back(ze / e);
  }
  return out;
}

}  // namespace keta

// Non-linear Newmark integration with hysteretic springs — see nldyn.h.
#include "keta/nldyn.h"

#include <cmath>

namespace keta {

ElastoPlastic::ElastoPlastic(double k, double fy, double alpha) : k_(k), fy_(fy), alpha_(alpha) {
  require(k > 0 && fy > 0 && alpha >= 0 && alpha < 1, "elasto-plastic spring: k, fy > 0, 0 <= alpha < 1");
}

void ElastoPlastic::state(double d, double& f, double& kt, double& dp, double& back) const {
  // elastic predictor on the elastic part k (1 - alpha) in parallel with alpha k
  const double ke = k_ * (1 - alpha_), kh = alpha_ * k_;
  dp = dp_;
  back = back_;
  double fe = ke * (d - dp);  // elasto-plastic branch force
  const double xi = fe - back;
  const double fyb = fy_ * (1 - alpha_);
  if (std::fabs(xi) > fyb) {
    // return mapping of the EPP branch with kinematic hardening of the back
    // force carried by the parallel spring: the parallel spring alpha k is
    // the hardening, the EPP branch keeps |fe| <= fy (1 - alpha)
    const double s = xi > 0 ? 1.0 : -1.0;
    fe = back + s * fyb;
    dp = d - fe / ke;
    kt = kh;
  } else {
    kt = k_;
  }
  f = fe + kh * d;
}

void ElastoPlastic::trial(double d, double& f, double& kt) const {
  double dp, back;
  state(d, f, kt, dp, back);
}

void ElastoPlastic::commit(double d) {
  double f, kt, dp, back;
  state(d, f, kt, dp, back);
  dp_ = dp;
  back_ = back;
}

Takeda::Takeda(double My, double phi_y, double Mu, double phi_u, double alpha)
    : My_(My), phiy_(phi_y), Mu_(Mu), phiu_(phi_u), alpha_(alpha) {
  require(My > 0 && phi_y > 0 && Mu > My && phi_u > phi_y && alpha >= 0,
          "takeda: 0 < My < Mu, 0 < phi_y < phi_u, alpha >= 0");
  k0_ = My / phi_y;
  kh_ = (Mu - My) / (phi_u - phi_y);
}

double Takeda::skeleton(double d, double& kt) const {
  const double a = std::fabs(d), s = d < 0 ? -1.0 : 1.0;
  if (a <= phiy_) {
    kt = k0_;
    return k0_ * d;
  }
  kt = kh_;
  return s * (My_ + kh_ * (a - phiy_));
}

// Walk the state s from (s.d, s.f) to the deformation d along the rules.
void Takeda::walk(State& s, double d, double& f, double& kt) const {
  if (d == s.d) {
    f = s.f;
    double k;
    if (s.on_skeleton) skeleton(s.d == 0 ? s.last_dir * 1e-300 : s.d, k);
    else k = k0_;
    kt = k;
    return;
  }
  const int dir = d > s.d ? 1 : -1;
  const int side = dir > 0 ? 0 : 1;
  const double sg = dir;
  // reversal: the departure point becomes an inner apex if off the skeleton
  if (s.last_dir != 0 && dir != s.last_dir && !s.on_skeleton && s.f * (-sg) > 0)
    s.apex[1 - side].push_back({s.d, s.f});
  s.last_dir = dir;
  for (int guard = 0; guard < 64; ++guard) {
    // 1. unloading toward zero force
    if (s.f * sg < 0) {
      const double ku = k0_ * std::pow(phiy_ / std::max(s.dmax, phiy_), alpha_);
      const double d0 = s.d - s.f / ku;  // zero-force deformation
      s.on_skeleton = false;
      if ((d - d0) * sg <= 0) {
        f = s.f + ku * (d - s.d);
        kt = ku;
        s.f = f;
        s.d = d;
        return;
      }
      s.d = d0;
      s.f = 0;
      continue;
    }
    // 2. on the skeleton beyond the peak: follow it
    if (s.on_skeleton) {
      double k;
      const double fs = skeleton(d, k);
      f = fs;
      kt = k;
      s.d = d;
      s.f = fs;
      s.dmax = std::max(s.dmax, std::fabs(d));
      if (std::fabs(d) > phiy_) {
        s.yielded[side] = true;
        s.peak[side] = {d, fs};
      }
      return;
    }
    // 3. reloading toward the next target: inner apex, then outer peak / yield
    std::vector<Pt>& ap = s.apex[side];
    while (!ap.empty() && (ap.back().d - s.d) * sg <= 0) ap.pop_back();  // already passed
    Pt tgt;
    bool to_skeleton = false;
    if (!ap.empty()) {
      tgt = ap.back();
    } else {
      tgt = s.yielded[side] ? s.peak[side] : Pt{sg * phiy_, sg * My_};
      to_skeleton = true;
    }
    const double k = (tgt.f - s.f) / (tgt.d - s.d);
    if ((d - tgt.d) * sg <= 0) {
      f = s.f + k * (d - s.d);
      kt = k;
      s.f = f;
      s.d = d;
      return;
    }
    s.d = tgt.d;
    s.f = tgt.f;
    if (to_skeleton) {
      s.on_skeleton = true;
      ap.clear();
    } else {
      ap.pop_back();
    }
  }
  fail("takeda: too many branch changes in one increment");
}

void Takeda::trial(double d, double& f, double& kt) const {
  State s = st_;
  walk(s, d, f, kt);
}

void Takeda::commit(double d) {
  double f, kt;
  walk(st_, d, f, kt);
}

double hysteresis_to_force(Hysteresis& h, double dc, double ft) {
  double fc, kc;
  h.trial(dc, fc, kc);
  if (fc == ft) {
    h.commit(dc);
    return dc;
  }
  const double dir = ft > fc ? 1.0 : -1.0;
  double step = std::max(std::fabs(ft - fc) / std::max(kc, 1e-300), 1e-12);
  double lo = dc, hi = dc;
  for (int i = 0; i < 200; ++i) {
    hi = dc + dir * step;
    double f, k;
    h.trial(hi, f, k);
    if ((f - ft) * dir >= 0) break;
    lo = hi;
    step *= 2;
  }
  for (int i = 0; i < 300 && std::fabs(hi - lo) > 1e-15 * std::max(1.0, std::fabs(hi)); ++i) {
    const double m = 0.5 * (lo + hi);
    double f, k;
    h.trial(m, f, k);
    ((f - ft) * dir < 0 ? lo : hi) = m;
  }
  const double d = 0.5 * (lo + hi);
  h.commit(d);
  return d;
}

TfpbBackbone tfpb_backbone(const TfpbParams& p) {
  require(p.mu1 < p.mu2 && p.mu2 < p.mu3, "tfpb: mu1 < mu2 < mu3 required");
  require(p.Re1 > 0 && p.Re1 < p.Re2 && p.Re2 <= p.Re3, "tfpb: Re1 < Re2 <= Re3 required");
  require(p.d1 > (p.mu3 - p.mu1) * p.Re1, "tfpb: d1 > (mu3 - mu1) Re1 required");
  TfpbBackbone b;
  b.F[0] = p.mu1;
  b.F[1] = p.mu2;
  b.F[2] = p.mu3;
  b.F[3] = p.mu2 + p.d2 / p.Re2;
  b.F[4] = p.mu3 + p.d3 / p.Re3;
  b.F[5] = p.mu1 + p.d1 / p.Re1 + p.d3 / p.Re3;
  const double flex[5] = {2 * p.Re1, p.Re1 + p.Re2, p.Re2 + p.Re3, p.Re1 + p.Re3, 2 * p.Re1};  // 1/k
  b.u[0] = 0;
  for (int i = 1; i < 6; ++i) b.u[i] = b.u[i - 1] + (b.F[i] - b.F[i - 1]) * flex[i - 1];
  return b;
}

double tfpb_displacement(const TfpbParams& p, double F) {
  const TfpbBackbone b = tfpb_backbone(p);
  const double a = std::fabs(F), s = F < 0 ? -1.0 : 1.0;
  require(a <= b.F[5] * (1 + 1e-12), "tfpb: force beyond the bearing capacity F6");
  if (a <= b.F[0]) return 0.0;
  for (int i = 1; i < 6; ++i)
    if (a <= b.F[i]) return s * (b.u[i - 1] + (a - b.F[i - 1]) / (b.F[i] - b.F[i - 1]) * (b.u[i] - b.u[i - 1]));
  return s * b.u[5];
}

NlHistory newmark_nonlinear(const NlSystem& sys, const std::function<VecX(double)>& p,
                            const NewmarkNlOptions& opt, const VecX& u0, const VecX& v0) {
  const int n = static_cast<int>(sys.M.rows());
  require(n > 0 && n <= 2000, "newmark_nonlinear: 1..2000 DOFs (dense)");
  require(sys.C.rows() == n && sys.K.rows() == n, "newmark_nonlinear: matrix sizes");
  require(opt.dt > 0 && opt.steps > 0, "newmark_nonlinear: dt and steps must be positive");
  std::vector<std::unique_ptr<Hysteresis>> law;
  for (const NlSpring& s : sys.springs) {
    require(s.law != nullptr && s.i >= 0 && s.i < n && s.j < n, "newmark_nonlinear: bad spring");
    law.push_back(s.law->clone());
  }
  auto spring_forces = [&](const VecX& u, VecX& fs, MatX* kt, std::vector<double>* fsp) {
    fs.setZero(n);
    if (kt) kt->setZero(n, n);
    if (fsp) fsp->assign(sys.springs.size(), 0.0);
    for (std::size_t k = 0; k < sys.springs.size(); ++k) {
      const NlSpring& s = sys.springs[k];
      const double d = u(s.i) - (s.j >= 0 ? u(s.j) : 0.0);
      double f, t;
      law[k]->trial(d, f, t);
      if (fsp) (*fsp)[k] = f;
      fs(s.i) += f;
      if (s.j >= 0) fs(s.j) -= f;
      if (kt) {
        (*kt)(s.i, s.i) += t;
        if (s.j >= 0) {
          (*kt)(s.j, s.j) += t;
          (*kt)(s.i, s.j) -= t;
          (*kt)(s.j, s.i) -= t;
        }
      }
    }
  };
  NlHistory h;
  VecX u = u0.size() ? u0 : VecX::Zero(n), v = v0.size() ? v0 : VecX::Zero(n);
  VecX fs;
  spring_forces(u, fs, nullptr, nullptr);
  // initial acceleration from equilibrium
  VecX a = sys.M.ldlt().solve(p(0.0) - sys.C * v - sys.K * u - fs);
  std::vector<double> fsp;
  spring_forces(u, fs, nullptr, &fsp);
  h.t.push_back(0);
  h.u.push_back(u);
  h.v.push_back(v);
  h.a.push_back(a);
  h.spring_force.push_back(fsp);
  const double dt = opt.dt, g = opt.gamma, b = opt.beta;
  for (int step = 1; step <= opt.steps; ++step) {
    const double t = step * dt;
    const VecX pt = p(t);
    VecX un = u + dt * v + dt * dt * (0.5 - b) * a;  // predictor (a(n+1) = 0)
    auto kin = [&](const VecX& x, VecX& an, VecX& vn) {
      an = (x - u - dt * v - dt * dt * (0.5 - b) * a) / (b * dt * dt);
      vn = v + dt * ((1 - g) * a + g * an);
    };
    VecX an, vn;
    MatX kt;
    bool done = false;
    for (int it = 0; it < opt.max_iterations; ++it) {
      kin(un, an, vn);
      spring_forces(un, fs, &kt, nullptr);
      const VecX r = pt - sys.M * an - sys.C * vn - sys.K * un - fs;
      const double scale = pt.lpNorm<Eigen::Infinity>() + (sys.M * an).lpNorm<Eigen::Infinity>() +
                           (sys.K * un + fs).lpNorm<Eigen::Infinity>() + 1e-300;
      ++h.newton_iterations;
      if (r.lpNorm<Eigen::Infinity>() <= opt.tol * scale) {
        done = true;
        break;
      }
      const MatX J = sys.M / (b * dt * dt) + sys.C * (g / (b * dt)) + sys.K + kt;
      un += J.fullPivLu().solve(r);
    }
    require(done, "newmark_nonlinear: Newton did not converge at t = " + std::to_string(t));
    kin(un, an, vn);
    spring_forces(un, fs, nullptr, &fsp);
    for (std::size_t k = 0; k < sys.springs.size(); ++k) {
      const NlSpring& s = sys.springs[k];
      law[k]->commit(un(s.i) - (s.j >= 0 ? un(s.j) : 0.0));
    }
    u = un;
    v = vn;
    a = an;
    h.t.push_back(t);
    h.u.push_back(u);
    h.v.push_back(v);
    h.a.push_back(a);
    h.spring_force.push_back(fsp);
  }
  return h;
}

}  // namespace keta

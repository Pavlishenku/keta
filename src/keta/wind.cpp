#include "keta/wind.h"

#include <complex>

namespace keta {

// ---------------- profile ----------------

double WindProfile::vm(double z) const {
  double ze = std::max(z, zmin);
  return kr() * std::log(ze / z0) * co * vb;
}
double WindProfile::iv(double z) const {
  double ze = std::max(z, zmin);
  return kI / (co * std::log(ze / z0));
}
double WindProfile::qp(double z) const {
  double v = vm(z);
  return 0.5 * rho * v * v * (1.0 + 7.0 * iv(z));
}

LoadCase mean_wind_case(const Model& model, const WindProfile& wind, const AeroSection& aero,
                        const Vec3& wind_dir, Id case_id, const std::string& name) {
  Vec3 d = wind_dir;
  d.z() = 0;
  require(d.norm() > 1e-9, "wind: direction must have a horizontal component");
  d.normalize();
  Vec3 lift(0, 0, 1);
  LoadCase lc;
  lc.id = case_id;
  lc.name = name;
  for (const Bar& b : model.bars) {
    if (b.kind == BarKind::Cable) continue;
    Vec3 x1 = model.node(b.ni).x, x2 = model.node(b.nj).x;
    // yaw factor: only the component of the flow perpendicular to the bar loads it
    Vec3 ex = (x2 - x1).normalized();
    double sin2 = 1.0 - std::pow(ex.dot(d), 2);
    auto q_at = [&](double z) { return wind.qp(z) * sin2; };
    double q1 = q_at(x1.z()), q2 = q_at(x2.z());
    Vec3 w1 = q1 * aero.B * (aero.cd * d + aero.cl * lift);
    Vec3 w2 = q2 * aero.B * (aero.cd * d + aero.cl * lift);
    lc.loads.push_back(Load::linear(b.id, w1, w2, LoadFrame::Global));
    if (aero.cm != 0) {
      // distributed pitching moment lumped consistently at the ends, about
      // lift x drag = z x d: positive nose-up (the upwind edge rises)
      double L = (x2 - x1).norm();
      Vec3 m_axis = lift.cross(d);
      double m1 = q1 * aero.cm * aero.B * aero.B, m2 = q2 * aero.cm * aero.B * aero.B;
      lc.loads.push_back(
          Load::point(b.id, 0.0, Vec3::Zero(), (L * (m1 / 3 + m2 / 6)) * m_axis));
      lc.loads.push_back(
          Load::point(b.id, L, Vec3::Zero(), (L * (m1 / 6 + m2 / 3)) * m_axis));
    }
  }
  return lc;
}

// ---------------- vortex ----------------

double vortex_critical_speed(double n1, double bwidth, double st) {
  require(n1 > 0 && bwidth > 0 && st > 0, "vortex: bad input");
  return n1 * bwidth / st;
}
double scruton_number(double m, double zeta, double rho, double bwidth) {
  return 4.0 * kPi * m * zeta / (rho * bwidth * bwidth);
}
double vortex_amplitude_ratio(double st, double sc, double K, double Kw, double clat) {
  require(st > 0 && sc > 0, "vortex: bad input");
  return K * Kw * clat / (st * st * sc);
}

// ---------------- buffeting ----------------

double u_spectrum(double n, const BuffetingOptions& o) {
  require(n > 0, "spectrum: n must be positive");
  double su = o.Iu * o.U;
  double fl = n * o.Lu / o.U;
  if (o.von_karman)
    return su * su / n * 4.0 * fl / std::pow(1.0 + 70.8 * fl * fl, 5.0 / 6.0);
  // EN 1991-1-4 Annex B (B.2), the Kaimal-type form on the integral length
  // scale: n S / sigma^2 = 6.8 f_L / (1 + 10.2 f_L)^(5/3), whose integral over
  // dn is sigma^2 exactly (6.8 / 10.2 x 3/2 = 1). The former 200 f / (1 +
  // 50 f / 3)^(5/3) was Kaimal's u*-normalised shape and integrated to 18 sigma^2.
  return su * su / n * 6.8 * fl / std::pow(1.0 + 10.2 * fl, 5.0 / 3.0);
}

BuffetingResults spectral_buffeting(const Model& model, const ModalResults& modal,
                                    const AeroSection& aero, const Vec3& wind_dir,
                                    const DampingSpec& damping, const BuffetingOptions& opt) {
  Vec3 d = wind_dir;
  d.z() = 0;
  require(d.norm() > 1e-9, "buffeting: bad wind direction");
  d.normalize();
  model.validate();
  Assembled asm_ = assemble_elements(model);

  // Strip data per bar: midpoint, length, per-mode modal ordinate (drag dir).
  struct Strip {
    Vec3 mid;
    double L;
    std::vector<double> phi;  // per mode
  };
  std::vector<Strip> strips;
  for (const Bar& b : model.bars) {
    if (b.kind == BarKind::Cable) continue;
    Strip s;
    Vec3 x1 = model.node(b.ni).x, x2 = model.node(b.nj).x;
    s.mid = 0.5 * (x1 + x2);
    s.L = (x2 - x1).norm();
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (std::size_t i = 0; i < modal.omegas.size(); ++i) {
      Vec3 phi_mid = Vec3::Zero();
      for (int c = 0; c < 3; ++c) {
        int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nj, c))];
        double v1 = f1 >= 0 ? modal.modes_free[i](f1) : 0.0;
        double v2 = f2 >= 0 ? modal.modes_free[i](f2) : 0.0;
        phi_mid(c) = 0.5 * (v1 + v2);
      }
      s.phi.push_back(phi_mid.dot(d));
    }
    strips.push_back(std::move(s));
  }
  require(!strips.empty(), "buffeting: no bars");

  // Quasi-steady drag admittance: dF = rho U C_D B u(t) per unit length.
  double c_force = opt.rho * opt.U * aero.cd * aero.B;

  // Log frequency grid around the modal band.
  double f_lo = modal.frequencies_hz.front() * 0.05;
  double f_hi = modal.frequencies_hz.back() * 5.0;
  std::vector<double> grid;
  for (int i = 0; i <= opt.n_freq; ++i)
    grid.push_back(f_lo * std::pow(f_hi / f_lo, static_cast<double>(i) / opt.n_freq));

  BuffetingResults out;
  std::size_t nm = modal.omegas.size();
  for (std::size_t i = 0; i < nm; ++i) {
    double om_i = modal.omegas[i];
    double z = damping.zeta_for(om_i, i);
    require(z > 0, "buffeting: every mode needs positive damping");
    // integrate response PSD
    double var = 0, var_n2 = 0, prev_n = grid[0], prev_s = 0, prev_sn2 = 0;
    for (std::size_t g = 0; g < grid.size(); ++g) {
      double n = grid[g];
      double su = u_spectrum(n, opt);
      // joint acceptance with Davenport coherence
      double J = 0;
      for (const Strip& sa : strips)
        for (const Strip& sb : strips) {
          double dist = (sa.mid - sb.mid).norm();
          double coh = std::exp(-opt.coherence_c * n * dist / opt.U);
          J += sa.phi[i] * sb.phi[i] * sa.L * sb.L * coh;
        }
      double SQ = c_force * c_force * su * J;  // modal force PSD
      double om = 2 * kPi * n;
      std::complex<double> H(om_i * om_i - om * om, 2 * z * om_i * om);
      double H2 = 1.0 / std::norm(H);
      double Sq = H2 * SQ;
      if (g > 0) {
        var += 0.5 * (prev_s + Sq) * (n - prev_n);
        var_n2 += 0.5 * (prev_sn2 + Sq * n * n) * (n - prev_n);
      }
      prev_n = n;
      prev_s = Sq;
      prev_sn2 = Sq * n * n;
    }
    BuffetingModeResult mr;
    mr.sigma_q = std::sqrt(var);
    mr.nu = var > 0 ? std::sqrt(var_n2 / var) : 0.0;
    double lt = std::log(std::max(1.0000001, mr.nu * opt.gust_duration));
    mr.peak_factor = std::sqrt(2 * lt) + 0.5772156649 / std::sqrt(2 * lt);
    out.modes.push_back(mr);
  }

  // Nodal RMS field: SRSS over modes (well-separated assumption documented).
  std::size_t nn = modal.modes.front().size();
  out.rms_disp.assign(nn, Vec6::Zero());
  out.peak_disp.assign(nn, Vec6::Zero());
  for (std::size_t ni = 0; ni < nn; ++ni)
    for (int c = 0; c < 6; ++c) {
      double s2 = 0, p2 = 0;
      for (std::size_t i = 0; i < nm; ++i) {
        int fi = asm_.dofs.free_ix[static_cast<std::size_t>(
            asm_.dofs.full(static_cast<int>(ni), c))];
        double phi = fi >= 0 ? modal.modes_free[i](fi) : 0.0;
        double v = out.modes[i].sigma_q * phi;
        s2 += v * v;
        double pv = v * out.modes[i].peak_factor;
        p2 += pv * pv;
      }
      out.rms_disp[ni](c) = std::sqrt(s2);
      out.peak_disp[ni](c) = std::sqrt(p2);
    }
  return out;
}

// ---------------- flutter ----------------

std::complex<double> theodorsen_c(double k) {
  require(k > 0, "theodorsen: k must be positive");
  std::complex<double> ik(0.0, 1.0 / k);
  return 1.0 - 0.165 / (1.0 - 0.0455 * ik) - 0.335 / (1.0 - 0.3 * ik);
}

namespace {

// Branch damping/frequency at speed U by frozen-k iteration on the quadratic
// eigenproblem det(p^2 M + p C + K) = 0, aero from Theodorsen at k = om b / U.
struct Branch {
  double omega, zeta;
};

std::vector<Branch> flutter_branches(const FlutterSection& s, double U) {
  using C = std::complex<double>;
  double rho = s.rho, b = s.b, a = s.a;
  auto eig_at = [&](double om_guess) {
    double k = std::max(1e-4, om_guess * b / U);
    C Ck = theodorsen_c(k);
    // Theodorsen aero (h down, alpha nose-up, L up, M about EA nose-up):
    // L = pi rho b^2 (hdd + U ad - b a add) + 2 pi rho U b C (hd + U a + b(1/2 - a) ad)
    // M = pi rho b^2 (b a hdd - U b (1/2 - a) ad - b^2(1/8 + a^2) add)
    //     + 2 pi rho U b^2 (a + 1/2) C (hd + U a + b(1/2 - a) ad)
    double pb2 = kPi * rho * b * b;
    C two = 2.0 * kPi * rho * U * b * Ck;
    // Mass (coefficients of hdd, add) entering  Ms qdd + ... + L-row / -M-row:
    Eigen::Matrix2cd Ma, Ca, Ka;
    Ma(0, 0) = pb2;
    Ma(0, 1) = -pb2 * b * a;
    Ma(1, 0) = -pb2 * b * a;
    Ma(1, 1) = pb2 * b * b * (0.125 + a * a);
    Ca(0, 0) = two;
    Ca(0, 1) = pb2 * U + two * b * (0.5 - a);
    Ca(1, 0) = -two * b * (a + 0.5);
    Ca(1, 1) = pb2 * U * b * (0.5 - a) - two * b * (a + 0.5) * b * (0.5 - a);
    Ka(0, 0) = 0;
    Ka(0, 1) = two * U;
    Ka(1, 0) = 0;
    Ka(1, 1) = -two * U * b * (a + 0.5);
    // structural
    Eigen::Matrix2cd Ms = Eigen::Matrix2cd::Zero(), Cs = Eigen::Matrix2cd::Zero(),
                     Ks = Eigen::Matrix2cd::Zero();
    Ms(0, 0) = s.m;
    Ms(0, 1) = s.S;
    Ms(1, 0) = s.S;
    Ms(1, 1) = s.I;
    Cs(0, 0) = 2 * s.zeta_h * s.omega_h * s.m;
    Cs(1, 1) = 2 * s.zeta_a * s.omega_a * s.I;
    Ks(0, 0) = s.m * s.omega_h * s.omega_h;
    Ks(1, 1) = s.I * s.omega_a * s.omega_a;
    Eigen::Matrix2cd M = Ms + Ma, Cm = Cs + Ca, K = Ks + Ka;
    // linearize: [0 I; -M^-1 K  -M^-1 C]
    Eigen::Matrix2cd Minv = M.inverse();
    Eigen::Matrix4cd A = Eigen::Matrix4cd::Zero();
    A.block<2, 2>(0, 2) = Eigen::Matrix2cd::Identity();
    A.block<2, 2>(2, 0) = -Minv * K;
    A.block<2, 2>(2, 2) = -Minv * Cm;
    Eigen::ComplexEigenSolver<Eigen::Matrix4cd> es(A);
    return es.eigenvalues();
  };

  std::vector<Branch> out;
  for (double om0 : {s.omega_h, s.omega_a}) {
    double om = om0;
    C p(0, om0);
    for (int it = 0; it < 60; ++it) {
      auto ev = eig_at(om);
      // pick the eigenvalue with Im > 0 closest to the current guess
      double best = 1e300;
      C pick(0, om);
      for (int i = 0; i < 4; ++i) {
        C pi = ev(i);
        if (pi.imag() <= 0) continue;
        double dist = std::fabs(pi.imag() - om);
        if (dist < best) {
          best = dist;
          pick = pi;
        }
      }
      double om_new = pick.imag();
      p = pick;
      if (std::fabs(om_new - om) < 1e-9 * om0) {
        om = om_new;
        break;
      }
      om = om_new;
    }
    Branch br;
    br.omega = om;
    br.zeta = -p.real() / std::abs(p);
    out.push_back(br);
  }
  return out;
}

}  // namespace

FlutterResult flutter_speed(const FlutterSection& s, double u_min, double u_max, int n_sweep) {
  require(s.m > 0 && s.I > 0 && s.b > 0 && s.omega_h > 0 && s.omega_a > 0, "flutter: bad section");
  FlutterResult out;
  auto min_zeta = [&](double U) {
    auto br = flutter_branches(s, U);
    double z = 1e300;
    double om = 0;
    for (const Branch& b : br)
      if (b.zeta < z) {
        z = b.zeta;
        om = b.omega;
      }
    return std::make_pair(z, om);
  };
  double u_prev = u_min;
  auto [z_prev, om_prev] = min_zeta(u_min);
  for (int i = 1; i <= n_sweep; ++i) {
    double U = u_min + (u_max - u_min) * i / n_sweep;
    auto [z, om] = min_zeta(U);
    if (z_prev >= 0 && z < 0) {
      // bisect the crossing
      double lo = u_prev, hi = U;
      for (int it = 0; it < 60; ++it) {
        double mid = 0.5 * (lo + hi);
        auto [zm, omm] = min_zeta(mid);
        (zm >= 0 ? lo : hi) = mid;
        om = omm;
      }
      out.U_flutter = 0.5 * (lo + hi);
      out.omega_flutter = om;
      out.found = true;
      return out;
    }
    u_prev = U;
    z_prev = z;
    om_prev = om;
  }
  return out;  // no crossing in the sweep
}

}  // namespace keta

// Force-based spatial fiber element and the bar basic system, shared by the
// fiber frame driver (fiber_frame.cpp) and the shell plasticity driver
// (shell_plastic.cpp, nonlinear bars). Internal header: implementation
// detail of those drivers, not a public API.
#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "keta/assembly.h"
#include "keta/fiber_frame.h"
#include "keta/fiber_hyst.h"

namespace keta::fiberel {

using Vec5 = Eigen::Matrix<double, 5, 1>;
using Mat5 = Eigen::Matrix<double, 5, 5>;
using Mat35 = Eigen::Matrix<double, 3, 5>;
using Mat612 = Eigen::Matrix<double, 6, 12>;

// span-load statics of the basic system (defined below)
struct SpanParticular {
  std::vector<Vec3> sp;           // per abscissa, per unit load factor
  Vec12 p_ss = Vec12::Zero();     // basic-system reactions on the member, local
  // torsion: T(x) = q_T + T_p(x); tor_int = int_0^L T_p dx (the elastic
  // twist it produces is tor_int / GJ)
  double tor_int = 0;
};
inline Mat612 basic_transform6(double L);
template <class LoadsUpto>
SpanParticular span_particular_core(double L, const std::vector<double>& xs, const std::vector<char>* incl,
                                    LoadsUpto&& loads_upto, double torque_lm = 0);
inline void release_statics(SpanParticular& sp, const std::array<bool, 12>& rl, double L,
                            const std::vector<double>& xs);
// Force interpolation b(x): rows (N, My, Mz) on q = (N, Myi, Myj, Mzi, Mzj).
inline Mat35 b_basic(double x, double L) {
  Mat35 b = Mat35::Zero();
  const double xi = x / L;
  b(0, 0) = 1.0;
  b(1, 1) = 1.0 - xi;
  b(1, 2) = -xi;
  b(2, 3) = 1.0 - xi;
  b(2, 4) = -xi;
  return b;
}

// ---------------- spatial force-based fiber element ----------------

struct FiberElement3D {
  const HystFiberSection* sec = nullptr;
  std::vector<const HystFiberSection*> secs;  // section per point (tapered member)
  std::vector<Mat3> k0s;                      // their sup-tangent stiffnesses
  // span-load PARTICULAR resultants (N, My, Mz) at the points: s(x) = b(x) q
  // + sp(x), the statics of the basic system under the span loads (exact)
  std::vector<Vec3> sp;
  // released basic components (q_c = 0, v_c free): 0 N, 1 My_i, 2 My_j,
  // 3 Mz_i, 4 Mz_j, 5 T; SHEAR releases (Vy = 0: Mz_i + Mz_j = 0; Vz = 0:
  // My_i + My_j = 0). The admissible basic forces are q = Cq Q (Cq an
  // orthonormal basis of the constraints' null space, built by init), the
  // compatibility Cq^T (v - v(q)) = 0 — a released deformation follows.
  std::array<bool, 6> rel{};
  bool shear_y = false, shear_z = false;
  std::array<bool, 12> rel_local{};  // the bar's released local DOFs (set_releases)
  Eigen::MatrixXd Cq;          // 5 x m
  Mat5 Pq = Mat5::Identity();  // Cq Cq^T
  // imposed strains (thermal, prestrain): mechanical strain e - e0 in the
  // section, the total e in the compatibility; per point, current values
  std::vector<Vec3> e0;
  // a basic deformation shift (dislocations, the twist of span torques):
  // compatibility with v - vshift
  Vec6 vshift = Vec6::Zero();
  // torsion carried OUTSIDE the element (an exact linear torsion block, e.g.
  // on a torsional bed): q_T = 0 here
  bool tor_external = false;
  std::vector<double> breaks;  // integration breakpoints (abscissae), set before init
  std::vector<char> incl;      // per point: right limit at a breakpoint (point loads at x included)
  std::vector<int> sub_of;     // per point: its sub-interval of the composite rule
  std::vector<std::pair<double, double>> subs;  // the sub-intervals [a, b]
  double L = 0;
  double GJ_L = 0;             // elastic Saint-Venant torsion stiffness GJ/L
  std::vector<double> xg, wg;  // Lobatto abscissae [m] and weights * L [m]
  Mat3 k0 = Mat3::Zero();      // section sup-tangent stiffness
  double floor = 1e-6;
  // round-off amplitude of the target deformations (stretch [m], angles
  // [rad]): corotational kinematics carry eps |X| cancellation noise, below
  // which a relative tolerance on v is meaningless (0: exact targets)
  double v_noise_len = 0, v_noise_ang = 0;
  int max_iter = 30;
  double tol = 1e-10;

  Vec5 q = Vec5::Zero();       // fiber basic forces (warm start)
  double qT = 0;               // torsion (closed form)
  Mat6 kb = Mat6::Identity();  // 6x6 basic tangent
  Mat5 kb_ref = Mat5::Identity();  // elastic 5x5 basic stiffness (metric)
  std::vector<Vec3> d, s;      // (eps0, ky, kz) / (N, My, Mz) at the points
  std::vector<Mat3> fsk;       // section flexibilities at the last evaluation
  std::vector<std::vector<HystState>> committed;  // material memory per point

  // Force interpolation b(x): rows (N, My, Mz) on q = (N, Myi, Myj, Mzi, Mzj).
  Mat35 b_of(double x) const { return b_basic(x, L); }

  void init(const HystFiberSection& section, double length, double gj, int n_ip,
            double tangent_floor, int element_iter, double element_tol) {
    sec = &section;
    L = length;
    GJ_L = gj / length;
    floor = tangent_floor;
    max_iter = element_iter;
    tol = element_tol;
    std::vector<double> x01, w01;
    gauss_lobatto01(n_ip, x01, w01);
    xg.clear();
    wg.clear();
    // composite Lobatto over the sub-intervals between the breakpoints (span
    // load points / segment ends: the internal force field has a kink, a jump
    // of its curvature, or — under a point load — a JUMP there). A shared end
    // point is kept TWICE, the left and the right limit, each with its own
    // weight and material state (incl: the right limit, loads at x included)
    incl.clear();
    sub_of.clear();
    subs.clear();
    std::vector<double> cuts{0.0};
    for (double c : breaks)
      if (c > 1e-12 * L && c < L * (1 - 1e-12)) cuts.push_back(c);
    cuts.push_back(L);
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    for (std::size_t sgi = 0; sgi + 1 < cuts.size(); ++sgi) {
      const double a = cuts[sgi], h = cuts[sgi + 1] - a;
      subs.emplace_back(a, cuts[sgi + 1]);
      for (std::size_t i = 0; i < x01.size(); ++i) {
        sub_of.push_back(static_cast<int>(sgi));
        const bool last = i + 1 == x01.size();
        xg.push_back(i == 0 ? a : (last ? cuts[sgi + 1] : a + x01[i] * h));
        wg.push_back(w01[i] * h);
        incl.push_back(sgi > 0 && i == 0 ? 1 : 0);
      }
    }
    if (secs.size() != xg.size()) secs.assign(xg.size(), sec);
    k0s.clear();
    committed.clear();
    for (const HystFiberSection* sk : secs) {
      k0s.push_back(sk->sup_stiffness());
      Eigen::FullPivLU<Mat3> lu(k0s.back());
      require(lu.isInvertible(), "fiber_frame: degenerate fiber section (fibers on a line?)");
      committed.push_back(sk->virgin_states());
    }
    k0 = k0s.front();
    d.assign(xg.size(), Vec3::Zero());
    s.assign(xg.size(), Vec3::Zero());
    fsk.assign(xg.size(), Mat3::Zero());
    if (sp.size() != xg.size()) sp.assign(xg.size(), Vec3::Zero());
    if (e0.size() != xg.size()) e0.assign(xg.size(), Vec3::Zero());
    build_constraints();
    kb_ref = basic5_ref();
    kb.setZero();
    kb.block<5, 5>(0, 0) = free_inverse(kb_ref.inverse());
    kb(5, 5) = torsion_stiffness();
  }
  // releases from the bar's local DOF flags (before init)
  void set_releases(const std::array<bool, 12>& rl) {
    rel_local = rl;
    rel = {};
    rel[0] = rl[0] || rl[6];
    rel[1] = rl[4];
    rel[2] = rl[10];
    rel[3] = rl[5];
    rel[4] = rl[11];
    rel[5] = rl[3] || rl[9];
    shear_y = rl[1] || rl[7];
    shear_z = rl[2] || rl[8];
  }
  void build_constraints() {
    std::vector<Vec5> rows;
    for (int c = 0; c < 5; ++c)
      if (rel[static_cast<std::size_t>(c)]) {
        Vec5 r = Vec5::Zero();
        r(c) = 1;
        rows.push_back(r);
      }
    if (shear_y) rows.push_back((Vec5() << 0, 0, 0, 1, 1).finished());
    if (shear_z) rows.push_back((Vec5() << 0, 1, 1, 0, 0).finished());
    if (rows.empty()) {
      Cq = Eigen::MatrixXd::Identity(5, 5);
    } else {
      Eigen::MatrixXd A(static_cast<Eigen::Index>(rows.size()), 5);
      for (std::size_t i = 0; i < rows.size(); ++i) A.row(static_cast<Eigen::Index>(i)) = rows[i].transpose();
      Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
      Eigen::MatrixXd N = lu.kernel();
      if (lu.rank() == 5) N.resize(5, 0);
      // orthonormal basis
      Eigen::HouseholderQR<Eigen::MatrixXd> qr(N);
      Cq = N.cols() == 0 ? Eigen::MatrixXd(5, 0)
                         : Eigen::MatrixXd(qr.householderQ() * Eigen::MatrixXd::Identity(5, N.cols()));
    }
    Pq = Cq * Cq.transpose();
  }
  // a section per point (tapered member): set before init
  void set_sections(const std::vector<const HystFiberSection*>& per_point) { secs = per_point; }

  // the free-block inverse of a 5x5 basic flexibility, zero on the released
  // components (their deformation is free, their force zero)
  // (on the admissible forces q = Cq Q: Cq (Cq^T fb Cq)^-1 Cq^T)
  Mat5 free_inverse(const Mat5& fb) const {
    if (Cq.cols() == 0) return Mat5::Zero();
    const Eigen::MatrixXd fQ = Cq.transpose() * fb * Cq;
    Eigen::FullPivLU<Eigen::MatrixXd> lu(fQ);
    require(lu.isInvertible(), "fiber_frame: singular element flexibility");
    return Cq * lu.inverse() * Cq.transpose();
  }
  double torsion_stiffness() const { return (rel[5] || tor_external) ? 0.0 : GJ_L; }
  // elastic reference basic stiffness (metric), from the per-point sections
  Mat5 basic5_ref() const {
    Mat5 fb = Mat5::Zero();
    for (std::size_t k = 0; k < xg.size(); ++k) {
      Eigen::FullPivLU<Mat3> lu(k0s[k]);
      const Mat3 fs = lu.inverse();
      const Mat35 b = b_of(xg[k]);
      fb += wg[k] * (b.transpose() * fs * b);
    }
    Eigen::FullPivLU<Mat5> lub(fb);
    require(lub.isInvertible(), "fiber_frame: singular basic flexibility");
    return lub.inverse();
  }

  // One clean evaluation of the stations at the current q. Returns false when
  // a section target is outside its capacity branch (typed Error from the
  // section Newton); fills v_cur / fb / d / s on success.
  bool eval_stations(Vec5& v_cur, Mat5& fb) {
    v_cur.setZero();
    fb.setZero();
    try {
      for (std::size_t k = 0; k < xg.size(); ++k) {
        Mat35 b = b_of(xg[k]);
        Vec3 sk = b * q + sp[k];
        if (!sp_bed.empty()) sk += sp_bed[k];
        // Warm start from the previous strain plane at this station, with a
        // section tolerance one decade TIGHTER than the element metric: a
        // warm-started section Newton returns as soon as it meets its
        // tolerance (no quadratic overshoot), so the section error must not
        // saturate the element residual.
        Vec3 e = secs[k]->solve_state(sk, committed[k], d[k], 100, std::max(1e-13, 0.1 * tol));
        d[k] = e;
        s[k] = sk;
        Mat3 ks = secs[k]->tangent(e, committed[k]) + floor * k0s[k];
        Eigen::FullPivLU<Mat3> lu(ks);
        require(lu.isInvertible(), "fiber_frame: singular section tangent");
        Mat3 fs = lu.inverse();
        fsk[k] = fs;
        v_cur += wg[k] * (b.transpose() * (e + e0[k]));
        fb += wg[k] * (b.transpose() * fs * b);
      }
    } catch (const NonConvergence&) {
      throw;
    } catch (const Error&) {
      return false;
    }
    return true;
  }

  // Element state determination for a target basic deformation v (6). The 5
  // fiber modes iterate (Spacone); torsion is elastic closed form. Two
  // robustness devices for history-carrying sections:
  //   * a trial q outside the capacity of the CURRENT committed state is
  //     bisected back towards the last cleanly evaluated q — and when the
  //     warm-start q itself is infeasible (states advanced since last step),
  //     the iteration restarts from q = 0 (always feasible);
  //   * a q update that INCREASES the residual energy is halved (damped
  //     iteration — branch flips of a cracked section make the plain
  //     iteration oscillate).
  void determine(const Vec6& v_in) {
    const Vec6 v = v_in - vshift;
    qT = torsion_stiffness() * v(5);
    const Vec5 vt = Pq * v.head<5>();  // the admissible part of the target
    q = Pq * q;
    double e_v = vt.dot(kb_ref * vt);
    double e_floor = kb_ref(0, 0) * v_noise_len * v_noise_len;
    for (int c = 1; c < 5; ++c) e_floor += kb_ref(c, c) * v_noise_ang * v_noise_ang;
    Vec5 q_good = q;
    bool good_valid = false;
    double e_r_prev = std::numeric_limits<double>::max();
    int rescue = 0, damped = 0, stall = 0, it = 0;
    // best-so-far snapshot for the graceful floor acceptance below
    double best_rel = std::numeric_limits<double>::max();
    Vec5 best_q = q;
    Mat6 best_kb = kb;
    std::vector<Vec3> best_d, best_s;
    std::vector<Mat3> best_fs;
    while (true) {
      Vec5 v_cur;
      Mat5 fb;
      if (!eval_stations(v_cur, fb)) {
        if (++rescue > 60)
          throw NonConvergence("fiber_frame: section target outside capacity in element");
        if (!good_valid) {
          if (q.norm() == 0)
            throw NonConvergence("fiber_frame: section infeasible at zero element force");
          q.setZero();  // elastic restart: zero force is always feasible
        } else if ((q - q_good).norm() == 0) {
          throw NonConvergence("fiber_frame: section target outside capacity in element");
        } else {
          q = 0.5 * (q + q_good);
        }
        continue;
      }
      Mat5 kb5;
      try {
        kb5 = free_inverse(fb);  // the plain inverse without releases
      } catch (const Error&) {
        throw NonConvergence("fiber_frame: singular element flexibility");
      }
      kb.setZero();
      kb.block<5, 5>(0, 0) = kb5;
      kb(5, 5) = torsion_stiffness();
      const Vec5 res = Pq * (vt - v_cur);  // released: v follows, q stays admissible
      double e_r = res.dot(kb_ref * res);
      double rel = e_v > 0 ? std::sqrt(std::fabs(e_r / e_v))
                           : std::sqrt(std::fabs(e_r));
      // (the section strains are only as accurate as the section solve: its
      // absolute floor at small forces bounds the reachable compatibility)
      const Vec5 vsec = section_noise();
      const double e_sec = vsec.dot(kb_ref.diagonal().cwiseProduct(vsec));
      bool ok = e_v > 0 ? (rel < tol || e_r <= e_floor + e_sec) : (std::fabs(e_r) < 1e-30);
      if (ok) return;
      // Damped iteration FIRST, outside the stall accounting: a q update that
      // overshoots a branch flip (crack closing) is walked back half-way
      // towards the last accepted q until the residual drops below the
      // pre-overshoot level — that walk is progress, not stagnation.
      if (e_r > e_r_prev && good_valid && damped < 20 && (q - q_good).norm() > 0) {
        ++damped;
        q = 0.5 * (q + q_good);
        continue;
      }
      if (rel < 0.99 * best_rel) {
        best_rel = rel;
        best_q = q;
        best_kb = kb;
        best_d = d;
        best_s = s;
        best_fs = fsk;
        stall = 0;
      } else {
        ++stall;
      }
      // GRACEFUL FLOOR: the residual can stagnate at the solver-precision
      // floor (section tolerance x floored-flexibility amplification) a few
      // decades above `tol` at deeply plastified states. When the iteration
      // stops improving there, accept the best state instead of failing —
      // the global Newton judges equilibrium in force space anyway. States
      // that are genuinely non-equilibrable keep a rel of order 1 and throw.
      if (stall >= 6 || it >= max_iter) {
        // (best_d/best_s are always filled when best_rel is finite — the first
        // pass improves on the DBL_MAX sentinel — but stay defensive.)
        if (best_rel < 1e3 * tol && !best_d.empty()) {
          q = best_q;
          kb = best_kb;
          d = best_d;
          s = best_s;
          fsk = best_fs;
          return;
        }
        throw NonConvergence("fiber_frame: element state determination did not converge");
      }
      ++it;
      e_r_prev = e_r;
      q_good = q;
      good_valid = true;
      q += kb5 * res;
      if (!q.allFinite()) throw NonConvergence("fiber_frame: element forces diverged");
    }
  }

  // The deformation noise of the section solves at the current stations:
  // each section strain is solved to tol_s relative to max(|s|, its scale
  // 1e-4 k0) (HystFiberSection::solve_state), i.e. to de = 10 tol_s max(|s|,
  // scale) / k0 — summed through the compatibility integral.
  Vec3 strain_noise(std::size_t k) const {
    const double tol_s = std::max(1e-13, 0.1 * tol);
    Vec3 de;
    for (int c = 0; c < 3; ++c) {
      const double kc = k0s[k](c, c);
      de(c) = 10 * tol_s * std::max(std::fabs(s[k](c)), std::max(1.0, 1e-4 * kc)) / kc;
    }
    return de;
  }
  Vec5 section_noise() const {
    Vec5 vn = Vec5::Zero();
    for (std::size_t k = 0; k < xg.size(); ++k) vn += wg[k] * (b_of(xg[k]).cwiseAbs().transpose() * strain_noise(k));
    return vn;
  }

  // Compliance of the converged state to a change of the particular field:
  // with v held, a perturbation dsp of the span-load resultants changes the
  // basic forces by dq = -kb5 G(dsp), G(dsp) = sum_k w_k b_k^T fs_k dsp_k —
  // i.e. like a basic deformation shift G(dsp).
  Vec5 sp_compliance(const std::vector<Vec3>& dsp) const {
    Vec5 g = Vec5::Zero();
    for (std::size_t k = 0; k < xg.size(); ++k) g += wg[k] * (b_of(xg[k]).transpose() * (fsk[k] * dsp[k]));
    return Pq * g;
  }


  // ---------------- Winkler bed (linear kinematics) ----------------
  //
  // A bed of moduli (kx, ky, kz) [N/m per m] reacts to the member's own
  // displacement field: the chord interpolation of the end translations plus
  // the deflection relative to the chord, integrated EXACTLY from the
  // section strains (piecewise-polynomial interpolation through the points of
  // each sub-interval: delta'' = curvature, delta(0) = delta(L) = 0; the
  // axial delta' = eps0 - mean). Its pressure p = -k u is a span load: its
  // particular resultants and basic-system reactions come from the same
  // statics as the applied loads (exact integration of the interpolated
  // pressure). The coupled element equations — compatibility in q and the bed
  // law at the points — are solved by Newton to round-off, the element
  // tangent being their exact linearisation (a collocation of the beam on
  // elastic foundation: it converges with fiber_points, it is not the
  // closed-form BEF solution).
  Vec3 kbed = Vec3::Zero();
  Eigen::MatrixXd Sop, Rop, Dop, Cop;  // see init_bed
  Eigen::VectorXd Pbed;                 // bed pressures at the points (px, py, pz), warm start
  std::vector<Vec3> sp_bed;             // their particular resultants
  bool has_bed() const { return kbed.squaredNorm() > 0; }

  // Gauss-Legendre on [0, 1] (Newton on P_n, exact to degree 2n - 1).
  static void gauss_legendre01(int n, std::vector<double>& x, std::vector<double>& w) {
    x.assign(static_cast<std::size_t>(n), 0.0);
    w.assign(static_cast<std::size_t>(n), 0.0);
    for (int i = 0; i < n; ++i) {
      double t = std::cos(3.14159265358979323846 * (i + 0.75) / (n + 0.5));
      double dp = 1;
      for (int it = 0; it < 100; ++it) {
        double p0 = 1, p1 = t;
        for (int k = 2; k <= n; ++k) {
          const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
          p0 = p1;
          p1 = p2;
        }
        dp = n * (t * p1 - p0) / (t * t - 1);
        const double dt = p1 / dp;
        t -= dt;
        if (std::fabs(dt) < 1e-16) break;
      }
      x[static_cast<std::size_t>(i)] = 0.5 * (1 - t);
      w[static_cast<std::size_t>(i)] = 1.0 / ((1 - t * t) * dp * dp);
    }
  }

  // The bed operators, after init (the abscissae final):
  //   Sop (3np x 3np): particular resultants (N, My, Mz) at the points per
  //     unit pressure (px, py, pz) at the points;
  //   Rop (12 x 3np): the basic-system reactions per unit pressure;
  //   Dop (3np x 3np): displacements relative to the chord (dx, dy, dz) at
  //     the points per section strain (eps0, ky, kz) at the points;
  //   Cop (3np x 12): the chord displacements at the points from the local
  //     end displacements.
  // Curvature-to-deflection signs from the basic system (ty = -w', tz = v'):
  // the chord-relative end slopes of delta must be the basic rotations
  // v = int b^T e of the same strains (ty_i = int (1 - xi) ky = -dz'(0),
  // tz_i = int (1 - xi) kz = dy'(0)), hence dz'' = ky, dy'' = -kz.
  void init_bed(const Vec3& k) {
    kbed = k;
    const std::size_t np = xg.size();
    const int n3 = static_cast<int>(3 * np);
    Sop = Eigen::MatrixXd::Zero(n3, n3);
    Rop = Eigen::MatrixXd::Zero(12, n3);
    Dop = Eigen::MatrixXd::Zero(n3, n3);
    Cop = Eigen::MatrixXd::Zero(n3, 12);
    Pbed = Eigen::VectorXd::Zero(n3);
    sp_bed.assign(np, Vec3::Zero());
    if (!has_bed()) return;
    std::vector<std::vector<std::size_t>> pts(subs.size());
    for (std::size_t j = 0; j < np; ++j) pts[static_cast<std::size_t>(sub_of[j])].push_back(j);
    auto lag = [&](std::size_t j, double t) {  // Lagrange basis of point j on its sub-interval
      double v = 1;
      for (std::size_t m : pts[static_cast<std::size_t>(sub_of[j])])
        if (m != j) v *= (t - xg[m]) / (xg[j] - xg[m]);
      return v;
    };
    std::vector<double> gx, gw;
    gauss_legendre01(static_cast<int>(pts.front().size()) + 2, gx, gw);
    auto integ = [&](std::size_t j, double a, double c, auto&& f) {  // int_a^c f l_j
      double acc = 0;
      if (c <= a) return acc;
      for (std::size_t g = 0; g < gx.size(); ++g) {
        const double t = a + gx[g] * (c - a);
        acc += gw[g] * (c - a) * f(t) * lag(j, t);
      }
      return acc;
    };
    for (std::size_t j = 0; j < np; ++j) {
      const double sa = subs[static_cast<std::size_t>(sub_of[j])].first, sb = subs[static_cast<std::size_t>(sub_of[j])].second;
      for (int c = 0; c < 3; ++c) {
        auto loads_upto = [&](double x, double x0, Vec3& Fq, Vec3& Mq, bool = false) {
          const double cc = std::min(sb, x);
          const double F0 = integ(j, sa, cc, [](double) { return 1.0; });
          const double F1 = integ(j, sa, cc, [&](double t) { return t - x0; });
          Vec3 e = Vec3::Zero();
          e(c) = 1;
          Fq = F0 * e;
          Mq = F1 * Vec3(0, -e.z(), e.y());  // int (s - x0) e_x x e l_j
        };
        SpanParticular spp = span_particular_core(L, xg, &incl, loads_upto);
        release_statics(spp, rel_local, L, xg);
        for (std::size_t k = 0; k < np; ++k)
          Sop.block<3, 1>(static_cast<Eigen::Index>(3 * k), static_cast<Eigen::Index>(3 * j + c)) = spp.sp[k];
        Rop.col(static_cast<Eigen::Index>(3 * j + c)) = spp.p_ss;
      }
    }
    for (std::size_t i = 0; i < np; ++i) {
      const double x = xg[i], xi = x / L;
      for (std::size_t j = 0; j < np; ++j) {
        const double sa = subs[static_cast<std::size_t>(sub_of[j])].first, sb = subs[static_cast<std::size_t>(sub_of[j])].second;
        const double c = std::min(sb, x);
        // delta = int_0^x (x - t) kappa - x/L int_0^L (L - t) kappa
        const double A = integ(j, sa, c, [&](double t) { return x - t; }) - xi * integ(j, sa, sb, [&](double t) { return L - t; });
        // axial: int_0^x eps - x/L int_0^L eps
        const double Ax = integ(j, sa, c, [](double) { return 1.0; }) - xi * integ(j, sa, sb, [](double) { return 1.0; });
        Dop(static_cast<Eigen::Index>(3 * i), static_cast<Eigen::Index>(3 * j)) = Ax;
        Dop(static_cast<Eigen::Index>(3 * i + 1), static_cast<Eigen::Index>(3 * j + 2)) = -A;  // dy'' = -kz
        Dop(static_cast<Eigen::Index>(3 * i + 2), static_cast<Eigen::Index>(3 * j + 1)) = A;   // dz'' = ky
      }
      for (int c = 0; c < 3; ++c) {
        Cop(static_cast<Eigen::Index>(3 * i + c), c) = 1 - xi;
        Cop(static_cast<Eigen::Index>(3 * i + c), 6 + c) = xi;
      }
    }
  }

  // The coupled element equations with the bed, for a target basic
  // deformation v_in and the bed law r2 = P + g0 + M (D E) at the points:
  // M the bed stiffness in the element's current axes (diag(k) in linear
  // kinematics; E^T R0^T diag(k) R0 E in large displacements), g0 the
  // pressure of the chord displacement field, D E the deflection relative to
  // the chord from the total strains. Solved by Newton to round-off; the
  // converged Jacobian (bed_J), the deflection D E (bed_delta) and the
  // operators for sensitivities stay for the caller's tangent.
  Eigen::MatrixXd bed_J, bed_Wb, bed_Fs, bed_MD;
  Eigen::VectorXd bed_delta;
  Eigen::VectorXd bed_pnoise;  // the bed pressures' round-off (64 eps of the bed law's magnitudes)
  void solve_bed(const Vec6& v_in, const Mat3& Mloc, const Eigen::VectorXd& g0) {
    const Vec6 v = v_in - vshift;
    qT = torsion_stiffness() * v(5);
    const int nf = static_cast<int>(Cq.cols());
    q = Pq * q;
    const std::size_t np = xg.size();
    const int n3 = static_cast<int>(3 * np);
    Eigen::MatrixXd Mblk = Eigen::MatrixXd::Zero(n3, n3);
    for (std::size_t k = 0; k < np; ++k) Mblk.block<3, 3>(static_cast<Eigen::Index>(3 * k), static_cast<Eigen::Index>(3 * k)) = Mloc;
    const Eigen::MatrixXd Mabs = Mblk.cwiseAbs();
    const Vec5 vt = Pq * v.head<5>();
    const double e_v = vt.dot(kb_ref * vt);
    Eigen::MatrixXd& J = bed_J;
    J = Eigen::MatrixXd::Zero(nf + n3, nf + n3);
    Eigen::VectorXd r(nf + n3);
    Eigen::MatrixXd Bst = Eigen::MatrixXd::Zero(n3, 5);
    Eigen::MatrixXd& Fs = bed_Fs;
    Fs = Eigen::MatrixXd::Zero(n3, n3);
    Eigen::MatrixXd& Wb = bed_Wb;
    Wb = Eigen::MatrixXd::Zero(5, n3);
    bed_MD = Mblk * Dop;
    const Eigen::MatrixXd& MD = bed_MD;
    Mat5 fb;
    Vec5 vabs = Vec5::Zero();  // sum w |b^T| |e|: the round-off scale of v_cur
    Eigen::VectorXd r2abs;     // the round-off scale of the bed residual
    const Eigen::MatrixXd CqT = Cq.transpose();
    auto evaluate = [&]() -> bool {
      for (std::size_t k = 0; k < np; ++k) sp_bed[k] = Sop.middleRows(static_cast<Eigen::Index>(3 * k), 3) * Pbed;
      Vec5 v_cur;
      if (!eval_stations(v_cur, fb)) return false;
      Eigen::VectorXd E(n3);  // total strains (imposed included): the deflection
      for (std::size_t k = 0; k < np; ++k) {
        const Eigen::Index o = static_cast<Eigen::Index>(3 * k);
        E.segment<3>(o) = d[k] + e0[k];
        Fs.block<3, 3>(o, o) = fsk[k];
        Bst.middleRows(o, 3) = b_of(xg[k]);
        Wb.middleCols(o, 3) = wg[k] * b_of(xg[k]).transpose();
      }
      vabs = Wb.cwiseAbs() * E.cwiseAbs();
      const Eigen::MatrixXd FsS = Fs * Sop, FsB = Fs * Bst;
      const Eigen::MatrixXd G = Wb * FsS;  // d v_cur / dP
      r.head(nf) = CqT * (vt - v_cur);
      J.topLeftCorner(nf, nf) = -CqT * fb * Cq;
      J.block(0, nf, nf, n3) = -CqT * G;
      J.block(nf, 0, n3, nf) = MD * FsB * Cq;
      J.block(nf, nf, n3, n3) = Eigen::MatrixXd::Identity(n3, n3) + MD * FsS;
      bed_delta = Dop * E;
      r.tail(n3) = Pbed + g0 + Mblk * bed_delta;
      // (+ the section solves' strain noise through the deflection)
      Eigen::VectorXd en(n3);
      for (std::size_t k = 0; k < np; ++k) en.segment<3>(static_cast<Eigen::Index>(3 * k)) = strain_noise(k);
      r2abs = Pbed.cwiseAbs() + g0.cwiseAbs() + Mabs * (Dop.cwiseAbs() * (E.cwiseAbs() + en / (64 * std::numeric_limits<double>::epsilon())));
      return true;
    };
    // residual measures above their round-off: compatibility relative to
    // the target deformation (energy norm), the bed law relative to the
    // pressure scale; each component's noise (64 eps of the magnitudes it is
    // summed from, and the section solves' floor) is discounted, so a
    // residual at round-off measures zero (section strains are solved to
    // 1e-13: compatibility is not asked finer than 1e-12)
    const double tol1 = std::max(tol, 1e-12);
    auto measures = [&](double& m1, double& m2) {
      const double eps64 = 64 * std::numeric_limits<double>::epsilon();
      const Vec5 rf = Cq * r.head(nf);
      const Vec5 vsec = section_noise();
      Vec5 r1 = Vec5::Zero();
      for (int c = 0; c < 5; ++c) r1(c) = std::max(0.0, std::fabs(rf(c)) - eps64 * (vabs(c) + std::fabs(vt(c))) - vsec(c));
      const double e_s = std::max({e_v, vabs.dot(kb_ref.diagonal().cwiseProduct(vabs)), 1e-300});
      m1 = std::sqrt(r1.dot(kb_ref * r1) / e_s);
      const double sc2 = std::max({g0.lpNorm<Eigen::Infinity>(), Pbed.lpNorm<Eigen::Infinity>(), 1e-300});
      m2 = 0;
      for (int i = 0; i < n3; ++i) m2 = std::max(m2, std::max(0.0, std::fabs(r(nf + i)) - eps64 * r2abs(i)) / sc2);
    };
    if (!evaluate()) throw NonConvergence("fiber bar on a bed: section infeasible at the warm start");
    for (int it = 0;; ++it) {
      double m1, m2;
      measures(m1, m2);
      if (m1 < tol1 && m2 < tol) break;
      if (it >= 60) throw NonConvergence("fiber bar on a bed: element state determination did not converge");
      const Eigen::FullPivLU<Eigen::MatrixXd> lu(J);
      if (!lu.isInvertible()) throw NonConvergence("fiber bar on a bed: singular element system");
      const Eigen::VectorXd dz = -lu.solve(r);
      const double merit0 = std::max(m1 / tol1, m2 / tol);
      const Vec5 q0 = q;
      const Eigen::VectorXd P0 = Pbed;
      double a = 1;
      bool ok = false;
      for (int back = 0; back < 30 && !ok; ++back, a *= 0.5) {
        q = q0 + a * (Cq * dz.head(nf));
        Pbed = P0 + a * dz.tail(n3);
        if (!evaluate()) continue;
        double t1, t2;
        measures(t1, t2);
        ok = std::max(t1 / tol1, t2 / tol) < merit0;
      }
      if (!ok) throw NonConvergence("fiber bar on a bed: element state determination failed");
    }
    kb.setZero();
    kb.block<5, 5>(0, 0) = free_inverse(fb);
    kb(5, 5) = torsion_stiffness();
    bed_pnoise = 64 * std::numeric_limits<double>::epsilon() * r2abs;
  }
  // Sensitivities at the converged state: dz = (dQ, dP) for perturbations of
  // the residual sources, columns of rhs = d(r)/d(source) (dz = -J^-1 rhs)
  Eigen::MatrixXd bed_sensitivity(const Eigen::MatrixXd& rhs) const {
    return -Eigen::FullPivLU<Eigen::MatrixXd>(bed_J).solve(rhs);
  }
  // d(r)/d(sp) for perturbations dsp (3np x c) of the particular field
  Eigen::MatrixXd bed_dr_dsp(const Eigen::MatrixXd& dsp) const {
    const int nf = static_cast<int>(Cq.cols());
    Eigen::MatrixXd out(nf + dsp.rows(), dsp.cols());
    out.topRows(nf) = -Cq.transpose() * (bed_Wb * (bed_Fs * dsp));
    out.bottomRows(dsp.rows()) = bed_MD * (bed_Fs * dsp);
    return out;
  }

  // State determination with the bed, LINEAR kinematics: local end
  // displacements ul (12) -> local end forces f (basic forces + bed
  // reactions; the applied span loads' p_ss are the caller's) and the exact
  // tangent K = df/dul.
  void determine_bed(const Vec12& ul, Vec12& f, Mat12& K) {
    const Mat612 a6 = basic_transform6(L);
    const std::size_t np = xg.size();
    const int n3 = static_cast<int>(3 * np);
    const int nf = static_cast<int>(Cq.cols());
    Eigen::VectorXd kf(n3);
    for (std::size_t k = 0; k < np; ++k) kf.segment<3>(static_cast<Eigen::Index>(3 * k)) = kbed;
    solve_bed(a6 * ul, kbed.asDiagonal().toDenseMatrix(), kf.cwiseProduct(Cop * ul));
    // tangent: J dz = -(dr/dul) dul, dr1/dul = Cq^T a6, dr2/dul = kf Cop
    Eigen::MatrixXd rhs(nf + n3, 12);
    rhs.topRows(nf) = Cq.transpose() * a6.topRows<5>();
    rhs.bottomRows(n3) = kf.asDiagonal() * Cop;
    const Eigen::MatrixXd dz = bed_sensitivity(rhs);
    const Eigen::Matrix<double, 5, 12> dq = Cq * dz.topRows(nf);
    K = a6.topRows<5>().transpose() * dq + Rop * dz.bottomRows(n3);
    K += torsion_stiffness() * a6.row(5).transpose() * a6.row(5);
    f = a6.transpose() * q6() + Rop * Pbed;
  }

  Vec6 q6() const {
    Vec6 o;
    o << q(0), q(1), q(2), q(3), q(4), qT;
    return o;
  }

  // Roll the committed material states forward to the last evaluated strains.
  void commit_states() {
    for (std::size_t k = 0; k < xg.size(); ++k) secs[k]->commit(d[k], committed[k]);
  }
};

// ---------------- kinematics ----------------

// Linear basic transform a6 (6x12) on LOCAL end displacements
// (uxi uyi uzi rxi ryi rzi | uxj uyj uzj rxj ryj rzj) — see the header for
// the sign conventions (ry = -w', rz = +v').
inline Mat612 basic_transform6(double L) {
  Mat612 a = Mat612::Zero();
  a(0, 0) = -1;
  a(0, 6) = 1;
  a(1, 2) = -1.0 / L;
  a(1, 4) = 1.0;
  a(1, 8) = 1.0 / L;
  a(2, 2) = -1.0 / L;
  a(2, 8) = 1.0 / L;
  a(2, 10) = 1.0;
  a(3, 1) = 1.0 / L;
  a(3, 5) = 1.0;
  a(3, 7) = -1.0 / L;
  a(4, 1) = 1.0 / L;
  a(4, 7) = -1.0 / L;
  a(4, 11) = 1.0;
  a(5, 3) = -1.0;
  a(5, 9) = 1.0;
  return a;
}

// Exact basic stiffness of an elastic bar from its condensed local matrix:
// kb = (a a^t)^-1 a k_local a^t (a a^t)^-1 (k_local = a^t kb a exactly, since
// null(k_local) = rigid modes = null(a)). Carries Timoshenko shear and end
// releases through to the corotational mode.
inline Mat6 elastic_basic6(const Mat12& k_local, double L) {
  Mat612 a = basic_transform6(L);
  Mat6 M = a * a.transpose();
  Eigen::FullPivLU<Mat6> lu(M);
  require(lu.isInvertible(), "fiber_frame: degenerate basic transform");
  Mat6 Mi = lu.inverse();
  return Mi * a * k_local * a.transpose() * Mi;
}


// ---------------- span loads on a force-based member ----------------
//
// The basic system (simply supported in both planes, axial and torsion held
// at end i) under the bar's span loads: its reactions p_ss (local end forces
// ON the member) and, at abscissae x, the PARTICULAR resultants (N, My, Mz)
// of the section — so that s(x) = b(x) q + lam s_p(x) is the exact internal
// force field and the member end forces are a6^T q + lam p_ss. Linear
// distributed forces / moments and point loads, closed-form statics (the
// moment of a linear load is integrated exactly by Simpson's rule). The
// section convention is calibrated once against b(x) q (signs of N, My, Mz
// relative to the cut-face equilibrium) and checked consistent. Torsional
// span moments are a typed rejection (the element's torsion is one constant
// basic force).
// The basic-system statics of span loads given by their resultants:
// loads_upto(x, x0, Fq, Mq, at_x) returns the resultant force Fq and moment
// Mq about x0 of the loads on [0, x] (point loads AT x included iff at_x).
template <class LoadsUpto>
SpanParticular span_particular_core(double L, const std::vector<double>& xs, const std::vector<char>* incl,
                                    LoadsUpto&& loads_upto, double torque_lm) {
  auto ex_cross = [](const Vec3& F) { return Vec3(0, -F.z(), F.y()); };  // e_x x F
  SpanParticular out;
  Vec3 Q, Mo;
  loads_upto(L * (1 + 1e-15) + 1e-300, 0.0, Q, Mo);
  const Vec3 Rj(0, -Mo.z() / L, Mo.y() / L);
  const Vec3 Ri = -Q - Rj;
  out.p_ss.segment<3>(0) = Ri;
  out.p_ss.segment<3>(6) = Rj;
  // torsion held at end i: its reaction -(total torque); the particular
  // torque T_p(x) = -(Mi_x + torque of the loads on [0, x]) (cut-face
  // convention of q_T: a6^T e_T puts -1 at rx_i), integrated over the bar:
  // int_0^L T_p = -(L Mi_x + int (L - s) m_x ds + sum (L - a) M_x)
  out.p_ss(3) = -Mo.x();
  {
    const Mat612 a6t = basic_transform6(L);
    const double sgT = -a6t(5, 3);  // T = sgT * (-(Mi_x + ...)) with Mi_x = a6(5,3) for q_T = 1
    out.tor_int = sgT * -(L * out.p_ss(3) + torque_lm);
  }
  // raw cut resultants (R_x, Mr_y, Mr_z) at x for end-i actions (F_i, M_i)
  auto raw = [&](double x, const Vec3& Fi, const Vec3& Mi, bool with_loads, bool at_x = false) {
    Vec3 Fq = Vec3::Zero(), Mq = Vec3::Zero();
    if (with_loads) loads_upto(x, x, Fq, Mq, at_x);
    const Vec3 R = -(Fi + Fq);
    const Vec3 Mr = -(Mi + (0.0 - x) * ex_cross(Fi) + Mq);
    return Vec3(R.x(), Mr.y(), Mr.z());
  };
  // calibration against b(x) q on the five fiber modes
  const Mat612 a6 = basic_transform6(L);
  Vec3 sg = Vec3::Zero();
  const double xc = 0.37 * L;
  for (int k = 0; k < 5; ++k) {
    const Vec12 f = a6.transpose().col(k);
    const Vec3 r = raw(xc, f.segment<3>(0), f.segment<3>(3), false);
    Vec5 ek = Vec5::Zero();
    ek(k) = 1;
    const double xi = xc / L;
    const Vec3 bq((k == 0) ? 1.0 : 0.0, k == 1 ? 1 - xi : (k == 2 ? -xi : 0.0), k == 3 ? 1 - xi : (k == 4 ? -xi : 0.0));
    for (int c = 0; c < 3; ++c)
      if (std::fabs(bq(c)) > 1e-12) {
        const double sgn = r(c) / bq(c);
        require(std::fabs(std::fabs(sgn) - 1) < 1e-9, "fiber bar: span-load statics calibration failed");
        if (sg(c) == 0) sg(c) = sgn > 0 ? 1.0 : -1.0;
        require(sg(c) * sgn > 0, "fiber bar: inconsistent span-load sign convention");
      }
  }
  for (std::size_t k = 0; k < xs.size(); ++k) {
    const bool at_x = incl != nullptr && k < incl->size() && (*incl)[k] != 0;
    const Vec3 r = raw(xs[k], Ri, Vec3::Zero(), true, at_x);
    out.sp.push_back(r.cwiseProduct(sg));
  }
  return out;
}

// incl (optional, per abscissa): the right limit at a point-load abscissa
// (loads AT x included); otherwise the left limit (loads strictly before x).
inline SpanParticular span_particular(const BarSpanLoads& sl, double L, const std::vector<double>& xs,
                                      const std::vector<char>* incl = nullptr) {
  // (imposed strains / dislocations carry no load: the caller's — see
  // imposed_fields)
  double torque_lm = 0;  // int (L - s) m_x ds + sum (L - a) M_x
  require_linear_span(sl.segs, "fibre element");
  for (const DistSeg& g : sl.segs) {
    if (g.b <= g.a) continue;
    const double h = g.b - g.a, mid = 0.5 * (g.a + g.b);
    torque_lm += h / 6 * ((L - g.a) * g.ma.x() + 4 * (L - mid) * 0.5 * (g.ma.x() + g.mb.x()) + (L - g.b) * g.mb.x());
  }
  for (const PointAt& pt : sl.points) torque_lm += (L - pt.a) * pt.M.x();
  auto ex_cross = [](const Vec3& F) { return Vec3(0, -F.z(), F.y()); };  // e_x x F
  // resultant force and moment about the point x0 of the loads on [0, x]
  auto loads_upto = [&](double x, double x0, Vec3& Fq, Vec3& Mq, bool at_x = false) {
    Fq.setZero();
    Mq.setZero();
    for (const DistSeg& g : sl.segs) {
      const double a = g.a, b = std::min(g.b, x);
      if (b <= a) continue;
      auto q_at = [&](double t) { return Vec3(g.qa + (g.qb - g.qa) * ((t - g.a) / (g.b - g.a))); };
      auto m_at = [&](double t) { return Vec3(g.ma + (g.mb - g.ma) * ((t - g.a) / (g.b - g.a))); };
      const double mid = 0.5 * (a + b), h = b - a;
      const Vec3 qa = q_at(a), qm = q_at(mid), qb = q_at(b);
      Fq += h / 6 * (qa + 4 * qm + qb);
      // integral of (s - x0) e_x x q(s): quadratic integrand, Simpson exact
      Mq += h / 6 * ((a - x0) * ex_cross(qa) + 4 * (mid - x0) * ex_cross(qm) + (b - x0) * ex_cross(qb));
      Mq += h / 6 * (m_at(a) + 4 * m_at(mid) + m_at(b));
    }
    for (const PointAt& pt : sl.points) {
      if (at_x ? pt.a > x : pt.a >= x) continue;
      Fq += pt.F;
      Mq += (pt.a - x0) * ex_cross(pt.F) + pt.M;
    }
  };
  return span_particular_core(L, xs, incl, loads_upto, torque_lm);
}

// The span-load statics on a member with END RELEASES (local DOF flags):
// the basic system's reactions moved so that every released end action
// vanishes — a fixed basic-force shift q_c, b(x) q_c added to the particular
// field, a6^T q_c to the reactions, L q_c(T) to the torsion integral. A load
// the released member cannot carry (a mechanism) is a typed rejection.
inline void release_statics(SpanParticular& sp, const std::array<bool, 12>& rl, double L,
                            const std::vector<double>& xs) {
  std::vector<int> R;
  for (int i = 0; i < 12; ++i)
    if (rl[static_cast<std::size_t>(i)]) R.push_back(i);
  if (R.empty()) return;
  const Mat612 a6 = basic_transform6(L);
  Eigen::MatrixXd M(static_cast<Eigen::Index>(R.size()), 6);
  Eigen::VectorXd rhs(static_cast<Eigen::Index>(R.size()));
  for (std::size_t i = 0; i < R.size(); ++i) {
    M.row(static_cast<Eigen::Index>(i)) = a6.col(R[i]).transpose();
    rhs(static_cast<Eigen::Index>(i)) = -sp.p_ss(R[i]);
  }
  const Eigen::VectorXd qc = M.completeOrthogonalDecomposition().solve(rhs);
  require((M * qc - rhs).norm() <= 1e-9 * (rhs.norm() + sp.p_ss.norm()) + 1e-300,
          "fiber bar: the released member cannot carry its span load (mechanism)");
  const Vec6 q6 = qc;
  sp.p_ss += a6.transpose() * q6;
  for (std::size_t k = 0; k < xs.size() && k < sp.sp.size(); ++k) sp.sp[k] += b_basic(xs[k], L) * q6.head<5>();
  sp.tor_int += L * q6(5);
}

// Imposed strains of a span record on a force-based member: the strain
// field (eps0, kappa_y, kappa_z) at the points (thermal, prestrain, linear
// variations — keta's curvature sign is the fiber section's: theta jumps by
// -kappa) and the basic deformation of the concentrated dislocations,
// v_k = b(a)^T (dux, theta_y, theta_z) plus the chord rotations of the
// transverse jumps (tz = -dvy / L, ty = +dvz / L at both ends) and the twist
// jump.
inline void imposed_fields(const BarSpanLoads& sl, double L, const std::vector<double>& xs, std::vector<Vec3>& e0,
                           Vec6& vkink) {
  e0.assign(xs.size(), Vec3::Zero());
  for (std::size_t k = 0; k < xs.size(); ++k)
    e0[k] = Vec3(sl.eps0 + (sl.deps0 + sl.c2eps0 * xs[k]) * xs[k], sl.kappa_y + (sl.dkappa_y + sl.c2kappa_y * xs[k]) * xs[k],
                 sl.kappa_z + (sl.dkappa_z + sl.c2kappa_z * xs[k]) * xs[k]);
  vkink.setZero();
  for (const KinkAt& kk : sl.kinks) {
    vkink.head<5>() += b_basic(kk.a, L).transpose() * Vec3(kk.dux, kk.theta_y, kk.theta_z);
    vkink(3) += -kk.dvy / L;
    vkink(4) += -kk.dvy / L;
    vkink(1) += kk.dvz / L;
    vkink(2) += kk.dvz / L;
    vkink(5) += kk.dthx;
  }
}

}  // namespace keta::fiberel

#include "keta/wind_synthesis.h"

#include <cmath>

namespace keta {

using C = std::complex<double>;

void fft(std::vector<C>& a, int sign) {
  const std::size_t n = a.size();
  require(n >= 1 && (n & (n - 1)) == 0, "fft: size must be a power of two");
  for (std::size_t i = 1, j = 0; i < n; ++i) {
    std::size_t bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (std::size_t len = 2; len <= n; len <<= 1) {
    const double ang = sign * 2 * kPi / static_cast<double>(len);
    const C wl(std::cos(ang), std::sin(ang));
    for (std::size_t i = 0; i < n; i += len) {
      C w(1, 0);
      for (std::size_t k = 0; k < len / 2; ++k) {
        const C u = a[i + k], v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        w *= wl;
      }
    }
  }
}

WindSeries synthesize_wind(const TurbulenceModel& m, const std::vector<WindPoint>& pts, double n_max, int n_freq,
                           std::uint64_t seed) {
  const std::size_t N = pts.size();
  require(N >= 1 && n_max > 0 && n_freq >= 2, "wind synthesis: points, n_max > 0, n_freq >= 2");
  const double dn = n_max / n_freq;
  // bins: n_ml = ((l - 1) N + m) dn / N, l = 1..n_freq, m = 1..N; nt > 2 x the top bin
  const std::size_t top = static_cast<std::size_t>(n_freq) * N;
  std::size_t nt = 1;
  while (nt < 4 * top) nt <<= 1;
  WindSeries out;
  out.nt = static_cast<int>(nt);
  out.dt = 1.0 / (dn / static_cast<double>(N)) / static_cast<double>(nt);  // T / nt, T = N / dn
  out.series.assign(N, {std::vector<double>(nt), std::vector<double>(nt)});
  std::uint64_t st = seed ? seed : 0x9e3779b97f4a7c15ULL;
  auto uni = [&]() {  // splitmix64: platform-independent
    st += 0x9e3779b97f4a7c15ULL;
    std::uint64_t z = st;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z ^= z >> 31;
    return static_cast<double>(z >> 11) * 0x1.0p-53;
  };
  for (int c = 0; c < 2; ++c) {
    const int comp = c == 0 ? 0 : 2;
    std::vector<std::vector<C>> X(N, std::vector<C>(nt, C(0, 0)));
    Eigen::MatrixXd cov = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
    for (int l = 1; l <= n_freq; ++l)
      for (std::size_t mm = 1; mm <= N; ++mm) {
        const std::size_t bin = static_cast<std::size_t>(l - 1) * N + mm;
        const double n = static_cast<double>(bin) * dn / static_cast<double>(N);
        // the component's cross-spectral matrix at n_ml and its Cholesky factor
        Eigen::MatrixXd S(static_cast<Eigen::Index>(N), static_cast<Eigen::Index>(N));
        for (std::size_t j = 0; j < N; ++j)
          for (std::size_t k = j; k < N; ++k) {
            const WindPoint &a = pts[j], &b = pts[k];
            const double sa = a.sigma[comp] > 0 ? std::sqrt(auto_spectrum(m.form, comp, n, a.U, a.sigma[comp], a.L[comp])) : 0;
            const double sb = b.sigma[comp] > 0 ? std::sqrt(auto_spectrum(m.form, comp, n, b.U, b.sigma[comp], b.L[comp])) : 0;
            const double coh = j == k ? 1.0 : davenport_coherence(n, a.y - b.y, a.z - b.z, a.U, b.U, m.C[comp][0], m.C[comp][1]);
            S(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) = S(static_cast<Eigen::Index>(k), static_cast<Eigen::Index>(j)) = sa * sb * coh;
          }
        const double jitter = 1e-14 * std::max(S.diagonal().maxCoeff(), 1e-300);
        const Eigen::LLT<Eigen::MatrixXd> llt(S + jitter * Eigen::MatrixXd::Identity(S.rows(), S.cols()));
        require(llt.info() == Eigen::Success, "wind synthesis: cross-spectral matrix not positive definite");
        const Eigen::MatrixXd H = llt.matrixL();
        const double phi = 2 * kPi * uni();
        const Eigen::Index col = static_cast<Eigen::Index>(mm - 1);
        for (std::size_t j = mm - 1; j < N; ++j) {
          const double h = H(static_cast<Eigen::Index>(j), col);  // real: theta = 0 (real coherence)
          X[j][bin] += std::sqrt(2 * dn) * h * C(std::cos(phi), std::sin(phi));
        }
        for (std::size_t j = mm - 1; j < N; ++j)
          for (std::size_t k = mm - 1; k < N; ++k)
            cov(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(k)) +=
                dn * H(static_cast<Eigen::Index>(j), col) * H(static_cast<Eigen::Index>(k), col);
      }
    for (std::size_t j = 0; j < N; ++j) {
      fft(X[j], +1);  // sum_bin X e^{+i 2 pi bin k / nt}
      for (std::size_t k = 0; k < nt; ++k) out.series[j][static_cast<std::size_t>(c)][k] = X[j][k].real();
    }
    out.cov[static_cast<std::size_t>(c)] = cov;
  }
  return out;
}

Eigen::MatrixXd modal_time_response(const ModalSystem& sys, const Eigen::MatrixXd& Q, double dt) {
  const Eigen::Index m = static_cast<Eigen::Index>(sys.omega.size());
  require(m >= 1 && Q.rows() == m && Q.cols() >= 2 && dt > 0, "modal time response: m rows of forces, dt > 0");
  Eigen::MatrixXd K = Eigen::MatrixXd::Zero(m, m), Cm = K;
  for (Eigen::Index i = 0; i < m; ++i) {
    const double w = sys.omega[static_cast<std::size_t>(i)];
    K(i, i) = w * w;
    Cm(i, i) = 2 * sys.zeta[static_cast<std::size_t>(i)] * w;
  }
  if (sys.K_ae.size()) K += sys.K_ae;
  if (sys.C_ae.size()) Cm += sys.C_ae;
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
  // average acceleration: Keff = K + 2/dt C + 4/dt^2 M
  const Eigen::PartialPivLU<Eigen::MatrixXd> lu(K + 2 / dt * Cm + 4 / (dt * dt) * I);
  Eigen::MatrixXd q = Eigen::MatrixXd::Zero(m, Q.cols());
  Eigen::VectorXd u = Eigen::VectorXd::Zero(m), v = u, a = Q.col(0) - K * u - Cm * v;
  for (Eigen::Index k = 1; k < Q.cols(); ++k) {
    const Eigen::VectorXd rhs = Q.col(k) + (4 / (dt * dt)) * u + (4 / dt) * v + a + Cm * ((2 / dt) * u + v);
    const Eigen::VectorXd un = lu.solve(rhs);
    const Eigen::VectorXd vn = 2 / dt * (un - u) - v;
    const Eigen::VectorXd an = 4 / (dt * dt) * (un - u) - 4 / dt * v - a;
    u = un;
    v = vn;
    a = an;
    q.col(k) = u;
  }
  return q;
}

TemporalBuffeting temporal_buffeting(const Model& model, const ModalResults& modal, const BuffetingSpec& spec,
                                     const std::vector<BuffetingRequest>& requests, double n_max, int n_freq,
                                     std::uint64_t seed, int periods) {
  require(periods >= 2, "temporal buffeting: 2 periods at least (the first is discarded)");
  const std::vector<Strip> strips = buffeting_strips(model, modal, spec);
  const ModalSystem sys = buffeting_system(modal, strips, spec);
  std::vector<WindPoint> pts;
  for (const Strip& s : strips) pts.push_back(s.wp);
  const WindSeries w = synthesize_wind(spec.turb, pts, n_max, n_freq, seed);
  const Eigen::Index m = static_cast<Eigen::Index>(modal.omegas.size());
  const std::size_t nt = static_cast<std::size_t>(w.nt);
  // modal forces of one period: Q_i = sum_s L_s phi_i . A_s (u, w)
  Eigen::MatrixXd Q1 = Eigen::MatrixXd::Zero(m, static_cast<Eigen::Index>(nt));
  const AeroTable::Value c = spec.aero.at(spec.alpha0);
  for (std::size_t s = 0; s < strips.size(); ++s) {
    const Strip& st = strips[s];
    const double U = st.wp.U, q = 0.5 * spec.turb.profile.rho * U * U, B = spec.B;
    Eigen::Matrix<double, 3, 2> A;
    A << q * B * 2 * c.cd / U, q * B * (c.dcd - c.cl) / U, q * B * 2 * c.cl / U, q * B * (c.dcl + c.cd) / U,
        q * B * B * 2 * c.cm / U, q * B * B * c.dcm / U;
    for (Eigen::Index i = 0; i < m; ++i) {
      const Vec3& ph = st.phi[static_cast<std::size_t>(i)];
      const double gu = st.L * ph.dot(A.col(0)), gw = st.L * ph.dot(A.col(1));
      for (std::size_t t = 0; t < nt; ++t)
        Q1(i, static_cast<Eigen::Index>(t)) += gu * w.series[s][0][t] + gw * w.series[s][1][t];
    }
  }
  Eigen::MatrixXd Q(m, static_cast<Eigen::Index>(nt) * periods);
  for (int p = 0; p < periods; ++p) Q.middleCols(static_cast<Eigen::Index>(nt) * p, static_cast<Eigen::Index>(nt)) = Q1;
  const Eigen::MatrixXd qt = modal_time_response(sys, Q, w.dt);
  const Assembled asm_ = assemble_elements(model);
  TemporalBuffeting out;
  out.dt = w.dt;
  for (const BuffetingRequest& r : requests) {
    require(r.comp >= 0 && r.comp < 6 && model.node_index(r.node) >= 0, "temporal buffeting: request node / component");
    const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(model.node_index(r.node), r.comp))];
    Eigen::VectorXd psi = Eigen::VectorXd::Zero(m);
    if (f >= 0)
      for (Eigen::Index i = 0; i < m; ++i) psi(i) = modal.modes_free[static_cast<std::size_t>(i)](f);
    TemporalResponse tr;
    const Eigen::Index k0 = static_cast<Eigen::Index>(nt);
    tr.series.reserve(static_cast<std::size_t>(qt.cols() - k0));
    double mean = 0;
    for (Eigen::Index k = k0; k < qt.cols(); ++k) {
      tr.series.push_back(psi.dot(qt.col(k)));
      mean += tr.series.back();
    }
    mean /= static_cast<double>(tr.series.size());
    double var = 0;
    tr.max = tr.min = tr.series.front();
    for (double v : tr.series) {
      var += (v - mean) * (v - mean);
      tr.max = std::max(tr.max, v);
      tr.min = std::min(tr.min, v);
    }
    tr.sigma = std::sqrt(var / static_cast<double>(tr.series.size()));
    out.responses.push_back(std::move(tr));
  }
  return out;
}

}  // namespace keta

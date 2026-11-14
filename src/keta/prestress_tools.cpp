#include "keta/prestress_tools.h"

#include <algorithm>
#include <cmath>

#include "keta/prestress.h"

namespace keta {

// ------------------------------------------------------------ lump-sum losses
double ip2_lump_shrinkage(double rh) { return rh < 65 ? 4e-4 : 3e-4; }
double bpel_lump_shrinkage(double rh) { return rh < 65 ? 3e-4 : 2e-4; }
double bpel_mu0(double rho1000) {
  if (rho1000 <= 2.5) return 0.43;
  if (rho1000 <= 8.0) return 0.35;
  return 0.30;
}

namespace {
void check_lump(double sigma_pi, double fprg, double Ep, double rh, const LumpSumOptions& o) {
  require(sigma_pi > 0 && fprg > 0 && Ep > 0, "lump-sum losses: sigma_pi, fprg and Ep > 0");
  require(rh >= 0 && rh <= 100, "lump-sum losses: relative humidity in [0, 100] %");
  require(o.n > 0, "lump-sum losses: n = Ea / Eb > 0");
}
}  // namespace

LumpSumLoss ip2_lump_sum(double sigma_b, double sigma_pi, double fprg, double Ep, double rh, double final_relaxation,
                         const LumpSumOptions& o) {
  check_lump(sigma_pi, fprg, Ep, rh, o);
  require(final_relaxation >= 0, "IP2 lump-sum losses: the final relaxation is >= 0");
  LumpSumLoss r;
  r.elastic = 0.5 * o.n * sigma_b;
  r.shrinkage = (o.eps_r >= 0 ? o.eps_r : ip2_lump_shrinkage(rh)) * Ep;
  r.creep = 2 * sigma_b * o.n;
  r.relaxation = final_relaxation;
  const double drop = r.shrinkage + r.creep, margin = sigma_pi - 0.55 * fprg;
  r.viscoelastic = margin > 0 && drop < margin ? drop + r.relaxation * (1 - drop / margin) : drop;
  return r;
}

LumpSumLoss bpel_lump_sum(double sigma_b, double sigma_pi, double fprg, double Ep, double rh, double rho1000,
                          const LumpSumOptions& o) {
  check_lump(sigma_pi, fprg, Ep, rh, o);
  require(rho1000 >= 0, "BPEL lump-sum losses: rho1000 >= 0");
  LumpSumLoss r;
  r.elastic = 0.5 * o.n * sigma_b;
  r.shrinkage = (o.eps_r >= 0 ? o.eps_r : bpel_lump_shrinkage(rh)) * Ep;
  r.creep = 2.5 * sigma_b * o.n;
  const double mu0 = o.mu0 >= 0 ? o.mu0 : bpel_mu0(rho1000);
  r.relaxation = 0.06 * rho1000 * std::max(sigma_pi / fprg - mu0, 0.0) * sigma_pi;
  r.viscoelastic = r.shrinkage + r.creep + 5.0 / 6.0 * r.relaxation;
  return r;
}

std::optional<StressPlane> plane_fit_stress(const std::vector<std::array<double, 3>>& points) {
  const std::size_t n = points.size();
  if (n < 3) return std::nullopt;
  double my = 0, mz = 0, ms = 0;
  for (const auto& p : points) {
    my += p[0];
    mz += p[1];
    ms += p[2];
  }
  my /= double(n);
  mz /= double(n);
  ms /= double(n);
  // centred normal equations [Syy Syz; Syz Szz] [b c] = [Sys Szs]
  double syy = 0, szz = 0, syz = 0, sys = 0, szs = 0;
  for (const auto& p : points) {
    const double y = p[0] - my, z = p[1] - mz, s = p[2] - ms;
    syy += y * y;
    szz += z * z;
    syz += y * z;
    sys += y * s;
    szs += z * s;
  }
  const double det = syy * szz - syz * syz;
  if (!(det > 1e-12 * syy * szz) || syy <= 0 || szz <= 0) return std::nullopt;
  StressPlane pl;
  pl.b = (sys * szz - szs * syz) / det;
  pl.c = (szs * syy - sys * syz) / det;
  pl.a = ms - pl.b * my - pl.c * mz;
  return pl;
}

double StressTable::at(double x) const {
  require(!points.empty(), "stress table: no point");
  if (x <= points.front().first) return points.front().second;
  if (x >= points.back().first) return points.back().second;
  const auto it = std::upper_bound(points.begin(), points.end(), x,
                                   [](double v, const std::pair<double, double>& p) { return v < p.first; });
  const auto& [x1, v1] = *it;
  const auto& [x0, v0] = *(it - 1);
  return v0 + (v1 - v0) * (x - x0) / (x1 - x0);
}

// ------------------------------------------------------------ trace transforms
namespace {
// coefficients of c(L - x) for the cubic c(x) = sum c_k x^k
std::array<double, 4> reflected(const std::array<double, 4>& c, double L) {
  static const double C[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
  std::array<double, 4> out{};
  for (int j = 0; j < 4; ++j) {
    double v = 0;
    for (int k = j; k < 4; ++k) v += c[static_cast<std::size_t>(k)] * C[k][j] * std::pow(L, k - j);
    out[static_cast<std::size_t>(j)] = (j % 2 ? -1.0 : 1.0) * v;
  }
  return out;
}
PolySegment translated_seg(PolySegment g, double dx) {
  g.s0 += dx;
  g.s1 += dx;
  if (g.exact == 1) g.p[0].x() += dx;
  if (g.exact == 2) g.C.x() += dx;
  return g;
}
PolySegment mirrored_seg(const PolySegment& g, double p) {
  PolySegment m = g;
  const double L = g.s1 - g.s0;
  m.s0 = 2 * p - g.s1;
  m.s1 = 2 * p - g.s0;
  m.cy = reflected(g.cy, L);
  m.cz = reflected(g.cz, L);
  const Vec3 flip(-1, 1, 1);
  if (g.exact == 1) {
    // q(u) = X(1 - u): q_j = (-1)^j sum_{k >= j} C(k, j) p_k; then s -> 2p - s
    static const double C[4][4] = {{1, 0, 0, 0}, {1, 1, 0, 0}, {1, 2, 1, 0}, {1, 3, 3, 1}};
    for (int j = 0; j < 4; ++j) {
      Vec3 q = Vec3::Zero();
      for (int k = j; k < 4; ++k) q += C[k][j] * g.p[static_cast<std::size_t>(k)];
      m.p[static_cast<std::size_t>(j)] = ((j % 2) ? -1.0 : 1.0) * q.cwiseProduct(flip);
    }
    m.p[0].x() += 2 * p;
  } else if (g.exact == 2) {
    // X'(phi) = M X(th - phi): U' = L (cos th U - sin th D), D' = -L (sin th U + cos th D)
    const double c = std::cos(g.th), s = std::sin(g.th);
    m.C = g.C.cwiseProduct(flip);
    m.C.x() += 2 * p;
    m.U = (c * g.U - s * g.D).cwiseProduct(flip);
    m.D = (-(s * g.U + c * g.D)).cwiseProduct(flip);
  }
  return m;
}
}  // namespace

std::vector<PolySegment> profile_translated(const std::vector<PolySegment>& profile, double dx) {
  std::vector<PolySegment> out;
  out.reserve(profile.size());
  for (const PolySegment& g : profile) out.push_back(translated_seg(g, dx));
  return out;
}

std::vector<PolySegment> profile_mirrored(const std::vector<PolySegment>& profile, double pivot) {
  std::vector<PolySegment> out;
  out.reserve(profile.size());
  for (auto it = profile.rbegin(); it != profile.rend(); ++it) out.push_back(mirrored_seg(*it, pivot));
  return out;
}

std::vector<PolySegment> profile_symmetric_completion(const std::vector<PolySegment>& half, double pivot) {
  require(!half.empty(), "symmetric completion: empty profile");
  const double span = std::max(std::abs(half.back().s1 - half.front().s0), 1.0);
  require(std::abs(half.back().s1 - pivot) <= 1e-9 * span, "symmetric completion: the half trace must end at the pivot");
  std::vector<PolySegment> out = half;
  out.back().s1 = pivot;
  const std::vector<PolySegment> m = profile_mirrored(half, pivot);
  out.insert(out.end(), m.begin(), m.end());
  out[half.size()].s0 = pivot;
  return out;
}

namespace {
double chain_end(const Model& model, const Tendon& t) { return t.s_end < 0 ? tendon_chain_length(model, t) : t.s_end; }
void check_extent(double a, double b, double L, const char* what) {
  const double tol = 1e-9 * std::max(L, 1.0);
  require(a >= -tol && b <= L + tol && b > a, std::string(what) + ": the transformed tendon leaves its chain");
}
}  // namespace

Tendon tendon_mirrored(const Model& model, const Tendon& t, double pivot) {
  const double L = tendon_chain_length(model, t), e = chain_end(model, t);
  Tendon m = t;
  m.profile = profile_mirrored(t.profile, pivot);
  m.s_begin = 2 * pivot - e;
  m.s_end = 2 * pivot - t.s_begin;
  check_extent(m.s_begin, m.s_end, L, "tendon mirror");
  if (std::abs(m.s_begin) <= 1e-12 * std::max(L, 1.0)) m.s_begin = 0;
  if (std::abs(m.s_end - L) <= 1e-12 * std::max(L, 1.0)) m.s_end = -1;
  m.force_profile.clear();
  for (auto it = t.force_profile.rbegin(); it != t.force_profile.rend(); ++it)
    m.force_profile.emplace_back(2 * pivot - it->first, it->second);
  m.free_runs.clear();
  for (auto it = t.free_runs.rbegin(); it != t.free_runs.rend(); ++it)
    m.free_runs.emplace_back(2 * pivot - it->second, 2 * pivot - it->first);
  if (!t.jack_both_ends) m.jack_from_end = !t.jack_from_end;
  else m.release_end_first = !t.release_end_first;
  return m;
}

Tendon tendon_translated(const Model& model, const Tendon& t, double dx) {
  const double L = tendon_chain_length(model, t), e = chain_end(model, t);
  Tendon m = t;
  m.profile = profile_translated(t.profile, dx);
  m.s_begin = t.s_begin + dx;
  m.s_end = e + dx;
  check_extent(m.s_begin, m.s_end, L, "tendon translation");
  if (std::abs(m.s_begin) <= 1e-12 * std::max(L, 1.0)) m.s_begin = 0;
  if (std::abs(m.s_end - L) <= 1e-12 * std::max(L, 1.0)) m.s_end = -1;
  for (auto& [s, P] : m.force_profile) s += dx;
  for (auto& [a, b] : m.free_runs) {
    a += dx;
    b += dx;
  }
  return m;
}

Tendon tendon_on_bars(const Model& model, const Tendon& t, const std::vector<Id>& bars, const std::vector<bool>& flips) {
  require(!bars.empty(), "tendon copy: no bar");
  require(flips.empty() || flips.size() == bars.size(), "tendon copy: one flip per bar");
  Tendon m = t;
  m.bars = bars;
  m.flips = flips;
  const double L0 = tendon_chain_length(model, t), L1 = tendon_chain_length(model, m);
  require(std::abs(L1 - L0) <= 1e-9 * std::max(L0, 1.0), "tendon copy: the new chain has another length");
  return m;
}

}  // namespace keta

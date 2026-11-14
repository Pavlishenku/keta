// Prestress tools:
// regulatory lump-sum deferred losses (IP2, BPEL), the final concrete stress
// field that drives them, and exact transforms of tendon traces.
//
// --- Lump-sum (forfaitaire) deferred losses ---------------------------------
// The French lump-sum rules (IP2, BPEL), from a USER-PREDICTED final
// concrete stress sigma_b at tendon level (compression positive), all
// stresses in Pa:
//   elastic loss            3 sigma_b  (= n/2 sigma_b, non-simultaneity +
//                                       elastic shortening, n = Ea/Eb = 6)
//   shrinkage               d_r  = eps_r Ep
//   creep                   d_fl = 2 sigma_b n (IP2) | 2.5 sigma_b n (BPEL)
//   relaxation              IP2: d_p, the final pure relaxation (given);
//                           BPEL: d_p = 6/100 rho1000 (sigma_pi/fprg - mu0)
//                           sigma_pi, clamped at 0
//   viscoelastic total      IP2: d_r + d_fl + d_p (1 - (d_r + d_fl) /
//                           (sigma_pi - 0.55 fprg)) if d_r + d_fl < sigma_pi -
//                           0.55 fprg, else d_r + d_fl (the relaxation is
//                           eaten by the stress drop);
//                           BPEL: d_r + d_fl + 5/6 d_p.
// Shrinkage strain eps_r: the rule as printed gives IP2 4e-4 for H >= 65 % and
// 3e-4 below, BPEL 3e-4 for H >= 65 % and 2e-4 below — i.e. MORE shrinkage
// in MORE humid air, which contradicts the IP2 shrinkage law
// (eps_c(H) = (-8.68e-3 H^2 + 0.545 H + 32.4) 1e-5, decreasing for H > 31 %)
// and the physics. keta takes the physically consistent reading (the larger
// value in the drier climate: IP2 4e-4 below 65 %, 3e-4 above; BPEL 3e-4 /
// 2e-4) and lets the caller impose eps_r (LumpSumOptions::eps_r) to
// reproduce a computation made with the printed rule.
// mu0 (BPEL) by rho1000 band: 0.43 for rho1000 <= 2.5 %,
// 0.35 for 2.5 < rho1000 <= 8 %, 0.30 above; LumpSumOptions::mu0 overrides it.
//
// Final stress field: uniform, tabulated sigma_b(x) (linear, clamped,
// O(log n) by bisection), or the least-squares plane sigma = a + b y + c z
// through section stress points (centred normal equations, O(n); a
// degenerate — collinear — cloud returns no plane).
//
// --- Trace transforms (exact, O(segments)) ----------------------------------
// translated: e_t(s) = e(s - dx); mirrored about the chain abscissa p:
// e_m(s) = e(2p - s). Every segment kind is transformed exactly — the cubic
// in s is re-expanded (binomial), the parametric cubic (exact = 1) is
// reparametrised u -> 1 - u, the circular arc (exact = 2) has its centre
// reflected and its frame (U, D) rotated so the arc runs the other way.
// symmetric completion: a half trace ending at p followed by its mirror;
// the result is continuous at p and C1 there iff e'(p) = 0 (a non-zero slope
// leaves a kink 2 e'(p), which is kept: a deviator can sit at mid-span).
// Tendon-level copies carry the anchors, the authored force profile, the free
// runs and the active anchor along.
#pragma once
#include <array>
#include <optional>
#include <utility>
#include <vector>

#include "keta/model.h"

namespace keta {

struct LumpSumLoss {
  double elastic = 0, shrinkage = 0, creep = 0, relaxation = 0, viscoelastic = 0;  // Pa
  double total() const { return elastic + viscoelastic; }
};
struct LumpSumOptions {
  double n = 6;       // Ea / Eb, the regulatory 6
  double eps_r = -1;  // < 0: the climate value (see above)
  double mu0 = -1;    // BPEL; < 0: by rho1000 band
};
double ip2_lump_shrinkage(double rh);   // eps_r for the relative humidity rh [%]
double bpel_lump_shrinkage(double rh);
double bpel_mu0(double rho1000);        // rho1000 in %
// sigma_b: final concrete stress at tendon level (compression > 0) [Pa];
// sigma_pi: initial tendon stress; fprg: guaranteed ultimate strength (Rg).
LumpSumLoss ip2_lump_sum(double sigma_b, double sigma_pi, double fprg, double Ep, double rh, double final_relaxation,
                         const LumpSumOptions& o = {});
LumpSumLoss bpel_lump_sum(double sigma_b, double sigma_pi, double fprg, double Ep, double rh, double rho1000,
                          const LumpSumOptions& o = {});

struct StressPlane {
  double a = 0, b = 0, c = 0;
  double at(double y, double z) const { return a + b * y + c * z; }
};
// least-squares plane through (y, z, sigma) points; nullopt when fewer than 3
// points or collinear (relative tolerance 1e-12 on the centred moments)
std::optional<StressPlane> plane_fit_stress(const std::vector<std::array<double, 3>>& points);
// sigma_b(x) tabulated: strictly increasing x, linear interpolation, clamped
struct StressTable {
  std::vector<std::pair<double, double>> points;
  double at(double x) const;
};

std::vector<PolySegment> profile_translated(const std::vector<PolySegment>& profile, double dx);
std::vector<PolySegment> profile_mirrored(const std::vector<PolySegment>& profile, double pivot);
// `half` must end at `pivot` (its last segment's s1, to 1e-9 of the span)
std::vector<PolySegment> profile_symmetric_completion(const std::vector<PolySegment>& half, double pivot);

// the tendon mirrored about the chain abscissa `pivot` on its own chain (the
// mirrored extent [2p - s_end, 2p - s_begin] must lie on the chain); a
// single-end tendon gets its active anchor at the mirrored end
Tendon tendon_mirrored(const Model& model, const Tendon& t, double pivot);
// the tendon shifted by dx along its own chain (the shifted extent must lie on it)
Tendon tendon_translated(const Model& model, const Tendon& t, double dx);
// the same trace and force on another chain of bars of the same length
// (relative tolerance 1e-9); flips as in Tendon::flips
Tendon tendon_on_bars(const Model& model, const Tendon& t, const std::vector<Id>& bars, const std::vector<bool>& flips = {});

}  // namespace keta

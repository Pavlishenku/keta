// Section drying and fiber shrinkage — see drying.h.
#include "keta/drying.h"

#include <cmath>

#include <Eigen/SparseLU>

namespace keta {

double DryingSection::diffusivity(double hh, double* dD) const {
  if (alpha0 >= 1) {
    if (dD) *dD = 0;
    return D1;
  }
  const double x = (1 - hh) / (1 - hc), xn = std::pow(std::max(x, 0.0), n);
  const double g = 1 / (1 + xn);
  if (dD) {
    // d/dh x^n = n x^(n-1) (-1 / (1 - hc))
    const double dxn = x > 0 ? n * std::pow(x, n - 1) * (-1 / (1 - hc)) : 0.0;
    *dD = D1 * (1 - alpha0) * (-g * g * dxn);
  }
  return D1 * (alpha0 + (1 - alpha0) * g);
}

double DryingField::content() const {
  double s = 0;
  for (double v : hum) s += v;
  return s;
}

DryingField drying_init(const DryingSection& s) {
  require(s.b > 0 && s.h > 0 && s.ny >= 1 && s.nz >= 1, "drying: a b x h section on an ny x nz grid");
  require(s.D1 > 0 && s.alpha0 > 0 && s.alpha0 <= 1 && s.hc > 0 && s.hc < 1 && s.n > 0 && s.surface >= 0 &&
              s.h_env >= 0 && s.h_env <= 1 && s.h0 > 0 && s.h0 <= 1 && s.substeps >= 1,
          "drying: D1 > 0, 0 < alpha0 <= 1, 0 < hc < 1, n > 0, surface >= 0, humidities in [0, 1], substeps >= 1");
  DryingField f;
  f.hum.assign(static_cast<std::size_t>(s.ny * s.nz), s.h0);
  return f;
}

namespace {
// one backward-Euler step of h_sub from f; false (f untouched) when Newton
// does not converge; the outflow over the step in *out
bool be_step(const DryingSection& s, DryingField& f, double h_sub, double* out) {
  const int ny = s.ny, nz = s.nz, nc = ny * nz;
  const double dy = s.b / ny, dz = s.h / nz, A = dy * dz;
  auto id = [&](int i, int j) { return i * nz + j; };
  const VecX h_old = Eigen::Map<const VecX>(f.hum.data(), nc);
  VecX hn = h_old;
  // residual R(h) = A (h - h_old) / dt + the outgoing fluxes
  auto assemble = [&](const VecX& hv, VecX& R, SpMat* J, double* out_flux) {
    R = A * (hv - h_old) / h_sub;
    std::vector<Triplet> t;
    if (J)
      for (int c = 0; c < nc; ++c) t.emplace_back(c, c, A / h_sub);
    double outf = 0;
    auto inner = [&](int a, int bcell, double len, double dist) {
      double dDm;
      const double hm = 0.5 * (hv(a) + hv(bcell)), D = s.diffusivity(hm, &dDm);
      const double q = D * len / dist * (hv(a) - hv(bcell));  // from a to b
      R(a) += q;
      R(bcell) -= q;
      if (J) {
        const double dqa = len / dist * (D + 0.5 * dDm * (hv(a) - hv(bcell)));
        const double dqb = len / dist * (-D + 0.5 * dDm * (hv(a) - hv(bcell)));
        t.emplace_back(a, a, dqa);
        t.emplace_back(a, bcell, dqb);
        t.emplace_back(bcell, a, -dqa);
        t.emplace_back(bcell, bcell, -dqb);
      }
    };
    auto face = [&](int a, double len, double half) {
      double dD;
      const double D = s.diffusivity(hv(a), &dD);
      // resistance: the half cell half / D, plus the surface 1 / surface
      const double r = half / D + (s.surface > 0 ? 1 / s.surface : 0.0);
      const double q = len * (hv(a) - s.h_env) / r;
      R(a) += q;
      outf += q;
      if (J) {
        const double dr = -half * dD / (D * D);
        t.emplace_back(a, a, len / r - len * (hv(a) - s.h_env) * dr / (r * r));
      }
    };
    for (int i = 0; i < ny; ++i)
      for (int j = 0; j < nz; ++j) {
        if (i + 1 < ny) inner(id(i, j), id(i + 1, j), dz, dy);
        if (j + 1 < nz) inner(id(i, j), id(i, j + 1), dy, dz);
        if (i == 0 && s.exposed[0]) face(id(i, j), dz, dy / 2);
        if (i == ny - 1 && s.exposed[1]) face(id(i, j), dz, dy / 2);
        if (j == 0 && s.exposed[2]) face(id(i, j), dy, dz / 2);
        if (j == nz - 1 && s.exposed[3]) face(id(i, j), dy, dz / 2);
      }
    if (J) {
      J->resize(nc, nc);
      J->setFromTriplets(t.begin(), t.end());
    }
    if (out_flux) *out_flux = outf;
  };
  const double scale = A / h_sub;
  bool conv = false;
  VecX R;
  assemble(hn, R, nullptr, nullptr);
  for (int it = 0; it < 50; ++it) {
    const double rn = R.cwiseAbs().maxCoeff();
    if (rn <= 1e-14 * scale) {
      conv = true;
      break;
    }
    SpMat J;
    VecX Rj;
    assemble(hn, Rj, &J, nullptr);
    J.makeCompressed();
    Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
    lu.compute(J);
    if (lu.info() != Eigen::Success) return false;
    const VecX dh = lu.solve(Rj);
    if (!dh.allFinite()) return false;
    // backtracking on |R| (the steep front of D(h))
    double a = 1;
    VecX ht, Rt;
    bool ok = false;
    for (int ls = 0; ls < 30; ++ls, a *= 0.5) {
      ht = hn - a * dh;
      assemble(ht, Rt, nullptr, nullptr);
      if (Rt.allFinite() && Rt.cwiseAbs().maxCoeff() < rn) {
        ok = true;
        break;
      }
    }
    if (!ok) {
      if (dh.cwiseAbs().maxCoeff() <= 1e-14) {
        conv = true;
        break;
      }
      return false;
    }
    hn = ht;
    R = Rt;
    if (a == 1 && dh.cwiseAbs().maxCoeff() <= 1e-15) {
      conv = true;
      break;
    }
  }
  if (!conv) return false;
  double outf = 0;
  assemble(hn, R, nullptr, &outf);
  *out = outf * h_sub;
  for (int c = 0; c < nc; ++c) f.hum[static_cast<std::size_t>(c)] = hn(c);
  return true;
}
// the step, halved recursively (deterministically) where Newton fails
double adaptive(const DryingSection& s, DryingField& f, double dt, int depth) {
  double out = 0;
  if (be_step(s, f, dt, &out)) return out;
  require(depth < 16, "drying: the diffusion step did not converge even halved 16 times");
  const double a = adaptive(s, f, dt / 2, depth + 1);
  return a + adaptive(s, f, dt / 2, depth + 1);
}
}  // namespace

double drying_step(const DryingSection& s, DryingField& f, double dt) {
  require(dt >= 0, "drying: dt >= 0");
  if (dt == 0) return 0;
  double out = 0;
  for (int sub = 0; sub < s.substeps; ++sub) out += adaptive(s, f, dt / s.substeps, 0);
  return out;
}

}  // namespace keta

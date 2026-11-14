// Drying of a concrete section and the shrinkage of each of its fibers.
//
// The pore relative humidity h(y, z, t) of a rectangular section b x h,
// discretised as the ny x nz grid of its fibers (HystFiberSection::rect
// order: fiber i nz + j at y = -b/2 + (i + 1/2) b/ny, z = -h/2 + (j + 1/2)
// h/nz), obeys the nonlinear diffusion of Bazant & Najjar (1972):
//   dh/dt = div(D(h) grad h),
//   D(h) = D1 (alpha0 + (1 - alpha0) / (1 + ((1 - h) / (1 - hc))^n)),
// from h0 (saturated: 1) at the bar's activation, the exposed faces in the
// ambient h_env — exactly (surface = 0) or through a surface factor
// (flux = surface (h_face - h_env)). Finite volumes on the grid (two-point
// fluxes, D at the mean humidity of the two cells; at an exposed face the
// half cell in series with the surface resistance), BACKWARD EULER in time,
// each step solved by Newton with the exact Jacobian (sparse LU: D(h)
// makes it non-symmetric) and backtracking, to round-off; a step whose
// Newton fails (the steep front of D(h)) is halved, recursively and
// deterministically. The scheme conserves moisture
// exactly: the change of the cells' content equals the time step times
// the boundary fluxes.
//
// Each fiber shrinks freely by eps_sh = eps_inf (k(h) - k(h0)), k(h) =
// 1 - h^3 (Bazant's humidity function; eps_inf < 0): the section's drying
// gradient gives the self-equilibrated eigenstresses, its asymmetry a
// curvature (keta/fiber_frame.h, FiberFrame3DOptions::drying).
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta {

struct DryingSection {
  double b = 0, h = 0;  // the rectangle (fiber grid ny x nz along y x z)
  int ny = 0, nz = 0;
  // exposed faces: y = -b/2, y = +b/2, z = -h/2, z = +h/2
  bool exposed[4] = {false, false, false, false};
  double h_env = 0.6, h0 = 1.0;
  double D1 = 1e-5;  // [m^2/day] (10 mm^2/day)
  double alpha0 = 0.05, hc = 0.75, n = 15;  // alpha0 = 1: constant D = D1
  double surface = 0;                        // [m/day]; 0: h = h_env on the face
  double eps_inf = -500e-6;                  // ultimate shrinkage at h = 0 (k = 1)
  int substeps = 1;                          // diffusion steps per time step
  double diffusivity(double hh, double* dD = nullptr) const;
  double shrinkage(double hh) const { return eps_inf * ((1 - hh * hh * hh) - (1 - h0 * h0 * h0)); }
};

struct DryingField {
  std::vector<double> hum;  // per cell (fiber order)
  double content() const;   // sum of the cells' h (times their area: the moisture, in h x m^2 per unit area ratio)
};
DryingField drying_init(const DryingSection& s);
// one backward-Euler step of dt days (s.substeps sub-steps); returns the
// moisture that left through the exposed faces (h x m^2) — for the balance
double drying_step(const DryingSection& s, DryingField& f, double dt);

}  // namespace keta

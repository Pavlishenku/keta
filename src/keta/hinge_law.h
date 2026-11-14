// Moment-rotation hinge laws depending on the axial force (a plastic hinge
// with an axial-force dependent moment-rotation relation; FEMA / ASCE 41
// backbones): a family of multilinear laws
// (theta, M) tabulated at axial forces N_0 < N_1 < ..., every law with the
// same number of points; the law at N is the point-by-point linear
// interpolation between the two bracketing laws (abscissas and moments),
// descending (softening) branches allowed. Outside [N_0, N_last] the law is
// not defined (typed rejection: no extrapolation rule is assumed).
//
// This is the LAW layer: a structural solve through a softening branch
// needs displacement control past the peak; the convex link solver
// (link_nl.h) keeps rejecting non-monotone laws.
#pragma once
#include "keta/model.h"

namespace keta {

struct HingeLawFamily {
  std::vector<double> N;      // ascending axial forces
  std::vector<LinkLaw> laws;  // one per N, same point count
  void validate() const;
  LinkLaw at(double n) const;  // the interpolated law
  double moment(double theta, double n) const { return at(n).force(theta); }
};

}  // namespace keta

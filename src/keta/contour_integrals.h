// Contour integrals of a uniaxial law over a polygonal section under a plane
// of strain — the integrator of the service states (keta/section_service.h)
// and of the resistance engine (keta/ec2_section.h).
//
// Under a plane eps(y, z) = e0 + ky z - kz y the stress of a material depends
// on one coordinate, s along the strain gradient; with t across it, Green's
// theorem turns every surface integral int f(s) t^k dA into the contour
// integral -oint f(s) t^(k+1) / (k + 1) ds over the rings (outer boundaries
// counter-clockwise, voids clockwise) less the same integral round each disc
// taken out. Along an edge the strain is linear and the integrand smooth
// between the kinks of the law: each piece is integrated by Gauss-Legendre
// 16 points after the change of variable xi^2 (3 - 2 xi), flat at both ends —
// exact for the polynomial laws, at round-off for an analytic one, and it
// absorbs the fractional power of a parabola at the start of its plateau.
// No mesh, no clipping.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "keta/fiber.h"

namespace keta {
using Vec2 = Eigen::Vector2d;  // section-plane point (y, z) [m]
}

namespace keta::contour {

// A uniaxial law as the contour integrals read it: stress and tangent modulus
// at a strain, and the strains where the law has a kink (a piece of boundary
// is cut there, the integrand being smooth in between). Without `tangent`
// (an empty function) the response carries the resultants alone, K = 0.
struct Law {
  std::function<double(double)> sigma, tangent;
  std::vector<double> kinks;
};
// A disc taken out of the rings.
struct Disc {
  Vec2 centre = Vec2::Zero();
  double area = 0;
};
// Resultants (N, My, Mz) — My = int sigma z, Mz = -int sigma y, about the
// origin — and their tangent to the plane (e0, ky, kz).
struct Response {
  Vec3 R = Vec3::Zero();
  Mat3 K = Mat3::Zero();
};
// A concrete law of keta/fiber.h as the contour integrals read it, with the
// strains where it changes branch: the parabola-rectangle (concrete_pr), the
// law of EN 1992-1-1 §3.1.5 (sargin), the parabola-rectangle with a softening
// tension branch (concrete_tension_softening) or a line (elastic). Another
// kind is a typed error, prefixed by `who`.
Law concrete_law(const UniaxialLaw& law, const std::string& who);

// The response of `law` over the rings less the discs under the plane e =
// (e0, ky, kz).
Response response(const std::vector<std::vector<Vec2>>& rings, const std::vector<Disc>& discs, const Law& law,
                  const Vec3& e);

}  // namespace keta::contour

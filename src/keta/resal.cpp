// Résal effect — see resal.h.
#include "keta/resal.h"

#include <cmath>

namespace keta {

ResalFactors resal_factors(const ResalSection& s) {
  for (double v : {s.A, s.Iy, s.Iz}) require(std::isfinite(v) && v > 0, "resal: A, Iy, Iz finite and > 0");
  ResalFactors k;
  for (const ResalChord& c : s.chords) {
    for (double v : {c.y, c.z, c.dy_dx, c.dz_dx}) require(std::isfinite(v), "resal: non-finite chord data");
    require(std::isfinite(c.area) && c.area > 0, "resal: chord area finite and > 0");
    // chord force per unit of (N, My, Mz)
    const Vec3 force = c.area * Vec3(1.0 / s.A, c.z / s.Iy, -c.y / s.Iz);
    k.y += c.dy_dx * force;
    k.z += c.dz_dx * force;
  }
  return k;
}

WebShear resal_web_shear(const ResalFactors& k, double N, double Vy, double Vz, double My, double Mz) {
  for (double v : {N, Vy, Vz, My, Mz}) require(std::isfinite(v), "resal: non-finite section force");
  const Vec3 forces(N, My, Mz);
  WebShear w;
  w.relief_y = k.y.dot(forces);
  w.relief_z = k.z.dot(forces);
  w.Vy = Vy - w.relief_y;
  w.Vz = Vz - w.relief_z;
  return w;
}

}  // namespace keta

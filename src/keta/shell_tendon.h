// Prestressing tendons in slabs and shells (post-tensioned flat slabs,
// decks modelled in facets): the tendon polyline at its true position inside
// the thickness, the friction losses, and its EQUIVALENT LOADS on the
// concrete facets:
//   * the anchors: P t at the jacked start (into the member), -P t at the
//     end;
//   * each vertex: the deviation force P_out t_out - P_in t_in (towards the
//     inside of the bend; P_out / P_in the forces on either side — friction
//     drops the force across the deviation);
//   * each straight run: the wobble friction d(P t)/ds = P'(s) t,
//     P(s) = P_v exp(-k (s - s_v)) — integrated against the facets' linear
//     shape functions piece by piece (the run cut at the facet edges), Gauss
//     8 points on each piece (an analytic integrand: round-off).
// The force profile: P(s) = P0 exp(-mu alpha(s) - k s) from the jacked end,
// alpha the cumulated deviation angle (both ends jacked: the larger of the
// two). Every force acts at the tendon point; it is carried to the facet
// that holds the point's orthogonal projection, as the force plus the moment
// r x F of its eccentricity r over the mid-surface, both split between the
// facet's nodes by the linear shape functions at the projected point. The
// system is SELF-EQUILIBRATED exactly (telescoping sum), so the section
// forces of a statically determinate member equal the tendon force and its
// eccentricity, whatever the mesh.
#pragma once
#include "keta/shell.h"

namespace keta {

// Polyline of a tendon from a plan trace and a vertical profile: `plan` the
// (x, y) trace points (straight between them), `profile` the control points
// (s along the plan trace, z) — consecutive extremes (`flat` = zero slope:
// high and low points) joined by two parabolas reversing at `lambda` of the
// stretch, an extreme and a free point by one parabola with its vertex at the
// extreme, two free points by a straight line — sampled `per_segment` chords
// per profile stretch (the deviation forces of the chords converge on the
// curvature pressure P / R as the sampling refines).
struct ProfilePoint {
  double s = 0, z = 0;
  bool flat = false;
  double lambda = 0.5;  // the reversal of the stretch that STARTS here (two extremes)
};
std::vector<Vec3> tendon_polyline(const std::vector<Eigen::Vector2d>& plan,
                                  const std::vector<ProfilePoint>& profile, int per_segment = 16);

// The force at each vertex of the tendon: (before, after) the vertex — equal
// without friction, dropping across a deviation with mu.
std::vector<std::pair<double, double>> shell_tendon_forces(const ShellTendon& t);

// The equivalent loads of a tendon times `factor` on the facets: (index into
// `shells`, 18 global nodal forces in the element DOF order).
std::vector<std::pair<std::size_t, Vec18>> shell_tendon_loads(const Model& model,
                                                              const std::vector<ShellElem>& shells,
                                                              const ShellTendon& t, double factor);

// The mean of the tendon force along its length (the reference of a loss
// spread over the profile in proportion to the force).
double shell_tendon_mean_force(const ShellTendon& t);

// The mean, over the tendon's length, of the strain of the concrete along
// the tendon at its position (bonded: the tendon's own strain change) under
// the node displacements `node_disp` (model node order): the facet strain
// at the tendon point (Allman membrane + z kappa, z its height over the
// mid-surface) along the tendon's in-plane direction, Gauss 8 points on each
// facet piece of each run (the strain is linear there: exact).
double shell_tendon_mean_strain(const Model& model, const std::vector<ShellElem>& shells, const ShellTendon& t,
                                const std::vector<Vec6>& node_disp);

}  // namespace keta

// Influence surfaces on an arbitrary cloud of deck nodes.
//
// A skewed slab, a fan-shaped ramp or an irregular grillage has no
// rectangular node grid (build_node_grid refuses it). Here the plan
// positions (x, y) of the chosen nodes are triangulated — the Delaunay
// triangulation of the point set, from keta's robust mesher (mesh2d: convex
// hull as the outline, every node a vertex, no refinement, no added node)
// — and eta is interpolated LINEARLY on each triangle (barycentric): exact
// at the nodes, exact for any field linear over a triangle, zero outside
// the hull. The ordinates are the unit-load sweep of the influence line
// (one factorisation, the adjoint sweep). A point is located through a
// uniform bucket grid of the triangles (O(1) on average; mesh2::locate
// rebuilds its adjacency at every call). The favourable-zone integrals are
// exact on each triangle (the zero line of a linear field is straight: the
// triangle is clipped).
#pragma once
#include <vector>

#include "keta/influence_surface.h"
#include "keta/mesh2d.h"

namespace keta {

class ScatteredSurface {
 public:
  // from plan positions and their values (>= 3 points, not collinear, no duplicates)
  ScatteredSurface(const std::vector<mesh2::Vec2>& points, const std::vector<double>& values);
  double at(double x, double y) const;  // linear on the triangle, 0 outside the hull
  double integral() const;
  double integral_pos() const;
  double integral_neg() const;
  const mesh2::TriMesh& mesh() const { return mesh_; }
  double xmin() const { return lo_.x(); }
  double xmax() const { return hi_.x(); }
  double ymin() const { return lo_.y(); }
  double ymax() const { return hi_.y(); }

 private:
  mesh2::TriMesh mesh_;
  std::vector<double> eta_;
  mesh2::Vec2 lo_, hi_;
  int gx_ = 1, gy_ = 1;
  double cw_ = 1, ch_ = 1;
  std::vector<std::vector<int>> cells_;  // triangles whose bounding box meets the cell
  double clipped(int sign) const;
};

// Unit loads of direction `dir` on each node of `nodes` (any layout), one response.
ScatteredSurface influence_surface_scattered(const Model& model, const std::vector<Id>& nodes, const Vec3& dir,
                                             const ResponseKind& response);

double wheel_set_value(const ScatteredSurface& s, const WheelSet& ws, double rx, double ry);
// the reference swept over the bounding box of the nodes
SurfaceExtremes wheel_set_extreme(const ScatteredSurface& s, const WheelSet& ws, int steps_x = 60, int steps_y = 60);

}  // namespace keta

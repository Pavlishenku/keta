// Resultant of the internal forces across a section cut of a mixed bar +
// facet-shell model — the "design element" integration of a T-beam
// modelled as a beam plus plate elements, done EXACTLY from the element nodal forces rather than from
// smoothed plate resultants:
//
//   the cut plane passes through mesh nodes (the mesh must conform: an
//   element with nodes strictly on both sides is a typed rejection); the
//   force that the material on the POSITIVE side (n . (x - p) > 0) exerts on
//   the negative side is the sum, over the elements of the positive side
//   that touch the cut, of minus their internal nodal forces at the cut
//   nodes: bars -T^t f_end (member end forces k u - f0, span loads and
//   rigid arms included), facet shells -(k_e u - f_e) with f_e the element's
//   own consistent load (pressure, self weight, thermal initial strain of
//   the result's load case). By equilibrium of each element this is the
//   exact FE section force — moments AND shears: integrating the plate
//   fields over the width and adding the beam forces converges to it. Only
//   NODAL loads applied at the cut nodes themselves are attributed to
//   neither side (no lever arm about axes in the cut plane; a shear cut
//   through a loaded node is ambiguous by nature).
//
// SectionCutter is built once per (model, result, cut direction): node
// coordinates along the direction sorted, node -> element adjacency; every
// cut is then a binary search plus the elements touching it: O(log n + k).
// An optional strip |t . (x - p)| <= half_width (t a direction in the cut
// plane) restricts the integration to an effective width.
#pragma once
#include <vector>

#include "keta/statics.h"

namespace keta {

struct CutResultant {
  Vec3 F = Vec3::Zero();  // force of the positive side on the negative side (global)
  Vec3 M = Vec3::Zero();  // its moment about the cut point (global)
  int n_bars = 0, n_shells = 0;
};

class SectionCutter {
 public:
  // `dir` the cut normal (normalised internally); tol the coordinate
  // tolerance [m] for nodes "on" the cut.
  SectionCutter(const Model& model, const CaseResult& result, const Vec3& dir, double tol = 1e-9);
  // Cut through `point` (its coordinate along dir selects the plane);
  // optional strip: only elements whose cut nodes satisfy |t . (x - point)|
  // <= half_width (t normalised, orthogonal to dir) contribute.
  CutResultant cut(const Vec3& point, const Vec3& strip_dir = Vec3::Zero(), double half_width = -1) const;

 private:
  const Model& model_;
  Vec3 n_;
  double tol_;
  std::vector<Vec3> x_;                       // node coordinates (dense)
  std::vector<std::pair<double, int>> along_;  // (n . x, node) sorted
  struct Elem {
    std::vector<int> nodes;          // dense node indices
    std::vector<Vec6> f;             // internal nodal forces (global), per node
  };
  std::vector<Elem> elems_;
  std::vector<std::vector<int>> node_elems_;  // node -> elements
};

}  // namespace keta

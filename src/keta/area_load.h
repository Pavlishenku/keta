// Load distribution zones: a uniform area load on a
// panel bounded by bars, with no slab in the model, handed to the bars by
// the classical tributary rules — exact as geometric rules (they are the
// definition of the distribution, not an approximation of a slab):
//   * TwoWay: every point of the panel goes to its NEAREST edge (for a convex
//     panel the straight skeleton of the 45-degree yield lines: triangles and
//     trapezoids of a rectangle). Edge k receives q h(s), h(s) the depth of
//     its region along the inward normal at the abscissa s;
//   * OneWay: strips along the span direction d; a strip's load goes half to
//     each edge it ends on: edge k receives q L(s) |d x t_k| / 2, L the
//     panel chord along d through the edge point.
// h(s) is piecewise linear between the projections of the region's (or the
// panel's) vertices: the bars carry it as exact linear loads between those
// breakpoints. Sum of the bar loads = q A, exactly.
#pragma once
#include "keta/model.h"

namespace keta {

struct AreaLoadZone {
  std::vector<Vec3> poly;   // coplanar CONVEX panel (either orientation)
  Vec3 q = Vec3::Zero();    // force per unit panel area [N/m^2], global vector
  enum class Mode { TwoWay, OneWay } mode = Mode::TwoWay;
  Vec3 span = Vec3::UnitX();  // OneWay: span direction (projected on the panel plane)
  std::vector<Id> bars;     // the carrying bars (empty: every bar lying on an edge)
};

// The bar loads of the zone (LoadFrame::Global linear loads), to add to a case.
// An edge stretch with load but no bar on it is a typed error.
std::vector<Load> area_load_to_bars(const Model& model, const AreaLoadZone& zone);

}  // namespace keta

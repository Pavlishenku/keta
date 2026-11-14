// Linearly varying edge tractions on 2-D continuum elements (G15): the
// traction goes linearly from `ta` at the edge's first corner to `tb` at its
// second (straight edge, mid node at the middle: linear in the position);
// exact consistent nodal forces, per radian for an axisymmetric element.
#pragma once
#include "keta/continuum.h"

namespace keta {

std::vector<Load> continuum_edge_load_linear(const Model& model, Id element, int edge, const Vec3& ta, const Vec3& tb);

}  // namespace keta

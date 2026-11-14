// Continuum elements (G15): plane stress, plane strain, axisymmetric and
// solid isoparametric elements on the model's nodes, assembled with the
// frames and shells (keta/assembly.h) — their nodes carry the translations
// they use (plane and axisymmetric: ux, uy in the global XY plane, the
// axisymmetric radius r = x and axis y; solid: ux, uy, uz), the unused DOFs
// being auto-fixed like those of a truss.
//
// Elements: triangles T3 / T6 and quadrilaterals Q4 / Q8 (serendipity),
// tetrahedron T4, hexahedra H8 / H20 (serendipity).
// Node order: corners counter-clockwise (2-D) / bottom then top faces
// (H8), then the mid-side nodes edge by edge (T6: 01 12 20; Q8: 01 12 23
// 30; H20: the four bottom edges, the four top edges, the four verticals).
// Stiffness by Gauss rules: T3 1 point, T6 3 (axisymmetric: 3 and 7), T4
// 1 — exact on straight-sided elements; Q4 2x2, Q8 2x2 (the reduced rule),
// H8 2x2x2, H20 3x3x3.
// Consistent mass by rules exact for affine elements (T3 3, T6 7, Q4 2x2,
// Q8 3x3, T4 4, H8 2x2x2, H20 3x3x3). Axisymmetric quantities are PER RADIAN (forces, masses,
// stiffness): a full-circumference load is 2 pi times the per-radian one.
// Plane elements carry the thickness t (plane strain: any t, usually 1).
//
// Stresses (Pa): plane (sx, sy, txy, sz) — sz = nu (sx + sy) in plane
// strain, 0 in plane stress; axisymmetric (sr, sz, stheta, trz); solid (sx,
// sy, sz, txy, tyz, tzx). Recovered from the displacements at any natural
// point (continuum_stress); the nodal average of the element values is
// left to the caller.
#pragma once
#include "keta/model.h"

namespace keta {

struct ContinuumKernel {
  int dim = 2;                // 2 (plane, axisymmetric) or 3 (solid)
  std::vector<int> nodes;     // dense node indices
  MatX K, M;                  // (dim n) x (dim n), DOF order node by node
};

// shape functions and natural derivatives of an n-node element of dimension
// `dim` at the natural point xi (exposed for direct use)
void continuum_shape(int dim, int n, const double* xi, VecX& N, MatX& dN);
// the element's coordinates in its own frame (x, y[, z]) — plane elements in
// the global XY plane
ContinuumKernel continuum_kernel(const Model& model, const Continuum& c);
// stress at the natural point xi from the model's node displacements
VecX continuum_stress(const Model& model, const Continuum& c, const std::vector<Vec6>& node_disp,
                      const std::array<double, 3>& xi);
// natural coordinates of the element's nodes (for nodal stresses)
std::vector<std::array<double, 3>> continuum_node_naturals(int dim, int n);
// consistent nodal forces of a traction (force per unit area, global
// vector; per radian of circumference for an axisymmetric element: the
// traction integrated over r dl) on the element side `face` (2-D: the edge
// from corner face to corner face + 1, mid node included; solid: face 0..5
// of a hexahedron — bottom, top, then y = -1, x = +1, y = +1, x = -1 — or
// 0..3 of a tetrahedron opposite... see continuum.cpp); returned as nodal
// loads to add to a case
std::vector<Load> continuum_face_load(const Model& model, Id element, int face, const Vec3& traction);
// the same for a pressure p along the inward normal (positive = pushing
// into the element)
std::vector<Load> continuum_face_pressure(const Model& model, Id element, int face, double p);

}  // namespace keta

// 2-D constrained Delaunay triangulation with Ruppert-style refinement —
// the shared meshing seam of the kernel:
//   * slab surface meshes for the flat-shell module (outline + holes +
//     internal constraint segments, one target edge length, 20-34 deg
//     minimum-angle policy);
//   * triangulations of arbitrary holed cross-sections for the exact
//     section engine (Tri6 warping/torsion FE).
//
// Contract:
//   * BYTE-DETERMINISTIC: same input -> identical node order and identical
//     triangle list, run to run (no address-based ordering, no randomness).
//   * Every input vertex appears verbatim among the output nodes; every
//     constraint segment is covered by a chain of mesh edges (split points
//     allowed); triangles inside holes are dropped.
//   * Typed failures (keta::fail) on self-intersecting outlines, duplicate
//     points closer than 1e-12 * bbox, zero-area outlines, or a refinement
//     that cannot reach the angle policy within the size budget — a partial
//     mesh is returned with `refinement_complete = false` ONLY when the
//     caller opts in (never silently).
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta::mesh2 {

using Vec2 = Eigen::Vector2d;

struct MeshInput {
  std::vector<Vec2> outline;                    // CCW closed polygon (no repeat of first point)
  std::vector<std::vector<Vec2>> holes;         // CW or CCW; interiority decides
  std::vector<std::pair<Vec2, Vec2>> segments;  // internal constraint segments
  std::vector<Vec2> points;                     // extra vertices to honor verbatim
  double target_h = 0;        // target edge length; 0 = derive from outline (~bbox/10)
  double min_angle_deg = 30;  // internal quality policy (20..34 supported)
  bool allow_partial = false; // true: return refinement_complete=false instead of failing
  // false: no constraint pre-split and no Ruppert refinement — the
  // constrained Delaunay triangulation of the input vertices only (with a
  // convex outline: the Delaunay triangulation of the point set)
  bool refine = true;
};

struct TriMesh {
  std::vector<Vec2> nodes;
  std::vector<std::array<int, 3>> tris;   // CCW, indices into nodes
  std::vector<std::array<int, 2>> boundary_edges;   // outer + hole boundaries
  std::vector<std::array<int, 2>> constraint_edges; // covering MeshInput::segments
  bool refinement_complete = true;
  double area() const;
};

TriMesh triangulate(const MeshInput& in);

// Locate the triangle containing p (barycentric-tolerant on edges); -1 if
// outside. Deterministic walk from a fixed seed.
int locate(const TriMesh& m, const Vec2& p);

}  // namespace keta::mesh2

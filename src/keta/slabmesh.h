// Slab surface meshing + conforming bar<->shell MIXING — the lowering pass
// that turns a descriptive `Slab` record into the `ShellTri` facets the engine
// assembles, welded into the existing structure. A `Slab` is lowered here
// exactly as a `Tendon` is lowered into per-bar equivalent loads: descriptive
// record in, per-element engine units out.
//
// ===========================================================================
// WHAT `mesh_slab` DOES (in this order, all of it deterministic)
// ===========================================================================
//  1. PLANE FIT. The slab plane is the plane through `outline[0]` with the
//     Newell normal of the outline ring. Every outline / hole vertex must lie
//     within `weld_tol` of that plane (typed failure otherwise — the raise
//     step must be able to weld corner nodes back onto pre-existing nodes, so
//     the meshing planarity gate IS the welding tolerance; `Model::validate`'s
//     looser 1e-6*diag gate only guards the record itself). The in-plane frame
//     is `shell.h::plane_frame(normal)` — the ONE shared geometric definition,
//     so the mesh frame and every facet frame derive from the same seed rule.
//     The projected outline is CCW in that frame by construction (the normal
//     comes from the ring itself); the shell normals therefore follow the
//     authored winding of the outline.
//  2. IN-PLANE BAR DETECTION. A `BarKind::Frame` bar whose BOTH end nodes lie
//     within `weld_tol` of the plane, whose endpoints project inside the
//     domain (boundary allowed) and whose midpoint is STRICTLY interior is a
//     bar "in the slab". A bar running exactly along the boundary (midpoint on
//     a ring within `weld_tol`) is EXCLUDED from mixing — documented, not an
//     error. Truss / cable bars are never mixed. Detection order =
//     `model.bars` order.
//  3. CONFLICT REJECTION. Two detected bars whose segments intersect anywhere
//     but at a shared endpoint (proper crossing, T-junction onto a bar
//     interior, collinear overlap) are a typed, actionable error naming both
//     bars: silently welding them together at a mesh vertex would change the
//     structure's connectivity.
//  4. IMPOSED POINTS. Every PRE-EXISTING model node lying within `weld_tol` of
//     the plane and strictly inside the domain (farther than `weld_tol` from
//     every ring) is passed to the mesher verbatim, so e.g. a column top under
//     the slab becomes a mesh vertex and welds back onto its own node. Points
//     within `weld_tol` of an already-collected input vertex are skipped
//     (dedup). Node order = `model.nodes` order.
//  5. TRIANGULATION. `mesh2::triangulate` with target_h = `Slab::mesh_size`,
//     the detected bars as constraint segments, `min_angle_deg` /
//     `allow_partial` from the options. Mesher failures (self-intersection,
//     a bar chord leaving the domain through a hole, ...) are re-raised with
//     the slab id attached. `refinement_complete == false` (reachable only
//     under `allow_partial`) is surfaced as a STRUCTURED WARNING in the
//     result, never silently.
//  6. NODE RAISE + WELD. Mesh nodes are raised to 3-D (exactly in the plane)
//     in mesher order; a raised node reuses the nearest PRE-EXISTING model
//     node within `weld_tol` (ties broken by lower node index; nodes created
//     by this very call are never candidates), otherwise a new node is
//     appended. Two mesh nodes claiming the same model node is a typed error
//     (the welding tolerance exceeds the mesh spacing).
//  7. SHELL CREATION. One `ShellTri` per mesh triangle, in mesher order,
//     carrying the slab's thickness / material / id.
//  8. AUTO-SPLIT + LOAD RE-SLICING. Each detected bar is replaced IN PLACE
//     (same position in `model.bars`, fresh ids) by sibling bars spanning the
//     consecutive mesh vertices along its axis, so the beam and the shells
//     share their nodes — the conforming mixing that makes an eccentric
//     (`Bar::ecc`) girder act compositely with the slab. Material, section,
//     kind, roll and the Winkler bed are copied onto every segment; end
//     releases are kept ONLY at the true ends; `ecc` is interpolated linearly
//     along the axis (exact for the constant offsets of composite girders).
//     Existing bar loads on a split bar are RE-SLICED: UDL / linear loads into
//     per-segment tranches (linear intensities re-interpolated, resultant
//     preserved exactly — the tranche abscissas partition the original range
//     bit-exactly), a concentrated load onto the containing segment, thermal /
//     prestrain copied per segment (uniform per length). A bar that CANNOT be
//     re-sliced cleanly is a typed, actionable error BEFORE any mutation: a
//     tendon chain through it, an ST1 `CARA VAR` exact kernel on it, or a
//     modal-mass entry targeting it.
//
// The model is mutated only after every check has passed (checks 1-5 raise
// before any node/shell/bar is touched; the split rejections of step 8 are
// evaluated before step 6 runs).
//
// ===========================================================================
// UNITS / TOLERANCES / DETERMINISM
// ===========================================================================
// SI: m throughout. ONE tolerance, `SlabMeshOptions::weld_tol` [m], governs
// plane membership (outline planarity, bar/node detection) and node welding;
// default 1e-6 m. The mesher's own relative tolerances (duplicate points at
// 1e-12 * bbox, collinearity of split points) are documented in mesh2d.h.
//
// BYTE-DETERMINISTIC: same model + same options => identical mesh, identical
// new-node coordinates and ids, identical shell and segment ids, run to run.
// Entity creation only ever iterates std::vector state in index order — never
// an unordered container.
//
// Documented limitations (v1):
//   * a bar whose chord DIPS OUT of the domain between its endpoints (concave
//     outline, or crossing a hole) fails inside the mesher ("constraint
//     segment crosses a domain boundary"), with the slab id but not the bar id
//     attached;
//   * a boundary-coincident bar is excluded from mixing (step 2): its nodes
//     still weld onto boundary mesh vertices when coincident, but the mesh
//     does not conform to it between vertices;
//   * detection uses the bar CHORD (node-to-node line); `Bar::ecc` offsets do
//     not move the detection line (the offset girder is exactly the composite
//     use case: its NODE line lies in the slab).
#pragma once
#include <string>
#include <vector>

#include "keta/mesh2d.h"
#include "keta/model.h"

namespace keta {

struct SlabMeshOptions {
  // Welding tolerance [m]: plane membership (outline planarity gate, in-plane
  // bar/node detection) AND node reuse radius. Must be well below the mesh
  // size; two mesh nodes within weld_tol of one model node is a typed error.
  double weld_tol = 1e-6;
  double min_angle_deg = 30.0;  // mesher quality policy (20..34)
  // true: an incomplete refinement is returned as a structured warning
  // (refinement_complete = false) instead of a typed failure.
  bool allow_partial = false;
};

// Structured warning — surfaced in the result, NEVER silent.
struct SlabMeshWarning {
  enum class Kind { RefinementIncomplete };
  Kind kind = Kind::RefinementIncomplete;
  std::string message;
};

// One detected in-plane bar and what became of it. When the mesh put no vertex
// strictly inside the axis the bar is left untouched and `segments` is just
// {original}; otherwise the original id is GONE from the model and `segments`
// holds the sibling bar ids in start->end order.
struct SplitBarRecord {
  Id original = 0;
  std::vector<Id> segments;     // >= 1 bar ids, along the axis
  std::vector<Id> chain_nodes;  // model node ids along the axis, size = segments.size() + 1
  std::vector<double> s;        // cumulated abscissa [m] of each chain node (s.front() = 0)
};

struct SlabMeshResult {
  Id slab = 0;
  Vec3 origin = Vec3::Zero();  // plane origin (= outline[0])
  Mat3 frame = Mat3::Identity();  // rows = plane x, y, z in global coords (plane_frame)
  mesh2::TriMesh mesh;         // the raw 2-D mesh in that frame (probing)
  std::vector<Id> nodes;       // model node id per mesh node, mesher order
  std::size_t welded_nodes = 0;  // how many of `nodes` reused a pre-existing node
  std::vector<Id> shells;      // created ShellTri ids, mesher triangle order
  std::vector<SplitBarRecord> split_bars;  // detected bars, model order
  bool refinement_complete = true;
  std::vector<SlabMeshWarning> warnings;
};

// Mesh the Slab record `slab` of `model` (add_slab first), weld the mesh into
// the model, create its ShellTris and auto-split the in-plane bars. Typed
// failures per the module header; the slab must not already carry shells.
SlabMeshResult mesh_slab(Model& model, Id slab, const SlabMeshOptions& opt = {});

}  // namespace keta

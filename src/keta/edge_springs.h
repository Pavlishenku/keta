// Springs along a line of shell (or bar) nodes — a facet edge on an elastic
// support, a semi-rigid joint between two slabs — as elastic node-to-node
// links (model.h ElasticLink), so every analysis that takes links takes them
// (and their non-linear laws). The stiffness per unit length k (6 x 6, in the
// frame of each segment: x along it, z towards `up`, y = z x x) is LUMPED to
// the nodes: each node takes half of its two segments, k L / 2 rotated to
// the global axes — exact for a relative displacement uniform along the
// line, the discretisation of the line spring otherwise.
// `other`: the nodes facing `line` one by one (a duplicated line: a joint);
// empty: the ground (a clamped node is created at each line node).
// Hinge lines: duplicate the nodes and couple the translations by rigid
// links (model.h RigidLink, mask {1,1,1,0,0,0}), plus a rotational edge
// spring for a semi-rigid hinge.
#pragma once
#include "keta/model.h"

namespace keta {

// A STRAIGHT EDGE for Allman facets: the drilling rotations of a line of
// nodes tied together (rigid links on the rotation about the facets' normal,
// which must be a global axis), so an edge whose nodes are imposed — a
// support line, a loaded edge — stays straight between them (Allman's edge
// displacement carries L/8 (theta_j - theta_i) on the normal). Returns the
// rigid link's id.
Id add_straight_edge(Model& model, const std::vector<Id>& line);

std::vector<Id> add_edge_springs(Model& model, const std::vector<Id>& line, const std::vector<Id>& other,
                                 const Mat6& k_per_length, const Vec3& up = Vec3::UnitZ());

}  // namespace keta

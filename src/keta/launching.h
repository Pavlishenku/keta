#pragma once
// Incremental launching along a constant-curvature path: straight, circular
// or helical.
//
// A girder can be pushed along its own axis only if that axis is a curve of
// constant curvature and torsion — a straight line, a circle in plan or a
// helix (a circle in plan on a constant grade): the launch is then a rigid
// screw motion that maps the curve onto itself and keeps the direction of
// gravity, so every position is congruent to the girder in its reference
// place with the bearings moved along it. The mesh is immutable and the
// SUPPORTS move: a pier at ground abscissa P bears, at advance `adv`, on the
// girder point of path abscissa
//     sigma = P - adv + Lg       (from the girder's first chain node),
// i.e. adv = 0 puts the girder's front end (sigma = Lg) at ground 0. The
// bearing stands at that exact abscissa: a contact between two chain nodes
// cuts the chain bar there (the same fraction of its chord as of the path
// interval) in the model solved for that position; a contact within 1e-7 of
// a bar length from a node bears on the node. The cut bar must be a
// prismatic frame bar on its node line, outside any tendon (typed refusal).
//
// Each engaged bearing holds the girder vertically and laterally in the
// frame of the path at its contact (x the plan tangent, y the plan normal,
// z global Z: a SkewSupport); the first engaged pier of the list also holds
// the longitudinal direction and the twist about the plan tangent (the
// launching jacks). A position with fewer than two piers under the girder
// is reported unstable and not solved.
#include <vector>

#include "keta/model.h"

namespace keta {

struct LaunchPath {
  enum class Kind { Straight, Circular, Helical };
  Kind kind = Kind::Straight;
  Vec3 origin = Vec3::Zero();  // the point of abscissa 0
  double heading = 0.0;        // plan angle of the tangent at 0, from +X towards +Y [rad]
  double radius = 0.0;         // Circular / Helical: plan radius [m], > 0 turning left (towards +Y at heading 0)
  double grade = 0.0;          // Straight / Helical: rise per plan length dz/ds_plan (Circular: 0)
};

// Position and unit tangent of the path at 3-D arc length s (s may be
// negative: the curve before its origin).
struct LaunchPoint {
  Vec3 x = Vec3::Zero();
  Vec3 t = Vec3::UnitX();
};
LaunchPoint launch_path_point(const LaunchPath& path, double s);

struct LaunchOptions {
  int stations_per_bar = 5;
  // Largest distance of a girder node from the path [m]: a girder off its
  // path cannot slide along it (typed refusal).
  double path_tol = 1e-6;
};

struct LaunchPosition {
  double advance = 0;
  bool stable = false;            // two engaged bearings or more: solved
  // per pier: path abscissa of its contact from the girder's first chain
  // node [m] (NaN: off the girder). Piers whose contacts coincide share
  // one bearing; the first of the list carries its reaction and guide.
  std::vector<double> contact_s;
  std::vector<double> reaction;   // per pier: vertical reaction [N], > 0 bearing, < 0 holding down
  std::vector<double> guide;      // per pier: lateral guide force along the plan normal [N]
  std::vector<std::size_t> uplift;  // the piers whose vertical reaction is negative (the girder lifts off)
  // free lengths of the girder beyond its outermost engaged piers [m], from
  // the exact contact abscissae: ahead of the front-most one (the
  // cantilever of the launching nose) and behind the rear-most one
  double front_overhang = 0, rear_overhang = 0;
  // per pier: girder moment My at the contact section [N m] (0 off the girder)
  std::vector<double> moment;
};

struct LaunchPathResults {
  std::vector<LaunchPosition> positions;
  // envelopes of the station moments My over the stable positions (bars
  // of the girder model in order, `stations_per_bar` stations each, on the
  // girder's own bars whatever the cuts)
  std::vector<double> m_min, m_max;
};

// `girder`: the girder model in its reference place (its own supports,
// cases and combinations are ignored), `chain`: its bars in launching order
// (from the rear end to the front end), on `path`. Self-weight per position.
LaunchPathResults launch_on_path(const Model& girder, const std::vector<Id>& chain, const LaunchPath& path,
                                 const std::vector<double>& pier_abscissae, const std::vector<double>& advances,
                                 const LaunchOptions& opt = {});

}  // namespace keta

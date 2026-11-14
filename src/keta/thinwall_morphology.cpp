// Wall graph from the contour of a section — see "Morphology" in thinwall.h.
#include <algorithm>
#include <cmath>
#include <limits>

#include "keta/thinwall.h"

namespace keta {

namespace {

// Two faces of a wall whose lengths differ by more than this share of the
// longer do not face each other: the wall has no midline to speak of. (The
// thicknesses of its two ends may differ freely: a tapered wall carries them.)
constexpr double kLargestFaceMismatch = 0.5;
constexpr double kDefaultSnapShare = 1e-6;  // of the extent of the contour

double contour_extent(const std::vector<Vec2>& contour) {
  require(!contour.empty(), "thinwall morphology: empty contour");
  Vec2 lo = contour.front(), hi = contour.front();
  for (const Vec2& p : contour) {
    require(std::isfinite(p.x()) && std::isfinite(p.y()), "thinwall morphology: non-finite contour point");
    lo = lo.cwiseMin(p);
    hi = hi.cwiseMax(p);
  }
  const double extent = (hi - lo).maxCoeff();
  require(extent > 0, "thinwall morphology: the contour has no extent");
  return extent;
}

const Vec2& contour_point(const std::vector<Vec2>& contour, int index, const std::string& owner) {
  require(index >= 0 && index < static_cast<int>(contour.size()),
          "thinwall morphology: " + owner + " references contour point " + std::to_string(index) + ", out of range");
  return contour[static_cast<std::size_t>(index)];
}

Vec2 left_normal(const Vec2& unit) { return Vec2(-unit.y(), unit.x()); }

void require_facing_faces(double length_a, double length_b, const std::string& owner) {
  require(std::fabs(length_a - length_b) <= kLargestFaceMismatch * std::max(length_a, length_b),
          "thinwall morphology: " + owner + " has no determinate midline (its faces " + std::to_string(length_a) +
              " and " + std::to_string(length_b) + " long differ by more than half)");
}

WallSegment wall_between(int a, int b, double t_a, double t_b) {
  WallSegment w{a, b, t_a};
  if (t_b != t_a) w.t_b = t_b;
  return w;
}

}  // namespace

ThinWall thinwall_from_quads(const std::vector<Vec2>& contour, const std::vector<std::array<int, 4>>& walls,
                             double snap) {
  require(!walls.empty(), "thinwall morphology: no wall");
  require(std::isfinite(snap) && snap >= 0, "thinwall morphology: snap tolerance >= 0");
  const double tol = snap > 0 ? snap : kDefaultSnapShare * contour_extent(contour);
  ThinWall out;
  auto node_at = [&](const Vec2& p) {
    for (std::size_t i = 0; i < out.nodes.size(); ++i)
      if ((out.nodes[i] - p).norm() <= tol) return static_cast<int>(i);
    out.nodes.push_back(p);
    return static_cast<int>(out.nodes.size()) - 1;
  };
  for (std::size_t k = 0; k < walls.size(); ++k) {
    const std::string owner = "wall " + std::to_string(k);
    const Vec2 p = contour_point(contour, walls[k][0], owner), q = contour_point(contour, walls[k][1], owner);
    const Vec2 r = contour_point(contour, walls[k][2], owner), s = contour_point(contour, walls[k][3], owner);
    const Vec2 start = 0.5 * (p + s), end = 0.5 * (q + r);
    const double length = (end - start).norm(), face_pq = (q - p).norm(), face_sr = (r - s).norm();
    require(length > tol && face_pq > tol && face_sr > tol,
            "thinwall morphology: " + owner + " has no determinate midline (degenerate length)");
    require_facing_faces(face_pq, face_sr, owner);
    const Vec2 normal = left_normal((end - start) / length);
    const double t_start = std::fabs((p - s).dot(normal)), t_end = std::fabs((q - r).dot(normal));
    require(t_start > 0 && t_end > 0, "thinwall morphology: " + owner + " has no determinate midline (no thickness)");
    const int a = node_at(start), b = node_at(end);
    require(a != b, "thinwall morphology: " + owner + " has no determinate midline (its ends merge)");
    out.segs.push_back(wall_between(a, b, t_start, t_end));
  }
  out.validate();
  return out;
}

ThinWall thinwall_from_nodes(const std::vector<Vec2>& contour, const std::vector<std::vector<int>>& nodes,
                             const std::vector<std::array<int, 2>>& branches) {
  require(!nodes.empty(), "thinwall morphology: no node");
  require(!branches.empty(), "thinwall morphology: no branch");
  const double tol = kDefaultSnapShare * contour_extent(contour);

  // step 1: the centroid of the contour points of each node
  std::vector<Vec2> centroid;
  for (std::size_t n = 0; n < nodes.size(); ++n) {
    const std::string owner = "node " + std::to_string(n);
    require(nodes[n].size() >= 2, "thinwall morphology: " + owner + " lists fewer than two contour points");
    Vec2 sum = Vec2::Zero();
    for (int index : nodes[n]) sum += contour_point(contour, index, owner);
    centroid.push_back(sum / static_cast<double>(nodes[n].size()));
  }

  // step 2: thickness and midline point of each branch end, across its axis
  struct Midline {
    Vec2 point, normal;
  };
  std::vector<std::vector<Midline>> meeting(nodes.size());
  ThinWall out;
  for (std::size_t k = 0; k < branches.size(); ++k) {
    const std::string owner = "branch " + std::to_string(k);
    const int a = branches[k][0], b = branches[k][1];
    require(a >= 0 && b >= 0 && a < static_cast<int>(nodes.size()) && b < static_cast<int>(nodes.size()),
            "thinwall morphology: " + owner + " references a missing node");
    require(a != b, "thinwall morphology: " + owner + " joins a node to itself");
    const Vec2 axis = centroid[static_cast<std::size_t>(b)] - centroid[static_cast<std::size_t>(a)];
    require(axis.norm() > tol, "thinwall morphology: " + owner + " has no determinate midline (coincident nodes)");
    const Vec2 across = left_normal(axis.normalized());
    struct End {
      double thickness;
      Vec2 midline_point;
    };
    auto end_at = [&](int node) {
      const Vec2& c = centroid[static_cast<std::size_t>(node)];
      double lo = std::numeric_limits<double>::infinity(), hi = -lo;
      for (int index : nodes[static_cast<std::size_t>(node)]) {
        const double offset = (contour[static_cast<std::size_t>(index)] - c).dot(across);
        lo = std::min(lo, offset);
        hi = std::max(hi, offset);
      }
      return End{hi - lo, c + 0.5 * (hi + lo) * across};
    };
    const End at_a = end_at(a), at_b = end_at(b);
    require(at_a.thickness > tol && at_b.thickness > tol,
            "thinwall morphology: " + owner + " has no determinate midline (no thickness)");
    const Vec2 along = at_b.midline_point - at_a.midline_point;
    require(along.norm() > tol, "thinwall morphology: " + owner + " has no determinate midline (its ends merge)");
    const Vec2 normal = left_normal(along.normalized());
    meeting[static_cast<std::size_t>(a)].push_back({at_a.midline_point, normal});
    meeting[static_cast<std::size_t>(b)].push_back({at_b.midline_point, normal});
    out.segs.push_back(wall_between(a, b, at_a.thickness, at_b.thickness));
  }

  // step 3: each node at the least-squares intersection of its midlines,
  // (sum n n^T) x = sum n (n . p)
  constexpr double kParallel = 1e-12;
  for (std::size_t n = 0; n < nodes.size(); ++n) {
    require(!meeting[n].empty(), "thinwall morphology: node " + std::to_string(n) + " carries no branch");
    Eigen::Matrix2d normal_matrix = Eigen::Matrix2d::Zero();
    Vec2 rhs = Vec2::Zero(), mean = Vec2::Zero();
    for (const Midline& line : meeting[n]) {
      normal_matrix += line.normal * line.normal.transpose();
      rhs += line.normal * line.normal.dot(line.point);
      mean += line.point / static_cast<double>(meeting[n].size());
    }
    const bool crossing = std::fabs(normal_matrix.determinant()) > kParallel * static_cast<double>(meeting[n].size());
    out.nodes.push_back(crossing ? Vec2(normal_matrix.inverse() * rhs) : mean);
  }
  out.validate();
  return out;
}

}  // namespace keta

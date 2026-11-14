#include "keta/edge_springs.h"

#include <algorithm>
#include <cmath>

#include <unordered_set>

namespace keta {

Id add_straight_edge(Model& model, const std::vector<Id>& line) {
  require(line.size() >= 2, "straight edge: a line needs two nodes at least");
  model.index();
  // the normal of the facets on the line
  Vec3 n = Vec3::Zero();
  const std::unordered_set<Id> on_line(line.begin(), line.end());
  for (const ShellTri& sh : model.shells) {
    int hits = 0;
    for (Id v : sh.n) hits += on_line.count(v) > 0;
    if (hits == 0) continue;
    const Vec3 a = model.node(sh.n[0]).x, b = model.node(sh.n[1]).x, c = model.node(sh.n[2]).x;
    Vec3 m = (b - a).cross(c - a).normalized();
    if (n.squaredNorm() > 0 && m.dot(n) < 0) m = -m;
    if (n.squaredNorm() == 0) n = m;
    require((m - n).norm() < 1e-9, "straight edge: the facets along the line are not coplanar");
  }
  require(n.squaredNorm() > 0, "straight edge: no facet touches the line");
  int axis = -1;
  for (int c = 0; c < 3; ++c)
    if (std::fabs(std::fabs(n(c)) - 1.0) < 1e-12) axis = c;
  require(axis >= 0, "straight edge: the facets' normal must be a global axis (the rigid links' masks are global)");
  std::array<bool, 6> mask{};
  mask[static_cast<std::size_t>(3 + axis)] = true;
  return model.add_rigid_link(line.front(), std::vector<Id>(line.begin() + 1, line.end()), mask);
}

std::vector<Id> add_edge_springs(Model& model, const std::vector<Id>& line, const std::vector<Id>& other,
                                 const Mat6& k, const Vec3& up) {
  require(line.size() >= 2, "edge springs: a line needs two nodes at least");
  require(other.empty() || other.size() == line.size(), "edge springs: one facing node per line node");
  require(k.allFinite(), "edge springs: non-finite stiffness");
  model.index();
  std::vector<Mat6> kn(line.size(), Mat6::Zero());
  for (std::size_t i = 0; i + 1 < line.size(); ++i) {
    const Vec3 a = model.node(line[i]).x, b = model.node(line[i + 1]).x;
    const double L = (b - a).norm();
    require(L > 0, "edge springs: two coincident consecutive nodes");
    const Vec3 ex = (b - a) / L;
    Vec3 ey = up.cross(ex);
    require(ey.norm() > 1e-12, "edge springs: `up` is along the line");
    ey.normalize();
    const Vec3 ez = ex.cross(ey);
    Mat6 R = Mat6::Zero();
    Mat3 r;
    r.row(0) = ex;
    r.row(1) = ey;
    r.row(2) = ez;
    R.block<3, 3>(0, 0) = r;
    R.block<3, 3>(3, 3) = r;
    const Mat6 kg = R.transpose() * k * R * (L / 2);
    kn[i] += kg;
    kn[i + 1] += kg;
  }
  std::vector<Id> ids;
  for (std::size_t i = 0; i < line.size(); ++i) {
    Id g;
    if (other.empty()) {
      const Vec3 x = model.node(line[i]).x;
      g = model.add_node(x.x(), x.y(), x.z());
      model.add_support(g, {true, true, true, true, true, true});
    } else {
      g = other[i];
    }
    ElasticLink l;
    l.ni = g;
    l.nj = line[i];
    l.k = 0.5 * (kn[i] + kn[i].transpose());
    l.axes = Mat3::Identity();
    const Id id = model.add_elastic_link(g, line[i], l.k, 0.5);
    model.links.back().axes = Mat3::Identity();
    ids.push_back(id);
  }
  return ids;
}

}  // namespace keta

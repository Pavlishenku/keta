// Section-cut resultants — see section_cut.h.
#include "keta/section_cut.h"

#include <algorithm>
#include <cmath>

#include "keta/assembly.h"

namespace keta {

SectionCutter::SectionCutter(const Model& model, const CaseResult& r, const Vec3& dir, double tol)
    : model_(model), n_(dir.normalized()), tol_(tol) {
  require(dir.norm() > 0, "section cut: zero direction");
  const Assembled asm_ = assemble_elements(model);
  const std::size_t nn = model.nodes.size();
  require(r.node_disp.size() == nn, "section cut: result does not match the model");
  x_.resize(nn);
  for (std::size_t i = 0; i < nn; ++i) {
    x_[i] = model.nodes[i].x;
    along_.emplace_back(n_.dot(x_[i]), static_cast<int>(i));
  }
  std::sort(along_.begin(), along_.end());
  node_elems_.assign(nn, {});
  // bars: global nodal internal forces T^t f_end (forces of the nodes on the member)
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) continue;
    require(e < r.bars.size(), "section cut: bar results missing");
    const Vec12 fg = asm_.elements[e].T().transpose() * r.bars[e].end_forces_local;
    Elem el;
    el.nodes = {model.node_index(b.ni), model.node_index(b.nj)};
    el.f = {fg.head<6>(), fg.tail<6>()};
    for (int k : el.nodes) node_elems_[static_cast<std::size_t>(k)].push_back(static_cast<int>(elems_.size()));
    elems_.push_back(el);
  }
  // facet shells: k_e u - f_e, f_e the element's own consistent load share
  // (pressure, self weight, thermal) of the result's load case
  std::vector<Vec18> shell_f(asm_.shells.size(), Vec18::Zero());
  if (!asm_.shells.empty()) {
    const LoadCase* lc = nullptr;
    for (const LoadCase& c : model.cases)
      if (c.id == r.id && c.name == r.name) lc = &c;
    require(lc != nullptr,
            "section cut: the result is not a load case of the model (a shell's own load share is needed)");
    shell_f = lower_case(model, asm_, *lc).shell_f;
  }
  for (const ShellElem& se : asm_.shells) {
    Vec18 u;
    for (int k = 0; k < 3; ++k) u.segment<6>(6 * k) = r.node_disp[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(k)])];
    const Vec18 f = se.el.ke() * u - shell_f[se.index];
    Elem el;
    el.nodes = {se.nodes[0], se.nodes[1], se.nodes[2]};
    el.f = {f.segment<6>(0), f.segment<6>(6), f.segment<6>(12)};
    for (int k : el.nodes) node_elems_[static_cast<std::size_t>(k)].push_back(static_cast<int>(elems_.size()));
    elems_.push_back(el);
  }
}

CutResultant SectionCutter::cut(const Vec3& p, const Vec3& strip_dir, double half_width) const {
  CutResultant out;
  const double s0 = n_.dot(p);
  auto lo = std::lower_bound(along_.begin(), along_.end(), std::make_pair(s0 - tol_, -1));
  auto hi = std::upper_bound(along_.begin(), along_.end(), std::make_pair(s0 + tol_, 1 << 30));
  require(lo != hi, "section cut: no mesh node on the cut plane");
  Vec3 t = Vec3::Zero();
  if (half_width >= 0) {
    require(strip_dir.norm() > 0, "section cut: strip direction is zero");
    t = (strip_dir - n_.dot(strip_dir) * n_).normalized();
  }
  std::vector<int> seen;
  for (auto it = lo; it != hi; ++it) {
    const int node = it->second;
    if (half_width >= 0 && std::fabs(t.dot(x_[static_cast<std::size_t>(node)] - p)) > half_width + tol_) continue;
    for (int ei : node_elems_[static_cast<std::size_t>(node)]) {
      const Elem& el = elems_[static_cast<std::size_t>(ei)];
      bool pos = false, neg = false;
      for (int k : el.nodes) {
        const double s = n_.dot(x_[static_cast<std::size_t>(k)]) - s0;
        pos |= s > tol_;
        neg |= s < -tol_;
      }
      require(!(pos && neg), "section cut: an element straddles the cut (the mesh must conform to it)");
      if (!pos) continue;
      // contribution of this element at this cut node
      for (std::size_t k = 0; k < el.nodes.size(); ++k) {
        if (el.nodes[k] != node) continue;
        const Vec6& f = el.f[k];
        const Vec3 F = -f.head<3>(), Mn = -f.tail<3>();
        out.F += F;
        out.M += Mn + (x_[static_cast<std::size_t>(node)] - p).cross(F);
      }
      if (std::find(seen.begin(), seen.end(), ei) == seen.end()) {
        seen.push_back(ei);
        (el.nodes.size() == 2 ? out.n_bars : out.n_shells)++;
      }
    }
  }
  return out;
}

}  // namespace keta

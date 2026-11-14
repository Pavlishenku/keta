#include "keta/composite.h"

#include <algorithm>

#include "keta/element.h"

namespace keta {

CompositeState composite_state(const Model& model, Id bar, const std::vector<double>& E_parts) {
  const auto it = model.composite.find(bar);
  require(it != model.composite.end(), "composite: bar " + std::to_string(bar) + " has no parts");
  const std::vector<SectionPart>& parts = it->second.parts;
  require(E_parts.size() == parts.size(), "composite: one modulus per part");
  CompositeState s;
  s.E = E_parts;
  double mass = 0;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (E_parts[i] <= 0) continue;
    const Material& m = model.material(parts[i].material);
    if (s.E_ref == 0) {
      s.E_ref = E_parts[i];
      s.nu_ref = m.nu;
    }
  }
  require(s.E_ref > 0, "composite: bar " + std::to_string(bar) + " has no active part");
  const double G_ref = s.E_ref / (2 * (1 + s.nu_ref));
  double A = 0, Sz = 0;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (E_parts[i] <= 0) continue;
    const double n = E_parts[i] / s.E_ref;
    A += n * parts[i].props.A;
    Sz += n * parts[i].props.A * parts[i].z;
    mass += model.material(parts[i].material).rho * parts[i].props.A;
  }
  s.zc = Sz / A;
  SectionProps& p = s.props;
  p.A = A;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (E_parts[i] <= 0) continue;
    const SectionPart& q = parts[i];
    const Material& m = model.material(q.material);
    const double n = E_parts[i] / s.E_ref;
    const double g = E_parts[i] / (2 * (1 + m.nu)) / G_ref;
    const double d = q.z - s.zc;
    p.Iy += n * (q.props.Iy + q.props.A * d * d);
    p.Iz += n * q.props.Iz;
    p.J += g * q.props.J;
    p.Avy += g * q.props.Avy;
    p.Avz += g * q.props.Avz;
  }
  s.rho = mass / A;
  return s;
}

void lower_composite_bar(Model& w, Id bar, const CompositeState& s) {
  Bar* b = nullptr;
  for (Bar& x : w.bars)
    if (x.id == bar) b = &x;
  require(b != nullptr, "composite: bar " + std::to_string(bar) + " not in the working model");
  Id mid = 0, sid = 0;
  for (const Material& m : w.materials) mid = std::max(mid, m.id);
  for (const Section& c : w.sections) sid = std::max(sid, c.id);
  Material m = w.material(b->material);
  m.id = mid + 1;
  m.E = s.E_ref;
  m.nu = s.nu_ref;
  m.rho = s.rho;
  m.concrete.reset();
  m.concrete_law.reset();
  w.materials.push_back(m);
  Section c;
  c.id = sid + 1;
  c.name = "composite " + std::to_string(bar);
  c.shape = ShapeGeneric{};
  c.props = s.props;
  w.sections.push_back(c);
  b->material = m.id;
  b->section = c.id;
  const Mat3 R = bar_rotation(w.node(b->ni).x, w.node(b->nj).x, b->roll);
  const Vec3 arm = R.row(2).transpose() * s.zc;
  b->ecc = {arm, arm};
  w.invalidate_index();
}

Model composite_lowered(const Model& model) {
  Model w = model;
  for (const auto& [bid, cb] : model.composite) {
    std::vector<double> E;
    for (const SectionPart& q : cb.parts) E.push_back(model.material(q.material).E);
    lower_composite_bar(w, bid, composite_state(model, bid, E));
  }
  w.composite.clear();
  w.invalidate_index();
  return w;
}

}  // namespace keta

// Prestress losses at the level of a section — see section_losses.h.
#include "keta/section_losses.h"

#include <algorithm>
#include <limits>
#include <string>

#include "keta/creep.h"
#include "keta/profiles.h"
#include "keta/section_exact.h"

namespace keta::ec2 {

namespace {

// Concrete section: the rings (outer boundaries counter-clockwise, holes
// clockwise) by Green's integrals of each ring, a hole counted negative, less
// the displaced discs.
PolygonInertia concrete_section(const std::vector<std::vector<Vec2>>& rings, const std::vector<DisplacedDisc>& discs) {
  double A = 0, Sy = 0, Sz = 0, Iy = 0, Iz = 0, Iyz = 0;  // about the origin
  for (const auto& ring : rings) {
    std::vector<Vec3> outline;
    for (const Vec2& p : ring) outline.emplace_back(p.x(), p.y(), 0.0);
    const PolygonInertia r = polygon_inertia(outline);
    const double sign = polygon_area(ring) > 0 ? 1.0 : -1.0, a = sign * r.A;
    A += a;
    Sy += a * r.zc;
    Sz += a * r.yc;
    Iy += sign * r.Iy + a * r.zc * r.zc;
    Iz += sign * r.Iz + a * r.yc * r.yc;
    Iyz += sign * r.Iyz + a * r.yc * r.zc;
  }
  for (const DisplacedDisc& d : discs) {
    const double y = d.centre.x(), z = d.centre.y(), own = d.area * d.area / (4 * kPi);  // about a diameter
    A -= d.area;
    Sy -= d.area * z;
    Sz -= d.area * y;
    Iy -= d.area * z * z + own;
    Iz -= d.area * y * y + own;
    Iyz -= d.area * y * z;
  }
  require(A > 0, "section losses: the rings enclose no concrete");
  PolygonInertia g;
  g.A = A;
  g.yc = Sz / A;
  g.zc = Sy / A;
  g.Iy = Iy - A * g.zc * g.zc;
  g.Iz = Iz - A * g.yc * g.yc;
  g.Iyz = Iyz - A * g.yc * g.zc;
  require(g.Iy * g.Iz - g.Iyz * g.Iyz > 0, "section losses: degenerate concrete section");
  return g;
}

// Stress at (y, z) of the linear section g under (N, My, Mz) at the origin:
// N / A + a (z - zc) + b (y - yc), the gradients from the centroidal moments
// My_c = a Iy + b Iyz and Mz_c = -(a Iyz + b Iz).
Vec2 stress_gradients(const PolygonInertia& g, const Vec3& action) {
  const double N = action(0), My = action(1) - N * g.zc, Mz = action(2) + N * g.yc;
  const double det = g.Iy * g.Iz - g.Iyz * g.Iyz;
  return Vec2((My * g.Iz + Mz * g.Iyz) / det, -(Mz * g.Iy + My * g.Iyz) / det);
}
double plane_stress(const PolygonInertia& g, const Vec3& action, double y, double z) {
  const Vec2 ab = stress_gradients(g, action);
  return action(0) / g.A + ab(0) * (z - g.zc) + ab(1) * (y - g.yc);
}

// Distance from p to the boundaries of the rings, and whether p is inside
// more outlines than voids.
double distance_to_rings(const std::vector<std::vector<Vec2>>& rings, const Vec2& p) {
  double nearest = std::numeric_limits<double>::infinity();
  for (const auto& ring : rings)
    for (std::size_t i = 0, j = ring.size() - 1; i < ring.size(); j = i++) {
      const Vec2 edge = ring[i] - ring[j];
      const double length2 = edge.squaredNorm();
      const double t = length2 > 0 ? std::clamp((p - ring[j]).dot(edge) / length2, 0.0, 1.0) : 0.0;
      nearest = std::min(nearest, (p - ring[j] - t * edge).norm());
    }
  return nearest;
}
bool inside_rings(const std::vector<std::vector<Vec2>>& rings, const Vec2& p) {
  int outlines = 0;
  for (const auto& ring : rings)
    if (point_in_polygon(ring, p)) outlines += polygon_area(ring) > 0 ? 1 : -1;
  return outlines > 0;
}
// A length below which a point is on a boundary: 1e-9 of the size of the rings.
double boundary_tolerance(const std::vector<std::vector<Vec2>>& rings) {
  constexpr double kOnBoundary = 1e-9, kInf = std::numeric_limits<double>::infinity();
  Vec2 low(kInf, kInf), high(-kInf, -kInf);
  for (const auto& ring : rings)
    for (const Vec2& v : ring) {
      low = low.cwiseMin(v);
      high = high.cwiseMax(v);
    }
  return rings.empty() ? 0.0 : kOnBoundary * (high - low).norm();
}
double radius(const DisplacedDisc& d) { return std::sqrt(d.area / kPi); }

// A point reinforcement of the elastic section: a tendon, or a passive bar (a
// bonded steel without initial strain).
struct Steel {
  double y = 0, z = 0, A = 0, E = 0, eps0 = 0;
  bool bonded = true;
};

// A unit tensile force at the steel, as an action at the origin.
Vec3 unit_force_at(const Steel& t) { return Vec3(1.0, t.z, -t.y); }

// The section data shared by the service state and the losses: the concrete
// of the rings and its steels — the tendons, then the bars counted.
struct ElasticSection {
  const RcPolySection& s;
  std::vector<DisplacedDisc> discs;
  PolygonInertia concrete;
  double Ecm = 0;
  double modular = 0;  // Ep / Ecm
  std::vector<Steel> steels;

  ElasticSection(const RcPolySection& section, double concrete_modulus, const PassiveSteel& bars = {})
      : s(section), Ecm(concrete_modulus) {
    require(std::isfinite(Ecm) && Ecm > 0, "section losses: Ecm > 0");
    require(s.tendons.empty() || (std::isfinite(s.strand.E) && s.strand.E > 0), "section losses: tendon modulus > 0");
    require(std::isfinite(bars.omega) && bars.omega >= 0, "section losses: bar scale omega >= 0");
    discs = displaced_concrete(s, bars);
    concrete = concrete_section(s.rings, discs);
    modular = s.strand.E / Ecm;
    for (const RcPolySection::Tendon& t : s.tendons) steels.push_back({t.y, t.z, t.Ap, s.strand.E, t.eps_p0, t.bonded});
    if (bars.omega == 0) return;
    require(std::isfinite(s.steel.Es) && s.steel.Es > 0, "section losses: bar modulus > 0");
    for (const RcPolySection::Bar& b : s.bars) {
      require(std::isfinite(b.w) && b.w >= 0, "section losses: bar weight >= 0");
      require(std::isfinite(b.eps0), "section losses: finite initial strain of a bar");
      steels.push_back({b.y, b.z, bars.omega * b.w, s.steel.Es, b.eps0, true});
    }
  }
  // The steels are the tendons alone unless bars were asked for.
  Eigen::Index size() const { return static_cast<Eigen::Index>(steels.size()); }
  const Steel& steel(Eigen::Index i) const { return steels[static_cast<std::size_t>(i)]; }
  const RcPolySection::Tendon& tendon(Eigen::Index i) const { return s.tendons[static_cast<std::size_t>(i)]; }

  // Concrete stress at steel i per unit of stress in each steel j (the row
  // of f_ij A_j), and under an action (N, My, Mz).
  Eigen::MatrixXd stress_per_steel_stress() const {
    Eigen::MatrixXd f(size(), size());
    for (Eigen::Index i = 0; i < size(); ++i)
      for (Eigen::Index j = 0; j < size(); ++j)
        f(i, j) = steel(j).A * plane_stress(concrete, unit_force_at(steel(j)), steel(i).y, steel(i).z);
    return f;
  }
  Eigen::VectorXd stress_under(const Vec3& action) const {
    Eigen::VectorXd sigma(size());
    for (Eigen::Index i = 0; i < size(); ++i) sigma(i) = plane_stress(concrete, action, steel(i).y, steel(i).z);
    return sigma;
  }

  // Steel stresses at equilibrium under `sustained`: the steel forces load
  // the concrete, sigma_c = sigma_c(sustained) - F sigma, and a bonded steel
  // follows it, sigma = E eps0 + (E / Ecm) sigma_c.
  Eigen::VectorXd steel_stresses(const Vec3& sustained) const {
    if (size() == 0) return Eigen::VectorXd(0);  // plain concrete
    const Eigen::MatrixXd f = stress_per_steel_stress();
    const Eigen::VectorXd sustained_stress = stress_under(sustained);
    Eigen::MatrixXd lhs = Eigen::MatrixXd::Identity(size(), size());
    Eigen::VectorXd rhs(size());
    for (Eigen::Index i = 0; i < size(); ++i) {
      rhs(i) = steel(i).E * steel(i).eps0;
      if (!steel(i).bonded) continue;  // a frozen force
      const double ratio = steel(i).E / Ecm;
      lhs.row(i) += ratio * f.row(i);
      rhs(i) += ratio * sustained_stress(i);
    }
    return lhs.partialPivLu().solve(rhs);
  }
  // What the concrete carries: `sustained` less the steel forces.
  Vec3 action_on_concrete(const Eigen::VectorXd& steel_stress, const Vec3& sustained) const {
    Vec3 action = sustained;
    for (Eigen::Index i = 0; i < size(); ++i) action -= steel(i).A * steel_stress(i) * unit_force_at(steel(i));
    return action;
  }
  Eigen::VectorXd concrete_stresses(const Eigen::VectorXd& steel_stress, const Vec3& sustained) const {
    return stress_under(action_on_concrete(steel_stress, sustained));
  }
};

std::vector<double> to_vector(const Eigen::VectorXd& v) { return std::vector<double>(v.data(), v.data() + v.size()); }

}  // namespace

std::vector<DisplacedDisc> displaced_concrete(const RcPolySection& s, const PassiveSteel& bars) {
  std::vector<DisplacedDisc> discs;
  for (const RcPolySection::Tendon& t : s.tendons) {
    require(std::isfinite(t.displaced_area) && t.displaced_area >= 0 && std::isfinite(t.displaced_dy) &&
                std::isfinite(t.displaced_dz),
            "displaced concrete: finite offset and area >= 0");
    if (t.displaced_area > 0) discs.push_back({Vec2(t.y + t.displaced_dy, t.z + t.displaced_dz), t.displaced_area});
  }
  if (bars.displace_concrete && bars.omega > 0)
    for (const RcPolySection::Bar& b : s.bars) {
      require(std::isfinite(b.w) && b.w >= 0, "displaced concrete: bar weight >= 0");
      if (b.w > 0) discs.push_back({Vec2(b.y, b.z), bars.omega * b.w});
    }
  const double tol = boundary_tolerance(s.rings);
  for (std::size_t i = 0; i < discs.size(); ++i) {
    const std::string which = "displaced concrete: disc " + std::to_string(i + 1);
    require(discs[i].centre.allFinite() && inside_rings(s.rings, discs[i].centre) &&
                distance_to_rings(s.rings, discs[i].centre) >= radius(discs[i]) - tol,
            which + " leaves the concrete of the rings");
    for (std::size_t j = 0; j < i; ++j)
      require((discs[i].centre - discs[j].centre).norm() >= radius(discs[i]) + radius(discs[j]) - tol,
              which + " overlaps disc " + std::to_string(j + 1));
  }
  return discs;
}

bool in_concrete(const RcPolySection& s, const std::vector<DisplacedDisc>& discs, const Vec2& p) {
  const double tol = boundary_tolerance(s.rings);
  for (const DisplacedDisc& d : discs)
    if ((p - d.centre).norm() < radius(d) - tol) return false;
  return distance_to_rings(s.rings, p) <= tol || inside_rings(s.rings, p);
}

PrestressServiceState prestress_service_state(const RcPolySection& s, double Ecm, const Vec3& sustained,
                                              const std::vector<Vec2>& points, const PassiveSteel& bars) {
  const ElasticSection e(s, Ecm, bars);
  const Eigen::VectorXd steel = e.steel_stresses(sustained);
  const Eigen::Index tendons = static_cast<Eigen::Index>(s.tendons.size());
  PrestressServiceState state;
  state.tendon_stress = to_vector(steel.head(tendons));
  state.bar_stress = to_vector(steel.tail(steel.size() - tendons));
  // the plane of the concrete stresses: at the tendons, then at the vertices
  // of the rings and at the points asked for
  const Vec3 action = e.action_on_concrete(steel, sustained);
  state.concrete_stress = to_vector(e.stress_under(action).head(tendons));
  constexpr double kInf = std::numeric_limits<double>::infinity();
  state.fibre_max.stress = -kInf;
  state.fibre_min.stress = kInf;
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      const FibreStress fibre{v.x(), v.y(), plane_stress(e.concrete, action, v.x(), v.y())};
      if (fibre.stress > state.fibre_max.stress) state.fibre_max = fibre;
      if (fibre.stress < state.fibre_min.stress) state.fibre_min = fibre;
    }
  const Vec2 ab = stress_gradients(e.concrete, action);
  state.plane = {(action(0) / e.concrete.A - ab(0) * e.concrete.zc - ab(1) * e.concrete.yc) / Ecm, ab(0) / Ecm,
                 -ab(1) / Ecm};
  state.decompressed = state.fibre_max.stress > 0;
  for (const Vec2& p : points) {
    require(p.allFinite() && in_concrete(s, e.discs, p),
            "prestress service state: point (" + std::to_string(p.x()) + ", " + std::to_string(p.y()) +
                ") is outside the concrete");
    state.point_stress.push_back(plane_stress(e.concrete, action, p.x(), p.y()));
  }
  return state;
}

SectionLosses section_losses(const RcPolySection& s, const LossState& state) {
  require(!s.tendons.empty(), "section losses: the section carries no tendon");
  const ElasticSection e(s, state.Ecm, state.bars);
  const Eigen::Index n = e.size(), tendons = static_cast<Eigen::Index>(s.tendons.size());
  require(std::isfinite(state.phi) && state.phi >= 0, "section losses: creep coefficient >= 0");
  require(std::isfinite(state.eps_cs), "section losses: finite shrinkage strain");
  require(state.relaxation.empty() || static_cast<Eigen::Index>(state.relaxation.size()) == tendons,
          "section losses: one relaxation loss per tendon, or none");
  for (double r : state.relaxation) require(std::isfinite(r) && r >= 0, "section losses: relaxation loss >= 0");
  require(state.ageing >= 0 && state.ageing <= 1 && state.relaxation_share >= 0 && state.relaxation_share <= 1,
          "section losses: ageing coefficient and relaxation share in [0, 1]");

  const Vec3 none = Vec3::Zero();
  const Eigen::VectorXd alone = e.concrete_stresses(e.steel_stresses(none), none);
  const Eigen::VectorXd stress = e.steel_stresses(state.sustained);
  const Eigen::VectorXd sigma_qp = e.concrete_stresses(stress, state.sustained);

  // (1 + (E_i / Ecm)(1 + chi phi) F) dsig = causes, one column per cause
  const Eigen::MatrixXd f = e.stress_per_steel_stress();
  Eigen::MatrixXd coupling = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd causes = Eigen::MatrixXd::Zero(n, 3);
  for (Eigen::Index i = 0; i < n; ++i) {
    const double ratio = e.steel(i).E / state.Ecm;
    coupling.row(i) += ratio * (1.0 + state.ageing * state.phi) * f.row(i);
    causes(i, 0) = -ratio * state.phi * sigma_qp(i);
    causes(i, 1) = -state.eps_cs * e.steel(i).E;
    if (i < tendons && !state.relaxation.empty())
      causes(i, 2) = state.relaxation_share * state.relaxation[static_cast<std::size_t>(i)];
  }
  const Eigen::MatrixXd parts = coupling.partialPivLu().solve(causes);
  const Eigen::VectorXd total = parts.rowwise().sum();
  // the concrete after the losses, on its short-term modulus: the stress the
  // decompression of each tendon has to undo
  const Eigen::VectorXd sigma_after = sigma_qp + f * total;

  SectionLosses out;
  for (Eigen::Index i = 0; i < n; ++i) {
    TendonLoss l;
    l.stress = stress(i);
    l.creep = parts(i, 0);
    l.shrinkage = parts(i, 1);
    l.relaxation = parts(i, 2);
    l.total = total(i);
    // the steel after its loss: its strain while the concrete at its level,
    // at the stress the losses left there, is unstressed (a bar is bonded)
    const bool bonded = i >= tendons || e.tendon(i).bonded;
    l.eps_p0_after = (stress(i) - total(i)) / e.steel(i).E - (bonded ? sigma_after(i) / state.Ecm : 0.0);
    if (i >= tendons) {
      out.bars.push_back(l);
      continue;
    }
    l.elastic_shortening = e.modular * std::fabs(alone(i));
    out.tendons.push_back(l);
  }
  // what the steels lose loads the concrete: the plane of its stress change
  Vec3 relief = Vec3::Zero();
  for (Eigen::Index i = 0; i < n; ++i) relief += e.steel(i).A * total(i) * unit_force_at(e.steel(i));
  const Vec2 ab = stress_gradients(e.concrete, relief);
  out.concrete_change = {relief(0) / e.concrete.A - ab(0) * e.concrete.zc - ab(1) * e.concrete.yc, ab(0), ab(1)};
  return out;
}

std::vector<TendonLoss> prestress_losses(const RcPolySection& s, const LossState& state) {
  return section_losses(s, state).tendons;
}

double sequential_tensioning_factor(int n_tendons) {
  require(n_tendons >= 1, "sequential tensioning: at least one tendon");
  return (n_tendons - 1) / (2.0 * n_tendons);
}

RcPolySection after_losses(const RcPolySection& s, const std::vector<TendonLoss>& losses) {
  require(losses.size() == s.tendons.size(), "after_losses: one loss per tendon");
  RcPolySection out = s;
  for (std::size_t i = 0; i < losses.size(); ++i) out.tendons[i].eps_p0 = losses[i].eps_p0_after;
  return out;
}

RcPolySection after_losses(const RcPolySection& s, const SectionLosses& losses) {
  RcPolySection out = after_losses(s, losses.tendons);
  if (losses.bars.empty()) return out;  // the bars were left out of the losses
  require(losses.bars.size() == s.bars.size(), "after_losses: one entry per bar");
  for (std::size_t i = 0; i < losses.bars.size(); ++i) out.bars[i].eps0 = losses.bars[i].eps_p0_after;
  return out;
}

std::vector<LossAtAge> prestress_loss_history(const RcPolySection& s, const LossLaws& laws,
                                              const std::vector<double>& ages) {
  constexpr double kHoursPerDay = 24.0;
  LossState state;
  state.Ecm = laws.Ecm > 0 ? laws.Ecm : ec2_ecm(laws.concrete.fck);
  state.sustained = laws.sustained;
  state.bars = laws.bars;
  const std::vector<double> sigma_pi = prestress_service_state(s, state.Ecm, laws.sustained, {}, laws.bars).tendon_stress;
  std::vector<LossAtAge> out;
  for (double age : ages) {
    require(std::isfinite(age) && age >= laws.t0, "loss history: ages at or after the transfer t0");
    state.phi = ec2_phi(age, laws.t0, laws.concrete);
    state.eps_cs = ec2_shrinkage(age, laws.concrete) - ec2_shrinkage(laws.t0, laws.concrete);
    state.relaxation.clear();
    for (double sigma : sigma_pi)
      state.relaxation.push_back(sigma * ec2_relaxation_ratio((age - laws.t0) * kHoursPerDay, sigma, laws.fpk,
                                                              laws.relaxation_class, laws.rho1000));
    const SectionLosses l = section_losses(s, state);
    out.push_back({age, l.tendons, l.bars});
  }
  return out;
}

}  // namespace keta::ec2

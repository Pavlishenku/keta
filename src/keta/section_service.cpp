// Service states beyond the uncracked one — see section_service.h.
#include "keta/section_service.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>

namespace keta::ec2 {

namespace {

using contour::Law;
using contour::Response;

struct Steel {
  double y = 0, z = 0, A = 0, eps0 = 0;
  bool bonded = true;
  const Law* law = nullptr;
};

struct Section {
  const std::vector<std::vector<Vec2>>& rings;
  std::vector<DisplacedDisc> discs;
  std::vector<Steel> steels;
  double size = 0;  // diagonal of the bounding box of the rings
};

Vec3 lever(const Steel& t) { return Vec3(1.0, t.z, -t.y); }

Response response(const Section& sec, const Law& concrete, const Vec3& e) {
  Response out = contour::response(sec.rings, sec.discs, concrete, e);
  for (const Steel& t : sec.steels) {
    const Vec3 g = lever(t);
    const double eps = t.eps0 + (t.bonded ? e.dot(g) : 0.0);
    out.R += t.A * t.law->sigma(eps) * g;
    if (t.bonded) out.K += t.A * t.law->tangent(eps) * g * g.transpose();
  }
  return out;
}

// The plane whose resultants are `target`, by a damped Newton iteration from
// the plane of the section on its initial moduli (`initial`: the linear law
// of the concrete at its initial modulus). Moments are compared to forces
// through the size of the section. A tangent that the iteration cannot use
// (nothing compressed nor bonded to answer the residual) is completed by a
// share of the initial stiffness.
// `response_at(e)` is the response of the section under the plane e, `first`
// its response at the zero plane on the initial moduli, `size` the size of the
// section and `steel_force` the largest force a steel carries at the zero
// plane (the scale of a prestress).
template <class ResponseAt>
Vec3 equilibrium(const ResponseAt& response_at, const Response& first, double size, double steel_force,
                 const Vec3& target, const std::string& who) {
  constexpr int kMaxIterations = 200, kMaxHalvings = 40;
  constexpr double kTolerance = 1e-13, kRoundOff = 1e-10;
  const Vec3 scale(1.0, 1.0 / size, 1.0 / size);
  const Mat3 K0 = scale.asDiagonal() * first.K * scale.asDiagonal();
  const Eigen::FullPivLU<Mat3> lu0(K0);
  require(lu0.isInvertible(), who + ": degenerate section");
  const double force =
      std::max({std::fabs(target(0)), (std::fabs(target(1)) + std::fabs(target(2))) / size, steel_force});
  if (force == 0) return Vec3::Zero();  // nothing acts on the section
  const auto residual = [&](const Response& r) { return Vec3((r.R - target).cwiseProduct(scale)); };
  const auto norm = [&](const Vec3& r) { return r.norm() / force; };

  Vec3 e = scale.cwiseProduct(lu0.solve(-residual(first)));
  Response now = response_at(e);
  double share = 0;
  for (int it = 0; it < kMaxIterations; ++it) {
    const Vec3 r = residual(now);
    if (norm(r) <= kTolerance) return e;
    const Mat3 K = scale.asDiagonal() * now.K * scale.asDiagonal() + share * K0;
    const Eigen::FullPivLU<Mat3> lu(K);
    bool stepped = false;
    if (lu.isInvertible()) {
      const Vec3 step = scale.cwiseProduct(lu.solve(-r));
      double t = 1;
      for (int h = 0; h < kMaxHalvings && !stepped; ++h, t /= 2) {
        const Vec3 trial = e + t * step;
        const Response there = response_at(trial);
        if (!(norm(residual(there)) < norm(r))) continue;
        e = trial;
        now = there;
        stepped = true;
      }
    }
    if (stepped) {
      share = share > 1e-8 ? share / 10 : 0.0;
      continue;
    }
    if (norm(r) <= kRoundOff) return e;  // the residual sits on its round-off
    share = std::max(10 * share, 1e-6);
    require(share <= 1e3, who + ": the section cannot carry the action");
  }
  require(false, who + ": no equilibrium found");
  return e;
}

// The state under the plane of equilibrium: the stresses of the steels, of
// the concrete at the tendons, at the vertices and at the points.
PrestressServiceState service_state(const RcPolySection& s, const Law& concrete, const Law& initial, const Law& tendon,
                                    const Law& bar, const Vec3& sustained, const std::vector<Vec2>& points,
                                    const PassiveSteel& bars, const std::string& who) {
  require(sustained.allFinite(), who + ": finite sustained action");
  require(std::isfinite(bars.omega) && bars.omega >= 0, who + ": bar scale omega >= 0");
  Section sec{s.rings, displaced_concrete(s, bars), {}, 0.0};
  constexpr double kInf = std::numeric_limits<double>::infinity();
  Vec2 low(kInf, kInf), high(-kInf, -kInf);
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      require(v.allFinite(), who + ": finite ring coordinates");
      low = low.cwiseMin(v);
      high = high.cwiseMax(v);
    }
  require(!s.rings.empty() && (high - low).norm() > 0, who + ": the rings enclose no concrete");
  sec.size = (high - low).norm();
  for (const RcPolySection::Tendon& t : s.tendons) {
    require(std::isfinite(t.y) && std::isfinite(t.z) && std::isfinite(t.eps_p0) && std::isfinite(t.Ap) && t.Ap >= 0,
            who + ": finite tendon data and Ap >= 0");
    sec.steels.push_back({t.y, t.z, t.Ap, t.eps_p0, t.bonded, &tendon});
  }
  const std::size_t tendons = sec.steels.size();
  if (bars.omega > 0)
    for (const RcPolySection::Bar& b : s.bars) {
      require(std::isfinite(b.y) && std::isfinite(b.z) && std::isfinite(b.eps0) && std::isfinite(b.w) && b.w >= 0,
              who + ": finite bar data and weight >= 0");
      sec.steels.push_back({b.y, b.z, bars.omega * b.w, b.eps0, true, &bar});
    }

  double steel_force = 0;
  for (const Steel& t : sec.steels) steel_force = std::max(steel_force, std::fabs(t.A * t.law->sigma(t.eps0)));
  const Vec3 e = equilibrium([&](const Vec3& plane) { return response(sec, concrete, plane); },
                             response(sec, initial, Vec3::Zero()), sec.size, steel_force, sustained, who);
  PrestressServiceState state;
  state.plane = {e(0), e(1), e(2)};
  for (std::size_t i = 0; i < sec.steels.size(); ++i) {
    const Steel& t = sec.steels[i];
    const double strain = state.plane.at(t.y, t.z);
    const double stress = t.law->sigma(t.eps0 + (t.bonded ? strain : 0.0));
    if (i < tendons) {
      state.tendon_stress.push_back(stress);
      state.concrete_stress.push_back(concrete.sigma(strain));
    } else {
      state.bar_stress.push_back(stress);
    }
  }
  double eps_max = -kInf, eps_min = kInf;
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) {
      const double strain = state.plane.at(v.x(), v.y());
      if (strain > eps_max) {
        eps_max = strain;
        state.fibre_max = {v.x(), v.y(), concrete.sigma(strain)};
      }
      if (strain < eps_min) {
        eps_min = strain;
        state.fibre_min = {v.x(), v.y(), concrete.sigma(strain)};
      }
    }
  state.decompressed = eps_max > 0;
  for (const Vec2& p : points) {
    require(p.allFinite() && in_concrete(s, sec.discs, p),
            who + ": point (" + std::to_string(p.x()) + ", " + std::to_string(p.y()) + ") is outside the concrete");
    state.point_stress.push_back(concrete.sigma(state.plane.at(p.x(), p.y())));
  }
  return state;
}

Law linear_law(double E) {
  return {[E](double eps) { return E * eps; }, [E](double) { return E; }, {}};
}

Law law_of(const UniaxialLaw& l) {
  return {[l](double eps) { return l.sigma(eps); }, [l](double eps) { return l.tangent(eps); }, {}};
}

// The bilinear curve of the bars (EN 1992-1-1 Fig. 3.8): Es eps to fyd, then
// the inclined branch to ftd at eps_ud, flat beyond; mirrored in compression.
Law bar_law(const SteelDesign& d) {
  const double eyd = d.fyd / d.Es, slope = d.eps_ud > eyd ? (d.ftd - d.fyd) / (d.eps_ud - eyd) : 0.0;
  return {[d, eyd, slope](double eps) {
            const double a = std::fabs(eps);
            const double sig = a <= eyd ? d.Es * a : d.fyd + slope * (std::min(a, d.eps_ud) - eyd);
            return eps < 0 ? -sig : sig;
          },
          [d, eyd, slope](double eps) {
            const double a = std::fabs(eps);
            return a <= eyd ? d.Es : a < d.eps_ud ? slope : 0.0;
          },
          {}};
}

}  // namespace

PrestressServiceState cracked_service_state(const RcPolySection& s, double Ecm, const Vec3& sustained,
                                            const std::vector<Vec2>& points, const PassiveSteel& bars) {
  const std::string who = "cracked service state";
  require(std::isfinite(Ecm) && Ecm > 0, who + ": Ecm > 0");
  require(s.tendons.empty() || (std::isfinite(s.strand.E) && s.strand.E > 0), who + ": tendon modulus > 0");
  require(bars.omega == 0 || (std::isfinite(s.steel.Es) && s.steel.Es > 0), who + ": bar modulus > 0");
  const Law concrete{[Ecm](double eps) { return eps < 0 ? Ecm * eps : 0.0; },
                     [Ecm](double eps) { return eps < 0 ? Ecm : 0.0; },
                     {0.0}};
  return service_state(s, concrete, linear_law(Ecm), linear_law(s.strand.E), linear_law(s.steel.Es), sustained,
                       points, bars, who);
}

PrestressServiceState nonlinear_service_state(const RcPolySection& s, const UniaxialLaw& concrete,
                                              const Vec3& sustained, const std::vector<Vec2>& points,
                                              const PassiveSteel& bars) {
  const std::string who = "nonlinear service state";
  const Law law = contour::concrete_law(concrete, who);
  const double E0 = concrete.sup_tangent();
  require(std::isfinite(E0) && E0 > 0, who + ": concrete law without initial modulus");
  using K = UniaxialLaw::Kind;
  require(s.tendons.empty() || s.strand.kind == K::Elastic || s.strand.kind == K::SteelBilinear ||
              s.strand.kind == K::SteelHardening || s.strand.kind == K::Strand,
          who + ": the tendon law must be elastic, steel, steel_hardening or strand");
  if (bars.omega != 0)
    require(std::isfinite(s.steel.Es) && s.steel.Es > 0 && std::isfinite(s.steel.fyd) && s.steel.fyd > 0 &&
                std::isfinite(s.steel.ftd) && s.steel.ftd >= s.steel.fyd && std::isfinite(s.steel.eps_ud),
            who + ": bar law with Es > 0 and ftd >= fyd > 0");
  const PrestressServiceState state =
      service_state(s, law, linear_law(E0), law_of(s.strand), bar_law(s.steel), sustained, points, bars, who);

  // a service state stays within the strain bounds of its laws
  constexpr double kSlack = 1 + 1e-12;
  double low = std::numeric_limits<double>::infinity();
  for (const auto& ring : s.rings)
    for (const Vec2& v : ring) low = std::min(low, state.plane.at(v.x(), v.y()));
  require(low >= concrete.eps_min() * kSlack, who + ": concrete strained past the bound of its law");
  for (const RcPolySection::Tendon& t : s.tendons) {
    const double eps = t.eps_p0 + (t.bonded ? state.plane.at(t.y, t.z) : 0.0);
    require(eps <= s.strand.eps_max() * kSlack && eps >= s.strand.eps_min() * kSlack,
            who + ": tendon strained past the bound of its law");
  }
  if (bars.omega > 0)
    for (const RcPolySection::Bar& b : s.bars)
      require(std::fabs(b.eps0 + state.plane.at(b.y, b.z)) <= s.steel.limit() * kSlack,
              who + ": bar strained past the bound of its law");
  return state;
}

namespace {

// A law of a part of a composite section with its kinks.
Law part_law(const UniaxialLaw& l, const std::string& who) {
  using K = UniaxialLaw::Kind;
  if (l.kind != K::SteelBilinear && l.kind != K::SteelHardening) return contour::concrete_law(l, who);
  const double ey = l.fy / l.E, end = l.kind == K::SteelBilinear ? l.eps_u : l.eps_ud;
  std::vector<double> kinks = {-ey, ey};
  if (end > 0) kinks.insert(kinks.end(), {-end, end});
  return {[l](double eps) { return l.sigma(eps); }, [l](double eps) { return l.tangent(eps); }, std::move(kinks)};
}

struct Composite {
  const CompositeSection& s;
  std::vector<Law> laws, initial;  // of the parts: their law, and the line of its initial modulus
  std::vector<Law> steel_laws;
  double size = 0;

  Composite(const CompositeSection& section, const std::string& who) : s(section) {
    constexpr double kInf = std::numeric_limits<double>::infinity();
    Vec2 low(kInf, kInf), high(-kInf, -kInf);
    for (const MaterialPart& part : s.parts) {
      laws.push_back(part_law(part.law, who));
      initial.push_back(linear_law(part.law.sup_tangent()));
      require(!part.rings.empty(), who + ": a part without ring");
      for (const auto& ring : part.rings)
        for (const Vec2& v : ring) {
          require(v.allFinite(), who + ": finite ring coordinates");
          low = low.cwiseMin(v);
          high = high.cwiseMax(v);
        }
    }
    size = (high - low).norm();
    require(size > 0, who + ": no part, or rings that enclose nothing");
    using K = UniaxialLaw::Kind;
    for (const PointSteel& t : s.steels) {
      require(std::isfinite(t.y) && std::isfinite(t.z) && std::isfinite(t.eps0) && std::isfinite(t.A) && t.A >= 0,
              who + ": finite steel data and area >= 0");
      require(t.law.kind == K::Elastic || t.law.kind == K::SteelBilinear || t.law.kind == K::SteelHardening ||
                  t.law.kind == K::Strand,
              who + ": the law of a point steel must be elastic, steel, steel_hardening or strand");
      steel_laws.push_back(law_of(t.law));
    }
  }
  // the response of part i alone, on its law or on its initial modulus
  Response part(std::size_t i, const Vec3& e, bool on_initial = false) const {
    return contour::response(s.parts[i].rings, {}, on_initial ? initial[i] : laws[i], e);
  }
  Response response(const Vec3& e, bool on_initial = false) const {
    Response out;
    for (std::size_t i = 0; i < s.parts.size(); ++i) {
      const Response r = part(i, e, on_initial);
      out.R += r.R;
      out.K += r.K;
    }
    for (std::size_t i = 0; i < s.steels.size(); ++i) {
      const PointSteel& t = s.steels[i];
      const Vec3 g(1.0, t.z, -t.y);
      const double eps = t.eps0 + e.dot(g);
      out.R += t.A * steel_laws[i].sigma(eps) * g;
      out.K += t.A * steel_laws[i].tangent(eps) * g * g.transpose();
    }
    return out;
  }
};

}  // namespace

Vec3 composite_resultants(const CompositeSection& s, const StrainPlane& plane) {
  const std::string who = "composite section";
  require(std::isfinite(plane.e0) && std::isfinite(plane.ky) && std::isfinite(plane.kz), who + ": finite plane");
  return Composite(s, who).response(Vec3(plane.e0, plane.ky, plane.kz)).R;
}

CompositeState composite_state(const CompositeSection& s, const Vec3& action) {
  const std::string who = "composite state";
  require(action.allFinite(), who + ": finite action");
  const Composite c(s, who);
  double steel_force = 0;
  for (std::size_t i = 0; i < s.steels.size(); ++i)
    steel_force = std::max(steel_force, std::fabs(s.steels[i].A * c.steel_laws[i].sigma(s.steels[i].eps0)));
  const Vec3 e = equilibrium([&](const Vec3& plane) { return c.response(plane); }, c.response(Vec3::Zero(), true),
                             c.size, steel_force, action, who);
  CompositeState state;
  state.plane = {e(0), e(1), e(2)};
  // a state stays within the strain bounds of its laws
  constexpr double kSlack = 1 + 1e-12;
  for (std::size_t i = 0; i < s.parts.size(); ++i) {
    const UniaxialLaw& law = s.parts[i].law;
    for (const auto& ring : s.parts[i].rings)
      for (const Vec2& v : ring) {
        const double eps = state.plane.at(v.x(), v.y());
        require(eps >= law.eps_min() * kSlack && eps <= law.eps_max() * kSlack,
                who + ": part " + std::to_string(i + 1) + " strained past the bound of its law");
      }
    state.part_resultants.push_back(c.part(i, e).R);
  }
  for (std::size_t i = 0; i < s.steels.size(); ++i) {
    const PointSteel& t = s.steels[i];
    const double eps = t.eps0 + state.plane.at(t.y, t.z);
    require(eps >= t.law.eps_min() * kSlack && eps <= t.law.eps_max() * kSlack,
            who + ": steel " + std::to_string(i + 1) + " strained past the bound of its law");
    state.steel_stress.push_back(c.steel_laws[i].sigma(eps));
  }
  return state;
}

}  // namespace keta::ec2

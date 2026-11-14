#include "keta/footbridge.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <string>

#include "keta/assembly.h"
#include "keta/element.h"

namespace keta {

namespace {

bool horizontal_lateral(PedestrianDirection dir) { return dir == PedestrianDirection::Lateral; }

// Trapezoid: 0 below a, linear to 1 on [a, b], 1 on [b, c], linear to 0 on [c, e], 0 above.
double trapezoid(double f, double a, double b, double c, double e) {
  if (f <= a || f >= e) return 0.0;
  if (f < b) return (f - a) / (b - a);
  if (f <= c) return 1.0;
  return (e - f) / (e - c);
}

void require_class(FootbridgeClass cls) {
  const int c = static_cast<int>(cls);
  require(c >= 1 && c <= 4, "footbridge class: I, II, III or IV");
}

}  // namespace

double setra_equivalent_pedestrians(double n, double xi, bool very_dense) {
  require(finite(n) && n >= 0, "setra_equivalent_pedestrians: the number of pedestrians n must be >= 0");
  if (very_dense) return 1.85 * std::sqrt(n);
  require(finite(xi) && xi >= 0, "setra_equivalent_pedestrians: the damping ratio xi must be >= 0");
  return 10.8 * std::sqrt(xi * n);
}

int setra_frequency_range(double f, PedestrianDirection dir) {
  require(finite(f) && f >= 0, "setra_frequency_range: the frequency must be >= 0 Hz");
  // (r1 lower, r1 upper, r2 lower, r2 upper, r3 upper)
  const bool lat = horizontal_lateral(dir);
  const double a = lat ? 0.3 : 1.0, b = lat ? 0.5 : 1.7, c = lat ? 1.1 : 2.1, d = lat ? 1.3 : 2.6,
               e = lat ? 2.5 : 5.0;
  if (f >= b && f <= c) return 1;
  if ((f >= a && f < b) || (f > c && f <= d)) return 2;
  if (f > d && f <= e) return 3;
  return 4;
}

int setra_comfort_range(double accel, PedestrianDirection dir) {
  require(finite(accel) && accel >= 0, "setra_comfort_range: the acceleration must be >= 0");
  const bool vert = dir == PedestrianDirection::Vertical;
  const double r1 = vert ? 0.5 : 0.15, r2 = vert ? 1.0 : 0.3, r3 = vert ? 2.5 : 0.8;
  if (accel <= r1) return 1;
  if (accel <= r2) return 2;
  if (accel <= r3) return 3;
  return 4;
}

double setra_psi(double f, PedestrianDirection dir, int harmonic) {
  require(finite(f) && f >= 0, "setra_psi: the frequency must be >= 0 Hz");
  require(harmonic == 1 || harmonic == 2, "setra_psi: harmonic 1 or 2");
  const bool lat = horizontal_lateral(dir);
  if (harmonic == 1)
    return lat ? trapezoid(f, 0.3, 0.5, 1.1, 1.3) : trapezoid(f, 1.0, 1.7, 2.1, 2.6);
  return lat ? trapezoid(f, 1.3, 1.7, 2.1, 2.5) : trapezoid(f, 2.6, 3.4, 4.2, 5.0);
}

double setra_pedestrian_force(PedestrianDirection dir, int harmonic) {
  require(harmonic == 1 || harmonic == 2, "setra_pedestrian_force: harmonic 1 or 2");
  switch (dir) {
    case PedestrianDirection::Vertical: return harmonic == 1 ? 280.0 : 70.0;
    case PedestrianDirection::Longitudinal: return harmonic == 1 ? 140.0 : 35.0;
    case PedestrianDirection::Lateral: return harmonic == 1 ? 35.0 : 7.0;
  }
  fail("setra_pedestrian_force: unknown direction");
}

double setra_crowd_density(FootbridgeClass cls, int load_case) {
  require_class(cls);
  const std::string who = "setra_crowd_density: load case " + std::to_string(load_case) + " does not apply to class " +
                          std::to_string(static_cast<int>(cls));
  switch (load_case) {
    case 1:
      if (cls == FootbridgeClass::III) return 0.5;
      if (cls == FootbridgeClass::II) return 0.8;
      fail(who + " (case 1: classes III and II)");
    case 2:
      if (cls == FootbridgeClass::I) return 1.0;
      fail(who + " (case 2: class I)");
    case 3:
      if (cls == FootbridgeClass::II) return 0.8;
      if (cls == FootbridgeClass::I) return 1.0;
      fail(who + " (case 3: classes II and I)");
    default: fail("setra_crowd_density: load case 1, 2 or 3");
  }
}

std::vector<int> setra_load_cases(FootbridgeClass cls, int range) {
  require_class(cls);
  require(range >= 1 && range <= 4, "setra_load_cases: frequency range 1 to 4");
  switch (cls) {
    case FootbridgeClass::III: return range == 1 ? std::vector<int>{1} : std::vector<int>{};
    case FootbridgeClass::II:
      if (range <= 2) return {1};
      if (range == 3) return {3};
      return {};
    case FootbridgeClass::I:
      if (range <= 2) return {2};
      if (range == 3) return {3};
      return {};
    case FootbridgeClass::IV: return {};
  }
  return {};
}

SetraSurfaceLoad setra_surface_load(FootbridgeClass cls, int load_case, PedestrianDirection dir, double f_hz,
                                    double xi, double area) {
  require(finite(area) && area > 0, "setra_surface_load: the loaded area must be > 0");
  SetraSurfaceLoad s;
  s.density = setra_crowd_density(cls, load_case);
  s.pedestrians = s.density * area;
  // random frequencies and phases (10.8 sqrt(xi N)) for case 1 and the class
  // II case 3; random phases only (1.85 sqrt(N)) for case 2 and the class I case 3
  const bool very_dense = load_case == 2 || (load_case == 3 && cls == FootbridgeClass::I);
  s.equivalent = setra_equivalent_pedestrians(s.pedestrians, xi, very_dense);
  const int harmonic = load_case == 3 ? 2 : 1;
  s.force = setra_pedestrian_force(dir, harmonic);
  s.psi = setra_psi(f_hz, dir, harmonic);
  s.amplitude = s.density * s.force * s.equivalent / s.pedestrians * s.psi;
  s.load_frequency = f_hz;
  return s;
}

LoadCase pedestrian_mode_load_case(const Model& model, const ModalResults& modal, int mode,
                                   const PedestrianModeLoad& load, Id case_id) {
  require(mode >= 0 && static_cast<std::size_t>(mode) < modal.modes.size(),
          "pedestrian_mode_load_case: mode " + std::to_string(mode) + " out of range (" +
              std::to_string(modal.modes.size()) + " modes)");
  require(load.comp >= 0 && load.comp <= 2, "pedestrian_mode_load_case: comp 0, 1 or 2 (global x, y, z)");
  require(finite(load.q), "pedestrian_mode_load_case: the amplitude q must be finite");
  const std::vector<Vec6>& shape = modal.modes[static_cast<std::size_t>(mode)];
  require(shape.size() == model.nodes.size(),
          "pedestrian_mode_load_case: the modal results were not solved on this model (" +
              std::to_string(shape.size()) + " node shapes, " + std::to_string(model.nodes.size()) + " nodes)");
  std::vector<const Bar*> bars;
  if (load.bars.empty()) {
    for (const Bar& b : model.bars)
      if (b.kind == BarKind::Frame) bars.push_back(&b);
  } else {
    std::set<Id> seen;
    for (Id id : load.bars) {
      require(seen.insert(id).second, "pedestrian_mode_load_case: bar " + std::to_string(id) + " listed twice");
      bars.push_back(&model.bar(id));
    }
  }
  require(!bars.empty(), "pedestrian_mode_load_case: no bar to load");
  const Assembled a = assemble_elements(model);
  LoadCase lc;
  lc.id = case_id;
  lc.name = "pedestrians, mode " + std::to_string(mode + 1);
  for (const Bar* b : bars) {
    const BarElement& el = a.elements[static_cast<std::size_t>(b - model.bars.data())];
    Vec12 un;
    un.head<6>() = shape[static_cast<std::size_t>(model.node_index(b->ni))];
    un.tail<6>() = shape[static_cast<std::size_t>(model.node_index(b->nj))];
    const Vec12 ul = el.local_displacements(un, Vec12::Zero());
    const double L = el.L();
    // the element field along xi = x / L: axial linear, transverse Hermite
    // cubics (v with theta_z, w with -theta_y), as polynomials c0 + c1 xi + ...
    using Poly = std::array<double, 4>;
    auto hermite = [L](double y1, double t1, double y2, double t2) -> Poly {
      // H1 = 1 - 3x^2 + 2x^3, H2 = x - 2x^2 + x^3, H3 = 3x^2 - 2x^3, H4 = -x^2 + x^3
      return {y1, L * t1, -3 * y1 - 2 * L * t1 + 3 * y2 - L * t2, 2 * y1 + L * t1 - 2 * y2 + L * t2};
    };
    const Poly pu = {ul(0), ul(6) - ul(0), 0, 0};
    const Poly pv = hermite(ul(1), ul(5), ul(7), ul(11));
    const Poly pw = hermite(ul(2), -ul(4), ul(8), -ul(10));
    const Mat3& R = el.R();  // local <- global: global comp = R(0,c) u + R(1,c) v + R(2,c) w
    Poly p{};
    for (int k = 0; k < 4; ++k) p[k] = R(0, load.comp) * pu[k] + R(1, load.comp) * pv[k] + R(2, load.comp) * pw[k];
    auto val = [&p](double x) { return p[0] + x * (p[1] + x * (p[2] + x * p[3])); };
    // the cubic is monotonic between the roots of its derivative: one bisection per piece
    std::vector<double> cuts = {0.0};
    {
      std::vector<double> crit;
      const double qa = 3 * p[3], qb = 2 * p[2], qc = p[1];
      if (std::abs(qa) > 0) {
        const double disc = qb * qb - 4 * qa * qc;
        if (disc > 0) {
          const double sq = std::sqrt(disc);
          crit = {(-qb - sq) / (2 * qa), (-qb + sq) / (2 * qa)};
        }
      } else if (std::abs(qb) > 0) {
        crit = {-qc / qb};
      }
      std::vector<double> knots = {0.0};
      for (double c : crit)
        if (c > 0 && c < 1) knots.push_back(c);
      std::sort(knots.begin(), knots.end());
      knots.push_back(1.0);
      const double scale = std::max({std::abs(p[0]), std::abs(p[1]), std::abs(p[2]), std::abs(p[3])});
      for (std::size_t k = 0; k + 1 < knots.size(); ++k) {
        double lo = knots[k], hi = knots[k + 1];
        double flo = val(lo), fhi = val(hi);
        if (!(flo * fhi < 0) || std::abs(flo) <= 1e-14 * scale || std::abs(fhi) <= 1e-14 * scale) continue;
        for (int it = 0; it < 200 && hi - lo > 1e-15; ++it) {
          const double mid = 0.5 * (lo + hi), fm = val(mid);
          if ((fm < 0) == (flo < 0)) {
            lo = mid;
            flo = fm;
          } else {
            hi = mid;
          }
        }
        cuts.push_back(0.5 * (lo + hi));
      }
    }
    cuts.push_back(1.0);
    for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
      if (cuts[k + 1] <= cuts[k]) continue;
      const double phi = val(0.5 * (cuts[k] + cuts[k + 1]));
      const double s = load.sign_cancelled ? (phi > 0 ? 1.0 : 0.0) : (phi > 0 ? 1.0 : (phi < 0 ? -1.0 : 0.0));
      if (s == 0.0) continue;
      Vec3 w = Vec3::Zero();
      w(load.comp) = s * load.q;
      const bool whole = cuts.size() == 2;
      lc.loads.push_back(whole ? Load::udl(b->id, w) : Load::udl(b->id, w, LoadFrame::Global, cuts[k] * L, cuts[k + 1] * L));
    }
  }
  return lc;
}

std::function<VecX(int)> pedestrian_modal_force(const Model& model, const ModalResults& modal, int mode,
                                                const PedestrianModeLoad& load, double f_hz, double dt) {
  require(finite(f_hz) && f_hz >= 0, "pedestrian_modal_force: the frequency must be >= 0 Hz");
  require(finite(dt) && dt > 0, "pedestrian_modal_force: the time step must be > 0");
  const LoadCase lc = pedestrian_mode_load_case(model, modal, mode, load);
  const Assembled a = assemble_elements(model);
  const VecX f = lower_case(model, a, lc).rhs;
  const double w = 2.0 * kPi * f_hz;
  return [f, w, dt](int step) -> VecX { return f * std::cos(w * dt * step); };
}

double fourier_factor(double t, double a0, const std::vector<double>& a_cos, const std::vector<double>& b_sin,
                      double period) {
  require(finite(period) && period > 0, "fourier_factor: the period must be > 0");
  const double w = 2.0 * kPi / period;
  double s = a0;
  for (std::size_t k = 0; k < a_cos.size(); ++k) s += a_cos[k] * std::cos(static_cast<double>(k + 1) * w * t);
  for (std::size_t k = 0; k < b_sin.size(); ++k) s += b_sin[k] * std::sin(static_cast<double>(k + 1) * w * t);
  return s;
}

std::vector<double> fourier_load_factors(double a0, const std::vector<double>& a_cos,
                                         const std::vector<double>& b_sin, double period, double dt, int nsteps) {
  require(finite(dt) && dt > 0, "fourier_load_factors: the time step must be > 0");
  require(nsteps >= 0, "fourier_load_factors: nsteps must be >= 0");
  std::vector<double> out;
  out.reserve(static_cast<std::size_t>(nsteps) + 1);
  for (int n = 0; n <= nsteps; ++n) out.push_back(fourier_factor(n * dt, a0, a_cos, b_sin, period));
  return out;
}

std::vector<double> high_frequency_filtered_damping(const std::vector<double>& freqs_hz,
                                                    const std::vector<double>& xi, double f_c, double d_c) {
  require(freqs_hz.size() == xi.size(), "high_frequency_filtered_damping: " + std::to_string(freqs_hz.size()) +
                                            " frequencies but " + std::to_string(xi.size()) + " damping ratios");
  require(finite(f_c) && finite(d_c) && d_c >= 0, "high_frequency_filtered_damping: finite f_c and d_c >= 0");
  std::vector<double> out(xi);
  for (std::size_t i = 0; i < out.size(); ++i)
    if (freqs_hz[i] > f_c) out[i] += d_c;
  return out;
}

}  // namespace keta

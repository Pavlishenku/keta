// keta Python bindings — the `keta._keta` extension behind `import keta`.
//
// The scripting surface: class and method names, keyword arguments, defaults,
// auto ids, tuple shapes and the KeyError / ValueError split are the public
// API that scripts depend on. Everything is SI (m, N, Pa, kg, rad); results
// come out as plain tuples and lists.
//
//   import keta as k
//   m = k.Model()
//   a, b = m.add_node(0, 0, 0), m.add_node(6, 0, 0)
//   beam = m.add_beam(a, b, section=m.add_section("IPE300"), material=m.steel())
//   m.fix(a)
//   m.add_point_load(b, fz=-10_000)
//   print(m.analyze().case(1).displacement(b))
//
// Mapping onto the keta model:
//   * cables are kept in their own id space and appended as Cable-kind bars
//     after the frame bars when a solver model is built, with one inert
//     placeholder material / section (exactly Assemble's lowering);
//   * the nonlinear, arc-length and ruin drivers use the corotational element
//     in its Sin mode, with the same settings as the Nonlinear / Arc Length
//     nodes;
//   * a modal `method` maps onto keta's dense / Lanczos switch.
// Conventions (documented in python/keta/__init__.py):
//   * default gravity is the standard 9.80665 m/s^2;
//   * catalog profiles carry their root fillets;
//   * Model.to_json / from_json use keta's model JSON;
//   * validate() reports the first problem keta's validation finds.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <charconv>
#include <map>
#include <set>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <sstream>

#include "keta/assembly.h"
#include "keta/cable.h"
#include "keta/catalog.h"
#include "keta/dynamics.h"
#include "keta/footbridge.h"
#include "keta/surch.h"
#include "keta/tangent.h"
#include "keta/accidental.h"
#include "keta/suspension.h"
#include "keta/rotating.h"
#include "keta/fiber_hyst.h"
#include "keta/fiber_frame.h"
#include "keta/supernodal.h"
#include "keta/hinge_law.h"
#include "keta/conduction.h"
#include "keta/continuum_loads.h"
#include "keta/continuum_plastic.h"
#include "keta/substructure.h"
#include "keta/eigsolve.h"
#include "keta/io.h"
#include "keta/launching.h"
#include "keta/nonlinear.h"
#include "keta/script.h"
#include "keta/sections.h"
#include "keta/ec2_section.h"
#include "keta/section_losses.h"
#include "keta/section_service.h"
#include "keta/section_sls.h"
#include "keta/wood_armer.h"
#include "keta/soil_springs.h"
#include "keta/track_bridge.h"
#include "keta/cable_tuning.h"
#include "keta/rail_dynamics.h"
#include "keta/nldyn_model.h"
#include "keta/traffic_surface.h"
#include "keta/traffic_groups.h"
#include "keta/combination.h"
#include "keta/ec4.h"
#include "keta/link_nl.h"
#include "keta/area_load.h"
#include "keta/statics.h"
#include "keta/phasing.h"
#include "keta/traffic.h"
#include "keta/traffic_fatigue.h"
#include "keta/traffic_multi.h"
#include "keta/traffic_transverse.h"
#include "keta/intensity_law.h"
#include "keta/influence_scatter.h"
#include "keta/rail_fatigue.h"
#include "keta/rail_actions.h"
#include "keta/wind.h"
#include "keta/wind_turbulence.h"
#include "keta/wind_buffeting.h"
#include "keta/wind_equivalent.h"
#include "keta/wind_flutter.h"
#include "keta/wind_synthesis.h"
#include "keta/fire.h"
#include "keta/firechain.h"
#include "keta/springmass.h"
#include "keta/soilzones.h"
#include "keta/warping.h"
#include "keta/hereditary.h"
#include "keta/creep.h"
#include "keta/rc.h"
#include "keta/eig_ext.h"
#include "keta/st1.h"
#include "keta/st1_run.h"
#include "keta/tapered.h"
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include "keta/matlaws.h"
#include "keta/section_exact.h"
#include "keta/resal.h"
#include "keta/thinwall.h"
#include "keta/combos_ec0.h"
#include "keta/envelope.h"
#include "keta/prestress.h"
#include "keta/prestress_tools.h"
#include "keta/ec3.h"

namespace py = pybind11;
using namespace keta;

namespace {

using Dof6 = std::tuple<double, double, double, double, double, double>;
using NodeDof6 = std::pair<Id, Dof6>;
using StationForces = std::tuple<double, double, double, double, double, double, double>;
// The memoryless design laws by name (uniaxial_path, fiber_frame).
std::optional<UniaxialLaw> design_law_of(const std::string& k, const std::function<double(const char*)>& g,
                                         const std::function<bool(const char*)>& has) {
  auto opt = [&](const char* key, double dflt) { return has(key) ? g(key) : dflt; };
  if (k == "sargin")
    return has("fck") ? UniaxialLaw::sargin_ec2(g("fck")) : UniaxialLaw::sargin(g("fcm"), g("Ecm"), g("eps_c1"), g("eps_cu1"));
  if (k == "steel_hardening") return UniaxialLaw::steel_hardening(g("E"), g("fy"), g("k"), g("eps_uk"), g("eps_ud"));
  if (k == "strand")
    return UniaxialLaw::strand(g("Ep"), g("fp01k"), g("fpk"), g("eps_uk"), opt("eps_ud", 0.0), opt("inclined", 1.0) != 0,
                               opt("gamma_s", 1.15));
  if (k == "concrete_tension_softening")
    return UniaxialLaw::concrete_tension_softening(g("fcd"), opt("eps_c2", 2e-3), opt("eps_cu2", 3.5e-3), opt("n", 2.0),
                                                   g("Et"), g("ft"), g("Gf"), g("lc"));
  if (k == "bfup_t3")
    return UniaxialLaw::bfup_t3(g("fck"), g("fcm"), g("Ecm"), g("fctm"), g("K"), g("fctk_el"), g("fctfk"),
                                g("eps_u_lim"), opt("gamma_c", 1.5), opt("gamma_cf", 1.3), opt("alpha_cc", 0.85));
  return std::nullopt;
}

using Vec3T = std::tuple<double, double, double>;

Dof6 dof6(const Vec6& d) { return {d(0), d(1), d(2), d(3), d(4), d(5)}; }

// Shortest round-trip decimal form without exponent: section names.
std::string num(double x) {
  if (std::isfinite(x) && std::trunc(x) == x && std::fabs(x) < 1e15) return std::to_string(static_cast<long long>(x));
  char buf[64];
  auto r = std::to_chars(buf, buf + sizeof buf, x);
  return std::string(buf, r.ptr);
}

std::string fmt(const char* f, double x) {
  char buf[96];
  std::snprintf(buf, sizeof buf, f, x);
  return buf;
}

// Run an engine call, mapping its failures onto ValueError.
template <class F>
auto engine(F&& f) -> decltype(f()) {
  try {
    return f();
  } catch (const py::error_already_set&) {
    throw;
  } catch (const py::key_error&) {
    throw;
  } catch (const py::value_error&) {
    throw;
  } catch (const std::exception& e) {
    throw py::value_error(e.what());
  }
}

std::size_t node_pos(const Model& m, Id node) {
  for (std::size_t i = 0; i < m.nodes.size(); ++i)
    if (m.nodes[i].id == node) return i;
  throw py::key_error("no node " + std::to_string(node));
}

std::size_t check_comp(std::size_t c) {
  if (c >= 6) throw py::value_error("component must be 0..6 = [dx, dy, dz, rx, ry, rz]");
  return c;
}

double max_translation(const std::vector<Vec6>& d) {
  double m = 0;
  for (const Vec6& v : d) m = std::max(m, v.head<3>().norm());
  return m;
}

std::vector<NodeDof6> reactions_of(const CaseResult& c) {
  std::vector<NodeDof6> out;
  for (const Reaction& r : c.reactions) out.emplace_back(r.node, dof6(r.R));
  return out;
}

Vec3T sum_reactions_of(const CaseResult& c) {
  Vec3 s = Vec3::Zero();
  for (const Reaction& r : c.reactions) s += r.R.head<3>();
  return {s.x(), s.y(), s.z()};
}

const BarResult& bar_of(const CaseResult& c, Id bar) {
  for (const BarResult& b : c.bars)
    if (b.bar == bar) return b;
  throw py::key_error("no bar " + std::to_string(bar));
}

// (xi = x / L, N, Vy, Vz, Mt, My, Mz) per station.
std::vector<StationForces> stations_of(const BarResult& b, double L = -1.0) {
  if (L < 0) L = b.stations.empty() ? 0.0 : b.stations.back().x;
  std::vector<StationForces> out;
  for (const Station& s : b.stations) out.emplace_back(L > 0 ? s.x / L : 0.0, s.N, s.Vy, s.Vz, s.T, s.My, s.Mz);
  return out;
}

const Bar& bar_of_model(const Model& m, Id bar) {
  for (const Bar& b : m.bars)
    if (b.id == bar) return b;
  throw py::key_error("no bar " + std::to_string(bar));
}
// Length of a bar between its (eccentric) ends: its stations run along it
// from the origin.
double bar_length(const Model& m, Id bar) {
  const Bar& b = bar_of_model(m, bar);
  return ((m.node(b.nj).x + b.ecc[1]) - (m.node(b.ni).x + b.ecc[0])).norm();
}

std::vector<double> vec12(const Vec12& v) { return std::vector<double>(v.data(), v.data() + 12); }

// Normal stress N/A + My z/Iy - Mz y/Iz of a station at the fibre (y, z); the
// terms of a zero property are skipped.
double normal_stress(const Station& s, const SectionProps& p, double y, double z) {
  double sig = 0;
  if (p.A != 0) sig += s.N / p.A;
  if (p.Iy != 0) sig += s.My * z / p.Iy;
  if (p.Iz != 0) sig -= s.Mz * y / p.Iz;
  return sig;
}

double trapz(const std::vector<double>& x, const std::vector<double>& y) {
  double s = 0;
  for (std::size_t i = 1; i < x.size(); ++i) s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return s;
}

LoadFrame parse_frame(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  if (s == "global") return LoadFrame::Global;
  if (s == "local") return LoadFrame::Local;
  if (s == "projected") return LoadFrame::GlobalProjected;
  throw py::value_error("frame must be 'global' | 'local' | 'projected', got '" + s + "'");
}

DampingSpec parse_damping(std::optional<double> zeta, std::optional<std::vector<double>> zetas,
                          std::optional<std::pair<double, double>> rayleigh, double dflt) {
  const int given = zeta.has_value() + zetas.has_value() + rayleigh.has_value();
  if (given > 1)
    throw py::value_error(
        "give at most one damping kwarg: zeta= (one ratio), zetas= (per mode) or rayleigh=(alpha, beta)");
  DampingSpec d;
  if (zetas) {
    d.kind = DampingSpec::Kind::PerMode;
    d.per_mode = *zetas;
  } else if (rayleigh) {
    d.kind = DampingSpec::Kind::Rayleigh;
    d.alpha = rayleigh->first;
    d.beta = rayleigh->second;
  } else {
    d.kind = DampingSpec::Kind::Modal;
    d.zeta = zeta.value_or(dflt);
  }
  return d;
}

// Every node's free-DOF displacement expanded to 6 components.
std::vector<Vec6> expand(const DofMap& dofs, const VecX& u) {
  std::vector<Vec6> out(static_cast<std::size_t>(dofs.n_nodes), Vec6::Zero());
  for (int ni = 0; ni < dofs.n_nodes; ++ni)
    for (int c = 0; c < 6; ++c) {
      const int f = dofs.free_of(ni, c);
      if (f >= 0) out[static_cast<std::size_t>(ni)](c) = u(f);
    }
  return out;
}

using ModelPtr = std::shared_ptr<const Model>;

// ============================ the builder ============================

struct CableRec {
  Id id = 0, n1 = 0, n2 = 0;
  double ea = 0, w = 0;
  int route = 0;  // 0 = l0, 1 = horizontal tension h, 2 = end tension
  double value = 0;
};

class PyModel {
 public:
  Model m;
  std::vector<CableRec> cables;
  // Global activation (set_active): the bars / support nodes the analyses
  // see; nullopt = all of them.
  std::optional<std::set<Id>> active_bars, active_supports;
  Id next_node = 1, next_bar = 1, next_cable = 1, next_pmass = 1, next_section = 1, next_material = 1, next_combo = 1;

  // --- geometry ---
  Id add_node(double x, double y, double z) {
    Node n;
    n.id = next_node++;
    n.x = Vec3(x, y, z);
    m.nodes.push_back(n);
    m.invalidate_index();
    return n.id;
  }

  // --- materials ---
  Id push_material(Material mat) {
    mat.id = next_material++;
    m.materials.push_back(mat);
    m.invalidate_index();
    return mat.id;
  }
  Id steel() {
    Material s;
    s.name = "Steel";
    s.E = 210e9;
    s.nu = 0.3;
    s.rho = 7850;
    s.alpha = 1.2e-5;
    return push_material(s);
  }
  Id concrete() {
    Material c;
    c.name = "Concrete";
    c.E = 33e9;
    c.nu = 0.2;
    c.rho = 2500;
    c.alpha = 1.0e-5;
    return push_material(c);
  }
  Id add_material(double e, double nu, std::optional<double> rho, std::optional<double> alpha,
                  std::optional<std::string> name, std::optional<double> g, std::optional<double> damping = {}) {
    if (damping && !(*damping >= 0 && *damping < 1)) throw py::value_error("material: damping ratio in [0, 1)");
    Material x;
    x.damping = damping.value_or(0.0);
    x.E = e;
    x.nu = nu;
    x.G_given = g.value_or(0.0);
    x.rho = rho.value_or(0.0);      // no density: no mass, no self-weight
    x.alpha = alpha.value_or(0.0);  // keta's "no alpha"
    x.name = name.value_or("");
    return push_material(x);
  }

  // --- sections ---
  Id push_section(Section s) {
    engine([&] {
      s.finalize();
      return 0;
    });
    s.id = next_section++;
    m.sections.push_back(std::move(s));
    m.invalidate_index();
    return m.sections.back().id;
  }
  Id add_section(const std::string& designation) {
    std::optional<Section> s = catalog_section(designation);
    if (!s) throw py::key_error("unknown profile '" + designation + "'");
    s->name = designation;
    return push_section(*s);
  }
  Id add_section_rect(double b, double h) {
    Section s;
    s.name = "Rect " + num(b) + "x" + num(h);
    s.shape = ShapeRect{b, h};
    return push_section(s);
  }
  Id add_section_tube(double d, double t) {
    Section s;
    s.name = "CHS " + num(d) + "x" + num(t);
    s.shape = ShapeTube{d, t};
    return push_section(s);
  }
  Id add_section_i(double h, double b, double tw, double tf) {
    Section s;
    s.name = "I " + num(h) + "x" + num(b);
    s.shape = ShapeI{h, b, tw, tf, 0.0};
    return push_section(s);
  }
  Id add_section_props(double sx, double iy, double iz, double ix, std::optional<double> sy,
                       std::optional<double> sz, std::optional<std::string> name, double iw = 0.0,
                       double beta_x = 0.0) {
    Section s;
    s.name = name.value_or("");
    s.shape = ShapeGeneric{};
    s.props.A = sx;
    s.props.Iy = iy;
    s.props.Iz = iz;
    s.props.J = ix;
    s.props.Avy = sy.value_or(0.0);  // none: rigid in shear
    s.props.Avz = sz.value_or(0.0);
    s.props.Iw = iw;
    s.props.beta_x = beta_x;
    s.id = next_section++;
    m.sections.push_back(std::move(s));
    m.invalidate_index();
    return m.sections.back().id;
  }

  // --- elements ---
  Id push_bar(Id n1, Id n2, Id section, Id material, BarKind kind) {
    Bar b;
    b.id = next_bar++;
    b.ni = n1;
    b.nj = n2;
    b.section = section;
    b.material = material;
    b.kind = kind;
    m.bars.push_back(b);
    m.resolve_shear(m.bars.back());
    m.invalidate_index();
    return b.id;
  }
  Id add_beam(Id n1, Id n2, Id section, Id material, double beta, std::optional<Vec3T> ecc,
               std::optional<Vec3T> ecc_end) {
    const Id id = push_bar(n1, n2, section, material, BarKind::Frame);
    Bar& b = m.bars.back();
    b.roll = beta;
    if (ecc) {
      const Vec3 e(std::get<0>(*ecc), std::get<1>(*ecc), std::get<2>(*ecc));
      b.ecc = {e, e};
    }
    if (ecc_end) b.ecc[1] = Vec3(std::get<0>(*ecc_end), std::get<1>(*ecc_end), std::get<2>(*ecc_end));
    return id;
  }
  Id add_truss(Id n1, Id n2, Id section, Id material) { return push_bar(n1, n2, section, material, BarKind::Truss); }
  Id add_cable(Id n1, Id n2, double ea, double w, std::optional<double> l0, std::optional<double> h,
               std::optional<double> end_tension) {
    const int given = l0.has_value() + h.has_value() + end_tension.has_value();
    if (given == 0) throw py::value_error("add_cable needs exactly one of l0=, h= or end_tension= (got none)");
    if (given > 1)
      throw py::value_error("add_cable takes exactly one of l0=, h= or end_tension= (got more than one)");
    CableRec c;
    c.id = next_cable++;
    c.n1 = n1;
    c.n2 = n2;
    c.ea = ea;
    c.w = w;
    c.route = l0 ? 0 : h ? 1 : 2;
    c.value = l0 ? *l0 : h ? *h : *end_tension;
    cables.push_back(c);
    return c.id;
  }
  Id add_continuum(const std::string& kind, const std::vector<Id>& nodes, Id material, double t) {
    Continuum::Kind k = Continuum::Kind::PlaneStress;
    if (kind == "plane_strain") k = Continuum::Kind::PlaneStrain;
    else if (kind == "axisymmetric") k = Continuum::Kind::Axisymmetric;
    else if (kind == "solid") k = Continuum::Kind::Solid;
    else if (kind != "plane_stress")
      throw py::value_error("add_continuum: kind is 'plane_stress', 'plane_strain', 'axisymmetric' or 'solid'");
    return engine([&] { return m.add_continuum(k, nodes, material, t); });
  }
  Id add_shell(Id n1, Id n2, Id n3, double t, Id material, Id slab, double offset, std::optional<std::array<double, 3>> t_corners,
               double bed, std::optional<py::dict> ortho, std::optional<py::dict> laminate) {
    ShellTri sh;
    Id mx = 0;
    for (const ShellTri& x : m.shells) mx = std::max(mx, x.id);
    sh.id = mx + 1;
    sh.n = {n1, n2, n3};
    sh.t = t;
    sh.material = material;
    sh.slab = slab;
    sh.enrich.ecc = offset;
    sh.enrich.kbed = bed;
    if (t_corners) sh.enrich.tn = *t_corners;
    auto num = [](const py::dict& d, const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
    if (ortho) {
      ShellOrtho o;
      o.Ex = num(*ortho, "Ex", 0);
      o.Ey = num(*ortho, "Ey", 0);
      o.nu_xy = num(*ortho, "nu_xy", 0);
      o.Gxy = num(*ortho, "Gxy", 0);
      o.Gxz = num(*ortho, "Gxz", 0);
      o.Gyz = num(*ortho, "Gyz", 0);
      o.angle = num(*ortho, "angle", 0);
      sh.enrich.ortho = o;
    }
    if (laminate) {
      ShellLaminate lm;
      for (const py::handle& h : (*laminate)["plies"].cast<py::list>()) {
        const py::dict d = py::reinterpret_borrow<py::dict>(h);
        ShellPly p;
        p.t = num(d, "t", 0);
        p.E1 = num(d, "E1", 0);
        p.E2 = num(d, "E2", 0);
        p.nu12 = num(d, "nu12", 0);
        p.G12 = num(d, "G12", 0);
        p.G13 = num(d, "G13", 0);
        p.G23 = num(d, "G23", 0);
        p.angle = num(d, "angle", 0);
        lm.plies.push_back(p);
      }
      if (laminate->contains("z_bottom")) {
        lm.centred = false;
        lm.z_bottom = (*laminate)["z_bottom"].cast<double>();
      }
      sh.enrich.laminate = lm;
    }
    (void)m.node(n1);
    (void)m.node(n2);
    (void)m.node(n3);
    (void)m.material(material);
    if (!(t > 0)) throw py::value_error("add_shell: thickness must be > 0");
    m.shells.push_back(sh);
    return sh.id;
  }
  void add_surface_pressure(Id slab, double qx, double qy, double qz, const std::string& frame, Id cs) {
    SurfaceFrame f = SurfaceFrame::Global;
    if (frame == "normal") f = SurfaceFrame::LocalNormal;
    else if (frame == "projected") f = SurfaceFrame::Projected;
    else if (frame != "global") throw py::value_error("add_surface_pressure: frame 'global', 'normal' or 'projected'");
    case_mut(cs).loads.push_back(Load::surface_pressure(slab, Vec3(qx, qy, qz), f));
  }
  void add_continuum_traction(Id element, int face, std::array<double, 3> tr, Id cs) {
    for (const Load& l : engine([&] { return continuum_face_load(m, element, face, Vec3(tr[0], tr[1], tr[2])); }))
      case_mut(cs).loads.push_back(l);
  }
  void add_continuum_pressure(Id element, int face, double p, Id cs) {
    for (const Load& l : engine([&] { return continuum_face_pressure(m, element, face, p); })) case_mut(cs).loads.push_back(l);
  }
  void add_continuum_thermal(Id element, const std::map<Id, double>& T, double T0, Id cs) {
    for (const Load& l : engine([&] { return continuum_thermal_load(m, element, T, T0); })) case_mut(cs).loads.push_back(l);
  }
  Id add_point_mass(Id node, double mass, std::optional<std::array<double, 3>> inertia) {
    PointMass p;
    p.node = node;
    p.mass = mass;
    if (inertia) p.inertia = Vec3((*inertia)[0], (*inertia)[1], (*inertia)[2]);
    m.point_masses.push_back(p);
    return next_pmass++;
  }
  Bar& bar(Id id) {
    for (Bar& b : m.bars)
      if (b.id == id) return b;
    throw py::key_error("no bar " + std::to_string(id));
  }
  void set_releases(Id id, bool ry_or, bool rz_or, bool ry_ex, bool rz_ex, const std::string& origin,
                    const std::string& end) {
    Bar& b = bar(id);
    b.releases = {};
    b.releases[4] = ry_or;
    b.releases[5] = rz_or;
    b.releases[10] = ry_ex;
    b.releases[11] = rz_ex;
    static const char* names[6] = {"dx", "dy", "dz", "rx", "ry", "rz"};
    for (int e = 0; e < 2; ++e) {
      std::istringstream in(e == 0 ? origin : end);
      std::string w;
      while (in >> w) {
        for (char& c : w) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        int k = -1;
        for (int i = 0; i < 6; ++i)
          if (w == names[i]) k = i;
        if (k < 0) throw py::value_error("set_releases: unknown component '" + w + "' (dx dy dz rx ry rz)");
        b.releases[static_cast<std::size_t>(6 * e + k)] = true;
      }
    }
  }
  void set_soil(Id id, double kfx, double kfy, double kfz, double kmx) {
    Bar& b = bar(id);
    b.soil = SoilProps{kfx, kfy, kfz, kmx};
  }

  // --- supports ---
  void fix(Id node) { add_support(node, true, true, true, true, true, true, {}, {}, {}, {}, {}, {}); }
  void pin(Id node) { add_support(node, true, true, true, false, false, false, {}, {}, {}, {}, {}, {}); }
  void add_support(Id node, bool dx, bool dy, bool dz, bool rx, bool ry, bool rz, std::optional<double> kx,
                   std::optional<double> ky, std::optional<double> kz, std::optional<double> krx,
                   std::optional<double> kry, std::optional<double> krz) {
    Support s;
    s.node = node;
    s.fixed = {dx, dy, dz, rx, ry, rz};
    const std::array<std::optional<double>, 6> k{kx, ky, kz, krx, kry, krz};
    for (int i = 0; i < 6; ++i)
      if (k[static_cast<std::size_t>(i)] && !s.fixed[static_cast<std::size_t>(i)])
        s.spring[static_cast<std::size_t>(i)] = *k[static_cast<std::size_t>(i)];  // a blocked DOF stays rigid
    m.supports.push_back(s);
  }

  // --- general constraints and skew supports ---
  Id add_constraint(const std::vector<std::tuple<Id, int, double>>& terms) {
    std::vector<ConstraintTerm> t;
    for (const auto& [n, c, k] : terms) t.push_back({n, c, k});
    return m.add_constraint(t);
  }
  Id add_skew_support(Id node, const std::array<std::array<double, 3>, 3>& axes, bool dx, bool dy, bool dz, bool rx,
                      bool ry, bool rz) {
    Mat3 A;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) A(i, j) = axes[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    return m.add_skew_support(node, A, {dx, dy, dz, rx, ry, rz});
  }

  // --- node-to-node links ---
  Id add_rigid_link(Id master, const std::vector<Id>& slaves, bool dx, bool dy, bool dz, bool rx, bool ry,
                    bool rz) {
    return m.add_rigid_link(master, slaves, {dx, dy, dz, rx, ry, rz});
  }
  Id add_elastic_link(Id ni, Id nj, double kx, double ky, double kz, double krx, double kry, double krz,
                      double s, double roll) {
    const Id id = m.add_elastic_link(ni, nj, ElasticLink::diagonal(kx, ky, kz, krx, kry, krz), s);
    m.links.back().roll = roll;
    return id;
  }
  void set_link_law(Id link, int comp, const std::string& kind, double k, double param, double kp,
                    std::optional<std::vector<std::pair<double, double>>> points) {
    if (comp < 0 || comp > 5) throw py::value_error("link law: comp must be 0..5");
    for (ElasticLink& l : m.links)
      if (l.id == link) {
        LinkLaw law;
        if (kind == "tension")
          law = LinkLaw::tension_only(k);
        else if (kind == "compression")
          law = LinkLaw::compression_only(k);
        else if (kind == "gap")
          law = LinkLaw::gap(param, k);
        else if (kind == "bilinear")
          law = LinkLaw::bilinear(k, param, kp);
        else if (kind == "points") {
          if (!points) throw py::value_error("link law 'points': give points=[(d, F), ...]");
          law.pts = *points;
        } else {
          throw py::value_error("link law: kind 'tension', 'compression', 'gap', 'bilinear' or 'points', got '" +
                                kind + "'");
        }
        l.law[static_cast<std::size_t>(comp)] = law;
        l.k(comp, comp) = 0;
        return;
      }
    throw py::key_error("no elastic link " + std::to_string(link));
  }

  // --- tuning (keta/cable_tuning.h) ---
  Id add_unit_prestress_case(Id bar, double P) {
    return engine([&] { return keta::add_unit_prestress_case(m, bar, P); });
  }

  // --- soil springs (keta/soil_springs.h) ---
  Id add_ground_spring(Id node, Vec3T dir, const std::vector<std::pair<double, double>>& points) {
    LinkLaw law;
    law.pts = points;
    return soil::add_ground_spring(m, node, Vec3(std::get<0>(dir), std::get<1>(dir), std::get<2>(dir)), law);
  }
  std::vector<Id> add_earth_pressure(const std::vector<Id>& wall, Vec3T away, Vec3T up, const py::list& layers,
                                     double width, const std::string& distribution) {
    std::vector<soil::EarthLayer> ls;
    for (const py::handle& h : layers) {
      const py::dict d = py::reinterpret_borrow<py::dict>(h);
      auto get = [&](const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
      soil::EarthLayer l;
      l.top = get("top", 0);
      l.bottom = get("bottom", 0);
      l.sv_top = get("sv_top", 0);
      l.gamma = get("gamma", 0);
      l.Ka = get("Ka", 0);
      l.K0 = get("K0", 0);
      l.Kp = get("Kp", 0);
      l.c = get("c", 0);
      l.Kac = get("Kac", 0);
      l.Kpc = get("Kpc", 0);
      l.kh_A = get("kh_A", 0);
      l.kh_B = get("kh_B", 0);
      l.kh_kappa = get("kh_kappa", 1);
      ls.push_back(l);
    }
    soil::EarthDistribution dist;
    if (distribution == "consistent")
      dist = soil::EarthDistribution::Consistent;
    else if (distribution == "lumped")
      dist = soil::EarthDistribution::Lumped;
    else
      throw py::value_error("add_earth_pressure: distribution 'consistent' or 'lumped', got '" + distribution + "'");
    auto v = [](Vec3T t) { return Vec3(std::get<0>(t), std::get<1>(t), std::get<2>(t)); };
    return engine([&] { return soil::add_earth_pressure(m, wall, v(away), v(up), ls, width, dist); });
  }

  // --- loads ---
  LoadCase& case_mut(Id id) {
    for (LoadCase& c : m.cases)
      if (c.id == id) return c;
    LoadCase c;
    c.id = id;
    c.name = "LC" + std::to_string(id);
    m.cases.push_back(c);
    std::sort(m.cases.begin(), m.cases.end(), [](const LoadCase& a, const LoadCase& b) { return a.id < b.id; });
    m.invalidate_index();
    for (LoadCase& x : m.cases)
      if (x.id == id) return x;
    throw std::logic_error("case insertion failed");
  }
  void add_point_load(Id node, double fx, double fy, double fz, double mx, double my, double mz, Id cs) {
    case_mut(cs).loads.push_back(Load::nodal(node, Vec3(fx, fy, fz), Vec3(mx, my, mz)));
  }
  static std::vector<Vec3> pts3(const std::vector<std::tuple<double, double, double>>& p) {
    std::vector<Vec3> o;
    for (const auto& [x, y, z] : p) o.emplace_back(x, y, z);
    return o;
  }
  int add_area_load(const std::vector<std::tuple<double, double, double>>& poly, double qx, double qy, double qz,
                    const std::string& mode, std::tuple<double, double, double> span, Id cs) {
    AreaLoadZone z;
    z.poly = pts3(poly);
    z.q = Vec3(qx, qy, qz);
    if (mode == "one_way")
      z.mode = AreaLoadZone::Mode::OneWay;
    else if (mode != "two_way")
      throw py::value_error("area load: mode 'two_way' or 'one_way', got '" + mode + "'");
    z.span = Vec3(std::get<0>(span), std::get<1>(span), std::get<2>(span));
    const std::vector<Load> ls = engine([&] { return area_load_to_bars(m, z); });
    for (const Load& l : ls) case_mut(cs).loads.push_back(l);
    return static_cast<int>(ls.size());
  }
  Id add_shell_tendon(const std::vector<std::tuple<double, double, double>>& pts, double P0, double mu, double k,
                      const std::string& jack, Id slab) {
    ShellTendon t;
    t.id = m.shell_tendons.empty() ? 1 : m.shell_tendons.back().id + 1;
    t.pts = pts3(pts);
    t.P0 = P0;
    t.mu = mu;
    t.k = k;
    if (jack == "end")
      t.jack = ShellTendon::Jack::End;
    else if (jack == "both")
      t.jack = ShellTendon::Jack::Both;
    else if (jack != "start")
      throw py::value_error("shell tendon: jack 'start', 'end' or 'both', got '" + jack + "'");
    t.slab = slab;
    m.shell_tendons.push_back(t);
    return t.id;
  }
  void add_shell_prestress(Id tendon, double factor, Id cs) {
    case_mut(cs).loads.push_back(Load::shell_prestress(tendon, factor));
  }
  void add_polygon_pressure(const std::vector<std::tuple<double, double, double>>& poly, double qx, double qy,
                            double qz, const std::vector<std::vector<std::tuple<double, double, double>>>& holes,
                            Id slab, Id cs) {
    std::vector<std::vector<Vec3>> h;
    for (const auto& r : holes) h.push_back(pts3(r));
    case_mut(cs).loads.push_back(Load::polygon_pressure(pts3(poly), Vec3(qx, qy, qz), h, slab));
  }
  void add_hydrostatic(double gamma, std::array<double, 3> point, Id slab, double side, std::array<double, 3> down, Id cs) {
    case_mut(cs).loads.push_back(
        Load::hydrostatic(gamma, Vec3(point[0], point[1], point[2]), slab, side, Vec3(down[0], down[1], down[2])));
  }
  void add_udl(Id b, double wx, double wy, double wz, Id cs, const std::string& frame, double mx, double my,
               double mz, std::optional<double> x1, std::optional<double> x2, bool absolute) {
    const LoadFrame fr = parse_frame(frame);
    (void)bar(b);
    Load l = Load::udl(b, Vec3(wx, wy, wz), fr);
    l.M = Vec3(mx, my, mz);
    if (x1 || x2) {
      l.a = x1.value_or(0.0);
      l.b = x2.value_or(absolute ? -1.0 : 1.0);
      l.relative = !absolute;
    }
    case_mut(cs).loads.push_back(l);
  }
  void add_linear_load(Id b, Vec3T w1, Vec3T w2, double x1, double x2, Id cs, const std::string& frame,
                       std::optional<Vec3T> m1, std::optional<Vec3T> m2, bool absolute) {
    const LoadFrame fr = parse_frame(frame);
    (void)bar(b);
    Load l = Load::linear(b, Vec3(std::get<0>(w1), std::get<1>(w1), std::get<2>(w1)),
                          Vec3(std::get<0>(w2), std::get<1>(w2), std::get<2>(w2)), fr, x1, x2);
    if (m1) l.M = Vec3(std::get<0>(*m1), std::get<1>(*m1), std::get<2>(*m1));
    if (m2) l.M2 = Vec3(std::get<0>(*m2), std::get<1>(*m2), std::get<2>(*m2));
    l.relative = !absolute;
    case_mut(cs).loads.push_back(l);
  }
  void add_parabolic_load(Id b, Vec3T w1, Vec3T wmid, Vec3T w2, double x1, double x2, Id cs,
                          const std::string& frame, bool absolute) {
    const LoadFrame fr = parse_frame(frame);
    (void)bar(b);
    auto v = [](const Vec3T& w) { return Vec3(std::get<0>(w), std::get<1>(w), std::get<2>(w)); };
    Load l = Load::parabolic(b, v(w1), v(wmid), v(w2), fr, x1, x2);
    l.relative = !absolute;
    case_mut(cs).loads.push_back(l);
  }
  void add_bar_point_load(Id b, double at, double fx, double fy, double fz, Id cs, const std::string& frame,
                          double mx, double my, double mz, bool absolute) {
    const LoadFrame fr = parse_frame(frame);
    (void)bar(b);
    Load l = Load::point(b, at, Vec3(fx, fy, fz), Vec3(mx, my, mz), fr);
    l.relative = !absolute;
    case_mut(cs).loads.push_back(l);
  }
  void add_thermal(Id b, double dt, double gy, double gz, Id cs) {
    case_mut(cs).loads.push_back(Load::thermal(b, dt, gy, gz));
  }
  void add_support_disp(Id node, std::optional<double> dx, std::optional<double> dy, std::optional<double> dz,
                        std::optional<double> rx, std::optional<double> ry, std::optional<double> rz, Id cs) {
    case_mut(cs).loads.push_back(Load::support_disp(node, {dx, dy, dz, rx, ry, rz}));
  }
  void add_self_weight(Id cs, double factor, std::optional<Vec3T> g) {
    const Vec3 v = g ? Vec3(std::get<0>(*g), std::get<1>(*g), std::get<2>(*g)) * factor
                     : Vec3(0, 0, -kGravity * factor);
    case_mut(cs).loads.push_back(Load::self_weight(v));
  }
  void name_case(Id cs, const std::string& name) { case_mut(cs).name = name; }
  Id add_combination(const std::vector<std::pair<double, Id>>& terms, std::optional<std::string> name) {
    Combination c;
    c.id = next_combo++;
    c.name = name.value_or("");
    for (const auto& [f, id] : terms) c.terms.emplace_back(id, f);
    m.combos.push_back(c);
    m.invalidate_index();
    return c.id;
  }

  // The model the solvers see: the frame plus the cables as Cable-kind bars
  // (ids after the largest bar id) on one inert placeholder material and
  // section — Assemble's lowering.
  ModelPtr solver_model() const {
    auto out = std::make_shared<Model>(m);
    if (active_bars) {
      std::erase_if(out->bars, [&](const Bar& b) { return !active_bars->count(b.id); });
      std::erase_if(out->bar_masses, [&](const BarMass& b) { return !active_bars->count(b.bar); });
      for (LoadCase& c : out->cases)
        std::erase_if(c.loads, [&](const Load& l) { return l.bar != 0 && !active_bars->count(l.bar); });
    }
    if (active_supports)
      std::erase_if(out->supports, [&](const Support& s) { return !active_supports->count(s.node); });
    if (!cables.empty()) {
      Id max_mat = 0, max_sec = 0, max_bar = 0;
      for (const Material& x : out->materials) max_mat = std::max(max_mat, x.id);
      for (const Section& x : out->sections) max_sec = std::max(max_sec, x.id);
      for (const Bar& x : out->bars) max_bar = std::max(max_bar, x.id);
      Material pm;
      pm.id = max_mat + 1;
      pm.name = "cable (placeholder)";
      pm.E = 1.0;
      pm.nu = 0.3;
      pm.rho = 0.0;
      out->materials.push_back(pm);
      Section ps;
      ps.id = max_sec + 1;
      ps.name = "cable (placeholder)";
      ps.shape = ShapeRect{1.0, 1.0};
      ps.finalize();
      out->sections.push_back(ps);
      for (std::size_t k = 0; k < cables.size(); ++k) {
        const CableRec& c = cables[k];
        Bar b;
        b.id = max_bar + 1 + static_cast<Id>(k);
        b.ni = c.n1;
        b.nj = c.n2;
        b.material = pm.id;
        b.section = ps.id;
        b.kind = BarKind::Cable;
        CableProps cp;
        cp.ea = c.ea;
        cp.w = c.w;
        if (c.route == 0) cp.unstressed_length = c.value;
        if (c.route == 1) cp.pretension = c.value;
        if (c.route == 2) cp.end_tension = c.value;
        b.cable = cp;
        out->bars.push_back(b);
        out->resolve_shear(out->bars.back());
      }
    }
    out->invalidate_index();
    return out;
  }
  Id cable_bar_id(Id cable) const {
    Id max_bar = 0;
    for (const Bar& x : m.bars) max_bar = std::max(max_bar, x.id);
    for (std::size_t k = 0; k < cables.size(); ++k)
      if (cables[k].id == cable) return max_bar + 1 + static_cast<Id>(k);
    throw py::key_error("no cable " + std::to_string(cable));
  }

  // Wrap a bare keta Model as a builder, resuming every auto id past the
  // largest present one; Cable-kind bars return to the cable list and the
  // placeholder material / section are dropped.
  static PyModel from_model(const Model& in) {
    PyModel p;
    p.m = in;
    Id ph_mat = -1, ph_sec = -1;
    std::vector<Bar> frame;
    for (const Bar& b : in.bars) {
      if (b.kind != BarKind::Cable || !b.cable) {
        frame.push_back(b);
        continue;
      }
      CableRec c;
      c.n1 = b.ni;
      c.n2 = b.nj;
      c.ea = b.cable->ea;
      c.w = b.cable->w.value_or(0.0);
      if (b.cable->unstressed_length > 0) {
        c.route = 0;
        c.value = b.cable->unstressed_length;
      } else if (b.cable->end_tension > 0) {
        c.route = 2;
        c.value = b.cable->end_tension;
      } else {
        c.route = 1;
        c.value = b.cable->pretension;
      }
      c.id = static_cast<Id>(p.cables.size()) + 1;
      p.cables.push_back(c);
      ph_mat = b.material;
      ph_sec = b.section;
    }
    p.m.bars = frame;
    std::erase_if(p.m.materials, [&](const Material& x) { return x.id == ph_mat && x.name == "cable (placeholder)"; });
    std::erase_if(p.m.sections, [&](const Section& x) { return x.id == ph_sec && x.name == "cable (placeholder)"; });
    auto next = [](auto first, auto last, auto key) {
      Id mx = 0;
      for (auto it = first; it != last; ++it) mx = std::max(mx, key(*it));
      return mx + 1;
    };
    p.next_node = next(p.m.nodes.begin(), p.m.nodes.end(), [](const Node& x) { return x.id; });
    p.next_bar = next(p.m.bars.begin(), p.m.bars.end(), [](const Bar& x) { return x.id; });
    p.next_cable = static_cast<Id>(p.cables.size()) + 1;
    p.next_pmass = static_cast<Id>(p.m.point_masses.size()) + 1;
    p.next_section = next(p.m.sections.begin(), p.m.sections.end(), [](const Section& x) { return x.id; });
    p.next_material = next(p.m.materials.begin(), p.m.materials.end(), [](const Material& x) { return x.id; });
    p.next_combo = next(p.m.combos.begin(), p.m.combos.end(), [](const Combination& x) { return x.id; });
    p.m.invalidate_index();
    return p;
  }

  Id case_or_first(std::optional<Id> c) const {
    if (c) return *c;
    if (m.cases.empty()) throw py::value_error("the model has no load cases");
    return m.cases.front().id;
  }

  std::vector<std::string> validate() const {
    try {
      solver_model()->validate();
    } catch (const std::exception& e) {
      return {e.what()};
    }
    return {};
  }

  std::string repr() const {
    return "<keta.Model: " + std::to_string(m.nodes.size()) + " nodes, " + std::to_string(m.bars.size()) +
           " bars, " + std::to_string(m.cases.size()) + " cases, " + std::to_string(m.combos.size()) + " combos>";
  }
};

// ============================ results ============================

struct PyResults {
  ModelPtr model;
  std::shared_ptr<const StaticResults> r;
};

struct PyCaseResult {
  ModelPtr model;
  std::shared_ptr<const StaticResults> r;
  std::size_t idx = 0;
  bool combo = false;
  const CaseResult& get() const { return combo ? r->combos[idx] : r->cases[idx]; }
};

// A single-case result (second order, cable statics) shares the accessors.
struct SingleCase {
  ModelPtr model;
  std::shared_ptr<const CaseResult> c;
};

template <class Cls, class Get>
void bind_case_accessors(Cls& cls, Get get, const char* end_forces_name) {
  cls.def("displacement",
          [get](const typename Cls::type& self, Id node) {
            const auto& [model, c] = get(self);
            return dof6(c.node_disp.at(node_pos(*model, node)));
          })
      .def("reactions", [get](const typename Cls::type& self) { return reactions_of(get(self).second); })
      .def("constraint_forces", [get](const typename Cls::type& self) { return get(self).second.constraint_forces; })
      .def("continuum_stress",
           [get](const typename Cls::type& self, Id element, std::array<double, 3> xi) {
             const auto& [model, c] = get(self);
             for (const Continuum& e : model->continua)
               if (e.id == element) {
                 const VecX st = engine([&] { return continuum_stress(*model, e, c.node_disp, xi); });
                 return std::vector<double>(st.data(), st.data() + st.size());
               }
             throw py::key_error("no continuum element " + std::to_string(element));
           },
           py::arg("element"), py::arg("xi") = std::array<double, 3>{0, 0, 0},
           "Stress of a continuum element at the natural point xi: plane (sx, sy, txy, sz), axisymmetric (sr, sz, "
           "stheta, trz), solid (sx, sy, sz, txy, tyz, tzx).")
      .def("skew_reactions",
           [get](const typename Cls::type& self) {
             std::vector<std::tuple<Id, Dof6, Dof6>> out;
             for (const auto& sr : get(self).second.skew_reactions)
               out.emplace_back(sr.node, dof6(sr.local), dof6(sr.global));
             return out;
           })
      .def("bar_forces",
           [get](const typename Cls::type& self, Id bar) {
             const auto& [model, c] = get(self);
             return stations_of(bar_of(c, bar), bar_length(*model, bar));
           })
      .def("bar_deflections",
           [get](const typename Cls::type& self, Id bar) {
             const auto& [model, c] = get(self);
             const BarResult& b = bar_of(c, bar);
             const double L = bar_length(*model, bar);
             std::vector<std::array<double, 7>> out;
             for (const Station& s : b.stations)
               out.push_back({L > 0 ? s.x / L : 0.0, s.defl(0), s.defl(1), s.defl(2), s.rot(0), s.rot(1), s.rot(2)});
             return out;
           },
           py::arg("bar"),
           "Displacements at the study sections, local frame: [(x/L, u, v, w, theta_x, theta_y, theta_z)] [m, rad] "
           "(theta_y = -dw/dx).")
      .def("soil_pressure",
           [get](const typename Cls::type& self, Id bar) {
             const auto& [model, c] = get(self);
             const BarResult& b = bar_of(c, bar);
             if (b.soil_p.empty()) throw py::value_error("bar " + std::to_string(bar) + " has no elastic foundation");
             const double L = bar_length(*model, bar);
             std::vector<std::array<double, 5>> out;
             for (std::size_t k = 0; k < b.stations.size() && k < b.soil_p.size(); ++k)
               out.push_back({L > 0 ? b.stations[k].x / L : 0.0, b.soil_p[k](0), b.soil_p[k](1), b.soil_p[k](2),
                              k < b.soil_m.size() ? b.soil_m[k] : 0.0});
             return out;
           },
           py::arg("bar"),
           "Soil reaction at the study sections of a bar on an elastic foundation (set_soil), local frame: "
           "[(x/L, px, py, pz, mx)] [N/m, N m/m], p = k u.")
      .def("stresses",
           [get](const typename Cls::type& self, Id bar) {
             const auto& [model, c] = get(self);
             const BarResult& b = bar_of(c, bar);
             const double L = bar_length(*model, bar);
             const Bar& br = bar_of_model(*model, bar);
             std::vector<std::tuple<double, std::string, double, double, double>> out;
             for (const Station& s : b.stations)
               for (const StressAtPoint& q : engine([&] { return station_stresses(*model, br, s); }))
                 out.emplace_back(L > 0 ? s.x / L : 0.0, q.point, q.sigma, q.tau, q.von_mises);
             return out;
           },
           py::arg("bar"),
           "Stresses at the named stress points of the bar's section, at every study section: [(x/L, point, sigma, "
           "tau, von_mises)] [Pa]; sigma = N/A + My z/Iy - Mz y/Iz, tau from the shear forces over the shear areas.")
      .def("shell_forces",
           [get](const typename Cls::type& self, Id shell) {
             const auto& [model, c] = get(self);
             for (const ShellResult& s : c.shells)
               if (s.shell == shell) {
                 const ShellResultants& r = s.r;
                 py::dict d;
                 d["n"] = std::array<double, 3>{r.nx, r.ny, r.nxy};
                 d["m"] = std::array<double, 3>{r.mx, r.my, r.mxy};
                 d["q"] = std::array<double, 2>{r.qx, r.qy};
                 d["n_principal"] = std::array<double, 3>{r.n1, r.n2, r.n_theta};
                 d["m_principal"] = std::array<double, 3>{r.m1, r.m2, r.m_theta};
                 return d;
               }
             throw py::key_error("no shell " + std::to_string(shell) + " in this result");
           },
           py::arg("shell"),
           "Resultants of a facet shell at its centroid, local frame: n = (nx, ny, nxy) [N/m], m = (mx, my, mxy) "
           "[N m/m], q = (qx, qy) [N/m], principal values with their angle [rad].")
      .def(end_forces_name,
           [get](const typename Cls::type& self, Id bar) { return vec12(bar_of(get(self).second, bar).end_forces_local); })
      .def("link_force",
           [get](const typename Cls::type& self, Id link) {
             for (const LinkResult& l : get(self).second.links)
               if (l.link == link) return dof6(l.force);
             throw py::index_error("no elastic link " + std::to_string(link));
           })
      .def("max_displacement", [get](const typename Cls::type& self) { return max_translation(get(self).second.node_disp); })
      .def("sum_reactions", [get](const typename Cls::type& self) { return sum_reactions_of(get(self).second); });
}

struct PyLinkLaw {
  ModelPtr model;
  std::shared_ptr<const LinkLawResult> r;
};
struct PyBuckling {
  ModelPtr model;
  std::shared_ptr<const BucklingResults> r;
};
struct PySecondOrder {
  ModelPtr model;
  std::shared_ptr<const SecondOrderResults> r;
};
// Dense index of a load case given its number.
std::size_t case_index_of(const Model& m, Id id) {
  for (std::size_t i = 0; i < m.cases.size(); ++i)
    if (m.cases[i].id == id) return i;
  throw py::key_error("no load case " + std::to_string(id));
}

struct PyNonlinear {
  ModelPtr model;
  std::shared_ptr<const CorotResults> r;
  int substeps = 1;
  std::vector<std::pair<Id, Id>> cable_bars;  // (Python cable id, solver bar id)
  Id case_id = 0;                             // appended: the case solved (tangent())
};
struct PyArcTrace {
  ModelPtr model;
  std::shared_ptr<const ArcTrace> r;
  std::shared_ptr<const DofMap> dofs;
};
struct PyModal {
  ModelPtr model;
  std::shared_ptr<const ModalResults> r;
};
struct PyCableStatic {
  ModelPtr model;
  std::shared_ptr<const CableStaticResults> r;
  std::vector<std::pair<Id, std::size_t>> cables;  // (cable id, model bar index)
  std::vector<double> l0;                          // per model bar
  Id case_id = 0;                                  // appended: the case solved (tangent())
};
struct PyHarmonic {
  ModelPtr model;
  std::shared_ptr<const HarmonicResults> r;
  std::vector<Id> bar_ids;
};
struct PyPsd {
  ModelPtr model;
  std::shared_ptr<const PsdResults> r;
  std::shared_ptr<const HarmonicResults> h;  // the unit-pattern FRF behind the stress PSDs
  std::vector<double> input;                 // S_alpha at each sweep frequency
  std::vector<Id> bar_ids;
};
struct PySpectrum {
  ModelPtr model;
  std::shared_ptr<const SpectrumResults> r;
  int dir = 0;
  double total_mass = 0;
};

std::vector<Id> node_ids_of(const Model& m) {
  std::vector<Id> v;
  for (const Node& n : m.nodes) v.push_back(n.id);
  return v;
}
std::vector<Id> bar_ids_of(const Model& m) {
  std::vector<Id> v;
  for (const Bar& b : m.bars) v.push_back(b.id);
  return v;
}

std::size_t bar_pos(const std::vector<Id>& ids, Id bar) {
  for (std::size_t i = 0; i < ids.size(); ++i)
    if (ids[i] == bar) return i;
  throw py::key_error("no bar " + std::to_string(bar));
}

std::size_t check_fi(std::size_t fi, std::size_t n) {
  if (fi >= n) throw py::key_error("no sweep frequency index " + std::to_string(fi));
  return fi;
}

const SectionProps& props_of(const Model& m, Id bar) {
  const Bar& b = m.bar(bar);
  return m.section(b.section).props;
}

double stress_amp(const HarmonicResults& h, const Model& m, const std::vector<Id>& bar_ids, std::size_t fi, Id bar,
                  std::size_t station, double y, double z) {
  const std::size_t bi = bar_pos(bar_ids, bar);
  const auto& re = h.bars_re.at(fi).at(bi).stations;
  if (station >= re.size()) throw py::key_error("no station " + std::to_string(station) + " on bar " + std::to_string(bar));
  const SectionProps& p = props_of(m, bar);
  const double a = normal_stress(re[station], p, y, z);
  const double b = normal_stress(h.bars_im[fi][bi].stations[station], p, y, z);
  return std::hypot(a, b);
}

// The mode shape list of a buckling / modal result.
std::vector<NodeDof6> mode_rows(const Model& m, const std::vector<std::vector<Vec6>>& modes, std::size_t k,
                                const char* what) {
  if (k >= modes.size()) throw py::key_error(std::string("no ") + what + std::to_string(k));
  std::vector<NodeDof6> out;
  for (std::size_t i = 0; i < m.nodes.size(); ++i) out.emplace_back(m.nodes[i].id, dof6(modes[k][i]));
  return out;
}

#include "py_phasing.inl"
#include "py_traffic.inl"
#include "py_wind.inl"
#include "py_fire.inl"
#include "py_st1.inl"
#include "py_stability.inl"
#include "py_creep.inl"
#include "py_rc.inl"
#include "py_dynamics.inl"
#include "py_modeling.inl"
#include "py_cables.inl"

// keta.uniaxial_path. A named function, like fiber_frame_py: MSVC 14.51 stops with an internal
// compiler error on the lambda form.
std::vector<std::pair<double, double>> uniaxial_path_py(std::string law, const py::dict& p,
                                                        std::vector<double> path, bool stress_driven) {
  auto g = [&](const char* k) {
    if (!p.contains(k)) throw py::value_error(std::string("uniaxial_path: missing parameter ") + k);
    return p[k].cast<double>();
  };
  HystLaw l;
  if (law == "steel_mixed")
    l = HystLaw::steel_mixed(g("E"), g("fy"), g("H"), g("iso"));
  else if (law == "asym_von_mises")
    l = HystLaw::asym_von_mises(g("E"), g("ft0"), g("Ht"), g("fc0"), g("Hc"));
  else if (law == "damage_simo_ju")
    l = HystLaw::damage_simo_ju(g("E"), g("r0"), g("A"), g("B"));
  else if (law == "plastic_damage_oliver")
    l = HystLaw::plastic_damage_oliver(g("E"), g("fy"), g("h"), g("r0"), g("A"), g("n"));
  else if (const auto dl = engine([&] {
             return design_law_of(law, g, [&](const char* key) { return p.contains(key); });
           }))
    l = HystLaw::wrap(*dl);
  else
    throw py::value_error("uniaxial_path: unknown law \"" + law + "\"");
  return engine([&] {
    std::vector<std::pair<double, double>> out;
    HystState st;
    for (double v : path) {
      double eps = v, sig = 0, Et = 0;
      if (stress_driven) eps = l.strain_at(st, v);
      l.response(st, eps, sig, Et);
      st = l.update(st, eps);
      out.emplace_back(eps, sig);
    }
    return out;
  });
}

// Model.fiber_frame. A named function rather than a lambda inside the module body: MSVC 14.51 stops
// with an internal compiler error on the lambda form.
py::dict fiber_frame_py(const PyModel& self, const py::dict& sections, std::optional<Id> cs, int load_steps, bool corotational,
             std::optional<py::dict> control, const py::list& stages) {
  ModelPtr model = self.solver_model();
  auto law_of = [](const py::dict& d) {
    const std::string k = d["kind"].cast<std::string>();
    auto g = [&](const char* key) { return d[key].cast<double>(); };
    HystLaw l;
    if (k == "steel_kinematic") l = HystLaw::steel_kinematic(g("E"), g("fy"), g("b"));
    else if (k == "steel_mixed") l = HystLaw::steel_mixed(g("E"), g("fy"), g("H"), g("iso"));
    else if (k == "concrete_kj") l = HystLaw::concrete_kj(g("fcd"));
    else if (k == "elastic") l = HystLaw::wrap(UniaxialLaw::elastic(g("E")));
    else if (const auto dl = engine([&] {
               return design_law_of(k, g, [&](const char* key) { return d.contains(key); });
             }))
      l = HystLaw::wrap(*dl);
    else throw py::value_error("fiber_frame: law kind 'steel_kinematic', 'steel_mixed', 'concrete_kj', "
                               "'elastic', 'sargin', 'steel_hardening', 'strand', "
                               "'concrete_tension_softening' or 'bfup_t3'");
    if (d.contains("creep")) {  // dict(a=[1/Pa], tau=[days], fck=0)
      const py::dict c = d["creep"].cast<py::dict>();
      l.creep_a = c["a"].cast<std::vector<double>>();
      l.creep_tau = c["tau"].cast<std::vector<double>>();
      if (c.contains("fck")) l.creep_fck = c["fck"].cast<double>();
      if (l.creep_a.size() != l.creep_tau.size() || l.creep_a.size() > 6)
        throw py::value_error("fiber_frame: creep a and tau paired, at most 6 units");
    }
    if (d.contains("shrink")) l.shrink = d["shrink"].cast<std::vector<std::pair<double, double>>>();  // [(age, eps)]
    return l;
  };
  FiberFrame3DOptions o;
  for (const auto& item : sections) {
    const py::handle kb = item.first;
    const py::dict d = item.second.cast<py::dict>();
    const auto rect = d["rect"].cast<std::tuple<double, double, int, int>>();
    HystFiberSection hs = HystFiberSection::rect(std::get<0>(rect), std::get<1>(rect), std::get<2>(rect), std::get<3>(rect),
                                                 law_of(d["law"].cast<py::dict>()));
    if (d.contains("rebars"))
      for (const py::handle& h : d["rebars"].cast<py::list>()) {
        const py::dict rb = py::reinterpret_borrow<py::dict>(h);
        hs.add_rebar(rb["y"].cast<double>(), rb["z"].cast<double>(), rb["A"].cast<double>(), law_of(rb["law"].cast<py::dict>()));
      }
    o.sections.emplace(kb.cast<Id>(), hs);
    if (d.contains("drying")) {  // dict(faces=['ymin','ymax','zmin','zmax'], h_env, h0, D1, alpha0, hc, n, surface, eps_inf, substeps)
      const py::dict dd = d["drying"].cast<py::dict>();
      DryingSection ds;
      ds.b = std::get<0>(rect);
      ds.h = std::get<1>(rect);
      ds.ny = std::get<2>(rect);
      ds.nz = std::get<3>(rect);
      for (const std::string& fc : dd["faces"].cast<std::vector<std::string>>()) {
        const int k = fc == "ymin" ? 0 : fc == "ymax" ? 1 : fc == "zmin" ? 2 : fc == "zmax" ? 3 : -1;
        if (k < 0) throw py::value_error("fiber_frame: drying faces among 'ymin', 'ymax', 'zmin', 'zmax'");
        ds.exposed[k] = true;
      }
      auto gd = [&](const char* key, double& v) {
        if (dd.contains(key)) v = dd[key].cast<double>();
      };
      gd("h_env", ds.h_env);
      gd("h0", ds.h0);
      gd("D1", ds.D1);
      gd("alpha0", ds.alpha0);
      gd("hc", ds.hc);
      gd("n", ds.n);
      gd("surface", ds.surface);
      gd("eps_inf", ds.eps_inf);
      if (dd.contains("substeps")) ds.substeps = dd["substeps"].cast<int>();
      o.drying[kb.cast<Id>()] = ds;
    }
  }
  o.load_steps = load_steps;
  o.corotational = corotational;
  if (control) {
    o.control = FiberFrame3DOptions::Control::Displacement;
    o.ctrl_node = (*control)["node"].cast<Id>();
    o.ctrl_comp = (*control)["comp"].cast<int>();
    o.ctrl_targets = (*control)["targets"].cast<std::vector<double>>();
    if (control->contains("steps_per_leg")) o.steps_per_leg = (*control)["steps_per_leg"].cast<int>();
  }
  for (const py::handle& h : stages) {
    const py::dict d = py::reinterpret_borrow<py::dict>(h);
    FiberFrame3DOptions::Stage st;
    if (d.contains("activate")) st.activate = d["activate"].cast<std::vector<Id>>();
    if (d.contains("case")) st.case_id = d["case"].cast<Id>();
    if (d.contains("load_steps")) st.load_steps = d["load_steps"].cast<int>();
    if (d.contains("duration")) st.duration = d["duration"].cast<double>();
    if (d.contains("time_steps")) st.time_steps = d["time_steps"].cast<int>();
    if (d.contains("activation_age")) st.activation_age = d["activation_age"].cast<double>();
    o.stages.push_back(st);
  }
  const Id id = o.stages.empty() ? self.case_or_first(cs) : 0;
  const FiberFrame3DResult r = engine([&] { return solve_fiber_frame3d(*model, id, o); });
  py::dict out;
  py::dict disp;
  for (std::size_t i = 0; i < model->nodes.size(); ++i) disp[py::int_(model->nodes[i].id)] = dof6(r.node_disp[i]);
  out["node_disp"] = disp;
  py::list tr;
  for (const FiberTracePoint& p : r.trace) tr.append(std::make_pair(p.lambda, p.u_ctrl));
  out["trace"] = tr;
  out["stage_end"] = r.stage_end;
  out["time"] = r.trace_time;
  py::dict hum;
  for (const auto& item : r.humidity) hum[py::int_(item.first)] = item.second;
  out["humidity"] = hum;
  out["lambda"] = r.lambda;
  return out;
}

}  // namespace

// The module body, cut into several functions: MSVC 14.51 stops with an internal compiler error
// on a single function of this size. Called in this order (a class is registered before the
// functions whose signatures name it).
void bind_module_model(py::module_& mod) {
  // ---------------- Model ----------------
  py::class_<PyModel> model_cls(mod, "Model", "The model builder. All add_* methods hand back integer ids.");
  model_cls.def(py::init<>());
  model_cls.def("add_node", &PyModel::add_node, py::arg("x"), py::arg("y"), py::arg("z"),
           "Add a node at (x, y, z) in metres. Returns its id.");
  model_cls.def("steel", &PyModel::steel, "Structural steel preset (E = 210 GPa, nu = 0.3, rho = 7850 kg/m^3).");
  model_cls.def("concrete", &PyModel::concrete, "Concrete preset (E = 33 GPa, nu = 0.2, rho = 2500 kg/m^3).");
  model_cls.def("add_material", &PyModel::add_material, py::arg("e"), py::arg("nu"), py::arg("rho") = py::none(),
           py::arg("alpha_temp") = py::none(), py::arg("name") = py::none(), py::arg("g") = py::none(),
           py::arg("damping") = py::none(),
           "Custom material. `e` in Pa, `rho` in kg/m^3 (None: massless), `alpha_temp` in 1/K, `g` the shear "
           "modulus of the bars in Pa when it is not E / 2(1 + nu) (None: isotropic; orthotropic timber, E/G ~ 16). "
           "A material with `g` cannot be used by shells or continuum elements. `damping`: the damping ratio of the "
           "material, read by strain_energy_damping (None: 0).");
  model_cls.def("add_section", &PyModel::add_section, py::arg("designation"),
           "European catalog steel profile by designation (\"IPE300\", \"HEB200\", ...).");
  model_cls.def("add_section_rect", &PyModel::add_section_rect, py::arg("b"), py::arg("h"),
           "Solid rectangle, `b` along local y, `h` along local z (metres).");
  model_cls.def("add_section_tube", &PyModel::add_section_tube, py::arg("d"), py::arg("t"),
           "Circular hollow section (outer diameter `d`, wall `t`).");
  model_cls.def("add_section_i", &PyModel::add_section_i, py::arg("h"), py::arg("b"), py::arg("tw"), py::arg("tf"),
           "Welded / custom doubly-symmetric I-section.");
  model_cls.def("add_section_props", &PyModel::add_section_props, py::arg("sx"), py::arg("iy"), py::arg("iz"), py::arg("ix"),
           py::arg("sy") = py::none(), py::arg("sz") = py::none(), py::arg("name") = py::none(),
           py::arg("iw") = 0.0, py::arg("beta_x") = 0.0,
           "Raw section properties: area sx, second moments iy / iz, torsion ix, shear areas sy / sz (None: rigid), "
           "warping constant iw [m^6] and Wagner coefficient beta_x [m] (read by warping / lateral_buckling).");
  model_cls.def("add_beam", &PyModel::add_beam, py::arg("n1"), py::arg("n2"), py::arg("section"), py::arg("material"),
           py::arg("beta") = 0.0, py::arg("ecc") = py::none(), py::arg("ecc_end") = py::none(),
           "A frame element n1 -> n2; `beta` rolls the section [rad]; `ecc` = (ex, ey, ez) rigid-arm offset (global) "
           "of both ends, `ecc_end` that of the end n2 when it differs.");
  model_cls.def("add_truss", &PyModel::add_truss, py::arg("n1"), py::arg("n2"), py::arg("section"), py::arg("material"),
           "An axial-only truss bar.");
  model_cls.def("add_cable", &PyModel::add_cable, py::arg("n1"), py::arg("n2"), py::arg("ea"), py::arg("w") = 0.0,
           py::arg("l0") = py::none(), py::arg("h") = py::none(), py::arg("end_tension") = py::none(),
           "An exact elastic-catenary cable: EA (N), weight w (N/m) and exactly one of l0= (unstrained length), "
           "h= (horizontal tension) or end_tension= (tension at end 2). Solve with cable_static().");
  model_cls.def("add_point_mass", &PyModel::add_point_mass, py::arg("node"), py::arg("mass"), py::arg("inertia") = py::none(),
           "Lump a mass (kg) and optional rotary inertia [Ixx, Iyy, Izz] (kg.m^2) onto a node.");
  model_cls.def("set_releases", &PyModel::set_releases, py::arg("bar"), py::arg("ry_or") = false, py::arg("rz_or") = false,
           py::arg("ry_ex") = false, py::arg("rz_ex") = false, py::arg("origin") = "", py::arg("end") = "",
           "Release components at the ends of a beam, local frame: the booleans release the bending rotations "
           "(\"or\" = origin, \"ex\" = end); origin= / end= take any of 'dx dy dz rx ry rz' (e.g. end='dx rx'). "
           "A release that leaves a mechanism is a typed error at the analysis.");
  model_cls.def("set_soil", &PyModel::set_soil, py::arg("bar"), py::arg("kfx") = 0.0, py::arg("kfy") = 0.0,
           py::arg("kfz") = 0.0, py::arg("kmx") = 0.0, "Winkler bed in the bar-local frame (N/m per m, N.m/m per rad).");
  model_cls.def("fix", &PyModel::fix, py::arg("node"), "Clamp: all six DOFs blocked.");
  model_cls.def("pin", &PyModel::pin, py::arg("node"), "Pin: translations blocked, rotations free.");
  model_cls.def("add_support", &PyModel::add_support, py::arg("node"), py::arg("dx") = false, py::arg("dy") = false,
           py::arg("dz") = false, py::arg("rx") = false, py::arg("ry") = false, py::arg("rz") = false,
           py::arg("kx") = py::none(), py::arg("ky") = py::none(), py::arg("kz") = py::none(),
           py::arg("krx") = py::none(), py::arg("kry") = py::none(), py::arg("krz") = py::none(),
           "Block the named DOFs; a spring k* acts on a DOF that is not blocked.");
  model_cls.def("add_rigid_link", &PyModel::add_rigid_link, py::arg("master"), py::arg("slaves"), py::arg("dx") = true,
           py::arg("dy") = true, py::arg("dz") = true, py::arg("rx") = true, py::arg("ry") = true,
           py::arg("rz") = true,
           "Rigid master-slave link (exact elimination): each slave's coupled global components follow "
           "the master's rigid-body motion. Linear static, modal, spectrum, buckling and phasing analyses.");
  model_cls.def("add_elastic_link", &PyModel::add_elastic_link, py::arg("ni"), py::arg("nj"), py::arg("kx") = 0.0,
           py::arg("ky") = 0.0, py::arg("kz") = 0.0, py::arg("krx") = 0.0, py::arg("kry") = 0.0,
           py::arg("krz") = 0.0, py::arg("s") = 0.5, py::arg("roll") = 0.0,
           "Elastic link between two nodes: springs in the link frame (x along ni -> nj, global for "
           "coincident nodes) at the point ni + s (nj - ni), reached by rigid arms.");
  model_cls.def("set_link_law", &PyModel::set_link_law, py::arg("link"), py::arg("comp"), py::arg("kind"),
           py::arg("k") = 0.0, py::arg("param") = 0.0, py::arg("kp") = 0.0, py::arg("points") = py::none(),
           "Non-linear law on one link-frame component: 'tension' / 'compression' (stiffness k on one side "
           "only), 'gap' (clearance param either side, then k), 'bilinear' (k up to |F| = param, then kp), "
           "'points' (non-decreasing [(d, F), ...]; F(0) may differ from 0: a force present at rest). Solved by "
           "analyze_link_laws.");
  model_cls.def("add_unit_prestress_case", &PyModel::add_unit_prestress_case, py::arg("bar"), py::arg("P") = 1.0,
           "A load case tensioning `bar` by P where it can shorten freely (imposed strain -P / EA): the unit "
           "tensile load case of cable tuning. Returns the case id.");
  model_cls.def(
          "tune_linear",
          [](const PyModel& self, const std::vector<Id>& base, const py::list& unknowns, const py::list& targets) {
            auto num = [](const py::dict& d, const char* k, double dflt) {
              return d.contains(k) ? d[k].cast<double>() : dflt;
            };
            std::vector<TuningUnknown> us;
            for (const py::handle& h : unknowns) {
              TuningUnknown u;
              if (py::isinstance<py::int_>(h)) {
                u.case_id = h.cast<Id>();
              } else {
                const py::dict d = py::reinterpret_borrow<py::dict>(h);
                u.case_id = d["case"].cast<Id>();
                u.lower = num(d, "lower", -std::numeric_limits<double>::infinity());
                u.upper = num(d, "upper", std::numeric_limits<double>::infinity());
                u.cost = num(d, "cost", 0);
                u.ref = num(d, "ref", 0);
              }
              us.push_back(u);
            }
            std::vector<TuningTarget> ts;
            for (const py::handle& h : targets) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              TuningTarget t;
              const std::string kind = d.contains("kind") ? d["kind"].cast<std::string>() : "displacement";
              if (kind == "displacement")
                t.kind = TuningTarget::Kind::Displacement;
              else if (kind == "reaction")
                t.kind = TuningTarget::Kind::Reaction;
              else if (kind == "bar_force")
                t.kind = TuningTarget::Kind::BarForce;
              else if (kind == "link_force")
                t.kind = TuningTarget::Kind::LinkForce;
              else
                throw py::value_error("tune_linear: target kind 'displacement', 'reaction', 'bar_force' or "
                                      "'link_force', got '" + kind + "'");
              const std::string mode = d.contains("mode") ? d["mode"].cast<std::string>() : "equal";
              if (mode == "equal")
                t.mode = TuningTarget::Mode::Equal;
              else if (mode == "at_most")
                t.mode = TuningTarget::Mode::AtMost;
              else if (mode == "at_least")
                t.mode = TuningTarget::Mode::AtLeast;
              else if (mode == "fit")
                t.mode = TuningTarget::Mode::Fit;
              else
                throw py::value_error("tune_linear: mode 'equal', 'at_most', 'at_least' or 'fit', got '" + mode + "'");
              t.id = d["id"].cast<Id>();
              t.comp = d.contains("comp") ? d["comp"].cast<int>() : 0;
              t.station = d.contains("station") ? d["station"].cast<int>() : 0;
              t.value = num(d, "value", 0);
              t.weight = num(d, "weight", 1);
              ts.push_back(t);
            }
            ModelPtr model = self.solver_model();
            const TuningResult r = engine([&] { return tune_linear(*model, base, us, ts); });
            py::dict o;
            o["x"] = r.x;
            o["achieved"] = r.achieved;
            o["base"] = r.base;
            return o;
          },
          py::arg("base"), py::arg("unknowns"), py::arg("targets"),
          "Linear tuning: the factors of the `unknowns` cases (ids, or dict(case, lower, upper, cost, ref)) that "
          "bring the `targets` (dict(kind, id, comp, value, mode, weight, station)) of base + sum x_j case_j to "
          "their values — exact, optimised (costs), best fit ('fit'), with bounds; an exact convex QP.");
  model_cls.def(
          "tune_cables",
          [](const PyModel& self, Id cs, const py::list& unknowns, const py::list& targets, double tol, int max_iter,
             int substeps, bool large_displacements) {
            auto num = [](const py::dict& d, const char* k, double dflt) {
              return d.contains(k) ? d[k].cast<double>() : dflt;
            };
            std::vector<TuningUnknown> us;
            for (const py::handle& h : unknowns) {
              TuningUnknown u;
              if (py::isinstance<py::int_>(h)) {
                u.case_id = self.cable_bar_id(h.cast<Id>());
              } else {
                const py::dict d = py::reinterpret_borrow<py::dict>(h);
                if (d.contains("cable") == d.contains("bar"))
                  throw py::value_error("tune_cables: an unknown is dict(cable=id) or dict(bar=id, ...)");
                if (d.contains("cable")) {
                  u.case_id = self.cable_bar_id(d["cable"].cast<Id>());
                } else {
                  u.case_id = d["bar"].cast<Id>();
                  u.kind = TuningUnknown::Kind::BarShortening;
                }
                u.lower = num(d, "lower", -std::numeric_limits<double>::infinity());
                u.upper = num(d, "upper", std::numeric_limits<double>::infinity());
                u.cost = num(d, "cost", 0);
                u.ref = num(d, "ref", 0);
              }
              us.push_back(u);
            }
            std::vector<TuningTarget> ts;
            for (const py::handle& h : targets) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              TuningTarget t;
              const std::string kind = d.contains("kind") ? d["kind"].cast<std::string>() : "displacement";
              t.id = d["id"].cast<Id>();
              if (kind == "displacement")
                t.kind = TuningTarget::Kind::Displacement;
              else if (kind == "reaction")
                t.kind = TuningTarget::Kind::Reaction;
              else if (kind == "bar_force")
                t.kind = TuningTarget::Kind::BarForce;
              else if (kind == "cable_tension") {
                t.kind = TuningTarget::Kind::CableTension;
                t.id = self.cable_bar_id(t.id);
              } else
                throw py::value_error("tune_cables: target kind 'displacement', 'reaction', 'bar_force' or "
                                      "'cable_tension', got '" + kind + "'");
              const std::string mode = d.contains("mode") ? d["mode"].cast<std::string>() : "equal";
              if (mode == "equal")
                t.mode = TuningTarget::Mode::Equal;
              else if (mode == "at_most")
                t.mode = TuningTarget::Mode::AtMost;
              else if (mode == "at_least")
                t.mode = TuningTarget::Mode::AtLeast;
              else if (mode == "fit")
                t.mode = TuningTarget::Mode::Fit;
              else
                throw py::value_error("tune_cables: mode 'equal', 'at_most', 'at_least' or 'fit', got '" + mode + "'");
              t.comp = d.contains("comp") ? d["comp"].cast<int>() : 0;
              t.station = d.contains("station") ? d["station"].cast<int>() : 0;
              t.value = num(d, "value", 0);
              t.weight = num(d, "weight", 1);
              ts.push_back(t);
            }
            ModelPtr model = self.solver_model();
            CableTuningOptions o;
            o.tol = tol;
            o.max_iter = max_iter;
            o.statics.substeps = std::max(substeps, 1);
            o.substeps = std::max(substeps, 1);
            o.large_displacements = large_displacements;
            const TuningResult r = engine([&] { return tune_cables(*model, cs, us, ts, o); });
            py::dict out;
            out["x"] = r.x;
            out["achieved"] = r.achieved;
            out["iterations"] = r.iterations;
            py::dict disp;
            for (std::size_t i = 0; i < model->nodes.size() && i < r.state.node_disp.size(); ++i)
              disp[py::int_(model->nodes[i].id)] = dof6(r.state.node_disp[i]);
            out["displacement"] = disp;
            return out;
          },
          py::arg("case"), py::arg("unknowns"), py::arg("targets"), py::arg("tol") = 1e-9, py::arg("max_iter") = 50,
          py::arg("substeps") = 1, py::arg("large_displacements") = false,
          "General form finding: the unstrained lengths of cables (cable ids, or dict(cable=id, lower, upper, cost, "
          "ref)) and the lack of fit of bars (dict(bar=id, ...): fabricated short by x, an imposed strain -x / L) "
          "that bring the `targets` (dict(kind 'displacement' | 'reaction' | 'bar_force' | 'cable_tension', id, "
          "comp, value, mode 'equal' | 'at_most' | 'at_least' | 'fit', weight, station)) of case `case` to their "
          "values, on the exact cable statics (or large displacements). Dependent equal targets are a typed error "
          "naming them. Returns dict(x, achieved, iterations, displacement {node: 6-tuple}).");
  model_cls.def("add_constraint", &PyModel::add_constraint, py::arg("terms"),
           "General linear constraint sum coef u(node, comp) = 0, terms [(node, comp 0..5, coef)]; the first "
           "term's DOF is eliminated (exact master-slave reduction).");
  model_cls.def("add_skew_support", &PyModel::add_skew_support, py::arg("node"), py::arg("axes"), py::arg("dx") = false,
           py::arg("dy") = false, py::arg("dz") = false, py::arg("rx") = false, py::arg("ry") = false,
           py::arg("rz") = false,
           "Support in any frame: axes = rows local x, y, z (orthonormal); blocks the named LOCAL components. "
           "Reactions in CaseResult.skew_reactions() (local, global).");
  model_cls.def("add_ground_spring", &PyModel::add_ground_spring, py::arg("node"), py::arg("direction"), py::arg("points"),
           "Zero-length spring from a new fully fixed node to `node`, its x axis along `direction`, with the law "
           "points [(d, F), ...] (d = node movement along direction, F > 0 tension). Returns the link id.");
  model_cls.def("add_earth_pressure", &PyModel::add_earth_pressure, py::arg("wall"), py::arg("away"), py::arg("up"),
           py::arg("layers"), py::arg("width") = 1.0, py::arg("distribution") = "consistent",
           "Trilinear earth-pressure springs (at rest at zero displacement, active / passive plateaus, kh = kh_A + "
           "kh_B (top - y)^kh_kappa) on a chain of wall nodes; `away` = the active direction (away from the soil); "
           "layers = [dict(top, bottom, sv_top, gamma, Ka, K0, Kp, c, Kac, Kpc, kh_A, kh_B, kh_kappa)]; "
           "distribution 'consistent' or 'lumped' (nodal pressure times tributary length). Returns the link ids "
           "per node (0 = none).");
  model_cls.def(
          "analyze_link_laws",
          [](const PyModel& self, std::optional<Id> cs, int stations) {
            ModelPtr model = self.solver_model();
            const Id id = self.case_or_first(cs);
            LinkLawOptions o;
            o.stations_per_bar = stations;
            auto r = engine([&] { return std::make_shared<const LinkLawResult>(solve_static_link_laws(*model, id, o)); });
            return PyLinkLaw{model, r};
          },
          py::arg("case") = py::none(), py::arg("stations") = 11,
          "Static analysis of one case with the non-linear link laws (exact, finite Newton on the convex energy).");
  model_cls.def(
          "continuum_nonlinear",
          [](const PyModel& self, Id cs, std::vector<double> lambdas, int substeps, const std::map<Id, std::vector<double>>& plastic,
             bool large, std::vector<double> times, const py::dict& frictional, bool stop_at_limit,
             const std::map<Id, double>& temperature) {
            ModelPtr model = self.solver_model();
            ContinuumPlasticOptions o;
            o.lambdas = lambdas;
            o.substeps = substeps;
            o.large_displacements = large;
            o.times = times;
            for (const auto& [mid, v] : plastic) {
              if (v.size() != 2 && v.size() != 5 && v.size() != 6)
                throw py::value_error(
                    "continuum_nonlinear: plastic = {material: (sigma_y, H) or (sigma_y, H, A, n, m[, Q])}");
              J2Material j{v[0], v[1]};
              if (v.size() >= 5) {
                j.creep_A = v[2];
                j.creep_n = v[3];
                j.creep_m = v[4];
              }
              if (v.size() == 6) j.creep_Q = v[5];
              o.plastic[mid] = j;
            }
            for (const auto& [k, v] : frictional) {
              const py::dict d = py::cast<py::dict>(v);
              FrictionalMaterial fm;
              const std::string kind = d.contains("kind") ? py::cast<std::string>(d["kind"]) : "mohr_coulomb";
              if (kind == "mohr_coulomb") fm.kind = FrictionalMaterial::Kind::MohrCoulomb;
              else if (kind == "drucker_prager") fm.kind = FrictionalMaterial::Kind::DruckerPrager;
              else throw py::value_error("continuum_nonlinear: frictional kind mohr_coulomb | drucker_prager");
              fm.c = d.contains("c") ? py::cast<double>(d["c"]) : 0.0;
              fm.phi = d.contains("phi") ? py::cast<double>(d["phi"]) : 0.0;
              fm.psi = d.contains("psi") ? py::cast<double>(d["psi"]) : fm.phi;
              fm.H = d.contains("H") ? py::cast<double>(d["H"]) : 0.0;
              const std::string fit = d.contains("fit") ? py::cast<std::string>(d["fit"]) : "outer";
              if (fit == "outer") fm.fit = FrictionalMaterial::DpFit::Outer;
              else if (fit == "inner") fm.fit = FrictionalMaterial::DpFit::Inner;
              else if (fit == "plane_strain") fm.fit = FrictionalMaterial::DpFit::PlaneStrain;
              else throw py::value_error("continuum_nonlinear: fit outer | inner | plane_strain");
              o.frictional[py::cast<Id>(k)] = fm;
            }
            o.stop_at_limit = stop_at_limit;
            o.temperature = temperature;
            const ContinuumPlasticResult r = engine([&] { return solve_continuum_plastic(*model, cs, o); });
            py::list steps;
            for (const auto& snap : r.node_disp) {
              py::dict d;
              for (std::size_t i = 0; i < model->nodes.size(); ++i) d[py::int_(model->nodes[i].id)] = dof6(snap[i]);
              steps.append(d);
            }
            py::dict out;
            out["steps"] = steps;
            out["gauss"] = r.gauss;
            out["iterations"] = r.iterations;
            py::list reac;
            for (const auto& snap : r.reactions) {
              py::dict dd;
              for (std::size_t i = 0; i < model->nodes.size(); ++i)
                if (snap[i].squaredNorm() > 0) dd[py::int_(model->nodes[i].id)] = dof6(snap[i]);
              reac.append(dd);
            }
            out["reactions"] = reac;
            out["lambdas"] = r.lambdas_done;
            out["limit_reached"] = r.limit_reached;
            return out;
          },
          py::arg("case") = 1, py::arg("lambdas") = std::vector<double>{1.0}, py::arg("substeps") = 4,
          py::arg("plastic") = std::map<Id, std::vector<double>>{}, py::arg("large_displacements") = false,
          py::arg("times") = std::vector<double>{}, py::arg("frictional") = py::dict(), py::arg("stop_at_limit") = false,
          py::arg("temperature") = std::map<Id, double>{},
          "Nonlinear continuum analysis (keta/continuum_plastic.h): von Mises plasticity, plastic = {material: (sigma_y, "
          "H)} or (sigma_y, H, A, n, m[, Q]) with time-hardening creep eps_c = A q^n t^m [exp(-Q / T), T interpolated from "
          "temperature = {node: T}] (then `times`: the time of each load factor), frictional = {material: dict(kind='mohr_coulomb'|'drucker_prager', c, phi, psi, H, fit='outer'|'inner'|"
          "'plane_strain')} (angles in rad), the case times each load factor of `lambdas` (prescribed support displacements "
          "too); stop_at_limit returns the converged factors; or total Lagrangian large displacements (elastic). Returns "
          "reactions and lambdas too, "
          "steps [{node: (ux..rz)}], gauss stresses (+ equivalent plastic strain) per element.");
  model_cls.def(
          "fiber_frame", &fiber_frame_py,
          py::arg("sections"), py::arg("case") = py::none(), py::arg("load_steps") = 10, py::arg("corotational") = false,
          py::arg("control") = py::none(), py::arg("stages") = py::list(),
          "Force-based fiber frame (keta/fiber_frame.h): sections = {bar: dict(rect=(b, h, ny, nz), law=dict(kind, ...), "
          "rebars=[dict(y, z, A, law)])}; load control, or control=dict(node, comp, targets, steps_per_leg) for "
          "displacement control; stages = [dict(activate=[bars], case, load_steps, duration=days, time_steps)] for "
          "staged construction with creep between phases (a law's creep=dict(a=[1/Pa], tau=[days], fck): Kelvin "
          "chain, EC2 nonlinear beyond 0.45 fck; shrink=[(age, eps)] imposed from the activation, activation_age on "
          "the stage; a section's drying=dict(faces, h_env, D1, alpha0, hc, n, surface, eps_inf): Bazant-Najjar "
          "diffusion on the fiber grid, each fiber's shrinkage from its humidity). Returns node_disp, humidity, the trace (lambda, u_ctrl), its times and the "
          "stage ends.");
  model_cls.def(
          "link_laws_path",
          [](const PyModel& self, Id cs, Id ctrl_node, int ctrl_comp, std::vector<double> targets, int steps,
             std::optional<Id> dead_case, const py::list& families, std::optional<py::dict> arc) {
            ModelPtr model = self.solver_model();
            LinkPathOptions o;
            if (arc) {
              const py::dict& a = *arc;
              o.arc_steps = a.contains("steps") ? a["steps"].cast<int>() : 100;
              o.arc_length = a["length"].cast<double>();
              o.arc_psi = a.contains("psi") ? a["psi"].cast<double>() : 0.0;
              if (a.contains("stop")) o.arc_stop_ctrl = a["stop"].cast<double>();
              o.arc_first_sign = a.contains("sign") ? a["sign"].cast<int>() : 1;
            }
            o.ctrl_node = ctrl_node;
            o.ctrl_comp = ctrl_comp;
            o.targets = targets;
            o.steps_per_leg = steps;
            if (dead_case) o.dead_case = *dead_case;
            for (const py::handle& h : families) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              auto fam = std::make_shared<HingeLawFamily>();
              fam->N = d["N"].cast<std::vector<double>>();
              for (const auto& pts : d["laws"].cast<std::vector<std::vector<std::pair<double, double>>>>()) {
                LinkLaw l;
                l.pts = pts;
                fam->laws.push_back(l);
              }
              LinkHingeFamily lf;
              lf.link = d["link"].cast<Id>();
              lf.comp = d.contains("comp") ? d["comp"].cast<int>() : 5;
              lf.axial_comp = d.contains("axial_comp") ? d["axial_comp"].cast<int>() : 0;
              lf.family = fam;
              o.families.push_back(lf);
            }
            const LinkPathResult r = engine([&] { return solve_link_laws_path(*model, cs, o); });
            py::list tr;
            for (const LinkPathPoint& p : r.trace) {
              py::dict d;
              d["lambda"] = p.lambda;
              d["u"] = p.u_ctrl;
              d["d"] = p.d;
              d["F"] = p.F;
              tr.append(d);
            }
            return tr;
          },
          py::arg("case"), py::arg("node"), py::arg("comp"), py::arg("targets"), py::arg("steps") = 20,
          py::arg("dead_case") = py::none(), py::arg("families") = py::list(), py::arg("arc") = py::none(),
          "Static path through SOFTENING link laws under displacement control of (node, comp) along `targets`, or "
          "with arc = dict(length, steps=100, psi=0, stop=u, sign=1) traced event to event by arc length (snap-backs; "
          "targets then unused): "
          "families = [dict(link, comp=5, axial_comp=0, N=[..], laws=[[(theta, M), ..], ..])] (M-theta following the "
          "link's axial force); exact equilibrium at every point (keta/link_nl.h). Returns [dict(lambda, u, d, F)].");
  model_cls.def("add_area_load", &PyModel::add_area_load, py::arg("polygon"), py::arg("qx") = 0.0, py::arg("qy") = 0.0,
           py::arg("qz") = 0.0, py::arg("mode") = "two_way", py::arg("span") = std::make_tuple(1.0, 0.0, 0.0),
           py::arg("case") = 1,
           "Load distribution zone: a uniform area load (N/m^2) on a convex panel bounded by bars, handed to "
           "them by the 'two_way' (nearest edge, 45-degree lines) or 'one_way' (strips along `span`) rule. "
           "Returns the number of bar loads added.");
  model_cls.def("add_shell_tendon", &PyModel::add_shell_tendon, py::arg("points"), py::arg("P0"), py::arg("mu") = 0.0,
           py::arg("k") = 0.0, py::arg("jack") = "start", py::arg("slab") = 0,
           "Tendon in the shells: polyline at its true position (inside the thickness), jacking force P0, "
           "friction mu [1/rad] and wobble k [1/m], jacked 'start' / 'end' / 'both'. Returns its id.");
  model_cls.def("add_shell_prestress", &PyModel::add_shell_prestress, py::arg("tendon"), py::arg("factor") = 1.0,
           py::arg("case") = 1, "The tendon's equivalent loads (anchors, deviations, friction) in a case.");
  model_cls.def("add_polygon_pressure", &PyModel::add_polygon_pressure, py::arg("polygon"), py::arg("qx") = 0.0,
           py::arg("qy") = 0.0, py::arg("qz") = 0.0,
           py::arg("holes") = std::vector<std::vector<std::tuple<double, double, double>>>{}, py::arg("slab") = 0,
           py::arg("case") = 1,
           "Free polygonal pressure (N/m^2 of plan area, global vector) on the shells it covers, openings "
           "excluded; slab 0 = every shell.");
  model_cls.def("add_continuum", &PyModel::add_continuum, py::arg("kind"), py::arg("nodes"), py::arg("material"), py::arg("t") = 1.0,
           "Continuum element (keta/continuum.h): 'plane_stress' / 'plane_strain' / 'axisymmetric' (T3, Q4, T6, Q8 in "
           "the XY plane; axisymmetric r = x, axis y, per radian) or 'solid' (T4, H8, H20); returns its id.");
  model_cls.def("add_continuum_traction", &PyModel::add_continuum_traction, py::arg("element"), py::arg("face"),
           py::arg("traction"), py::arg("case") = 1,
           "Uniform traction (global vector, force per area) on an edge / hexahedron face: exact consistent nodal loads.");
  model_cls.def("add_continuum_pressure", &PyModel::add_continuum_pressure, py::arg("element"), py::arg("face"), py::arg("p"),
           py::arg("case") = 1, "Uniform pressure (positive pushing into the element) on an edge / hexahedron face.");
  model_cls.def("add_continuum_thermal", &PyModel::add_continuum_thermal, py::arg("element"), py::arg("T"), py::arg("T0") = 0.0,
           py::arg("case") = 1,
           "Thermal load of a continuum element (keta/conduction.h): the nodal temperatures T = {node: T} (e.g. from "
           "conduction()), reference T0, the material's alpha; consistent nodal forces added to the case.");
  model_cls.def(
          "conduction",
          [](const PyModel& self, const py::dict& materials, const std::map<Id, double>& fixed, const py::list& faces,
             const std::map<Id, double>& sources, const py::list& links, std::vector<double> times, int substeps,
             double initial, const std::map<Id, double>& initial_nodes) {
            ModelPtr model = self.solver_model();
            ConductionProblem p;
            for (const auto& [k, v] : materials) {
              ConductionMaterial cm;
              if (py::isinstance<py::dict>(v)) {
                const py::dict d = py::cast<py::dict>(v);
                const double kk = d.contains("k") ? py::cast<double>(d["k"]) : 1.0;
                cm.kx = d.contains("kx") ? py::cast<double>(d["kx"]) : kk;
                cm.ky = d.contains("ky") ? py::cast<double>(d["ky"]) : kk;
                cm.kz = d.contains("kz") ? py::cast<double>(d["kz"]) : kk;
                cm.k1 = d.contains("k1") ? py::cast<double>(d["k1"]) : 0.0;
                cm.capacity = d.contains("capacity") ? py::cast<double>(d["capacity"]) : 0.0;
              } else {
                cm.kx = cm.ky = cm.kz = py::cast<double>(v);
              }
              p.materials[py::cast<Id>(k)] = cm;
            }
            p.fixed = fixed;
            for (const py::handle& h : faces) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              ConductionProblem::Face f;
              f.element = py::cast<Id>(d["element"]);
              f.face = py::cast<int>(d["face"]);
              f.q = d.contains("q") ? py::cast<double>(d["q"]) : 0.0;
              f.h = d.contains("h") ? py::cast<double>(d["h"]) : 0.0;
              f.u_env = d.contains("u_env") ? py::cast<double>(d["u_env"]) : 0.0;
              p.faces.push_back(f);
            }
            p.sources = sources;
            for (const py::handle& h : links) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              ConductionProblem::Link l;
              l.a = py::cast<Id>(d["a"]);
              l.b = d.contains("b") ? py::cast<Id>(d["b"]) : 0;
              l.k = py::cast<double>(d["k"]);
              l.u_env = d.contains("u_env") ? py::cast<double>(d["u_env"]) : 0.0;
              p.links.push_back(l);
            }
            p.u_initial = initial;
            p.initial = initial_nodes;
            ConductionOptions o;
            o.times = times;
            o.substeps = substeps;
            const ConductionResult r = engine([&] { return solve_conduction(*model, p, o); });
            py::dict out;
            out["u"] = r.u;
            out["reactions"] = r.reactions;
            out["iterations"] = r.iterations;
            return out;
          },
          py::arg("materials"), py::arg("fixed") = std::map<Id, double>{}, py::arg("faces") = py::list(),
          py::arg("sources") = std::map<Id, double>{}, py::arg("links") = py::list(), py::arg("times") = std::vector<double>{},
          py::arg("substeps") = 1, py::arg("initial") = 0.0, py::arg("initial_nodes") = std::map<Id, double>{},
          "Steady / transient scalar conduction (heat, seepage) on the continuum elements (keta/conduction.h): "
          "materials {mat: k or dict(k | kx, ky, kz, k1, capacity)}, fixed {node: u}, faces [dict(element, face, q, h, "
          "u_env)] (flux in, convection), sources {node: Q}, links [dict(a, b (0 = environment), k, u_env)]; times "
          "empty = steady, else backward Euler from initial / initial_nodes. Returns u and reactions per step "
          "({node: value}) and the Newton iterations.");
  model_cls.def("add_shell", &PyModel::add_shell, py::arg("n1"), py::arg("n2"), py::arg("n3"), py::arg("t"), py::arg("material"),
           py::arg("slab") = 1, py::arg("offset") = 0.0, py::arg("t_corners") = py::none(), py::arg("bed") = 0.0,
           py::arg("ortho") = py::none(), py::arg("laminate") = py::none(),
           "Flat facet shell (keta/shell.h) on three nodes: thickness, mid-surface offset, corner thicknesses, Winkler bed, "
           "orthotropy dict(Ex, Ey, nu_xy, Gxy, Gxz, Gyz, angle) or laminate dict(plies=[dict(t, E1, E2, nu12, G12, G13, G23, "
           "angle)], z_bottom); returns its id.");
  model_cls.def("add_surface_pressure", &PyModel::add_surface_pressure, py::arg("slab"), py::arg("qx") = 0.0, py::arg("qy") = 0.0,
           py::arg("qz") = 0.0, py::arg("frame") = "global", py::arg("case") = 1,
           "Uniform pressure on the shells of `slab`: 'global' (force/area vector), 'normal' (qz along each facet normal) or "
           "'projected' (per plan area).");
  model_cls.def("add_hydrostatic", &PyModel::add_hydrostatic, py::arg("gamma"), py::arg("surface_point"), py::arg("slab") = 0,
           py::arg("side") = 1.0, py::arg("down") = std::array<double, 3>{0, 0, -1}, py::arg("case") = 1,
           "Hydrostatic pressure gamma x depth below the free surface through surface_point (depth along `down`, "
           "zero above) on the shells of `slab` (0: all), pushing along side x the facet normal; exact consistent "
           "loads, facets cut by the free surface integrated over their wet part.");
  model_cls.def("add_point_load", &PyModel::add_point_load, py::arg("node"), py::arg("fx") = 0.0, py::arg("fy") = 0.0,
           py::arg("fz") = 0.0, py::arg("mx") = 0.0, py::arg("my") = 0.0, py::arg("mz") = 0.0, py::arg("case") = 1,
           "Concentrated force / moment at a node (global frame, N and N.m).");
  model_cls.def("add_udl", &PyModel::add_udl, py::arg("bar"), py::arg("wx") = 0.0, py::arg("wy") = 0.0, py::arg("wz") = 0.0,
           py::arg("case") = 1, py::arg("frame") = "global", py::arg("mx") = 0.0, py::arg("my") = 0.0,
           py::arg("mz") = 0.0, py::arg("x1") = py::none(), py::arg("x2") = py::none(), py::arg("absolute") = false,
           "Uniform load on a bar (N/m): 'global' | 'local' | 'projected'; mx, my, mz = distributed couple [N m/m]; "
           "x1, x2 = loaded range (relative in [0, 1], or [m] from the origin with absolute=True; default: the "
           "whole bar).");
  model_cls.def("add_linear_load", &PyModel::add_linear_load, py::arg("bar"), py::arg("w1"), py::arg("w2"),
           py::arg("x1") = 0.0, py::arg("x2") = 1.0, py::arg("case") = 1, py::arg("frame") = "global",
           py::arg("m1") = py::none(), py::arg("m2") = py::none(), py::arg("absolute") = false,
           "Load ramping from w1 at x1 to w2 at x2 (relative abscissae in [0, 1], or [m] with absolute=True); "
           "m1, m2 = distributed couple (mx, my, mz) [N m/m] at x1 and x2.");
  model_cls.def("add_parabolic_load", &PyModel::add_parabolic_load, py::arg("bar"), py::arg("w1"), py::arg("wmid"),
           py::arg("w2"), py::arg("x1") = 0.0, py::arg("x2") = 1.0, py::arg("case") = 1, py::arg("frame") = "global",
           py::arg("absolute") = false,
           "Load along the parabola through w1 at x1, wmid at (x1 + x2) / 2 and w2 at x2 (each (wx, wy, wz) [N/m]; "
           "relative abscissae in [0, 1], or [m] with absolute=True). Exact equivalent loads, internal forces and "
           "deflections in the linear statics; analyses that carry linear span loads only refuse it by name.");
  model_cls.def("add_bar_point_load", &PyModel::add_bar_point_load, py::arg("bar"), py::arg("at"), py::arg("fx") = 0.0,
           py::arg("fy") = 0.0, py::arg("fz") = 0.0, py::arg("case") = 1, py::arg("frame") = "global",
           py::arg("mx") = 0.0, py::arg("my") = 0.0, py::arg("mz") = 0.0, py::arg("absolute") = false,
           "Concentrated force (fx, fy, fz) and couple (mx, my, mz) on a bar at abscissa `at`: relative in [0, 1], "
           "or [m] from the origin with absolute=True.");
  model_cls.def("add_thermal", &PyModel::add_thermal, py::arg("bar"), py::arg("dt"), py::arg("gy") = 0.0, py::arg("gz") = 0.0,
           py::arg("case") = 1, "Uniform temperature change dt (K) and gradients gy, gz (K/m).");
  model_cls.def("add_support_disp", &PyModel::add_support_disp, py::arg("node"), py::arg("dx") = py::none(),
           py::arg("dy") = py::none(), py::arg("dz") = py::none(), py::arg("rx") = py::none(), py::arg("ry") = py::none(),
           py::arg("rz") = py::none(), py::arg("case") = 1,
           "Imposed support displacement (settlement) on rigidly blocked DOFs (m, rad).");
  model_cls.def("add_self_weight", &PyModel::add_self_weight, py::arg("case") = 1, py::arg("factor") = 1.0,
           py::arg("g") = py::none(), "Self-weight: g = (0, 0, -9.80665 factor), or factor * g when given.");
  model_cls.def("name_case", &PyModel::name_case, py::arg("case"), py::arg("name"), "Name a load case.");
  model_cls.def("add_combination", &PyModel::add_combination, py::arg("terms"), py::arg("name") = py::none(),
           "A linear combination: `terms` = [(factor, case_id), ...].");
  model_cls.def("validate", &PyModel::validate, "Validate; returns the problems found (empty = sound model).");
  model_cls.def(
          "read_results",
          [](const PyModel& self, const std::string& path) {
            ModelPtr model = self.solver_model();
            auto r = engine([&] {
              StaticResults loaded = load_results_file(path).results;
              // the accessors address the results through this model: same nodes and bars, in order
              auto matches = [&](const CaseResult& c) {
                if (c.node_disp.size() != model->nodes.size() || c.bars.size() != model->bars.size()) return false;
                for (std::size_t i = 0; i < c.bars.size(); ++i)
                  if (c.bars[i].bar != model->bars[i].id) return false;
                return true;
              };
              for (const auto* group : {&loaded.cases, &loaded.combos})
                for (const CaseResult& c : *group)
                  require(matches(c), "read_results: '" + path + "' does not belong to this model (case " +
                                          std::to_string(c.id) + " has other nodes or bars)");
              return std::make_shared<const StaticResults>(std::move(loaded));
            });
            return PyResults{model, r};
          },
          py::arg("path"),
          "The results stored by Results.save, bound to this model and read back bit for bit. ValueError when "
          "the document is not a keta result store, is of a newer schema, or was computed on other nodes or bars.");
  model_cls.def(
          "analyze",
          [](const PyModel& self, int stations, std::optional<std::map<Id, std::vector<double>>> stations_abs) {
            ModelPtr model = self.solver_model();
            StaticOptions o;
            o.stations_per_bar = stations;
            std::map<Id, std::vector<double>> sa;
            if (stations_abs) {
              for (auto [bar, xs] : *stations_abs) {
                (void)bar_of_model(*model, bar);
                std::sort(xs.begin(), xs.end());
                xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
                const double L = bar_length(*model, bar);
                for (double x : xs)
                  if (!(x >= 0 && x <= L * (1 + 1e-12)))
                    throw py::value_error("analyze: station abscissa " + num(x) + " m outside bar " +
                                          std::to_string(bar) + " (length " + num(L) + " m)");
                sa[bar] = xs;
              }
              o.stations_abs = &sa;
            }
            auto r = engine([&] { return std::make_shared<const StaticResults>(solve_static(*model, o)); });
            return PyResults{model, r};
          },
          py::arg("stations") = 11, py::arg("stations_abs") = py::none(),
          "Linear static analysis of every case and combination. stations = evenly spaced study sections per bar; "
          "stations_abs = {bar: [x, ...]} study sections at given abscissae [m] from the bar origin (those bars only).");
  model_cls.def(
          "buckling",
          [](const PyModel& self, std::optional<Id> cs, int n_modes) {
            ModelPtr model = self.solver_model();
            BucklingOptions o;
            o.case_id = cs.value_or(0);
            o.n_modes = n_modes;
            auto r = engine([&] { return std::make_shared<const BucklingResults>(solve_buckling(*model, o)); });
            return PyBuckling{model, r};
          },
          py::arg("case") = py::none(), py::arg("n_modes") = 6,
          "Linear buckling of a reference case (default: the first); factors scale exactly that case.");
  model_cls.def(
          "apply_imperfection",
          [](const PyModel& self, const PyBuckling& b, std::size_t mode, double scale) {
            if (mode >= b.r->modes.size())
              throw py::value_error("imperfection: buckling mode index " + std::to_string(mode) + " out of range (" +
                                    std::to_string(b.r->modes.size()) + " mode(s) available)");
            if (!std::isfinite(scale)) throw py::value_error("imperfection: non-finite imperfection scale " + num(scale));
            if (b.model->nodes.size() != self.m.nodes.size())
              throw py::value_error("imperfection: the buckling result belongs to another model");
            const Model imp = engine([&] { return apply_imperfection(self.m, b.r->modes[mode], scale); });
            PyModel out = self;
            for (std::size_t i = 0; i < out.m.nodes.size(); ++i) out.m.nodes[i].x = imp.nodes[i].x;
            out.m.invalidate_index();
            return out;
          },
          py::arg("buckling"), py::arg("mode") = 0, py::arg("scale") = 1.0,
          "A NEW model whose nodes are offset by scale x the translations of a unit-crest buckling mode.");
  model_cls.def(
          "apply_imperfection_shape",
          [](const PyModel& self, const std::vector<Dof6>& shape, const std::vector<Id>& node_ids, double scale) {
            if (shape.size() != node_ids.size())
              throw py::value_error("imperfection: imperfection shape has " + std::to_string(shape.size()) +
                                    " rows but node_ids has " + std::to_string(node_ids.size()) + " (must match)");
            if (!std::isfinite(scale)) throw py::value_error("imperfection: non-finite imperfection scale " + num(scale));
            std::vector<Vec6> rows(self.m.nodes.size(), Vec6::Zero());
            for (std::size_t k = 0; k < node_ids.size(); ++k) {
              std::size_t i = self.m.nodes.size();
              for (std::size_t j = 0; j < self.m.nodes.size(); ++j)
                if (self.m.nodes[j].id == node_ids[k]) i = j;
              if (i == self.m.nodes.size())
                throw py::value_error("imperfection: imperfection references node " + std::to_string(node_ids[k]) +
                                      ", absent from the model");
              const auto& s = shape[k];
              rows[i] << std::get<0>(s), std::get<1>(s), std::get<2>(s), std::get<3>(s), std::get<4>(s), std::get<5>(s);
            }
            const Model imp = engine([&] { return apply_imperfection(self.m, rows, scale); });
            PyModel out = self;
            for (std::size_t i = 0; i < out.m.nodes.size(); ++i) out.m.nodes[i].x = imp.nodes[i].x;
            out.m.invalidate_index();
            return out;
          },
          py::arg("shape"), py::arg("node_ids"), py::arg("scale") = 1.0,
          "A NEW model offset by scale x an explicit per-node shape [(dx, dy, dz, rx, ry, rz), ...].");
}

void bind_module_model_2(py::module_& mod) {
  auto model_cls = py::reinterpret_borrow<py::class_<PyModel>>(mod.attr("Model"));
  model_cls.def(
          "second_order",
          [](const PyModel& self, std::optional<Id> cs, int stations, int max_iterations, bool exact_stability) {
            ModelPtr model = self.solver_model();
            const Id id = self.case_or_first(cs);
            SecondOrderOptions o;
            o.stations_per_bar = stations;
            o.max_iterations = max_iterations;
            o.exact_stability = exact_stability;
            auto r = engine([&] { return std::make_shared<const SecondOrderResults>(solve_second_order(*model, id, o)); });
            return PySecondOrder{model, r};
          },
          py::arg("case") = py::none(), py::arg("stations") = 11, py::arg("max_iterations") = 30,
          py::arg("exact_stability") = false,
          "Second-order (P-Delta) static analysis of one case (default: the first). exact_stability: the member "
          "tangent is the closed-form beam-column solution (stability functions), exact with one element per member; "
          "default: the consistent geometric stiffness (mesh the members).");
  model_cls.def(
          "nonlinear",
          [](const PyModel& self, std::optional<Id> cs, int substeps, int max_iterations) {
            ModelPtr model = self.solver_model();
            const Id id = self.case_or_first(cs);
            CorotOptions o;
            o.measure = CorotRotationMeasure::Sin;  // sine-form EICR rotations
            o.tol_rel = 1e-10;
            o.substeps = std::max(substeps, 1);
            o.max_iterations = max_iterations;
            auto r = engine([&] { return std::make_shared<const CorotResults>(solve_corotational(*model, id, o)); });
            std::vector<std::pair<Id, Id>> cb;
            for (const auto& c : self.cables) cb.emplace_back(c.id, self.cable_bar_id(c.id));
            return PyNonlinear{model, r, o.substeps, cb, id};
          },
          py::arg("case") = py::none(), py::arg("substeps") = 1, py::arg("max_iterations") = 30,
          "Geometrically nonlinear (co-rotational) static analysis of one case.");
  model_cls.def(
          "arc_length",
          [](const PyModel& self, std::optional<Id> cs, double arc_length, int max_steps, int desired_iterations,
             std::optional<double> stop_lambda, std::optional<double> stop_disp) {
            ModelPtr model = self.solver_model();
            const Id id = self.case_or_first(cs);
            ArcLengthOptions o;
            o.measure = CorotRotationMeasure::Sin;
            o.dl0 = arc_length;
            o.max_steps = max_steps;
            o.max_iterations = 20;
            o.max_cuts = 8;
            o.tol_rel = 1e-10;
            o.desired_iterations = desired_iterations;
            o.stop_lambda = stop_lambda.value_or(0.0);
            o.stop_disp_inf = stop_disp.value_or(0.0);
            auto r = engine([&] { return std::make_shared<const ArcTrace>(arc_length_trace(*model, id, o)); });
            auto dofs = engine([&] { return std::make_shared<const DofMap>(assemble_elements(*model).dofs); });
            return PyArcTrace{model, r, dofs};
          },
          py::arg("case") = py::none(), py::arg("arc_length") = 1.0, py::arg("max_steps") = 100,
          py::arg("desired_iterations") = 0, py::arg("stop_lambda") = py::none(), py::arg("stop_disp") = py::none(),
          "Arc-length continuation through limit points; `arc_length` is the displacement-increment norm per step.");
  model_cls.def(
          "ruin",
          [](const PyModel& self, std::optional<Id> cs) {
            ModelPtr model = self.solver_model();
            const Id id = self.case_or_first(cs);
            RuinOptions o;
            o.corot.measure = CorotRotationMeasure::Sin;
            const RuinResult r = engine([&] { return ruin_load(*model, id, o); });
            std::string oracle;
            if (r.oracle == "indefinite")
              oracle = fmt("tangent indefinite at lambda = %.6f", r.lambda_upper);
            else if (r.oracle == "singular")
              oracle = fmt("singular tangent at lambda = %.6f", r.lambda_upper);
            else if (r.oracle == "snap")
              oracle = fmt("snap-through: tangent indefinite between two equilibria at lambda = %.6f",
                           r.lambda_upper);
            else
              oracle = fmt("Newton diverged at lambda = %.6f", r.lambda_upper);
            return std::make_tuple(r.lambda_lower, r.lambda_upper, oracle);
          },
          py::arg("case") = py::none(), "First-failure load bracket (lambda_lower, lambda_upper, oracle).");
  model_cls.def(
          "influence_surface",
          [](const PyModel& self, const std::vector<Id>& nodes, int nx, int ny, const py::dict& response,
             std::array<double, 3> dir) {
            ModelPtr model = self.solver_model();
            const std::string kind = response.contains("kind") ? response["kind"].cast<std::string>() : "displacement";
            ResponseKind rk;
            if (kind == "displacement")
              rk = ResponseKind::displacement(response["node"].cast<Id>(), response["comp"].cast<int>());
            else if (kind == "reaction")
              rk = ResponseKind::reaction(response["node"].cast<Id>(), response["comp"].cast<int>());
            else if (kind == "bar_force")
              rk = ResponseKind::bar_force(response["bar"].cast<Id>(), response["x"].cast<double>(),
                                           response["comp"].cast<int>());
            else
              throw py::value_error("influence_surface: response kind 'displacement', 'reaction' or 'bar_force'");
            return engine([&] {
              const NodeGrid g = build_node_grid(*model, nodes, nx, ny);
              return influence_surface(*model, g, Vec3(dir[0], dir[1], dir[2]), rk);
            });
          },
          py::arg("nodes"), py::arg("nx"), py::arg("ny"), py::arg("response"), py::arg("dir") = std::array<double, 3>{0, 0, -1},
          "Influence surface of one response (dict kind 'displacement' / 'reaction' (node, comp) or 'bar_force' (bar, x, "
          "comp)) for a unit load of direction `dir` on a plan-rectangular node grid (row-major `nodes`, nx x ny); one "
          "factorization (keta/influence_surface.h).");
  model_cls.def(
          "lm1_surface",
          [](const PyModel& self, const std::vector<Id>& nodes, int nx, int ny, const py::dict& response, double y0,
             double y1, bool minimise, double step_y, const std::string& factors) {
            ModelPtr model = self.solver_model();
            const std::string kind = response.contains("kind") ? response["kind"].cast<std::string>() : "displacement";
            ResponseKind rk;
            if (kind == "displacement")
              rk = ResponseKind::displacement(response["node"].cast<Id>(), response["comp"].cast<int>());
            else if (kind == "reaction")
              rk = ResponseKind::reaction(response["node"].cast<Id>(), response["comp"].cast<int>());
            else if (kind == "bar_force")
              rk = ResponseKind::bar_force(response["bar"].cast<Id>(), response["x"].cast<double>(),
                                           response["comp"].cast<int>());
            else
              throw py::value_error("lm1_surface: response kind 'displacement', 'reaction' or 'bar_force'");
            Lm1SurfaceOptions o;
            o.y0 = y0;
            o.y1 = y1;
            o.minimise = minimise;
            o.step_y = step_y;
            if (factors == "class2")
              o.factors = Lm1Factors::class2();
            else if (factors != "class1")
              throw py::value_error("lm1_surface: factors 'class1' or 'class2'");
            const Lm1SurfaceResult r = engine([&] {
              const NodeGrid g = build_node_grid(*model, nodes, nx, ny);
              const InfluenceSurface srf = influence_surface(*model, g, Vec3(0, 0, -1), rk);
              return lm1_surface_extreme(srf, o);
            });
            py::dict out;
            out["value"] = r.value;
            out["remaining"] = r.remaining;
            py::list lanes;
            for (const Lm1Lane& l : r.lanes) {
              py::dict d;
              d["rank"] = l.rank;
              d["y_start"] = l.y_start;
              d["tandem"] = l.tandem;
              d["ts_x"] = l.ts_x;
              d["ts"] = l.ts;
              d["udl"] = l.udl;
              lanes.append(d);
            }
            out["lanes"] = lanes;
            return out;
          },
          py::arg("nodes"), py::arg("nx"), py::arg("ny"), py::arg("response"), py::arg("y0"), py::arg("y1"),
          py::arg("minimise") = false, py::arg("step_y") = 0.05, py::arg("factors") = "class1",
          "EN 1991-2 LM1 (lanes, tandems, UDL on unfavourable zones) placed for the extreme of one response on the "
          "influence surface of a plan-rectangular node grid (row-major `nodes`, nx x ny), carriageway y0..y1 "
          "(keta/traffic_surface.h). Returns the value and the lanes.");
  model_cls.def(
          "nonlinear_time_history",
          [](const PyModel& self, const py::list& devices, double dt, int steps, const py::list& cases,
             std::optional<int> ground_dir, std::optional<std::pair<std::vector<double>, std::vector<double>>> ground,
             double rayleigh_a, double rayleigh_b, const std::vector<std::pair<Id, int>>& observe) {
            // a piecewise-linear time function from (times, values)
            auto table = [](std::vector<double> ts, std::vector<double> vs) {
              if (ts.size() != vs.size() || ts.empty()) throw py::value_error("time function: times and values of one length");
              return std::function<double(double)>([ts, vs](double t) {
                if (t <= ts.front()) return vs.front();
                if (t >= ts.back()) return vs.back();
                const std::size_t k = static_cast<std::size_t>(std::upper_bound(ts.begin(), ts.end(), t) - ts.begin());
                const double a = (t - ts[k - 1]) / (ts[k] - ts[k - 1]);
                return vs[k - 1] + a * (vs[k] - vs[k - 1]);
              });
            };
            std::vector<LinkDevice> devs;
            for (const py::handle& h : devices) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              auto num = [&](const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
              LinkDevice dv;
              dv.link = d["link"].cast<Id>();
              dv.comp = d.contains("comp") ? d["comp"].cast<int>() : 0;
              const std::string law = d.contains("law") ? d["law"].cast<std::string>() : "none";
              if (law == "lrb")
                dv.law = lrb_law(num("k1", 0), num("Qd", 0), num("k2", 0));
              else if (law == "fps")
                dv.law = fps_law(num("W", 0), num("R", 0), num("mu", 0), num("dy", 1e-3));
              else if (law == "bilinear")
                dv.law = std::make_shared<ElastoPlastic>(num("k", 0), num("fy", 0), num("alpha", 0));
              else if (law != "none")
                throw py::value_error("device law 'lrb', 'fps', 'bilinear' or 'none', got '" + law + "'");
              dv.C = num("C", 0);
              dv.alpha = num("alpha_c", 1);
              dv.v_lin = num("v_lin", 1e-4);
              devs.push_back(dv);
            }
            NlExcitation ex;
            for (const py::handle& h : cases) {
              const auto c = h.cast<std::tuple<Id, std::vector<double>, std::vector<double>>>();
              ex.cases.emplace_back(std::get<0>(c), table(std::get<1>(c), std::get<2>(c)));
            }
            if (ground_dir && ground) {
              ex.ground_dir = *ground_dir;
              ex.ground_accel = table(ground->first, ground->second);
            }
            ModelNlDynOptions o;
            o.dt = dt;
            o.steps = steps;
            o.rayleigh_a = rayleigh_a;
            o.rayleigh_b = rayleigh_b;
            o.observe = observe;
            ModelPtr model = self.solver_model();
            const ModelNlDynResult r = engine([&] { return newmark_model(*model, devs, ex, o); });
            py::dict out;
            out["t"] = r.t;
            out["u"] = r.u;
            out["device_force"] = r.device_force;
            out["device_deformation"] = r.device_deform;
            out["device_work"] = r.device_work;
            out["input_work"] = r.input_work;
            out["kinetic"] = r.kinetic;
            out["strain"] = r.strain;
            out["rayleigh_work"] = r.rayleigh_work;
            return out;
          },
          py::arg("devices"), py::arg("dt"), py::arg("steps"), py::arg("cases") = py::list(),
          py::arg("ground_dir") = py::none(), py::arg("ground") = py::none(), py::arg("rayleigh_a") = 0.0,
          py::arg("rayleigh_b") = 0.0, py::arg("observe") = std::vector<std::pair<Id, int>>{},
          "Nonlinear time history with devices on elastic links: devices = [dict(link, comp, law='lrb' (k1, Qd, k2) "
          "| 'fps' (W, R, mu, dy) | 'bilinear' (k, fy, alpha) | 'none', C, alpha_c, v_lin)], cases = [(case, "
          "times, values)], ground acceleration (times, values) along ground_dir; Newmark + Newton, sparse "
          "(keta/nldyn_model.h). Returns histories and the energy balance terms.");
  model_cls.def(
          "large_displacement_time_history",
          [](const PyModel& self, Id cs, double dt, int steps, std::optional<std::pair<std::vector<double>, std::vector<double>>> f,
             const std::vector<std::pair<Id, int>>& observe, bool lumped, double rayleigh_a, double rayleigh_b,
             const std::map<Id, std::array<double, 6>>& initial_velocity, double rho_inf,
             const PyNonlinear* initial_state, std::optional<int> rupture_step, const std::vector<Id>& rupture_bars,
             const std::vector<Id>& rupture_cables) {
            ModelPtr model = self.solver_model();
            std::function<double(double)> lam = [](double) { return 1.0; };
            if (f) {
              const auto ts = f->first, vs = f->second;
              if (ts.empty() || ts.size() != vs.size()) throw py::value_error("large_displacement_time_history: f = (times, values)");
              lam = [ts, vs](double t) {
                if (t <= ts.front()) return vs.front();
                if (t >= ts.back()) return vs.back();
                std::size_t k = 1;
                while (ts[k] < t) ++k;
                const double a = (t - ts[k - 1]) / (ts[k] - ts[k - 1]);
                return vs[k - 1] + a * (vs[k] - vs[k - 1]);
              };
            }
            CorotDynOptions o;
            o.dt = dt;
            o.steps = steps;
            o.observe = observe;
            o.lumped_mass = lumped;
            o.rayleigh_a = rayleigh_a;
            o.rayleigh_b = rayleigh_b;
            for (const auto& [node, v] : initial_velocity) o.initial_velocity[node] = Eigen::Map<const Vec6>(v.data());
            o.rho_inf = rho_inf;
            if (initial_state) {
              if (initial_state->r->node_disp.size() != model->nodes.size())
                throw py::value_error("large_displacement_time_history: initial_state belongs to another model");
              o.initial_displacement = initial_state->r->node_disp;
            }
            std::vector<Id> cut = rupture_bars;
            for (Id c : rupture_cables) cut.push_back(self.cable_bar_id(c));
            if (!rupture_step && !cut.empty())
              throw py::value_error("large_displacement_time_history: rupture_bars / rupture_cables need rupture_step");
            const CorotDynResult r = engine([&] {
              return rupture_step ? corot_rupture_history(*model, cs, lam, o, *rupture_step, cut)
                                  : corot_newmark(*model, cs, lam, o);
            });
            py::dict out;
            out["t"] = r.t;
            out["u"] = r.u;
            out["iterations"] = r.iterations;
            py::dict fin;
            for (std::size_t i = 0; i < model->nodes.size(); ++i) fin[py::int_(model->nodes[i].id)] = dof6(r.node_disp[i]);
            out["final"] = fin;
            return out;
          },
          py::arg("case") = 1, py::arg("dt") = 0.01, py::arg("steps") = 100, py::arg("f") = py::none(),
          py::arg("observe") = std::vector<std::pair<Id, int>>{}, py::arg("lumped") = true, py::arg("rayleigh_a") = 0.0,
          py::arg("rayleigh_b") = 0.0, py::arg("initial_velocity") = std::map<Id, std::array<double, 6>>{},
          py::arg("rho_inf") = -1.0, py::arg("initial_state") = py::none(), py::arg("rupture_step") = py::none(),
          py::arg("rupture_bars") = std::vector<Id>{}, py::arg("rupture_cables") = std::vector<Id>{},
          "Time history in large displacements (corotational Newmark, keta/nonlinear.h corot_newmark), catenary "
          "cables included (their mass lumped at the ends): the case times f(t) = (times, values) (default 1), from "
          "rest, from initial_state (a Model.nonlinear(...) result: the motion starts hanging in that equilibrium) "
          "and the initial velocities {node: (vx, vy, vz, wx, wy, wz)}; rho_inf in [0, 1] = generalized-alpha "
          "(Chung-Hulbert) with that high-frequency spectral radius, -1 = average-acceleration Newmark. "
          "rupture_step with rupture_bars / rupture_cables: those members vanish at t = rupture_step x dt and the "
          "motion goes on from the state reached. Returns t, the observed (node, dof) histories, the Newton "
          "iterations and the final displacements {node: 6-tuple}.");
  model_cls.def(
          "staged_large_displacements",
          [](const PyModel& self, const py::list& stages) {
            ModelPtr model = self.solver_model();
            std::vector<CorotStage> st;
            for (const py::handle& h : stages) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              CorotStage c;
              if (d.contains("activate")) c.activate = d["activate"].cast<std::vector<Id>>();
              if (d.contains("case")) c.case_id = d["case"].cast<Id>();
              if (d.contains("cable_l0"))
                for (const auto& [k, v] : d["cable_l0"].cast<std::map<Id, double>>()) c.cable_l0.emplace_back(k, v);
              st.push_back(std::move(c));
            }
            CorotOptions o;
            o.tol_rel = 1e-10;
            const CorotStagedResults r = engine([&] { return solve_corotational_staged(*model, st, o); });
            py::list out;
            for (const auto& snap : r.node_disp) {
              py::dict d;
              for (std::size_t i = 0; i < model->nodes.size(); ++i) d[py::int_(model->nodes[i].id)] = dof6(snap[i]);
              out.append(d);
            }
            return out;
          },
          py::arg("stages"),
          "Staged construction in large displacements (keta/nonlinear.h solve_corotational_staged): stages = "
          "[dict(activate=[bar ids], case=id (its loads added), cable_l0={cable: L0})]; new nodes placed stress-free on "
          "the deformed geometry. Returns, per stage, {node: (ux, uy, uz, rx, ry, rz)}.");
  model_cls.def(
          "rail_speed_sweep",
          [](const PyModel& self, const std::vector<Id>& path, const py::object& train, const std::vector<double>& speeds,
             const std::vector<std::pair<Id, int>>& observe, double zeta, int n_modes, double free_time,
             int steps_per_period, double max_travel) {
            ModelPtr model = self.solver_model();
            RailSweepOptions o;
            o.speeds = speeds;
            o.observe = observe;
            o.damping.zeta = zeta;
            o.free_time = free_time;
            o.steps_per_period = steps_per_period;
            o.max_travel = max_travel;
            const Vec3 down(0, 0, -1);
            ModalOptions mo;
            mo.n_modes = n_modes;
            const RailSweepResult r = engine([&]() -> RailSweepResult {
              const ModalResults modes = solve_modal(*model, mo);
              if (py::isinstance<py::str>(train)) {
                const std::string t = train.cast<std::string>();
                require(t.rfind("hslm_a", 0) == 0, "rail_speed_sweep: train 'hslm_a1'..'hslm_a10' or a list");
                return rail_speed_sweep(*model, modes, path, rail_hslm_a(std::stoi(t.substr(6))), down, o);
              }
              const py::list l = train.cast<py::list>();
              bool sprung = false;
              for (const py::handle& h : l) sprung |= py::isinstance<py::dict>(h);
              if (!sprung) {
                ConvoyPattern c;
                for (const py::handle& h : l) c.axles.push_back(h.cast<std::pair<double, double>>());
                return rail_speed_sweep(*model, modes, path, c, down, o);
              }
              VehicleTrain vt;
              for (const py::handle& h : l) {
                const py::dict d = py::reinterpret_borrow<py::dict>(h);
                SprungAxle a;
                a.offset = d["offset"].cast<double>();
                a.m_wheel = d["m_wheel"].cast<double>();
                a.m_body = d["m_body"].cast<double>();
                a.k = d["k"].cast<double>();
                a.c = d.contains("c") ? d["c"].cast<double>() : 0.0;
                vt.axles.push_back(a);
              }
              return rail_speed_sweep_vbi(*model, modes, path, vt, down, o);
            });
            py::dict out;
            py::list pts;
            for (const RailSweepPoint& p : r.points) {
              py::dict d;
              d["speed"] = p.speed;
              d["u_max"] = p.u_max;
              d["u_min"] = p.u_min;
              d["a_max"] = p.a_max;
              pts.append(d);
            }
            out["points"] = pts;
            out["static_max"] = r.static_max;
            out["static_min"] = r.static_min;
            out["dynamic_factor"] = r.dynamic_factor;
            out["critical_speed"] = r.critical_speed;
            return out;
          },
          py::arg("path"), py::arg("train"), py::arg("speeds"), py::arg("observe"), py::arg("zeta") = 0.02,
          py::arg("n_modes") = 10, py::arg("free_time") = 0.0, py::arg("steps_per_period") = 20,
          py::arg("max_travel") = 0.5,
          "Railway speed sweep (vertical loads, downward): the train — 'hslm_a1'..'hslm_a10', [(offset <= 0, load "
          "N)], or sprung axles [dict(offset, m_wheel, m_body, k, c)] (vehicle-bridge interaction) — crossing the bar "
          "path at each speed [m/s]; envelopes of the observed (node, dof) displacements and accelerations, static "
          "envelope, dynamic factor and critical speed. Exact modal moving loads (keta/rail_dynamics.h).");
  model_cls.def(
          "modal",
          [](const PyModel& self, int n_modes, std::string method, std::optional<double> f_max,
             std::optional<double> mass_ratio, int increment) {
            for (char& c : method) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            ModalExtOptions o;
            o.n_modes = n_modes;
            if (method == "dense")
              o.dense_threshold = std::numeric_limits<int>::max();
            else if (method == "lanczos")
              o.dense_threshold = 0;
            else if (method != "auto")
              throw py::value_error("modal: unknown method \"" + method + "\"; use \"auto\", \"dense\" or \"lanczos\"");
            if (f_max && mass_ratio) throw py::value_error("modal: give f_max= or mass_ratio=, not both");
            ModelPtr model = self.solver_model();
            if (f_max) {
              if (!(*f_max > 0)) throw py::value_error("modal: f_max must be > 0 [Hz]");
              auto fr = engine([&] { return solve_modal_frequency(*model, *f_max, o); });
              if (!fr.complete) PyErr_WarnEx(PyExc_RuntimeWarning, fr.warning.c_str(), 1);
              return PyModal{model, std::make_shared<const ModalResults>(std::move(fr.modal))};
            }
            if (mass_ratio) {
              if (!(*mass_ratio > 0 && *mass_ratio <= 1)) throw py::value_error("modal: mass_ratio in (0, 1]");
              if (increment < 1) throw py::value_error("modal: increment must be >= 1");
              auto rr = engine([&] { return solve_modal_ratio(*model, 100.0 * *mass_ratio, increment, o); });
              if (!rr.target_reached) PyErr_WarnEx(PyExc_RuntimeWarning, rr.warning.c_str(), 1);
              return PyModal{model, std::make_shared<const ModalResults>(std::move(rr.modal))};
            }
            auto r = engine([&] { return std::make_shared<const ModalResults>(solve_modal(*model, o)); });
            if (!r->warning.empty()) PyErr_WarnEx(PyExc_RuntimeWarning, r->warning.c_str(), 1);
            return PyModal{model, r};
          },
          py::arg("n_modes") = 6, py::arg("method") = "auto", py::arg("f_max") = py::none(),
          py::arg("mass_ratio") = py::none(), py::arg("increment") = 10,
          "Free-vibration frequencies and mode shapes: n_modes modes; or every mode below f_max [Hz] (Sturm count, "
          "checked); or enough modes for every translational direction carrying mass to reach mass_ratio (0.9 = 90 %) "
          "of its total mass, the count grown by `increment`. The lowest modes are certified complete by a Sturm count "
          "(repeated frequencies included); a free body (no supports) is solved with a spectral shift and its rigid "
          "modes come first, at 0 Hz. A count or a target not reached raises a RuntimeWarning.");
  model_cls.def(
          "modal_substructured",
          [](const PyModel& self, const py::list& subs, std::vector<Id> masters, int n_modes) {
            std::vector<Substructure> ss;
            for (const py::handle& h : subs) {
              const py::dict d = py::reinterpret_borrow<py::dict>(h);
              Substructure s;
              if (d.contains("bars")) s.bars = d["bars"].cast<std::vector<Id>>();
              if (d.contains("shells")) s.shells = d["shells"].cast<std::vector<Id>>();
              if (d.contains("n_modes")) s.n_modes = d["n_modes"].cast<int>();
              ss.push_back(std::move(s));
            }
            ModalOptions o;
            o.n_modes = n_modes;
            ModelPtr model = self.solver_model();
            auto r = engine([&] { return std::make_shared<const SubstructuredModal>(solve_modal_substructured(*model, ss, masters, o)); });
            py::dict out;
            out["modal"] = PyModal{model, std::make_shared<const ModalResults>(r->modal)};
            out["n_boundary"] = r->n_boundary;
            out["n_reduced"] = r->n_reduced;
            py::list fi;
            for (const SubstructureInfo& i : r->subs) fi.append(i.fixed_interface_hz);
            out["fixed_interface_hz"] = fi;
            return out;
          },
          py::arg("substructures"), py::arg("masters") = std::vector<Id>{}, py::arg("n_modes") = 6,
          "Modal analysis by superelements: substructures [dict(bars=[...], shells=[...], n_modes=0)] condensed "
          "on their boundary (Guyan, n_modes = 0) plus fixed-interface modes (Craig-Bampton); masters = nodes kept. "
          "Returns dict(modal, n_boundary, n_reduced, fixed_interface_hz) (keta/substructure.h).");
  model_cls.def(
          "modal_rotating",
          [](const PyModel& self, double omega, std::array<double, 3> axis, std::array<double, 3> point, int n_modes,
             bool stiffening, bool softening, bool lumped) {
            Spin sp;
            sp.omega = omega;
            sp.axis = Vec3(axis[0], axis[1], axis[2]);
            sp.point = Vec3(point[0], point[1], point[2]);
            sp.stiffening = stiffening;
            sp.softening = softening;
            ModalOptions o;
            o.n_modes = n_modes;
            o.lumped_mass = lumped;
            ModelPtr model = self.solver_model();
            auto r = engine([&] { return std::make_shared<const ModalResults>(solve_modal_rotating(*model, sp, o).modal); });
            return PyModal{model, r};
          },
          py::arg("omega"), py::arg("axis") = std::array<double, 3>{1, 0, 0}, py::arg("point") = std::array<double, 3>{0, 0, 0},
          py::arg("n_modes") = 6, py::arg("stiffening") = true, py::arg("softening") = true, py::arg("lumped") = false,
          "Modes of the structure spinning at omega [rad/s] about the axis (point, direction): exact centrifugal "
          "prestress and spin softening -Omega^2 M_perp (keta/rotating.h).");
  model_cls.def(
          "harmonic",
          [](const PyModel& self, std::size_t cs, std::vector<double> freqs, std::optional<double> zeta,
             std::optional<std::vector<double>> zetas, std::optional<std::pair<double, double>> rayleigh, int n_modes,
             int stations_per_bar, std::optional<Id> case_id) {
            if (freqs.empty()) throw py::value_error("harmonic: freqs must be a non-empty frequency sweep (Hz)");
            const DampingSpec d = parse_damping(zeta, zetas, rayleigh, 0.02);
            ModelPtr model = self.solver_model();
            if (case_id) cs = case_index_of(*model, *case_id);
            if (cs >= model->cases.size())
              throw py::value_error("harmonic: load-case index " + std::to_string(cs) + " is out of range (the model has " +
                                    std::to_string(model->cases.size()) + " cases)");
            const Id id = model->cases[cs].id;
            auto r = engine([&] {
              return std::make_shared<const HarmonicResults>(solve_harmonic(*model, id, freqs, d, n_modes, stations_per_bar));
            });
            return PyHarmonic{model, r, bar_ids_of(*model)};
          },
          py::arg("case") = 0, py::arg("freqs") = std::vector<double>{}, py::arg("zeta") = py::none(),
          py::arg("zetas") = py::none(), py::arg("rayleigh") = py::none(), py::arg("n_modes") = 0,
          py::arg("stations_per_bar") = 11, py::arg("case_id") = py::none(),
          "Steady-state harmonic response of one case over a frequency sweep (Hz): `case` is the DENSE 0-based "
          "index, `case_id` (preferred) the load-case number.");
  model_cls.def(
          "psd",
          [](const PyModel& self, std::size_t cs, std::vector<double> freqs, std::optional<double> zeta,
             std::optional<std::vector<double>> zetas, std::optional<std::pair<double, double>> rayleigh, int n_modes,
             int stations_per_bar, std::optional<double> psd, std::optional<std::vector<std::pair<double, double>>> table,
             std::optional<Id> case_id) {
            if (freqs.empty()) throw py::value_error("psd: freqs must be a non-empty frequency sweep (Hz)");
            if (psd && table)
              throw py::value_error("psd: give at most one input-PSD kwarg: psd= (white level) or psd_table=");
            const DampingSpec d = parse_damping(zeta, zetas, rayleigh, 0.02);
            ModelPtr model = self.solver_model();
            if (case_id) cs = case_index_of(*model, *case_id);
            if (cs >= model->cases.size())
              throw py::value_error("psd: load-case index " + std::to_string(cs) + " is out of range (the model has " +
                                    std::to_string(model->cases.size()) + " cases)");
            // Input PSD: the white level, or the table linearly interpolated
            // and clamped at both ends.
            std::vector<double> s(freqs.size(), psd.value_or(1.0));
            if (table) {
              const auto& t = *table;
              for (std::size_t i = 0; i < freqs.size(); ++i) {
                const double f = freqs[i];
                if (t.empty()) {
                  s[i] = 0;
                } else if (f <= t.front().first) {
                  s[i] = t.front().second;
                } else if (f >= t.back().first) {
                  s[i] = t.back().second;
                } else {
                  for (std::size_t k = 1; k < t.size(); ++k)
                    if (f <= t[k].first) {
                      const double u = (f - t[k - 1].first) / (t[k].first - t[k - 1].first);
                      s[i] = t[k - 1].second + u * (t[k].second - t[k - 1].second);
                      break;
                    }
                }
              }
            }
            const Id id = model->cases[cs].id;
            auto r = engine([&] {
              return std::make_shared<const PsdResults>(solve_psd(*model, id, freqs, s, d, n_modes, stations_per_bar));
            });
            auto h = engine([&] {
              return std::make_shared<const HarmonicResults>(solve_harmonic(*model, id, freqs, d, n_modes, stations_per_bar));
            });
            return PyPsd{model, r, h, s, bar_ids_of(*model)};
          },
          py::arg("case") = 0, py::arg("freqs") = std::vector<double>{}, py::arg("zeta") = py::none(),
          py::arg("zetas") = py::none(), py::arg("rayleigh") = py::none(), py::arg("n_modes") = 0,
          py::arg("stations_per_bar") = 11, py::arg("psd") = py::none(), py::arg("psd_table") = py::none(),
          py::arg("case_id") = py::none(),
          "Random response to a fully correlated pattern F(t) = alpha(t) P with input PSD S_alpha; `case` is the "
          "DENSE 0-based index, `case_id` (preferred) the load-case number.");
  model_cls.def(
          "response_spectrum",
          [](const PyModel& self, std::vector<std::pair<double, double>> spectrum, int dir, std::optional<double> zeta,
             std::optional<std::vector<double>> zetas, std::optional<std::pair<double, double>> rayleigh,
             std::string rule, int n_modes, double table_zeta, std::optional<Vec3T> direction, double duration) {
            const std::string bad = "bad analysis options: response spectrum: ";
            if (spectrum.empty()) throw py::value_error(bad + "the spectrum table needs at least one (period, Sa) point");
            for (const auto& [t, sa] : spectrum) {
              if (!std::isfinite(t) || t < 0) throw py::value_error(bad + "every period must be finite and non-negative");
              if (!std::isfinite(sa) || sa < 0) throw py::value_error(bad + "every Sa ordinate must be finite and non-negative");
            }
            for (std::size_t i = 1; i < spectrum.size(); ++i)
              if (spectrum[i].first <= spectrum[i - 1].first)
                throw py::value_error(bad + "periods must be strictly increasing");
            const DampingSpec d = parse_damping(zeta, zetas, rayleigh, 0.05);
            for (char& c : rule) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            SpectrumOptions o;
            if (rule == "srss")
              o.rule = ComboRule::Srss;
            else if (rule == "cqc")
              o.rule = ComboRule::Cqc;
            else if (rule == "10%" || rule == "tenpercent")
              o.rule = ComboRule::TenPercent;
            else if (rule == "abs")
              o.rule = ComboRule::Abs;
            else if (rule == "nrl")
              o.rule = ComboRule::Nrl;
            else if (rule == "double_sum") {
              o.rule = ComboRule::DoubleSum;
              o.duration = duration;
            } else
              throw py::value_error("rule must be 'srss', 'cqc', '10%', 'double_sum', 'abs' or 'nrl', got '" + rule + "'");
            if (dir < 0 || dir > 5)
              throw py::value_error("excitation direction " + std::to_string(dir) + " out of range (0..5)");
            ModelPtr model = self.solver_model();
            ResponseSpectrum sp;
            sp.kind = ResponseSpectrum::Kind::Table;
            sp.table = spectrum;
            sp.table_zeta = table_zeta;
            o.dir = dir;
            if (direction) o.dir_vector << std::get<0>(*direction), std::get<1>(*direction), std::get<2>(*direction), 0, 0, 0;
            const double total = engine([&] {
              ModalOptions mo;
              mo.n_modes = 1;
              return solve_modal(*model, mo).total_mass[static_cast<std::size_t>(dir)];
            });
            auto r = engine([&] {
              o.n_modes = n_modes > 0 ? n_modes : assemble_elements(*model, true).dofs.n_free;
              return std::make_shared<const SpectrumResults>(solve_spectrum(*model, sp, d, o));
            });
            return PySpectrum{model, r, dir, total};
          },
          py::arg("spectrum") = std::vector<std::pair<double, double>>{}, py::arg("dir") = 0,
          py::arg("zeta") = py::none(), py::arg("zetas") = py::none(), py::arg("rayleigh") = py::none(),
          py::arg("rule") = "srss", py::arg("n_modes") = 0, py::arg("table_zeta") = 0.0,
          py::arg("direction") = py::none(), py::arg("duration") = 20.0,
          "Response-spectrum analysis: [(period_s, Sa), ...] in one direction 0..5 (or along `direction`, a "
          "vector of translations); modal rule 'srss', 'cqc', '10%', 'double_sum' (strong-motion duration "
          "`duration` [s]), 'abs' or 'nrl'; table_zeta > 0: the table's damping, corrected to the modal "
          "damping by the EN 1998-1 eta.");
  model_cls.def(
          "cable_static",
          [](const PyModel& self, std::optional<Id> cs, int stations, int substeps) {
            if (self.cables.empty()) throw py::value_error("cable_static: the model carries no cable");
            ModelPtr model = self.solver_model();
            Id id = 0;
            if (!cs && model->cases.empty()) {  // cables under their own weight alone: an empty case
              auto mm = std::make_shared<Model>(*model);
              LoadCase lc;
              lc.id = 1;
              lc.name = "cables";
              mm->cases.push_back(lc);
              mm->invalidate_index();
              model = mm;
              id = 1;
            } else {
              id = self.case_or_first(cs);
            }
            CableStaticOptions o;
            o.stations_per_bar = stations;
            o.substeps = std::max(substeps, 1);
            auto r = engine([&] { return std::make_shared<const CableStaticResults>(solve_cable_static(*model, id, o)); });
            PyCableStatic out{model, r, {}, engine([&] { return resolve_cable_l0(*model); }), id};
            for (const CableRec& c : self.cables)
              out.cables.emplace_back(c.id, static_cast<std::size_t>(model->bar_index(self.cable_bar_id(c.id))));
            return out;
          },
          py::arg("case") = py::none(), py::arg("stations") = 11, py::arg("substeps") = 1,
          "Coupled elastic-catenary cables + linear frame statics (non-superposable).");
  model_cls.def(
          "export_script",
          [](const PyModel& self, const PyResults* results, int stations) {
            ScriptOptions o;
            o.stations_per_bar = stations;
            return engine([&] { return export_script(*self.solver_model(), results ? results->r.get() : nullptr, o); });
          },
          py::arg("results") = py::none(), py::arg("stations") = 11,
          "The model as a Python script re-solving it (with `results`: asserting the same numbers).");
  model_cls.def("to_json", [](const PyModel& self) { return model_to_json(*self.solver_model()).dump(1); },
           "The model as keta model JSON.");
  model_cls.def_static(
          "from_json",
          [](const std::string& text) {
            return engine([&] { return PyModel::from_model(model_from_json(json::Value::parse(text))); });
          },
          py::arg("text"), "Rebuild a model from to_json output.");
  model_cls.def_property_readonly("node_count", [](const PyModel& s) { return s.m.nodes.size(); });
  model_cls.def_property_readonly("bar_count", [](const PyModel& s) { return s.m.bars.size(); });
  model_cls.def("__repr__", &PyModel::repr);
}

void bind_module_statics(py::module_& mod) {
  // ---------------- static results ----------------
  py::class_<PyResults>(mod, "Results", "All analysis results (every case and combination).")
      .def(
          "save",
          [](const PyResults& s, const std::string& path, const std::string& note) {
            engine([&] {
              save_json_file(path, results_to_json(*s.r, {"", note}));
              return 0;
            });
          },
          py::arg("path"), py::arg("note") = "",
          "Write the results to a JSON store (schema version, provenance with this note, sign convention, every "
          "case and combination); Model.read_results reads it back bit for bit.")
      .def(
          "case",
          [](const PyResults& s, Id id) {
            for (std::size_t i = 0; i < s.r->cases.size(); ++i)
              if (s.r->cases[i].id == id) return PyCaseResult{s.model, s.r, i, false};
            throw py::key_error("no load case " + std::to_string(id));
          },
          py::arg("id"))
      .def(
          "combo",
          [](const PyResults& s, Id id) {
            for (std::size_t i = 0; i < s.r->combos.size(); ++i)
              if (s.r->combos[i].id == id) return PyCaseResult{s.model, s.r, i, true};
            throw py::key_error("no combination " + std::to_string(id));
          },
          py::arg("id"))
      .def_property_readonly("case_ids",
                             [](const PyResults& s) {
                               std::vector<Id> v;
                               for (const auto& c : s.r->cases) v.push_back(c.id);
                               return v;
                             })
      .def_property_readonly("combo_ids",
                             [](const PyResults& s) {
                               std::vector<Id> v;
                               for (const auto& c : s.r->combos) v.push_back(c.id);
                               return v;
                             })
      .def("__repr__", [](const PyResults& s) {
        return "<keta.Results: " + std::to_string(s.r->cases.size()) + " case(s), " + std::to_string(s.r->combos.size()) +
               " combo(s)>";
      });

  {
    py::class_<PyCaseResult> cr(mod, "CaseResult", "One case's (or combination's) results.");
    bind_case_accessors(cr, [](const PyCaseResult& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, s.get()); },
                        "end_forces");
    cr.def("__repr__", [](const PyCaseResult& s) {
      const CaseResult& c = s.get();
      return "<keta.CaseResult '" + (c.name.empty() ? std::to_string(c.id) : c.name) + "': " +
             std::to_string(c.node_disp.size()) + " nodes, " + std::to_string(c.bars.size()) + " bars>";
    });
  }

  // ---------------- buckling ----------------
  py::class_<PyBuckling>(mod, "Buckling", "Linear buckling: critical load factors and unit-crest modes.")
      .def_property_readonly("factors", [](const PyBuckling& s) { return s.r->load_factors; })
      .def_property_readonly("n_modes", [](const PyBuckling& s) { return s.r->modes.size(); })
      .def_property_readonly("node_ids", [](const PyBuckling& s) { return node_ids_of(*s.model); })
      .def(
          "mode_displacement",
          [](const PyBuckling& s, std::size_t mode, Id node) {
            if (mode >= s.r->modes.size()) throw py::key_error("no buckling mode " + std::to_string(mode));
            return dof6(s.r->modes[mode][node_pos(*s.model, node)]);
          },
          py::arg("mode"), py::arg("node"))
      .def("mode", [](const PyBuckling& s, std::size_t k) { return mode_rows(*s.model, s.r->modes, k, "buckling mode "); },
           py::arg("mode"))
      .def("__repr__", [](const PyBuckling& s) {
        if (s.r->load_factors.empty()) return std::string("<keta.Buckling: no instability (0 modes)>");
        return fmt("<keta.Buckling: lambda_1 = %.6f", s.r->load_factors[0]) + ", " + std::to_string(s.r->modes.size()) +
               " mode(s)>";
      });

  // ---------------- second order ----------------
  {
    py::class_<PySecondOrder> so(mod, "SecondOrder", "Second-order (P-Delta) results of one case.");
    so.def_property_readonly("iterations", [](const PySecondOrder& s) { return s.r->iterations; });
    bind_case_accessors(
        so, [](const PySecondOrder& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, s.r->result); },
        "bar_end_forces");
    so.def("__repr__", [](const PySecondOrder& s) {
      const CaseResult& c = s.r->result;
      return "<keta.SecondOrder '" + (c.name.empty() ? std::to_string(c.id) : c.name) + "': " +
             std::to_string(s.r->iterations) + " iteration(s), " + std::to_string(c.node_disp.size()) + " nodes, " +
             std::to_string(c.bars.size()) + " bars>";
    });
  }

  {
    py::class_<PyLinkLaw> ll(mod, "LinkLawResult", "Static state of one case with non-linear link laws.");
    ll.def_property_readonly("iterations", [](const PyLinkLaw& s) { return s.r->iterations; });
    bind_case_accessors(
        ll, [](const PyLinkLaw& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, s.r->result); },
        "bar_end_forces");
  }

  // ---------------- nonlinear ----------------
  py::class_<PyNonlinear>(mod, "Nonlinear", "Converged large-displacement state of one case.")
      .def_property_readonly("iterations", [](const PyNonlinear& s) { return s.r->iterations; })
      .def(
          "cable_tension",
          [](const PyNonlinear& s, Id cable) {
            for (const auto& [cid, bid] : s.cable_bars)
              if (cid == cable) {
                const auto& t = s.r->cable_T[static_cast<std::size_t>(s.model->bar_index(bid))];
                return std::make_tuple(t[0], t[1]);
              }
            throw py::key_error("no cable " + std::to_string(cable));
          },
          py::arg("cable"), "End tensions (end 1, end 2) of a cable in the deformed state (exact catenary).")
      .def_property_readonly("substeps", [](const PyNonlinear& s) { return s.substeps; })
      .def(
          "displacement", [](const PyNonlinear& s, Id node) { return dof6(s.r->node_disp.at(node_pos(*s.model, node))); },
          py::arg("node"))
      .def("reactions",
           [](const PyNonlinear& s) {
             std::vector<NodeDof6> out;
             for (const Reaction& r : s.r->reactions) out.emplace_back(r.node, dof6(r.R));
             return out;
           })
      .def(
          "bar_end_forces",
          [](const PyNonlinear& s, std::size_t bar_index) {
            if (bar_index >= s.r->bar_forces.size()) throw py::key_error("no bar at index " + std::to_string(bar_index));
            return vec12(s.r->bar_forces[bar_index]);
          },
          py::arg("bar_index"), "Member end forces of the bar at 0-based index `bar_index` (see end_forces).")
      .def(
          "end_forces",
          [](const PyNonlinear& s, Id bar) {
            const int ix = s.model->bar_index(bar);
            if (ix < 0 || static_cast<std::size_t>(ix) >= s.r->bar_forces.size())
              throw py::key_error("no bar " + std::to_string(bar));
            return vec12(s.r->bar_forces[static_cast<std::size_t>(ix)]);
          },
          py::arg("bar"), "Member end forces of bar `bar` (its id, as every other accessor).")
      .def("max_displacement", [](const PyNonlinear& s) { return max_translation(s.r->node_disp); })
      .def("__repr__", [](const PyNonlinear& s) {
        return "<keta.Nonlinear: " + std::to_string(s.substeps) + " substep(s), " + std::to_string(s.r->iterations) +
               " iteration(s), " + std::to_string(s.model->nodes.size()) + " nodes>";
      });

  // ---------------- arc length ----------------
  py::class_<PyArcTrace>(mod, "ArcTrace", "An arc-length equilibrium path.")
      .def_property_readonly("n_points", [](const PyArcTrace& s) { return s.r->points.size(); })
      .def("lambdas",
           [](const PyArcTrace& s) {
             std::vector<double> v;
             for (const ArcPoint& p : s.r->points) v.push_back(p.lambda);
             return v;
           })
      .def(
          "n_neg",
          [](const PyArcTrace& s, std::size_t k) {
            if (k >= s.r->points.size()) throw py::key_error("no path point " + std::to_string(k));
            return s.r->points[k].n_neg;
          },
          py::arg("point"))
      .def("crossings",
           [](const PyArcTrace& s) {
             std::vector<std::size_t> v;
             for (std::size_t k = 0; k < s.r->points.size(); ++k)
               if (s.r->points[k].crossing) v.push_back(k);
             return v;
           })
      .def(
          "node_disp",
          [](const PyArcTrace& s, std::size_t k, Id node) {
            if (k >= s.r->points.size()) throw py::key_error("no path point " + std::to_string(k));
            return dof6(expand(*s.dofs, s.r->points[k].u).at(node_pos(*s.model, node)));
          },
          py::arg("point"), py::arg("node"))
      .def(
          "max_disp",
          [](const PyArcTrace& s, std::size_t k) {
            if (k >= s.r->points.size()) throw py::key_error("no path point " + std::to_string(k));
            return max_translation(expand(*s.dofs, s.r->points[k].u));
          },
          py::arg("point"))
      .def("__repr__", [](const PyArcTrace& s) {
        double lo = 0, hi = 0;
        bool first = true;
        std::size_t crossings = 0;
        for (const ArcPoint& p : s.r->points) {
          lo = first ? p.lambda : std::min(lo, p.lambda);
          hi = first ? p.lambda : std::max(hi, p.lambda);
          first = false;
          crossings += p.crossing;
        }
        return "<keta.ArcTrace: " + std::to_string(s.r->points.size()) + " points, lambda in [" + fmt("%.4f", lo) +
               ", " + fmt("%.4f", hi) + "], " + std::to_string(crossings) + " crossing(s)>";
      });

  // ---------------- modal ----------------
  {
    auto ext = [](const SurfaceExtremes& e) {
      py::dict d;
      d["max"] = e.max;
      d["min"] = e.min;
      d["at_max"] = std::make_pair(e.rx_max, e.ry_max);
      d["at_min"] = std::make_pair(e.rx_min, e.ry_min);
      return d;
    };
    auto wheels = [](const std::vector<std::array<double, 3>>& w) {
      WheelSet ws;
      ws.wheels = w;
      return ws;
    };
    py::class_<InfluenceSurface>(mod, "InfluenceSurface", "Influence surface on a plan-rectangular grid (bilinear).")
        .def_readonly("x", &InfluenceSurface::x)
        .def_readonly("y", &InfluenceSurface::y)
        .def_property_readonly("eta",
                               [](const InfluenceSurface& s) {
                                 std::vector<std::vector<double>> rows(static_cast<std::size_t>(s.ny));
                                 for (int j = 0; j < s.ny; ++j)
                                   for (int i = 0; i < s.nx; ++i) rows[static_cast<std::size_t>(j)].push_back(s.node_eta(i, j));
                                 return rows;
                               })
        .def("at", &InfluenceSurface::at, py::arg("x"), py::arg("y"), "Bilinear ordinate (0 outside the grid).")
        .def("integral", &InfluenceSurface::integral)
        .def("integral_pos", &InfluenceSurface::integral_pos)
        .def("integral_neg", &InfluenceSurface::integral_neg)
        .def("integral_box", &InfluenceSurface::integral_box, py::arg("x0"), py::arg("x1"), py::arg("y0"), py::arg("y1"),
             "Exact integral over a box (clipped to the grid).")
        .def("point_extreme", &InfluenceSurface::point_extreme, "(max, min) of a single unit load.")
        .def(
            "patch_extreme",
            [ext](const InfluenceSurface& s, double a, double b, double p, int sx, int sy) {
              return ext(engine([&] { return patch_extreme(s, a, b, p, sx, sy); }));
            },
            py::arg("a"), py::arg("b"), py::arg("p"), py::arg("steps_x") = 60, py::arg("steps_y") = 60)
        .def(
            "udl_extreme", [ext](const InfluenceSurface& s, double p) { return ext(udl_surface_extreme(s, p)); }, py::arg("p"),
            "A free UDL on the favourable zones: p x integral_pos / integral_neg.")
        .def(
            "wheel_set_extreme",
            [ext, wheels](const InfluenceSurface& s, const std::vector<std::array<double, 3>>& w, int sx, int sy) {
              return ext(engine([&] { return wheel_set_extreme(s, wheels(w), sx, sy); }));
            },
            py::arg("wheels"), py::arg("steps_x") = 60, py::arg("steps_y") = 60,
            "Slide a wheel set [(dx, dy, P)] over the grid (wheels off the deck count zero).")
        .def(
            "vehicle_extreme",
            [](const InfluenceSurface& s, const std::vector<std::array<double, 3>>& w, double y0, double y1, bool both,
               bool minimise) {
              VehicleBand b;
              b.y0 = y0;
              b.y1 = y1;
              b.both_directions = both;
              b.minimise = minimise;
              const VehiclePlacement p = engine([&] { return vehicle_surface_extreme(s, w, b); });
              py::dict d;
              d["value"] = p.value;
              d["x"] = p.x;
              d["y"] = p.y;
              d["reversed"] = p.reversed;
              return d;
            },
            py::arg("wheels"), py::arg("y0") = -1e300, py::arg("y1") = 1e300, py::arg("both_directions") = true,
            py::arg("minimise") = false,
            "EXACT extreme of a vehicle [(dx, dy, P)] whose wheels stay within [y0, y1] (keta/traffic_groups.h).")
        .def(
            "load_group",
            [](const InfluenceSurface& s, const std::string& group, double y0, double y1, bool minimise, const py::object& factors,
               const std::vector<std::pair<double, double>>& footways, const std::vector<std::array<double, 3>>& vehicle,
               double clear, double step_y, const std::vector<InfluenceSurface>& concomitant) {
              GroupSpec g;
              Lm1SurfaceOptions& o = g.lm1;
              o.y0 = y0;
              o.y1 = y1;
              o.minimise = minimise;
              o.step_y = step_y;
              if (py::isinstance<py::str>(factors)) {
                const std::string f = factors.cast<std::string>();
                if (f == "class2") o.factors = Lm1Factors::class2();
                else if (f != "class1") throw py::value_error("load_group: factors 'class1', 'class2' or dict(aQ, aq1, aq_other)");
              } else {
                const py::dict d = factors.cast<py::dict>();
                const auto aQ = d["aQ"].cast<std::array<double, 3>>();
                for (int k = 0; k < 3; ++k) o.factors.aQ[k] = aQ[static_cast<std::size_t>(k)];
                o.factors.aq1 = d["aq1"].cast<double>();
                o.factors.aq_other = d["aq_other"].cast<double>();
              }
              if (group == "gr1a") {
                g.group = LoadGroup::Gr1a;
                for (const auto& [a, b] : footways) g.footways.push_back({a, b});
              } else if (group == "gr4") {
                g.group = LoadGroup::Gr4;
                g.y0 = footways.empty() ? y0 : std::min(y0, footways.front().first);
                g.y1 = footways.empty() ? y1 : std::max(y1, footways.back().second);
              } else if (group == "gr5") {
                g.group = LoadGroup::Gr5;
                g.vehicle = vehicle;
                g.band.y0 = y0;
                g.band.y1 = y1;
                g.band.minimise = minimise;
                g.clear = clear;
              } else {
                throw py::value_error("load_group: 'gr1a', 'gr4' or 'gr5'");
              }
              const GroupResult r = engine([&] { return group_extreme(s, g); });
              py::dict d;
              d["value"] = r.value;
              d["lm1"] = r.lm1.value;
              d["footways"] = r.footways;
              d["vehicle"] = r.vehicle.value;
              d["vehicle_at"] = std::make_tuple(r.vehicle.x, r.vehicle.y, r.vehicle.reversed);
              std::vector<double> cc;
              for (const InfluenceSurface& t : concomitant) cc.push_back(engine([&] { return group_concomitant(t, s, r, g); }));
              d["concomitant"] = cc;
              return d;
            },
            py::arg("group"), py::arg("y0"), py::arg("y1"), py::arg("minimise") = false, py::arg("factors") = py::str("class1"),
            py::arg("footways") = std::vector<std::pair<double, double>>{}, py::arg("vehicle") = std::vector<std::array<double, 3>>{},
            py::arg("clear") = 0.0, py::arg("step_y") = 0.05, py::arg("concomitant") = std::vector<InfluenceSurface>{},
            "EN 1991-2 load group on the surface (carriageway y0..y1): 'gr1a' (LM1 + footways [(y0, y1)] at 3 kN/m2), "
            "'gr4' (crowd 5 kN/m2), 'gr5' (special vehicle [(dx, dy, P)] + frequent LM1 kept `clear` m beyond its axles, "
            "exact branch and bound); factors 'class1', 'class2' or dict(aQ=[..], aq1, aq_other) for other annexes; "
            "concomitant = [surfaces of other responses on the same grid]: their values under the arrangement found.")
        .def(
            "wheel_set_value",
            [wheels](const InfluenceSurface& s, const std::vector<std::array<double, 3>>& w, double rx, double ry) {
              return wheel_set_value(s, wheels(w), rx, ry);
            },
            py::arg("wheels"), py::arg("rx"), py::arg("ry"));
  }
  py::class_<PyModal>(mod, "Modal", "Natural frequencies and mode shapes.")
      .def_property_readonly("frequencies", [](const PyModal& s) { return s.r->frequencies_hz; })
      .def_property_readonly("omegas", [](const PyModal& s) { return s.r->omegas; })
      .def_property_readonly("periods",
                             [](const PyModal& s) {
                               std::vector<double> v;
                               for (double f : s.r->frequencies_hz)
                                 v.push_back(f > 0 ? 1.0 / f : std::numeric_limits<double>::infinity());
                               return v;
                             })
      .def_property_readonly("n_modes", [](const PyModal& s) { return s.r->modes.size(); })
      .def_property_readonly("node_ids", [](const PyModal& s) { return node_ids_of(*s.model); })
      .def(
          "mode_displacement",
          [](const PyModal& s, std::size_t mode, Id node) {
            if (mode >= s.r->modes.size()) throw py::key_error("no mode " + std::to_string(mode));
            return dof6(s.r->modes[mode][node_pos(*s.model, node)]);
          },
          py::arg("mode"), py::arg("node"))
      .def("mode", [](const PyModal& s, std::size_t k) { return mode_rows(*s.model, s.r->modes, k, "mode "); },
           py::arg("mode"))
      .def("participation",
           [](const PyModal& s) {
             std::vector<Dof6> v;
             for (const auto& g : s.r->participation) v.push_back({g[0], g[1], g[2], g[3], g[4], g[5]});
             return v;
           })
      .def("effective_mass",
           [](const PyModal& s) {
             std::vector<Dof6> v;
             for (const auto& e : s.r->effective_mass) v.push_back({e[0], e[1], e[2], e[3], e[4], e[5]});
             return v;
           })
      .def("total_mass",
           [](const PyModal& s) {
             const auto& t = s.r->total_mass;
             return Dof6{t[0], t[1], t[2], t[3], t[4], t[5]};
           })
      .def("mass_ratio",
           [](const PyModal& s) {
             const auto& t = s.r->total_mass;
             std::vector<Dof6> v;
             for (const auto& e : s.r->effective_mass) {
               auto r = [&](int d) {
                 return t[static_cast<std::size_t>(d)] != 0 ? e[static_cast<std::size_t>(d)] / t[static_cast<std::size_t>(d)]
                                                            : 0.0;
               };
               v.push_back({r(0), r(1), r(2), r(3), r(4), r(5)});
             }
             return v;
           })
      .def("__repr__", [](const PyModal& s) {
        if (s.r->frequencies_hz.empty()) return std::string("<keta.Modal: no modes>");
        return fmt("<keta.Modal: f_1 = %.2f Hz", s.r->frequencies_hz[0]) + ", " + std::to_string(s.r->modes.size()) +
               " mode(s)>";
      });

  // ---------------- cable statics ----------------
  {
    py::class_<PyCableStatic> cs(mod, "CableStatic", "Cable + frame equilibrium (non-superposable).");
    auto cable_of = [](const PyCableStatic& s, Id id) {
      for (const auto& [cid, bi] : s.cables)
        if (cid == id) return bi;
      throw py::key_error("no cable " + std::to_string(id));
    };
    // The cable's catenary at its converged (deformed) anchors.
    auto catenary = [](const PyCableStatic& s, std::size_t bi) {
      const Model& m = *s.model;
      const Bar& b = m.bars[bi];
      const std::size_t i1 = static_cast<std::size_t>(m.node_index(b.ni)), i2 = static_cast<std::size_t>(m.node_index(b.nj));
      const Vec3 p1 = m.nodes[i1].x + s.r->result.node_disp[i1].head<3>();
      const Vec3 p2 = m.nodes[i2].x + s.r->result.node_disp[i2].head<3>();
      const auto [ea, w] = cable_ea_w(m, b);
      const Vec3 d = p2 - p1;
      const double dx = std::hypot(d.x(), d.y());
      Catenary2D c = engine([&] { return solve_catenary(dx, d.z(), s.l0[bi], ea, w); });
      return std::make_tuple(c, p1, p2);
    };
    cs.def_property_readonly("iterations", [](const PyCableStatic& s) { return s.r->iterations; })
        .def_property_readonly("cable_ids",
                               [](const PyCableStatic& s) {
                                 std::vector<Id> v;
                                 for (const auto& c : s.cables) v.push_back(c.first);
                                 return v;
                               })
        .def(
            "cable",
            [cable_of, catenary](const PyCableStatic& s, Id id) {
              const std::size_t bi = cable_of(s, id);
              const auto [c, p1, p2] = catenary(s, bi);
              return std::make_tuple(c.H, c.V1, c.V2(), std::max(c.T1(), c.T2()), s.l0[bi]);
            },
            py::arg("id"), "(H, V1, V2, T_max, L0) of a cable at the converged state.")
        .def(
            "cable_shape",
            [cable_of, catenary](const PyCableStatic& s, Id id) {
              const std::size_t bi = cable_of(s, id);
              const auto [c, p1, p2] = catenary(s, bi);
              Vec3 eh(p2.x() - p1.x(), p2.y() - p1.y(), 0.0);
              eh = eh.norm() > 0 ? Vec3(eh / eh.norm()) : Vec3(1, 0, 0);
              std::vector<Vec3T> out;
              for (int k = 0; k <= 20; ++k) {
                Catenary2D sub = c;
                sub.L0 = c.L0 * k / 20.0;
                double lx = 0, lz = 0;
                if (k > 0) catenary_projections(sub, lx, lz);
                const Vec3 p = p1 + lx * eh + Vec3(0, 0, lz);
                out.emplace_back(p.x(), p.y(), p.z());
              }
              return out;
            },
            py::arg("id"), "21 points along the converged cable centre-line, end 1 to end 2.");
    bind_case_accessors(
        cs, [](const PyCableStatic& s) { return std::pair<ModelPtr, const CaseResult&>(s.model, s.r->result); },
        "bar_end_forces");
    cs.def("__repr__", [](const PyCableStatic& s) {
      if (s.cables.empty()) return std::string("<keta.CableStatic: no cables>");
      const std::size_t bi = s.cables[0].second;
      return "<keta.CableStatic: " + std::to_string(s.cables.size()) + " cable(s), H_1 = " +
             fmt("%.1f", s.r->cable_H[bi]) + " N, T_max_1 = " + fmt("%.1f", s.r->cable_Tmax[bi]) + " N, " +
             std::to_string(s.r->iterations) + " iteration(s)>";
    });
  }
}

void bind_module_dynamics(py::module_& mod) {
  // ---------------- harmonic ----------------
  py::class_<PyHarmonic>(mod, "Harmonic", "Steady-state frequency response by modal superposition.")
      .def("freqs", [](const PyHarmonic& s) { return s.r->freqs_hz; })
      .def_property_readonly("omegas", [](const PyHarmonic& s) { return s.r->omegas; })
      .def_property_readonly("zetas", [](const PyHarmonic& s) { return s.r->zetas; })
      .def_property_readonly("n_modes", [](const PyHarmonic& s) { return s.r->omegas.size(); })
      .def_property_readonly("node_ids", [](const PyHarmonic& s) { return node_ids_of(*s.model); })
      .def_property_readonly("bar_ids", [](const PyHarmonic& s) { return s.bar_ids; })
      .def(
          "amplitude",
          [](const PyHarmonic& s, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            std::vector<double> v;
            for (std::size_t fi = 0; fi < s.r->freqs_hz.size(); ++fi) v.push_back(s.r->disp_amp[fi][ni](static_cast<int>(comp)));
            return v;
          },
          py::arg("node"), py::arg("comp"))
      .def(
          "amplitude_at",
          [](const PyHarmonic& s, std::size_t fi, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            return s.r->disp_amp[check_fi(fi, s.r->freqs_hz.size())][ni](static_cast<int>(comp));
          },
          py::arg("fi"), py::arg("node"), py::arg("comp"))
      .def(
          "real",
          [](const PyHarmonic& s, std::size_t fi, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            return s.r->disp_re[check_fi(fi, s.r->freqs_hz.size())][ni](static_cast<int>(comp));
          },
          py::arg("fi"), py::arg("node"), py::arg("comp"))
      .def(
          "imag",
          [](const PyHarmonic& s, std::size_t fi, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            return s.r->disp_im[check_fi(fi, s.r->freqs_hz.size())][ni](static_cast<int>(comp));
          },
          py::arg("fi"), py::arg("node"), py::arg("comp"))
      .def(
          "phase",
          [](const PyHarmonic& s, std::size_t fi, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            fi = check_fi(fi, s.r->freqs_hz.size());
            return std::atan2(s.r->disp_im[fi][ni](static_cast<int>(comp)), s.r->disp_re[fi][ni](static_cast<int>(comp)));
          },
          py::arg("fi"), py::arg("node"), py::arg("comp"))
      .def(
          "peak",
          [](const PyHarmonic& s, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            std::pair<double, double> best{0.0, -std::numeric_limits<double>::infinity()};
            for (std::size_t fi = 0; fi < s.r->freqs_hz.size(); ++fi) {
              const double a = s.r->disp_amp[fi][ni](static_cast<int>(comp));
              if (a > best.second) best = {s.r->freqs_hz[fi], a};
            }
            return best;
          },
          py::arg("node"), py::arg("comp"), "(frequency_hz, amplitude) of the largest response (first maximum).")
      .def(
          "stress_amplitude",
          [](const PyHarmonic& s, std::size_t fi, Id bar, std::size_t station, double y, double z) {
            return stress_amp(*s.r, *s.model, s.bar_ids, check_fi(fi, s.r->freqs_hz.size()), bar, station, y, z);
          },
          py::arg("fi"), py::arg("bar"), py::arg("station"), py::arg("y"), py::arg("z"),
          "|sigma| of the complex normal stress N/A + My z/Iy - Mz y/Iz at fibre (y, z).")
      .def("__repr__", [](const PyHarmonic& s) {
        const auto& f = s.r->freqs_hz;
        return "<keta.Harmonic: " + std::to_string(f.size()) + " sweep freq(s) in [" + fmt("%.3f", f.empty() ? 0 : f.front()) +
               ", " + fmt("%.3f", f.empty() ? 0 : f.back()) + "] Hz, " + std::to_string(s.r->omegas.size()) + " mode(s)>";
      });

  // ---------------- PSD ----------------
  py::class_<PyPsd>(mod, "Psd", "Random (PSD) response.")
      .def("freqs", [](const PyPsd& s) { return s.r->freqs_hz; })
      .def("input_psd", [](const PyPsd& s) { return s.input; })
      .def_property_readonly("n_modes", [](const PyPsd& s) { return s.h->omegas.size(); })
      .def_property_readonly("node_ids", [](const PyPsd& s) { return node_ids_of(*s.model); })
      .def_property_readonly("bar_ids", [](const PyPsd& s) { return s.bar_ids; })
      .def(
          "disp_psd",
          [](const PyPsd& s, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            std::vector<double> v;
            for (const auto& f : s.r->disp_psd) v.push_back(f[ni](static_cast<int>(comp)));
            return v;
          },
          py::arg("node"), py::arg("comp"))
      .def(
          "disp_psd_at",
          [](const PyPsd& s, std::size_t fi, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            return s.r->disp_psd[check_fi(fi, s.r->freqs_hz.size())][ni](static_cast<int>(comp));
          },
          py::arg("fi"), py::arg("node"), py::arg("comp"))
      .def(
          "peak_disp_psd",
          [](const PyPsd& s, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            std::pair<double, double> best{0.0, -std::numeric_limits<double>::infinity()};
            for (std::size_t fi = 0; fi < s.r->freqs_hz.size(); ++fi) {
              const double v = s.r->disp_psd[fi][ni](static_cast<int>(comp));
              if (v > best.second) best = {s.r->freqs_hz[fi], v};
            }
            return best;
          },
          py::arg("node"), py::arg("comp"))
      .def(
          "disp_rms",
          [](const PyPsd& s, Id node, std::size_t comp) {
            const std::size_t ni = node_pos(*s.model, node);
            comp = check_comp(comp);
            std::vector<double> v;
            for (const auto& f : s.r->disp_psd) v.push_back(f[ni](static_cast<int>(comp)));
            return std::sqrt(std::max(0.0, trapz(s.r->freqs_hz, v)));
          },
          py::arg("node"), py::arg("comp"), "sqrt of the trapezoidal integral of the displacement PSD.")
      .def(
          "stress_psd",
          [](const PyPsd& s, std::size_t fi, Id bar, std::size_t station, double y, double z) {
            fi = check_fi(fi, s.r->freqs_hz.size());
            const double a = stress_amp(*s.h, *s.model, s.bar_ids, fi, bar, station, y, z);
            return a * a * s.input[fi];
          },
          py::arg("fi"), py::arg("bar"), py::arg("station"), py::arg("y"), py::arg("z"))
      .def(
          "stress_rms",
          [](const PyPsd& s, Id bar, std::size_t station, double y, double z) {
            std::vector<double> v;
            for (std::size_t fi = 0; fi < s.r->freqs_hz.size(); ++fi) {
              const double a = stress_amp(*s.h, *s.model, s.bar_ids, fi, bar, station, y, z);
              v.push_back(a * a * s.input[fi]);
            }
            return std::sqrt(std::max(0.0, trapz(s.r->freqs_hz, v)));
          },
          py::arg("bar"), py::arg("station"), py::arg("y"), py::arg("z"))
      .def("harmonic", [](const PyPsd& s) { return PyHarmonic{s.model, s.h, s.bar_ids}; },
           "The unit-pattern harmonic response behind the PSDs.")
      .def("__repr__", [](const PyPsd& s) {
        const auto& f = s.r->freqs_hz;
        return "<keta.Psd: " + std::to_string(f.size()) + " sweep freq(s) in [" + fmt("%.3f", f.empty() ? 0 : f.front()) +
               ", " + fmt("%.3f", f.empty() ? 0 : f.back()) + "] Hz, " + std::to_string(s.h->omegas.size()) + " mode(s)>";
      });

  // ---------------- spectrum ----------------
  py::class_<PySpectrum>(mod, "Spectrum", "Response-spectrum results in one direction.")
      .def_property_readonly("dir", [](const PySpectrum& s) { return s.dir; })
      .def_property_readonly("n_modes", [](const PySpectrum& s) { return s.r->frequencies_hz.size(); })
      .def_property_readonly("node_ids", [](const PySpectrum& s) { return node_ids_of(*s.model); })
      .def("frequencies", [](const PySpectrum& s) { return s.r->frequencies_hz; })
      .def("omegas",
           [](const PySpectrum& s) {
             std::vector<double> v;
             for (double f : s.r->frequencies_hz) v.push_back(2 * kPi * f);
             return v;
           })
      .def("zetas", [](const PySpectrum& s) { return s.r->zetas; })
      .def("participation", [](const PySpectrum& s) { return s.r->participation; })
      .def("effective_mass",
           [](const PySpectrum& s) {
             std::vector<double> v;
             for (double g : s.r->participation) v.push_back(g * g);
             return v;
           })
      .def("total_mass", [](const PySpectrum& s) { return s.total_mass; })
      .def(
          "combined_disp", [](const PySpectrum& s, Id node) { return dof6(s.r->combined_disp.at(node_pos(*s.model, node))); },
          py::arg("node"))
      .def(
          "modal_disp",
          [](const PySpectrum& s, std::size_t mode, Id node) {
            if (mode >= s.r->modal_disp.size()) throw py::key_error("no mode " + std::to_string(mode));
            return dof6(s.r->modal_disp[mode][node_pos(*s.model, node)]);
          },
          py::arg("mode"), py::arg("node"))
      .def("base_shear", [](const PySpectrum& s) { return s.r->base_shear; })
      .def("modal_base_shear", [](const PySpectrum& s) { return s.r->modal_base_shear; })
      .def("__repr__", [](const PySpectrum& s) {
        return "<keta.Spectrum: dir " + std::to_string(s.dir) + ", " + std::to_string(s.r->frequencies_hz.size()) +
               " mode(s), base shear = " + fmt("%.1f", s.r->base_shear) + " N>";
      });

  mod.def(
      "section_props",
      [](const std::string& designation) {
        std::optional<Section> s = catalog_section(designation);
        if (!s) throw py::key_error("unknown profile '" + designation + "'");
        std::map<std::string, double> m{{"A", s->props.A}, {"Iy", s->props.Iy}, {"Iz", s->props.Iz}, {"J", s->props.J}};
        if (s->props.Avy > 0) m["Asy"] = s->props.Avy;
        if (s->props.Avz > 0) m["Asz"] = s->props.Avz;
        if (s->props.Iw > 0) m["Iw"] = s->props.Iw;
        return m;
      },
      py::arg("designation"), "Catalog section properties: A, Iy, Iz, J, Asy, Asz, Iw (SI).");
}

void bind_module_ec2(py::module_& mod) {
  // ---------------- EN 1992 design drivers ----------------
  mod.def(
      "rc_design_section",
      [](const std::vector<std::vector<std::pair<double, double>>>& rings,
         const std::vector<std::tuple<double, double, double>>& bars, double fck, double fcd, double N, double My,
         double Mz, double fyd, double ftd, double eps_ud, double Es, double omega_max) {
        ec2::RcPolySection s;
        for (const auto& r : rings) {
          std::vector<Vec2> ring;
          for (const auto& [y, z] : r) ring.emplace_back(y, z);
          s.rings.push_back(ring);
        }
        for (const auto& [y, z, w] : bars) s.bars.push_back({y, z, w});
        s.fck = fck;
        s.fcd = fcd;
        s.steel.fyd = fyd;
        s.steel.ftd = ftd;
        s.steel.eps_ud = eps_ud;
        s.steel.Es = Es;
        const ec2::SectionDesign r = engine([&] { return ec2::design_section(s, N, My, Mz, omega_max); });
        py::dict d;
        d["feasible"] = r.feasible;
        d["As"] = r.As;
        d["omega"] = r.omega;
        d["plane"] = std::make_tuple(r.plane.e0, r.plane.ky, r.plane.kz);
        d["resultants"] = std::make_tuple(r.R(0), r.R(1), r.R(2));
        return d;
      },
      py::arg("rings"), py::arg("bars"), py::arg("fck"), py::arg("fcd"), py::arg("N") = 0.0, py::arg("My") = 0.0,
      py::arg("Mz") = 0.0, py::arg("fyd") = 500e6 / 1.15, py::arg("ftd") = 525e6 / 1.15, py::arg("eps_ud") = 0.025,
      py::arg("Es") = 200e9, py::arg("omega_max") = 1.0,
      "EN 1992-1-1 required reinforcement of a polygonal section (rings in (y, z): outer CCW, holes CW) under "
      "N + My + Mz, bars [(y, z, weight)]: As_i = omega weight_i, the smallest omega (keta/ec2_section.h).");
  mod.def(
      "rc_shell_reinforcement",
      [](double nx, double ny, double nxy, double mx, double my, double mxy, double qx, double qy, double d,
         double z, double fck, double fcd, double fyd, double fywd, double cot_theta, double rho_l) {
        ec2::ShellDesignInput in;
        in.d = d;
        in.z = z;
        in.fck = fck;
        in.fcd = fcd;
        in.fyd = fyd;
        in.fywd = fywd;
        in.cot_theta = cot_theta;
        in.rho_l = rho_l;
        const ec2::ShellReinforcement r =
            engine([&] { return ec2::design_shell_point(nx, ny, nxy, mx, my, mxy, qx, qy, in); });
        py::dict o;
        o["asx_top"] = r.asx_top;
        o["asy_top"] = r.asy_top;
        o["asx_bot"] = r.asx_bot;
        o["asy_bot"] = r.asy_bot;
        o["asw"] = r.asw;
        o["v0"] = r.v0;
        o["vrdc"] = r.vrdc;
        return o;
      },
      py::arg("nx") = 0.0, py::arg("ny") = 0.0, py::arg("nxy") = 0.0, py::arg("mx") = 0.0, py::arg("my") = 0.0,
      py::arg("mxy") = 0.0, py::arg("qx") = 0.0, py::arg("qy") = 0.0, py::arg("d") = 0.0, py::arg("z") = 0.0,
      py::arg("fck") = 30e6, py::arg("fcd") = 0.85 * 30e6 / 1.5, py::arg("fyd") = 500e6 / 1.15,
      py::arg("fywd") = 500e6 / 1.15, py::arg("cot_theta") = 1.0, py::arg("rho_l") = 0.0,
      "Shell layer reinforcement at a point (sandwich + Baumann / Annex F, shear by 6.2a / 6.8); m > 0 "
      "stretches the face +z ('top').");
  {
    auto pts = [](const LinkLaw& l) { return l.pts; };
    mod.def("soil_ka_rankine", &soil::ka_rankine, py::arg("phi"));
    mod.def("soil_kp_rankine", &soil::kp_rankine, py::arg("phi"));
    mod.def("soil_k0_jaky", &soil::k0_jaky, py::arg("phi"));
    mod.def("soil_ka_coulomb", &soil::ka_coulomb, py::arg("phi"), py::arg("delta"));
    mod.def("soil_kp_coulomb", &soil::kp_coulomb, py::arg("phi"), py::arg("delta"));
    mod.def(
        "soil_py_soft_clay",
        [pts](double c, double gamma, double X, double D, double eps50, double J) {
          return pts(engine([&] { return soil::py_soft_clay(c, gamma, X, D, eps50, J); }));
        },
        py::arg("c"), py::arg("gamma"), py::arg("X"), py::arg("D"), py::arg("eps50"), py::arg("J") = 0.5,
        "API RP 2A soft-clay p-y law per unit length (Matlock, static): points [(y, p)].");
    mod.def(
        "soil_py_stiff_clay",
        [pts](double c, double gamma, double X, double D, double eps50, int n) {
          return pts(engine([&] { return soil::py_stiff_clay(c, gamma, X, D, eps50, n); }));
        },
        py::arg("c"), py::arg("gamma"), py::arg("X"), py::arg("D"), py::arg("eps50"), py::arg("n") = 16,
        "Welch-Reese stiff-clay p-y law per unit length (above water, static).");
    mod.def(
        "soil_py_sand",
        [pts](double phi, double gamma, double X, double D, double kk, bool cyclic, int n) {
          return pts(engine([&] { return soil::py_sand(phi, gamma, X, D, kk, cyclic, n); }));
        },
        py::arg("phi"), py::arg("gamma"), py::arg("X"), py::arg("D"), py::arg("k"), py::arg("cyclic") = false,
        py::arg("n") = 16, "API RP 2A sand p-y law per unit length (tanh).");
    mod.def(
        "soil_tz_clay", [pts](double t, double D) { return pts(engine([&] { return soil::tz_clay(t, D); })); },
        py::arg("tmax"), py::arg("D"), "API t-z law of clay (peak kept, no softening) per unit length.");
    mod.def(
        "soil_tz_sand", [pts](double t) { return pts(engine([&] { return soil::tz_sand(t); })); }, py::arg("tmax"),
        "API t-z law of sand per unit length.");
    mod.def(
        "soil_qz", [pts](double Qp, double D) { return pts(engine([&] { return soil::qz_tip(Qp, D); })); },
        py::arg("Qp"), py::arg("D"), "API q-z tip law (spring axis pointing up; compression only).");
    mod.def(
        "soil_earth_law",
        [pts](double P0, double Pa, double Pp, double K) {
          return pts(engine([&] { return soil::earth_pressure_law(P0, Pa, Pp, K); }));
        },
        py::arg("P0"), py::arg("Pa"), py::arg("Pp"), py::arg("K"),
        "Trilinear earth-pressure law: F(d) = -clamp(P0 - K d, Pa, Pp), d = movement away from the soil.");
  }
  mod.def(
      "track_bridge_uic774",
      [](const py::list& decks, double rail_A, double rail_I, double rail_dT, double k_unloaded, double r_unloaded,
         double k_loaded, double r_loaded, const py::dict& train, double element, double emb_left, double emb_right,
         double rail_E, double rail_alpha, int steps) {
        TrackBridgeInput in;
        in.rail_A = rail_A;
        in.rail_I = rail_I;
        in.rail_dT = rail_dT;
        in.rail_E = rail_E;
        in.rail_alpha = rail_alpha;
        in.k_unloaded = k_unloaded;
        in.r_unloaded = r_unloaded;
        in.k_loaded = k_loaded;
        in.r_loaded = r_loaded;
        in.element = element;
        in.emb_left = emb_left;
        in.emb_right = emb_right;
        for (const py::handle& h : decks) {
          const py::dict d = py::reinterpret_borrow<py::dict>(h);
          auto get = [&](const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
          TrackDeck td;
          td.length = get("length", 0);
          td.A = get("A", 0);
          td.I = get("I", 0);
          td.E = get("E", 210e9);
          td.alpha = get("alpha", 1e-5);
          td.e_track = get("e_track", 0);
          td.e_bearing = get("e_bearing", 0);
          td.k_left = get("k_left", 0);
          td.k_right = get("k_right", 0);
          td.dT = get("dT", 0);
          in.decks.push_back(td);
        }
        auto tget = [&](const char* k, double dflt) { return train.contains(k) ? train[k].cast<double>() : dflt; };
        TrainPass t;
        t.length = tget("length", 300);
        t.horizontal = tget("horizontal", 20e3);
        t.vertical = tget("vertical", 80e3);
        t.forward = train.contains("forward") ? train["forward"].cast<bool>() : true;
        t.front_first = tget("front_first", 0);
        t.front_last = tget("front_last", 0);
        t.positions = train.contains("positions") ? train["positions"].cast<int>() : 16;
        const TrackBridgeResult r = engine([&] { return track_bridge_uic774(in, t, steps); });
        auto state = [](const TrackState& s) {
          py::dict o;
          o["rail_min"] = s.rail_min;
          o["rail_max"] = s.rail_max;
          o["bearing"] = s.bearing;
          o["rel_deck"] = s.rel_deck;
          o["rel_ground"] = s.rel_ground;
          return o;
        };
        py::dict o;
        o["thermal"] = state(r.thermal);
        o["envelope"] = state(r.envelope);
        py::list ps;
        for (const TrackState& s : r.positions) ps.append(state(s));
        o["positions"] = ps;
        return o;
      },
      py::arg("decks"), py::arg("rail_A"), py::arg("rail_I"), py::arg("rail_dT"), py::arg("k_unloaded"),
      py::arg("r_unloaded"), py::arg("k_loaded"), py::arg("r_loaded"), py::arg("train"), py::arg("element") = 2.0,
      py::arg("emb_left") = 300.0, py::arg("emb_right") = 300.0, py::arg("rail_E") = 210e9,
      py::arg("rail_alpha") = 1.2e-5, py::arg("steps") = 10,
      "UIC 774-3 track-structure interaction: thermal state with the unloaded ballast, then every train position "
      "from it with the loaded ballast under the train. decks = [dict(length, A, I, E, alpha, e_track, e_bearing, "
      "k_left, k_right, dT)], train = dict(length, horizontal, vertical, forward, front_first, front_last, "
      "positions). Returns the rail stresses, bearing forces and relative displacements (thermal, per position, "
      "envelope).");
  mod.def(
      "uniaxial_path", &uniaxial_path_py,
      py::arg("law"), py::arg("params"), py::arg("path"), py::arg("stress_driven") = false,
      "Uniaxial history: law 'steel_mixed' (E, fy, H, iso), 'asym_von_mises' (E, ft0, Ht, fc0, Hc), "
      "'damage_simo_ju' (E, r0, A, B) or 'plastic_damage_oliver' (E, fy, h, r0, A, n) driven along the strains "
      "(or stresses) of `path`, one committed step each; returns [(strain, stress)] (keta/fiber_hyst.h). The "
      "memoryless design laws (keta/fiber.h): 'sargin' (fck, or fcm, Ecm, eps_c1, eps_cu1: EN 1992-1-1 Expr. "
      "3.14), 'steel_hardening' (E, fy, k, eps_uk, eps_ud: Fig. 3.8), 'strand' (Ep, fp01k, fpk, eps_uk, eps_ud, "
      "inclined=1, gamma_s=1.15: §3.3.6), 'concrete_tension_softening' (fcd, Et, ft, Gf, lc, eps_c2, eps_cu2, n: "
      "crack band) and 'bfup_t3' (fck, fcm, Ecm, fctm, K, fctk_el, fctfk, eps_u_lim, gamma_c=1.5, gamma_cf=1.3: "
      "NF P18-710).");
  mod.def(
      "hinge_law_at",
      [](std::vector<double> N, std::vector<std::vector<std::pair<double, double>>> laws, double n) {
        HingeLawFamily f;
        f.N = N;
        for (auto& pts : laws) {
          LinkLaw l;
          l.pts = pts;
          f.laws.push_back(l);
        }
        return engine([&] { return f.at(n).pts; });
      },
      py::arg("N"), py::arg("laws"), py::arg("n"),
      "M-theta hinge law at the axial force n, interpolated point by point between the laws [(theta, M)] tabulated "
      "at the ascending forces N (softening allowed; keta/hinge_law.h).");
  mod.def(
      "annex_a_vehicle",
      [](const std::string& name, const std::vector<std::pair<double, double>>& layout) {
        return engine([&] {
          const auto lines = annex_a_axle_lines(name);
          return layout.empty() ? Wheels{} : axle_lines_to_wheels(lines, layout);
        });
      },
      py::arg("name"), py::arg("layout"),
      "EN 1991-2 Annex A special vehicle '600/150' ... '3600/240/240' as wheels [(dx, dy, P [N])]: the axle lines "
      "(1.50 m, 12 m between double bogies) spread over the transverse layout [(dy, fraction)] the caller gives.");
  mod.def(
      "lm2_wheels", [](double beta_q) { return engine([&] { return lm2_wheels(beta_q); }); }, py::arg("beta_q") = 1.0,
      "EN 1991-2 LM2: beta_Q x 400 kN on two wheels 2.00 m apart.");
  mod.def(
      "combine_extreme",
      [](const std::vector<double>& effects, const py::list& actions, const std::string& state, bool maximise,
         const py::dict& rules) {
        std::vector<CombAction> A;
        Id k = 1;
        for (const py::handle& h : actions) {
          const py::dict d = py::reinterpret_borrow<py::dict>(h);
          CombAction a;
          a.case_id = d.contains("case") ? d["case"].cast<Id>() : k;
          ++k;
          a.name = d.contains("name") ? d["name"].cast<std::string>() : "";
          const std::string kd = d.contains("kind") ? d["kind"].cast<std::string>() : "variable";
          a.kind = kd == "permanent" ? ActionKind::Permanent : kd == "settlement" ? ActionKind::Settlement : ActionKind::Variable;
          if (kd != "permanent" && kd != "settlement" && kd != "variable")
            throw py::value_error("combine_extreme: kind 'permanent', 'settlement' or 'variable'");
          const bool var = a.kind == ActionKind::Variable;
          a.gamma_sup = d.contains("gamma_sup") ? d["gamma_sup"].cast<double>() : (var ? 1.5 : -1.0);
          a.gamma_inf = d.contains("gamma_inf") ? d["gamma_inf"].cast<double>() : (var ? 0.0 : -1.0);
          a.psi0 = d.contains("psi0") ? d["psi0"].cast<double>() : 0.7;
          a.psi1 = d.contains("psi1") ? d["psi1"].cast<double>() : 0.5;
          a.psi2 = d.contains("psi2") ? d["psi2"].cast<double>() : 0.3;
          a.exclusive = d.contains("exclusive") ? d["exclusive"].cast<std::string>() : "";
          a.may_lead = d.contains("may_lead") ? d["may_lead"].cast<bool>() : true;
          A.push_back(a);
        }
        CombinationRules r;
        auto num = [&](const char* key, double& v) {
          if (rules.contains(key)) v = rules[key].cast<double>();
        };
        num("gamma_G_sup", r.gamma_G_sup);
        num("gamma_G_inf", r.gamma_G_inf);
        num("gamma_set_sup", r.gamma_set_sup);
        num("gamma_set_inf", r.gamma_set_inf);
        num("xi", r.xi);
        if (rules.contains("eq610ab")) r.eq610ab = rules["eq610ab"].cast<bool>();
        LimitState ls = LimitState::ULS;
        if (state == "characteristic") ls = LimitState::SLS_characteristic;
        else if (state == "frequent") ls = LimitState::SLS_frequent;
        else if (state == "quasi_permanent") ls = LimitState::SLS_quasi_permanent;
        else if (state != "uls") throw py::value_error("combine_extreme: state 'uls', 'characteristic', 'frequent', 'quasi_permanent'");
        const CombinationExtreme c = engine([&] { return combination_extreme(effects, A, r, ls, maximise); });
        py::dict out;
        out["value"] = c.value;
        out["label"] = c.label;
        std::vector<std::pair<Id, double>> t;
        for (const CombTerm& x : c.terms) t.emplace_back(x.case_id, x.factor);
        out["terms"] = t;
        return out;
      },
      py::arg("effects"), py::arg("actions"), py::arg("state") = "uls", py::arg("maximise") = true, py::arg("rules") = py::dict(),
      "EN 1990 combination as an EXACT extreme of one response: effects per action (unit factor), actions = [dict(case, "
      "name, kind='permanent'|'settlement'|'variable', gamma_sup, gamma_inf, psi0, psi1, psi2, exclusive='tag', "
      "may_lead)], rules = the annex's dict(gamma_G_sup, gamma_G_inf, gamma_set_sup, gamma_set_inf, xi, eq610ab); "
      "returns value, label and the acting terms (keta/combination.h).");
  mod.def(
      "response_spectrum_of",
      [](const std::vector<double>& accel, double dt, const std::vector<double>& periods, double zeta) {
        return engine([&] { return response_spectrum_of(accel, dt, periods, zeta); });
      },
      py::arg("accel"), py::arg("dt"), py::arg("periods"), py::arg("zeta") = 0.05,
      "Pseudo-acceleration response spectrum of an acceleration history (exact oscillator steps).");
  mod.def(
      "wood_armer",
      [](double mx, double my, double mxy, double alpha) {
        const WoodArmerMoments w = engine([&] { return wood_armer(mx, my, mxy, alpha); });
        py::dict o;
        o["mx_top"] = w.mx_top;
        o["my_top"] = w.my_top;
        o["mx_bot"] = w.mx_bot;
        o["my_bot"] = w.my_bot;
        return o;
      },
      py::arg("mx"), py::arg("my"), py::arg("mxy"), py::arg("alpha") = kPi / 2,
      "Wood-Armer design moments (x bars and y' bars at alpha from x, rad): top >= 0, bottom <= 0; m > 0 "
      "stretches the face +z. Combine the moments first, then design.");
  mod.def(
      "wood_armer_assess",
      [](double mx, double my, double mxy, double mx_top, double my_top, double mx_bot, double my_bot, double alpha) {
        const WoodArmerAssessment a = engine(
            [&] { return wood_armer_assess(mx, my, mxy, SlabCapacity{mx_top, my_top, mx_bot, my_bot}, alpha); });
        py::dict o;
        o["u"] = a.u;
        o["u_top"] = a.u_top;
        o["u_bot"] = a.u_bot;
        o["theta_top"] = a.theta_top;
        o["theta_bot"] = a.theta_bot;
        return o;
      },
      py::arg("mx"), py::arg("my"), py::arg("mxy"), py::arg("mx_top"), py::arg("my_top"), py::arg("mx_bot"),
      py::arg("my_bot"), py::arg("alpha") = kPi / 2,
      "Utilisation of a slab point against its four layer capacities (all >= 0): the worst demand / capacity "
      "ratio over every yield line, exact; theta = critical yield-line normal.");
  mod.def(
      "clark_nielsen",
      [](double nx, double ny, double nxy, double mx, double my, double mxy, double z_top, double z_bot,
         double alpha) {
        const ClarkNielsenForces f =
            engine([&] { return clark_nielsen(nx, ny, nxy, mx, my, mxy, z_top, z_bot, alpha); });
        py::dict o;
        o["nx_top"] = f.nx_top;
        o["ny_top"] = f.ny_top;
        o["nx_bot"] = f.nx_bot;
        o["ny_bot"] = f.ny_bot;
        o["fc_top"] = f.fc_top;
        o["fc_bot"] = f.fc_bot;
        return o;
      },
      py::arg("nx") = 0.0, py::arg("ny") = 0.0, py::arg("nxy") = 0.0, py::arg("mx") = 0.0, py::arg("my") = 0.0,
      py::arg("mxy") = 0.0, py::arg("z_top") = 0.0, py::arg("z_bot") = 0.0, py::arg("alpha") = kPi / 2,
      "Clark-Nielsen steel forces (>= 0) and concrete principal forces (<= 0) of the two sandwich layers at "
      "z_top > z_bot from the mid-surface.");
}

void bind_module_ec4(py::module_& mod) {
  // ---------------- EN 1994 ----------------
  auto steel_of = [](const py::dict& d) {
    ec4::SteelI s;
    auto get = [&](const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
    s.h = get("h", 0);
    s.b = get("b", 0);
    s.tw = get("tw", 0);
    s.tf = get("tf", 0);
    s.r = get("r", 0);
    s.A = get("A", 0);
    s.Wpl = get("Wpl", 0);
    s.I = get("I", 0);
    s.fy = get("fy", 355e6);
    s.gM0 = get("gM0", 1.0);
    s.E = get("E", 210e9);
    return s;
  };
  auto slab_of = [](const py::dict& d) {
    ec4::Slab s;
    auto get = [&](const char* k, double dflt) { return d.contains(k) ? d[k].cast<double>() : dflt; };
    s.hs = get("hs", 0);
    s.hp = get("hp", 0);
    s.fck = get("fck", 30e6);
    s.gamma_c = get("gamma_c", 1.5);
    return s;
  };
  mod.def(
      "ec4_sagging",
      [=](const py::dict& steel, const py::dict& slab, double beff, double Nc) {
        const ec4::SaggingResult r = engine([&] { return ec4::sagging(steel_of(steel), slab_of(slab), beff, Nc); });
        py::dict o;
        o["Mpl"] = r.Mpl;
        o["MRd"] = r.MRd;
        o["eta"] = r.eta;
        o["x"] = r.x;
        o["Ncf"] = r.Ncf;
        return o;
      },
      py::arg("steel"), py::arg("slab"), py::arg("beff"), py::arg("Nc") = 1e300,
      "EN 1994-1-1 6.2.1.2 sagging plastic moment (partial connection 6.1 when Nc < Ncf).");
  mod.def(
      "ec4_hogging",
      [=](const py::dict& steel, double As, double fsk, double zs, double gamma_s) {
        ec4::Rebar r;
        r.As = As;
        r.fsk = fsk;
        r.zs = zs;
        r.gamma_s = gamma_s;
        const ec4::HoggingResult h = engine([&] { return ec4::hogging(steel_of(steel), r); });
        py::dict o;
        o["Mpl"] = h.Mpl;
        o["x"] = h.x;
        o["alpha"] = h.alpha;
        o["psi"] = h.psi;
        o["web_class"] = h.web_class;
        o["flange_class"] = h.flange_class;
        o["plastic_valid"] = h.plastic_valid();
        return o;
      },
      py::arg("steel"), py::arg("As"), py::arg("fsk"), py::arg("zs"), py::arg("gamma_s") = 1.15,
      "EN 1994-1-1 6.2.1.2 hogging plastic moment with the classes of the compressed web and flange.");
  mod.def(
      "ec4_deflection",
      [=](const py::dict& steel, const py::dict& slab, double L, double beff, double g_a, double g_c, double g_s,
          double q, bool propped, double n0, double phi_t, double phi_sh, double eps_cs) {
        ec4::DeflectionInput in;
        in.L = L;
        in.beff = beff;
        in.g_a = g_a;
        in.g_c = g_c;
        in.g_s = g_s;
        in.q = q;
        in.propped = propped;
        in.n0 = n0;
        in.phi_t = phi_t;
        in.phi_sh = phi_sh;
        in.eps_cs = eps_cs;
        const ec4::DeflectionResult r =
            engine([&] { return ec4::deflection_simply_supported(steel_of(steel), slab_of(slab), in); });
        py::dict o;
        o["construction"] = r.d_construction;
        o["permanent"] = r.d_perm;
        o["variable"] = r.d_var;
        o["shrinkage"] = r.d_sh;
        o["total"] = r.d_total;
        o["nL"] = r.nL;
        return o;
      },
      py::arg("steel"), py::arg("slab"), py::arg("L"), py::arg("beff"), py::arg("g_a") = 0.0, py::arg("g_c") = 0.0,
      py::arg("g_s") = 0.0, py::arg("q") = 0.0, py::arg("propped") = false, py::arg("n0") = 210.0 / 33.0, py::arg("phi_t") = 0.0,
      py::arg("phi_sh") = 0.0, py::arg("eps_cs") = 0.0,
      "EN 1994-1-1 SLS deflection of a simply supported composite beam with its construction stage.");
  mod.def(
      "ec4_eta_min",
      [](double fy, double Le, double ratio) { return engine([&] { return ec4::eta_min(fy, Le, ratio); }); },
      py::arg("fy"), py::arg("Le"), py::arg("Abf_over_Atf") = 1.0,
      "EN 1994-1-1 6.6.1.2 minimum degree of shear connection.");

  mod.def(
      "set_shear_areas",
      [](const std::string& convention) {
        if (convention == "energy")
          set_default_shear_areas(ShearAreas::Energy);
        else if (convention == "shape_factor")
          set_default_shear_areas(ShearAreas::ShapeFactor);
        else
          throw py::value_error("shear areas: 'energy' or 'shape_factor', got '" + convention + "'");
      },
      py::arg("convention"),
      "Shear areas of the sections given by their shape (rectangle, circle, tube, I, box, T, polygon — the "
      "catalogue profiles of add_section included, with their root fillets) for the beams added FROM NOW ON, "
      "process-wide: 'energy' "
      "(default) = the energy-consistent areas of Saint-Venant flexure for the material's Poisson ratio; "
      "'shape_factor' = the classical shape factors 5/6 (rectangle), 0.9 (circle). A section given by its "
      "properties (add_section_props) keeps the areas it was given.");
  mod.def("shear_areas",
          [] { return std::string(default_shear_areas() == ShearAreas::Energy ? "energy" : "shape_factor"); },
          "The current shear-area convention ('energy' or 'shape_factor').");
  mod.def("set_solver_threads", [](int n) { set_solver_threads(n); }, py::arg("n"),
          "Threads of the supernodal factorization of large systems (1: sequential, the default; 0: all); the "
          "results are bit-identical whatever the count.");
  mod.def("solver_threads", [] { return solver_threads(); }, "The current thread count of the factorization.");
}

PYBIND11_MODULE(_keta, mod) {
  mod.doc() = "keta — bar-structure analysis kernel (C++), Python scripting surface";

  bind_module_model(mod);
  bind_module_model_2(mod);
  bind_module_statics(mod);
  bind_module_dynamics(mod);
  bind_module_ec2(mod);
  bind_module_ec4(mod);
  bind_phasing(mod);
  bind_traffic(mod);
  bind_wind(mod);
  bind_fire(mod);
  bind_st1(mod);
  bind_stability(mod);
  bind_creep(mod);
  bind_rc(mod);
  bind_dynamics(mod);
  bind_modeling(mod);
  bind_cables(mod);
  mod.def(
      "results_info",
      [](const std::string& path) {
        const StoredResults s = engine([&] { return load_results_file(path); });
        py::dict d;
        d["schema_version"] = s.schema_version;
        d["generator"] = s.provenance.generator;
        d["note"] = s.provenance.note;
        d["sign_convention"] = s.sign_convention;
        d["cases"] = s.results.cases.size();
        d["combinations"] = s.results.combos.size();
        return d;
      },
      py::arg("path"),
      "Header of a result store written by Results.save: dict(schema_version, generator, note, sign_convention, "
      "cases, combinations).");
  mod.attr("__version__") = KETA_VERSION;
}

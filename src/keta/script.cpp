// Export to script (see script.h). The blocker list names what the model
// carries beyond what the Python builder can express.
#include "keta/script.h"

#include <charconv>
#include <cmath>
#include <cstdio>
#include <sstream>

namespace keta {

namespace {

// Shortest round-trip literal Python reads as a float ("0.005", "2.1e+11",
// "6.0").
std::string f(double x) {
  char buf[64];
  auto r = std::to_chars(buf, buf + sizeof buf, x);
  std::string s(buf, r.ptr);
  if (s.find_first_of(".eEni") == std::string::npos) s += ".0";
  return s;
}

std::string id(Id i) { return std::to_string(i); }

// A double-quoted Python literal, pure ASCII (non-ASCII UTF-8 as \x / \u / \U).
std::string py_str(const std::string& s) {
  std::string out = "\"";
  for (std::size_t i = 0; i < s.size();) {
    const unsigned char c = static_cast<unsigned char>(s[i]);
    if (c < 0x80) {
      switch (c) {
        case '\\': out += "\\\\"; break;
        case '"': out += "\\\""; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default: out += static_cast<char>(c);
      }
      ++i;
      continue;
    }
    // Decode one UTF-8 sequence (invalid bytes pass as \xHH).
    int len = c >= 0xF0 ? 4 : c >= 0xE0 ? 3 : c >= 0xC0 ? 2 : 1;
    std::uint32_t cp = len == 1 ? c : c & (0xFF >> (len + 1));
    bool ok = len > 1 && i + static_cast<std::size_t>(len) <= s.size();
    for (int k = 1; ok && k < len; ++k) {
      const unsigned char cc = static_cast<unsigned char>(s[i + static_cast<std::size_t>(k)]);
      ok = (cc & 0xC0) == 0x80;
      cp = (cp << 6) | (cc & 0x3F);
    }
    if (!ok) {
      cp = c;
      len = 1;
    }
    char buf[16];
    if (cp <= 0xFF)
      std::snprintf(buf, sizeof buf, "\\x%02x", static_cast<unsigned>(cp));
    else if (cp <= 0xFFFF)
      std::snprintf(buf, sizeof buf, "\\u%04x", static_cast<unsigned>(cp));
    else
      std::snprintf(buf, sizeof buf, "\\U%08x", static_cast<unsigned>(cp));
    out += buf;
    i += static_cast<std::size_t>(len);
  }
  return out + "\"";
}

const char* frame_str(LoadFrame fr) {
  switch (fr) {
    case LoadFrame::Global: return "global";
    case LoadFrame::Local: return "local";
    case LoadFrame::GlobalProjected: return "projected";
  }
  return "global";
}

std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (const std::string& x : v) s += (s.empty() ? "" : ", ") + x;
  return s;
}

std::string feat(std::size_t n, const char* sing, const char* plur, const char* why) {
  return std::to_string(n) + " " + (n == 1 ? sing : plur) + " (" + why + ")";
}

bool is_bar_load(Load::Kind k) {
  return k == Load::Kind::BarUdl || k == Load::Kind::BarLinear || k == Load::Kind::BarParabolic ||
         k == Load::Kind::BarPoint;
}

std::vector<std::string> blockers(const Model& m) {
  std::vector<std::string> out;
  std::size_t n = 0;
  for (const Bar& b : m.bars) n += b.kind == BarKind::Cable;
  if (n)
    out.push_back(feat(n, "elastic-catenary cable", "elastic-catenary cables",
                       "linear analyze() rejects cables; use the separate cable_static solver"));
  if (!m.tendons.empty())
    out.push_back(feat(m.tendons.size(), "prestress tendon", "prestress tendons", "no keta Python surface"));
  n = m.slabs.size() + m.shells.size();
  for (const LoadCase& c : m.cases)
    for (const Load& l : c.loads)
      n += l.kind == Load::Kind::SurfacePressure || l.kind == Load::Kind::SurfaceThermal;
  if (n)
    out.push_back(feat(n, "facet slab / shell / surface load", "facet slabs / shells / surface loads",
                       "no keta Python surface"));
  if (!m.tapered.empty())
    out.push_back(feat(m.tapered.size(), "tapered bar", "tapered bars", "no keta Python surface"));
  n = 0;
  for (const LoadCase& c : m.cases)
    for (const Load& l : c.loads)
      n += l.kind == Load::Kind::Prestrain || l.kind == Load::Kind::PrestrainPoint ||
           l.kind == Load::Kind::PrestrainLinear || l.kind == Load::Kind::PrestrainParabolic;
  if (n) out.push_back(feat(n, "prestrain load", "prestrain loads", "no keta Python surface"));
  n = m.nodal_masses.size() + m.bar_masses.size();
  if (n) out.push_back(feat(n, "nodal / bar mass record", "nodal / bar mass records", "no keta Python surface"));
  n = 0;
  for (const Support& s : m.supports) {
    bool odd = s.rot != 0 || s.decol_axis >= 0 || s.spring6.has_value();
    for (int k : s.oneside) odd = odd || k != 0;
    n += odd;
  }
  if (n)
    out.push_back(feat(n, "rotated / unilateral / coupled support", "rotated / unilateral / coupled supports",
                       "the Python builder has rigid and diagonal-spring supports only"));
  n = 0;
  for (const Bar& b : m.bars)
    for (int k : {0, 1, 2, 3, 6, 7, 8, 9}) {
      if (b.releases[static_cast<std::size_t>(k)]) {
        ++n;
        break;
      }
    }
  if (n)
    out.push_back(feat(n, "bar with releases other than ry / rz", "bars with releases other than ry / rz",
                       "set_releases covers the bending rotations only"));
  n = 0;
  for (const Bar& b : m.bars) n += (b.ecc[0] - b.ecc[1]).squaredNorm() > 0;
  if (n)
    out.push_back(feat(n, "bar with unequal end eccentricities", "bars with unequal end eccentricities",
                       "add_beam(ecc=) offsets both ends alike"));
  n = 0;
  for (const LoadCase& c : m.cases)
    for (const Load& l : c.loads) n += is_bar_load(l.kind) && l.M.squaredNorm() > 0;
  if (n) out.push_back(feat(n, "span moment load", "span moment loads", "no keta Python surface"));

  // Every numeric the emitter would bake.
  std::vector<double> v;
  auto add3 = [&](const Vec3& x) { v.insert(v.end(), {x.x(), x.y(), x.z()}); };
  for (const Node& nd : m.nodes) add3(nd.x);
  for (const Material& x : m.materials) v.insert(v.end(), {x.E, x.nu, x.rho, x.alpha});
  for (const Section& s : m.sections) v.insert(v.end(), {s.props.A, s.props.Iy, s.props.Iz, s.props.J, s.props.Avy, s.props.Avz});
  for (const Bar& b : m.bars) {
    v.push_back(b.roll);
    add3(b.ecc[0]);
    v.insert(v.end(), {b.soil.kfx, b.soil.kfy, b.soil.kfz, b.soil.kmx});
  }
  for (const PointMass& p : m.point_masses) {
    v.push_back(p.mass);
    add3(p.inertia);
  }
  for (const Support& s : m.supports) v.insert(v.end(), s.spring.begin(), s.spring.end());
  for (const LoadCase& c : m.cases)
    for (const Load& l : c.loads) {
      add3(l.F);
      add3(l.M);
      add3(l.F2);
      add3(l.g);
      v.insert(v.end(), {l.a, l.b, l.dt, l.gy, l.gz});
      for (const auto& d : l.disp)
        if (d) v.push_back(*d);
    }
  for (const Combination& c : m.combos)
    for (const auto& t : c.terms) v.push_back(t.second);
  n = 0;
  for (double x : v) n += !std::isfinite(x);
  if (n) out.push_back(std::to_string(n) + " non-finite numeric value(s) (NaN or infinity has no faithful re-solve)");
  return out;
}

// A bar load's relative range [x1, x2] on its bar.
std::pair<double, double> rel_range(const Model& m, const Load& l) {
  if (l.relative) return {l.a, l.b < 0 ? 1.0 : l.b};
  const double L = m.bar_chord_length(m.bar(l.bar));
  return {L > 0 ? l.a / L : 0.0, l.b < 0 ? 1.0 : (L > 0 ? l.b / L : 1.0)};
}

void emit_load(std::ostream& s, const Model& m, const Load& l, Id cs) {
  const std::string tail = "case=" + id(cs);
  switch (l.kind) {
    case Load::Kind::NodalForce: {
      std::vector<std::string> kw;
      const char* names[] = {"fx", "fy", "fz", "mx", "my", "mz"};
      for (int k = 0; k < 6; ++k) {
        const double x = k < 3 ? l.F(k) : l.M(k - 3);
        if (x != 0) kw.push_back(std::string(names[k]) + "=" + f(x));
      }
      kw.push_back(tail);
      s << "m.add_point_load(n" << l.node << ", " << join(kw) << ")\n";
      break;
    }
    case Load::Kind::BarParabolic: {
      const auto [x1, x2] = rel_range(m, l);
      auto t3 = [&](const Vec3& w) { return "(" + f(w.x()) + ", " + f(w.y()) + ", " + f(w.z()) + ")"; };
      s << "m.add_parabolic_load(b" << l.bar << ", " << t3(l.F) << ", " << t3(l.Fmid) << ", " << t3(l.F2)
        << ", x1=" << f(x1) << ", x2=" << f(x2) << ", " << tail << ", frame=" << py_str(frame_str(l.frame))
        << ")\n";
      break;
    }
    case Load::Kind::BarUdl:
    case Load::Kind::BarLinear: {
      const auto [x1, x2] = rel_range(m, l);
      const Vec3 w2 = l.kind == Load::Kind::BarUdl ? l.F : l.F2;
      if (l.kind == Load::Kind::BarUdl && x1 == 0 && x2 == 1) {
        std::vector<std::string> kw;
        const char* names[] = {"wx", "wy", "wz"};
        for (int k = 0; k < 3; ++k)
          if (l.F(k) != 0) kw.push_back(std::string(names[k]) + "=" + f(l.F(k)));
        kw.push_back(tail);
        kw.push_back(std::string("frame=") + py_str(frame_str(l.frame)));
        s << "m.add_udl(b" << l.bar << ", " << join(kw) << ")\n";
      } else {
        s << "m.add_linear_load(b" << l.bar << ", (" << f(l.F.x()) << ", " << f(l.F.y()) << ", " << f(l.F.z())
          << "), (" << f(w2.x()) << ", " << f(w2.y()) << ", " << f(w2.z()) << "), x1=" << f(x1) << ", x2=" << f(x2)
          << ", " << tail << ", frame=" << py_str(frame_str(l.frame)) << ")\n";
      }
      break;
    }
    case Load::Kind::BarPoint: {
      const double at = rel_range(m, l).first;
      std::vector<std::string> kw;
      const char* names[] = {"fx", "fy", "fz"};
      for (int k = 0; k < 3; ++k)
        if (l.F(k) != 0) kw.push_back(std::string(names[k]) + "=" + f(l.F(k)));
      kw.push_back(tail);
      kw.push_back(std::string("frame=") + py_str(frame_str(l.frame)));
      s << "m.add_bar_point_load(b" << l.bar << ", " << f(at) << ", " << join(kw) << ")\n";
      break;
    }
    case Load::Kind::BarThermal: {
      std::vector<std::string> kw;
      if (l.gy != 0) kw.push_back("gy=" + f(l.gy));
      if (l.gz != 0) kw.push_back("gz=" + f(l.gz));
      kw.push_back(tail);
      s << "m.add_thermal(b" << l.bar << ", " << f(l.dt) << ", " << join(kw) << ")\n";
      break;
    }
    case Load::Kind::SelfWeight:
      s << "m.add_self_weight(" << tail << ", g=(" << f(l.g.x()) << ", " << f(l.g.y()) << ", " << f(l.g.z()) << "))\n";
      break;
    case Load::Kind::SupportDisp: {
      std::vector<std::string> kw;
      const char* names[] = {"dx", "dy", "dz", "rx", "ry", "rz"};
      for (std::size_t k = 0; k < 6; ++k)
        if (l.disp[k]) kw.push_back(std::string(names[k]) + "=" + f(*l.disp[k]));
      kw.push_back(tail);
      s << "m.add_support_disp(n" << l.node << ", " << join(kw) << ")\n";
      break;
    }
    default:
      throw std::logic_error("export_script: a blocked load reached emission");
  }
}

void case_asserts(std::ostream& s, const CaseResult& cr, const std::string& accessor) {
  s << "cr = " << accessor << "\n";
  Vec3 sum = Vec3::Zero();
  double scale = 0;
  for (const Reaction& r : cr.reactions) {
    sum += r.R.head<3>();
    for (int k = 0; k < 3; ++k) scale = std::max(scale, std::fabs(r.R(k)));
  }
  const std::string tol = f(std::max(1e-6 * scale, 1e-30));
  s << "sr = cr.sum_reactions()\n";
  for (int k = 0; k < 3; ++k)
    s << "assert abs(sr[" << k << "] - " << f(sum(k)) << ") <= " << tol << ", sr[" << k << "]\n";
  double md = 0;
  for (const Vec6& d : cr.node_disp) md = std::max(md, d.head<3>().norm());
  if (md != 0) {
    s << "md = cr.max_displacement()\n";
    s << "assert abs(md - " << f(md) << ") <= " << f(1e-6 * md) << ", md\n";
  }
  // Spot forces: the (bar, station) of the global max |My| and |N| (first wins).
  struct Best {
    Id bar = 0;
    std::size_t idx = 0;
    double v = 0;
    bool set = false;
  } my, nn;
  for (const BarResult& b : cr.bars)
    for (std::size_t i = 0; i < b.stations.size(); ++i) {
      const Station& st = b.stations[i];
      if (!my.set || std::fabs(st.My) > std::fabs(my.v)) my = {b.bar, i, st.My, true};
      if (!nn.set || std::fabs(st.N) > std::fabs(nn.v)) nn = {b.bar, i, st.N, true};
    }
  auto spot = [&](const Best& b, int col) {
    if (!b.set || std::fabs(b.v) == 0) return;
    const std::string at = "bf[" + std::to_string(b.idx) + "][" + std::to_string(col) + "]";
    s << "bf = cr.bar_forces(b" << b.bar << ")\n";
    s << "assert abs(" << at << " - " << f(b.v) << ") <= " << f(std::max(1e-6 * std::fabs(b.v), 1e-30)) << ", " << at
      << "\n";
  };
  spot(my, 5);
  spot(nn, 1);
}

}  // namespace

std::string export_script(const Model& m, const StaticResults* results, const ScriptOptions& opt) {
  const std::vector<std::string> feats = blockers(m);
  if (!feats.empty()) {
    std::string msg = "cannot export to Python script: " + std::to_string(feats.size()) +
                      " unsupported model feature(s)\n";
    for (const std::string& x : feats) msg += "  - " + x + "\n";
    ScriptError e(msg);
    e.features = feats;
    throw e;
  }
  if (results)
    for (const auto* list : {&results->cases, &results->combos})
      for (const CaseResult& c : *list)
        for (const BarResult& b : c.bars)
          if (static_cast<int>(b.stations.size()) != opt.stations_per_bar)
            throw ScriptError("cannot export to Python script: results station-count mismatch — options request " +
                              std::to_string(opt.stations_per_bar) + " station(s) per bar but the results carry " +
                              std::to_string(b.stations.size()) +
                              "; re-solve and export with the same stations_per_bar");

  std::ostringstream s;
#ifdef KETA_VERSION
  s << "# Generated by keta v" << KETA_VERSION << "\n";
#else
  s << "# Generated by keta\n";
#endif
  s << "#\n"
       "# TWO FORMATS, ONE MODEL:\n"
       "#   * a parametric document is the LIVING DOCUMENT (sliders, rules, node graph).\n"
       "#   * THIS script is a FROZEN MODEL-LEVEL SNAPSHOT: the solved model as currently\n"
       "#     parameterized, with every slider and rule baked into concrete nodes/bars/loads.\n"
       "#     It re-solves headlessly to the same numbers; it does NOT reproduce the graph.\n"
       "#\n"
       "# Units are SI base units throughout: m, N, Pa, kg, rad.\n"
       "#\n"
       "# Sections are emitted as RAW numeric properties via m.add_section_props(...):\n"
       "#   sx (area), iy, iz, ix (torsion J), and optional sy/sz shear areas (omitted:\n"
       "#   rigid in shear). The raw props are numerically complete for this script.\n"
       "\n"
       "import keta as k\n"
       "\n"
       "m = k.Model()\n";

  s << "\n# --- nodes ---\n";
  for (const Node& n : m.nodes)
    s << "n" << n.id << " = m.add_node(" << f(n.x.x()) << ", " << f(n.x.y()) << ", " << f(n.x.z()) << ")\n";

  s << "\n# --- materials ---\n";
  for (const Material& x : m.materials) {
    std::vector<std::string> kw;
    if (x.rho != 0) kw.push_back("rho=" + f(x.rho));
    if (x.alpha != 0) kw.push_back("alpha_temp=" + f(x.alpha));
    if (!x.name.empty()) kw.push_back("name=" + py_str(x.name));
    if (x.G_given > 0) kw.push_back("g=" + f(x.G_given));
    s << "mat" << x.id << " = m.add_material(" << f(x.E) << ", " << f(x.nu) << (kw.empty() ? std::string() : ", " + join(kw))
      << ")\n";
  }

  s << "\n# --- sections (raw props: sx, iy, iz, ix [+ sy, sz shear areas]) ---\n";
  for (const Section& x : m.sections) {
    const SectionProps& p = x.props;
    std::vector<std::string> kw;
    if (p.Avy > 0) kw.push_back("sy=" + f(p.Avy));
    if (p.Avz > 0) kw.push_back("sz=" + f(p.Avz));
    if (!x.name.empty()) kw.push_back("name=" + py_str(x.name));
    s << "s" << x.id << " = m.add_section_props(" << f(p.A) << ", " << f(p.Iy) << ", " << f(p.Iz) << ", " << f(p.J)
      << (kw.empty() ? std::string() : ", " + join(kw)) << ")\n";
  }

  s << "\n# --- bars (beams + trusses) ---\n";
  for (const Bar& b : m.bars) {
    if (b.kind == BarKind::Truss) {
      s << "b" << b.id << " = m.add_truss(n" << b.ni << ", n" << b.nj << ", section=s" << b.section << ", material=mat"
        << b.material << ")\n";
      continue;
    }
    s << "b" << b.id << " = m.add_beam(n" << b.ni << ", n" << b.nj << ", section=s" << b.section << ", material=mat"
      << b.material;
    if (b.roll != 0) s << ", beta=" << f(b.roll);
    if (b.has_ecc()) s << ", ecc=(" << f(b.ecc[0].x()) << ", " << f(b.ecc[0].y()) << ", " << f(b.ecc[0].z()) << ")";
    s << ")\n";
  }
  bool any_rel = false;
  for (const Bar& b : m.bars) {
    std::vector<std::string> kw;
    if (b.releases[4]) kw.push_back("ry_or=True");
    if (b.releases[5]) kw.push_back("rz_or=True");
    if (b.releases[10]) kw.push_back("ry_ex=True");
    if (b.releases[11]) kw.push_back("rz_ex=True");
    if (kw.empty()) continue;
    if (!any_rel) s << "\n# --- end releases (or = origin, ex = extremity) ---\n";
    any_rel = true;
    s << "m.set_releases(b" << b.id << ", " << join(kw) << ")\n";
  }
  bool any_soil = false;
  for (const Bar& b : m.bars) {
    if (!b.soil.any()) continue;
    if (!any_soil) s << "\n# --- Winkler foundation beds (bar-local bed moduli) ---\n";
    any_soil = true;
    std::vector<std::string> kw;
    const std::pair<const char*, double> mod[] = {
        {"kfx", b.soil.kfx}, {"kfy", b.soil.kfy}, {"kfz", b.soil.kfz}, {"kmx", b.soil.kmx}};
    for (const auto& [name, v] : mod)
      if (v != 0) kw.push_back(std::string(name) + "=" + f(v));
    s << "m.set_soil(b" << b.id << ", " << join(kw) << ")\n";
  }

  if (!m.point_masses.empty()) {
    s << "\n# --- point masses ---\n";
    for (const PointMass& p : m.point_masses) {
      s << "m.add_point_mass(n" << p.node << ", " << f(p.mass);
      if (p.inertia.squaredNorm() > 0)
        s << ", inertia=[" << f(p.inertia.x()) << ", " << f(p.inertia.y()) << ", " << f(p.inertia.z()) << "]";
      s << ")\n";
    }
  }

  bool any_sup = false;
  for (const Support& sp : m.supports) {
    if (!sp.any()) continue;  // blocks nothing: not a support (keta's validation rejects it)
    if (!any_sup) s << "\n# --- supports ---\n";
    any_sup = true;
    bool springs = false;
    for (int k = 0; k < 6; ++k) springs = springs || (sp.spring[static_cast<std::size_t>(k)] != 0 && !sp.fixed[static_cast<std::size_t>(k)]);
    const std::array<bool, 6> all{true, true, true, true, true, true}, pin{true, true, true, false, false, false};
    if (!springs && sp.fixed == all) {
      s << "m.fix(n" << sp.node << ")\n";
    } else if (!springs && sp.fixed == pin) {
      s << "m.pin(n" << sp.node << ")\n";
    } else {
      std::vector<std::string> kw;
      const char* names[] = {"dx", "dy", "dz", "rx", "ry", "rz"};
      const char* knames[] = {"kx", "ky", "kz", "krx", "kry", "krz"};
      for (std::size_t k = 0; k < 6; ++k)
        if (sp.fixed[k]) kw.push_back(std::string(names[k]) + "=True");
      for (std::size_t k = 0; k < 6; ++k)
        if (sp.spring[k] != 0 && !sp.fixed[k]) kw.push_back(std::string(knames[k]) + "=" + f(sp.spring[k]));
      s << "m.add_support(n" << sp.node << (kw.empty() ? std::string() : ", " + join(kw)) << ")\n";
    }
  }

  s << "\n# --- loads (grouped by load case) ---\n";
  for (const LoadCase& c : m.cases) {
    s << "\n# load case " << c.id << "\n";
    if (!c.name.empty()) s << "m.name_case(" << c.id << ", " << py_str(c.name) << ")\n";
    for (const Load& l : c.loads) emit_load(s, m, l, c.id);
    // An empty unnamed case is otherwise never created on the Python side.
    if (c.name.empty() && c.loads.empty()) s << "m.name_case(" << c.id << ", \"LC" << c.id << "\")\n";
  }

  if (!m.combos.empty()) {
    s << "\n# --- combinations ---\n";
    for (const Combination& c : m.combos) {
      std::vector<std::string> terms;
      for (const auto& [cs, fac] : c.terms) terms.push_back("(" + f(fac) + ", " + id(cs) + ")");
      s << "c" << c.id << " = m.add_combination([" << join(terms) << "]"
        << (c.name.empty() ? std::string() : ", name=" + py_str(c.name)) << ")\n";
    }
  }

  s << "\n# --- validate & solve ---\n"
       "problems = m.validate()\n"
       "assert not problems, problems\n"
       "r = m.analyze(stations="
    << opt.stations_per_bar << ")\n";

  if (results) {
    s << "\n# --- assertions: the frozen numbers this model solved to ---\n"
         "# Each tolerance is 1e-6 * S, S the meaningful scale of that quantity in\n"
         "# this case (max |reaction| for the sums, the value itself for a peak),\n"
         "# never a bare unit constant.\n";
    for (const CaseResult& c : results->cases) {
      s << "\n# case " << c.id << "\n";
      case_asserts(s, c, "r.case(" + id(c.id) + ")");
    }
    for (const CaseResult& c : results->combos) {
      s << "\n# combination " << c.id << "\n";
      case_asserts(s, c, "r.combo(c" + id(c.id) + ")");
    }
    s << "\nprint(\"OK - model re-solved with matching numbers\")\n";
  }
  return s.str();
}

}  // namespace keta

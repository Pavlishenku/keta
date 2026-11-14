// keta_bench: timings of the main drivers on large deterministic models
// (performance review).
//
// Usage:
//   keta_bench [scale [threads]]          frame + deck (static, 20 modes)
//   keta_bench <scenario> [sizes...]      one scenario, a size sweep
//
//   frame    : 3-D building frame (nx x ny bays, nz storeys), 8 load cases,
//              static (stations 11) and 20 modes
//   deck     : facet-shell deck n x m on longitudinal girders (bars), 8 cases,
//              static and 20 modes
//   phasing  : span-by-span continuous concrete girder, size = spans (8 bars
//              each): one cast + tension per span then a long wait; creep along
//              the bar, exact bonded losses, local losses (the ST1 options)
//   phasing_nb : the same girder, 8 spans, size = bars per span
//   cable    : cable-stayed deck, size = stays per side (semi-fan, catenary
//              stays at an end tension), coupled cable + frame Newton
//   il       : influence line of a midspan moment along the central girder of
//              a 3-span grillage deck, size = bars per girder (5 girders)
//   is       : influence surface of an interior support reaction on the same
//              grillage (every node a load position), size = bars per girder
// Quadratic-loop hunt (size = bars per girder of the 5-girder grillage unless
// stated):
//   diag     : model diagnostics (coincident nodes, ...), no mechanism probe
//   uplift   : every girder node on a one-sided vertical spring (active set)
//   shellval : Model::validate of n x 8 facet shells with one ShellImposed
//              load per shell (the creep fields of a shell phasing)
//   static   : linear statics, 8 cases
//   envelope : envelope of 32 static cases
//   convoy   : convoy extreme on the influence line, one reference position
//              per station, 2 axles + 1 segment
//   st1      : ST1 import of a continuous beam of n bars (3 cases)
//   mesh     : constrained Delaunay + refinement of a 10 x 10 m plate with two
//              holes and a constraint segment, size = edges per side (h = 10/n)
//
// Each sweep prints the time of the SOLVE per size (the model construction is
// not timed) and the empirical exponent p = d log t / d log n between
// consecutive sizes (t ~ n^p).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

#include "keta/cable.h"
#include "keta/eigsolve.h"
#include "keta/influence.h"
#include "keta/influence_surface.h"
#include "keta/phasing.h"
#include "keta/statics.h"
#include "keta/supernodal.h"
#include "keta/diagnostics.h"
#include "keta/envelope.h"
#include "keta/st1.h"
#include "keta/mesh2d.h"

using namespace keta;

namespace {
double now() {
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
Model frame(int nx, int ny, int nz) {
  Model m;
  Material s;
  s.E = 210e9;
  s.nu = 0.3;
  s.rho = 7850;
  const Id mi = m.add_material(s);
  Section c;
  c.props.A = 0.02;
  c.props.Iy = c.props.Iz = 3e-4;
  c.props.J = 1e-5;
  const Id si = m.add_section(c);
  auto id = [&](int i, int j, int k) { return static_cast<Id>(1 + i + (nx + 1) * (j + (ny + 1) * k)); };
  for (int k = 0; k <= nz; ++k)
    for (int j = 0; j <= ny; ++j)
      for (int i = 0; i <= nx; ++i) m.add_node(6.0 * i, 6.0 * j, 3.5 * k);
  for (int k = 0; k <= nz; ++k)
    for (int j = 0; j <= ny; ++j)
      for (int i = 0; i <= nx; ++i) {
        if (k == 0) {
          m.add_support(id(i, j, k), {true, true, true, true, true, true});
          continue;
        }
        m.add_bar(id(i, j, k - 1), id(i, j, k), mi, si);
        if (i > 0) m.add_bar(id(i - 1, j, k), id(i, j, k), mi, si);
        if (j > 0) m.add_bar(id(i, j - 1, k), id(i, j, k), mi, si);
      }
  for (int lc = 0; lc < 8; ++lc) {
    std::vector<Load> ls;
    for (int j = 0; j <= ny; ++j)
      for (int k = 1; k <= nz; ++k) ls.push_back(Load::nodal(id(lc % (nx + 1), j, k), Vec3(1e4 * (lc + 1), 2e3, -5e3)));
    m.add_case("c" + std::to_string(lc), ls);
  }
  return m;
}
Model deck(int n, int mm) {
  Model m;
  Material s;
  s.E = 35e9;
  s.nu = 0.2;
  s.rho = 2500;
  const Id mi = m.add_material(s);
  Section g;
  g.props.A = 0.8;
  g.props.Iy = g.props.Iz = 0.4;
  g.props.J = 0.05;
  const Id si = m.add_section(g);
  const double L = 60, B = 12;
  auto id = [&](int i, int j) { return static_cast<Id>(1 + i + (n + 1) * j); };
  for (int j = 0; j <= mm; ++j)
    for (int i = 0; i <= n; ++i) m.add_node(L * i / n, B * j / mm, 0);
  for (int j = 0; j < mm; ++j)
    for (int i = 0; i < n; ++i) {
      m.add_shell({id(i, j), id(i + 1, j), id(i + 1, j + 1)}, 0.25, mi, 1);
      m.add_shell({id(i, j), id(i + 1, j + 1), id(i, j + 1)}, 0.25, mi, 1);
    }
  for (int j = 0; j <= mm; j += std::max(1, mm / 4))
    for (int i = 0; i < n; ++i) m.add_bar(id(i, j), id(i + 1, j), mi, si);
  for (int j = 0; j <= mm; ++j) {
    m.add_support(id(0, j), {true, true, true, false, false, false});
    m.add_support(id(n, j), {false, true, true, false, false, false});
  }
  for (int lc = 0; lc < 8; ++lc) {
    std::vector<Load> ls{Load::surface_pressure(1, Vec3(0, 0, -1e3 * (lc + 1)))};
    ls.push_back(Load::nodal(id(n / 2, mm / 2), Vec3(0, 0, -1e5)));
    m.add_case("c" + std::to_string(lc), ls);
  }
  return m;
}
void run(const char* name, const Model& m) {
  double t0 = now();
  StaticOptions so;
  const StaticResults r = solve_static(m, so);
  double t1 = now();
  ModalOptions mo;
  mo.n_modes = 20;
  const ModalResults md = solve_modal(m, mo);
  double t2 = now();
  std::printf("%-6s nodes %7zu bars %7zu shells %7zu | static (8 cases) %8.3f s | modal (20) %8.3f s | f1 %.4f Hz, u %.3e\n",
              name, m.nodes.size(), m.bars.size(), m.shells.size(), t1 - t0, t2 - t1, md.frequencies_hz[0],
              r.cases[0].node_disp.back()(0));
}

// ---------------- staged construction ----------------

// Continuous prestressed concrete girder of `spans` spans of 30 m, `nb` bars
// per span, built span by span: phase k casts span k (age 7 days), loads it
// with its self weight and tensions its tendon (parabolic, -0.4 m at midspan),
// 14 days apart; a last phase waits to 10 000 days.
struct Staged {
  Model m;
  PhasePlan plan;
};
Staged staged_girder(int spans, int nb) {
  Staged s;
  Model& m = s.m;
  Material conc;
  conc.E = 36e9;
  conc.rho = 0;
  ConcreteTime ct;
  ct.fck = 40e6;
  ct.rh = 70;
  ct.h0 = 0.5;
  conc.concrete = ct;
  m.add_material(conc);
  Section sec;
  sec.shape = ShapeRect{1.0, 2.0};
  m.add_section(sec);
  const double Ls = 30.0;
  std::vector<Id> nodes;
  for (int i = 0; i <= spans * nb; ++i) nodes.push_back(m.add_node(Ls * i / nb, 0, 0));
  for (int i = 0; i < spans * nb; ++i)
    m.add_bar(nodes[static_cast<std::size_t>(i)], nodes[static_cast<std::size_t>(i) + 1], 1, 1);
  for (int k = 0; k <= spans; ++k) {
    const Id n = nodes[static_cast<std::size_t>(k * nb)];
    if (k == 0)
      m.add_support(n, {true, true, true, true, false, false});
    else
      m.add_support(n, {false, true, true, false, false, false});
  }
  for (int k = 0; k < spans; ++k) {
    Tendon t;
    t.id = static_cast<Id>(k + 1);
    for (int i = 0; i < nb; ++i) t.bars.push_back(static_cast<Id>(1 + k * nb + i));
    PolySegment seg;  // e(s) = -0.4 * 4 s (L - s) / L^2
    seg.s0 = 0;
    seg.s1 = Ls;
    seg.cz = {0.0, -1.6 / Ls, 1.6 / (Ls * Ls), 0.0};
    t.profile = {seg};
    t.jack_force = 8e6;
    t.losses.mu = 0.19;
    t.losses.k = 0.01;
    t.losses.slip = 0.006;
    t.losses.Ep = 195e9;
    t.losses.Ap = 5.6e-3;
    t.losses.fpk = 1860e6;
    m.tendons.push_back(t);
    std::vector<Load> g;
    for (Id b : t.bars) g.push_back(Load::udl(b, Vec3(0, 0, -50e3)));
    m.add_case("g" + std::to_string(k), g);
  }
  m.validate();
  for (int k = 0; k < spans; ++k) {
    Phase p;
    p.name = "span " + std::to_string(k);
    p.time = 28.0 + 14.0 * k;
    for (int i = 0; i < nb; ++i) {
      const Id b = static_cast<Id>(1 + k * nb + i);
      p.activate_bars.push_back(b);
      p.activate_age[b] = 7.0;
    }
    p.cases = {static_cast<Id>(k + 1)};
    p.tension = {{static_cast<Id>(k + 1), 1.0}};
    s.plan.phases.push_back(p);
  }
  Phase w;
  w.name = "long term";
  w.time = 10000;
  s.plan.phases.push_back(w);
  return s;
}
PhasingOptions staged_options() {
  PhasingOptions o;
  o.creep_along_bar = true;
  o.local_losses = true;
  o.exact_losses = true;
  o.born_deformed = true;
  o.creep_engine = CreepEngine::Exact;
  return o;
}

// ---------------- cable-stayed deck ----------------

// Symmetric semi-fan cable-stayed deck: 2 n stays at 10 m along the deck
// (2 deck bars between anchors), anchored on the pylon top 2 m apart, each
// stay at the end tension carrying its 10 m of deck load; deck pinned at both
// ends, pylon clamped 20 m below the deck.
Model cable_stayed(int n) {
  Model m;
  Material conc;
  conc.E = 36e9;
  conc.rho = 0;
  const Id mc = m.add_material(conc);
  Material steel;
  steel.E = 195e9;
  steel.rho = 7850;
  const Id ms = m.add_material(steel);
  Section box;
  box.props.A = 6.0;
  box.props.Iy = 4.0;
  box.props.Iz = 60.0;
  box.props.J = 10.0;
  const Id sb = m.add_section(box);
  Section pyl;
  pyl.props.A = 20.0;
  pyl.props.Iy = pyl.props.Iz = 60.0;
  pyl.props.J = 80.0;
  const Id sp = m.add_section(pyl);
  Section stay;
  stay.props.A = 6e-3;
  stay.props.Iy = stay.props.Iz = 1e-8;
  stay.props.J = 1e-8;
  const Id sc = m.add_section(stay);
  const double a = 10.0, q = 150e3;
  const int nd = 2 * (2 * n + 1);  // deck bars: anchors at +/- a, ..., +/- n a, 2 bars each
  const double L = a * (2 * n + 1);
  std::vector<Id> deck;
  for (int i = 0; i <= nd; ++i) deck.push_back(m.add_node(-L / 2 + L * i / nd, 0, 0));
  for (int i = 0; i < nd; ++i) m.add_bar(deck[static_cast<std::size_t>(i)], deck[static_cast<std::size_t>(i) + 1], mc, sb);
  m.add_support(deck.front(), {true, true, true, true, false, false});
  m.add_support(deck.back(), {false, true, true, true, false, false});
  // pylon: base -20 m, deck crossing, anchor zone from 0.5 L / 2 up
  const Id base = m.add_node(0, 0.5, -20);
  m.add_support(base, {true, true, true, true, true, true});
  const Id cross = m.add_node(0, 0.5, 0);
  m.add_bar(base, cross, mc, sp);
  const double z0 = 0.25 * L;
  Id below = cross;
  std::vector<Id> anchor;
  for (int k = 0; k < n; ++k) {
    const Id t = m.add_node(0, 0.5, z0 + 2.0 * k);
    m.add_bar(below, t, mc, sp);
    anchor.push_back(t);
    below = t;
  }
  // deck node at x = i a (i = +/-1..n): index (L/2 + i a) / (L / nd)
  for (int side : {-1, 1})
    for (int i = 1; i <= n; ++i) {
      const std::size_t di = static_cast<std::size_t>(std::lround((L / 2 + side * i * a) / (L / nd)));
      const Id top = anchor[static_cast<std::size_t>(n - i)];  // the outer stays anchor highest
      const Vec3 d = m.nodes[static_cast<std::size_t>(m.node_index(top))].x - m.nodes[static_cast<std::size_t>(m.node_index(deck[di]))].x;
      const double sin_t = d(2) / d.norm();
      const Id c = m.add_bar(deck[di], top, ms, sc, BarKind::Cable);
      CableProps cp;
      cp.end_tension = q * a / sin_t;
      m.bars[static_cast<std::size_t>(m.bar_index(c))].cable = cp;
    }
  std::vector<Load> ls;
  for (int i = 0; i < nd; ++i) ls.push_back(Load::udl(static_cast<Id>(i + 1), Vec3(0, 0, -q)));
  m.add_case("q", ls);
  m.validate();
  return m;
}

// ---------------- grillage deck (influence) ----------------

// 3 continuous spans of 30 m, `ng` girders 3 m apart, `nx` bars per girder
// (nx a multiple of 3), cross beams at every node; vertical supports at the
// span ends of every girder. Node (i, j) = 1 + i + (nx + 1) j.
Model grillage(int nx, int ng) {
  Model m;
  Material c;
  c.E = 36e9;
  c.nu = 0.2;
  const Id mi = m.add_material(c);
  Section g;
  g.props.A = 1.2;
  g.props.Iy = 0.5;
  g.props.Iz = 0.6;
  g.props.J = 0.1;
  const Id sg = m.add_section(g);
  Section t;
  t.props.A = 0.3;
  t.props.Iy = 0.02;
  t.props.Iz = 0.05;
  t.props.J = 0.01;
  const Id st = m.add_section(t);
  const double L = 90.0, b = 3.0;
  auto id = [&](int i, int j) { return static_cast<Id>(1 + i + (nx + 1) * j); };
  for (int j = 0; j < ng; ++j)
    for (int i = 0; i <= nx; ++i) m.add_node(L * i / nx, b * j, 0);
  for (int j = 0; j < ng; ++j)
    for (int i = 0; i < nx; ++i) m.add_bar(id(i, j), id(i + 1, j), mi, sg);
  for (int j = 1; j < ng; ++j)
    for (int i = 0; i <= nx; ++i) m.add_bar(id(i, j - 1), id(i, j), mi, st);
  for (int j = 0; j < ng; ++j)
    for (int k = 0; k <= 3; ++k) {
      const int i = k * nx / 3;
      if (i == 0)
        m.add_support(id(i, j), {true, j == 0, true, false, false, false});
      else
        m.add_support(id(i, j), {false, false, true, false, false, false});
    }
  m.add_case("none", {});
  m.validate();
  return m;
}

// one-sided vertical springs on every girder node (plus the rigid lines)
Model grillage_uplift(int nx, int ng) {
  Model m = grillage(nx, ng);
  auto id = [&](int i, int j) { return static_cast<Id>(1 + i + (nx + 1) * j); };
  for (int j = 0; j < ng; ++j)
    for (int i = 1; i < nx; ++i) {
      if (i % (nx / 3) == 0) continue;
      Support sp;
      sp.node = id(i, j);
      sp.spring[2] = 2e7;
      sp.oneside[2] = 1;
      m.supports.push_back(sp);
    }
  std::vector<Load> ls;
  for (int j = 0; j < ng; ++j) {
    ls.push_back(Load::nodal(id(nx / 6, j), Vec3(0, 0, -2e6)));
    ls.push_back(Load::nodal(id(5 * nx / 6, j), Vec3(0, 0, 1e6)));
  }
  m.cases.clear();
  m.add_case("lift", ls);
  m.validate();
  return m;
}
// n x 8 facet shells and one ShellImposed load per shell
Model shells_imposed(int n) {
  Model m;
  Material c;
  c.E = 35e9;
  c.nu = 0.2;
  m.add_material(c);
  const int mm = 4;
  auto id = [&](int i, int j) { return static_cast<Id>(1 + i + (n + 1) * j); };
  for (int j = 0; j <= mm; ++j)
    for (int i = 0; i <= n; ++i) m.add_node(1.0 * i, 1.0 * j, 0);
  for (int j = 0; j < mm; ++j)
    for (int i = 0; i < n; ++i) {
      m.add_shell({id(i, j), id(i + 1, j), id(i + 1, j + 1)}, 0.25, 1, 1);
      m.add_shell({id(i, j), id(i + 1, j + 1), id(i, j + 1)}, 0.25, 1, 1);
    }
  for (int j = 0; j <= mm; ++j) m.add_support(id(0, j), {true, true, true, true, true, true});
  std::vector<Load> ls;
  Eigen::Matrix<double, 18, 1> d = Eigen::Matrix<double, 18, 1>::Constant(1e-5);
  for (const ShellTri& sh : m.shells) ls.push_back(Load::shell_imposed(sh.id, d));
  m.add_case("creep", ls);
  return m;
}
std::string st1_beam(int n) {
  std::string t = "OPTION PLANE\nTITRE 'poutre'\nNOEUD\n";
  for (int i = 0; i <= n; ++i) t += std::to_string(i + 1) + " " + std::to_string(0.5 * i) + " 0.0\n";
  for (int i = 0; i <= n; i += std::max(1, n / 8)) t += "APPUI " + std::to_string(i + 1) + (i == 0 ? " DX DY\n" : " DY\n");
  t += "BARRE\n";
  for (int i = 0; i < n; ++i) t += std::to_string(i + 1) + " " + std::to_string(i + 1) + " " + std::to_string(i + 2) + "\n";
  t += "CARA 1 a " + std::to_string(n) + " SX 0.5 IZ 0.05\nCONS TOUT E 3e6 NU 0.2 RO 2.5\n";
  t += "CHARG 1 'pp'\n   POIDS PROPRE TOUT\nFIN\n";
  t += "CHARG 2 'q'\n   BARRE\n   1 a " + std::to_string(n) + " UNI FY -2.0\nFIN\n";
  t += "CHARG 3 'P'\n";
  for (int i = 1; i < n; i += 3) t += "   NOEUD " + std::to_string(i + 1) + " FY -1.0\n";
  t += "FIN\nEXEC CHARG TOUT\n";
  return t;
}

// ---------------- sweeps ----------------

// job(n, t): builds the model, solves it, sets t to the solve time [s] and
// returns a summary line
using Job = std::function<std::string(int, double&)>;
void sweep(const char* name, const std::vector<int>& sizes, const Job& job) {
  double t_prev = 0;
  int n_prev = 0;
  for (int n : sizes) {
    double t = 0;
    std::string info;
    try {
      info = job(n, t);
    } catch (const std::exception& e) {
      std::printf("%-10s n %6d | failed: %s\n", name, n, e.what());
      std::fflush(stdout);
      t_prev = 0;
      continue;
    }
    std::printf("%-10s n %6d | %9.3f s", name, n, t);
    if (n_prev > 0 && t_prev > 0 && t > 0)
      std::printf(" | p %5.2f", std::log(t / t_prev) / std::log(static_cast<double>(n) / n_prev));
    else
      std::printf(" |        ");
    std::printf(" | %s\n", info.c_str());
    std::fflush(stdout);
    t_prev = t;
    n_prev = n;
  }
}
std::string fmt(const char* f, double a, double b = 0, double c = 0, double d = 0) {
  char buf[256];
  std::snprintf(buf, sizeof buf, f, a, b, c, d);
  return buf;
}
std::vector<int> sizes_or(int argc, char** argv, std::vector<int> def) {
  if (argc <= 2) return def;
  std::vector<int> s;
  for (int i = 2; i < argc; ++i) s.push_back(std::atoi(argv[i]));
  return s;
}
}  // namespace

int main(int argc, char** argv) {
  if (argc > 1 && !(std::isdigit(static_cast<unsigned char>(argv[1][0])) || argv[1][0] == '.')) {
    const std::string sc = argv[1];
    std::printf("solver threads: %d\n", solver_threads());
    if (sc == "phasing" || sc == "phasing_nb") {
      const bool by_spans = sc == "phasing";
      sweep(sc.c_str(), sizes_or(argc, argv, by_spans ? std::vector<int>{2, 4, 8, 16} : std::vector<int>{4, 8, 16, 32}),
            [&](int n, double& t) {
              const Staged s = by_spans ? staged_girder(n, 8) : staged_girder(8, n);
              const double t0 = now();
              const PhasedResults r = solve_phased(s.m, s.plan, staged_options());
              t = now() - t0;
              const CaseResult& last = r.phases.back().cumulated;
              double wmin = 0;
              for (const auto& u : last.node_disp) wmin = std::min(wmin, u(2));
              return fmt("bars %.0f phases %.0f | min uz %.4e m", static_cast<double>(s.m.bars.size()),
                         static_cast<double>(r.phases.size()), wmin);
            });
    } else if (sc == "cable") {
      sweep("cable", sizes_or(argc, argv, {4, 8, 16, 32}), [](int n, double& t) {
        const Model m = cable_stayed(n);
        const double t0 = now();
        const CableStaticResults r = solve_cable_static(m, 1);
        t = now() - t0;
        double hmax = 0;
        for (double h : r.cable_H) hmax = std::max(hmax, h);
        return fmt("bars %.0f stays %.0f | Newton %.0f | H max %.4e N", static_cast<double>(m.bars.size()), 2.0 * n,
                   r.iterations, hmax);
      });
    } else if (sc == "il" || sc == "is") {
      const bool line = sc == "il";
      const int ng = 5;
      sweep(sc.c_str(), sizes_or(argc, argv, {60, 120, 240, 480}), [&](int nx, double& t) {
        const Model m = grillage(nx, ng);
        const int jc = ng / 2;
        auto id = [&](int i, int j) { return static_cast<Id>(1 + i + (nx + 1) * j); };
        if (line) {
          std::vector<Id> path_bars;
          for (int i = 0; i < nx; ++i) path_bars.push_back(static_cast<Id>(1 + jc * nx + i));
          const double t0 = now();
          const LoadPath path = build_path(m, path_bars, 11);
          // midspan of the central span, girder jc: bar starting at node 3 nx / 6 ... = nx/2
          const Id bar = static_cast<Id>(1 + jc * nx + nx / 2);
          const InfluenceLine il =
              influence_line(m, path, Vec3(0, 0, -1), ResponseKind::bar_force(bar, 0.0, 4));
          t = now() - t0;
          const auto ex = il.point_extreme();
          return fmt("dofs %.0f positions %.0f | eta max %.4e min %.4e", 6.0 * m.nodes.size(),
                     static_cast<double>(il.s.size()), ex.first, ex.second);
        }
        std::vector<Id> nodes;
        for (int j = 0; j < ng; ++j)
          for (int i = 0; i <= nx; ++i) nodes.push_back(id(i, j));
        const double t0 = now();
        const NodeGrid grid = build_node_grid(m, nodes, nx + 1, ng);
        const InfluenceSurface s =
            influence_surface(m, grid, Vec3(0, 0, -1), ResponseKind::reaction(id(nx / 3, jc), 2));
        t = now() - t0;
        const auto ex = s.point_extreme();
        return fmt("dofs %.0f positions %.0f | eta max %.4e min %.4e", 6.0 * m.nodes.size(),
                   static_cast<double>(s.eta.size()), ex.first, ex.second);
      });
    } else if (sc == "diag" || sc == "uplift" || sc == "static" || sc == "envelope" || sc == "convoy") {
      sweep(sc.c_str(), sizes_or(argc, argv, {240, 480, 960, 1920, 3840}), [&](int nx, double& t) {
        const int ng = 5;
        if (sc == "diag") {
          const Model m = grillage(nx, ng);
          const double t0 = now();
          const auto d = diagnose(m, false);
          t = now() - t0;
          return fmt("nodes %.0f | diagnostics %.0f", static_cast<double>(m.nodes.size()), static_cast<double>(d.size()));
        }
        if (sc == "uplift") {
          const Model m = grillage_uplift(nx, ng);
          const double t0 = now();
          const StaticResults r = solve_static(m, StaticOptions{});
          t = now() - t0;
          double lifted = 0;
          for (const auto& re : r.cases[0].reactions) lifted += re.R(2) == 0.0;
          return fmt("supports %.0f | zero reactions %.0f", static_cast<double>(m.supports.size()), lifted);
        }
        Model m = grillage(nx, ng);
        m.cases.clear();
        const int nc = sc == "envelope" ? 32 : 8;
        for (int c = 0; c < nc; ++c) {
          std::vector<Load> ls;
          for (int j = 0; j < ng; ++j) ls.push_back(Load::nodal(static_cast<Id>(1 + (c * nx / nc) + (nx + 1) * j), Vec3(0, 0, -1e5)));
          m.add_case("c" + std::to_string(c), ls);
        }
        if (sc == "static") {
          const double t0 = now();
          const StaticResults r = solve_static(m, StaticOptions{});
          t = now() - t0;
          return fmt("dofs %.0f | uz %.4e", 6.0 * m.nodes.size(), r.cases[0].node_disp[static_cast<std::size_t>(nx / 2)](2));
        }
        if (sc == "envelope") {
          const StaticResults r = solve_static(m, StaticOptions{});
          std::vector<const CaseResult*> src;
          for (const auto& c : r.cases) src.push_back(&c);
          const double t0 = now();
          const Envelope e = envelope(src);
          t = now() - t0;
          (void)e;
          return fmt("cases %.0f bars %.0f", static_cast<double>(nc), static_cast<double>(m.bars.size()));
        }
        // convoy
        std::vector<Id> path_bars;
        for (int i = 0; i < nx; ++i) path_bars.push_back(static_cast<Id>(1 + (ng / 2) * nx + i));
        const LoadPath path = build_path(m, path_bars, 11);
        const InfluenceLine il = influence_line(m, path, Vec3(0, 0, -1),
                                                ResponseKind::bar_force(static_cast<Id>(1 + (ng / 2) * nx + nx / 2), 0.0, 4));
        ConvoyPattern cv;
        cv.axles = {{0.0, 1.2e5}, {1.35, 1.2e5}};
        cv.segments = {{3.0, 20.0, 1e5}};
        const double t0 = now();
        const Extremes ex = convoy_extreme(il, cv, static_cast<int>(il.s.size()));
        t = now() - t0;
        return fmt("positions %.0f | max %.4e min %.4e", static_cast<double>(il.s.size()), ex.max, ex.min);
      });
    } else if (sc == "shellval") {
      sweep("shellval", sizes_or(argc, argv, {250, 500, 1000, 2000, 4000}), [](int n, double& t) {
        const Model m = shells_imposed(n);
        const double t0 = now();
        m.validate();
        t = now() - t0;
        return fmt("shells %.0f loads %.0f", static_cast<double>(m.shells.size()), static_cast<double>(m.cases[0].loads.size()));
      });
    } else if (sc == "mesh") {
      sweep("mesh", sizes_or(argc, argv, {25, 50, 100, 200, 400}), [](int n, double& t) {
        mesh2::MeshInput in;
        in.outline = {{0, 0}, {10, 0}, {10, 10}, {0, 10}};
        for (const auto& [cx, cy] : {std::pair{3.0, 3.0}, std::pair{7.0, 6.5}}) {
          std::vector<mesh2::Vec2> h;
          for (int k = 0; k < 48; ++k) h.push_back({cx + std::cos(-2 * M_PI * k / 48), cy + std::sin(-2 * M_PI * k / 48)});
          in.holes.push_back(h);
        }
        in.segments = {{{1.0, 8.5}, {9.0, 9.0}}};
        in.target_h = 10.0 / n;
        const double t0 = now();
        const mesh2::TriMesh tm = mesh2::triangulate(in);
        t = now() - t0;
        return fmt("nodes %.0f triangles %.0f | area %.6f", static_cast<double>(tm.nodes.size()),
                   static_cast<double>(tm.tris.size()), tm.area());
      });
    } else if (sc == "st1") {
      sweep("st1", sizes_or(argc, argv, {500, 1000, 2000, 4000, 8000}), [](int n, double& t) {
        const std::string text = st1_beam(n);
        const double t0 = now();
        const St1Import imp = import_st1(text);
        t = now() - t0;
        return fmt("bars %.0f cases %.0f | script %.0f bytes", static_cast<double>(imp.model.bars.size()),
                   static_cast<double>(imp.model.cases.size()), static_cast<double>(text.size()));
      });
    } else {
      std::fprintf(stderr, "unknown scenario '%s'\n", sc.c_str());
      return 2;
    }
    return 0;
  }
  const double s = argc > 1 ? std::atof(argv[1]) : 1.0;
  if (argc > 2) set_solver_threads(std::atoi(argv[2]));  // threads of the factorization
  std::printf("solver threads: %d\n", solver_threads());
  run("frame", frame(static_cast<int>(10 * s), static_cast<int>(10 * s), static_cast<int>(20 * s)));
  run("deck", deck(static_cast<int>(120 * s), static_cast<int>(24 * s)));
  return 0;
}

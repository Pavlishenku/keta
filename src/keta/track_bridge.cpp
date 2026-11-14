#include "keta/track_bridge.h"

#include <algorithm>
#include <cmath>

namespace keta {

namespace {
constexpr double kRigid = 1e12;  // vertical interaction stiffness per unit length [N/m/m]
constexpr double kBearingZ = 1e14;

struct Spring {
  Id link = 0;
  double x = 0, L = 0;
  bool deck = false;
};

struct Built {
  Model m;
  std::vector<Id> rail_bars;
  std::vector<double> xr;
  std::vector<Spring> springs;
  std::vector<Id> bearings;  // 2 per deck
  Id thermal = 0;
  std::vector<LinkHistory::Component> comps;
};

Built build(const TrackBridgeInput& in) {
  require(in.element > 0 && in.emb_left > 0 && in.emb_right > 0, "track_bridge: element and embankments > 0");
  require(in.rail_A > 0 && in.rail_I > 0 && in.rail_E > 0, "track_bridge: rail A, I, E > 0");
  require(in.k_unloaded > 0 && in.r_unloaded > 0 && in.k_loaded > 0 && in.r_loaded > 0,
          "track_bridge: ballast stiffness and resistance > 0");
  require(!in.decks.empty(), "track_bridge: at least one deck");
  Built b;
  Model& m = b.m;
  Material rail;
  rail.E = in.rail_E;
  rail.nu = 0.3;
  rail.alpha = in.rail_alpha;
  const Id mr = m.add_material(rail);
  Section sr;
  sr.props.A = in.rail_A;
  sr.props.Iy = sr.props.Iz = in.rail_I;
  sr.props.J = in.rail_I;
  const Id sir = m.add_section(sr);
  // breakpoints: embankment, decks, embankment, each cut in equal elements
  std::vector<double> xs{0.0};
  auto cut = [&](double a, double L) {
    const int n = std::max(1, static_cast<int>(std::ceil(L / in.element - 1e-9)));
    for (int i = 1; i <= n; ++i) xs.push_back(a + L * i / n);
  };
  cut(0, in.emb_left);
  std::vector<double> d0;  // deck starts
  double x = in.emb_left;
  for (const TrackDeck& d : in.decks) {
    require(d.length > 0 && d.A > 0 && d.I > 0 && d.E > 0 && d.e_track >= 0 && d.e_bearing >= 0,
            "track_bridge: deck length, A, I, E > 0, eccentricities >= 0");
    d0.push_back(x);
    cut(x, d.length);
    x += d.length;
  }
  cut(x, in.emb_right);
  b.xr = xs;
  const std::size_t n = xs.size() - 1;
  std::vector<Id> r;
  for (double xi : xs) r.push_back(m.add_node(xi, 0, 0));
  for (std::size_t i = 0; i < n; ++i) b.rail_bars.push_back(m.add_bar(r[i], r[i + 1], mr, sir));
  for (std::size_t i = 0; i <= n; ++i)
    m.add_support(r[i], {i == 0 || i == n, true, false, true, false, true});
  auto near = [](double a, double c) { return std::fabs(a - c) < 1e-9 * std::max(1.0, std::fabs(c)); };
  const PlasticLink unl{in.k_unloaded, in.r_unloaded, 0.0};
  auto spring = [&](Id ground_or_deck, Id rail_node, double L, double xi, bool deck) {
    const Id l = m.add_elastic_link(ground_or_deck, rail_node, ElasticLink::diagonal(0, 0, kRigid * L, 0, 0, 0));
    m.links.back().axes = Mat3::Identity();
    b.comps.push_back({l, 0, {unl.k * L, unl.fy * L, 0.0}});
    b.springs.push_back({l, xi, L, deck});
  };
  // decks
  std::vector<Load> th;
  if (in.rail_dT != 0)
    for (Id bar : b.rail_bars) th.push_back(Load::thermal(bar, in.rail_dT));
  std::vector<std::vector<std::pair<double, Id>>> tops(in.decks.size());  // (x, top node)
  for (std::size_t k = 0; k < in.decks.size(); ++k) {
    const TrackDeck& d = in.decks[k];
    Material mat;
    mat.E = d.E;
    mat.nu = 0.3;
    mat.alpha = d.alpha;
    const Id md = m.add_material(mat);
    Section s;
    s.props.A = d.A;
    s.props.Iy = s.props.Iz = d.I;
    s.props.J = d.I;
    const Id sd = m.add_section(s);
    std::vector<Id> na;
    for (double xi : xs)
      if (xi >= d0[k] - 1e-9 && xi <= d0[k] + d.length + 1e-9) {
        const Id nn = m.add_node(xi, 0, -d.e_track);
        m.add_support(nn, {false, true, false, true, false, true});
        na.push_back(nn);
        if (d.e_track > 0) {
          const Id t = m.add_node(xi, 0, 0);
          m.add_rigid_link(nn, {t});
          tops[k].push_back({xi, t});
        } else {
          tops[k].push_back({xi, nn});
        }
      }
    for (std::size_t i = 0; i + 1 < na.size(); ++i) {
      const Id bar = m.add_bar(na[i], na[i + 1], md, sd);
      if (d.dT != 0) th.push_back(Load::thermal(bar, d.dT));
    }
    for (int end = 0; end < 2; ++end) {
      const Id master = end == 0 ? na.front() : na.back();
      const Vec3 p = m.node(master).x;
      const Id bn = m.add_node(p.x(), 0, -d.e_track - d.e_bearing);
      if (d.e_bearing > 0) m.add_rigid_link(master, {bn});
      const Id g = m.add_node(p.x(), 0, -d.e_track - d.e_bearing);
      m.add_support(g, {true, true, true, true, true, true});
      const Id target = d.e_bearing > 0 ? bn : master;
      b.bearings.push_back(
          m.add_elastic_link(g, target, ElasticLink::diagonal(end == 0 ? d.k_left : d.k_right, 0, kBearingZ, 0, 0, 0)));
      m.links.back().axes = Mat3::Identity();
    }
  }
  // interaction springs, node by node: each half of the tributary length
  // goes to the support under it (embankment or the deck spanning it)
  auto deck_at = [&](double a, double c) -> int {  // the deck holding the segment [a, c]
    for (std::size_t k = 0; k < in.decks.size(); ++k)
      if (a >= d0[k] - 1e-9 && c <= d0[k] + in.decks[k].length + 1e-9) return static_cast<int>(k);
    return -1;
  };
  for (std::size_t i = 0; i <= n; ++i) {
    double Lg = 0;
    for (int side : {-1, 1}) {
      const long o = static_cast<long>(i) + side;
      if (o < 0 || o > static_cast<long>(n)) continue;
      const double a = std::min(xs[i], xs[static_cast<std::size_t>(o)]),
                   c = std::max(xs[i], xs[static_cast<std::size_t>(o)]);
      const double half = (c - a) / 2;
      const int k = deck_at(a, c);
      if (k < 0) {
        Lg += half;
        continue;
      }
      Id top = 0;
      for (const auto& [xt, t] : tops[static_cast<std::size_t>(k)])
        if (near(xt, xs[i])) top = t;
      require(top != 0, "track_bridge: internal error (deck top node)");
      spring(top, r[i], half, xs[i], true);
    }
    if (Lg > 0) {
      const Id g = m.add_node(xs[i], 0, 0);
      m.add_support(g, {true, true, true, true, true, true});
      spring(g, r[i], Lg, xs[i], false);
    }
  }
  b.thermal = m.add_case("thermal", th);
  return b;
}

TrackState measure(const Built& b, const TrackBridgeInput& in, const LinkHistory& h) {
  TrackState s;
  const CaseResult c = h.result(2);
  s.rail_min = 1e300;
  s.rail_max = -1e300;
  for (Id bar : b.rail_bars)
    for (const Station& st : c.bars[static_cast<std::size_t>(b.m.bar_index(bar))].stations) {
      s.rail_min = std::min(s.rail_min, st.N / in.rail_A);
      s.rail_max = std::max(s.rail_max, st.N / in.rail_A);
    }
  for (Id l : b.bearings)
    for (const LinkResult& lr : c.links)
      if (lr.link == l) s.bearing.push_back(lr.force(0));
  for (const Spring& sp : b.springs) {
    const double d = std::fabs(h.deformation(sp.link, 0));
    (sp.deck ? s.rel_deck : s.rel_ground) = std::max(sp.deck ? s.rel_deck : s.rel_ground, d);
  }
  return s;
}
}  // namespace

TrackBridgeResult track_bridge_uic774(const TrackBridgeInput& in, const TrainPass& train, int steps) {
  require(train.length > 0 && train.positions >= 1, "track_bridge: train length > 0, positions >= 1");
  Built b = build(in);
  // the train cases (one per position), defined before the history is built
  std::vector<std::pair<double, double>> spans;
  std::vector<Id> cases;
  for (int p = 0; p < train.positions; ++p) {
    const double f = train.positions == 1 ? train.front_first
                                          : train.front_first + (train.front_last - train.front_first) * p /
                                                                    (train.positions - 1);
    const double lo = train.forward ? f - train.length : f, hi = train.forward ? f : f + train.length;
    spans.emplace_back(lo, hi);
    std::vector<Load> L;
    for (std::size_t i = 0; i + 1 < b.xr.size(); ++i) {
      const double a = std::max(lo, b.xr[i]), c = std::min(hi, b.xr[i + 1]);
      if (c - a <= 1e-9) continue;
      L.push_back(Load::udl(b.rail_bars[i], Vec3(train.horizontal, 0, -train.vertical), LoadFrame::Global,
                            a - b.xr[i], c - b.xr[i]));
    }
    cases.push_back(b.m.add_case("train " + std::to_string(p), L));
  }
  LinkHistory base(b.m, b.comps);
  base.apply(b.thermal, steps);
  TrackBridgeResult out;
  out.thermal = measure(b, in, base);
  out.envelope.rail_min = 1e300;
  out.envelope.rail_max = -1e300;
  out.envelope.bearing.assign(2 * in.decks.size(), 0.0);
  for (std::size_t p = 0; p < cases.size(); ++p) {
    LinkHistory h = base;
    const auto [lo, hi] = spans[p];
    // the loaded resistance on the springs whose tributary length lies under the train
    for (const Spring& sp : b.springs)
      if (sp.x >= lo - 1e-9 && sp.x <= hi + 1e-9)
        h.set_law(sp.link, 0, {in.k_loaded * sp.L, in.r_loaded * sp.L, 0.0});
    h.apply(cases[p], steps);
    const TrackState s = measure(b, in, h);
    out.envelope.rail_min = std::min(out.envelope.rail_min, s.rail_min);
    out.envelope.rail_max = std::max(out.envelope.rail_max, s.rail_max);
    out.envelope.rel_deck = std::max(out.envelope.rel_deck, s.rel_deck);
    out.envelope.rel_ground = std::max(out.envelope.rel_ground, s.rel_ground);
    for (std::size_t j = 0; j < s.bearing.size(); ++j)
      out.envelope.bearing[j] = std::max(out.envelope.bearing[j], std::fabs(s.bearing[j]));
    out.positions.push_back(s);
  }
  return out;
}

}  // namespace keta

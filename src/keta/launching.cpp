#include "keta/launching.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>

#include "keta/influence.h"
#include "keta/statics.h"

namespace keta {

namespace {

constexpr double kPi = 3.14159265358979323846;
// A contact closer to a chain node than this fraction of its bar bears on
// the node. A cut at a distance d (fraction of the bar) loses about
// 1e-16 / d of the reactions to round-off (the short piece is a near-rigid
// element), bearing on the node misplaces them by about d: measured on a
// 20 m girder of 4 m bars, both stay near 1e-7 of the weight at 1e-7.
constexpr double kContactSnap = 1e-7;

void check_path(const LaunchPath& p) {
  require(std::isfinite(p.heading) && std::isfinite(p.radius) && std::isfinite(p.grade) &&
              std::isfinite(p.origin.x()) && std::isfinite(p.origin.y()) && std::isfinite(p.origin.z()),
          "launch path: non-finite data");
  switch (p.kind) {
    case LaunchPath::Kind::Straight: break;
    case LaunchPath::Kind::Circular:
      require(p.radius != 0.0, "launch path: a circular path needs a non-zero radius");
      require(p.grade == 0.0, "launch path: a circular path is level (a grade makes it helical)");
      break;
    case LaunchPath::Kind::Helical:
      require(p.radius != 0.0, "launch path: a helical path needs a non-zero radius");
      break;
  }
}

// 3-D arc length per plan length
double stretch(const LaunchPath& p) { return std::sqrt(1.0 + p.grade * p.grade); }

// centre of the plan circle
Vec3 centre(const LaunchPath& p) {
  return p.origin + p.radius * Vec3(-std::sin(p.heading), std::cos(p.heading), 0.0);
}

// The path abscissa of a point on it and its distance to the path; `near`:
// the plan angle (curved paths) or abscissa to stay closest to.
double abscissa_of(const LaunchPath& p, const Vec3& x, double near_s) {
  const double k = stretch(p);
  if (p.kind == LaunchPath::Kind::Straight) {
    const Vec3 tp(std::cos(p.heading), std::sin(p.heading), 0.0);
    return (x - p.origin).dot(tp) * k;
  }
  // angle of (x - C) in plan = heading + phi - sgn(R) pi/2, phi = s_plan / R
  const Vec3 r = x - centre(p);
  const double sgn = p.radius > 0 ? 1.0 : -1.0;
  double phi = std::atan2(r.y(), r.x()) - p.heading + sgn * kPi / 2;
  // the branch closest to the expected angle (unwraps a helix turn by turn)
  const double phi_near = near_s / k / p.radius;
  phi += 2 * kPi * std::round((phi_near - phi) / (2 * kPi));
  return p.radius * phi * k;
}

}  // namespace

LaunchPoint launch_path_point(const LaunchPath& p, double s) {
  check_path(p);
  const double k = stretch(p), sp = s / k;
  LaunchPoint out;
  if (p.kind == LaunchPath::Kind::Straight) {
    const Vec3 tp(std::cos(p.heading), std::sin(p.heading), 0.0);
    out.x = p.origin + sp * tp + Vec3(0, 0, p.grade * sp);
    out.t = (tp + Vec3(0, 0, p.grade)) / k;
    return out;
  }
  const double a = p.heading + sp / p.radius;
  out.x = centre(p) - p.radius * Vec3(-std::sin(a), std::cos(a), 0.0) + Vec3(0, 0, p.grade * sp);
  out.t = Vec3(std::cos(a), std::sin(a), p.grade) / k;
  return out;
}

namespace {

// A bearing of one position: its contact abscissa on the girder (from the
// first chain node), the chain interval holding it, the node it lands on in
// the working model, the first pier of the list that bears there.
struct Bearing {
  double s = 0;
  std::size_t interval = 0;  // chain bar index k: sig[k] <= s <= sig[k + 1]
  bool on_node = false;      // within the snap distance of chain node `interval` (or the last one)
  Id node = 0;
  std::size_t pier = 0;
};

// Split the chain bar `bar` of `w` at the local abscissae `cuts` (from its
// node i, ascending, strictly inside): the bar keeps its id from node i to
// the first cut, new bars carry the rest; each piece keeps the material,
// section, roll, bed and composite parts of the bar, the end releases go to
// the end pieces. Returns the ids of the pieces in order from node i and
// the new nodes at the cuts.
void split_bar(Model& w, Id bar, const std::vector<double>& cuts, std::vector<Id>& pieces, std::vector<Id>& cut_nodes) {
  std::size_t ix = 0;
  while (w.bars[ix].id != bar) ++ix;
  const Bar orig = w.bars[ix];
  require(orig.kind == BarKind::Frame && !orig.has_ecc() && w.tapered_of(bar) == nullptr,
          "launch: a pier bears inside bar " + std::to_string(bar) +
              ", which is not a prismatic frame bar on its node line (variable section, rigid arms or "
              "cable): put a node at the contact");
  for (const Tendon& t : w.tendons)
    require(std::find(t.bars.begin(), t.bars.end(), bar) == t.bars.end(),
            "launch: a pier bears inside bar " + std::to_string(bar) + ", run by a tendon: put a node at the contact");
  const Vec3 xi = w.node(orig.ni).x, xj = w.node(orig.nj).x;
  const double L = (xj - xi).norm();
  Id next = 0;
  for (const Bar& b : w.bars) next = std::max(next, b.id);
  cut_nodes.clear();
  for (double a : cuts) {
    const Vec3 x = xi + (a / L) * (xj - xi);
    cut_nodes.push_back(w.add_node(x.x(), x.y(), x.z()));
  }
  pieces.assign(1, bar);
  Releases none{};
  Bar& first = w.bars[ix];
  first.nj = cut_nodes.front();
  for (int c = 6; c < 12; ++c) first.releases[static_cast<std::size_t>(c)] = false;
  for (std::size_t q = 0; q < cut_nodes.size(); ++q) {
    Bar piece = orig;
    piece.id = ++next;
    piece.ni = cut_nodes[q];
    piece.nj = q + 1 < cut_nodes.size() ? cut_nodes[q + 1] : orig.nj;
    piece.releases = none;
    if (q + 1 == cut_nodes.size())
      for (int c = 6; c < 12; ++c) piece.releases[static_cast<std::size_t>(c)] = orig.releases[static_cast<std::size_t>(c)];
    w.bars.push_back(piece);
    if (auto it = w.composite.find(bar); it != w.composite.end()) w.composite[piece.id] = it->second;
    pieces.push_back(piece.id);
  }
  w.invalidate_index();
}

}  // namespace

LaunchPathResults launch_on_path(const Model& girder, const std::vector<Id>& chain, const LaunchPath& path,
                                 const std::vector<double>& pier_abscissae, const std::vector<double>& advances,
                                 const LaunchOptions& opt) {
  check_path(path);
  require(!chain.empty() && !pier_abscissae.empty() && !advances.empty(), "launch: empty input");
  require(opt.stations_per_bar >= 2 && opt.path_tol > 0, "launch: stations_per_bar >= 2 and path_tol > 0");
  // the chain nodes in launching order, their path abscissae (the girder
  // must lie on its path: it slides along it)
  const LoadPath lp = build_path(girder, chain, 2);
  std::vector<Id> nodes;
  for (const PathStation& st : lp.stations) {
    const Bar& b = girder.bar(st.bar);
    nodes.push_back(st.x_local < girder.bar_chord_length(b) / 2 ? b.ni : b.nj);
  }
  require(nodes.size() == chain.size() + 1, "launch: the chain must be a simple sequence of bars");
  std::vector<double> sig(nodes.size());
  {
    // the first node: on a helix its height fixes the turn
    const Vec3 x0 = girder.node(nodes.front()).x;
    double s_prev = 0.0;
    if (path.kind == LaunchPath::Kind::Helical && path.grade != 0.0)
      s_prev = (x0.z() - path.origin.z()) / path.grade * stretch(path);
    for (std::size_t k = 0; k < nodes.size(); ++k) {
      const Vec3 x = girder.node(nodes[k]).x;
      const double s = abscissa_of(path, x, s_prev);
      const double off = (launch_path_point(path, s).x - x).norm();
      require(off <= opt.path_tol, "launch: node " + std::to_string(nodes[k]) + " lies " + std::to_string(off) +
                                       " m off the launch path (a girder slides along its own path only)");
      sig[k] = s;
      s_prev = s;
    }
  }
  const double s0 = sig.front();
  for (double& s : sig) s -= s0;
  for (std::size_t k = 1; k < sig.size(); ++k)
    require(sig[k] > sig[k - 1], "launch: the chain must run along the path direction, rear end first");
  const double Lg = sig.back();
  const std::size_t nb = chain.size();
  const int nst = opt.stations_per_bar;

  LaunchPathResults out;
  bool first_env = true;
  for (double adv : advances) {
    LaunchPosition pos;
    pos.advance = adv;
    const std::size_t np = pier_abscissae.size();
    pos.contact_s.assign(np, std::numeric_limits<double>::quiet_NaN());
    pos.reaction.assign(np, 0.0);
    pos.guide.assign(np, 0.0);
    pos.moment.assign(np, 0.0);
    // the bearings: one per distinct contact abscissa (piers closer than the
    // snap distance share one, the first of the list carries it)
    std::vector<Bearing> bearings;
    std::vector<int> bearing_of(np, -1);
    for (std::size_t i = 0; i < np; ++i) {
      const double s = pier_abscissae[i] - adv + Lg;
      if (s < -1e-9 || s > Lg + 1e-9) continue;
      Bearing b;
      b.s = std::clamp(s, 0.0, Lg);
      b.pier = i;
      b.interval = static_cast<std::size_t>(std::upper_bound(sig.begin(), sig.end(), b.s) - sig.begin());
      b.interval = std::min(b.interval == 0 ? 0 : b.interval - 1, nb - 1);
      const std::size_t k = b.interval;
      const double tol = kContactSnap * (sig[k + 1] - sig[k]);
      if (b.s - sig[k] <= tol) {
        b.on_node = true;
        b.s = sig[k];
      } else if (sig[k + 1] - b.s <= tol) {
        b.on_node = true;
        b.s = sig[k + 1];
        b.interval = k + 1;  // the node index
      }
      pos.contact_s[i] = b.s;
      int same = -1;
      for (std::size_t q = 0; q < bearings.size(); ++q)
        if (std::fabs(bearings[q].s - b.s) <= tol) same = static_cast<int>(q);
      if (same < 0) {
        bearing_of[i] = static_cast<int>(bearings.size());
        bearings.push_back(b);
      } else {
        bearing_of[i] = same;
        pos.contact_s[i] = bearings[static_cast<std::size_t>(same)].s;
      }
    }
    pos.stable = bearings.size() >= 2;
    if (!bearings.empty()) {
      double smin = Lg, smax = 0;
      for (std::size_t i = 0; i < np; ++i)
        if (bearing_of[i] >= 0) {
          smin = std::min(smin, pos.contact_s[i]);
          smax = std::max(smax, pos.contact_s[i]);
        }
      pos.front_overhang = Lg - smax;
      pos.rear_overhang = smin;
    }
    if (!pos.stable) {
      out.positions.push_back(std::move(pos));
      continue;
    }
    Model w = girder;
    w.supports.clear();
    w.skew_supports.clear();
    w.cases.clear();
    w.combos.clear();
    // a node at every contact between two chain nodes: the chain bar is cut
    // there (several cuts per bar if several piers bear inside it)
    std::map<Id, std::vector<double>> stations;  // the pieces' station abscissae
    std::map<Id, std::vector<std::pair<Id, std::size_t>>> station_of;  // original bar: (piece, station) per station
    for (std::size_t k = 0; k < nb; ++k) {
      std::vector<std::size_t> inside;
      for (std::size_t q = 0; q < bearings.size(); ++q)
        if (!bearings[q].on_node && bearings[q].interval == k) inside.push_back(q);
      if (inside.empty()) continue;
      const Bar& b = girder.bar(chain[k]);
      const bool forward = b.ni == nodes[k];
      const double L = girder.bar_chord_length(b);
      // local abscissa (from node i) of each cut: the same fraction of the
      // chord as of the path interval
      std::vector<std::pair<double, std::size_t>> cut;
      for (std::size_t q : inside) {
        const double t = (bearings[q].s - sig[k]) / (sig[k + 1] - sig[k]);
        cut.emplace_back((forward ? t : 1.0 - t) * L, q);
      }
      std::sort(cut.begin(), cut.end());
      std::vector<double> a;
      for (const auto& c : cut) a.push_back(c.first);
      std::vector<Id> pieces, cut_nodes;
      split_bar(w, chain[k], a, pieces, cut_nodes);
      for (std::size_t c = 0; c < cut.size(); ++c) bearings[cut[c].second].node = cut_nodes[c];
      // the original stations, each read on the piece holding it; both ends
      // of every piece (the contacts' moments)
      std::vector<double> from(1, 0.0);
      for (double x : a) from.push_back(x);
      from.push_back(L);
      for (std::size_t p = 0; p < pieces.size(); ++p) stations[pieces[p]] = {0.0, from[p + 1] - from[p]};
      std::vector<std::pair<Id, double>> where;
      for (int j = 0; j < nst; ++j) {
        const double x = L * j / (nst - 1);
        std::size_t p = 0;
        while (p + 1 < pieces.size() && x > from[p + 1]) ++p;
        const double xl = std::clamp(x - from[p], 0.0, from[p + 1] - from[p]);
        where.emplace_back(pieces[p], xl);
        stations[pieces[p]].push_back(xl);
      }
      for (auto& [id, xs] : stations) {
        std::sort(xs.begin(), xs.end());
        xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
      }
      auto& so = station_of[chain[k]];
      for (const auto& [piece, xl] : where) {
        const std::vector<double>& xs = stations[piece];
        so.emplace_back(piece, static_cast<std::size_t>(std::lower_bound(xs.begin(), xs.end(), xl) - xs.begin()));
      }
    }
    for (Bearing& b : bearings)
      if (b.on_node) b.node = nodes[b.interval];
    // one bearing per contact, in the path frame there; bearings are in
    // pier order, the first engaged pier carries the jacks
    for (std::size_t q = 0; q < bearings.size(); ++q) {
      const Vec3 t = launch_path_point(path, bearings[q].s + s0).t;
      const Vec3 tp = Vec3(t.x(), t.y(), 0.0).normalized();
      Mat3 axes;
      axes.row(0) = tp.transpose();
      axes.row(1) = Vec3(-tp.y(), tp.x(), 0.0).transpose();
      axes.row(2) = Vec3(0, 0, 1).transpose();
      std::array<bool, 6> fix{};
      fix[1] = fix[2] = true;            // lateral guide, vertical bearing
      if (q == 0) fix[0] = fix[3] = true;  // the jacks: longitudinal and twist
      w.add_skew_support(bearings[q].node, axes, fix);
    }
    w.cases.push_back({1, "launch", {Load::self_weight()}});
    StaticOptions so;
    so.stations_per_bar = nst;
    so.stations_abs = &stations;
    const StaticResults r = solve_static(w, so);
    const CaseResult& cr = r.cases.front();
    std::map<Id, const BarResult*> res;
    for (const BarResult& br : cr.bars) res[br.bar] = &br;
    // My at a contact: on a piece of the chain starting (else ending) there
    std::set<Id> chain_pieces(chain.begin(), chain.end());
    for (const auto& [id, xs] : stations) chain_pieces.insert(id);
    auto moment_at = [&](Id n) {
      for (const Bar& b : w.bars)
        if (chain_pieces.count(b.id) && b.ni == n) return res.at(b.id)->stations.front().My;
      for (const Bar& b : w.bars)
        if (chain_pieces.count(b.id) && b.nj == n) return res.at(b.id)->stations.back().My;
      return 0.0;
    };
    double W = 0;
    for (std::size_t q = 0; q < bearings.size(); ++q) {
      const Vec6& loc = cr.skew_reactions[q].local;
      pos.reaction[bearings[q].pier] = loc(2);
      pos.guide[bearings[q].pier] = loc(1);
      W += loc(2);
    }
    for (std::size_t i = 0; i < np; ++i)
      if (bearing_of[i] >= 0) pos.moment[i] = moment_at(bearings[static_cast<std::size_t>(bearing_of[i])].node);
    for (std::size_t i = 0; i < np; ++i)
      if (pos.reaction[i] < -1e-9 * std::fabs(W)) pos.uplift.push_back(i);
    // the envelope on the girder's own station grid
    std::size_t k = 0;
    for (const Bar& b : girder.bars) {
      const auto split = station_of.find(b.id);
      for (int j = 0; j < nst; ++j, ++k) {
        double My = 0;
        if (split == station_of.end()) {
          My = res.at(b.id)->stations[static_cast<std::size_t>(j)].My;
        } else {
          const auto& [piece, ix] = split->second[static_cast<std::size_t>(j)];
          My = res.at(piece)->stations[ix].My;
        }
        if (first_env) {
          out.m_min.push_back(My);
          out.m_max.push_back(My);
        } else {
          out.m_min[k] = std::min(out.m_min[k], My);
          out.m_max[k] = std::max(out.m_max[k], My);
        }
      }
    }
    first_env = false;
    out.positions.push_back(std::move(pos));
  }
  return out;
}

}  // namespace keta

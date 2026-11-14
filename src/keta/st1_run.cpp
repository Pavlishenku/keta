#include "keta/st1_run.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <set>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include "keta/st1.h"
#include "keta/bef.h"
#include "keta/st1_domains.h"
#include "keta/st1_edit.h"

namespace keta::st1 {
namespace {

// `APPUI ... DECOL NON` as parsed (cancels a DECOL declared before), apart
// from 0 = "no DECOL on this line".
constexpr int kDecolOff = 2;

// id -> the FIRST bar of that id in m.bars (what a std::find_if over m.bars
// returns), for per-bar lookups inside result loops.
std::unordered_map<Id, const Bar*> bars_by_id(const Model& m) {
  std::unordered_map<Id, const Bar*> out;
  out.reserve(m.bars.size());
  for (const Bar& b : m.bars) out.emplace(b.id, &b);
  return out;
}

// Position of bar `id` in m.bars through the model index, checked against
// the slot (-1: absent, or an index out of step with the vector).
int model_has_bar(const Model& m, Id id) {
  try {
    const int ix = m.bar_index(id);
    return m.bars[static_cast<std::size_t>(ix)].id == id ? ix : -1;
  } catch (const Error&) {
    return -1;
  }
}

std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

bool is_integer(const std::string& s) {
  if (s.empty()) return false;
  std::size_t k = (s[0] == '-' || s[0] == '+') ? 1 : 0;
  if (k == s.size()) return false;
  for (; k < s.size(); ++k)
    if (!std::isdigit(static_cast<unsigned char>(s[k]))) return false;
  return true;
}

// Ids are lexically integer tokens. The interpreter formats integer-valued
// expressions as plain integers ("nseg+1" -> "5"), so a decimal-shaped token
// ("1.0", "0.3") is always a value — never an id. This is what disambiguates
// 'CHARG 1 1.0' (id 1, coefficient 1.0) and 'BETA 3 0.3'.
bool is_id_tok(const Tok& t) { return !t.str && is_integer(t.s); }

// A LIST ITEM is always an UNSIGNED integer: ST1 numbers nodes, bars, supports,
// loadings… from 1. A signed token is therefore a value, never an item — which
// is what separates the coefficient run of `CHARG 5 -30 20` (source 5, box
// [-30, 20], manual p.234) from a three-element list. The `PAS` step and the
// range bound keep `is_id_tok`, since a descending range legitimately steps by
// a negative number.
bool is_uid_tok(const Tok& t) {
  return !t.str && is_integer(t.s) && t.s[0] != '-' && t.s[0] != '+';
}

bool is_num_tok(const Tok& t) { return !t.str && is_st1_number(t.s); }
double num(const Tok& t) { return st1_num(t.s); }

std::string fmt(double v) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

// ST1 `CARA ...` characteristic keywords -> CaraProps, one keyword + value at a
// time (manual p.59-60). `set` records which of the 13 fields were named so a
// redefinition only overwrites what it mentions.
bool cara_keyword(const std::string& key, double v, CaraProps& p, std::array<bool, 13>& set) {
  if (key == "SX") { p.sx = v; set[0] = true; }
  else if (key == "SY") { p.sy = v; set[1] = true; }
  else if (key == "SZ") { p.sz = v; set[2] = true; }
  else if (key == "IX") { p.ix = v; set[3] = true; }
  else if (key == "IY") { p.iy = v; set[4] = true; }
  else if (key == "IZ") { p.iz = v; set[5] = true; }
  else if (key == "IPOL") { p.ipol = v; set[6] = true; }
  else if (key == "VY") { p.vy = v; set[7] = true; }
  else if (key == "WY") { p.wy = v; set[8] = true; }
  else if (key == "VZ") { p.vz = v; set[9] = true; }
  else if (key == "WZ") { p.wz = v; set[10] = true; }
  else if (key == "EXTY") { p.exty = v; set[11] = true; }
  else if (key == "EXTZ") { p.extz = v; set[12] = true; }
  else return false;
  return true;
}

// Merge `src` onto `dst`, field by field, honouring the `set` mask.
void cara_merge(CaraProps& dst, const CaraProps& src, const std::array<bool, 13>& set) {
  if (set[0]) dst.sx = src.sx;
  if (set[1]) dst.sy = src.sy;
  if (set[2]) dst.sz = src.sz;
  if (set[3]) dst.ix = src.ix;
  if (set[4]) dst.iy = src.iy;
  if (set[5]) dst.iz = src.iz;
  if (set[6]) dst.ipol = src.ipol;
  if (set[7]) dst.vy = src.vy;
  if (set[8]) dst.wy = src.wy;
  if (set[9]) dst.vz = src.vz;
  if (set[10]) dst.wz = src.wz;
  if (set[11]) dst.exty = src.exty;
  if (set[12]) dst.extz = src.extz;
}

// ST1 characteristics -> keta section properties (tapered.h's own mapping:
// A <- SX, Avy <- SY, Avz <- SZ, J <- IX, Iy <- IY, Iz <- IZ).
// Same completion as the constant CARA (build_model): a PLANE deck gives IZ
// only, a GRILL deck IY only — the missing inertia and IX default like there
// (arc: `CARA VAR LIN Y 101 / OR SX .. IZ ..` has no IY).
SectionProps cara_to_props(const CaraProps& c) {
  SectionProps p;
  p.A = c.sx;
  p.Avy = c.sy;
  p.Avz = c.sz;
  p.Iz = c.iz > 0 ? c.iz : c.iy;
  p.Iy = c.iy > 0 ? c.iy : c.iz;
  p.J = c.ix > 0 ? c.ix : std::max(p.Iy, p.Iz);
  return p;
}

// ST1 GRILL orders a 3-DOF support / mass half-matrix as (rx, dz, ry) while the
// PLANE order is (dx, dy, rz) — the same lower triangle on keta DOFs (3, 2, 4)
// instead of (0, 1, 5) (manual p.51, p.71). Permute the PLANE-built matrix.
Mat6 to_grill_frame(const Mat6& plane) {
  static const int kPlane[3] = {0, 1, 5};
  static const int kGrill[3] = {3, 2, 4};
  Mat6 out = Mat6::Zero();
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) out(kGrill[a], kGrill[b]) = plane(kPlane[a], kPlane[b]);
  return out;
}

// A load scaled by `c`: every magnitude it carries (a COMB weight applied to
// the loads themselves, for the nonlinear solves).
Load scaled_load(Load l, double c) {
  l.F *= c;
  l.M *= c;
  l.F2 *= c;
  l.M2 *= c;
  l.dt *= c;
  l.gy *= c;
  l.gz *= c;
  l.g *= c;
  for (auto& d : l.disp)
    if (d) *d *= c;
  l.strain *= c;
  l.kappa_y *= c;
  l.kappa_z *= c;
  l.strain2 *= c;
  l.kappa_y2 *= c;
  l.kappa_z2 *= c;
  l.dvy *= c;
  l.dvz *= c;
  l.dthx *= c;
  l.grad *= c;
  return l;
}

}  // namespace

bool CaraProps::operator<(const CaraProps& o) const {
  return std::tie(sx, sy, sz, ix, iy, iz, ipol, vy, wy, vz, wz, exty, extz) <
         std::tie(o.sx, o.sy, o.sz, o.ix, o.iy, o.iz, o.ipol, o.vy, o.wy, o.vz, o.wz,
                  o.exty, o.extz);
}
bool CaraProps::operator==(const CaraProps& o) const { return !(*this < o) && !(o < *this); }
bool ConsProps::operator<(const ConsProps& o) const {
  return std::tie(E, nu, ro, alpha, g_dyn, amortissement, mat, mat_origin) <
         std::tie(o.E, o.nu, o.ro, o.alpha, o.g_dyn, o.amortissement, o.mat, o.mat_origin);
}

const CaseResult* ResultsStore::find(const std::string& kind, Id id) const {
  auto pick = [&](const std::map<Id, CaseResult>& m) -> const CaseResult* {
    auto it = m.find(id);
    return it == m.end() ? nullptr : &it->second;
  };
  if (kind == "CHARG") return pick(charg);
  if (kind == "COMB") return pick(comb);
  if (kind == "ETAT") return pick(etat);
  if (kind == "PHASAGE") return pick(phasage);
  if (kind == "SPECTRE") return pick(spectre);
  if (kind == "INCENDIE") return pick(incendie);
  if (kind == "PHASE_STEP") return pick(phase_step);
  return nullptr;
}

void Output::line(const std::string& s) {
  if (to_file && file && file->good())
    *file << s << '\n';
  else {
    console += s + "\n";
    if (on_console) on_console(s);
  }
}

// ================================================================= runtime --
void St1Runtime::stamp_cara(const std::string& kind, Id id) {
  if (!cara_snap_ || *cara_snap_ != cara_)
    cara_snap_ = std::make_shared<const std::map<Id, CaraProps>>(cara_);
  result_cara_[{kind, id}] = cara_snap_;
  const std::vector<Id> act = all_appui_ids();
  result_appuis_[{kind, id}] = std::set<Id>(act.begin(), act.end());
}

std::set<Id> St1Runtime::appuis_of(const std::string& kind, Id id) const {
  auto it = result_appuis_.find({kind, id});
  if (it != result_appuis_.end()) return it->second;
  const std::vector<Id> act = all_appui_ids();
  return std::set<Id>(act.begin(), act.end());
}

std::map<Id, SupportDef> St1Runtime::appui_defs_of(const std::string& kind, Id id) const {
  std::map<Id, SupportDef> d = appui_;
  auto it = result_appui_defs_.find({kind, id});
  if (it != result_appui_defs_.end())
    for (const auto& [a, s] : it->second) d[a] = s;
  return d;
}

std::map<Id, Vec6> St1Runtime::appui_reactions(const CaseResult& r, const std::set<Id>& active,
                                               bool zero_at_held_node,
                                               const std::map<Id, SupportDef>* defs) const {
  const std::map<Id, SupportDef>& appui_ = defs ? *defs : this->appui_;
  std::map<Id, Vec6> out;
  std::map<Id, std::array<bool, 6>> taken;  // node -> components already given
  for (Id aid : appui_order_) {
    if (!active.count(aid)) continue;
    auto ad = appui_.find(aid);
    if (ad == appui_.end()) continue;
    const SupportDef& s = ad->second;
    Vec6 R = Vec6::Zero();
    const Reaction* re = nullptr;
    for (const Reaction& x : r.reactions)
      if (x.node == s.node) re = &x;
    std::array<bool, 6>& tk = taken[s.node];
    for (std::size_t c = 0; c < 6; ++c) {
      const bool holds = s.fix[c] || s.spring[c] != 0 || !s.spring_lt.empty();
      if (!holds || tk[c]) continue;
      tk[c] = true;
      if (re) R(static_cast<int>(c)) = re->R(static_cast<int>(c));
    }
    out[aid] = R;
  }
  // ST1 keeps its reactions by node: an inactive APPUI on a node that another
  // active APPUI holds reads 0 there (e.g. an RZ support beside another
  // one on the same node enters the ENV with Min Mz = 0 from the loadings
  // where it was off), while one alone on its node is absent.
  if (zero_at_held_node)
    for (Id aid : appui_order_) {
      if (active.count(aid) || out.count(aid)) continue;
      auto ad = appui_.find(aid);
      if (ad != appui_.end() && taken.count(ad->second.node)) out[aid] = Vec6::Zero();
    }
  return out;
}

const std::map<Id, CaraProps>& St1Runtime::cara_of(const std::string& kind, Id id) const {
  auto it = result_cara_.find({kind, id});
  return it != result_cara_.end() && it->second ? *it->second : cara_;
}

std::vector<double> St1Runtime::contr_fibres(const std::map<Id, CaraProps>& cara, Id bid,
                                             const Station& s) const {
  std::vector<double> v(option_ == 2 ? 4 : 2, 0.0);
  auto it = cara.find(bid);
  if (it == cara.end()) return v;
  CaraProps cp = it->second;
  // CARA VAR (p.59): the section at the station — SX and inertias by the
  // tapered law, the fibre distances interpolated like the dimension that
  // varies (linear under LIN, parabolic through OR / MI / EX under PARA).
  if (auto vr = cara_var_.find(bid); vr != cara_var_.end() && vr->second.has_or && vr->second.has_ex) {
    const CaraVarDef& d = vr->second;
    double L = 0;
    if (auto rb = bars_.find(bid); rb != bars_.end())
      L = (nodes_.at(rb->second.n2) - nodes_.at(rb->second.n1)).norm();
    const double t = L > 0 ? std::clamp(s.x / L, 0.0, 1.0) : 0.0;
    TaperedVar var;
    var.law = d.law;
    var.orig = cara_to_props(d.orig);
    var.mid = cara_to_props(d.mid);
    var.ext = cara_to_props(d.ext);
    const SectionProps sp = tapered_props_at(var, t);
    auto interp = [&](double a, double m, double b) {
      if (d.law == TaperLaw::Para && d.has_mi)  // Lagrange through 0, 1/2, 1
        return a * (1 - t) * (1 - 2 * t) + 4 * m * t * (1 - t) + b * t * (2 * t - 1);
      return a + (b - a) * t;
    };
    cp.sx = sp.A;
    cp.iy = sp.Iy;
    cp.iz = sp.Iz;
    cp.vy = interp(d.orig.vy, d.mid.vy, d.ext.vy);
    cp.wy = interp(d.orig.wy, d.mid.wy, d.ext.wy);
    cp.vz = interp(d.orig.vz, d.mid.vz, d.ext.vz);
    cp.wz = interp(d.orig.wz, d.mid.wz, d.ext.wz);
  }
  if (cp.sx <= 0 && option_ != 1) return v;  // (GRILL: bending alone, no SX needed)
  auto sigma = [&](double y, double z) {
    double sg = cp.sx > 0 ? -s.N / cp.sx : 0.0;
    if (cp.iy > 0) sg += -s.My * z / cp.iy;
    if (cp.iz > 0) sg += s.Mz * y / cp.iz;
    return sg;
  };
  if (option_ == 0) v = {sigma(cp.vy, 0), sigma(-cp.wy, 0)};
  else if (option_ == 1) v = {sigma(0, cp.vz), sigma(0, -cp.wz)};
  else v = {sigma(cp.vy, cp.vz), sigma(-cp.wy, cp.vz), sigma(cp.vy, -cp.wz), sigma(-cp.wy, -cp.wz)};
  return v;
}


St1Runtime::St1Runtime(Env& env, const St1Options& opt) : env_(env), opt_(opt) {}
St1Runtime::~St1Runtime() = default;

void St1Runtime::warn(const std::string& s) {
  warnings_.push_back(s);
  if (out_.message_level >= 2) out_.line("### AVERTISSEMENT ### " + s);
}

std::vector<Id> St1Runtime::all_bar_ids() const {
  std::vector<Id> v;
  for (Id b : bar_order_)
    if (!selective_activation_ || active_bars_.count(b)) v.push_back(b);
  return v;
}
std::vector<Id> St1Runtime::all_node_ids() const {
  std::vector<Id> v;
  for (Id n : node_order_) v.push_back(n);
  return v;
}
double St1Runtime::bar_length(Id bar) const {
  auto rb = bars_.find(bar);
  if (rb == bars_.end()) return 0.0;
  auto n1 = nodes_.find(rb->second.n1), n2 = nodes_.find(rb->second.n2);
  return n1 != nodes_.end() && n2 != nodes_.end() ? (n2->second - n1->second).norm() : 0.0;
}
std::vector<Id> St1Runtime::all_appui_ids() const {
  std::vector<Id> v;
  for (Id a : appui_order_)
    if (!selective_activation_ || active_appuis_.count(a)) v.push_back(a);
  return v;
}

std::vector<Id> St1Runtime::id_list(const Statement& st, std::size_t& i,
                                    const std::vector<Id>& all, bool comma_only) const {
  std::vector<Id> out;
  bool expect_item = true;  // comma_only: an item is due (start or after ',')
  while (i < st.size()) {
    if (st[i].s == ",") {
      ++i;
      expect_item = true;
      continue;
    }
    if (comma_only && !expect_item) break;
    expect_item = false;
    std::string u = st[i].str ? "" : upper(st[i].s);
    if (u == "TOUT" || u == "TOUTES") {
      out = all;
      ++i;
      continue;
    }
    if (!is_uid_tok(st[i])) break;
    Id start = static_cast<Id>(num(st[i]));
    ++i;
    // range "i a j [pas p]" — commas may separate items around it
    std::size_t save = i;
    if (i < st.size() && upper(st[i].s) == "A") {
      ++i;
      if (i < st.size() && is_id_tok(st[i])) {
        Id end = static_cast<Id>(num(st[i]));
        ++i;
        Id step = 1;
        if (i + 1 < st.size() && upper(st[i].s) == "PAS" && is_id_tok(st[i + 1])) {
          step = static_cast<Id>(num(st[i + 1]));
          i += 2;
        }
        require(step != 0, "ST1: PAS nul dans une liste");
        if (step > 0)
          for (Id k = start; k <= end; k += step) out.push_back(k);
        else
          for (Id k = start; k >= end; k += step) out.push_back(k);
        continue;
      }
      i = save;  // 'A' belonged to something else
    }
    out.push_back(start);
  }
  return out;
}

// ------------------------------------------------------------------ nodes --
void St1Runtime::cmd_noeud(const Statement& st) {
  if (st.size() == 1) {
    mode_ = Mode::Nodes;
    return;
  }
  // The id list and the positional coordinates are both numeric runs; the
  // boundary is resolved by lexical shape: a token with '.' or an exponent
  // is a coordinate, an integer literal before the first such token is an
  // id. When EVERYTHING is integer-shaped, the trailing 2 (PLANE/GRILL) or
  // 3 (SPATIALE, when at least 4 numbers) values are the coordinates.
  bool has_kw = false;
  for (std::size_t k = 1; k < st.size(); ++k) {
    std::string u = st[k].str ? "" : upper(st[k].s);
    if (u == "X" || u == "Y" || u == "Z") has_kw = true;
  }
  Vec3 c = Vec3::Zero();
  std::vector<Id> ids;
  if (has_kw) {
    std::size_t i = 1;
    ids = id_list(st, i, all_node_ids());
    require(!ids.empty(), "ST1: NOEUD sans identifiant");
    while (i + 1 < st.size()) {
      std::string k = upper(st[i].s);
      if (!is_num_tok(st[i + 1])) break;
      double v = num(st[i + 1]);
      if (k == "X") c.x() = v;
      else if (k == "Y") c.y() = v;
      else if (k == "Z") c.z() = v;
      else break;
      i += 2;
    }
  } else {
    // Positional form: the COORDINATES are the trailing numeric tokens
    // (3 in SPATIALE, 2 otherwise), the id list is everything before. Walk
    // backwards, never crossing a range keyword ('A'/'PAS' binds its
    // following number to the list), and give one token back if the id part
    // would end up empty (SPATIALE decks may omit z).
    std::size_t arity = option_ == 2 ? 3 : 2;
    std::size_t split = st.size();
    std::size_t got = 0;
    while (split > 1 && got < arity) {
      const Tok& t = st[split - 1];
      if (t.str || !is_num_tok(t)) break;
      if (split >= 2) {
        std::string prev = st[split - 2].str ? "" : upper(st[split - 2].s);
        if (prev == "A" || prev == "PAS") break;  // range endpoint stays in the list
      }
      --split;
      ++got;
    }
    bool list_empty = split <= 1 || !(is_id_tok(st[1]) || upper(st[1].s) == "TOUT" ||
                                      upper(st[1].s) == "TOUTES");
    // If peeling consumed the whole numeric run, hand the first one back as
    // the node id.
    if (split == 1 && got > 2) {
      ++split;
      --got;
    } else if (list_empty && got > 2 && is_id_tok(st[split])) {
      ++split;
      --got;
    }
    require(got >= 2, "ST1: coordonnees de noeud incompletes");
    Statement listpart(st.begin(), st.begin() + static_cast<long>(split));
    std::size_t i = 1;
    ids = id_list(listpart, i, all_node_ids());
    require(!ids.empty(), "ST1: NOEUD sans identifiant");
    std::vector<double> coords;
    for (std::size_t k = split; k < st.size(); ++k)
      if (is_num_tok(st[k])) coords.push_back(num(st[k]));
    require(coords.size() >= 2, "ST1: coordonnees de noeud incompletes");
    c = Vec3(coords[0], coords[1], coords.size() > 2 ? coords[2] : 0.0);
  }
  for (Id n : ids) {
    if (!nodes_.count(n)) node_order_.push_back(n);
    nodes_[n] = c;
  }
  invalidate();
}

void St1Runtime::cmd_gener(const Statement& st) {
  require(st.size() >= 3 && is_num_tok(st[1]), "ST1: GENER invalide");
  int n = static_cast<int>(num(st[1]));
  std::string what = upper(st[2].s);
  std::size_t i = 3;
  auto grab = [&](double& v0, double& pas) {
    require(i < st.size() && is_num_tok(st[i]), "ST1: GENER — valeur attendue");
    v0 = num(st[i++]);
    pas = (i < st.size() && is_num_tok(st[i])) ? num(st[i++]) : 0.0;
  };
  // The ID keyword is optional: numbers right after NOEUD / BARRE are the
  // first id and its step (cable_et_phasage: `gener 10 barre 3 2 de 300 100 a
  // 301 100` makes bars 3, 5, ..., 21).
  auto leading_id = [&](double& j, double& pas0) {
    if (i < st.size() && is_num_tok(st[i])) grab(j, pas0);
  };
  if (what == "NOEUD") {
    double j = 1, pas0 = 1, x0 = 0, px = 0, y0 = 0, py = 0, z0 = 0, pz = 0;
    leading_id(j, pas0);
    while (i < st.size()) {
      std::string key = upper(st[i++].s);
      if (key == "ID") grab(j, pas0);
      else if (key == "X") grab(x0, px);
      else if (key == "Y") grab(y0, py);
      else if (key == "Z") grab(z0, pz);
      else warn("ST1: GENER NOEUD — mot '" + key + "' ignore");
    }
    if (pas0 == 0) pas0 = 1;
    for (int k = 0; k < n; ++k) {
      Id id = static_cast<Id>(j + k * pas0);
      Vec3 c(x0 + k * px, y0 + k * py, z0 + k * pz);
      if (!nodes_.count(id)) node_order_.push_back(id);
      nodes_[id] = c;
    }
  } else if (what == "BARRE") {
    double j = 1, pas0 = 1, n1 = 0, p1 = 0, n2 = 0, p2 = 0;
    leading_id(j, pas0);
    while (i < st.size()) {
      std::string key = upper(st[i++].s);
      if (key == "ID") grab(j, pas0);
      else if (key == "DE") grab(n1, p1);
      else if (key == "A") grab(n2, p2);
      else warn("ST1: GENER BARRE — mot '" + key + "' ignore");
    }
    if (pas0 == 0) pas0 = 1;
    for (int k = 0; k < n; ++k) {
      Id id = static_cast<Id>(j + k * pas0);
      RawBar b{id, static_cast<Id>(n1 + k * p1), static_cast<Id>(n2 + k * p2)};
      if (!bars_.count(id)) bar_order_.push_back(id);
      bars_[id] = b;
    }
  } else {
    warn("ST1: GENER '" + what + "' non reconnu");
  }
  invalidate();
}

void St1Runtime::cmd_barre(const Statement& st) {
  if (st.size() == 1) {
    mode_ = Mode::Bars;
    return;
  }
  // 'BARRE liste DE j1 A j2' is unambiguous; 'BARRE liste j1 j2' peels the
  // trailing two integers as the node pair (the id list is greedy otherwise).
  bool has_de = false;
  for (std::size_t k = 1; k < st.size(); ++k)
    if (!st[k].str && upper(st[k].s) == "DE") has_de = true;
  Id n1 = 0, n2 = 0;
  std::vector<Id> ids;
  if (has_de) {
    std::size_t i = 1;
    ids = id_list(st, i, all_bar_ids());
    require(!ids.empty(), "ST1: BARRE sans identifiant");
    require(i + 3 < st.size() && upper(st[i].s) == "DE" && upper(st[i + 2].s) == "A",
            "ST1: BARRE DE j1 A j2");
    n1 = static_cast<Id>(num(st[i + 1]));
    n2 = static_cast<Id>(num(st[i + 3]));
  } else {
    std::vector<std::size_t> numpos;
    for (std::size_t k = 1; k < st.size(); ++k)
      if (is_num_tok(st[k])) numpos.push_back(k);
    require(numpos.size() >= 3, "ST1: BARRE — noeuds manquants");
    std::size_t split = numpos[numpos.size() - 2];
    Statement listpart(st.begin(), st.begin() + static_cast<long>(split));
    std::size_t i = 1;
    ids = id_list(listpart, i, all_bar_ids());
    require(!ids.empty(), "ST1: BARRE sans identifiant");
    n1 = static_cast<Id>(num(st[numpos[numpos.size() - 2]]));
    n2 = static_cast<Id>(num(st[numpos[numpos.size() - 1]]));
  }
  for (Id b : ids) {
    if (!bars_.count(b)) bar_order_.push_back(b);
    bars_[b] = {b, n1, n2};
  }
  invalidate();
}

// ---------------------------------------------------------------- support --
void St1Runtime::cmd_appui(const Statement& st) {
  if (st.size() == 1) {  // block instruction (manual p.47, p.49)
    mode_ = Mode::Appui;
    return;
  }
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_node_ids());
  require(!ids.empty(), "ST1: APPUI sans liste");
  // Optional NOEUD liste2 (support ids != node ids)
  std::vector<Id> node_ids = ids;
  if (i < st.size() && upper(st[i].s) == "NOEUD") {
    ++i;
    node_ids = id_list(st, i, all_node_ids());
    require(node_ids.size() == ids.size(),
            "ST1: APPUI — listes appuis/noeuds de longueurs differentes");
  }
  int decol = 0, decol_axis = -1;
  double rot = 0;
  std::array<bool, 6> fix{};
  std::vector<double> el;
  bool el_di = false, has_el = false;
  while (i < st.size()) {
    std::string d = st[i].str ? "" : upper(st[i].s);
    if (d == "DX") fix[0] = true;
    else if (d == "DY") fix[1] = true;
    else if (d == "DZ") fix[2] = true;
    else if (d == "RX") fix[3] = true;
    else if (d == "RY") fix[4] = true;
    else if (d == "RZ") fix[5] = true;
    else if (d == "DECOL") {
      decol = +1;
      // DECOL <NON, <X,Y,Z> <POS,NEG>>
      if (i + 1 < st.size()) {
        std::string a = upper(st[i + 1].s);
        if (a == "NON") {
          decol = kDecolOff;
          ++i;
        } else if (a == "X" || a == "Y" || a == "Z") {
          decol_axis = a == "X" ? 0 : (a == "Y" ? 1 : 2);
          ++i;
          if (i + 1 < st.size()) {
            std::string sgn = upper(st[i + 1].s);
            if (sgn == "POS") {
              decol = +1;
              ++i;
            } else if (sgn == "NEG") {
              decol = -1;
              ++i;
            }
          }
        }
      }
    } else if (d == "POS" && decol != 0 && decol != kDecolOff) decol = +1;
    else if (d == "NEG" && decol != 0 && decol != kDecolOff) decol = -1;
    else if (d == "ROT" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      rot = num(st[i + 1]);
      ++i;
    } else if (d == "EL") {
      has_el = true;
      if (i + 1 < st.size() && upper(st[i + 1].s) == "DI") {
        el_di = true;
        ++i;
      }
      ++i;
      while (i < st.size() && (is_num_tok(st[i]) || st[i].s == ",")) {
        if (st[i].s != ",") el.push_back(num(st[i]));
        ++i;
      }
      continue;  // EL consumes to end of numeric run
    } else {
      warn("ST1: APPUI — mot '" + st[i].s + "' ignore");
    }
    ++i;
  }
  // `EL` given over several statements: the real decks lay the 21-term (or
  // 6-term) lower triangle out row by row. Hold the description until the
  // count is complete; `appui_stmt` absorbs the all-numeric continuations.
  if (has_el) {
    bool plane = option_ != 2;
    std::size_t need = el_di ? (plane ? 3u : 6u) : (plane ? 6u : 21u);
    if (el.size() < need) {
      appui_el_open_ = true;
      appui_el_ids_ = ids;
      appui_el_nodes_ = node_ids;
      appui_el_fix_ = fix;
      appui_el_decol_ = decol;
      appui_el_axis_ = decol_axis;
      appui_el_rot_ = rot;
      appui_el_vals_ = el;
      appui_el_di_ = el_di;
      appui_el_need_ = need;
      mode_ = Mode::Appui;
      return;
    }
  }
  apply_appui(ids, node_ids, fix, decol, decol_axis, rot, el, el_di, has_el);
}

void St1Runtime::apply_appui(const std::vector<Id>& ids, const std::vector<Id>& node_ids,
                             const std::array<bool, 6>& fix, int decol, int decol_axis,
                             double rot, const std::vector<double>& el, bool el_di,
                             bool has_el) {
  for (std::size_t k = 0; k < ids.size(); ++k) {
    Id a = ids[k];
    if (!appui_.count(a)) appui_order_.push_back(a);
    SupportDef& s = appui_[a];
    s.node = node_ids[k];
    // A redescribed APPUI (manual p.42): the last description holds, the
    // values it does not redefine keep their previous ones — a new set of
    // blocked DOFs replaces the old one, a new EL the old springs.
    if (std::any_of(fix.begin(), fix.end(), [](bool f) { return f; })) s.fix = fix;
    s.rot = rot;
    if (decol == kDecolOff) {  // DECOL NON
      s.decol = 0;
    } else if (decol != 0) {
      s.decol = decol;
      s.decol_axis = decol_axis;
    }
    if (has_el) {
      // `EL DI` on the structure frame stays on the cheap per-DOF diagonal
      // path (bit-identical); everything else — a full lower triangle, or any
      // matrix given in a ROT alpha frame — becomes a `Support::spring6`
      // built by springmass.h at model-build time (manual p.49-51).
      if (el_di && rot == 0) {
        std::array<double, 6> ks{};
        if (option_ == 0) {  // PLANE order dx, dy, rz
          if (el.size() > 0) ks[0] = el[0];
          if (el.size() > 1) ks[1] = el[1];
          if (el.size() > 2) ks[5] = el[2];
        } else if (option_ == 1) {  // GRILL order rx, dz, ry (manual p.51)
          if (el.size() > 0) ks[3] = el[0];
          if (el.size() > 1) ks[2] = el[1];
          if (el.size() > 2) ks[4] = el[2];
        } else {
          for (std::size_t c = 0; c < 6 && c < el.size(); ++c) ks[c] = el[c];
        }
        s.spring = ks;
        s.spring_lt.clear();
      } else {
        s.spring = {};
        s.spring_lt = el;
        s.spring_di = el_di;
      }
    }
  }
  invalidate();
}

void St1Runtime::cmd_exc(const Statement& st) {
  if (st.size() == 1) {  // block instruction: per-bar lines follow (exemple1/3)
    mode_ = Mode::Exc;
    return;
  }
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_bar_ids());
  std::array<Vec3, 2> e{Vec3::Zero(), Vec3::Zero()};
  while (i < st.size()) {
    std::string end = upper(st[i++].s);
    int slot = end == "OR" ? 0 : (end == "EX" ? 1 : -1);
    if (slot < 0) {
      warn("ST1: EXC — mot '" + end + "' ignore");
      continue;
    }
    // `(X) x (Y) y (Z) z`: a label sends the next value to its component;
    // unlabelled values fill x, y, z in order from the current one
    // (`EX Y -1` is a pure y offset, `EX 0 -1` the same, positionally).
    Vec3 vv = Vec3::Zero();
    int comp = 0;
    while (i < st.size()) {
      std::string k = st[i].str ? "" : upper(st[i].s);
      if (k == "X" || k == "Y" || k == "Z") {
        comp = k[0] - 'X';
        ++i;
        continue;
      }
      if (!is_num_tok(st[i])) break;
      const double x = num(st[i++]);
      if (comp < 3) vv(comp++) = x;
      else warn("ST1: EXC — composante en trop ignoree");
    }
    e[static_cast<std::size_t>(slot)] = vv;
  }
  for (Id b : ids) exc_[b] = e;
  invalidate();
}

void St1Runtime::cmd_art(const Statement& st) {
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_bar_ids());
  Releases rel{};
  int base = 0;
  while (i < st.size()) {
    std::string t = upper(st[i++].s);
    if (t == "OR") base = 0;
    else if (t == "EX") base = 6;
    else if (t == "RX") rel[static_cast<std::size_t>(base + 3)] = true;
    else if (t == "RY") rel[static_cast<std::size_t>(base + 4)] = true;
    else if (t == "RZ") rel[static_cast<std::size_t>(base + 5)] = true;
    else warn("ST1: ART — mot '" + t + "' ignore");
  }
  // Redeclaring ART without a DOF re-blocks it (manual p.55): last wins.
  for (Id b : ids) art_[b] = rel;
  invalidate();
}

void St1Runtime::cmd_beta(const Statement& st) {
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_bar_ids());
  require(i < st.size() && is_num_tok(st[i]), "ST1: BETA attend un angle");
  double b = num(st[i]);
  for (Id id : ids) beta_[id] = b;
  invalidate();
}

void St1Runtime::cmd_cara(const Statement& st) {
  if (st.size() == 1) {  // block instruction (manual p.59)
    cara_pending_bars_.clear();
    mode_ = Mode::Cara;
    return;
  }
  // CARA VAR <LIN Y|LIN Z|PARA> / CARA PSE / CARA SPECIAL open their own
  // FIN-less block modes (manual p.59-60 for the simple AND block forms,
  // p.210 for SPECIAL).
  if (st.size() >= 2 && !st[1].str) {
    std::string w1 = upper(st[1].s);
    if (w1 == "VAR") {
      std::size_t i = 2;
      require(i < st.size(), "ST1: CARA VAR <LIN Y, LIN Z, PARA>");
      std::string k = upper(st[i++].s);
      if (k == "LIN") {
        require(i < st.size(), "ST1: CARA VAR LIN <Y, Z>");
        std::string ax = upper(st[i++].s);
        require(ax == "Y" || ax == "Z", "ST1: CARA VAR LIN attend Y ou Z (manuel p.59)");
        cara_block_law_ = ax == "Y" ? TaperLaw::LinY : TaperLaw::LinZ;
      } else if (k == "PARA") {
        cara_block_law_ = TaperLaw::Para;
      } else {
        fail("ST1: CARA VAR '" + k + "' — attendu LIN Y, LIN Z ou PARA (manuel p.59)");
      }
      cara_block_bars_ = id_list(st, i, all_bar_ids());
      for (Id b : cara_block_bars_) {
        CaraVarDef d;
        d.law = cara_block_law_;
        cara_var_[b] = d;
      }
      mode_ = Mode::CaraVar;
      invalidate();
      return;
    }
    if (w1 == "PSE") {
      if (st.size() == 2) {  // block instruction: the bar lists follow
        cara_block_bars_.clear();
        mode_ = Mode::CaraPse;
        return;
      }
      mode_ = Mode::CaraPse;
      Statement rest(st.begin() + 2, st.end());
      cara_pse_stmt(rest);
      return;
    }
    if (w1 == "SPECIAL") {
      require(st.size() >= 3 && is_id_tok(st[2]), "ST1: CARA SPECIAL j (manuel p.210)");
      cara_special_bar_ = static_cast<Id>(num(st[2]));
      cara_special_[cara_special_bar_] = CaraSpecialDef{};
      cara_special_section_ = 0;
      cara_special_nums_.clear();
      mode_ = Mode::CaraSpecial;
      invalidate();
      return;
    }
  }
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_bar_ids());
  // `CARA liste` with NOTHING after the list: an ST1 command runs until the
  // next keyword, so the characteristics arrive on the FOLLOWING line(s)
  // (OA108 writes `cara 8001 a 8343` / `sx 1e-5 ix 1e-5 ...`). Open the
  // block with the list pending instead of storing empty properties.
  if (i >= st.size() && !ids.empty()) {
    cara_pending_bars_ = ids;
    mode_ = Mode::Cara;
    return;
  }
  CaraProps base;
  std::array<bool, 13> set{};
  while (i + 1 < st.size()) {
    std::string key = upper(st[i].s);
    if (!is_num_tok(st[i + 1])) {
      warn("ST1: CARA — mot '" + st[i].s + "' ignore");
      ++i;
      continue;
    }
    if (!cara_keyword(key, num(st[i + 1]), base, set))
      warn("ST1: CARA — cle '" + key + "' ignoree");
    i += 2;
  }
  for (Id b : ids) {
    CaraProps p = cara_.count(b) ? cara_[b] : CaraProps{};
    cara_merge(p, base, set);
    cara_[b] = p;
  }
  invalidate();
}

// Content lines of a bare `CARA` block (manual p.59, "INSTRUCTION DE BLOC"):
//   liste SX sx (SY sy) … (EXTZ extz)
// A line carrying ONLY keywords continues the previous line's bar list, which
// is how the real decks wrap a long characteristic list over two lines.
void St1Runtime::cara_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  std::size_t i = 0;
  if (is_uid_tok(st[0]) || head == "TOUT" || head == "TOUTES") {
    cara_pending_bars_ = id_list(st, i, all_bar_ids());
    require(!cara_pending_bars_.empty(), "ST1: CARA — liste de barres vide");
  } else {
    CaraProps probe;
    std::array<bool, 13> probe_set{};
    bool is_kw = st.size() >= 2 && is_num_tok(st[1]) &&
                 cara_keyword(head, num(st[1]), probe, probe_set);
    if (!is_kw || cara_pending_bars_.empty()) {
      mode_ = Mode::Top;
      execute(st);
      return;
    }
  }
  CaraProps base;
  std::array<bool, 13> set{};
  while (i + 1 < st.size()) {
    if (!is_num_tok(st[i + 1]) || !cara_keyword(upper(st[i].s), num(st[i + 1]), base, set)) {
      warn("ST1: CARA — mot '" + st[i].s + "' ignore");
      ++i;
      continue;
    }
    i += 2;
  }
  for (Id b : cara_pending_bars_) {
    CaraProps p = cara_.count(b) ? cara_.at(b) : CaraProps{};
    cara_merge(p, base, set);
    cara_[b] = p;
  }
  invalidate();
}

// Content lines of a bare `CONS` block: `liste E e (NU nu) (RO ro) …`.
void St1Runtime::cons_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  if (!is_uid_tok(st[0]) && head != "TOUT" && head != "TOUTES") {
    mode_ = Mode::Top;
    execute(st);
    return;
  }
  Statement full;
  full.push_back({"CONS", false});
  for (const Tok& t : st) full.push_back(t);
  cmd_cons(full);
}

// Content lines of a bare `APPUI` block, plus the all-numeric continuations of
// an `EL` half-matrix spread over several statements (manual p.49-50).
void St1Runtime::appui_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (appui_el_open_) {
    bool all_num = true;
    for (const Tok& t : st) all_num = all_num && (is_num_tok(t) || t.s == ",");
    if (all_num) {
      for (const Tok& t : st)
        if (is_num_tok(t)) appui_el_vals_.push_back(num(t));
      if (appui_el_vals_.size() >= appui_el_need_) {
        apply_appui(appui_el_ids_, appui_el_nodes_, appui_el_fix_, appui_el_decol_,
                    appui_el_axis_, appui_el_rot_, appui_el_vals_, appui_el_di_, true);
        appui_el_open_ = false;
        invalidate();
      }
      return;
    }
    fail("ST1: APPUI EL — demi-matrice incomplete (" +
         std::to_string(appui_el_vals_.size()) + " termes sur " +
         std::to_string(appui_el_need_) + ", manuel p.49)");
  }
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  if (!is_uid_tok(st[0]) && head != "TOUT" && head != "TOUTES") {
    mode_ = Mode::Top;
    execute(st);
    return;
  }
  Statement full;
  full.push_back({"APPUI", false});
  for (const Tok& t : st) full.push_back(t);
  cmd_appui(full);
}

// MINIMAL FALLBACK material reader. `MAT` / `MATERIAU` is a DOMAIN family
// (st1_domains.cpp gets first refusal on every statement, see `domain_handler`)
// — this reader only runs when no domain extension is installed, and only
// understands the elementary constants a `CONS liste MAT i` needs (E, NU, RO,
// TEMP, G_DYN, AMORTISSEMENT). Every other line of the block is skipped with a
// warning, so a deck whose material is a plain set of constants still solves
// headless, and nothing is ever silently invented.
void St1Runtime::mat_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  ConsProps& p = mat_props_[open_mat_];
  if (st.size() >= 2 && is_num_tok(st[1])) {
    double v = num(st[1]);
    if (head == "E") { p.E = v; return; }
    if (head == "NU") { p.nu = v; return; }
    if (head == "RO") { p.ro = v; return; }
    if (head == "TEMP") { p.alpha = v; return; }
    if (head == "G_DYN") { p.g_dyn = v; return; }
    if (head == "AMORTISSEMENT") { p.amortissement = static_cast<int>(v); return; }
  }
  warn("ST1: MAT " + std::to_string(open_mat_) + " — ligne '" + st[0].s +
       "' hors du lecteur de repli (extension domaine requise)");
}

// -------------------------------------------------------------- CARA VAR --
// OR / MI / EX control lines of a `CARA VAR` block (manual p.59, p.61-64).
// A leading identifier list retargets the block (block instruction form); any
// other head closes the block and is re-dispatched.
void St1Runtime::cara_var_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  if (head == "OR" || head == "MI" || head == "EX") {
    require(!cara_block_bars_.empty(), "ST1: CARA VAR — ligne OR/MI/EX sans liste de barres");
    CaraProps base;
    std::array<bool, 13> set{};
    std::size_t i = 1;
    while (i + 1 < st.size()) {
      if (!is_num_tok(st[i + 1]) || !cara_keyword(upper(st[i].s), num(st[i + 1]), base, set)) {
        warn("ST1: CARA VAR — mot '" + st[i].s + "' ignore");
        ++i;
        continue;
      }
      i += 2;
    }
    for (Id b : cara_block_bars_) {
      CaraVarDef& d = cara_var_[b];
      if (head == "OR") {
        cara_merge(d.orig, base, set);
        d.has_or = true;
        // The OR station also serves as the bar's nominal CARA (fibres for the
        // stress output, section used by mass / Kg — see model.h TaperedBar).
        CaraProps p = cara_.count(b) ? cara_.at(b) : CaraProps{};
        cara_merge(p, base, set);
        cara_[b] = p;
      } else if (head == "MI") {
        cara_merge(d.mid, base, set);
        d.has_mi = true;
      } else {
        cara_merge(d.ext, base, set);
        d.has_ex = true;
      }
    }
    invalidate();
    return;
  }
  if (is_id_tok(st[0]) || head == "TOUT" || head == "TOUTES") {
    std::size_t i = 0;
    cara_block_bars_ = id_list(st, i, all_bar_ids());
    for (Id b : cara_block_bars_) {
      CaraVarDef d;
      d.law = cara_block_law_;
      cara_var_[b] = d;
    }
    invalidate();
    return;
  }
  mode_ = Mode::Top;
  execute(st);
}

// -------------------------------------------------------------- CARA PSE --
// ZONE / PLASTIQUE lines of a `CARA PSE` block (manual p.59-61, soilzones.h).
void St1Runtime::cara_pse_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    return;
  }
  auto moduli = [&](std::size_t i, double& kfx, double& kfy, double& kfz, double& kmx,
                    double& xl, bool& has_xl) {
    while (i < st.size()) {
      std::string k = upper(st[i].s);
      if (i + 1 >= st.size() || !is_num_tok(st[i + 1])) {
        warn("ST1: CARA PSE — mot '" + st[i].s + "' ignore");
        ++i;
        continue;
      }
      double v = num(st[i + 1]);
      if (k == "KFX") kfx = v;
      else if (k == "KFY") kfy = v;
      else if (k == "KFZ") kfz = v;
      else if (k == "KMX") kmx = v;
      else if (k == "XL") { xl = v; has_xl = true; }
      else warn("ST1: CARA PSE — cle '" + k + "' ignoree");
      i += 2;
    }
  };
  if (head == "ZONE") {
    require(!cara_block_bars_.empty(), "ST1: CARA PSE ZONE sans liste de barres");
    std::size_t i = 1;
    if (i < st.size() && is_id_tok(st[i])) ++i;  // the zone number is positional
    SoilZone z;
    double xl = 0;
    bool has_xl = false;
    moduli(i, z.kfx, z.kfy, z.kfz, z.kmx, xl, has_xl);
    z.xl_end = xl;  // ignored on the LAST zone (manual p.61)
    for (Id b : cara_block_bars_) soil_zones_[b].zones.push_back(z);
    invalidate();
    return;
  }
  if (head == "PLASTIQUE") {
    require(!cara_block_bars_.empty(), "ST1: CARA PSE PLASTIQUE sans liste de barres");
    SoilPlastic pl;
    std::size_t i = 1;
    while (i < st.size()) {
      std::string k = upper(st[i].s);
      if (k == "POS") { pl.pos = true; ++i; continue; }
      if (k == "NEG") { pl.neg = true; ++i; continue; }
      break;
    }
    double xl = 0;
    bool has_xl = false;
    moduli(i, pl.ffx, pl.ffy, pl.ffz, pl.fmx, xl, has_xl);
    for (Id b : cara_block_bars_) soil_zones_[b].plastic = pl;
    invalidate();
    return;
  }
  if (is_uid_tok(st[0]) || head == "TOUT" || head == "TOUTES") {
    // A bar list (+ its ordinary characteristics), optionally followed on the
    // SAME statement by the first ZONE / PLASTIQUE line — both shapes appear
    // in the real decks. Re-declaring a bar REPLACES its previous soil
    // description entirely (manual p.42: the last description wins).
    std::size_t i = 0;
    cara_block_bars_ = id_list(st, i, all_bar_ids());
    require(!cara_block_bars_.empty(), "ST1: CARA PSE — liste de barres vide");
    CaraProps base;
    std::array<bool, 13> set{};
    while (i + 1 < st.size() && is_num_tok(st[i + 1])) {
      if (!cara_keyword(upper(st[i].s), num(st[i + 1]), base, set)) break;
      i += 2;
    }
    for (Id b : cara_block_bars_) {
      CaraProps p = cara_.count(b) ? cara_.at(b) : CaraProps{};
      cara_merge(p, base, set);
      cara_[b] = p;
      soil_zones_[b] = BarSoil{};
    }
    invalidate();
    if (i < st.size()) {
      Statement rest(st.begin() + static_cast<long>(i), st.end());
      cara_pse_stmt(rest);
    }
    return;
  }
  mode_ = Mode::Top;
  execute(st);
}

// ---------------------------------------------------------- CARA SPECIAL --
// `RIGIDITE CARA` / `EFFORT` sections of a `CARA SPECIAL j` block, both given
// as bare positional numbers (manual p.210):
//   OPTION PLANE    : sx sy iz          | fx fy mz
//   OPTION SPATIALE : sx sy sz ix iy iz ip | fx fy fz mx my mz
void St1Runtime::cara_special_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  auto flush = [&]() {
    if (cara_special_section_ == 0) return;
    CaraSpecialDef& d = cara_special_[cara_special_bar_];
    const std::vector<double>& v = cara_special_nums_;
    if (cara_special_section_ == 1) {
      if (option_ == 2) {
        require(v.size() >= 6, "ST1: CARA SPECIAL RIGIDITE CARA — 6 ou 7 valeurs en SPATIALE");
        d.cara.sx = v[0];
        d.cara.sy = v[1];
        d.cara.sz = v[2];
        d.cara.ix = v[3];
        d.cara.iy = v[4];
        d.cara.iz = v[5];
        if (v.size() > 6) d.cara.ipol = v[6];
      } else {
        require(v.size() >= 3, "ST1: CARA SPECIAL RIGIDITE CARA — 3 valeurs en PLANE/GRILL");
        d.cara.sx = v[0];
        d.cara.sy = v[1];
        d.cara.iz = v[2];
      }
      d.has_cara = true;
      CaraProps p = d.cara;
      cara_[cara_special_bar_] = p;
    } else {
      if (option_ == 2) {
        require(v.size() >= 6, "ST1: CARA SPECIAL EFFORT — 6 valeurs en SPATIALE");
        for (std::size_t c = 0; c < 6; ++c) d.lim.cap[c] = std::fabs(v[c]);
      } else {
        require(v.size() >= 3, "ST1: CARA SPECIAL EFFORT — 3 valeurs en PLANE/GRILL");
        d.lim.cap[0] = std::fabs(v[0]);
        d.lim.cap[1] = std::fabs(v[1]);
        d.lim.cap[5] = std::fabs(v[2]);
      }
      d.lim.validate();
      d.has_effort = true;
    }
    cara_special_nums_.clear();
    cara_special_section_ = 0;
  };
  if (head == "FIN") {
    flush();
    mode_ = Mode::Top;
    invalidate();
    return;
  }
  if (head == "RIGIDITE") {
    flush();
    cara_special_section_ = 1;
    for (std::size_t i = 1; i < st.size(); ++i)
      if (is_num_tok(st[i])) cara_special_nums_.push_back(num(st[i]));
    return;
  }
  if (head == "EFFORT") {
    flush();
    cara_special_section_ = 2;
    for (std::size_t i = 1; i < st.size(); ++i)
      if (is_num_tok(st[i])) cara_special_nums_.push_back(num(st[i]));
    return;
  }
  bool all_num = true;
  for (const Tok& t : st) all_num = all_num && is_num_tok(t);
  if (all_num && cara_special_section_ != 0) {
    for (const Tok& t : st) cara_special_nums_.push_back(num(t));
    return;
  }
  flush();
  mode_ = Mode::Top;
  execute(st);
}

void St1Runtime::cmd_cons(const Statement& st) {
  if (st.size() == 1) {  // block instruction
    mode_ = Mode::Cons;
    return;
  }
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, all_bar_ids());
  // A redescribed bar (manual p.42) keeps the constants this line does not
  // redefine; `CONS ... MAT i` starts the bar afresh on that material.
  std::vector<std::function<void(ConsProps&)>> sets;
  int mat = 0;
  bool any = false;
  while (i < st.size()) {
    std::string key = upper(st[i].s);
    if (key == "MAT" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      mat = static_cast<int>(num(st[i + 1]));
      i += 2;
      any = true;
      continue;
    }
    if (i + 1 >= st.size() || !is_num_tok(st[i + 1])) {
      warn("ST1: CONS — mot '" + st[i].s + "' ignore");
      ++i;
      continue;
    }
    double v = num(st[i + 1]);
    if (key == "E") sets.push_back([v](ConsProps& c) { c.E = v; });
    else if (key == "NU") sets.push_back([v](ConsProps& c) { c.nu = v; });
    else if (key == "RO") sets.push_back([v](ConsProps& c) { c.ro = v; });
    else if (key == "TEMP") sets.push_back([v](ConsProps& c) { c.alpha = v; });
    else if (key == "G_DYN") sets.push_back([v](ConsProps& c) { c.g_dyn = v; });
    else if (key == "AMORTISSEMENT")
      sets.push_back([v](ConsProps& c) { c.amortissement = static_cast<int>(v); });
    else warn("ST1: CONS — cle '" + key + "' ignoree");
    i += 2;
    any = true;
  }
  require(any, "ST1: CONS sans constantes");
  for (Id b : ids) {
    ConsProps c;
    if (mat != 0) c.mat = mat;
    else if (auto it = cons_.find(b); it != cons_.end()) c = it->second;
    for (const auto& f : sets) f(c);
    cons_[b] = c;
  }
  invalidate();
}

// ------------------------------------------------------------------ CHARG --
void St1Runtime::cmd_charg_open(const Statement& st) {
  ChargDef def;
  std::size_t i = 1;
  if (i < st.size() && is_id_tok(st[i])) {
    def.id = static_cast<Id>(num(st[i]));
    ++i;
  } else {
    def.id = 1;  // manual p.106: "numero du chargement (par defaut i=1)"
  }
  if (i < st.size() && !st[i].str && upper(st[i].s) == "COMB") {
    def.is_comb = true;
    ++i;
  }
  if (i < st.size() && st[i].str) def.name = st[i].s;
  if (def.name.empty()) def.name = "CHARG " + std::to_string(def.id);
  if (!chargs_.count(def.id)) charg_order_.push_back(def.id);
  chargs_[def.id] = def;
  open_charg_ = &chargs_[def.id];
  mode_ = Mode::Charg;
}

void St1Runtime::charg_line(std::vector<Load>& loads, Mode& sub, const Statement& st) {
  ChargDef tmp;
  tmp.loads = std::move(loads);
  ChargDef* const save_open = open_charg_;
  const Mode save_mode = mode_;
  const bool save_sel = selective_activation_;
  auto restore = [&] {
    loads = std::move(tmp.loads);
    open_charg_ = save_open;
    mode_ = save_mode;
    selective_activation_ = save_sel;
  };
  open_charg_ = &tmp;
  selective_activation_ = false;
  mode_ = sub == Mode::Top ? Mode::Charg : sub;
  // a keyword line closes the open sub-list (as the top-level dispatcher)
  const std::string head = st[0].str ? "" : upper(st[0].s);
  if (mode_ != Mode::Charg && !is_uid_tok(st[0]) && head != "TOUT") mode_ = Mode::Charg;
  try {
    charg_stmt(st);
  } catch (...) {
    restore();
    throw;
  }
  if (!tmp.cable_actions.empty())
    warn("ST1: PHASAGE CHARG CABLE — utiliser TENDRE / DETENDRE dans le phasage");
  sub = mode_;
  restore();
}

void St1Runtime::charg_stmt(const Statement& st) {
  ChargDef& c = *open_charg_;
  std::string head = st[0].str ? "" : upper(st[0].s);

  // Sub-mode switches
  if (head == "FIN") {
    open_charg_ = nullptr;
    mode_ = Mode::Top;
    return;
  }
  if (head == "NOEUD" && st.size() == 1) {
    mode_ = Mode::ChargNodes;
    return;
  }
  if (head == "BARRE" && st.size() == 1) {
    mode_ = Mode::ChargBars;
    return;
  }
  // Simple form on one statement: `BARRE liste <UNI,CON,LIN> …` (manual p.106).
  if (head == "BARRE") {
    Statement rest(st.begin() + 1, st.end());
    Mode save = mode_;
    mode_ = Mode::ChargBars;
    charg_stmt(rest);
    mode_ = save == Mode::ChargBars ? Mode::ChargBars : Mode::Charg;
    return;
  }

  if (c.is_comb) {
    // CHARG i COMB body: 'CHARG liste (coef)' lines.
    require(head == "CHARG", "ST1: CHARG COMB attend des lignes 'CHARG liste coef'");
    std::size_t i = 1;
    std::vector<Id> ids;
    for (Id cid : charg_order_) ids.push_back(cid);
    std::vector<Id> list = id_list(st, i, ids);
    double coef = (i < st.size() && is_num_tok(st[i])) ? num(st[i]) : 1.0;
    for (Id cid : list) c.parts.emplace_back(cid, coef);
    return;
  }

  if (head == "POIDS") {
    // POIDS PROPRE (liste) | POIDS PROPRE <X,Y,Z> (POS,NEG) liste
    require(option_ != 1, "ST1: POIDS PROPRE interdit en GRILL (manuel p.104)");
    std::size_t i = 1;
    if (i < st.size() && upper(st[i].s) == "PROPRE") ++i;
    Vec3 dir = Vec3::Zero();
    if (option_ == 0)
      dir.y() = -1;
    else
      dir.z() = -1;
    if (i < st.size()) {
      std::string a = upper(st[i].s);
      double sgn = -1;
      int axis = -1;
      if (a == "X") axis = 0;
      else if (a == "Y") axis = 1;
      else if (a == "Z") axis = 2;
      if (axis >= 0) {
        ++i;
        if (i < st.size() && upper(st[i].s) == "POS") {
          sgn = +1;
          ++i;
        } else if (i < st.size() && upper(st[i].s) == "NEG") {
          sgn = -1;
          ++i;
        }
        dir = Vec3::Zero();
        dir(axis) = sgn;
      }
    }
    std::vector<Id> list = id_list(st, i, all_bar_ids());
    for (Load& l : self_weight_loads(list, dir)) c.loads.push_back(std::move(l));
    return;
  }

  if (head == "NOEUD") {
    std::size_t i = 1;
    std::vector<Id> ids = id_list(st, i, all_node_ids());
    Vec3 F = Vec3::Zero(), M = Vec3::Zero();
    while (i + 1 < st.size()) {
      std::string key = upper(st[i].s);
      if (!is_num_tok(st[i + 1])) break;
      double v = num(st[i + 1]);
      if (key == "FX") F.x() += v;
      else if (key == "FY") F.y() += v;
      else if (key == "FZ") F.z() += v;
      else if (key == "MX") M.x() += v;
      else if (key == "MY") M.y() += v;
      else if (key == "MZ") M.z() += v;
      else break;
      i += 2;
    }
    if (i < st.size())
      warn("ST1: CHARG NOEUD — mot '" + st[i].s + "' ignore (avec la suite de la ligne)");
    for (Id n : ids) c.loads.push_back(Load::nodal(n, F, M));
    return;
  }

  if (head == "APPUI" || mode_ == Mode::ChargAppui) {
    // Imposed support displacements (global frame, rotations in radians);
    // `APPUI` alone opens a block of `liste DX dx …` lines.
    if (head == "APPUI" && st.size() == 1) {
      mode_ = Mode::ChargAppui;
      return;
    }
    std::size_t i = head == "APPUI" ? 1 : 0;
    std::vector<Id> ids = id_list(st, i, all_node_ids());
    std::array<std::optional<double>, 6> d{};
    while (i + 1 < st.size()) {
      std::string key = upper(st[i].s);
      if (!is_num_tok(st[i + 1])) break;
      double v = num(st[i + 1]);
      if (key == "DX") d[0] = v;
      else if (key == "DY") d[1] = v;
      else if (key == "DZ") d[2] = v;
      else if (key == "RX") d[3] = v;
      else if (key == "RY") d[4] = v;
      else if (key == "RZ") d[5] = v;
      else break;
      i += 2;
    }
    if (i < st.size())
      warn("ST1: CHARG APPUI — mot '" + st[i].s + "' ignore (avec la suite de la ligne)");
    for (Id a : ids) {
      Id n = appui_.count(a) ? appui_.at(a).node : a;
      c.loads.push_back(Load::support_disp(n, d));
    }
    return;
  }

  if (head == "TEMP" || mode_ == Mode::ChargTemp) {
    std::size_t i = head == "TEMP" ? 1 : 0;
    if (head == "TEMP" && st.size() == 1) {
      mode_ = Mode::ChargTemp;
      return;
    }
    std::vector<Id> ids = id_list(st, i, all_bar_ids());
    if (ids.empty()) ids = all_bar_ids();
    double dt = 0, dty = 0, dtz = 0;
    while (i + 1 < st.size()) {
      std::string key = upper(st[i].s);
      if (!is_num_tok(st[i + 1])) break;
      double v = num(st[i + 1]);
      if (key == "UNI") dt = v;
      else if (key == "GY") dty = v;
      else if (key == "GZ") dtz = v;
      else break;
      i += 2;
    }
    if (i < st.size())
      warn("ST1: CHARG TEMP — mot '" + st[i].s + "' ignore (avec la suite de la ligne)");
    // Fibre height VY + WY (WY defaults to VY), from the given properties.
    auto height = [](double v, double w) { return v + (w > 0 ? w : v); };
    for (Id b : ids) {
      double gy = 0, gz = 0;
      auto var = cara_var_.find(b);
      if ((dty != 0 || dtz != 0) && var != cara_var_.end() && var->second.has_or &&
          var->second.has_ex) {
        // CARA VAR: the fibres follow the section from OR to EX (linearly),
        // so the curvature alpha dT / h(x) varies along the bar. It is
        // lowered as the LINEAR curvature with the same zeroth and first
        // moments over the bar — exact end rotation and deflection of a
        // statically determinate member; the hyperstatic block sees the
        // projection of 1/h on linear fields.
        const CaraProps& o = var->second.orig;
        const CaraProps& e = var->second.ext;
        require(cons_.count(b) && cons_.at(b).alpha != 0,
                "ST1: TEMP gradient — coefficient TEMP (CONS) requis barre " + std::to_string(b));
        const double al = cons_.at(b).alpha;
        const CaraProps& ce = e.vy > 0 || e.vz > 0 ? e : o;
        // Linear field a + b t (t = x / L) with the moments of 1 / h(t).
        auto inv_h_linear = [](double h0, double h1, double& k0, double& k1) {
          const double d = h1 - h0;
          if (std::fabs(d) < 1e-9 * h0) {
            k0 = k1 = 1 / h0;
            return;
          }
          const double m0 = std::log(h1 / h0) / d, m1 = (1 - h0 * m0) / d;
          const double bb = 12 * (m1 - m0 / 2), aa = m0 - bb / 2;
          k0 = aa;
          k1 = aa + bb;
        };
        double kz[2] = {0, 0}, ky[2] = {0, 0};
        if (dty != 0) {
          double h0 = height(o.vy, o.wy), h1 = height(ce.vy, ce.wy);
          require(h0 > 0 && h1 > 0, "ST1: TEMP GY — VY/WY requis barre " + std::to_string(b));
          inv_h_linear(h0, h1, kz[0], kz[1]);
          kz[0] *= al * dty;
          kz[1] *= al * dty;
        }
        if (dtz != 0) {
          double h0 = height(o.vz, o.wz), h1 = height(ce.vz, ce.wz);
          require(h0 > 0 && h1 > 0, "ST1: TEMP GZ — VZ/WZ requis barre " + std::to_string(b));
          inv_h_linear(h0, h1, ky[0], ky[1]);
          ky[0] *= -al * dtz;
          ky[1] *= -al * dtz;
        }
        if (dt != 0) c.loads.push_back(Load::thermal(b, dt));
        c.loads.push_back(Load::prestrain_linear(b, 0.0, ky[0], kz[0], 0.0, ky[1], kz[1]));
        continue;
      }
      if (dty != 0 || dtz != 0) {
        require(cara_.count(b),
                "ST1: TEMP gradient — fibres CARA manquantes barre " + std::to_string(b));
        const CaraProps& cp = cara_.at(b);
        if (dty != 0) {
          double h = height(cp.vy, cp.wy);
          require(h > 0, "ST1: TEMP GY — VY/WY requis barre " + std::to_string(b));
          gy = dty / h;
        }
        if (dtz != 0) {
          double h = height(cp.vz, cp.wz);
          require(h > 0, "ST1: TEMP GZ — VZ/WZ requis barre " + std::to_string(b));
          gz = dtz / h;
        }
      }
      c.loads.push_back(Load::thermal(b, dt, gy, gz));
    }
    return;
  }

  if (head == "DEFOR" || mode_ == Mode::ChargDefor) {
    std::size_t i = head == "DEFOR" ? 1 : 0;
    if (head == "DEFOR" && st.size() == 1) {
      mode_ = Mode::ChargDefor;
      return;
    }
    std::vector<Id> ids = id_list(st, i, all_bar_ids());
    if (ids.empty()) ids = all_bar_ids();
    require(i < st.size(), "ST1: DEFOR incomplet");
    std::string kind = upper(st[i++].s);
    // Local frame of a bar (GLO components are turned into it).
    auto frame_of = [&](Id b) {
      const RawBar& rb = bars_.at(b);
      return bar_rotation(nodes_.at(rb.n1), nodes_.at(rb.n2), beta_.count(b) ? beta_.at(b) : 0.0);
    };
    // The distributed generalized strains keta's Prestrain does not carry —
    // shear DY DZ (transverse displacement gradients), and every component
    // over a part [x1, x2] of the bar (DEFOR LIN XL) — values ga -> gb:
    // lowered as `kDeforPts` concentrated dislocations at the midpoints of
    // equal sub-segments (a transverse slip or a twist jump sets the same end
    // forces wherever it sits, so the bar forces are exact and the
    // deflection converges with the count). The twist RX over the whole bar
    // is the exact imposed twist rate kappa_x.
    constexpr int kDeforPts = 16;
    auto lumped = [&](Id b, double x1, double x2, const std::array<double, 6>& ga,
                      const std::array<double, 6>& gb, bool all) {
      const double h = (x2 - x1) / kDeforPts;
      for (int k = 0; k < kDeforPts; ++k) {
        const double f = (k + 0.5) / kDeforPts, x = x1 + (k + 0.5) * h;
        std::array<double, 6> g{};
        for (int c = 0; c < 6; ++c) g[c] = (ga[c] + f * (gb[c] - ga[c])) * h;
        // rotation jumps are the opposite of keta's kinks (as DEFOR CON)
        Load l = Load::prestrain_point(b, x, all ? g[0] : 0.0, all ? -g[4] : 0.0,
                                       all ? -g[5] : 0.0);
        l.dvy = g[1];
        l.dvz = g[2];
        l.dthx = g[3];
        if (l.strain != 0 || l.kappa_y != 0 || l.kappa_z != 0 || l.dvy != 0 || l.dvz != 0 ||
            l.dthx != 0)
          c.loads.push_back(l);
      }
    };
    static const std::string comps[] = {"DX", "DY", "DZ", "RX", "RY", "RZ"};
    // GLO: the displacement and rotation gradient vectors into the bar axes
    auto to_local = [&](Id b, std::array<double, 6> g) {
      const Mat3 R = frame_of(b);
      const Vec3 d = R * Vec3(g[0], g[1], g[2]), r = R * Vec3(g[3], g[4], g[5]);
      return std::array<double, 6>{d.x(), d.y(), d.z(), r.x(), r.y(), r.z()};
    };
    if (kind == "UNI") {
      // `DEFOR liste UNI (DX dx) (DY dy) ... (RZ rz) (GLO,LOC)` (p.120):
      // DX the axial strain, RY RZ the rotation gradients (the opposite of
      // keta's curvatures), RX the twist rate (keta's kappa_x, same sign),
      // DY DZ shear (lowered as dislocations)
      std::array<double, 6> g{};
      bool glo = false;
      if (i < st.size() && is_num_tok(st[i])) {
        g[0] = num(st[i]);
      } else {
        while (i < st.size()) {
          const std::string key = st[i].str ? "" : upper(st[i].s);
          if (key == "GLO") glo = true;
          else if (key == "LOC") glo = false;
          else if (const auto* k = std::find(std::begin(comps), std::end(comps), key);
                   k != std::end(comps) && i + 1 < st.size() && is_num_tok(st[i + 1])) {
            g[static_cast<std::size_t>(k - std::begin(comps))] = num(st[i + 1]);
            ++i;
          } else {
            warn("ST1: DEFOR UNI — mot '" + st[i].s + "' ignore");
          }
          ++i;
        }
      }
      for (Id b : ids) {
        const std::array<double, 6> gl = glo ? to_local(b, g) : g;
        Load l;
        l.kind = Load::Kind::Prestrain;
        l.bar = b;
        l.strain = gl[0];
        l.kappa_y = -gl[4];  // ST1: imposed rotation gradient,
        l.kappa_z = -gl[5];  // the opposite of keta's kappa
        l.kappa_x = gl[3];   // d(theta_x)/dx, the same measure
        if (l.strain != 0 || l.kappa_y != 0 || l.kappa_z != 0 || l.kappa_x != 0) c.loads.push_back(l);
        if (gl[1] != 0 || gl[2] != 0) {
          std::array<double, 6> shear = gl;
          shear[3] = 0.0;
          lumped(b, 0.0, bar_length(b), shear, shear, false);
        }
      }
    } else if (kind == "CON") {
      // `DEFOR liste CON (XL x0) (REL/ABS) DX dux RY thy RZ thz`: concentrated
      // dislocation at abscissa x0 — DX an axial gap [m], RY/RZ rotation kinks
      // [rad]. Documented hypothesis: the manual's DEFOR chapter gives no
      // worked CON/LIN example, so the syntax mirrors the other BARRE loads
      // (XL/REL/ABS abscissa markers, then component/value pairs).
      // (manual p.119: DX DY DZ jumps of the displacement, RX RY RZ of the
      // rotation, in the bar axes unless GLO; the markers may come anywhere)
      double x0 = 0;
      bool relative = false, glo = false;
      double d[6] = {0, 0, 0, 0, 0, 0};
      while (i < st.size()) {
        std::string t = st[i].str ? "" : upper(st[i].s);
        if (t == "REL") relative = true;
        else if (t == "ABS") relative = false;
        else if (t == "GLO") glo = true;
        else if (t == "LOC") glo = false;
        else if ((t == "XL" || t == "X") && i + 1 < st.size() && is_num_tok(st[i + 1])) {
          x0 = num(st[i + 1]);
          ++i;
        } else if (const auto* k = std::find(std::begin(comps), std::end(comps), t);
                   k != std::end(comps) && i + 1 < st.size() && is_num_tok(st[i + 1])) {
          d[k - std::begin(comps)] = num(st[i + 1]);
          ++i;
        } else {
          warn("ST1: DEFOR CON — mot '" + st[i].s + "' ignore");
        }
        ++i;
      }
      for (Id b : ids) {
        // the rotation jumps are the opposite of keta's kinks (ST1's
        // RZ 1 is the influence line of Mz)
        Load l = Load::prestrain_point(b, x0, d[0], -d[4], -d[5]);
        l.dvy = d[1];
        l.dvz = d[2];
        l.dthx = d[3];
        l.relative = relative;
        if (glo) l.frame = LoadFrame::Global;
        c.loads.push_back(l);
      }
    } else if (kind == "LIN") {
      // `DEFOR liste LIN DX e_i e_j RY ky_i ky_j RZ kz_i kz_j`: strain /
      // curvature fields varying linearly from end i to end j (two values per
      // component, mirroring `BARRE ... LIN FY va vb`; hypothesis documented
      // as for CON above).
      // (p.121: DX DY DZ RX RY RZ, two values each, GLO/LOC; an `XL x1 x2
      // (REL,ABS)` range lowers every component as dislocations over it)
      std::array<double, 6> ga{}, gb{};
      double r1 = 0, r2 = -1;
      bool rrel = false, glo = false;
      while (i < st.size()) {
        const std::string t = st[i].str ? "" : upper(st[i].s);
        if (t == "REL") rrel = true;
        else if (t == "ABS") rrel = false;
        else if (t == "GLO") glo = true;
        else if (t == "LOC") glo = false;
        else if ((t == "XL" || t == "X") && i + 2 < st.size() && is_num_tok(st[i + 1]) &&
                 is_num_tok(st[i + 2])) {
          r1 = num(st[i + 1]);
          r2 = num(st[i + 2]);
          i += 2;
        } else if (const auto* k = std::find(std::begin(comps), std::end(comps), t);
                   k != std::end(comps) && i + 2 < st.size() && is_num_tok(st[i + 1]) &&
                   is_num_tok(st[i + 2])) {
          const std::size_t c6 = static_cast<std::size_t>(k - std::begin(comps));
          ga[c6] = num(st[i + 1]);
          gb[c6] = num(st[i + 2]);
          i += 2;
        } else {
          warn("ST1: DEFOR LIN — mot '" + st[i].s + "' ignore");
        }
        ++i;
      }
      for (Id b : ids) {
        const double L = bar_length(b), x1 = rrel ? r1 * L : r1,
                     x2 = r2 < 0 ? L : (rrel ? r2 * L : r2);
        const std::array<double, 6> la = glo ? to_local(b, ga) : ga,
                                    lb = glo ? to_local(b, gb) : gb;
        const bool whole = x1 <= 1e-9 * (1 + L) && std::fabs(x2 - L) <= 1e-6 * (1 + L);
        if (whole) {
          if (la[0] != 0 || lb[0] != 0 || la[3] != 0 || lb[3] != 0 || la[4] != 0 || lb[4] != 0 || la[5] != 0 ||
              lb[5] != 0) {
            Load l = Load::prestrain_linear(b, la[0], -la[4], -la[5], lb[0], -lb[4], -lb[5]);
            l.kappa_x = la[3];
            l.kappa_x2 = lb[3];
            c.loads.push_back(l);
          }
          if (la[1] != 0 || lb[1] != 0 || la[2] != 0 || lb[2] != 0) {
            std::array<double, 6> sa = la, sb = lb;
            sa[3] = sb[3] = 0.0;
            lumped(b, 0.0, L, sa, sb, false);
          }
        } else {
          lumped(b, x1, x2, la, lb, true);
        }
      }
    } else {
      warn("ST1: DEFOR " + kind + " — genre non reconnu (UNI/CON/LIN tenus)");
    }
    return;
  }

  if (head == "CABLE") {
    ChargDef::CableAction a;
    std::size_t i = 1;
    a.cables = id_list(st, i, {});
    while (i < st.size()) {
      std::string k = upper(st[i].s);
      if (k == "TENSION") a.tension = true;
      else if (k == "PERTE") a.perte = true;
      else if (k == "INST") a.perte_inst = true;
      else if (k == "RETRAIT" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
        a.retrait = num(st[++i]);
      } else if (k == "UNI" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
        a.uni = num(st[++i]);
      } else if (k == "RELAX") {
        a.relax = true;
        if (i + 1 < st.size() && is_num_tok(st[i + 1])) a.relax_cf = num(st[++i]);
      } else {
        warn("ST1: CHARG CABLE — mot '" + st[i].s + "' ignore");
      }
      ++i;
    }
    if (!a.tension && !a.perte) a.tension = true;
    c.cable_actions.push_back(a);
    // Lowered into equivalent loads by the prestress domain at EXEC CHARG.
    return;
  }

  // ChargBars / ChargNodes content lines
  if (mode_ == Mode::ChargBars) {
    std::size_t i = 0;
    std::vector<Id> ids = id_list(st, i, all_bar_ids());
    require(!ids.empty() && i < st.size(), "ST1: ligne de charge de barre invalide");
    std::string bad;
    if (!st1_bar_load_line(st, i, ids, c.loads, &bad))
      warn("ST1: charge de barre '" + st[i].s + "' non reconnue");
    else if (!bad.empty())
      warn("ST1: charge de barre — mot '" + bad + "' ignore (avec la suite de la ligne)");
    return;
  }
  if (mode_ == Mode::ChargNodes) {
    std::size_t i = 0;
    std::vector<Id> ids = id_list(st, i, all_node_ids());
    Vec3 F = Vec3::Zero(), M = Vec3::Zero();
    while (i + 1 < st.size()) {
      std::string key = upper(st[i].s);
      if (!is_num_tok(st[i + 1])) break;
      double v = num(st[i + 1]);
      if (key == "FX") F.x() += v;
      else if (key == "FY") F.y() += v;
      else if (key == "FZ") F.z() += v;
      else if (key == "MX") M.x() += v;
      else if (key == "MY") M.y() += v;
      else if (key == "MZ") M.z() += v;
      else break;
      i += 2;
    }
    if (i < st.size())
      warn("ST1: CHARG NOEUD — mot '" + st[i].s + "' ignore (avec la suite de la ligne)");
    for (Id n : ids) c.loads.push_back(Load::nodal(n, F, M));
    return;
  }
  warn("ST1: instruction de chargement '" + st[0].s + "' non reconnue");
}

bool st1_bar_load_line(const Statement& st, std::size_t i, const std::vector<Id>& ids,
                       std::vector<Load>& out, std::string* ignored) {
  require(i < st.size(), "ST1: ligne de charge de barre invalide");
  const std::string kind = upper(st[i++].s);
  if (kind != "UNI" && kind != "CON" && kind != "LIN") return false;
  const int nv = kind == "LIN" ? 2 : 1;
  LoadFrame frame = LoadFrame::Global;
  bool relative = false, has_x = false;
  double x1 = 0, x2 = -1;
  Vec3 Fa = Vec3::Zero(), Fb = Vec3::Zero(), Ma = Vec3::Zero(), Mb = Vec3::Zero();
  static const std::string comps[] = {"FX", "FY", "FZ", "MX", "MY", "MZ"};
  while (i < st.size()) {
    const std::string t = st[i].str ? "" : upper(st[i].s);
    if (t == "GLO") frame = LoadFrame::Global;
    else if (t == "LOC") frame = LoadFrame::Local;
    else if (t == "PRO") frame = LoadFrame::GlobalProjected;
    else if (t == "REL") relative = true;
    else if (t == "ABS") relative = false;
    else if ((t == "XL" || t == "X") && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      x1 = num(st[++i]);
      if (i + 1 < st.size() && is_num_tok(st[i + 1])) x2 = num(st[++i]);
      has_x = true;
    } else {
      const auto* c = std::find(std::begin(comps), std::end(comps), t);
      if (c == std::end(comps) || i + nv >= st.size() + 0 || !is_num_tok(st[i + 1]) ||
          (nv == 2 && !is_num_tok(st[i + 2]))) {
        if (ignored) *ignored = st[i].s;
        break;
      }
      const int k = static_cast<int>(c - std::begin(comps));
      const double va = num(st[i + 1]), vb = nv == 2 ? num(st[i + 2]) : va;
      (k < 3 ? Fa : Ma)(k % 3) += va;
      (k < 3 ? Fb : Mb)(k % 3) += vb;
      i += nv;
    }
    ++i;
  }
  for (Id b : ids) {
    Load l;
    if (kind == "UNI") {
      l = Load::udl(b, Fa, frame);
      l.M = Ma;  // distributed couple (no shear)
    } else if (kind == "CON") {
      l = Load::point(b, has_x ? x1 : 0.0, Fa, Ma, frame);
    } else {
      l = Load::linear(b, Fa, Fb, frame, has_x ? x1 : 0.0, has_x ? x2 : -1.0);
      l.M = Ma;
      l.M2 = Mb;
    }
    l.relative = relative;
    out.push_back(l);
  }
  return true;
}

// ------------------------------------------------------------------- COMB --
void St1Runtime::cmd_comb(const Statement& st) {
  require(st.size() >= 2 && is_id_tok(st[1]), "ST1: COMB invalide");
  CombDef cb;
  cb.id = static_cast<Id>(num(st[1]));
  std::size_t i = 2;
  if (i < st.size() && st[i].str) cb.name = st[i++].s;
  if (cb.name.empty()) cb.name = "COMB " + std::to_string(cb.id);
  // Inline compat form: c1 CHARG i1 c2 CHARG i2 ...
  bool inline_terms = false;
  while (i + 2 < st.size() + 0 && i < st.size() && is_num_tok(st[i])) {
    double coef = num(st[i]);
    if (i + 2 >= st.size() || upper(st[i + 1].s) != "CHARG") break;
    cb.terms.push_back({"CHARG", static_cast<Id>(num(st[i + 2])), coef});
    i += 3;
    inline_terms = true;
  }
  if (inline_terms) {
    combs_[cb.id] = cb;
    compute_comb(cb);
    return;
  }
  open_comb_ = cb;
  comb_open_ = true;
  mode_ = Mode::Comb;
}

void St1Runtime::comb_stmt(const Statement& st) {
  std::string head = upper(st[0].s);
  if (head == "FIN") {
    combs_[open_comb_.id] = open_comb_;
    compute_comb(open_comb_);
    comb_open_ = false;
    mode_ = Mode::Top;
    return;
  }
  if (head == "CHARGE") head = "CHARG";  // attested truncation-style alias (OA108)
  static const std::set<std::string> kinds = {"CHARG",   "COMB",     "ETAT",
                                             "PHASAGE", "SPECTRE",  "INCENDIE"};
  require(kinds.count(head), "ST1: COMB — source '" + st[0].s + "' inconnue");
  std::size_t i = 1;
  std::vector<Id> ids = id_list(st, i, {});
  double coef = (i < st.size() && is_num_tok(st[i])) ? num(st[i]) : 1.0;
  for (Id id : ids) open_comb_.terms.push_back({head, id, coef});
}

void St1Runtime::compute_comb(const CombDef& cb) {
  // Immediately calculated and stored (manual p.230) — sources must exist.
  std::vector<const CaseResult*> srcs;
  std::vector<double> coefs;
  for (const auto& t : cb.terms) {
    const CaseResult* r = store_.find(t.kind, t.id);
    require(r, "ST1: COMB " + std::to_string(cb.id) + " — " + t.kind + " " +
                   std::to_string(t.id) + " non calcule (EXEC prealable requis)");
    srcs.push_back(r);
    coefs.push_back(t.coef);
  }
  if (srcs.empty()) {
    warn("ST1: COMB " + std::to_string(cb.id) + " vide");
    return;
  }
  // A COMB adds stored RESULTS (p.230), nonlinear ones too: ST1 combines an
  // INCENDIE with a loading — the sum is ST1's, not a
  // nonlinear solve of the summed loads.
  std::vector<CaseResult> lin;
  lin.reserve(srcs.size());
  // A spectral response holds magnitudes; ST1 adds them as PRINTED, i.e.
  // positive in its edition of the efforts (N, T, M printed as -station
  // value): e.g. a COMB = G + |S| (N 183.17 = 77 + 107). Reactions are
  // printed as they are stored.
  for (std::size_t k = 0; k < srcs.size(); ++k)
    if (cb.terms[k].kind == "SPECTRE") {
      lin.push_back(*srcs[k]);
      for (BarResult& b : lin.back().bars)
        for (Station& st : b.stations) {
          st.N = -st.N;
          st.Vy = -st.Vy;
          st.Vz = -st.Vz;
          st.T = -st.T;
          st.My = -st.My;
          st.Mz = -st.Mz;
        }
      srcs[k] = &lin.back();
    }
  for (const CaseResult*& r : srcs)
    if (r->non_superposable) {
      lin.push_back(*r);
      lin.back().non_superposable = false;
      r = &lin.back();
    }
  if (!lin.empty())
    warn("ST1: COMB " + std::to_string(cb.id) +
         " — resultats non lineaires superposes tels quels (comme ST1)");
  CaseResult r = combine_results(srcs, coefs, cb.id, cb.name);
  store_.comb[cb.id] = std::move(r);
  stamp_cara("COMB", cb.id);
}

// ------------------------------------------------- ENV evaluation (p.231) --

namespace {

// Merge of `CaseResult`s taken ENTITY BY ENTITY (bar id, reaction node id)
// instead of by position — the ENV counterpart of `combine_results`.
//
// WHY. A staged deck (`ACTIVER` / `DESACTIVER`, e.g. pseudo.st1: the
// launching of a bridge) envelopes results computed under DIFFERENT activation
// states. Two `CaseResult`s of the same ST1 description are then NOT
// positionally comparable:
//   * `bars` holds only the ACTIVE bars, so its length changes when the nose
//     is dismantled (`DESACTIVER BARRE bec`) — `combine_results` rejects that
//     outright ("combine: shape mismatch");
//   * `reactions` follows `build_model`'s support order, i.e. the ACTIVE
//     `APPUI` nodes first and the option-auto-blocked nodes after, so the same
//     node sits at a different index from one phase to the next — a positional
//     merge would silently mix the reaction of one node with another's;
//   * `node_disp` is the only array that keeps its layout (the node set does
//     not depend on the activation).
//
// WHAT. ST1 envelopes the STRUCTURE, not one phase's layout (manual p.231):
// a phase in which an entity does not exist simply carries no term for it. The
// merge therefore runs over the UNION of the entities, and a source that does
// not carry an entity is SKIPPED — never zero-filled, which would drag a real
// extremum towards 0 and claim the phase held a null force there.
//
// `how` selects the two ST1 readings of an ENV block:
//   Sum       — `ENV i COMB` (§13.4, p.234): the sources are ADDITIVE terms;
//   Min / Max — plain `ENV i` (§13.2, p.231): the sources are ALTERNATIVES and
//               the envelope is their pointwise extreme.
//
// Each entity is templated from the FIRST source that carries it (station
// abscissae, soil-vector lengths). When the number of study sections of one bar
// differs between sources (an `ETUDE` re-declared between phases) the common
// leading grid is used and the bar id is reported in `clipped` so the caller
// can warn.
enum class EntityMerge { Sum, Min, Max };

CaseResult merge_by_entity(const std::vector<const CaseResult*>& srcs,
                           const std::vector<double>& factors, EntityMerge how, Id id,
                           const std::string& name, std::vector<Id>* clipped) {
  require(!srcs.empty() && srcs.size() == factors.size(), "ST1: ENV — sources incoherentes");
  for (const CaseResult* s : srcs)
    require(!s->non_superposable,
            "ST1: ENV — source '" + s->name + "' non superposable (elements speciaux / PSE)");

  // the bars already reported clipped (membership in O(1), the list keeps its order)
  std::unordered_set<Id> clipped_set;
  if (clipped) clipped_set.insert(clipped->begin(), clipped->end());
  CaseResult out;
  out.id = id;
  out.name = name;

  std::size_t nn = 0;
  for (const CaseResult* s : srcs) nn = std::max(nn, s->node_disp.size());
  out.node_disp.assign(nn, Vec6::Zero());
  std::vector<char> node_seen(nn, 0);

  // Union layouts, in first-seen order: deterministic run to run.
  std::map<Id, std::size_t> rindex, bindex;
  for (const CaseResult* s : srcs) {
    for (const Reaction& re : s->reactions)
      if (!rindex.count(re.node)) {
        rindex[re.node] = out.reactions.size();
        Reaction r;
        r.node = re.node;
        out.reactions.push_back(r);
      }
    for (const BarResult& b : s->bars) {
      auto it = bindex.find(b.bar);
      if (it == bindex.end()) {
        bindex[b.bar] = out.bars.size();
        BarResult t = b;  // template: abscissae + soil-vector lengths
        t.end_forces_local.setZero();
        for (Station& st : t.stations) {
          st.N = st.Vy = st.Vz = st.T = st.My = st.Mz = 0;
          st.defl.setZero();
          st.rot.setZero();
        }
        for (Vec3& p : t.soil_p) p.setZero();
        for (double& m : t.soil_m) m = 0;
        out.bars.push_back(std::move(t));
        continue;
      }
      BarResult& t = out.bars[it->second];
      if (b.stations.size() != t.stations.size()) {
        if (clipped && clipped_set.insert(b.bar).second) clipped->push_back(b.bar);
        if (b.stations.size() < t.stations.size()) {
          t.stations.resize(b.stations.size());
          if (t.soil_p.size() > t.stations.size()) t.soil_p.resize(t.stations.size());
          if (t.soil_m.size() > t.stations.size()) t.soil_m.resize(t.stations.size());
        }
      }
    }
  }
  std::vector<char> reac_seen(out.reactions.size(), 0), bar_seen(out.bars.size(), 0);

  // One scalar step of the merge. `first` = this source is the first one that
  // carries the entity, so it SETS the value instead of merging into it.
  auto acc = [&](double& dst, double v, bool first) {
    require(finite(v), "ST1: ENV — valeur non finie dans une source");
    if (how == EntityMerge::Sum) dst += v;
    else if (first) dst = v;
    else if (how == EntityMerge::Max) dst = std::max(dst, v);
    else dst = std::min(dst, v);
  };

  for (std::size_t s = 0; s < srcs.size(); ++s) {
    const CaseResult& c = *srcs[s];
    double f = factors[s];
    for (std::size_t i = 0; i < c.node_disp.size() && i < out.node_disp.size(); ++i) {
      bool first = !node_seen[i];
      for (int k = 0; k < 6; ++k) acc(out.node_disp[i](k), f * c.node_disp[i](k), first);
      node_seen[i] = 1;
    }
    for (const Reaction& re : c.reactions) {
      std::size_t ix = rindex.at(re.node);
      bool first = !reac_seen[ix];
      for (int k = 0; k < 6; ++k) acc(out.reactions[ix].R(k), f * re.R(k), first);
      reac_seen[ix] = 1;
    }
    for (const BarResult& b : c.bars) {
      std::size_t ix = bindex.at(b.bar);
      BarResult& d = out.bars[ix];
      bool first = !bar_seen[ix];
      for (int k = 0; k < 12; ++k)
        acc(d.end_forces_local(k), f * b.end_forces_local(k), first);
      for (std::size_t k = 0; k < d.stations.size() && k < b.stations.size(); ++k) {
        Station& a = d.stations[k];
        const Station& q = b.stations[k];
        acc(a.N, f * q.N, first);
        acc(a.Vy, f * q.Vy, first);
        acc(a.Vz, f * q.Vz, first);
        acc(a.T, f * q.T, first);
        acc(a.My, f * q.My, first);
        acc(a.Mz, f * q.Mz, first);
        for (int u = 0; u < 3; ++u) acc(a.defl(u), f * q.defl(u), first);
        for (int u = 0; u < 3; ++u) acc(a.rot(u), f * q.rot(u), first);
      }
      for (std::size_t k = 0; k < d.soil_p.size() && k < b.soil_p.size(); ++k)
        for (int u = 0; u < 3; ++u) acc(d.soil_p[k](u), f * b.soil_p[k](u), first);
      for (std::size_t k = 0; k < d.soil_m.size() && k < b.soil_m.size(); ++k)
        acc(d.soil_m[k], f * b.soil_m[k], first);
      bar_seen[ix] = 1;
    }
  }
  return out;
}

}  // namespace

// Field-by-field lo/hi bounds of a SURCH: `run_surch` stores, per studied
// section, the two governing CONCOMITANT states (surch.h `SurchExtreme::state`).
// A SURCH is a SET of possible states, not an additive action, so it cannot be
// handed to `env_run` as a single CaseResult: keta envelopes ALL its stored
// states into two pseudo-cases and feeds them as the alternative set of one
// `EnvSource` (etude.h's documented ENV-of-ENV mechanism). DEVIATION, documented:
// the pair is GLOBAL (one lo/hi for the whole SURCH) while ST1 keeps a governing
// placement per section, so the companion state of a section other than the
// governing one is the envelope of the stored states, not a replayed placement.
const Envelope& St1Runtime::surch_bounds(Id id) {
  auto cached = surch_bounds_.find(id);
  if (cached != surch_bounds_.end()) return cached->second;
  // ST1 edition (p.137): the per-section extremes of `run_surch_envelope`,
  // laid on a result of the EXEC CHARG shape (same stations, same reactions)
  // so ENV / ENV COMB fold them field by field with the loadings.
  auto ev = surch_envelopes.find(id);
  if (ev != surch_envelopes.end()) {
    build_model();
    Model mm = model_;
    mm.cases.clear();
    mm.combos.clear();
    LoadCase lc;
    lc.id = id;
    lc.name = "SURCH";
    mm.cases.push_back(lc);
    CaseResult tmpl = solve_static(mm, recovery_options()).cases.front();
    for (Vec6& u : tmpl.node_disp) u.setZero();
    Envelope e;
    e.lo.value = tmpl;
    e.hi.value = tmpl;
    auto station_at = [&](CaseResult& cr, Id bar, double x) -> Station* {
      BarResult& br = cr.bars[static_cast<std::size_t>(mm.bar_index(bar))];
      Station* best = nullptr;
      for (Station& st : br.stations)
        if (!best || std::fabs(st.x - x) < std::fabs(best->x - x)) best = &st;
      return best;
    };
    for (int side = 0; side < 2; ++side) {
      CaseResult& cr = side == 0 ? e.lo.value : e.hi.value;
      for (const SurchEnvRow& r : ev->second.efforts) {
        Station* st = station_at(cr, r.id, r.x);
        if (!st) continue;
        const std::array<double, 6>& v = side == 0 ? r.lo : r.hi;
        st->N = v[0];
        st->Vy = v[1];
        st->Vz = v[2];
        st->T = v[3];
        st->My = v[4];
        st->Mz = v[5];
      }
      for (const SurchEnvRow& r : ev->second.deplas) {
        Station* st = station_at(cr, r.id, r.x);
        if (!st) continue;
        const std::array<double, 6>& v = side == 0 ? r.lo : r.hi;
        const Bar& b = mm.bar(r.id);
        Mat3 R = bar_rotation(mm.node(b.ni).x, mm.node(b.nj).x, b.roll);
        st->defl = R * Vec3(v[0], v[1], v[2]);
        st->rot = R * Vec3(v[3], v[4], v[5]);
        // PRESS: the soil pressures of the row (SurchEnvSpec::depla_fibres)
        const std::vector<double>& pf = side == 0 ? r.fibre_lo : r.fibre_hi;
        BarResult& br = cr.bars[static_cast<std::size_t>(mm.bar_index(r.id))];
        const std::size_t k = static_cast<std::size_t>(st - br.stations.data());
        if (!pf.empty() && k < br.soil_p.size()) {
          Vec3 q = Vec3::Zero();
          double mx = 0;
          if (option_ == 0) q = Vec3(pf[0], pf.size() > 1 ? pf[1] : 0.0, 0.0);
          else if (option_ == 1) {
            q = Vec3(0.0, 0.0, pf[0]);
            mx = pf.size() > 1 ? pf[1] : 0.0;
          } else if (pf.size() >= 4) {
            q = Vec3(pf[0], pf[1], pf[2]);
            mx = pf[3];
          }
          br.soil_p[k] = q;
          if (k < br.soil_m.size()) br.soil_m[k] = mx;
        }
      }
      for (const SurchEnvRow& r : ev->second.reactions)
        for (Reaction& re : cr.reactions)
          if (re.node == r.id)
            for (int c = 0; c < 6; ++c)
              re.R(c) = side == 0 ? r.lo[static_cast<std::size_t>(c)] : r.hi[static_cast<std::size_t>(c)];
    }
    surch_bounds_[id] = std::move(e);
    return surch_bounds_.at(id);
  }
  auto it = surch_results.find(id);
  require(it != surch_results.end(),
          "ST1: SURCH " + std::to_string(id) + " non calculee (EXEC SURCH prealable requis)");
  std::vector<const CaseResult*> states;
  for (const SurchTargetResult& t : it->second.sections) {
    if (t.lo.has_state) states.push_back(&t.lo.state);
    if (t.hi.has_state) states.push_back(&t.hi.state);
  }
  for (const SurchTargetResult& t : it->second.reactions) {
    if (t.lo.has_state) states.push_back(&t.lo.state);
    if (t.hi.has_state) states.push_back(&t.hi.state);
  }
  require(!states.empty(), "ST1: SURCH " + std::to_string(id) +
                               " — aucun etat concomitant stocke (option concomitance requise)");
  surch_bounds_[id] = envelope(states);
  return surch_bounds_.at(id);
}

void St1Runtime::compute_env(const EnvDef& def) {
  if (def.sources.empty()) {
    warn("ST1: ENV " + std::to_string(def.id) + " sans source, ignoree");
    return;
  }
  build_model();

  // ST1 has TWO envelope commands and they do NOT compose their sources the
  // same way:
  //
  //   * plain `ENV i` (§13.2, p.231) — the sources are ALTERNATIVES. The
  //     manual's own decks say so unambiguously: exemple2.st1 writes
  //     `ENV 1 'enveloppe de la poussee des terres : coef. : 0.25 et 0.50'`
  //     over `CHARG 2 .25` + `CHARG 2 .50` (one earth-pressure case taken at
  //     its two extreme coefficients, NOT 0.75x it), exemple6.st1 envelopes
  //     four ULS combinations under `ENV 1000`, and pseudo.st1 accumulates the
  //     launching phases by re-declaring `ENV 1` with itself plus the current
  //     phase. Consistently, p.231 states an envelope "contient donc deux fois
  //     plus de resultats qu'un cas de charge" — a fixed-coefficient SUM would
  //     hold a single state.
  //   * `ENV i COMB` (§13.4, p.234, "COMBINAISONS ENVELOPPES") — the sources
  //     are ADDITIVE terms, each carrying a [cfmin, cfmax] box; the extremes
  //     are taken at the corners of that box.
  //
  // Resolution is shared: a plain result contributes its own CaseResult, a
  // SURCH or a nested ENV contributes its two bounds as ALTERNATIVES (an
  // envelope is a set of states, never a sum — etude.h `EnvSource::alts`).
  std::vector<EnvSource> sources;
  // Corner sets used for the plain field-by-field bounds (see below).
  std::vector<std::vector<std::pair<const CaseResult*, double>>> corners;
  // Plain ENV only: the flat alternative set and the coefficient-scaled copies
  // it points at (they must outlive `sources`, hence the reserve — no
  // reallocation may move them).
  std::vector<std::pair<const CaseResult*, double>> plain_alts;
  std::vector<CaseResult> plain_scaled;
  // CONTR envelope (see EnvResult::contr): every alternative with the sections
  // it was computed with; a nested ENV brings its own stress bounds. One group
  // per source under ENV COMB (its corners), a single group otherwise.
  struct ContrAlt {
    const CaseResult* r = nullptr;
    const std::map<Id, CaraProps>* cara = nullptr;
    const std::map<Id, std::vector<EnvResult::ContrBound>>* bounds = nullptr;  // ENV / SURCH
    // CONTRG bounds of an ENV / SURCH source (CONTRG id -> bar -> bounds)
    const std::map<Id, std::map<Id, std::vector<EnvResult::ContrBound>>>* gbounds = nullptr;
    bool cst = true;  // the CONTRG constant applies (loadings, not surcharges: p.99)
    double c = 1.0;
  };
  std::vector<std::vector<ContrAlt>> contr_groups(def.comb ? 0 : 1);
  // Reactions per APPUI, likewise: every alternative with the APPUI its
  // source was computed with.
  using AppuiMap = std::map<Id, Vec6>;
  std::vector<std::vector<std::pair<AppuiMap, double>>> appui_groups(def.comb ? 0 : 1);
  for (const EnvDef::Source& s : def.sources) {
    // Plain ENV: one coefficient (cmin == cmax). ENV COMB: the [cfmin, cfmax]
    // box, extremes taken at its corners (p.234).
    double cmin = s.cmin, cmax = s.cmax;
    std::vector<const CaseResult*> alts;
    std::vector<ContrAlt> calts;
    // A source whose result is not in the database is DROPPED with a loud
    // warning rather than aborting the session: an envelope routinely mixes
    // families the core does not own (SURCH, SPECTRE, HISTORIQUE…), and a deck
    // must stay runnable while their engine is wired. COMB stays strict
    // (manual p.230) — only ENV is lenient, and never silently.
    if (s.kind == "SURCH") {
      if (!surch_results.count(s.id)) {
        warn("ST1: ENV " + std::to_string(def.id) + " — SURCH " + std::to_string(s.id) +
             " non calculee, source ignoree");
        continue;
      }
      const Envelope& e = surch_bounds(s.id);
      alts = {&e.lo.value, &e.hi.value};
      if (auto sc = surch_contr.find(s.id); sc != surch_contr.end()) {
        ContrAlt ca;
        ca.bounds = &sc->second;
        calts = {ca};
      } else {
        calts = {{&e.lo.value, &cara_of("SURCH", s.id)}, {&e.hi.value, &cara_of("SURCH", s.id)}};
      }
      auto sg = surch_contrg.find(s.id);
      for (ContrAlt& ca : calts) {
        ca.cst = false;
        if (sg != surch_contrg.end()) ca.gbounds = &sg->second;
      }
    } else if (s.kind == "HISTORIQUE") {
      auto it = historique_env.find(s.id);
      if (it == historique_env.end()) {
        warn("ST1: ENV " + std::to_string(def.id) + " — HISTORIQUE " + std::to_string(s.id) +
             " non calcule, source ignoree");
        continue;
      }
      // the time series when kept: every step is a concomitant state (the
      // field bounds mix components of different instants)
      if (auto ser = historique_series.find(s.id);
          ser != historique_series.end() && !ser->second.states.empty()) {
        alts.clear();
        calts.clear();
        for (const CaseResult& st : ser->second.states) {
          alts.push_back(&st);
          calts.push_back({&st, &cara_});
        }
      } else {
        alts = {&it->second.first, &it->second.second};
        calts = {{alts[0], &cara_}, {alts[1], &cara_}};
      }
    } else if (s.kind == "ENV") {
      auto it = env_full.find(s.id);
      if (it == env_full.end()) {
        warn("ST1: ENV " + std::to_string(def.id) + " — ENV " + std::to_string(s.id) +
             " non calculee, source ignoree");
        continue;
      }
      alts = {&it->second.bounds.lo.value, &it->second.bounds.hi.value};
      ContrAlt ca;
      ca.bounds = &it->second.contr;
      ca.gbounds = &it->second.contrg_b;
      calts = {ca};
    } else {
      const CaseResult* r = store_.find(s.kind, s.id);
      if (!r) {
        warn("ST1: ENV " + std::to_string(def.id) + " — " + s.kind + " " +
             std::to_string(s.id) + " non calcule, source ignoree");
        continue;
      }
      alts = {r};
      calts = {{r, &cara_of(s.kind, s.id)}};
    }
    {
      std::vector<AppuiMap> amaps;
      if (s.kind == "ENV") {
        const EnvResult& ne = env_full.at(s.id);
        AppuiMap lo, hi;
        for (const auto& [aid, b] : ne.appui_bounds) {
          lo[aid] = b.first;
          hi[aid] = b.second;
        }
        amaps = {lo, hi};
      } else {
        const std::set<Id> act = appuis_of(s.kind, s.id);
        const std::map<Id, SupportDef> defs = appui_defs_of(s.kind, s.id);
        for (const CaseResult* a : alts) amaps.push_back(appui_reactions(*a, act, true, &defs));
      }
      auto& ag = def.comb ? appui_groups.emplace_back() : appui_groups.front();
      for (const AppuiMap& am : amaps) {
        ag.push_back({am, cmax});
        if (cmin != cmax) ag.push_back({am, cmin});
      }
    }
    {
      std::vector<ContrAlt>& g = def.comb ? contr_groups.emplace_back() : contr_groups.front();
      for (ContrAlt ca : calts) {
        ca.c = cmax;
        g.push_back(ca);
        if (cmin != cmax) {
          ca.c = cmin;
          g.push_back(ca);
        }
      }
    }
    if (def.comb) {
      EnvSource es;
      es.cmin = cmin;
      es.cmax = cmax;
      if (alts.size() == 1) es.r = alts[0];
      else es.alts = alts;
      sources.push_back(es);

      std::vector<std::pair<const CaseResult*, double>> cs;
      for (const CaseResult* a : alts) {
        cs.emplace_back(a, cmax);
        if (cmin != cmax) cs.emplace_back(a, cmin);
      }
      corners.push_back(std::move(cs));
    } else {
      // Plain ENV: every alternative of every source joins ONE alternative set,
      // already scaled by the source's coefficient. §13.2 gives ONE coefficient
      // per source; a deck that writes two anyway gets both as alternatives
      // rather than a silently dropped one.
      for (const CaseResult* a : alts) {
        plain_alts.emplace_back(a, cmax);
        if (cmin != cmax) plain_alts.emplace_back(a, cmin);
      }
    }
  }
  if (!def.comb) {
    // One `EnvSource` holding the whole alternative set: `pick_corner`
    // (etude.h) then returns the single most favourable alternative for the
    // privileged run, which is exactly ST1's plain envelope.
    plain_scaled.reserve(plain_alts.size());
    for (const auto& [r, c] : plain_alts)
      plain_scaled.push_back(combine_results({r}, {c}, 0, "alternative"));
    if (!plain_scaled.empty()) {  // every source dropped: the guard below warns
      EnvSource es;
      es.cmin = es.cmax = 1.0;
      es.skip_null = true;
      for (const CaseResult& r : plain_scaled) es.alts.push_back(&r);
      for (const auto& [r, c] : plain_alts) es.raw.push_back(r);
      sources.push_back(es);
    }
  }
  if (sources.empty()) {
    warn("ST1: ENV " + std::to_string(def.id) + " — aucune source calculee, enveloppe ignoree");
    return;
  }

  EnvResult out;
  out.id = def.id;
  out.name = def.name;
  out.has_priv = def.has_priv;
  out.what = def.what;
  out.comp = def.comp;
  out.contrg = def.contrg;
  out.fiber = def.fiber;

  // --- field-by-field bounds -------------------------------------------
  // Sources may come from DIFFERENT activation states (a launching deck
  // re-envelopes phase after phase), so every merge below runs by entity id
  // rather than by position — see `merge_by_entity`.
  //
  // Section deformations are enveloped component by component in GLOBAL
  // axes, as ST1 lists them (and as the SURCH editions already are): the
  // merged copies are rotated to global first, the bounds back to the local
  // frame of the CaseResult convention afterwards (RESU rotates them again).
  std::map<Id, Mat3> rot_cache;
  const std::unordered_map<Id, const Bar*> model_bar = bars_by_id(model_);  // model_ is fixed from here
  auto bar_rot = [&](Id bar, Mat3& R) {
    auto c = rot_cache.find(bar);
    if (c != rot_cache.end()) {
      R = c->second;
      return true;
    }
    if (auto it = model_bar.find(bar); it != model_bar.end()) {
      const Bar* mb = it->second;
      R = bar_rotation(model_.node(mb->ni).x, model_.node(mb->nj).x, mb->roll);
    } else {
      auto rb = bars_.find(bar);
      if (rb == bars_.end()) return false;
      R = bar_rotation(nodes_.at(rb->second.n1), nodes_.at(rb->second.n2), 0.0);
    }
    rot_cache[bar] = R;
    return true;
  };
  auto deform_frame = [&](CaseResult& r, bool to_global) {
    for (BarResult& b : r.bars) {
      Mat3 R;
      if (!bar_rot(b.bar, R)) continue;
      const Mat3 M = to_global ? Mat3(R.transpose()) : R;
      for (Station& st : b.stations) {
        st.defl = M * st.defl;
        st.rot = M * st.rot;
      }
    }
  };
  std::vector<Id> clipped;
  if (def.comb) {
    // ENV COMB: the extremes of a SUM of independent terms are the sum of the
    // per-term extremes, so each source is first reduced to its own (lo, hi)
    // pseudo-case over its corner set and the reduced cases are then added.
    std::vector<CaseResult> per_lo, per_hi;
    for (const auto& cs : corners) {
      std::vector<CaseResult> scaled;
      scaled.reserve(cs.size());
      for (const auto& [r, c] : cs) scaled.push_back(combine_results({r}, {c}, 0, "corner"));
      for (CaseResult& r : scaled) deform_frame(r, true);
      std::vector<const CaseResult*> ptrs;
      std::vector<double> ones(scaled.size(), 1.0);
      for (const CaseResult& r : scaled) ptrs.push_back(&r);
      per_lo.push_back(
          merge_by_entity(ptrs, ones, EntityMerge::Min, 0, "corner min", &clipped));
      per_hi.push_back(
          merge_by_entity(ptrs, ones, EntityMerge::Max, 0, "corner max", &clipped));
    }
    std::vector<const CaseResult*> lo_ptr, hi_ptr;
    for (const CaseResult& r : per_lo) lo_ptr.push_back(&r);
    for (const CaseResult& r : per_hi) hi_ptr.push_back(&r);
    std::vector<double> ones(lo_ptr.size(), 1.0);
    out.bounds.lo.value =
        merge_by_entity(lo_ptr, ones, EntityMerge::Sum, def.id, def.name + " min", &clipped);
    out.bounds.hi.value =
        merge_by_entity(hi_ptr, ones, EntityMerge::Sum, def.id, def.name + " max", &clipped);
  } else {
    // Plain ENV: the pointwise extremes of the alternative set.
    std::vector<CaseResult> glob = plain_scaled;
    for (CaseResult& r : glob) deform_frame(r, true);
    std::vector<const CaseResult*> ptrs;
    for (const CaseResult& r : glob) ptrs.push_back(&r);
    std::vector<double> ones(ptrs.size(), 1.0);
    out.bounds.lo.value =
        merge_by_entity(ptrs, ones, EntityMerge::Min, def.id, def.name + " min", &clipped);
    out.bounds.hi.value =
        merge_by_entity(ptrs, ones, EntityMerge::Max, def.id, def.name + " max", &clipped);
  }
  deform_frame(out.bounds.lo.value, false);
  deform_frame(out.bounds.hi.value, false);
  // An ENV lists every DECLARED support, a loading only its active ones: a
  // support no source carries reads 0 (pushing deck `pseudo`: the bearings
  // inactive in the first enveloped phase keep a 0 minimum, those active in
  // it their true one, since later inactive phases add nothing).
  for (CaseResult* cr : {&out.bounds.lo.value, &out.bounds.hi.value}) {
    std::set<Id> have;
    for (const Reaction& re : cr->reactions) have.insert(re.node);
    for (const auto& [aid, ad] : appui_) {
      (void)aid;
      if (have.insert(ad.node).second) cr->reactions.push_back({ad.node, Vec6::Zero()});
    }
  }
  // --- reactions per APPUI ------------------------------------------------
  // Per group the componentwise extremes over the alternatives that carry
  // the APPUI (one inactive in a loading is skipped), the groups summed (ENV
  // COMB; a group without it adds nothing). An ENV holds every declared
  // APPUI, 0 where no source had it active, and a nested ENV passes those 0
  // on — ST1's pseudo: appui 15 Min 993.51 (skipped while inactive), appui
  // 16 Min 0 (inactive in every source of the ENV it extends).
  {
    std::set<Id> keys(appui_order_.begin(), appui_order_.end());
    for (const auto& g : appui_groups)
      for (const auto& [am, c] : g)
        for (const auto& [aid, v] : am) keys.insert(aid);
    for (Id aid : keys) {
      Vec6 lo = Vec6::Zero(), hi = Vec6::Zero();
      for (const auto& g : appui_groups) {
        if (g.empty()) continue;
        Vec6 glo = Vec6::Constant(1e300), ghi = Vec6::Constant(-1e300);
        bool any = false;
        for (const auto& [am, c] : g) {
          auto it = am.find(aid);
          if (it == am.end()) continue;
          const Vec6 v = c * it->second;
          glo = glo.cwiseMin(v);
          ghi = ghi.cwiseMax(v);
          any = true;
        }
        if (!any) continue;
        lo += glo;
        hi += ghi;
      }
      out.appui_bounds[aid] = {lo, hi};  // never active: 0
    }
  }
  // --- CONTR bounds -----------------------------------------------------
  // On the grid of the effort bounds: per group the fibre-wise extremes over
  // its alternatives, the groups summed (ENV COMB) — a source without the
  // station adds nothing.
  {
    auto near = [](double a, double b) { return std::fabs(a - b) <= 1e-6 * (1.0 + std::fabs(a)); };
    // the fibre rows of one alternative at (bar, x): one row, or two (nested ENV)
    auto rows_of = [&](const ContrAlt& ca, Id bar, double x) {
      std::vector<std::vector<double>> rows;
      auto scaled = [&](std::vector<double> v) {
        for (double& t : v) t *= ca.c;
        return v;
      };
      if (ca.bounds) {
        auto b = ca.bounds->find(bar);
        if (b == ca.bounds->end()) return rows;
        for (const EnvResult::ContrBound& cb : b->second)
          if (near(cb.x, x)) {
            rows.push_back(scaled(cb.lo));
            rows.push_back(scaled(cb.hi));
            break;
          }
        return rows;
      }
      for (const BarResult& br : ca.r->bars) {
        if (br.bar != bar) continue;
        for (const Station& st : br.stations)
          if (near(st.x, x)) {
            rows.push_back(scaled(contr_fibres(*ca.cara, bar, st)));
            break;
          }
        break;
      }
      return rows;
    };
    for (const BarResult& br : out.bounds.lo.value.bars) {
      std::vector<EnvResult::ContrBound>& cbs = out.contr[br.bar];
      for (const Station& st : br.stations) {
        EnvResult::ContrBound cb;
        cb.x = st.x;
        bool first_group = true;
        for (const std::vector<ContrAlt>& g : contr_groups) {
          std::vector<double> glo, ghi;
          for (const ContrAlt& ca : g)
            for (const std::vector<double>& row : rows_of(ca, br.bar, st.x)) {
              if (glo.empty()) {
                glo = ghi = row;
                continue;
              }
              for (std::size_t k = 0; k < row.size() && k < glo.size(); ++k) {
                glo[k] = std::min(glo[k], row[k]);
                ghi[k] = std::max(ghi[k], row[k]);
              }
            }
          if (glo.empty()) continue;
          if (first_group) {
            cb.lo = glo;
            cb.hi = ghi;
            first_group = false;
          } else {
            for (std::size_t k = 0; k < glo.size() && k < cb.lo.size(); ++k) {
              cb.lo[k] += glo[k];
              cb.hi[k] += ghi[k];
            }
          }
        }
        if (!cb.lo.empty()) cbs.push_back(std::move(cb));
      }
    }
  }
  // --- CONTRG bounds ----------------------------------------------------
  // Each CONTRG on the grid of the effort bounds, as the CONTR above: per
  // group the extremes over the alternatives (a loading with its constant, a
  // nested ENV / SURCH with its own bounds), the groups summed.
  for (const auto& [cid, cd] : contrgs) {
    auto near = [](double a, double b) { return std::fabs(a - b) <= 1e-6 * (1.0 + std::fabs(a)); };
    auto values_of = [&](const ContrAlt& ca, Id bar, double x, double L) {
      std::vector<double> v;
      if (ca.gbounds) {
        auto g = ca.gbounds->find(cid);
        if (g == ca.gbounds->end()) return v;
        auto b = g->second.find(bar);
        if (b == g->second.end()) return v;
        for (const EnvResult::ContrBound& cb : b->second)
          if (near(cb.x, x) && !cb.lo.empty()) {
            v = {ca.c * cb.lo[0], ca.c * cb.hi[0]};
            break;
          }
        return v;
      }
      if (!ca.r) return v;
      for (const BarResult& br : ca.r->bars) {
        if (br.bar != bar) continue;
        for (const Station& st : br.stations)
          if (near(st.x, x)) {
            const double xx = L > 0 ? std::clamp(st.x / L, 0.0, 1.0) : 0.0;
            v = {ca.c * contrg_at(cd, bar, xx, st, ca.cst)};
            break;
          }
        break;
      }
      return v;
    };
    for (const BarResult& br : out.bounds.lo.value.bars) {
      if (!cd.defined_for(br.bar)) continue;
      const double L = bar_length(br.bar);
      std::vector<EnvResult::ContrBound>& cbs = out.contrg_b[cid][br.bar];
      for (const Station& st : br.stations) {
        double lo = 0, hi = 0;
        bool any = false;
        for (const std::vector<ContrAlt>& g : contr_groups) {
          double glo = 1e300, ghi = -1e300;
          for (const ContrAlt& ca : g)
            for (double v : values_of(ca, br.bar, st.x, L)) {
              glo = std::min(glo, v);
              ghi = std::max(ghi, v);
            }
          if (glo > ghi) continue;
          lo += glo;
          hi += ghi;
          any = true;
        }
        if (any) cbs.push_back({st.x, {lo}, {hi}});
      }
    }
  }
  if (!clipped.empty()) {
    std::string list;
    for (Id b : clipped) list += (list.empty() ? "" : ",") + std::to_string(b);
    warn("ST1: ENV " + std::to_string(def.id) +
         " — sections d'etude differentes entre sources sur la/les barre(s) " + list +
         "; enveloppe sur la grille commune la plus courte");
  }

  // --- privileged run ---------------------------------------------------
  if (def.has_priv) {
    EnvSpec spec;
    spec.what = def.what;
    spec.comp = def.comp;
    spec.fiber = def.fiber;
    spec.sigma_st1_sign = true;  // ST1 lists normal stresses compression-positive (p.66)
    if (def.what == EnvWhat::Contrg) {
      auto it = contrgs.find(def.contrg);
      require(it != contrgs.end(),
              "ST1: ENV CONTRG " + std::to_string(def.contrg) + " non definie");
      spec.contrg = &it->second;
    }
    if (def.what == EnvWhat::Depla) {
      // A DEPLA envelope has no section: its extremes ARE the field-by-field
      // bounds above (exact, per node), only the cross-component concomitance
      // is lost — documented deviation from ST1 §13.3.1.
      warn("ST1: ENV " + std::to_string(def.id) +
           " DEPLA — extremes nodaux exacts, concomitance inter-composantes non tenue");
    } else {
      out.points = env_run_sections(model_, spec, sources);
      // One representative full run for the `env_results` seam store.
      if (!out.points.empty()) {
        spec.bar = out.points.front().bar;
        spec.station = out.points.front().station;
        if (def.has_reac) {
          spec.has_reac = true;
          spec.reac_comp = def.reac_comp;
          out.reac_comp = def.reac_comp;
          // the APPUI named, else every APPUI carrying a reaction
          std::vector<Id> aps;
          if (def.reac_appui != 0) aps.push_back(def.reac_appui);
          else
            for (Id a : all_appui_ids()) aps.push_back(a);
          const CaseResult& shape =
              sources.front().r ? *sources.front().r : *sources.front().alts.front();
          for (Id a : aps) {
            const Id node = appui_.count(a) ? appui_.at(a).node : a;
            if (std::none_of(shape.reactions.begin(), shape.reactions.end(),
                             [&](const Reaction& r) { return r.node == node; }))
              continue;
            spec.reac_node = node;
            const EnvRun rr = env_run(model_, spec, sources);
            out.reac_concom[a] = {rr.reac_state_min, rr.reac_state_max};
            if (!out.has_run) {
              out.run = rr;
              out.has_run = true;
            }
          }
        }
        if (!out.has_run) {
          spec.has_reac = false;
          out.run = env_run(model_, spec, sources);
          out.has_run = true;
        }
      }
    }
  }
  if (out.has_run) env_results[def.id] = out.run;
  env_full[def.id] = std::move(out);
  stamp_cara("ENV", def.id);
  out_.line("ENV " + std::to_string(def.id) + " : " + std::to_string(def.sources.size()) +
            " source(s) enveloppees");
}

// ------------------------------------------------------------------ ETUDE --
// `ETUDE (SUP, SEUL) (EFFORT) (DEPLA)` + `liste1 SE liste2 (PAS lpas) (ABS,REL)`
// lines (manual p.97-98). The block is stored as an `EtudeSpec` (etude.h) and
// replayed onto an `EtudePlan`, so ST1's accumulation rules live in the module
// that documents them.
void St1Runtime::cmd_etude(const Statement& st) {
  // Compat: `ETUDE PAS n` alone is the old global station-count hint.
  if (st.size() <= 3)
    for (std::size_t i = 1; i + 1 < st.size(); ++i)
      if (upper(st[i].s) == "PAS" && is_num_tok(st[i + 1])) {
        etude_pas_hint = num(st[i + 1]);
        return;
      }
  open_etude_ = EtudeSpec{};
  std::size_t i = 1;
  for (; i < st.size(); ++i) {
    std::string k = st[i].str ? "" : upper(st[i].s);
    if (k == "SUP") open_etude_.mode = EtudeMode::Sup;
    else if (k == "SEUL") open_etude_.mode = EtudeMode::Seul;
    else if (k == "EFFORT") open_etude_.effort = true;
    else if (k == "DEPLA") open_etude_.depla = true;
    else break;  // the header stops at the first content token
  }
  etude_open_ = true;
  mode_ = Mode::Etude;
  // Single-statement form: header AND content line on the same statement
  // (`ETUDE EFFORT DEPLA TOUT SE 0 A 1 PAS 0.1`).
  if (i < st.size()) {
    Statement rest(st.begin() + static_cast<long>(i), st.end());
    etude_stmt(rest);
  }
}

void St1Runtime::etude_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    if (etude_open_ && !open_etude_.lines.empty()) etude_.blocks.push_back(open_etude_);
    etude_open_ = false;
    etude_cache_valid_ = false;
    mode_ = Mode::Top;
    return;
  }
  std::size_t i = 0;
  bool all = head == "TOUT" || head == "TOUTES";
  std::vector<Id> bars = id_list(st, i, all_bar_ids());
  if (bars.empty() && !all) {
    // Not an ETUDE line: close the block and re-dispatch.
    if (etude_open_ && !open_etude_.lines.empty()) etude_.blocks.push_back(open_etude_);
    etude_open_ = false;
    etude_cache_valid_ = false;
    mode_ = Mode::Top;
    execute(st);
    return;
  }
  EtudeLine line;
  line.all = all;
  if (!all) line.bars = bars;
  while (i < st.size()) {
    std::string k = st[i].str ? "" : upper(st[i].s);
    if (k == "SE") {
      ++i;
      // `SE v1 (A v2) (, v3 …)`: the manual's own `SE 0. a 1. PAS .1` form is a
      // RANGE, which etude.h reads as the two-value abscissa list + step.
      while (i < st.size()) {
        if (st[i].s == "," || (!st[i].str && upper(st[i].s) == "A")) { ++i; continue; }
        if (!is_num_tok(st[i])) break;
        line.abscissas.push_back(num(st[i++]));
      }
      continue;
    }
    if (k == "PAS" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      line.step = num(st[i + 1]);
      i += 2;
      continue;
    }
    if (k == "REL") { line.absolute = false; ++i; continue; }
    if (k == "ABS") { line.absolute = true; ++i; continue; }
    warn("ST1: ETUDE — mot '" + st[i].s + "' ignore");
    ++i;
  }
  require(!line.abscissas.empty() || line.step > 0, "ST1: ETUDE — 'SE liste2' requis");
  open_etude_.lines.push_back(std::move(line));
  etude_cache_valid_ = false;
}

// Lowering of the whole ETUDE declaration to per-bar abscissas [m].
std::map<Id, std::vector<double>> St1Runtime::etude_stations(StudyFamily which) {
  if (!etude_.any() && which == StudyFamily::Any) return {};
  build_model();
  return etude_stations_for(model_, which);
}

std::map<Id, std::vector<double>> St1Runtime::etude_stations_for(const Model& m,
                                                                 StudyFamily which) const {
  if (!etude_.any() && which == StudyFamily::Any) return {};
  EtudePlan plan;
  // ST1 default: an EFFORT section at each end of every bar (p.97). A first
  // block in SUP mode adds to it; a Default block re-installs it; SEUL drops it.
  plan.add_default(m);
  for (const EtudeSpec& spec : etude_.blocks) plan.add(m, spec);
  return plan.stations_abs(which);
}

StaticOptions St1Runtime::recovery_options() {
  StaticOptions sopt;
  if (etude_pas_hint > 0) sopt.stations_per_bar = std::max(2, static_cast<int>(etude_pas_hint) + 1);
  sopt.stations_abs = etude_stations_ref();
  return sopt;
}

const std::map<Id, std::vector<double>>* St1Runtime::etude_stations_ref() {
  if (!etude_.any()) return nullptr;
  if (!etude_cache_valid_) {
    etude_cache_ = etude_stations();
    etude_cache_valid_ = true;
  }
  return &etude_cache_;
}

// ----------------------------------------------------------- MASSE PROPRE --
// `MASSE PROPRE (…) FIN` (manual p.69-71). Both the simple form (everything on
// the MASSE PROPRE line) and the block form are accepted; the content lines are
// handled by `masse_stmt`.
void St1Runtime::cmd_masse(const Statement& st) {
  masse_pending_nodes_.clear();
  masse_pending_bars_.clear();
  masse_exc_or_ = Vec3::Zero();
  masse_exc_ex_ = Vec3::Zero();
  mode_ = Mode::Masse;
  std::size_t i = 1;
  if (i < st.size() && !st[i].str && upper(st[i].s) == "PROPRE") ++i;
  if (i < st.size()) {
    Statement rest(st.begin() + static_cast<long>(i), st.end());
    masse_stmt(rest);
  }
}

void St1Runtime::masse_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    mode_ = Mode::Top;
    invalidate();
    return;
  }
  std::size_t i = 0;
  if (head == "NOEUD") {
    i = 1;
    masse_pending_nodes_ = id_list(st, i, all_node_ids());
    masse_pending_bars_.clear();
    require(!masse_pending_nodes_.empty(), "ST1: MASSE PROPRE NOEUD sans liste");
  } else if (head == "BARRE") {
    i = 1;
    masse_pending_bars_ = id_list(st, i, all_bar_ids());
    masse_pending_nodes_.clear();
    masse_exc_or_ = Vec3::Zero();
    masse_exc_ex_ = Vec3::Zero();
    require(!masse_pending_bars_.empty(), "ST1: MASSE PROPRE BARRE sans liste");
    // Optional `EXC (OR (X) x (Y) y (Z) z) (EX …)` — global-frame offsets.
    if (i < st.size() && upper(st[i].s) == "EXC") {
      ++i;
      while (i < st.size()) {
        std::string end = upper(st[i].s);
        if (end != "OR" && end != "EX") break;
        ++i;
        std::vector<double> v;
        while (i < st.size()) {
          std::string k = st[i].str ? "" : upper(st[i].s);
          if (k == "X" || k == "Y" || k == "Z") { ++i; continue; }
          if (!is_num_tok(st[i])) break;
          v.push_back(num(st[i++]));
        }
        Vec3 e = Vec3::Zero();
        for (std::size_t c = 0; c < 3 && c < v.size(); ++c) e(static_cast<int>(c)) = v[c];
        (end == "OR" ? masse_exc_or_ : masse_exc_ex_) = e;
      }
    }
  }
  masse_.selective = true;
  // `BARRE listeb (MASSE (DI))`: these bars carry their automatic rho*A mass
  auto mark_auto = [&]() {
    for (Id b : masse_pending_bars_) {
      if (std::any_of(masse_.bars.begin(), masse_.bars.end(), [&](const BarMass& x) {
            return x.bar == b && x.kind == BarMassKind::Auto;
          }))
        continue;
      BarMass bm;
      bm.bar = b;
      bm.kind = BarMassKind::Auto;
      masse_.bars.push_back(bm);
    }
  };
  if (i >= st.size()) {
    // list-only line: the default diagonal masses unless a MASSE line follows
    if (!masse_pending_bars_.empty()) {
      masse_.lumped_default = true;
      mark_auto();
      invalidate();
    }
    return;
  }
  // `MASSE` is the manual's keyword on every form (p.69-71:
  // `NOEUD liste MASSE (SUP) (DI) …`, `BARRE listeb <MASSE (DI), MASSE NUL>`,
  // `BARRE listeb (EXC …) MASSE (SUP) <UNI, LIN, PARABOL> …`) and SUP and DI
  // may both follow it, in that order. Real decks also elide it and write the
  // form keyword straight after the list (`NOEUD 2 DI 1000 1000 1000`,
  // `BARRE 1 NUL`); that elision is accepted here because the word can then
  // only be a mass form. Anything else stays the typed error it was.
  static const std::set<std::string> kMasseForms = {"SUP", "DI",  "NUL",
                                                    "UNI", "LIN", "PARABOL"};
  if (upper(st[i].s) == "MASSE") ++i;
  else
    require(kMasseForms.count(upper(st[i].s)) > 0,
            "ST1: MASSE PROPRE — mot '" + st[i].s + "' inattendu (MASSE attendu)");
  bool sup = false, di = false, nul = false;
  BarMassKind line_kind = BarMassKind::Auto;
  while (i < st.size()) {
    std::string k = upper(st[i].s);
    if (k == "SUP") { sup = true; ++i; continue; }
    if (k == "DI") { di = true; ++i; continue; }
    if (k == "NUL") { nul = true; ++i; continue; }
    if (k == "UNI") { line_kind = BarMassKind::Uni; ++i; continue; }
    if (k == "LIN") { line_kind = BarMassKind::Lin; ++i; continue; }
    if (k == "PARABOL") { line_kind = BarMassKind::Parabol; ++i; continue; }
    break;
  }
  std::vector<double> vals;
  for (; i < st.size(); ++i)
    if (is_num_tok(st[i])) vals.push_back(num(st[i]));

  if (!masse_pending_nodes_.empty()) {
    bool plane = option_ != 2;
    // `MASSE DI` takes 3 diagonal terms in PLANE / GRILL and 6 in SPATIALE,
    // the last three being the ROTARY inertias (p.70-71, whose own SPATIALE
    // example writes `1. 2. 3. 0. 0. 0.`). A deck that stops after the
    // translational masses means exactly that state — no rotary inertia — so
    // the tail is completed with zeros instead of rejecting the line. The
    // lower-triangle form is NOT padded: a partial triangle is ambiguous.
    if (di && !vals.empty() && vals.size() < (plane ? 3u : 6u)) {
      warn("ST1: MASSE PROPRE NOEUD — " + std::to_string(vals.size()) +
           " termes diagonaux au lieu de " + std::to_string(plane ? 3 : 6) +
           "; termes d'inertie de rotation completes a zero (manuel p.70-71)");
      vals.resize(plane ? 3u : 6u, 0.0);
    }
    Mat6 mm = di ? mass_from_diagonal(vals, plane) : mass_from_lower_triangle(vals, plane);
    if (option_ == 1) mm = to_grill_frame(mm);
    for (Id n : masse_pending_nodes_) {
      NodalMass nm;
      nm.node = n;
      nm.m = mm;
      nm.add = sup;  // ST1 SUP: accumulate instead of replacing (p.69)
      masse_.nodal.push_back(nm);
    }
    invalidate();
    return;
  }
  require(!masse_pending_bars_.empty(), "ST1: MASSE PROPRE — ni NOEUD ni BARRE");
  for (Id b : masse_pending_bars_) {
    if (nul) {
      BarMass bm;
      bm.bar = b;
      bm.kind = BarMassKind::Nul;
      masse_.bars.push_back(bm);
      continue;
    }
    if (line_kind == BarMassKind::Auto) {
      // `BARRE liste MASSE (DI)`: keep the automatic rho*A mass; DI selects the
      // diagonal form, MASSE alone the consistent one, bar by bar (both may
      // be mixed in one structure).
      if (di) masse_.lumped = true;
      else masse_.consistent = true;
      auto it = std::find_if(masse_.bars.begin(), masse_.bars.end(), [&](const BarMass& x) {
        return x.bar == b && x.kind == BarMassKind::Auto;
      });
      if (it == masse_.bars.end()) {
        BarMass bm;
        bm.bar = b;
        bm.kind = BarMassKind::Auto;
        masse_.bars.push_back(bm);
        it = std::prev(masse_.bars.end());
      }
      it->form = di ? BarMass::Form::Diagonal : BarMass::Form::Consistent;
      continue;
    }
    std::size_t need = line_kind == BarMassKind::Uni ? 1 : (line_kind == BarMassKind::Lin ? 2 : 3);
    require(vals.size() >= need, "ST1: MASSE PROPRE — masses lineiques manquantes");
    // Without SUP the entry REPLACES the line masses already declared on the
    // bar (manual p.69); SUP accumulates.
    if (!sup)
      std::erase_if(masse_.bars, [&](const BarMass& x) {
        return x.bar == b && x.kind != BarMassKind::Nul;
      });
    BarMass bm;
    bm.bar = b;
    bm.kind = line_kind;
    bm.m1 = vals[0];
    if (need > 1) bm.m2 = vals[1];
    if (need > 2) bm.m3 = vals[2];
    bm.ecc_or = masse_exc_or_;
    bm.ecc_ex = masse_exc_ex_;
    masse_.bars.push_back(bm);
  }
  invalidate();
}

// ----------------------------------------------------------------- CONTRG --
namespace {

// `(N cfn) (TY cfty) (TZ cftz) (MX cfmx) (MY cfmy) (MZ cfmz) (CST cst)` —
// manual p.99. Returns the index just past the coefficient run.
std::size_t parse_contrg_coeffs(const Statement& st, std::size_t i, ContrgCoeffs& c,
                                std::vector<std::string>& unknown) {
  while (i + 1 < st.size() && is_st1_number(st[i + 1].s) && !st[i + 1].str) {
    std::string k = upper(st[i].s);
    double v = st1_num(st[i + 1].s);
    if (k == "N") c.cfn = v;
    else if (k == "TY") c.cfty = v;
    else if (k == "TZ") c.cftz = v;
    else if (k == "MX") c.cfmx = v;
    else if (k == "MY") c.cfmy = v;
    else if (k == "MZ") c.cfmz = v;
    else if (k == "CST") c.cst = v;
    else {
      unknown.push_back(st[i].s);
      break;
    }
    i += 2;
  }
  return i;
}

}  // namespace

// `CONTRG i ('Titre') … FIN` (manual p.99-101). The coefficients are authored in
// ST1's listing convention (N/MY/MZ compression-positive, p.66), which is what
// `ForceConvention::St1` means in etude.h.
void St1Runtime::cmd_contrg(const Statement& st) {
  require(st.size() >= 2 && is_id_tok(st[1]), "ST1: CONTRG i (manuel p.99)");
  open_contrg_ = ContrgDef{};
  open_contrg_.id = static_cast<Id>(num(st[1]));
  open_contrg_.convention = ForceConvention::St1;
  for (std::size_t i = 2; i < st.size(); ++i)
    if (st[i].str) open_contrg_.title = st[i].s;
  contrg_kind_ = ContrgKind::Cte;
  contrg_pending_bars_.clear();
  contrg_or_ = contrg_mi_ = contrg_ex_ = ContrgCoeffs{};
  contrg_open_ = true;
  mode_ = Mode::Contrg;
}

void St1Runtime::contrg_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  auto flush = [&]() {
    if (contrg_pending_bars_.empty()) return;
    switch (contrg_kind_) {
      case ContrgKind::Cte:
        open_contrg_.set_cte(contrg_pending_bars_, contrg_or_);
        break;
      case ContrgKind::VarLin:
        open_contrg_.set_var_lin(contrg_pending_bars_, contrg_or_, contrg_ex_);
        break;
      case ContrgKind::VarPara:
        open_contrg_.set_var_para(contrg_pending_bars_, contrg_or_, contrg_mi_, contrg_ex_);
        break;
    }
    contrg_pending_bars_.clear();
    contrg_or_ = contrg_mi_ = contrg_ex_ = ContrgCoeffs{};
  };
  if (head == "FIN") {
    flush();
    contrgs[open_contrg_.id] = open_contrg_;
    contrg_open_ = false;
    mode_ = Mode::Top;
    return;
  }
  if (head == "CTE") {
    flush();
    contrg_kind_ = ContrgKind::Cte;
    return;
  }
  if (head == "VAR") {
    flush();
    require(st.size() >= 2, "ST1: CONTRG VAR <LIN, PARA>");
    std::string k = upper(st[1].s);
    require(k == "LIN" || k == "PARA", "ST1: CONTRG VAR attend LIN ou PARA (manuel p.99)");
    contrg_kind_ = k == "LIN" ? ContrgKind::VarLin : ContrgKind::VarPara;
    return;
  }
  std::vector<std::string> unknown;
  if (head == "OR" || head == "MI" || head == "EX") {
    require(!contrg_pending_bars_.empty(), "ST1: CONTRG — ligne OR/MI/EX sans liste de barres");
    ContrgCoeffs& slot = head == "OR" ? contrg_or_ : (head == "MI" ? contrg_mi_ : contrg_ex_);
    parse_contrg_coeffs(st, 1, slot, unknown);
    for (const std::string& u : unknown) warn("ST1: CONTRG — cle '" + u + "' ignoree");
    return;
  }
  if (is_id_tok(st[0]) || head == "TOUT" || head == "TOUTES") {
    flush();
    std::size_t i = 0;
    contrg_pending_bars_ = id_list(st, i, all_bar_ids());
    require(!contrg_pending_bars_.empty(), "ST1: CONTRG — liste de barres vide");
    // `liste OR …` (VAR forms) or `liste N cfn …` (CTE form) on the same line.
    if (i < st.size() && !st[i].str) {
      std::string k = upper(st[i].s);
      if (k == "OR" || k == "MI" || k == "EX") {
        ContrgCoeffs& slot = k == "OR" ? contrg_or_ : (k == "MI" ? contrg_mi_ : contrg_ex_);
        parse_contrg_coeffs(st, i + 1, slot, unknown);
      } else {
        parse_contrg_coeffs(st, i, contrg_or_, unknown);
      }
    }
    for (const std::string& u : unknown) warn("ST1: CONTRG — cle '" + u + "' ignoree");
    return;
  }
  flush();
  contrgs[open_contrg_.id] = open_contrg_;
  contrg_open_ = false;
  mode_ = Mode::Top;
  execute(st);
}

// -------------------------------------------------------------------- ENV --
// `ENV i (COMB) ('Titre') … FIN` (manual p.231 / p.234).
void St1Runtime::cmd_env(const Statement& st) {
  require(st.size() >= 2 && is_id_tok(st[1]), "ST1: ENV i (manuel p.231)");
  open_env_ = EnvDef{};
  open_env_.id = static_cast<Id>(num(st[1]));
  std::size_t i = 2;
  if (i < st.size() && !st[i].str && upper(st[i].s) == "COMB") {
    open_env_.comb = true;
    ++i;
  }
  for (; i < st.size(); ++i)
    if (st[i].str) open_env_.name = st[i].s;
  if (open_env_.name.empty()) open_env_.name = "ENV " + std::to_string(open_env_.id);
  env_cur_kind_.clear();
  env_open_ = true;
  mode_ = Mode::Env;
}

void St1Runtime::env_stmt(const Statement& st) {
  std::string head = st[0].str ? "" : upper(st[0].s);
  if (head == "FIN") {
    envs[open_env_.id] = open_env_;
    env_open_ = false;
    mode_ = Mode::Top;
    compute_env(envs.at(open_env_.id));
    return;
  }
  // --- privileged component (p.231) ---
  static const std::vector<std::string> kEffort = {"N", "TY", "TZ", "MX", "MY", "MZ"};
  static const std::vector<std::string> kDepla = {"DX", "DY", "DZ", "RX", "RY", "RZ"};
  static const std::vector<std::string> kPress = {"PX", "PY", "PZ", "MX"};
  static const std::vector<std::string> kFiber = {"VY",   "WY",   "VZ",   "WZ",
                                                  "VYVZ", "VYWZ", "WYVZ", "WYWZ"};
  auto index_of = [](const std::vector<std::string>& v, const std::string& s) {
    for (std::size_t k = 0; k < v.size(); ++k)
      if (v[k] == s) return static_cast<int>(k);
    return -1;
  };
  if (head == "EFFORT" || head == "DEPLA" || head == "PRESS" || head == "CONTR" ||
      head == "CONTRG") {
    require(st.size() >= 2, "ST1: ENV " + head + " — composante attendue");
    std::string c = upper(st[1].s);
    open_env_.has_priv = true;
    if (head == "EFFORT") {
      open_env_.what = EnvWhat::Effort;
      open_env_.comp = index_of(kEffort, c);
    } else if (head == "DEPLA") {
      open_env_.what = EnvWhat::Depla;
      open_env_.comp = index_of(kDepla, c);
    } else if (head == "PRESS") {
      open_env_.what = EnvWhat::Press;
      open_env_.comp = index_of(kPress, c);
    } else if (head == "CONTR") {
      open_env_.what = EnvWhat::Contr;
      int f = index_of(kFiber, c);
      require(f >= 0, "ST1: ENV CONTR — fibre '" + c + "' inconnue (manuel p.231)");
      open_env_.fiber = static_cast<ContrFiber>(f);
      open_env_.comp = 0;
    } else {
      open_env_.what = EnvWhat::Contrg;
      require(is_id_tok(st[1]), "ST1: ENV CONTRG ic — numero attendu");
      open_env_.contrg = static_cast<Id>(num(st[1]));
      open_env_.comp = 0;
    }
    require(open_env_.what == EnvWhat::Contr || open_env_.what == EnvWhat::Contrg ||
                open_env_.comp >= 0,
            "ST1: ENV " + head + " — composante '" + c + "' inconnue (manuel p.231)");
    return;
  }
  if (head == "REAC") {
    require(st.size() >= 2, "ST1: ENV REAC — composante attendue");
    static const std::vector<std::string> kReac = {"FX", "FY", "FZ", "MX", "MY", "MZ"};
    std::string c = upper(st[1].s);
    int k = index_of(kReac, c);
    require(k >= 0, "ST1: ENV REAC — composante '" + c + "' inconnue");
    open_env_.has_reac = true;
    open_env_.reac_comp = k;
    std::size_t i = 2;
    std::vector<Id> ids = id_list(st, i, all_appui_ids());
    if (!ids.empty()) open_env_.reac_appui = ids.front();
    return;
  }
  // --- source lines ---
  static const std::set<std::string> kKinds = {"CHARG",   "SURCH",   "COMB",       "ETAT",
                                               "PHASAGE", "SPECTRE", "HISTORIQUE", "INCENDIE",
                                               "ENV"};
  std::size_t i = 1;
  if (!kKinds.count(head)) {
    // Sub-list form: `CHARG` alone opens a run of `id (coef)` lines.
    require(!env_cur_kind_.empty() && (is_uid_tok(st[0]) || head == "TOUT"),
            "ST1: ENV — source '" + st[0].s + "' inconnue (manuel p.231)");
    head = env_cur_kind_;
    i = 0;
  } else {
    env_cur_kind_ = head;
  }
  std::vector<Id> ids = id_list(st, i, {}, /*comma_only=*/true);
  if (ids.empty()) return;  // `CHARG` alone: the ids come on the following lines
  double c1 = 1.0, c2 = 1.0;
  bool has1 = false;
  if (i < st.size() && is_num_tok(st[i])) {
    c1 = num(st[i++]);
    has1 = true;
  }
  if (i < st.size() && is_num_tok(st[i])) c2 = num(st[i++]);
  else c2 = c1;  // "cfmax = cfmin par defaut" (p.234)
  if (!has1) c1 = c2 = 1.0;
  for (Id id : ids) {
    EnvDef::Source s;
    s.kind = head;
    s.id = id;
    s.cmin = std::min(c1, c2);
    s.cmax = std::max(c1, c2);
    open_env_.sources.push_back(s);
  }
}

// ---------------------------------------------------------------- ACTIVER --
void St1Runtime::cmd_activer(const Statement& st, bool activate) {
  std::size_t i = 1;
  if (i < st.size() && (upper(st[i].s) == "TOUT" || upper(st[i].s) == "TOUTES")) {
    if (activate) {
      for (Id b : bar_order_) active_bars_.insert(b);
      for (Id a : appui_order_) active_appuis_.insert(a);
      selective_activation_ = true;  // explicit full activation
    } else {
      active_bars_.clear();
      active_appuis_.clear();
      selective_activation_ = true;
    }
    invalidate();
    return;
  }
  // Leaving the default state (everything active): "L'utilisation des mots
  // cles ACTIVER ou DESACTIVER supprime toute activation par defaut" (p.99-
  // 100) — the selective state starts EMPTY: `ACTIVER BARRE 1 a 5` alone
  // keeps those bars only (and no support), a first DESACTIVER leaves
  // nothing active.
  auto go_selective = [&] {
    if (selective_activation_) return;
    active_bars_.clear();
    active_appuis_.clear();
    selective_activation_ = true;
  };
  while (i < st.size()) {
    std::string k = upper(st[i].s);
    if (k == "BARRE") {
      ++i;
      std::vector<Id> ids = id_list(st, i, [&] {
        std::vector<Id> v(bar_order_.begin(), bar_order_.end());
        return v;
      }());
      go_selective();
      for (Id b : ids)
        if (activate) active_bars_.insert(b);
        else active_bars_.erase(b);
    } else if (k == "APPUI") {
      ++i;
      std::vector<Id> ids = id_list(st, i, [&] {
        std::vector<Id> v(appui_order_.begin(), appui_order_.end());
        return v;
      }());
      go_selective();
      for (Id a : ids)
        if (activate) active_appuis_.insert(a);
        else active_appuis_.erase(a);
    } else {
      warn("ST1: ACTIVER/DESACTIVER — mot '" + st[i].s + "' ignore");
      ++i;
    }
  }
  invalidate();
}

// ------------------------------------------------------------------- EXEC --
void St1Runtime::cmd_exec(const Statement& st) {
  std::string what = st.size() > 1 ? upper(st[1].s) : "CHARG";
  if (what == "CHARG") {
    exec_charg(st);
  } else if (what == "MODE") {
    exec_mode(st);
  } else if (what == "FLAMBEMENT") {
    exec_flambement(st);
  } else if (domain_handler && domain_handler(st)) {
    // EXEC SURCH / PHASAGE / SPECTRE REPONSE / HISTORIQUE / INCENDIE.
  } else {
    warn("ST1: EXEC " + what + " — integration moteur en cours");
  }
  mode_ = Mode::Top;
}

void St1Runtime::exec_charg(const Statement& st) {
  exec_charg_requested = true;
  std::size_t i = 2;
  bool pse = false;
  if (i < st.size() && !st[i].str && upper(st[i].s) == "PSE") {
    pse = true;
    ++i;
  }
  std::vector<Id> all;
  for (Id c : charg_order_)
    if (!chargs_.at(c).is_comb) all.push_back(c);
  std::vector<Id> requested = id_list(st, i, all);
  if (requested.empty()) requested = all;
  // A combined loading (CHARG i COMB, p.125) is the superposition of its
  // parts: `EXEC CHARG i` solves the parts too, as auxiliary cases that are
  // dropped from the store afterwards unless they were already there.
  std::vector<Id> list;
  std::set<Id> aux;
  // Under PSE the soil makes the response nonlinear: a combined loading is
  // then ONE loading, the weighted sum of its parts' loads,
  // solved as such rather than superposed.
  std::vector<Id> pse_combs;
  std::function<void(Id, bool)> want = [&](Id cid, bool asked) {
    require(chargs_.count(cid),
            "ST1: EXEC CHARG — chargement " + std::to_string(cid) + " inconnu");
    const ChargDef& c = chargs_.at(cid);
    if (c.is_comb && pse && asked) {
      if (std::find(pse_combs.begin(), pse_combs.end(), cid) == pse_combs.end())
        pse_combs.push_back(cid);
      return;
    }
    if (c.is_comb) {
      for (const auto& part : c.parts) want(part.first, false);
      return;
    }
    if (std::find(list.begin(), list.end(), cid) != list.end()) {
      if (asked) aux.erase(cid);
      return;
    }
    list.push_back(cid);
    if (!asked && !store_.charg.count(cid)) aux.insert(cid);
  };
  for (Id cid : requested) want(cid, true);

  build_model();
  // Assemble the requested cases into the snapshot and solve them together
  // (multi-RHS on one factorization).
  model_.cases.clear();
  for (Id cid : list) {
    require(chargs_.count(cid),
            "ST1: EXEC CHARG — chargement " + std::to_string(cid) + " inconnu");
    const ChargDef& c = chargs_.at(cid);
    LoadCase lc;
    lc.id = c.id;
    lc.name = c.name;
    lc.loads = lower_loads(c.loads);
    model_.cases.push_back(lc);
  }
  for (Id cid : pse_combs) {
    const ChargDef& c = chargs_.at(cid);
    LoadCase lc;
    lc.id = c.id;
    lc.name = c.name;
    std::function<void(Id, double)> gather = [&](Id pid, double coef) {
      require(chargs_.count(pid), "ST1: CHARG COMB " + std::to_string(cid) + " — chargement " +
                                      std::to_string(pid) + " inconnu");
      const ChargDef& p = chargs_.at(pid);
      if (p.is_comb) {
        for (const auto& [q, cq] : p.parts) gather(q, coef * cq);
        return;
      }
      for (const Load& l : lower_loads(p.loads)) lc.loads.push_back(scaled_load(l, coef));
    };
    gather(cid, 1.0);
    model_.cases.push_back(lc);
  }
  model_.combos.clear();
  // ST1 ETUDE (manual p.97-98): the study sections drive the recovery grid.
  StaticOptions sopt = recovery_options();
  std::string how;
  if (!cara_special_.empty()) {
    exec_charg_special(list);
    how = " (elements speciaux, secante)";
  } else if (!soil_zones_.empty()) {
    // `CARA PSE` bars carry the exact bedded block in model_ (bef.h). `EXEC
    // CHARG PSE` on a model with a unilateral / capped law runs the
    // pseudo-plastic active set (manual p.121) on the Hermite zoned bed
    // instead (the exact blocks removed so the soil is counted once); bars
    // that could not take the exact block always use the zoned bed.
    auto nonlinear = [](const BarSoil& bs) {
      if (!bs.plastic) return false;
      const SoilPlastic& p = *bs.plastic;
      return !p.bilateral() || p.ffx > 0 || p.ffy > 0 || p.ffz > 0 || p.fmx > 0;
    };
    bool any_nl = false;
    SoilZoneMap fallback;
    for (const auto& [bid, bs] : soil_zones_) {
      if (nonlinear(bs)) any_nl = true;
      if (!bef_bars_.count(bid)) fallback[bid] = bs;
    }
    // ST1 2.24 stops under EXEC CHARG PSE on a bar bedded on more than one
    // soil layer ("barre avec ressort non lineaire et plus d une couche de sol
    // non autorise"; the CARA PSE notes: one layer per bar on a
    // pseudo-plastic soil). keta computes it; say that ST1 would not.
    if (pse)
      for (const auto& [bid, bs] : soil_zones_)
        if (bs.zones.size() > 1)
          warn("ST1: EXEC CHARG PSE — barre " + std::to_string(bid) + " sur " +
               std::to_string(bs.zones.size()) + " couches de sol : ST1 2.24 refuse ce cas "
               "(une seule couche par barre en pseudo-plastique), keta le calcule");
    PseOptions po;
    po.sub_per_zone = 1;  // ST1 decides on the mean pressure of the bar (p.60)
    if (pse && any_nl) {
      Model mm = model_;
      for (Id b : bef_bars_) mm.tapered.erase(b);
      PseResults pr = solve_static_pse(mm, soil_zones_, sopt, po);
      for (const PseCaseResult& pc : pr.cases) {
        store_.charg[pc.result.id] = pc.result;
        pse_soil[pc.result.id] = pc.soil;
        stamp_cara("CHARG", pc.result.id);
      }
      how = " sur sol pseudo-plastique";
    } else if (!fallback.empty()) {
      PseResults pr = solve_static_pse(model_, fallback, sopt, po);
      for (const PseCaseResult& pc : pr.cases) {
        store_.charg[pc.result.id] = pc.result;
        pse_soil[pc.result.id] = pc.soil;
        stamp_cara("CHARG", pc.result.id);
      }
      how = " sur sol elastique";
    } else {
      StaticResults r = solve_static(model_, sopt);
      for (const CaseResult& cr : r.cases) {
        store_.charg[cr.id] = cr;
        pse_soil.erase(cr.id);
        stamp_cara("CHARG", cr.id);
      }
      how = " sur sol elastique (element exact)";
    }
  } else {
    if (pse)
      warn("ST1: EXEC CHARG PSE — aucune barre CARA PSE definie, calcul elastique classique");
    StaticResults r = solve_static(model_, sopt);
    for (const CaseResult& cr : r.cases) {
      store_.charg[cr.id] = cr;
      stamp_cara("CHARG", cr.id);
    }
  }

  // Combined loadings (CHARG COMB): linear superposition of the parts —
  // identical to a joint solve in linear statics (documented; the PSE
  // nonlinear path re-solves as one case at integration).
  for (Id cid : charg_order_) {
    const ChargDef& c = chargs_.at(cid);
    if (!c.is_comb) continue;
    if (std::count(pse_combs.begin(), pse_combs.end(), cid)) continue;  // solved as one
    std::vector<const CaseResult*> srcs;
    std::vector<double> coefs;
    bool ok = true;
    for (auto& [pid, coef] : c.parts) {
      auto it = store_.charg.find(pid);
      if (it == store_.charg.end()) {
        ok = false;
        break;
      }
      srcs.push_back(&it->second);
      coefs.push_back(coef);
    }
    if (!ok) {
      warn("ST1: CHARG COMB " + std::to_string(cid) + " — composantes non calculees");
      continue;
    }
    store_.charg[cid] = combine_results(srcs, coefs, cid, c.name);
    stamp_cara("CHARG", cid);
  }
  for (Id cid : aux) {
    store_.charg.erase(cid);
    pse_soil.erase(cid);
  }
  out_.line("EXEC CHARG : " + std::to_string(requested.size()) + " chargement(s) calcules" + how);
}

// ST1 `CARA SPECIAL` (manual p.210): a classic R.D.M. bar whose force
// components are capped. The pure per-element pieces live in tapered.h
// (`special_ratios` / `special_dof_scales` / `special_secant_stiffness`); the
// ITERATION is the runtime's job, and it is a SECANT fixed point: the member
// stiffness of a capped bar is replaced by S k S (S diagonal, a congruence, so
// symmetry and positive semi-definiteness survive) and the scales are
// accumulated until every capped component is within its threshold.
//
// Convergence is monitored on `special_overshoot`; a run that has not settled
// after `kMaxSweeps` is a loud warning with the residual overshoot, never a
// silently unconverged result. Loads stay elastic (the `TaperedBar::f0`
// callback is left empty on purpose — the prismatic span-load lowering is used).
void St1Runtime::exec_charg_special(const std::vector<Id>& cases) {
  constexpr int kMaxSweeps = 40;
  constexpr double kTol = 1e-6;
  std::map<Id, std::array<double, 12>> scales;
  std::map<Id, Mat12> elastic;
  std::unordered_set<Id> bar_ids;
  for (const Bar& b : model_.bars) bar_ids.insert(b.id);
  for (const auto& [bid, d] : cara_special_) {
    const bool present = bar_ids.count(bid) > 0;
    if (!present) continue;
    require(d.has_effort, "ST1: CARA SPECIAL " + std::to_string(bid) + " sans bloc EFFORT");
    const Bar& bar = model_.bar(bid);
    elastic[bid] = BarElement(model_, bar).k_local();
    scales[bid].fill(1.0);
  }
  StaticOptions sopt;
  sopt.stations_abs = etude_stations_ref();
  if (etude_pas_hint > 0) sopt.stations_per_bar = std::max(2, static_cast<int>(etude_pas_hint) + 1);

  double worst = 0;
  StaticResults r;
  for (int sweep = 0; sweep < kMaxSweeps; ++sweep) {
    model_.tapered.clear();
    for (const auto& [bid, s] : scales) {
      TaperedBar tb;
      tb.ke = special_secant_stiffness(elastic.at(bid), s);
      model_.tapered[bid] = std::move(tb);
    }
    r = solve_static(model_, sopt);
    worst = 0;
    bool changed = false;
    for (auto& [bid, s] : scales) {
      const SpecialLimits& lim = cara_special_.at(bid).lim;
      if (!lim.any()) continue;
      for (const CaseResult& cr : r.cases)
        for (const BarResult& br : cr.bars) {
          if (br.bar != bid) continue;
          worst = std::max(worst, special_overshoot(lim, br.end_forces_local));
          std::array<double, 12> sc = special_dof_scales(lim, br.end_forces_local);
          for (std::size_t d = 0; d < 12; ++d)
            if (sc[d] < 1.0 - 1e-15) {
              s[d] *= sc[d];
              changed = true;
            }
        }
    }
    if (worst <= kTol || !changed) break;
  }
  if (worst > kTol)
    warn("ST1: CARA SPECIAL — iteration secante non convergee (depassement relatif " +
         fmt(worst) + ", manuel p.210)");
  for (const CaseResult& cr : r.cases) {
    CaseResult out = cr;
    // Report the member end forces clamped onto their thresholds (p.210).
    for (BarResult& br : out.bars) {
      auto it = cara_special_.find(br.bar);
      if (it == cara_special_.end() || !it->second.lim.any()) continue;
      br.end_forces_local = special_clamp(it->second.lim, br.end_forces_local);
    }
    out.non_superposable = true;  // a capped element is not a linear system
    const Id oid = out.id;
    store_.charg[oid] = std::move(out);
    stamp_cara("CHARG", oid);
  }
  (void)cases;
  model_.tapered.clear();
  invalidate();
}

// `EXEC MODE (nb) (METHODE <1,2>) (<FREQUENCE fq, RATIO rt>) (AVANCE nm)`
// (manual p.194-195), driven through eig_ext.h: METHODE 1 is the deflated
// inverse iteration, METHODE 2 (the default) the subspace machinery, and the
// FREQUENCE / RATIO stopping criteria are the module's own — Sturm-counted for
// the cutoff, participation-cumulated for the ratio, both reporting a loud
// warning instead of silently returning an incomplete mode set.
void St1Runtime::exec_mode(const Statement& st) {
  exec_modal_requested = true;
  build_model();
  ModalExtOptions mo;
  mo.lumped_mass = masse_.diagonal();  // ST1 MASSE DI or the default form (p.69)
  mo.n_modes = 100;  // "par defaut au maximum 100 modes seront calcules" (p.198)
  std::size_t i = 2;
  double freq_cut = 0, ratio = 0;
  int avance = 10;
  while (i < st.size()) {
    std::string k = st[i].str ? "" : upper(st[i].s);
    if (is_num_tok(st[i]) && i == 2) mo.n_modes = static_cast<int>(num(st[i]));
    else if (k == "METHODE" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      mo.method = static_cast<int>(num(st[i + 1])) == 1 ? EigMethod::InverseIteration
                                                        : EigMethod::Subspace;
      ++i;
    } else if (k == "FREQUENCE" && i + 1 < st.size() && is_num_tok(st[i + 1])) freq_cut = num(st[++i]);
    else if (k == "RATIO" && i + 1 < st.size() && is_num_tok(st[i + 1])) ratio = num(st[++i]);
    else if (k == "AVANCE" && i + 1 < st.size() && is_num_tok(st[i + 1]))
      avance = static_cast<int>(num(st[++i]));
    ++i;
  }
  if (freq_cut > 0) {
    ModalFrequencyResults fr = solve_modal_frequency(model_, freq_cut, mo);
    if (!fr.warning.empty()) warn("ST1: EXEC MODE FREQUENCE — " + fr.warning);
    store_.modal = std::move(fr.modal);
  } else if (ratio > 0) {
    ModalRatioResults rr = solve_modal_ratio(model_, ratio, avance, mo);
    if (!rr.warning.empty()) warn("ST1: EXEC MODE RATIO — " + rr.warning);
    store_.modal = std::move(rr.modal);
  } else {
    store_.modal = solve_modal_method(model_, mo);
  }
  store_.has_modal = true;
  out_.line("EXEC MODE : " + std::to_string(store_.modal.frequencies_hz.size()) +
            " mode(s) calcules");
}

// `EXEC FLAMBEMENT (nb) (METHODE <1,2>) (CTE <CHARG,COMB,ETAT,PHASAGE> i)
//                  (VAR <CHARG,COMB,ETAT,PHASAGE> j)` (manual p.211-212).
//
// ST1 amplifies the VARIABLE axial state on top of a FIXED one:
//     [K - KG(Ni)] x = c [KG(Nj)] x                                   (p.212)
// which is exactly `solve_buckling_cte_var_fields` (eig_ext.h) on keta's
// tension-positive convention. Both states are read from the STORED results —
// the manual is explicit that they must have been computed beforehand — so any
// result kind (CHARG / COMB / ETAT / PHASAGE) works, not just a load case.
// Default mode count: 2 (p.211).
void St1Runtime::exec_flambement(const Statement& st) {
  exec_buckling_requested = true;
  build_model();
  int n_modes = 2;  // "par defaut 2 modes seront calcules" (p.211)
  BucklingCteVarOptions bo;
  std::size_t i = 2;
  Id var_case = 0, cte_case = 0;
  std::string var_kind = "CHARG", cte_kind = "CHARG";
  while (i < st.size()) {
    std::string k = st[i].str ? "" : upper(st[i].s);
    if (is_num_tok(st[i]) && i == 2) n_modes = std::max(1, static_cast<int>(num(st[i])));
    else if (k == "METHODE" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
      bo.method = static_cast<int>(num(st[i + 1])) == 1 ? EigMethod::InverseIteration
                                                        : EigMethod::Subspace;
      ++i;
    } else if ((k == "CTE" || k == "VAR") && i + 1 < st.size()) {
      std::string kind = "CHARG";
      std::size_t j = i + 1;
      if (!is_num_tok(st[j])) kind = upper(st[j++].s);
      require(j < st.size() && is_num_tok(st[j]),
              "ST1: EXEC FLAMBEMENT " + k + " — numero de resultat attendu");
      Id id = static_cast<Id>(num(st[j]));
      if (k == "CTE") { cte_kind = kind; cte_case = id; }
      else { var_kind = kind; var_case = id; }
      i = j;
    }
    ++i;
  }
  auto field_of = [&](const std::string& kind, Id id) -> std::vector<double> {
    if (id == 0) return {};  // "il suffit d'appliquer une charge 0" (p.212)
    const CaseResult* r = store_.find(kind, id);
    require(r, "ST1: EXEC FLAMBEMENT — " + kind + " " + std::to_string(id) +
                   " non calcule (les cas doivent avoir ete calcules au prealable, p.211)");
    return axial_field(*r);
  };
  std::vector<double> N_cte = field_of(cte_kind, cte_case);
  std::vector<double> N_var = field_of(var_kind, var_case);
  if (N_var.empty()) {
    // "Si on met une charge 0 sur la composante variable, ST1 considere un
    // effort normal uniforme de compression des barres" (p.212).
    bo.uniform_var_if_zero = true;
    N_var.assign(model_.bars.size(), 0.0);
  }
  // a VAR loading producing no axial force amplifies nothing: ST1 says
  // nothing and keeps the previous critical modes
  bool any_var = bo.uniform_var_if_zero;
  for (double v : N_var) any_var = any_var || v != 0.0;
  if (!any_var && store_.has_buckling) {
    warn("ST1: EXEC FLAMBEMENT — effort normal VAR nul : rien a amplifier, les modes "
         "critiques precedents sont conserves (comme ST1)");
    return;
  }
  BucklingCteVar b = solve_buckling_cte_var_fields(model_, N_cte, N_var, n_modes, bo);
  store_.buckling.load_factors = b.load_factors;
  store_.buckling.modes = b.modes;
  store_.buckling.residuals = b.residuals;
  store_.buckling.gen_stiffness = b.gen_stiffness;
  store_.buckling.participation = b.participation;
  store_.buckling.total_kg = b.total_kg;
  store_.has_buckling = true;
  out_.line("EXEC FLAMBEMENT : " + std::to_string(store_.buckling.load_factors.size()) +
            " mode(s) critiques");
}

// -------------------------------------------------------------------- GET --
// `GET <what> (MAX,MIN) …` (manual p.270-276). Every failure sets `$err = 1`
// and warns instead of aborting the run (p.270), which is what the ST1
// pseudo-programming layer tests with `SI ($err …)`.
void St1Runtime::cmd_get(const Statement& st) {
  env_.scalars["$err"] = 0.0;
  try {
    std::size_t i = 1;
    require(st.size() >= 3, "ST1: GET incomplet");
    std::string what = upper(st[i++].s);
    Id contrg_id = 0;
    if (what == "CONTRG") {
      require(i < st.size() && is_uid_tok(st[i]), "ST1: GET CONTRG k BARRE …");
      contrg_id = static_cast<Id>(num(st[i++]));
    }
    int side = 0;  // 0 = a single stored state, +1 = MAX, -1 = MIN
    if (i < st.size() && !st[i].str) {
      std::string k = upper(st[i].s);
      if (k == "MAX") { side = +1; ++i; }
      else if (k == "MIN") { side = -1; ++i; }
    }

    auto set6 = [&](const char* const names[6], const Vec6& v, double sign = 1.0) {
      for (int c = 0; c < 6; ++c) env_.scalars[names[c]] = sign * v(c);
    };
    static const char* kD[6] = {"$dx", "$dy", "$dz", "$rx", "$ry", "$rz"};
    static const char* kF[6] = {"$fx", "$fy", "$fz", "$mx", "$my", "$mz"};

    // The stored state a `<kind> j` reference denotes. An ENV or a SURCH is a
    // PAIR of states: MAX / MIN selects one of its two bounds (p.271-274).
    auto pick = [&](const std::string& kind, Id id) -> const CaseResult* {
      if (kind == "ENV") {
        auto it = env_full.find(id);
        require(it != env_full.end(), "ST1: GET — ENV " + std::to_string(id) + " non calculee");
        require(side != 0, "ST1: GET … ENV j — MAX ou MIN requis (manuel p.271)");
        return side < 0 ? &it->second.bounds.lo.value : &it->second.bounds.hi.value;
      }
      if (kind == "SURCH") {
        const Envelope& e = surch_bounds(id);
        require(side != 0, "ST1: GET … SURCH j — MAX ou MIN requis (manuel p.272)");
        return side < 0 ? &e.lo.value : &e.hi.value;
      }
      const CaseResult* r = store_.find(kind, id);
      require(r, "ST1: GET — resultat " + kind + " " + std::to_string(id) + " absent");
      return r;
    };

    if (what == "ALPHA") {
      require(i < st.size() && upper(st[i].s) == "FLAMBEMENT", "ST1: GET ALPHA FLAMBEMENT i");
      ++i;
      std::size_t idx = 1;
      if (i < st.size() && is_uid_tok(st[i])) idx = static_cast<std::size_t>(num(st[i]));
      require(store_.has_buckling && idx >= 1 && idx <= store_.buckling.load_factors.size(),
              "ST1: GET ALPHA — mode critique absent");
      env_.scalars["$alpha"] = store_.buckling.load_factors[idx - 1];
      return;
    }

    if (what == "TENSION") {
      // `GET TENSION CABLE i SEC n <CHARG,ETAT,PHASAGE> j` (manual p.276). GET
      // is a core command but the tendon results belong to the prestress
      // domain: the runtime only reads `cable_tension_store`, which the domain
      // extension fills. An empty store is a typed $err, never a silent zero.
      require(i < st.size() && upper(st[i].s) == "CABLE", "ST1: GET TENSION CABLE i SEC n …");
      ++i;
      require(i < st.size() && is_uid_tok(st[i]), "ST1: GET TENSION CABLE — numero attendu");
      Id cid = static_cast<Id>(num(st[i++]));
      require(i + 1 < st.size() && upper(st[i].s) == "SEC", "ST1: GET TENSION — SEC n requis");
      ++i;
      std::size_t sec = static_cast<std::size_t>(num(st[i++]));
      auto it = cable_tension_store.find(cid);
      require(it != cable_tension_store.end() && !it->second.empty(),
              "ST1: GET TENSION CABLE " + std::to_string(cid) +
                  " — aucune tension stockee (module precontrainte requis)");
      require(sec >= 1 && sec <= it->second.size(), "ST1: GET TENSION — section hors bornes");
      const std::array<double, 7>& row = it->second[sec - 1];
      env_.scalars["$nsec"] = static_cast<double>(it->second.size());
      env_.scalars["$xl"] = row[1];
      env_.scalars["$sigp"] = row[2];
      env_.scalars["$dsigp"] = row[3];
      env_.scalars["$x"] = row[4];
      env_.scalars["$y"] = row[5];
      env_.scalars["$z"] = row[6];
      return;
    }

    if (what == "DEPLA" && i < st.size() && upper(st[i].s) == "NOEUD") {
      ++i;
      Id nid = static_cast<Id>(num(st[i++]));
      std::string kind = upper(st[i++].s);
      Id cid = static_cast<Id>(num(st[i++]));
      const CaseResult* r = pick(kind, cid);
      std::size_t ni = model_node_index(nid);
      require(ni < r->node_disp.size(), "ST1: GET DEPLA — noeud hors du resultat");
      set6(kD, r->node_disp[ni]);
      return;
    }
    if (what == "REAC") {
      require(i < st.size() && upper(st[i].s) == "APPUI", "ST1: GET REAC APPUI i ...");
      ++i;
      Id aid = static_cast<Id>(num(st[i++]));
      std::string kind = upper(st[i++].s);
      Id cid = static_cast<Id>(num(st[i++]));
      const CaseResult* r = pick(kind, cid);
      Id nid = appui_.count(aid) ? appui_.at(aid).node : aid;
      for (const Reaction& re : r->reactions)
        if (re.node == nid) {
          set6(kF, re.R);
          return;
        }
      fail("ST1: GET REAC — appui " + std::to_string(aid) + " sans reaction");
    }
    if (what == "EFFORT" || what == "DEPLA" || what == "CONTR" || what == "PRESS" ||
        what == "CONTRG") {
      require(i < st.size() && upper(st[i].s) == "BARRE", "ST1: GET ... BARRE i SEC n ...");
      ++i;
      Id bid = static_cast<Id>(num(st[i++]));
      require(i + 1 < st.size() && upper(st[i].s) == "SEC", "ST1: GET — SEC n requis");
      ++i;
      std::size_t sec = static_cast<std::size_t>(num(st[i++]));
      std::string kind = upper(st[i++].s);
      Id cid = static_cast<Id>(num(st[i++]));
      const CaseResult* r = pick(kind, cid);
      // the result lists the bars in model order: read the bar's own slot
      // when it holds that bar, else the scan (the last match)
      const BarResult* br = nullptr;
      if (const int bx = model_has_bar(model_, bid); bx >= 0 && static_cast<std::size_t>(bx) < r->bars.size() &&
                                                   r->bars[static_cast<std::size_t>(bx)].bar == bid)
        br = &r->bars[static_cast<std::size_t>(bx)];
      else
        for (const BarResult& b : r->bars)
          if (b.bar == bid) br = &b;
      require(br, "ST1: GET — barre " + std::to_string(bid) + " absente du resultat");
      require(sec >= 1 && sec <= br->stations.size(), "ST1: GET — section hors bornes");
      Station s = br->stations[sec - 1];
      // An ENV declaring a privileged component carries the CONCOMITANT state
      // of its extreme at every study section (p.232): prefer it over the
      // field-by-field bound.
      if (kind == "ENV") {
        auto it = env_full.find(cid);
        if (it != env_full.end() && it->second.has_priv)
          for (const EnvPoint& p : it->second.points)
            if (p.bar == bid && p.station == sec - 1) s = side < 0 ? p.st_min : p.st_max;
      }
      env_.scalars["$xl"] = s.x;
      env_.scalars["$nsec"] = static_cast<double>(br->stations.size());
      if (what == "EFFORT") {
        // ST1 lists N positive in compression (manual p.60): flip keta's
        // tension-positive N at this boundary only.
        env_.scalars["$n"] = -s.N;
        env_.scalars["$ty"] = s.Vy;
        env_.scalars["$tz"] = s.Vz;
        env_.scalars["$mx"] = s.T;
        env_.scalars["$my"] = s.My;
        env_.scalars["$mz"] = s.Mz;
      } else if (what == "DEPLA") {
        env_.scalars["$dx"] = s.defl.x();
        env_.scalars["$dy"] = s.defl.y();
        env_.scalars["$dz"] = s.defl.z();
        env_.scalars["$rx"] = s.rot.x();
        env_.scalars["$ry"] = s.rot.y();
        env_.scalars["$rz"] = s.rot.z();
      } else if (what == "CONTR") {
        // Fibre stresses via the section fibre ordinates (compression > 0).
        require(cara_.count(bid), "ST1: GET CONTR — CARA absent");
        const CaraProps& cp = cara_.at(bid);
        double A = cp.sx, iy = cp.iy, iz = cp.iz;
        require(A > 0, "ST1: GET CONTR — SX nul");
        auto sigma = [&](double y, double z) {
          double sg = -s.N / A;  // compression-positive
          if (iy > 0) sg += -s.My * z / iy;
          if (iz > 0) sg += s.Mz * y / iz;
          return sg;
        };
        env_.scalars["$vy"] = sigma(cp.vy, 0);
        env_.scalars["$wy"] = sigma(-cp.wy, 0);
        env_.scalars["$vz"] = sigma(0, cp.vz);
        env_.scalars["$wz"] = sigma(0, -cp.wz);
        env_.scalars["$vyvz"] = sigma(cp.vy, cp.vz);
        env_.scalars["$vywz"] = sigma(cp.vy, -cp.wz);
        env_.scalars["$wyvz"] = sigma(-cp.wy, cp.vz);
        env_.scalars["$wywz"] = sigma(-cp.wy, -cp.wz);
      } else if (what == "CONTRG") {
        auto it = contrgs.find(contrg_id);
        require(it != contrgs.end(),
                "ST1: GET CONTRG " + std::to_string(contrg_id) + " — non definie");
        require(it->second.defined_for(bid),
                "ST1: GET CONTRG — non definie sur la barre " + std::to_string(bid) +
                    " (manuel p.101)");
        double L = bar_axis_length(model_, model_.bar(bid));
        // The constant term applies to loadings, not to surcharges (p.99).
        bool cst = kind != "SURCH";
        env_.scalars["$contrg"] = contrg_at(it->second, bid, L > 0 ? s.x / L : 0.0, s, cst);
      } else {  // PRESS — soil pressures of a bar on elastic ground (p.274)
        Vec3 p = Vec3::Zero();
        double mx = 0;
        bool found = false;
        if (sec - 1 < br->soil_p.size()) {
          p = br->soil_p[sec - 1];
          mx = br->soil_m[sec - 1];
          found = true;
        } else {
          auto it = pse_soil.find(cid);
          if (it != pse_soil.end())
            for (const PseBarSoil& bs : it->second)
              if (bs.bar == bid && sec - 1 < bs.soil_p.size()) {
                p = bs.soil_p[sec - 1];
                mx = bs.soil_m[sec - 1];
                found = true;
              }
        }
        require(found, "ST1: GET PRESS — barre " + std::to_string(bid) +
                           " sans sol elastique (CARA PSE requis)");
        env_.scalars["$px"] = p.x();
        env_.scalars["$py"] = p.y();
        env_.scalars["$pz"] = p.z();
        env_.scalars["$mx"] = mx;
      }
      return;
    }
    fail("ST1: GET '" + what + "' non reconnu");
  } catch (const Error& e) {
    env_.scalars["$err"] = 1.0;
    warn(std::string("ST1: GET en erreur — ") + e.what());
  }
}

std::size_t St1Runtime::model_node_index(Id nid) const {
  for (std::size_t k = 0; k < model_.nodes.size(); ++k)
    if (model_.nodes[k].id == nid) return k;
  fail("ST1: noeud " + std::to_string(nid) + " inconnu");
}

// ----------------------------------------------------------------- output --
void St1Runtime::cmd_sortie(const Statement& st) {
  require(st.size() >= 2, "ST1: SORTIE incomplet");
  if (st[1].str) {
    bool suite = false;
    for (std::size_t i = 2; i < st.size(); ++i) {
      std::string k = upper(st[i].s);
      if (k == "SUITE") suite = true;
      // NON PAGINE: pagination is a print-layout concern — recorded, no-op.
    }
    out_.file = std::make_unique<std::ofstream>(
        st[1].s, suite ? std::ios::app : std::ios::trunc);
    require(out_.file->good(), "ST1: SORTIE — impossible d'ouvrir '" + st[1].s + "'");
    out_.last_file_path = st[1].s;
    out_.to_file = true;
    return;
  }
  std::string k = upper(st[1].s);
  if (k == "CONSOLE") out_.to_file = false;
  else if (k == "FICHIER") {
    require(!out_.last_file_path.empty(), "ST1: SORTIE FICHIER sans fichier prealable");
    out_.file = std::make_unique<std::ofstream>(out_.last_file_path, std::ios::app);
    out_.to_file = true;
  } else warn("ST1: SORTIE — forme non reconnue");
}

void St1Runtime::cmd_projet(const Statement& st, bool save) {
  // ENREGISTRER (PROJET) 'nom' — writes a summary of the model (see
  // save_model_json), not a base that PROJET could read back.
  std::string name;
  for (std::size_t i = 1; i < st.size(); ++i)
    if (st[i].str) name = st[i].s;
  if (save) {
    require(!name.empty(), "ST1: ENREGISTRER — nom de projet requis");
    build_model();
    save_model_json(name);
    out_.line("PROJET '" + name + "' enregistre (" + name + ".json)");
  } else {
    if (name.empty()) {
      out_.line("PROJET COURANT : " + (project_.empty() ? "(aucun)" : project_));
      return;
    }
    warn("ST1: PROJET '" + name + "' — relecture de base non tenue (recharger le .st1)");
  }
  project_ = name;
}

void St1Runtime::save_model_json(const std::string& name) {
  // A summary only (project name, node and bar counts). The model JSON of
  // keta/io.h does not carry what a deck can build — bars with their own
  // kernel (variable sections, capped elements), soil zones, the phasing —
  // so writing it here would save a model that differs from the deck.
  std::ofstream f(name + ".json");
  require(f.good(), "ST1: ENREGISTRER — impossible d'ecrire '" + name + ".json'");
  f << "{ \"projet\": \"" << name << "\", \"noeuds\": " << nodes_.size()
    << ", \"barres\": " << bars_.size() << " }\n";
}

void St1Runtime::cmd_lister(const Statement& st) {
  std::string what = st.size() > 1 && !st[1].str ? upper(st[1].s) : "GEOM";
  auto& o = out_;
  // Optional id list after the family keyword (manual p.236).
  std::size_t li = 2;
  std::vector<Id> want = st.size() > 2 ? id_list(st, li, {}) : std::vector<Id>{};
  auto wanted = [&](Id id) { return want.empty() || std::count(want.begin(), want.end(), id); };
  if (what == "PROJET") {
    o.line("LISTER PROJET : " + (project_.empty() ? "(aucun)" : project_));
    return;
  }
  if (what == "GEOM" || what == "NOEUD") {
    o.line("NOEUDS (" + std::to_string(nodes_.size()) + ")");
    for (Id n : node_order_) {
      if (what == "NOEUD" && !wanted(n)) continue;
      const Vec3& c = nodes_.at(n);
      o.line("  " + std::to_string(n) + "  " + fmt(c.x()) + "  " + fmt(c.y()) + "  " +
             fmt(c.z()));
    }
    if (what == "NOEUD") return;
  }
  if (what == "GEOM" || what == "BARRE") {
    o.line("BARRES (" + std::to_string(bars_.size()) + ")");
    for (Id b : bar_order_) {
      if (what == "BARRE" && !wanted(b)) continue;
      const RawBar& rb = bars_.at(b);
      o.line("  " + std::to_string(b) + "  DE " + std::to_string(rb.n1) + " A " +
             std::to_string(rb.n2));
    }
    if (what == "BARRE") return;
  }
  if (what == "GEOM" || what == "APPUI") {
    o.line("APPUIS (" + std::to_string(appui_.size()) + ")");
    for (Id a : appui_order_) {
      if (what == "APPUI" && !wanted(a)) continue;
      const SupportDef& s = appui_.at(a);
      std::string ddl;
      const char* names[6] = {"DX", "DY", "DZ", "RX", "RY", "RZ"};
      for (int c = 0; c < 6; ++c)
        if (s.fix[static_cast<std::size_t>(c)]) ddl += std::string(" ") + names[c];
      if (!s.spring_lt.empty()) ddl += s.spring_di ? " EL DI" : " EL";
      for (int c = 0; c < 6; ++c)
        if (s.spring[static_cast<std::size_t>(c)] != 0)
          ddl += " K" + std::string(names[c]) + " " + fmt(s.spring[static_cast<std::size_t>(c)]);
      if (s.rot != 0) ddl += " ROT " + fmt(s.rot);
      if (s.decol != 0) ddl += std::string(" DECOL") + (s.decol > 0 ? " POS" : " NEG");
      o.line("  " + std::to_string(a) + " NOEUD " + std::to_string(s.node) + ddl);
    }
    return;
  }
  if (what == "CARA") {
    for (auto& [b, p] : cara_) {
      if (!wanted(b)) continue;
      std::string line = "CARA " + std::to_string(b) + " SX " + fmt(p.sx) + " IY " + fmt(p.iy) +
                         " IZ " + fmt(p.iz) + " IX " + fmt(p.ix) + " SY " + fmt(p.sy) + " SZ " +
                         fmt(p.sz);
      auto v = cara_var_.find(b);
      if (v != cara_var_.end())
        line += std::string(" VAR ") + (v->second.law == TaperLaw::Para
                                            ? "PARA"
                                            : (v->second.law == TaperLaw::LinY ? "LIN Y" : "LIN Z"));
      auto sz = soil_zones_.find(b);
      if (sz != soil_zones_.end())
        line += " PSE " + std::to_string(sz->second.zones.size()) + " zone(s)" +
                (sz->second.plastic && sz->second.plastic->nonlinear() ? " PLASTIQUE" : "");
      if (cara_special_.count(b)) line += " SPECIAL";
      o.line(line);
      if (v != cara_var_.end()) {
        const CaraVarDef& d = v->second;
        o.line("    OR SX " + fmt(d.orig.sx) + " IY " + fmt(d.orig.iy) + " IZ " + fmt(d.orig.iz));
        if (d.has_mi)
          o.line("    MI SX " + fmt(d.mid.sx) + " IY " + fmt(d.mid.iy) + " IZ " + fmt(d.mid.iz));
        o.line("    EX SX " + fmt(d.ext.sx) + " IY " + fmt(d.ext.iy) + " IZ " + fmt(d.ext.iz));
      }
      if (sz != soil_zones_.end()) {
        for (const SoilZone& z : sz->second.zones)
          o.line("    ZONE KFX " + fmt(z.kfx) + " KFY " + fmt(z.kfy) + " KFZ " + fmt(z.kfz) +
                 " KMX " + fmt(z.kmx) + " XL " + fmt(z.xl_end));
        if (sz->second.plastic) {
          const SoilPlastic& pl = *sz->second.plastic;
          o.line(std::string("    PLASTIQUE") + (pl.pos ? " POS" : "") + (pl.neg ? " NEG" : "") +
                 " KFX " + fmt(pl.ffx) + " KFY " + fmt(pl.ffy) + " KFZ " + fmt(pl.ffz) + " KMX " +
                 fmt(pl.fmx));
        }
      }
      if (cara_special_.count(b)) {
        const SpecialLimits& lim = cara_special_.at(b).lim;
        o.line("    EFFORT " + fmt(lim.cap[0]) + " " + fmt(lim.cap[1]) + " " + fmt(lim.cap[2]) +
               " " + fmt(lim.cap[3]) + " " + fmt(lim.cap[4]) + " " + fmt(lim.cap[5]));
      }
    }
    return;
  }
  if (what == "CONTRG") {
    for (const auto& [id, def] : contrgs) {
      if (!wanted(id)) continue;
      o.line("CONTRG " + std::to_string(id) + " '" + def.title + "' (" +
             std::to_string(def.bars.size()) + " barre(s))");
      for (const auto& [bid, rule] : def.bars) {
        const char* kind = rule.kind == ContrgKind::Cte
                               ? "CTE"
                               : (rule.kind == ContrgKind::VarLin ? "VAR LIN" : "VAR PARA");
        o.line("  BARRE " + std::to_string(bid) + " " + kind + " N " + fmt(rule.or_.cfn) +
               " TY " + fmt(rule.or_.cfty) + " TZ " + fmt(rule.or_.cftz) + " MX " +
               fmt(rule.or_.cfmx) + " MY " + fmt(rule.or_.cfmy) + " MZ " + fmt(rule.or_.cfmz) +
               " CST " + fmt(rule.or_.cst));
      }
    }
    return;
  }
  if (what == "MASSE") {
    o.line("MASSES PROPRES");
    for (const NodalMass& nm : masse_.nodal) {
      std::string d;
      for (int c = 0; c < 6; ++c) d += " " + fmt(nm.m(c, c));
      o.line(std::string("  NOEUD ") + std::to_string(nm.node) + (nm.add ? " SUP" : "") +
             " DIAG" + d);
    }
    for (const BarMass& bm : masse_.bars) {
      if (bm.kind == BarMassKind::Auto) {
        o.line("  BARRE " + std::to_string(bm.bar) + " MASSE AUTOMATIQUE");
        continue;
      }
      const char* k = bm.kind == BarMassKind::Nul
                          ? "NUL"
                          : (bm.kind == BarMassKind::Uni
                                 ? "UNI"
                                 : (bm.kind == BarMassKind::Lin ? "LIN" : "PARABOL"));
      o.line("  BARRE " + std::to_string(bm.bar) + " MASSE " + k + " " + fmt(bm.m1) + " " +
             fmt(bm.m2) + " " + fmt(bm.m3));
    }
    o.line(std::string("  MASSES AUTOMATIQUES : ") + (masse_.diagonal() ? "DIAGONALES" : "CONSISTANTES"));
    return;
  }
  if (what == "ENV") {
    for (const auto& [id, def] : envs) {
      if (!wanted(id)) continue;
      o.line("ENV " + std::to_string(id) + (def.comb ? " COMB '" : " '") + def.name + "' (" +
             std::to_string(def.sources.size()) + " source(s))" +
             (env_full.count(id) ? " calculee" : " non calculee"));
    }
    return;
  }
  if (what == "CHARG" || what == "CHARGE") {
    for (Id c : charg_order_) {
      if (!wanted(c)) continue;
      o.line("CHARG " + std::to_string(c) + " '" + chargs_.at(c).name + "' (" +
             std::to_string(chargs_.at(c).loads.size()) + " charges)");
    }
    return;
  }
  if (what == "ETUDE") {
    const std::map<Id, std::vector<double>>* m = etude_stations_ref();
    if (!m) {
      o.line("ETUDE : sections par defaut (extremites de barre)");
      return;
    }
    for (const auto& [bid, xs] : *m) {
      if (!wanted(bid)) continue;
      std::string line = "ETUDE BARRE " + std::to_string(bid) + " SE";
      for (double x : xs) line += " " + fmt(x);
      o.line(line);
    }
    return;
  }
  if (what == "ACTIVATION") {
    o.line(selective_activation_ ? "ACTIVATION SELECTIVE" : "ACTIVATION PAR DEFAUT (TOUT)");
    std::string bl = "  BARRES ACTIVES :";
    for (Id b : all_bar_ids()) bl += " " + std::to_string(b);
    o.line(bl);
    std::string al = "  APPUIS ACTIFS  :";
    for (Id a : all_appui_ids()) al += " " + std::to_string(a);
    o.line(al);
    return;
  }
  warn("ST1: LISTER " + what + " — famille listee a l'integration");
}

// ------------------------------------------------------------------- RESU --
void St1Runtime::cmd_resu_line(const Statement& st) {
  std::string head = upper(st[0].s);
  if (head == "FIN") {
    // ST1 edition: a line naming a result kind opens a NEW request once the
    // current one carries table options of its own (Exemple6: `spectre appui`
    // then `spectre barre effort depla` print the support blocks, then the
    // bar blocks); kind lines without options (`comb tout` / `charg tout` /
    // ... / `barre tout effort`) keep sharing the options that follow.
    static const std::set<std::string> kKinds = {"CHARG", "CHARGE", "COMB", "ENV", "SURCH",
                                                 "ETAT", "PHASAGE", "SPECTRE", "INCENDIE",
                                                 "HISTORIQUE", "MODE", "FLAMBEMENT"};
    static const std::set<std::string> kTables = {"NOEUD", "APPUI", "BARRE", "EFFORT", "DEPLA",
                                                  "CONTR", "CONTRG", "PRESS", "CABLE", "FREQUENCE"};
    std::vector<std::vector<std::vector<std::string>>> groups(1);
    bool tables = false;
    for (const auto& line : resu_req_) {
      const bool kind_line = !line.empty() && kKinds.count(line[0]);
      if (opt_.st1_edition && kind_line && tables && !groups.back().empty()) {
        groups.emplace_back();
        tables = false;
      }
      groups.back().push_back(line);
      for (std::size_t k = kind_line ? 1 : 0; k < line.size(); ++k)
        if (kTables.count(line[k])) tables = true;
    }
    const auto all = resu_req_;
    for (const auto& g : groups) {
      resu_req_ = g;
      render_resu();
    }
    resu_req_ = all;
    resu_req_.clear();
    mode_ = Mode::Top;
    return;
  }
  // A `LISTER …` line inside a RESU block is executed on the spot (the real
  // decks do exactly that to interleave data recalls with the result edition).
  if (head == "LISTER") {
    cmd_lister(st);
    return;
  }
  // Accumulate request tokens (kinds, lists, component filters).
  std::vector<std::string> toks;
  for (const Tok& t : st) toks.push_back(t.str ? t.s : upper(t.s));
  resu_req_.push_back(toks);
}

namespace {

// A parsed RESU request: which stored results to edit, narrowed by the
// NOEUD / APPUI / BARRE lists and the component filters (manual p.240-241).
// The whole block is one request: ST1 explicitly allows every keyword on a
// single line, so the tokens are scanned as ONE keyword-driven stream.
struct ResuReq {
  std::map<std::string, std::vector<Id>> kinds;  // CHARG/COMB/ENV/SURCH/... -> ids (empty = all)
  std::vector<Id> nodes, appuis, bars;
  bool has_nodes = false, has_appuis = false, has_bars = false;
  bool appui_loc = false;  // APPUI LOC: reactions in each support's ROT frame
  // CABLE … XS (REL,ABS) listexs: the cable rows at these curvilinear
  // abscissae (REL: fractions of the cable length) instead of every point
  std::vector<double> xs;
  bool has_xs = false, xs_rel = false;
  bool effort = false, depla = false, contr = false, press = false;
  std::vector<Id> contrg;
  bool contrg_all = false;
  bool mode_freq = false, mode_depla = false;
  std::vector<Id> modes;
  bool any_comp() const { return effort || depla || contr || press || contrg_all || !contrg.empty(); }
};

}  // namespace

void St1Runtime::render_resu() {
  // model_ is not rebuilt while the results are edited: one id -> bar table
  // for the per-bar (and per-station) lookups of the editions
  const std::unordered_map<Id, const Bar*> model_bar = bars_by_id(model_);
  auto& o = out_;
  // ---- one keyword-driven pass over every token of the block ----
  ResuReq req;
  std::string cur;  // the keyword the following integers belong to
  static const std::set<std::string> kResultKinds = {
      "CHARG", "CHARGE", "COMB", "ENV", "SURCH", "ETAT", "PHASAGE", "SPECTRE",
      "INCENDIE", "HISTORIQUE", "MODE", "FLAMBEMENT", "CABLE"};
  auto canonical = [](std::string k) { return k == "CHARGE" ? "CHARG" : k; };
  for (const auto& line : resu_req_) {
    for (std::size_t k = 0; k < line.size(); ++k) {
      const std::string& t = line[k];
      if (t == "," ) continue;
      if (cur == "XS") {  // REL / ABS, reals, `a A b PAS p` ranges, commas
        if (t == "REL" || t == "ABS") {
          req.xs_rel = t == "REL";
          continue;
        }
        if (is_st1_number(t)) {
          req.xs.push_back(st1_num(t));
          continue;
        }
        if (t == "A" && k + 1 < line.size() && is_st1_number(line[k + 1]) && !req.xs.empty()) {
          const double a = req.xs.back(), b = st1_num(line[k + 1]);
          double step = b > a ? b - a : 0;
          std::size_t k2 = k + 1;
          if (k + 3 < line.size() && line[k + 2] == "PAS" && is_st1_number(line[k + 3])) {
            step = st1_num(line[k + 3]);
            k2 = k + 3;
          }
          if (step > 0)
            for (double x = a + step; x <= b + 1e-9 * (1 + std::fabs(b)); x += step)
              req.xs.push_back(x);
          k = k2;
          continue;
        }
      }
      if (is_integer(t) && t[0] != '-') {
        Id id = static_cast<Id>(std::stoll(t));
        if (cur == "NOEUD") req.nodes.push_back(id);
        else if (cur == "APPUI") req.appuis.push_back(id);
        else if (cur == "BARRE") req.bars.push_back(id);
        else if (cur == "CONTRG") req.contrg.push_back(id);
        else if (cur == "MODE" || cur == "FLAMBEMENT") req.modes.push_back(id);
        else if (!cur.empty() && kResultKinds.count(cur)) req.kinds[canonical(cur)].push_back(id);
        continue;
      }
      if (t == "LOC" && cur == "APPUI") { req.appui_loc = true; continue; }
      if (t == "XS" && cur == "CABLE") {
        cur = t;
        req.has_xs = true;
        continue;
      }
      if (t == "A" || t == "TOUT" || t == "TOUTES" || t == "PAS" || t == "LOC" ||
          t == "REL" || t == "ABS" || t == "XS")
        continue;  // range / scope words: the default is already "everything"
      if (kResultKinds.count(t)) {
        // `MODE` is NOT always a request for the modal edition: the RESU
        // grammar (manual p.240) nests it inside two other families —
        // `(FLAMBEMENT (listem) (MODE) (DEPLA))` and
        // `(SPECTRE (liste) (MODE (listem)) …)` — where it is a SUBOPTION of
        // the family already being read (p.241: "La commande FLAMBEMENT MODE
        // fournit les facteurs d'amplification critiques des modes de
        // flambement élastiques calculés"). Reading it as a new kind there
        // asked for modal results the deck never computed (`RESU MODE — EXEC
        // MODE prealable requis` on exemple8.st1). Keep routing the following
        // mode numbers to `req.modes`, but do not materialize a MODE request.
        if (t == "MODE" && (cur == "FLAMBEMENT" || cur == "SPECTRE")) {
          cur = t;
          continue;
        }
        cur = t;
        req.kinds[canonical(t)];  // materialize with an empty (= all) id list
        continue;
      }
      if (t == "NOEUD") { cur = t; req.has_nodes = true; continue; }
      if (t == "APPUI") { cur = t; req.has_appuis = true; continue; }
      if (t == "BARRE") { cur = t; req.has_bars = true; continue; }
      if (t == "EFFORT") { req.effort = true; continue; }
      if (t == "DEPLA") {
        if (cur == "MODE" || cur == "FLAMBEMENT") req.mode_depla = true;
        else req.depla = true;
        continue;
      }
      if (t == "CONTR") { req.contr = true; continue; }
      if (t == "PRESS") { req.press = true; continue; }
      if (t == "CONTRG") { cur = t; req.contrg_all = true; continue; }
      if (t == "FREQUENCE") { req.mode_freq = true; continue; }
      if (t == "VITESSE" || t == "ACCELERATION") continue;  // historique: domain families
      warn("ST1: RESU — mot '" + t + "' ignore");
    }
  }
  bool comps = req.any_comp();
  auto bar_wanted = [&](Id b) {
    return !req.has_bars || req.bars.empty() || std::count(req.bars.begin(), req.bars.end(), b);
  };
  auto node_wanted = [&](Id n) {
    return !req.has_nodes || req.nodes.empty() || std::count(req.nodes.begin(), req.nodes.end(), n);
  };
  // APPUI LOC (manual p.240): a reaction in the frame of its support, i.e.
  // turned by the APPUI's ROT alpha about Z (forces and moments alike).
  auto appui_frame = [&](Id aid, Vec6 R) {
    auto ad = appui_.find(aid);
    if (!req.appui_loc || ad == appui_.end() || ad->second.rot == 0) return R;
    const double c = std::cos(ad->second.rot), sn = std::sin(ad->second.rot);
    for (int k : {0, 3}) {
      const double x = R(k), y = R(k + 1);
      R(k) = c * x + sn * y;
      R(k + 1) = -sn * x + c * y;
    }
    return R;
  };
  auto appui_wanted = [&](Id n) {
    if (!req.has_appuis || req.appuis.empty()) return true;
    for (Id a : req.appuis) {
      Id nid = appui_.count(a) ? appui_.at(a).node : a;
      if (nid == n) return true;
    }
    return false;
  };

  auto fibre_line = [&](Id bid, const Station& s) {
    auto it = cara_.find(bid);
    if (it == cara_.end() || it->second.sx <= 0) return std::string();
    const CaraProps& cp = it->second;
    // ST1 lists normal stresses COMPRESSION-POSITIVE (manual p.66-68):
    // sigma = N/SX + MY*z/IY - MZ*y/IZ with N compression-positive.
    auto sigma = [&](double y, double z) {
      double sg = -s.N / cp.sx;
      if (cp.iy > 0) sg += -s.My * z / cp.iy;
      if (cp.iz > 0) sg += s.Mz * y / cp.iz;
      return sg;
    };
    if (option_ == 0)
      return "  VY " + fmt(sigma(cp.vy, 0)) + " WY " + fmt(sigma(-cp.wy, 0));
    if (option_ == 1)
      return "  VZ " + fmt(sigma(0, cp.vz)) + " WZ " + fmt(sigma(0, -cp.wz));
    return "  VYVZ " + fmt(sigma(cp.vy, cp.vz)) + " VYWZ " + fmt(sigma(cp.vy, -cp.wz)) +
           " WYVZ " + fmt(sigma(-cp.wy, cp.vz)) + " WYWZ " + fmt(sigma(-cp.wy, -cp.wz));
  };

  auto tabulate = [&](const CaseResult& r, const std::string& tag,
                      const std::vector<PseBarSoil>* soil) {
    o.line("== " + tag + " '" + r.name + "' ==");
    // `NOEUD` asks for the node displacements whatever filters the bar line
    // carries (exemple2: `noeud` + `barre effort` lists both, like ST1).
    if (!comps || req.depla || req.has_nodes) {
      o.line("  DEPLACEMENTS DES NOEUDS");
      for (std::size_t k = 0; k < r.node_disp.size(); ++k) {
        Id nid = model_.nodes[k].id;
        if (!node_wanted(nid)) continue;
        const Vec6& d = r.node_disp[k];
        o.line("    NOEUD " + std::to_string(nid) + "  " + fmt(d(0)) + " " + fmt(d(1)) + " " +
               fmt(d(2)) + " " + fmt(d(3)) + " " + fmt(d(4)) + " " + fmt(d(5)));
      }
    }
    if (!comps) {
      o.line("  REACTIONS D'APPUI");
      for (const Reaction& re : r.reactions) {
        if (!appui_wanted(re.node)) continue;
        o.line("    NOEUD " + std::to_string(re.node) + "  " + fmt(re.R(0)) + " " +
               fmt(re.R(1)) + " " + fmt(re.R(2)) + " " + fmt(re.R(3)) + " " + fmt(re.R(4)) +
               " " + fmt(re.R(5)));
      }
    } else if (req.has_appuis) {
      o.line("  REACTIONS D'APPUI");
      for (const Reaction& re : r.reactions) {
        if (!appui_wanted(re.node)) continue;
        o.line("    NOEUD " + std::to_string(re.node) + "  " + fmt(re.R(0)) + " " +
               fmt(re.R(1)) + " " + fmt(re.R(2)) + " " + fmt(re.R(3)) + " " + fmt(re.R(4)) +
               " " + fmt(re.R(5)));
      }
    }
    if (!comps || req.effort) {
      o.line("  EFFORTS (N positif en compression)");
      for (const BarResult& b : r.bars) {
        if (!bar_wanted(b.bar)) continue;
        for (const Station& s : b.stations)
          o.line("    BARRE " + std::to_string(b.bar) + " XL " + fmt(s.x) + "  N " + fmt(-s.N) +
                 " TY " + fmt(s.Vy) + " TZ " + fmt(s.Vz) + " MX " + fmt(s.T) + " MY " +
                 fmt(s.My) + " MZ " + fmt(s.Mz));
      }
    }
    if (comps && req.depla && req.has_bars) {
      // Global axes, as ST1 prints them (a section at x = 0 reads its node).
      o.line("  DEFORMEES DES SECTIONS D'ETUDE (repere global)");
      for (const BarResult& b : r.bars) {
        if (!bar_wanted(b.bar)) continue;
        // The result may come from an earlier activation state (pushing,
        // phasing): a bar absent from the current model reads its raw geometry.
        Mat3 Rt;
        if (auto it = model_bar.find(b.bar); it != model_bar.end()) {
          const Bar* mb = it->second;
          Rt = bar_rotation(model_.node(mb->ni).x, model_.node(mb->nj).x, mb->roll).transpose();
        } else {
          auto rb = bars_.find(b.bar);
          if (rb == bars_.end()) continue;
          Rt = bar_rotation(nodes_.at(rb->second.n1), nodes_.at(rb->second.n2), 0.0).transpose();
        }
        for (const Station& s : b.stations) {
          Vec3 d = Rt * s.defl, w = Rt * s.rot;
          o.line("    BARRE " + std::to_string(b.bar) + " XL " + fmt(s.x) + "  DX " + fmt(d.x()) +
                 " DY " + fmt(d.y()) + " DZ " + fmt(d.z()) + " RX " + fmt(w.x()) + " RY " +
                 fmt(w.y()) + " RZ " + fmt(w.z()));
        }
      }
    }
    if (req.contr) {
      o.line("  CONTRAINTES NORMALES (positives en compression)");
      for (const BarResult& b : r.bars) {
        if (!bar_wanted(b.bar)) continue;
        for (const Station& s : b.stations) {
          std::string f = fibre_line(b.bar, s);
          if (f.empty()) continue;
          o.line("    BARRE " + std::to_string(b.bar) + " XL " + fmt(s.x) + f);
        }
      }
    }
    if (req.contrg_all || !req.contrg.empty()) {
      for (const auto& [cid, def] : contrgs) {
        if (!req.contrg.empty() && !std::count(req.contrg.begin(), req.contrg.end(), cid))
          continue;
        o.line("  CONTRAINTE GENERALISEE " + std::to_string(cid) + " '" + def.title + "'");
        for (const BarResult& b : r.bars) {
          if (!bar_wanted(b.bar) || !def.defined_for(b.bar)) continue;
          double L = bar_axis_length(model_, model_.bar(b.bar));
          for (const Station& s : b.stations)
            o.line("    BARRE " + std::to_string(b.bar) + " XL " + fmt(s.x) + "  CONTRG " +
                   fmt(contrg_at(def, b.bar, L > 0 ? s.x / L : 0.0, s)));
        }
      }
    }
    if (req.press) {
      o.line("  PRESSIONS DU SOL");
      for (const BarResult& b : r.bars) {
        if (!bar_wanted(b.bar)) continue;
        const std::vector<Vec3>* p = &b.soil_p;
        const std::vector<double>* mm = &b.soil_m;
        if (p->empty() && soil)
          for (const PseBarSoil& bs : *soil)
            if (bs.bar == b.bar) {
              p = &bs.soil_p;
              mm = &bs.soil_m;
            }
        for (std::size_t k = 0; k < p->size() && k < b.stations.size(); ++k)
          o.line("    BARRE " + std::to_string(b.bar) + " XL " + fmt(b.stations[k].x) + "  PX " +
                 fmt((*p)[k].x()) + " PY " + fmt((*p)[k].y()) + " PZ " + fmt((*p)[k].z()) +
                 " MX " + fmt(k < mm->size() ? (*mm)[k] : 0.0));
      }
    }
  };

  // ---- ST1 edition (St1Options::st1_edition, keta/st1_edit.h) ----
  // The tables as ST1 prints them: Fortran formats, the components of the
  // analysis option, ST1's sign conventions (N compression-positive, Ty Tz Mx
  // My Mz opposite to keta's stations; displacements and reactions in global
  // axes, keta's own signs), only the tables the RESU block asks for.
  std::vector<int> cols;  // components of the analysis option (DOF 0..5)
  if (option_ == 0) cols = {0, 1, 5};
  else if (option_ == 1) cols = {2, 3, 4};
  else cols = {0, 1, 2, 3, 4, 5};
  static const char* const kDofName[6] = {"Dx", "Dy", "Dz", "Rx", "Ry", "Rz"};
  static const char* const kReacName[6] = {"Fx", "Fy", "Fz", "Mx", "My", "Mz"};
  static const char* const kEffName[6] = {"N", "Ty", "Tz", "Mx", "My", "Mz"};
  // `CABLE` in a RESU block asks for the cable tables of its results (ST1
  // edition); asked alone, it is all the block prints.
  const bool want_cable = req.kinds.count("CABLE") > 0;
  const bool none = !req.has_nodes && !req.has_appuis && !req.has_bars && !comps && !want_cable;
  const bool want_nodes = req.has_nodes || none;
  const bool want_appuis = req.has_appuis || none;
  const bool want_effort = req.effort || none || (req.has_bars && !comps);
  const bool want_depla = req.depla;
  // a column name starting at `pos` (one blank at least before it)
  auto at = [](std::string s, std::size_t pos, const std::string& name) {
    if (s.size() < pos) s.append(pos - s.size(), ' ');
    else s += ' ';
    return s + name;
  };
  auto eff_st1 = [](const Station& s, int c) {
    switch (c) {
      case 0: return -s.N;
      case 1: return -s.Vy;
      case 2: return -s.Vz;
      case 3: return -s.T;
      case 4: return -s.My;
      default: return -s.Mz;
    }
  };
  // Station displacements in global axes (a bar absent from the current
  // model, e.g. a phased one, reads its raw geometry).
  auto depl_st1 = [&](Id bar, const Station& s) {
    Mat3 Rt = Mat3::Identity();
    if (auto it = model_bar.find(bar); it != model_bar.end()) {
      const Bar* mb = it->second;
      Rt = bar_rotation(model_.node(mb->ni).x, model_.node(mb->nj).x, mb->roll).transpose();
    } else if (auto rb = bars_.find(bar); rb != bars_.end()) {
      Rt = bar_rotation(nodes_.at(rb->second.n1), nodes_.at(rb->second.n2), 0.0).transpose();
    }
    const Vec3 d = Rt * s.defl, w = Rt * s.rot;
    return std::array<double, 6>{d.x(), d.y(), d.z(), w.x(), w.y(), w.z()};
  };
  // CONTR (p.66-68): normal stresses at the section fibres, compression-positive.
  std::vector<std::string> contr_names;
  if (option_ == 0) contr_names = {"Ctr Vy", "Ctr Wy"};
  else if (option_ == 1) contr_names = {"Ctr Vz", "Ctr Wz"};
  else contr_names = {"Ctr VyVz", "Ctr WyVz", "Ctr VyWz", "Ctr WyWz"};
  // the soil of the PSE loading being edited (PRESS), else null
  const std::vector<PseBarSoil>* st1_soil = nullptr;
  // PRESS (p.274): the soil pressures of the option, in the bar axes.
  std::vector<std::string> press_names;
  if (option_ == 0) press_names = {"Px", "Py"};
  else if (option_ == 1) press_names = {"Pz", "Mx"};
  else press_names = {"Px", "Py", "Pz", "Mx"};
  // the soil pressures at station k of a bar, when it lies on elastic ground
  // (BarResult::soil_p, else the PSE soil of the loading being edited)
  auto press_at = [&](const BarResult& b, std::size_t k) -> std::optional<std::vector<double>> {
    const std::vector<Vec3>* p = &b.soil_p;
    const std::vector<double>* mm = &b.soil_m;
    if (p->empty() && st1_soil)
      for (const PseBarSoil& bs : *st1_soil)
        if (bs.bar == b.bar) {
          p = &bs.soil_p;
          mm = &bs.soil_m;
        }
    if (k >= p->size()) return std::nullopt;
    const Vec3& q = (*p)[k];
    const double mx = k < mm->size() ? (*mm)[k] : 0.0;
    if (option_ == 0) return std::vector<double>{q.x(), q.y()};
    if (option_ == 1) return std::vector<double>{q.z(), mx};
    return std::vector<double>{q.x(), q.y(), q.z(), mx};
  };
  // The sections of the result being edited (stamp_cara): ST1 computes its
  // stresses with the CARA in force at EXEC time.
  const std::map<Id, CaraProps>* res_cara = &cara_;
  auto use_cara_of = [&](const std::string& kind, Id id) { res_cara = &cara_of(kind, id); };
  // the stress bounds of the ENV being edited (EnvResult::contr), else null
  const std::map<Id, std::vector<EnvResult::ContrBound>>* env_contr = nullptr;
  // the CONTRG bounds of the ENV / SURCH being edited, else null; whether the
  // CONTRG constant applies (loadings, not surcharges: p.99)
  const std::map<Id, std::map<Id, std::vector<EnvResult::ContrBound>>>* env_contrg = nullptr;
  bool contrg_cst = true;
  // the reactions per APPUI of the result being edited (lo, hi), else null
  const std::map<Id, std::pair<Vec6, Vec6>>* appui_lohi = nullptr;
  std::string appui_head;  // a line before the APPUI table (ENV REAC concomitance)
  // ENV with a privileged EFFORT: the efforts and displacements edited are the
  // states concomitant with the extrema at each study section (p.231)
  const EnvResult* env_concom = nullptr;
  std::map<Id, std::pair<Vec6, Vec6>> appui_tmp;
  // the APPUI descriptions of the result being edited (MODIFIER APPUI moves
  // one to another node during a phasage)
  std::map<Id, SupportDef> appui_print = appui_;
  auto appui_single = [&](const CaseResult& r, const std::string& kind, Id id) {
    appui_tmp.clear();
    appui_print = appui_defs_of(kind, id);
    for (const auto& [aid, v] : appui_reactions(r, appuis_of(kind, id), false, &appui_print))
      appui_tmp[aid] = {v, v};
    appui_lohi = &appui_tmp;
  };
  // the cable tables of the result being edited (cable_tables), else null;
  // `cable_var`: with the Var.Tension column (construction states)
  const std::vector<CableTable>* cable_tabs = nullptr;
  bool cable_var = false;
  auto contr_st1 = [&](Id bid, const Station& s) { return contr_fibres(*res_cara, bid, s); };
  // The APPUI number ST1 prints in front of a support node (the smallest one
  // defined there).
  auto appui_of = [&](Id node) {
    Id best = 0;
    for (const auto& [aid, sd] : appui_)
      if (sd.node == node && (best == 0 || aid < best)) best = aid;
    return best;
  };
  // Round-off of an exact zero (a shear at mid-span, ~1e-16 of the moments)
  // prints as ST1's 0.0000: values below 1e-11 of the largest of their row.
  auto denoise = [](std::vector<double> v) {
    double big = 0;
    for (double x : v) big = std::max(big, std::fabs(x));
    for (double& x : v)
      if (std::fabs(x) < 1e-11 * big) x = 0.0;
    return v;
  };
  // One ST1 block: a result (`hi` null) or an envelope pair printed as Min /
  // Max lines (component-wise bounds in ST1's own signs). Each table ends on a
  // blank line, the RESU command on one more.
  bool st1_blocks = false;
  auto st1_block = [&](const std::string& label, Id id, const std::string& title,
                       const CaseResult& lo, const CaseResult* hi) {
    const bool pair = hi != nullptr;
    o.line(" --- RESULTATS " + label + " :" + st1::fortran_i(id, 9) + " ---");
    if (!title.empty()) {
      o.line(" ");
      o.line(" Titre : " + title);
    }
    o.line("  ");
    auto bounds = [&](double a, double b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
    if (want_nodes) {
      std::string h = "    Noeud";
      for (std::size_t i = 0; i < cols.size(); ++i)
        h = at(h, 17 + 12 * i, kDofName[cols[i]]);
      o.line(h);
      for (std::size_t k = 0; k < lo.node_disp.size(); ++k) {
        const Id nid = model_.nodes.size() > k ? model_.nodes[k].id : static_cast<Id>(k + 1);
        if (!node_wanted(nid)) continue;
        std::vector<double> va, vb;
        for (int c : cols) {
          va.push_back(lo.node_disp[k](c));
          vb.push_back(pair ? hi->node_disp[k](c) : 0.0);
        }
        if (!pair) {
          std::string ln = st1::fortran_i(nid, 11) + "    ";
          for (double v : denoise(va)) ln += " " + st1::fortran_g(v, 11, 5);
          o.line(ln);
        } else {
          std::vector<double> mn(va.size()), mx(va.size());
          for (std::size_t c = 0; c < va.size(); ++c)
            std::tie(mn[c], mx[c]) = bounds(va[c], vb[c]);
          std::string l1 = st1::fortran_i(nid, 11) + " Min", l2 = std::string(11, ' ') + " Max";
          for (double v : denoise(mn)) l1 += " " + st1::fortran_g(v, 11, 5);
          for (double v : denoise(mx)) l2 += " " + st1::fortran_g(v, 11, 5);
          o.line(l1);
          o.line(l2);
        }
      }
      o.line("");
    }
    if (want_appuis && appui_lohi) {
      // per declared APPUI (ST1 lists those of the result, an ENV the union
      // over its sources), in APPUI order
      if (!appui_head.empty()) o.line(appui_head);
      std::string h = "    Appui    Noeud";
      for (std::size_t i = 0; i < cols.size(); ++i)
        h = at(h, 26 + 12 * i, kReacName[cols[i]]);
      o.line(h);
      for (const auto& [aid, b] : *appui_lohi) {
        auto ad = appui_print.find(aid);
        if (ad == appui_print.end() || !appui_wanted(ad->second.node)) continue;
        const std::string lead = st1::fortran_i(aid, 11) + st1::fortran_i(ad->second.node, 11);
        // (an envelope's bounds are per global component: turned as a pair
        // they bound the local components only approximately)
        const Vec6 b1 = appui_frame(aid, b.first), b2 = appui_frame(aid, b.second);
        std::vector<double> mn, mx;
        for (int c : cols) {
          mn.push_back(std::min(b1(c), b2(c)));
          mx.push_back(std::max(b1(c), b2(c)));
        }
        if (!pair) {
          std::string ln = lead + "    ";
          for (double v : denoise(mn)) ln += " " + st1::fortran_g(v, 11, 5);
          o.line(ln);
          continue;
        }
        std::string l1 = lead + " Min", l2 = std::string(22, ' ') + " Max";
        for (double v : denoise(mn)) l1 += " " + st1::fortran_g(v, 11, 5);
        for (double v : denoise(mx)) l2 += " " + st1::fortran_g(v, 11, 5);
        o.line(l1);
        o.line(l2);
      }
      o.line("");
    } else if (want_appuis) {
      std::string h = "    Appui    Noeud";
      for (std::size_t i = 0; i < cols.size(); ++i)
        h = at(h, 26 + 12 * i, kReacName[cols[i]]);
      o.line(h);
      std::vector<std::pair<Id, std::size_t>> rows;  // (appui, reaction index)
      for (std::size_t r = 0; r < lo.reactions.size(); ++r)
        if (appui_wanted(lo.reactions[r].node)) rows.push_back({appui_of(lo.reactions[r].node), r});
      std::stable_sort(rows.begin(), rows.end(),
                       [](const auto& a, const auto& b) { return a.first < b.first; });
      for (const auto& [aid, r] : rows) {
        const Reaction& ra = lo.reactions[r];
        const std::string lead = st1::fortran_i(aid, 11) + st1::fortran_i(ra.node, 11);
        const Vec6 Ra = appui_frame(aid, ra.R);
        std::vector<double> va;
        for (int c : cols) va.push_back(Ra(c));
        if (!pair) {
          std::string ln = lead + "    ";
          for (double v : denoise(va)) ln += " " + st1::fortran_g(v, 11, 5);
          o.line(ln);
          continue;
        }
        const Reaction* rb = nullptr;
        for (const Reaction& x : hi->reactions)
          if (x.node == ra.node) rb = &x;
        std::vector<double> mn(va.size()), mx(va.size());
        for (std::size_t c = 0; c < va.size(); ++c)
          std::tie(mn[c], mx[c]) = bounds(va[c], rb ? appui_frame(aid, rb->R)(cols[c]) : va[c]);
        std::string l1 = lead + " Min", l2 = std::string(22, ' ') + " Max";
        for (double v : denoise(mn)) l1 += " " + st1::fortran_g(v, 11, 5);
        for (double v : denoise(mx)) l2 += " " + st1::fortran_g(v, 11, 5);
        o.line(l1);
        o.line(l2);
      }
      o.line("");
    }
    // a bar table: `vals(bar, station)` gives the printed values of a station
    // `fix` may overwrite the bounds of an envelope pair (CONTR of an ENV)
    // `only` (when set) restricts the bars (a CONTRG table: the bars it is
    // defined on); `tab_bar` / `tab_k` name the row being built for `vals`.
    using FixFn = std::function<void(Id, double, std::vector<double>&, std::vector<double>&)>;
    const BarResult* tab_bar = nullptr;
    std::size_t tab_k = 0;
    auto bar_table = [&](const std::string& header, auto&& vals, const FixFn& fix = nullptr,
                         const std::function<bool(Id)>& only = nullptr) {
      o.line(header);
      for (std::size_t j = 0; j < lo.bars.size(); ++j) {
        const BarResult& ba = lo.bars[j];
        if (!bar_wanted(ba.bar) || ba.stations.empty() || (only && !only(ba.bar))) continue;
        const BarResult* bb = pair && j < hi->bars.size() ? &hi->bars[j] : nullptr;
        for (std::size_t k = 0; k < ba.stations.size(); ++k) {
          tab_bar = &ba;
          tab_k = k;
          const std::vector<double> va = vals(ba.bar, ba.stations[k]);
          const std::string lead = st1::fortran_i(ba.bar, 11) + st1::fortran_f(ba.stations[k].x, 11, 3);
          if (!pair) {
            std::string ln = lead + "    ";
            for (double v : denoise(va)) ln += st1::fortran_g(v, 13, 5);
            o.line(ln);
            continue;
          }
          tab_bar = bb;
          const std::vector<double> vb =
              bb && k < bb->stations.size() ? vals(bb->bar, bb->stations[k]) : va;
          std::vector<double> mn(va.size()), mx(va.size());
          for (std::size_t c = 0; c < va.size(); ++c)
            std::tie(mn[c], mx[c]) = bounds(va[c], c < vb.size() ? vb[c] : va[c]);
          if (fix) fix(ba.bar, ba.stations[k].x, mn, mx);
          std::string l1 = lead + " Min", l2 = std::string(22, ' ') + " Max";
          for (double v : denoise(mn)) l1 += st1::fortran_g(v, 13, 5);
          for (double v : denoise(mx)) l2 += st1::fortran_g(v, 13, 5);
          o.line(l1);
          o.line(l2);
        }
        o.line("");
      }
    };
    // the stresses follow the efforts on their rows in PLANE / GRILL; in
    // SPATIALE ST1 edits their four corners as a table of their own after
    // the efforts
    const bool contr_inline = req.contr && option_ != 2;
    FixFn contr_fix = req.contr && env_contr ? FixFn([&](Id bar, double x, std::vector<double>& mn,
                                                         std::vector<double>& mx) {
      auto b = env_contr->find(bar);
      if (b == env_contr->end()) return;
      for (const EnvResult::ContrBound& cb : b->second) {
        if (std::fabs(cb.x - x) > 1e-6 * (1.0 + std::fabs(x))) continue;
        const std::size_t off = mn.size() - cb.lo.size();
        for (std::size_t k = 0; k < cb.lo.size(); ++k) {
          mn[off + k] = cb.lo[k];
          mx[off + k] = cb.hi[k];
        }
        return;
      }
    }) : FixFn();
    // the concomitant tables of an ENV with a privileged EFFORT / CONTR /
    // CONTRG (p.231): per study section the state at the min of the
    // component, then at its max, under the line naming the component
    std::string concom_line;
    bool concom_swap = false;  // ST1 lists every effort negated: its min is keta's max
    if (env_concom) {
      static const char* const kFib[8] = {"Vy", "Wy", "Vz", "Wz", "VyVz", "VyWz", "WyVz", "WyWz"};
      if (env_concom->what == EnvWhat::Effort) {
        concom_line = " ==> Resultats concomitants des extrema de l'effort " +
                      std::string(kEffName[std::clamp(env_concom->comp, 0, 5)]);
        concom_swap = true;
      } else if (env_concom->what == EnvWhat::Contr) {
        concom_line = " ==> Resultats concomitants des extrema de la contrainte Ctr " +
                      std::string(kFib[std::clamp(static_cast<int>(env_concom->fiber), 0, 7)]);
      } else {
        concom_line = " ==> Resultats concomitants des extrema de la contrainte generalisee :" +
                      st1::fortran_i(env_concom->contrg, 10);
      }
    }
    auto concom_table = [&](const std::string& header, auto&& vals) {
      o.line(concom_line);
      o.line(header);
      Id last = 0;
      for (const EnvPoint& p : env_concom->points) {
        if (!bar_wanted(p.bar)) continue;
        if (last != 0 && p.bar != last) o.line("");
        last = p.bar;
        const EnvPoint* pp = &p;
        std::string l1 = st1::fortran_i(p.bar, 11) + st1::fortran_f(p.x, 11, 3) + " Min";
        std::string l2 = std::string(22, ' ') + " Max";
        for (double v : denoise(vals(*pp, concom_swap ? p.st_max : p.st_min)))
          l1 += st1::fortran_g(v, 13, 5);
        for (double v : denoise(vals(*pp, concom_swap ? p.st_min : p.st_max)))
          l2 += st1::fortran_g(v, 13, 5);
        o.line(l1);
        o.line(l2);
      }
      if (last != 0) o.line("");
    };
    if (want_effort && env_concom) {
      std::string h = "    Barre   Abscisse";
      for (std::size_t i = 0; i < cols.size(); ++i) h = at(h, 26 + 13 * i, kEffName[cols[i]]);
      if (req.contr && option_ != 2) {
        h += "         " + contr_names[0];
        for (std::size_t i = 1; i < contr_names.size(); ++i) h += "      " + contr_names[i];
        h += "   ";
      } else {
        h += "      ";
      }
      concom_table(h, [&](const EnvPoint& p, const Station& s) {
        std::vector<double> v;
        for (int c : cols) v.push_back(eff_st1(s, c));
        if (req.contr && option_ != 2)
          for (double x : contr_st1(p.bar, s)) v.push_back(x);
        return v;
      });
      if (req.contr && option_ == 2) {
        std::string hc = "    Barre   Abscisse   " + contr_names[0];
        for (std::size_t i = 1; i < contr_names.size(); ++i) hc += "    " + contr_names[i];
        concom_table(hc + "  ", [&](const EnvPoint& p, const Station& s) { return contr_st1(p.bar, s); });
      }
    } else if (want_effort) {
      std::string h = "    Barre   Abscisse";
      for (std::size_t i = 0; i < cols.size(); ++i) h = at(h, 26 + 13 * i, kEffName[cols[i]]);
      if (contr_inline) {
        h += "         " + contr_names[0];
        for (std::size_t i = 1; i < contr_names.size(); ++i) h += "      " + contr_names[i];
        h += "   ";
      } else {
        h += "      ";
      }
      bar_table(h, [&](Id bar, const Station& s) {
        std::vector<double> v;
        for (int c : cols) v.push_back(eff_st1(s, c));
        if (contr_inline)
          for (double x : contr_st1(bar, s)) v.push_back(x);
        return v;
      }, contr_inline ? contr_fix : FixFn());
    }
    if (req.contr && !contr_inline) {
      std::string h = "    Barre   Abscisse   " + contr_names[0];
      for (std::size_t i = 1; i < contr_names.size(); ++i) h += "    " + contr_names[i];
      bar_table(h + "  ", [&](Id bar, const Station& s) { return contr_st1(bar, s); }, contr_fix);
    }
    if (want_depla && env_concom) {
      std::string h = "    Barre   Abscisse";
      for (std::size_t i = 0; i < cols.size(); ++i) h = at(h, 26 + 13 * i, kDofName[cols[i]]);
      concom_table(h + "      ", [&](const EnvPoint& p, const Station& s) {
        const std::array<double, 6> d = depl_st1(p.bar, s);
        std::vector<double> v;
        for (int c : cols) v.push_back(d[static_cast<std::size_t>(c)]);
        return v;
      });
    } else if (want_depla) {
      // PRESS (p.274): the soil pressures follow the displacements on the
      // bars on elastic ground (Exemple11: Px Py, Py = k·Dy), the other rows
      // stop after the displacements.
      std::string h = "    Barre   Abscisse";
      for (std::size_t i = 0; i < cols.size(); ++i) h = at(h, 26 + 13 * i, kDofName[cols[i]]);
      if (req.press) {
        std::size_t pos = 26 + 13 * cols.size() + 2;
        for (std::size_t i = 0; i < press_names.size(); ++i, pos += 12) h = at(h, pos, press_names[i]);
        h += "   ";
      } else {
        h += "      ";
      }
      bar_table(h, [&](Id bar, const Station& s) {
        const std::array<double, 6> d = depl_st1(bar, s);
        std::vector<double> v;
        for (int c : cols) v.push_back(d[static_cast<std::size_t>(c)]);
        if (req.press && tab_bar)
          if (auto p = press_at(*tab_bar, tab_k))
            for (double x : *p) v.push_back(x);
        return v;
      });
    }
    if (req.contrg_all || !req.contrg.empty()) {
      // CONTRG (p.99): one table per generalized stress, on the bars it is
      // defined on — laid out as the CONTRG table of ST1's ENV (Exemple7).
      // the CONTRG listed, else all of them (`CONTRG` alone)
      std::vector<Id> cgs;
      for (Id cid : req.contrg.empty() ? std::vector<Id>{} : req.contrg)
        if (contrgs.count(cid)) cgs.push_back(cid);
        else warn("ST1: RESU CONTRG " + std::to_string(cid) + " — non definie");
      if (req.contrg.empty())
        for (const auto& [cid, cd] : contrgs) cgs.push_back(cid);
      // ST1: a title line per generalized stress, then ONE table with a
      // column each, on the bars where one of them is defined
      if (!cgs.empty()) {
        for (Id cid : cgs)
          o.line(" Contrainte" + st1::fortran_i(cid, 10) + " " + contrgs.at(cid).title);
        o.line(" ");
        std::string h = "    Barre Abscisse" + st1::fortran_i(cgs.front(), 16);
        for (std::size_t k = 1; k < cgs.size(); ++k) h += st1::fortran_i(cgs[k], 12);
        if (env_concom) {
          concom_table(h + " ", [&](const EnvPoint& p, const Station& s) {
            const double L = bar_length(p.bar);
            const double xx = L > 0 ? std::clamp(s.x / L, 0.0, 1.0) : 0.0;
            std::vector<double> v;
            for (Id cid : cgs) {
              const ContrgDef& cd = contrgs.at(cid);
              v.push_back(cd.defined_for(p.bar) ? contrg_at(cd, p.bar, xx, s, contrg_cst) : 0.0);
            }
            return v;
          });
        } else
        bar_table(h + " ",
                  [&](Id bar, const Station& s) {
                    const double L = bar_length(bar);
                    const double xx = L > 0 ? std::clamp(s.x / L, 0.0, 1.0) : 0.0;
                    std::vector<double> v;
                    for (Id cid : cgs) {
                      const ContrgDef& cd = contrgs.at(cid);
                      v.push_back(cd.defined_for(bar) ? contrg_at(cd, bar, xx, s, contrg_cst) : 0.0);
                    }
                    return v;
                  },
                  env_contrg ? FixFn([&](Id bar, double x, std::vector<double>& mn,
                                         std::vector<double>& mx) {
                    for (std::size_t k = 0; k < cgs.size(); ++k) {
                      auto g = env_contrg->find(cgs[k]);
                      if (g == env_contrg->end()) continue;
                      auto b = g->second.find(bar);
                      if (b == g->second.end()) continue;
                      for (const EnvResult::ContrBound& cb : b->second)
                        if (std::fabs(cb.x - x) <= 1e-6 * (1.0 + std::fabs(x))) {
                          mn[k] = cb.lo[0];
                          mx[k] = cb.hi[0];
                          break;
                        }
                    }
                  })
                             : FixFn(),
                  [&](Id bar) {
                    return std::any_of(cgs.begin(), cgs.end(),
                                       [&](Id cid) { return contrgs.at(cid).defined_for(bar); });
                  });
      }
    }
    if (want_cable && cable_tabs && !cable_tabs->empty()) {
      o.line(cable_var ? "  Cable       X           Y           Z          Tension      Var.Tension"
                       : "  Cable       X           Y           Z          Tension");
      for (const CableTable& ct : *cable_tabs) {
        // XS: the rows at the asked curvilinear abscissae, interpolated
        // between the cable points (the length is the polyline's)
        std::vector<CableRow> rows = ct.rows;
        if (req.has_xs && ct.rows.size() >= 2) {
          std::vector<double> sc{0.0};
          for (std::size_t k = 1; k < ct.rows.size(); ++k)
            sc.push_back(sc.back() + std::hypot(ct.rows[k].x - ct.rows[k - 1].x,
                                                ct.rows[k].y - ct.rows[k - 1].y,
                                                ct.rows[k].z - ct.rows[k - 1].z));
          const double len = sc.back();
          rows.clear();
          for (double x : req.xs) {
            const double sx = req.xs_rel ? x * len : x;
            if (sx < -1e-9 || sx > len * (1 + 1e-9)) continue;
            std::size_t k = 1;
            while (k + 1 < sc.size() && sc[k] < sx) ++k;
            const double t = sc[k] > sc[k - 1] ? std::clamp((sx - sc[k - 1]) / (sc[k] - sc[k - 1]), 0.0, 1.0) : 0.0;
            const CableRow &a = ct.rows[k - 1], &b = ct.rows[k];
            CableRow r;
            r.x = a.x + t * (b.x - a.x);
            r.y = a.y + t * (b.y - a.y);
            r.z = a.z + t * (b.z - a.z);
            r.sigma = a.sigma + t * (b.sigma - a.sigma);
            r.dsigma = a.dsigma + t * (b.dsigma - a.dsigma);
            rows.push_back(r);
          }
        }
        for (const CableRow& r : rows) {
          std::string ln = st1::fortran_i(ct.cable, 7);
          for (double v : {r.x, r.y, r.z, r.sigma}) ln += st1::fortran_g(v, 14, 5);
          if (cable_var) ln += st1::fortran_g(r.dsigma, 14, 5);
          o.line(ln);
        }
        o.line(" ");
      }
      o.line("");
    }
    st1_blocks = true;
  };
  // The title ST1 prints, or none when the result carries keta's default name.
  auto st1_title = [](const std::string& name, const std::string& kind, Id id) {
    return (name.empty() || name == kind + " " + std::to_string(id)) ? std::string() : name;
  };
  static const std::map<std::string, std::string> kSt1Label = {
      {"CHARG", "DU CHARGEMENT"},  {"COMB", "DE LA COMBINAISON"}, {"ETAT", "DE L'ETAT"},
      {"PHASAGE", "DU PHASAGE"},   {"SURCH", "DE LA SURCHARGE"},  {"ENV", "DE L'ENVELOPPE"},
      {"INCENDIE", "CHARGE INCENDIE"}};
  const bool st1_ed = opt_.st1_edition;

  bool any = false;
  for (const auto& [kind, ids_in] : req.kinds) {
    if (kind == "MODE" && st1_ed) {
      // MODE PROPRE (p.240): period, participation along the sum of the
      // translation directions of the option and along each of them, then
      // the mass-normalized shape (Masse Modale Generalisee 1).
      require(store_.has_modal, "ST1: RESU MODE — EXEC MODE prealable requis");
      const ModalResults& mr = store_.modal;
      const std::vector<int> dirs = option_ == 0   ? std::vector<int>{0, 1}
                                    : option_ == 1 ? std::vector<int>{2}
                                                   : std::vector<int>{0, 1, 2};
      double mass_d = 0;
      for (int d : dirs) mass_d += mr.total_mass[static_cast<std::size_t>(d)];
      auto pad = [](std::string s, std::size_t n) {
        if (s.size() < n) s.append(n - s.size(), ' ');
        return s;
      };
      auto pct = [](double a, double b) { return b > 0 ? 100.0 * a / b : 0.0; };
      double cum_d = 0;
      std::array<double, 6> cum{};
      static const char* const kAxis[3] = {"X", "Y", "Z"};
      for (std::size_t k = 0; k < mr.frequencies_hz.size(); ++k) {
        double g = 0;
        for (int d : dirs) g += k < mr.participation.size() ? mr.participation[k][static_cast<std::size_t>(d)] : 0.0;
        // ST1 orients a mode by its participation factor (positive), and a
        // round-off factor (a mode orthogonal to the excitation) reads 0.
        const bool null_g = std::fabs(g) <= 1e-9 * std::sqrt(std::max(mass_d, 1e-300));
        if (null_g) g = 0;
        const double orient = g < 0 ? -1.0 : 1.0;
        g *= orient;
        cum_d += g * g;
        for (int d : dirs)
          if (k < mr.effective_mass.size())
            cum[static_cast<std::size_t>(d)] += mr.effective_mass[k][static_cast<std::size_t>(d)];
        if (!req.modes.empty() &&
            !std::count(req.modes.begin(), req.modes.end(), static_cast<Id>(k + 1)))
          continue;
        const double f = mr.frequencies_hz[k];
        o.line(" --- RESULTATS DU MODE PROPRE :" + st1::fortran_i(static_cast<long long>(k + 1), 9) + " ---");
        o.line("  ");
        if (req.mode_freq || !req.mode_depla) {
          o.line(pad(pad("  Periode", 30) + st1::fortran_g(f > 0 ? 1.0 / f : 0.0, 11, 5) + " sec", 54) +
                 "Frequence    :  " + st1::fortran_g(f, 12, 5) + " Hz");
          o.line(pad(pad("  Facteur Participation Mod", 30) + st1::fortran_g(g, 11, 5), 54) +
                 "Masse Modale :  " + st1::fortran_g(g * g, 12, 5));
          o.line(pad("  Masse Modale Generalisee", 30) + st1::fortran_g(1.0, 11, 5));
          o.line(pad(pad("  Taux Participation Modale", 30) + st1::fortran_g(pct(g * g, mass_d), 11, 5) + "  %", 54) +
                 "Cumul        :  " + st1::fortran_g(pct(cum_d, mass_d), 12, 5) + "  %");
          o.line(" ");
          for (int d : dirs) {
            const std::size_t du = static_cast<std::size_t>(d);
            const double md = mr.total_mass[du];
            double eff = k < mr.effective_mass.size() ? mr.effective_mass[k][du] : 0.0;
            if (eff <= 1e-18 * std::max(md, 1e-300)) eff = 0;
            o.line(std::string("                ---> Direction :  ") + kAxis[d]);
            o.line("            Taux Participation Modale   " + st1::fortran_g(pct(eff, md), 11, 5) +
                   " %   Cumul : " + st1::fortran_g(pct(cum[du], md), 11, 5) +
                   " %   Masse mobile cumulee : " + st1::fortran_g(cum[du], 11, 5));
          }
          o.line(" ");
        }
        if (req.mode_depla && k < mr.modes.size()) {
          std::string h = "  Noeud";
          for (std::size_t i = 0; i < cols.size(); ++i) {
            std::string nm = kDofName[cols[i]];
            for (char& ch : nm) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            h = at(h, 14 + 13 * i, nm);
          }
          o.line(h);
          // unit-crest shape x its crest in the M-normalized one
          const double c = k < mr.modes_free.size() && mr.modes_free[k].size() > 0
                               ? mr.modes_free[k].cwiseAbs().maxCoeff()
                               : 1.0;
          for (std::size_t n = 0; n < mr.modes[k].size(); ++n) {
            const Id nid = model_.nodes.size() > n ? model_.nodes[n].id : static_cast<Id>(n + 1);
            if (!node_wanted(nid)) continue;
            std::string ln = st1::fortran_i(nid, 10);
            std::vector<double> v;
            for (int col : cols) v.push_back(orient * c * mr.modes[k][n](col));
            for (double x : denoise(v)) ln += st1::fortran_g(x, 13, 5);
            o.line(ln);
          }
        }
        o.line("");
      }
      st1_blocks = true;
      any = true;
      continue;
    }
    if (kind == "FLAMBEMENT" && st1_ed) {
      // MODE CRITIQUE (p.241, p.215): the critical coefficient, the modal
      // quantities of the dynamics (p.200) with the mass replaced by the
      // geometric stiffness KG of the variable load — participation on the
      // unit vector of every free DOF, Rigidite Modale = gamma^2, rate and
      // cumulated rate over R^T KG R — and the shape normalized to
      // x^T KG x = 1 (Exemple8 as ST1's, up to the sign of each mode).
      require(store_.has_buckling, "ST1: RESU FLAMBEMENT — EXEC FLAMBEMENT requis");
      const BucklingResults& br = store_.buckling;
      auto pad = [](std::string s, std::size_t n) {
        if (s.size() < n) s.append(n - s.size(), ' ');
        return s;
      };
      double cum = 0;
      for (std::size_t k = 0; k < br.load_factors.size(); ++k) {
        const double gam = k < br.participation.size() ? br.participation[k] : 0.0;
        const double rate = br.total_kg > 0 ? 100.0 * gam * gam / br.total_kg : 0.0;
        cum += rate;
        if (!req.modes.empty() &&
            !std::count(req.modes.begin(), req.modes.end(), static_cast<Id>(k + 1)))
          continue;
        o.line(" --- RESULTATS DU MODE CRITIQUE :" + st1::fortran_i(static_cast<long long>(k + 1), 5) + " ---");
        o.line("  ");
        o.line("  Coeff critique (val propre)  :" + st1::fortran_g(br.load_factors[k], 14, 5) + "    ");
        o.line(pad("  Facteur Participation Modale :" + st1::fortran_g(gam, 14, 5), 59) +
               "Rigidite Modale :" + st1::fortran_g(gam * gam, 14, 5));
        o.line("  Rigidite Modale Generalisee   " + st1::fortran_g(1.0, 14, 5));
        o.line(pad("  Taux Participation Modale    :" + st1::fortran_g(rate, 14, 5) + "  %", 62) +
               "Cumul        :" + st1::fortran_g(cum, 14, 5) + "  %");
        o.line(" ");
        if (req.mode_depla && k < br.modes.size()) {
          std::string h = "  Noeud";
          for (std::size_t i = 0; i < cols.size(); ++i) {
            std::string nm = kDofName[cols[i]];
            for (char& ch : nm) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
            h = at(h, 14 + 13 * i, nm);
          }
          o.line(h);
          const double gk = k < br.gen_stiffness.size() ? br.gen_stiffness[k] : 0.0;
          const double scale = gk > 0 ? 1.0 / std::sqrt(gk) : 1.0;
          for (std::size_t n = 0; n < br.modes[k].size(); ++n) {
            const Id nid = model_.nodes.size() > n ? model_.nodes[n].id : static_cast<Id>(n + 1);
            if (!node_wanted(nid)) continue;
            std::string ln = st1::fortran_i(nid, 10);
            std::vector<double> v;
            for (int col : cols) v.push_back(scale * br.modes[k][n](col));
            for (double x : denoise(v)) ln += st1::fortran_g(x, 14, 5);
            o.line(ln);
          }
        }
        o.line("");
      }
      st1_blocks = true;
      any = true;
      continue;
    }
    if (kind == "SPECTRE" && st1_ed) {
      // REPONSE SPECTRALE (p.240): the combined (CQC / SRSS) magnitudes, in
      // G13.4, a blank line after every support and every bar.
      std::vector<Id> ids = ids_in;
      if (ids.empty()) ids = spectre_last;
      if (ids.empty())
        for (const auto& [sid, r] : store_.spectre) ids.push_back(sid);
      static const char* const kSpecReac[6] = {"Rx", "Ry", "Rz", "Mx", "My", "Mz"};
      for (Id id : ids) {
        auto it = store_.spectre.find(id);
        if (it == store_.spectre.end()) {
          warn("ST1: RESU — SPECTRE " + std::to_string(id) + " non calcule");
          continue;
        }
        const CaseResult& r = it->second;
        o.line(" --- RESULTATS REPONSE SPECTRALE :" + st1::fortran_i(id, 9) + " ---");
        o.line("  ");
        auto right = [](std::string s, std::size_t end, const std::string& name) {
          const std::size_t pos = end > name.size() ? end - name.size() : 0;
          if (s.size() < pos) s.append(pos - s.size(), ' ');
          else s += ' ';
          return s + name;
        };
        if (want_appuis) {
          std::string h = "   Appui";
          for (std::size_t i = 0; i < cols.size(); ++i) h = right(h, 16 + 13 * i, kSpecReac[cols[i]]);
          o.line(h);
          for (const auto& [aid, ad] : appui_) {
            if (!appui_wanted(ad.node)) continue;
            Vec6 R = Vec6::Zero();
            for (const Reaction& re : r.reactions)
              if (re.node == ad.node) R = re.R;
            std::string ln = st1::fortran_i(aid, 10);
            for (int col : cols) ln += st1::fortran_g(std::fabs(R(col)), 13, 4);
            o.line(ln);
            o.line(" ");
          }
        }
        auto bar_rows = [&](const std::string& h0, auto&& names, auto&& vals) {
          std::string h = h0;
          for (std::size_t i = 0; i < cols.size(); ++i) h = right(h, 26 + 13 * i, names(cols[i]));
          o.line(h);
          for (const BarResult& ba : r.bars) {
            if (!bar_wanted(ba.bar) || ba.stations.empty()) continue;
            for (const Station& stn : ba.stations) {
              std::string ln = st1::fortran_i(ba.bar, 10) + st1::fortran_f(stn.x, 8, 3) + "    ";
              for (double v : vals(ba.bar, stn)) ln += st1::fortran_g(std::fabs(v), 13, 4);
              o.line(ln);
            }
            o.line(" ");
          }
        };
        if (want_effort)
          bar_rows("   Barre  Abscisse", [&](int c) { return std::string(kEffName[c]); },
                   [&](Id, const Station& stn) {
                     std::vector<double> v;
                     for (int c : cols) v.push_back(eff_st1(stn, c));
                     return v;
                   });
        if (want_depla)
          bar_rows("   Barre  Abscisse", [&](int c) { return std::string(kDofName[c]); },
                   [&](Id bar, const Station& stn) {
                     const std::array<double, 6> d = depl_st1(bar, stn);
                     std::vector<double> v;
                     for (int c : cols) v.push_back(d[static_cast<std::size_t>(c)]);
                     return v;
                   });
      }
      st1_blocks = true;
      any = true;
      continue;
    }
    if (kind == "HISTORIQUE" && st1_ed) {
      // HISTORIQUE (p.240-241): per support its reactions step by step, per
      // bar and study section the efforts, then the generalized stresses
      // asked; `historique i` alone prints the block header only (Exemple7).
      std::vector<Id> ids = ids_in;
      if (ids.empty())
        for (const auto& [hid, hs] : historique_series) ids.push_back(hid);
      auto pad = [](std::string s, std::size_t n) {
        if (s.size() < n) s.append(n - s.size(), ' ');
        return s;
      };
      for (Id id : ids) {
        auto it = historique_series.find(id);
        if (it == historique_series.end() || it->second.states.empty()) {
          warn("ST1: RESU HISTORIQUE " + std::to_string(id) +
               " — serie temporelle non calculee (EXEC HISTORIQUE, TEMPS OUI)");
          continue;
        }
        const HistorySeries& hs = it->second;
        o.line(" --- RESULTATS DE L HISTORIQUE :" + st1::fortran_i(id, 9) + " --- " + hs.title);
        o.line("  ");
        if (req.has_appuis) {
          for (const auto& [aid, ad] : appui_) {
            if (!appui_wanted(ad.node)) continue;
            o.line("  Appui :" + st1::fortran_i(aid, 10));
            std::string h = "     Temps";
            for (std::size_t c = 0; c < cols.size(); ++c) {
              std::string nm = kReacName[cols[c]];
              for (char& ch : nm) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
              h = pad(h, 35 + 13 * c) + nm;
            }
            o.line(h);
            for (std::size_t s = 0; s < hs.states.size(); ++s) {
              Vec6 R = Vec6::Zero();
              for (const Reaction& re : hs.states[s].reactions)
                if (re.node == ad.node) R = re.R;
              std::vector<double> v;
              for (int c : cols) v.push_back(R(c));
              std::string ln = pad("   " + st1::fortran_g(hs.times[s], 11, 5), 30);
              for (double x : denoise(v)) ln += st1::fortran_g(x, 13, 5);
              o.line(ln);
            }
            o.line("");
          }
        }
        const bool efforts = req.has_bars && (req.effort || !comps);
        std::vector<Id> cgs = req.contrg;
        if (cgs.empty() && req.contrg_all)
          for (const auto& [cid, cd] : contrgs) cgs.push_back(cid);
        if (req.has_bars && (efforts || !cgs.empty())) {
          o.line("           Sollicitation  :  K*DU+C*VU+M*AU ");
          const CaseResult& s0 = hs.states.front();
          for (std::size_t j = 0; j < s0.bars.size(); ++j) {
            const Id bar = s0.bars[j].bar;
            if (!bar_wanted(bar) || s0.bars[j].stations.empty()) continue;
            o.line("  Barre :" + st1::fortran_i(bar, 10));
            double L = 0;
            if (auto rb = bars_.find(bar); rb != bars_.end())
              L = (nodes_.at(rb->second.n2) - nodes_.at(rb->second.n1)).norm();
            const std::size_t ns = s0.bars[j].stations.size();
            auto lead = [&](std::size_t s, double x) {
              return "   " + st1::fortran_g(hs.times[s], 11, 5) + "   " + st1::fortran_g(x, 11, 5);
            };
            if (efforts)
              for (std::size_t k = 0; k < ns; ++k) {
                std::string h = "         Temps   Abscisse";
                for (std::size_t c = 0; c < cols.size(); ++c) h = pad(h, 35 + 13 * c) + kEffName[cols[c]];
                o.line(h + "   ");
                for (std::size_t s = 0; s < hs.states.size(); ++s) {
                  const Station& st = hs.states[s].bars[j].stations[k];
                  std::vector<double> v;  // (Mz: the stored series' ST1 sign)
                  for (int c : cols) v.push_back(eff_st1(st, c));
                  std::string ln = pad(lead(s, st.x), 34);
                  for (double x : denoise(v)) ln += st1::fortran_g(x, 13, 5);
                  o.line(ln);
                }
                o.line("");
              }
            for (Id cid : cgs) {
              auto cg = contrgs.find(cid);
              if (cg == contrgs.end() || !cg->second.defined_for(bar)) continue;
              for (std::size_t k = 0; k < ns; ++k) {
                o.line("   Temps  Abscisse   contrainte gen" + st1::fortran_i(cid, 12));
                for (std::size_t s = 0; s < hs.states.size(); ++s) {
                  const Station& st = hs.states[s].bars[j].stations[k];
                  const double xx = L > 0 ? std::clamp(st.x / L, 0.0, 1.0) : 0.0;
                  o.line(pad(lead(s, st.x), 39) +
                         st1::fortran_g(contrg_at(cg->second, bar, xx, st), 13, 5));
                }
                o.line("");
              }
            }
          }
        }
        o.line("");
      }
      st1_blocks = true;
      any = true;
      continue;
    }
    if (kind == "MODE") {
      require(store_.has_modal, "ST1: RESU MODE — EXEC MODE prealable requis");
      o.line("== MODES PROPRES ==");
      for (std::size_t k = 0; k < store_.modal.frequencies_hz.size(); ++k) {
        if (!req.modes.empty() &&
            !std::count(req.modes.begin(), req.modes.end(), static_cast<Id>(k + 1)))
          continue;
        o.line("  MODE " + std::to_string(k + 1) + "  F " +
               fmt(store_.modal.frequencies_hz[k]) + " Hz");
      }
      any = true;
      continue;
    }
    if (kind == "FLAMBEMENT") {
      require(store_.has_buckling, "ST1: RESU FLAMBEMENT — EXEC FLAMBEMENT requis");
      o.line("== MODES CRITIQUES ==");
      for (std::size_t k = 0; k < store_.buckling.load_factors.size(); ++k) {
        if (!req.modes.empty() &&
            !std::count(req.modes.begin(), req.modes.end(), static_cast<Id>(k + 1)))
          continue;
        o.line("  MODE " + std::to_string(k + 1) + "  ALPHA " +
               fmt(store_.buckling.load_factors[k]));
        if (req.mode_depla && k < store_.buckling.modes.size())
          for (std::size_t n = 0; n < store_.buckling.modes[k].size(); ++n)
            o.line("    NOEUD " + std::to_string(model_.nodes[n].id) + "  " +
                   fmt(store_.buckling.modes[k][n](0)) + " " +
                   fmt(store_.buckling.modes[k][n](1)) + " " +
                   fmt(store_.buckling.modes[k][n](2)));
      }
      any = true;
      continue;
    }
    if (kind == "ENV") {
      std::vector<Id> ids = ids_in;
      if (ids.empty())
        for (const auto& [id, r] : env_full) ids.push_back(id);
      for (Id id : ids) {
        auto it = env_full.find(id);
        if (it == env_full.end()) {
          warn("ST1: RESU — ENV " + std::to_string(id) + " non calculee");
          continue;
        }
        const EnvResult& er = it->second;
        if (st1_ed) {
          use_cara_of("ENV", id);
          env_contr = &er.contr;
          env_contrg = &er.contrg_b;
          appui_print = appui_;
          appui_lohi = &er.appui_bounds;
          if (er.has_priv && !er.points.empty() &&
              (er.what == EnvWhat::Effort || er.what == EnvWhat::Contr ||
               er.what == EnvWhat::Contrg))
            env_concom = &er;
          if (!er.reac_concom.empty()) {  // ENV … REAC comp (p.231)
            static const char* const kR[6] = {"Fx", "Fy", "Fz", "Mx", "My", "Mz"};
            appui_lohi = &er.reac_concom;
            appui_head = " ==> Resultats concomitants des extrema de la reaction " +
                         std::string(kR[std::clamp(er.reac_comp, 0, 5)]);
          }
          st1_block(kSt1Label.at("ENV"), id, st1_title(er.name, "ENV", id), er.bounds.lo.value,
                    &er.bounds.hi.value);
          env_contr = nullptr;
          env_contrg = nullptr;
          appui_lohi = nullptr;
          appui_head.clear();
          env_concom = nullptr;
          any = true;
          continue;
        }
        o.line("== ENV " + std::to_string(id) + " '" + er.name + "' ==");
        if (er.has_priv)
          o.line("  ==> Resultats concomitants des extremes de la composante privilegiee");
        tabulate(er.bounds.hi.value, "ENV " + std::to_string(id) + " MAX", nullptr);
        tabulate(er.bounds.lo.value, "ENV " + std::to_string(id) + " MIN", nullptr);
        any = true;
      }
      continue;
    }
    if (kind == "CABLE" && st1_ed) continue;  // the tables go with their results
    if (kind == "SURCH" || kind == "HISTORIQUE" || kind == "CABLE") {
      // SURCH / HISTORIQUE / CABLE editions belong to the domain extension;
      // what the core can honour is the SURCH envelope pair.
      if (kind != "SURCH") {
        warn("ST1: RESU " + kind + " — edition assuree par l'extension de domaine");
        continue;
      }
      std::vector<Id> ids = ids_in;
      if (ids.empty())
        for (const auto& [id, r] : surch_results) ids.push_back(id);
      for (Id id : ids) {
        if (!surch_results.count(id)) {
          warn("ST1: RESU — SURCH " + std::to_string(id) + " non calculee");
          continue;
        }
        const Envelope& e = surch_bounds(id);
        if (st1_ed) {
          auto t = surch_titles.find(id);
          use_cara_of("SURCH", id);
          auto sc = surch_contr.find(id);
          env_contr = sc != surch_contr.end() ? &sc->second : nullptr;
          auto sg = surch_contrg.find(id);
          env_contrg = sg != surch_contrg.end() ? &sg->second : nullptr;
          contrg_cst = false;
          {
            const std::set<Id> act = appuis_of("SURCH", id);
            const std::map<Id, Vec6> lo = appui_reactions(e.lo.value, act),
                                     hi = appui_reactions(e.hi.value, act);
            appui_tmp.clear();
            for (const auto& [aid, v] : lo) appui_tmp[aid] = {v, hi.count(aid) ? hi.at(aid) : v};
            appui_lohi = &appui_tmp;
          }
          st1_block(kSt1Label.at("SURCH"), id, t != surch_titles.end() ? t->second : std::string(),
                    e.lo.value, &e.hi.value);
          env_contr = nullptr;
          env_contrg = nullptr;
          contrg_cst = true;
          appui_lohi = nullptr;
          any = true;
          continue;
        }
        tabulate(e.hi.value, "SURCH " + std::to_string(id) + " MAX", nullptr);
        tabulate(e.lo.value, "SURCH " + std::to_string(id) + " MIN", nullptr);
        any = true;
      }
      continue;
    }
    const std::map<Id, CaseResult>* m = nullptr;
    if (kind == "CHARG") m = &store_.charg;
    else if (kind == "COMB") m = &store_.comb;
    else if (kind == "ETAT") m = &store_.etat;
    else if (kind == "PHASAGE") m = &store_.phasage;
    else if (kind == "SPECTRE") m = &store_.spectre;
    else if (kind == "INCENDIE") m = &store_.incendie;
    if (!m) continue;
    std::vector<Id> ids = ids_in;
    if (ids.empty())
      for (auto& [id, r] : *m) ids.push_back(id);
    for (Id id : ids) {
      auto it = m->find(id);
      if (it == m->end()) {
        warn("ST1: RESU — " + kind + " " + std::to_string(id) + " non calcule");
        continue;
      }
      const std::vector<PseBarSoil>* soil = nullptr;
      auto ps = pse_soil.find(id);
      if (kind == "CHARG" && ps != pse_soil.end()) soil = &ps->second;
      if (st1_ed && kSt1Label.count(kind)) {
        use_cara_of(kind, id);
        auto ct = cable_tables.find({kind, id});
        cable_tabs = ct != cable_tables.end() ? &ct->second : nullptr;
        cable_var = kind != "CHARG";
        appui_single(it->second, kind, id);
        st1_soil = soil;
        st1_block(kSt1Label.at(kind), id, st1_title(it->second.name, kind, id), it->second, nullptr);
        st1_soil = nullptr;
        cable_tabs = nullptr;
        appui_lohi = nullptr;
      }
      else
        tabulate(it->second, kind + " " + std::to_string(id), soil);
      any = true;
    }
  }
  if (!any && req.kinds.empty()) {
    // RESU with no result kind: everything in the loading store.
    for (auto& [id, r] : store_.charg) {
      const std::vector<PseBarSoil>* soil = nullptr;
      auto ps = pse_soil.find(id);
      if (ps != pse_soil.end()) soil = &ps->second;
      if (st1_ed) {
        use_cara_of("CHARG", id);
        auto ct = cable_tables.find({"CHARG", id});
        cable_tabs = ct != cable_tables.end() ? &ct->second : nullptr;
        cable_var = false;
        appui_single(r, "CHARG", id);
        st1_soil = soil;
        st1_block(kSt1Label.at("CHARG"), id, st1_title(r.name, "CHARG", id), r, nullptr);
        st1_soil = nullptr;
        cable_tabs = nullptr;
        appui_lohi = nullptr;
      }
      else
        tabulate(r, "CHARG " + std::to_string(id), soil);
    }
  }
  if (st1_blocks) o.line("");
}

// -------------------------------------------------------------- the model --
// ST1's default G_DYN (RO / G_DYN is the mass density; manual, CONS).
constexpr double kSt1DefaultGDyn = 9.81;

void St1Runtime::build_model() {
  if (built_) return;
  Model m;
  for (Id n : node_order_) m.nodes.push_back({n, nodes_.at(n)});

  std::map<CaraProps, Id> sec_ids;
  std::map<ConsProps, Id> mat_ids;
  for (Id bid : bar_order_) {
    if (selective_activation_ && !active_bars_.count(bid)) continue;
    const RawBar& rb = bars_.at(bid);
    require(cara_.count(bid), "ST1: barre " + std::to_string(bid) + " sans CARA");
    require(cons_.count(bid), "ST1: barre " + std::to_string(bid) + " sans CONS");
    const CaraProps& cp = cara_.at(bid);
    ConsProps mp = cons_.at(bid);
    if (mp.mat != 0) {
      // `CONS liste MAT i` (manual p.72): resolve the referenced material.
      auto mt = mat_props_.find(mp.mat);
      require(mt != mat_props_.end(),
              "ST1: barre " + std::to_string(bid) + " — MAT " + std::to_string(mp.mat) +
                  " non defini (bloc MAT/MATERIAU requis avant CONS)");
      const ConsProps& src = mt->second;
      if (mp.E == 0) mp.E = src.E;
      if (mp.ro == 0) mp.ro = src.ro;
      if (mp.alpha == 0) mp.alpha = src.alpha;
      if (mp.g_dyn == 0) mp.g_dyn = src.g_dyn;
      if (mp.nu == 0.3 && src.nu != 0.3) mp.nu = src.nu;
      if (mp.amortissement == 0) mp.amortissement = src.amortissement;
    }
    require(mp.E > 0, "ST1: barre " + std::to_string(bid) + " — module E non defini");
    if (!sec_ids.count(cp)) {
      Section s;
      s.name = "CARA" + std::to_string(sec_ids.size() + 1);
      s.shape = ShapeGeneric{};
      // GRILL carries no axial force (its DX DY RZ are blocked): SX is not
      // required there, a unit area stands in for it
      s.props.A = cp.sx > 0 || option_ != 1 ? cp.sx : 1.0;
      s.props.Avy = cp.sy;
      s.props.Avz = cp.sz;
      s.props.Iz = cp.iz > 0 ? cp.iz : cp.iy;
      s.props.Iy = cp.iy > 0 ? cp.iy : cp.iz;
      s.props.J = cp.ix > 0 ? cp.ix : std::max(s.props.Iy, s.props.Iz);
      s.props.Ip = cp.ipol;  // CARA IPOL: rotational mass inertia (p.60)
      // the CARA fibres VY WY VZ WZ as the section's four corners (an ENV
      // with privileged CONTR reads its fibres there)
      if (cp.vy > 0 || cp.wy > 0 || cp.vz > 0 || cp.wz > 0)
        for (double y : {cp.vy, -cp.wy})
          for (double z : {cp.vz, -cp.wz}) {
            StressPoint sp;
            sp.y = y;
            sp.z = z;
            s.stress_points.push_back(sp);
          }
      sec_ids[cp] = m.add_section(s);
    }
    if (!mat_ids.count(mp)) {
      Material mt;
      mt.name = "CONS" + std::to_string(mat_ids.size() + 1);
      mt.E = mp.E;
      mt.nu = mp.nu;
      // RO is a weight density: the mass density is RO / G_DYN, the
      // self-weight force RO x A is restored by lower_loads.
      const double gd = mp.g_dyn > 0 ? mp.g_dyn : kSt1DefaultGDyn;
      mt.rho = mp.ro / gd;
      mt.alpha = mp.alpha != 0 ? mp.alpha : 1.2e-5;
      mat_ids[mp] = m.add_material(mt);
    }
    Bar b;
    b.id = bid;
    b.ni = rb.n1;
    b.nj = rb.n2;
    b.material = mat_ids[mp];
    b.section = sec_ids[cp];
    if (exc_.count(bid)) b.ecc = exc_.at(bid);
    if (art_.count(bid)) b.releases = art_.at(bid);
    if (beta_.count(bid)) b.roll = beta_.at(bid);
    m.bars.push_back(b);
  }

  // Supports: user APPUI + option auto-blocked DOFs on every node.
  //
  // ST1 lets SEVERAL supports sit on the SAME node — that is the documented
  // point of numbering supports apart from nodes (manual p.48, p.51: "elle
  // permet aussi la definition de deux types d'appui sur un meme noeud").
  // keta carries one `Support` per node, so the descriptions are MERGED here:
  // blocked DOFs are unioned, diagonal springs and full 6x6 matrices are
  // summed (both are additive stiffness paths), and the DECOL switch of the
  // last description that carries one wins.
  std::array<bool, 6> auto_block{};
  if (option_ == 0) auto_block = {false, false, true, true, true, false};
  else if (option_ == 1) auto_block = {true, true, false, false, false, true};
  std::vector<Id> support_order;
  std::map<Id, Support> merged;
  for (Id aid : appui_order_) {
    if (selective_activation_ && !active_appuis_.count(aid)) continue;
    const SupportDef& s = appui_.at(aid);
    std::array<bool, 6> fix = auto_block;
    for (int c = 0; c < 6; ++c)
      fix[static_cast<std::size_t>(c)] =
          fix[static_cast<std::size_t>(c)] || s.fix[static_cast<std::size_t>(c)];
    // ST1 `APPUI ... (ROT alpha) EL m1 … m21` / `EL DI d1 … dj`: the FULL
    // symmetric support matrix, built by springmass.h from the row-by-row
    // lower triangle (manual p.49-51) and kept verbatim in the support frame —
    // `Support::rot` carries alpha and assembly rotates it to the global frame.
    std::optional<Mat6> k6;
    if (!s.spring_lt.empty()) {
      bool plane = option_ != 2;
      Mat6 kk = s.spring_di ? spring_from_diagonal(s.spring_lt, plane)
                            : spring_from_lower_triangle(s.spring_lt, plane);
      if (option_ == 1) kk = to_grill_frame(kk);  // GRILL order rx, dz, ry
      k6 = kk;
    }
    bool any = false;
    for (bool f : fix) any |= f;
    for (double k : s.spring) any |= (k != 0);
    if (k6 && k6->squaredNorm() > 0) any = true;
    if (!any) continue;
    if (!merged.count(s.node)) {
      support_order.push_back(s.node);
      merged[s.node].node = s.node;
    }
    Support& sup = merged[s.node];
    for (int c = 0; c < 6; ++c) {
      sup.fixed[static_cast<std::size_t>(c)] =
          sup.fixed[static_cast<std::size_t>(c)] || fix[static_cast<std::size_t>(c)];
      sup.spring[static_cast<std::size_t>(c)] += s.spring[static_cast<std::size_t>(c)];
    }
    if (k6 && k6->squaredNorm() > 0) {
      require(!sup.spring6 || sup.rot == s.rot,
              "ST1: noeud " + std::to_string(s.node) +
                  " — deux appuis elastiques de rotations ROT differentes ne sont pas cumulables");
      sup.spring6 = sup.spring6 ? Mat6(*sup.spring6 + *k6) : *k6;
      sup.rot = s.rot;
    }
    if (s.decol != 0) {
      // `DECOL <X,Y,Z> <POS,NEG>`: the criterion bears on ONE global axis and
      // releases the WHOLE support (manual p.47, p.49) — springmass.h's
      // all-or-nothing switch. A bare `DECOL` names no axis; keta then keeps
      // the per-DOF one-sided release of the blocked DOFs (documented).
      if (s.decol_axis >= 0) {
        sup.decol_axis = s.decol_axis;
        sup.decol_sign = s.decol;
      } else {
        for (int c = 0; c < 6; ++c)
          if (s.fix[static_cast<std::size_t>(c)])
            sup.oneside[static_cast<std::size_t>(c)] = s.decol;
      }
    }
  }
  for (Id nid : support_order) m.supports.push_back(merged.at(nid));
  if (option_ != 2) {
    for (const Node& n : m.nodes) {
      if (merged.count(n.id)) continue;
      m.add_support(n.id, auto_block, {});
      m.supports.back().implicit = true;
    }
  }
  // ST1 MASSE PROPRE (manual p.69-71): nodal mass matrices + per-bar line
  // masses / MASSE NUL, applied verbatim (springmass.h scatters them).
  for (const NodalMass& nm : masse_.nodal) {
    if (!nodes_.count(nm.node)) {
      warn("ST1: MASSE PROPRE — noeud " + std::to_string(nm.node) + " inconnu, ignore");
      continue;
    }
    m.nodal_masses.push_back(nm);
  }
  m.bar_mass_selective = masse_.selective;
  // the bars of the model (complete from here on): id -> FIRST position, for
  // the per-bar side tables below (a scan of m.bars per entry was O(n^2))
  std::unordered_map<Id, std::size_t> bar_pos;
  for (std::size_t k = 0; k < m.bars.size(); ++k) bar_pos.emplace(m.bars[k].id, k);
  for (const BarMass& bm : masse_.bars) {
    const bool active = bar_pos.count(bm.bar) > 0;
    if (!active) continue;  // MASSE on a de-activated bar: out of scope
    m.bar_masses.push_back(bm);
  }

  // ST1 CARA VAR (manual p.59-64): exact tapered member kernels installed on
  // the model side-table (see `TaperedBar` in model.h). Built here because the
  // element length and the material are only known once the model is complete.
  for (const auto& [bid, d] : cara_var_) {
    const bool active = bar_pos.count(bid) > 0;
    if (!active) continue;
    require(d.has_or && d.has_ex,
            "ST1: CARA VAR barre " + std::to_string(bid) + " — stations OR et EX requises");
    require(d.law != TaperLaw::Para || d.has_mi,
            "ST1: CARA VAR PARA barre " + std::to_string(bid) + " — station MI requise");
    TaperedVar var;
    var.law = d.law;
    var.orig = cara_to_props(d.orig);
    var.mid = cara_to_props(d.mid);
    var.ext = cara_to_props(d.ext);
    install_tapered(m, bid, var);
  }

  // ST1 CARA PSE (manual p.60-61): the exact beam on elastic foundation of
  // the note RDM §6 (bef.h), the formulation ST1 itself uses, on the model
  // side-table so that every analysis (EXEC CHARG, SURCH influence lines,
  // modes...) sees the soil. The pseudo-plastic driver (EXEC CHARG PSE)
  // keeps its own zoned bed.
  bef_bars_.clear();
  for (const auto& [bid, bs] : soil_zones_) {
    const auto pos = bar_pos.find(bid);
    if (pos == bar_pos.end() || bs.zones.empty()) continue;
    const auto it = m.bars.begin() + static_cast<std::ptrdiff_t>(pos->second);
    std::string why;
    if (m.tapered.count(bid)) why = "CARA VAR";
    for (bool r : it->releases)
      if (r) why = "liberation";
    if (!why.empty()) {
      warn("ST1: CARA PSE barre " + std::to_string(bid) + " (" + why +
           ") — lit de Hermite par zones au lieu de l'element exact sur sol elastique");
      continue;
    }
    std::vector<BefZone> zs;
    for (const SoilZone& z : bs.zones) zs.push_back({z.xl_end, z.kfx, z.kfy, z.kfz, z.kmx});
    install_bef(m, bid, zs);
    bef_bars_.insert(bid);
  }

  m.validate();
  model_ = std::move(m);
  built_ = true;
}

std::vector<Load> St1Runtime::self_weight_loads(const std::vector<Id>& list, const Vec3& dir) {
  std::vector<Load> out;
  // No list (or every active bar): the model-wide self-weight. A partial list
  // is one self-weight per bar, tagged with its bar: RO and SX are resolved
  // when the case is lowered (lower_loads) — a RO given in a MAT block only
  // reaches the CONS table when the model is built, after this CHARG is read.
  const std::vector<Id> active = all_bar_ids();
  if (list.empty() || std::set<Id>(list.begin(), list.end()) ==
                          std::set<Id>(active.begin(), active.end())) {
    out.push_back(Load::self_weight(dir * opt_.gravity));
    return out;
  }
  for (Id b : list) {
    require(cara_.count(b) && cons_.count(b),
            "ST1: POIDS PROPRE — CARA/CONS manquants pour la barre " + std::to_string(b));
    Load l = Load::self_weight(dir * opt_.gravity);
    l.bar = b;
    out.push_back(l);
  }
  return out;
}

Model& St1Runtime::snapshot() {
  build_model();
  return model_;
}

Model St1Runtime::snapshot_exact_bars(const std::set<Id>& bars) {
  const unsigned long gen = description_gen_;  // a snapshot describes nothing new
  const std::set<Id> saved = active_bars_;
  const bool saved_sel = selective_activation_;
  active_bars_.clear();
  for (Id b : bars)
    if (bars_.count(b)) active_bars_.insert(b);
  selective_activation_ = true;
  invalidate();
  build_model();
  Model m = model_;
  active_bars_ = saved;
  selective_activation_ = saved_sel;
  invalidate();
  build_model();  // lower_loads reads the materials of the current model
  description_gen_ = gen;
  return m;
}

std::vector<Load> St1Runtime::lower_loads(const std::vector<Load>& loads, const Model* target) {
  build_model();
  const Model& tm = target ? *target : model_;
  // G_DYN of a bar: its CONS (or MAT) value, ST1's 9.81 by default — read off
  // the description, so a target model built under another activation (a
  // PHASAGE base) lowers with its own bars.
  auto gdyn = [&](Id bid) {
    ConsProps cp = cons_.at(bid);
    if (cp.g_dyn == 0 && cp.mat != 0 && mat_props_.count(cp.mat)) cp.g_dyn = mat_props_.at(cp.mat).g_dyn;
    return cp.g_dyn > 0 ? cp.g_dyn : kSt1DefaultGDyn;
  };
  std::set<double> gds;
  for (const Bar& b : tm.bars) gds.insert(gdyn(b.id));
  std::vector<Load> out;
  out.reserve(loads.size());
  std::set<Id> present;
  for (const Bar& b : tm.bars) present.insert(b.id);
  for (const Load& l : loads) {
    // a bar load lands on the bars ACTIVE when the loading is executed
    // (p.112: "TOUT … agira sur toutes les barres activees au moment de
    // l'execution du chargement")
    const bool on_bar = l.bar != 0 && l.kind != Load::Kind::NodalForce &&
                        l.kind != Load::Kind::SupportDisp;
    if (on_bar && !present.count(l.bar)) continue;
    if (l.kind != Load::Kind::SelfWeight) {
      out.push_back(l);
      continue;
    }
    if (l.bar != 0) {
      // POIDS PROPRE of ONE listed bar (self_weight_loads): a span load
      // RO x SX(x) in the direction l.g (gravity factor included).
      const Id b = l.bar;
      ConsProps cp = cons_.at(b);
      if (cp.ro == 0 && cp.mat != 0 && mat_props_.count(cp.mat)) cp.ro = mat_props_.at(cp.mat).ro;
      auto v = cara_var_.find(b);
      if (v != cara_var_.end() && v->second.has_or && v->second.has_ex) {
        // CARA VAR: RO SX(x), linear under both LIN laws (exact), the
        // endpoint-linear approximation under PARA (as the model-wide
        // self-weight of build_model) — never the OR area along the bar.
        TaperedVar var;
        var.law = v->second.law;
        var.orig = cara_to_props(v->second.orig);
        var.mid = cara_to_props(v->second.mid);
        var.ext = cara_to_props(v->second.ext);
        out.push_back(Load::linear(b, l.g * (cp.ro * tapered_props_at(var, 0.0).A),
                                   l.g * (cp.ro * tapered_props_at(var, 1.0).A), LoadFrame::Global));
        continue;
      }
      out.push_back(Load::udl(b, l.g * (cp.ro * cara_.at(b).sx), LoadFrame::Global));
      continue;
    }
    if (gds.size() <= 1) {  // one G_DYN: rho * (G_DYN g) = RO g
      Load w = l;
      w.g *= gds.empty() ? 1.0 : *gds.begin();
      out.push_back(w);
      continue;
    }
    // Materials disagree on G_DYN: each bar weighs RO A = rho G_DYN A.
    for (const Bar& b : tm.bars) {
      require(!tm.tapered.count(b.id),
              "ST1: POIDS PROPRE — barre CARA VAR " + std::to_string(b.id) +
                  " dans un modele dont les materiaux ont des G_DYN differents (non tenu)");
      const double rhoA = tm.material(b.material).rho * gdyn(b.id) * tm.section(b.section).props.A;
      out.push_back(Load::udl(b.id, rhoA * l.g, LoadFrame::Global));
    }
  }
  return out;
}

// --------------------------------------------------------------- dispatch --
void St1Runtime::execute(const Statement& st) {
  if (finished || st.empty()) return;
  std::string head = st[0].str ? "" : upper(st[0].s);
  // Real ST1 accepts a few spellings of the same keyword; normalize the ones
  // the decks actually use so nothing is ever reported "non reconnue".
  if (head == "CHARGE") head = "CHARG";
  if (head == "TOUTES") head = "TOUT";

  // Block-content modes.
  if (mode_ == Mode::Nodes || mode_ == Mode::Bars) {
    if (is_uid_tok(st[0]) || head == "TOUT") {
      Statement full;
      full.push_back({mode_ == Mode::Nodes ? "NOEUD" : "BARRE", false});
      for (const Tok& t : st) full.push_back(t);
      if (mode_ == Mode::Nodes) cmd_noeud(full);
      else cmd_barre(full);
      return;
    }
    mode_ = Mode::Top;  // fall through to top-level dispatch
  }
  if (mode_ == Mode::Charg || mode_ == Mode::ChargBars || mode_ == Mode::ChargNodes ||
      mode_ == Mode::ChargDefor || mode_ == Mode::ChargTemp || mode_ == Mode::ChargAppui) {
    // Keywords that end sub-list modes but stay inside the CHARG block.
    if (mode_ != Mode::Charg && !is_uid_tok(st[0]) && head != "TOUT") {
      if (head == "NOEUD" && st.size() == 1) {
        mode_ = Mode::ChargNodes;
        return;
      }
      if (head == "BARRE" && st.size() == 1) {
        mode_ = Mode::ChargBars;
        return;
      }
      mode_ = Mode::Charg;
    }
    // A top-level command that can never be a load instruction closes a
    // CHARG block left without FIN: ST1 accepts `CHARG 20 '…' / TEMP TOUT
    // UNI +32 / EXEC CHARG 20` (the arc reference deck) and reports CHARG 20.
    static const std::set<std::string> kCloses = {"EXEC", "CHARG", "ENV",  "COMB",
                                                  "SURCH", "RESU", "SORTIE"};
    if (open_charg_ && !open_charg_->is_comb && kCloses.count(head)) {
      open_charg_ = nullptr;
      mode_ = Mode::Top;
      execute(st);
      return;
    }
    charg_stmt(st);
    return;
  }
  if (mode_ == Mode::Comb) {
    comb_stmt(st);
    return;
  }
  if (mode_ == Mode::Etude) {
    etude_stmt(st);
    return;
  }
  if (mode_ == Mode::Cara) {
    cara_stmt(st);
    return;
  }
  if (mode_ == Mode::Exc) {
    if (is_uid_tok(st[0]) || (!st[0].str && (upper(st[0].s) == "TOUT" ||
                                             upper(st[0].s) == "TOUTES"))) {
      Statement full;
      full.push_back({"EXC", false});
      for (const Tok& t : st) full.push_back(t);
      // Re-enter without disturbing the open block.
      Mode saved = mode_;
      mode_ = Mode::Top;
      cmd_exc(full);
      mode_ = saved;
      return;
    }
    mode_ = Mode::Top;  // any keyword line closes the block; re-dispatch
    if (!(st.size() == 1 && !st[0].str && upper(st[0].s) == "FIN")) execute(st);
    return;
  }
  if (mode_ == Mode::Cons) {
    cons_stmt(st);
    return;
  }
  if (mode_ == Mode::Appui) {
    appui_stmt(st);
    return;
  }
  if (mode_ == Mode::CaraVar) {
    cara_var_stmt(st);
    return;
  }
  if (mode_ == Mode::CaraPse) {
    cara_pse_stmt(st);
    return;
  }
  if (mode_ == Mode::CaraSpecial) {
    cara_special_stmt(st);
    return;
  }
  if (mode_ == Mode::Masse) {
    masse_stmt(st);
    return;
  }
  if (mode_ == Mode::Contrg) {
    contrg_stmt(st);
    return;
  }
  if (mode_ == Mode::Env) {
    env_stmt(st);
    return;
  }
  if (mode_ == Mode::Mat) {
    mat_stmt(st);
    return;
  }
  if (mode_ == Mode::Pagination) {
    if (head == "FIN") mode_ = Mode::Top;
    return;  // layout parameters recorded as no-ops (headless)
  }
  if (mode_ == Mode::Resu) {
    cmd_resu_line(st);
    return;
  }

  // Domain extension: gets first refusal on everything the core block modes
  // did not claim (its own FIN-blocks stay coherent because it is consulted
  // on every statement while one is open).
  if (domain_handler && domain_handler(st)) return;

  // Pending unwired block: swallow its body so its content lines never leak
  // into the top-level dispatch (a TABLIER's 'BARRE listeb' must not create
  // bars). FIN closes it.
  if (!pending_block_.empty()) {
    if (head == "FIN") pending_block_.clear();
    return;
  }

  // Top level.
  if (head == "OPTION") {
    std::string o = st.size() > 1 ? upper(st[1].s) : "SPATIALE";
    option_ = o == "PLANE" ? 0 : (o == "GRILL" ? 1 : 2);
    invalidate();
  } else if (head == "TITRE") {
    title_ = st.size() > 1 && st[1].str ? st[1].s : "";
  } else if (head == "NOEUD") cmd_noeud(st);
  else if (head == "GENER") cmd_gener(st);
  else if (head == "BARRE") cmd_barre(st);
  else if (head == "APPUI") cmd_appui(st);
  else if (head == "EXC") cmd_exc(st);
  else if (head == "ART") cmd_art(st);
  else if (head == "BETA") cmd_beta(st);
  else if (head == "CARA") cmd_cara(st);
  else if (head == "CONS") cmd_cons(st);
  else if (head == "CHARG") cmd_charg_open(st);
  else if (head == "COMB") cmd_comb(st);
  else if (head == "ETUDE") cmd_etude(st);
  else if (head == "CONTRG") cmd_contrg(st);
  else if (head == "ENV") cmd_env(st);
  else if (head == "MASSE") cmd_masse(st);
  else if (head == "ACTIVER") cmd_activer(st, true);
  else if (head == "DESACTIVER") cmd_activer(st, false);
  else if (head == "EXEC") cmd_exec(st);
  else if (head == "GET") cmd_get(st);
  else if (head == "SORTIE") cmd_sortie(st);
  else if (head == "MESSAGE") {
    std::string lvl = st.size() > 1 ? upper(st[1].s) : "TOUT";
    out_.message_level = lvl == "AUCUN" ? 0 : (lvl == "REDUIT" ? 1 : 2);
  } else if (head == "PAGINATION") {
    mode_ = Mode::Pagination;
  } else if (head == "ENREGISTRER") cmd_projet(st, true);
  else if (head == "PROJET") cmd_projet(st, false);
  else if (head == "LISTER") cmd_lister(st);
  else if (head == "RESU") {
    resu_req_.clear();
    mode_ = Mode::Resu;
    if (st.size() > 1) {
      Statement rest(st.begin() + 1, st.end());
      cmd_resu_line(rest);
    }
  } else if (head == "QUITTER" || head == "QUIT") {
    finished = true;
  } else if (head == "FIN") {
    mode_ = Mode::Top;  // stray FIN
  } else if (head == "MAT" || head == "MATERIAU") {
    // Fallback material reader — the domain extension had first refusal above
    // (see `mat_stmt`).
    require(st.size() >= 2 && is_uid_tok(st[1]), "ST1: MAT i ('titre')");
    open_mat_ = static_cast<int>(num(st[1]));
    mat_props_[open_mat_] = ConsProps{};
    mode_ = Mode::Mat;
  } else if (head == "PARA_STRUCTURE" || head == "MX_ETUDE") {
    // Memory sizing hints: keta allocates dynamically — accepted, no-op.
  } else if (head == "PREC" || head == "CABLE" || head == "TABLIER" || head == "SURCH" ||
             head == "EUROCODE" || head == "ALG" || head == "CONV" || head == "CONVOI" ||
             head == "PHASAGE" || head == "AMORTISSEMENT" || head == "SPECTRE" ||
             head == "ACCELEROGRAMME" || head == "ACCELEROGRAM" || head == "DYNAMIQUE" ||
             head == "EXPOSITION" || head == "SECTION" || head == "ACIER" ||
             head == "ARMATURE" || head == "FEU" || head == "INCENDIE" ||
             head == "ANALYSE" || head == "DESS") {
    warn("ST1: " + head + " — famille reconnue, cablage moteur a l'integration");
    pending_block_ = head;
  } else {
    warn("ST1: instruction '" + st[0].s + "' non reconnue");
  }
}

// ---------------------------------------------------------------- session --
St1Session::St1Session(const St1Options& opt) {
  runtime = std::make_unique<St1Runtime>(env, opt);
  interp = std::make_unique<Interp>(env, [this](const Statement& st) {
    const Interp::Where at = interp->where();
    if (on_statement) on_statement(st, at, false);
    runtime->execute(st);
    if (warning_where.size() < runtime->warnings().size())
      warning_where.resize(runtime->warnings().size(), at);
    if (on_statement) on_statement(st, at, true);
  });
  interp->set_echo([this](const std::string& l) { runtime->output().line(l); });
  // Full-language sessions carry the domain extension (MAT, PREC/CABLE,
  // TABLIER/SURCH, PHASAGE, dynamics, fire, DESS) by default.
  install_st1_domains(*runtime);
}

void St1Session::run(const std::string& text, const std::string& file) {
  try {
    interp->run(text, file);
  } catch (...) {
    // Warnings raised by the statement that failed belong to it.
    const std::vector<Interp::Where>& e = interp->error_where();
    if (warning_where.size() < runtime->warnings().size())
      warning_where.resize(runtime->warnings().size(), e.empty() ? Interp::Where{} : e[0]);
    throw;
  }
}

void ErreurLog::header(const std::string& deck) {
  sink_(" keta - execution du jeu ST1 " + deck + "\n" +
        "                      RAPPEL DES MESSAGES D'ERREUR\n" +
        " " + std::string(79, '-') + "\n \n");
}

void ErreurLog::step(const std::string& text) {
  sink_(" \n - " + text + " - \n \n");
  line_ = -1;  // the next message restates its line
}

void ErreurLog::elapsed(double seconds) {
  char buf[64];
  std::snprintf(buf, sizeof buf, "  Temps de calcul : %.2f s\n", seconds);
  sink_(buf);
}

void ErreurLog::message(const Interp::Where& at, bool error, const std::string& text) {
  std::string out;
  if (at.file != file_ || !file_known_) {
    out += " --- Fichier : >" + at.file + "\n";
    file_ = at.file;
    file_known_ = true;
    line_ = -1;
  }
  if (at.line != line_ && at.line > 0) {
    char buf[48];
    std::snprintf(buf, sizeof buf, " --- Ligne : %5d ---\n", at.line);
    out += buf;
  }
  line_ = at.line;
  std::string m = text.rfind("ST1: ", 0) == 0 ? text.substr(5) : text;
  for (char& c : m)
    if (c == '\n' || c == '\r') c = ' ';
  out += std::string(" ") + (error ? "erreur" : "attention") + " : " + m + "\n";
  sink_(out);
}

}  // namespace keta::st1

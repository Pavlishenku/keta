#include "keta/st1_domains.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <unordered_map>
#include <unordered_set>

#include "keta/assembly.h"
#include "keta/creep.h"
#include "keta/element.h"
#include "keta/etude.h"
#include "keta/influence.h"
#include "keta/prestress.h"

namespace keta::st1 {
namespace {

// --------------------------------------------------------------- lexical --
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
bool is_id_tok(const Tok& t) { return !t.str && is_integer(t.s); }
bool is_num_tok(const Tok& t) { return !t.str && is_st1_number(t.s); }
double num(const Tok& t) { return st1_num(t.s); }

// Upper-cased keyword at position i ("" when out of range or a quoted string).
std::string kw(const Statement& st, std::size_t i) {
  if (i >= st.size() || st[i].str) return std::string();
  return upper(st[i].s);
}
bool has_kw(const Statement& st, const char* k) {
  for (std::size_t i = 0; i < st.size(); ++i)
    if (kw(st, i) == k) return true;
  return false;
}
std::string first_string(const Statement& st, std::size_t from) {
  for (std::size_t i = from; i < st.size(); ++i)
    if (st[i].str) return st[i].s;
  return std::string();
}
std::string fmt(double v) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "%.6g", v);
  return buf;
}

// Consume a run of numbers (commas tolerated) starting at i.
std::vector<double> num_run(const Statement& st, std::size_t& i) {
  std::vector<double> v;
  while (i < st.size()) {
    if (st[i].s == ",") {
      ++i;
      continue;
    }
    if (!is_num_tok(st[i])) break;
    v.push_back(num(st[i]));
    ++i;
  }
  return v;
}

// --------------------------------------------------------------- mapping --
StrengthUnit unit_from(const std::string& k, bool& found) {
  found = true;
  if (k == "MPA") return StrengthUnit::MPa;
  if (k == "TM2") return StrengthUnit::TM2;
  if (k == "KNM2") return StrengthUnit::KNM2;
  found = false;
  return StrengthUnit::MPa;
}

bool law_from(const std::string& k, MatLaw::Code& out) {
  if (k == "BPEL") out = MatLaw::Code::BPEL;
  else if (k == "BHP1") out = MatLaw::Code::BHP1;
  else if (k == "BHP2") out = MatLaw::Code::BHP2;
  else if (k == "BPEL_BL") out = MatLaw::Code::BPEL_BL;
  else if (k == "CEB") out = MatLaw::Code::CEB;
  else if (k == "EC2") out = MatLaw::Code::EC2;
  else if (k == "EC2_BHP1") out = MatLaw::Code::EC2_BHP1;
  else if (k == "EC2_BHP2") out = MatLaw::Code::EC2_BHP2;
  else if (k == "FIB_2010" || k == "FIB2010") out = MatLaw::Code::FIB2010;
  else if (k == "ACI_B3") out = MatLaw::Code::ACI_B3;
  else if (k == "EC2_2027") out = MatLaw::Code::EC2_2027;
  else return false;
  return true;
}

bool cement_from(const std::string& k, CementClass& out) {
  // (the EC2 2027 classes are written CS / CN / CR in the manual's examples)
  if (k == "S" || k == "SL" || k == "CS") out = CementClass::S;
  else if (k == "N" || k == "CN") out = CementClass::N;
  else if (k == "R" || k == "CR") out = CementClass::R;
  else if (k == "RS") out = CementClass::RS;
  else return false;
  return true;
}

// EFFORT component names -> keta station order (N Vy Vz T My Mz).
int effort_comp(const std::string& k) {
  if (k == "N") return 0;
  if (k == "TY") return 1;
  if (k == "TZ") return 2;
  if (k == "MX") return 3;
  if (k == "MY") return 4;
  if (k == "MZ") return 5;
  return -1;
}
int dof_comp(const std::string& k) {
  if (k == "DX" || k == "FX") return 0;
  if (k == "DY" || k == "FY") return 1;
  if (k == "DZ" || k == "FZ") return 2;
  if (k == "RX" || k == "MX") return 3;
  if (k == "RY" || k == "MY") return 4;
  if (k == "RZ" || k == "MZ") return 5;
  return -1;
}

// EN 1991-2 Annex A special vehicles, ST1 keyword `VS_<total>_<par essieu>`
// (manual p.167, e.g. `VS_1800_200` = 1800 kN carried by 200 kN axles).  The
// annex fixes the axle spacing at 1,50 m and the vehicle width at 3,00 m; the
// two wheel lines are placed at +-1,00 m.  The weights are handed over in the
// deck's OWN force unit (`SurchDef::tonne_force` is 1 for ST1), so a kN deck
// gets kN axles — exactly what the keyword states.
bool vs_generic(const std::string& name, ConvoyDef& out) {
  if (name.rfind("VS_", 0) != 0) return false;
  std::size_t sep = name.find('_', 3);
  if (sep == std::string::npos) return false;
  std::string a = name.substr(3, sep - 3), b = name.substr(sep + 1);
  if (a.empty() || b.empty()) return false;
  for (char c : a + b)
    if (!std::isdigit(static_cast<unsigned char>(c))) return false;
  double total = std::stod(a), per_axle = std::stod(b);
  if (per_axle <= 0 || total < per_axle) return false;
  int n = static_cast<int>(std::llround(total / per_axle));
  if (n < 1 || n > 40) return false;
  out = ConvoyDef{};
  out.name = name;
  out.source = "EN 1991-2 annexe A (vehicule special), ST1 p.167 — entraxe 1,50 m, "
               "largeur 3,00 m, deux files de roues a +-1,00 m (reconstruction)";
  const std::vector<double> wheels = {-1.0, 1.0};
  for (int k = 0; k < n; ++k) {
    ConvoyAxle ax;
    ax.xl = 1.5 * k;
    ax.weight_t = per_axle;
    ax.yl = wheels;
    ax.impact_l = 0.15;
    ax.impact_t = 0.35;
    out.axles.push_back(ax);
  }
  out.longueur = 1.5 * (n - 1) + 1.5;
  out.largeur = 3.0;
  out.max_cam = 1;
  out.max_file = 1;
  out.trans = TransLaw::None;
  return true;
}

// EN 1991-2 §4.6.4 fatigue load model 3 (ST1 p.177, `LMF3` / `LMF3B`): four
// 120 kN axles at 0 / 1,20 / 7,20 / 8,40 m, wheels at +-1,00 m (contact
// 0,40 x 0,40 m), centred on its lane. LMF3B adds a second vehicle in the
// same lane with 36 kN axles, its centre 40 m behind (EN 1991-2 §4.6.4(3),
// the minimum spacing). Like the VS (p.167), the weights are handed over in
// the deck's own force unit.
bool lmf3_convoy(const std::string& name, ConvoyDef& out) {
  if (name != "LMF3" && name != "LMF3B") return false;
  out = ConvoyDef{};
  out.name = name;
  out.source = "EN 1991-2 4.6.4 (LMF3), ST1 p.177 — reconstruction";
  const std::vector<double> wheels = {-1.0, 1.0};
  auto vehicle = [&](double x0, double w) {
    for (double x : {0.0, 1.2, 7.2, 8.4}) {
      ConvoyAxle ax;
      ax.xl = x0 + x;
      ax.weight_t = w;
      ax.yl = wheels;
      ax.impact_l = 0.4;
      ax.impact_t = 0.4;
      out.axles.push_back(ax);
    }
  };
  vehicle(0.0, 120.0);
  if (name == "LMF3B") vehicle(40.0, 36.0);
  out.longueur = out.axles.back().xl + 0.4;
  out.largeur = 3.0;
  out.max_cam = 1;
  out.max_file = 1;
  out.axis_ecc = 0.0;  // centred on its lane (p.177)
  out.trans = TransLaw::None;
  return true;
}

// French national annex earth-moving plant (ST1 p.178, `CH_23M3` / `CH_34M3`).
// The manual gives them as figures only, lost to the text extraction: the
// three-axle layouts below are a RECONSTRUCTION (confidence LOW) fitted on
// ST1 2.24 envelopes (reactions, Ty and Mz of a
// 20 + 25 m two-span beam; residual 0.35 % / 0.55 % of the peak, against
// 1.5 % for the best two-axle layout). Weights in the deck's own force unit
// (e.g. kN), as the VS (p.167).
bool chantier_convoy(const std::string& name, ConvoyDef& out) {
  if (name != "CH_23M3" && name != "CH_34M3") return false;
  const bool big = name == "CH_34M3";
  out = ConvoyDef{};
  out.name = name;
  out.source = "NF EN 1991-2/NA (engin de terrassement), ST1 p.178 — reconstruction "
               "ajustee sur ST1 2.24 (couverture D05), confiance BASSE";
  const std::vector<double> wheels = {-1.2, 1.2};
  const std::array<std::pair<double, double>, 3> axles =
      big ? std::array<std::pair<double, double>, 3>{{{0.0, 470.0}, {1.0, 304.0}, {10.0, 728.0}}}
          : std::array<std::pair<double, double>, 3>{{{0.0, 363.0}, {0.5, 351.0}, {8.0, 491.0}}};
  for (const auto& [x, w] : axles) {
    ConvoyAxle ax;
    ax.xl = x;
    ax.weight_t = w;
    ax.yl = wheels;
    ax.impact_l = 0.5;
    ax.impact_t = 0.8;
    out.axles.push_back(ax);
  }
  out.longueur = axles.back().first + 0.5;
  out.largeur = 3.5;
  out.max_cam = 1;
  out.max_file = 1;
  out.axis_ecc = 0.0;
  out.trans = TransLaw::None;
  return true;
}

// Expand a free-DOF vector into per-node Vec6 (blocked DOFs read zero).
std::vector<Vec6> expand_free(const DofMap& d, const VecX& x) {
  std::vector<Vec6> out(static_cast<std::size_t>(d.n_nodes), Vec6::Zero());
  for (int n = 0; n < d.n_nodes; ++n)
    for (int c = 0; c < 6; ++c) {
      int f = d.free_of(n, c);
      if (f >= 0) out[static_cast<std::size_t>(n)](c) = x(f);
    }
  return out;
}

}  // namespace

// ===========================================================================
// SPECTRE — the two forms of the acceleration spectrum (manual p.198-200)
// ===========================================================================
ResponseSpectrum SpectrumDef::build(int dir) const {
  ResponseSpectrum s;
  if (!eurocode) {
    const std::vector<std::pair<double, double>>& t = dir == 0 ? ax : (dir == 1 ? ay : az);
    require(!t.empty(), "SPECTRE " + std::to_string(id) +
                            ": aucune table ACCELERATION pour la direction demandee (p.198)");
    s.kind = ResponseSpectrum::Kind::Table;
    s.table = t;
    for (auto& p : s.table) p.second *= pond;
    return s;
  }
  // EN 1998-1 §3.2.2.2/§3.2.2.5. ST1's `ST` (topography) multiplies ag, `S` is
  // the soil factor; `ALPHA`/`ETA` are the plateau amplification and the
  // damping correction, both built into keta's EC8 branch (2.5 and
  // sqrt(10/(5+100 zeta))) — a deviation is warned about by the parser.
  const Ec8Branch& b = ec8[static_cast<std::size_t>(std::clamp(dir, 0, 2))];
  require(b.set, "SPECTRE " + std::to_string(id) +
                     ": aucune branche EUROCODE pour la direction demandee "
                     "(HORIZONTAL / VERTICAL, p.199)");
  double a_g = b.ag * b.st_topo * pond;
  // ALPHA (2.5, vertical 3) is the plateau; ETA the damping correction,
  // 1.00 by default (p.199): ST1 does not derive it from the damping.
  if (reduit) {
    ResponseSpectrum r = ResponseSpectrum::ec8_design(a_g, b.s_soil, b.tb, b.tc, b.td, b.q);
    r.plateau = b.alpha;
    r.beta_floor = b.beta0;
    return r;
  }
  ResponseSpectrum e = ResponseSpectrum::ec8_elastic(a_g, b.s_soil, b.tb, b.tc, b.td);
  e.beta_floor = b.beta0;
  e.plateau = b.alpha;
  e.eta_fixed = b.eta > 0 ? b.eta : 1.0;
  return e;
}

// ===========================================================================
// Installation
// ===========================================================================
void install_st1_domains(St1Runtime& rt) {
  auto st = std::make_shared<St1Domains>(rt);
  rt.domain_state = st;
  St1Domains* raw = st.get();
  rt.domain_handler = [raw](const Statement& s) { return raw->handle(s); };
}

St1Domains* st1_domains(St1Runtime& rt) {
  return static_cast<St1Domains*>(rt.domain_state.get());
}

// ===========================================================================
// Dispatch
// ===========================================================================
bool St1Domains::handle(const Statement& st) {
  if (st.empty()) return false;
  const std::string head = st[0].str ? std::string() : upper(st[0].s);

  if (in_dess_) {
    dess_stmt(st, head);
    return true;
  }
  if (block_ != Block::None) {
    block_stmt(st, head);
    return true;
  }
  if (head == "EXEC") return exec_stmt(st);
  if (head == "GET" && kw(st, 1) == "TENSION") {
    get_tension(st);
    return true;
  }
  if (head == "LISTER" && lister(st, kw(st, 1))) return true;
  return open_block(st, head);
}

Model& St1Domains::model() {
  sync_materials();
  Model& m = rt_->snapshot();
  // MASSE PROPRE (p.69-71) is a description, not a builder map of the core:
  // re-apply it on every snapshot (the core rebuilds on any invalidation).
  if (masse.declared) {
    std::set<Id> bars;
    for (const Bar& b : m.bars) bars.insert(b.id);
    std::set<Id> nodes;
    for (const Node& n : m.nodes) nodes.insert(n.id);
    m.bar_masses.clear();
    for (const BarMass& bm : masse.bars)
      if (bars.count(bm.bar)) m.bar_masses.push_back(bm);
    m.nodal_masses.clear();
    for (const NodalMass& nm : masse.nodes)
      if (nodes.count(nm.node)) m.nodal_masses.push_back(nm);
  }
  return m;
}

// ===========================================================================
// MAT -> keta Material constants (manual p.72-73 CONS ... MAT i)
// ===========================================================================
// `E FEU EC2` assigns the tangent modulus at the origin at 20 degC of the
// EN 1992-1-2 3.2.2.1 law, E = 3 fck / (2 eps_c1,20) (manual p.219), not Ecm.
static double fire_cold_modulus(double fck) { return 1.5 * fck / concrete_eps_c1(20.0); }

void St1Domains::sync_materials() {
  bool changed = false;
  // A fire bar carries no CONS/CARA (manual p.224: "affectation au materiau
  // implicite"): synthesise both from its SECTION so the core snapshot builds.
  for (const auto& [bid, fb] : fire_bars) {
    auto fs = fire_sections.find(fb.section);
    if (fs == fire_sections.end()) continue;
    const FireSection& sec = fs->second.section;
    double A = 0, iy = 0, iz = 0;
    if (sec.kind == FireSection::Kind::Perimeter) {
      // Shoelace area + second moments of the polygon (y, z) about its centroid.
      const std::vector<Vec3>& p = sec.polygon;
      double a2 = 0, cy = 0, cz = 0;
      for (std::size_t k = 0; k + 0 < p.size(); ++k) {
        const Vec3& q0 = p[k];
        const Vec3& q1 = p[(k + 1) % p.size()];
        double cr = q0.x() * q1.y() - q1.x() * q0.y();
        a2 += cr;
        cy += (q0.x() + q1.x()) * cr;
        cz += (q0.y() + q1.y()) * cr;
      }
      A = std::fabs(a2) / 2;
      if (A > 0) {
        cy /= 3 * a2;
        cz /= 3 * a2;
        double ixx = 0, iyy = 0;
        for (std::size_t k = 0; k < p.size(); ++k) {
          double y0 = p[k].x() - cy, z0 = p[k].y() - cz;
          double y1 = p[(k + 1) % p.size()].x() - cy, z1 = p[(k + 1) % p.size()].y() - cz;
          double cr = y0 * z1 - y1 * z0;
          ixx += (z0 * z0 + z0 * z1 + z1 * z1) * cr;
          iyy += (y0 * y0 + y0 * y1 + y1 * y1) * cr;
        }
        iy = std::fabs(ixx) / 12;
        iz = std::fabs(iyy) / 12;
      }
    } else {
      // Thickness along local Y (SUP = VY, manual p220-223): Iz = w t^3 / 12
      // carries the bending of the slab.
      double t = sec.thickness, w = sec.width;
      A = t * w;
      iz = w * t * t * t / 12.0;
      iy = t * w * w * w / 12.0;
    }
    if (A <= 0) continue;
    if (!rt_->cara_.count(bid)) {
      CaraProps cp;
      cp.sx = A;
      // ST1 gives a SECTION-defined bar the gross area as its shear sections
      // (Exemple9: the cold static deflection of a 1 m slab over 3.6 m reads
      // ST1's -0.3096 mm only with SY = SX; the clamped-beam moments do not
      // depend on it, which is why the forces already matched).
      cp.sy = A;
      cp.sz = A;
      cp.iy = iy;
      cp.iz = iz;
      cp.ix = iy + iz;
      cp.vy = sec.thickness / 2;
      cp.wy = sec.thickness / 2;
      cp.vz = sec.width / 2;
      cp.wz = sec.width / 2;
      rt_->cara_[bid] = cp;
      changed = true;
    }
    if (!rt_->cons_.count(bid)) {
      ConsProps c;
      c.mat = fs->second.mat;
      auto mit = mats.find(fs->second.mat);
      c.E = (mit != mats.end() && mit->second.fire)
                ? fire_cold_modulus(mit->second.fire_mat.fck) / mit->second.stress_scale()
                : 35e9;
      c.nu = 0.2;
      c.ro = 0;
      rt_->cons_[bid] = c;
      changed = true;
    }
  }

  for (auto& [bid, cp] : rt_->cons_) {
    (void)bid;
    if (cp.mat == 0) continue;
    auto it = mats.find(cp.mat);
    if (it == mats.end()) {
      warn("ST1: CONS MAT " + std::to_string(cp.mat) + " — materiau non defini (p.73)");
      continue;
    }
    const MatDef& d = it->second;
    // Back to the deck's unit system (MatDef::stress_scale) — a Direct
    // 'E value' law is already in user units and must NOT be rescaled.
    double us = d.law.code == MatLaw::Code::Direct && !d.fire ? 1.0 : d.stress_scale();
    double E = (d.fire ? fire_cold_modulus(d.fire_mat.fck) : d.law.E(28.0)) / us;
    if (cp.E != E) changed = true;
    cp.E = E;
    if (d.has_nu) cp.nu = d.nu;
    if (d.has_ro) cp.ro = d.ro;
    if (d.g_dyn != 0) cp.g_dyn = d.g_dyn;
    if (d.has_alpha) cp.alpha = d.alpha;
    if (d.amortissement != 0) cp.amortissement = d.amortissement;
    // Materialized: the core's build must not demand a MAT registry it does
    // not own (its typed 'MAT non defini' check remains the backstop for a
    // CONS MAT i with no MAT block, or with no domain installed). The MAT
    // identity survives in `bar_mat` for the phasing creep path and in
    // `mats` for LISTER.
    bar_mat[bid] = cp.mat;
    cp.mat_origin = cp.mat;
    cp.mat = 0;
    changed = true;
  }
  if (changed) rt_->invalidate();
}

// ===========================================================================
// Block openers
// ===========================================================================
bool St1Domains::open_block(const Statement& st, const std::string& head) {
  auto id_at = [&](std::size_t i, Id dflt) {
    return (i < st.size() && is_id_tok(st[i])) ? static_cast<Id>(num(st[i])) : dflt;
  };

  if (head == "MAT" || head == "MATERIAU") {
    MatDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    mats[d.id] = d;
    rt_->touch_description();
    block_ = Block::Mat;
    block_id_ = d.id;
    sub_ = 0;
    // A one-line `MAT i E ...` form is legal too: replay the non-title tail.
    std::size_t k = 2;
    while (k < st.size() && st[k].str) ++k;
    if (k < st.size()) {
      Statement tail(st.begin() + static_cast<long>(k), st.end());
      mat_stmt(tail, kw(tail, 0));
    }
    return true;
  }
  // MASSE PROPRE: owned by the core runtime (st1_run) — never claimed here.
  if (head == "PREC") {
    PrecDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    precs[d.id] = d;
    rt_->touch_description();
    block_ = Block::Prec;
    block_id_ = d.id;
    return true;
  }
  if (head == "CABLE") {
    CableDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    cables[d.id] = d;
    rt_->touch_description();
    block_ = Block::Cable;
    block_id_ = d.id;
    sub_ = 0;
    return true;
  }
  if (head == "TABLIER") {
    has_tablier = true;
    block_ = Block::Tablier;
    sub_ = 0;
    return true;
  }
  if (head == "SURCH") {
    SurchDef d;
    d.tonne_force = 1.0;  // p.123: the French catalogue stays in TONNES; POND converts
    d.ec = eurocode;
    Id id = id_at(1, 1);
    surchs[id] = d;
    rt_->surch_titles[id] = first_string(st, 1);
    block_ = Block::Surch;
    block_id_ = id;
    return true;
  }
  if (head == "ALG") {
    AlgDef d;
    d.name = first_string(st, 1);
    Id id = id_at(1, 1);
    algs[id] = d;
    block_ = Block::Alg;
    block_id_ = id;
    sub_ = 0;
    return true;
  }
  if (head == "CONV" || head == "CONVOI") {
    ConvoyDef d;
    d.name = first_string(st, 1);
    d.source = "ST1 CONV generalise (p.160-162)";
    Id id = id_at(1, 1);
    convs[id] = d;
    block_ = Block::Conv;
    block_id_ = id;
    sub_ = 0;
    return true;
  }
  if (head == "EUROCODE") {
    block_ = Block::Eurocode;
    return true;
  }
  if (head == "PHASAGE") {
    PhasageDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    phasages[d.id] = d;
    phas_id_ = d.id;
    phas_cur_ = Phase{};
    phas_appui_on_.clear();
    phas_appui_off_.clear();
    phas_appui_mods_.clear();
    phas_appui_defs_.clear();
    phas_env_on_.clear();
    phas_kind_ = 0;
    phas_time_ = 0;
    phas_date_step_ = false;
    phas_dirty_ = false;
    phas_charg_ = false;
    block_ = Block::Phasage;
    block_id_ = d.id;
    return true;
  }
  if (head == "AMORTISSEMENT") {
    DampingDef d;
    d.id = id_at(1, 1);
    dampings[d.id] = d;
    block_ = Block::Amortissement;
    block_id_ = d.id;
    return true;
  }
  if (head == "SPECTRE") {
    Id id = id_at(1, 1);
    if (has_kw(st, "REPONSE")) {
      SpectreResponseDef d;
      d.id = id;
      spectre_reps[id] = d;
      block_ = Block::SpectreRep;
    } else {
      SpectrumDef d;
      d.id = id;
      spectra[id] = d;
      block_ = Block::Spectre;
      sub_ = 0;
    }
    block_id_ = id;
    return true;
  }
  if (head == "ACCELEROGRAMME" || head == "ACCELEROGRAM") {
    AccelDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    accels[d.id] = d;
    block_ = Block::Accel;
    block_id_ = d.id;
    sub_ = 0;
    return true;
  }
  if (head == "DYNAMIQUE") {
    std::string what = kw(st, 1);
    Id id = id_at(2, 1);
    if (what == "CONVOI") {
      DynConvoiDef d;
      d.id = id;
      dyn_convois[id] = d;
      block_ = Block::DynConvoi;
    } else if (what == "CHARGE") {
      DynChargeDef d;
      d.id = id;
      dyn_charges[id] = d;
      block_ = Block::DynCharge;
    } else {
      warn("ST1: DYNAMIQUE " + what + " — forme inconnue (p.205-209)");
      return false;
    }
    block_id_ = id;
    return true;
  }
  if (head == "EXPOSITION") {
    // `EXPOSITION (FEU) i` (p.216). The FEU word is optional in the reference
    // decks. When the line carries the description too it is a ONE-LINER
    // (`EXPOSITION FEU 1 TEMP 20`) and no FIN follows; a bare header opens the
    // block.
    std::size_t k = 1;
    if (kw(st, k) == "FEU") ++k;
    Id id = id_at(k, 1);
    if (k < st.size() && is_id_tok(st[k])) ++k;
    exposures[id] = Exposure{};
    block_id_ = id;
    if (k < st.size()) {
      Block save = block_;
      block_ = Block::Exposition;
      Statement tail(st.begin() + static_cast<long>(k), st.end());
      exposition_stmt(tail, kw(tail, 0));
      block_ = save;  // one-liner: nothing stays open
      return true;
    }
    block_ = Block::Exposition;
    return true;
  }
  if (head == "SECTION") {
    FireSectionDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    d.expo_ids = {0, 0, 0, 0};
    fire_sections[d.id] = d;
    block_ = Block::FireSection;
    block_id_ = d.id;
    sub_ = 0;
    return true;
  }
  if (head == "ACIER") {
    Id id = id_at(1, 1);
    aciers[id] = RebarSteel{};
    block_ = Block::Acier;
    block_id_ = id;
    return true;
  }
  if (head == "ARMATURE") {
    ArmatureDef d;
    d.id = id_at(1, 1);
    armatures[d.id] = d;
    block_ = Block::Armature;
    block_id_ = d.id;
    return true;
  }
  if (head == "CARA") {
    // The FIRE form only (manual p.224): `CARA ib SECTION i (ARMATURE j)`, or
    // the block header `CARA ib` once a fire SECTION exists. Everything else is
    // the R.D.M. CARA of p.59-68 and belongs to the core.
    const bool listed = st.size() >= 2 && (is_id_tok(st[1]) || kw(st, 1) == "TOUT");  // (TOUT too)
    bool inline_form = st.size() >= 3 && listed && has_kw(st, "SECTION");
    bool header_form = st.size() == 2 && is_id_tok(st[1]) && !fire_sections.empty();
    if (!inline_form && !header_form) return false;
    std::size_t i = 1;
    cara_ids_ = rt_->id_list(st, i, rt_->all_bar_ids());
    require(!cara_ids_.empty(), "ST1: CARA (feu) sans identifiant de barre (p.224)");
    for (Id b : cara_ids_) fire_bars[b] = FireBarDef{b, 0, 0};
    block_id_ = cara_ids_.front();
    if (inline_form) {
      Statement tail(st.begin() + static_cast<long>(i), st.end());
      cara_feu_stmt(tail, kw(tail, 0));
      // Synthesize the CARA/CONS of the fire bars EAGERLY: a following
      // `CHARG ... POIDS PROPRE liste` lowers per-bar weights at parse time
      // and must find them (exemple9.st1).
      sync_materials();
      return true;  // consumed in place; a trailing FIN is harmless at top level
    }
    block_ = Block::CaraFeu;
    return true;
  }
  if (head == "FEU") {
    FeuDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    feux[d.id] = d;
    block_ = Block::Feu;
    block_id_ = d.id;
    return true;
  }
  if (head == "INCENDIE") {
    IncendieDef d;
    d.id = id_at(1, 1);
    d.title = first_string(st, 1);
    incendies[d.id] = d;
    block_ = Block::Incendie;
    block_id_ = d.id;
    return true;
  }
  if (head == "ANALYSE" && kw(st, 1) == "FEU") {
    block_ = Block::AnalyseFeu;
    return true;
  }
  if (head == "DESS") {
    in_dess_ = true;
    return true;
  }
  return false;
}

void St1Domains::close_block() {
  block_ = Block::None;
  block_id_ = 0;
  sub_ = 0;
}

void St1Domains::block_stmt(const Statement& st, const std::string& head) {
  if (block_ == Block::Phasage) {
    phasage_stmt(st, head);
    return;
  }
  if (head == "FIN") {
    close_block();
    return;
  }
  switch (block_) {
    // (a MAT / PREC / CABLE line describes: the generation advances)
    case Block::Mat: rt_->touch_description(); mat_stmt(st, head); break;
    case Block::Prec: rt_->touch_description(); prec_stmt(st, head); break;
    case Block::Cable: rt_->touch_description(); cable_stmt(st, head); break;
    case Block::Tablier: tablier_stmt(st, head); break;
    case Block::Surch: surch_stmt(st, head); break;
    case Block::Alg: alg_stmt(st, head); break;
    case Block::Conv: conv_stmt(st, head); break;
    case Block::Eurocode: eurocode_stmt(st, head); break;
    case Block::Amortissement: amortissement_stmt(st, head); break;
    case Block::Spectre: spectre_stmt(st, head); break;
    case Block::SpectreRep: spectre_rep_stmt(st, head); break;
    case Block::Accel: accel_stmt(st, head); break;
    case Block::DynConvoi: dyn_convoi_stmt(st, head); break;
    case Block::DynCharge: dyn_charge_stmt(st, head); break;
    case Block::Masse: masse_stmt(st, head); break;
    case Block::Exposition: exposition_stmt(st, head); break;
    case Block::FireSection: fire_section_stmt(st, head); break;
    case Block::Acier: acier_stmt(st, head); break;
    case Block::Armature: armature_stmt(st, head); break;
    case Block::CaraFeu: cara_feu_stmt(st, head); break;
    case Block::Feu: feu_stmt(st, head); break;
    case Block::Incendie: incendie_stmt(st, head); break;
    case Block::AnalyseFeu: analyse_feu_stmt(st, head); break;
    default: break;
  }
}

// ===========================================================================
// MAT (p.74-88) + the fire branch `E FEU EC2` (p.214-215)
// ===========================================================================
void St1Domains::mat_stmt(const Statement& st, const std::string& head) {
  (void)head;
  MatDef& d = mats[block_id_];
  std::size_t i = 0;
  auto val = [&](double dflt) {
    if (i + 1 < st.size() && is_num_tok(st[i + 1])) {
      double v = num(st[i + 1]);
      i += 2;
      return v;
    }
    ++i;
    return dflt;
  };
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k.empty()) {
      ++i;
      continue;
    }
    if (k == "E") {
      ++i;
      if (kw(st, i) == "FEU") {
        // E FEU EC2 (<MPA,TM2,KNM2>) FCK fck (TRACTION <OUI,NON>)
        d.fire = true;
        ++i;
        if (kw(st, i) == "EC2") ++i;
        bool ok = false;
        StrengthUnit u = unit_from(kw(st, i), ok);
        if (ok) {
          d.unit = u;
          ++i;
        }
        if (kw(st, i) == "FCK" || kw(st, i) == "FC28") {
          ++i;
          require(i < st.size() && is_num_tok(st[i]), "ST1: MAT E FEU — FCK attend une valeur");
          d.fire_mat.fck = strength_to_pa(num(st[i]), StrengthUnit::MPa);  // always MPa (p.75)
          ++i;
        }
        if (kw(st, i) == "TRACTION") {
          d.fire_mat.tension = kw(st, i + 1) == "OUI";
          i += 2;
        }
        continue;
      }
      if (i < st.size() && is_num_tok(st[i])) {
        d.law.code = MatLaw::Code::Direct;
        d.law.e_direct = num(st[i]);
        ++i;
        continue;
      }
      MatLaw::Code code;
      if (!law_from(kw(st, i), code)) {
        warn("ST1: MAT " + std::to_string(d.id) + " — loi E '" + kw(st, i) + "' inconnue (p.75-86)");
        ++i;
        continue;
      }
      d.law.code = code;
      ++i;
      bool ok = false;
      StrengthUnit u = unit_from(kw(st, i), ok);
      if (ok) {
        d.unit = u;
        ++i;
      }
      if (kw(st, i) == "FC28" || kw(st, i) == "FCK") {
        ++i;
        require(i < st.size() && is_num_tok(st[i]), "ST1: MAT E — FC28/FCK attend une valeur");
        // fc28 / fck are always in MPa; the unit tag only says in which unit
        // the MODULUS is expressed (manual p.76).
        d.law.fc = strength_to_pa(num(st[i]), StrengthUnit::MPa);
        ++i;
      }
      continue;
    }
    if (k == "RETRAIT") {
      d.has_retrait = true;
      ++i;
      MatLaw::Code c;
      if (law_from(kw(st, i), c)) {
        d.code_retrait = c;
        ++i;
      }
      if (kw(st, i) == "ER") {
        ++i;
        require(i < st.size() && is_num_tok(st[i]), "ST1: MAT RETRAIT ER attend une valeur");
        d.law.er = std::fabs(num(st[i]));
        ++i;
      }
      continue;
    }
    if (k == "FLUAGE") {
      d.has_fluage = true;
      ++i;
      MatLaw::Code c;
      if (law_from(kw(st, i), c)) {
        d.code_fluage = c;
        ++i;
      }
      while (i < st.size()) {
        std::string t = kw(st, i);
        if (t == "TRACTION") {
          d.creep_traction = kw(st, i + 1) != "NON";
          i += 2;
        } else if (t == "LINEAIRE") {
          d.creep_linear = kw(st, i + 1) != "NON";
          i += 2;
        } else break;
      }
      continue;
    }
    if (k == "AVANCE") {
      ++i;
      while (i + 1 < st.size() && is_num_tok(st[i + 1])) {
        std::string t = kw(st, i);
        double v = num(st[i + 1]);
        // (the short names are aliases, p.91)
        AvanceFactors& av = d.law.avance;
        if (t == "K_MODULE" || t == "KM") av.k_module = v;
        else if (t == "K_RETRAIT_ENDOGENE" || t == "KRE") av.k_retrait_endogene = v;
        else if (t == "K_RETRAIT_DESSICATION" || t == "KRD") av.k_retrait_dessication = v;
        else if (t == "K_FLUAGE_ENDOGENE" || t == "KFE") av.k_fluage_endogene = v;
        else if (t == "K_FLUAGE_DESSICATION" || t == "KFD") av.k_fluage_dessication = v;
        else if (t == "K_FLUAGE_CINETIQUE" || t == "KFC") av.k_fluage_cinetique = v;
        else if (t == "K_RETRAIT_ENDOGENE_CINETIQUE" || t == "KREC") av.k_retrait_endogene_cin = v;
        else if (t == "K_RETRAIT_DESSICATION_CINETIQUE" || t == "K_RETRAIT_DESSIC_CINETIQUE" ||
                 t == "KRDC")
          av.k_retrait_dessication_cin = v;
        else if (t == "K_FLUAGE_ENDOGENE_CINETIQUE" || t == "KFEC") av.k_fluage_endogene_cin = v;
        else if (t == "K_FLUAGE_DESSICATION_CINETIQUE" || t == "K_FLUAGE_DESSIC_CINETIQUE" ||
                 t == "KFDC")
          av.k_fluage_dessication_cin = v;
        else break;
        i += 2;
      }
      continue;
    }
    if (k == "SILICE") {  // EC2_2027 creep with silica fume (p.90), NON by default
      d.law.silica = kw(st, i + 1) == "OUI";
      i += 2;
      continue;
    }
    if (k == "CIMENT") {
      ++i;
      CementClass cc;
      if (i < st.size() && is_num_tok(st[i])) {
        d.law.b3_cement_type = static_cast<int>(num(st[i]));  // ACI_B3: 1/2/3
        ++i;
      } else if (cement_from(kw(st, i), cc)) {
        d.law.cement = cc;
        ++i;
      } else {
        warn("ST1: MAT CIMENT '" + kw(st, i) + "' inconnu (p.78-86)");
        ++i;
      }
      continue;
    }
    if (k == "CONDUCTIVITE") {
      ++i;
      if (kw(st, i) == "EC2") {
        ++i;
        d.fire_mat.thermal.cond = ThermalMaterial::Conductivity::Ec2;
        if (kw(st, i) == "AN") {
          d.fire_mat.thermal.lambda_alpha = kw(st, i + 1) == "BHP" ? 0.0 : 1.0;
          i += 2;
        } else if (i < st.size() && is_num_tok(st[i])) {
          d.fire_mat.thermal.lambda_alpha = num(st[i]);
          ++i;
        }
      } else if (i < st.size() && is_num_tok(st[i])) {
        d.fire_mat.thermal.cond = ThermalMaterial::Conductivity::Constant;
        d.fire_mat.thermal.lambda_const = num(st[i]);
        ++i;
      }
      continue;
    }
    if (k == "GRANULAT") {
      d.fire_mat.thermal.aggregate =
          kw(st, i + 1) == "CALCAIRE" ? Aggregate::Calcareous : Aggregate::Siliceous;
      i += 2;
      continue;
    }
    if (k == "RHO_FEU") {
      ++i;
      if (i < st.size() && is_num_tok(st[i])) {
        d.fire_mat.thermal.rho20 = num(st[i]);
        ++i;
      }
      if (kw(st, i) == "CTE") {
        d.fire_mat.thermal.rho_law = ThermalMaterial::DensityLaw::Constant;
        ++i;
      } else if (kw(st, i) == "EC2") ++i;
      continue;
    }
    // Plain `KEY value` parameters.
    if (k == "NU") { d.nu = val(d.nu); d.has_nu = true; continue; }
    if (k == "RO") { d.ro = val(d.ro); d.has_ro = true; continue; }
    if (k == "G_DYN") { d.g_dyn = val(d.g_dyn); continue; }
    if (k == "TEMP") { d.alpha = val(d.alpha); d.has_alpha = true; continue; }
    if (k == "AMORTISSEMENT") { d.amortissement = static_cast<int>(val(0)); continue; }
    if (k == "RH" || k == "ROH") { d.law.rh = val(d.law.rh); continue; }
    if (k == "H" || k == "RM_EC") { d.law.h_mm = val(d.law.h_mm); continue; }
    if (k == "RM") {
      d.law.rm_cm = val(d.law.rm_cm);
      d.law.h_mm = 20.0 * d.law.rm_cm;  // matlaws.h: h_mm = 20 rm_cm, same section
      continue;
    }
    if (k == "ROS") { d.law.ros = val(d.law.ros); continue; }
    if (k == "ROBS") { d.law.robs = val(d.law.robs); continue; }
    if (k == "EA") { d.law.ea = val(d.law.ea); continue; }
    if (k == "K_FORME") { d.law.k_forme = val(d.law.k_forme); continue; }
    if (k == "K_CURE") { d.law.k_cure = val(d.law.k_cure); continue; }
    if (k == "TEMPS_CURE") { d.law.ts_dry = val(d.law.ts_dry); continue; }
    if (k == "QUANTITE_EAU") { d.law.water = val(d.law.water); continue; }
    if (k == "QUANTITE_CIMENT") { d.law.cement_qty = val(d.law.cement_qty); continue; }
    if (k == "QUANTITE_AGREGAT") { d.law.aggregate = val(d.law.aggregate); continue; }
    if (k == "UU") { d.fire_mat.thermal.uu = val(0); continue; }
    if (k == "EMISSIVITE") { d.fire_mat.thermal.emissivity = val(0.7); continue; }
    if (k == "CHALEUR_SPECIFIQUE_MASSIQUE") { d.fire_mat.thermal.c_const = val(0); continue; }
    if (k == "T_THERM") {
      // T_THERM FCMP fcmp (LAMBDA l) DTP dt1,dt2.. THETA t1,t2.. (EC2, p.80)
      // or FCP fcp (LAMBDA l) DT tpa THETA tmax (BPEL, p.75)
      d.therm = true;
      ++i;
      auto list = [&](std::vector<double>& out) {
        out.clear();
        while (i < st.size() && (is_num_tok(st[i]) || st[i].s == ",")) {
          if (st[i].s != ",") out.push_back(num(st[i]));
          ++i;
        }
      };
      while (i < st.size()) {
        const std::string s = kw(st, i);
        if (s == "FCMP" || s == "FCP") {
          d.therm_fcmp = num(st[i + 1]);
          i += 2;
        } else if (s == "LAMBDA") {
          d.therm_lambda = num(st[i + 1]);
          i += 2;
        } else if (s == "DTP" || s == "DT") {
          ++i;
          list(d.therm_dt);
        } else if (s == "THETA") {
          ++i;
          list(d.therm_theta);
        } else {
          break;
        }
      }
      continue;
    }
    warn("ST1: MAT " + std::to_string(d.id) + " — mot '" + k + "' ignore");
    ++i;
  }
}

double MatDef::therm_end() const {
  double tp = 0;
  for (double dt : therm_dt) tp += dt;
  return tp;
}

double MatDef::therm_dtheta_max() const {
  double tmax = 20.0;
  for (double th : therm_theta) tmax = std::max(tmax, th);
  return tmax - 20.0;
}

// EN 1992-1-1 eq. 10.2 (ST1 annex B4.6): t_eb = t_pa 1.14^(theta_max - 20),
// t_pa = sum (theta_i - 20) dt_i / (theta_max - 20).
double MatDef::therm_teb() const {
  const double dmax = therm_dtheta_max();
  if (!therm || dmax <= 0) return 0.0;
  double s = 0;
  for (std::size_t k = 0; k < therm_dt.size() && k < therm_theta.size(); ++k)
    s += (therm_theta[k] - 20.0) * therm_dt[k];
  return s / dmax * std::pow(1.14, dmax);
}

// ===========================================================================
// PREC (p.89-90)
// ===========================================================================
void St1Domains::prec_stmt(const Statement& st, const std::string& head) {
  (void)head;
  PrecDef& d = precs[block_id_];
  std::size_t i = 0;
  auto val = [&](double dflt) {
    if (i + 1 < st.size() && is_num_tok(st[i + 1])) {
      double v = num(st[i + 1]);
      i += 2;
      return v;
    }
    ++i;
    return dflt;
  };
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k.empty()) {
      ++i;
      continue;
    }
    if (k == "SECTION") { d.ap = val(d.ap); continue; }
    if (k == "TENSION") { d.sigma = val(d.sigma); continue; }
    if (k == "E") { d.ep = val(d.ep); continue; }
    if (k == "RG" || k == "FPK") { d.fpk = val(d.fpk); continue; }
    if (k == "R1000") { d.rho1000 = val(d.rho1000); continue; }
    if (k == "RECUL") { d.recul = val(d.recul); continue; }
    if (k == "SEUIL") { d.seuil = val(0.0); continue; }
    if (k == "LOI") {
      const std::string law = kw(st, i + 1);
      d.loi_ec2_2027 = law == "EC2_2027";
      d.loi = law == "BPEL" ? 1 : (law == "EC2" ? 2 : (law == "EC2_2027" ? 3 : 0));
      if (d.loi == 0) warn("ST1: PREC LOI — loi '" + law + "' inconnue (BPEL, EC2, EC2_2027)");
      i += 2;
      if (kw(st, i) == "KP") d.kp = val(d.kp);
      continue;
    }
    if (k == "KP") { d.kp = val(d.kp); continue; }
    if (k == "PERTE") {
      ++i;
      while (i < st.size()) {
        std::string t = kw(st, i);
        if (t == "F") { d.f = val(0); }
        else if (t == "PHI") { d.phi = val(0); }
        else if (t == "MU") { d.f = val(0); d.eurocode_form = true; }
        else if (t == "KPHI") { d.kphi = val(0); d.eurocode_form = true; }
        else if (t == "INST") { d.perte_inst = val(0); }
        else break;
      }
      continue;
    }
    if (k == "PRETENSION") {
      d.pretension = true;
      ++i;
      while (i < st.size()) {
        std::string t = kw(st, i);
        if (t == "LINEAIRE") { d.pret_law = 1; ++i; }
        else if (t == "PARABOL") { d.pret_law = 2; ++i; }
        else if (t == "BILINEAIRE") { d.pret_law = 3; ++i; }
        else if (t == "LG_SCEL") {
          d.lg_scel1 = val(0);
          if (d.pret_law == 3) {
            --i;  // val() consumed the keyword slot: read lg2 right after lg1
            d.lg_scel2 = val(d.lg_scel1);
          }
        } else if (t == "COEF") { d.pret_alpha = val(1); }
        else break;
      }
      require(d.pret_law != 3 || (d.lg_scel2 > d.lg_scel1 && d.lg_scel1 > 0 &&
                                  d.pret_alpha >= 0 && d.pret_alpha <= 1),
              "ST1: PREC " + std::to_string(d.id) +
                  " PRETENSION BILINEAIRE — LG_SCEL lg1 lg2 (0 < lg1 < lg2), COEF 0..1 (p.93)");
      continue;
    }
    warn("ST1: PREC " + std::to_string(d.id) + " — mot '" + k + "' ignore");
    ++i;
  }
}

// ===========================================================================
// CABLE (p.91-94)
// ===========================================================================
void St1Domains::cable_stmt(const Statement& st, const std::string& head) {
  CableDef& d = cables[block_id_];
  std::size_t i = 0;
  if (sub_ == 0) {
    if (head == "PREC") {
      d.prec = static_cast<Id>(num(st[1]));
      for (std::size_t k = 2; k < st.size(); ++k) {
        if (kw(st, k) == "NON") d.non_simul = true;
        if (kw(st, k) == "SIMUL" && k + 1 < st.size() && is_num_tok(st[k + 1]))
          d.simul = static_cast<int>(num(st[k + 1]));
      }
      return;
    }
    if (head == "TENSION") {
      // OR / EX / OR EX / EX OR (p.95): the order of a pair is the order of
      // release of the two anchors.
      std::size_t por = 0, pex = 0;
      for (std::size_t k = 1; k < st.size(); ++k) {
        if (kw(st, k) == "OR" && !por) por = k;
        if (kw(st, k) == "EX" && !pex) pex = k;
      }
      d.anchor = (por && pex) ? (pex < por ? 3 : 2) : (pex ? 1 : 0);
      return;
    }
    if (head == "BARRE") {
      i = 1;
      d.bars = rt_->id_list(st, i, rt_->all_bar_ids());
      return;
    }
    if (head == "PAS_CABLE") {
      d.pas_cable = num(st[1]);
      return;
    }
    if (head == "LG_GAINE") {
      for (std::size_t k = 1; k + 1 < st.size(); ++k) {
        if (kw(st, k) == "OR") d.lg_gaine_or = num(st[k + 1]);
        if (kw(st, k) == "EX") d.lg_gaine_ex = num(st[k + 1]);
      }
      return;  // pretension: no bond over these lengths (build_tendon)
    }
    if (head == "TRACE") {
      sub_ = 1;
      i = 1;
      if (kw(st, i) == "INTERIEUR") ++i;
      else if (kw(st, i) == "EXTERIEUR") {
        d.exterieur = true;
        ++i;
      }
      if (i >= st.size()) return;
      // a point may follow on the same line
    } else {
      warn("ST1: CABLE " + std::to_string(d.id) + " — mot '" + head + "' ignore");
      return;
    }
  }
  // TRACE point line.
  TracePoint tp;
  int axis = 0;
  bool any = false;
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (st[i].s == ",") {
      ++i;
      continue;
    }
    if (k == "X" || k == "Y" || k == "Z") {
      axis = k == "X" ? 0 : (k == "Y" ? 1 : 2);
      ++i;
      require(i < st.size() && is_num_tok(st[i]), "ST1: CABLE TRACE — coordonnee attendue");
      tp.p(axis) = num(st[i]);
      ++axis;
      any = true;
      ++i;
      continue;
    }
    if (is_num_tok(st[i])) {
      require(axis < 3, "ST1: CABLE TRACE — plus de trois coordonnees");
      tp.p(axis) = num(st[i]);
      ++axis;
      any = true;
      ++i;
      continue;
    }
    if (k == "ALIGNE") { tp.aligne = true; ++i; continue; }
    if (k == "PENTE") { tp.pente = num(st[i + 1]); i += 2; continue; }
    if (k == "GIS") { tp.gis = num(st[i + 1]); i += 2; continue; }
    if (k == "RAYON") {  // a pole: circular fillet of that radius (p.98)
      tp.rayon = num(st[i + 1]);
      i += 2;
      continue;
    }
    if (k == "INTERIEUR") { ++i; continue; }
    if (k == "EXTERIEUR") { d.exterieur = true; ++i; continue; }
    warn("ST1: CABLE TRACE — mot '" + k + "' ignore");
    ++i;
  }
  if (any) d.trace.push_back(tp);
}

// ===========================================================================
// Tendon construction (CABLE trace -> keta Tendon)
// ===========================================================================
std::vector<St1Runtime::CableRow> St1Domains::cable_rows(
    const Model& m, const Tendon& t, double pas, const std::function<double(double)>& sigma_at,
    double sigma0, const std::function<double(double)>& ref) const {
  std::vector<St1Runtime::CableRow> rows;
  struct Frame {
    double s0 = 0, L = 0;
    Vec3 xi, xj;
    Mat3 R;
  };
  std::vector<Frame> fr;
  double s0 = 0;
  for (Id bid : t.bars) {
    const Bar& b = m.bar(bid);
    Frame f;
    f.s0 = s0;
    f.xi = m.node(b.ni).x;
    f.xj = m.node(b.nj).x;
    f.L = (f.xj - f.xi).norm();
    f.R = bar_rotation(f.xi, f.xj, b.roll);
    s0 += f.L;
    fr.push_back(f);
  }
  if (fr.empty() || s0 <= 0) return rows;
  auto pos = [&](double s) {
    std::size_t bi = 0;
    for (std::size_t j = 0; j < fr.size(); ++j)
      if (s >= fr[j].s0 - 1e-12) bi = j;
    const Frame& f = fr[bi];
    const double u = f.L > 0 ? std::clamp((s - f.s0) / f.L, 0.0, 1.0) : 0.0;
    return Vec3(f.xi + u * (f.xj - f.xi) + tendon_ey(t, s) * f.R.row(1).transpose() +
                tendon_ez(t, s) * f.R.row(2).transpose());
  };
  const double sb = std::clamp(t.s_begin, 0.0, s0);
  const double se = t.s_end < 0 ? s0 : std::clamp(t.s_end, sb, s0);
  const double step = 4.0 * (pas > 0 ? pas : 0.25);
  auto row = [&](double s) {
    const Vec3 p = pos(s);
    St1Runtime::CableRow r;
    r.x = p.x();
    r.y = p.y();
    r.z = p.z();
    r.sigma = sigma_at(s);
    r.dsigma = r.sigma - (ref ? ref(s) : sigma0);  // Var.Tension
    rows.push_back(r);
  };
  // arc length along the 3-D trace, sampled finely; rows on its grid
  const int n = std::max(8, static_cast<int>(std::ceil((se - sb) / std::min(0.05, step / 20))));
  double arc = 0, next = 0;
  Vec3 prev = pos(sb);
  row(sb);
  next = step;
  for (int i = 1; i <= n; ++i) {
    const double sa = sb + (se - sb) * (i - 1) / n, sc = sb + (se - sb) * i / n;
    const Vec3 cur = pos(sc);
    const double d = (cur - prev).norm();
    while (d > 0 && arc + d >= next - 1e-9 && next < arc + d - 1e-9 * step) {
      row(sa + (sc - sa) * (next - arc) / d);
      next += step;
    }
    arc += d;
    prev = cur;
  }
  row(se);
  return rows;
}

Tendon St1Domains::build_tendon(const Model& m, const CableDef& c, const PrecDef& p,
                                bool with_losses) {
  Tendon t;
  t.id = c.id;
  t.name = c.title.empty() ? ("CABLE " + std::to_string(c.id)) : c.title;
  t.exact_angle = true;  // ST1's external method: the force along the cable
  if (c.bars.empty())
    for (const Bar& b : m.bars) t.bars.push_back(b.id);
  else
    for (Id b : c.bars) t.bars.push_back(b);
  require(!t.bars.empty(), "ST1: CABLE " + std::to_string(c.id) + " — aucune barre support");

  // Chain geometry: cumulated chord length + per-bar local frame.
  struct Seg {
    Vec3 xi, xj;
    Mat3 R;
    double s0 = 0, L = 0;
  };
  std::vector<Seg> segs;
  double s = 0;
  Id prev_nj = 0;
  for (Id bid : t.bars) {
    const Bar& b = m.bar(bid);
    if (prev_nj != 0 && b.ni != prev_nj)
      warn("ST1: CABLE " + std::to_string(c.id) + " — la chaine de barres n'est pas continue "
           "a la barre " + std::to_string(bid) + " (p.92 : suite geometriquement continue)");
    prev_nj = b.nj;
    Seg g;
    g.xi = m.node(b.ni).x;
    g.xj = m.node(b.nj).x;
    g.L = (g.xj - g.xi).norm();
    g.R = bar_rotation(g.xi, g.xj, b.roll);
    g.s0 = s;
    s += g.L;
    segs.push_back(g);
  }
  double Lchain = s;

  // Project every trace point onto the chain: (s, e_y, e_z) in the bar frame.
  struct Pt {
    double s = 0, ey = 0, ez = 0;
    bool aligne = false;
    std::optional<double> gis, pente, rayon;
  };
  std::vector<Pt> pts;
  for (const TracePoint& tp : c.trace) {
    double best = 1e300;
    Pt out;
    for (std::size_t gi = 0; gi < segs.size(); ++gi) {
      const Seg& g = segs[gi];
      if (g.L <= 0) continue;
      Vec3 d = g.xj - g.xi;
      const double raw = (tp.p - g.xi).dot(d) / d.squaredNorm();
      double u = std::clamp(raw, 0.0, 1.0);
      double dist = (tp.p - (g.xi + u * d)).norm();
      if (dist < best) {
        best = dist;
        // a point beyond the chain ends (Fontenoy's renfort cables start
        // 0.6 m before the first bar) keeps its own abscissa on the
        // prolonged end bar, so the trace keeps its true length
        if ((gi == 0 && raw < 0) || (gi + 1 == segs.size() && raw > 1)) u = raw;
        Vec3 off = tp.p - (g.xi + u * d);
        out.s = g.s0 + u * g.L;
        out.ey = g.R.row(1).dot(off);
        out.ez = g.R.row(2).dot(off);
      }
    }
    out.aligne = tp.aligne;
    out.gis = tp.gis;
    out.pente = tp.pente;
    out.rayon = tp.rayon;
    pts.push_back(out);
  }
  for (std::size_t k = 1; k < pts.size(); ++k)
    if (pts[k].s < pts[k - 1].s - 1e-9) {
      warn("ST1: CABLE " + std::to_string(c.id) + " — points de passage non ordonnes le long "
           "des barres (point " + std::to_string(k + 1) + ") : keta les trie, alors que ST1 suit "
           "l'ordre de saisie et fait rebrousser le cable (deviation et frottement en plus)");
      break;
    }
  std::stable_sort(pts.begin(), pts.end(), [](const Pt& a, const Pt& b) { return a.s < b.s; });
  // Drop duplicates in s (a repeated abscissa would produce a null segment).
  std::vector<Pt> q;
  for (const Pt& x : pts)
    if (q.empty() || x.s > q.back().s + 1e-9) q.push_back(x);

  if (q.empty()) {
    PolySegment seg;
    seg.s0 = 0;
    seg.s1 = Lchain;
    t.profile.push_back(seg);  // concentric tendon
  } else if (q.size() == 1) {
    PolySegment seg;
    seg.s0 = 0;
    seg.s1 = Lchain;
    seg.cy = {q[0].ey, 0, 0, 0};
    seg.cz = {q[0].ez, 0, 0, 0};
    t.profile.push_back(seg);
  } else {
    // The trace ends ARE the anchors (p.92: the cable runs between its first
    // and last points): a post-tensioned cable may start or stop inside its
    // bars (Fontenoy: 1106 anchored 4 m into the span). The profile is still
    // extended flat to the chain ends below, carrying no force there.
    // Pretension keeps its bonded-length law over the whole chain.
    if (!p.pretension) {
      const double tol = 1e-6 * std::max(1.0, Lchain);
      if (q.front().s > tol) t.s_begin = q.front().s;
      if (q.back().s < Lchain - tol) t.s_end = q.back().s;
    }
    // The trace is extended flat to the chain ends (constant segments added
    // around it below). The extensions are fictitious: they never steer the
    // slopes of the trace itself, whose ends keep their own tangents.
    const double s_first = q.front().s, s_last = q.back().s;
    if (s_first > 1e-9) {
      PolySegment seg;
      seg.s0 = 0;
      seg.s1 = s_first;
      seg.cy = {q.front().ey, 0, 0, 0};
      seg.cz = {q.front().ez, 0, 0, 0};
      t.profile.push_back(seg);
    }
    std::size_t n = q.size();
    // ALIGNE on a point (p.97): the segment from THIS point to the NEXT one is
    // straight, and the curved neighbours run tangent to it ("X 0 Z 0 ALIGNE /
    // X 1 Z -0.5" is "practically equivalent" to "X 1 Z -0.5 PENTE -0.5").
    // An ALIGNE on a point that already ENDS an alignment is void: the next
    // run is a curve tangent to both straights (ST1: ALIGNE on
    // the three first points of 0 / 12.5 / 27.5 / 40 leaves 12.5-27.5 a
    // Hermite cubic with the slopes of its neighbours, Mz = 1.39 P at
    // midspan where the polyline would give P).
    std::vector<char> aligned(n, 0);
    for (std::size_t k = 0; k + 1 < n; ++k)
      aligned[k] = q[k].aligne && !(k > 0 && aligned[k - 1]);
    auto straight = [&](std::size_t k) { return k + 1 < n && aligned[k]; };
    // TRACE by POLES (p.98, any RAYON given): the broken line through the
    // points, each pole with a RAYON filleted by the circular arc of that
    // radius tangent to its two legs (the others stay sharp) — built exactly
    // below, the arcs sampled as Hermite pieces exact in position and slope.
    bool poles = false;
    for (const Pt& x : q) poles = poles || (x.rayon && *x.rayon > 0);
    // Only the points WITH a RAYON are poles: the segments
    // touching a pole are straight legs, the other points stay points of the
    // trace's spline, tangent to the legs next to them (e.g. a cubic
    // with a free end from (0, 0.3) to (13.75, -1) tangent to the flat leg,
    // y(0.99) = 0.1598, ST1 0.15960). The fillet arcs are put on below.
    auto is_pole = [&](std::size_t k) { return q[k].rayon && *q[k].rayon > 0; };
    if (poles)
      for (std::size_t k = 0; k + 1 < n; ++k)
        if (is_pole(k) || is_pole(k + 1)) aligned[k] = 1;
    // Node slopes: prescribed by PENTE (the angle with the horizontal: the
    // in-plane eccentricity in PLANE, the vertical one in SPATIALE) / GIS
    // (the plan angle), or by a straight neighbour; the others make the trace
    // a C2 cubic spline, with a free (zero-curvature) end where the end slope
    // is not given — in ST1 a free end leaves at -0.2 between
    // (0, 0) and (6, -0.8 PENTE 0), not on the chord.
    std::vector<std::pair<double, double>> pole_lines;  // straight legs (s range)
    {
      // ST1's trace: per segment a PARAMETRIC cubic, the three
      // coordinates (s, e_y, e_z) Hermite cubics of sigma, the cumulative
      // chord of the trace points; node derivatives dP/dsigma:
      //  - a straight neighbour (ALIGNE): the unit chord vector;
      //  - PENTE p and GIS g (angles, rad, p.103): the unit vector (cos p cos g,
      //    cos p sin g, sin p); GIS alone: (cos g, sin g) in plan, unit;
      //  - PENTE alone: lambda (cos p, ., sin p), lambda the x-derivative
      //    of the C2 spline x(sigma) solved first with that node free;
      //  - the free components: C2 at an interior node, zero curvature at an
      //    end (0.13 / 0.6 / 0.27 mm from ST1 on three cables, where the cubic
      //    in s with slope = PENTE left 5.5, 38 and 14 mm; <= 0.9 mm with end
      //    slopes, the print floor of ST1's broken line).
      // Each parametric cubic is then sampled into Hermite pieces in s, exact
      // in position and slope.
      const bool plane = rt_->option_ == 0;
      const int iv = plane ? 1 : 2, ih = plane ? 2 : 1;  // vertical / transverse
      auto P = [&](std::size_t k) { return Vec3(q[k].s, q[k].ey, q[k].ez); };
      std::vector<double> sig(n, 0.0), h(n > 1 ? n - 1 : 0);
      for (std::size_t k = 0; k + 1 < n; ++k) {
        h[k] = (P(k + 1) - P(k)).norm();
        sig[k + 1] = sig[k] + h[k];
      }
      std::vector<Vec3> M(n, Vec3::Zero());
      std::vector<std::array<char, 3>> fixed(n, {0, 0, 0});
      // C2 at a free interior node, zero curvature at a free end: one
      // equation per free derivative of coordinate c3, solved densely (a
      // trace has a handful of points).
      auto solve = [&](int c3) {
        const std::size_t c = static_cast<std::size_t>(c3);
        std::vector<std::size_t> idx(n, n);
        std::size_t nf = 0;
        for (std::size_t k = 0; k < n; ++k)
          if (!fixed[k][c]) idx[k] = nf++;
        if (nf == 0) return;
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(nf), static_cast<Eigen::Index>(nf));
        Eigen::VectorXd r = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(nf));
        auto add = [&](std::size_t row, std::size_t k, double coef) {
          if (fixed[k][c]) r(static_cast<Eigen::Index>(row)) -= coef * M[k][c3];
          else A(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(idx[k])) += coef;
        };
        for (std::size_t k = 0; k < n; ++k) {
          if (fixed[k][c]) continue;
          const std::size_t row = idx[k];
          // P'' just right of k (segment k) and just left of k (segment k-1)
          if (k + 1 < n) {  // +[6 d/h^2 - (4 m_k + 2 m_k+1)/h]
            const double d = P(k + 1)[c3] - P(k)[c3];
            r(static_cast<Eigen::Index>(row)) -= 6 * d / (h[k] * h[k]);
            add(row, k, -4 / h[k]);
            add(row, k + 1, -2 / h[k]);
          }
          if (k > 0) {  // -[-6 d/h^2 + (2 m_k-1 + 4 m_k)/h]
            const double d = P(k)[c3] - P(k - 1)[c3];
            r(static_cast<Eigen::Index>(row)) -= 6 * d / (h[k - 1] * h[k - 1]);
            add(row, k - 1, -2 / h[k - 1]);
            add(row, k, -4 / h[k - 1]);
          }
        }
        const Eigen::VectorXd x = A.fullPivLu().solve(r);
        for (std::size_t k = 0; k < n; ++k)
          if (!fixed[k][c]) M[k][c3] = x(static_cast<Eigen::Index>(idx[k]));
      };
      // the angles of a node (PLANE: global, the eccentricity in the bar axes)
      std::vector<double> pv(n, 0.0), gv(n, 0.0);
      std::vector<char> p_only(n, 0);
      for (std::size_t k = 0; k < n; ++k) {
        const bool gis = q[k].gis && !plane;
        if (q[k].pente) {
          pv[k] = *q[k].pente;
          if (plane)
            for (const Seg& gs : segs)
              if (q[k].s >= gs.s0 - 1e-9 && q[k].s <= gs.s0 + gs.L + 1e-9) {
                pv[k] -= std::atan2(gs.R(0, 1), gs.R(0, 0));
                break;
              }
        }
        if (gis) gv[k] = *q[k].gis;
        // A straight neighbour sets the node (between two straights the node
        // derivative is unused: each straight segment keeps its own direction).
        std::size_t j = n;
        if (k > 0 && straight(k - 1)) j = k - 1;
        else if (straight(k)) j = k;
        if (j < n) {
          M[k] = (P(j + 1) - P(j)) / h[j];
          fixed[k] = {1, 1, 1};
        } else if (q[k].pente && gis) {
          M[k][0] = std::cos(pv[k]) * std::cos(gv[k]);
          M[k][ih] = std::cos(pv[k]) * std::sin(gv[k]);
          M[k][iv] = std::sin(pv[k]);
          fixed[k] = {1, 1, 1};
        } else if (gis) {  // GIS alone: the plan direction, unit
          M[k][0] = std::cos(gv[k]);
          M[k][ih] = std::sin(gv[k]);
          fixed[k][0] = fixed[k][static_cast<std::size_t>(ih)] = 1;
        } else if (q[k].pente) {
          p_only[k] = 1;
        }
      }
      // PENTE alone: its magnitude lambda is the x-derivative of the
      // C2 spline x(sigma) solved first with that node free; the node then
      // takes lambda (cos p, ., sin p), the transverse component left free.
      solve(0);
      for (std::size_t k = 0; k < n; ++k)
        if (p_only[k]) {
          const double lam = M[k][0];
          M[k][0] = lam * std::cos(pv[k]);
          M[k][iv] = lam * std::sin(pv[k]);
          fixed[k][0] = fixed[k][static_cast<std::size_t>(iv)] = 1;
        }
      for (int c3 = 0; c3 < 3; ++c3) solve(c3);
      for (std::size_t k = 0; k + 1 < n; ++k) {
        const Vec3 P0 = P(k), P1 = P(k + 1);
        Vec3 M0 = M[k], M1 = M[k + 1];
        if (straight(k)) M0 = M1 = (P1 - P0) / h[k];
        // the parametric cubic and its derivative at u in [0, 1]
        auto at = [&](double u, Vec3& x, Vec3& dx) {
          const double u2 = u * u, u3 = u2 * u;
          x = (2 * u3 - 3 * u2 + 1) * P0 + (u3 - 2 * u2 + u) * h[k] * M0 + (-2 * u3 + 3 * u2) * P1 +
              (u3 - u2) * h[k] * M1;
          dx = ((6 * u2 - 6 * u) * P0 + (3 * u2 - 4 * u + 1) * h[k] * M0 + (-6 * u2 + 6 * u) * P1 +
                (3 * u2 - 2 * u) * h[k] * M1);
        };
        // The parametric cubic itself, read exactly (PolySegment::exact = 1),
        // when its abscissa advances along it; else Hermite pieces in s.
        bool monotone = !straight(k);
        for (int j = 0; j <= 64 && monotone; ++j) {
          Vec3 xu, du;
          at(j / 64.0, xu, du);
          monotone = du[0] > 1e-9 * h[k];
        }
        if (monotone) {
          Vec3 xa, da, xb, db;
          at(0.0, xa, da);
          at(1.0, xb, db);
          const double hs = q[k + 1].s - q[k].s;
          auto hermite = [&](double v0, double v1, double m0, double m1) {
            std::array<double, 4> c{};
            c[0] = v0;
            c[1] = m0;
            c[2] = (3 * (v1 - v0) / hs - 2 * m0 - m1) / hs;
            c[3] = (2 * (v0 - v1) / hs + m0 + m1) / (hs * hs);
            return c;
          };
          PolySegment seg;
          seg.s0 = q[k].s;
          seg.s1 = q[k + 1].s;
          seg.cy = hermite(xa[1], xb[1], da[1] / da[0], db[1] / db[0]);
          seg.cz = hermite(xa[2], xb[2], da[2] / da[0], db[2] / db[0]);
          seg.exact = 1;
          // X(u) = P0 h00 + h M0 h10 + P1 h01 + h M1 h11 in powers of u
          const Vec3 m0 = h[k] * M0, m1 = h[k] * M1;
          seg.p = {P0, m0, -3 * P0 - 2 * m0 + 3 * P1 - m1, 2 * P0 + m0 - 2 * P1 + m1};
          seg.p[0][0] = q[k].s;  // the abscissae exactly those of the points
          t.profile.push_back(seg);
          continue;
        }
        const int pieces = straight(k) ? 1 : std::max(4, static_cast<int>(std::ceil(h[k] / 2.0)));
        Vec3 xa, da;
        at(0.0, xa, da);
        for (int pc = 1; pc <= pieces; ++pc) {
          Vec3 xb, db;
          at(static_cast<double>(pc) / pieces, xb, db);
          const double hs = xb[0] - xa[0];
          if (hs > 1e-12) {
            auto hermite = [&](double v0, double v1, double m0, double m1) {
              std::array<double, 4> c{};
              c[0] = v0;
              c[1] = m0;
              c[2] = (3 * (v1 - v0) / hs - 2 * m0 - m1) / hs;
              c[3] = (2 * (v0 - v1) / hs + m0 + m1) / (hs * hs);
              return c;
            };
            const double ma_y = da[1] / da[0], mb_y = db[1] / db[0];
            const double ma_z = da[2] / da[0], mb_z = db[2] / db[0];
            PolySegment seg;
            seg.s0 = pc == 1 ? q[k].s : xa[0];
            seg.s1 = pc == pieces ? q[k + 1].s : xb[0];
            seg.cy = hermite(xa[1], xb[1], ma_y, mb_y);
            seg.cz = hermite(xa[2], xb[2], ma_z, mb_z);
            t.profile.push_back(seg);
          }
          xa = xb;
          da = db;
        }
      }
    }
    // RAYON poles: each corner between its two straight legs filleted by the
    // circular arc of that radius tangent to both (p.98), read exactly
    // (PolySegment::exact = 2); the legs are trimmed to the tangent points.
    if (poles) {
      auto Pt3 = [&](std::size_t k) { return Vec3(q[k].s, q[k].ey, q[k].ez); };
      auto leg_poly = [](PolySegment& g, const Vec3& from, const Vec3& d) {
        const double sy = d.y() / d.x(), sz = d.z() / d.x();
        g.cy = {from.y() + sy * (g.s0 - from.x()), sy, 0, 0};
        g.cz = {from.z() + sz * (g.s0 - from.x()), sz, 0, 0};
      };
      for (std::size_t k = 1; k + 1 < n; ++k) {
        if (!is_pole(k)) continue;
        const Vec3 d1 = (Pt3(k) - Pt3(k - 1)).normalized(), d2 = (Pt3(k + 1) - Pt3(k)).normalized();
        const double cth = std::clamp(d1.dot(d2), -1.0, 1.0), th = std::acos(cth);
        if (th < 1e-9) continue;
        double T = *q[k].rayon * std::tan(0.5 * th);
        const double lmax = 0.5 * std::min((Pt3(k) - Pt3(k - 1)).norm(), (Pt3(k + 1) - Pt3(k)).norm());
        if (T > lmax) {
          warn("ST1: CABLE " + std::to_string(c.id) + " — RAYON " + fmt(*q[k].rayon) +
               " trop grand au pole " + std::to_string(k + 1) +
               " : raccordement reduit a la demi-longueur des alignements");
          T = lmax;
        }
        require(d1.x() > 1e-9 && d2.x() > 1e-9,
                "ST1: CABLE " + std::to_string(c.id) + " — raccordement RAYON sur un trace qui revient en arriere");
        const double Rr = T / std::tan(0.5 * th);
        const Vec3 A = Pt3(k) - T * d1, B = Pt3(k) + T * d2;
        Vec3 u = d2 - cth * d1;
        u.normalize();
        std::size_t ib = t.profile.size(), ia = t.profile.size();
        for (std::size_t g = 0; g < t.profile.size(); ++g) {
          if (std::fabs(t.profile[g].s1 - q[k].s) < 1e-9) ib = g;
          if (std::fabs(t.profile[g].s0 - q[k].s) < 1e-9) ia = g;
        }
        require(ib < t.profile.size() && ia < t.profile.size(),
                "ST1: CABLE " + std::to_string(c.id) + " — alignements du pole introuvables");
        t.profile[ib].s1 = A.x();
        leg_poly(t.profile[ib], Pt3(k - 1), d1);
        t.profile[ia].s0 = B.x();
        leg_poly(t.profile[ia], Pt3(k), d2);
        PolySegment arc;
        arc.s0 = A.x();
        arc.s1 = B.x();
        const double hs = arc.s1 - arc.s0;
        auto herm = [&](double v0, double v1, double m0, double m1) {
          std::array<double, 4> c4{};
          c4[0] = v0;
          c4[1] = m0;
          c4[2] = (3 * (v1 - v0) / hs - 2 * m0 - m1) / hs;
          c4[3] = (2 * (v0 - v1) / hs + m0 + m1) / (hs * hs);
          return c4;
        };
        arc.cy = herm(A.y(), B.y(), d1.y() / d1.x(), d2.y() / d2.x());
        arc.cz = herm(A.z(), B.z(), d1.z() / d1.x(), d2.z() / d2.x());
        arc.exact = 2;
        arc.C = A + Rr * u;
        arc.U = u;
        arc.D = d1;
        arc.R = Rr;
        arc.th = th;
        t.profile.insert(t.profile.begin() + static_cast<long>(ia), arc);
      }
      // the straight legs (EXTERIEUR: no rubbing there)
      for (const PolySegment& g : t.profile)
        if (g.exact == 0)
          for (std::size_t k = 0; k + 1 < n; ++k)
            if (straight(k) && g.s0 >= q[k].s - 1e-9 && g.s1 <= q[k + 1].s + 1e-9) {
              pole_lines.push_back({g.s0, g.s1});
              break;
            }
    }
    // TRACE EXTERIEUR (p.96): "un cable exterieur n'est pas lie a la
    // structure sur les troncons rectilignes" — on a straight run it does not
    // rub (no wobble term), the tension stays constant there. A run is known
    // straight when ALIGNE says so, or when both its ends carry the same full
    // tangent (PENTE and GIS) along the chord: ST1 keeps the wobble on a run
    // whose ends give PENTE alone (Fontenoy: 20101, PENTE 0 GIS 0, constant N
    // along bar 2; 20112, PENTE 0 only, the friction / draw-in profile).
    if (c.exterieur && poles)  // the straight legs between the fillets
      for (const auto& [l0, l1] : pole_lines)
        if (l1 - l0 > 1e-9) t.free_runs.push_back({l0, l1});
    if (c.exterieur && !poles)
      for (std::size_t k = 0; k + 1 < n; ++k) {
        const double h = q[k + 1].s - q[k].s;
        bool free_run = straight(k);
        if (!free_run && q[k].pente && q[k].gis && q[k + 1].pente && q[k + 1].gis && h > 0) {
          const double cy = (q[k + 1].ey - q[k].ey) / h, cz = (q[k + 1].ez - q[k].ez) / h;
          const double tol = 1e-4;
          free_run = std::fabs(std::tan(*q[k].pente) - cz) < tol &&
                     std::fabs(std::tan(*q[k + 1].pente) - cz) < tol &&
                     std::fabs(std::tan(*q[k].gis) - cy) < tol && std::fabs(std::tan(*q[k + 1].gis) - cy) < tol;
        }
        if (free_run) t.free_runs.push_back({q[k].s, q[k + 1].s});
      }
    if (s_last < Lchain - 1e-9) {
      PolySegment seg;
      seg.s0 = s_last;
      seg.s1 = Lchain;
      seg.cy = {q.back().ey, 0, 0, 0};
      seg.cz = {q.back().ez, 0, 0, 0};
      t.profile.push_back(seg);
    }
  }

  // SIMUL n: the ST1 cable stands for n identical elementary tendons (p.95).
  double n_simul = std::max(1, c.simul);
  t.jack_force = n_simul * p.jack_force();
  t.jack_both_ends = c.anchor >= 2;
  t.jack_from_end = c.anchor == 1;
  t.release_end_first = c.anchor == 3;
  t.losses.Ep = p.ep;
  t.losses.Ap = n_simul * p.ap;
  // NON SIMUL n: jacked one after the other (annex B1.1.3, phasing only)
  if (c.non_simul) t.losses.non_simul = static_cast<int>(n_simul);
  if (p.pretension) {
    // Pretension: no anchor, the force is transmitted by bond from each
    // bonded end (after the LG_GAINE debonded lengths) along the PRETENSION
    // law; ST1 tables show exactly these shapes (Exemple5). The plateau stays
    // `jack_force`, the reference of the fractions and stresses.
    double L = tendon_chain_length(m, t);
    double g0 = c.lg_gaine_or, g1 = L - c.lg_gaine_ex;
    std::vector<double> ds = {0.0};
    if (p.pret_law == 1) ds.push_back(p.lg_scel1);
    if (p.pret_law == 3) ds.insert(ds.end(), {p.lg_scel1, p.lg_scel2});
    // PARABOL: the parabola sampled by 256 chords (force error <= 2^-16 / 2
    // of the plateau, below any printed digit)
    if (p.pret_law == 2)
      for (int k = 1; k <= 256; ++k) ds.push_back(p.lg_scel1 * k / 256.0);
    std::vector<double> ss = {0.0, L};
    for (double d : ds) {
      ss.push_back(g0 + d);
      ss.push_back(g1 - d);
      if (d == 0) {  // the force jumps from 0 at a bonded end without law
        ss.push_back(g0 + 1e-9 * L);
        ss.push_back(g1 - 1e-9 * L);
      }
    }
    std::sort(ss.begin(), ss.end());
    t.force_profile.clear();
    for (double s : ss) {
      if (s < 0 || s > L) continue;
      if (!t.force_profile.empty() && s <= t.force_profile.back().first + 1e-12 * L) continue;
      double phi = std::min(p.transmission(s - g0), p.transmission(g1 - s));
      t.force_profile.push_back({s, phi * t.jack_force});
    }
  }
  t.losses.rho1000 = p.rho1000;
  // Deferred-loss coupling in the phasing creep loop: same conventions as
  // `deferred_fraction` below — FPK from the deck when given, else the usual
  // sigma_0 = 0.75 fpk jacking assumption; R1000 <= 2.5 % is the EC2 class-2
  // (low relaxation) curve.
  t.losses.fpk = p.fpk > 0 ? p.fpk : (p.sigma > 0 ? p.sigma / 0.75 : 0.0);
  t.losses.relaxation_class = p.rho1000 <= 2.5 ? 2 : 1;
  if (with_losses) {
    // BPEL : P = P0 exp(-(f alpha + phi x))  ->  mu = f, k = phi / f.
    // Eurocode : PERTE MU f KPHI kphi with kphi f = phi  ->  mu = f, k = kphi.
    t.losses.mu = p.f;
    t.losses.k = p.eurocode_form ? p.kphi : (p.f > 0 ? p.phi / p.f : 0.0);
    if (!p.eurocode_form && p.f <= 0) t.losses.phi_abs = p.phi;  // F 0 PHI phi
    t.losses.slip = p.recul;
  } else {
    t.losses.mu = 0;
    t.losses.k = 0;
    t.losses.slip = 0;
  }
  return t;
}

double PrecDef::transmission(double d) const {
  if (d <= 0) return 0.0;
  if (pret_law == 0 || lg_scel1 <= 0) return 1.0;
  if (pret_law == 3) {
    if (d < lg_scel1) return pret_alpha * d / lg_scel1;
    if (d < lg_scel2) return pret_alpha + (1 - pret_alpha) * (d - lg_scel1) / (lg_scel2 - lg_scel1);
    return 1.0;
  }
  double u = std::min(1.0, d / lg_scel1);
  return pret_law == 2 ? 1 - (1 - u) * (1 - u) : u;
}

double St1Domains::deferred_loss(const PrecDef& p, const ChargDef::CableAction& a,
                                 double sigma) const {
  double dsig = 0;
  if (a.retrait != 0) dsig += std::fabs(a.retrait) * p.ep;  // RETRAIT er (a strain)
  if (a.uni != 0) dsig += a.uni;                            // UNI dsig (a stress, > 0 = loss)
  if (a.relax) {
    // EC2_2027 outside a phasing: the BPEL rule with mu0 = 0.35 (annex B2.2)
    const double fprg = p.fpk > 0 ? p.fpk : p.sigma / 0.75;
    const double mu0 = p.loi_ec2_2027 ? 0.35
                                      : (p.rho1000 <= 2.5 ? 0.43 : (p.rho1000 < 8.0 ? 0.35 : 0.30));
    dsig += a.relax_cf * 0.06 * p.rho1000 * std::max(0.0, sigma / fprg - mu0) * sigma;
  }
  return dsig;
}

double St1Domains::deferred_fraction(const PrecDef& p, const ChargDef::CableAction& a) const {
  if (p.sigma <= 0) return 1.0;
  double dsig = 0;
  if (a.retrait != 0) dsig += std::fabs(a.retrait) * p.ep;  // RETRAIT er (a strain)
  if (a.uni != 0) dsig += a.uni;                            // UNI dsig (a stress, > 0 = loss)
  // (PREC PERTE INST, a pretension-only loss, travels with the transmission
  // deficit in lower_cable_actions)
  if (a.relax) {
    // Outside a phasing ST1 uses the BPEL91 3.3,23 rule for BPEL, EC2 and
    // EC2_2027 (annex B2.2): dsig = 6/100 rho1000 (sigma_pi/fprg - mu0)
    // sigma_pi, mu0 = 0.43 (rho1000 <= 2.5, TBR), 0.35 (< 8), 0.30 (RN), and
    // 0.35 ("valeur moyenne BPEL") for EC2_2027, whose KP only enters the
    // phased law (e.g. 0.18079 MN); weighted by cf.
    double fprg = p.fpk > 0 ? p.fpk : p.sigma / 0.75;
    double mu0 = p.loi_ec2_2027 ? 0.35
                                : (p.rho1000 <= 2.5 ? 0.43 : (p.rho1000 < 8.0 ? 0.35 : 0.30));
    double r = 0.06 * p.rho1000 * std::max(0.0, p.sigma / fprg - mu0);
    dsig += a.relax_cf * r * p.sigma;
  }
  return std::max(0.0, 1.0 - dsig / p.sigma);
}

// CHARG ... CABLE ... (p.118-120): lower the stored actions into the load list
// of their CHARG, so the core's own EXEC CHARG solves the enriched case.
void St1Domains::lower_cable_actions() {
  bool any = false;
  for (const auto& [cid, cd] : rt_->chargs_)
    if (!cd.cable_actions.empty() && !cable_lowered_.count(cid)) any = true;
  if (!any) return;
  Model& m = model();
  for (auto& [cid, cd] : rt_->chargs_) {
    if (cd.cable_actions.empty() || cable_lowered_.count(cid)) continue;
    for (const ChargDef::CableAction& a : cd.cable_actions) {
      std::vector<Id> ids = a.cables;
      if (ids.empty())
        for (const auto& [k, v] : cables) {
          (void)v;
          ids.push_back(k);
        }
      for (Id cab : ids) {
        auto it = cables.find(cab);
        if (it == cables.end()) {
          warn("ST1: CHARG CABLE — cable " + std::to_string(cab) + " non defini (p.91)");
          continue;
        }
        auto pit = precs.find(it->second.prec);
        if (pit == precs.end()) {
          warn("ST1: CABLE " + std::to_string(cab) + " — PREC " +
               std::to_string(it->second.prec) + " non defini (p.89)");
          continue;
        }
        Tendon t = build_tendon(m, it->second, pit->second, a.perte && a.perte_inst);
        // PERTE RETRAIT / UNI / RELAX (p.120): ST1 takes the deferred losses
        // point by point, the relaxation on the LOCAL stress after the
        // instantaneous losses (e.g. 193 to 206 MPa along a cable
        // whose jacking stress alone gives 215) — the tendon then carries
        // that profile.
        double frac = 1.0;
        if (a.perte && (a.retrait != 0 || a.uni != 0 || a.relax) && pit->second.sigma > 0) {
          const ForceProfile fpt = tendon_force_profile(m, t);
          const ForceProfile fpi =
              tendon_force_profile(m, build_tendon(m, it->second, pit->second, true));
          const double Ap = t.losses.Ap > 0 ? t.losses.Ap : pit->second.ap;  // all its units
          std::vector<std::pair<double, double>> prof;
          for (std::size_t k = 0; k < fpt.s.size(); ++k) {
            // the instantaneous force at this node: the second node of a
            // slope break is its RIGHT limit (the same trace, the same grid)
            const double Pi = fpi.s == fpt.s ? fpi.P[k]
                              : k > 0 && fpt.s[k] <= fpt.s[k - 1] ? fpi.at(fpt.s[k] + 1e-9)
                                                                   : fpi.at(fpt.s[k]);
            const double sig_i = Ap > 0 ? Pi / Ap : pit->second.sigma;
            // a pretension strand loses in proportion to the force its
            // transmission length has built up (none at a bonded end: Exemple5
            // CHARG 104, ST1 N = 0 at x = 0, -6 on the plateau)
            const double w = pit->second.pretension && t.jack_force > 0
                                 ? std::clamp(Pi / t.jack_force, 0.0, 1.0)
                                 : 1.0;
            prof.push_back({fpt.s[k], fpt.P[k] - w * Ap * deferred_loss(pit->second, a, sig_i)});
          }
          t.force_profile = std::move(prof);
        } else if (a.perte) {
          frac = deferred_fraction(pit->second, a);
        }
        tendons[cab] = t;
        double mean_L = tendon_chain_length(m, t) / static_cast<double>(t.bars.size());
        int nsub = std::max(2, static_cast<int>(std::ceil(mean_L / std::max(1e-6, it->second.pas_cable))));
        nsub = std::min(nsub, 64);
        // TENSION + PERTE: the tension after losses. PERTE alone: only the
        // effect of the losses, i.e. that minus the loss-free tension (p.122).
        std::vector<Load> ls = tendon_loads(m, t, frac, nsub);
        const Tendon t0 = build_tendon(m, it->second, pit->second, false);
        if (!a.tension) {
          std::vector<Load> l0 = tendon_loads(m, t0, -1.0, nsub);
          ls.insert(ls.end(), l0.begin(), l0.end());
        }
        // Pretension, `CABLE i (TENSION) PERTE INST`: ST1 adds the transmission
        // deficit P0 - P(s) of the bonded ends with the sign of the prestress
        // itself, the PREC PERTE INST loss Ap pinst taken off uniformly, debonded
        // runs included (Exemple5 CHARG 107: N = 300, 150, 5.19, 0 at 0, 1, 2,
        // 3 m, the complement of CHARG 101; on a pretension case TENSION PERTE INST is
        // the uniform (sigma - pinst) Ap, PERTE INST alone 828 under the
        // sheath and -12 on the plateau). Reproduced as edited.
        const bool deficit = a.perte_inst && pit->second.pretension;
        const double inst_force = t.losses.Ap * pit->second.perte_inst;
        if (deficit) {
          Tendon td = t0;
          for (auto& [sp, P] : td.force_profile) P = td.jack_force - P - inst_force;
          std::vector<Load> ld = tendon_loads(m, td, 1.0, nsub);
          ls.insert(ls.end(), ld.begin(), ld.end());
        }
        cd.loads.insert(cd.loads.end(), ls.begin(), ls.end());
        // Curvilinear tension table (GET TENSION CABLE, p.276).
        ForceProfile fp = tendon_force_profile(m, t);
        CableTension ct;
        ct.sigma0 = pit->second.sigma;
        double s0 = 0;
        std::vector<std::pair<double, std::pair<Vec3, Mat3>>> frames;
        for (Id bid : t.bars) {
          const Bar& b = m.bar(bid);
          frames.push_back({s0, {m.node(b.ni).x, bar_rotation(m.node(b.ni).x, m.node(b.nj).x, b.roll)}});
          s0 += m.bar_chord_length(b);
        }
        for (std::size_t k = 0; k < fp.s.size(); ++k) {
          double sk = fp.s[k];
          // a slope break is a double node of the profile (its two limits):
          // one section in the table, its left limit
          if (k > 0 && sk <= fp.s[k - 1]) continue;
          ct.s.push_back(sk);
          ct.P.push_back(frac * fp.P[k]);
          ct.sigma.push_back(t.losses.Ap > 0 ? frac * fp.P[k] / t.losses.Ap : 0.0);
          // Global position of the tendon point.
          std::size_t bi = 0;
          for (std::size_t j = 0; j < frames.size(); ++j)
            if (sk >= frames[j].first - 1e-12) bi = j;
          const Bar& b = m.bar(t.bars[bi]);
          Vec3 xi = m.node(b.ni).x, xj = m.node(b.nj).x;
          double Lb = (xj - xi).norm();
          double u = Lb > 0 ? std::clamp((sk - frames[bi].first) / Lb, 0.0, 1.0) : 0.0;
          const Mat3& R = frames[bi].second.second;
          Vec3 pos = xi + u * (xj - xi) + tendon_ey(t, sk) * R.row(1).transpose() +
                     tendon_ez(t, sk) * R.row(2).transpose();
          ct.pos.push_back(pos);
        }
        cable_tension[cab] = ct;
        if (t.losses.Ap > 0) {
          // TENSION: the tension after the losses asked; PERTE alone: the
          // stress the losses take (Exemple5 CHARG 104, RETRAIT 1e-4: -2000
          // on the plateau), with the pretension transmission deficit.
          const double Ap = t.losses.Ap;
          // the debonded runs of a pretension strand carry nothing in the
          // table, whatever the bar efforts say
          const double g0 = it->second.lg_gaine_or;
          const double g1 = tendon_chain_length(m, t) - it->second.lg_gaine_ex;
          const bool pret = pit->second.pretension;
          auto sheathed = [pret, g0, g1](double sk) {
            return pret && (sk < g0 - 1e-9 || sk > g1 + 1e-9);
          };
          const ForceProfile fp0 = tendon_force_profile(m, t0);
          const double P0 = t.jack_force;
          const bool tension = a.tension;
          std::function<double(double)> sig = [&, deficit, P0, tension](double sk) {
            if (sheathed(sk)) return 0.0;
            const double p0 = fp0.at(sk);
            return (frac * fp.at(sk) - (tension ? 0.0 : p0) +
                    (deficit ? P0 - p0 - inst_force : 0.0)) /
                   Ap;
          };
          St1Runtime::CableTable tab;
          tab.cable = cab;
          tab.rows = cable_rows(m, t, it->second.pas_cable, sig, pit->second.sigma);
          auto& v = rt_->cable_tables[{"CHARG", cid}];
          std::erase_if(v, [&](const St1Runtime::CableTable& x) { return x.cable == cab; });
          v.push_back(std::move(tab));
        }
      }
    }
    cable_lowered_.insert(cid);
  }
}

// ===========================================================================
// TABLIER (p.126-135)
// ===========================================================================
void St1Domains::tablier_stmt(const Statement& st, const std::string& head) {
  // `KBARRE listekb KTRANS listek` from word `i` (the KBARRE keyword).
  auto kbarre_line = [&](std::size_t i) {
    RepTransBlock b;
    ++i;
    b.bars = rt_->id_list(st, i, rt_->all_bar_ids());
    if (kw(st, i) == "KTRANS") {
      ++i;
      b.k = num_run(st, i);
    }
    tablier.rep_trans.push_back(b);
  };
  if (sub_ == 1) {
    // REP_TRANS body: `KBARRE listekb KTRANS listek`
    if (head == "KBARRE") {
      kbarre_line(0);
      return;
    }
    sub_ = 0;  // any other keyword leaves the REP_TRANS body
  }
  auto zone_line = [&](const Statement& s, std::size_t i) {
    int iz = static_cast<int>(num(s[i]));
    require(iz >= 1, "ST1: TABLIER ZONE_TRANS — numero de zone >= 1 (p.127)");
    if (static_cast<int>(tablier.zones.size()) < iz)
      tablier.zones.resize(static_cast<std::size_t>(iz));
    TransZone& z = tablier.zones[static_cast<std::size_t>(iz) - 1];
    for (std::size_t k = i + 1; k + 1 < s.size(); ++k) {
      if (kw(s, k) == "LARG") z.width = num(s[k + 1]);
      if (kw(s, k) == "VOIE") z.nvoies = static_cast<int>(num(s[k + 1]));
    }
  };
  if (sub_ == 2) {
    // ZONE_TRANS body: `izone LARG larg (VOIE nvoie)` lines.
    if (!st.empty() && is_id_tok(st[0]) && has_kw(st, "LARG")) {
      zone_line(st, 0);
      return;
    }
    sub_ = 0;  // any other keyword leaves the ZONE_TRANS body
  }
  if (head == "CLASSE") {
    std::size_t i = 1;
    if (kw(st, i) == "TRAFIC") {
      tablier.trafic_class = true;
      ++i;
      if (i < st.size() && is_num_tok(st[i]))
        tablier.classe_trafic = static_cast<int>(num(st[i]));
      return;
    }
    if (i < st.size() && is_num_tok(st[i])) tablier.classe = static_cast<int>(num(st[i]));
    return;
  }
  if (head == "ZONE_TRANS") {
    if (st.size() > 1 && is_id_tok(st[1])) zone_line(st, 1);
    else sub_ = 2;  // the zones follow on their own lines
    return;
  }
  if (head == "PAS_LIGNE") {
    tablier.pas_ligne = num(st[1]);
    return;
  }
  if (head == "LDIF") {
    tablier.ldif = num(st[1]);
    return;
  }
  if (head == "BARRE" || head == "BARRES") {
    std::size_t i = 1;
    tablier.bars = rt_->id_list(st, i, rt_->all_bar_ids());
    return;
  }
  if (head == "REP_TRANS") {
    sub_ = 1;
    // `REP_TRANS KBARRE ... KTRANS ...` on one line (arc deck): its first
    // block sits on the keyword line itself.
    if (kw(st, 1) == "KBARRE") kbarre_line(1);
    return;
  }
  warn("ST1: TABLIER — mot '" + head + "' ignore");
}

// ===========================================================================
// SURCH (p.136-163, p.166)
// ===========================================================================
void St1Domains::surch_stmt(const Statement& st, const std::string& head) {
  (void)head;  // SURCH re-reads the whole line from word 0 (head == st.words[0])
  SurchDef& d = surchs[block_id_];
  std::size_t i = 0;

  auto set_conv = [&](const std::string& name) {
    d.load = SurchLoad::Convoy;
    d.convoy = name;
  };
  static const std::vector<std::string> kCatalog = st1_convoy_names();
  auto in_catalog = [&](const std::string& n) {
    return std::find(kCatalog.begin(), kCatalog.end(), n) != kCatalog.end();
  };

  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k.empty() || st[i].s == ",") {
      ++i;
      continue;
    }
    if (k == "EFFORT") {
      d.priv = SurchPriv::Effort;
      int c = effort_comp(kw(st, i + 1));
      if (c >= 0) {
        d.priv_comp = c;
        d.conc_effort = c;
        ++i;
      }
      ++i;
      continue;
    }
    if (k == "DEPLA") {
      d.priv = SurchPriv::Depla;
      int c = dof_comp(kw(st, i + 1));
      if (c >= 0) {
        d.priv_comp = c;
        d.conc_depla = c;
        ++i;
      }
      ++i;
      continue;
    }
    if (k == "REAC") {
      // the concomitance of the reactions (p.141), studied beside the
      // sections' privileged component — not a replacement for it
      int c = dof_comp(kw(st, i + 1));
      if (c >= 0) {
        d.conc_reac = c;
        ++i;
      }
      ++i;
      continue;
    }
    if (k == "CONTR" || k == "CONTRG") {
      // the extremes of a fibre stress / generalized stress, with their
      // concomitant efforts and fibres (p.137): run_surch_envelope's
      // privileged fibre
      ++i;
      std::string what;
      if (i < st.size() && (k == "CONTRG" ? is_num_tok(st[i]) : !is_num_tok(st[i]))) {
        what = k == "CONTRG" ? std::to_string(static_cast<long>(num(st[i]))) : kw(st, i);
        ++i;
      }
      surch_priv_fibre[block_id_] = {k, what};
      continue;
    }
    if (k == "PRESS") {  // the extremes of a soil pressure (p.137, p.274)
      ++i;
      std::string what;
      if (i < st.size() && !is_num_tok(st[i]) && kw(st, i) != "REAC") what = kw(st, i++);
      surch_priv_fibre[block_id_] = {k, what};
      continue;
    }
    if (k == "ZONE") {
      ++i;
      while (i < st.size() && (is_num_tok(st[i]) || st[i].s == ",")) {
        if (st[i].s != ",") d.zones.push_back(static_cast<int>(num(st[i])));
        ++i;
      }
      continue;
    }
    if (k == "POND") { d.pond = num(st[i + 1]); i += 2; continue; }
    if (k == "PAS_CONV") { d.pas_conv = num(st[i + 1]); i += 2; continue; }
    if (k == "SENS") { d.sens = static_cast<int>(num(st[i + 1])); i += 2; continue; }
    if (k == "SEUIL") { d.seuil = num(st[i + 1]); i += 2; continue; }
    if (k == "PONDVS") { d.pondvs = num(st[i + 1]); i += 2; continue; }
    if (k == "LGEXCLU") { d.lgexclu = num(st[i + 1]); i += 2; continue; }
    // ---- load families ----
    if (k == "AL") {
      ++i;
      std::string a = kw(st, i), b = kw(st, i + 1);
      if (a == "FREIN") { d.load = SurchLoad::AlFrein; ++i; }
      else if (a == "PIETON") { d.load = SurchLoad::AlPieton; ++i; }
      else if (a == "PLANCHER") {
        d.load = (b == "FREIN") ? SurchLoad::AlPlancherFrein : SurchLoad::AlPlancher;
        i += (b == "FREIN") ? 2 : 1;
      } else d.load = SurchLoad::Al;
      continue;
    }
    if (k == "TROT") { d.load = SurchLoad::Trot; ++i; continue; }
    if (k == "TROT_EC") { d.load = SurchLoad::TrotEc; ++i; continue; }
    if (k == "TROT_LM1") { d.load = SurchLoad::TrotLm1; ++i; continue; }
    if (k == "LMF3" || k == "LMF3B" || k == "CH_23M3" || k == "CH_34M3") {
      // fatigue vehicles (p.177), earth-moving plant (p.178)
      ConvoyDef g;
      if (!lmf3_convoy(k, g)) chantier_convoy(k, g);
      d.load = SurchLoad::Convoy;
      d.conv = g;
      d.convoy.clear();
      ++i;
      continue;
    }
    if (k == "ALG") {
      ++i;
      Id j = (i < st.size() && is_id_tok(st[i])) ? static_cast<Id>(num(st[i])) : 1;
      auto it = algs.find(j);
      require(it != algs.end(), "ST1: SURCH ALG " + std::to_string(j) + " non defini (p.156)");
      d.load = SurchLoad::Alg;
      d.alg = it->second;
      ++i;
      continue;
    }
    if (k == "CONV" || k == "CONVOI") {
      ++i;
      Id j = (i < st.size() && is_id_tok(st[i])) ? static_cast<Id>(num(st[i])) : 1;
      auto it = convs.find(j);
      require(it != convs.end(), "ST1: SURCH CONV " + std::to_string(j) + " non defini (p.160)");
      d.load = SurchLoad::Convoy;
      d.conv = it->second;
      d.convoy.clear();
      ++i;
      continue;
    }
    if (k == "LM1") {
      // Inside an LM3 block, `LM1 ...` is the ACCOMPANYING model (p.167).
      if (d.load != SurchLoad::Lm3) d.load = SurchLoad::Lm1;
      ++i;
      continue;
    }
    if (k == "LM2") { d.load = SurchLoad::Lm2; ++i; continue; }
    if (k == "LM4") {
      d.load = kw(st, i + 1) == "VAR" ? SurchLoad::Lm4Var : SurchLoad::Lm4Cte;
      i += 2;
      continue;
    }
    if (k == "LM3") {
      d.load = SurchLoad::Lm3;
      ++i;
      std::string vs = kw(st, i);
      if (i < st.size() && is_id_tok(st[i])) {
        // `LM3 nj`: a user vehicle (CONV generalise, p.160 / p.167)
        const Id j = static_cast<Id>(num(st[i]));
        auto it = convs.find(j);
        require(it != convs.end(), "ST1: SURCH LM3 " + std::to_string(j) + " — CONV non defini");
        d.conv = it->second;
        d.convoy.clear();
        ++i;
      } else if (!vs.empty() && !is_num_tok(st[i])) {
        if (in_catalog(vs)) {
          d.convoy = vs;
          ++i;
        } else {
          ConvoyDef g;
          if (vs_generic(vs, g)) {
            d.conv = g;
            d.convoy.clear();
            ++i;
          }
        }
      }
      continue;
    }
    // LM3 transverse options (p.171).
    if (k == "ZONEVS") {
      if (i + 1 < st.size() && is_num_tok(st[i + 1])) d.zonevs = static_cast<int>(num(st[i + 1]));
      i += 2;
      continue;
    }
    if (k == "EXC" || k == "DEBORDVOIE") {
      double v = (i + 1 < st.size() && is_num_tok(st[i + 1])) ? num(st[i + 1]) : 0.0;
      if (k == "EXC") {
        d.vs_exc = v;
        d.has_exc = true;
      } else {
        d.vs_debv = v;
        d.has_debv = true;
      }
      i += 2;
      continue;
    }
    if (k == "YVS") {
      ++i;
      std::vector<double> v = num_run(st, i);
      if (v.size() > 0) d.yvs1 = v[0];
      if (v.size() > 1) d.yvs2 = v[1];
      continue;
    }
    if (k == "CORRELE") {
      ++i;
      std::vector<double> v = num_run(st, i);
      if (v.size() > 0) d.correle_j = static_cast<int>(v[0]);
      if (v.size() > 1) d.correle_k = static_cast<int>(v[1]);
      continue;
    }
    if (k == "ABS" || k == "REL") {
      d.yvs_abs = k == "ABS";
      ++i;
      continue;
    }
    // inside an LM3 block they set the ACCOMPANYING LM1 (p.167: `LM1 <CARA,
    // (FREQ), COMB, NUL>`)
    if (k == "CARA") { if (d.load == SurchLoad::Lm3) d.lm3_lm1 = 2; else d.freq = false; ++i; continue; }
    if (k == "FREQ") { if (d.load == SurchLoad::Lm3) d.lm3_lm1 = 1; else d.freq = true; ++i; continue; }
    if (k == "COMB" && d.load == SurchLoad::Lm3) { d.lm3_lm1 = 3; ++i; continue; }
    if (k == "CENTRE") { d.lm1_libre = false; ++i; continue; }
    if (k == "LIBRE") { d.lm1_libre = true; ++i; continue; }
    if (k == "NUL") { d.lm3_lm1 = 0; ++i; continue; }
    if (k == "BC" && kw(st, i + 1) == "FREIN") { set_conv("BC_FREIN"); i += 2; continue; }
    if (in_catalog(k)) { set_conv(k); ++i; continue; }
    warn("ST1: SURCH — mot '" + k + "' ignore (p.136-167)");
    ++i;
  }
}

// ===========================================================================
// ALG (p.156-159)
// ===========================================================================
void St1Domains::alg_stmt(const Statement& st, const std::string& head) {
  AlgDef& d = algs[block_id_];
  if (head == "DENS") {
    sub_ = 1;
    std::size_t i = 1;
    std::vector<double> v = num_run(st, i);
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) d.dens.push_back({v[k], v[k + 1]});
    return;
  }
  if (head == "DIR") {
    sub_ = 0;
    std::size_t i = 1;
    d.dir_local = kw(st, i) == "LOC";
    if (kw(st, i) == "LOC" || kw(st, i) == "GLO") ++i;
    Vec3 f = Vec3::Zero();
    bool moment = false;
    while (i + 1 < st.size()) {
      std::string k = kw(st, i);
      if (!is_num_tok(st[i + 1])) break;
      int c = dof_comp(k);
      if (k == "MX" || k == "MY" || k == "MZ") moment = true;
      else if (c >= 0) f(c) = num(st[i + 1]);
      i += 2;
    }
    if (moment)
      warn("ST1: ALG DIR — les composantes MOMENT ne sont pas portees par le moteur "
           "de lignes d'influence (charges de FORCE unitaires seulement, cf. surch.h)");
    // the coefficients multiply the density (p.160): the norm is kept, the
    // engine folds it in (SurchCore::dir_scale)
    if (f.norm() > 0) d.dir = f;
    return;
  }
  if (head == "LARG") {
    sub_ = 0;
    for (std::size_t k = 1; k + 1 < st.size(); ++k)
      if (kw(st, k) == "VOIE") d.larg_voie = num(st[k + 1]);
    return;
  }
  if (head == "COEF") {
    sub_ = 0;
    std::size_t i = 1;
    if (kw(st, i) == "TRANS") ++i;
    d.coef_trans = num_run(st, i);
    return;
  }
  if (sub_ == 1 && !st.empty() && is_num_tok(st[0])) {
    std::size_t i = 0;
    std::vector<double> v = num_run(st, i);
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) d.dens.push_back({v[k], v[k + 1]});
    return;
  }
  warn("ST1: ALG — mot '" + head + "' ignore");
}

// ===========================================================================
// CONV (p.160-162)
// ===========================================================================
void St1Domains::conv_stmt(const Statement& st, const std::string& head) {
  ConvoyDef& d = convs[block_id_];
  if (head == "ESSIEU") {
    sub_ = 1;
    if (st.size() > 1) {
      Statement tail(st.begin() + 1, st.end());
      conv_stmt(tail, kw(tail, 0));
    }
    return;
  }
  if (head == "DIR") {
    sub_ = 0;
    std::size_t i = 1;
    bool loc = kw(st, i) == "LOC";
    if (loc || kw(st, i) == "GLO") ++i;
    bool moment = false;
    Vec3 f = Vec3::Zero();
    while (i + 1 < st.size()) {
      std::string k = kw(st, i);
      if (!is_num_tok(st[i + 1])) break;
      if (k == "MX" || k == "MY" || k == "MZ") moment = true;
      else if (int c = dof_comp(k); c >= 0 && c < 3) f(c) = num(st[i + 1]);
      if (k == "FX" && loc) d.longitudinal = true;
      i += 2;
    }
    // the coefficients multiply the load (p.164)
    if (loc && (f.y() != 0 || f.z() != 0))
      warn("ST1: CONV DIR LOC — seule la composante FX (axe ox local) est portee");
    if (f.norm() > 0) d.dir_coef = f;
    if (moment)
      warn("ST1: CONV DIR — les composantes MOMENT ne sont pas portees par le moteur "
           "de lignes d'influence (cf. surch.h : contourner par REP_TRANS)");
    return;
  }
  if (head == "MAX_CAM") { d.max_cam = static_cast<int>(num(st[1])); sub_ = 0; return; }
  if (head == "MAX_FILE") { d.max_file = static_cast<int>(num(st[1])); sub_ = 0; return; }
  if (head == "LONG") { d.longueur = num(st[1]); sub_ = 0; return; }
  if (head == "LARG") { d.largeur = num(st[1]); sub_ = 0; return; }
  if (head == "COEF" || head == "COEFF") {
    sub_ = 0;
    std::size_t i = 1;
    if (kw(st, i) == "TRANS") ++i;
    d.coef_trans = num_run(st, i);
    d.trans = d.coef_trans.empty() ? TransLaw::None : TransLaw::Table;
    return;
  }
  if (sub_ == 1) {
    // `j XL xl POIDS poids (YL listeyl) (IMPACT impactl impactt)`
    ConvoyAxle a;
    std::size_t i = 0;
    if (i < st.size() && is_id_tok(st[i]) && kw(st, i + 1) == "XL") ++i;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "XL") { a.xl = num(st[i + 1]); i += 2; continue; }
      if (k == "POIDS") { a.weight_t = num(st[i + 1]); i += 2; continue; }
      if (k == "YL") {
        ++i;
        a.yl = num_run(st, i);
        continue;
      }
      if (k == "IMPACT") {
        ++i;
        std::vector<double> v = num_run(st, i);
        if (v.size() > 0) a.impact_l = v[0];
        if (v.size() > 1) a.impact_t = v[1];
        continue;
      }
      ++i;
    }
    d.axles.push_back(a);
    return;
  }
  warn("ST1: CONV — mot '" + head + "' ignore");
}

// ===========================================================================
// EUROCODE (p.165)
// ===========================================================================
void St1Domains::eurocode_stmt(const Statement& st, const std::string& head) {
  std::size_t i = 0;
  if (head == "PSI") {
    // `PSI FREQ (TS) (UDL)`, `PSI COMB (TS) (UDL)` (the LM3 accompanying
    // values, p.169)
    i = 1;
    const bool comb = kw(st, i) == "COMB";
    if (kw(st, i) == "FREQ" || comb) ++i;
    while (i + 1 < st.size()) {
      std::string k = kw(st, i);
      if (!is_num_tok(st[i + 1])) break;
      if (k == "TS") (comb ? eurocode.psi_comb_ts : eurocode.psi_freq_ts) = num(st[i + 1]);
      else if (k == "UDL") (comb ? eurocode.psi_comb_udl : eurocode.psi_freq_udl) = num(st[i + 1]);
      else if (k == "LM2" && !comb) eurocode.psi_freq_lm2 = num(st[i + 1]);
      else warn("ST1: EUROCODE PSI — mot '" + k + "' ignore (p.169)");
      i += 2;
    }
    return;
  }
  if (head == "POIDS") {
    i = 1;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "TS") {
        ++i;
        std::vector<double> v = num_run(st, i);
        for (std::size_t j = 0; j < 3 && j < v.size(); ++j) eurocode.Qk[j] = v[j];
        continue;
      }
      if (k == "UDL") {
        ++i;
        std::vector<double> v = num_run(st, i);
        for (std::size_t j = 0; j < 4 && j < v.size(); ++j) eurocode.qk[j] = v[j];
        continue;
      }
      ++i;
    }
    return;
  }
  if (head == "CLASSE" || head == "AJUSTEMENT") {
    // `CLASSE n` (optionally followed, on the same or the next line, by
    // `AJUSTEMENT TS aQ1,aQi UDL aq1,aqi(,aqr) LM2 beta`).
    if (head == "CLASSE") {
      ec_cur_class_ = static_cast<int>(num(st[1]));
      if (!ec_class_.count(ec_cur_class_)) {
        EurocodeParams base = eurocode;
        if (ec_cur_class_ >= 1 && ec_cur_class_ <= 6) {
          EurocodeParams t = EurocodeParams::traffic_class(ec_cur_class_);
          base.aQ1 = t.aQ1;
          base.aQi = t.aQi;
          base.aq1 = t.aq1;
          base.aqi = t.aqi;
          base.aqr = t.aqr;
          base.beta_lm2 = t.beta_lm2;
        }
        ec_class_[ec_cur_class_] = base;
      }
      i = 2;
    } else {
      require(ec_cur_class_ > 0, "ST1: EUROCODE AJUSTEMENT — CLASSE n prealable (p.170)");
      i = 1;
    }
    EurocodeParams& e = ec_class_[ec_cur_class_];
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "AJUSTEMENT") { ++i; continue; }
      if (k == "TS") {
        ++i;
        std::vector<double> v = num_run(st, i);
        if (v.size() > 0) e.aQ1 = v[0];
        if (v.size() > 1) e.aQi = v[1];
        continue;
      }
      if (k == "UDL") {
        ++i;
        std::vector<double> v = num_run(st, i);
        if (v.size() > 0) e.aq1 = v[0];
        if (v.size() > 1) e.aqi = v[1];
        if (v.size() > 2) e.aqr = v[2];
        continue;
      }
      if (k == "LM2") { e.beta_lm2 = num(st[i + 1]); i += 2; continue; }
      ++i;
    }
    return;
  }
  warn("ST1: EUROCODE — mot '" + head + "' ignore (p.165)");
}

// ===========================================================================
// PHASAGE (p.175-192)
// ===========================================================================
Support St1Domains::support_of(Id aid) const {
  Support sp;
  auto it = rt_->appui_.find(aid);
  require(it != rt_->appui_.end(), "ST1: PHASAGE — appui " + std::to_string(aid) + " non defini");
  auto mod = phas_appui_defs_.find(aid);
  const SupportDef& s = mod != phas_appui_defs_.end() ? mod->second : it->second;
  sp.node = s.node;
  std::array<bool, 6> autob{};
  if (rt_->option_ == 0) autob = {false, false, true, true, true, false};
  else if (rt_->option_ == 1) autob = {true, true, false, false, false, true};
  for (int c = 0; c < 6; ++c)
    sp.fixed[static_cast<std::size_t>(c)] =
        autob[static_cast<std::size_t>(c)] || s.fix[static_cast<std::size_t>(c)];
  sp.spring = s.spring;
  sp.rot = s.rot;
  if (s.decol != 0) {
    sp.decol_axis = s.decol_axis;
    sp.decol_sign = s.decol;
  }
  return sp;
}

void St1Domains::flush_phase(bool essential) {
  if (!phas_dirty_) return;
  PhasageDef& def = phasages[phas_id_];
  phas_cur_.time = phas_time_;
  if (phas_cur_.name.empty())
    phas_cur_.name = "PHASE " + std::to_string(def.plan.phases.size() + 1);
  def.plan.phases.push_back(phas_cur_);
  def.appui_ops.push_back({phas_appui_on_, phas_appui_off_});
  def.appui_mods.push_back(phas_appui_mods_);
  def.phase_envs.push_back(essential ? phas_env_on_ : std::set<Id>{});
  phas_cur_ = Phase{};
  phas_appui_on_.clear();
  phas_appui_off_.clear();
  phas_appui_mods_.clear();
  phas_kind_ = 0;
  phas_dirty_ = false;
}

void St1Domains::phase_step(int kind) {
  // (the step is computed, so an ENV of the phasage sees it: an ENV takes
  // the state before a later ACTIVER APPUI of the same ETAT)
  if (phas_date_step_ || (phas_dirty_ && kind < phas_kind_)) flush_phase(true);
  phas_date_step_ = false;
  phas_kind_ = std::max(phas_kind_, kind);
}

void St1Domains::phasage_stmt(const Statement& st, const std::string& head) {
  PhasageDef& def = phasages[phas_id_];
  if (phas_charg_) {
    phas_charg_stmt(st, head);
    return;
  }
  if (head == "FIN") {
    flush_phase(true);
    close_block();
    return;
  }
  if (head == "ET") {
    // `ET` chains a command onto the previous one without an intermediate
    // computation (p.180): keta computes at dates and states only, so the
    // chained command is simply the next action of the same phase.
    Statement rest(st.begin() + 1, st.end());
    phas_et_ = true;
    if (!rest.empty()) phasage_stmt(rest, kw(rest, 0));
    phas_et_ = false;
    return;
  }
  if (head == "CHARG") {
    phase_step(4);
    phas_charg_ = true;
    phas_charg_mode_ = St1Runtime::Mode::Top;
    phas_case_ = LoadCase{};
    phas_case_.id = 9000 + static_cast<Id>(def.cases.size()) + 1;
    phas_case_.name = first_string(st, 1);
    if (phas_case_.name.empty())
      phas_case_.name = "PHASAGE " + std::to_string(def.id) + " CHARG";
    // `CHARG POIDS PROPRE liste` (Fontenoy): the rest of the line is the
    // first instruction of the block.
    if (st.size() > 1 && !st[1].str) {
      Statement rest(st.begin() + 1, st.end());
      phas_charg_stmt(rest, kw(rest, 0));
    }
    return;
  }
  if (head == "DATE") {
    require(st.size() >= 2 && is_num_tok(st[1]), "ST1: PHASAGE DATE attend une date (p.181)");
    double d = num(st[1]);
    require(d >= phas_time_ - 1e-12,
            "ST1: PHASAGE DATE — la chronologie doit etre croissante (p.181)");
    flush_phase(true);
    phas_dirty_ = d > phas_time_ + 1e-12;  // a new date is itself a computed step (creep interval)
    phas_date_step_ = phas_dirty_;
    phas_time_ = d;
    return;
  }
  if (head == "ETAT") {
    Id e = st.size() > 1 && is_id_tok(st[1]) ? static_cast<Id>(num(st[1])) : 1;
    flush_phase(true);
    if (def.plan.phases.empty()) {
      warn("ST1: PHASAGE ETAT " + std::to_string(e) + " — aucune phase avant l'etat (p.190)");
      return;
    }
    def.etats.push_back({def.plan.phases.size() - 1, e});
    def.etat_names.push_back(first_string(st, 2));
    return;
  }
  if (head == "ACTIVER" || head == "DESACTIVER") {
    bool on = head == "ACTIVER";
    std::size_t i = 1;
    std::vector<Id> line_bars;
    bool mode_abs = false;  // MODE ABS: back to the undeformed geometry (p.187)
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "BARRE") {
        ++i;
        phase_step(on ? 1 : 2);
        std::vector<Id> ids = rt_->id_list(st, i, rt_->defined_bar_ids());
        for (Id b : ids) {
          if (on) phas_cur_.activate_bars.push_back(b);
          else phas_cur_.deactivate_bars.push_back(b);
        }
        line_bars.insert(line_bars.end(), ids.begin(), ids.end());
        phas_dirty_ = true;
        continue;
      }
      if (k == "APPUI") {
        ++i;
        phase_step(on ? 1 : 2);
        std::vector<Id> ids = rt_->id_list(st, i, rt_->all_appui_ids());
        for (Id a : ids) {
          if (on) phas_cur_.add_supports.push_back(support_of(a));
          else phas_cur_.remove_support_nodes.push_back(support_of(a).node);
          (on ? phas_appui_on_ : phas_appui_off_).push_back(a);
        }
        phas_dirty_ = true;
        continue;
      }
      if (k == "ENV") {
        // (DES)ACTIVER ENV liste (p.184): suspend / resume the extremes
        ++i;
        std::vector<Id> all;
        for (const EnvDef& e : def.envs) all.push_back(e.id);
        for (Id e : rt_->id_list(st, i, all)) {
          if (on) phas_env_on_.insert(e);
          else phas_env_on_.erase(e);
        }
        continue;
      }
      if (k == "AGE") {
        // Concrete age at activation (p.187), applied to the bars this
        // ACTIVER line has listed so far: it sets the ages of the creep law
        // and of the modulus E_j. Negative: cast that many days LATER (a
        // pretension bed, Exemple5 AGE -2: no shrinkage, no creep before).
        ++i;
        if (i < st.size() && is_num_tok(st[i])) {
          double age = num(st[i]);
          if (on && age != 0)
            for (Id b : line_bars) phas_cur_.activate_age[b] = age;
          ++i;
        }
        continue;
      }
      if (k == "MODE") {
        ++i;
        mode_abs = kw(st, i) == "ABS";
        ++i;
        continue;
      }
      if (k == "ET") { ++i; continue; }
      ++i;
    }
    // MODE REL (the default): the new nodes continue the deformed structure.
    if (on && !mode_abs)
      for (Id b : line_bars) phas_cur_.activate_rel.insert(b);
    return;
  }
  if (head == "VERINER") {
    // VERINER APPUI liste (p.186): the active support is brought back to a
    // zero TOTAL displacement on its blocked components — from the state the
    // previous commands left, a computed step of its own (an ENV sees the
    // support moved before it is jacked).
    if (phas_dirty_) flush_phase(true);
    phase_step(3);
    std::size_t i = 1;
    if (kw(st, i) == "APPUI") ++i;
    std::vector<Id> ids = rt_->id_list(st, i, rt_->all_appui_ids());
    for (Id a : ids) phas_cur_.jack_support_nodes.push_back(support_of(a).node);
    phas_dirty_ = true;
    return;
  }
  if (head == "MODIFIER") {
    // MODIFIER APPUI liste (NOEUD liste2) (DX ... RZ) (p.195): the APPUI takes
    // a new node and / or new blocked components for the rest of the
    // phasage. On another node: the old one releases its reaction, the new
    // one is set on the deformed structure; on the same node: the freed
    // components release theirs, the new ones start from the deformed state.
    phase_step(3);
    std::size_t i = 1;
    if (kw(st, i) == "APPUI") ++i;
    std::vector<Id> ids = rt_->id_list(st, i, rt_->all_appui_ids());
    std::vector<Id> nodes;
    if (kw(st, i) == "NOEUD") {
      ++i;
      nodes = rt_->id_list(st, i, {});
      require(nodes.size() == ids.size(),
              "ST1: PHASAGE MODIFIER APPUI — listes appuis/noeuds de longueurs differentes");
    }
    std::array<bool, 6> fix{};
    bool any_fix = false;
    for (; i < st.size(); ++i) {
      static const char* names[6] = {"DX", "DY", "DZ", "RX", "RY", "RZ"};
      std::string k = kw(st, i);
      bool known = false;
      for (int c = 0; c < 6; ++c)
        if (k == names[c]) {
          fix[static_cast<std::size_t>(c)] = true;
          any_fix = known = true;
        }
      if (!known && k != "ET") warn("ST1: PHASAGE MODIFIER APPUI — mot '" + k + "' ignore");
    }
    for (std::size_t k = 0; k < ids.size(); ++k) {
      Id a = ids[k];
      const Support before = support_of(a);
      auto cur = phas_appui_defs_.find(a);
      SupportDef d = cur != phas_appui_defs_.end() ? cur->second : rt_->appui_.at(a);
      if (!nodes.empty()) d.node = nodes[k];
      if (any_fix) d.fix = fix;
      phas_appui_defs_[a] = d;
      phas_appui_mods_[a] = d;
      const Support after = support_of(a);
      if (after.node != before.node) {
        phas_cur_.add_supports.push_back(after);
        phas_cur_.remove_support_nodes.push_back(before.node);
      } else {
        phas_cur_.modify_supports.push_back(after);
      }
    }
    phas_dirty_ = true;
    return;
  }
  if (head == "SUITE") {
    // SUITE PHASAGE liste (p.182): this phasage starts from the end state of
    // the listed ones — their phases are replayed ahead of its own at EXEC.
    std::size_t i = 1;
    if (kw(st, i) == "PHASAGE") ++i;
    std::vector<Id> all;
    for (const auto& [k, v] : phasages) {
      (void)v;
      if (k != phas_id_) all.push_back(k);
    }
    for (Id p : rt_->id_list(st, i, all)) {
      require(p != phas_id_ && phasages.count(p),
              "ST1: SUITE PHASAGE " + std::to_string(p) + " — phasage non defini (p.182)");
      def.suite.push_back(p);
    }
    return;
  }
  if (head == "ENV") {
    // ENV v ('Titre') (<EFFORT c, DEPLA c, CONTR f, CONTRG c, PRESS c>) (REAC c)
    // (p.183): the extremes over the essential steps (ETAT, DATE) of the
    // phasage while the envelope is active — from its declaration on.
    require(st.size() >= 2 && is_id_tok(st[1]), "ST1: PHASAGE ENV v (p.183)");
    EnvDef e;
    e.id = static_cast<Id>(num(st[1]));
    e.name = first_string(st, 2);
    if (e.name.empty()) e.name = "ENV " + std::to_string(e.id);
    static const std::vector<std::string> kEffort = {"N", "TY", "TZ", "MX", "MY", "MZ"};
    static const std::vector<std::string> kDepla = {"DX", "DY", "DZ", "RX", "RY", "RZ"};
    static const std::vector<std::string> kPress = {"PX", "PY", "PZ", "MX"};
    static const std::vector<std::string> kFiber = {"VY",   "WY",   "VZ",   "WZ",
                                                    "VYVZ", "VYWZ", "WYVZ", "WYWZ"};
    static const std::vector<std::string> kReac = {"FX", "FY", "FZ", "MX", "MY", "MZ"};
    auto index_of = [](const std::vector<std::string>& v, const std::string& x) {
      for (std::size_t k = 0; k < v.size(); ++k)
        if (v[k] == x) return static_cast<int>(k);
      return -1;
    };
    for (std::size_t i = 2; i + 1 < st.size(); ++i) {
      const std::string k = kw(st, i), c = kw(st, i + 1);
      if (k == "EFFORT" || k == "DEPLA" || k == "PRESS" || k == "CONTR") {
        e.has_priv = true;
        e.what = k == "EFFORT" ? EnvWhat::Effort
                 : k == "DEPLA" ? EnvWhat::Depla
                 : k == "PRESS" ? EnvWhat::Press
                                : EnvWhat::Contr;
        if (k == "CONTR") {
          int f = index_of(kFiber, c);
          require(f >= 0, "ST1: PHASAGE ENV CONTR — fibre '" + c + "' inconnue (p.183)");
          e.fiber = static_cast<ContrFiber>(f);
        } else {
          e.comp = index_of(k == "EFFORT" ? kEffort : k == "DEPLA" ? kDepla : kPress, c);
          require(e.comp >= 0, "ST1: PHASAGE ENV " + k + " — composante '" + c + "' inconnue");
        }
        ++i;
      } else if (k == "CONTRG" && is_id_tok(st[i + 1])) {
        e.has_priv = true;
        e.what = EnvWhat::Contrg;
        e.contrg = static_cast<Id>(num(st[i + 1]));
        ++i;
      } else if (k == "REAC") {
        e.has_reac = true;
        e.reac_comp = index_of(kReac, c);
        require(e.reac_comp >= 0, "ST1: PHASAGE ENV REAC — composante '" + c + "' inconnue");
        ++i;
      }
    }
    def.envs.push_back(e);
    phas_env_on_.insert(e.id);
    return;
  }
  if (head == "TENDRE" || head == "DETENDRE" || head == "INJECTER" ||
      head == "TENDRE_BANC" || head == "RELACHER_BANC") {
    // (grouting goes with the support changes, before the loads of the phase)
    phase_step(head == "INJECTER" ? 3 : 5);
    std::size_t i = 1;
    if (kw(st, i) == "CABLE" || kw(st, i) == "CABLES") ++i;
    std::vector<Id> ids = rt_->id_list(st, i, {});
    // TENDRE CABLE liste (NON) (INJECTE) (p.191): injected by default
    bool unbonded = false;
    for (std::size_t k = i; k < st.size(); ++k)
      if (kw(st, k) == "NON" && kw(st, k + 1) == "INJECTE") unbonded = true;
    if (ids.empty())
      for (const auto& [k, v] : cables) {
        (void)v;
        ids.push_back(k);
      }
    if (head == "INJECTER") {
      for (Id c : ids) phas_cur_.inject.push_back(c);
      phas_dirty_ = true;
      return;
    }
    if (head == "TENDRE_BANC" || head == "RELACHER_BANC") {
      // pretension on a bed (p.193-194): tensioned against the bed, then
      // transmitted by bond at release (keta/phasing.h)
      for (Id c : ids) {
        TensionAction ta{c, 1.0};
        ta.mode = head == "TENDRE_BANC" ? TensionAction::Mode::Bed : TensionAction::Mode::Release;
        phas_cur_.tension.push_back(ta);
      }
      phas_dirty_ = true;
      return;
    }
    // DETENDRE: the negative target fraction is the full-release marker that
    // keta/phasing.h turns into the spring-back of the current effective force.
    double frac = head == "DETENDRE" ? -1.0 : 1.0;
    // a line chained by ET joins the tensioning group before it (no
    // computation between them: simultaneous, K = 0 among them)
    bool joins = false;
    if (phas_et_)
      for (const TensionAction& ta : phas_cur_.tension)
        joins = joins || ta.group == phas_tension_group_;
    if (!joins) ++phas_tension_group_;
    for (Id c : ids) {
      TensionAction ta{c, frac};
      ta.unbonded = head == "TENDRE" && unbonded;
      ta.group = phas_tension_group_;
      phas_cur_.tension.push_back(ta);
    }
    phas_dirty_ = true;
    return;
  }
  warn("ST1: PHASAGE — instruction '" + head + "' non reconnue (p.175-192)");
}

// The chapter-7 load syntax, re-parsed inside a PHASAGE (p.186).
void St1Domains::phas_charg_stmt(const Statement& st, const std::string& head) {
  PhasageDef& def = phasages[phas_id_];
  auto close_charg = [&] {
    def.cases.push_back(phas_case_);
    phas_cur_.cases.push_back(phas_case_.id);
    phas_dirty_ = true;
    phas_charg_ = false;
  };
  // The one-line form `CHARG (definition)` (p.190) needs no FIN: a phase
  // command, or `FIN PHASAGE`, closes the loading implicitly (the block form,
  // or Fontenoy's `CHARG POIDS PROPRE ... FIN`, closes on its FIN).
  static const std::set<std::string> kPhaseCmd = {
      "DATE",    "ETAT",     "ACTIVER",     "DESACTIVER",    "VERINER", "MODIFIER",
      "TENDRE",  "DETENDRE", "INJECTER",    "TENDRE_BANC",   "RELACHER_BANC",
      "SUITE",   "ENV",      "CHARG"};
  if (kPhaseCmd.count(head) || (head == "FIN" && kw(st, 1) == "PHASAGE")) {
    close_charg();
    phasage_stmt(st, head);
    return;
  }
  if (head == "FIN") {
    close_charg();
    return;
  }
  if (head == "ET") {
    // `et TENDRE CABLE ...` inside the CHARG block (Fontenoy): a phase-level
    // command chained onto the loading (p.180), not a load line.
    Statement rest(st.begin() + 1, st.end());
    if (rest.empty()) return;
    phas_charg_ = false;
    phas_et_ = true;
    phasage_stmt(rest, kw(rest, 0));
    phas_et_ = false;
    phas_charg_ = true;
    return;
  }
  // every load line: the CHARG grammar of the runtime (one parser for both)
  rt_->charg_line(phas_case_.loads, phas_charg_mode_, st);
}

// ===========================================================================
// Dynamics blocks
// ===========================================================================
void St1Domains::amortissement_stmt(const Statement& st, const std::string& head) {
  DampingDef& d = dampings[block_id_];
  if (head == "AUTOMATIQUE") {
    d.spec.kind = DampingSpec::Kind::Rayleigh;
    for (std::size_t k = 1; k + 1 < st.size(); ++k) {
      if (kw(st, k) == "MASSE") d.spec.alpha = num(st[k + 1]);
      if (kw(st, k) == "RIGIDITE") d.spec.beta = num(st[k + 1]);
    }
    return;
  }
  if (head == "MODE" && kw(st, 1) == "TOUT") {  // every mode alike — the
    // default of the modes a `MODE j` line does not name
    if (d.spec.kind == DampingSpec::Kind::Rayleigh) d.spec.kind = DampingSpec::Kind::Modal;
    d.spec.zeta = 0;
    for (std::size_t k = 2; k + 1 < st.size(); ++k)
      if (kw(st, k) == "KSI") d.spec.zeta = num(st[k + 1]);
    return;
  }
  if (head == "MODE") {
    if (d.spec.kind != DampingSpec::Kind::PerMode) {
      if (d.spec.kind != DampingSpec::Kind::Modal) d.spec.zeta = 0;  // ksij = 0 by default
      d.spec.kind = DampingSpec::Kind::PerMode;
    }
    // `MODE liste KSI ksi` (manual p.201): a list with commas and `a`
    // intervals ("MODE 2 a 3 KSI 3"), every mode named takes the KSI.
    std::size_t i = 1;
    while (i < st.size()) {
      std::vector<std::size_t> js;
      while (i < st.size() && kw(st, i) != "KSI") {
        if (is_id_tok(st[i])) {
          js.push_back(static_cast<std::size_t>(num(st[i])));
        } else if (kw(st, i) == "A" && !js.empty() && i + 1 < st.size() && is_id_tok(st[i + 1])) {
          const std::size_t hi = static_cast<std::size_t>(num(st[i + 1]));
          for (std::size_t j = js.back() + 1; j <= hi; ++j) js.push_back(j);
          ++i;
        }
        ++i;
      }
      double ksi = 0;
      if (kw(st, i) == "KSI" && i + 1 < st.size() && is_num_tok(st[i + 1])) {
        ksi = num(st[i + 1]);
        i += 2;
      } else {
        ++i;
      }
      for (std::size_t j : js) {
        if (j == 0) continue;
        if (d.spec.per_mode.size() < j) d.spec.per_mode.resize(j, std::nan(""));
        d.spec.per_mode[j - 1] = ksi;
      }
    }
    return;
  }
  if (head == "ENERGIE") {
    // EN 1998-2 4.1.3 (p.200): the KSI of the bars of this set, each mode
    // taking the strain-energy weighted mean of its bars (damping_for)
    d.spec.kind = DampingSpec::Kind::Modal;
    d.energie = true;
    for (std::size_t k = 1; k + 1 < st.size(); ++k)
      if (kw(st, k) == "KSI") d.spec.zeta = num(st[k + 1]);
    return;
  }
  if (head == "POND") {
    d.pond = num(st[1]);
    return;
  }
  warn("ST1: AMORTISSEMENT — mot '" + head + "' ignore (p.196)");
}

DampingSpec St1Domains::damping_for() const {
  if (dampings.empty()) {
    DampingSpec d;
    d.kind = DampingSpec::Kind::Modal;
    d.zeta = 0.0;
    return d;
  }
  std::set<int> used;
  for (const auto& [b, c] : rt_->cons_) {
    (void)b;
    if (c.amortissement != 0) used.insert(c.amortissement);
  }
  Id pick = dampings.begin()->first;
  if (used.size() == 1) pick = static_cast<Id>(*used.begin());
  // ENERGIE, or several sets referenced by CONS / MAT (p.200): each mode
  // takes the mean of its bars' KSI weighted by their modal strain energy
  // (EN 1998-2 eq. 4.1: xi = sum xi_b E_b / sum E_b), a bar without a set
  // taking the default one's
  bool energie = used.size() > 1;
  for (int u : used)
    if (auto du = dampings.find(static_cast<Id>(u)); du != dampings.end() && du->second.energie)
      energie = true;
  if (dampings.count(pick) && dampings.at(pick).energie) energie = true;
  if (energie && rt_->results().has_modal) {
    rt_->build_model();
    const Model& m = rt_->model_;
    const ModalResults& mo = rt_->results().modal;
    const Assembled asm_ = assemble_elements(m);
    auto zeta_of = [&](Id set) {
      auto d = dampings.find(set);
      if (d == dampings.end()) d = dampings.find(pick);
      if (d == dampings.end()) return 0.0;
      const DampingSpec& sp = d->second.spec;
      return (sp.kind == DampingSpec::Kind::Rayleigh ? 0.0 : sp.zeta) * d->second.pond;
    };
    DampingSpec s;
    s.kind = DampingSpec::Kind::PerMode;
    s.per_mode.assign(mo.modes.size(), 0.0);
    for (std::size_t j = 0; j < mo.modes.size(); ++j) {
      double num = 0, den = 0;
      for (std::size_t e = 0; e < m.bars.size(); ++e) {
        const Bar& b = m.bars[e];
        const int ni = m.node_index(b.ni), nj = m.node_index(b.nj);
        Vec12 u;
        u.head<6>() = mo.modes[j][static_cast<std::size_t>(ni)];
        u.tail<6>() = mo.modes[j][static_cast<std::size_t>(nj)];
        const Vec12 ul = asm_.elements[e].local_displacements(u, Vec12::Zero());
        const double eb = 0.5 * ul.dot(asm_.elements[e].end_forces_local(ul, Vec12::Zero()));
        int set = 0;
        if (auto c = rt_->cons_.find(b.id); c != rt_->cons_.end()) set = c->second.amortissement;
        num += zeta_of(set != 0 ? static_cast<Id>(set) : pick) * eb;
        den += eb;
      }
      s.per_mode[j] = den > 0 ? num / den : zeta_of(pick);
    }
    return s;
  }
  if (energie)
    rt_->warn("ST1: AMORTISSEMENT ENERGIE / plusieurs jeux — EXEC MODE requis pour la ponderation "
              "energetique : le jeu " + std::to_string(pick) + " est retenu");
  auto it = dampings.find(pick);
  if (it == dampings.end()) return DampingSpec{};
  DampingSpec s = it->second.spec;
  if (s.kind == DampingSpec::Kind::Rayleigh) {
    s.alpha *= it->second.pond;
    s.beta *= it->second.pond;
  } else if (s.kind == DampingSpec::Kind::Modal) {
    s.zeta *= it->second.pond;
  } else {
    // a mode the block does not name takes the MODE TOUT value, else 0
    // (p.196: "Par defaut, ksij=0") — the list covers any mode count
    if (s.per_mode.size() < 10000) s.per_mode.resize(10000, std::nan(""));
    for (double& z : s.per_mode) z = (std::isnan(z) ? s.zeta : z) * it->second.pond;
  }
  return s;
}

void St1Domains::spectre_stmt(const Statement& st, const std::string& head) {
  SpectrumDef& d = spectra[block_id_];
  if (head == "ACCELERATION") {
    sub_ = 1;
    if (st.size() > 1) {
      Statement tail(st.begin() + 1, st.end());
      spectre_stmt(tail, kw(tail, 0));
    }
    return;
  }
  if (head == "POND") { d.pond = num(st[1]); return; }
  // a BETA line is always the angle (EUROCODE REDUIT's beta0 sits on its own
  // HORIZONTAL / VERTICAL line, p.202)
  if (head == "BETA") { d.beta = num(st[1]); return; }
  if (head == "EUROCODE" || ((head == "HORIZONTAL" || head == "VERTICAL") && sub_ == 2)) {
    // ONE `SPECTRE isp` may carry a HORIZONTAL line AND a VERTICAL line
    // (p.199), each with its own ag / S / ST / periods — hence one EC8 branch
    // per global axis.  `EUROCODE (REDUIT)` alone just opens the axis block.
    d.eurocode = true;
    sub_ = 2;
    std::size_t i = head == "EUROCODE" ? 1 : 0;
    if (kw(st, i) == "REDUIT") { d.reduit = true; ++i; }
    std::vector<int> axes;
    if (kw(st, i) == "HORIZONTAL") {
      ++i;
      if (kw(st, i) == "X") { axes = {0}; ++i; }
      else if (kw(st, i) == "Y") { axes = {1}; ++i; }
      else axes = {0, 1};
    } else if (kw(st, i) == "VERTICAL") {
      axes = {2};
      ++i;
    }
    if (axes.empty()) {
      if (i >= st.size()) return;  // the bare `EUROCODE` header
      axes = {0, 1};
    }
    Ec8Branch b;
    b.set = true;
    b.alpha = axes[0] == 2 ? 3.0 : 2.5;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "PERIODE") {
        ++i;
        std::vector<double> v = num_run(st, i);
        if (v.size() > 0) b.tb = v[0];
        if (v.size() > 1) b.tc = v[1];
        if (v.size() > 2) b.td = v[2];
        if (v.size() > 3) b.te = v[3];
        continue;
      }
      if (i + 1 < st.size() && is_num_tok(st[i + 1])) {
        double v = num(st[i + 1]);
        if (k == "AX") b.ag = v;
        else if (k == "S") b.s_soil = v;
        else if (k == "ST") b.st_topo = v;
        else if (k == "ALPHA") b.alpha = v;
        else if (k == "ETA") b.eta = v;
        else if (k == "Q") b.q = v;
        else if (k == "GAMMA") b.gamma = v;
        else if (k == "BETA") b.beta0 = v;
        i += 2;
        continue;
      }
      ++i;
    }
    for (int a : axes) d.ec8[static_cast<std::size_t>(a)] = b;
    return;
  }
  if (sub_ == 1 || head == "PERIODE") {
    // `PERIODE tj AX axj AY ayj (AZ azj)`
    std::size_t i = 0;
    double T = 0;
    bool haveT = false;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "PERIODE") {
        T = num(st[i + 1]);
        haveT = true;
        i += 2;
        continue;
      }
      if (k == "AX" && haveT) { d.ax.push_back({T, num(st[i + 1])}); i += 2; continue; }
      if (k == "AY" && haveT) { d.ay.push_back({T, num(st[i + 1])}); i += 2; continue; }
      if (k == "AZ" && haveT) { d.az.push_back({T, num(st[i + 1])}); i += 2; continue; }
      ++i;
    }
    return;
  }
  warn("ST1: SPECTRE — mot '" + head + "' ignore (p.198)");
}

void St1Domains::spectre_rep_stmt(const Statement& st, const std::string& head) {
  SpectreResponseDef& d = spectre_reps[block_id_];
  if (head == "SPECTRE") { d.spectrum = static_cast<Id>(num(st[1])); return; }
  if (head == "SRSS") { d.rule = ComboRule::Srss; d.linear = false; return; }
  if (head == "CQC") { d.rule = ComboRule::Cqc; d.linear = false; return; }
  if (head == "LIN") { d.linear = true; return; }
  if (head == "POND") {
    for (std::size_t k = 1; k + 1 < st.size(); ++k) {
      std::string c = kw(st, k);
      if (c == "X") d.pond[0] = num(st[k + 1]);
      else if (c == "Y") d.pond[1] = num(st[k + 1]);
      else if (c == "Z") d.pond[2] = num(st[k + 1]);
    }
    return;
  }
  warn("ST1: SPECTRE REPONSE — mot '" + head + "' ignore (p.201)");
}

void St1Domains::accel_stmt(const Statement& st, const std::string& head) {
  AccelDef& d = accels[block_id_];
  if (head == "METHODE") {
    d.methode = static_cast<int>(num(st[1]));
    for (std::size_t k = 2; k + 1 < st.size(); ++k)
      if (kw(st, k) == "ALPHA") d.hht_alpha = num(st[k + 1]);
    return;
  }
  if (head == "PAS_TEMPS") { d.dt = num(st[1]); return; }
  if (head == "POND") { d.pond = num(st[1]); return; }
  if (head == "BETA") { d.beta = num(st[1]); return; }
  if (head == "ACCELERATION") {
    sub_ = 1;
    if (st.size() > 1) {
      Statement tail(st.begin() + 1, st.end());
      accel_stmt(tail, kw(tail, 0));
    }
    return;
  }
  if (head == "EFFORT" || head == "DEPLA" || head == "CONTR" || head == "CONTRG" ||
      head == "PRESS" || head == "REAC") {
    warn("ST1: ACCELEROGRAMME — la composante privilegiee '" + head +
         "' est enregistree mais l'historique keta enveloppe TOUS les champs (p.203)");
    return;
  }
  if (sub_ == 1) {
    std::size_t i = 0;
    std::array<double, 3> a{{0, 0, 0}};
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "AX") { a[0] = num(st[i + 1]); i += 2; continue; }
      if (k == "AY") { a[1] = num(st[i + 1]); i += 2; continue; }
      if (k == "AZ") { a[2] = num(st[i + 1]); i += 2; continue; }
      ++i;
    }
    d.accel.push_back(a);
    return;
  }
  warn("ST1: ACCELEROGRAMME — mot '" + head + "' ignore (p.203)");
}

void St1Domains::dyn_convoi_stmt(const Statement& st, const std::string& head) {
  DynConvoiDef& d = dyn_convois[block_id_];
  if (head == "METHODE") {
    d.methode = static_cast<int>(num(st[1]));
    for (std::size_t k = 2; k + 1 < st.size(); ++k)
      if (kw(st, k) == "ALPHA") d.hht_alpha = num(st[k + 1]);
    return;
  }
  if (head == "CONVOI" || head == "CONV") { d.convoi = static_cast<Id>(num(st[1])); return; }
  if (head == "PAS_TEMPS") { d.dt = num(st[1]); return; }
  if (head == "VITESSE") { d.speed = num(st[1]); return; }
  if (head == "TEMPS") { d.duration = num(st[1]); return; }
  if (head == "POND") { d.pond = num(st[1]); return; }
  if (head == "EFFORT" || head == "DEPLA" || head == "CONTR" || head == "CONTRG" ||
      head == "PRESS" || head == "REAC")
    return;
  warn("ST1: DYNAMIQUE CONVOI — mot '" + head + "' ignore (p.205)");
}

void St1Domains::dyn_charge_stmt(const Statement& st, const std::string& head) {
  DynChargeDef& d = dyn_charges[block_id_];
  if (head == "METHODE") {
    d.methode = static_cast<int>(num(st[1]));
    for (std::size_t k = 2; k + 1 < st.size(); ++k)
      if (kw(st, k) == "ALPHA") d.hht_alpha = num(st[k + 1]);
    return;
  }
  if (head == "PAS_TEMPS") { d.dt = num(st[1]); return; }
  if (head == "TEMPS") { d.duration = num(st[1]); return; }
  if (head == "MODE") { d.mode = static_cast<int>(num(st[1])); return; }
  if (head == "FREQUENCE") { d.f0 = num(st[1]); return; }
  if (head == "PHI") {
    std::size_t i = 1;
    d.phi = num_run(st, i);
    return;
  }
  if (head == "COEF") {
    std::size_t i = 1;
    d.coef = num_run(st, i);
    return;
  }
  if (head == "BARRE") {
    DynChargeDef::BarLine bl;
    std::size_t i = 1;
    bl.bars = rt_->id_list(st, i, rt_->all_bar_ids());
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "GLO") { bl.frame = LoadFrame::Global; ++i; continue; }
      if (k == "LOC") { bl.frame = LoadFrame::Local; ++i; continue; }
      if (i + 1 < st.size() && is_num_tok(st[i + 1])) {
        int c = dof_comp(k);
        if (k == "FX" || k == "FY" || k == "FZ") bl.F(c) = num(st[i + 1]);
        else if (k == "MX" || k == "MY" || k == "MZ") bl.M(c - 3) = num(st[i + 1]);
        i += 2;
        continue;
      }
      ++i;
    }
    d.lines.push_back(bl);
    return;
  }
  if (head == "EFFORT" || head == "DEPLA" || head == "CONTR" || head == "CONTRG" ||
      head == "PRESS" || head == "REAC")
    return;
  warn("ST1: DYNAMIQUE CHARGE — mot '" + head + "' ignore (p.207)");
}

// ===========================================================================
// MASSE PROPRE (p.69-71)
// ===========================================================================
// `MASSE (SUP)` = full symmetric matrix (6 terms PLANE / 21 SPATIALE),
// `MASSE (SUP) DI` = the diagonal (3 / 6), `MASSE NUL` drops a bar's automatic
// rho*A mass, and `UNI / LIN / PARABOL` add a line mass on a bar.  Masses are
// in MASS units (the deck's RO / G_DYN pair is what the core converts).
void St1Domains::masse_stmt(const Statement& st, const std::string& head) {
  if (head == "PROPRE") {
    if (st.size() > 1) {
      Statement tail(st.begin() + 1, st.end());
      masse_stmt(tail, kw(tail, 0));
    }
    return;
  }
  bool plane = rt_->option_ == 0;
  if (head == "NOEUD") {
    std::size_t i = 1;
    std::vector<Id> ids = rt_->id_list(st, i, rt_->all_node_ids());
    bool sup = false, di = false;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "MASSE") { ++i; continue; }
      if (k == "SUP") { sup = true; ++i; continue; }
      if (k == "DI") { di = true; ++i; continue; }
      break;
    }
    std::vector<double> v = num_run(st, i);
    if (v.empty()) return;
    Mat6 mm = Mat6::Zero();
    if (di) {
      // PLANE order dx dy rz, SPATIALE dx dy dz rx ry rz.
      static const int kPlane[3] = {0, 1, 5};
      for (std::size_t c = 0; c < v.size(); ++c) {
        int dof = plane ? (c < 3 ? kPlane[c] : -1) : (c < 6 ? static_cast<int>(c) : -1);
        if (dof >= 0) mm(dof, dof) = v[c];
      }
    } else {
      // Lower triangle, row by row.
      std::size_t nd = plane ? 3 : 6;
      static const int kPlane[3] = {0, 1, 5};
      std::size_t k = 0;
      for (std::size_t r = 0; r < nd; ++r)
        for (std::size_t c = 0; c <= r; ++c, ++k) {
          if (k >= v.size()) break;
          int dr = plane ? kPlane[r] : static_cast<int>(r);
          int dc = plane ? kPlane[c] : static_cast<int>(c);
          mm(dr, dc) = v[k];
          mm(dc, dr) = v[k];
        }
    }
    for (Id n : ids) {
      NodalMass nm;
      nm.node = n;
      nm.m = mm;
      nm.add = sup;
      masse.nodes.push_back(nm);
    }
    return;
  }
  if (head == "BARRE" || head == "BARRES") {
    std::size_t i = 1;
    std::vector<Id> ids = rt_->id_list(st, i, rt_->all_bar_ids());
    BarMassKind kind = BarMassKind::Auto;
    Vec3 ecc_or = Vec3::Zero(), ecc_ex = Vec3::Zero();
    std::vector<double> v;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "MASSE" || k == "SUP" || k == "DI") { ++i; continue; }
      if (k == "NUL") { kind = BarMassKind::Nul; ++i; continue; }
      if (k == "UNI" || k == "LIN" || k == "PARABOL") {
        kind = k == "UNI" ? BarMassKind::Uni
                          : (k == "LIN" ? BarMassKind::Lin : BarMassKind::Parabol);
        ++i;
        v = num_run(st, i);
        continue;
      }
      if (k == "EXC") {
        ++i;
        int slot = -1, axis = 0;
        while (i < st.size()) {
          std::string t = kw(st, i);
          if (t == "OR") { slot = 0; axis = 0; ++i; continue; }
          if (t == "EX") { slot = 1; axis = 0; ++i; continue; }
          if (t == "X" || t == "Y" || t == "Z") {
            axis = t == "X" ? 0 : (t == "Y" ? 1 : 2);
            ++i;
            continue;
          }
          if (is_num_tok(st[i]) && slot >= 0 && axis < 3) {
            (slot == 0 ? ecc_or : ecc_ex)(axis) = num(st[i]);
            ++axis;
            ++i;
            continue;
          }
          break;
        }
        continue;
      }
      ++i;
    }
    if (kind == BarMassKind::Auto) return;  // explicit no-op (p.69)
    for (Id b : ids) {
      BarMass bm;
      bm.bar = b;
      bm.kind = kind;
      if (v.size() > 0) bm.m1 = v[0];
      if (v.size() > 1) bm.m2 = v[1];
      if (v.size() > 2) bm.m3 = v[2];
      bm.ecc_or = ecc_or;
      bm.ecc_ex = ecc_ex;
      masse.bars.push_back(bm);
    }
    return;
  }
  warn("ST1: MASSE PROPRE — mot '" + head + "' ignore (p.69)");
}

// ===========================================================================
// Fire chain (p.213-228)
// ===========================================================================
void St1Domains::exposition_stmt(const Statement& st, const std::string& head) {
  Exposure& e = exposures[block_id_];
  if (head == "CN") { e.kind = Exposure::Kind::Standard; return; }
  if (head == "HC") { e.kind = Exposure::Kind::Hydrocarbon; return; }
  if (head == "HCM") { e.kind = Exposure::Kind::Hcm; return; }
  if (head == "ADIABATIQUE") { e.kind = Exposure::Kind::Adiabatic; return; }
  if (head == "TEMP") {
    e.kind = Exposure::Kind::ConstantGas;
    e.temp = st.size() > 1 && is_num_tok(st[1]) ? num(st[1]) : 20.0;
    return;
  }
  if (head == "LOI_GAZ" || head == "LOI_PAROI") {
    e.kind = head == "LOI_GAZ" ? Exposure::Kind::GasTable : Exposure::Kind::WallTable;
    std::size_t i = 1;
    std::vector<double> v = num_run(st, i);
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) e.table.push_back({v[k], v[k + 1]});
    return;
  }
  if (head == "PROTECTION") { e.protection = num(st[1]); return; }
  if (head == "CONVECTION") {
    if (kw(st, 1) == "EC1") e.user_convection = false;
    else if (st.size() > 1 && is_num_tok(st[1])) {
      e.user_convection = true;
      e.alpha_c_user = num(st[1]);
    }
    return;
  }
  if (is_num_tok(st[0]) && (e.kind == Exposure::Kind::GasTable || e.kind == Exposure::Kind::WallTable)) {
    std::size_t i = 0;
    std::vector<double> v = num_run(st, i);
    for (std::size_t k = 0; k + 1 < v.size(); k += 2) e.table.push_back({v[k], v[k + 1]});
    return;
  }
  warn("ST1: EXPOSITION FEU — mot '" + head + "' ignore (p.216)");
}

void St1Domains::fire_section_stmt(const Statement& st, const std::string& head) {
  FireSectionDef& d = fire_sections[block_id_];
  if (head == "DALLE" || head == "RECTANGLE" || head == "PERIMETRE") {
    d.section.kind = head == "DALLE" ? FireSection::Kind::Slab
                                     : (head == "RECTANGLE" ? FireSection::Kind::Rectangle
                                                            : FireSection::Kind::Perimeter);
    if (head == "PERIMETRE") sub_ = 1;
    for (std::size_t k = 1; k + 1 < st.size(); ++k) {
      if (kw(st, k) == "EPAISSEUR") d.section.thickness = num(st[k + 1]);
      if (kw(st, k) == "LARGEUR") d.section.width = num(st[k + 1]);
    }
    return;
  }
  if (head == "EXPOSITION") {
    std::size_t i = 1;
    if (kw(st, i) == "FEU") ++i;
    while (i + 1 < st.size()) {
      std::string k = kw(st, i);
      if (!is_num_tok(st[i + 1])) break;
      Id e = static_cast<Id>(num(st[i + 1]));
      if (k == "SUP" || k == "VY") d.expo_ids[0] = e;
      else if (k == "INF" || k == "WY") d.expo_ids[1] = e;
      else if (k == "VZ") d.expo_ids[2] = e;
      else if (k == "WZ") d.expo_ids[3] = e;
      else break;
      i += 2;
    }
    return;
  }
  if (head == "MAT" || head == "MATERIAU") { d.mat = static_cast<Id>(num(st[1])); return; }
  if (head == "MAILLE") { d.section.maille = num(st[1]); return; }
  if (head == "EPAISSEUR") { d.section.thickness = num(st[1]); return; }
  if (head == "LARGEUR") { d.section.width = num(st[1]); return; }
  if (sub_ == 1 && (head == "Y" || is_num_tok(st[0]))) {
    // PERIMETRE vertex: `Y y Z z EXPOSITION expo`
    double y = 0, z = 0;
    Id e = 0;
    std::size_t i = 0;
    int axis = 0;
    while (i < st.size()) {
      std::string k = kw(st, i);
      if (k == "Y") { y = num(st[i + 1]); axis = 1; i += 2; continue; }
      if (k == "Z") { z = num(st[i + 1]); axis = 2; i += 2; continue; }
      if (k == "EXPOSITION") { e = static_cast<Id>(num(st[i + 1])); i += 2; continue; }
      if (is_num_tok(st[i])) {
        if (axis == 0) { y = num(st[i]); axis = 1; }
        else if (axis == 1) { z = num(st[i]); axis = 2; }
        ++i;
        continue;
      }
      ++i;
    }
    d.section.polygon.push_back(Vec3(y, z, 0));
    d.section.seg_exposure.push_back(static_cast<int>(e));
    return;
  }
  warn("ST1: SECTION (feu) — mot '" + head + "' ignore (p.218)");
}

void St1Domains::acier_stmt(const Statement& st, const std::string& head) {
  RebarSteel& s = aciers[block_id_];
  if (head == "DIAMETRE") { s.diameter = num(st[1]); return; }
  if (head == "SECTION") { s.area = num(st[1]); return; }
  if (head == "CLASSE") { s.cls = kw(st, 1) == "X" ? SteelClass::X : SteelClass::N; return; }
  if (head == "LAMINE" || head == "LAMINES") {  // both spellings in the manual
    s.forming = kw(st, 1) == "CHAUD" ? SteelForming::HotRolled : SteelForming::ColdWorked;
    return;
  }
  if (head == "E") { s.E = num(st[1]); return; }
  if (head == "FYK") { s.fyk = num(st[1]); return; }
  warn("ST1: ACIER — mot '" + head + "' ignore (p.223)");
}

void St1Domains::armature_stmt(const Statement& st, const std::string& head) {
  ArmatureDef& d = armatures[block_id_];
  if (head != "ACIER") {
    warn("ST1: ARMATURE — mot '" + head + "' ignore (p.223)");
    return;
  }
  RebarLayer l;
  Id steel = 0;
  std::size_t i = 0;
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k == "ACIER") { steel = static_cast<Id>(num(st[i + 1])); i += 2; continue; }
    if (k == "NOMBRE") { l.count = static_cast<int>(num(st[i + 1])); i += 2; continue; }
    if (k == "ENROBAGE") {
      ++i;
      std::string p = kw(st, i);
      l.placement = p == "INF" ? RebarLayer::Placement::CoverInf : RebarLayer::Placement::CoverSup;
      ++i;
      if (i < st.size() && is_num_tok(st[i])) { l.cover = num(st[i]); ++i; }
      continue;
    }
    if (k == "LOCAL") {
      l.placement = RebarLayer::Placement::Local;
      ++i;
      continue;
    }
    if (k == "Y") { l.y = num(st[i + 1]); i += 2; continue; }
    if (k == "Z") { l.z = num(st[i + 1]); i += 2; continue; }
    ++i;
  }
  auto it = aciers.find(steel);
  require(it != aciers.end(), "ST1: ARMATURE — ACIER " + std::to_string(steel) + " non defini");
  l.steel = it->second;
  d.layers.push_back(l);
  d.steels.push_back(steel);
}

void St1Domains::cara_feu_stmt(const Statement& st, const std::string& head) {
  (void)head;
  Id sec = 0, arm = 0;
  std::size_t i = 0;
  while (i + 1 < st.size()) {
    std::string k = kw(st, i);
    if (k == "SECTION") { sec = static_cast<Id>(num(st[i + 1])); i += 2; continue; }
    if (k == "ARMATURE" || k == "ARMA") { arm = static_cast<Id>(num(st[i + 1])); i += 2; continue; }
    ++i;
  }
  std::vector<Id> ids = cara_ids_.empty() ? std::vector<Id>{block_id_} : cara_ids_;
  for (Id b : ids) {
    FireBarDef& d = fire_bars[b];
    d.bar = b;
    if (sec != 0) d.section = sec;
    if (arm != 0) d.armature = arm;
  }
}

void St1Domains::feu_stmt(const Statement& st, const std::string& head) {
  FeuDef& d = feux[block_id_];
  if (head == "DUREE") { d.duree = num(st[1]); return; }
  if (head == "TEMP") { d.temp = num(st[1]); return; }
  if (head == "MAILLE") { d.maille = num(st[1]); return; }
  warn("ST1: FEU — mot '" + head + "' ignore (p.225)");
}

void St1Domains::incendie_stmt(const Statement& st, const std::string& head) {
  IncendieDef& d = incendies[block_id_];
  if (head != "FEU") {
    warn("ST1: INCENDIE — mot '" + head + "' ignore (p.226)");
    return;
  }
  std::size_t i = 1;
  if (i < st.size() && is_num_tok(st[i])) { d.feu = static_cast<Id>(num(st[i])); ++i; }
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k == "LOCAL") { d.local = true; ++i; continue; }
    if (k == "CHARG") { d.charg = static_cast<Id>(num(st[i + 1])); i += 2; continue; }
    ++i;
  }
}

void St1Domains::analyse_feu_stmt(const Statement& st, const std::string& head) {
  if (head == "NEWTON_RAPHSON") {
    for (std::size_t k = 1; k + 1 < st.size(); ++k) {
      if (kw(st, k) == "MAX_ITER") fire_opt.max_iter = static_cast<int>(num(st[k + 1]));
      if (kw(st, k) == "PERIODE") fire_opt.period = static_cast<int>(num(st[k + 1]));
    }
    return;
  }
  if (head == "Y") { fire_opt.centroid_var = kw(st, 1) != "CTE"; return; }
  if (head == "CARA") { fire_opt.keep_hot = kw(st, 1) == "GARDER"; return; }
  if (head == "SPATIALE") return;  // recorded: keta reads OPTION for the model kind
  warn("ST1: ANALYSE FEU — mot '" + head + "' ignore (p.226)");
}

// ===========================================================================
// EXEC
// ===========================================================================
bool St1Domains::exec_stmt(const Statement& st) {
  std::string what = kw(st, 1);
  if (what == "CHARG") {
    model();
    lower_cable_actions();
    return false;  // the core's exec_charg solves the enriched cases
  }
  if (what == "MODE" || what == "FLAMBEMENT") {
    model();  // materials + MASSE PROPRE folded into the snapshot the core reuses
    return false;
  }
  if (what == "SURCH") { exec_surch(st); return true; }
  if (what == "PHASAGE") { exec_phasage(st); return true; }
  if (what == "SPECTRE") { exec_spectre_reponse(st); return true; }
  if (what == "HISTORIQUE") { exec_historique(st); return true; }
  if (what == "INCENDIE") { exec_incendie(st); return true; }
  return false;
}

EurocodeParams St1Domains::eurocode_for_class(int cls) const {
  EurocodeParams e = eurocode;  // loads (POIDS) and PSI factors
  auto it = ec_class_.find(cls);
  EurocodeParams a = it != ec_class_.end() ? it->second : EurocodeParams::traffic_class(cls);
  e.aQ1 = a.aQ1;
  e.aQi = a.aQi;
  e.aq1 = a.aq1;
  e.aqi = a.aqi;
  e.aqr = a.aqr;
  e.beta_lm2 = a.beta_lm2;
  return e;
}

void St1Domains::exec_surch(const Statement& st) {
  require(has_tablier, "ST1: EXEC SURCH — aucun TABLIER defini (p.126)");
  Model& m = model();
  tablier.validate();
  std::size_t i = 2;
  std::vector<Id> all;
  for (const auto& [k, v] : surchs) {
    (void)v;
    all.push_back(k);
  }
  std::vector<Id> list = rt_->id_list(st, i, all);
  if (list.empty()) list = all;

  // Study sections: the ETUDE domain restricted to the deck chain (p.174).
  // The core replays its ETUDE blocks through the etude module and hands the
  // resolved per-bar metre abscissas back (`etude_stations()`).
  SurchOptions opt;
  opt.stations_per_bar = 3;
  if (rt_->option_ == 0) {  // PLANE: the frame lives in x-y, y up (exemple1/2/4)
    opt.down = Vec3(0, -1, 0);
    opt.reac_comp = 1;
  } else {
    // GRILL / SPATIALE: "effet de l'excentrement des voies" (p.128).
    opt.eccentric_lanes = true;
  }
  if (!tablier.rep_trans_uniform())
    warn("ST1: REP_TRANS — profils differents selon les barres du tablier : le profil de la "
         "barre etudiee (de la premiere barre TABLIER hors tablier) vaut pour tout le chemin "
         "(approximation documentee, surch.h)");
  std::set<Id> deck(tablier.bars.begin(), tablier.bars.end());
  for (const auto& [b, xs] : rt_->etude_stations()) {
    if (!deck.count(b)) continue;
    for (double x : xs) opt.sections.push_back({b, x});
  }
  // ST1 edition (p.137, p.141): every EFFORT / DEPLA study section of every
  // bar and every support, the components of the analysis option.
  SurchEnvSpec base;
  for (const auto& [b, xs] : rt_->etude_stations(StudyFamily::Effort))
    for (double x : xs) base.effort_sections.push_back({b, x});
  for (const auto& [b, xs] : rt_->etude_stations(StudyFamily::Depla))
    for (double x : xs) base.depla_sections.push_back({b, x});
  {
    std::set<Id> seen;
    for (const Support& s : m.supports)
      if (seen.insert(s.node).second) base.support_nodes.push_back(s.node);
  }
  const std::vector<int> comps = rt_->option_ == 0   ? std::vector<int>{0, 1, 5}
                                 : rt_->option_ == 1 ? std::vector<int>{2, 3, 4}
                                                     : std::vector<int>{0, 1, 2, 3, 4, 5};
  base.effort_comps = base.depla_comps = base.reac_comps = comps;
  // CONTR (ST1 edition): the fibres of each EFFORT section as combinations of
  // the station components, read off the fibre stress of unit components with
  // the sections in force now. The CONTRG defined on the bar follow, each on
  // its own influence line too (without the constant, p.99): `n_contr` and
  // `sec_contrg` tell them apart.
  std::vector<std::size_t> n_contr;
  std::vector<std::vector<Id>> sec_contrg;
  if (rt_->opt_.st1_edition) {
    for (const SurchSection& sec : base.effort_sections) {
      std::vector<std::array<double, 6>> fs;
      for (int c = 0; c < 6; ++c) {
        Station unit;
        unit.x = sec.x;
        double* comp[6] = {&unit.N, &unit.Vy, &unit.Vz, &unit.T, &unit.My, &unit.Mz};
        *comp[c] = 1.0;
        const std::vector<double> v = rt_->contr_fibres(rt_->cara_, sec.bar, unit);
        if (fs.empty()) fs.assign(v.size(), std::array<double, 6>{});
        for (std::size_t k = 0; k < v.size(); ++k) fs[k][static_cast<std::size_t>(c)] = v[k];
      }
      bool any = false;
      for (const auto& f : fs)
        for (double a : f) any = any || a != 0.0;
      if (!any) fs.clear();
      n_contr.push_back(fs.size());
      std::vector<Id> cids;
      const double L = rt_->bar_length(sec.bar);
      for (const auto& [cid, cd] : rt_->contrgs) {
        if (!cd.defined_for(sec.bar)) continue;
        std::array<double, 6> g{};
        for (int c = 0; c < 6; ++c) {
          Station unit;
          unit.x = sec.x;
          double* comp[6] = {&unit.N, &unit.Vy, &unit.Vz, &unit.T, &unit.My, &unit.Mz};
          *comp[c] = 1.0;
          g[static_cast<std::size_t>(c)] =
              contrg_at(cd, sec.bar, L > 0 ? std::clamp(sec.x / L, 0.0, 1.0) : 0.0, unit, false);
        }
        fs.push_back(g);
        cids.push_back(cid);
      }
      sec_contrg.push_back(cids);
      base.effort_fibres.push_back(fs);
    }
  }
  // PRESS (p.274): the soil pressures of a DEPLA section on elastic ground,
  // p = k x the local displacement (Py = ky Dy, Exemple11), as combinations
  // of the global components: Px Py (PLANE), Pz Mx (GRILL), Px Py Pz Mx.
  for (const SurchSection& sec : base.depla_sections) {
    std::vector<std::array<double, 6>> fs;
    const Bar& b = m.bar(sec.bar);
    // the bed of the section: the bar's own, else its CARA PSE zone (the
    // last zone runs to the bar end)
    SoilProps soil = b.soil;
    if (!soil.any())
      if (auto sz = rt_->soil_zones_.find(sec.bar); sz != rt_->soil_zones_.end()) {
        const auto& zs = sz->second.zones;
        for (std::size_t z = 0; z < zs.size(); ++z)
          if (z + 1 == zs.size() || sec.x <= zs[z].xl_end + 1e-9) {
            soil = {zs[z].kfx, zs[z].kfy, zs[z].kfz, zs[z].kmx};
            break;
          }
      }
    if (soil.any()) {
      const Mat3 R = bar_rotation(m.node(b.ni).x, m.node(b.nj).x, b.roll);
      auto trans = [&](int c, double k) {
        std::array<double, 6> a{};
        for (int j = 0; j < 3; ++j) a[static_cast<std::size_t>(j)] = k * R(c, j);
        return a;
      };
      std::array<double, 6> mx{};
      for (int j = 0; j < 3; ++j) mx[static_cast<std::size_t>(3 + j)] = soil.kmx * R(0, j);
      if (rt_->option_ == 0) fs = {trans(0, soil.kfx), trans(1, soil.kfy)};
      else if (rt_->option_ == 1) fs = {trans(2, soil.kfz), mx};
      else fs = {trans(0, soil.kfx), trans(1, soil.kfy), trans(2, soil.kfz), mx};
    }
    base.depla_fibres.push_back(fs);
  }
  auto env_spec = [&](const SurchDef& d, Id id) {
    SurchEnvSpec sp = base;
    sp.effort_priv = d.conc_effort;
    sp.depla_priv = d.conc_depla;
    sp.reac_priv = d.conc_reac;
    // CONTR <fibre> / CONTRG j privileged: its index among each section's
    // fibres (the CONTR ones first, in contr_fibres' order, then the CONTRG)
    if (auto pf = surch_priv_fibre.find(id); pf != surch_priv_fibre.end()) {
      const auto& [kind, what] = pf->second;
      int contr_k = -1;
      if (kind == "CONTR") {
        static const std::vector<std::string> plane = {"VY", "WY"}, grill = {"VZ", "WZ"},
                                              spat = {"VYVZ", "WYVZ", "VYWZ", "WYWZ"};
        const auto& names = rt_->option_ == 0 ? plane : (rt_->option_ == 1 ? grill : spat);
        auto it = std::find(names.begin(), names.end(), what);
        contr_k = it != names.end() ? static_cast<int>(it - names.begin()) : 0;
        if (it == names.end())
          warn("ST1: SURCH " + std::to_string(id) + " CONTR '" + what +
               "' — fibre inconnue, la premiere est retenue");
      }
      if (kind == "PRESS") {
        static const std::vector<std::string> plane = {"PX", "PY"}, grill = {"PZ", "MX"},
                                              spat = {"PX", "PY", "PZ", "MX"};
        const auto& names = rt_->option_ == 0 ? plane : (rt_->option_ == 1 ? grill : spat);
        auto it = std::find(names.begin(), names.end(), what);
        const int pk = it != names.end() ? static_cast<int>(it - names.begin()) : 0;
        if (it == names.end())
          warn("ST1: SURCH " + std::to_string(id) + " PRESS '" + what +
               "' — composante inconnue, la premiere est retenue");
        sp.depla_priv_fibre.assign(sp.depla_fibres.size(), -1);
        for (std::size_t r = 0; r < sp.depla_fibres.size(); ++r)
          if (pk < static_cast<int>(sp.depla_fibres[r].size())) sp.depla_priv_fibre[r] = pk;
        return sp;
      }
      sp.effort_priv_fibre.assign(sp.effort_fibres.size(), -1);
      for (std::size_t r = 0; r < sp.effort_fibres.size(); ++r) {
        if (kind == "CONTR") {
          if (r < n_contr.size() && contr_k < static_cast<int>(n_contr[r]))
            sp.effort_priv_fibre[r] = contr_k;
        } else {
          const Id j = static_cast<Id>(std::stol(what.empty() ? "0" : what));
          for (std::size_t g = 0; r < sec_contrg.size() && g < sec_contrg[r].size(); ++g)
            if (sec_contrg[r][g] == j)
              sp.effort_priv_fibre[r] = static_cast<int>(n_contr[r] + g);
        }
      }
      if (!rt_->opt_.st1_edition || sp.effort_fibres.empty())
        warn("ST1: SURCH " + std::to_string(id) + " — composante privilegiee " + kind +
             " sans fibres d'etude : l'etude est menee sur les efforts");
    }
    return sp;
  };
  for (Id id : list) {
    auto it = surchs.find(id);
    if (it == surchs.end()) {
      warn("ST1: EXEC SURCH — surcharge " + std::to_string(id) + " non definie");
      continue;
    }
    SurchDef d = it->second;
    // TABLIER CLASSE TRAFIC icltra (default 2, p.170) selects the adjustment.
    d.ec = (d.load == SurchLoad::Lm1 || d.load == SurchLoad::Lm2 || d.load == SurchLoad::Lm3 ||
            d.load == SurchLoad::Lm4Cte || d.load == SurchLoad::Lm4Var ||
            d.load == SurchLoad::TrotEc || d.load == SurchLoad::TrotLm1)
               ? eurocode_for_class(tablier.eurocode_class())
               : eurocode;
    // PLANE: the 3-D default privileged component (My / ry) is identically
    // zero; the plane frame bends about z and deflects along y.
    if (rt_->option_ == 0 && d.priv_comp == 4)
      d.priv_comp = d.priv == SurchPriv::Effort ? 5 : 1;
    d.pas_conv_grid_only = rt_->opt_.pas_conv_grid_only;  // ST1's own scan (p.147)
    d.validate();
    // A SURCH without ZONE: the load per metre of width on the deck axis, a
    // single convoy, no transverse profile ("A(l) sur 1 m sans profil";
    // ALG: the density x 1 m, the axle x 1)
    Tablier tab = tablier;
    if (d.zones.empty()) {
      tab.zones.clear();
      tab.rep_trans.clear();
    }
    rt_->surch_results[id] = run_surch(m, tab, d, opt);
    rt_->stamp_cara("SURCH", id);
    rt_->surch_envelopes[id] = run_surch_envelope(m, tab, d, env_spec(d, id), opt);
    {
      auto& sc = rt_->surch_contr[id];
      auto& sg = rt_->surch_contrg[id];
      sc.clear();
      sg.clear();
      const std::vector<SurchEnvRow>& rows = rt_->surch_envelopes[id].efforts;
      for (std::size_t r = 0; r < rows.size() && r < n_contr.size(); ++r) {
        const SurchEnvRow& row = rows[r];
        const std::size_t nc = n_contr[r];
        if (row.fibre_lo.size() < nc + sec_contrg[r].size()) continue;
        if (nc > 0)
          sc[row.id].push_back({row.x,
                                std::vector<double>(row.fibre_lo.begin(), row.fibre_lo.begin() + nc),
                                std::vector<double>(row.fibre_hi.begin(), row.fibre_hi.begin() + nc)});
        for (std::size_t g = 0; g < sec_contrg[r].size(); ++g)
          sg[sec_contrg[r][g]][row.id].push_back(
              {row.x, {row.fibre_lo[nc + g]}, {row.fibre_hi[nc + g]}});
      }
      if (sc.empty()) rt_->surch_contr.erase(id);
      if (sg.empty()) rt_->surch_contrg.erase(id);
    }
    rt_->surch_bounds_.erase(id);
    // Laid on the model the SURCH ran on, now: a later ENV may come after a
    // DESACTIVER TOUT (a PHASAGE deck), when the current model has no bars.
    rt_->surch_bounds(id);
  }
  rt_->output().line("EXEC SURCH : " + std::to_string(list.size()) + " surcharge(s) calculees");
}

namespace {
// Field-by-field equality of what a chain of phasages is made of: the
// condition for a `SUITE PHASAGE` to continue a kept run instead of
// replaying its prefix (never a tolerance — the continuation must be the
// phasage that ran, to the bit).
bool same_load(const Load& a, const Load& b) {
  return a.kind == b.kind && a.node == b.node && a.bar == b.bar && a.F == b.F && a.M == b.M &&
         a.F2 == b.F2 && a.frame == b.frame && a.a == b.a && a.b == b.b && a.relative == b.relative &&
         a.dt == b.dt && a.gy == b.gy && a.gz == b.gz && a.g == b.g && a.disp == b.disp &&
         a.strain == b.strain && a.kappa_y == b.kappa_y && a.kappa_z == b.kappa_z &&
         a.slab == b.slab && a.sframe == b.sframe && a.grad == b.grad && a.strain2 == b.strain2 &&
         a.kappa_y2 == b.kappa_y2 && a.kappa_z2 == b.kappa_z2 && a.M2 == b.M2 && a.dvy == b.dvy &&
         a.dvz == b.dvz && a.dthx == b.dthx && a.poly == b.poly && a.holes == b.holes && a.proj == b.proj &&
         a.hydro_gamma == b.hydro_gamma && a.hydro_side == b.hydro_side && a.hydro_point == b.hydro_point &&
         a.Fmid == b.Fmid && a.kappa_x == b.kappa_x && a.kappa_x2 == b.kappa_x2 && a.strain_m == b.strain_m &&
         a.kappa_y_m == b.kappa_y_m && a.kappa_z_m == b.kappa_z_m && a.kappa_x_m == b.kappa_x_m;
}
bool same_supports(const std::vector<Support>& a, const std::vector<Support>& b) {
  if (a.size() != b.size()) return false;
  for (std::size_t k = 0; k < a.size(); ++k)
    if (!equal_support(a[k], b[k])) return false;
  return true;
}
bool same_action(const TensionAction& a, const TensionAction& b) {
  return a.tendon == b.tendon && a.fraction == b.fraction && a.mode == b.mode &&
         a.unbonded == b.unbonded && a.group == b.group;
}
bool same_phase(const Phase& a, const Phase& b) {
  if (a.name != b.name || a.time != b.time || a.activate_bars != b.activate_bars ||
      a.deactivate_bars != b.deactivate_bars || !same_supports(a.add_supports, b.add_supports) ||
      a.remove_support_nodes != b.remove_support_nodes || a.cases != b.cases ||
      a.tension.size() != b.tension.size() || a.activate_age != b.activate_age ||
      a.activate_rel != b.activate_rel || !same_supports(a.modify_supports, b.modify_supports) ||
      a.jack_support_nodes != b.jack_support_nodes || a.inject != b.inject)
    return false;
  for (std::size_t k = 0; k < a.tension.size(); ++k)
    if (!same_action(a.tension[k], b.tension[k])) return false;
  return true;
}
bool same_support_def(const SupportDef& a, const SupportDef& b) {
  return a.fix == b.fix && a.spring == b.spring && a.spring_lt == b.spring_lt &&
         a.spring_di == b.spring_di && a.rot == b.rot && a.decol == b.decol &&
         a.decol_axis == b.decol_axis && a.node == b.node;
}
bool same_support_defs(const std::map<Id, SupportDef>& a, const std::map<Id, SupportDef>& b) {
  if (a.size() != b.size()) return false;
  for (auto ia = a.begin(), ib = b.begin(); ia != a.end(); ++ia, ++ib)
    if (ia->first != ib->first || !same_support_def(ia->second, ib->second)) return false;
  return true;
}
}  // namespace

void St1Domains::exec_phasage(const Statement& st) {
  std::size_t i = 2;
  // EXEC PHASAGE PSE (p.192): the CARA PSE bars with a PLASTIQUE law iterate
  // on their soil states (PhasingOptions::pse); the bars that could not take
  // the exact bed (CARA VAR, releases) stay elastic.
  const bool pse = kw(st, i) == "PSE";
  if (pse) ++i;
  std::vector<Id> all;
  for (const auto& [k, v] : phasages) {
    (void)v;
    all.push_back(k);
  }
  std::vector<Id> list = rt_->id_list(st, i, all);
  if (list.empty()) list = all;

  for (Id id : list) {
    auto it = phasages.find(id);
    if (it == phasages.end()) {
      warn("ST1: EXEC PHASAGE — phasage " + std::to_string(id) + " non defini");
      continue;
    }
    PhasageDef& def = it->second;
    require(!def.plan.phases.empty(),
            "ST1: PHASAGE " + std::to_string(id) + " — aucune phase decrite (p.175)");
    // SUITE PHASAGE (p.182): the continued phasages are replayed ahead of this
    // one (recursively), so it starts from their end state. Their CHARG blocks
    // are renumbered past every id already taken.
    PhasePlan plan;
    std::vector<LoadCase> plan_cases;
    std::vector<std::pair<std::vector<Id>, std::vector<Id>>> plan_appui_ops;
    std::set<Id> visiting;
    using AppuiOp = std::pair<std::vector<Id>, std::vector<Id>>;
    using AppuiMods = std::map<Id, SupportDef>;
    using Seq = std::vector<std::tuple<Phase, AppuiOp, AppuiMods>>;
    std::vector<AppuiMods> plan_appui_mods;
    // A phasage continuing SEVERAL ones (`SUITE PHASAGE 1,2`, e.g.
    // two cantilevers each built from date 0) starts from their histories
    // MERGED by date — parallel works on different bars — then runs its own.
    std::function<Seq(Id)> chain = [&](Id pid) {
      const PhasageDef& pd = phasages.at(pid);
      require(visiting.insert(pid).second,
              "ST1: SUITE PHASAGE — reprise circulaire du phasage " + std::to_string(pid));
      std::vector<Seq> parents;
      for (Id parent : pd.suite) parents.push_back(chain(parent));
      visiting.erase(pid);
      Seq out;
      std::vector<std::size_t> at(parents.size(), 0);
      for (;;) {
        std::size_t best = parents.size();
        for (std::size_t k = 0; k < parents.size(); ++k)
          if (at[k] < parents[k].size() &&
              (best == parents.size() ||
               std::get<0>(parents[k][at[k]]).time <
                   std::get<0>(parents[best][at[best]]).time - 1e-12))
            best = k;
        if (best == parents.size()) break;
        out.push_back(parents[best][at[best]++]);
      }
      Id offset = 0;
      for (const LoadCase& lc : plan_cases) offset = std::max(offset, lc.id);
      std::map<Id, Id> remap;
      for (LoadCase lc : pd.cases) {
        remap[lc.id] = lc.id + offset;
        lc.id += offset;
        plan_cases.push_back(std::move(lc));
      }
      for (std::size_t k = 0; k < pd.plan.phases.size(); ++k) {
        Phase ph = pd.plan.phases[k];
        for (Id& c : ph.cases) c = remap.at(c);
        if (!out.empty())
          require(ph.time >= std::get<0>(out.back()).time - 1e-12,
                  "ST1: SUITE PHASAGE — la date " + fmt(ph.time) + " du phasage " +
                      std::to_string(pid) + " precede la fin du phasage repris (p.181)");
        out.emplace_back(std::move(ph), k < pd.appui_ops.size() ? pd.appui_ops[k] : AppuiOp{},
                         k < pd.appui_mods.size() ? pd.appui_mods[k] : AppuiMods{});
      }
      return out;
    };
    const Seq seq = chain(id);
    for (const auto& [ph, op, mods] : seq) {
      plan.phases.push_back(ph);
      plan_appui_ops.push_back(op);
      plan_appui_mods.push_back(mods);
    }
    const std::size_t own_first = plan.phases.size() - def.plan.phases.size();
    // SUITE PHASAGE of ONE phasage whose run is kept: its chain is the prefix
    // of this one, to be continued from its end state rather than replayed
    // (a phasage continuing several merges their histories: replayed).
    std::shared_ptr<PhasageRun> parent;
    if (def.suite.size() == 1)
      if (auto pit = phasage_runs.find(def.suite.front()); pit != phasage_runs.end())
        parent = pit->second;
    const std::size_t prefix = parent ? parent->seq.size() : 0;
    // Base: the bars the phases cast — a phasage starts from an empty
    // structure, nothing is active by default (p.180, p.186); a bar of the
    // deck the phases never activate takes no part in it.
    std::set<Id> cast;
    for (const Phase& p : plan.phases) cast.insert(p.activate_bars.begin(), p.activate_bars.end());
    sync_materials();  // MAT blocks owned by this extension (what model() does first)
    // the run kept for a later SUITE PHASAGE: the model and everything the
    // options point at live as long as the session does
    auto run = std::make_shared<PhasageRun>();
    run->model = rt_->snapshot_exact_bars(cast);
    run->description_gen = rt_->description_gen();
    run->seq = seq;
    Model& m = run->model;
    m.cases.clear();
    m.combos.clear();
    for (LoadCase lc : plan_cases) {
      lc.loads = rt_->lower_loads(lc.loads, &m);
      m.cases.push_back(std::move(lc));
    }

    // Tendons of every cable a TENDRE/DETENDRE references.
    std::set<Id> want;
    for (const Phase& p : plan.phases)
      for (const TensionAction& t : p.tension) want.insert(t.tendon);
    m.tendons.clear();
    for (Id c : want) {
      auto cit = cables.find(c);
      if (cit == cables.end()) {
        warn("ST1: PHASAGE TENDRE CABLE " + std::to_string(c) + " — cable non defini");
        continue;
      }
      auto pit = precs.find(cit->second.prec);
      if (pit == precs.end()) {
        warn("ST1: CABLE " + std::to_string(c) + " — PREC non defini");
        continue;
      }
      m.tendons.push_back(build_tendon(m, cit->second, pit->second, true));
    }

    // Creep: keta/phasing.h integrates the time-dependent law of each MAT
    // declaring FLUAGE (and RETRAIT): the EC2 one through Material::concrete,
    // every other one (BPEL, BHP, CEB, fib, B3 ...) through its own MatLaw
    // (Material::concrete_law), as ST1 does (manual annex B3-B4).
    std::set<Id> bpel_family_mats;  // model materials whose cables relax by BPEL
    // bars that creep by ST1's incremental method: the EC2 family
    // (EC2, CEB, EC2_2027: no tension creep) and BPEL (tension creep); BHP,
    // fib and B3 keep the superposition until a case identifies theirs
    std::set<Id>& st1_creep_bars = run->st1_creep_bars;
    std::set<Id>& st1_tension_creep_bars = run->st1_tension_creep_bars;
    std::unordered_map<Id, Material*> mat_of;  // the model's materials by id
    for (Material& mm : m.materials) mat_of.emplace(mm.id, &mm);
    for (const Bar& b : m.bars) {
      auto cit = rt_->cons_.find(b.id);
      int mid = 0;
      if (auto bm = bar_mat.find(b.id); bm != bar_mat.end()) mid = bm->second;
      else if (cit != rt_->cons_.end()) mid = cit->second.mat;
      if (cit == rt_->cons_.end() || mid == 0) continue;
      auto mit = mats.find(mid);
      // RETRAIT alone (no FLUAGE) still shortens the bar
      if (mit == mats.end() || (!mit->second.has_fluage && !mit->second.has_retrait)) continue;
      const MatDef& md = mit->second;
      const MatLaw& law = md.law;
      using C = MatLaw::Code;
      if (law.code == C::BPEL || law.code == C::BPEL_BL || law.code == C::BHP1 ||
          law.code == C::BHP2)
        bpel_family_mats.insert(b.material);
      {
        const C cf = md.code_fluage ? *md.code_fluage : law.code;
        if (cf == C::EC2 || cf == C::CEB || cf == C::EC2_2027) st1_creep_bars.insert(b.id);
        if (cf == C::BPEL || cf == C::BPEL_BL) {
          st1_creep_bars.insert(b.id);
          st1_tension_creep_bars.insert(b.id);
        }
      }
      // Every law through its MatLaw (manual annex B4): the EC2 family too, so
      // phi is referenced to Ec = 1.05 Ecm (§3.1.4(2), annex B4.6) and the
      // AVANCE factors and RETRAIT ER apply.
      auto cl = std::make_shared<ConcreteLaw>();
      const double scale = md.stress_scale();  // law moduli in Pa -> deck units
      const bool shrinks = md.has_retrait;
      // the creep and the shrinkage of the law each one names (FLUAGE <loi>,
      // RETRAIT <loi>), E's by default
      MatLaw lawf = law, laws = law;
      if (md.code_fluage) lawf.code = *md.code_fluage;
      if (md.code_retrait) laws.code = *md.code_retrait;
      // Drying runs from casting (ts = 0) in ST1 for every law but EC2 and
      // ACI_B3 (whose TEMPS_CURE sets it): the BPEL family by annex B4.3.1,
      // CEB / fib / BHP / EC2_BHP / EC2_2027 by ST1's results (CEB
      // and fib shrinkage to 4 digits at every ETAT, EC2_2027 16 % -> 3 %).
      const double ts = laws.code == C::EC2 || laws.code == C::ACI_B3 ? laws.ts_dry : 0.0;
      const bool creeps = md.has_fluage;
      cl->phi = [lawf, creeps](double t, double t0) { return creeps ? lawf.phi(t, t0) : 0.0; };
      cl->E_ref = [lawf, scale](double t0) { return lawf.E_ref(t0) / scale; };
      if (md.therm && laws.code == C::EC2) {
        // heat cure (annex B4.6): no autogenous part, drying from the cure's end
        ConcreteTime ct;
        ct.fck = law.fc;
        ct.rh = law.rh;
        ct.h0 = law.h_mm / 1000.0;
        ct.cement_s =
            law.cement == CementClass::R ? 0.20 : (law.cement == CementClass::S ? 0.38 : 0.25);
        ct.t_shrink0 = md.therm_end();
        ct.autogenous = false;
        const double krd = law.avance.k_retrait_dessication;
        cl->shrinkage = [ct, krd, shrinks](double t) {
          return shrinks ? krd * ec2_shrinkage(t, ct) : 0.0;
        };
      } else if (md.therm && (laws.code == C::BPEL || laws.code == C::BHP1 ||
                              laws.code == C::BHP2 || laws.code == C::BPEL_BL)) {
        // heat cure of the BPEL family (BPEL 99 annexe 6, manual annex
        // B4.3): the shrinkage law and the hardening factor of creep run on
        // the age shifted by teb = tpa 1.08^(theta_max - 20), tpa the cure
        // time-temperature integral over (theta_max - 20); f(j - j0) is kept
        const double dmax = md.therm_dtheta_max();
        double tpa = 0;
        for (std::size_t k = 0; k < md.therm_dt.size() && k < md.therm_theta.size(); ++k)
          tpa += (md.therm_theta[k] - 20.0) * md.therm_dt[k];
        const double teb = dmax > 0 ? tpa / dmax * std::pow(1.08, dmax) : 0.0;
        cl->shrinkage = [laws, shrinks, ts, teb](double t) {
          return shrinks ? laws.eps_sh(t + teb, ts) : 0.0;
        };
        cl->phi = [lawf, creeps, teb](double t, double t0) {
          return creeps ? lawf.phi(t + teb, t0 + teb) : 0.0;
        };
        cl->age_shift = teb;
      } else {
        cl->shrinkage = [laws, shrinks, ts](double t) { return shrinks ? laws.eps_sh(t, ts) : 0.0; };
      }
      cl->creep_tension = md.creep_traction;
      if (!md.creep_linear)
        cl->nonlinear_fcm = [law, scale](double t0) { return law.fc_mean(t0) / scale; };
      mat_of.at(b.material)->concrete_law = cl;
    }
    // A cable relaxes by the law of its concrete (manual p.93): BPEL 91 with
    // the equivalent-time method for the BPEL family.
    for (Tendon& t : m.tendons)
      if (!t.bars.empty() && bpel_family_mats.count(m.bar(t.bars.front()).material))
        t.losses.relaxation_bpel = true;
    // ... unless its PREC imposes one: `LOI <BPEL, EC2>` (p.93)
    for (Tendon& t : m.tendons) {
      auto cit = cables.find(t.id);
      if (cit == cables.end()) continue;
      auto pit = precs.find(cit->second.prec);
      if (pit == precs.end()) continue;
      if (pit->second.loi == 1) t.losses.relaxation_bpel = true;
      else if (pit->second.loi == 2) t.losses.relaxation_bpel = false;
      else if (pit->second.loi == 3)
        warn("ST1: PREC " + std::to_string(cit->second.prec) +
             " LOI EC2_2027 — en phasage, keta applique la relaxation EC2 (loi en puissance "
             "tenue hors phasage seulement)");
    }
    // Heat-cured pretension (annex B4.6): the relaxation clock gains t_eb and
    // the strand loses Ep alpha_b (theta_max - 20)(1 - lambda) at release.
    for (Tendon& t : m.tendons) {
      auto cit = cables.find(t.id);
      if (cit == cables.end() || t.bars.empty()) continue;
      auto pit = precs.find(cit->second.prec);
      if (pit == precs.end() || !pit->second.pretension) continue;
      const Id b0 = t.bars.front();
      int mid = 0;
      if (auto bm = bar_mat.find(b0); bm != bar_mat.end()) mid = bm->second;
      else if (auto cc = rt_->cons_.find(b0); cc != rt_->cons_.end()) mid = cc->second.mat;
      auto mit = mats.find(mid);
      if (mit == mats.end() || !mit->second.therm) continue;
      const MatDef& md = mit->second;
      t.losses.relax_teb_days = md.therm_teb();
      t.losses.thermal_loss =
          pit->second.ep * md.alpha * md.therm_dtheta_max() * (1.0 - md.therm_lambda);
    }
    m.validate();

    // Rate-type creep: ST1's incremental method telescopes over any
    // subdivision of the dates, which the single AAEM step does not (a deck
    // dated every few days would creep several times over).
    // ST1 creeps every fibre of every study section (annex B3.1): along the
    // bar, at the ETUDE abscissae the results are edited at.
    PhasingOptions& po = run->po;
    // Each increment with the modulus E_j of the concrete's age on its date
    // (annex B4: E_j from fcm(j)). After a heat cure fcm(j) is fcmp up to its
    // end tp, then runs to fcm at 28 days (EN 1992-1-1 §10.3.1.1): Exemple5,
    // released at j = tp, deflects as E(fcmp) = 0.91 Ecm.
    po.modulus_at = [&](Id bar, double age) {
      int mid = 0;
      if (auto bm = bar_mat.find(bar); bm != bar_mat.end()) mid = bm->second;
      else if (auto cc = rt_->cons_.find(bar); cc != rt_->cons_.end()) mid = cc->second.mat;
      auto mit = mats.find(mid);
      const Material& mm = m.material(m.bar(bar).material);
      if (mit == mats.end() || mit->second.fire || mit->second.law.code == MatLaw::Code::Direct)
        return mm.E;
      const MatDef& md = mit->second;
      using C = MatLaw::Code;
      const bool ec2 = md.law.code == C::EC2 || md.law.code == C::EC2_BHP1 ||
                       md.law.code == C::EC2_BHP2;
      // BPEL: fcj stops at fc28 beyond 28 days (a step at 100 days
      // jumps as Ei28, 37 619 MPa, not 39 258)
      const bool bpel28 = md.law.code == C::BPEL || md.law.code == C::BPEL_BL;
      double E = md.law.E(bpel28 ? std::min(age, 28.0) : age);
      if (md.therm && ec2 && md.therm_fcmp > 0) {
        const double tp = md.therm_end(), fcm = md.law.fc / 1e6 + 8.0;
        if (age < 28.0) {
          const double f = age <= tp ? md.therm_fcmp
                                     : md.therm_fcmp + (fcm - md.therm_fcmp) /
                                                           std::log10(28.0 - tp + 1.0) *
                                                           std::log10(age - tp + 1.0);
          E = md.law.E(28.0) * std::pow(f / fcm, 0.3);
        }
      }
      // BPEL family (BPEL 99 annexe 6, manual annex B4.3): fcj = fcp up to
      // the end of the cure tp, then a log interpolation to fc28 at 28 days;
      // Eij = 11000 fcj^(1/3). A pretension bed is released at the end of the cure as fcp.
      const bool bpel = md.law.code == C::BPEL || md.law.code == C::BHP1 ||
                        md.law.code == C::BHP2 || md.law.code == C::BPEL_BL;
      if (md.therm && bpel && md.therm_fcmp > 0 && age < 28.0) {
        const double tp = md.therm_end(), fc28 = md.law.fc / 1e6;
        const double f = age <= tp ? md.therm_fcmp
                                   : md.therm_fcmp + (fc28 - md.therm_fcmp) /
                                                         std::log10(28.0 - tp + 1.0) *
                                                         std::log10(age - tp + 1.0);
        E = md.law.E(28.0) * std::cbrt(f / fc28);
      }
      return E / md.stress_scale();
    };
    // ST1's own incremental creep method (CreepEngine::St1), unless the session
    // asks for the EN 1992-1-1 superposition (St1Options::creep_superposition)
    po.creep_engine = rt_->opt_.creep_superposition ? CreepEngine::Exact : CreepEngine::St1;
    po.exact_losses = rt_->opt_.exact_losses;
    po.st1_creep_bars = &st1_creep_bars;
    po.st1_tension_creep_bars = &st1_tension_creep_bars;
    po.st1_fibre_layout = rt_->option_;  // 0 PLANE, 1 GRILL, 2 SPATIALE
    po.creep_along_bar = true;
    po.local_losses = true;  // ST1 loses each point of a cable its own strain
    po.born_deformed = true;
    po.threads = rt_->opt_.threads;  // station recovery spread over the bars
    if (rt_->etude_pas_hint > 0)
      po.stations_per_bar = std::max(2, static_cast<int>(rt_->etude_pas_hint) + 1);
    run->sections = rt_->etude_stations_for(m);
    if (!run->sections.empty()) po.stations_abs = &run->sections;
    SoilZoneMap& pse_soil = run->pse_soil;  // the bars of `m` carrying the exact bed
    if (pse)
      for (const auto& [bid, bs] : rt_->soil_zones_) {
        if (!bs.plastic || !bs.plastic->nonlinear()) continue;
        if (m.tapered.count(bid) && !rt_->cara_var_.count(bid)) pse_soil[bid] = bs;
        else
          warn("ST1: EXEC PHASAGE PSE — barre " + std::to_string(bid) +
               " sans element exact sur sol elastique : son sol reste elastique");
      }
    if (!pse_soil.empty()) po.pse = &pse_soil;
    run->pse = pse;
    run->cases = m.cases;
    // Continue the kept run when this chain IS its chain plus own phases: the
    // description at the same generation, the same phases, the same CHARG
    // blocks on the bars it knew, the engine accepting the new model (keta/
    // phasing.h PhasingSession::fork — same footprint, same options). Any
    // difference replays the chain from its first date, as before: exact
    // either way, only the time differs.
    std::string why;  // why the kept run is not continued (empty: it is)
    if (parent) {
      if (parent->description_gen != run->description_gen)
        why = "description modifiee depuis son execution";
      else if (parent->pse != pse)
        why = "options de phasage differentes";
      else if (parent->cases.size() > m.cases.size())
        why = "chargements differents";
      for (std::size_t k = 0; why.empty() && k < prefix; ++k) {
        if (!same_phase(std::get<0>(parent->seq[k]), std::get<0>(seq[k])) ||
            std::get<1>(parent->seq[k]) != std::get<1>(seq[k]) ||
            !same_support_defs(std::get<2>(parent->seq[k]), std::get<2>(seq[k])))
          why = "phasage repris redefini";
      }
      if (why.empty()) {
        // its CHARG blocks lowered on its own bars: the same loads, those on
        // the bars this chain adds set aside
        std::unordered_set<Id> known;
        for (const Bar& b : parent->model.bars) known.insert(b.id);
        for (std::size_t k = 0; why.empty() && k < parent->cases.size(); ++k) {
          const LoadCase& a = parent->cases[k];
          const LoadCase& b = m.cases[k];
          if (a.id != b.id) {
            why = "chargements differents";
            break;
          }
          std::size_t q = 0;
          for (const Load& l : b.loads) {
            const bool on_bar = l.bar != 0 && l.kind != Load::Kind::NodalForce &&
                                l.kind != Load::Kind::SupportDisp;
            if (on_bar && !known.count(l.bar)) continue;
            if (q >= a.loads.size() || !same_load(a.loads[q], l)) {
              why = "chargements differents";
              break;
            }
            ++q;
          }
          if (why.empty() && q != a.loads.size()) why = "chargements differents";
        }
      }
      if (why.empty()) why = parent->session->fork_reject(m, po);
    }
    // the phases this run computes itself: all of them, or the own ones
    std::size_t res_offset = 0;
    if (parent && why.empty()) {
      run->session = parent->session->fork(m, po);
      PhasePlan rest;
      rest.phases.assign(plan.phases.begin() + static_cast<std::ptrdiff_t>(prefix), plan.phases.end());
      require(!rest.phases.empty(), "ST1: PHASAGE " + std::to_string(id) + " — aucune phase propre");
      if (std::getenv("KETA_TIMING"))
        std::fprintf(stderr, "PHASAGE %lld: reprise du phasage %lld a la date %g (%zu phases non rejouees)\n",
                     static_cast<long long>(id), static_cast<long long>(def.suite.front()),
                     run->session->time(), prefix);
      run->session->run(rest);
      run->resumed = true;
      res_offset = prefix;
    } else {
      if (parent && std::getenv("KETA_TIMING"))
        std::fprintf(stderr, "PHASAGE %lld: phasage %lld rejoue depuis la date 0 — %s\n",
                     static_cast<long long>(id), static_cast<long long>(def.suite.front()), why.c_str());
      // bars no phase activates are active from the start (none here: the
      // base holds exactly the cast bars) — the convention of solve_phased
      std::set<Id> from_start;
      for (const Bar& b : m.bars)
        if (!cast.count(b.id)) from_start.insert(b.id);
      run->session = std::make_unique<PhasingSession>(m, po, from_start);
      run->session->run(plan);
    }
    PhasedResults pr;
    pr.phases = run->session->take_results();
    pr.kelvin_fits = run->session->kelvin_fits();
    pr.cache_stats = run->session->cache_stats();
    phasage_runs[id] = run;
    require(!pr.phases.empty(), "ST1: EXEC PHASAGE — aucun resultat");
    // the result of chain phase `pi` (null when it belongs to the continued run)
    auto phase_at = [&](std::size_t pi) -> const PhaseResult* {
      return pi >= res_offset && pi - res_offset < pr.phases.size() ? &pr.phases[pi - res_offset] : nullptr;
    };
    // A construction state is the PERMANENT part of an ST1 combination (ENV i
    // COMB ... PHASAGE j / ETAT e, manual p.231-234): taken whole, the
    // variable actions superposed on it — so the stored states combine.
    CaseResult final_ = pr.phases.back().cumulated;
    final_.non_superposable = false;
    final_.id = id;
    final_.name = def.title.empty() ? ("PHASAGE " + std::to_string(id)) : def.title;
    rt_->results().phasage[id] = final_;
    rt_->stamp_cara("PHASAGE", id);
    // the APPUI active in each phase: the top-level ones, then the phases'
    std::vector<std::set<Id>> phase_appuis(plan.phases.size());
    // and their description (MODIFIER APPUI), when one was redefined
    std::vector<std::map<Id, SupportDef>> phase_appui_defs(plan.phases.size());
    {
      const std::vector<Id> top = rt_->all_appui_ids();
      std::set<Id> act(top.begin(), top.end());
      std::map<Id, SupportDef> defs;
      for (std::size_t k = 0; k < plan.phases.size(); ++k) {
        if (k < plan_appui_ops.size()) {
          act.insert(plan_appui_ops[k].first.begin(), plan_appui_ops[k].first.end());
          for (Id a : plan_appui_ops[k].second) act.erase(a);
        }
        if (k < plan_appui_mods.size())
          for (const auto& [a, d] : plan_appui_mods[k]) defs[a] = d;
        phase_appuis[k] = act;
        phase_appui_defs[k] = defs;
      }
    }
    if (!phase_appuis.empty()) {
      rt_->result_appuis_[{"PHASAGE", id}] = phase_appuis.back();
      rt_->result_appui_defs_[{"PHASAGE", id}] = phase_appui_defs.back();
    }
    // Cable tables of the states: the live cables, their friction / draw-in
    // profile scaled by the effective force of the state.
    // the local loss field of a tendon (TendonPhaseForce::local_loss) at s
    auto interp_loss = [](const std::vector<std::pair<double, double>>& f, double s) {
      if (f.empty()) return 0.0;
      if (s <= f.front().first) return f.front().second;
      if (s >= f.back().first) return f.back().second;
      auto it = std::upper_bound(f.begin(), f.end(), s,
                                 [](double v, const auto& p) { return v < p.first; });
      const auto& b = *it;
      const auto& a = *(it - 1);
      const double u = (s - a.first) / (b.first - a.first);
      return a.second + u * (b.second - a.second);
    };
    // the friction / draw-in profile of each cable: the session's (computed
    // once for the chain, inherited by its continuations); the tables of
    // every state scale it
    auto phase_cables = [&](const PhaseResult& ph) {
      std::vector<St1Runtime::CableTable> tabs;
      for (const TendonPhaseForce& tf : ph.tendons) {
        if (tf.jacked_to <= 0) continue;
        const Tendon* ct = &m.tendon(tf.tendon);  // the phase lists the model's tendons
        auto cd = cables.find(tf.tendon);
        if (cd == cables.end() || ct->losses.Ap <= 0 || ct->jack_force <= 0) continue;
        auto pd = precs.find(cd->second.prec);
        const double sigma0 = pd != precs.end() ? pd->second.sigma : 0.0;
        const ForceProfile& fp = run->session->force_profile(ct->id);
        const double ratio = tf.effective / ct->jack_force, Ap = ct->losses.Ap;
        St1Runtime::CableTable tab;
        tab.cable = tf.tendon;
        // pretension: the strand on its bed is at sigma0 over its bonded
        // length (0 where LG_GAINE sheathes it, ST1 Exemple5), the reference
        // of Var.Tension before and after the release
        const bool pret = pd != precs.end() && pd->second.pretension;
        const double Lc = tendon_chain_length(m, *ct);
        auto bonded = [&](double sk) {
          return sk >= cd->second.lg_gaine_or - 1e-9 &&
                 sk <= Lc - cd->second.lg_gaine_ex + 1e-9;
        };
        std::function<double(double)> ref;
        if (pret) ref = [&](double sk) { return bonded(sk) ? sigma0 : 0.0; };
        if (tf.bed)
          tab.rows = cable_rows(m, *ct, cd->second.pas_cable,
                                [&](double sk) { return bonded(sk) ? ratio * sigma0 : 0.0; },
                                sigma0, ref);
        else
          tab.rows = cable_rows(m, *ct, cd->second.pas_cable,
                                [&](double sk) {
                                  const double P = fp.at(sk);
                                  if (P <= 0) return 0.0;  // not bonded / outside the anchors
                                  return (ratio * P - interp_loss(tf.local_loss, sk)) / Ap;
                                },
                                sigma0, ref);
        tabs.push_back(std::move(tab));
      }
      std::sort(tabs.begin(), tabs.end(),
                [](const auto& a, const auto& b) { return a.cable < b.cable; });
      return tabs;
    };
    rt_->cable_tables[{"PHASAGE", id}] = phase_cables(pr.phases.back());
    for (std::size_t k = 0; k < def.etats.size(); ++k) {
      std::size_t pi = own_first + def.etats[k].first;
      Id eid = def.etats[k].second;
      const PhaseResult* pres = phase_at(pi);
      if (!pres) continue;
      CaseResult e = pres->cumulated;
      e.non_superposable = false;
      e.id = eid;
      e.name = def.etat_names[k].empty() ? ("ETAT " + std::to_string(eid)) : def.etat_names[k];
      rt_->results().etat[eid] = e;
      rt_->stamp_cara("ETAT", eid);
      rt_->result_appuis_[{"ETAT", eid}] = phase_appuis[pi];
      rt_->result_appui_defs_[{"ETAT", eid}] = phase_appui_defs[pi];
      rt_->cable_tables[{"ETAT", eid}] = phase_cables(*pres);
    }
    // ENV of the phasage (p.183): the essential steps where each is active
    if (!def.envs.empty()) {
      // the bars active after each phase (never activated in the plan: from
      // the start); an inactive bar is out of the envelope at that step,
      // an inactive APPUI reads 0 there
      std::vector<std::set<Id>> active_after(plan.phases.size());
      {
        std::set<Id> ever, act;
        for (const Phase& ph : plan.phases) ever.insert(ph.activate_bars.begin(), ph.activate_bars.end());
        for (const Bar& b : m.bars)
          if (!ever.count(b.id)) act.insert(b.id);
        for (std::size_t k = 0; k < plan.phases.size(); ++k) {
          act.insert(plan.phases[k].activate_bars.begin(), plan.phases[k].activate_bars.end());
          for (Id b : plan.phases[k].deactivate_bars) act.erase(b);
          active_after[k] = act;
        }
      }
      // every declared APPUI (an inactive one reads 0), DESACTIVER TOUT or not
      const std::vector<Id>& all_appuis = rt_->appui_order_;
      std::map<Id, std::vector<Id>> env_steps;
      for (std::size_t k = 0; k < def.phase_envs.size(); ++k) {
        const std::size_t pi = own_first + k;
        const PhaseResult* pres = def.phase_envs[k].empty() ? nullptr : phase_at(pi);
        if (!pres) continue;
        const Id uid = id * 100000 + static_cast<Id>(k) + 1;
        CaseResult r = pres->cumulated;
        std::erase_if(r.bars, [&](const BarResult& b) { return !active_after[pi].count(b.bar); });
        // a support node no active APPUI holds any more reacts 0 (its
        // accumulated reaction was released when it was removed)
        {
          std::set<Id> held;
          for (Id a : phase_appuis[pi]) {
            auto d = phase_appui_defs[pi].find(a);
            held.insert(d != phase_appui_defs[pi].end() ? d->second.node : rt_->appui_.at(a).node);
          }
          for (Reaction& re : r.reactions)
            if (!held.count(re.node)) re.R.setZero();
        }
        r.non_superposable = false;
        r.id = uid;
        r.name = pres->name;
        rt_->results().phase_step[uid] = r;
        rt_->stamp_cara("PHASE_STEP", uid);
        rt_->result_appuis_[{"PHASE_STEP", uid}] =
            std::set<Id>(all_appuis.begin(), all_appuis.end());
        rt_->result_appui_defs_[{"PHASE_STEP", uid}] = phase_appui_defs[pi];
        for (Id e : def.phase_envs[k]) env_steps[e].push_back(uid);
      }
      for (EnvDef e : def.envs) {
        for (Id uid : env_steps[e.id]) e.sources.push_back({"PHASE_STEP", uid, 1.0, 1.0});
        if (e.sources.empty()) {
          warn("ST1: PHASAGE " + std::to_string(id) + " ENV " + std::to_string(e.id) +
               " — aucune etape de calcul pendant que l'enveloppe est active");
          continue;
        }
        rt_->envs[e.id] = e;
        rt_->compute_env(rt_->envs.at(e.id));
      }
    }
    rt_->output().line("EXEC PHASAGE " + std::to_string(id) + " : " +
                       std::to_string(plan.phases.size()) + " phase(s), " +
                       std::to_string(def.etats.size()) + " etat(s)");
  }
}

void St1Domains::exec_spectre_reponse(const Statement& st) {
  require(kw(st, 2) == "REPONSE", "ST1: EXEC SPECTRE REPONSE liste (p.202)");
  require(rt_->results().has_modal,
          "ST1: EXEC SPECTRE REPONSE — EXEC MODE prealable requis (p.202)");
  Model& m = model();
  std::size_t i = 3;
  std::vector<Id> all;
  for (const auto& [k, v] : spectre_reps) {
    (void)v;
    all.push_back(k);
  }
  std::vector<Id> list = rt_->id_list(st, i, all);
  if (list.empty()) list = all;
  rt_->spectre_last.clear();
  DampingSpec damp = damping_for();
  int nmodes = static_cast<int>(rt_->results().modal.frequencies_hz.size());

  for (Id id : list) {
    auto it = spectre_reps.find(id);
    if (it == spectre_reps.end()) {
      warn("ST1: EXEC SPECTRE REPONSE — reponse " + std::to_string(id) + " non definie");
      continue;
    }
    const SpectreResponseDef& d = it->second;
    auto sit = spectra.find(d.spectrum);
    require(sit != spectra.end(),
            "ST1: SPECTRE REPONSE " + std::to_string(id) + " — SPECTRE " +
                std::to_string(d.spectrum) + " non defini (p.201)");
    ComboRule rule = d.rule;
    if (rule == ComboRule::Cqc) {
      bool zero = damp.kind == DampingSpec::Kind::Modal && damp.zeta == 0;
      if (zero) {
        warn("ST1: SPECTRE REPONSE CQC sans amortissement — degenere en SRSS (p.242)");
        rule = ComboRule::Srss;
      }
    }
    // The POND directions excite TOGETHER (ST1, Exemple6: POND X 1 Y 1 and
    // X 1 Y -1 give Mz 22.45 and 4468 at the root): each mode's signed peak
    // response is summed over the weighted directions, THEN the modes are
    // combined by the rule (LIN: summed). Bars and reactions are recovered
    // mode by mode from the modal fields (load-free static recovery).
    const int nst = rt_->etude_pas_hint > 0 ? std::max(2, static_cast<int>(rt_->etude_pas_hint) + 1)
                                            : 11;
    std::vector<SpectrumResults> per_dir;
    std::vector<double> wts;
    // SPECTRE BETA (p.202, SPATIALE only): the spectrum's local Ox makes the
    // angle beta with the global X in the horizontal plane. A mode's signed
    // peak is linear in the excitation direction, so the local X (Y) branch
    // excites global X and Y with weights (cos, sin) ((-sin, cos)).
    const double beta = rt_->option_ == 2 ? sit->second.beta : 0.0;
    for (int k = 0; k < 3; ++k) {
      if (d.pond[static_cast<std::size_t>(k)] == 0) continue;
      // EUROCODE HORIZONTAL / VERTICAL: in PLANE the frame's Y axis is the
      // vertical one (Exemple6: the VERTICAL branch on POND Y)
      const int branch = sit->second.eurocode && rt_->option_ == 0 && k == 1 ? 2 : k;
      std::vector<std::pair<int, double>> globals = {{k, 1.0}};
      if (beta != 0 && k < 2)
        globals = k == 0 ? std::vector<std::pair<int, double>>{{0, std::cos(beta)}, {1, std::sin(beta)}}
                         : std::vector<std::pair<int, double>>{{0, -std::sin(beta)}, {1, std::cos(beta)}};
      for (const auto& [g, c] : globals) {
        if (c == 0) continue;
        SpectrumOptions so;
        so.dir = g;
        so.rule = rule;
        so.n_modes = std::max(1, nmodes);
        so.stations = nst;
        so.stations_abs = rt_->etude_stations_ref();  // the static ones' sections
        per_dir.push_back(solve_spectrum(m, sit->second.build(branch), damp, so));
        wts.push_back(c * d.pond[static_cast<std::size_t>(k)]);
      }
    }
    require(!per_dir.empty(), "ST1: SPECTRE REPONSE " + std::to_string(id) +
                                  " — aucune direction pondere (POND X/Y/Z, p.201)");
    const std::size_t nm = per_dir.front().modal_disp.size();
    const std::vector<double>& om = rt_->results().modal.omegas;
    const std::vector<double>& zt = per_dir.front().zetas;
    // the signed per-mode responses, directions summed
    std::vector<std::vector<Vec6>> fld(nm, std::vector<Vec6>(m.nodes.size(), Vec6::Zero()));
    std::vector<std::vector<BarResult>> bars(nm);
    std::vector<std::vector<Reaction>> reac(nm);
    auto scale_bars = [](std::vector<BarResult>& v, double w) {
      for (BarResult& br : v)
        for (Station& q : br.stations) {
          q.N *= w; q.Vy *= w; q.Vz *= w; q.T *= w; q.My *= w; q.Mz *= w;
          q.defl *= w; q.rot *= w;
        }
    };
    for (std::size_t i = 0; i < nm; ++i)
      for (std::size_t k = 0; k < per_dir.size(); ++k) {
        const SpectrumResults& sr = per_dir[k];
        const double w = wts[k];
        for (std::size_t n = 0; n < fld[i].size() && n < sr.modal_disp[i].size(); ++n)
          fld[i][n] += w * sr.modal_disp[i][n];
        if (i < sr.modal_bars.size()) {
          if (bars[i].empty()) {
            bars[i] = sr.modal_bars[i];
            scale_bars(bars[i], w);
          } else {
            for (std::size_t e = 0; e < bars[i].size(); ++e)
              for (std::size_t q = 0; q < bars[i][e].stations.size(); ++q) {
                Station& x = bars[i][e].stations[q];
                const Station& y = sr.modal_bars[i][e].stations[q];
                x.N += w * y.N; x.Vy += w * y.Vy; x.Vz += w * y.Vz; x.T += w * y.T;
                x.My += w * y.My; x.Mz += w * y.Mz; x.defl += w * y.defl; x.rot += w * y.rot;
              }
          }
        }
        if (i < sr.modal_reactions.size()) {
          if (reac[i].empty()) {
            reac[i] = sr.modal_reactions[i];
            for (Reaction& r : reac[i]) r.R *= w;
          } else {
            for (std::size_t r = 0; r < reac[i].size() && r < sr.modal_reactions[i].size(); ++r)
              reac[i][r].R += w * sr.modal_reactions[i][r].R;
          }
        }
      }
    // the modal combination of one quantity
    auto combine = [&](auto value_of) {
      if (d.linear) {
        double s = 0;
        for (std::size_t i = 0; i < nm; ++i) s += value_of(i);
        return s;
      }
      double s = 0;
      for (std::size_t i = 0; i < nm; ++i)
        for (std::size_t j = 0; j < nm; ++j) {
          if (i != j && rule == ComboRule::Srss) continue;
          double rho = 1.0;
          if (i != j) {
            if (rule == ComboRule::Cqc) {
              rho = cqc_rho(om[i], om[j], zt[i], zt[j]);
            } else if (rule == ComboRule::DoubleSum) {
              rho = double_sum_rho(om[i], om[j], zt[i], zt[j], 20.0);
            } else {  // TenPercent: closely spaced modes added absolutely
              const double lo = std::min(om[i], om[j]), hi = std::max(om[i], om[j]);
              s += (hi - lo) <= 0.10 * lo ? std::fabs(value_of(i) * value_of(j)) : 0.0;
              continue;
            }
          }
          s += rho * value_of(i) * value_of(j);
        }
      return std::sqrt(std::max(0.0, s));
    };
    std::vector<Vec6> comb(m.nodes.size(), Vec6::Zero());
    for (std::size_t n = 0; n < comb.size(); ++n)
      for (int c = 0; c < 6; ++c) comb[n](c) = combine([&](std::size_t i) { return fld[i][n](c); });
    CaseResult cr;
    cr.id = id;
    cr.name = "SPECTRE " + std::to_string(id);
    cr.node_disp = comb;
    if (nm > 0 && !bars.front().empty()) {
      cr.bars = bars.front();
      for (std::size_t e = 0; e < cr.bars.size(); ++e)
        for (std::size_t q = 0; q < cr.bars[e].stations.size(); ++q) {
          Station& st = cr.bars[e].stations[q];
          auto one = [&](auto get) {
            return combine([&](std::size_t i) { return get(bars[i][e].stations[q]); });
          };
          st.N = one([](const Station& x) { return x.N; });
          st.Vy = one([](const Station& x) { return x.Vy; });
          st.Vz = one([](const Station& x) { return x.Vz; });
          st.T = one([](const Station& x) { return x.T; });
          st.My = one([](const Station& x) { return x.My; });
          st.Mz = one([](const Station& x) { return x.Mz; });
          for (int c = 0; c < 3; ++c) {
            st.defl(c) = one([c](const Station& x) { return x.defl(c); });
            st.rot(c) = one([c](const Station& x) { return x.rot(c); });
          }
        }
    }
    if (nm > 0 && !reac.front().empty()) {
      // the restrained DOFs only: a free DOF of a support node reads K u =
      // its inertia force, not a reaction (ST1: 0 on an APPUI without DOF)
      cr.reactions = reac.front();
      std::unordered_map<Id, std::vector<const Support*>> supports_of;  // per node, in order
      for (const Support& sp : m.supports) supports_of[sp.node].push_back(&sp);
      for (std::size_t r = 0; r < cr.reactions.size(); ++r) {
        std::array<bool, 6> held{};
        for (const Support* sp : supports_of[cr.reactions[r].node])
          for (std::size_t c = 0; c < 6; ++c)
            held[c] = held[c] || sp->fixed[c] || sp->spring[c] != 0;
        for (int c = 0; c < 6; ++c)
          cr.reactions[r].R(c) =
              held[static_cast<std::size_t>(c)]
                  ? combine([&](std::size_t i) { return r < reac[i].size() ? reac[i][r].R(c) : 0.0; })
                  : 0.0;
      }
    }
    // ST1 treats the spectral response as an ordinary LOAD CASE and tells the
    // user to fold the +/- effect through an ENV with +1/-1 coefficients
    // (manual p.201) — so the stored case must remain envelope-usable. The
    // combined field is NOT physically superposable across distinct dynamic
    // cases; the ENV bounding (per-corner scaling) is the sanctioned use.
    cr.non_superposable = false;
    rt_->results().spectre[id] = cr;
    if (std::find(rt_->spectre_last.begin(), rt_->spectre_last.end(), id) == rt_->spectre_last.end())
      rt_->spectre_last.push_back(id);
  }
}

void St1Domains::exec_historique(const Statement& st) {
  Model& m = model();
  std::size_t i = 2;
  bool keep_time = true;
  int family = 0;  // 0 accelerogramme, 1 convoi, 2 charge
  while (i < st.size()) {
    std::string k = kw(st, i);
    if (k == "TEMPS") {
      keep_time = kw(st, i + 1) != "NON";
      i += 2;
      continue;
    }
    if (k == "CONVOI") { family = 1; ++i; continue; }
    if (k == "CHARGE") { family = 2; ++i; continue; }
    break;
  }
  std::vector<Id> all;
  if (family == 0)
    for (const auto& [k, v] : accels) { (void)v; all.push_back(k); }
  else if (family == 1)
    for (const auto& [k, v] : dyn_convois) { (void)v; all.push_back(k); }
  else
    for (const auto& [k, v] : dyn_charges) { (void)v; all.push_back(k); }
  std::vector<Id> list = rt_->id_list(st, i, all);
  if (list.empty()) list = all;

  DampingSpec damp = damping_for();
  bool modal_route = damp.kind != DampingSpec::Kind::Rayleigh &&
                     !(damp.kind == DampingSpec::Kind::Modal && damp.zeta == 0);
  if (modal_route)
    require(rt_->results().has_modal,
            "ST1: EXEC HISTORIQUE — un amortissement MODE/ENERGIE impose la "
            "superposition modale : EXEC MODE prealable requis (p.204)");

  Assembled asm_ = assemble_elements(m);
  const DofMap& dofs = asm_.dofs;
  // the states of the time series, at the bar ends: "il n'y a pas d'analyse
  // temporelle des sections d'etude [...] efforts et deplacements pour toutes
  // les barres en debut et fin" (p.208)
  const CaseLoads cl0 = lower_case(m, asm_, LoadCase{});
  StaticOptions sopt;
  sopt.stations_per_bar = 2;

  auto run = [&](Id id, const std::string& name, int n_steps, double dt,
                 const std::function<VecX(int)>& force, const IntegrationScheme& sch) {
    TimeHistoryResults th = modal_route
                                ? time_history_modal(m, rt_->results().modal, n_steps, dt, force,
                                                     damp, sch)
                                : time_integrate(m, n_steps, dt, force, damp, sch,
                                                 rt_->masse_.diagonal());
    HistoryResult h;
    h.id = id;
    h.name = name;
    h.steps = n_steps;
    h.dt = dt;
    h.keep_time = keep_time;
    h.max.id = id;
    h.max.name = name + " MAX";
    h.min.id = id;
    h.min.name = name + " MIN";
    h.max.node_disp.assign(m.nodes.size(), Vec6::Zero());
    h.min.node_disp.assign(m.nodes.size(), Vec6::Zero());
    for (std::size_t s = 0; s < th.u.size(); ++s) {
      std::vector<Vec6> u = expand_free(dofs, th.u[s]);
      double peak = 0;
      for (std::size_t n = 0; n < u.size(); ++n) {
        for (int c = 0; c < 6; ++c) {
          h.max.node_disp[n](c) = std::max(h.max.node_disp[n](c), u[n](c));
          h.min.node_disp[n](c) = std::min(h.min.node_disp[n](c), u[n](c));
        }
        peak = std::max(peak, u[n].head<3>().cwiseAbs().maxCoeff());
      }
      h.peak_disp.push_back(peak);
    }
    // Bounds of the time history: envelope-usable by the core's ENV
    // (per-corner scaling), like the spectral cases — manual p.231.
    h.max.non_superposable = false;
    h.min.non_superposable = false;
    historique_env[id] = {h.max, h.min};
    historique[id] = h;
    St1Runtime::HistorySeries& ser = rt_->historique_series[id];
    ser = St1Runtime::HistorySeries{};
    if (keep_time) {
      // ST1 keeps n states for an n-step history (Exemple7: 17 states, t = 0
      // to 0.16 s, for 18 accelerogram values; each matches keta's state of
      // the same time): the final state is not part of its series.
      const std::size_t ns = th.u.size() > 1 ? th.u.size() - 1 : th.u.size();
      ser.times.reserve(ns);
      ser.states.reserve(ns);
      for (std::size_t s = 0; s < ns; ++s) {
        ser.times.push_back(static_cast<double>(s) * dt);
        CaseResult cr = recover_case(m, asm_, cl0, th.u[s], sopt, id, name);
        cr.non_superposable = false;
        // reactions on the restrained DOFs only: a free DOF of a support
        // node carries its inertia / damping force, not a reaction (ST1:
        // APPUI 2 of Exemple7 reads MZ alone)
        for (Reaction& re : cr.reactions) {
          std::array<bool, 6> held{};
          for (const Support& sp : m.supports)
            if (sp.node == re.node)
              for (std::size_t c = 0; c < 6; ++c)
                held[c] = held[c] || sp.fixed[c] || sp.spring[c] != 0;
          for (std::size_t c = 0; c < 6; ++c)
            if (!held[c]) re.R(static_cast<int>(c)) = 0;
        }
        // ST1's time-series states carry Mz with the sign OPPOSITE to its
        // static editions (Exemple7: dMz/dx = +Ty in RESU HISTORIQUE and in
        // the ENV built on it, -Ty in Exemple3's CHARG; N, Ty, Tz, My keep
        // the static signs): the stored states follow ST1, so its tables and
        // the concomitant states of an ENV read alike.
        for (BarResult& br : cr.bars)
          for (Station& q : br.stations) q.Mz = -q.Mz;
        ser.states.push_back(std::move(cr));
      }
    }
    rt_->historique_env[id] = {h.min, h.max};  // core seam: {min, max}
    rt_->output().line("EXEC HISTORIQUE " + std::to_string(id) + " : " +
                       std::to_string(n_steps) + " pas de " + fmt(dt) + " s");
  };

  auto scheme_of = [](int methode, double alpha) {
    if (methode == 2) return IntegrationScheme::newmark_linear();
    if (methode == 3) return IntegrationScheme::hht(alpha);
    return IntegrationScheme::newmark_average();
  };

  for (Id id : list) {
    if (family == 0) {
      auto it = accels.find(id);
      if (it == accels.end()) { warn("ST1: EXEC HISTORIQUE — ACCELEROGRAMME " +
                                     std::to_string(id) + " non defini"); continue; }
      const AccelDef& d = it->second;
      require(d.dt > 0, "ST1: ACCELEROGRAMME " + std::to_string(id) +
                            " — PAS_TEMPS obligatoire (p.203)");
      require(d.accel.size() >= 2, "ST1: ACCELEROGRAMME " + std::to_string(id) +
                                       " — table ACCELERATION vide (p.204)");
      int n = static_cast<int>(d.accel.size()) - 1;
      std::vector<std::function<VecX(int)>> fs;
      // BETA (p.207, SPATIALE): AX AY are given along the local axes whose Ox
      // makes the angle beta with the global X in the horizontal plane
      const double cb = std::cos(rt_->option_ == 2 ? d.beta : 0.0),
                   sb = std::sin(rt_->option_ == 2 ? d.beta : 0.0);
      for (int c = 0; c < 3; ++c) {
        std::vector<double> a;
        bool any = false;
        for (const auto& r : d.accel) {
          const double g = c == 0 ? cb * r[0] - sb * r[1] : (c == 1 ? sb * r[0] + cb * r[1] : r[2]);
          a.push_back(d.pond * g);
          if (g != 0) any = true;
        }
        // the mass form of MASSE PROPRE (p.69), as EXEC MODE uses it
        if (any) fs.push_back(seismic_force(m, c, a, rt_->masse_.diagonal()));
      }
      require(!fs.empty(), "ST1: ACCELEROGRAMME " + std::to_string(id) +
                               " — accelerations toutes nulles");
      auto force = [fs](int s) {
        VecX v = fs[0](s);
        for (std::size_t k = 1; k < fs.size(); ++k) v += fs[k](s);
        return v;
      };
      run(id, "HISTORIQUE " + std::to_string(id), n, d.dt, force,
          scheme_of(d.methode, d.hht_alpha));
      rt_->historique_series[id].title = d.title;
      continue;
    }
    if (family == 1) {
      auto it = dyn_convois.find(id);
      if (it == dyn_convois.end()) { warn("ST1: EXEC HISTORIQUE CONVOI " +
                                          std::to_string(id) + " non defini"); continue; }
      const DynConvoiDef& d = it->second;
      require(d.dt > 0 && d.duration > 0 && d.speed != 0,
              "ST1: DYNAMIQUE CONVOI " + std::to_string(id) +
                  " — PAS_TEMPS / VITESSE / TEMPS obligatoires (p.205)");
      auto cit = convs.find(d.convoi);
      require(cit != convs.end(), "ST1: DYNAMIQUE CONVOI — CONVOI " +
                                      std::to_string(d.convoi) + " non defini (p.160)");
      std::vector<Id> bars = has_tablier && !tablier.bars.empty() ? tablier.bars
                                                                  : rt_->all_bar_ids();
      LoadPath path = build_path(m, bars, 11);
      ConvoyPattern pat = convoy_pattern(cit->second, d.pond, has_tablier ? tablier.ldif : 0.0);
      Vec3 dir = rt_->option_ == 0 ? Vec3(0, -1, 0) : Vec3(0, 0, -1);
      int n = std::max(1, static_cast<int>(d.duration / d.dt));
      auto force = moving_convoy_force(m, path, pat, dir, d.speed, d.dt);
      warn("ST1: DYNAMIQUE CONVOI (p.205) — un seul camion, centre sur les barres du "
           "tablier ; le positionnement transversal (ZONE / YL) est ignore, conforme "
           "au manuel");
      run(id, "CONVOI " + std::to_string(id), n, d.dt, force, scheme_of(d.methode, d.hht_alpha));
      continue;
    }
    auto it = dyn_charges.find(id);
    if (it == dyn_charges.end()) { warn("ST1: EXEC HISTORIQUE CHARGE " +
                                        std::to_string(id) + " non defini"); continue; }
    const DynChargeDef& d = it->second;
    require(d.dt > 0 && d.duration > 0,
            "ST1: DYNAMIQUE CHARGE " + std::to_string(id) +
                " — PAS_TEMPS / TEMPS obligatoires (p.207)");
    require(rt_->results().has_modal,
            "ST1: DYNAMIQUE CHARGE — le signe de la deformee du mode " +
                std::to_string(d.mode) + " est requis : EXEC MODE prealable (p.208)");
    const ModalResults& mo = rt_->results().modal;
    std::size_t mi = static_cast<std::size_t>(std::max(1, d.mode)) - 1;
    require(mi < mo.modes.size(), "ST1: DYNAMIQUE CHARGE MODE " + std::to_string(d.mode) +
                                      " — mode absent de EXEC MODE");
    double om0 = d.f0 > 0 ? 2 * kPi * d.f0 : mo.omegas[mi];
    // One lowered RHS per BARRE line, signed by the mode shape (p.208).
    std::vector<VecX> rhs;
    std::vector<double> phi, coef;
    for (std::size_t li = 0; li < d.lines.size(); ++li) {
      const DynChargeDef::BarLine& bl = d.lines[li];
      LoadCase lc;
      lc.id = 1;
      lc.name = "dyn";
      for (Id b : bl.bars) {
        // the sign of the work of the line load on the mode, forces and
        // (distributed) moments alike (BARRE ... MX MY MZ, p.211): its
        // consistent nodal loads against the mode's nodal displacements
        double s = 0;
        {
          LoadCase one;
          one.id = 1;
          Load lb = Load::udl(b, bl.F, bl.frame);
          lb.M = bl.M;
          one.loads.push_back(lb);
          const CaseLoads c1 = lower_case(m, asm_, one);
          const DofMap& dm = asm_.dofs;
          for (std::size_t nidx = 0; nidx < mo.modes[mi].size(); ++nidx)
            for (int c = 0; c < 6; ++c) {
              const int fi = dm.free_ix[static_cast<std::size_t>(dm.full(static_cast<int>(nidx), c))];
              if (fi >= 0) s += c1.rhs(fi) * mo.modes[mi][nidx](c);
            }
        }
        double sign = s >= 0 ? 1.0 : -1.0;
        if (bl.F.norm() == 0 && bl.M.norm() == 0) sign = 1.0;
        Load l = Load::udl(b, sign * bl.F, bl.frame);
        l.M = sign * bl.M;
        lc.loads.push_back(l);
      }
      CaseLoads cl = lower_case(m, asm_, lc);
      rhs.push_back(cl.rhs);
      phi.push_back(li < d.phi.size() ? d.phi[li] : (d.phi.empty() ? 0.0 : d.phi.back()));
      coef.push_back(li < d.coef.size() ? d.coef[li] : (d.coef.empty() ? 1.0 : d.coef.back()));
    }
    require(!rhs.empty(), "ST1: DYNAMIQUE CHARGE — aucune ligne BARRE (p.209)");
    double dt = d.dt;
    auto force = [rhs, phi, coef, om0, dt](int s) {
      VecX v = VecX::Zero(rhs[0].size());
      double t = s * dt;
      for (std::size_t k = 0; k < rhs.size(); ++k)
        v += (coef[k] * std::sin(om0 * t + phi[k])) * rhs[k];
      return v;
    };
    int n = std::max(1, static_cast<int>(d.duration / d.dt));
    run(id, "CHARGE " + std::to_string(id), n, d.dt, force, scheme_of(d.methode, d.hht_alpha));
  }
}

void St1Domains::exec_incendie(const Statement& st) {
  require(!fire_bars.empty(), "ST1: EXEC INCENDIE — aucune barre au feu "
                              "(CARA ib SECTION i, p.224)");
  require(rt_->option_ != 1,
          "ST1: EXEC INCENDIE — le calcul au feu n'est pas disponible en option GRILL (p.213)");
  Model base = model();
  std::size_t i = 2;
  std::vector<Id> all;
  for (const auto& [k, v] : incendies) { (void)v; all.push_back(k); }
  std::vector<Id> list = rt_->id_list(st, i, all);
  if (list.empty()) list = all;

  // Dense exposure table + ST1 id -> index map.
  std::vector<Exposure> ex;
  std::map<Id, int> ex_ix;
  for (const auto& [k, v] : exposures) {
    ex_ix[k] = static_cast<int>(ex.size());
    ex.push_back(v);
  }
  auto ix_of = [&](Id id) {
    auto it = ex_ix.find(id);
    return it == ex_ix.end() ? -1 : it->second;
  };

  for (Id id : list) {
    auto iit = incendies.find(id);
    if (iit == incendies.end()) {
      warn("ST1: EXEC INCENDIE — incendie " + std::to_string(id) + " non defini");
      continue;
    }
    const IncendieDef& inc = iit->second;
    auto fit = feux.find(inc.feu);
    require(fit != feux.end(),
            "ST1: INCENDIE " + std::to_string(id) + " — FEU " + std::to_string(inc.feu) +
                " non defini (p.225)");
    const FeuDef& feu = fit->second;

    Model m = base;
    m.cases.clear();
    for (Id cid : rt_->charg_order_) {
      const ChargDef& c = rt_->chargs_.at(cid);
      if (c.is_comb) continue;
      LoadCase lc;
      lc.id = c.id;
      lc.name = c.name;
      lc.loads = rt_->lower_loads(c.loads);
      m.cases.push_back(lc);
    }
    // A CHARG COMB preload (p.225) is superposed into one real case.
    if (inc.charg != 0 && rt_->chargs_.count(inc.charg) && rt_->chargs_.at(inc.charg).is_comb) {
      const ChargDef& c = rt_->chargs_.at(inc.charg);
      LoadCase lc;
      lc.id = c.id;
      lc.name = c.name;
      for (const auto& [pid, coef] : c.parts) {
        auto pc = rt_->chargs_.find(pid);
        if (pc == rt_->chargs_.end()) continue;
        for (Load l : rt_->lower_loads(pc->second.loads)) {
          l.F *= coef;
          l.M *= coef;
          l.F2 *= coef;
          l.dt *= coef;
          l.gy *= coef;
          l.gz *= coef;
          l.g *= coef;
          l.strain *= coef;
          l.kappa_y *= coef;
          l.kappa_z *= coef;
          lc.loads.push_back(l);
        }
      }
      m.cases.push_back(lc);
    }
    m.validate();

    FireProblem prob;
    prob.exposures = ex;
    prob.duration = feu.duree;
    prob.t_init = feu.temp;
    prob.local = inc.local;
    prob.preload_case = inc.charg;
    prob.option = rt_->option_ == 0 ? ModelOption::Plane : ModelOption::Spatial;
    for (const auto& [bid, fb] : fire_bars) {
      auto sit = fire_sections.find(fb.section);
      require(sit != fire_sections.end(),
              "ST1: CARA " + std::to_string(bid) + " SECTION " + std::to_string(fb.section) +
                  " — section feu non definie (p.218)");
      FireBarSpec spec;
      spec.bar = bid;
      spec.section = sit->second.section;
      spec.section.expo_vy = ix_of(sit->second.expo_ids[0]);
      spec.section.expo_wy = ix_of(sit->second.expo_ids[1]);
      spec.section.expo_vz = ix_of(sit->second.expo_ids[2]);
      spec.section.expo_wz = ix_of(sit->second.expo_ids[3]);
      for (int& sg : spec.section.seg_exposure) sg = ix_of(static_cast<Id>(sg));
      if (feu.maille > 0) spec.section.maille = feu.maille;
      auto mit = mats.find(sit->second.mat);
      require(mit != mats.end() && mit->second.fire,
              "ST1: SECTION " + std::to_string(sit->second.id) + " MAT " +
                  std::to_string(sit->second.mat) +
                  " — materiau feu (MAT ... E FEU EC2) non defini (p.214)");
      spec.material = mit->second.fire_mat;
      // firechain's hot stress amplitudes are kc(theta)*fck: hand it fck in
      // the DECK's unit system so the hot EA/EI swap consistently against the
      // user-unit cold model (the thermal side stays SI — temperatures only).
      spec.material.fck /= mit->second.stress_scale();
      if (fb.armature != 0) {
        auto ait = armatures.find(fb.armature);
        require(ait != armatures.end(), "ST1: CARA " + std::to_string(bid) + " ARMATURE " +
                                            std::to_string(fb.armature) + " non definie (p.223)");
        spec.rebars = ait->second.layers;
      }
      // `SECTION i DALLE` (p.218): the slab's heat problem is 1-D through the
      // thickness (SUP/INF exposed, lateral faces adiabatic). heat2d meshes a
      // Slab as ONE cell column, which leaves the hot EIz identically zero —
      // so lower a DALLE WITH a width onto the strictly equivalent RECTANGLE:
      // same SUP/INF exposures, lateral faces unexposed (= adiabatic, i.e.
      // the exact 1-D field extruded over LARGEUR), correct layered EIz.
      if (spec.section.kind == FireSection::Kind::Slab) {
        require(spec.section.width > 0,
                "ST1: SECTION " + std::to_string(sit->second.id) + " DALLE (p.218) — le "
                "calcul au feu d'une dalle demande sa LARGEUR (equivalent RECTANGLE a "
                "faces laterales adiabatiques)");
        spec.section.kind = FireSection::Kind::Rectangle;
        spec.section.expo_vz = -1;
        spec.section.expo_wz = -1;
      }
      spec.section.validate();
      prob.bars.push_back(spec);
    }
    FireOutcome out = solve_fire(m, prob, fire_opt);
    fire_outcomes[id] = out;
    CaseResult cr = out.result;
    // Y VAR: ST1 prints the member forces about the hot centroid
    // (yc, zc) of each exposed bar, M_c = M_0 - (0, yc, zc) x (N, Vy, Vz)
    // (a column head free in rotation: Mz = N yc; in space: Mx = -Vz yc).
    if (!fire_opt.fold_centroid_ecc)
      for (const HotSectionState& hs : out.states)
        for (BarResult& br : cr.bars)
          if (br.bar == hs.bar)
            for (Station& stn : br.stations) {
              stn.T += -hs.yc * stn.Vz + hs.zc * stn.Vy;
              stn.My -= hs.zc * stn.N;
              stn.Mz += hs.yc * stn.N;
            }
    cr.id = id;
    cr.name = inc.title.empty() ? ("INCENDIE " + std::to_string(id)) : inc.title;
    cr.non_superposable = true;
    rt_->results().incendie[id] = cr;
    rt_->output().line("EXEC INCENDIE " + std::to_string(id) + " : " +
                       std::to_string(out.iterations) + " iteration(s), duree " +
                       fmt(feu.duree) + " min");
    // `RESU INCENDIE i BARRE b CARA CHAUD (TEMP)` (p.227) has no slot in the
    // core's ResultsStore: render the hot-section table here.
    for (const HotSectionState& hs : out.states) {
      std::string ln = "INCENDIE " + std::to_string(id) + " BARRE " +
                       std::to_string(hs.bar) + " CARA CHAUD  EA " + fmt(hs.EA) + "  EIY " +
                       fmt(hs.EIy) + "  EIZ " + fmt(hs.EIz) + "  YC " + fmt(hs.yc) + "  ZC " +
                       fmt(hs.zc) + "  TEMP MAX " + fmt(hs.theta_max) + "  MOY " +
                       fmt(hs.theta_mean) + "  ARM " + fmt(hs.theta_rebar_max);
      hot_cara_lines.push_back(ln);
      rt_->output().line(ln);
    }
  }
}

// ===========================================================================
// GET TENSION CABLE (p.276)
// ===========================================================================
void St1Domains::get_tension(const Statement& st) {
  Env& env = rt_->env_;
  env.scalars["$err"] = 0.0;
  try {
    require(kw(st, 2) == "CABLE", "ST1: GET TENSION CABLE i SEC n <CHARG,ETAT,PHASAGE> j");
    Id cab = static_cast<Id>(num(st[3]));
    auto it = cable_tension.find(cab);
    require(it != cable_tension.end(),
            "ST1: GET TENSION CABLE " + std::to_string(cab) +
                " — aucune table de tension (EXEC CHARG avec une action CHARG CABLE requis)");
    const CableTension& t = it->second;
    std::size_t sec = 1;
    for (std::size_t i = 4; i + 1 < st.size(); ++i)
      if (kw(st, i) == "SEC" && is_num_tok(st[i + 1]))
        sec = static_cast<std::size_t>(num(st[i + 1]));
    require(sec >= 1 && sec <= t.s.size(), "ST1: GET TENSION CABLE — section hors bornes");
    std::size_t k = sec - 1;
    env.scalars["$nsec"] = static_cast<double>(t.s.size());
    env.scalars["$xl"] = t.s[k];
    env.scalars["$sigp"] = t.sigma[k];
    env.scalars["$dsigp"] = t.sigma0 - t.sigma[k];
    env.scalars["$x"] = t.pos[k].x();
    env.scalars["$y"] = t.pos[k].y();
    env.scalars["$z"] = t.pos[k].z();
  } catch (const Error& e) {
    env.scalars["$err"] = 1.0;
    warn(std::string("ST1: GET TENSION en erreur — ") + e.what());
  }
}

// ===========================================================================
// LISTER (p.236-239)
// ===========================================================================
bool St1Domains::lister(const Statement& st, const std::string& what) {
  Output& o = rt_->output();
  if (what == "PREC") {
    for (const auto& [id, p] : precs)
      o.line("PREC " + std::to_string(id) + " '" + p.title + "' SECTION " + fmt(p.ap) +
             " TENSION " + fmt(p.sigma) + " E " + fmt(p.ep) + " RECUL " + fmt(p.recul) +
             " F " + fmt(p.f) + " PHI " + fmt(p.eurocode_form ? p.kphi : p.phi));
    return true;
  }
  if (what == "CABLE") {
    for (const auto& [id, c] : cables)
      o.line("CABLE " + std::to_string(id) + " '" + c.title + "' PREC " +
             std::to_string(c.prec) + " BARRE " + std::to_string(c.bars.size()) +
             " TRACE " + std::to_string(c.trace.size()) + " points PAS_CABLE " +
             fmt(c.pas_cable));
    return true;
  }
  if (what == "CARA" && kw(st, 2) == "CABLE") {
    for (const auto& [id, t] : cable_tension) {
      o.line("CARA CABLE " + std::to_string(id) + " (" + std::to_string(t.s.size()) +
             " sections)");
      for (std::size_t k = 0; k < t.s.size(); ++k)
        o.line("  XS " + fmt(t.s[k]) + "  P " + fmt(t.P[k]) + "  SIGP " + fmt(t.sigma[k]));
    }
    return true;
  }
  if (what == "TABLIER") {
    if (!has_tablier) {
      o.line("TABLIER : (aucun)");
      return true;
    }
    o.line("TABLIER CLASSE " + std::to_string(tablier.bridge_class()) +
           (tablier.trafic_class ? " (TRAFIC)" : "") + " LARGEUR " + fmt(tablier.width()) +
           " PAS_LIGNE " + fmt(tablier.pas_ligne) + " LDIF " + fmt(tablier.ldif) +
           " BARRE " + std::to_string(tablier.bars.size()));
    for (std::size_t k = 0; k < tablier.zones.size(); ++k)
      o.line("  ZONE_TRANS " + std::to_string(k + 1) + " LARG " + fmt(tablier.zones[k].width) +
             " VOIE " + std::to_string(tablier.zones[k].nvoies));
    return true;
  }
  if (what == "ALG") {
    for (const auto& [id, a] : algs) {
      o.line("ALG " + std::to_string(id) + " '" + a.name + "' LARG VOIE " + fmt(a.larg_voie));
      for (const auto& [l, dd] : a.dens) o.line("  DENS " + fmt(l) + "  " + fmt(dd));
    }
    return true;
  }
  if (what == "CONV" || what == "CONVOI") {
    for (const auto& [id, c] : convs)
      o.line("CONV " + std::to_string(id) + " '" + c.name + "' " +
             std::to_string(c.axles.size()) + " essieu(x), total " + fmt(c.total_t()) +
             " t, LONG " + fmt(c.longueur) + " LARG " + fmt(c.largeur) + " MAX_CAM " +
             std::to_string(c.max_cam));
    return true;
  }
  if (what == "SURCH") {
    for (const auto& [id, s] : surchs)
      o.line("SURCH " + std::to_string(id) + " charge " + std::to_string(static_cast<int>(s.load)) +
             (s.convoy.empty() ? "" : " '" + s.convoy + "'") + " POND " + fmt(s.pond) +
             " PAS_CONV " + fmt(s.pas_conv) + " SENS " + std::to_string(s.sens) + " SEUIL " +
             fmt(s.seuil));
    return true;
  }
  if (what == "PHASAGE") {
    for (const auto& [id, p] : phasages) {
      o.line("PHASAGE " + std::to_string(id) + " '" + p.title + "' " +
             std::to_string(p.plan.phases.size()) + " phase(s), " +
             std::to_string(p.etats.size()) + " etat(s)");
      for (const Phase& ph : p.plan.phases)
        o.line("  " + ph.name + " DATE " + fmt(ph.time) + " ACTIVER " +
               std::to_string(ph.activate_bars.size()) + " barre(s), CHARG " +
               std::to_string(ph.cases.size()));
    }
    return true;
  }
  if (what == "AMORTISSEMENT") {
    for (const auto& [id, d] : dampings) {
      std::string k = d.spec.kind == DampingSpec::Kind::Rayleigh
                          ? ("RAYLEIGH MASSE " + fmt(d.spec.alpha) + " RIGIDITE " + fmt(d.spec.beta))
                          : (d.spec.kind == DampingSpec::Kind::Modal
                                 ? ("ENERGIE KSI " + fmt(d.spec.zeta))
                                 : ("MODE (" + std::to_string(d.spec.per_mode.size()) + " valeurs)"));
      o.line("AMORTISSEMENT " + std::to_string(id) + " " + k + " POND " + fmt(d.pond));
    }
    return true;
  }
  if (what == "SPECTRE") {
    for (const auto& [id, s] : spectra) {
      // The EN 1998-1 parameters live per axis (`Ec8Branch`), so the listing
      // shows the first declared branch, named by its axis (p.199).
      std::size_t bi = 0;
      while (bi < s.ec8.size() && !s.ec8[bi].set) ++bi;
      const bool has_branch = bi < s.ec8.size();
      const Ec8Branch& b = has_branch ? s.ec8[bi] : s.ec8[0];
      static const char* kAxis[3] = {"X", "Y", "Z"};
      o.line("SPECTRE " + std::to_string(id) + (s.eurocode ? " EUROCODE" : " ACCELERATION") +
             (s.reduit ? " REDUIT" : "") + " POND " + fmt(s.pond) +
             (s.eurocode ? (has_branch ? (std::string(" A") + kAxis[bi] + " " + fmt(b.ag) + " S " +
                                          fmt(b.s_soil) + " PERIODE " + fmt(b.tb) + " " +
                                          fmt(b.tc) + " " + fmt(b.td))
                                       : std::string(" (aucune branche)"))
                         : (" " + std::to_string(s.ax.size()) + " periode(s)")));
    }
    for (const auto& [id, r] : spectre_reps)
      o.line("SPECTRE " + std::to_string(id) + " REPONSE SPECTRE " + std::to_string(r.spectrum) +
             (r.linear ? " LIN" : (r.rule == ComboRule::Cqc ? " CQC" : " SRSS")) + " POND X " +
             fmt(r.pond[0]) + " Y " + fmt(r.pond[1]) + " Z " + fmt(r.pond[2]));
    return true;
  }
  if (what == "HISTORIQUE") {
    for (const auto& [id, h] : historique)
      o.line("HISTORIQUE " + std::to_string(id) + " '" + h.name + "' " +
             std::to_string(h.steps) + " pas de " + fmt(h.dt) + " s");
    for (const auto& [id, a] : accels)
      o.line("ACCELEROGRAMME " + std::to_string(id) + " PAS_TEMPS " + fmt(a.dt) + " " +
             std::to_string(a.accel.size()) + " point(s)");
    return true;
  }
  if (what == "SECTION") {
    for (const auto& [id, s] : fire_sections)
      o.line("SECTION " + std::to_string(id) + " '" + s.title + "' " +
             (s.section.kind == FireSection::Kind::Slab
                  ? "DALLE"
                  : (s.section.kind == FireSection::Kind::Rectangle ? "RECTANGLE" : "PERIMETRE")) +
             " EPAISSEUR " + fmt(s.section.thickness) + " LARGEUR " + fmt(s.section.width) +
             " MAILLE " + fmt(s.section.maille) + " MAT " + std::to_string(s.mat));
    return true;
  }
  if (what == "MATERIAU" || what == "MAT") {
    for (const auto& [id, d] : mats) {
      if (d.fire)
        o.line("MAT " + std::to_string(id) + " '" + d.title + "' E FEU EC2 FCK " +
               fmt(d.fire_mat.fck) + " TRACTION " + (d.fire_mat.tension ? "OUI" : "NON"));
      else
        o.line("MAT " + std::to_string(id) + " '" + d.title + "' loi " +
               std::to_string(static_cast<int>(d.law.code)) + " FC " + fmt(d.law.fc) +
               " E(28) " + fmt(d.law.E(28.0)) + " RETRAIT " + (d.has_retrait ? "OUI" : "NON") +
               " FLUAGE " + (d.has_fluage ? "OUI" : "NON"));
    }
    return true;
  }
  if (what == "EXPOSITION") {
    for (const auto& [id, e] : exposures)
      o.line("EXPOSITION FEU " + std::to_string(id) + " type " +
             std::to_string(static_cast<int>(e.kind)) + " PROTECTION " + fmt(e.protection));
    return true;
  }
  if (what == "ACIER") {
    for (const auto& [id, s] : aciers)
      o.line("ACIER " + std::to_string(id) + " DIAMETRE " + fmt(s.diameter) + " SECTION " +
             fmt(s.bar_area()) + " E " + fmt(s.E) + " FYK " + fmt(s.fyk));
    return true;
  }
  if (what == "ARMATURE") {
    for (const auto& [id, a] : armatures)
      o.line("ARMATURE " + std::to_string(id) + " " + std::to_string(a.layers.size()) + " lit(s)");
    return true;
  }
  return false;
}

// ===========================================================================
// Graphics (p.243-258)
// ===========================================================================
std::vector<const CaseResult*> St1Domains::draw_sources(const Statement& st,
                                                        std::size_t i) const {
  std::vector<const CaseResult*> out;
  std::string kind = kw(st, i);
  if (kind.empty()) return out;
  ++i;
  std::vector<Id> ids = rt_->id_list(const_cast<Statement&>(st), i, {});
  if (kind == "ENV") {
    for (Id id : ids) {
      auto it = rt_->env_results.find(id);
      if (it == rt_->env_results.end()) {
        rt_->warn("ST1: DESS ENV " + std::to_string(id) +
                  " — aucune enveloppe en base (ENV n'est pas cable par cette extension)");
        continue;
      }
      out.push_back(&it->second.state_max);
      out.push_back(&it->second.state_min);
    }
    return out;
  }
  if (kind == "SURCH") {
    for (Id id : ids) {
      auto it = rt_->surch_results.find(id);
      if (it == rt_->surch_results.end()) {
        rt_->warn("ST1: DESS SURCH " + std::to_string(id) + " — EXEC SURCH prealable requis");
        continue;
      }
      for (const SurchTargetResult& t : it->second.sections) {
        if (t.hi.has_state) out.push_back(&t.hi.state);
        if (t.lo.has_state) out.push_back(&t.lo.state);
        break;  // one governing pair is enough for a diagram
      }
    }
    return out;
  }
  for (Id id : ids) {
    const CaseResult* r = rt_->results().find(kind, id);
    if (!r) {
      rt_->warn("ST1: DESS — resultat " + kind + " " + std::to_string(id) + " absent");
      continue;
    }
    out.push_back(r);
  }
  return out;
}

void St1Domains::dess_stmt(const Statement& st, const std::string& head) {
  if (head == "RETOUR") {
    in_dess_ = false;
    return;
  }
  if (head == "ECRAN") {
    draw.screen = screen_from_name(st.size() > 1 ? st[1].s : "w1024_768");
    return;
  }
  if (head == "TRA") {
    if (st.size() > 1 && st[1].str) {
      // keta extension: `TRA 'chemin.ext'` fixes the file SAUV writes to.
      last_drawing_path = st[1].s;
      std::string p = upper(last_drawing_path);
      draw.plotter = p.size() > 4 && p.substr(p.size() - 4) == ".DXF" ? Plotter::Dxf : Plotter::Post;
    } else if (st.size() > 1) {
      draw.plotter = plotter_from_name(st[1].s);
    }
    for (std::size_t k = 1; k < st.size(); ++k) {
      if (kw(st, k) == "COULEUR") draw.monochrome = false;
      if (kw(st, k) == "MONOCHROME") draw.monochrome = true;
    }
    return;
  }
  if (head == "SAUV") {
    std::string path = st.size() > 1 && st[1].str ? st[1].s : last_drawing_path;
    if (path.empty()) {
      last_drawing_path = sauv(draw, ".");
      rt_->output().line("SAUV : " + last_drawing_path);
    } else {
      save_drawing(draw, path);
      last_drawing_path = path;
      rt_->output().line("SAUV : " + path);
    }
    return;
  }
  if (head == "PRINT") {
    warn("ST1: DESS PRINT (p.245) — l'impression Windows n'a pas d'equivalent "
         "headless : utiliser SAUV");
    return;
  }
  if (head == "CLS") {
    cls(draw);
    return;
  }
  if (head == "VUE") {
    Vec3 v = Vec3::Zero();
    int axis = 0;
    for (std::size_t i = 1; i < st.size(); ++i) {
      std::string k = kw(st, i);
      if (k == "X") { axis = 0; continue; }
      if (k == "Y") { axis = 1; continue; }
      if (k == "Z") { axis = 2; continue; }
      if (is_num_tok(st[i]) && axis < 3) v(axis++) = num(st[i]);
    }
    vue(draw, v);
    return;
  }
  if (head == "ZOOM") {
    std::size_t i = 1;
    std::vector<double> v = num_run(st, i);
    if (v.size() >= 3) zoom(draw, v[0], v[1], v[2]);
    else if (v.size() >= 1) zoom(draw, v[0]);
    return;
  }
  if (head == "SELEC") {
    std::size_t i = 1;
    std::string k = kw(st, i);
    if (k == "TOUT") {
      selec(draw, Selection::everything());
      return;
    }
    ++i;
    if (k == "BARRE") selec(draw, Selection::of_bars(rt_->id_list(st, i, rt_->all_bar_ids())));
    else if (k == "NOEUD") selec(draw, Selection::of_nodes(rt_->id_list(st, i, rt_->all_node_ids())));
    else warn("ST1: DESS SELEC — mot '" + k + "' ignore");
    return;
  }
  if (head == "CENTR") {
    std::string k = kw(st, 1);
    if (k.empty() || k == "GEOM") {
      if (k == "GEOM") centr_geom(draw, model());
      else centr(draw, model());
      return;
    }
    double xmax = 0;
    for (std::size_t i = 1; i + 1 < st.size(); ++i)
      if (kw(st, i) == "MAX") xmax = num(st[i + 1]);
    std::string fam = k;
    if (fam == "DEPLA" || fam == "EFFORT" || fam == "CONTR" || fam == "CONTRG" ||
        fam == "PRESS") {
      if (xmax > 0) centr_max(draw, family_from_name(fam), xmax);
      return;
    }
    warn("ST1: DESS CENTR — famille '" + k + "' inconnue (p.248)");
    return;
  }
  if (head == "COUL") {
    std::string fam = kw(st, 1);
    std::string col = kw(st, 2);
    coul(draw, family_from_name(fam), color_from_name(col));
    return;
  }
  if (head == "GEOM") {
    draw_geom(draw, model());
    return;
  }
  if (head == "NUM") {
    NumOptions no;
    bool any = has_kw(st, "NOEUD") || has_kw(st, "BARRE");
    if (any) {
      no.nodes = has_kw(st, "NOEUD");
      no.bars = has_kw(st, "BARRE");
    }
    draw_num(draw, model(), no);
    return;
  }
  if (head == "DEPLA" || head == "EFFORT" || head == "CONTR" || head == "PRESS" ||
      head == "CONTRG") {
    Model& m = model();
    std::size_t i = 1;
    std::string comp;
    if (head == "CONTRG") {
      warn("ST1: DESS CONTRG (p.253) — aucune contrainte generalisee n'est stockee "
           "par le coeur ST1 : le trace est ignore");
      return;
    }
    if (i < st.size() && !st[i].str) {
      std::string k = kw(st, i);
      static const std::set<std::string> kinds = {"CHARG", "COMB",     "ENV",    "SURCH",
                                                  "ETAT",  "PHASAGE",  "SPECTRE", "INCENDIE"};
      if (!kinds.count(k)) {
        comp = k;
        ++i;
      }
    }
    // an optional trailing colour
    Color col = Color::Noir;
    bool has_col = false;
    Statement body;
    for (std::size_t k = 0; k < st.size(); ++k) {
      if (k >= i) {
        std::string t = kw(st, k);
        bool is_color = false;
        try {
          Color c = color_from_name(t);
          col = c;
          is_color = true;
        } catch (const Error&) {
        }
        if (is_color) {
          has_col = true;
          continue;
        }
      }
      body.push_back(st[k]);
    }
    Family fam = head == "DEPLA" ? Family::Depla
                                 : (head == "EFFORT" ? Family::Effort
                                                     : (head == "CONTR" ? Family::Contr : Family::Press));
    if (has_col) coul(draw, fam, col);
    std::vector<const CaseResult*> srcs = draw_sources(body, i);
    for (const CaseResult* r : srcs) {
      if (head == "DEPLA") draw_depla(draw, m, *r, disp_comp_from_name(comp.empty() ? "ALL" : comp));
      else if (head == "EFFORT") draw_effort(draw, m, *r, force_comp_from_name(comp));
      else if (head == "CONTR") draw_contr(draw, m, *r, fiber_from_name(comp));
      else draw_press(draw, m, *r, press_comp_from_name(comp));
    }
    return;
  }
  warn("ST1: DESS — commande '" + head + "' non reconnue (p.243-258)");
}

}  // namespace keta::st1

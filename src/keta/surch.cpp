#include "keta/surch.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <list>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>

namespace keta {

namespace {
// Runs fn(i) for i in [0, n) on the hardware threads (the per-target SURCH
// evaluations are independent and read-only on the shared core); the first
// exception is rethrown on the caller's thread. Results are written by
// index, so the outcome is identical to the sequential loop.
void parallel_for_index(std::size_t n, const std::function<void(std::size_t)>& fn) {
  std::size_t nt = std::max<std::size_t>(1, std::thread::hardware_concurrency());
  nt = std::min(nt, n);
  if (nt <= 1) {
    for (std::size_t i = 0; i < n; ++i) fn(i);
    return;
  }
  std::atomic<std::size_t> next{0};
  std::atomic<bool> stop{false};
  std::exception_ptr err;
  std::mutex mu;
  std::vector<std::thread> pool;
  for (std::size_t t = 0; t < nt; ++t)
    pool.emplace_back([&] {
      for (std::size_t i; !stop && (i = next++) < n;) {
        try {
          fn(i);
        } catch (...) {
          std::lock_guard<std::mutex> lk(mu);
          if (!err) err = std::current_exception();
          stop = true;
        }
      }
    });
  for (std::thread& th : pool) th.join();
  if (err) std::rethrow_exception(err);
}
}  // namespace

// ===========================================================================
// 1. TABLIER (manuel p.126-135)
// ===========================================================================

int st1_lane_count(double width, bool eurocode) {
  require(width >= 0 && finite(width), "TABLIER: largeur de zone non finie ou negative");
  if (eurocode) {
    // EN 1991-2 §4.2.3 tel que tabulé p.164.
    if (width < 5.40) return width > 0 ? 1 : 0;
    if (width < 6.00) return 2;
    return static_cast<int>(width / 3.0);
  }
  // Fasc. 61 titre II art. 2, tel que tabulé p.126.
  if (width < 5.00) return width > 0 ? 1 : 0;
  if (width <= 6.00) return 2;
  return static_cast<int>(width / 3.0);
}

double Tablier::width() const {
  double w = 0;
  for (const TransZone& z : zones) w += z.width;
  return w;
}

std::vector<int> Tablier::resolved_zones(const std::vector<int>& ids) const {
  std::vector<int> out;
  if (zones.empty()) return out;
  if (ids.empty()) {
    for (std::size_t i = 0; i < zones.size(); ++i)
      if (zones[i].width > 0) out.push_back(static_cast<int>(i) + 1);
    return out;
  }
  for (int id : ids) {
    require(id >= 1 && id <= static_cast<int>(zones.size()),
            "SURCH ZONE " + std::to_string(id) + " : zone transversale non definie");
    if (zones[static_cast<std::size_t>(id) - 1].width > 0) out.push_back(id);
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

double Tablier::loaded_width(const std::vector<int>& ids) const {
  double w = 0;
  for (int id : resolved_zones(ids)) w += zones[static_cast<std::size_t>(id) - 1].width;
  return w;
}

std::vector<DeckLane> Tablier::lanes(const std::vector<int>& ids) const {
  return lanes(ids, trafic_class);
}

std::vector<DeckLane> Tablier::lanes(const std::vector<int>& ids, bool eurocode) const {
  const bool trafic_class = eurocode;  // the lane rules of the family asked for
  std::vector<DeckLane> out;
  std::vector<int> zs = resolved_zones(ids);
  if (zs.empty()) return out;
  const double half = width() / 2.0;
  // Abscisse transversale du bord gauche de chaque zone (le tablier est centré
  // sur l'axe des barres, zone 1 à gauche — p.128).
  std::vector<double> left(zones.size(), 0.0);
  double acc = -half;
  for (std::size_t i = 0; i < zones.size(); ++i) {
    left[i] = acc;
    acc += zones[i].width;
  }
  for (int id : zs) {
    const TransZone& z = zones[static_cast<std::size_t>(id) - 1];
    int n = z.nvoies > 0 ? z.nvoies : st1_lane_count(z.width, trafic_class);
    if (n <= 0) continue;
    // Fasc.61 : la « largeur réelle de la voie » partage la zone (a2 = v0/v).
    // Eurocode : voies de 3,00 m (ou L/2 entre 5,40 et 6,00), le reste est
    // une aire résiduelle.
    double lw = z.width / n;
    if (trafic_class) lw = (z.width >= 5.40 && z.width < 6.00) ? z.width / 2 : 3.0;
    for (int j = 0; j < n; ++j) {
      DeckLane l;
      l.width = lw;
      l.y = left[static_cast<std::size_t>(id) - 1] + (j + 0.5) * lw;
      l.zone = id;
      out.push_back(l);
    }
  }
  return out;
}

namespace {
// Profile of one bar (nullptr = coefficient 1: not in any listekb, p.134).
const RepTransBlock* rep_block_of(const Tablier& t, Id bar) {
  const RepTransBlock* blk = nullptr;
  for (const RepTransBlock& b : t.rep_trans)
    for (Id id : b.bars)
      if (id == bar) blk = &b;  // a later KBARRE line redefines the bar
  return blk;
}
const std::vector<double> kUnitProfile = {1.0};
const std::vector<double>& rep_profile(const RepTransBlock* b) { return b ? b->k : kUnitProfile; }
}  // namespace

bool Tablier::rep_trans_uniform() const {
  if (rep_trans.empty() || bars.empty()) return true;
  const std::vector<double>& k0 = rep_profile(rep_block_of(*this, bars.front()));
  for (Id b : bars)
    if (rep_profile(rep_block_of(*this, b)) != k0) return false;
  return true;
}

double Tablier::rep_coef(Id study_bar, double y, double half_width) const {
  if (rep_trans.empty()) return 1.0;
  const bool on_deck = std::find(bars.begin(), bars.end(), study_bar) != bars.end();
  const RepTransBlock* blk = nullptr;
  if (!bars.empty() && (!on_deck || rep_trans_uniform()))
    blk = rep_block_of(*this, bars.front());  // the deck's (uniform) profile
  else
    blk = rep_block_of(*this, study_bar);
  if (blk == nullptr) return 1.0;  // p.134 : coefficient 1 hors de listekb
  const std::vector<double>& k = blk->k;
  require(!k.empty(), "REP_TRANS: liste de coefficients vide");
  if (k.size() == 1) return k[0];
  const double half = width() / 2.0;
  require(half > 0, "REP_TRANS: ZONE_TRANS doit etre defini (p.130)");
  auto at = [&](double yy) {
    double t = (std::clamp(yy, -half, half) + half) / (2 * half);  // 0..1
    double u = t * static_cast<double>(k.size() - 1);
    std::size_t i = static_cast<std::size_t>(std::floor(u));
    if (i + 1 >= k.size()) return k.back();
    double f = u - static_cast<double>(i);
    return k[i] + f * (k[i + 1] - k[i]);
  };
  if (half_width <= 0) return at(y);
  // Impact diffusé : moyenne du coefficient sur l'empreinte, intégrale exacte
  // de la ligne brisée (trapèzes entre ses nœuds, constante hors du tablier).
  double a = y - half_width, b = y + half_width;
  std::vector<double> cuts{a, b};
  for (std::size_t i = 0; i < k.size(); ++i) {
    const double yi = -half + 2 * half * static_cast<double>(i) / static_cast<double>(k.size() - 1);
    if (yi > a && yi < b) cuts.push_back(yi);
  }
  std::sort(cuts.begin(), cuts.end());
  double sum = 0;
  for (std::size_t i = 1; i < cuts.size(); ++i) sum += 0.5 * (cuts[i] - cuts[i - 1]) * (at(cuts[i - 1]) + at(cuts[i]));
  return sum / (b - a);
}

double Tablier::rep_extreme(Id study_bar, double lo, double hi, double half_width,
                            bool maximise) const {
  if (hi < lo) std::swap(lo, hi);
  std::vector<double> c{lo, hi};
  const RepTransBlock* blk = nullptr;
  const bool on_deck = std::find(bars.begin(), bars.end(), study_bar) != bars.end();
  if (!rep_trans.empty() && !bars.empty() && (!on_deck || rep_trans_uniform()))
    blk = rep_block_of(*this, bars.front());
  else if (!rep_trans.empty())
    blk = rep_block_of(*this, study_bar);
  if (blk && blk->k.size() > 1 && width() > 0) {
    const double half = width() / 2;
    const std::size_t n = blk->k.size();
    for (std::size_t i = 0; i < n; ++i) {
      const double yi = -half + 2 * half * static_cast<double>(i) / static_cast<double>(n - 1);
      for (double y : {yi - half_width, yi + half_width, yi})
        if (y > lo && y < hi) c.push_back(y);
    }
  }
  std::sort(c.begin(), c.end());
  // the vertex of each quadratic piece (3-point fit)
  const std::size_t m = c.size();
  for (std::size_t i = 0; i + 1 < m; ++i) {
    const double a = c[i], b = c[i + 1];
    if (b - a <= 1e-12) continue;
    const double fa = rep_coef(study_bar, a, half_width), fm = rep_coef(study_bar, 0.5 * (a + b), half_width),
                 fb = rep_coef(study_bar, b, half_width);
    const double den = fa - 2 * fm + fb;
    if (std::fabs(den) > 1e-15 * (std::fabs(fa) + std::fabs(fb) + 1)) {
      const double t = 0.5 * (fa - fb) / den + 0.5;  // vertex in [0, 1]
      if (t > 0 && t < 1) c.push_back(a + t * (b - a));
    }
  }
  bool first = true;
  double best = 0;
  for (double y : c) {
    const double v = rep_coef(study_bar, y, half_width);
    if (first || (maximise ? v > best : v < best)) best = v;
    first = false;
  }
  return best;
}

void Tablier::validate() const {
  require(!bars.empty(), "TABLIER: la liste BARRE est vide");
  require(pas_ligne > 0 && finite(pas_ligne), "TABLIER PAS_LIGNE doit etre > 0");
  require(ldif >= 0 && finite(ldif), "TABLIER LDIF doit etre >= 0");
  require(eurocode_class() >= 1 && eurocode_class() <= 6,
          "TABLIER CLASSE TRAFIC: classe Eurocode 1..6");
  require(fasc61_class() >= 1 && fasc61_class() <= 3,
          "TABLIER CLASSE: classe Fasc.61 1..3 (art. 3)");
  for (const TransZone& z : zones) {
    require(finite(z.width) && z.width >= 0, "ZONE_TRANS: largeur non finie ou negative");
    require(z.nvoies >= 0, "ZONE_TRANS VOIE: nombre de voies negatif");
  }
  require(rep_trans.empty() || !zones.empty(),
          "REP_TRANS ne peut etre defini sans ZONE_TRANS (p.130)");
  for (const RepTransBlock& b : rep_trans) {
    require(!b.k.empty(), "REP_TRANS KTRANS: liste vide");
    for (double v : b.k) require(finite(v), "REP_TRANS KTRANS: coefficient non fini");
  }
}

// ===========================================================================
// 2. TABLES ET DENSITÉS RÉGLEMENTAIRES (p.138-146)
// ===========================================================================

namespace {
// Table [classe][nombre de voies chargées 1..5], 0 = combinaison interdite.
constexpr double kA1[3][5] = {
    {1.00, 1.00, 0.90, 0.75, 0.70},  // classe 1
    {1.00, 0.90, 0.00, 0.00, 0.00},  // classe 2
    {0.90, 0.80, 0.00, 0.00, 0.00},  // classe 3
};
constexpr double kBc[3][5] = {
    {1.20, 1.10, 0.95, 0.80, 0.70},
    {1.00, 1.00, 0.00, 0.00, 0.00},
    {1.00, 0.80, 0.00, 0.00, 0.00},
};
constexpr double kBt[3][3] = {
    {1.00, 1.00, 0.00},
    {0.90, 0.90, 0.00},
    {0.00, 0.00, 0.00},
};
constexpr double kV0[3] = {3.50, 3.00, 2.75};

double table_lookup(const double* row, int n_cols, int n_lanes, const char* what) {
  require(n_lanes >= 1, std::string(what) + ": au moins une voie chargee");
  int i = std::min(n_lanes, n_cols) - 1;
  return row[i];
}
}  // namespace

double st1_a1(int cls, int n_lanes) {
  require(cls >= 1 && cls <= 3, "a1: classe 1..3 (Fasc.61 art. 3)");
  return table_lookup(kA1[cls - 1], 5, n_lanes, "a1");
}
double st1_v0(int cls) {
  require(cls >= 1 && cls <= 3, "v0: classe 1..3");
  return kV0[cls - 1];
}
double st1_bc(int cls, int n_lanes) {
  require(cls >= 1 && cls <= 3, "bc: classe 1..3");
  return table_lookup(kBc[cls - 1], 5, n_lanes, "bc");
}
double st1_bt(int cls, int n_lanes) {
  require(cls >= 1 && cls <= 3, "bt: classe 1..3");
  return table_lookup(kBt[cls - 1], 3, n_lanes, "bt");
}

double st1_al(double L) {
  require(L >= 0 && finite(L), "A(l): longueur chargee non finie ou negative");
  return 0.23 + 36.0 / (12.0 + L);
}
double st1_al_frein(double L, double S) {
  require(S >= 0 && finite(S), "A(l) FREIN: surface chargee non finie ou negative");
  return st1_al(L) / (20.0 + 0.0035 * S);
}
double st1_al_plancher(double L) {
  require(L >= 0 && finite(L), "A(l) PLANCHER: longueur chargee non finie");
  return std::max(0.0, 0.4 - 0.0002 * L);
}
double st1_al_plancher_frein(double L, double S) {
  require(S >= 0 && finite(S), "A(l) PLANCHER FREIN: surface chargee non finie");
  return st1_al_plancher(L) / (20.0 + 0.0035 * S);
}
double st1_trot() { return 0.150; }
double st1_al_pieton(double L) {
  require(L >= 0 && finite(L), "A(l) PIETON: longueur chargee non finie");
  return 0.2 + 15.0 / (50.0 + L);
}

double AlgDef::density(double L) const {
  require(!dens.empty(), "ALG: table DENS vide");
  if (dens.size() == 1) return dens[0].second;
  if (L <= dens.front().first) return dens.front().second;
  if (L >= dens.back().first) return dens.back().second;
  for (std::size_t i = 1; i < dens.size(); ++i) {
    if (L <= dens[i].first) {
      double x0 = dens[i - 1].first, x1 = dens[i].first;
      double y0 = dens[i - 1].second, y1 = dens[i].second;
      double t = (x1 > x0) ? (L - x0) / (x1 - x0) : 0.0;
      return y0 + t * (y1 - y0);
    }
  }
  return dens.back().second;
}

void AlgDef::validate() const {
  require(!dens.empty(), "ALG: la table DENS est obligatoire");
  double prev = -1e300;
  for (const auto& [l, d] : dens) {
    require(finite(l) && finite(d) && l >= 0 && d >= 0, "ALG DENS: couple non fini ou negatif");
    require(l > prev, "ALG DENS: longueurs chargees non strictement croissantes");
    prev = l;
  }
  require(larg_voie >= 0 && finite(larg_voie), "ALG LARG VOIE: valeur non finie ou negative");
  require(finite(dir), "ALG DIR: direction non finie");
  for (double c : coef_trans) require(finite(c), "ALG COEF TRANS: coefficient non fini");
}

EurocodeParams EurocodeParams::scaled_forces(double f) const {
  require(finite(f) && f > 0, "EurocodeParams: facteur de force invalide");
  EurocodeParams p = *this;
  for (double& v : p.Qk) v *= f;
  for (double& v : p.qk) v *= f;
  p.qrk *= f;
  p.force_scale *= f;
  return p;
}

EurocodeParams EurocodeParams::traffic_class(int cls) {
  EurocodeParams p;
  switch (cls) {
    // p.164 : αQ1, αQi, αq1, αqi, αqr, QLM2 (= beta).
    case 1: p.aQ1 = 1.0;  p.aQi = 1.0;  p.aq1 = 1.00; p.aqi = 1.20; p.aqr = 1.20; p.beta_lm2 = 0.8; break;
    case 2: p.aQ1 = 0.9;  p.aQi = 0.8;  p.aq1 = 0.70; p.aqi = 1.00; p.aqr = 1.00; p.beta_lm2 = 0.8; break;
    case 4: p.aQ1 = 0.32; p.aQi = 0.48; p.aq1 = 0.58; p.aqi = 2.09; p.aqr = 2.09; p.beta_lm2 = 0.56; break;
    case 5: p.aQ1 = 0.20; p.aQi = 0.30; p.aq1 = 0.36; p.aqi = 1.30; p.aqr = 1.30; p.beta_lm2 = 0.36; break;
    case 6: p.aQ1 = 0.06; p.aQi = 0.09; p.aq1 = 0.21; p.aqi = 0.76; p.aqr = 0.76; p.beta_lm2 = 0.15; break;
    // p.168 : la classe 3 est tabulée (0,8 ; 0,5 ; 0,5 ; 1 ; 1 ; 0,8)
    case 3: p.aQ1 = 0.8;  p.aQi = 0.5;  p.aq1 = 0.50; p.aqi = 1.00; p.aqr = 1.00; p.beta_lm2 = 0.8; break;
    default: fail("EUROCODE CLASSE TRAFIC: classe 1 a 6");
  }
  return p;
}

// ===========================================================================
// 3. CATALOGUE DE CONVOIS (p.145-155)
// ===========================================================================

double ConvoyDef::total_t() const {
  double t = 0;
  for (const ConvoyAxle& a : axles) t += a.weight_t;
  return t;
}

ConvoyDef ConvoyDef::mirrored() const {
  ConvoyDef m = *this;
  m.name += "/SENS2";
  for (ConvoyAxle& a : m.axles) a.xl = longueur - a.xl;
  std::reverse(m.axles.begin(), m.axles.end());
  return m;
}

void ConvoyDef::validate() const {
  require(!axles.empty(), "CONV: aucun essieu defini (ESSIEU)");
  require(max_cam >= 1, "CONV MAX_CAM: au moins un camion");
  require(max_file >= 0, "CONV MAX_FILE: valeur negative");
  require(finite(longueur) && longueur > 0, "CONV LONG: encombrement longitudinal non fini");
  require(finite(largeur) && largeur >= 0, "CONV LARG: encombrement transversal non fini");
  for (const ConvoyAxle& a : axles) {
    require(finite(a.xl) && finite(a.weight_t), "CONV ESSIEU: XL/POIDS non fini");
    require(a.weight_t >= 0, "CONV ESSIEU POIDS: le poids doit etre positif (p.160)");
    require(a.impact_l >= 0 && a.impact_t >= 0 && finite(a.impact_l) && finite(a.impact_t),
            "CONV ESSIEU IMPACT: dimension non finie ou negative");
    for (double y : a.yl) require(finite(y), "CONV ESSIEU YL: position non finie");
  }
  for (double c : coef_trans) require(finite(c), "CONV COEF TRANS: coefficient non fini");
}

namespace {

ConvoyAxle mk_axle(double xl, double w, double il = 0, double it = 0,
                   std::vector<double> yl = {}) {
  ConvoyAxle a;
  a.xl = xl;
  a.weight_t = w;
  a.impact_l = il;
  a.impact_t = it;
  a.yl = std::move(yl);
  return a;
}

// Train d'essieux régulier : n essieux de poids w, pas constant.
void push_regular(ConvoyDef& c, int n, double w, double pitch, double x0, double il, double it,
                  const std::vector<double>& yl) {
  for (int i = 0; i < n; ++i)
    c.axles.push_back(mk_axle(x0 + i * pitch, w, il, it, yl));
}

std::string upper(const std::string& s) {
  std::string o;
  for (char ch : s)
    if (ch != ' ' && ch != '\t') o.push_back(static_cast<char>(std::toupper(ch)));
  return o;
}

// --- catalogue ---------------------------------------------------------
//
// Chaque convoi est cité avec la page du manuel ST1 v24 dont il est tiré et
// le niveau de CONFIANCE de la reconstruction, l'OCR des figures du manuel
// perdant la mise en page des cotes :
//   HAUTE   — le manuel donne l'équivalent CONV complet (Bc p.162, Mc120 p.162)
//             ou les cotes sont sans ambiguïté et bouclent arithmétiquement ;
//   MOYENNE — les cotes bouclent (somme des entraxes = longueur totale, somme
//             des poids = tonnage annoncé) mais l'ORDRE des essieux est déduit ;
//   BASSE   — seul le tonnage total est certain (il est dans le mot-clé), le
//             schéma d'entraxes est reconstruit et signalé comme tel.
ConvoyDef build_convoy(const std::string& n) {
  ConvoyDef c;
  c.name = n;
  const std::vector<double> w2 = {-1.0, 1.0};  // roues à +/- 1,00 m (voie de 2,00 m)

  if (n == "BC") {
    // p.145 + équivalent CONV p.162 (CONFIANCE HAUTE) : deux camions de 30 t
    // par file, essieux 6 t / 12 t / 12 t aux abscisses 2,25 / 6,75 / 8,25 m
    // depuis l'avant, encombrement 10,50 x 2,50 m, roues à +/- 1,00 m.
    c.source = "ST1 p.145 + CONV equivalent p.162; Fasc.61 II art. 5.2 (confiance HAUTE)";
    c.axles = {mk_axle(2.25, 6, 0.20, 0.20, w2), mk_axle(6.75, 12, 0.25, 0.25, w2),
               mk_axle(8.25, 12, 0.25, 0.25, w2)};
    c.longueur = 10.50;
    c.largeur = 2.50;
    c.max_cam = 2;
    c.trans = TransLaw::Fasc61Bc;
  } else if (n == "BC_FREIN") {
    // p.145 : « un groupe de trois forces dirigees suivant l'axe ox local » —
    // un seul camion Bc, forces longitudinales (CONFIANCE HAUTE).
    c.source = "ST1 p.145; Fasc.61 II art. 6 (confiance HAUTE)";
    c.axles = {mk_axle(2.25, 6, 0.20, 0.20, w2), mk_axle(6.75, 12, 0.25, 0.25, w2),
               mk_axle(8.25, 12, 0.25, 0.25, w2)};
    c.longueur = 10.50;
    c.largeur = 2.50;
    c.max_cam = 1;
    c.max_file = 1;
    c.longitudinal = true;
    // the braking of ONE truck, without the bc coefficient (Fasc.61 II
    // art. 6; ST1 2.24: Fx = 30 t on a class-1 deck)
    c.trans = TransLaw::None;
  } else if (n == "BT") {
    // p.146 : tandem de 2 essieux de 16 t a 1,35 m, roues a +/- 1,00 m,
    // impact 0,25 x 0,60 m ; au plus 2 tandems (table bt) (CONFIANCE HAUTE).
    c.source = "ST1 p.146; Fasc.61 II art. 5.4 (confiance HAUTE)";
    c.axles = {mk_axle(0.675, 16, 0.25, 0.60, w2), mk_axle(2.025, 16, 0.25, 0.60, w2)};
    c.longueur = 2.70;
    c.largeur = 3.00;  // p.150 transverse schema: two tandems of 3,00 m
    c.max_cam = 1;
    c.max_file = 2;
    c.trans = TransLaw::Fasc61Bt;
  } else if (n == "MC80") {
    // p.146 : char a 2 chenilles de 4,90 m, 72 t, largeur hors-tout 3,65 m,
    // chenille de 0,85 m => entraxe 2,80 m ; « dix forces separees au minimum
    // par 30,50 m » = au plus 10 chars par file, pas mini 4,90 + 30,50 m.
    // Les deux chenilles etant identiques longitudinalement, l'empreinte est
    // portee par un unique « essieu » reparti de 6,10/4,90 m (CONFIANCE MOYENNE).
    c.source = "ST1 p.146; Fasc.61 II art. 9.21 (confiance MOYENNE)";
    c.axles = {mk_axle(2.45, 72, 4.90, 0.85, {-1.40, 1.40})};
    c.longueur = 4.90 + 30.50;
    c.largeur = 3.65;
    c.max_cam = 10;
    c.max_file = 1;
  } else if (n == "MC120") {
    // p.147 + equivalent CONV p.162 (CONFIANCE HAUTE) : 110 t sur deux
    // chenilles de 6,10 x 1,00 m, entraxe 3,30 m, largeur 4,30 m,
    // LONG 36,60 m (= 6,10 + 30,50), MAX_CAM 10, MAX_FILE 1.
    c.source = "ST1 p.147 + CONV equivalent p.162; Fasc.61 II art. 9.21 (confiance HAUTE)";
    c.axles = {mk_axle(3.05, 110, 6.10, 1.00, {-1.65, 1.65})};
    c.longueur = 36.60;
    c.largeur = 4.30;
    c.max_cam = 10;
    c.max_file = 1;
  } else if (n == "ME80") {
    // p.146 : deux essieux de 22 t a 1,50 m, impact long. 0,12 m, largeur
    // hors-tout 3,50 m (CONFIANCE MOYENNE).
    c.source = "ST1 p.146; Fasc.61 II art. 9.21 (confiance MOYENNE)";
    c.axles = {mk_axle(0.06, 22, 0.12, 3.50, {0.0}), mk_axle(1.56, 22, 0.12, 3.50, {0.0})};
    c.longueur = 1.62;
    c.largeur = 3.50;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "ME120") {
    // p.147 : deux essieux de 33 t a 1,80 m, impact long. 0,15 m, largeur
    // hors-tout 4,00 m (CONFIANCE MOYENNE).
    c.source = "ST1 p.147; Fasc.61 II art. 9.21 (confiance MOYENNE)";
    c.axles = {mk_axle(0.075, 33, 0.15, 4.00, {0.0}), mk_axle(1.875, 33, 0.15, 4.00, {0.0})};
    c.longueur = 1.95;
    c.largeur = 4.00;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_F61_D") {
    // p.147 : convoi D, 280 t = 2 x 140 t repartis sur 11 m, entre-deux 8 m,
    // largeur 3,30 m ; « circule seul sans autre ponderation » (p.143).
    c.source = "ST1 p.143/p.147; Fasc.61 II convoi D (confiance HAUTE)";
    c.axles = {mk_axle(5.5, 140, 11.0, 3.30, {0.0}), mk_axle(24.5, 140, 11.0, 3.30, {0.0})};
    c.longueur = 30.0;
    c.largeur = 3.30;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_F61_E") {
    // p.147 : convoi E, 400 t = 2 x 200 t sur 15 m, entre-deux 18 m, l. 3,30 m.
    c.source = "ST1 p.143/p.147; Fasc.61 II convoi E (confiance HAUTE)";
    c.axles = {mk_axle(7.5, 200, 15.0, 3.30, {0.0}), mk_axle(40.5, 200, 15.0, 3.30, {0.0})};
    c.longueur = 48.0;
    c.largeur = 3.30;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_1946_3CL") {
    // p.154 : camion de 3e classe 1946, 70 t = 2 empreintes de 35 t sur 1,00 m
    // separees de 1,60 m ; 10 camions par file, entraxe 40 m (CONFIANCE MOYENNE
    // sur les cotes transversales, HAUTE sur le tonnage et l'entraxe).
    c.source = "ST1 p.143/p.154; Fasc.61 de 1946, 3e classe (confiance MOYENNE)";
    c.axles = {mk_axle(0.5, 35, 1.0, 2.50, {0.0}), mk_axle(3.1, 35, 1.0, 2.50, {0.0})};
    c.longueur = 40.0;
    c.largeur = 3.50;
    c.max_cam = 10;
    c.max_file = 1;
  } else if (n == "VS_1946_4CL") {
    // p.154 : camion de 4e classe 1946, 100 t = 2 x 50 t sur 1,00 m separees
    // de 1,80 m ; 10 camions par file, entraxe 60 m.
    c.source = "ST1 p.143/p.154; Fasc.61 de 1946, 4e classe (confiance MOYENNE)";
    c.axles = {mk_axle(0.5, 50, 1.0, 2.50, {0.0}), mk_axle(3.3, 50, 1.0, 2.50, {0.0})};
    c.longueur = 60.0;
    c.largeur = 4.50;
    c.max_cam = 10;
    c.max_file = 1;
  } else if (n == "CV_C1") {
    // Guide CEREMA 2022 annexe 3 fig. A3.2 (Setra 1982), modelisation
    // simplifiee : 6 t, 12 t a 2,75 m, 12 t a 1,35 m, puis 64 t uniformement
    // repartis sur 7,75 m a 4,225 m du dernier essieu (2,75 + 1,35 + 4,225 +
    // 7,75 = 16,075 m), largeur 3,20 m.
    c.source = "ST1 p.148; guide CEREMA 2022 fig. A3.2 (confiance HAUTE)";
    c.axles = {mk_axle(0.0, 6, 0, 3.20, {0.0}), mk_axle(2.75, 12, 0, 3.20, {0.0}),
               mk_axle(4.10, 12, 0, 3.20, {0.0}), mk_axle(12.20, 64.0, 7.75, 3.20, {0.0})};
    c.longueur = 16.075;
    c.largeur = 3.20;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_C2") {
    // fig. A3.3 : 7 t, 12,75 t a 2,75 m, 12,75 t a 1,35 m, puis 87,5 t sur
    // 7,75 m a 2,725 m (2,75 + 1,35 + 2,725 + 7,75 = 14,575 m), l. 3,20 m.
    c.source = "ST1 p.148; guide CEREMA 2022 fig. A3.3 (confiance HAUTE)";
    c.axles = {mk_axle(0.0, 7, 0, 3.20, {0.0}), mk_axle(2.75, 12.75, 0, 3.20, {0.0}),
               mk_axle(4.10, 12.75, 0, 3.20, {0.0}), mk_axle(10.70, 87.5, 7.75, 3.20, {0.0})};
    c.longueur = 14.575;
    c.largeur = 3.20;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_D2F1") {
    c.source = "ST1 p.149; R/EG.3 1983 convoi D2F1 (confiance HAUTE)";
    c.axles = {mk_axle(10.85, 245, 21.70, 3.20, {0.0})};
    c.longueur = 21.70;
    c.largeur = 3.20;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_D3F1") {
    c.source = "ST1 p.149; R/EG.3 1983 convoi D3F1 (confiance HAUTE)";
    c.axles = {mk_axle(8.525, 250, 17.05, 5.15, {0.0})};
    c.longueur = 17.05;
    c.largeur = 5.15;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_D3F2") {
    // p.149 : 2 x 125 t sur 9,30 m, entre-deux 14,70 m (9,30+14,70+9,30 = 33,30).
    c.source = "ST1 p.149; R/EG.3 1983 convoi D3F2 (confiance HAUTE)";
    c.axles = {mk_axle(4.65, 125, 9.30, 5.15, {0.0}), mk_axle(28.65, 125, 9.30, 5.15, {0.0})};
    c.longueur = 33.30;
    c.largeur = 5.15;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_E2F1") {
    c.source = "ST1 p.150; R/EG.3 1983 convoi E2F1 (confiance HAUTE)";
    c.axles = {mk_axle(15.50, 350, 31.00, 3.20, {0.0})};
    c.longueur = 31.00;
    c.largeur = 3.20;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_E3F1") {
    c.source = "ST1 p.150; R/EG.3 1983 convoi E3F1 (confiance HAUTE)";
    c.axles = {mk_axle(13.175, 400, 26.35, 5.15, {0.0})};
    c.longueur = 26.35;
    c.largeur = 5.15;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "CV_E3F2") {
    // p.150 : 2 x 188 t sur 12,40 m, entre-deux 11,60 m (12,40+11,60+12,40 = 36,40).
    c.source = "ST1 p.150; R/EG.3 1983 convoi E3F2 (confiance HAUTE)";
    c.axles = {mk_axle(6.20, 188, 12.40, 5.15, {0.0}), mk_axle(30.20, 188, 12.40, 5.15, {0.0})};
    c.longueur = 36.40;
    c.largeur = 5.15;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_72_OAN" || n == "VS_94_OAN" || n == "VS_120_OAN") {
    // Guide CEREMA 2022 « Transports exceptionnels », fig. 13, 18 et 23
    // (ouvrages a construire) : essieux de 12 t au pas de 1,35 m, le 94 t
    // commencant par un essieu de 10 t ; impact 0,40 x 0,40 m, entraxe
    // transversal des roues 2,25 m, encombrement 2,65 m.
    c.source = "ST1 p.151-153; guide CEREMA 2022 fig. 13/18/23 (confiance HAUTE)";
    const std::vector<double> wr = {-1.125, 1.125};
    if (n == "VS_72_OAN") {
      push_regular(c, 6, 12, 1.35, 0.0, 0.40, 0.40, wr);
    } else if (n == "VS_94_OAN") {
      c.axles.push_back(mk_axle(0.0, 10, 0.40, 0.40, wr));
      push_regular(c, 7, 12, 1.35, 1.35, 0.40, 0.40, wr);
    } else {
      push_regular(c, 10, 12, 1.35, 0.0, 0.40, 0.40, wr);
    }
    c.longueur = c.axles.back().xl + 0.40;
    c.largeur = 2.65;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_72_OAE" || n == "VS_94_OAE" || n == "VS_120_OAE") {
    // Guide CEREMA 2022, fig. 14, 19 et 24 (ouvrages existants) : semi-remorque,
    // essieu avant, deux essieux moteurs a 1,35 m, puis la remorque au pas de
    // 1,35 m. 72 t : 6 / 9 9 / 4 x 12 (3,20 1,35 4,50) ; 94 t : 10 / 12 12 /
    // 5 x 12 (3,20 1,35 5,05) ; 120 t : 12 / 12 12 / 7 x 12 (3,20 1,35 5,07).
    // Entraxe transversal des roues 2,00 m (modelisation simplifiee du guide),
    // encombrement 3,20 m, porte-a-faux 1,30 m avant et 2,50 m arriere.
    c.source = "ST1 p.151-153; guide CEREMA 2022 fig. 14/19/24 (confiance HAUTE)";
    const std::vector<double> wr = {-1.0, 1.0};
    double front = 6, drive = 9, gap = 4.50;
    int trailer = 4;
    if (n == "VS_94_OAE") {
      front = 10, drive = 12, gap = 5.05, trailer = 5;
    } else if (n == "VS_120_OAE") {
      front = 12, drive = 12, gap = 5.07, trailer = 7;
    }
    c.axles.push_back(mk_axle(0.0, front, 0.40, 0.40, wr));
    c.axles.push_back(mk_axle(3.20, drive, 0.40, 0.40, wr));
    c.axles.push_back(mk_axle(4.55, drive, 0.40, 0.40, wr));
    push_regular(c, trailer, 12, 1.35, 4.55 + gap, 0.40, 0.40, wr);
    c.longueur = c.axles.back().xl + 0.40;
    c.largeur = 3.20;
    c.max_cam = 1;
    c.max_file = 1;
  } else if (n == "VS_48_GRUE" || n == "VS_72_GRUE" || n == "VS_94_GRUE" ||
             n == "VS_108_GRUE") {
    // Guide CEREMA 2022, fig. 9, 14, 19 et 24 : grues automotrices, impact
    // 0,40 x 0,40 m. 48 t : 4 x 12 (1,65 2,03 1,65), encombrement 2,75 m,
    // roues a 2,35 m ; 72 t : 6 x 12 (1,75 2,20 1,75 2,20 1,75) ; 94 t :
    // 8 x 11,75 (1,60 2,80 1,60 1,60 2,80 1,60 1,60) ; 108 t : 9 x 12
    // (1,60 1,60 2,80 1,60 1,60 2,80 1,60 1,60) ; encombrement 3,00 m, roues
    // a 2,55 m.
    c.source = "ST1 p.150-153; guide CEREMA 2022 fig. 9/14/19/24 (confiance HAUTE)";
    std::vector<double> wr = {-1.275, 1.275};
    std::vector<double> gaps;
    double w = 12;
    if (n == "VS_48_GRUE") {
      gaps = {1.65, 2.03, 1.65};
      wr = {-1.175, 1.175};
      c.largeur = 2.75;
    } else if (n == "VS_72_GRUE") {
      gaps = {1.75, 2.20, 1.75, 2.20, 1.75};
      c.largeur = 3.00;
    } else if (n == "VS_94_GRUE") {
      gaps = {1.60, 2.80, 1.60, 1.60, 2.80, 1.60, 1.60};
      w = 11.75;
      c.largeur = 3.00;
    } else {
      gaps = {1.60, 1.60, 2.80, 1.60, 1.60, 2.80, 1.60, 1.60};
      c.largeur = 3.00;
    }
    double x = 0;
    c.axles.push_back(mk_axle(x, w, 0.40, 0.40, wr));
    for (double g : gaps) {
      x += g;
      c.axles.push_back(mk_axle(x, w, 0.40, 0.40, wr));
    }
    c.longueur = x + 0.40;
    c.max_cam = 1;
    c.max_file = 1;
  } else {
    fail("SURCH: convoi inconnu '" + n + "' (voir st1_convoy_names())");
  }
  c.validate();
  return c;
}

}  // namespace

std::vector<std::string> st1_convoy_names() {
  return {"BC",         "BC_FREIN",   "BT",          "MC80",       "ME80",
          "MC120",      "ME120",      "VS_F61_D",    "VS_F61_E",   "VS_1946_3CL",
          "VS_1946_4CL","CV_C1",      "CV_C2",       "CV_D2F1",    "CV_D3F1",
          "CV_D3F2",    "CV_E2F1",    "CV_E3F1",     "CV_E3F2",    "VS_48_GRUE",
          "VS_72_OAN",  "VS_72_GRUE", "VS_72_OAE",   "VS_94_OAN",  "VS_94_GRUE",
          "VS_94_OAE",  "VS_120_OAN", "VS_108_GRUE", "VS_120_OAE"};
}

ConvoyDef st1_convoy(const std::string& name) {
  ConvoyDef c = build_convoy(upper(name));
  if (c.name.rfind("CV_", 0) == 0) c.axis_ecc = 0.30;  // R/EG.3 : dans l'axe de la zone
  return c;
}

ConvoyPattern convoy_pattern(const ConvoyDef& c, double tonne_force, double ldif) {
  require(finite(tonne_force) && tonne_force > 0, "convoi: facteur tonne->force invalide");
  require(ldif >= 0 && finite(ldif), "convoi: LDIF invalide");
  ConvoyPattern p;
  for (const ConvoyAxle& a : c.axles) {
    double P = a.weight_t * tonne_force;
    double Le = a.impact_l + 2 * ldif;  // diffusion longitudinale (p.126)
    if (Le > 1e-9)
      p.segments.push_back({a.xl - Le / 2, Le, P});
    else
      p.axles.emplace_back(a.xl, P);
  }
  return p;
}

// ===========================================================================
// 4. GÉOMÉTRIE DU CHEMIN : abscisse curviligne <-> (barre, abscisse locale)
// ===========================================================================

namespace {

struct DeckSpan {
  Id bar = 0;
  double s0 = 0, L = 0;
  bool flip = false;  // la barre est parcourue de j vers i
  double x_of(double s) const { return flip ? L - (s - s0) : (s - s0); }
};

// Même logique de chaînage que build_path (influence.cpp) : l'orientation de la
// première barre est fixée par la seconde, puis on suit le nœud de tête.
std::vector<DeckSpan> deck_spans(const Model& m, const std::vector<Id>& bars) {
  require(!bars.empty(), "TABLIER: liste de barres vide");
  std::vector<DeckSpan> out;
  double s = 0;
  Id head = 0;
  for (std::size_t i = 0; i < bars.size(); ++i) {
    const Bar& b = m.bar(bars[i]);
    DeckSpan sp;
    sp.bar = b.id;
    sp.L = m.bar_chord_length(b);
    sp.s0 = s;
    if (i == 0) {
      if (bars.size() > 1) {
        const Bar& b2 = m.bar(bars[1]);
        if (b.nj == b2.ni || b.nj == b2.nj) {
          sp.flip = false;
          head = b.nj;
        } else if (b.ni == b2.ni || b.ni == b2.nj) {
          sp.flip = true;
          head = b.ni;
        } else {
          fail("TABLIER: les barres " + std::to_string(b.id) + " et " +
               std::to_string(b2.id) + " ne sont pas contigues");
        }
      } else {
        head = b.nj;
      }
    } else {
      if (b.ni == head) {
        sp.flip = false;
        head = b.nj;
      } else if (b.nj == head) {
        sp.flip = true;
        head = b.ni;
      } else {
        fail("TABLIER: la barre " + std::to_string(b.id) + " ne prolonge pas la chaine");
      }
    }
    s += sp.L;
    out.push_back(sp);
  }
  return out;
}

const DeckSpan* span_at(const std::vector<DeckSpan>& sp, double s) {
  for (const DeckSpan& d : sp)
    if (s >= d.s0 - 1e-9 && s <= d.s0 + d.L + 1e-9) return &d;
  return nullptr;
}

// Tangente unitaire (repère global) du chemin à l'abscisse s.
Vec3 path_tangent(const Model& m, const DeckSpan& d) {
  const Bar& b = m.bar(d.bar);
  Vec3 t = m.node(b.nj).x - m.node(b.ni).x;
  double L = t.norm();
  require(L > 0, "TABLIER: barre degeneree sur le chemin de roulement");
  t /= L;
  return d.flip ? Vec3(-t) : t;
}

// Couple unitaire de l'excentrement à la travée `d` : une charge P d à
// l'abscisse transversale y (positive à gauche du parcours) vaut une force
// P d sur l'axe plus le couple P y (n_gauche x d).
Vec3 ecc_torque_dir(const Model& m, const DeckSpan& d, const Vec3& dir) {
  Vec3 t = path_tangent(m, d);
  Vec3 nl = (-dir).cross(t);
  double n = nl.norm();
  if (n < 1e-12) return Vec3::Zero();
  return (nl / n).cross(dir);
}

// Charge ponctuelle à l'abscisse s.
void add_point(const Model& m, const std::vector<DeckSpan>& sp, double s, double P,
               const Vec3& dir, bool local, std::vector<Load>& out) {
  const DeckSpan* d = span_at(sp, s);
  if (d == nullptr || P == 0) return;
  double x = std::clamp(d->x_of(s), 0.0, d->L);
  if (local)
    out.push_back(Load::point(d->bar, x, Vec3(d->flip ? -P : P, 0, 0), Vec3::Zero(),
                              LoadFrame::Local));
  else
    out.push_back(Load::point(d->bar, x, P * dir));
  (void)m;
}

// Charge uniforme de densité q [N/m] sur [s0, s1], découpée aux joints de barre.
void add_udl(const Model& m, const std::vector<DeckSpan>& sp, double s0, double s1, double q,
             const Vec3& dir, bool local, std::vector<Load>& out) {
  if (s1 <= s0 || q == 0) return;
  for (const DeckSpan& d : sp) {
    double a = std::max(s0, d.s0), b = std::min(s1, d.s0 + d.L);
    if (b <= a + 1e-12) continue;
    double xa = d.x_of(a), xb = d.x_of(b);
    if (xa > xb) std::swap(xa, xb);
    xa = std::clamp(xa, 0.0, d.L);
    xb = std::clamp(xb, 0.0, d.L);
    if (xb <= xa + 1e-12) continue;
    if (local)
      out.push_back(Load::udl(d.bar, Vec3(d.flip ? -q : q, 0, 0), LoadFrame::Local, xa, xb));
    else
      out.push_back(Load::udl(d.bar, q * dir, LoadFrame::Global, xa, xb));
  }
  (void)m;
}

}  // namespace

// ===========================================================================
// 5. LIGNES D'INFLUENCE GROUPÉES
// ===========================================================================
//
// influence_line() (influence.h) assemble et factorise à CHAQUE appel et
// n'extrait qu'une réponse. Une étude SURCH demande des dizaines de réponses
// sur le MÊME modèle : on factorise donc une fois et on extrait toutes les
// réponses de chaque colonne du balayage. L'extraction scalaire est celle de
// influence.cpp, partagée telle quelle.

namespace {

// L'extraction scalaire elle-même est PARTAGÉE avec la ligne d'influence
// (`extract_response`, influence.h) : une seule définition, un seul jeu de
// conventions de signe pour les deux voies.

// Balayage unitaire groupé : une factorisation, `dirs.size()` balayages.
// Chaque réponse est AFFINE en (déplacements, charge) : valeur = g.u + h(charge),
// g creux (au plus les DDL que l'extracteur lit, `response_dofs`) obtenu une
// fois par sondes unitaires, h non nul seulement quand la charge unitaire est
// posée sur la barre étudiée (effort / déplacement de section) ou sur une
// barre incidente au nœud (réaction). Une ordonnée coûte donc O(|g|) au lieu
// d'une restitution complète : indispensable pour étudier toutes les
// composantes de toutes les sections (édition ST1, p.137).
// Lignes d'influence par Müller-Breslau / Maxwell-Betti. La ligne d'une
// réponse R est la déformée de la structure sous la sollicitation duale de R :
// le long d'une barre de portique prismatique (Euler-Bernoulli ou Timoshenko,
// libérations et excentrements compris) c'est un polynôme de degré <= 3 de
// l'abscisse de la charge — la partie globale (déplacements adjoints
// interpolés par les fonctions de forme exactes) partout, plus, sur la barre
// de la section étudiée, le terme local à extrémités bloquées, cubique de
// part et d'autre de la section (saut ou coude du théorème). Chaque ligne est
// donc EXACTEMENT une cubique d'Hermite par morceau de barre du tracé (coupé à
// la section de la réponse), obtenue de 4 ordonnées intérieures par morceau :
// quelques sondes par barre au lieu d'un échantillonnage fin. `ok` faux (et
// rien calculé) quand une barre du tracé sort de ce cadre (CARA VAR, sol,
// treillis, câble) ou quand le champ de directions varie dans une barre.
std::vector<InfluenceLine> il_batch_exact(const Model& model, const LoadPath& path,
                                          const std::vector<Vec3>& station_dir,
                                          const std::vector<ResponseKind>& responses,
                                          const std::vector<Vec3>* torque_dir,
                                          std::vector<InfluenceLine>* out_torque, bool& ok) {
  ok = false;
  const std::size_t ns = path.stations.size();
  if (ns < 2 || station_dir.size() != ns) return {};
  // the runs: consecutive stations on one bar
  struct Run {
    Id bar;
    std::size_t k0, k1;
    double s0, s1, x0, x1;
  };
  std::vector<Run> runs;
  for (std::size_t k = 0; k < ns;) {
    std::size_t j = k;
    while (j + 1 < ns && path.stations[j + 1].bar == path.stations[k].bar) ++j;
    const PathStation &a = path.stations[k], &b = path.stations[j];
    if (b.s > a.s) runs.push_back({a.bar, k, j, a.s, b.s, a.x_local, b.x_local});
    k = j + 1;
  }
  if (runs.empty()) return {};
  for (std::size_t r = 0; r + 1 < runs.size(); ++r)
    if (std::fabs(runs[r + 1].s0 - runs[r].s1) > 1e-9 * std::max(1.0, path.length)) return {};
  for (const Run& r : runs) {
    const Bar& b = model.bar(r.bar);
    if (b.kind != BarKind::Frame || model.tapered_of(b.id) || b.soil.any()) return {};
    for (std::size_t k = r.k0; k <= r.k1; ++k) {
      if ((station_dir[k] - station_dir[r.k0]).lpNorm<Eigen::Infinity>() > 1e-12) return {};
      if (torque_dir && ((*torque_dir)[k] - (*torque_dir)[r.k0]).lpNorm<Eigen::Infinity>() > 1e-12) return {};
    }
  }
  // runs joined at a node: same load direction(s) on both sides
  std::vector<char> run_join(runs.size(), 0);
  for (std::size_t r = 1; r < runs.size(); ++r) {
    const std::size_t a = runs[r - 1].k1, b = runs[r].k0;
    bool same = (station_dir[a] - station_dir[b]).lpNorm<Eigen::Infinity>() <= 1e-12;
    if (torque_dir) same = same && ((*torque_dir)[a] - (*torque_dir)[b]).lpNorm<Eigen::Infinity>() <= 1e-12;
    run_join[r] = same;
  }
  model.validate();
  Assembled asm_ = assemble_elements(model);
  SpMat K = assemble_stiffness(model, asm_);
  FactorCache fac;
  fac.factorize(K);
  const DofMap& dm = asm_.dofs;
  const CaseLoads empty = empty_case_loads(asm_);
  std::vector<std::vector<std::pair<int, double>>> grad(responses.size());
  {
    VecX probe = VecX::Zero(dm.n_free);
    for (std::size_t j = 0; j < responses.size(); ++j)
      for (int f : response_dofs(model, asm_, responses[j])) {
        int fi = dm.free_ix[static_cast<std::size_t>(f)];
        if (fi < 0) continue;
        probe(fi) = 1;
        double gv = extract_response(model, asm_, empty, probe, responses[j]);
        probe(fi) = 0;
        if (gv != 0.0) grad[j].emplace_back(fi, gv);
      }
  }
  // the responses whose ordinate has a local term when the load is on bar b
  // node -> the bars ending on it, in model order (each bar once), built on
  // the first reaction response: a scan of every bar per reaction was O(R x B)
  std::unordered_map<Id, std::vector<Id>> incident;
  bool incident_built = false;
  auto bars_at_node = [&](Id node) -> const std::vector<Id>& {
    if (!incident_built) {
      for (const Bar& b : model.bars) {
        incident[b.ni].push_back(b.id);
        if (b.nj != b.ni) incident[b.nj].push_back(b.id);
      }
      incident_built = true;
    }
    static const std::vector<Id> none;
    auto it = incident.find(node);
    return it == incident.end() ? none : it->second;
  };
  std::map<Id, std::vector<std::size_t>> local_of_bar;
  for (std::size_t j = 0; j < responses.size(); ++j) {
    const ResponseKind& rk = responses[j];
    switch (rk.what) {
      case ResponseKind::What::Displacement: break;
      case ResponseKind::What::BarForce:
      case ResponseKind::What::BarDisplacement: local_of_bar[rk.bar].push_back(j); break;
      case ResponseKind::What::Reaction:
        for (Id bid : bars_at_node(rk.node)) local_of_bar[bid].push_back(j);
        break;
    }
  }
  const VecX zero = VecX::Zero(dm.n_free);
  // 4 interior points of a piece (any 4 distinct points fix a cubic; inside,
  // a load never sits on the section or a bar end)
  static const double tq[4] = {0.125, 0.375, 0.625, 0.875};
  // Lagrange basis at tq: coefficients of the cubic through 4 values
  Eigen::Matrix4d Vinv;
  {
    Eigen::Matrix4d V;
    for (int r = 0; r < 4; ++r)
      for (int c = 0; c < 4; ++c) V(r, c) = std::pow(tq[r], c);
    Vinv = V.inverse();
  }
  auto fit = [&](const double v[4]) {
    Eigen::Vector4d y(v[0], v[1], v[2], v[3]);
    return Eigen::Vector4d(Vinv * y);  // p(t) = c0 + c1 t + c2 t^2 + c3 t^3
  };
  auto peval = [](const Eigen::Vector4d& c, double t) { return c(0) + t * (c(1) + t * (c(2) + t * c(3))); };
  auto pder = [](const Eigen::Vector4d& c, double t) { return c(1) + t * (2 * c(2) + t * 3 * c(3)); };

  auto sweep = [&](const std::vector<Vec3>& field, bool moment) {
    // the unit action at (bar, x) in direction d, lowered
    auto unit_case = [&](Id bar, double x, const Vec3& d) {
      LoadCase lc;
      lc.id = 1;
      lc.name = "unit";
      lc.loads.push_back(moment ? Load::point(bar, x, Vec3::Zero(), d) : Load::point(bar, x, d));
      return lower_case(model, asm_, lc);
    };
    const std::size_t nr = runs.size(), np = 4 * nr;
    // global part at the shared probes: one multi-RHS solve
    MatX rhs(dm.n_free, static_cast<Eigen::Index>(np));
    std::vector<CaseLoads> pcl;
    pcl.reserve(np);
    for (std::size_t r = 0; r < nr; ++r)
      for (int q = 0; q < 4; ++q) {
        const Run& R = runs[r];
        const double x = R.x0 + tq[q] * (R.x1 - R.x0);
        pcl.push_back(unit_case(R.bar, x, field[R.k0]));
        rhs.col(static_cast<Eigen::Index>(4 * r + static_cast<std::size_t>(q))) = pcl.back().rhs;
      }
    const MatX U = fac.solve(rhs);
    const MatX UT = U.transpose();
    // the local term at sample points of a split piece: per (bar, abscissa)
    std::map<std::pair<Id, double>, CaseLoads> split_cl;
    auto split_case = [&](Id bar, double x, const Vec3& d) -> const CaseLoads& {
      auto it = split_cl.find({bar, x});
      if (it == split_cl.end()) it = split_cl.emplace(std::make_pair(bar, x), unit_case(bar, x, d)).first;
      return it->second;
    };
    std::vector<InfluenceLine> res(responses.size());
    std::vector<double> g(np);
    for (std::size_t j = 0; j < responses.size(); ++j) {
      std::fill(g.begin(), g.end(), 0.0);
      for (const auto& [fi, gv] : grad[j]) {
        const double* col = UT.col(fi).data();
        for (std::size_t p = 0; p < np; ++p) g[p] += gv * col[p];
      }
      const ResponseKind& rk = responses[j];
      // a section at a node: the line jumps there, no join with that run
      const bool sec_kind = rk.what == ResponseKind::What::BarForce || rk.what == ResponseKind::What::BarDisplacement;
      auto joins = [&](std::size_t r) {
        return r > 0 && run_join[r] && !(sec_kind && (runs[r].bar == rk.bar || runs[r - 1].bar == rk.bar));
      };
      // pieces: (s_a, s_b, cubic in t over the piece)
      // join: the piece continues the previous one at a shared station (the
      // structure depends only on the runs and the section, never on values,
      // so lines of one section and both sweeps stay index-aligned)
      // A jump at the section spans delta = 1e-7 L on the far side of the
      // section (refine_path: a load at a == x counts on the bar's origin
      // side): the piece there starts (or ends) delta away, its cubic
      // evaluated at u0 / u1 inside its own parameter range.
      struct Piece {
        double sa, sb;
        Eigen::Vector4d c;
        bool join;
        double u0 = 0, u1 = 1;
        Id bar = 0;
        double xa = 0, xb = 0;  // local abscissae at sa, sb
      };
      std::vector<Piece> pieces;
      for (std::size_t r = 0; r < nr; ++r) {
        const Run& R = runs[r];
        const double Lbar = model.bar_chord_length(model.bar(R.bar));
        double gv4[4] = {g[4 * r], g[4 * r + 1], g[4 * r + 2], g[4 * r + 3]};
        auto lit = local_of_bar.find(R.bar);
        const bool local = lit != local_of_bar.end() &&
                           std::find(lit->second.begin(), lit->second.end(), j) != lit->second.end();
        const bool sec = local && (rk.what == ResponseKind::What::BarForce ||
                                   rk.what == ResponseKind::What::BarDisplacement) && rk.bar == R.bar;
        // the section on this run, strictly inside: a split
        double ts = -1;
        if (sec && std::fabs(R.x1 - R.x0) > 0) {
          const double t = (rk.x - R.x0) / (R.x1 - R.x0);
          if (t > 1e-12 && t < 1 - 1e-12) ts = t;
        }
        if (ts < 0) {
          double v[4];
          for (int q = 0; q < 4; ++q) {
            v[q] = gv4[q];
            if (local) {
              const CaseLoads& cl = pcl[4 * r + static_cast<std::size_t>(q)];
              if (rk.what == ResponseKind::What::BarDisplacement) {
                std::vector<double> dv;
                std::vector<const ResponseKind*> one{&rk};
                extract_bar_displacements(model, asm_, cl, zero, rk.bar, one, dv);
                v[q] += dv[0];
              } else {
                v[q] += extract_response(model, asm_, cl, zero, rk);
              }
            }
          }
          Piece P{R.s0, R.s1, fit(v), joins(r)};
          P.bar = R.bar;
          P.xa = R.x0;
          P.xb = R.x1;
          if (sec && std::fabs(rk.x) <= 1e-9 * Lbar) {  // section at the bar origin
            const bool flip = R.x1 < R.x0;
            const double du = 1e-7 * Lbar / (R.s1 - R.s0);
            if (flip)
              P.u1 = 1 - du;
            else
              P.u0 = du;
          }
          pieces.push_back(P);
          continue;
        }
        // split at the section: the global cubic of the run, the local term
        // sampled inside each side
        const Eigen::Vector4d gc = fit(gv4);
        for (int side = 0; side < 2; ++side) {
          const double ta = side == 0 ? 0.0 : ts, tb = side == 0 ? ts : 1.0;
          double v[4];
          for (int q = 0; q < 4; ++q) {
            const double t = ta + tq[q] * (tb - ta);
            const double x = R.x0 + t * (R.x1 - R.x0);
            const CaseLoads& cl = split_case(R.bar, x, field[R.k0]);
            double loc;
            if (rk.what == ResponseKind::What::BarDisplacement) {
              std::vector<double> dv;
              std::vector<const ResponseKind*> one{&rk};
              extract_bar_displacements(model, asm_, cl, zero, rk.bar, one, dv);
              loc = dv[0];
            } else {
              loc = extract_response(model, asm_, cl, zero, rk);
            }
            v[q] = peval(gc, t) + loc;
          }
          Piece P{R.s0 + ta * (R.s1 - R.s0), R.s0 + tb * (R.s1 - R.s0), fit(v), side == 0 && joins(r)};
          P.bar = R.bar;
          P.xa = R.x0 + ta * (R.x1 - R.x0);
          P.xb = R.x0 + tb * (R.x1 - R.x0);
          const bool flip = R.x1 < R.x0;
          const double du = 1e-7 * Lbar / (P.sb - P.sa);
          if (!flip && side == 1) P.u0 = du;   // x + delta: the increasing-x limit
          if (flip && side == 0) P.u1 = 1 - du;
          pieces.push_back(P);
        }
      }
      // a section at the bar END x = L: the jump spans delta on the side of
      // the next bar (the load at a == L counts on the origin side)
      if (sec_kind) {
        const double Ls = model.bar_chord_length(model.bar(rk.bar));
        if (std::fabs(rk.x - Ls) <= 1e-9 * Ls)
          for (std::size_t q = 1; q < pieces.size(); ++q) {
            Piece &A = pieces[q - 1], &B = pieces[q];
            if (B.join || A.bar == B.bar) continue;
            if (A.bar == rk.bar && std::fabs(A.xb - Ls) <= 1e-9 * Ls)
              B.u0 = 1e-7 * Ls / (B.sb - B.sa);
            else if (B.bar == rk.bar && std::fabs(B.xa - Ls) <= 1e-9 * Ls)
              A.u1 = 1 - 1e-7 * Ls / (A.sb - A.sa);
          }
      }
      // the Hermite line through the pieces (a jump: a double station)
      InfluenceLine& il = res[j];
      for (std::size_t q = 0; q < pieces.size(); ++q) {
        const Piece& P = pieces[q];
        const double h = P.sb - P.sa;
        const double va = peval(P.c, P.u0), vb = peval(P.c, P.u1);
        const double da = pder(P.c, P.u0) / h, db = pder(P.c, P.u1) / h;
        const double sa = P.sa + P.u0 * h, sb = P.sa + P.u1 * h;
        if (P.join && !il.s.empty()) {
          // a node with one load direction on both sides: the line is
          // continuous there, the two ends agree to round-off
          il.eta.back() = 0.5 * (il.eta.back() + va);
          il.dr.back() = da;
        } else {
          il.s.push_back(sa);
          il.eta.push_back(va);
          il.dl.push_back(da);
          il.dr.push_back(da);
        }
        il.s.push_back(sb);
        il.eta.push_back(vb);
        il.dl.push_back(db);
        il.dr.push_back(db);
      }
      for (double e : il.eta) require(std::isfinite(e), "SURCH: ordonnee de ligne d'influence non finie");
    }
    return res;
  };
  std::vector<InfluenceLine> out = sweep(station_dir, false);
  if (torque_dir && out_torque) *out_torque = sweep(*torque_dir, true);
  ok = true;
  return out;
}

std::vector<InfluenceLine> il_batch(const Model& model, const LoadPath& path,
                                    const std::vector<Vec3>& station_dir,
                                    const std::vector<ResponseKind>& responses,
                                    const std::vector<Vec3>* torque_dir = nullptr,
                                    std::vector<InfluenceLine>* out_torque = nullptr) {
  const std::size_t ns = path.stations.size();
  require(station_dir.size() == ns, "SURCH: champ de directions incoherent");
  require(!responses.empty(), "SURCH: aucune reponse demandee");
  if (!std::getenv("KETA_IL_SAMPLED")) {  // Müller-Breslau: exact pieces
    bool ok = false;
    std::vector<InfluenceLine> ex =
        il_batch_exact(model, path, station_dir, responses, torque_dir, out_torque, ok);
    if (ok) return ex;
  }
  model.validate();
  Assembled asm_ = assemble_elements(model);
  SpMat K = assemble_stiffness(model, asm_);
  FactorCache fac;
  fac.factorize(K);
  const DofMap& dm = asm_.dofs;

  // Gradients g (indices libres, coefficients) par sondes unitaires.
  const CaseLoads empty = empty_case_loads(asm_);
  std::vector<std::vector<std::pair<int, double>>> grad(responses.size());
  {
    VecX probe = VecX::Zero(dm.n_free);
    for (std::size_t j = 0; j < responses.size(); ++j)
      for (int f : response_dofs(model, asm_, responses[j])) {
        int fi = dm.free_ix[static_cast<std::size_t>(f)];
        if (fi < 0) continue;
        probe(fi) = 1;
        double gv = extract_response(model, asm_, empty, probe, responses[j]);
        probe(fi) = 0;
        if (gv != 0.0) grad[j].emplace_back(fi, gv);
      }
  }
  // Réponses dont la charge posée sur la barre `b` porte un terme local.
  // node -> the bars ending on it, in model order (each bar once), built on
  // the first reaction response: a scan of every bar per reaction was O(R x B)
  std::unordered_map<Id, std::vector<Id>> incident;
  bool incident_built = false;
  auto bars_at_node = [&](Id node) -> const std::vector<Id>& {
    if (!incident_built) {
      for (const Bar& b : model.bars) {
        incident[b.ni].push_back(b.id);
        if (b.nj != b.ni) incident[b.nj].push_back(b.id);
      }
      incident_built = true;
    }
    static const std::vector<Id> none;
    auto it = incident.find(node);
    return it == incident.end() ? none : it->second;
  };
  std::map<Id, std::vector<std::size_t>> local_of_bar;
  for (std::size_t j = 0; j < responses.size(); ++j) {
    const ResponseKind& r = responses[j];
    switch (r.what) {
      case ResponseKind::What::Displacement: break;
      case ResponseKind::What::BarForce:
      case ResponseKind::What::BarDisplacement: local_of_bar[r.bar].push_back(j); break;
      case ResponseKind::What::Reaction:
        for (Id bid : bars_at_node(r.node)) local_of_bar[bid].push_back(j);
        break;
    }
  }

  // Probe points: the path stations, then 4 points inside each interval (at
  // e and 2e from both ends, on the bar that carries it) whose ordinates give
  // the exact slopes of the cubic line (influence.h, InfluenceLine::dl/dr).
  struct Probe {
    Id bar = 0;
    double x = 0;
    std::size_t field_k = 0;  // the station whose load direction it takes
  };
  std::vector<Probe> probes;
  for (std::size_t k = 0; k < ns; ++k) probes.push_back({path.stations[k].bar, path.stations[k].x_local, k});
  struct Iv {
    bool ok = false;
    double e = 0;
    std::size_t first = 0;
  };
  std::vector<Iv> ivs(ns > 0 ? ns - 1 : 0);
  {
    auto bar_end = [&](const PathStation& st) {
      const double Lb = model.bar_chord_length(model.bar(st.bar));
      return st.x_local <= 1e-9 * Lb || st.x_local >= Lb * (1 - 1e-9);
    };
    for (std::size_t k = 0; k + 1 < ns; ++k) {
      const PathStation &a = path.stations[k], &b = path.stations[k + 1];
      const double h = b.s - a.s;
      if (h <= 0) continue;
      const PathStation& ref = a.bar == b.bar ? a : (bar_end(a) ? b : a);
      const std::size_t fk = &ref == &a ? k : k + 1;
      double dir_s = 0;
      for (std::size_t q = 0; q < ns; ++q)
        if (path.stations[q].bar == ref.bar && path.stations[q].s != ref.s) {
          dir_s = (path.stations[q].x_local - ref.x_local) / (path.stations[q].s - ref.s) > 0 ? 1.0 : -1.0;
          break;
        }
      if (dir_s == 0) continue;
      ivs[k].ok = true;
      ivs[k].e = 1e-3 * h;
      ivs[k].first = probes.size();
      for (double sp : {a.s + ivs[k].e, a.s + 2 * ivs[k].e, b.s - ivs[k].e, b.s - 2 * ivs[k].e})
        probes.push_back({ref.bar, ref.x_local + dir_s * (sp - ref.s), fk});
    }
  }
  const std::size_t np = probes.size();

  // Un champ de charges unitaires (forces, ou couples pour l'excentrement)
  // balayé le long du chemin. Direction constante -> un seul balayage ;
  // sinon les 3 axes globaux (la réponse étant linéaire en la direction).
  const VecX zero = VecX::Zero(dm.n_free);
  auto sweep_field = [&](const std::vector<Vec3>& field, bool moment) {
    std::vector<std::vector<double>> val(responses.size(), std::vector<double>(np, 0.0));
    bool uniform = true;
    for (std::size_t k = 1; k < ns; ++k)
      if ((field[k] - field[0]).lpNorm<Eigen::Infinity>() > 1e-12) uniform = false;
    std::vector<Vec3> sweeps;
    if (uniform)
      sweeps = {field.empty() ? Vec3(0, 0, -1) : field[0]};
    else
      sweeps = {Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    for (std::size_t d = 0; d < sweeps.size(); ++d) {
      auto unit_case = [&](std::size_t p) {
        LoadCase lc;
        lc.id = 1;
        lc.name = "unit";
        const Probe& pr = probes[p];
        lc.loads.push_back(moment ? Load::point(pr.bar, pr.x, Vec3::Zero(), sweeps[d])
                                  : Load::point(pr.bar, pr.x, sweeps[d]));
        return lower_case(model, asm_, lc);
      };
      MatX rhs(dm.n_free, static_cast<Eigen::Index>(np));
      for (std::size_t p = 0; p < np; ++p) rhs.col(static_cast<Eigen::Index>(p)) = unit_case(p).rhs;
      MatX U = fac.solve(rhs);
      // the global part, response by response over all the probes at once:
      // U transposed so that one DOF across the probes is contiguous (the
      // same sums in the same order as probe by probe)
      {
        const MatX UT = U.transpose();
        std::vector<double> wv(np);
        for (std::size_t p = 0; p < np; ++p)
          wv[p] = uniform ? 1.0 : field[probes[p].field_k](static_cast<Eigen::Index>(d));
        std::vector<double> v(np);
        for (std::size_t j = 0; j < responses.size(); ++j) {
          std::fill(v.begin(), v.end(), 0.0);
          for (const auto& [fi, gv] : grad[j]) {
            const double* col = UT.col(fi).data();
            for (std::size_t p = 0; p < np; ++p) v[p] += gv * col[p];
          }
          double* out = val[j].data();
          for (std::size_t p = 0; p < np; ++p)
            if (wv[p] != 0.0) out[p] += wv[p] * v[p];
        }
      }
      for (std::size_t p = 0; p < np; ++p) {
        double w = uniform ? 1.0 : field[probes[p].field_k](static_cast<Eigen::Index>(d));
        if (w == 0.0) continue;
        auto loc = local_of_bar.find(probes[p].bar);
        if (loc != local_of_bar.end()) {
          CaseLoads cl = unit_case(p);
          // the bar displacements of the probe's bar from ONE recovery
          std::vector<const ResponseKind*> disp;
          std::vector<std::size_t> disp_j;
          for (std::size_t j : loc->second) {
            if (responses[j].what == ResponseKind::What::BarDisplacement && responses[j].bar == probes[p].bar) {
              disp.push_back(&responses[j]);
              disp_j.push_back(j);
              continue;
            }
            val[j][p] += w * extract_response(model, asm_, cl, zero, responses[j]);
          }
          if (!disp.empty()) {
            std::vector<double> dv;
            extract_bar_displacements(model, asm_, cl, zero, probes[p].bar, disp, dv);
            for (std::size_t q = 0; q < disp.size(); ++q) val[disp_j[q]][p] += w * dv[q];
          }
        }
      }
    }
    std::vector<InfluenceLine> res(responses.size());
    for (std::size_t j = 0; j < responses.size(); ++j) {
      InfluenceLine& il = res[j];
      il.s.assign(ns, 0.0);
      il.eta.assign(ns, 0.0);
      il.dl.assign(ns, 0.0);
      il.dr.assign(ns, 0.0);
      for (std::size_t k = 0; k < ns; ++k) {
        il.s[k] = path.stations[k].s;
        il.eta[k] = val[j][k];
      }
      for (std::size_t k = 0; k + 1 < ns; ++k) {
        if (!ivs[k].ok) continue;
        const double e = ivs[k].e;
        const std::size_t b0 = ivs[k].first;
        const double f1 = val[j][b0], f2 = val[j][b0 + 1], g1 = val[j][b0 + 2], g2 = val[j][b0 + 3];
        const double fa = 2 * f1 - f2, gb = 2 * g1 - g2;
        il.dr[k] = (-3 * fa + 4 * f1 - f2) / (2 * e);
        il.dl[k + 1] = (3 * gb - 4 * g1 + g2) / (2 * e);
      }
    }
    return res;
  };
  std::vector<InfluenceLine> out = sweep_field(station_dir, false);
  if (torque_dir && out_torque) {
    require(torque_dir->size() == ns, "SURCH: champ de couples incoherent");
    *out_torque = sweep_field(*torque_dir, true);
    for (const InfluenceLine& il : *out_torque)
      for (double e : il.eta) require(std::isfinite(e), "SURCH: ordonnee de couple non finie");
  }
  for (const InfluenceLine& il : out)
    for (double e : il.eta) require(std::isfinite(e), "SURCH: ordonnee de ligne d'influence non finie");
  return out;
}

}  // namespace

// ===========================================================================
// 6. ZONES DE SIGNE ET LONGUEUR CHARGÉE
// ===========================================================================

namespace {
int cubic_roots01(double c0, double c1, double c2, double c3, double out[3]);
}

std::vector<SignZone> il_sign_zones(const InfluenceLine& il, double seuil) {
  require(seuil >= 0 && seuil < 1, "SURCH SEUIL: doit etre dans [0, 1[");
  std::vector<SignZone> zones;
  const std::size_t n = il.eta.size();
  if (n < 2) return zones;
  // the extreme ordinates, the interior extrema of a cubic line included
  const auto [emax, emin] = il.point_extreme();
  double ymax = std::max(std::fabs(emax), std::fabs(emin));
  for (double e : il.eta) ymax = std::max(ymax, std::fabs(e));
  if (ymax <= 0) return zones;
  const double thr = ymax * seuil;

  // p.140 : un point appartient à une zone positive si y > ymax·seuil, à une
  // zone négative si y < -ymax·seuil. Une zone est donc une suite maximale
  // de tels points, bornée aux franchissements (interpolés) du seuil : les
  // extrémités en deçà du seuil ne sont pas chargées, et un creux sous le
  // seuil sépare deux zones de même signe, chargeables séparément (Exemple3,
  // tranchant au droit d'un hauban). Un point d'ordonnée nulle (appui, saut
  // d'un effort tranchant à une jonction) n'appartient à aucune zone (p.141).
  // Avec SEUIL 0 les bornes sont les zéros exacts.
  // The line is walked interval by interval: inside each, the crossings of
  // +thr and -thr (roots of the exact cubic, or of the chord) cut it into
  // pieces of constant class (+1 above thr, -1 below -thr, 0 between); a
  // zone is a maximal chain of pieces of one class, broken at any crossing
  // and at any station whose ordinate (either side of a double station)
  // falls back to the threshold. A sign zone may thus lie wholly between
  // two stations — an exact line has a handful per bar.
  const bool cub = il.cubic();
  int zs = 0;  // class of the open zone (0: none)
  double z0 = 0;
  auto close = [&](double at) {
    if (zs != 0 && at > z0) {
      SignZone z;
      z.s0 = z0;
      z.s1 = at;
      z.sign = zs;
      z.integral = il.integral(z0, at);
      zones.push_back(z);
    }
    zs = 0;
  };
  for (std::size_t i = 0; i + 1 < n; ++i) {
    // the station itself: a zone does not run through an ordinate within the
    // threshold (support, zero of a shear jump)
    if (zs != 0 && !(zs * il.eta[i] > thr)) close(il.s[i]);
    const double a = il.s[i], h = il.s[i + 1] - a;
    if (!(h > 0)) continue;
    const double e0 = il.eta[i], e1 = il.eta[i + 1];
    double c[4] = {e0, e1 - e0, 0, 0};
    if (cub) {
      const double m0 = h * il.dr[i], m1 = h * il.dl[i + 1];
      c[1] = m0;
      c[2] = -3 * e0 - 2 * m0 + 3 * e1 - m1;
      c[3] = 2 * e0 + m0 - 2 * e1 + m1;
    }
    double ts[8];
    int nt = 0;
    ts[nt++] = 0.0;
    for (double lv : {thr, -thr}) {
      double r[3];
      const int nr = cubic_roots01(c[0] - lv, c[1], c[2], c[3], r);
      for (int k = 0; k < nr; ++k) ts[nt++] = r[k];
      if (thr == 0) break;
    }
    std::sort(ts + 1, ts + nt);
    ts[nt++] = 1.0;
    auto p = [&](double t) { return c[0] + t * (c[1] + t * (c[2] + t * c[3])); };
    for (int k = 0; k + 1 < nt; ++k) {
      const double ta = ts[k], tb = ts[k + 1];
      if (!(tb > ta)) continue;
      const double v = p(0.5 * (ta + tb));
      const int cls = v > thr ? 1 : (v < -thr ? -1 : 0);
      const double sa = a + ta * h;
      if (k > 0 || cls != zs) close(sa);  // a crossing always cuts
      if (cls != 0 && zs == 0) {
        zs = cls;
        z0 = sa;
      }
    }
  }
  close(il.s.back());
  return zones;
}

namespace {

// Choix des zones favorables sous densité DÉPENDANT de la longueur chargée.
// q étant décroissante, « tout charger » n'est pas toujours l'optimum : on
// énumère les sous-ensembles (<= 2^12) et on retient le meilleur ; au-delà,
// repli glouton par |contribution| / longueur décroissante.
struct LoadedChoice {
  std::vector<SignZone> zones;
  double length = 0;
  double integral = 0;  // somme signée des intégrales retenues
  double value = 0;     // q(L) * integral (avant facteur transversal)
};

LoadedChoice pick_zones(const std::vector<SignZone>& all, double gain,
                        const std::function<double(double)>& q, bool maximise) {
  // contribution unitaire d'une zone = gain * integrale
  std::vector<SignZone> cand;
  for (const SignZone& z : all) {
    double c = gain * z.integral;
    if (maximise ? (c > 0) : (c < 0)) cand.push_back(z);
  }
  LoadedChoice best;
  if (cand.empty()) return best;
  auto score = [&](const std::vector<SignZone>& sel) {
    double L = 0, I = 0;
    for (const SignZone& z : sel) {
      L += z.s1 - z.s0;
      I += z.integral;
    }
    if (L <= 0) return LoadedChoice{};
    LoadedChoice lc;
    lc.zones = sel;
    lc.length = L;
    lc.integral = I;
    lc.value = q(L) * gain * I;
    return lc;
  };
  if (cand.size() <= 12) {
    const std::uint32_t nmask = 1u << cand.size();
    bool first = true;
    for (std::uint32_t mask = 1; mask < nmask; ++mask) {
      std::vector<SignZone> sel;
      for (std::size_t i = 0; i < cand.size(); ++i)
        if (mask & (1u << i)) sel.push_back(cand[i]);
      LoadedChoice lc = score(sel);
      if (lc.length <= 0) continue;
      if (first || (maximise ? lc.value > best.value : lc.value < best.value)) {
        best = lc;
        first = false;
      }
    }
  } else {
    std::vector<SignZone> sorted = cand;
    std::sort(sorted.begin(), sorted.end(), [](const SignZone& a, const SignZone& b) {
      double la = std::max(1e-12, a.s1 - a.s0), lb = std::max(1e-12, b.s1 - b.s0);
      return std::fabs(a.integral) / la > std::fabs(b.integral) / lb;
    });
    bool first = true;
    for (std::size_t k = 1; k <= sorted.size(); ++k) {
      std::vector<SignZone> sel(sorted.begin(), sorted.begin() + static_cast<long>(k));
      LoadedChoice lc = score(sel);
      if (lc.length <= 0) continue;
      if (first || (maximise ? lc.value > best.value : lc.value < best.value)) {
        best = lc;
        first = false;
      }
    }
  }
  return best;
}

}  // namespace

// ===========================================================================
// 7. PLACEMENT DES CONVOIS
// ===========================================================================

namespace {

// The distributed segments integrate through the line's prefix sums (O(log n)
// per segment instead of a sweep of the whole line: the placement evaluates
// this at thousands of candidate positions for every target — 4 min 40 s of
// cable_et_phasage's 4 min 44 s went into that sweep).
double pattern_effect(const InfluenceLine& il, const IlPrefix& pre, const ConvoyPattern& p, double r) {
  double E = 0;
  for (const auto& [a, P] : p.axles) E += P * il.at(r + a);
  for (const auto& seg : p.segments)
    if (seg[1] > 0) E += seg[2] / seg[1] * pre.integral(r + seg[0], r + seg[0] + seg[1]);
  return E;
}

// Positions candidates : grille de pas `step` UNION les points de rupture
// exacts {s_k - offset}. La LI étant affine par morceaux, l'effet d'un train
// d'essieux PONCTUELS est affine par morceaux de sommets exactement dans cet
// ensemble : le balayage est donc exact.
// `grid_only` (balayage de ST1, p.147) : la seule grille de pas `step`
// EXACTEMENT, depuis l'entrée du convoi sur le chemin (dernier essieu à
// l'origine) jusqu'à sa sortie — l'extrême trouvé est alors celui de la grille,
// égal ou inférieur à l'extrême exact.
std::vector<double> candidate_positions(const InfluenceLine& il, const ConvoyPattern& p,
                                        double step, bool grid_only = false, bool with_grid = true) {
  double amin = 0, amax = 0;
  std::vector<double> offs;
  for (const auto& [a, P] : p.axles) {
    (void)P;
    amin = std::min(amin, a);
    amax = std::max(amax, a);
    offs.push_back(a);
  }
  for (const auto& seg : p.segments) {
    amin = std::min(amin, seg[0]);
    amax = std::max(amax, seg[0] + seg[1]);
    offs.push_back(seg[0]);
    offs.push_back(seg[0] + seg[1]);
  }
  const double r0 = il.s.front() - amax, r1 = il.s.back() - amin;
  std::vector<double> pos;
  require(step > 0 && finite(step), "SURCH PAS_CONV: le pas doit etre > 0");
  if (grid_only) {
    for (long k = 0;; ++k) {
      const double r = r0 + static_cast<double>(k) * step;
      if (r > r1 + 1e-9) break;
      pos.push_back(r);
    }
    return pos;
  }
  int n = static_cast<int>(std::ceil((r1 - r0) / step));
  n = std::max(n, 1);
  if (with_grid)
    for (int k = 0; k <= n; ++k) pos.push_back(r0 + (r1 - r0) * k / n);
  else
    pos.insert(pos.end(), {r0, r1});
  for (double s : il.s)
    for (double a : offs) {
      double r = s - a;
      if (r >= r0 - 1e-9 && r <= r1 + 1e-9) pos.push_back(r);
    }
  std::sort(pos.begin(), pos.end());
  pos.erase(std::unique(pos.begin(), pos.end(),
                        [](double a, double b) { return std::fabs(a - b) < 1e-10; }),
            pos.end());
  return pos;
}

// Real roots in (0, 1) of c0 + c1 t + c2 t^2 + c3 t^3: the monotone pieces
// between the roots of its derivative, each bisected where it changes sign.
// At most 3 roots, written to `out`; returns their count (no allocation).
int cubic_roots01(double c0, double c1, double c2, double c3, double out[3]) {
  auto f = [&](double t) { return c0 + t * (c1 + t * (c2 + t * c3)); };
  double cut[4] = {0.0, 1.0, 0.0, 0.0};
  int nc = 2;
  // f' = c1 + 2 c2 t + 3 c3 t^2
  const double A = 3 * c3, B = 2 * c2, C = c1;
  if (std::fabs(A) > 1e-300) {
    const double disc = B * B - 4 * A * C;
    if (disc >= 0) {
      const double sq = std::sqrt(disc), q = -0.5 * (B + (B >= 0 ? sq : -sq));
      for (double t : {q / A, q != 0 ? C / q : -B / (2 * A)})
        if (t > 0 && t < 1) cut[nc++] = t;
    }
  } else if (std::fabs(B) > 1e-300) {
    const double t = -C / B;
    if (t > 0 && t < 1) cut[nc++] = t;
  }
  std::sort(cut, cut + nc);
  int n = 0;
  for (int i = 0; i + 1 < nc; ++i) {
    double lo = cut[i], hi = cut[i + 1];
    double flo = f(lo), fhi = f(hi);
    if (flo == 0 && lo > 0) out[n++] = lo;
    if (flo * fhi >= 0) continue;
    for (int it = 0; it < 200 && hi - lo > 1e-15; ++it) {
      const double mid = 0.5 * (lo + hi), fm = f(mid);
      if ((fm < 0) == (flo < 0)) {
        lo = mid;
        flo = fm;
      } else {
        hi = mid;
      }
    }
    out[n++] = 0.5 * (lo + hi);
  }
  return n;
}

// The line as polynomials in d = x - s_i per interval: eta = sum a[k] d^k and
// its antiderivative F = cum_i + sum A[k] d^k (the Hermite cubic exactly).
struct IlPoly {
  std::vector<std::array<double, 4>> a;
  std::vector<std::array<double, 5>> A;
  std::vector<double> cum;
  void build(const InfluenceLine& il) {
    const std::size_t n = il.s.size();
    a.assign(n, {0, 0, 0, 0});
    A.assign(n, {0, 0, 0, 0, 0});
    cum.assign(n, 0.0);
    for (std::size_t i = 0; i + 1 < n; ++i) {
      const double h = il.s[i + 1] - il.s[i];
      if (h > 0) {
        const double e0 = il.eta[i], e1 = il.eta[i + 1];
        double c[4] = {e0, e1 - e0, 0, 0};
        if (il.cubic()) {
          const double m0 = h * il.dr[i], m1 = h * il.dl[i + 1];
          c[1] = m0;
          c[2] = -3 * e0 - 2 * m0 + 3 * e1 - m1;
          c[3] = 2 * e0 + m0 - 2 * e1 + m1;
        }
        double hk = 1;
        for (int k = 0; k < 4; ++k) {
          a[i][static_cast<std::size_t>(k)] = c[k] / hk;
          A[i][static_cast<std::size_t>(k) + 1] = c[k] / ((k + 1) * hk);
          hk *= h;
        }
        double tot = 0, hp = h;
        for (int k = 1; k <= 4; ++k) {
          tot += A[i][static_cast<std::size_t>(k)] * hp;
          hp *= h;
        }
        cum[i + 1] = cum[i] + tot;
      } else {
        cum[i + 1] = cum[i];
      }
      A[i][0] = cum[i];
    }
  }
};

// The effect of a convoy at every SORTED position, exactly pattern_effect's
// arithmetic (same terms, same order), each term read by a pointer that only
// moves forward instead of a binary search per position.
void effects_sorted(const InfluenceLine& il, const IlPrefix& pre, const ConvoyPattern& p,
                    const std::vector<double>& pos, std::vector<double>& E) {
  const std::size_t m = pos.size(), ns = il.s.size();
  E.assign(m, 0.0);
  if (ns == 0) return;
  const auto& sv = il.s;
  // the interval of x (s[i] <= x < s[i+1]) moving forward from i
  auto advance = [&](std::size_t& i, double x) {
    while (i + 2 < ns && sv[i + 1] <= x) ++i;
  };
  for (const auto& [a, P] : p.axles) {
    std::size_t i = 0;
    for (std::size_t k = 0; k < m; ++k) {
      const double x = pos[k] + a;
      double v;
      if (x <= sv.front()) v = x < sv.front() ? 0.0 : il.eta.front();
      else if (x >= sv.back()) v = x > sv.back() ? 0.0 : il.eta.back();
      else {
        advance(i, x);
        v = il.on(i, x);
      }
      E[k] += P * v;
    }
  }
  for (const auto& seg : p.segments) {
    if (!(seg[1] > 0)) continue;
    const double q = seg[2] / seg[1];
    std::size_t ia = 0, ib = 0;
    auto prefix = [&](std::size_t& i, double x) {
      if (ns < 2 || x <= sv.front()) return 0.0;
      if (x >= sv.back()) return pre.total();
      advance(i, x);
      return pre.prefix_in(i, x);
    };
    for (std::size_t k = 0; k < m; ++k) {
      const double xa = pos[k] + seg[0], xb = xa + seg[1];
      const double fb = prefix(ib, xb), fa = prefix(ia, xa);
      E[k] += q * (xb <= xa ? 0.0 : fb - fa);
    }
  }
}

// What a placement needs that does not depend on the sign sought: the
// candidates (grid / breakpoints), their effects, and for a cubic line the
// effect polynomial of every interval between breakpoints, for a single truck
// and for files of c trucks at the minimum pitch (MAX_CAM). Kept per thread:
// the max and the min of one line, and the same line met again, reuse it.
struct PlacementPrep {
  // key
  std::vector<double> s, eta, dl, dr;
  std::vector<std::pair<double, double>> axles;
  std::vector<std::array<double, 3>> segs;
  int max_cam = 0;
  double pitch = 0, step = 0;
  bool grid_only = false, valid = false;
  // data
  std::vector<double> pos, E;
  struct Region {  // one truck: the effect polynomial and its bounds
    double ra, H;
    double e[5];
    double up, dn;  // upper bounds of +E and -E on the interval
    double B;       // cheap bound of |E| on the interval
    bool built;     // e / up / dn computed
  };
  std::vector<Region> regions;
  double gmax = 0, gmin = 0;  // upper bounds of +E and -E over the whole line
  std::vector<double> mabs, apre;  // max |eta| per interval, prefix of mabs h
  bool ip_built = false;
  bool lazy = false;              // effects at the candidates computed on demand
  std::shared_ptr<struct PlacementGeom> geom;  // the geometry it was prepared on
  std::vector<double> cub;        // per candidate: bound of |E| (its intervals)
  std::vector<std::size_t> order; // candidates by decreasing bound
  std::vector<char> have;         // E[k] computed
  // scratch
  IlPoly ip;
  std::vector<double> br, extra, Eextra, allpos, allE;
  std::vector<int> prev;
  std::vector<double> f;
  std::vector<int> take;
};
// a few preparations per thread: the max and the min of a target come in
// separate sweeps over its lanes and directions
thread_local std::list<PlacementPrep> tl_preps;  // stable addresses

bool same_key(const PlacementPrep& k, const InfluenceLine& il, const ConvoyPattern& p, int max_cam,
              double pitch, double step, bool grid_only) {
  return k.valid && k.max_cam == max_cam && k.pitch == pitch && k.step == step &&
         k.grid_only == grid_only && k.s == il.s && k.eta == il.eta && k.dl == il.dl &&
         k.dr == il.dr && k.axles == p.axles && k.segs == p.segments;
}

// What depends only on the line's ABSCISSAE and the pattern — identical for
// every target of one SURCH (one load path): the candidates and, per file of
// c trucks, the intervals between breakpoints with, for each term, the
// interval of the line it reads and its offset there. Kept per thread.
struct PlacementGeom {
  std::vector<double> s;
  std::vector<std::pair<double, double>> axles;
  std::vector<std::array<double, 3>> segs;
  int max_cam = 0;
  double pitch = 0, step = 0;
  bool grid_only = false, stationary = false, valid = false;
  std::vector<double> pos;
  struct Term {
    double o, w;
    bool F;
  };
  std::vector<std::vector<Term>> terms;  // per c - 1
  struct Region {
    double ra, H;
    int c;
    std::size_t first;  // c = 1: into idx / sh (terms[0].size() entries);
                        // c > 1: into r1 / off (c entries)
  };
  std::vector<Region> regions;  // the single-truck ones first, by ra
  std::size_t n1 = 0;           // their count
  std::vector<int> idx;     // interval read (-1 before the line, -2 after it)
  std::vector<double> sh;   // offset d = u + sh in that interval
  // a file of c trucks: truck j reads single-truck region r1, at u + off
  std::vector<std::size_t> r1;
  std::vector<double> off;
};
// a few geometries per thread: the patterns alternate (both directions,
// several convoys) over the same line abscissae
thread_local std::list<std::shared_ptr<PlacementGeom>> tl_geoms;

void prepare_geom(PlacementGeom& Gm, const InfluenceLine& il, const ConvoyPattern& p, int max_cam,
                  double pitch, double step, bool grid_only) {
  Gm.s = il.s;
  Gm.axles = p.axles;
  Gm.segs = p.segments;
  Gm.max_cam = max_cam;
  Gm.pitch = pitch;
  Gm.step = step;
  Gm.grid_only = grid_only;
  Gm.stationary = il.cubic() && !grid_only;
  Gm.valid = true;
  // On a cubic line the effect between breakpoints is a polynomial whose
  // extremes are its ends or its stationary points: the PAS_CONV grid adds
  // nothing there (kept for an affine line, and for --pas-conv).
  Gm.pos = candidate_positions(il, p, step, grid_only, /*with_grid=*/!Gm.stationary);
  Gm.terms.clear();
  Gm.regions.clear();
  Gm.idx.clear();
  Gm.sh.clear();
  Gm.r1.clear();
  Gm.off.clear();
  Gm.n1 = 0;
  if (!Gm.stationary || Gm.pos.size() < 2) return;
  const auto& sv = il.s;
  std::vector<PlacementGeom::Term> one;
  for (const auto& [a, w] : p.axles) one.push_back({a, w, false});
  for (const auto& seg : p.segments)
    if (seg[1] > 0) {
      const double q = seg[2] / seg[1];
      one.push_back({seg[0] + seg[1], q, true});
      one.push_back({seg[0], -q, true});
    }
  const double r0 = Gm.pos.front(), r1 = Gm.pos.back();
  const int cmax = std::max(1, std::min(max_cam, static_cast<int>((r1 - r0) / pitch) + 1));
  std::vector<double> br;
  std::vector<std::size_t> ptr;
  for (int c = 1; c <= cmax; ++c) {
    const double span = (c - 1) * pitch;
    if (r1 - span < r0) break;
    std::vector<PlacementGeom::Term> terms;
    for (int j = 0; j < c; ++j)
      for (const PlacementGeom::Term& t : one) terms.push_back({t.o + j * pitch, t.w, t.F});
    br.assign({r0, r1 - span});
    for (const PlacementGeom::Term& t : terms)
      for (double x : sv) {
        const double r = x - t.o;
        if (r > r0 && r < r1 - span) br.push_back(r);
      }
    std::sort(br.begin(), br.end());
    br.erase(std::unique(br.begin(), br.end(), [](double x, double y) { return std::fabs(x - y) < 1e-10; }),
             br.end());
    ptr.assign(terms.size(), 0);
    if (c == 1) Gm.n1 = 0;
    std::vector<std::size_t> rp(static_cast<std::size_t>(c), 0);
    for (std::size_t q = 0; q + 1 < br.size(); ++q) {
      const double ra = br[q], H = br[q + 1] - ra;
      if (H <= 1e-9) continue;
      const double rm = ra + 0.5 * H;
      if (c > 1) {
        // the file's polynomial is the sum of its trucks' single-truck ones
        Gm.regions.push_back({ra, H, c, Gm.r1.size()});
        for (int j = 0; j < c; ++j) {
          std::size_t& k = rp[static_cast<std::size_t>(j)];
          const double x = rm + j * pitch;
          while (k + 1 < Gm.n1 && Gm.regions[k + 1].ra <= x) ++k;
          Gm.r1.push_back(k);
          Gm.off.push_back(ra + j * pitch - Gm.regions[k].ra);
        }
        continue;
      }
      Gm.regions.push_back({ra, H, c, Gm.idx.size()});
      ++Gm.n1;
      for (std::size_t ti = 0; ti < terms.size(); ++ti) {
        const double xm = rm + terms[ti].o;
        if (xm <= sv.front()) {
          Gm.idx.push_back(-1);
          Gm.sh.push_back(0);
        } else if (xm >= sv.back()) {
          Gm.idx.push_back(-2);
          Gm.sh.push_back(0);
        } else {
          std::size_t& i = ptr[ti];  // xm increases with q: forward only
          while (i + 2 < sv.size() && sv[i + 1] <= xm) ++i;
          Gm.idx.push_back(static_cast<int>(i));
          Gm.sh.push_back(ra + terms[ti].o - sv[i]);
        }
      }
    }
    Gm.terms.push_back(std::move(terms));
  }
}

void prepare(PlacementPrep& P, const InfluenceLine& il, const IlPrefix& pre, const ConvoyPattern& p,
             int max_cam, double pitch, double step, bool grid_only) {
  P.s = il.s;
  P.eta = il.eta;
  P.dl = il.dl;
  P.dr = il.dr;
  P.axles = p.axles;
  P.segs = p.segments;
  P.max_cam = max_cam;
  P.pitch = pitch;
  P.step = step;
  P.grid_only = grid_only;
  P.valid = true;
  std::shared_ptr<PlacementGeom> found;
  for (const auto& gp : tl_geoms) {
    const PlacementGeom& g = *gp;
    if (g.valid && g.max_cam == max_cam && g.pitch == pitch && g.step == step &&
        g.grid_only == grid_only && g.stationary == (il.cubic() && !grid_only) && g.axles == p.axles &&
        g.segs == p.segments && g.s == il.s) {
      found = gp;
      break;
    }
  }
  if (!found) {
    if (tl_geoms.size() >= 16) tl_geoms.pop_front();
    found = std::make_shared<PlacementGeom>();
    prepare_geom(*found, il, p, max_cam, pitch, step, grid_only);
    tl_geoms.push_back(found);
  }
  P.geom = found;  // held by the preparation: never evicted under it
  PlacementGeom& Gm = *found;
  P.pos = Gm.pos;
  P.regions.clear();
  P.gmax = P.gmin = 0;
  P.lazy = !Gm.regions.empty();
  if (!P.lazy) {  // an affine line / the PAS_CONV grid: every candidate
    effects_sorted(il, pre, p, P.pos, P.E);
    for (double e : P.E) {
      P.gmax = std::max(P.gmax, e);
      P.gmin = std::max(P.gmin, -e);
    }
    return;
  }
  P.ip_built = false;
  // max |eta| on each interval (its ends and the interior extrema of the
  // cubic) and the prefix of max |eta| h: the cheap bounds of the effect
  const std::size_t ns = il.s.size();
  P.mabs.assign(ns, 0.0);
  P.apre.assign(ns + 1, 0.0);
  for (std::size_t i = 0; i + 1 < ns; ++i) {
    const double h = il.s[i + 1] - il.s[i];
    double m = std::max(std::fabs(il.eta[i]), std::fabs(il.eta[i + 1]));
    if (h > 0 && il.cubic()) {
      const double e0 = il.eta[i], e1 = il.eta[i + 1], m0 = il.dr[i] * h, m1 = il.dl[i + 1] * h;
      const double A = 3 * (2 * e0 + m0 - 2 * e1 + m1), B = 2 * (-3 * e0 - 2 * m0 + 3 * e1 - m1), C = m0;
      double ts[2];
      int nt = 0;
      if (std::fabs(A) > 1e-300) {
        const double disc = B * B - 4 * A * C;
        if (disc >= 0) {
          ts[nt++] = (-B + std::sqrt(disc)) / (2 * A);
          ts[nt++] = (-B - std::sqrt(disc)) / (2 * A);
        }
      } else if (std::fabs(B) > 1e-300) {
        ts[nt++] = -C / B;
      }
      for (int q = 0; q < nt; ++q)
        if (ts[q] > 0 && ts[q] < 1) m = std::max(m, std::fabs(il.on(i, il.s[i] + ts[q] * h)));
    }
    P.mabs[i] = m * (1 + 1e-12);
    P.apre[i + 1] = P.apre[i] + P.mabs[i] * std::max(h, 0.0);
  }
  P.apre[ns] = P.apre[ns - 1];
  const std::vector<PlacementGeom::Term>& terms = Gm.terms[0];
  P.regions.resize(Gm.n1);
  double gB = 0;
  for (std::size_t g = 0; g < Gm.n1; ++g) {
    const PlacementGeom::Region& G = Gm.regions[g];
    double B = 0;
    for (std::size_t ti = 0; ti < terms.size(); ++ti) {
      const PlacementGeom::Term& t = terms[ti];
      const int i = Gm.idx[G.first + ti];
      if (!t.F) {
        if (i >= 0) B += std::fabs(t.w) * P.mabs[static_cast<std::size_t>(i)];
        continue;
      }
      // a segment: its end term (+q F) then its start term (-q F): q times
      // the integral of |eta| over the intervals it can cover
      const int ie = i, is = Gm.idx[G.first + ti + 1];
      const std::size_t lo = is == -1 ? 0 : is == -2 ? ns : static_cast<std::size_t>(is);
      const std::size_t hi = ie == -1 ? 0 : ie == -2 ? ns : static_cast<std::size_t>(ie) + 1;
      if (hi > lo) B += std::fabs(t.w) * (P.apre[hi] - P.apre[lo]);
      ++ti;
    }
    PlacementPrep::Region& R = P.regions[g];
    R.ra = G.ra;
    R.H = G.H;
    R.B = B;
    R.built = false;
    gB = std::max(gB, B);
  }
  P.gmax = std::max(P.gmax, gB);
  P.gmin = std::max(P.gmin, gB);
  // each candidate's bound: the intervals it bounds (any covering one)
  const std::size_t m = P.pos.size();
  P.cub.assign(m, std::numeric_limits<double>::infinity());
  {
    std::size_t g = 0;
    for (std::size_t k = 0; k < m; ++k) {
      const double x = P.pos[k];
      while (g < Gm.n1 && Gm.regions[g].ra + Gm.regions[g].H < x - 1e-9) ++g;
      double b = -1;
      for (std::size_t q = g; q < Gm.n1 && Gm.regions[q].ra <= x + 1e-9; ++q)
        if (Gm.regions[q].ra + Gm.regions[q].H >= x - 1e-9) b = std::max(b, P.regions[q].B);
      if (b >= 0) P.cub[k] = b;
    }
  }
  P.E.assign(m, 0.0);
  P.have.assign(m, 0);
}

// The effect polynomial of single-truck interval g, built on demand.
void build_region(PlacementPrep& P, const PlacementGeom& Gm, const InfluenceLine& il, std::size_t g) {
  PlacementPrep::Region& R = P.regions[g];
  if (R.built) return;
  if (!P.ip_built) {
    P.ip.build(il);
    P.ip_built = true;
  }
  static const double binom[5][5] = {
      {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};
  const double cum_end = P.ip.cum.back();
  const PlacementGeom::Region& G = Gm.regions[g];
  const std::vector<PlacementGeom::Term>& terms = Gm.terms[0];
  for (double& x : R.e) x = 0;
  for (std::size_t ti = 0; ti < terms.size(); ++ti) {
    const PlacementGeom::Term& t = terms[ti];
    const int i = Gm.idx[G.first + ti];
    double al[5] = {0, 0, 0, 0, 0};
    double sh = 0;
    if (i == -1) continue;  // eta = 0, F = 0
    if (i == -2) {
      if (!t.F) continue;
      al[0] = cum_end;
    } else {
      const std::size_t ii = static_cast<std::size_t>(i);
      if (t.F)
        for (int k = 0; k < 5; ++k) al[k] = P.ip.A[ii][static_cast<std::size_t>(k)];
      else
        for (int k = 0; k < 4; ++k) al[k] = P.ip.a[ii][static_cast<std::size_t>(k)];
      sh = Gm.sh[G.first + ti];
    }
    const double shp[5] = {1, sh, sh * sh, sh * sh * sh, sh * sh * sh * sh};
    for (int jj = 0; jj < 5; ++jj)
      for (int k = jj; k < 5; ++k) R.e[jj] += t.w * al[k] * binom[k][jj] * shp[k - jj];
  }
  double up = R.e[0], dn = -R.e[0], mag = std::fabs(R.e[0]), hp = 1;
  for (int j = 1; j < 5; ++j) {
    hp *= G.H;
    up += std::max(R.e[j] * hp, 0.0);
    dn += std::max(-R.e[j] * hp, 0.0);
    mag += std::fabs(R.e[j] * hp);
  }
  R.up = up + 1e-12 * mag;  // round-off of the coefficients themselves
  R.dn = dn + 1e-12 * mag;
  R.built = true;
}

// The DP over candidate positions: at most K trucks, pairwise >= pitch apart,
// maximising the sum of `val`. Returns the best total and the chosen indices.
double best_trucks(const std::vector<double>& pos, const std::vector<double>& val, int K, double pitch,
                   PlacementPrep& P, std::vector<std::size_t>* chosen) {
  const std::size_t m = pos.size();
  std::vector<int>& prev = P.prev;
  prev.assign(m, -1);
  {
    std::size_t j = 0;
    for (std::size_t i = 0; i < m; ++i) {
      while (j < m && pos[j] <= pos[i] - pitch + 1e-12) ++j;
      prev[i] = static_cast<int>(j) - 1;
    }
  }
  const std::size_t KK = static_cast<std::size_t>(K), W = m + 1;
  P.f.assign((KK + 1) * W, 0.0);
  P.take.assign((KK + 1) * W, -1);
  auto F = [&](std::size_t k, std::size_t i) -> double& { return P.f[k * W + i]; };
  auto T = [&](std::size_t k, std::size_t i) -> int& { return P.take[k * W + i]; };
  for (std::size_t k = 1; k <= KK; ++k)
    for (std::size_t i = 0; i < m; ++i) {
      double skip = F(k, i);
      double base = (prev[i] >= 0) ? F(k - 1, static_cast<std::size_t>(prev[i]) + 1) : 0.0;
      double use = val[i] + base;
      if (use > skip) {
        F(k, i + 1) = use;
        T(k, i + 1) = static_cast<int>(i);
      } else {
        F(k, i + 1) = skip;
        T(k, i + 1) = -1;
      }
    }
  const double best = F(KK, m);
  if (chosen) {
    chosen->clear();
    std::size_t k = KK, i = m;
    while (k > 0 && i > 0) {
      if (T(k, i) >= 0) {
        std::size_t idx = static_cast<std::size_t>(T(k, i));
        chosen->push_back(idx);
        i = (prev[idx] >= 0) ? static_cast<std::size_t>(prev[idx]) + 1 : 0;
        --k;
      } else {
        --i;
      }
    }
  }
  return best;
}

}  // namespace

ConvoyPlacement convoy_placement(const InfluenceLine& il, const ConvoyPattern& pat, int max_cam,
                                 double pitch, double step, bool maximise, bool grid_only) {
  require(max_cam >= 1, "SURCH: MAX_CAM doit valoir au moins 1");
  require(pitch > 0 && finite(pitch), "SURCH: pas minimal entre camions invalide");
  ConvoyPlacement out;
  if (il.s.size() < 2) return out;
  const IlPrefix pre(il);
  PlacementPrep* hit = nullptr;
  for (auto it = tl_preps.begin(); it != tl_preps.end(); ++it)
    if (same_key(*it, il, pat, max_cam, pitch, step, grid_only)) {
      tl_preps.splice(tl_preps.end(), tl_preps, it);  // most recent last
      hit = &tl_preps.back();
      break;
    }
  if (!hit) {
    if (tl_preps.size() >= 24) tl_preps.pop_front();
    tl_preps.emplace_back();
    hit = &tl_preps.back();
    prepare(*hit, il, pre, pat, max_cam, pitch, step, grid_only);
  }
  PlacementPrep& P = *hit;
  if (P.pos.empty()) return out;
  const double sg = maximise ? 1.0 : -1.0;
  const int K = max_cam;
  const double G = std::max(maximise ? P.gmax : P.gmin, 0.0);  // no truck adds more
  // The breakpoint candidates in use and the DP on them: a lower bound of the
  // optimum. Lazily (cubic line): by decreasing bound, until no candidate
  // left can take part in a better placement (its bound plus K - 1 trucks
  // of at most G each).
  std::vector<double>& cpos = P.allpos;
  std::vector<double>& val = P.Eextra;
  double best0 = 0;
  if (!P.lazy) {
    cpos = P.pos;
    val.resize(P.E.size());
    for (std::size_t i = 0; i < P.E.size(); ++i) val[i] = sg * P.E[i];
    best0 = best_trucks(cpos, val, K, pitch, P, nullptr);
  } else {
    const std::size_t m = P.pos.size();
    auto effect = [&](std::size_t k) {
      if (!P.have[k]) {
        P.E[k] = pattern_effect(il, pre, pat, P.pos[k]);
        P.have[k] = 1;
      }
      return sg * P.E[k];
    };
    auto dp_on = [&](const std::vector<std::size_t>& ks) {  // ks ascending
      cpos.clear();
      val.clear();
      for (std::size_t k : ks) {
        cpos.push_back(P.pos[k]);
        val.push_back(effect(k));
      }
      return best_trucks(cpos, val, K, pitch, P, nullptr);
    };
    // the best-bounded few give a first optimum, then one pass keeps every
    // candidate whose bound can still beat it (an optimum only grows)
    std::vector<std::size_t> top(m);
    for (std::size_t k = 0; k < m; ++k) top[k] = k;
    const std::size_t nt = std::min<std::size_t>(m, 64);
    std::nth_element(top.begin(), top.begin() + static_cast<long>(nt - 1), top.end(),
                     [&](std::size_t a, std::size_t b) { return P.cub[a] > P.cub[b]; });
    top.resize(nt);
    std::sort(top.begin(), top.end());
    best0 = dp_on(top);
    std::vector<std::size_t> used;
    for (std::size_t k = 0; k < m; ++k)
      if (P.cub[k] + (K - 1) * G > best0) used.push_back(k);
    // keep the first ones too (their effects are known; harmless)
    std::vector<std::size_t> all;
    std::merge(used.begin(), used.end(), top.begin(), top.end(), std::back_inserter(all));
    all.erase(std::unique(all.begin(), all.end()), all.end());
    best0 = dp_on(all);
  }
  // the stationary points of the intervals that can still beat it (an exact
  // upper bound of the polynomial on [0, H]: the positive parts of its terms)
  std::vector<double>& extra = P.extra;
  extra.clear();
  static const double binom[5][5] = {
      {1, 0, 0, 0, 0}, {1, 1, 0, 0, 0}, {1, 2, 1, 0, 0}, {1, 3, 3, 1, 0}, {1, 4, 6, 4, 1}};
  auto roots_of = [&](const double* e, double ra, double H, int c) {
    double t[3];
    const int nr = cubic_roots01(e[1], 2 * e[2] * H, 3 * e[3] * H * H, 4 * e[4] * H * H * H, t);
    for (int q = 0; q < nr; ++q)
      for (int j = 0; j < c; ++j) extra.push_back(ra + t[q] * H + j * pitch);
  };
  const PlacementGeom& Gm = *P.geom;
  for (std::size_t g = 0; g < P.regions.size(); ++g) {  // one truck
    PlacementPrep::Region& R = P.regions[g];
    if (R.B <= 0 || R.B + (K - 1) * G <= best0) continue;  // cheap bound first
    build_region(P, Gm, il, g);
    const double up = maximise ? R.up : R.dn;
    if (up <= 0) continue;                       // a truck there adds nothing
    if (up + (K - 1) * G <= best0) continue;     // cannot beat the known optimum
    roots_of(R.e, R.ra, R.H, 1);
  }
  for (std::size_t g = Gm.n1; g < Gm.regions.size(); ++g) {  // files in contact
    const PlacementGeom::Region& F = Gm.regions[g];
    double ub = 0;  // the sum of its trucks' cheap bounds
    for (int j = 0; j < F.c; ++j) ub += P.regions[Gm.r1[F.first + static_cast<std::size_t>(j)]].B;
    if (ub <= 0 || ub + (K - F.c) * G <= best0) continue;
    double up = 0;  // the sum of its trucks' interval bounds
    for (int j = 0; j < F.c; ++j) {
      const std::size_t g1 = Gm.r1[F.first + static_cast<std::size_t>(j)];
      build_region(P, Gm, il, g1);
      const PlacementPrep::Region& R = P.regions[g1];
      up += maximise ? R.up : R.dn;
    }
    if (up <= 0 || up + (K - F.c) * G <= best0) continue;
    double e[5] = {0, 0, 0, 0, 0};  // sum of the trucks' polynomials, shifted
    for (int j = 0; j < F.c; ++j) {
      const std::size_t k = F.first + static_cast<std::size_t>(j);
      const PlacementPrep::Region& R = P.regions[Gm.r1[k]];
      const double o = Gm.off[k], op[5] = {1, o, o * o, o * o * o, o * o * o * o};
      for (int jj = 0; jj < 5; ++jj)
        for (int kk = jj; kk < 5; ++kk) e[jj] += R.e[kk] * binom[kk][jj] * op[kk - jj];
    }
    roots_of(e, F.ra, F.H, F.c);
    // a file in contact may also rest on a breakpoint: its partners
    for (int j = 1; j < F.c; ++j) {
      extra.push_back(F.ra + j * pitch);
      extra.push_back(F.ra + F.H + j * pitch);
    }
  }
  std::vector<double>* posv = &cpos;
  std::vector<double>* valv = &val;
  std::vector<double> mp, mv;
  if (!extra.empty()) {
    // the candidates in use and the stationary points, merged in order
    std::sort(extra.begin(), extra.end());
    mp.reserve(cpos.size() + extra.size());
    mv.reserve(cpos.size() + extra.size());
    std::size_t a = 0, b = 0;
    while (a < cpos.size() || b < extra.size()) {
      const bool from_c = b == extra.size() || (a < cpos.size() && cpos[a] <= extra[b]);
      const double x = from_c ? cpos[a] : extra[b];
      if (mp.empty() || std::fabs(x - mp.back()) >= 1e-10) {
        mp.push_back(x);
        mv.push_back(from_c ? val[a] : sg * pattern_effect(il, pre, pat, x));
      }
      (from_c ? a : b)++;
    }
    posv = &mp;
    valv = &mv;
  }
  std::vector<std::size_t> chosen;
  const double best = best_trucks(*posv, *valv, K, pitch, P, &chosen);
  out.value = maximise ? best : -best;
  for (std::size_t idx : chosen) out.positions.push_back((*posv)[idx]);
  std::sort(out.positions.begin(), out.positions.end());
  // Un camion d'effet nul n'apporte rien : on ne conserve que les utiles.
  std::vector<double> keep;
  for (double r : out.positions)
    if (std::fabs(pattern_effect(il, pre, pat, r)) > 0) keep.push_back(r);
  out.positions = keep;
  return out;
}

// ===========================================================================
// 8. FACTEUR TRANSVERSAL (p.125 « cumul de l'effet de chaque file »)
// ===========================================================================

namespace {

struct TransChoice {
  double factor = 1.0;  // sans dimension (convois) ou en m (charges réparties)
  int n_lanes = 1;
  double phys_width = 0;  // largeur physique chargée [m]
};

// Coefficients de répartition transversale des voies chargées, vus par la
// section d'étude, triés par valeur décroissante (les voies les plus efficaces
// sont chargées d'abord).
std::vector<double> lane_coefs(const Tablier& t, const std::vector<DeckLane>& lanes,
                               Id study_bar, double impact_half, bool maximise) {
  std::vector<double> k;
  for (const DeckLane& l : lanes) {
    double hw = impact_half > 0 ? std::min(impact_half, l.width / 2) : l.width / 2;
    k.push_back(t.rep_coef(study_bar, l.y, hw));
  }
  std::sort(k.begin(), k.end(),
            [&](double a, double b) { return maximise ? a > b : a < b; });
  return k;
}

// Convois : b(n) * somme des n meilleurs coefficients, n balayé de 1 à n_max.
TransChoice convoy_trans(const Tablier& t, const ConvoyDef& c, const std::vector<int>& zone_ids,
                         Id study_bar, bool maximise) {
  TransChoice out;
  std::vector<DeckLane> lanes = t.lanes(zone_ids, false);
  if (lanes.empty()) {
    // p.144 : sans ZONE_TRANS, ST1 se limite au positionnement longitudinal
    // d'UNE file, sans pondération transversale.
    out.factor = 1.0;
    out.n_lanes = 1;
    out.phys_width = c.largeur;
    return out;
  }
  double impact_half = c.largeur > 0 ? c.largeur / 2 : 0;
  // p.148 : « les charges roulantes, affectées à des zones spécifiques, sont
  // placées sur toute la largeur des zones attribuées » — des files de la
  // largeur du convoi, accolées depuis un bord de la zone (autant de files
  // que de voies dans la zone). Avec REP_TRANS, le bord retenu est le plus
  // défavorable ; sans, toutes les files ont le coefficient 1.
  auto file_coefs = [&](bool from_left) {
    std::vector<double> k;
    const double half = t.width() / 2;
    double zl = -half;
    for (std::size_t z = 0; z < t.zones.size(); ++z) {
      const double zw = t.zones[z].width;
      const int id = static_cast<int>(z) + 1;
      int nl = 0;
      for (const DeckLane& l : lanes)
        if (l.zone == id) ++nl;
      const double fw = c.largeur > 0 ? std::min(c.largeur, zw) : zw / std::max(nl, 1);
      for (int j = 0; j < nl; ++j) {
        if ((j + 1) * fw > zw + 1e-9) break;
        double y = from_left ? zl + (j + 0.5) * fw : zl + zw - (j + 0.5) * fw;
        k.push_back(t.rep_coef(study_bar, y, fw / 2));
      }
      zl += zw;
    }
    std::sort(k.begin(), k.end(), [&](double a, double b) { return maximise ? a > b : a < b; });
    return k;
  };
  std::vector<double> k = lane_coefs(t, lanes, study_bar, impact_half, maximise);
  if (c.largeur > 0) {
    std::vector<double> kl = file_coefs(true), kr = file_coefs(false);
    if (!kl.empty()) {
      double sl = 0, sr = 0;
      for (double v : kl) sl += v;
      for (double v : kr) sr += v;
      k = (maximise ? sl >= sr : sl <= sr) ? kl : kr;
    }
  }
  int n_max = static_cast<int>(k.size());
  if (c.max_file > 0) n_max = std::min(n_max, c.max_file);
  double sum = 0;
  bool first = true;
  for (int n = 1; n <= n_max; ++n) {
    sum += k[static_cast<std::size_t>(n) - 1];
    double b = 1.0;
    switch (c.trans) {
      case TransLaw::Fasc61Bc: b = st1_bc(t.fasc61_class(), n); break;
      case TransLaw::Fasc61Bt: b = st1_bt(t.fasc61_class(), n); break;
      case TransLaw::Table:
        b = c.coef_trans.empty()
                ? 1.0
                : c.coef_trans[std::min(static_cast<std::size_t>(n), c.coef_trans.size()) - 1];
        break;
      case TransLaw::None: b = 1.0; break;
    }
    if (b == 0) continue;  // combinaison interdite par la table (0 dans le manuel)
    double f = b * sum;
    if (first || (maximise ? f > out.factor : f < out.factor)) {
      out.factor = f;
      out.n_lanes = n;
      double w = 0;
      for (int j = 0; j < n; ++j) w += lanes[static_cast<std::size_t>(j)].width;
      out.phys_width = w;
      first = false;
    }
  }
  if (first) fail("SURCH: aucune combinaison de voies autorisee par la table (classe/voies)");
  return out;
}

}  // namespace

// ===========================================================================
// 9. ÉVALUATION D'UNE SURCHARGE SUR UNE LIGNE D'INFLUENCE
// ===========================================================================

namespace {

struct Evaluated {
  double value = 0;
  SurchPlacement pl;
};

// --- convois --------------------------------------------------------------

Evaluated eval_convoy_dir(const InfluenceLine& il, const Tablier& t, const SurchDef& d,
                          const ConvoyDef& c, const Vec3& dir, Id study_bar, bool maximise,
                          bool tmax) {
  Evaluated ev;
  TransChoice tr = convoy_trans(t, c, d.zones, study_bar, tmax);
  double scale = tr.factor * d.pond;
  // CONV DIR coefficients: the multiplier of the load along its direction
  if (c.longitudinal && c.dir_coef.x() != 0) scale *= c.dir_coef.x();
  else if (!c.longitudinal && c.dir_coef.norm() > 0) scale *= c.dir_coef.norm();
  ConvoyPattern pat = convoy_pattern(c, d.tonne_force, t.ldif);
  ConvoyPlacement cp = convoy_placement(il, pat, c.max_cam, c.longueur, d.pas_conv,
                                        maximise ? scale >= 0 : scale < 0, d.pas_conv_grid_only);
  ev.value = scale * cp.value;
  for (auto& [a, P] : pat.axles) P *= scale;
  for (auto& seg : pat.segments) seg[2] *= scale;
  PlacedGroup g;
  g.pattern = pat;
  g.positions = cp.positions;
  if (!g.positions.empty()) ev.pl.groups.push_back(g);
  ev.pl.dir = dir;
  ev.pl.dir_local = c.longitudinal;
  ev.pl.trans_factor = tr.factor;
  ev.pl.n_lanes = tr.n_lanes;
  return ev;
}

Evaluated eval_convoy(const InfluenceLine& il, const Tablier& t, const SurchDef& d,
                      const ConvoyDef& c, const Vec3& dir, Id study_bar, bool maximise,
                      bool tmax) {
  require(d.sens >= 1 && d.sens <= 3, "SURCH SENS: 1 (croissant), 2 (decroissant) ou 3");
  // ST1 SENS (p.146): 1 croissant = the convoy MIRRORED from its definition
  // order, 2 decroissant = as defined (with the opposite reading, ST1's
  // SENS 1 and SENS 2 results would swap)
  Evaluated best;
  bool first = true;
  if (d.sens == 2 || d.sens == 3) {
    best = eval_convoy_dir(il, t, d, c, dir, study_bar, maximise, tmax);
    first = false;
  }
  if (d.sens == 1 || d.sens == 3) {
    Evaluated m = eval_convoy_dir(il, t, d, c.mirrored(), dir, study_bar, maximise, tmax);
    m.pl.mirrored = true;
    if (first || (maximise ? m.value > best.value : m.value < best.value)) best = m;
  }
  return best;
}

// --- charges réparties entre zéros ---------------------------------------

// gain = facteur transversal en m (largeur pondérée) ; q(L) en N/m².
Evaluated eval_distributed(const InfluenceLine& il, const SurchDef& d, double gain,
                           const std::function<double(double)>& q, bool maximise,
                           const Vec3& dir, bool local, int n_lanes, double phys_width) {
  Evaluated ev;
  ev.pl.dir = dir;
  ev.pl.dir_local = local;
  ev.pl.trans_factor = gain;
  ev.pl.n_lanes = n_lanes;
  std::vector<SignZone> zones = il_sign_zones(il, d.seuil);
  LoadedChoice ch = pick_zones(zones, gain, q, maximise);
  if (ch.zones.empty()) return ev;
  double dens = q(ch.length);
  double line = dens * gain;
  ev.value = line * ch.integral;
  for (const SignZone& z : ch.zones) ev.pl.udl.push_back({z.s0, z.s1, line});
  ev.pl.loaded_length = ch.length;
  ev.pl.loaded_surface = ch.length * phys_width;
  ev.pl.density = dens;
  return ev;
}

}  // namespace

// ===========================================================================
// 10. RECONSTRUCTION DU CAS DE CHARGE (concomitance, p.137)
// ===========================================================================

namespace {

LoadCase placement_case(const Model& m, const std::vector<DeckSpan>& sp,
                        const SurchPlacement& pl, const std::string& name) {
  LoadCase lc;
  lc.id = 1;
  lc.name = name;
  // The couple of an eccentric load (GRILL / SPATIALE): point couples; a
  // distributed load's couple is lumped on 10 points per segment (the
  // concomitant replay of an eccentric UDL is thus discretized).
  auto add_couple = [&](double s, double Mt) {
    const DeckSpan* d = span_at(sp, s);
    if (d == nullptr || Mt == 0) return;
    Vec3 mv = ecc_torque_dir(m, *d, pl.dir);
    if (mv.isZero()) return;
    double x = std::clamp(d->x_of(s), 0.0, d->L);
    lc.loads.push_back(Load::point(d->bar, x, Vec3::Zero(), Mt * mv));
  };
  auto add_udl_couple = [&](double s0, double s1, double q, double e) {
    if (e == 0 || s1 <= s0) return;
    const int n = 10;
    for (int i = 0; i < n; ++i) add_couple(s0 + (s1 - s0) * (i + 0.5) / n, q * (s1 - s0) / n * e);
  };
  for (const PlacedGroup& g : pl.groups)
    for (double r : g.positions) {
      for (const auto& [a, P] : g.pattern.axles) {
        add_point(m, sp, r + a, P, pl.dir, pl.dir_local, lc.loads);
        if (g.ecc != 0 && !pl.dir_local) add_couple(r + a, P * g.ecc);
      }
      for (const auto& seg : g.pattern.segments)
        if (seg[1] > 0) {
          add_udl(m, sp, r + seg[0], r + seg[0] + seg[1], seg[2] / seg[1], pl.dir, pl.dir_local,
                  lc.loads);
          if (g.ecc != 0 && !pl.dir_local)
            add_udl_couple(r + seg[0], r + seg[0] + seg[1], seg[2] / seg[1], g.ecc);
        }
    }
  for (std::size_t i = 0; i < pl.udl.size(); ++i) {
    const auto& u = pl.udl[i];
    add_udl(m, sp, u[0], u[1], u[2], pl.dir, pl.dir_local, lc.loads);
    if (i < pl.udl_ecc.size() && !pl.dir_local) add_udl_couple(u[0], u[1], u[2], pl.udl_ecc[i]);
  }
  return lc;
}

}  // namespace

// ===========================================================================
// 11. RUNNER
// ===========================================================================

void SurchDef::validate() const {
  require(finite(pond), "SURCH POND: coefficient non fini");
  require(pas_conv > 0 && finite(pas_conv), "SURCH PAS_CONV: le pas doit etre > 0");
  require(sens >= 1 && sens <= 3, "SURCH SENS: 1, 2 ou 3");
  require(seuil >= 0 && seuil < 1, "SURCH SEUIL: doit etre dans [0, 1[");
  require(finite(tonne_force) && tonne_force > 0, "SURCH: facteur tonne->force invalide");
  require(priv_comp >= 0 && priv_comp <= 5, "SURCH: composante privilegiee 0..5");
  require(lgexclu >= 0 && finite(lgexclu), "SURCH LGEXCLU: longueur non finie ou negative");
  require(finite(pondvs), "SURCH PONDVS: coefficient non fini");
  require(lm3_lm1 >= 0 && lm3_lm1 <= 3,
          "SURCH LM3: accompagnement 0 (NUL), 1 (FREQ), 2 (CARA) ou 3 (COMB)");
  if (load == SurchLoad::Alg) alg.validate();
  if ((load == SurchLoad::Convoy || load == SurchLoad::Lm3) && convoy.empty())
    // LONG/LARG par défaut sont complétés par surch_convoy() avant validate().
    require(!conv.axles.empty(),
            "SURCH: convoi non defini (mot-cle du catalogue ou CONV generalise)");
}

namespace {

// Convoi effectif d'une surcharge roulante : catalogue ou CONV généralisé.
ConvoyDef surch_convoy(const SurchDef& d) {
  if (!d.convoy.empty()) return st1_convoy(d.convoy);
  ConvoyDef c = d.conv;
  if (c.longueur <= 0) {
    // LONG par défaut = max(xl + impact/2) (p.160).
    double L = 0;
    for (const ConvoyAxle& a : c.axles) L = std::max(L, a.xl + a.impact_l / 2);
    c.longueur = std::max(L, 1e-6);
  }
  if (c.coef_trans.empty() && c.trans == TransLaw::Table) c.trans = TransLaw::None;
  c.validate();
  return c;
}

struct LaneRank {
  double k = 0;
  double width = 0;
  double y = 0;  // lane axis (transverse abscissa, deck axis = 0)
};

// Voies triées par |coefficient| décroissant : la numérotation Eurocode la plus
// défavorable donne le rang 1 à la voie la plus efficace (EN 1991-2 §4.2.4).
std::vector<LaneRank> ranked_lanes(const Tablier& t, const SurchDef& d, Id bar, bool maximise) {
  std::vector<DeckLane> lanes = t.lanes(d.zones, true);
  std::vector<LaneRank> out;
  if (lanes.empty()) {
    LaneRank r;
    r.k = 1.0;
    r.width = 3.0;
    out.push_back(r);
    return out;
  }
  for (const DeckLane& l : lanes) {
    LaneRank r;
    // LM1 CENTRE : coefficient au centre de la voie ; LIBRE : meilleur point.
    if (d.lm1_libre) {
      r.k = t.rep_extreme(bar, l.y - l.width / 2, l.y + l.width / 2, 0.0, maximise);
    } else {
      r.k = t.rep_coef(bar, l.y);
    }
    r.width = l.width;
    r.y = l.y;
    out.push_back(r);
  }
  std::sort(out.begin(), out.end(), [&](const LaneRank& a, const LaneRank& b) {
    return maximise ? a.k > b.k : a.k < b.k;
  });
  return out;
}

// LM1 : TS (2 essieux à 1,20 m) cumulés sur les voies + UDL par voie.
struct Lm1Loads {
  double ts_axle = 0;  // charge d'un essieu du tandem cumulé [N]
  double udl = 0;      // densité linéique cumulée [N/m]
  // UDL of the lanes / residual strips whose REP_TRANS coefficient has the
  // OTHER sign: EN 1991-2 4.3.2(1)b loads the unfavourable parts of the
  // influence SURFACE k(y) IL(x), so they go on the zones of opposite sign
  // of the line (ranks continue after the main family's).
  double udl_other = 0;
  int n_lanes = 0;
  // Per loaded lane, in Eurocode rank order (LM3: the VS lanes are removed
  // from the LM1 inside the exclusion window).
  std::vector<double> ts_rank, udl_rank;
  std::vector<LaneRank> lanes;
};

Lm1Loads lm1_loads(const Tablier& t, const SurchDef& d, Id bar, bool maximise, double psi_ts,
                   double psi_udl) {
  Lm1Loads o;
  std::vector<LaneRank> lanes = ranked_lanes(t, d, bar, maximise);
  std::vector<std::size_t> other;  // lanes of the opposite coefficient sign
  const EurocodeParams& e = d.ec;
  double res_width = t.loaded_width(d.zones);
  // A lane (or residual strip) is loaded only where it is unfavourable
  // (EN 1991-2 §4.2.4-4.3.2): with REP_TRANS coefficients of mixed signs, a
  // lane of the wrong sign for the sought transverse extreme stays empty.
  auto useful = [&](double k) { return maximise ? k > 0 : k < 0; };
  int rank = 0;
  for (std::size_t i = 0; i < lanes.size(); ++i) {
    res_width -= lanes[i].width;
    if (!t.rep_trans.empty() && !useful(lanes[i].k)) {
      if (lanes[i].k != 0) other.push_back(i);
      continue;
    }
    ++rank;
    double ts = psi_ts * e.aQ(rank) * e.Q(rank) * lanes[i].k;
    double ud = psi_udl * e.aq(rank) * e.q(rank) * lanes[i].width * lanes[i].k;
    o.ts_axle += ts;
    o.udl += ud;
    o.ts_rank.push_back(ts);
    o.udl_rank.push_back(ud);
    o.lanes.push_back(lanes[i]);
  }
  // Residual areas: the strip of each loaded zone left of its lanes, with
  // ITS OWN transverse coefficient (it used to borrow lane 1's).
  if (res_width > 1e-9 && !t.zones.empty()) {
    const double half = t.width() / 2;
    std::vector<double> left(t.zones.size(), 0.0);
    double acc = -half;
    for (std::size_t z = 0; z < t.zones.size(); ++z) {
      left[z] = acc;
      acc += t.zones[z].width;
    }
    std::vector<DeckLane> dl = t.lanes(d.zones, true);
    for (int id : t.resolved_zones(d.zones)) {
      const TransZone& z = t.zones[static_cast<std::size_t>(id) - 1];
      double used = 0;
      for (const DeckLane& l : dl)
        if (l.zone == id) used += l.width;
      double w = z.width - used;
      if (w <= 1e-9) continue;
      double yc = left[static_cast<std::size_t>(id) - 1] + used + w / 2;
      double k = lanes.empty() ? 1.0 : t.rep_coef(bar, yc, w / 2);
      if (!t.rep_trans.empty() && !useful(k)) {
        o.udl_other += psi_udl * e.aqr * e.qrk * w * k;
        continue;
      }
      o.udl += psi_udl * e.aqr * e.qrk * w * k;
    }
  }
  int rk = rank;
  for (std::size_t i : other) {
    ++rk;
    o.udl_other += psi_udl * e.aq(rk) * e.q(rk) * lanes[i].width * lanes[i].k;
  }
  o.n_lanes = rank;
  return o;
}

}  // namespace


namespace {

// Contexte partagé par les deux points d'entrée : chemin de roulement,
// discrétisation de la ligne d'influence et champ de directions de charge.
struct SurchCore {
  std::vector<DeckSpan> spans;
  LoadPath path;
  std::vector<Vec3> dirs;     // direction unitaire (globale) station par station
  Vec3 gdir = Vec3(0, 0, -1); // direction retenue (ou (+/-1,0,0) en repère local)
  bool local = false;         // charges suivant ox local des barres
  double dir_scale = 1.0;     // norme de la direction ALG, repliée dans la densité
  // Excentrement des voies (GRILL / SPATIALE) : couple unitaire n_gauche x d
  // par station, n_gauche = (-d) x tangente.
  bool ecc = false;
  std::vector<Vec3> tdirs;
};


// Chemin de roulement aux discontinuités exactes. Une ligne d'influence
// d'effort SAUTE au droit de sa section (effort tranchant : +1/2 -> -1/2) et
// une jonction de barres porte deux limites (fin de la barre amont, début de
// la barre aval). Le chemin uniforme n'y garde qu'un point et lisse le saut
// sur un pas, ce qui fausse les zones de signe et la longueur chargée. On
// insère donc, à chaque jonction, le début de la barre aval, et à chaque
// section d'étude intérieure à une barre du tablier, ses deux limites
// (x et x + delta), en ordre de parcours.
void refine_path(const Model& model, const std::vector<DeckSpan>& spans,
                 const std::vector<SurchSection>& sections, LoadPath& path) {
  (void)model;
  // A bar's origin x = 0 lying at an end of the path is pulled delta inside
  // the bar: a unit load at a == x counts on the origin side (influence.cpp),
  // so a section at that origin would otherwise see an isolated spike — the
  // load "before" a section that nothing precedes. ST1 likewise excludes a
  // load sitting on the end node from the bar-end forces (free cantilever
  // tip, abutment).
  auto pull_in = [&](PathStation& st, double sign) {
    for (const DeckSpan& d : spans)
      if (d.bar == st.bar && st.x_local == 0.0) {
        const double delta = 1e-7 * d.L;
        st.x_local = delta;
        st.s += sign * delta;
      }
  };
  std::vector<PathStation> add;
  for (std::size_t i = 1; i < spans.size(); ++i) {
    const DeckSpan& d = spans[i];
    add.push_back({d.bar, d.flip ? d.L : 0.0, d.s0});
  }
  for (const SurchSection& sec : sections)
    for (const DeckSpan& d : spans) {
      if (d.bar != sec.bar) continue;
      const double delta = 1e-7 * d.L;
      // A section at the bar origin counts a load at a = 0 on its origin
      // side; the other limit is a load just inside the bar, at x = delta —
      // without that station the jump would be smeared over a grid step
      // (stay anchorages, Exemple3). At x = L the jump sits between the bar
      // end and the next bar's junction station, already on the path.
      if (sec.x <= 1e-9 * d.L) {
        add.push_back({d.bar, delta, d.flip ? d.s0 + d.L - delta : d.s0 + delta});
        continue;
      }
      if (sec.x >= d.L * (1 - 1e-9)) continue;
      double s = d.flip ? d.s0 + d.L - sec.x : d.s0 + sec.x;
      // a <= x est compté du côté de l'origine de la barre : (b, x) est la
      // limite côté x décroissant, (b, x + delta) la limite côté x croissant.
      if (d.flip) {
        add.push_back({d.bar, sec.x + delta, s - delta});
        add.push_back({d.bar, sec.x, s});
      } else {
        add.push_back({d.bar, sec.x, s});
        add.push_back({d.bar, sec.x + delta, s + delta});
      }
    }
  if (!path.stations.empty()) {
    pull_in(path.stations.front(), 1.0);
    pull_in(path.stations.back(), -1.0);
  }
  if (add.empty()) return;
  path.stations.insert(path.stations.end(), add.begin(), add.end());
  // Ordre de parcours ; à abscisse égale, la barre amont d'abord (ordre des
  // travées), puis l'abscisse locale dans le sens de parcours.
  std::map<Id, std::size_t> rank;
  std::map<Id, bool> flipped;
  for (std::size_t i = 0; i < spans.size(); ++i) {
    rank[spans[i].bar] = i;
    flipped[spans[i].bar] = spans[i].flip;
  }
  std::stable_sort(path.stations.begin(), path.stations.end(),
                   [&](const PathStation& a, const PathStation& b) {
                     if (a.s != b.s) return a.s < b.s;
                     if (a.bar != b.bar) return rank[a.bar] < rank[b.bar];
                     return flipped[a.bar] ? a.x_local > b.x_local : a.x_local < b.x_local;
                   });
  path.stations.erase(std::unique(path.stations.begin(), path.stations.end(),
                                  [](const PathStation& a, const PathStation& b) {
                                    return a.bar == b.bar && a.x_local == b.x_local;
                                  }),
                      path.stations.end());
}

SurchCore build_core(const Model& model, const Tablier& t, const SurchDef& d,
                     const SurchOptions& opt, const std::vector<SurchSection>& sections = {}) {
  SurchCore c;
  require(opt.down.allFinite() && std::fabs(opt.down.norm() - 1.0) < 1e-9,
          "SURCH: direction de pesanteur non unitaire");
  c.gdir = opt.down;
  c.spans = deck_spans(model, t.bars);
  double step = opt.il_step > 0 ? opt.il_step : t.pas_ligne;
  require(step > 0 && finite(step), "SURCH: pas de ligne d'influence invalide");
  double Lmax = 0;
  for (const DeckSpan& s : c.spans) Lmax = std::max(Lmax, s.L);
  int nper = std::max(3, static_cast<int>(std::ceil(Lmax / step)) + 1);
  c.path = build_path(model, t.bars, nper);
  refine_path(model, c.spans, sections.empty() ? opt.sections : sections, c.path);

  c.local = (d.load == SurchLoad::AlFrein || d.load == SurchLoad::AlPlancherFrein);
  if (c.local) c.gdir = Vec3(1, 0, 0);
  // ST1 2.24 applies AL PLANCHER FREIN against ox (e.g. a support reaction
  // Fx of +5.84 where AL FREIN gives -12.86), although the
  // manual (p.143) states both along ox: reproduced as ST1 runs it.
  if (d.load == SurchLoad::AlPlancherFrein) c.gdir = Vec3(-1, 0, 0);
  if (d.load == SurchLoad::Alg) {
    // no DIR: along the self-weight (p.160), the opt.down set above
    double nrm = d.alg.dir.norm();
    require(finite(nrm), "ALG DIR: direction non finie");
    if (nrm > 0) {
      c.dir_scale = nrm;
      c.gdir = d.alg.dir / nrm;
      c.local = d.alg.dir_local;
    }
    if (c.local) {
      require(std::fabs(std::fabs(c.gdir(0)) - 1.0) < 1e-9,
              "ALG DIR LOC: seule la composante FX (axe ox local) est portee par le moteur de "
              "lignes d'influence; utiliser DIR GLO, ou REP_TRANS pour un excentrement");
      c.gdir = Vec3(c.gdir(0) > 0 ? 1.0 : -1.0, 0, 0);
    }
  }
  if (d.load == SurchLoad::Convoy && !surch_convoy(d).longitudinal &&
      surch_convoy(d).dir_coef.norm() > 0)
    c.gdir = surch_convoy(d).dir_coef.normalized();  // CONV DIR GLO
  if (d.load == SurchLoad::Convoy && surch_convoy(d).longitudinal) {
    c.local = true;
    c.gdir = Vec3(1, 0, 0);
  }

  c.dirs.assign(c.path.stations.size(), c.gdir);
  if (opt.eccentric_lanes && !c.local && !t.zones.empty()) {
    c.ecc = true;
    c.tdirs.assign(c.path.stations.size(), Vec3::Zero());
    for (std::size_t k = 0; k < c.path.stations.size(); ++k) {
      const DeckSpan* sp = span_at(c.spans, c.path.stations[k].s);
      require(sp != nullptr, "SURCH: station hors du chemin de roulement");
      c.tdirs[k] = ecc_torque_dir(model, *sp, c.gdir);
    }
  }
  if (c.local) {
    for (std::size_t k = 0; k < c.path.stations.size(); ++k) {
      const DeckSpan* sp = span_at(c.spans, c.path.stations[k].s);
      require(sp != nullptr, "SURCH: station hors du chemin de roulement");
      c.dirs[k] = c.gdir(0) * path_tangent(model, *sp);
    }
  }
  return c;
}

// Évalue les deux extrêmes d'une réponse.
// `maximise` choisit la branche (max / min) ; `tmax` choisit le sens de la
// sélection TRANSVERSALE (voies les plus / moins efficaces). Les deux sont
// indépendants : avec des coefficients REP_TRANS de signes mêlés, le minimum
// d'un effet peut être atteint avec le facteur transversal MAXIMAL. On évalue
// donc les deux coins et on retient le meilleur (inutile sans REP_TRANS).
// GRILL / SPATIALE — étude transversale voie par voie (p.128-129) : chaque
// voie (charges réparties à a1/a2), zone (trottoirs) ou file de camions
// (charges roulantes, accolées depuis un bord de zone, p.148) a sa propre
// ligne d'influence eta_F + y eta_T, y étant son excentrement (positif à
// gauche du parcours). Chaque élément est placé indépendamment au plus
// défavorable, puis les éléments sont cumulés avec le coefficient
// réglementaire du nombre d'éléments chargés (a1, bc, bt, COEF TRANS), les
// plus défavorables d'abord ; les zones de trottoir se cumulent sans
// coefficient. Renvoie false pour les charges non concernées (freinage,
// Eurocode : étude scalaire inchangée).
bool eval_eccentric(const InfluenceLine& ilF, const InfluenceLine& ilT, const Tablier& t,
                    const SurchDef& d, const SurchCore& core, Id study_bar, bool maximise,
                    Evaluated& out) {
  // Also in PLANE with REP_TRANS (ilT = 0): each lane / file then carries its
  // own coefficient and is placed on its own unfavourable zones — the
  // influence SURFACE k(y) IL(x) of EN 1991-2 4.3.2(1)b / Fasc.61, with no
  // torque couple (`ecc` 0 below).
  if ((!core.ecc && t.rep_trans.empty()) || core.local) return false;
  const double ecc_on = core.ecc ? 1.0 : 0.0;
  auto lane_il = [&](double y_left) {
    InfluenceLine il = ilF;
    require(ilT.s.size() == il.s.size(), "SURCH: lignes flexion / torsion non alignees");
    for (std::size_t k = 0; k < il.eta.size(); ++k) il.eta[k] += y_left * ilT.eta[k];
    if (il.cubic() && ilT.cubic())  // the slopes combine as the ordinates
      for (std::size_t k = 0; k < il.eta.size(); ++k) {
        il.dl[k] += y_left * ilT.dl[k];
        il.dr[k] += y_left * ilT.dr[k];
      }
    else {
      il.dl.clear();
      il.dr.clear();
    }
    return il;
  };
  auto better = [&](double a, double b) { return maximise ? a > b : a < b; };
  const int cls = t.fasc61_class();
  const double tf = d.tonne_force;

  // Combine per-element results: `coef(n)` (0 = forbidden), elements sorted
  // most unfavourable first; `scale_placement` rescales an element's loads.
  struct Elem {
    double value = 0;
    Evaluated ev;
  };
  auto combine = [&](std::vector<Elem>& el, int n_max, const std::function<double(int)>& coef,
                     bool scale_udl) {
    std::sort(el.begin(), el.end(), [&](const Elem& a, const Elem& b) { return better(a.value, b.value); });
    n_max = std::min<int>(n_max, static_cast<int>(el.size()));
    bool found = false;
    double best = 0;
    int bn = 0;
    double sum = 0;
    for (int n = 1; n <= n_max; ++n) {
      sum += el[static_cast<std::size_t>(n) - 1].value;
      double c = coef(n);
      if (c == 0) continue;
      double tot = c * sum;
      if (!found || better(tot, best)) {
        best = tot;
        bn = n;
        found = true;
      }
    }
    out = Evaluated{};
    out.pl.dir = core.gdir;
    if (!found || !better(best, 0.0)) return;
    out.value = best;
    out.pl.n_lanes = bn;
    const double c = coef(bn);
    for (int i = 0; i < bn; ++i) {
      const Evaluated& e = el[static_cast<std::size_t>(i)].ev;
      for (PlacedGroup g : e.pl.groups) {
        for (auto& [a, P] : g.pattern.axles) P *= c;
        for (auto& seg : g.pattern.segments) seg[2] *= c;
        out.pl.groups.push_back(g);
      }
      for (std::size_t k = 0; k < e.pl.udl.size(); ++k) {
        auto u = e.pl.udl[k];
        if (scale_udl) u[2] *= c;
        out.pl.udl.push_back(u);
        out.pl.udl_ecc.push_back(k < e.pl.udl_ecc.size() ? e.pl.udl_ecc[k] : 0.0);
      }
    }
  };

  switch (d.load) {
    case SurchLoad::Convoy: {
      ConvoyDef c = surch_convoy(d);
      if (c.longitudinal) return false;
      std::vector<DeckLane> lanes = t.lanes(d.zones, false);
      if (lanes.empty()) return false;
      ConvoyPattern pat = convoy_pattern(c, tf, t.ldif);
      // One file of width `fw` centred at `y`, placed longitudinally on its
      // own lane line (files are independent longitudinally).
      auto file_elem = [&](double y, double fw) {
        double k = t.rep_coef(study_bar, y, fw / 2);
        double scale = k * d.pond;
        InfluenceLine il = lane_il(-y);
        Elem e;
        bool first = true;
        for (int dirn = 0; dirn < 2; ++dirn) {
          if (dirn == 0 && d.sens == 1) continue;  // SENS: see eval_convoy
          if (dirn == 1 && d.sens == 2) continue;
          ConvoyPattern p = dirn == 0 ? pat : pat.mirrored();
          ConvoyPlacement cp = convoy_placement(il, p, c.max_cam, c.longueur, d.pas_conv,
                                                maximise ? scale >= 0 : scale < 0,
                                                d.pas_conv_grid_only);
          double v = scale * cp.value;
          if (first || better(v, e.value)) {
            e.value = v;
            e.ev = Evaluated{};
            PlacedGroup g;
            g.pattern = p;
            for (auto& [a, P] : g.pattern.axles) P *= scale;
            for (auto& sg : g.pattern.segments) sg[2] *= scale;
            g.positions = cp.positions;
            g.ecc = -y * ecc_on;
            if (!g.positions.empty()) e.ev.pl.groups.push_back(g);
            first = false;
          }
        }
        return e;
      };
      // Within a zone, the files only need to fit side by side without
      // overlapping (p.148: "positionnés de la façon la plus défavorable"):
      // the extreme arrangements pack `kl` files from the left edge and `kr`
      // from the right edge, kl + kr <= files that fit. Two Bc files in a
      // 7 m carriageway can thus run along opposite kerbs — the checkerboard
      // ("damier") loading that governs torsion. Zones then combine by a
      // knapsack on the file count; the transverse coefficient applies to
      // the total.
      struct Pick {
        bool ok = false;
        double sum = 0;
        std::vector<const Elem*> el;
      };
      std::vector<std::vector<Elem>> store;  // per zone: left slots, then right slots
      std::vector<Pick> dp(1);
      dp[0].ok = true;
      double zl = -t.width() / 2;
      store.reserve(t.zones.size());
      for (std::size_t z = 0; z < t.zones.size(); ++z) {
        const double zw = t.zones[z].width;
        const int id = static_cast<int>(z) + 1;
        int nl = 0;
        for (const DeckLane& l : lanes)
          if (l.zone == id) ++nl;
        const double fw = c.largeur > 0 ? std::min(c.largeur, zw) : zw / std::max(nl, 1);
        int nfit = 0;
        while (nfit < nl && (nfit + 1) * fw <= zw + 1e-9) ++nfit;
        store.emplace_back();
        std::vector<Elem>& slots = store.back();
        if (c.axis_ecc >= 0 && nfit > 0) {
          // one file along the zone axis, +/- axis_ecc within the room left
          nfit = 1;
          const double e = std::min(c.axis_ecc, std::max(0.0, (zw - fw) / 2));
          slots.push_back(file_elem(zl + zw / 2 - e, fw));
          slots.push_back(file_elem(zl + zw / 2 + e, fw));
        } else {
          for (int j = 0; j < nfit; ++j) slots.push_back(file_elem(zl + (j + 0.5) * fw, fw));
          for (int j = 0; j < nfit; ++j) slots.push_back(file_elem(zl + zw - (j + 0.5) * fw, fw));
        }
        zl += zw;
        if (nfit == 0) continue;
        std::vector<Pick> next = dp;
        for (std::size_t n0 = 0; n0 < dp.size(); ++n0) {
          if (!dp[n0].ok) continue;
          for (int kl = 0; kl <= nfit; ++kl) {
            double sl = 0;
            for (int j = 0; j < kl; ++j) sl += slots[static_cast<std::size_t>(j)].value;
            for (int kr = (kl == 0 ? 1 : 0); kl + kr <= nfit; ++kr) {
              double sum = dp[n0].sum + sl;
              for (int j = 0; j < kr; ++j) sum += slots[static_cast<std::size_t>(nfit + j)].value;
              const std::size_t n = n0 + static_cast<std::size_t>(kl + kr);
              if (next.size() <= n) next.resize(n + 1);
              if (next[n].ok && !better(sum, next[n].sum)) continue;
              Pick pk;
              pk.ok = true;
              pk.sum = sum;
              pk.el = dp[n0].el;
              for (int j = 0; j < kl; ++j) pk.el.push_back(&slots[static_cast<std::size_t>(j)]);
              for (int j = 0; j < kr; ++j) pk.el.push_back(&slots[static_cast<std::size_t>(nfit + j)]);
              next[n] = std::move(pk);
            }
          }
        }
        dp = std::move(next);
      }
      auto coef = [&](int n) {
        switch (c.trans) {
          case TransLaw::Fasc61Bc: return st1_bc(cls, n);
          case TransLaw::Fasc61Bt: return st1_bt(cls, n);
          case TransLaw::Table:
            return c.coef_trans.empty()
                       ? 1.0
                       : c.coef_trans[std::min(static_cast<std::size_t>(n), c.coef_trans.size()) - 1];
          case TransLaw::None: return 1.0;
        }
        return 1.0;
      };
      int n_max = static_cast<int>(dp.size()) - 1;
      if (c.max_file > 0) n_max = std::min(n_max, c.max_file);
      int bn = 0;
      double bestv = 0;
      for (int n = 1; n <= n_max; ++n) {
        if (!dp[static_cast<std::size_t>(n)].ok) continue;
        const double cn = coef(n);
        if (cn == 0) continue;
        const double v = cn * dp[static_cast<std::size_t>(n)].sum;
        if (bn == 0 || better(v, bestv)) {
          bn = n;
          bestv = v;
        }
      }
      out = Evaluated{};
      out.pl.dir = core.gdir;
      if (bn == 0 || !better(bestv, 0.0)) return true;
      out.value = bestv;
      out.pl.n_lanes = bn;
      const double cb = coef(bn);
      for (const Elem* e : dp[static_cast<std::size_t>(bn)].el)
        for (PlacedGroup g : e->ev.pl.groups) {
          for (auto& [a, P] : g.pattern.axles) P *= cb;
          for (auto& seg : g.pattern.segments) seg[2] *= cb;
          out.pl.groups.push_back(g);
        }
      return true;
    }
    case SurchLoad::Al:
    case SurchLoad::AlPlancher:
    case SurchLoad::Trot:
    case SurchLoad::AlPieton:
    case SurchLoad::Alg: {
      std::function<double(double)> q;
      switch (d.load) {
        case SurchLoad::Al: q = [tf](double L) { return st1_al(L) * tf; }; break;
        case SurchLoad::AlPlancher: q = [tf](double L) { return st1_al_plancher(L) * tf; }; break;
        case SurchLoad::Trot: q = [tf](double) { return st1_trot() * tf; }; break;
        case SurchLoad::AlPieton: q = [tf](double L) { return st1_al_pieton(L) * tf; }; break;
        default: {
          const AlgDef* a = &d.alg;
          q = [a, tf](double L) { return a->density(L) * tf; };
        }
      }
      const bool per_lane = d.load == SurchLoad::Al || d.load == SurchLoad::AlPlancher ||
                            (d.load == SurchLoad::Alg && d.alg.larg_voie > 0);
      std::vector<Elem> el;
      if (per_lane) {
        std::vector<DeckLane> lanes = t.lanes(d.zones, false);
        if (lanes.empty()) return false;
        const double v0 = d.load == SurchLoad::Alg ? d.alg.larg_voie : st1_v0(cls);
        for (const DeckLane& l : lanes) {
          double k = t.rep_coef(study_bar, l.y, l.width / 2);
          double gain = v0 * k * d.pond * core.dir_scale;
          Elem e;
          e.ev = eval_distributed(lane_il(-l.y), d, gain, q, maximise, core.gdir, false, 1, l.width);
          e.value = e.ev.value;
          e.ev.pl.udl_ecc.assign(e.ev.pl.udl.size(), -l.y * ecc_on);
          el.push_back(e);
        }
        const bool a1_law = d.load == SurchLoad::Al;
        combine(el, static_cast<int>(el.size()),
                [&](int n) {
                  if (a1_law) return st1_a1(cls, n);
                  if (d.load == SurchLoad::Alg && !d.alg.coef_trans.empty())
                    return d.alg.coef_trans[std::min(static_cast<std::size_t>(n), d.alg.coef_trans.size()) - 1];
                  return 1.0;
                },
                true);
        return true;
      }
      // Trottoirs, A(l) piétons, ALG sans LARG VOIE : chaque zone chargée
      // seule si elle est défavorable, sans coefficient.
      std::vector<int> ids = t.resolved_zones(d.zones);
      if (ids.empty()) return false;
      out = Evaluated{};
      out.pl.dir = core.gdir;
      for (int id : ids) {
        double zl = -t.width() / 2;
        for (int j = 1; j < id; ++j) zl += t.zones[static_cast<std::size_t>(j) - 1].width;
        const double w = t.zones[static_cast<std::size_t>(id) - 1].width;
        const double y = zl + w / 2;
        double gain = w * t.rep_coef(study_bar, y, w / 2) * d.pond * core.dir_scale;
        Evaluated e = eval_distributed(lane_il(-y), d, gain, q, maximise, core.gdir, false, 1, w);
        if (!better(e.value, 0.0)) continue;
        out.value += e.value;
        for (const auto& u : e.pl.udl) {
          out.pl.udl.push_back(u);
          out.pl.udl_ecc.push_back(-y * ecc_on);
        }
        ++out.pl.n_lanes;
      }
      return true;
    }
    case SurchLoad::Lm1: {
      // EN 1991-2 4.2.4 / 4.3.2 in SPATIALE / GRILL: the notional lanes of
      // each loaded zone form one block of 3 m strips that sits against
      // either edge of the zone (the residual area takes the rest); each
      // strip has its own line eta_F + y eta_T and the lane numbers (TS and
      // UDL of rank 1, 2, ...) go to the most unfavourable strips first. A
      // strip that would relieve the effect stays empty. Without torsion
      // (ilT = 0) every strip reads the same line and the scalar study is
      // recovered.
      if (!core.ecc) return false;
      const double psi_ts = d.freq ? d.ec.psi_freq_ts : 1.0;
      const double psi_udl = d.freq ? d.ec.psi_freq_udl : 1.0;
      const EurocodeParams& e = d.ec;
      std::vector<int> zids = t.resolved_zones(d.zones);
      std::vector<DeckLane> all = t.lanes(d.zones, true);
      if (zids.empty()) return false;
      std::vector<double> zleft(t.zones.size(), 0.0);
      {
        double acc = -t.width() / 2;
        for (std::size_t z = 0; z < t.zones.size(); ++z) {
          zleft[z] = acc;
          acc += t.zones[z].width;
        }
      }
      ConvoyPattern unit_ts;
      unit_ts.axles = {{0.0, 1.0}, {1.2, 1.0}};
      // per strip (centre y, width w): its unit TS placement and unit UDL zones
      struct Strip {
        double y = 0, w = 0, k = 1;
        InfluenceLine il;
        ConvoyPlacement ts;
        double udl = 0;  // effect of a unit density on the favourable zones
        std::vector<SignZone> zones;
      };
      auto make_strip = [&](double y, double w) {
        Strip st;
        st.y = y;
        st.w = w;
        st.k = t.rep_coef(study_bar, y, w / 2);
        st.il = lane_il(-y);
        const bool up = maximise ? st.k >= 0 : st.k < 0;
        st.ts = convoy_placement(st.il, unit_ts, 1, 1e-6, d.pas_conv, up, d.pas_conv_grid_only);
        for (const SignZone& z : il_sign_zones(st.il, d.seuil))
          if (up ? z.integral > 0 : z.integral < 0) {
            st.udl += z.integral;
            st.zones.push_back(z);
          }
        return st;
      };
      // arrangements: per loaded zone, kl of its notional lanes packed from
      // the left edge and the other kr = n - kl from the right edge, the
      // residual area between them. EN 1991-2 4.2.4(2): the location and the
      // numbering of the lanes are free, so two lanes may run along opposite
      // kerbs with the residual strip between them (as in ST1: the lane-2 TS
      // and UDL then act at the kerb eccentricity) -- the same packing as
      // the Bc files of the Convoy case above. Without residual width the
      // splits coincide and only one is kept.
      struct ZoneLanes {
        double zl = 0, zw = 0, rw = 0;
        std::vector<double> ws;
        int nsplit = 1;
      };
      std::vector<ZoneLanes> zlanes;
      std::size_t ncomb = 1;
      for (int id : zids) {
        ZoneLanes z;
        z.zw = t.zones[static_cast<std::size_t>(id) - 1].width;
        z.zl = zleft[static_cast<std::size_t>(id) - 1];
        double used = 0;
        for (const DeckLane& l : all)
          if (l.zone == id) {
            z.ws.push_back(l.width);
            used += l.width;
          }
        z.rw = std::max(0.0, z.zw - used);
        z.nsplit = z.rw > 1e-9 ? static_cast<int>(z.ws.size()) + 1 : 1;
        ncomb *= static_cast<std::size_t>(z.nsplit);
        zlanes.push_back(z);
      }
      Evaluated best;
      bool have = false;
      for (std::size_t comb = 0; comb < ncomb; ++comb) {
        std::vector<Strip> lanes_s, resid;
        std::size_t rem = comb;
        for (const ZoneLanes& z : zlanes) {
          const int n = static_cast<int>(z.ws.size());
          const int kl = z.nsplit == 1 ? n : static_cast<int>(rem % static_cast<std::size_t>(z.nsplit));
          rem /= static_cast<std::size_t>(z.nsplit);
          double y0 = z.zl;
          for (int j = 0; j < kl; ++j) {
            const double w = z.ws[static_cast<std::size_t>(j)];
            lanes_s.push_back(make_strip(y0 + w / 2, w));
            y0 += w;
          }
          double y1 = z.zl + z.zw;
          for (int j = n - 1; j >= kl; --j) {
            const double w = z.ws[static_cast<std::size_t>(j)];
            y1 -= w;
            lanes_s.push_back(make_strip(y1 + w / 2, w));
          }
          if (z.rw > 1e-9) resid.push_back(make_strip((y0 + y1) / 2, z.rw));
        }
        {
          // lane numbers to the strips, most unfavourable first
          Evaluated ev;
          ev.pl.dir = core.gdir;
          std::vector<bool> taken(lanes_s.size(), false);
          int rank = 0;
          for (;;) {
            const int r = rank + 1;
            int bj = -1;
            double bv = 0;
            for (std::size_t j = 0; j < lanes_s.size(); ++j) {
              if (taken[j]) continue;
              const Strip& st = lanes_s[j];
              const double v = st.k * d.pond *
                               (psi_ts * e.aQ(r) * e.Q(r) * st.ts.value +
                                psi_udl * e.aq(r) * e.q(r) * st.w * st.udl);
              if (bj < 0 || better(v, bv)) {
                bj = static_cast<int>(j);
                bv = v;
              }
            }
            if (bj < 0 || !better(bv, 0.0)) break;
            taken[static_cast<std::size_t>(bj)] = true;
            rank = r;
            const Strip& st = lanes_s[static_cast<std::size_t>(bj)];
            ev.value += bv;
            const double P = st.k * d.pond * psi_ts * e.aQ(r) * e.Q(r);
            if (P != 0 && !st.ts.positions.empty()) {
              PlacedGroup g;
              g.pattern.axles = {{0.0, P}, {1.2, P}};
              g.positions = st.ts.positions;
              g.ecc = -st.y * ecc_on;
              ev.pl.groups.push_back(g);
            }
            const double line = st.k * d.pond * psi_udl * e.aq(r) * e.q(r) * st.w;
            for (const SignZone& z : st.zones) {
              ev.pl.udl.push_back({z.s0, z.s1, line});
              ev.pl.udl_ecc.push_back(-st.y * ecc_on);
            }
          }
          ev.pl.n_lanes = rank;
          for (const Strip& st : resid) {
            const double line = st.k * d.pond * psi_udl * e.aqr * e.qrk * st.w;
            const double v = line * st.udl;
            if (!better(v, 0.0)) continue;
            ev.value += v;
            for (const SignZone& z : st.zones) {
              ev.pl.udl.push_back({z.s0, z.s1, line});
              ev.pl.udl_ecc.push_back(-st.y * ecc_on);
            }
          }
          if (!have || better(ev.value, best.value)) {
            best = std::move(ev);
            have = true;
          }
        }
      }
      out = have ? best : Evaluated{};
      out.pl.dir = core.gdir;
      return true;
    }
    default: return false;
  }
}

void evaluate_target(const Tablier& t, const SurchDef& d, const SurchCore& core,
                     const InfluenceLine& il, Id study_bar, SurchTargetResult& out,
                     const InfluenceLine* ilT = nullptr) {
  InfluenceLine zero_t;
  if (!ilT && !t.rep_trans.empty() && !t.zones.empty() && !core.local) {
    zero_t = il;
    std::fill(zero_t.eta.begin(), zero_t.eta.end(), 0.0);
    std::fill(zero_t.dl.begin(), zero_t.dl.end(), 0.0);
    std::fill(zero_t.dr.begin(), zero_t.dr.end(), 0.0);
    ilT = &zero_t;  // PLANE + REP_TRANS: per-lane coefficients, no torque
  }
  auto eval_with = [&](bool maximise, bool tmax) -> Evaluated {
    if (ilT) {
      Evaluated e;
      if (eval_eccentric(il, *ilT, t, d, core, study_bar, maximise, e)) return e;
    }
    switch (d.load) {
      case SurchLoad::Convoy: {
        ConvoyDef c = surch_convoy(d);
        return eval_convoy(il, t, d, c, core.gdir, study_bar, maximise, tmax);
      }
      case SurchLoad::Al:
      case SurchLoad::AlFrein:
      case SurchLoad::AlPlancher:
      case SurchLoad::AlPlancherFrein:
      case SurchLoad::Trot:
      case SurchLoad::AlPieton:
      case SurchLoad::Alg: {
        std::vector<DeckLane> lanes = t.lanes(d.zones, false);
        double phys = t.loaded_width(d.zones);
        int cls = t.fasc61_class();
        double gain = 0;
        int n_used = static_cast<int>(lanes.size());
        // Les charges dont la densité est pondérée VOIE PAR VOIE (a1/a2, ou
        // ALG avec LARG VOIE) : la largeur pondérée vaut a1 * v0 * somme(k),
        // puisque a2 * largeur réelle de voie = v0 (p.138).
        bool per_lane = (d.load == SurchLoad::Al || d.load == SurchLoad::AlFrein ||
                         d.load == SurchLoad::AlPlancher ||
                         d.load == SurchLoad::AlPlancherFrein ||
                         (d.load == SurchLoad::Alg && d.alg.larg_voie > 0));
        if (lanes.empty()) {
          // p.127 : sans ZONE_TRANS, les charges réparties sont appliquées sur
          // 1,00 m de largeur, sans pondération transversale.
          gain = 1.0;
          n_used = 1;
          phys = 1.0;
        } else if (per_lane) {
          double v0 = (d.load == SurchLoad::Alg) ? d.alg.larg_voie : st1_v0(cls);
          std::vector<double> k = lane_coefs(t, lanes, study_bar, 0, tmax);
          bool a1_law = (d.load == SurchLoad::Al || d.load == SurchLoad::AlFrein);
          double sum = 0, best = 0, bw = 0;
          bool first = true;
          for (std::size_t i = 0; i < k.size(); ++i) {
            int n = static_cast<int>(i) + 1;
            sum += k[i];
            double a1 = 1.0;
            if (a1_law) a1 = st1_a1(cls, n);
            if (d.load == SurchLoad::Alg && !d.alg.coef_trans.empty())
              a1 = d.alg.coef_trans[std::min(static_cast<std::size_t>(n),
                                             d.alg.coef_trans.size()) -
                                    1];
            if (a1 == 0) continue;  // combinaison interdite (0 dans la table)
            double f = a1 * v0 * sum;
            double w = 0;
            for (int j = 0; j < n; ++j) w += lanes[static_cast<std::size_t>(j)].width;
            if (first || (tmax ? f > best : f < best)) {
              best = f;
              n_used = n;
              bw = w;
              first = false;
            }
          }
          require(!first, "SURCH: aucune combinaison de voies autorisee (table a1)");
          gain = best;
          phys = bw;
        } else {
          // TROT / AL PIETON / ALG sans LARG VOIE : toute la largeur chargée.
          for (int id : t.resolved_zones(d.zones)) {
            const TransZone& z = t.zones[static_cast<std::size_t>(id) - 1];
            double acc = -t.width() / 2;
            for (int j = 1; j < id; ++j) acc += t.zones[static_cast<std::size_t>(j) - 1].width;
            gain += z.width * t.rep_coef(study_bar, acc + z.width / 2, z.width / 2);
          }
        }
        gain *= d.pond * core.dir_scale;
        double tf = d.tonne_force;
        std::function<double(double)> q;
        switch (d.load) {
          case SurchLoad::Al: q = [tf](double L) { return st1_al(L) * tf; }; break;
          case SurchLoad::AlPlancher:
            q = [tf](double L) { return st1_al_plancher(L) * tf; };
            break;
          case SurchLoad::Trot: q = [tf](double) { return st1_trot() * tf; }; break;
          case SurchLoad::AlPieton: q = [tf](double L) { return st1_al_pieton(L) * tf; }; break;
          case SurchLoad::AlFrein: {
            double w = phys;
            q = [tf, w](double L) { return st1_al_frein(L, L * w) * tf; };
            break;
          }
          case SurchLoad::AlPlancherFrein: {
            double w = phys;
            q = [tf, w](double L) { return st1_al_plancher_frein(L, L * w) * tf; };
            break;
          }
          case SurchLoad::Alg: {
            const AlgDef* a = &d.alg;
            q = [a, tf](double L) { return a->density(L) * tf; };
            break;
          }
          default: fail("SURCH: loi de densite manquante");
        }
        return eval_distributed(il, d, gain, q, maximise, core.gdir, core.local, n_used, phys);
      }
      case SurchLoad::Lm4Cte:
      case SurchLoad::Lm4Var:
      case SurchLoad::TrotEc:
      case SurchLoad::TrotLm1: {
        double phys = t.loaded_width(d.zones);
        double gain = 0;
        if (t.zones.empty()) {
          gain = 1.0;
          phys = 1.0;
        } else {
          for (int id : t.resolved_zones(d.zones)) {
            const TransZone& z = t.zones[static_cast<std::size_t>(id) - 1];
            double acc = -t.width() / 2;
            for (int j = 1; j < id; ++j) acc += t.zones[static_cast<std::size_t>(j) - 1].width;
            gain += z.width * t.rep_coef(study_bar, acc + z.width / 2, z.width / 2);
          }
        }
        gain *= d.pond;
        std::function<double(double)> q;
        const double fs = d.ec.force_scale;  // N -> the deck's force unit
        if (d.load == SurchLoad::Lm4Cte)
          q = [fs](double) { return 5.0e3 * fs; };  // p.166 : qfk = 5 kN/m2
        else if (d.load == SurchLoad::Lm4Var)
          // p.166 : qfk = 2 + 120/(L+30) kN/m2, borné à [2,5 ; 5,0]
          q = [fs](double L) { return std::clamp(2.0 + 120.0 / (L + 30.0), 2.5, 5.0) * 1e3 * fs; };
        else if (d.load == SurchLoad::TrotEc)
          q = [fs](double) { return 5.0e3 * fs; };  // p.166 : trottoir Eurocode
        else
          q = [fs](double) { return 3.0e3 * fs; };  // p.166 : trottoir accompagnant le LM1
        return eval_distributed(il, d, gain, q, maximise, core.gdir, false,
                                static_cast<int>(t.lanes(d.zones, true).size()), phys);
      }
      case SurchLoad::Lm1: {
        double psi_ts = d.freq ? d.ec.psi_freq_ts : 1.0;
        double psi_udl = d.freq ? d.ec.psi_freq_udl : 1.0;
        Lm1Loads L = lm1_loads(t, d, study_bar, tmax, psi_ts, psi_udl);
        Evaluated ev;
        ev.pl.dir = core.gdir;
        ev.pl.n_lanes = L.n_lanes;
        ConvoyPattern ts;
        ts.axles = {{0.0, L.ts_axle * d.pond}, {1.2, L.ts_axle * d.pond}};
        ConvoyPlacement cp =
            convoy_placement(il, ts, 1, 1e-6, d.pas_conv, maximise, d.pas_conv_grid_only);
        // UDL : densité constante, donc pas de point fixe sur la longueur.
        // Each family of coefficient sign loads its own unfavourable zones.
        double line = L.udl * d.pond, line_o = L.udl_other * d.pond;
        std::vector<SignZone> zones = il_sign_zones(il, d.seuil);
        double udl_val = 0, len = 0;
        for (const SignZone& z : zones)
          for (double ln : {line, line_o}) {
            double c = ln * z.integral;
            if (ln != 0 && (maximise ? c > 0 : c < 0)) {
              udl_val += c;
              len += z.s1 - z.s0;
              ev.pl.udl.push_back({z.s0, z.s1, ln});
            }
          }
        ev.value = cp.value + udl_val;
        if (!cp.positions.empty()) {
          PlacedGroup g;
          g.pattern = ts;
          g.positions = cp.positions;
          ev.pl.groups.push_back(g);
        }
        ev.pl.loaded_length = len;
        ev.pl.density = line;
        return ev;
      }
      case SurchLoad::Lm2: {
        double psi = d.freq ? d.ec.psi_freq_lm2 : 1.0;
        // p.171 : le déplacement du LM2 est entièrement libre longitudinalement
        // ET transversalement -> coefficient transversal le plus défavorable.
        double kbest = 1.0;
        if (!t.zones.empty()) {
          double half = t.width() / 2;
          kbest = t.rep_extreme(study_bar, -half, half, 0.0, tmax);
        }
        ConvoyPattern ax;
        ax.axles = {{0.0, psi * d.ec.beta_lm2 * 400e3 * d.ec.force_scale * kbest * d.pond}};
        ConvoyPlacement cp =
            convoy_placement(il, ax, 1, 1e-6, d.pas_conv, maximise, d.pas_conv_grid_only);
        Evaluated ev;
        ev.value = cp.value;
        ev.pl.dir = core.gdir;
        ev.pl.trans_factor = kbest;
        ev.pl.n_lanes = 1;
        if (!cp.positions.empty()) {
          PlacedGroup g;
          g.pattern = ax;
          g.positions = cp.positions;
          ev.pl.groups.push_back(g);
        }
        return ev;
      }
      case SurchLoad::Lm3: {
        // Un SEUL véhicule VS (p.167) + LM1 d'accompagnement. Transversalement
        // (p.171) le VS roule dans la bande ZONEVS/YVS (ou sur la voie 1 du
        // LM1, EXC / DEBORDVOIE) au point le plus défavorable des
        // coefficients REP_TRANS ; le LM1 n'est exclu, sur l'emprise du VS
        // prolongée de LGEXCLU, que dans les voies que le VS occupe (CORRELE,
        // 1 ou 2 voies selon sa largeur — EN 1991-2 A.3) : les autres voies
        // restent chargées.
        ConvoyDef c = surch_convoy(d);
        double psi_ts = 0, psi_udl = 0;
        if (d.lm3_lm1 == 1) {
          psi_ts = d.ec.psi_freq_ts;
          psi_udl = d.ec.psi_freq_udl;
        } else if (d.lm3_lm1 == 2) {
          psi_ts = 1.0;
          psi_udl = 1.0;
        } else if (d.lm3_lm1 == 3) {
          psi_ts = d.ec.psi_comb_ts;
          psi_udl = d.ec.psi_comb_udl;
        }
        Lm1Loads L = lm1_loads(t, d, study_bar, tmax, psi_ts, psi_udl);
        const int nvs = c.largeur > 3.0 + 1e-9 ? 2 : 1;
        std::vector<int> vs_ranks = {d.correle_j};
        if (nvs == 2) vs_ranks.push_back(d.correle_k);
        // Coefficient transversal du VS.
        double k_vs = 1.0;
        if (!t.rep_trans.empty() && !t.zones.empty()) {
          const double w = c.largeur > 0 ? c.largeur : 3.0;
          std::vector<double> cands;
          if (d.zonevs > 0 && d.zonevs <= static_cast<int>(t.zones.size())) {
            double zl = -t.width() / 2;
            for (int z = 1; z < d.zonevs; ++z) zl += t.zones[static_cast<std::size_t>(z) - 1].width;
            const double zw = t.zones[static_cast<std::size_t>(d.zonevs) - 1].width;
            double y1 = d.yvs_abs ? d.yvs1 : d.yvs1 * zw;
            double y2 = d.yvs2 < 0 ? zw : (d.yvs_abs ? d.yvs2 : d.yvs2 * zw);
            double lo = zl + y1 + w / 2, hi = zl + y2 - w / 2;
            if (hi < lo) lo = hi = zl + 0.5 * (y1 + y2);
            // the exact extreme over the band, through its abscissa
            k_vs = t.rep_extreme(study_bar, lo, hi, w / 2, tmax);  // cands left empty
          } else if (!L.lanes.empty()) {
            double left = 1e300, right = -1e300;
            for (int r : vs_ranks)
              if (r >= 1 && r <= static_cast<int>(L.lanes.size())) {
                const LaneRank& ln = L.lanes[static_cast<std::size_t>(r) - 1];
                left = std::min(left, ln.y - ln.width / 2);
                right = std::max(right, ln.y + ln.width / 2);
              }
            if (left <= right) {
              double cc = 0.5 * (left + right);
              if (d.has_exc) cands = {cc + d.vs_exc, cc - d.vs_exc};
              else if (d.has_debv) cands = {left + d.vs_debv + w / 2, right - d.vs_debv - w / 2};
              else cands = {cc};
            }
          }
          bool kfirst = true;
          for (double y : cands) {
            double k = t.rep_coef(study_bar, y, w / 2);
            if (kfirst || (tmax ? k > k_vs : k < k_vs)) k_vs = k;
            kfirst = false;
          }
        }
        ConvoyPattern vs = convoy_pattern(c, d.tonne_force, t.ldif);
        for (auto& ax : vs.axles) ax.second *= d.pondvs * d.pond * k_vs;
        for (auto& seg : vs.segments) seg[2] *= d.pondvs * d.pond * k_vs;
        // LM1 hors fenêtre (toutes les voies) / dans la fenêtre (sans les
        // voies du VS).
        double ts_in = L.ts_axle, udl_in = L.udl;
        for (int r : vs_ranks)
          if (r >= 1 && r <= static_cast<int>(L.ts_rank.size())) {
            ts_in -= L.ts_rank[static_cast<std::size_t>(r) - 1];
            udl_in -= L.udl_rank[static_cast<std::size_t>(r) - 1];
          }
        ConvoyPattern ts_all, ts_oth, ts_vsl;
        ts_all.axles = {{0.0, L.ts_axle * d.pond}, {1.2, L.ts_axle * d.pond}};
        ts_oth.axles = {{0.0, ts_in * d.pond}, {1.2, ts_in * d.pond}};
        const double ts_vl = L.ts_axle - ts_in;  // the tandems of the VS's lanes
        ts_vsl.axles = {{0.0, ts_vl * d.pond}, {1.2, ts_vl * d.pond}};
        const double line_out = L.udl * d.pond, line_in = udl_in * d.pond;
        std::vector<SignZone> zones = il_sign_zones(il, d.seuil);
        std::vector<double> vs_pos =
            candidate_positions(il, vs, d.pas_conv, d.pas_conv_grid_only);
        std::vector<double> ts_grid =
            (L.ts_axle != 0) ? candidate_positions(il, ts_all, d.pas_conv, d.pas_conv_grid_only)
                             : std::vector<double>{};
        Evaluated best;
        best.pl.dir = core.gdir;
        best.pl.n_lanes = L.n_lanes;
        bool first = true;
        const IlPrefix pre(il);
        auto good = [&](double e) { return maximise ? e > 0 : e < 0; };
        for (double r : vs_pos) {
          double e_vs = pattern_effect(il, pre, vs, r);
          double w0 = r - d.lgexclu, w1 = r + c.longueur + d.lgexclu;
          double e_udl = 0;
          std::vector<std::array<double, 3>> segs;
          for (const SignZone& z : zones) {
            double e_z = 0;
            std::vector<std::array<double, 3>> zs;
            auto add = [&](double a, double b, double line) {
              if (b <= a || line == 0) return;
              e_z += line * pre.integral(a, b);
              zs.push_back({a, b, line});
            };
            add(z.s0, std::min(z.s1, w0), line_out);
            add(std::max(z.s0, w0), std::min(z.s1, w1), line_in);
            add(std::max(z.s0, w1), z.s1, line_out);
            if (!good(e_z)) continue;
            e_udl += e_z;
            segs.insert(segs.end(), zs.begin(), zs.end());
          }
          // Tandems (p.174 : « la position longitudinale de chacun est
          // indépendante ») : ceux des voies du VS hors de la fenêtre, ceux
          // des autres voies n'importe où, chacun à sa place la plus
          // défavorable (R9 1498 / ST1 1494, 1467 avec un seul
          // tandem commun).
          double e_ts = 0;
          std::vector<PlacedGroup> ts_groups;
          auto place_ts = [&](const ConvoyPattern& pat, double axle, bool outside) {
            if (ts_grid.empty() || axle == 0) return;
            bool tfirst = true;
            double bestts = 0, bestr = 0;
            for (double rt : ts_grid) {
              if (outside && rt + 1.2 >= w0 && rt <= w1) continue;
              double e = pattern_effect(il, pre, pat, rt);
              if (tfirst || (maximise ? e > bestts : e < bestts)) {
                bestts = e;
                bestr = rt;
                tfirst = false;
              }
            }
            if (tfirst || !good(bestts)) return;
            e_ts += bestts;
            PlacedGroup gt;
            gt.pattern = pat;
            gt.positions = {bestr};
            ts_groups.push_back(gt);
          };
          place_ts(ts_vsl, ts_vl, true);
          place_ts(ts_oth, ts_in, false);
          double tot = e_vs + e_udl + e_ts;
          if (first || (maximise ? tot > best.value : tot < best.value)) {
            best = Evaluated{};
            best.value = tot;
            best.pl.dir = core.gdir;
            best.pl.n_lanes = L.n_lanes;
            PlacedGroup gv;
            gv.pattern = vs;
            gv.positions = {r};
            best.pl.groups.push_back(gv);
            for (const PlacedGroup& g : ts_groups) best.pl.groups.push_back(g);
            best.pl.udl = segs;
            first = false;
          }
        }
        // Without the VS (a VS and its exclusion longer than the
        // bridge): the accompanying LM1 alone on every lane is a member of
        // the family.
        {
          double e_udl = 0;
          std::vector<std::array<double, 3>> segs;
          for (const SignZone& z : zones) {
            if (line_out == 0 || z.s1 <= z.s0) continue;
            const double e_z = line_out * pre.integral(z.s0, z.s1);
            if (!good(e_z)) continue;
            e_udl += e_z;
            segs.push_back({z.s0, z.s1, line_out});
          }
          double e_ts = 0;
          PlacedGroup gt;
          bool have_ts = false;
          if (!ts_grid.empty()) {
            bool tfirst = true;
            double bestts = 0, bestr = 0;
            for (double rt : ts_grid) {
              double e = pattern_effect(il, pre, ts_all, rt);
              if (tfirst || (maximise ? e > bestts : e < bestts)) {
                bestts = e;
                bestr = rt;
                tfirst = false;
              }
            }
            if (!tfirst && good(bestts)) {
              e_ts = bestts;
              gt.pattern = ts_all;
              gt.positions = {bestr};
              have_ts = true;
            }
          }
          const double tot = e_udl + e_ts;
          if (first || (maximise ? tot > best.value : tot < best.value)) {
            best = Evaluated{};
            best.value = tot;
            best.pl.dir = core.gdir;
            best.pl.n_lanes = L.n_lanes;
            if (have_ts) best.pl.groups.push_back(gt);
            best.pl.udl = segs;
          }
        }
        return best;
      }
    }
    fail("SURCH: type de charge non gere");
  };

  auto run = [&](bool maximise) -> Evaluated {
    Evaluated a = eval_with(maximise, true);
    if (t.rep_trans.empty()) return a;  // sans REP_TRANS les deux coins coïncident
    Evaluated b = eval_with(maximise, false);
    return (maximise ? (a.value >= b.value) : (a.value <= b.value)) ? a : b;
  };
  Evaluated hi = run(true), lo = run(false);
  out.hi.value = hi.value;
  out.hi.placement = hi.pl;
  out.lo.value = lo.value;
  out.lo.placement = lo.pl;
}

// Corps commun : évaluation + concomitance sur une liste de réponses.
void solve_targets(const Model& model, const Tablier& t, const SurchDef& d,
                   const SurchOptions& opt, const SurchCore& core,
                   const std::vector<SurchTarget>& targets, const std::vector<bool>& is_reac,
                   SurchResults& out) {
  std::vector<ResponseKind> responses;
  responses.reserve(targets.size());
  for (const SurchTarget& tg : targets) responses.push_back(tg.response);
  std::vector<InfluenceLine> ilT;
  const auto tm0 = std::chrono::steady_clock::now();
  std::vector<InfluenceLine> ils =
      il_batch(model, core.path, core.dirs, responses, core.ecc ? &core.tdirs : nullptr, &ilT);
  const auto tm1 = std::chrono::steady_clock::now();
  struct TimingReport {
    std::chrono::steady_clock::time_point a, b;
    std::size_t n;
    ~TimingReport() {
      if (std::getenv("KETA_TIMING"))
        std::fprintf(stderr, "SURCH run: %zu cibles, lignes %.2f s, placements %.2f s\n", n,
                     std::chrono::duration<double>(b - a).count(),
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - b).count());
    }
  } timing_report{tm0, tm1, targets.size()};

  std::unique_ptr<StaticContext> ctx;
  if (opt.concomitant) ctx = std::make_unique<StaticContext>(StaticContext::build(model));
  StaticOptions so;
  so.stations_per_bar = opt.out_stations;

  std::vector<SurchTargetResult> evals(targets.size());
  parallel_for_index(targets.size(), [&](std::size_t i) {
    SurchTargetResult& tr = evals[i];
    tr.target = targets[i];
    tr.il = ils[i];
    evaluate_target(t, d, core, tr.il, targets[i].bar, tr, core.ecc ? &ilT[i] : nullptr);
  });
  for (std::size_t i = 0; i < targets.size(); ++i) {
    SurchTargetResult tr = std::move(evals[i]);
    if (ctx) {
      for (SurchExtreme* ex : {&tr.lo, &tr.hi}) {
        if (ex->placement.empty()) continue;
        LoadCase lc = placement_case(model, core.spans, ex->placement, "SURCH");
        if (lc.loads.empty()) continue;
        ex->state = ctx->solve_case(lc, so);
        ex->has_state = true;
      }
    }
    if (is_reac[i])
      out.reactions.push_back(std::move(tr));
    else
      out.sections.push_back(std::move(tr));
  }
}

}  // namespace

SurchTargetResult surch_response(const Model& model, const Tablier& tablier, const SurchDef& def,
                                 const ResponseKind& response, const SurchOptions& opt) {
  tablier.validate();
  def.validate();
  model.index();
  // The studied section joins the path refinement (its jump stations).
  std::vector<SurchSection> secs = opt.sections;
  if (response.what == ResponseKind::What::BarForce ||
      response.what == ResponseKind::What::BarDisplacement) {
    SurchSection sec;
    sec.bar = response.bar;
    sec.x = response.x;
    secs.push_back(sec);
  }
  SurchCore core = build_core(model, tablier, def, opt, secs);
  SurchTarget tg;
  tg.response = response;
  tg.bar = response.what == ResponseKind::What::BarForce ? response.bar : 0;
  tg.x = response.x;
  tg.node = response.node;
  tg.comp = response.comp;
  SurchResults out;
  out.path = core.path;
  bool reac = response.what == ResponseKind::What::BarForce ? false : true;
  solve_targets(model, tablier, def, opt, core, {tg}, {reac}, out);
  return reac ? out.reactions.front() : out.sections.front();
}

SurchResults run_surch(const Model& model, const Tablier& tablier, const SurchDef& def,
                       const SurchOptions& opt) {
  tablier.validate();
  def.validate();
  model.index();
  SurchCore core = build_core(model, tablier, def, opt);

  SurchResults out;
  out.path = core.path;

  std::vector<SurchTarget> targets;
  std::vector<bool> is_reac;

  std::vector<SurchSection> secs = opt.sections;
  if (secs.empty()) {
    int n = std::max(2, opt.stations_per_bar);
    for (const DeckSpan& s : core.spans)
      for (int k = 0; k < n; ++k) secs.push_back({s.bar, s.L * k / (n - 1)});
  }
  for (const SurchSection& s : secs) {
    SurchTarget tg;
    switch (def.priv) {
      case SurchPriv::Effort:
        tg.response = ResponseKind::bar_force(s.bar, s.x, def.priv_comp);
        tg.bar = s.bar;
        tg.x = s.x;
        break;
      case SurchPriv::Depla: {
        // Divergence documentée : le manuel (p.137) n'autorise l'enveloppe des
        // déplacements qu'aux sections d'étude de barre ; keta la porte sur le
        // NŒUD le plus proche de la section (le moteur de lignes d'influence
        // n'extrait pas la flèche à une station).
        const Bar& b = model.bar(s.bar);
        Id nd = (s.x <= model.bar_chord_length(b) / 2) ? b.ni : b.nj;
        tg.response = ResponseKind::displacement(nd, def.priv_comp);
        tg.node = nd;
        // le coefficient REP_TRANS reste celui de la barre d'étude
        tg.bar = s.bar;
        break;
      }
      case SurchPriv::Reac:
        fail("SURCH: composante privilegiee REAC — la renseigner via SurchOptions::reactions");
    }
    tg.comp = def.priv_comp;
    targets.push_back(tg);
    is_reac.push_back(false);
  }

  std::vector<std::pair<Id, int>> reacs = opt.reactions;
  if (reacs.empty() && opt.auto_reactions) {
    for (const Support& s : model.supports) {
      std::size_t c = static_cast<std::size_t>(opt.reac_comp);
      require(opt.reac_comp >= 0 && opt.reac_comp <= 5, "SURCH: composante de reaction 0..5");
      if (s.fixed[c] || s.spring[c] != 0) reacs.emplace_back(s.node, opt.reac_comp);
    }
  }
  for (const auto& [nd, c] : reacs) {
    SurchTarget tg;
    tg.response = ResponseKind::reaction(nd, c);
    tg.node = nd;
    tg.comp = c;
    targets.push_back(tg);
    is_reac.push_back(true);
  }
  require(!targets.empty(), "SURCH: aucune section d'etude ni reaction a etudier");

  solve_targets(model, tablier, def, opt, core, targets, is_reac, out);
  return out;
}

SurchEnvelope run_surch_envelope(const Model& model, const Tablier& tablier, const SurchDef& def,
                                 const SurchEnvSpec& spec, const SurchOptions& opt) {
  tablier.validate();
  def.validate();
  model.index();
  for (const std::vector<int>* cs : {&spec.effort_comps, &spec.depla_comps, &spec.reac_comps})
    for (int c : *cs) require(c >= 0 && c <= 5, "SURCH: composante 0..5");
  std::vector<SurchSection> all_sections = spec.effort_sections;
  all_sections.insert(all_sections.end(), spec.depla_sections.begin(), spec.depla_sections.end());
  SurchCore core = build_core(model, tablier, def, opt, all_sections);

  SurchEnvelope env;
  for (const SurchSection& s : spec.effort_sections) env.efforts.push_back({s.bar, s.x, {}, {}});
  for (const SurchSection& s : spec.depla_sections) env.deplas.push_back({s.bar, s.x, {}, {}});
  for (Id n : spec.support_nodes) env.reactions.push_back({n, 0.0, {}, {}});

  // Cibles : (famille, ligne, composante) ; une famille privilégiée n'étudie
  // que sa composante privilégiée.
  struct Target {
    int fam;  // 0 efforts, 1 déplacements, 2 réactions
    std::size_t row;
    int comp;
  };
  std::vector<Target> targets;
  std::vector<ResponseKind> responses;
  auto family = [&](int fam, std::vector<SurchEnvRow>& rows, const std::vector<int>& comps,
                    int priv) {
    std::vector<int> cs = priv >= 0 ? std::vector<int>{priv} : comps;
    for (std::size_t r = 0; r < rows.size(); ++r)
      for (int c : cs) {
        targets.push_back({fam, r, c});
        const SurchEnvRow& row = rows[r];
        responses.push_back(fam == 0   ? ResponseKind::bar_force(row.id, row.x, c)
                            : fam == 1 ? ResponseKind::bar_displacement(row.id, row.x, c)
                                       : ResponseKind::reaction(row.id, c));
      }
  };
  bool fibre_priv = false;
  for (int k : spec.effort_priv_fibre) fibre_priv = fibre_priv || k >= 0;
  // a privileged fibre carries the whole effort family (its concomitants)
  if (!fibre_priv) family(0, env.efforts, spec.effort_comps, spec.effort_priv);
  bool depla_fibre_priv = false;
  for (int k : spec.depla_priv_fibre) depla_fibre_priv = depla_fibre_priv || k >= 0;
  if (!depla_fibre_priv) family(1, env.deplas, spec.depla_comps, spec.depla_priv);
  family(2, env.reactions, spec.reac_comps, spec.reac_priv);
  // Concomitance : état statique rejoué au placement gouvernant, restitué aux
  // sections d'étude exactes.
  const bool any_priv =
      spec.effort_priv >= 0 || spec.depla_priv >= 0 || spec.reac_priv >= 0 || fibre_priv ||
      depla_fibre_priv;
  std::map<Id, std::vector<double>> stations;
  for (const SurchEnvRow& r : env.efforts) stations[r.id].push_back(r.x);
  for (const SurchEnvRow& r : env.deplas) stations[r.id].push_back(r.x);
  for (auto& [b, xs] : stations) {
    std::sort(xs.begin(), xs.end());
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
  }
  std::unique_ptr<StaticContext> ctx;
  if (any_priv) ctx = std::make_unique<StaticContext>(StaticContext::build(model));
  StaticOptions so;
  so.stations_per_bar = 2;
  so.stations_abs = &stations;
  auto station_at = [&](const CaseResult& cr, Id bar, double x) -> const Station& {
    const BarResult& br = cr.bars[static_cast<std::size_t>(model.bar_index(bar))];
    std::size_t best = 0;
    for (std::size_t k = 1; k < br.stations.size(); ++k)
      if (std::fabs(br.stations[k].x - x) < std::fabs(br.stations[best].x - x)) best = k;
    return br.stations[best];
  };
  auto read_state = [&](int fam, const SurchEnvRow& row, const CaseResult& cr) {
    std::array<double, 6> v{};
    if (fam == 0) {
      const Station& st = station_at(cr, row.id, row.x);
      v = {st.N, st.Vy, st.Vz, st.T, st.My, st.Mz};
    } else if (fam == 1) {
      const Station& st = station_at(cr, row.id, row.x);
      const Bar& b = model.bar(row.id);
      Mat3 R = bar_rotation(model.node(b.ni).x, model.node(b.nj).x, b.roll);
      Vec3 d = R.transpose() * st.defl, r = R.transpose() * st.rot;
      v = {d(0), d(1), d(2), r(0), r(1), r(2)};
    } else {
      for (const Reaction& re : cr.reactions)
        if (re.node == row.id)
          for (int c = 0; c < 6; ++c) v[static_cast<std::size_t>(c)] = re.R(c);
    }
    return v;
  };

  // Fibres : les lignes d'influence des composantes qu'elles combinent
  // (famille 0 : efforts de section ; 1 : déplacements globaux, PRESS).
  auto fibre_family = [&](int fam) {
    std::vector<SurchEnvRow>& rows = fam == 0 ? env.efforts : env.deplas;
    const auto& coefs = fam == 0 ? spec.effort_fibres : spec.depla_fibres;
    const auto& privs = fam == 0 ? spec.effort_priv_fibre : spec.depla_priv_fibre;
    struct Fibre {
      std::size_t row, k;
    };
    std::vector<Fibre> fibres;
    std::vector<ResponseKind> fresp;
    std::map<std::pair<std::size_t, int>, std::size_t> fidx;  // (ligne, composante) -> fresp
    for (std::size_t r = 0; r < rows.size() && r < coefs.size(); ++r) {
      const auto& fs = coefs[r];
      rows[r].fibre_lo.assign(fs.size(), 0.0);
      rows[r].fibre_hi.assign(fs.size(), 0.0);
      for (std::size_t k = 0; k < fs.size(); ++k) {
        fibres.push_back({r, k});
        for (int c = 0; c < 6; ++c)
          if (fs[k][static_cast<std::size_t>(c)] != 0.0 && !fidx.count({r, c})) {
            fidx[{r, c}] = fresp.size();
            fresp.push_back(fam == 0 ? ResponseKind::bar_force(rows[r].id, rows[r].x, c)
                                     : ResponseKind::bar_displacement(rows[r].id, rows[r].x, c));
          }
      }
    }
    if (fibres.empty()) return;
    std::vector<InfluenceLine> filT;
    std::vector<InfluenceLine> fils =
        il_batch(model, core.path, core.dirs, fresp, core.ecc ? &core.tdirs : nullptr, &filT);
    std::vector<SurchTargetResult> fev(fibres.size());
    parallel_for_index(fibres.size(), [&](std::size_t i) {
      const Fibre& f = fibres[i];
      const std::array<double, 6>& coef = coefs[f.row][f.k];
      InfluenceLine il, ilt;
      bool first = true;
      for (int c = 0; c < 6; ++c) {
        const double a = coef[static_cast<std::size_t>(c)];
        if (a == 0.0) continue;
        const std::size_t j = fidx.at({f.row, c});
        // a linear combination of lines: the ordinates and the slopes alike
        auto start = [](InfluenceLine& o, const InfluenceLine& src) {
          o.s = src.s;
          o.eta.assign(src.eta.size(), 0.0);
          if (src.cubic()) {
            o.dl.assign(src.eta.size(), 0.0);
            o.dr.assign(src.eta.size(), 0.0);
          }
        };
        auto add = [](InfluenceLine& o, const InfluenceLine& src, double w) {
          require(src.s.size() == o.s.size(), "SURCH: lignes de fibre non alignees");
          for (std::size_t q = 0; q < o.eta.size(); ++q) o.eta[q] += w * src.eta[q];
          if (o.cubic() && src.cubic())
            for (std::size_t q = 0; q < o.eta.size(); ++q) {
              o.dl[q] += w * src.dl[q];
              o.dr[q] += w * src.dr[q];
            }
        };
        if (first) {
          start(il, fils[j]);
          if (core.ecc) start(ilt, filT[j]);
          first = false;
        }
        add(il, fils[j], a);
        if (core.ecc) add(ilt, filT[j], a);
      }
      if (first) return;  // fibre sans composante : 0
      SurchTargetResult& tr = fev[i];
      tr.il = il;
      evaluate_target(tablier, def, core, tr.il, rows[f.row].id, tr, core.ecc ? &ilt : nullptr);
    });
    for (std::size_t i = 0; i < fibres.size(); ++i) {
      SurchEnvRow& row = rows[fibres[i].row];
      row.fibre_lo[fibres[i].k] = fev[i].lo.value;
      row.fibre_hi[fibres[i].k] = fev[i].hi.value;
    }
    // The privileged fibre: the concomitant state at each of its extremes.
    for (std::size_t i = 0; i < fibres.size(); ++i) {
      const Fibre& f = fibres[i];
      if (f.row >= privs.size() || privs[f.row] != static_cast<int>(f.k)) continue;
      SurchEnvRow& row = rows[f.row];
      const auto& fs = coefs[f.row];
      std::vector<double> lo_f(fs.size(), 0.0), hi_f(fs.size(), 0.0);
      for (int side = 0; side < 2; ++side) {
        const SurchExtreme& ex = side == 0 ? fev[i].lo : fev[i].hi;
        std::array<double, 6>& out = side == 0 ? row.lo : row.hi;
        out.fill(0.0);
        if (ex.placement.empty()) continue;
        LoadCase lc = placement_case(model, core.spans, ex.placement, "SURCH");
        if (lc.loads.empty()) continue;
        out = read_state(fam, row, ctx->solve_case(lc, so));
        std::vector<double>& fo = side == 0 ? lo_f : hi_f;
        for (std::size_t k = 0; k < fs.size(); ++k)
          for (int c = 0; c < 6; ++c)
            fo[k] += fs[k][static_cast<std::size_t>(c)] * out[static_cast<std::size_t>(c)];
      }
      // the privileged fibre keeps its own extremes (the same values, without
      // the replay's round-off); the others take their concomitant values
      for (std::size_t k = 0; k < fs.size(); ++k)
        if (k != f.k) {
          row.fibre_lo[k] = lo_f[k];
          row.fibre_hi[k] = hi_f[k];
        }
    }
  };
  {
    const auto tf0 = std::chrono::steady_clock::now();
    fibre_family(0);
    fibre_family(1);
    if (std::getenv("KETA_TIMING"))
      std::fprintf(stderr, "SURCH env: fibres %.2f s\n",
                   std::chrono::duration<double>(std::chrono::steady_clock::now() - tf0).count());
  }
  if (targets.empty()) return env;

  std::vector<InfluenceLine> ilT;
  const auto tm0 = std::chrono::steady_clock::now();
  std::vector<InfluenceLine> ils =
      il_batch(model, core.path, core.dirs, responses, core.ecc ? &core.tdirs : nullptr, &ilT);
  const auto tm1 = std::chrono::steady_clock::now();
  struct TimingReport {
    std::chrono::steady_clock::time_point a, b;
    std::size_t n;
    ~TimingReport() {
      if (std::getenv("KETA_TIMING"))
        std::fprintf(stderr, "SURCH env: %zu cibles, lignes %.2f s, placements %.2f s\n", n,
                     std::chrono::duration<double>(b - a).count(),
                     std::chrono::duration<double>(std::chrono::steady_clock::now() - b).count());
    }
  } timing_report{tm0, tm1, targets.size()};

  std::vector<SurchTargetResult> evals(targets.size());
  parallel_for_index(targets.size(), [&](std::size_t i) {
    const Target& tg = targets[i];
    const std::vector<SurchEnvRow>& rows =
        tg.fam == 0 ? env.efforts : tg.fam == 1 ? env.deplas : env.reactions;
    const SurchEnvRow& row = rows[tg.row];
    SurchTargetResult& tr = evals[i];
    tr.il = ils[i];
    // REP_TRANS of a support: that of a deck bar incident to it (the
    // coefficients are attached to bars, p.130).
    Id study = row.id;
    if (tg.fam == 2) {
      study = 0;
      for (Id b : tablier.bars) {
        const Bar& bb = model.bar(b);
        if (bb.ni == row.id || bb.nj == row.id) {
          study = b;
          break;
        }
      }
    }
    evaluate_target(tablier, def, core, tr.il, study, tr, core.ecc ? &ilT[i] : nullptr);
  });
  for (std::size_t i = 0; i < targets.size(); ++i) {
    const Target& tg = targets[i];
    std::vector<SurchEnvRow>& rows = tg.fam == 0 ? env.efforts : tg.fam == 1 ? env.deplas : env.reactions;
    SurchEnvRow& row = rows[tg.row];
    SurchTargetResult tr = std::move(evals[i]);
    int priv = tg.fam == 0 ? spec.effort_priv : tg.fam == 1 ? spec.depla_priv : spec.reac_priv;
    if (priv < 0) {
      row.lo[static_cast<std::size_t>(tg.comp)] = tr.lo.value;
      row.hi[static_cast<std::size_t>(tg.comp)] = tr.hi.value;
      continue;
    }
    // a privileged component identically null (the Mz of a pinned
    // support): ST1 places nothing, the concomitants are null
    double il_peak = 0;
    for (double e : tr.il.eta) il_peak = std::max(il_peak, std::fabs(e));
    if (core.ecc)  // the lanes' eccentric torque line (the torsion of a deck)
      for (double e : ilT[i].eta) il_peak = std::max(il_peak, std::fabs(e));
    for (int side = 0; side < 2; ++side) {
      const SurchExtreme& ex = side == 0 ? tr.lo : tr.hi;
      std::array<double, 6>& out = side == 0 ? row.lo : row.hi;
      if (ex.placement.empty() || il_peak < 1e-12) {
        out.fill(0.0);
        continue;
      }
      LoadCase lc = placement_case(model, core.spans, ex.placement, "SURCH");
      if (lc.loads.empty()) {
        out.fill(0.0);
        continue;
      }
      out = read_state(tg.fam, row, ctx->solve_case(lc, so));
    }
  }
  return env;
}

}  // namespace keta

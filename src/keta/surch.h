// ST1 SURCH — moteur de charges d'exploitation (charges roulantes et charges
// réparties entre zéros de ligne d'influence), fidèle au manuel ST1 v24
// chapitre 8 « Charges d'exploitation » (p.122-174).
//
// ---------------------------------------------------------------------------
// MÉTHODOLOGIE (manuel ST1 p.125 « méthodes utilisées »)
// ---------------------------------------------------------------------------
// Pour chaque section d'étude et pour chaque effet (effort, déplacement,
// réaction), on construit la LIGNE D'INFLUENCE de l'effet (balayage multi-RHS
// d'une charge unitaire sur une factorisation persistée = Müller-Breslau par
// Maxwell-Betti, cf. `influence.h`), puis :
//
//   * CHARGES ROULANTES (Bc, Bt, Mc, convois exceptionnels, TS, VS) — on fait
//     glisser le convoi sur la ligne d'influence. Les positions candidates sont
//     la grille de pas `pas_conv` UNION les points de rupture exacts
//     {s_k − xl_j} (station k de la LI, essieu j) : la LI étant affine par
//     morceaux, l'effet d'un train d'essieux ponctuels est affine par morceaux
//     et son extrême est ATTEINT sur cette grille — le balayage est donc EXACT
//     pour les essieux ponctuels (approché au pas pour les impacts répartis,
//     dont l'effet est quadratique par morceaux).
//     Plusieurs camions par file (MAX_CAM) sont placés INDÉPENDAMMENT les uns
//     des autres « hormis le chevauchement » (p.145) : c'est un problème de
//     sélection sous contrainte d'écart minimal LONG, résolu EXACTEMENT par
//     programmation dynamique sur les positions triées (un camion défavorable
//     n'entre jamais : l'optimum « au plus n camions » l'exclut de lui-même).
//     SENS 1 = sens croissant, 2 = décroissant (convoi miroir), 3 = les deux.
//
//   * CHARGES RÉPARTIES ENTRE ZÉROS (A(l), trottoir, UDL, ALG) — on découpe la
//     LI en zones de signe avec le SEUIL du manuel (p.136 : yi appartient à une
//     zone positive si yi > ymax·seuil, négative si yi < −ymax·seuil, avec
//     ymax = max|yi|), et on ne charge que les zones favorables. La densité
//     dépendant de la LONGUEUR CHARGÉE (A(l) = 0,23 + 36/(12+L)), le choix des
//     zones et la densité sont couplés : on résout ce point fixe EXACTEMENT en
//     énumérant les sous-ensembles de zones favorables (≤ 2^12 ; au-delà, repli
//     glouton par densité d'effet décroissante) et en retenant celui qui
//     maximise q(L)·Σ∫η — q étant décroissante, charger « tout » n'est pas
//     toujours l'optimum, et l'énumération le tranche sans hypothèse.
//
//   * TRANSVERSALEMENT (p.125) : « étude indépendante de chaque file chargeable
//     et cumul de l'effet de chaque file avec pondération de dégressivité
//     transversale réglementaire ». On calcule donc un FACTEUR TRANSVERSAL
//     scalaire (Σ des coefficients de répartition REP_TRANS des voies chargées
//     × le coefficient réglementaire a1/bc/bt/COEF TRANS du nombre de voies
//     chargées), le nombre de voies chargées étant balayé de 1 à n pour retenir
//     le plus défavorable (les coefficients étant dégressifs, l'optimum n'est
//     pas toujours « toutes les voies »).
//
//   * CONCOMITANCE (p.137) : pour le placement gouvernant de la composante
//     privilégiée, on RECONSTRUIT le cas de charge réel (charges ponctuelles /
//     réparties sur les barres du tablier) et on RÉSOUT la statique — l'état
//     compagnon complet (efforts, déplacements, réactions, contraintes) est
//     donc exact et cohérent, jamais un mélange d'extrêmes.
//
// ---------------------------------------------------------------------------
// UNITÉS
// ---------------------------------------------------------------------------
// Le CATALOGUE réglementaire français est en TONNES (manuel p.123 : « Les
// charges prédéfinies du règlement français sont exprimées en tonnes »). Le
// moteur, lui, reste en unités de FORCE SI (N) : la conversion se fait par le
// paramètre `SurchDef::tonne_force` (par défaut kTonneForce = 1000·g =
// 9806,65 N/t ; l'import ST1 passera 9,81e3 pour reproduire le POND 9.81 du
// manuel, ou 1.0 pour rester en tonnes). Les charges Eurocode sont en kN et
// sont fournies directement en N dans `EurocodeParams` (p.163).
// Longueurs en m, densités surfaciques du catalogue en t/m², densités
// linéiques résultantes en N/m.
//
// ---------------------------------------------------------------------------
// LIMITES CONNUES (jamais silencieuses — rejet typé)
// ---------------------------------------------------------------------------
//  * Les composantes MOMENT de DIR GLO/LOC (ALG/CONV, ex. « DIR LOC MX 1.75 »)
//    ne sont pas portées : le moteur de lignes d'influence n'applique que des
//    charges de FORCE unitaires. Rejet typé, avec le contournement documenté
//    (excentrement à traiter par REP_TRANS).
//  * L'enveloppe des DÉPLACEMENTS porte sur les déplacements NODAUX (le manuel
//    ne l'autorise qu'aux sections d'étude de barre) — divergence documentée.
//  * L'étude transversale « n camions de front » (2ᵉ étape p.125, nécessaire à
//    la torsion) n'est pas faite : le facteur transversal est scalaire.
//  * LDIF diffuse les impacts longitudinalement (impact + 2·ldif) et
//    transversalement (moyenne du coefficient REP_TRANS sur l'impact diffusé,
//    borné à la largeur du tablier) ; conformément au manuel (p.126) la
//    diffusion ne s'applique JAMAIS aux charges réparties.
//  * LM3 : un seul véhicule VS (conforme p.167), accompagnement LM1 avec
//    longueur d'exclusion LGEXCLU appliquée à l'UDL ET au TS.
#pragma once
#include "keta/influence.h"

namespace keta {

// 1 tonne-force [N] — facteur de conversion par défaut du catalogue français.
inline constexpr double kTonneForce = 1000.0 * kGravity;  // 9806.65 N

// ===========================================================================
// TABLIER (manuel p.126-135)
// ===========================================================================

// Nombre de voies par défaut d'une zone transversale.
//   Fasc. 61 titre II (TABLIER sans CLASSE TRAFIC, p.126) :
//     0 <= larg < 5,00 -> 1 ; 5,00 <= larg <= 6,00 -> 2 ; sinon larg/3.
//   Eurocode (CLASSE TRAFIC, p.164) : < 5,40 -> 1 ; 5,40 <= L < 6,00 -> 2 ;
//     sinon Entier(L/3).
int st1_lane_count(double width, bool eurocode = false);

struct TransZone {
  double width = 0;  // LARG [m]
  int nvoies = 0;    // VOIE (0 = règle par défaut ci-dessus)
};

// REP_TRANS : coefficients de répartition transversale (type Courbon) attachés
// à une liste de barres CHARGÉES du tablier (`listekb` parmi les barres
// TABLIER, p.130 / p.134 : « coefficients pris égaux à 1 pour les barres de
// listeb non utilisées dans listekb »).
//   1 coefficient  -> constant sur tout le profil (simple pondération) ;
//   2 coefficients -> valeurs aux DEUX BORDS du profil, interpolation linéaire ;
//   n > 2          -> premier et dernier aux bords, les autres à intervalles
//                     constants (interpolation linéaire par morceaux).
// Les barres non citées ont un coefficient 1.
struct RepTransBlock {
  std::vector<Id> bars;
  std::vector<double> k;
};

struct DeckLane {
  double y = 0;      // abscisse transversale du centre de voie (axe des barres = 0)
  double width = 0;  // largeur réelle de la voie [m]
  int zone = 0;      // numéro (1-based) de la zone transversale porteuse
};

struct Tablier {
  bool trafic_class = false;  // CLASSE TRAFIC (Eurocode) au lieu de CLASSE (Fasc.61)
  int classe = 0;             // 0 = défaut (1 en Fasc.61, 2 en CLASSE TRAFIC)
  // `CLASSE TRAFIC n` when given on its own line next to a Fasc.61 `CLASSE n`
  // (p.130, p.170: a deck carrying both families, Exemple4).
  int classe_trafic = 0;
  std::vector<TransZone> zones;  // ZONE_TRANS, de gauche à droite
  double pas_ligne = 0.5;        // PAS_LIGNE [m]
  double ldif = 0.0;             // LDIF [m]
  std::vector<Id> bars;          // BARRE (chaîne géométriquement continue)
  std::vector<RepTransBlock> rep_trans;

  int bridge_class() const { return classe > 0 ? classe : (trafic_class ? 2 : 1); }
  // Class of each family: the Fasc.61 coefficients (a1, a2, bc, bt) and the
  // Eurocode adjustment are chosen independently.
  int fasc61_class() const {
    return (classe > 0 && !(trafic_class && classe_trafic == 0)) ? classe : 1;
  }
  int eurocode_class() const {
    return classe_trafic > 0 ? classe_trafic : (trafic_class && classe > 0 ? classe : 2);
  }
  double width() const;  // somme des largeurs de zones (0 si non défini)
  // Zones chargées (numéros 1-based ; vide = toutes les zones définies).
  std::vector<int> resolved_zones(const std::vector<int>& ids) const;
  double loaded_width(const std::vector<int>& ids) const;
  // Voies des zones chargées, de gauche à droite.
  std::vector<DeckLane> lanes(const std::vector<int>& ids) const;
  // Lanes under the rules of one family (Eurocode 3 m lanes, or Fasc.61
  // lanes sharing the zone), whatever the deck's CLASSE line.
  std::vector<DeckLane> lanes(const std::vector<int>& ids, bool eurocode) const;
  // Coefficient de répartition transversale vu par `study_bar` pour une charge
  // placée à l'abscisse transversale y (moyenné sur ±half_width si > 0).
  // The coefficients belong to the loaded deck bars: when every TABLIER bar
  // carries the same profile it applies to ANY studied response (arch,
  // hangers, reactions — arc deck); otherwise (one profile per span) the
  // profile of the studied deck bar is used for the whole path, a
  // non-deck study bar taking that of the first TABLIER bar
  // (`rep_trans_uniform()` false: documented approximation).
  double rep_coef(Id study_bar, double y, double half_width = 0) const;
  // The exact extreme of rep_coef(study_bar, y, half_width) over y in
  // [lo, hi] (the footprint mean of the broken line is piecewise quadratic in
  // y: its extremes sit at the ends, at the nodes +- half_width, or at the
  // vertex of a piece).
  double rep_extreme(Id study_bar, double lo, double hi, double half_width, bool maximise) const;
  bool rep_trans_uniform() const;
  void validate() const;
};

// ===========================================================================
// CATALOGUE DE CONVOIS (manuel p.145-155)
// ===========================================================================

// Loi de pondération transversale d'un convoi.
enum class TransLaw {
  None,       // aucune pondération (convois exceptionnels : « circule seul »)
  Fasc61Bc,   // table bc (p.145)
  Fasc61Bt,   // table bt (p.146)
  Table,      // COEF TRANS explicite (CONV généralisé, p.160)
};

// Un essieu (ou une empreinte répartie : impact_l > 0) d'un convoi.
// `xl` = position du CENTRE de l'empreinte par rapport à l'avant du camion
// (convention ST1 « ESSIEU j XL xl », p.160).
struct ConvoyAxle {
  double xl = 0;              // [m]
  double weight_t = 0;        // poids TOTAL de l'essieu [t]
  std::vector<double> yl;     // positions transversales des roues / chenilles [m]
  double impact_l = 0;        // dimension longitudinale de l'impact [m]
  double impact_t = 0;        // dimension transversale de l'impact [m]
};

struct ConvoyDef {
  std::string name;
  std::string source;  // citation : page du manuel + article réglementaire
  std::vector<ConvoyAxle> axles;
  double longueur = 0;  // LONG : encombrement / pas minimal en file [m]
  double largeur = 0;   // LARG : encombrement transversal [m]
  int max_cam = 1;      // MAX_CAM : camions par file
  int max_file = 0;     // MAX_FILE : files chargées (0 = nombre de voies chargées)
  TransLaw trans = TransLaw::None;
  std::vector<double> coef_trans;  // COEF TRANS (TransLaw::Table)
  bool longitudinal = false;       // forces suivant ox local (freinage)
  // Convois exceptionnels R/EG.3 : la file roule dans l'axe de sa zone, à
  // +/- axis_ecc près (0,30 m, relevé sur ST1 : Mx = 0,30 x poids pour C2 et
  // D2F1 sur une zone de 9,44 m), au lieu d'être accolée à un bord de zone.
  // < 0 : placement accolé aux bords (p.148).
  double axis_ecc = -1;
  // CONV DIR coefficients (p.164: "coefficients multiplicateurs de la charge
  // pour obtenir la composante suivant l'axe"): zero = the gravity direction.
  // GLO: the axle loads act along dir_coef, scaled by its norm; with
  // `longitudinal` (DIR LOC FX f) only dir_coef.x() is read, as the signed
  // multiplier of the force along the bars' ox.
  Vec3 dir_coef = Vec3::Zero();

  double total_t() const;
  ConvoyDef mirrored() const;  // SENS 2 : convoi retourné dans son encombrement
  void validate() const;
};

// Convoi du catalogue ST1 par mot-clé (BC, BT, MC120, CV_D3F2, VS_120_OAN…).
// Rejet typé si le nom est inconnu ; `st1_convoy_names()` liste le catalogue.
ConvoyDef st1_convoy(const std::string& name);
std::vector<std::string> st1_convoy_names();

// Coefficients réglementaires (tables du manuel).
double st1_a1(int bridge_class, int n_lanes);   // A(l), p.138
double st1_v0(int bridge_class);                // 3,50 / 3,00 / 2,75
double st1_bc(int bridge_class, int n_lanes);   // p.145
double st1_bt(int bridge_class, int n_lanes);   // p.146

// Densités du Fasc. 61 titre II, en TONNES/m² (p.138-139).
double st1_al(double loaded_length);
double st1_al_frein(double loaded_length, double loaded_surface);
double st1_al_plancher(double loaded_length);
double st1_al_plancher_frein(double loaded_length, double loaded_surface);
double st1_trot();  // 0,150 t/m²
double st1_al_pieton(double loaded_length);

// ===========================================================================
// ALG — A(l) généralisé (p.156-159)
// ===========================================================================

struct AlgDef {
  std::string name;
  // DENS : couples (longueur chargée, densité pour 1 m de largeur).
  // Interpolation linéaire ; hors intervalle, densité de la borne la plus proche.
  std::vector<std::pair<double, double>> dens;
  double larg_voie = 0;            // LARG VOIE v0 (0 = non défini)
  std::vector<double> coef_trans;  // COEF TRANS par nombre de voies chargées
  Vec3 dir = Vec3::Zero();         // DIR GLO/LOC FX FY FZ (zero = the gravity direction)
  bool dir_local = false;          // DIR LOC (repère local des barres)
  double density(double loaded_length) const;
  void validate() const;
};

// ===========================================================================
// EUROCODE (p.163-165)
// ===========================================================================

struct EurocodeParams {
  // Valeurs caractéristiques (p.163), en N et N/m².
  std::array<double, 3> Qk{{300e3, 200e3, 100e3}};  // POIDS TS, voies 1..3
  std::array<double, 4> qk{{9e3, 2.5e3, 2.5e3, 2.5e3}};  // POIDS UDL, voies 1..3 + suivantes
  double qrk = 2.5e3;                                    // aires résiduelles
  // Coefficients d'ajustement (CLASSE i AJUSTEMENT).
  double aQ1 = 1, aQi = 1, aq1 = 1, aqi = 1, aqr = 1, beta_lm2 = 1;
  // PSI FREQ (valeurs fréquentes, p.164).
  double psi_freq_ts = 0.75, psi_freq_udl = 0.4, psi_freq_lm2 = 0.75;
  // PSI COMB (accompanying values of the LM1 with an LM3, p.169-170).
  double psi_comb_ts = 0.5, psi_comb_udl = 0.3;
  // Unit of the hard-wired loads (LM2 400 kN, LM4 and sidewalk densities),
  // in N: `scaled_forces(f)` multiplies it with the tabulated ones.
  double force_scale = 1.0;

  // Classe de trafic 1, 2, 3 (p.168) et 4/5/6 (limitations de tonnage françaises).
  static EurocodeParams traffic_class(int cls);
  // The same parameters with every FORCE value (Qk, qk, qrk) multiplied by f:
  // ST1 states the predefined Eurocode loads in kN (p.127), i.e. f = 1e-3.
  EurocodeParams scaled_forces(double f) const;
  double aQ(int rank) const { return rank <= 1 ? aQ1 : aQi; }
  double aq(int rank) const { return rank <= 1 ? aq1 : aqi; }
  double Q(int rank) const { return rank >= 1 && rank <= 3 ? Qk[static_cast<std::size_t>(rank) - 1] : 0.0; }
  double q(int rank) const {
    return qk[static_cast<std::size_t>(std::min(std::max(rank, 1), 4)) - 1];
  }
};

// ===========================================================================
// SURCH (p.136, p.142, p.166)
// ===========================================================================

enum class SurchLoad {
  // charges réparties entre zéros de ligne d'influence (p.136)
  Al,
  AlFrein,
  AlPlancher,
  AlPlancherFrein,
  Trot,
  AlPieton,
  Alg,
  // charges roulantes (p.142)
  Convoy,
  // Eurocode (p.166)
  Lm1,
  Lm2,
  Lm3,
  Lm4Cte,
  Lm4Var,
  TrotEc,
  TrotLm1,
};

// Composante privilégiée (p.137) : l'étude cherche ses deux extrêmes et l'état
// concomitant complet.
enum class SurchPriv { Effort, Depla, Reac };

struct SurchDef {
  SurchLoad load = SurchLoad::Al;
  std::string convoy;  // mot-clé du catalogue (SurchLoad::Convoy / VS du LM3)
  ConvoyDef conv;      // CONV généralisé (utilisé si `convoy` est vide)
  AlgDef alg;          // ALG j
  EurocodeParams ec;

  bool freq = false;      // LM1/LM2 FREQ (sinon CARA)
  bool lm1_libre = false; // LM1 LIBRE (sinon CENTRE)
  int lm3_lm1 = 1;        // accompagnement du LM3 : 0 = NUL, 1 = FREQ, 2 = CARA, 3 = COMB
  double lgexclu = 10.0;  // LGEXCLU [m]
  double pondvs = 1.0;    // PONDVS
  // LM3 transverse placement of the VS (p.171): ZONEVS kvs (0 = not given),
  // YVS yvs1 yvs2 (ABS / REL on the zone width; yvs2 < 0 = whole zone),
  // EXC exc (offset from its lane axis), DEBORDVOIE debv (from the edge of
  // LM1 lane 1), CORRELE vj vk (the LM1 lanes it runs on, default 1 and 2).
  int zonevs = 0;
  double yvs1 = 0, yvs2 = -1;
  bool yvs_abs = false;
  double vs_exc = 0, vs_debv = 0;
  bool has_exc = false, has_debv = false;
  int correle_j = 1, correle_k = 2;

  std::vector<int> zones;  // ZONE (vide = toutes les zones du tablier)
  double pond = 1.0;       // POND
  double pas_conv = 0.5;   // PAS_CONV [m]
  // Convoys moved on the PAS_CONV grid ONLY (ST1's own scan, manual p.147),
  // instead of the grid plus the exact breakpoints of the influence line —
  // an execution choice (St1Options::pas_conv_grid_only), not a deck keyword.
  bool pas_conv_grid_only = false;
  int sens = 3;            // SENS 1 / 2 / 3
  double seuil = 0.01;     // SEUIL
  double tonne_force = kTonneForce;  // conversion catalogue (t -> N)

  SurchPriv priv = SurchPriv::Effort;
  int priv_comp = 4;  // EFFORT : 0..5 = N,Vy,Vz,T,My,Mz ; REAC/DEPLA : DDL 0..5
  // Composantes privilégiées DÉCLARÉES par famille (p.141 : EFFORT c, DEPLA c,
  // REAC c), -1 = non déclarée, c.-à-d. sans concomitance pour cette famille
  // dans l'édition ST1 (`run_surch_envelope`).
  int conc_effort = -1, conc_depla = -1, conc_reac = -1;
  void validate() const;
};

// ===========================================================================
// Résultats
// ===========================================================================

// Un groupe de véhicules identiques posés sur le chemin. Le motif porte déjà
// les charges en N, facteur transversal et POND inclus.
struct PlacedGroup {
  ConvoyPattern pattern;
  std::vector<double> positions;  // abscisses curvilignes de l'origine du motif
  double ecc = 0;  // excentrement transversal, positif à GAUCHE du sens de parcours [m]
};

// Le placement gouvernant, tel qu'il est rejoué pour la concomitance.
struct SurchPlacement {
  std::vector<PlacedGroup> groups;
  bool mirrored = false;  // convoi parcouru en sens décroissant (SENS 2)
  // Charges réparties : (s0, s1, densité linéique [N/m]) sur le chemin.
  std::vector<std::array<double, 3>> udl;
  std::vector<double> udl_ecc;  // excentrement de chaque segment de `udl` (vide = 0)
  Vec3 dir = Vec3(0, 0, -1);  // direction (globale) des charges
  bool dir_local = false;     // direction = ox local des barres
  double trans_factor = 1.0;  // facteur transversal retenu
  int n_lanes = 0;            // nombre de voies chargées retenues
  double loaded_length = 0;   // longueur chargée L convergée [m]
  double loaded_surface = 0;  // surface chargée S [m²]
  double density = 0;         // densité surfacique retenue [t/m² ou N/m²]
  bool empty() const { return groups.empty() && udl.empty(); }
};

struct SurchExtreme {
  double value = 0;
  SurchPlacement placement;
  CaseResult state;      // état concomitant complet (si opt.concomitant)
  bool has_state = false;
};

struct SurchTarget {
  ResponseKind response;
  Id bar = 0;      // section d'étude (SurchPriv::Effort)
  double x = 0;
  Id node = 0;     // SurchPriv::Reac / Depla
  int comp = 0;
};

struct SurchTargetResult {
  SurchTarget target;
  InfluenceLine il;
  SurchExtreme lo, hi;
};

struct SurchResults {
  std::vector<SurchTargetResult> sections;   // sections d'étude
  std::vector<SurchTargetResult> reactions;  // réactions d'appui
  LoadPath path;
};

struct SurchSection {
  Id bar = 0;
  double x = 0;
};

struct SurchOptions {
  // Sections d'étude (vide = `stations_per_bar` abscisses par barre du tablier).
  std::vector<SurchSection> sections;
  int stations_per_bar = 3;
  // Réactions étudiées (vide = tous les nœuds appuyés, composante `reac_comp`).
  std::vector<std::pair<Id, int>> reactions;
  int reac_comp = 2;
  // Direction of the gravity-type deck loads (unit, global). A ST1 PLANE
  // model is drawn in the x-y plane with y up: (0, -1, 0) and reac_comp 1.
  Vec3 down = Vec3(0, 0, -1);
  bool auto_reactions = true;
  // Pas de discrétisation de la ligne d'influence (0 = tablier.pas_ligne).
  double il_step = 0;
  bool concomitant = true;   // rejouer la statique au placement gouvernant
  // GRILL / SPATIALE (p.128) : « effet de l'excentrement des voies » — chaque
  // voie, file ou zone a sa ligne d'influence eta_F + y eta_T (eta_T : couple
  // unitaire n_gauche x d), étudiée séparément puis cumulée avec les
  // coefficients réglementaires. Sans effet en PLANE.
  bool eccentric_lanes = false;
  int out_stations = 11;     // stations_per_bar du re-calcul concomitant
};

SurchResults run_surch(const Model& model, const Tablier& tablier, const SurchDef& def,
                       const SurchOptions& opt = {});

// ===========================================================================
// ÉDITION ST1 D'UNE SURCHARGE (p.137, p.141) : enveloppe par section d'étude
// ===========================================================================
// « Les commandes EFFORT, DEPLA, CONTR, CONTRG j, PRESS ou REAC permettent de
// définir une composante privilégiée. L'étude de la surcharge consiste à
// chercher les deux extrêmes de cette composante et les valeurs concomitantes
// des autres composantes. Par défaut il n'y a pas de concomitance, c'est à
// dire que l'étude de la surcharge consiste à établir les extrêmes de toutes
// les composantes. »
//
// Pour chaque famille (efforts aux sections d'étude EFFORT, déplacements aux
// sections d'étude DEPLA, réactions aux appuis) :
//   * sans composante privilégiée : chaque composante a sa propre ligne
//     d'influence et ses deux extrêmes, indépendants (lo[c], hi[c]) ;
//   * avec composante privilégiée c : seule c est étudiée ; ses deux
//     placements gouvernants sont rejoués en statique et la ligne lo (hi)
//     porte TOUTES les composantes de cet état concomitant.
// Les sections d'étude sont celles de TOUTES les barres (le convoi roule sur
// le tablier, ses effets sont lus partout), comme dans les éditions de ST1.
struct SurchEnvSpec {
  std::vector<SurchSection> effort_sections;  // sections d'étude EFFORT
  std::vector<SurchSection> depla_sections;   // sections d'étude DEPLA
  std::vector<Id> support_nodes;              // appuis étudiés
  std::vector<int> effort_comps;  // composantes de station 0..5 = N Vy Vz T My Mz
  std::vector<int> depla_comps;   // composantes globales 0..5 = Dx Dy Dz Rx Ry Rz
  std::vector<int> reac_comps;    // composantes globales 0..5 = Fx Fy Fz Mx My Mz
  // Composante privilégiée de chaque famille (concomitance), -1 = aucune.
  int effort_priv = -1, depla_priv = -1, reac_priv = -1;
  // Contraintes normales (CONTR, p.66) : par section EFFORT (même ordre que
  // `effort_sections`, vide = aucune), les fibres comme combinaisons linéaires
  // des composantes de station (N Vy Vz T My Mz). Chaque fibre a sa ligne
  // d'influence — la même combinaison des leurs — et ses deux extrêmes : ce
  // sont les contraintes concomitantes, pas celles des efforts bornés.
  std::vector<std::vector<std::array<double, 6>>> effort_fibres;
  // Privileged FIBRE of each EFFORT section (SURCH CONTR / CONTRG j, p.137),
  // an index into its `effort_fibres` (-1 / absent = none): its two extremes
  // place the load, and the section's components and fibres are read off the
  // static state replayed there (concomitants) — no other effort target.
  std::vector<int> effort_priv_fibre;
  // The same on the DEPLA sections: combinations of the GLOBAL components
  // (Dx Dy Dz Rx Ry Rz) — the soil pressures of a bar on elastic ground
  // (PRESS, p.274: p = k x the local displacement) — and the privileged one
  // (SURCH PRESS <PX, PY, ...>).
  std::vector<std::vector<std::array<double, 6>>> depla_fibres;
  std::vector<int> depla_priv_fibre;
};

struct SurchEnvRow {
  Id id = 0;      // barre (efforts / déplacements) ou nœud (réactions)
  double x = 0;   // abscisse de la section [m] (barres)
  std::array<double, 6> lo{}, hi{};  // extrêmes (ou états concomitants) par composante
  std::vector<double> fibre_lo, fibre_hi;  // extrêmes des fibres (effort_fibres)
};

struct SurchEnvelope {
  std::vector<SurchEnvRow> efforts;
  std::vector<SurchEnvRow> deplas;
  std::vector<SurchEnvRow> reactions;
};

SurchEnvelope run_surch_envelope(const Model& model, const Tablier& tablier, const SurchDef& def,
                                 const SurchEnvSpec& spec, const SurchOptions& opt = {});

// Étude d'une seule réponse (brique unitaire, utile à l'import ST1).
SurchTargetResult surch_response(const Model& model, const Tablier& tablier,
                                 const SurchDef& def, const ResponseKind& response,
                                 const SurchOptions& opt = {});

// --- briques exposées (testables indépendamment) ---

// Zones de signe d'une ligne d'influence au sens du SEUIL (p.136).
struct SignZone {
  double s0 = 0, s1 = 0;   // bornes sur le chemin
  double integral = 0;     // ∫η ds sur la zone
  int sign = 0;            // +1 / -1
};
std::vector<SignZone> il_sign_zones(const InfluenceLine& il, double seuil);

// Placement optimal de `max_cam` véhicules (pas minimal `pitch`) sur la LI :
// programmation dynamique exacte sur les positions candidates.
struct ConvoyPlacement {
  double value = 0;
  std::vector<double> positions;  // abscisses de l'avant de chaque camion retenu
};
// `grid_only`: the positions of the PAS_CONV grid alone (see
// SurchDef::pas_conv_grid_only); otherwise the grid plus the exact breakpoints.
ConvoyPlacement convoy_placement(const InfluenceLine& il, const ConvoyPattern& pat,
                                 int max_cam, double pitch, double step, bool maximise,
                                 bool grid_only = false);

// Abaissement d'un ConvoyDef en motif longitudinal (impacts diffusés de LDIF).
ConvoyPattern convoy_pattern(const ConvoyDef& c, double tonne_force, double ldif);

}  // namespace keta

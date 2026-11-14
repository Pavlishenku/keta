// Prestress losses at the level of a section (EN 1992-1-1 §5.10.5.1 and
// §5.10.6), per tendon of an RcPolySection (keta/ec2_section.h).
//
// MODEL. Service state: the concrete is the section of the rings less the
// concrete its steels displace (area Ac, centroid, inertia tensor), linear
// elastic of modulus Ecm; the tendons
// are linear elastic of modulus Ep (the modulus of the section's tendon law).
// As in the resistance engine, eps_p0 is the initial strain Delta eps_p(0) of
// a tendon: its strain while the concrete at its level is unstressed. A
// bonded tendon then carries Ep (eps_p0 + eps_c), eps_c the concrete strain
// at its level — the homogenised section loaded by the tendon forces at
// decompression; an unbonded tendon keeps Ep eps_p0. An action (N, My, Mz) is
// given at the origin of the section axes, My = int sigma z, Mz = -int sigma y
// (keta's convention; My > 0 stretches the top fibre).
//
// ELASTIC SHORTENING (§5.10.5.1). The concrete stress sigma_c at a tendon
// under the prestress alone shortens a tendon that follows it by
// (Ep / Ecm) |sigma_c| — for a bonded tendon the gap between Ep eps_p0 and its
// stress. That is the whole loss of a pretensioned strand (it follows the
// transfer), and for n tendons post-tensioned one after the other its mean
// share j = (n - 1) / (2 n) (Expr. (5.44)).
//
// TIME-DEPENDENT LOSS. A tendon follows the concrete at its level. Between
// t0 and t that concrete shortens by the shrinkage and by the creep of the
// stress sigma_c,QP present under the quasi-permanent actions, and recovers
// the strain of the compression that the losses themselves remove, on the
// age-adjusted modulus Ecm / (1 + 0.8 phi); the relaxation enters reduced to
// 0.8 of its pure value. Written for tendon i (losses positive):
//
//   dsig_i + (Ep / Ecm)(1 + 0.8 phi) sum_j f_ij Ap_j dsig_j
//       = -eps_cs Ep + 0.8 dsig_pr,i - (Ep / Ecm) phi sigma_c,QP,i
//
// with f_ij the concrete stress at tendon i per unit force at tendon j,
// 1 / Ac + e_i . I^-1 . e_j. With one tendon this IS Expr. (5.46),
// f = (1 + Ac zcp^2 / Ic) / Ac. With several the same compatibility couples
// them: tendons at one point lose as the single tendon of their total area,
// and the loss of a tendon relieves the concrete at the others.
// Signs: eps_cs < 0 is a shortening (keta/creep.h), sigma_c < 0 a compression;
// a tensile sigma_c,QP gives a negative creep part (a gain).
// EN 1992-1-1 §5.10.6(2): the expression holds for bonded tendons with the
// local stresses, and for unbonded tendons with the mean stresses along the
// tendon — for those the caller passes mean values.
#pragma once
#include "keta/contour_integrals.h"
#include "keta/ec2_section.h"

namespace keta::ec2 {

// A fibre of the concrete and its stress [Pa, tension positive].
struct FibreStress {
  double y = 0, z = 0, stress = 0;
};

// Service state of the section under its prestress and `sustained` =
// (N, My, Mz): stress of each tendon and of the concrete at its level [Pa,
// tension positive], in the order of the tendons.
// The concrete is uncracked: its stress is the plane N_c / Ac + bending of the
// concrete section under the action left to it, `sustained` less
// the tendon forces — for bonded tendons, N / A + M z / I of the homogenised
// section (Ac + (Ep / Ecm) Ap) under `sustained` and the decompression forces
// Ap Ep eps_p0. `point_stress` is that stress at each of `points` = (y, z), in
// order; `fibre_max` and `fibre_min` are the extreme fibres, the most
// tensioned and the most compressed vertex of the rings (a plane over
// polygons is extreme at a vertex; along an edge at the extreme stress, one
// of its ends). A point outside the concrete — beyond the outlines or in a
// void — is a typed error.
// DISPLACED CONCRETE. The concrete is the rings less the disc each tendon
// takes away (RcPolySection::Tendon::displaced_area, its centre offset from
// the tendon): the area of an empty duct gives the net section, the tendon's
// own area Ac + (Ep / Ecm - 1) Ap. A disc counts for its area, its first
// moments and its own inertia area^2 / (4 pi) about a diameter. It must lie in
// the concrete of the rings and clear of the other discs (typed error); a
// point inside a disc is outside the concrete. A tendon without displaced
// area takes nothing away — the same void can be drawn as a ring instead.
// PASSIVE BARS. With `bars.omega` > 0 the bars of the section, of area
// omega w_i and modulus s.steel.Es, belong to the homogenised section: a bar
// follows the concrete at its level from an unstressed state, sigma_s =
// (Es / Ecm) sigma_c — a bonded tendon without initial strain and of another
// modulus. `bar_stress` is then the stress of each bar of the section, in
// order; it stays empty while the bars are left out. A bar with an initial
// strain eps0 (RcPolySection::Bar) carries Es eps0 + (Es / Ecm) sigma_c, as a
// bonded tendon does. The service state needs no tendon: with bars alone it
// is the uncracked state of a reinforced section.
struct PassiveSteel {
  double omega = 0;  // bar areas omega w_i; 0: the bars are left out of the section
  // a bar counted takes a disc of its own area of concrete away:
  // Ac + (Es / Ecm - 1) As instead of Ac + (Es / Ecm) As
  bool displace_concrete = false;
};
// The discs of concrete the steels of the section displace: one per tendon
// with a displaced area, then one per bar counted when the bars displace
// theirs. Typed errors: a negative or non-finite area or offset, a disc
// leaving the concrete of the rings or overlapping another.
using DisplacedDisc = contour::Disc;  // centre and area
std::vector<DisplacedDisc> displaced_concrete(const RcPolySection& s, const PassiveSteel& bars = {});
// Whether p lies in the concrete: in the rings (their boundaries included,
// within 1e-9 of the size of the section) and not inside a displaced disc.
bool in_concrete(const RcPolySection& s, const std::vector<DisplacedDisc>& discs, const Vec2& p);
struct PrestressServiceState {
  std::vector<double> tendon_stress, concrete_stress;
  // --- appended: the concrete away from the tendons ---
  std::vector<double> point_stress;
  FibreStress fibre_max, fibre_min;
  // --- appended: the passive bars ---
  std::vector<double> bar_stress;
  // --- appended: the strain plane ---
  // eps(y, z) = e0 + ky z - kz y of the concrete and of the steels bonded to
  // it; `decompressed`: some fibre of the rings is stretched (fibre_max in
  // tension here, where the uncracked concrete carries it).
  StrainPlane plane;
  bool decompressed = false;
};
PrestressServiceState prestress_service_state(const RcPolySection& s, double Ecm,
                                              const Vec3& sustained = Vec3::Zero(),
                                              const std::vector<Vec2>& points = {},
                                              const PassiveSteel& bars = {});

struct TendonLoss {
  double stress = 0;  // tendon stress at t0 under the prestress and the sustained action
  // (Ep / Ecm) |sigma_c| under the prestress alone: the pretensioning value,
  // an instantaneous loss that is NOT part of `total`
  double elastic_shortening = 0;
  double creep = 0, shrinkage = 0, relaxation = 0;  // parts of the time-dependent loss [Pa]
  double total = 0;                                 // their sum, Delta sigma_p,c+s+r
  double eps_p0_after = 0;  // initial strain Delta eps_p(0) of the tendon once `total` is lost
};

struct LossState {
  double Ecm = 0;
  double phi = 0;     // creep coefficient phi(t, t0)
  double eps_cs = 0;  // shrinkage strain between t0 and t, negative = shortening
  // pure relaxation loss Delta sigma_pr >= 0 of each tendon between t0 and t [Pa]; empty = none
  std::vector<double> relaxation;
  Vec3 sustained = Vec3::Zero();  // quasi-permanent (N, My, Mz) besides the prestress
  // --- appended: the passive bars and the coefficients of the method ---
  // With bars.omega > 0 the bars belong to the section at t0 (as in
  // prestress_service_state) and to the compatibility between t0 and t: a
  // bar follows the concrete at its level, so it restrains the creep and the
  // shrinkage, takes compression and unloads the concrete.
  PassiveSteel bars;
  double ageing = 0.8;            // ageing coefficient chi of the concrete (0.8 in Expr. (5.46))
  double relaxation_share = 0.8;  // share of the pure relaxation kept by the tendon (0.8 in Expr. (5.46))
};
// One entry per tendon of the section, in order. Typed errors: no tendon,
// Ecm <= 0, phi < 0, a relaxation list of another length or with a negative
// value, an ageing coefficient or a relaxation share outside [0, 1], a
// degenerate section.
std::vector<TendonLoss> prestress_losses(const RcPolySection& s, const LossState& state);

// PASSIVE BARS IN THE LOSSES. The compatibility above written for every
// steel of the section, tendon or bar, each with its own modulus E_i (the
// age-adjusted effective modulus method on the transformed section; Expr.
// (5.46) is its case of one tendon and no bar):
//
//   dsig_i + (E_i / Ecm)(1 + chi phi) sum_j f_ij A_j dsig_j
//       = -eps_cs E_i + rho dsig_pr,i - (E_i / Ecm) phi sigma_c,QP,i
//
// chi = `ageing`, rho = `relaxation_share`, dsig_pr = 0 for a bar. A bar,
// unstressed apart from the elastic state at t0, ends compressed by what it
// restrained. `tendons` is the result of prestress_losses; `bars` holds one
// entry per bar of the section when they are counted — `stress` its stress
// at t0, `creep`, `shrinkage`, `relaxation`, `total` the tension it loses
// (the compression it gains) by cause and `eps_p0_after` its initial strain
// afterwards (`elastic_shortening` unused);
// `concrete_change` is the change of the concrete stress between t0 and t,
// a plane over the section (positive where the compression decreases).
// `eps_p0_after` of a bar is its initial strain once the losses are taken
// (RcPolySection::Bar::eps0): its strain while the concrete at its level is
// unstressed. after_losses(s, SectionLosses) gives it to the bars as it gives
// the tendons theirs, so that the uncracked state of the returned section
// under the same sustained action IS the state after the losses — tendons and
// bars at their stresses less their losses, the concrete at its stress plus
// `concrete_change` — and that a later action (a live load, cracking, the
// resistance) starts from it.
struct StressPlane {
  double s0 = 0, sz = 0, sy = 0;  // sigma(y, z) = s0 + sz z + sy y [Pa]
  double at(double y, double z) const { return s0 + sz * z + sy * y; }
};
struct SectionLosses {
  std::vector<TendonLoss> tendons, bars;
  StressPlane concrete_change;
};
SectionLosses section_losses(const RcPolySection& s, const LossState& state);

// Mean share j of the elastic shortening lost by n tendons post-tensioned in
// succession, Expr. (5.44): (n - 1) / (2 n); a single tendon loses none.
double sequential_tensioning_factor(int n_tendons);

// The section once the losses are taken: each tendon at its `eps_p0_after`,
// ready for the resistance engine or for a later service state.
RcPolySection after_losses(const RcPolySection& s, const std::vector<TendonLoss>& losses);
// The same with the bars: each bar at the initial strain its entry of
// `losses.bars` carries (none when the bars were left out of the losses).
RcPolySection after_losses(const RcPolySection& s, const SectionLosses& losses);

// Losses at successive ages with the laws of EN 1992-1-1 (Annex B creep and
// shrinkage, §3.3.2 relaxation classes): each entry is the loss accumulated
// from t0, on the section as it stands at t0 — not an increment from the
// previous age. The shrinkage counts from t0 (a tendon does not see what the
// concrete shrank before it was stressed); the relaxation follows Expr.
// (3.28)-(3.30) at t - t0, without the 500 000 h cap of §3.3.2(8).
struct LossLaws {
  ConcreteTime concrete;
  double Ecm = 0;             // 0 = ec2_ecm(concrete.fck)
  double fpk = 1860e6;        // characteristic strength of the tendons (relaxation level sigma_pi / fpk)
  int relaxation_class = 2;   // §3.3.2: 1 ordinary, 2 low relaxation, 3 bars
  double rho1000 = 2.5;       // relaxation loss at 1000 h [%]
  double t0 = 28.0;           // age of the concrete at transfer [days]
  Vec3 sustained = Vec3::Zero();
  // --- appended: the passive bars in the compatibility (see LossState::bars)
  PassiveSteel bars;
};
struct LossAtAge {
  double age = 0;  // [days]
  std::vector<TendonLoss> losses;
  // --- appended: one entry per bar when the bars are counted (as SectionLosses::bars)
  std::vector<TendonLoss> bars;
};
std::vector<LossAtAge> prestress_loss_history(const RcPolySection& s, const LossLaws& laws,
                                              const std::vector<double>& ages);

}  // namespace keta::ec2

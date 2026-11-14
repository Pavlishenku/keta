// Service states of a concrete section whose concrete does not follow one
// line: the CRACKED state — the concrete linear in compression and without
// tension — and the NONLINEAR state, each material on its own law. The
// uncracked state is keta/section_losses.h (prestress_service_state).
//
// MODEL. Plane sections, eps(y, z) = e0 + ky z - kz y (tension positive). The
// concrete is the rings of the section less the discs its steels displace
// (displaced_concrete, keta/section_losses.h). A bonded tendon carries the
// strain eps_p0 + eps(y, z) from its initial strain eps_p0 (EN 1992-1-1
// §6.1(6), Fig. 6.1), an unbonded one keeps the stress of eps_p0; a bar counted
// (bars.omega > 0) carries eps(y, z). The plane is the one whose resultants
// (N, My, Mz) — My = int sigma z, Mz = -int sigma y, about the origin of the
// section axes — are `sustained`. The section needs no tendon: with bars
// alone it is a reinforced section.
//
// CRACKED STATE. sigma_c = Ecm eps where eps < 0, zero elsewhere; tendons and
// bars linear elastic (s.strand.E, s.steel.Es). Give the effective modulus
// Ecm / (1 + phi) for a sustained action. The potential of this state is
// convex: the plane is unique whenever the section can carry the action.
//
// NONLINEAR STATE. The concrete follows `concrete`, a law of keta/fiber.h:
// the parabola-rectangle (concrete_pr), the law of §3.1.5 for nonlinear
// analysis (sargin, Expr. (3.14)), the parabola-rectangle with a softening
// tension branch (concrete_tension_softening) or a line (elastic). The
// bonded tendons follow s.strand from eps_p0, the unbonded ones keep its
// stress at eps_p0, the bars follow the bilinear curve of s.steel (fyd, Es,
// ftd, eps_ud). The laws are read as given: a service state wants them with
// characteristic or mean strengths, without partial factor. A state where a
// vertex of the concrete, a bonded tendon or a bar passes the strain bound of
// its law is refused: this is a service state, not a resistance. A law that
// softens (the descending branch of Sargin past eps_c1, the tension
// softening) has no convex potential: the plane found is an equilibrium, not
// proven unique.
//
// INTEGRATION. The resultants of the concrete and their tangent are the
// contour integrals of keta/contour_integrals.h over the rings, less the same
// integrals round each displaced disc: no mesh, no clipping, exact for the
// linear and the parabolic laws, at round-off for Sargin. The plane comes
// from a damped Newton iteration started on the section at its initial moduli.
#pragma once
#include "keta/section_losses.h"

namespace keta::ec2 {

// The cracked state of the section under `sustained` = (N, My, Mz): the
// fields of PrestressServiceState, the concrete stresses being Ecm eps where
// the plane compresses and zero where it stretches. `fibre_max` and
// `fibre_min` are the vertices of the rings of largest and smallest strain,
// with the concrete stress there (zero at the most stretched vertex of a
// decompressed section). A section that the action leaves compressed
// everywhere is in its uncracked state, the same as prestress_service_state.
// Typed errors: Ecm <= 0, a point outside the concrete, an action the
// section cannot carry (no compressed concrete nor bonded steel to answer it).
PrestressServiceState cracked_service_state(const RcPolySection& s, double Ecm,
                                            const Vec3& sustained = Vec3::Zero(),
                                            const std::vector<Vec2>& points = {},
                                            const PassiveSteel& bars = {});

// The nonlinear state of the section under `sustained`, the concrete on
// `concrete`: the fields of PrestressServiceState as for the cracked state,
// the stresses read on the laws. Typed errors: a concrete law of another
// kind, a strain past the bound of a law, a point outside the concrete, an
// action for which no equilibrium is found.
PrestressServiceState nonlinear_service_state(const RcPolySection& s, const UniaxialLaw& concrete,
                                              const Vec3& sustained = Vec3::Zero(),
                                              const std::vector<Vec2>& points = {},
                                              const PassiveSteel& bars = {});

// COMPOSITE SECTION. Several materials, each on its own law over its own
// rings — a steel girder under a concrete slab, a concrete cast in two
// strengths, an encased profile — without any mesh: the resultants of each
// part are the contour integrals of its law over its rings (keta/
// contour_integrals.h). Plane sections: the parts share one plane eps(y, z) =
// e0 + ky z - kz y, as do the bonded point steels, each from its initial
// strain eps0. The rings of two parts must not overlap: a part lying inside
// another is drawn as a void of the outer one.
// Laws of a part: the concrete laws of the nonlinear state (concrete_pr,
// sargin, concrete_tension_softening), elastic, steel (bilinear with a
// plateau) and steel_hardening; of a point steel: elastic, steel,
// steel_hardening, strand.
struct MaterialPart {
  std::vector<std::vector<Vec2>> rings;  // outer boundaries counter-clockwise, voids clockwise
  UniaxialLaw law;
};
struct PointSteel {
  double y = 0, z = 0, A = 0;
  UniaxialLaw law;
  double eps0 = 0;  // initial strain: its strain while the plane is at zero
};
struct CompositeSection {
  std::vector<MaterialPart> parts;
  std::vector<PointSteel> steels;
};
// (N, My, Mz) of the section under a plane, about the origin: the sum of the
// parts and of the steels. Typed errors: no part, a law of another kind,
// non-finite data.
Vec3 composite_resultants(const CompositeSection& s, const StrainPlane& plane);
// The plane whose resultants are `action` = (N, My, Mz), by the damped
// Newton iteration of the service states started on the initial moduli;
// `part_resultants` is what each part carries and `steel_stress` the stress
// of each point steel. A state that strains a vertex of a part or a steel
// past the bound of its law is refused (crushing, rupture: a resistance, not
// a state). With a softening law the plane found is an equilibrium, not
// proven unique. Typed errors: those of composite_resultants, an action for
// which no equilibrium is found, a strain past a bound.
struct CompositeState {
  StrainPlane plane;
  std::vector<Vec3> part_resultants;
  std::vector<double> steel_stress;
};
CompositeState composite_state(const CompositeSection& s, const Vec3& action);

}  // namespace keta::ec2

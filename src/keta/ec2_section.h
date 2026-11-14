// EN 1992-1-1 ULS design DRIVERS on top of the formulas of ec2_design.h:
//
//  1. Required reinforcement of an ARBITRARY concrete section (polygon with
//     holes) under N + My + Mz, its bars laid out by the caller with relative
//     weights (As_i = omega w_i): the smallest omega whose resistance domain
//     holds (NEd, MyEd, MzEd). §6.1: plane sections, parabola-rectangle
//     concrete (§3.1.7(1) eq. 3.17, n / eps_c2 / eps_cu2 of Table 3.1 for
//     fck, no tension), the design steel law of `SteelDesign` (Figure 3.8,
//     inclined or flat top branch), and the ultimate strain planes of
//     §6.1(5) / Figure 6.1 — pivot A (the most tensioned bar at the steel
//     limit), pivot B (the most compressed fibre at -eps_cu2), pivot C (-eps_c2
//     at (1 - eps_c2 / eps_cu2) h from the most compressed fibre).
//     EXACT integrals: the concrete stress depends on one coordinate s
//     (along the strain gradient); Green's theorem turns its surface
//     integrals into contour integrals over the rings (keta/
//     contour_integrals.h), each edge cut at the kinks of the law and
//     integrated by Gauss-Legendre — exact for the parabola of exponent 2
//     and the bilinear relation, at round-off for a fractional exponent and
//     for Expr. (3.14). Bars are points (the concrete they displace is not
//     deducted, as the design charts do, unless
//     RcPolySection::deduct_displaced asks for it).
//     Root finding (strain plane at N = NEd, neutral-axis angle putting the
//     resisting moment along MEd, omega) by bracketed bisection to
//     round-off.
//
//  2. Shell (plate + membrane) reinforcement at a point: the sandwich model
//     — the outer layers carry n/2 -+ m/z (and nxy/2 -+ mxy/z) — with the
//     orthogonal reinforcement of each layer by Baumann / EN 1992-1-1 Annex F
//     (membrane_reinforcement), the transverse shear by §6.2.2 / §6.2.3
//     (v0 = sqrt(qx^2 + qy^2), VRd,c without shear reinforcement, else Asw /
//     s for cot theta). Field driver over a CaseResult's facet resultants.
//
//  3. Prestressed sections: tendons are points of area Ap on the
//     design law `strand` (keta/fiber.h UniaxialLaw::strand, §3.3.6 Fig.
//     3.10). A BONDED tendon carries its initial strain eps_p0 (§6.1(6),
//     Fig. 6.1 Delta eps_p(0): the tendon strain while the concrete around it
//     is unstrained) and follows the plane: sigma = law(eps_p0 + eps(y, z)),
//     its strain bound eps_ud entering pivot A as the extra strain
//     eps_ud - eps_p0. An UNBONDED tendon keeps the stress of its initial
//     strain plus the lump sum Delta sigma_p,ULS of §5.10.8(2) (100 MPa
//     recommended), capped by the top of its law, whatever the plane.
//
// Sign conventions: section axes (y, z), strain eps(y, z) = e0 + ky z - kz y
// (tension positive), My = int sigma z dA, Mz = -int sigma y dA (keta/fiber.h).
#pragma once
#include <optional>
#include <vector>

#include "keta/ec2_design.h"
#include "keta/fiber.h"
#include "keta/statics.h"

namespace keta {
using Vec2 = Eigen::Vector2d;  // section-plane point (y, z) [m]
}

namespace keta::ec2 {

struct RcPolySection {
  // rings in (y, z): outer boundaries counter-clockwise, holes clockwise
  std::vector<std::vector<Vec2>> rings;
  struct Bar {
    double y = 0, z = 0, w = 1;  // position and relative weight (area = omega w)
    // --- appended: initial strain ---
    // The strain of the bar while the concrete at its level is unstrained:
    // the bar carries steel(eps0 + eps(y, z)). Zero for a bar cast in an
    // unloaded section; negative once the bar has restrained the creep and
    // the shrinkage of the concrete (after_losses, keta/section_losses.h).
    double eps0 = 0;
  };
  std::vector<Bar> bars;
  double fck = 30e6, fcd = 0.85 * 30e6 / 1.5;  // fcd = alpha_cc fck / gamma_c (caller)
  SteelDesign steel;
  // --- appended: prestressing tendons ---
  struct Tendon {
    double y = 0, z = 0, Ap = 0;  // position and area [m^2] (not scaled by omega)
    double eps_p0 = 0;            // initial strain Delta eps_p(0) (Fig. 6.1)
    bool bonded = true;
    double dsigma_uls = 100e6;    // unbonded: Delta sigma_p,ULS (§5.10.8(2), recommended)
    // --- appended: the concrete the tendon takes away ---
    // A disc of concrete of this area [m^2], its centre offset by (dy, dz)
    // from the tendon (a tendon lies against the wall of its duct): the duct
    // while it is empty, the tendon's own area once it is grouted or for a
    // pretensioned strand. Read by the service states and the losses
    // (keta/section_losses.h); the resistance engine deducts no concrete
    // behind its steels, bars or tendons.
    double displaced_area = 0, displaced_dy = 0, displaced_dz = 0;
  };
  std::vector<Tendon> tendons;
  // design law of the tendons: Y1860S7 strand (Ep 195 GPa, fp0,1k 1640 MPa,
  // fpk 1860 MPa, eps_uk 3.5 %, eps_ud = 0.9 eps_uk, inclined top branch)
  UniaxialLaw strand = UniaxialLaw::strand(195e9, 1640e6, 1860e6, 0.035, 0.0315);
  // --- appended: the concrete on another law, the displaced concrete ---
  // Law of the concrete for the resistance engine. Unset: the
  // parabola-rectangle of Table 3.1 for fck at fcd. Set, it gives the pivots
  // their strains: crushing at the compressive bound of the law, pivot C at
  // the strain where the law reaches its largest compression. Kinds: concrete_pr with any exponent and
  // strains — the exponent 1 is the bilinear relation of §3.1.7(2), see
  // bilinear_concrete — and sargin (§3.1.5, read with design values as
  // §5.8.6(3) asks). A law that carries tension is refused: its tensile
  // strain bound would be a pivot of its own.
  std::optional<UniaxialLaw> concrete;
  // The resistance engine deducts the concrete its steels displace: the disc
  // of each tendon (displaced_area) and a disc of its own area omega w for
  // each bar (displaced_concrete of keta/section_losses.h). False: no
  // concrete is deducted behind the steels.
  bool deduct_displaced = false;
};

// The bilinear stress-strain relation of EN 1992-1-1 §3.1.7(2) for fck at
// fcd: a line to fcd at eps_c3, a plateau to eps_cu3 (Table 3.1: 1,75 and
// 3,5 per mille up to C50/60; above, 1,75 + 0,55 (fck - 50) / 40 and
// 2,6 + 35 ((90 - fck) / 100)^4 per mille).
UniaxialLaw bilinear_concrete(double fck, double fcd);

struct StrainPlane {
  double e0 = 0, ky = 0, kz = 0;
  double at(double y, double z) const { return e0 + ky * z - kz * y; }
};

// (N, My, Mz) of the section, bar areas omega w_i, under a strain plane.
Vec3 section_resultants(const RcPolySection& s, double omega, const StrainPlane& p);

// The ultimate strain plane of direction alpha (compression on the side of
// u = (cos alpha, sin alpha) in (y, z)) and domain parameter tau in [0, 3]
// (pivot A for [0, 1], B for [1, 2], C for [2, 3]); N decreases with tau.
// Pivot A turns about the first reinforcement to reach its bound: the bars
// at the steel limit, the bonded tendons with a strain bound at
// eps_ud - eps_p0 (the concrete edge at the steel limit without any).
StrainPlane ultimate_plane(const RcPolySection& s, double alpha, double tau);

struct SectionDesign {
  bool feasible = false;
  double omega = 0;  // bar area scale: As_i = omega w_i
  double As = 0;     // omega sum w_i
  double alpha = 0, tau = 0;
  StrainPlane plane;  // the ultimate plane carrying (NEd, MEd)
  Vec3 R = Vec3::Zero();
};
// Resistance moment along the direction theta of (My, Mz) at N = NEd: the
// ultimate plane whose N is NEd and whose moment points along theta.
// feasible = false when NEd is outside the axial capacity for this omega.
SectionDesign moment_resistance(const RcPolySection& s, double omega, double NEd, double theta);
// The smallest omega (<= omega_max) whose domain holds (NEd, MyEd, MzEd);
// feasible = false when even omega_max does not.
SectionDesign design_section(const RcPolySection& s, double NEd, double MyEd, double MzEd,
                             double omega_max = 1.0);

// ---------------------------------------------------- prestress
// Stress of a tendon of the section under a strain plane (the rule of point 3
// above); typed error for a bonded tendon strained past its law's bound or a
// law that is not a steel law (elastic, steel, steel_hardening, strand).
double tendon_stress(const RcPolySection& s, const RcPolySection::Tendon& t, const StrainPlane& p);
// The initial strain Delta eps_p(0) of a tendon whose stress is sigma_pm while
// the concrete at its level carries the stress sigma_cp (negative in
// compression, modulus Ecm): the strain of sigma_pm on the law's rising
// branch plus the concrete shortening -sigma_cp / Ecm the tendon follows once
// bonded (Fig. 6.1). sigma_cp = 0: the stress before elastic shortening.
double tendon_initial_strain(const UniaxialLaw& law, double sigma_pm, double sigma_cp = 0.0, double Ecm = 0.0);
// The isostatic action of the prestress on the section: minus the resultant
// of the tendon forces at their initial strain, (N, My, Mz) = (-P, -P z_p,
// +P y_p) summed over the tendons — what a frame analysis of a statically
// determinate member reports (keta/prestress.h: N = -P, My = -P e_z, Mz = +P e_y).
Vec3 prestress_isostatic_action(const RcPolySection& s);
// A total prestress effect at the section (from a frame analysis) split into
// its isostatic part and the hyperstatic (secondary) remainder total - iso.
struct PrestressDecomposition {
  Vec3 total = Vec3::Zero(), isostatic = Vec3::Zero(), hyperstatic = Vec3::Zero();
};
PrestressDecomposition decompose_prestress(const RcPolySection& s, const Vec3& total);

// N-M interaction diagram of direction alpha (compression on the side of
// (cos alpha, sin alpha)): the resultants (N, My, Mz) of the ultimate planes
// at n_points values of tau spread over [0, 3] (pivots A, B, C), N decreasing.
std::vector<Vec3> interaction_diagram(const RcPolySection& s, double omega, double alpha, int n_points = 61);

// ---------------------------------------------------- shell reinforcement
struct ShellDesignInput {
  double h = 0.2;          // thickness
  double d = 0;            // effective depth (shear check, and z = 0.9 d when z = 0)
  double z = 0;            // lever arm of the sandwich layers
  double fck = 30e6, fcd = 0.85 * 30e6 / 1.5, gamma_c = 1.5;
  double fyd = 500e6 / 1.15;   // stress of the layer steel (Annex F / Baumann)
  double fywd = 500e6 / 1.15;  // shear links
  double cot_theta = 1.0;      // strut inclination for Asw (eq. 6.8)
  double rho_l = 0;            // longitudinal ratio entering VRd,c (eq. 6.2a)
};
struct ShellReinforcement {
  // layer reinforcement [m^2/m]: "top" = the face +z of the facet frame
  double asx_top = 0, asy_top = 0, asx_bot = 0, asy_bot = 0;
  double asw = 0;           // shear links [m^2/m^2] (0 when VRd,c holds)
  double v0 = 0, vrdc = 0;  // principal shear and VRd,c [N/m]
};
// Resultants per unit width in keta's shell convention (keta/shell.h
// ShellResultants: m > 0 stretches the face +z): membrane nx, ny, nxy [N/m],
// moments mx, my, mxy [N m/m], shears qx, qy [N/m]. Layer forces
// n_top = n / 2 + m / z, n_bot = n / 2 - m / z (and the same for nxy, mxy).
ShellReinforcement design_shell_point(double nx, double ny, double nxy, double mx, double my, double mxy,
                                      double qx, double qy, const ShellDesignInput& in);

}  // namespace keta::ec2

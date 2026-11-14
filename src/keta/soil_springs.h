// Soil-structure interaction springs: the non-linear laws of earth pressure
// on walls and of soil resistance along piles, generated from the soil
// parameters and attached as zero-length links to an anchored ground node
// (solved by solve_static_link_laws, keta/link_nl.h).
//
// --- Ground springs ---
// add_ground_spring: a new node at the structure node's position, fully
// fixed, and a zero-length link between them whose x axis is `dir`; the law
// acts on that component (d = displacement of the structure node along dir,
// F > 0 = tension: the force on the structure is -F along dir).
//
// --- Earth pressure (trilinear earth-pressure springs) ---
// A soil layer: vertical effective stress sv(y) = sv_top + gamma (top - y)
// between its levels (y = position along `up`), at-rest p0 = K0 sv, active
// pa = max(0, Ka sv - Kac c), passive pp = Kp sv + Kpc c, and a modulus of
// subgrade reaction kh = A + B (top - y)^kappa. With d the movement of the
// wall AWAY from the soil (the active direction), the pressure is
// p(d) = clamp(p0 - kh d, pa, pp): p0 at rest (F(0) != 0), pa beyond
// (p0 - pa) / kh, pp beyond (pp - p0) / kh into the soil. The pressures are
// distributed to the wall nodes CONSISTENTLY: every characteristic level
// (P0, Pa, Pp, K of a node) is the integral of the field times the node's
// linear shape function along the chain, over the part of each adjacent
// segment inside the layer (Gauss-Legendre, exact for linear profiles; the
// active cut-off split at its zero), times the width. The sum of the plateau
// forces is then the exact thrust (e.g. Ka (q H + gamma H^2 / 2) for a
// uniform layer), the at-rest forces the exact at-rest thrust. A node at a
// layer interface gets the sum of both layers' laws. Mobilisation per node:
// (P0 - Pa) / K and (Pp - P0) / K.
// Coefficients: Rankine Ka = tan^2(45 - phi/2), Kp = tan^2(45 + phi/2);
// Jaky K0 = 1 - sin phi; Coulomb (vertical wall, horizontal ground, wall
// friction delta).
//
// --- Piles (API RP 2A-WSD, 21st ed., §6.7-6.8) ---
// Laws per unit length of pile (scale by the tributary length) or for the
// tip, as functions of the pile movement d along the spring axis:
//  * p-y soft clay (Matlock 1970), static: API's own piecewise-linear table
//    p/pu = 0, 0.23, 0.33, 0.50, 0.72, 1.00 at y/yc = 0, 0.1, 0.3, 1, 3, 8
//    (flat beyond), yc = 2.5 eps50 D, pu = min(3 c + gamma' X + J c X / D,
//    9 c) D; odd-symmetric;
//  * p-y stiff clay above water (Welch & Reese 1972), static:
//    p / pu = 0.5 (y / y50)^(1/4) up to 16 y50, pu beyond; y50 = 2.5 eps50 D;
//  * p-y sand: p = A pu tanh(k X y / (A pu)), A = max(0.9, 3 - 0.8 X / D)
//    (static; 0.9 cyclic), pu = min((C1 X + C2 D) gamma' X, C3 D gamma' X),
//    C1..C3 from phi' by the wedge formulas of Reese et al. (1974) that
//    API's figure plots;
//  * t-z (clay: t/tmax = 0.30, 0.50, 0.75, 0.90, 1.00 at z/D = 0.0016,
//    0.0031, 0.0057, 0.0080, 0.0100; sand: tmax at z = 2.54 mm), flat beyond
//    — API's post-peak softening to 0.7-0.9 tmax is NOT representable (a
//    decreasing law breaks the convexity of the static solve) and is
//    left out: the peak is kept;
//  * q-z (Q/Qp = 0.25, 0.50, 0.75, 0.90, 1.00 at z/D = 0.002, 0.013, 0.042,
//    0.073, 0.10), compression only, for a spring whose axis points UP (the
//    tip moving down = d < 0 = compression).
// The curves defined by a formula are sampled ON the curve and joined by
// chords: stiff clay in n equal force steps, sand in n equal steps of
// k X y / (A pu) up to 4 (flat beyond at tanh 4 = 0.99933 A pu; chord error
// below 0.77 (4/n)^2 / 8 of A pu: 0.6 % at n = 16); the tables are exact.
#pragma once
#include <vector>

#include "keta/model.h"

namespace keta::soil {

// ---------------------------------------------------------- coefficients
double ka_rankine(double phi);
double kp_rankine(double phi);
double k0_jaky(double phi);
// Coulomb, vertical wall, horizontal ground, wall friction delta (rad)
double ka_coulomb(double phi, double delta);
double kp_coulomb(double phi, double delta);

// ---------------------------------------------------------- ground springs
Id add_ground_spring(Model& m, Id node, const Vec3& dir, const LinkLaw& law);
// Pointwise sum of two laws (merged breakpoints).
LinkLaw law_sum(const LinkLaw& a, const LinkLaw& b);
LinkLaw law_scaled(const LinkLaw& a, double f);

// ---------------------------------------------------------- earth pressure
struct EarthLayer {
  double top = 0, bottom = 0;  // levels along `up` (top > bottom)
  double sv_top = 0;           // vertical effective stress at the top (surcharge + overburden) [Pa]
  double gamma = 0;            // effective unit weight [N/m^3]
  double Ka = 0, K0 = 0, Kp = 0;
  double c = 0, Kac = 0, Kpc = 0;  // cohesion terms: pa = Ka sv - Kac c, pp = Kp sv + Kpc c
  double kh_A = 0, kh_B = 0, kh_kappa = 1;  // kh = A + B (top - y)^kappa [N/m^3]
};
struct EarthPoint {
  double p0 = 0, pa = 0, pp = 0, kh = 0;
};
EarthPoint earth_pressure_at(const EarthLayer& l, double y);
// The law of one node from its integrated levels (forces) and stiffness.
LinkLaw earth_pressure_law(double P0, double Pa, double Pp, double K);
// Distribution of the pressure fields to the nodes: Consistent (integrals
// against the shape functions, above) or Lumped — the field at the node's
// level times the length of its two half-segments inside the layer (the
// node's level clamped into the layer), the usual practice for
// earth-pressure springs on an embedded wall. Both
// give the exact thrust of linear profiles when the layer interfaces are
// nodes.
enum class EarthDistribution { Consistent, Lumped };
// Springs on a chain of wall nodes (ordered), soil on the side opposite to
// `away` (the active direction); one link per loaded node. Returns the link
// ids, in chain order (0 for a node outside every layer).
std::vector<Id> add_earth_pressure(Model& m, const std::vector<Id>& wall, const Vec3& away, const Vec3& up,
                                   const std::vector<EarthLayer>& layers, double width,
                                   EarthDistribution dist = EarthDistribution::Consistent);

// ---------------------------------------------------------- piles (API RP 2A)
// tributary length of each node of a chain (half of each adjacent segment)
std::vector<double> tributary_lengths(const Model& m, const std::vector<Id>& chain);

// ultimate lateral resistance per unit length [N/m]
double pu_clay(double c, double gamma_eff, double X, double D, double J = 0.5);
struct SandCoefficients {
  double C1 = 0, C2 = 0, C3 = 0;
};
SandCoefficients api_sand_coefficients(double phi);
double pu_sand(double phi, double gamma_eff, double X, double D);

// p-y laws per unit length of pile
LinkLaw py_soft_clay(double c, double gamma_eff, double X, double D, double eps50, double J = 0.5);
LinkLaw py_stiff_clay(double c, double gamma_eff, double X, double D, double eps50, int n = 16);
LinkLaw py_sand(double phi, double gamma_eff, double X, double D, double k, bool cyclic = false, int n = 16);
// t-z per unit length (tmax_line = unit skin friction x perimeter), q-z of
// the tip (Qp = unit end bearing x tip area)
LinkLaw tz_clay(double tmax_line, double D);
LinkLaw tz_sand(double tmax_line);
LinkLaw qz_tip(double Qp, double D);
// API alpha of clay: f = alpha c, psi = c / sv' (alpha <= 1)
double api_alpha(double c, double sv_eff);

}  // namespace keta::soil

// Fiber section engine: point-quadrature fibers bound to uniaxial laws
// (EN 1992-1-1 §3.1.7 parabola-rectangle concrete, bilinear elastic-plastic
// steel, linear elastic), plane-section state determination by full 3x3
// Newton with a sup-tangent gradient fallback (plateaus zero the tangent),
// axial solve by bracketed bisection, and moment-curvature / capacity curves.
// Strain plane: eps(y, z) = eps0 + kappa_y z - kappa_z y  (tension-positive,
// matching sigma = N/A + My z/Iy - Mz y/Iz).
#pragma once
#include <limits>
#include <map>
#include <string>
#include "keta/model.h"

namespace keta {

struct UniaxialLaw {
  enum class Kind {
    Elastic, SteelBilinear, ConcretePR, SteelElliptic,
    // --- appended ---
    ConcreteTensionSoftening,  // parabola-rectangle + linear tension softening (crack band)
    SteelHardening,            // EN 1992-1-1 Fig. 3.8 design curve, inclined top branch to eps_ud
    Sargin,                    // EN 1992-1-1 §3.1.5 Expr. (3.14), compression only
    Strand,                    // EN 1992-1-1 §3.3.6 Fig. 3.10 prestressing steel, design curve
    Bfup                       // NF P18-710 UHPFRC: bilinear compression, class T3 tension (thick elements)
  };
  Kind kind = Kind::Elastic;
  double E = 210e9;
  double fy = 355e6;         // SteelBilinear yield (symmetric)
  double eps_u = 0;          // 0 = unlimited plateau
  // ConcretePR (compression negative strains):
  double fcd = 20e6;
  double eps_c2 = 0.002, eps_cu2 = 0.0035, n_exp = 2.0;
  // SteelElliptic (EN 199x-1-2 fire law): proportional limit + ellipse to 2%.
  double fp = 0;

  double sigma(double eps) const;
  double tangent(double eps) const;
  double sup_tangent() const;

  static UniaxialLaw elastic(double E);
  static UniaxialLaw steel(double E, double fy, double eps_u = 0);
  static UniaxialLaw concrete_pr(double fcd, double eps_c2 = 0.002, double eps_cu2 = 0.0035,
                                 double n = 2.0);
  static UniaxialLaw steel_elliptic(double E, double fp, double fy);

  // --- appended ---
  // ConcreteTensionSoftening: compression as ConcretePR; tension linear to
  // ft at eps_cr = ft / E_t, then linear softening to zero at
  // eps_tu = eps_cr + 2 Gf / (ft lc) (crack band: the dissipated energy per
  // unit crack area is Gf whatever the band length lc), zero beyond.
  // SteelHardening: E eps to fy at fy / E, then the inclined branch of slope
  // `hardening` up to eps_ud (the strain bound), mirrored in compression.
  // Sargin: sigma = -fcm (k eta - eta^2) / (1 + (k - 2) eta), eta = -eps /
  // eps_c1, k = 1.05 Ecm eps_c1 / fcm, to eps_cu1; no tension.
  // Strand: Ep eps to fpd = fp01k / gamma_s, then (inclined) to fpk / gamma_s
  // at eps_uk, bounded at eps_ud, or (horizontal) the plateau fpd; mirrored.
  // Bfup: compression Ecm eps to fcd at eps_c0d = fcd / Ecm, plateau to eps_cud;
  // tension Ecm eps to ft_el at ft_el / Ecm, then linear to ft_f at eps_lim, zero beyond.
  double ft = 0, E_t = 0, eps_tu = 0;     // tension branch (softening concrete)
  double hardening = 0, eps_ud = 0;       // SteelHardening, Strand (inclined): slope, strain bound
  double fcm = 0, eps_c1 = 0, eps_cu1 = 0, k_sargin = 0;  // Sargin
  double fpd = 0;                         // Strand plateau start (design)
  double eps_cud = 0, ft_el = 0, ft_f = 0, eps_lim = 0;  // Bfup (fcd, E shared)

  static UniaxialLaw concrete_tension_softening(double fcd, double eps_c2, double eps_cu2, double n, double E_t,
                                                double ft, double Gf, double lc);
  // design curve B(b) of Fig. 3.8: design yield fy (= fyk / gamma_s), ratio
  // k = (ft / fy)_k > 1, characteristic strain at maximum force eps_uk, design
  // bound eps_ud <= eps_uk
  static UniaxialLaw steel_hardening(double E, double fy, double k, double eps_uk, double eps_ud);
  static UniaxialLaw sargin(double fcm, double Ecm, double eps_c1, double eps_cu1);
  // Table 3.1 from fck: fcm = fck + 8 MPa, Ecm = 22 (fcm / 10)^0.3 GPa,
  // eps_c1 = 0.7 fcm^0.31 permil <= 2.8, eps_cu1 = 3.5 permil (fck < 50 MPa) else
  // 2.8 + 27 ((98 - fcm) / 100)^4 permil
  static UniaxialLaw sargin_ec2(double fck);
  // EN 1992-1-1 §3.3.6(7): Ep (195 GPa strands, 205 GPa bars), fp01k, fpk,
  // gamma_s (1.15), eps_uk, eps_ud (0.9 eps_uk recommended); inclined top
  // branch or horizontal plateau (no strain bound)
  static UniaxialLaw strand(double Ep, double fp01k, double fpk, double eps_uk, double eps_ud,
                            bool inclined = true, double gamma_s = 1.15);
  // NF P18-710 design (ULS) laws, class T3 thick element: fcd = alpha_cc fck /
  // gamma_c (alpha_cc = 0.85), eps_cud = (1 + 14 fctm / (K fcm)) eps_c0d
  // (Fig. 3.202); tension to fctk_el / gamma_cf at its elastic strain, then to
  // fctfk / (gamma_cf K) at eps_u_lim (Fig. 3.206).
  static UniaxialLaw bfup_t3(double fck, double fcm, double Ecm, double fctm, double K, double fctk_el,
                             double fctfk, double eps_u_lim, double gamma_c = 1.5, double gamma_cf = 1.3,
                             double alpha_cc = 0.85);

  // The strain bounds of the law (the capacity search of the section stops
  // there): eps_min < 0 compression, eps_max > 0 tension; +-infinity = none.
  double eps_min() const;
  double eps_max() const;
};

// Laws by name (a section refers to its materials by name).
class LawRegistry {
 public:
  void add(const std::string& name, const UniaxialLaw& law);  // typed error: empty or taken name
  const UniaxialLaw& get(const std::string& name) const;     // typed error: unknown name
  bool has(const std::string& name) const { return laws_.count(name) != 0; }
  std::vector<std::string> names() const;
 private:
  std::map<std::string, UniaxialLaw> laws_;
};

struct Fiber {
  double y = 0, z = 0, A = 0;
  int law = 0;
  // --- appended ---
  // Initial strain: the strain of the fibre while the plane is at rest at its
  // level, sigma = law(eps0 + eps(y, z)) — a passive bar after the
  // time-dependent losses of the section (compressed by the creep and the
  // shrinkage it restrained), a bar stressed before it was bonded. Read by
  // FiberSection (resultants, tangent, strain bounds on eps0 + eps);
  // HystFiberSection::wrap and fire_moment_ratio refuse a non-zero value
  // (the hysteretic section and the fire laws do not carry it).
  double eps0 = 0;
};

struct FiberSection {
  std::vector<Fiber> fibers;
  std::vector<UniaxialLaw> laws;

  double area() const;
  // Resultants of a strain plane (N, My, Mz) and the tangent 3x3.
  Vec3 resultants(const Vec3& e) const;  // e = (eps0, kappa_y, kappa_z)
  Mat3 tangent(const Vec3& e) const;

  // State determination: find e with resultants(e) = target, every fibre and
  // bonded tendon within the strain bounds of its law (a law gives no stress
  // past its bound: a plane out there can balance the target on the remaining
  // steels alone, and is no state). Newton from the plane of the section on
  // its initial moduli. Throws on unreachable targets (outside the capacity
  // polytope).
  Vec3 solve_state(const Vec3& target, int max_iter = 100, double tol = 1e-8) const;
  // Axial-only solve at fixed curvatures (bracketed bisection).
  double solve_axial(double n_target, double kappa_y, double kappa_z) const;
  // Moment-curvature about y at constant N.
  struct MKappaPoint {
    double kappa = 0, moment = 0, eps0 = 0;
  };
  std::vector<MKappaPoint> m_kappa(double n_target, const std::vector<double>& kappas) const;

  // Mesh builders.
  static FiberSection rect(double b, double h, int ny, int nz, const UniaxialLaw& law);
  void add_rebar(double y, double z, double A, const UniaxialLaw& law);
  // stress of fibre i under the plane e = (eps0, kappa_y, kappa_z)
  double fiber_stress(std::size_t i, const Vec3& e) const;
  // a bar with an initial strain (Fiber::eps0), within the bounds of its law
  void add_rebar(double y, double z, double A, const UniaxialLaw& law, double eps0);

  // --- appended: prestressing tendons ---
  // A bonded tendon follows the plane from its initial strain eps_p0 (EN
  // 1992-1-1 §6.1(6), Fig. 6.1: the tendon strain while the concrete around
  // it is unstrained): sigma = law(eps_p0 + eps(y, z)). An unbonded tendon
  // keeps the force of its initial strain whatever the plane (no overstress,
  // no stiffness). They enter the resultants, the tangent, the state and
  // axial solves, m_kappa and interaction_curve (bonded: within the law's
  // strain bounds); the hysteretic wrap and the fire drivers refuse them.
  struct Tendon {
    double y = 0, z = 0, A = 0;
    int law = 0;
    double eps_p0 = 0;
    bool bonded = true;
  };
  std::vector<Tendon> tendons;
  void add_tendon(double y, double z, double A, const UniaxialLaw& law, double eps_p0, bool bonded = true);
  // stress of tendon i under the plane e = (eps0, kappa_y, kappa_z)
  double tendon_stress(std::size_t i, const Vec3& e) const;
};

}  // namespace keta

namespace keta {
// The ultimate point of bending in the direction alpha (kappa_y = k cos a,
// kappa_z = k sin a) at constant N: the largest curvature whose axial
// equilibrium keeps every fibre and bonded tendon within the strain bounds of
// its law — the first limit strain reached, where a moment-curvature curve
// ends. `eps_steel_max` also bounds the bilinear and elliptic steel fibres,
// whose plateau has no end of its own. Typed errors: N out of reach at a
// vanishing curvature; no bound ever reached (no ultimate point).
struct UltimatePoint {
  double kappa = 0, eps0 = 0;
  double my = 0, mz = 0;
};
UltimatePoint ultimate_point(const FiberSection& s, double n_target, double alpha,
                             double eps_steel_max = std::numeric_limits<double>::infinity());

// N-M interaction curve of a fiber section: the resisting moments (my, mz) of
// the ultimate points of n_alpha bending directions spread over the turn.
struct InteractionPoint {
  double alpha = 0;
  double my = 0, mz = 0;
};
std::vector<InteractionPoint> interaction_curve(const FiberSection& s, double n_target,
                                                int n_alpha = 16, double eps_steel_max = 0.01);
}  // namespace keta

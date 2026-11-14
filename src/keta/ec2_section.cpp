#include "keta/ec2_section.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "keta/contour_integrals.h"
#include "keta/section_losses.h"

namespace keta::ec2 {

namespace {
constexpr double kPiS = 3.14159265358979323846;
struct PrP {
  double n, ec2, ecu2;
};
// EN 1992-1-1 Table 3.1 (fck in Pa)
PrP pr_of(double fck) {
  const double f = fck / 1e6;
  if (f <= 50) return {2.0, 2.0e-3, 3.5e-3};
  const double r = (90 - f) / 100;
  return {1.4 + 23.4 * std::pow(r, 4), (2.0 + 0.085 * std::pow(f - 50, 0.53)) * 1e-3,
          (2.6 + 35 * std::pow(r, 4)) * 1e-3};
}

// The strains of the concrete that set the ultimate planes: crushing (pivot
// B) and the strain of pivot C, where the law reaches its largest compression.
struct PivotStrains {
  double ecu = 0, ec = 0;
};
PivotStrains pivot_strains(const RcPolySection& s) {
  if (!s.concrete) {
    const PrP pr = pr_of(s.fck);
    return {pr.ecu2, pr.ec2};
  }
  using K = UniaxialLaw::Kind;
  const UniaxialLaw& l = *s.concrete;
  require(l.kind == K::ConcretePR || l.kind == K::Sargin,
          "ec2 section: the concrete law of the resistance engine must be concrete_pr or sargin");
  const PivotStrains p = l.kind == K::Sargin ? PivotStrains{l.eps_cu1, l.eps_c1} : PivotStrains{l.eps_cu2, l.eps_c2};
  require(std::isfinite(p.ecu) && std::isfinite(p.ec) && p.ec > 0 && p.ecu >= p.ec,
          "ec2 section: concrete law with 0 < peak strain <= crushing strain");
  return p;
}

// The concrete by the contour integrals of keta/contour_integrals.h: the law
// of the section (the parabola-rectangle of Table 3.1 for fck at fcd when
// unset) over the rings less the discs its steels displace. An ultimate
// plane never strains the concrete past crushing; a strain beyond it — an ulp
// at the pivot itself, or a plane a caller hands in — reads the stress of the
// bound, not the zero a law gives past it.
Vec3 concrete_by_contours(const RcPolySection& s, double omega, const StrainPlane& p) {
  const PivotStrains pivots = pivot_strains(s);
  const PrP pr = pr_of(s.fck);
  const UniaxialLaw law = s.concrete ? *s.concrete : UniaxialLaw::concrete_pr(s.fcd, pr.ec2, pr.ecu2, pr.n);
  contour::Law read = contour::concrete_law(law, "ec2 section");
  const double bound = -pivots.ecu * (1 - 1e-14);
  read.sigma = [law, bound](double eps) { return law.sigma(std::max(eps, bound)); };
  read.tangent = nullptr;  // the resistance engine reads resultants alone
  std::vector<contour::Disc> discs;
  if (s.deduct_displaced) discs = displaced_concrete(s, {omega, true});
  return contour::response(s.rings, discs, read, Vec3(p.e0, p.ky, p.kz)).R;
}
}  // namespace

UniaxialLaw bilinear_concrete(double fck, double fcd) {
  require(std::isfinite(fck) && fck > 0 && std::isfinite(fcd) && fcd > 0, "bilinear concrete: fck and fcd > 0");
  const double f = fck / 1e6;
  const double ec3 = f <= 50 ? 1.75e-3 : (1.75 + 0.55 * (f - 50) / 40) * 1e-3;
  const double ecu3 = f <= 50 ? 3.5e-3 : (2.6 + 35 * std::pow((90 - f) / 100, 4)) * 1e-3;
  return UniaxialLaw::concrete_pr(fcd, ec3, ecu3, 1.0);
}

Vec3 section_resultants(const RcPolySection& s, double omega, const StrainPlane& p) {
  const Vec3 c = concrete_by_contours(s, omega, p);  // (N, My, Mz)
  double N = c(0), Sy = -c(2), Sz = c(1);
  for (const auto& b : s.bars) {
    const double A = omega * b.w;
    if (A == 0) continue;
    const double sg = steel_stress(s.steel, b.eps0 + p.at(b.y, b.z));
    N += sg * A;
    Sy += sg * A * b.y;
    Sz += sg * A * b.z;
  }
  for (const auto& t : s.tendons) {
    const double sg = tendon_stress(s, t, p);
    N += sg * t.Ap;
    Sy += sg * t.Ap * t.y;
    Sz += sg * t.Ap * t.z;
  }
  return Vec3(N, Sz, -Sy);
}

namespace {
void require_steel_law(const UniaxialLaw& l) {
  using K = UniaxialLaw::Kind;
  require(l.kind == K::Elastic || l.kind == K::SteelBilinear || l.kind == K::SteelHardening || l.kind == K::Strand,
          "ec2 tendon: the law must be elastic, steel, steel_hardening or strand");
}
// the largest stress of a steel law (+inf for the elastic one)
double law_top(const UniaxialLaw& l) {
  if (l.kind == UniaxialLaw::Kind::Elastic) return std::numeric_limits<double>::infinity();
  const double em = l.eps_max();
  return l.sigma(std::isfinite(em) ? em : 1.0);
}
}  // namespace

double tendon_stress(const RcPolySection& s, const RcPolySection::Tendon& t, const StrainPlane& p) {
  require_steel_law(s.strand);
  require(std::isfinite(t.eps_p0) && std::isfinite(t.Ap) && t.Ap >= 0, "ec2 tendon: finite eps_p0 and Ap >= 0");
  if (!t.bonded) {
    require(std::isfinite(t.dsigma_uls) && t.dsigma_uls >= 0, "ec2 tendon: unbonded dsigma_uls finite and >= 0");
    return std::min(s.strand.sigma(t.eps_p0) + t.dsigma_uls, law_top(s.strand));
  }
  const double e = t.eps_p0 + p.at(t.y, t.z);
  require(e <= s.strand.eps_max() * (1 + 1e-12) && e >= s.strand.eps_min() * (1 + 1e-12),
          "ec2 tendon: strain beyond the bound of the tendon law");
  return s.strand.sigma(e);
}

double tendon_initial_strain(const UniaxialLaw& law, double sigma_pm, double sigma_cp, double Ecm) {
  require_steel_law(law);
  require(std::isfinite(sigma_pm) && sigma_pm >= 0 && sigma_pm < law_top(law),
          "tendon_initial_strain: 0 <= sigma_pm below the top of the law");
  require(std::isfinite(sigma_cp) && (sigma_cp == 0 || (std::isfinite(Ecm) && Ecm > 0)),
          "tendon_initial_strain: a concrete stress needs Ecm > 0");
  // the law rises strictly up to its top: bisection on [0, eps_max]
  double lo = 0, hi = std::isfinite(law.eps_max()) ? law.eps_max() : 1.0;
  if (law.kind == UniaxialLaw::Kind::Elastic) hi = sigma_pm / law.E * 2 + 1e-300;
  for (int i = 0; i < 200 && hi - lo > 0; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (mid == lo || mid == hi) break;
    (law.sigma(mid) < sigma_pm ? lo : hi) = mid;
  }
  const double e = law.sigma(hi) - sigma_pm < sigma_pm - law.sigma(lo) ? hi : lo;
  return e + (sigma_cp == 0 ? 0.0 : -sigma_cp / Ecm);
}

Vec3 prestress_isostatic_action(const RcPolySection& s) {
  Vec3 a = Vec3::Zero();
  for (const auto& t : s.tendons) {
    require_steel_law(s.strand);
    const double P = s.strand.sigma(t.eps_p0) * t.Ap;
    a += Vec3(-P, -P * t.z, P * t.y);
  }
  return a;
}

PrestressDecomposition decompose_prestress(const RcPolySection& s, const Vec3& total) {
  PrestressDecomposition d;
  d.total = total;
  d.isostatic = prestress_isostatic_action(s);
  d.hyperstatic = total - d.isostatic;
  return d;
}

std::vector<Vec3> interaction_diagram(const RcPolySection& s, double omega, double alpha, int n_points) {
  require(n_points >= 2, "interaction_diagram: n_points >= 2");
  std::vector<Vec3> out;
  out.reserve(static_cast<std::size_t>(n_points));
  for (int k = 0; k < n_points; ++k)
    out.push_back(section_resultants(s, omega, ultimate_plane(s, alpha, 3.0 * k / (n_points - 1))));
  return out;
}

StrainPlane ultimate_plane(const RcPolySection& s, double alpha, double tau) {
  const PivotStrains pivots = pivot_strains(s);
  // compression on the +u side: eps = e_top + g (s - s_top), g < 0 ... take
  // eps(s) = A + B s with the pivots
  const Vec2 u(std::cos(alpha), std::sin(alpha));
  double sc_max = -1e300, sc_min = 1e300;
  for (const auto& r : s.rings)
    for (const Vec2& p : r) {
      sc_max = std::max(sc_max, p.dot(u));
      sc_min = std::min(sc_min, p.dot(u));
    }
  // tension bounds of the reinforcement along u: (level s_i, strain bound);
  // the concrete edge at the steel limit without bars nor bounded tendons
  struct Lim {
    double s, lim;
  };
  std::vector<Lim> lims;
  for (const auto& b : s.bars) lims.push_back({Vec2(b.y, b.z).dot(u), s.steel.limit() - b.eps0});
  for (const auto& t : s.tendons)
    if (t.bonded && std::isfinite(s.strand.eps_max())) lims.push_back({Vec2(t.y, t.z).dot(u), s.strand.eps_max() - t.eps_p0});
  if (lims.empty()) lims.push_back({sc_min, s.steel.limit()});
  double e_a = lims.front().lim;
  for (const Lim& l : lims) e_a = std::min(e_a, l.lim);
  require(e_a > 0, "ec2: the initial strain of a bar or of a bonded tendon at or past its bound");
  const double h = sc_max - sc_min, ecu = pivots.ecu, ec2 = pivots.ec;
  // the steepest slope keeping every bound for a top strain e_top <= e_a:
  // eps(s_i) = e_top - slope_i (sc_max - s_i) <= lim_i
  auto slope_of = [&](double e_top) {
    bool any = false;
    double sl = 0;
    for (const Lim& l : lims) {
      const double d = sc_max - l.s;
      if (d <= 0) continue;
      const double v = (e_top - l.lim) / d;
      sl = any ? std::max(sl, v) : v;
      any = true;
    }
    return any ? sl : 0.0;
  };
  double e_top, slope;  // eps(s) = e_top + slope (s - sc_max)
  if (tau <= 1) {
    // pivot A: the first bound reached, eps(sc_max) from e_a to -ecu2
    e_top = e_a - tau * (e_a + ecu);
    slope = slope_of(e_top);
  } else if (tau <= 2) {
    // pivot B: eps(sc_max) = -ecu2, eps(sc_min) from its value at the end of
    // A to 0
    e_top = -ecu;
    const double slope_a = slope_of(-ecu);
    const double eb_a = -ecu - slope_a * h;  // eps at sc_min at the A/B joint
    const double eb = eb_a + (tau - 1) * (0.0 - eb_a);
    slope = (e_top - eb) / h;
  } else {
    // pivot C: -ec2 at sC = sc_max - (1 - ec2/ecu2) h, eps(sc_max) from -ecu2 to -ec2
    const double t = std::min(tau, 3.0) - 2;
    e_top = -ecu + t * (ecu - ec2);
    const double dc = (1 - ec2 / ecu) * h;
    slope = dc > 0 ? (e_top + ec2) / dc : 0.0;
  }
  // eps(p) = e_top + slope (p.u - sc_max) = e0 + ky z - kz y
  StrainPlane p;
  p.e0 = e_top - slope * sc_max;
  p.ky = slope * u.y();
  p.kz = -slope * u.x();
  return p;
}

SectionDesign moment_resistance(const RcPolySection& s, double omega, double NEd, double theta) {
  SectionDesign out;
  out.omega = omega;
  double wsum = 0;
  for (const auto& b : s.bars) wsum += b.w;
  out.As = omega * wsum;
  // the plane of direction alpha at N = NEd (N decreases with tau)
  auto at_N = [&](double alpha, double& tau, Vec3& R) {
    double lo = 0, hi = 3;
    const Vec3 Rlo = section_resultants(s, omega, ultimate_plane(s, alpha, lo));
    const Vec3 Rhi = section_resultants(s, omega, ultimate_plane(s, alpha, hi));
    if (NEd > Rlo(0) || NEd < Rhi(0)) return false;
    for (int it = 0; it < 200 && hi - lo > 1e-15; ++it) {
      const double mid = 0.5 * (lo + hi);
      (section_resultants(s, omega, ultimate_plane(s, alpha, mid))(0) > NEd ? lo : hi) = mid;
    }
    tau = 0.5 * (lo + hi);
    R = section_resultants(s, omega, ultimate_plane(s, alpha, tau));
    return true;
  };
  const Vec2 dir(std::cos(theta), std::sin(theta));  // (My, Mz) direction
  // the component of M_R across dir, as alpha turns: a sign change where
  // M_R points along dir
  auto cross = [&](double alpha, double& tau, Vec3& R, bool& ok) {
    ok = at_N(alpha, tau, R);
    if (!ok) return 0.0;
    return dir.x() * R(2) - dir.y() * R(1);
  };
  const int ns = 144;
  double best_a = 0, best_c = 0, tau = 0;
  Vec3 R;
  bool ok_prev = false;
  double a_prev = 0, c_prev = 0;
  bool found = false;
  for (int k = 0; k <= ns && !found; ++k) {
    const double a = 2 * kPiS * k / ns;
    bool ok;
    const double c = cross(a, tau, R, ok);
    if (!ok) return out;  // NEd outside the axial capacity
    const double along = dir.x() * R(1) + dir.y() * R(2);
    // a sample that is itself the root (a symmetric section loaded along an
    // axis: the crossing falls on alpha = 0): the sign test below would miss
    // it, the neighbours having the same sign on both sides
    if (along > 0 && std::fabs(c) <= 1e-12 * std::hypot(R(1), R(2))) {
      a_prev = best_a = a;
      found = true;
      break;
    }
    if (ok_prev && along > 0 && ((c_prev <= 0 && c >= 0) || (c_prev >= 0 && c <= 0))) {
      best_a = a;
      best_c = c;
      found = true;
      break;
    }
    ok_prev = along > 0;
    a_prev = a;
    c_prev = c;
  }
  if (!found) return out;
  double lo = a_prev, hi = best_a, clo = c_prev;
  (void)best_c;
  for (int it = 0; it < 200 && hi - lo > 1e-15; ++it) {
    const double mid = 0.5 * (lo + hi);
    bool ok;
    const double c = cross(mid, tau, R, ok);
    if ((c <= 0) == (clo <= 0)) {
      lo = mid;
      clo = c;
    } else {
      hi = mid;
    }
  }
  out.alpha = 0.5 * (lo + hi);
  bool ok;
  cross(out.alpha, out.tau, out.R, ok);
  out.plane = ultimate_plane(s, out.alpha, out.tau);
  out.feasible = ok;
  return out;
}

SectionDesign design_section(const RcPolySection& s, double NEd, double MyEd, double MzEd, double omega_max) {
  const double M = std::hypot(MyEd, MzEd);
  const double theta = std::atan2(MzEd, MyEd);
  // capacity margin along MEd: |M_R| - |M_Ed| (negative / infeasible below)
  auto holds = [&](double omega, SectionDesign& r) {
    r = moment_resistance(s, omega, NEd, M > 0 ? theta : 0.0);
    if (!r.feasible) return false;
    return M == 0 || std::hypot(r.R(1), r.R(2)) >= M;
  };
  SectionDesign r;
  if (holds(0.0, r)) return r;
  if (!holds(omega_max, r)) {
    r.feasible = false;
    return r;
  }
  double lo = 0, hi = omega_max;
  for (int it = 0; it < 200 && hi - lo > 1e-15 * omega_max; ++it) {
    const double mid = 0.5 * (lo + hi);
    SectionDesign t;
    (holds(mid, t) ? hi : lo) = mid;
  }
  holds(hi, r);
  return r;
}

// ---------------------------------------------------- shell reinforcement
ShellReinforcement design_shell_point(double nx, double ny, double nxy, double mx, double my, double mxy,
                                      double qx, double qy, const ShellDesignInput& in) {
  ShellReinforcement r;
  const double z = in.z > 0 ? in.z : 0.9 * in.d;
  require(z > 0, "ec2 shell design: lever arm z (or d) needed");
  // sandwich layers (keta: m > 0 stretches the face +z)
  const MembraneSteel top = membrane_reinforcement(nx / 2 + mx / z, ny / 2 + my / z, nxy / 2 + mxy / z);
  const MembraneSteel bot = membrane_reinforcement(nx / 2 - mx / z, ny / 2 - my / z, nxy / 2 - mxy / z);
  r.asx_top = top.Zx / in.fyd;
  r.asy_top = top.Zy / in.fyd;
  r.asx_bot = bot.Zx / in.fyd;
  r.asy_bot = bot.Zy / in.fyd;
  // transverse shear: principal v0, eq. 6.2a per unit width, else eq. 6.8
  r.v0 = std::hypot(qx, qy);
  if (in.d > 0) {
    VrdcInput vi;
    vi.bw = 1.0;
    vi.d = in.d;
    vi.Asl = in.rho_l * in.d;
    vi.fck = in.fck;
    vi.gamma_c = in.gamma_c;
    vi.fcd = in.fcd;
    r.vrdc = vrd_c_en(vi);
    if (r.v0 > r.vrdc) r.asw = asw_s(r.v0, z, in.fywd, in.cot_theta);
  }
  return r;
}

}  // namespace keta::ec2

#include "keta/ec4.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "keta/core.h"
#include "keta/ec3.h"

namespace keta::ec4 {

namespace {
constexpr double kPi = 3.14159265358979323846;
double eps_of(double fy) { return std::sqrt(235e6 / fy); }
}  // namespace

double stud_prd_solid(double d, double fu, double hsc, double fck, double Ecm, double gamma_v) {
  require(d > 0 && fu > 0 && hsc > 0 && fck > 0 && Ecm > 0 && gamma_v > 0, "ec4: invalid stud data");
  const double r = hsc / d;
  require(r >= 3.0, "ec4: §6.6.3.1(1) needs hsc / d >= 3");
  const double alpha = r > 4.0 ? 1.0 : 0.2 * (r + 1.0);                       // (6.20), (6.21)
  const double steel = 0.8 * fu * kPi * d * d / 4.0 / gamma_v;                  // (6.18)
  const double concrete = 0.29 * alpha * d * d * std::sqrt(fck * Ecm) / gamma_v;  // (6.19) (fck, Ecm in Pa: N/m^2)
  return std::min(steel, concrete);
}

double kt_transverse(int nr, double b0, double hp, double hsc, double t_sheet, bool through_welded) {
  require(nr >= 1 && b0 > 0 && hp > 0 && hsc > hp, "ec4: invalid rib data (hsc > hp needed)");
  const double kt = 0.7 / std::sqrt(static_cast<double>(nr)) * b0 / hp * (hsc / hp - 1.0);  // (6.23)
  // Table 6.2 (d <= 20 mm through-deck welded, d = 19 / 22 mm with holes)
  const bool thin = t_sheet <= 1.0e-3;
  double kmax = 0;
  if (nr == 1)
    kmax = through_welded ? (thin ? 0.85 : 1.0) : 0.75;
  else
    kmax = through_welded ? (thin ? 0.70 : 0.80) : 0.60;
  return std::min(kt, kmax);
}

double beff(double b0, double Le, double b1, double b2) {
  return b0 + std::min(Le / 8.0, b1) + std::min(Le / 8.0, b2);  // (5.3)
}

SaggingResult sagging(const SteelI& s, const Slab& sl, double be, double Nc) {
  SaggingResult r;
  const double fyd = s.fy / s.gM0, hc = sl.hc();
  r.Ncf = 0.85 * sl.fcd() * be * hc;
  r.Npla = s.A * fyd;
  r.Mpla = s.Wpl * fyd;
  if (r.Npla <= r.Ncf) {
    // PNA in the slab
    r.pna = PnaIn::slab;
    r.x = r.Npla / (0.85 * sl.fcd() * be);
    r.Mpl = r.Npla * (s.h / 2 + sl.hs - r.x / 2);
  } else if (r.Npla - r.Ncf <= 2 * s.b * s.tf * fyd) {
    // PNA in the top flange: Npla = Ncf + 2 Nac,f
    r.pna = PnaIn::flange;
    r.x = (r.Npla - r.Ncf) / (2 * s.b * fyd);
    const double Nacf = s.b * r.x * fyd;
    r.Mpl = r.Ncf * (s.h / 2 + sl.hs - hc / 2) + 2 * Nacf * (s.h / 2 - r.x / 2);
  } else {
    // PNA in the web (the web in tension below it carries Ncf / 2 tw fyd
    // extra on each side of the mid-depth)
    r.pna = PnaIn::web;
    r.x = s.h / 2 - r.Ncf / (2 * s.tw * fyd);
    r.Mpl = r.Mpla + r.Ncf * (s.h / 2 + sl.hs - hc / 2) - r.Ncf * r.Ncf / (4 * s.tw * fyd);
  }
  r.Nc = std::min(Nc, r.Ncf);
  r.eta = r.Nc / r.Ncf;
  r.MRd = r.Mpla + (r.Mpl - r.Mpla) * r.eta;  // (6.1)
  return r;
}

ElasticSection elastic_cracked(const SteelI& s, const Rebar& r) {
  require(s.A > 0 && s.I > 0 && r.As >= 0, "ec4: elastic_cracked needs the steel A and I");
  ElasticSection e;
  const double zr = s.h / 2 + r.zs;
  e.A = s.A + r.As;
  e.z = r.As * zr / e.A;
  e.I = s.I + s.A * e.z * e.z + r.As * (zr - e.z) * (zr - e.z);
  return e;
}

ElasticSection elastic_uncracked(const SteelI& s, const Slab& sl, double be, double n, const Rebar* r) {
  require(s.A > 0 && s.I > 0 && be > 0 && n > 0 && sl.hc() > 0, "ec4: invalid elastic section data");
  const double hc = sl.hc(), Ac = be * hc / n, Ic = be * hc * hc * hc / 12.0 / n;
  const double zc = s.h / 2 + sl.hp + hc / 2;
  const double As = r ? r->As : 0.0, zr = r ? s.h / 2 + r->zs : 0.0;
  ElasticSection e;
  e.A = s.A + Ac + As;
  e.z = (Ac * zc + As * zr) / e.A;
  e.I = s.I + s.A * e.z * e.z + Ic + Ac * (zc - e.z) * (zc - e.z) + As * (zr - e.z) * (zr - e.z);
  return e;
}

double n_long_term(double n0, double psiL, double phi_t) {
  require(n0 > 0 && psiL >= 0 && phi_t >= 0, "ec4: invalid modular ratio data");
  return n0 * (1.0 + psiL * phi_t);  // (5.6)
}

HoggingResult hogging(const SteelI& s, const Rebar& rb) {
  require(s.h > 0 && s.b > 0 && s.tw > 0 && s.tf > 0 && s.A > 0 && s.Wpl > 0 && s.fy > 0,
          "ec4: invalid steel section");
  require(rb.As >= 0 && rb.fsk > 0 && rb.gamma_s > 0, "ec4: invalid reinforcement");
  HoggingResult r;
  const double fyd = s.fy / s.gM0;
  r.Ns = rb.As * rb.fsd();
  r.Npla = s.A * fyd;
  r.Mpla = s.Wpl * fyd;
  const double lever = s.h / 2 + rb.zs;  // reinforcement above the steel mid-depth
  if (r.Ns >= r.Npla) {
    // the whole steel section in compression, the reinforcement balancing it
    // (not yielded): PNA at the reinforcement
    r.pna = PnaIn::slab;
    r.x = 0;
    r.Mpl = r.Npla * lever;
  } else if (r.Ns <= r.Npla - 2 * s.b * s.tf * fyd) {
    // PNA in the web: it rises from mid-depth by a = Ns / (2 tw fyd)
    r.pna = PnaIn::web;
    const double a = r.Ns / (2 * s.tw * fyd);
    r.x = s.h / 2 - a;
    r.Mpl = r.Mpla + r.Ns * lever - r.Ns * r.Ns / (4 * s.tw * fyd);
  } else {
    // PNA in the top flange: tension over a depth x, 2 b x fyd = Npla - Ns
    r.pna = PnaIn::flange;
    r.x = (r.Npla - r.Ns) / (2 * s.b * fyd);
    r.Mpl = r.Ns * lever + 2 * s.b * r.x * fyd * (s.h / 2 - r.x / 2);
  }
  // classes (EN 1993-1-1 Table 5.2): web of clear depth c = hw - 2 r, bottom
  // flange outstand (b - tw - 2 r) / 2 in uniform compression
  const double eps = eps_of(s.fy);
  const double zt = s.h / 2 - s.tf - s.r, zb = -zt, c = zt - zb;
  require(c > 0, "ec4: web clear depth c <= 0");
  const double zpna = r.pna == PnaIn::slab ? s.h / 2 : s.h / 2 - r.x;
  r.alpha = std::clamp((std::min(zt, zpna) - zb) / c, 0.0, 1.0);
  require(r.alpha > 0, "ec4: hogging web without compression");
  const bool has_I = s.I > 0;
  if (has_I) {
    const ElasticSection e = elastic_cracked(s, rb);
    r.psi = (zt - e.z) / (zb - e.z);
  } else {
    r.psi = -1.0;
  }
  r.web_class = ec3::class_internal(c / s.tw, eps, r.alpha, r.psi);
  if (!has_I && r.web_class > 2) r.web_class = 3;
  const double cf = (s.b - s.tw - 2 * s.r) / 2;
  r.flange_class = ec3::class_outstand(cf / s.tf, eps, 1.0, 1.0, true);
  return r;
}

DeflectionResult deflection_simply_supported(const SteelI& s, const Slab& sl, const DeflectionInput& in) {
  require(in.L > 0 && in.n0 > 0 && in.phi_t >= 0 && in.phi_sh >= 0 && in.eps_cs >= 0,
          "ec4: invalid deflection data");
  require(s.E > 0 && s.I > 0, "ec4: deflection needs the steel E and I");
  DeflectionResult r;
  const double L4 = in.L * in.L * in.L * in.L, Ea = s.E;
  auto udl = [&](double w, double I) { return 5.0 * w * L4 / (384.0 * Ea * I); };
  r.nL = n_long_term(in.n0, 1.1, in.phi_t);
  r.n_sh = n_long_term(in.n0, 0.55, in.phi_sh);
  r.perm = elastic_uncracked(s, sl, in.beff, r.nL);
  r.var = elastic_uncracked(s, sl, in.beff, in.n0);
  r.sh = elastic_uncracked(s, sl, in.beff, r.n_sh);
  const double gac = in.g_a + in.g_c;
  r.d_construction = in.propped ? 0.0 : udl(gac, s.I);
  r.d_perm = udl(in.g_s + (in.propped ? gac : 0.0), r.perm.I);
  r.d_var = udl(in.q, r.var.I);
  // primary effects of shrinkage: the restraint force of the free shortening
  // of the slab, applied at its centroid to the composite section
  const double hc = sl.hc(), zc = s.h / 2 + sl.hp + hc / 2;
  r.Nsh = in.eps_cs * Ea / r.n_sh * in.beff * hc;
  r.e_sh = zc - r.sh.z;
  r.Msh = r.Nsh * r.e_sh;
  r.d_sh = r.Msh * in.L * in.L / (8.0 * Ea * r.sh.I);
  r.d_total = r.d_construction + r.d_perm + r.d_var + r.d_sh;
  return r;
}

namespace {
double eta_bound(double v) { return std::max(v, 0.4); }
}  // namespace

double eta_min(double fy, double Le, double ratio) {
  require(fy > 0 && Le > 0, "ec4: invalid minimum-connection data");
  require(ratio >= 1.0 && ratio <= 3.0,
          "ec4: §6.6.1.2(1) covers 1 <= Abf / Atf <= 3 (bottom flange not smaller than the top one)");
  const double k = 355e6 / fy;
  const double equal = Le > 25.0 ? 1.0 : eta_bound(1.0 - k * (0.75 - 0.03 * Le));    // (6.12), (6.13)
  const double three = Le > 20.0 ? 1.0 : eta_bound(1.0 - k * (0.30 - 0.015 * Le));   // (6.14), (6.15)
  return equal + (three - equal) * (ratio - 1.0) / 2.0;
}

double eta_min_sheeting(double fy, double Le, double d, double hsc, double b0, double hp, int studs_per_rib) {
  require(fy > 0 && Le > 0, "ec4: invalid minimum-connection data");
  require(std::fabs(d - 0.019) < 1e-9 && hsc >= 0.076 - 1e-12,
          "ec4: §6.6.1.2(3)(a) needs studs d = 19 mm, hsc >= 76 mm");
  require(hp <= 0.060 + 1e-12 && b0 / hp >= 2.0, "ec4: §6.6.1.2(3)(e) needs b0 / hp >= 2, hp <= 60 mm");
  require(studs_per_rib == 1, "ec4: §6.6.1.2(3)(d) needs one stud per rib");
  const double k = 355e6 / fy;
  return Le > 25.0 ? 1.0 : eta_bound(1.0 - k * (1.0 - 0.04 * Le));  // (6.16), (6.17)
}

double cracked_length(double L1, double L2) {
  require(L1 > 0 && L2 > 0, "ec4: spans must be positive");
  require(std::min(L1, L2) / std::max(L1, L2) >= 0.6,
          "ec4: 5.4.2.3(3) needs adjacent spans with Lmin / Lmax >= 0.6 (else the cracked analysis of (2))");
  return 0.15;  // the fraction of each span, each side of the support
}

double uframe_ks(double alpha, double Ea, double I2, double a, double tw, double hs, double nu) {
  require(alpha > 0 && Ea > 0 && I2 > 0 && a > 0 && tw > 0 && hs > 0, "ec4: invalid U-frame data");
  const double k1 = alpha * Ea * I2 / a;                                   // (6.9)
  const double k2 = Ea * tw * tw * tw / (4.0 * (1.0 - nu * nu) * hs);     // (6.10)
  return k1 * k2 / (k1 + k2);                                              // (6.8)
}

double table61_max_depth(double fy, bool he, bool encased) {
  int g = -1;
  const double f = fy / 1e6;
  if (std::fabs(f - 235) < 1e-9) g = 0;
  if (std::fabs(f - 275) < 1e-9) g = 1;
  if (std::fabs(f - 355) < 1e-9) g = 2;
  if (std::fabs(f - 420) < 1e-9 || std::fabs(f - 460) < 1e-9) g = 3;
  require(g >= 0, "ec4: Table 6.1 covers S235, S275, S355, S420 and S460");
  static const double t[4][4] = {{600, 550, 400, 270},     // unencased IPE
                                 {800, 700, 650, 500},     // unencased HE
                                 {800, 700, 550, 420},     // partially encased IPE
                                 {1000, 1000, 850, 650}};  // partially encased HE
  return t[(encased ? 2 : 0) + (he ? 1 : 0)][g] * 1e-3;
}

double npl_rd(const CompositeColumn& c) {
  double N = 0;
  for (const ColumnRect& r : c.parts) {
    const double A = r.b * r.h;
    N += A * (r.kind == ColumnRect::Steel ? c.fyd : r.kind == ColumnRect::Bar ? c.fsd : c.fcd_eff);
  }
  return N;
}

std::pair<double, double> plastic_nm(const CompositeColumn& c, double z0) {
  // compression above z0 (stress -f), tension below (+f, concrete 0)
  double N = 0, M = 0;
  for (const ColumnRect& r : c.parts) {
    const double f = r.kind == ColumnRect::Steel ? c.fyd : r.kind == ColumnRect::Bar ? c.fsd : c.fcd_eff;
    const double lo = r.z - r.h / 2, hi = r.z + r.h / 2;
    // part above z0: compression
    const double a0 = std::max(lo, z0), a1 = hi;
    if (a1 > a0) {
      const double A = r.b * (a1 - a0), zc = 0.5 * (a0 + a1);
      N -= f * A;
      M += f * A * zc;  // compression above the axis: positive moment
    }
    if (r.kind == ColumnRect::Concrete) continue;
    const double b0 = lo, b1 = std::min(hi, z0);
    if (b1 > b0) {
      const double A = r.b * (b1 - b0), zc = 0.5 * (b0 + b1);
      N += f * A;
      M -= f * A * zc;
    }
  }
  return {N, M};
}

double mpl_n(const CompositeColumn& c, double N) {
  double zmin = 1e300, zmax = -1e300;
  for (const ColumnRect& r : c.parts) {
    zmin = std::min(zmin, r.z - r.h / 2);
    zmax = std::max(zmax, r.z + r.h / 2);
  }
  // N(z0) increases (less compression) as the neutral axis rises
  double lo = zmin, hi = zmax;
  const double Nlo = plastic_nm(c, lo).first, Nhi = plastic_nm(c, hi).first;
  require(N >= Nlo && N <= Nhi, "ec4: axial force outside the plastic resistance of the section");
  for (int it = 0; it < 200 && hi - lo > 1e-16 * (zmax - zmin); ++it) {
    const double mid = 0.5 * (lo + hi);
    (plastic_nm(c, mid).first < N ? lo : hi) = mid;
  }
  return plastic_nm(c, 0.5 * (lo + hi)).second;
}

double ei_eff(const CompositeColumn& c, double Ke, double Ke2, double phi_t, double NG_over_N) {
  const double Ec = c.Ecm / (1.0 + NG_over_N * phi_t);  // (6.41)
  // second moments about the section's y axis (z = 0)
  double Ia = 0, Is = 0, Ic = 0;
  for (const ColumnRect& r : c.parts) {
    const double I = r.b * r.h * r.h * r.h / 12.0 + r.b * r.h * r.z * r.z;
    (r.kind == ColumnRect::Steel ? Ia : r.kind == ColumnRect::Bar ? Is : Ic) += I;
  }
  return Ke2 * (c.Ea * Ia + c.Es * Is + Ke * Ec * Ic);
}

ColumnCheck column_check(const CompositeColumn& c, double L, double NEd, double MEd1, double r_end, ec3::Curve curve,
                         double d0, double alpha_M, double phi_t, double NG_over_N) {
  require(L > 0 && NEd < 0, "ec4: column check needs L > 0 and a compression NEd < 0");
  ColumnCheck k;
  constexpr double kPiC = 3.14159265358979323846;
  // characteristic plastic resistance (material factors 1) for lambda (6.39)
  k.Npl_Rd = npl_rd(c);
  k.Npl_Rk = k.Npl_Rd;  // the caller's design strengths; pass characteristic ones for Npl,Rk
  k.Ncr = kPiC * kPiC * ei_eff(c, 0.6, 1.0, phi_t, NG_over_N) / (L * L);
  k.lambda = std::sqrt(k.Npl_Rk / k.Ncr);
  k.chi = ec3::chi(k.lambda, ec3::imperfection(curve));
  const double N = -NEd;
  k.util_N = N / (k.chi * k.Npl_Rd);
  // second-order moments (6.43, 6.44): Ncr,eff of (EI)eff,II = 0.9 (Ea Ia +
  // Es Is + 0.5 Ecm,eff Ic); the end-moment term with beta, the member
  // imperfection N e0 with beta = 1
  const double Ncr2 = kPiC * kPiC * ei_eff(c, 0.5, 0.9, phi_t, NG_over_N) / (L * L);
  require(N < Ncr2, "ec4: NEd beyond the second-order critical force");
  const double beta = std::max(0.66 + 0.44 * r_end, 0.44);
  const double kb = std::max(beta / (1.0 - N / Ncr2), 1.0), k1 = 1.0 / (1.0 - N / Ncr2);
  k.MEd2 = kb * std::fabs(MEd1) + k1 * N * L / d0;
  // 6.7.3.6: mu_d = M_pl,N,Rd / M_pl,Rd
  const double Mpl = mpl_n(c, 0.0);
  k.mu_d = std::min(1.0, mpl_n(c, NEd) / Mpl);
  k.MRd = alpha_M * k.mu_d * Mpl;
  k.util_M = k.MEd2 / k.MRd;
  return k;
}

namespace {
// Plastic moment of a Tee (flange b x tf on top, web tw x hw below it, a
// stiffener area Ar at depth dr) about its plastic neutral axis.
double tee_plastic_moment(double b, double tf, double tw, double hw, double Ar, double dr, double fy, double fyr) {
  // area (force) above depth z
  const double Af = b * tf * fy, Aw = tw * hw * fy, Fr = Ar * fyr;
  const double total = Af + Aw + Fr;
  auto above = [&](double z) {
    double f = std::min(z, tf) / tf * Af;
    if (z > tf) f += std::min(z - tf, hw) / hw * Aw;
    if (z > dr) f += Fr;
    return f;
  };
  // bisection on the depth of the neutral axis (piecewise linear, exact to round-off)
  double lo = 0, hi = tf + hw;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    (above(mid) < total / 2 ? lo : hi) = mid;
  }
  const double z = 0.5 * (lo + hi);
  // first moments of the yield stresses about z: flange [0, tf] (width b),
  // web [tf, tf + hw] (width tw), stiffener at dr
  double M = 0;
  M += z <= tf ? b * fy * (z * z / 2 + (tf - z) * (tf - z) / 2) : b * fy * (z * tf - tf * tf / 2);
  if (z <= tf) {
    M += tw * hw * fy * (tf + hw / 2 - z);
  } else {
    const double a = z - tf, c = tf + hw - z;
    M += tw * fy * (a * a / 2 + c * c / 2);
  }
  M += Fr * std::fabs(dr - z);
  return M;
}
}  // namespace

OpeningResult opening(const BeamData& bd, const Opening& o, double MEd, double VEd, double VEd_edge) {
  const SteelI& s = bd.steel;
  const Slab& sl = bd.slab;
  OpeningResult r;
  const double fyd = s.fy / s.gM0, eps = eps_of(s.fy), hw = s.hw(), hc = sl.hc();
  require(o.ho > 0 && o.ho < hw, "ec4 P355: opening depth outside the web");
  r.ho_eq = o.circular ? 0.9 * o.ho : o.ho;     // §3.1.4
  r.lo_eq = o.circular ? 0.45 * o.ho : o.lo;
  const double Aflange = (s.A - hw * s.tw) / 2;  // flange + root (tabulated area)
  auto tee = [&](double hwT) {
    TeeProps t;
    t.hw = hwT;
    t.A = Aflange + hwT * s.tw;
    if (o.stiffened) t.Ar = o.tr * std::min(o.br, 10 * eps * o.tr);
    const double dr = s.tf + hwT - o.er;
    t.z = (s.b * s.tf * s.tf / 2 + hwT * s.tw * (s.tf + hwT / 2) + t.Ar * dr) / (s.b * s.tf + hwT * s.tw + t.Ar);
    t.NRd = t.A * fyd + t.Ar * o.fyr / s.gM0;  // (5)
    return t;
  };
  const double hwt = hw / 2 - r.ho_eq / 2 - o.eo, hwb = hw - hwt - r.ho_eq;
  require(hwt > 0 && hwb > 0, "ec4 P355: no Tee web left");
  r.top = tee(hwt);
  r.bottom = tee(hwb);
  r.heff = s.h - r.top.z - r.bottom.z;
  // global bending at the opening centre
  r.beff_o = o.x < bd.Le / 4 ? 3 * bd.Le / 16 + o.x / 4 : bd.Le / 4;  // (11)
  r.beff_o = std::min(r.beff_o / 2, bd.bl) + std::min(r.beff_o / 2, bd.br);
  r.NcRd1 = 0.85 * sl.fcd() * r.beff_o * hc;  // (6)
  const double ribs = o.x / bd.conn.s;
  r.nsc = (bd.conn.floor_studs ? std::floor(ribs) : ribs) * bd.conn.nr;
  r.NcRd2 = r.nsc * bd.conn.PRd;
  r.NcRd = std::min(r.NcRd1, r.NcRd2);
  if (r.NcRd < r.top.NRd && r.NcRd < r.bottom.NRd) {
    // PNA in the top Tee (§3.2.2): the slab at its resistance, zc = hc
    r.pna_in_top_tee = true;
    const double zc = hc;
    r.NbTEd = MEd / (r.heff + r.top.z + sl.hs - 0.5 * zc);  // (8)
    // the compressed top Tee web is limited to 10 eps tw (class 2)
    r.AtTcap = r.top.A - hwt * s.tw + s.tw * std::min(10 * eps * s.tw, hwt);
    r.NtTRd = r.AtTcap * fyd + r.top.Ar * o.fyr / s.gM0;
    r.MoRd = std::min(r.bottom.NRd, r.NtTRd + r.NcRd) * r.heff + r.NcRd * (r.top.z + sl.hs - 0.5 * hc);  // (9)
  } else {
    // PNA in the slab (§3.2.1): the bottom Tee at its resistance
    r.pna_in_top_tee = false;
    const double zc = std::min(hc, r.bottom.NRd / (0.85 * sl.fcd() * r.beff_o));
    r.NbTEd = MEd / (r.heff + r.top.z + sl.hs - 0.5 * zc);
    r.AtTcap = r.top.A;
    r.NtTRd = r.top.NRd;
    r.MoRd = r.bottom.NRd * (r.heff + r.top.z + sl.hs - 0.5 * zc);  // (7)
  }
  // shear (§3.3): the Tees of the REAL opening depth, EN 1993-1-1 6.2.6(3)
  // shear area of a rolled Tee, and the slab (EN 1992-1-1 6.2.2 form)
  {
    const double hvt = hw / 2 - o.ho / 2 - o.eo, hvb = hw - hvt - o.ho;
    auto av = [&](double hv) { return Aflange + hv * s.tw - s.b * s.tf + (2 * s.r + s.tw) * 0.5 * s.tf; };
    r.AvT = av(hvt);
    r.VtRd = av(hvt) * fyd / std::sqrt(3.0);
    r.VbRd = av(hvb) * fyd / std::sqrt(3.0);
    const double fck = sl.fck / 1e6;  // MPa
    const double CRdc = 0.18 / sl.gamma_c, k1 = 0.15;
    const double k = std::min(1 + std::sqrt(200.0 / (hc * 1e3)), 2.0);
    const double hs_eff = 0.75 * sl.hs, bw = s.b + 2 * hs_eff;
    const double Asl = kPi * bd.dl * bd.dl / 4 / bd.sl;  // per unit width
    const double rho = Asl / hc;
    r.sigma_cp = std::min(r.NbTEd / (r.beff_o * hc), 0.2 * sl.fck / sl.gamma_c);
    const double scp = r.sigma_cp / 1e6;
    r.VcRd_a = (CRdc * k * std::cbrt(100 * rho * fck) + k1 * scp) * 1e6 * bw * hc;          // (12)
    r.VcRd_min = (0.035 * std::pow(k, 1.5) * std::sqrt(fck) + k1 * scp) * 1e6 * bw * hc;    // (13)
    r.VcRd = std::max(r.VcRd_a, r.VcRd_min);
    r.VRd = r.VtRd + r.VbRd + r.VcRd;
  }
  // Vierendeel bending (§3.4)
  r.lo_eff = o.circular ? 0.7 * o.ho : o.lo;  // (15)
  r.mu = std::fabs(VEd) / r.VRd;
  auto class_limited = [&](double h) {
    if (o.circular || o.stiffened) return h;  // (circular: the equivalent Tee; stiffened: class 2 assumed)
    const double a = 32 * eps * s.tw / r.lo_eff;
    if (a >= 1) return h;
    return std::min(h, 10 * eps * s.tw / std::sqrt(1 - a * a));  // Table 3.2
  };
  r.hwt_v = class_limited(hwt);
  r.hwb_v = class_limited(hwb);
  // web thickness reduced by the shear (EN 1993-1-1 6.2.8(3)) above mu = 0.5
  const double rho_v = r.mu > 0.5 ? (2 * r.mu - 1) * (2 * r.mu - 1) : 0.0;
  const double tw_v = s.tw * (1 - rho_v);
  auto mtv = [&](double h, const TeeProps& t) {
    return tee_plastic_moment(s.b, s.tf, tw_v, h, t.Ar, s.tf + h - o.er, fyd, o.fyr / s.gM0);  // (19)
  };
  r.MtTVRd = mtv(r.hwt_v, r.top);
  r.MbTVRd = mtv(r.hwb_v, r.bottom);
  r.NtEd = std::max(r.NbTEd - r.NcRd, 0.0);
  // (20): M (1 - (N / NRd)^2) for a class-2 Tee; a STIFFENED Tee (no longer
  // a symmetric web outstand) interacts linearly, M (1 - N / NRd), as the
  // Appendix A worked example applies it
  auto inter = [&](double n) { return o.stiffened ? 1 - n : 1 - n * n; };
  r.MtTNVRd = std::max(r.MtTVRd * inter(r.NtEd / r.NtTRd), 0.0);
  r.MbTNVRd = std::max(r.MbTVRd * inter(r.NbTEd / r.bottom.NRd), 0.0);
  // local composite action (§3.4.6)
  const double so = r.lo_eq / bd.conn.s;
  r.nsco = (bd.conn.floor_studs ? std::floor(so) : so) * bd.conn.nr;
  r.dNcRd = r.nsco * bd.conn.PRd;
  if (o.circular)
    r.ko = 1.0;
  else
    r.ko = std::max(0.0, 1 - o.lo / ((o.stiffened && o.ko_eq61 ? 35.0 : 25.0) * (s.tf + r.hwt_v)));  // (25) / (61)
  r.MvcRd = r.dNcRd * (sl.hs + r.top.z - 0.5 * hc) * r.ko;  // (24)
  r.MVbRd = 2 * r.MtTNVRd + 2 * r.MbTNVRd + r.MvcRd;       // (14)
  r.MVbEd = std::fabs(VEd_edge) * r.lo_eq;
  return r;
}

WebPostResult web_post(const BeamData& bd, const WebPostInput& in) {
  const SteelI& s = bd.steel;
  const Slab& sl = bd.slab;
  WebPostResult r;
  const double fyd = s.fy / s.gM0;
  const double lever = in.heff + in.zt + sl.hs - 0.5 * sl.hc();
  r.VwpEd = std::fabs(in.VEd) * in.s / lever;                 // (29)
  r.dNcsRd = in.nscs * bd.conn.PRd;                            // (30)
  r.VwpRd = in.s0 * s.tw * fyd / std::sqrt(3.0);               // (36)
  r.MwpRd = in.s0 * in.s0 * s.tw / 6 * fyd;                    // (44)
  const double lambda1 = kPi * std::sqrt(s.E / s.fy);
  r.lambda = 2.5 * std::sqrt(in.s0 * in.s0 + in.ho * in.ho) / s.tw / lambda1;  // (41)
  r.chi = ec3::chi(r.lambda, ec3::imperfection(ec3::Curve::c));
  r.NwpRd = r.chi * in.s0 * s.tw * s.fy / s.gM1;               // (42)
  const double arm = (in.heff + sl.hs - 0.5 * sl.hc()) / (in.heff + 2 * in.eo + in.ho);
  r.VRdM = (2 * r.MwpRd / in.s + 4 * in.MbTNVRd / in.lo) * arm;            // (51)
  r.VRdN = (r.NwpRd * in.ho / in.s + 4 * in.MbTNVRd / in.lo) * arm;        // (53)
  return r;
}

}  // namespace keta::ec4

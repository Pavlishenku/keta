#include "keta/bef.h"

#include "keta/element.h"

#include <algorithm>
#include <cmath>
#include <memory>

namespace keta {

namespace {

// ---------------------------------------------------------------------------
// Transfer functions (note RDM ST1 §6)
// ---------------------------------------------------------------------------

// A decoupled field of the member: order 4 (bending, y'''' = -r y + s) or
// order 2 (axial / torsion, y'' = +r y + s), r = k / R >= 0.
struct FieldKind {
  int n = 4;
  double sgn = -1.0;
};
constexpr FieldKind kBending{4, -1.0};
constexpr FieldKind kSecond{2, +1.0};

// phi_p(x) = sum_m (sgn r)^m x^(p+n m) / (p+n m)!  for p >= 0 (f_p / g_p of
// the note); for p < 0, phi_p = sgn r phi_(p+n) (derivative of phi_0).
double phi(const FieldKind& f, int p, double x, double r) {
  if (p < 0) return f.sgn * r * phi(f, p + f.n, x, r);
  if (x <= 0) return p == 0 ? 1.0 : 0.0;
  double term = 1.0;
  for (int i = 1; i <= p; ++i) term *= x / i;
  double sum = term;
  const double rr = f.sgn * r;
  if (rr == 0.0) return sum;
  for (int m = 1; m < 200; ++m) {
    const int base = p + f.n * (m - 1);
    for (int i = 1; i <= f.n; ++i) term *= x / (base + i);
    term *= rr;
    sum += term;
    if (std::fabs(term) <= 1e-18 * std::fabs(sum)) break;
  }
  return sum;
}

// A particular-solution source inside a piece: amp * phi_(base - i)(xi - a)
// for the i-th derivative, active for xi >= a.
struct Src {
  double a = 0;
  double amp = 0;
  int base = 0;
  bool jump = false;  // a concentrated jump (side matters exactly at xi == a)
};

struct Piece {
  double xa = 0, xb = 0;
  double r = 0;   // k / R
  double k = 0;   // modulus (soil pressure recovery)
  double len() const { return xb - xa; }
};

// One field of the bar: its pieces, sources (bar coordinates) and imposed
// strain / curvature aux(x) = aux0 + aux1 x.
struct Field {
  FieldKind kind;
  double R = 1;  // EI / EA / GJ
  bool bending = true;
  std::vector<Piece> pieces;
};

struct Loads {
  // bar-coordinate distributed segments: s(x) linear on [a, b]
  struct Dist {
    double a, b, sa, sb;
  };
  std::vector<Dist> dist;
  struct Jump {
    double a, amp;
    int base;
  };
  std::vector<Jump> jumps;
  double aux0 = 0, aux1 = 0;
};

// Sources of one piece (local coordinate xi = x - xa).
std::vector<Src> piece_sources(const Field& f, std::size_t p, const Loads& ld) {
  const Piece& pc = f.pieces[p];
  const int n = f.kind.n;
  const bool last = p + 1 == f.pieces.size();
  std::vector<Src> out;
  for (const Loads::Dist& d : ld.dist) {
    double a = std::max(d.a, pc.xa), b = std::min(d.b, pc.xb);
    if (b <= a) continue;
    double s1 = d.b > d.a ? (d.sb - d.sa) / (d.b - d.a) : 0.0;
    double sa = d.sa + s1 * (a - d.a), sb = d.sa + s1 * (b - d.a);
    out.push_back({a - pc.xa, sa, n, false});
    if (s1 != 0.0) out.push_back({a - pc.xa, s1, n + 1, false});
    if (b < pc.xb) {
      out.push_back({b - pc.xa, -sb, n, false});
      if (s1 != 0.0) out.push_back({b - pc.xa, -s1, n + 1, false});
    }
  }
  for (const Loads::Jump& j : ld.jumps) {
    bool mine = (j.a >= pc.xa && j.a < pc.xb) || (last && j.a >= pc.xb);
    if (!mine) continue;
    out.push_back({std::min(j.a, pc.xb) - pc.xa, j.amp, j.base, true});
  }
  return out;
}

// i-th derivative of the particular solution at xi. `at_bar_end`: a
// concentrated source sitting exactly at the station is excluded (it acts on
// the joint — the statics convention of the end station).
double particular(const FieldKind& fk, const std::vector<Src>& srcs, int i, double xi, double r,
                  bool at_bar_end) {
  double v = 0;
  for (const Src& s : srcs) {
    if (xi < s.a) continue;
    if (s.jump && xi == s.a && at_bar_end) continue;
    v += s.amp * phi(fk, s.base - i, xi - s.a, r);
  }
  return v;
}

// Initial state s_0..s_(n-1) of a piece from its end displacements.
std::array<double, 4> initial_state(const Field& f, const Piece& pc, const std::vector<Src>& srcs,
                                    const double* dA, const double* dB) {
  const int n = f.kind.n, h = n / 2;
  const double l = pc.len();
  std::array<double, 4> s{};
  for (int j = 0; j < h; ++j) s[static_cast<std::size_t>(j)] = dA[j];
  Eigen::Matrix2d B = Eigen::Matrix2d::Identity();
  Eigen::Vector2d rhs = Eigen::Vector2d::Zero();
  for (int i = 0; i < h; ++i) {
    double v = dB[i] - particular(f.kind, srcs, i, l, pc.r, false);
    for (int j = 0; j < h; ++j) v -= dA[j] * phi(f.kind, j - i, l, pc.r);
    rhs(i) = v;
    for (int j = h; j < n; ++j) B(i, j - h) = phi(f.kind, j - i, l, pc.r);
  }
  Eigen::Vector2d x = h == 2 ? Eigen::Vector2d(B.lu().solve(rhs))
                             : Eigen::Vector2d(rhs(0) / B(0, 0), 0.0);
  for (int j = h; j < n; ++j) s[static_cast<std::size_t>(j)] = x(j - h);
  return s;
}

double deriv(const Field& f, const Piece& pc, const std::vector<Src>& srcs,
             const std::array<double, 4>& s, int k, double xi, bool at_bar_end) {
  double v = particular(f.kind, srcs, k, xi, pc.r, at_bar_end);
  for (int j = 0; j < f.kind.n; ++j) v += s[static_cast<std::size_t>(j)] * phi(f.kind, j - k, xi, pc.r);
  return v;
}

// Canonical end forces of a piece (DOF order y0, y0', y1, y1' for bending; y0,
// y1 for second order):
//   bending  [ R(y'''+aux')(0), -R(y''+aux)(0), -R(y'''+aux')(l), R(y''+aux)(l) ]
//   2nd ord. [ -R(y' - aux)(0),  R(y' - aux)(l) ]
Eigen::Vector4d piece_end_forces(const Field& f, const Piece& pc, const std::vector<Src>& srcs,
                                 const std::array<double, 4>& s, double aux0, double aux1) {
  const double l = pc.len();
  Eigen::Vector4d e = Eigen::Vector4d::Zero();
  if (f.kind.n == 4) {
    double ka = aux0 + aux1 * pc.xa, kb = aux0 + aux1 * pc.xb;
    // the start of the piece is its LEFT limit: a concentrated source at
    // xi = 0 (a load on a piece node) acts inside the piece — taken on both
    // ends, it would cancel out of the piece's fixed-end forces
    e(0) = f.R * (deriv(f, pc, srcs, s, 3, 0.0, true) + aux1);
    e(1) = -f.R * (deriv(f, pc, srcs, s, 2, 0.0, true) + ka);
    e(2) = -f.R * (deriv(f, pc, srcs, s, 3, l, false) + aux1);
    e(3) = f.R * (deriv(f, pc, srcs, s, 2, l, false) + kb);
  } else {
    double ea = aux0 + aux1 * pc.xa, eb = aux0 + aux1 * pc.xb;
    e(0) = -f.R * (deriv(f, pc, srcs, s, 1, 0.0, true) - ea);
    e(1) = f.R * (deriv(f, pc, srcs, s, 1, l, false) - eb);
  }
  return e;
}

// Exact piece stiffness (n x n in canonical DOF order).
Eigen::Matrix4d piece_stiffness(const Field& f, const Piece& pc) {
  const int n = f.kind.n, h = n / 2;
  Eigen::Matrix4d K = Eigen::Matrix4d::Zero();
  const std::vector<Src> none;
  for (int c = 0; c < n; ++c) {
    double dA[2] = {0, 0}, dB[2] = {0, 0};
    if (c < h)
      dA[c] = 1;
    else
      dB[c - h] = 1;
    auto s = initial_state(f, pc, none, dA, dB);
    K.col(c) = piece_end_forces(f, pc, none, s, 0.0, 0.0);
  }
  return K;
}

// Chain of pieces condensed on the two bar ends.
struct Chain {
  Field f;
  std::vector<Eigen::Matrix4d> kp;  // piece stiffness
  int h = 2, nn = 0;                // dofs per node, node count
  MatX Kii, Kie;                    // interior / interior-end blocks
  Eigen::PartialPivLU<MatX> lu;     // of Kii (when interior nodes exist)
  Eigen::Matrix4d Kc = Eigen::Matrix4d::Zero();  // condensed (end) stiffness

  int dof(int node, int c) const { return node * h + c; }

  void build() {
    const int n = f.kind.n;
    h = n / 2;
    nn = static_cast<int>(f.pieces.size()) + 1;
    MatX K = MatX::Zero(nn * h, nn * h);
    for (std::size_t p = 0; p < f.pieces.size(); ++p) {
      kp.push_back(piece_stiffness(f, f.pieces[p]));
      for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b) {
          int ga = dof(static_cast<int>(p) + a / h, a % h), gb = dof(static_cast<int>(p) + b / h, b % h);
          K(ga, gb) += kp.back()(a, b);
        }
    }
    const int ni = (nn - 2) * h;
    std::vector<int> ends;
    for (int c = 0; c < h; ++c) ends.push_back(dof(0, c));
    for (int c = 0; c < h; ++c) ends.push_back(dof(nn - 1, c));
    MatX Kee(2 * h, 2 * h);
    for (int a = 0; a < 2 * h; ++a)
      for (int b = 0; b < 2 * h; ++b) Kee(a, b) = K(ends[static_cast<std::size_t>(a)], ends[static_cast<std::size_t>(b)]);
    if (ni > 0) {
      Kii = K.block(h, h, ni, ni);
      Kie.resize(ni, 2 * h);
      for (int a = 0; a < ni; ++a)
        for (int b = 0; b < 2 * h; ++b) Kie(a, b) = K(h + a, ends[static_cast<std::size_t>(b)]);
      lu.compute(Kii);
      Kee -= Kie.transpose() * lu.solve(Kie);
    }
    Kc.setZero();
    Kc.topLeftCorner(2 * h, 2 * h) = Kee;
  }

  // Fixed-end forces of the loads on every piece, assembled on the chain.
  VecX fixed(const Loads& ld, std::vector<std::vector<Src>>& srcs) const {
    VecX F = VecX::Zero(nn * h);
    srcs.clear();
    const double zero[2] = {0, 0};
    for (std::size_t p = 0; p < f.pieces.size(); ++p) {
      srcs.push_back(piece_sources(f, p, ld));
      auto s = initial_state(f, f.pieces[p], srcs.back(), zero, zero);
      Eigen::Vector4d e = piece_end_forces(f, f.pieces[p], srcs.back(), s, ld.aux0, ld.aux1);
      for (int a = 0; a < f.kind.n; ++a) F(dof(static_cast<int>(p) + a / h, a % h)) += e(a);
    }
    return F;
  }

  // Condensed fixed-end force vector (end dofs).
  Eigen::Vector4d fixed_end(const Loads& ld) const {
    std::vector<std::vector<Src>> srcs;
    VecX F = fixed(ld, srcs);
    Eigen::Vector4d out = Eigen::Vector4d::Zero();
    for (int c = 0; c < h; ++c) {
      out(c) = F(dof(0, c));
      out(h + c) = F(dof(nn - 1, c));
    }
    if (nn > 2) {
      VecX Fi = F.segment(h, (nn - 2) * h);
      VecX corr = Kie.transpose() * lu.solve(Fi);
      for (int c = 0; c < 2 * h; ++c) out(c) -= corr(c);
    }
    return out;
  }

  // Node displacements of the chain for end displacements `ue` (canonical).
  VecX node_disp(const Loads& ld, const Eigen::Vector4d& ue, std::vector<std::vector<Src>>& srcs) const {
    VecX F = fixed(ld, srcs);
    VecX u = VecX::Zero(nn * h);
    for (int c = 0; c < h; ++c) {
      u(dof(0, c)) = ue(c);
      u(dof(nn - 1, c)) = ue(h + c);
    }
    if (nn > 2) {
      VecX Fi = F.segment(h, (nn - 2) * h);
      VecX uee(2 * h);
      for (int c = 0; c < 2 * h; ++c) uee(c) = ue(c);
      u.segment(h, (nn - 2) * h) = -lu.solve(Kie * uee + Fi);
    }
    return u;
  }
};

// The prepared member: four chains and their DOF maps.
struct BefModel {
  BefBar bar;
  Chain v, w, ax, tor;
};

std::vector<double> piece_nodes(const BefBar& b) {
  std::vector<double> nodes = {0.0};
  double x0 = 0;
  for (std::size_t z = 0; z < b.zones.size(); ++z) {
    const BefZone& zn = b.zones[z];
    double x1 = z + 1 == b.zones.size() ? b.L : std::min(std::max(zn.x1, x0), b.L);
    double len = x1 - x0;
    if (len <= 1e-12 * b.L) continue;
    double lmax = len;
    auto cap4 = [&](double k, double R) {
      if (k > 0 && R > 0) lmax = std::min(lmax, std::pow(R / k, 0.25));
    };
    auto cap2 = [&](double k, double R) {
      if (k > 0 && R > 0) lmax = std::min(lmax, std::sqrt(R / k));
    };
    cap4(zn.ky, b.EIz);
    cap4(zn.kz, b.EIy);
    cap2(zn.kx, b.EA);
    cap2(zn.kt, b.GJ);
    int m = std::max(1, static_cast<int>(std::ceil(len / lmax - 1e-9)));
    for (int i = 1; i <= m; ++i) nodes.push_back(x0 + len * i / m);
    x0 = x1;
    if (x0 >= b.L) break;
  }
  if (nodes.back() < b.L) nodes.push_back(b.L);
  return nodes;
}

double zone_modulus(const BefBar& b, double xm, int comp) {
  double x0 = 0;
  for (std::size_t z = 0; z < b.zones.size(); ++z) {
    double x1 = z + 1 == b.zones.size() ? b.L : b.zones[z].x1;
    if (xm <= x1 || z + 1 == b.zones.size()) {
      const BefZone& zn = b.zones[z];
      return comp == 0 ? zn.kx : comp == 1 ? zn.ky : comp == 2 ? zn.kz : zn.kt;
    }
    x0 = x1;
  }
  (void)x0;
  return 0;
}

std::shared_ptr<const BefModel> prepare(const BefBar& b) {
  require(b.L > 0 && std::isfinite(b.L), "BEF: longueur de barre invalide");
  require(b.EA > 0 && b.EIy > 0 && b.EIz > 0 && b.GJ > 0,
          "BEF: rigidites EA, EIy, EIz et GJ strictement positives requises");
  require(!b.zones.empty(), "BEF: au moins une zone de sol");
  auto m = std::make_shared<BefModel>();
  m->bar = b;
  std::vector<double> nodes = piece_nodes(b);
  auto make = [&](Chain& c, FieldKind k, double R, int comp) {
    c.f.kind = k;
    c.f.R = R;
    c.f.bending = k.n == 4;
    for (std::size_t i = 0; i + 1 < nodes.size(); ++i) {
      Piece p;
      p.xa = nodes[i];
      p.xb = nodes[i + 1];
      p.k = zone_modulus(b, 0.5 * (p.xa + p.xb), comp);
      p.r = p.k / R;
      c.f.pieces.push_back(p);
    }
    c.build();
  };
  make(m->v, kBending, b.EIz, 1);
  make(m->w, kBending, b.EIy, 2);
  make(m->ax, kSecond, b.EA, 0);
  make(m->tor, kSecond, b.GJ, 3);
  return m;
}

// Canonical loads of each field from a keta span record.
void field_loads(const BefBar& b, const BarSpanLoads& sp, Loads& lv, Loads& lw, Loads& la,
                 Loads& lt) {
  // concentrated dislocations (DEFOR CON, creep fields lowered as Gauss-point
  // dislocations): a jump of the field (base 0) or of its slope (base 1),
  // the homogeneous solution phi_base started at the abscissa — the same
  // sign as the prismatic recovery (v' jumps -theta_z, w' +theta_y)
  for (const KinkAt& k : sp.kinks) {
    if (k.theta_z != 0.0) lv.jumps.push_back({k.a, -k.theta_z, 1});
    if (k.dvy != 0.0) lv.jumps.push_back({k.a, k.dvy, 0});
    if (k.theta_y != 0.0) lw.jumps.push_back({k.a, k.theta_y, 1});
    if (k.dvz != 0.0) lw.jumps.push_back({k.a, k.dvz, 0});
    if (k.dux != 0.0) la.jumps.push_back({k.a, k.dux, 0});
    if (k.dthx != 0.0) lt.jumps.push_back({k.a, k.dthx, 0});
  }
  require_linear_span(sp.segs, "beam on elastic foundation");
  for (const DistSeg& d : sp.segs) {
    if (d.b <= d.a) continue;
    lv.dist.push_back({d.a, d.b, d.qa.y() / b.EIz, d.qb.y() / b.EIz});
    lw.dist.push_back({d.a, d.b, d.qa.z() / b.EIy, d.qb.z() / b.EIy});
    la.dist.push_back({d.a, d.b, -d.qa.x() / b.EA, -d.qb.x() / b.EA});
    // a distributed moment, by parts: a point moment field is equivalent to
    // the transverse load -m' plus the end forces m(b) at b, -m(a) at a (the
    // deflection is exact, the shear the equivalent one)
    if (d.ma.z() != 0 || d.mb.z() != 0) {
      const double dm = (d.mb.z() - d.ma.z()) / (d.b - d.a);
      lv.dist.push_back({d.a, d.b, -dm / b.EIz, -dm / b.EIz});
      lv.jumps.push_back({d.a, -d.ma.z() / b.EIz, 3});
      lv.jumps.push_back({d.b, d.mb.z() / b.EIz, 3});
    }
    if (d.ma.y() != 0 || d.mb.y() != 0) {  // theta_y = -dw/dx
      const double dm = (d.mb.y() - d.ma.y()) / (d.b - d.a);
      lw.dist.push_back({d.a, d.b, dm / b.EIy, dm / b.EIy});
      lw.jumps.push_back({d.a, d.ma.y() / b.EIy, 3});
      lw.jumps.push_back({d.b, -d.mb.y() / b.EIy, 3});
    }
    // a distributed torque: GJ th'' - k th = -m, the axial field's form
    if (d.ma.x() != 0 || d.mb.x() != 0) lt.dist.push_back({d.a, d.b, -d.ma.x() / b.GJ, -d.mb.x() / b.GJ});
  }
  require_no_imposed_twist(sp, "beam on elastic foundation");
  require_linear_strain(sp, "beam on elastic foundation");
  if (sp.deps0 != 0.0) la.dist.push_back({0.0, b.L, sp.deps0, sp.deps0});
  for (const PointAt& p : sp.points) {
    if (p.F.y() != 0.0) lv.jumps.push_back({p.a, p.F.y() / b.EIz, 3});
    if (p.M.z() != 0.0) lv.jumps.push_back({p.a, -p.M.z() / b.EIz, 2});
    if (p.F.z() != 0.0) lw.jumps.push_back({p.a, p.F.z() / b.EIy, 3});
    if (p.M.y() != 0.0) lw.jumps.push_back({p.a, p.M.y() / b.EIy, 2});
    if (p.F.x() != 0.0) la.jumps.push_back({p.a, -p.F.x() / b.EA, 1});
    if (p.M.x() != 0.0) lt.jumps.push_back({p.a, -p.M.x() / b.GJ, 1});
  }
  lv.aux0 = sp.kappa_z;
  lv.aux1 = sp.dkappa_z;
  lw.aux0 = -sp.kappa_y;
  lw.aux1 = -sp.dkappa_y;
  la.aux0 = sp.eps0;
  la.aux1 = sp.deps0;
}

// Canonical DOF maps (local index, sign) of the four fields.
struct DofMap4 {
  int idx[4];
  double sgn[4];
  int n;
};
constexpr DofMap4 kMapV{{1, 5, 7, 11}, {1, 1, 1, 1}, 4};
constexpr DofMap4 kMapW{{2, 4, 8, 10}, {1, -1, 1, -1}, 4};
constexpr DofMap4 kMapA{{0, 6, 0, 0}, {1, 1, 0, 0}, 2};
constexpr DofMap4 kMapT{{3, 9, 0, 0}, {1, 1, 0, 0}, 2};

Mat12 ke_of(const BefModel& m) {
  Mat12 K = Mat12::Zero();
  auto put = [&](const Chain& c, const DofMap4& d) {
    for (int a = 0; a < d.n; ++a)
      for (int b = 0; b < d.n; ++b) K(d.idx[a], d.idx[b]) += d.sgn[a] * d.sgn[b] * c.Kc(a, b);
  };
  put(m.v, kMapV);
  put(m.w, kMapW);
  put(m.ax, kMapA);
  put(m.tor, kMapT);
  return K;
}

Vec12 f0_of(const BefModel& m, const BarSpanLoads& sp) {
  Loads lv, lw, la, lt;
  field_loads(m.bar, sp, lv, lw, la, lt);
  Vec12 f0 = Vec12::Zero();
  auto put = [&](const Chain& c, const DofMap4& d, const Loads& ld) {
    Eigen::Vector4d fe = c.fixed_end(ld);
    for (int a = 0; a < d.n; ++a) f0(d.idx[a]) -= d.sgn[a] * fe(a);
  };
  put(m.v, kMapV, lv);
  put(m.w, kMapW, lw);
  put(m.ax, kMapA, la);
  put(m.tor, kMapT, lt);
  return f0;
}

// Field state (derivatives 0..3) at bar abscissa x.
struct FieldEval {
  const Chain* c;
  std::vector<std::vector<Src>> srcs;
  VecX u;
};

FieldEval eval_field(const Chain& c, const DofMap4& d, const Loads& ld, const Vec12& ul) {
  FieldEval fe;
  fe.c = &c;
  Eigen::Vector4d ue = Eigen::Vector4d::Zero();
  for (int a = 0; a < d.n; ++a) ue(a) = d.sgn[a] * ul(d.idx[a]);
  fe.u = c.node_disp(ld, ue, fe.srcs);
  return fe;
}

// Returns (y, y', y'', y''') and the piece modulus at x.
std::array<double, 5> state_at(const FieldEval& fe, double x, double L) {
  const Chain& c = *fe.c;
  const std::vector<Piece>& ps = c.f.pieces;
  std::size_t p = 0;
  while (p + 1 < ps.size() && x >= ps[p].xb) ++p;
  const Piece& pc = ps[p];
  const int h = c.h;
  double dA[2], dB[2];
  for (int j = 0; j < h; ++j) {
    dA[j] = fe.u(c.dof(static_cast<int>(p), j));
    dB[j] = fe.u(c.dof(static_cast<int>(p) + 1, j));
  }
  auto s = initial_state(c.f, pc, fe.srcs[p], dA, dB);
  double xi = std::min(std::max(x - pc.xa, 0.0), pc.len());
  bool at_end = x >= L;
  std::array<double, 5> out{};
  for (int k = 0; k < c.f.kind.n; ++k)
    out[static_cast<std::size_t>(k)] = deriv(c.f, pc, fe.srcs[p], s, k, xi, at_end);
  out[4] = pc.k;
  return out;
}

}  // namespace

Mat12 bef_ke(const BefBar& bar) { return ke_of(*prepare(bar)); }

Vec12 bef_f0(const BefBar& bar, const BarSpanLoads& sp) { return f0_of(*prepare(bar), sp); }

namespace {

void recover_with(const BefModel& m, const BarSpanLoads& sp, const Vec12& ul,
                  const std::vector<double>& xs, BarResult& out) {
  Loads lv, lw, la, lt;
  field_loads(m.bar, sp, lv, lw, la, lt);
  FieldEval ev = eval_field(m.v, kMapV, lv, ul);
  FieldEval ew = eval_field(m.w, kMapW, lw, ul);
  FieldEval ea = eval_field(m.ax, kMapA, la, ul);
  FieldEval et = eval_field(m.tor, kMapT, lt, ul);
  const BefBar& b = m.bar;
  out.stations.resize(xs.size());
  out.soil_p.assign(xs.size(), Vec3::Zero());
  out.soil_m.assign(xs.size(), 0.0);
  for (std::size_t k = 0; k < xs.size(); ++k) {
    double x = std::min(std::max(xs[k], 0.0), b.L);
    auto sv = state_at(ev, x, b.L), sw = state_at(ew, x, b.L);
    auto sa = state_at(ea, x, b.L), st = state_at(et, x, b.L);
    double kz = lv.aux0 + lv.aux1 * x, kyc = lw.aux0 + lw.aux1 * x;
    double eps = la.aux0 + la.aux1 * x;
    Station& s = out.stations[k];
    s.x = x;
    s.Mz = b.EIz * (sv[2] + kz);
    s.Vy = -b.EIz * (sv[3] + lv.aux1);
    s.My = -b.EIy * (sw[2] + kyc);
    s.Vz = -b.EIy * (sw[3] + lw.aux1);
    s.N = b.EA * (sa[1] - eps);
    s.T = b.GJ * st[1];
    s.defl = Vec3(sa[0], sv[0], sw[0]);
    s.rot = Vec3(st[0], -sw[1], sv[1]);
    out.soil_p[k] = Vec3(sa[4] * sa[0], sv[4] * sv[0], sw[4] * sw[0]);
    out.soil_m[k] = st[4] * st[0];
  }
}

}  // namespace

void install_bef(Model& model, Id bar_id, const std::vector<BefZone>& zones) {
  const Bar& b = model.bar(bar_id);
  for (bool r : b.releases)
    require(!r, "BEF: barre " + std::to_string(bar_id) + " sur sol elastique avec liberation non portee");
  require(!b.soil.any(), "BEF: barre " + std::to_string(bar_id) + " porte deja un lit Bar::soil");
  require(b.kind == BarKind::Frame, "BEF: seule une barre de portique peut reposer sur le sol");
  const Material& mat = model.material(b.material);
  const SectionProps& sec = model.section(b.section).props;
  BefBar bb;
  // Element axis length: between the OFFSET ends of an eccentric bar (the
  // rigid arms belong to the element transform, the bed to the axis).
  bb.L = BarElement(model, b).L();
  bb.EA = mat.E * sec.A;
  bb.EIy = mat.E * sec.Iy;
  bb.EIz = mat.E * sec.Iz;
  bb.GJ = mat.G() * sec.J;
  bb.zones = zones;
  std::shared_ptr<const BefModel> m = prepare(bb);
  TaperedBar tb;
  tb.ke = ke_of(*m);
  tb.f0 = [](const Load&) -> Vec12 {
    fail("BEF: charges abaissees par f0_span uniquement");
  };
  tb.f0_span = [m](const BarSpanLoads& sp) { return f0_of(*m, sp); };
  tb.recover = [m](const BarSpanLoads& sp, const Vec12& ul, const std::vector<double>& xs,
                   BarResult& out) { recover_with(*m, sp, ul, xs, out); };
  model.tapered[bar_id] = std::move(tb);
}

bool exact_bed_applies(const Model& model, const Bar& bar) {
  if (!bar.soil.any() || bar.kind != BarKind::Frame || model.tapered.count(bar.id)) return false;
  for (bool r : bar.releases)
    if (r) return false;
  const SectionProps& sec = model.section(bar.section).props;
  return sec.Avy == 0 && sec.Avz == 0;
}

int install_exact_beds(Model& model) {
  int converted = 0;
  for (std::size_t e = 0; e < model.bars.size(); ++e) {
    const Bar& b = model.bars[e];
    if (!exact_bed_applies(model, b)) continue;
    const SoilProps soil = b.soil;
    const double length = BarElement(model, b).L();
    model.bars[e].soil = SoilProps{};
    install_bef(model, b.id, {BefZone{length, soil.kfx, soil.kfy, soil.kfz, soil.kmx}});
    ++converted;
  }
  return converted;
}

}  // namespace keta

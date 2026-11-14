#include "keta/nonlinear.h"

#include <Eigen/SparseLU>
#include <limits>
#include <map>
#include <type_traits>

#include "keta/cable.h"
#include "keta/corot_kin.h"
#include "keta/jet.h"
#include "keta/springmass.h"

namespace keta {

// ---------------- rotation helpers ----------------

Mat3 rot_exp(const Vec3& th) {
  double t = th.norm();
  if (t < 1e-12) {
    Mat3 S;
    S << 0, -th.z(), th.y(), th.z(), 0, -th.x(), -th.y(), th.x(), 0;
    return Mat3::Identity() + S + 0.5 * S * S;
  }
  Vec3 a = th / t;
  Mat3 S;
  S << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
  return Mat3::Identity() + std::sin(t) * S + (1 - std::cos(t)) * S * S;
}

Vec3 rot_log(const Mat3& R) {
  double c = std::clamp((R.trace() - 1.0) / 2.0, -1.0, 1.0);
  double t = std::acos(c);
  Vec3 ax(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
  double s = std::sin(t);
  if (t < 1e-8) return 0.5 * ax;  // small angle
  if (t > kPi - 1e-6) {
    // near pi: use the symmetric part
    Mat3 A = 0.5 * (R + Mat3::Identity());
    Vec3 v(std::sqrt(std::max(0.0, A(0, 0))), std::sqrt(std::max(0.0, A(1, 1))),
           std::sqrt(std::max(0.0, A(2, 2))));
    // fix signs from the off-diagonal terms
    if (ax.norm() > 1e-12) {
      for (int i = 0; i < 3; ++i)
        if (ax(i) < 0) v(i) = -v(i);
    }
    return t * v.normalized();
  }
  return (t / (2 * s)) * ax;
}

namespace {
using namespace keta::corot;

// Left Jacobian of the exponential map: dexp so that spin = J(theta) dtheta.
Mat3 rot_left_jacobian(const Vec3& th) {
  double t = th.norm();
  Mat3 S;
  S << 0, -th.z(), th.y(), th.z(), 0, -th.x(), -th.y(), th.x(), 0;
  if (t < 1e-6) return Mat3::Identity() + 0.5 * S + S * S / 6.0;
  double c2 = (1 - std::cos(t)) / (t * t);
  double c3 = (t - std::sin(t)) / (t * t * t);
  return Mat3::Identity() + c2 * S + c3 * S * S;
}


// ---- node-to-node links in large displacements ----
// ELASTIC LINK: the spring point P reached from both nodes by rigid arms
// that turn with them, p_i = x_i + u_i + R_i (P - X_i); the relative
// displacement read in the MEAN frame R_m = R_i exp(psi / 2), psi =
// log(R_i^T R_j) (the relative rotation, the same components in R_i's and
// R_m's frames since exp(-psi/2) psi = psi), and the link axes A:
//   d = A R_m^T (p_j - p_i),   dtheta = A psi,   E = (d, dtheta)' k (d, dtheta) / 2
// — objective (a rigid motion of the pair strains nothing, at any angle), its
// linearization the linear link (link_element). Exact AD force and tangent.
struct CorotLink {
  int ni = 0, nj = 0;
  Vec3 ai = Vec3::Zero(), aj = Vec3::Zero();  // P - X_i, P - X_j
  Mat3 A = Mat3::Identity();                  // rows: link axes in global
  Mat6 k = Mat6::Zero();

  template <class T>
  T energy(const std::array<T, 12>& q) const {
    const M3t<T> Ri = exp_t<T>({q[3], q[4], q[5]}), Rj = exp_t<T>({q[9], q[10], q[11]});
    M3t<T> Rr;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        T v = Ri[static_cast<std::size_t>(a)] * Rj[static_cast<std::size_t>(b)];
        for (int c = 1; c < 3; ++c) v = v + Ri[static_cast<std::size_t>(3 * c + a)] * Rj[static_cast<std::size_t>(3 * c + b)];
        Rr[static_cast<std::size_t>(3 * a + b)] = v;
      }
    const V3t<T> psi = log_t(Rr);
    const M3t<T> Hh = exp_t<T>({psi[0] * 0.5, psi[1] * 0.5, psi[2] * 0.5});
    M3t<T> Rm;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        T v = Ri[static_cast<std::size_t>(3 * a)] * Hh[static_cast<std::size_t>(b)];
        for (int c = 1; c < 3; ++c) v = v + Ri[static_cast<std::size_t>(3 * a + c)] * Hh[static_cast<std::size_t>(3 * c + b)];
        Rm[static_cast<std::size_t>(3 * a + b)] = v;
      }
    V3t<T> dl;
    for (int c = 0; c < 3; ++c) {
      T v = q[static_cast<std::size_t>(6 + c)] - q[static_cast<std::size_t>(c)];
      for (int k2 = 0; k2 < 3; ++k2) {
        if (aj(k2) != 0.0) v = v + Rj[static_cast<std::size_t>(3 * c + k2)] * aj(k2);
        if (ai(k2) != 0.0) v = v - Ri[static_cast<std::size_t>(3 * c + k2)] * ai(k2);
      }
      dl[static_cast<std::size_t>(c)] = v - (aj(c) - ai(c));
    }
    std::array<T, 6> e;
    V3t<T> dgl;  // R_m^T (p_j - p_i)
    for (int a = 0; a < 3; ++a) {
      T dg = Rm[static_cast<std::size_t>(a)] * dl[0];
      for (int c = 1; c < 3; ++c) dg = dg + Rm[static_cast<std::size_t>(3 * c + a)] * dl[static_cast<std::size_t>(c)];
      dgl[static_cast<std::size_t>(a)] = dg;
    }
    for (int r = 0; r < 3; ++r) {
      T v = dgl[0] * A(r, 0), w = psi[0] * A(r, 0);
      for (int c = 1; c < 3; ++c) {
        v = v + dgl[static_cast<std::size_t>(c)] * A(r, c);
        w = w + psi[static_cast<std::size_t>(c)] * A(r, c);
      }
      e[static_cast<std::size_t>(r)] = v;
      e[static_cast<std::size_t>(3 + r)] = w;
    }
    T E = e[0] * 0.0;
    for (int a = 0; a < 6; ++a)
      for (int b = 0; b < 6; ++b)
        if (k(a, b) != 0.0) E = E + e[static_cast<std::size_t>(a)] * e[static_cast<std::size_t>(b)] * (0.5 * k(a, b));
    return E;
  }
  ad::Jet<12> jet(const Eigen::Matrix<double, 12, 1>& q) const {
    std::array<ad::Jet<12>, 12> x;
    for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Jet<12>::seed(q(i), i);
    return energy(x);
  }
};

// RIGID LINK: a slave DOF as an exact function of its master's — a
// translation u_s,c = u_m,c + ((R(theta_m) - I) r)_c, r = X_s - X_m (the
// slave rides on the master's body at any rotation), a rotation theta_s =
// theta_m (the same rotation). The solve runs on the independent DOFs q:
// u = g(q), residual J' r, tangent J' K J + sum_i r_i d2g_i/dq2 (the second
// variation of the constraint: the rigid arm's term at finite rotation).
struct CorotSlave {
  int fs = -1;    // free index of the slave DOF
  int sn = 0;     // slave node index
  int comp = 0;   // 0..5
  int m = 0;      // master node index
  Vec3 r = Vec3::Zero();
};

// A LINEAR slave of the large-displacement reduction: a translation (or a
// fully held rotation) eliminated by a linear constraint over independent
// translations — exact at any displacement, translations being additive
// (skew supports, translational general constraints).
struct CorotLinSlave {
  int fs = -1;                              // free index of the slave DOF
  std::vector<std::pair<int, double>> terms;  // (free index of an independent DOF, coefficient)
  // for the reactions: the slave (node, comp) and every term (supported ones
  // included), its imbalance carried to them (general constraints only)
  int sn = 0, comp = 0;
  std::vector<std::tuple<int, int, double>> all_terms;
};

// An exact elastic catenary between two nodes (keta/cable.h): its forces on
// the nodes and their translational tangent come from eval_cable at the
// current positions; its weight acts in full from the start.
struct CorotCable {
  std::size_t bar_ix = 0;
  int ni = 0, nj = 0;
  Vec3 X1 = Vec3::Zero(), X2 = Vec3::Zero();
  double L0 = 0, EA = 0, w = 0;
};

struct CorotSystem {
  const Model* model;
  Assembled asm_;
  std::vector<CorotElem> elems;
  std::vector<CorotCable> cables;
  std::vector<int> elem_of_bar;  // bar index -> elems index (-1: a cable)
  std::vector<CorotShell> shells;  // facet shells (EICR)
  double force_floor = 0;          // round-off floor of the element forces
  VecX Fext;      // reference nodal FORCES (free DOFs); moments handled apart
  std::vector<std::pair<int, Vec3>> ext_moments;  // (node index, spatial moment)
  std::vector<Vec3> nodal_forces;  // reference nodal FORCES on every node (blocked DOFs included)
  std::vector<CorotBarLoad> bar_loads;            // span loads, state-dependent
  VecX spring;    // per free-DOF spring stiffness
  // node-to-node links (see CorotLink / CorotSlave): the solve runs on the
  // independent DOFs q (n_red of them), u = full(q)
  std::vector<CorotLink> clinks;
  std::vector<CorotSlave> slaves;
  std::vector<CorotLinSlave> lin_slaves;
  std::vector<int> red_of_free, free_of_red;
  int n_red = 0;
  bool slave_loaded = false;  // a nodal force on a slave translation (state-dependent J' F)
  VecX Fref;                  // the reference load on q (constant unless slave_loaded)

  static CorotSystem build(const Model& model, Id case_id, CorotRotationMeasure measure, bool allow_cables = false) {
    model.validate();
    // The corotational residual/tangent carries springs as a DIAGONAL vector;
    // a coupled 6x6 elastic support (ST1 APPUI EL) is a typed reject rather
    // than a silently dropped stiffness.
    require(!has_full_spring(model),
            "corotational: full 6x6 elastic support (ST1 APPUI EL) is not supported — "
            "use per-DOF springs");
    CorotSystem s;
    s.model = &model;
    s.asm_ = assemble_elements(model, true);  // links: handled below (CorotLink, CorotSlave)
    std::vector<double> l0;
    for (const Bar& b : model.bars)
      if (b.kind == BarKind::Cable) {
        require(allow_cables, "corotational: cables are solved by solve_corotational (or solve_cable_static), "
                              "not by the path-following drivers");
        l0 = resolve_cable_l0(model);
        break;
      }
    s.elem_of_bar.assign(model.bars.size(), -1);
    for (std::size_t e = 0; e < model.bars.size(); ++e) {
      const Bar& b = model.bars[e];
      if (b.kind == BarKind::Cable) {
        CorotCable c;
        c.bar_ix = e;
        c.ni = model.node_index(b.ni);
        c.nj = model.node_index(b.nj);
        c.X1 = model.node(b.ni).x;
        c.X2 = model.node(b.nj).x;
        c.L0 = l0[e];
        std::tie(c.EA, c.w) = cable_ea_w(model, b);
        s.cables.push_back(c);
        const double size = std::max({c.X1.cwiseAbs().maxCoeff(), c.X2.cwiseAbs().maxCoeff(), c.L0});
        s.force_floor = std::max(s.force_floor, 64 * std::numeric_limits<double>::epsilon() * size * c.EA / c.L0);
        continue;
      }
      s.elem_of_bar[e] = static_cast<int>(s.elems.size());
      s.elems.push_back(make_corot_elem(model, s.asm_, e, measure));
      // the bar's round-off floor: chord and frame are differences of
      // positions (cancellation ~ eps |X|), the stretch noise eps |X| and the
      // angle noise eps |X| / l0 times the axial / rotational stiffnesses
      const CorotElem& ce = s.elems.back();
      const double size = std::max({ce.X1.cwiseAbs().maxCoeff(), ce.X2.cwiseAbs().maxCoeff(), ce.l0});
      const double krot = ce.Kbar.bottomRightCorner<6, 6>().cwiseAbs().rowwise().sum().maxCoeff();
      s.force_floor = std::max(s.force_floor, 64 * std::numeric_limits<double>::epsilon() * size *
                                                  std::max({std::fabs(ce.Kbar(0, 0)), krot / ce.l0, krot / (ce.l0 * ce.l0)}));
    }
    for (const ShellElem& se : s.asm_.shells) {
      s.shells.push_back(make_corot_shell(se));
      const CorotShell& cs = s.shells.back();
      double size = 0;
      for (const Vec3& x : cs.X) size = std::max(size, x.cwiseAbs().maxCoeff());
      s.force_floor = std::max(s.force_floor, 64 * std::numeric_limits<double>::epsilon() *
                                                  cs.Kl.cwiseAbs().rowwise().sum().maxCoeff() * size);
    }
    const LoadCase& lc = model.load_case(case_id);
    s.Fext = VecX::Zero(s.asm_.dofs.n_free);
    s.nodal_forces.assign(model.nodes.size(), Vec3::Zero());
    auto add_bar_load = [&](std::size_t e, const Load& l) {
      require(s.elem_of_bar[e] >= 0, "corotational: a span load on a cable (its weight is intrinsic)");
      s.bar_loads.push_back(make_corot_bar_load(model, s.asm_, e, static_cast<std::size_t>(s.elem_of_bar[e]), l));
    };
    for (const Load& l : lc.loads) {
      switch (l.kind) {
        case Load::Kind::NodalForce: {
          int ni = model.node_index(l.node);
          s.nodal_forces[static_cast<std::size_t>(ni)] += l.F;
          for (int c = 0; c < 3; ++c) {
            int fi = s.asm_.dofs.free_ix[static_cast<std::size_t>(s.asm_.dofs.full(ni, c))];
            if (fi >= 0) s.Fext(fi) += l.F(c);
          }
          // A spatially-fixed nodal moment is conjugate to the SPIN, not to the
          // rotation vector: its generalized force is J(theta)^t M,
          // state-dependent.
          if (l.M.squaredNorm() > 0) s.ext_moments.emplace_back(ni, l.M);
          break;
        }
        case Load::Kind::BarUdl:
        case Load::Kind::BarLinear:
        case Load::Kind::BarPoint: {
          require(l.frame != LoadFrame::GlobalProjected,
                  "corotational: projected (PRO) bar loads are not supported — use a "
                  "global or local intensity");
          std::size_t e = static_cast<std::size_t>(model.bar_index(l.bar));
          add_bar_load(e, l);
          break;
        }
        case Load::Kind::SelfWeight: {
          // Dead load: rho A(x) g per bar (dead, its exact fixed-end forces)
          // + point-mass forces.
          for (std::size_t e = 0; e < model.bars.size(); ++e) {
            const Bar& b = model.bars[e];
            const Material& mat = model.material(b.material);
            (void)mat;
            // rho A(x) g: the variable weight of a tapered bar included (a
            // cable carries its own weight intrinsically)
            if (s.elem_of_bar[e] < 0) continue;
            s.bar_loads.push_back(
                make_corot_self_weight(model, s.asm_, e, static_cast<std::size_t>(s.elem_of_bar[e]), l.g));
          }
          for (const PointMass& pm : model.point_masses) {
            int ni = model.node_index(pm.node);
            Vec3 F = pm.mass * l.g;
            s.nodal_forces[static_cast<std::size_t>(ni)] += F;
            for (int c = 0; c < 3; ++c) {
              int fi = s.asm_.dofs.free_ix[static_cast<std::size_t>(s.asm_.dofs.full(ni, c))];
              if (fi >= 0) s.Fext(fi) += F(c);
            }
          }
          break;
        }
        default:
          fail("corotational: unsupported load kind in case '" + lc.name +
               "' (nodal forces and bar UDL/linear/point/self-weight only)");
      }
    }
    s.spring = VecX::Zero(s.asm_.dofs.n_free);
    for (std::size_t f = 0; f < s.asm_.spring.size(); ++f) {
      int fi = s.asm_.dofs.free_ix[f];
      if (fi >= 0) s.spring(fi) += s.asm_.spring[f];
    }
    s.build_links();
    return s;
  }

  void build_links() {
    const Model& m = *model;
    const DofMap& d = asm_.dofs;
    for (const ElasticLink& l : m.links) {
      require(!l.has_law(), "large displacements: link " + std::to_string(l.id) +
                                " carries a non-linear law (solve_static_link_laws, linear kinematics)");
      const LinkElem le = link_element(m, l);
      CorotLink cl;
      cl.ni = le.ni;
      cl.nj = le.nj;
      const Vec3 xi = m.node(l.ni).x, xj = m.node(l.nj).x, P = xi + l.s * (xj - xi);
      cl.ai = P - xi;
      cl.aj = P - xj;
      cl.A = le.B.block<3, 3>(0, 6);
      cl.k = le.k;
      const double size = std::max({xi.cwiseAbs().maxCoeff(), xj.cwiseAbs().maxCoeff(), 1.0});
      force_floor = std::max(force_floor, 64 * std::numeric_limits<double>::epsilon() *
                                              le.K.cwiseAbs().rowwise().sum().maxCoeff() * size);
      clinks.push_back(cl);
    }
    std::vector<char> is_slave_node(m.nodes.size(), 0);
    for (const RigidLink& l : m.rigid_links)
      for (Id sid : l.slaves) is_slave_node[static_cast<std::size_t>(m.node_index(sid))] = 1;
    for (const RigidLink& l : m.rigid_links) {
      const int mi = m.node_index(l.master);
      require(!is_slave_node[static_cast<std::size_t>(mi)],
              "large displacements: rigid link " + std::to_string(l.id) +
                  ": its master is a slave of another link (chained links: tie every slave to the root master)");
      const int nrot = l.dofs[3] + l.dofs[4] + l.dofs[5];
      require(nrot == 0 || nrot == 3,
              "large displacements: rigid link " + std::to_string(l.id) +
                  " ties some rotation components only — a component of the rotation vector is not a hinge at "
                  "finite rotation (tie all three rotations, or none)");
      for (Id sid : l.slaves) {
        const int si = m.node_index(sid);
        for (int c = 0; c < 6; ++c) {
          if (!l.dofs[static_cast<std::size_t>(c)]) continue;
          CorotSlave sl;
          sl.fs = d.free_of(si, c);
          require(sl.fs >= 0, "large displacements: a slave DOF of rigid link " + std::to_string(l.id) + " is supported");
          sl.comp = c;
          sl.sn = si;
          sl.m = mi;
          sl.r = m.node(sid).x - m.node(l.master).x;
          slaves.push_back(sl);
        }
      }
    }
    // linear slaves: translational general constraints and skew supports
    std::vector<char> taken(static_cast<std::size_t>(d.n_free), 0);
    for (const CorotSlave& sl : slaves) taken[static_cast<std::size_t>(sl.fs)] = 1;
    auto add_lin = [&](int fs, const std::vector<std::pair<int, double>>& terms, const std::string& who) {
      require(fs >= 0, "large displacements: " + who + ": its eliminated DOF is supported");
      require(!taken[static_cast<std::size_t>(fs)], "large displacements: " + who + ": its eliminated DOF is already constrained");
      taken[static_cast<std::size_t>(fs)] = 1;
      CorotLinSlave ls;
      ls.fs = fs;
      ls.terms = terms;
      lin_slaves.push_back(ls);
    };
    for (const Constraint& c : m.constraints) {
      const std::string who = "constraint " + std::to_string(c.id);
      for (const ConstraintTerm& t : c.terms)
        require(t.comp < 3, "large displacements: " + who +
                                " involves a rotation — a linear relation of rotation-vector components is not "
                                "exact at finite rotation (translations only)");
      const ConstraintTerm& t0 = c.terms.front();
      std::vector<std::pair<int, double>> terms;
      for (std::size_t k = 1; k < c.terms.size(); ++k) {
        const int f = d.free_of(m.node_index(c.terms[k].node), c.terms[k].comp);
        if (f >= 0 && c.terms[k].coef != 0) terms.emplace_back(f, -c.terms[k].coef / t0.coef);
      }
      add_lin(d.free_of(m.node_index(t0.node), t0.comp), terms, who);
      CorotLinSlave& ls = lin_slaves.back();
      ls.sn = m.node_index(t0.node);
      ls.comp = t0.comp;
      for (std::size_t k = 1; k < c.terms.size(); ++k)
        ls.all_terms.emplace_back(m.node_index(c.terms[k].node), c.terms[k].comp, -c.terms[k].coef / t0.coef);
    }
    for (const SkewSupport& sk : m.skew_supports) {
      const std::string who = "skew support " + std::to_string(sk.id);
      const int ni = m.node_index(sk.node);
      const int nrot = sk.fixed[3] + sk.fixed[4] + sk.fixed[5];
      require(nrot == 0 || nrot == 3, "large displacements: " + who +
                                          " holds some rotations only — not a hinge at finite rotation (all or none)");
      for (int k = 0; k < 6; ++k) {
        if (!sk.fixed[static_cast<std::size_t>(k)]) continue;
        if (k >= 3) {  // all rotations held: each global rotation is zero
          add_lin(d.free_of(ni, k), {}, who);
          continue;
        }
        const Vec3 a = sk.axes.row(k).transpose();
        int best = -1;
        for (int c = 0; c < 3; ++c) {
          const int f = d.free_of(ni, c);
          if (f < 0 || taken[static_cast<std::size_t>(f)] || std::fabs(a(c)) < 1e-12) continue;
          if (best < 0 || std::fabs(a(c)) > std::fabs(a(best))) best = c;
        }
        require(best >= 0, "large displacements: " + who + ": no free translation left to eliminate");
        std::vector<std::pair<int, double>> terms;
        for (int c = 0; c < 3; ++c) {
          const int f = d.free_of(ni, c);
          if (c != best && f >= 0 && a(c) != 0) terms.emplace_back(f, -a(c) / a(best));
        }
        add_lin(d.free_of(ni, best), terms, who);
      }
    }
    for (const CorotLinSlave& ls : lin_slaves)
      for (const auto& [f, c] : ls.terms)
        require(!taken[static_cast<std::size_t>(f)],
                "large displacements: a constraint refers to another eliminated DOF (chains are not resolved here)");
    red_of_free.assign(static_cast<std::size_t>(d.n_free), 0);
    for (const CorotSlave& sl : slaves) red_of_free[static_cast<std::size_t>(sl.fs)] = -1;
    for (const CorotLinSlave& ls : lin_slaves) red_of_free[static_cast<std::size_t>(ls.fs)] = -1;
    free_of_red.clear();
    for (int f = 0; f < d.n_free; ++f)
      if (red_of_free[static_cast<std::size_t>(f)] == 0) {
        red_of_free[static_cast<std::size_t>(f)] = static_cast<int>(free_of_red.size());
        free_of_red.push_back(f);
      }
    n_red = static_cast<int>(free_of_red.size());
    for (const CorotSlave& sl : slaves)
      if (sl.comp < 3 && Fext(sl.fs) != 0) slave_loaded = true;
    Fref = linked() ? VecX(jacobian(VecX::Zero(d.n_free)).transpose() * Fext) : Fext;
  }

  bool linked() const { return !slaves.empty() || !lin_slaves.empty(); }
  int n() const { return linked() ? n_red : asm_.dofs.n_free; }
  VecX restrict(const VecX& u) const {
    if (!linked()) return u;
    VecX q(n_red);
    for (int k = 0; k < n_red; ++k) q(k) = u(free_of_red[static_cast<std::size_t>(k)]);
    return q;
  }
  // the master's 6 DOFs at u (blocked ones: 0)
  Vec6 master_dofs(const VecX& u, int mi) const {
    Vec6 v = Vec6::Zero();
    for (int c = 0; c < 6; ++c) {
      const int f = asm_.dofs.free_of(mi, c);
      if (f >= 0) v(c) = u(f);
    }
    return v;
  }
  // a translation slave: (R(theta) - I) r component c as a jet in theta
  static ad::Jet<3> arm_jet(const Vec3& th, const Vec3& r, int c) {
    const M3t<ad::Jet<3>> R = exp_t<ad::Jet<3>>({ad::Jet<3>::seed(th(0), 0), ad::Jet<3>::seed(th(1), 1), ad::Jet<3>::seed(th(2), 2)});
    ad::Jet<3> v = R[static_cast<std::size_t>(3 * c)] * r(0);
    for (int k = 1; k < 3; ++k) v = v + R[static_cast<std::size_t>(3 * c + k)] * r(k);
    return v + (-r(c));
  }
  VecX full(const VecX& q) const {
    if (!linked()) return q;
    VecX u = VecX::Zero(asm_.dofs.n_free);
    for (int k = 0; k < n_red; ++k) u(free_of_red[static_cast<std::size_t>(k)]) = q(k);
    for (const CorotSlave& sl : slaves) {
      const Vec6 mv = master_dofs(u, sl.m);
      u(sl.fs) = sl.comp < 3 ? mv(sl.comp) + arm_jet(mv.tail<3>(), sl.r, sl.comp).v : mv(sl.comp);
    }
    for (const CorotLinSlave& ls : lin_slaves) {
      double v = 0;
      for (const auto& [f, c] : ls.terms) v += c * u(f);
      u(ls.fs) = v;
    }
    return u;
  }
  // J = du/dq (n_free x n_red)
  SpMat jacobian(const VecX& u) const {
    std::vector<Triplet> t;
    for (int k = 0; k < n_red; ++k) t.emplace_back(free_of_red[static_cast<std::size_t>(k)], k, 1.0);
    for (const CorotSlave& sl : slaves) {
      auto put = [&](int c, double v) {
        const int f = asm_.dofs.free_of(sl.m, c);
        if (f >= 0 && v != 0) t.emplace_back(sl.fs, red_of_free[static_cast<std::size_t>(f)], v);
      };
      put(sl.comp, 1.0);
      if (sl.comp < 3) {
        const ad::Jet<3> a = arm_jet(master_dofs(u, sl.m).tail<3>(), sl.r, sl.comp);
        for (int k = 0; k < 3; ++k) put(3 + k, a.g[static_cast<std::size_t>(k)]);
      }
    }
    for (const CorotLinSlave& ls : lin_slaves)
      for (const auto& [f, c] : ls.terms) t.emplace_back(ls.fs, red_of_free[static_cast<std::size_t>(f)], c);
    SpMat J(asm_.dofs.n_free, n_red);
    J.setFromTriplets(t.begin(), t.end());
    return J;
  }
  // sum_i w_i d2g_i/dq2 over the translation slaves (w on the extended space)
  void add_constraint_hessian(const VecX& u, const VecX& w, std::vector<Triplet>& t) const {
    for (const CorotSlave& sl : slaves) {
      if (sl.comp >= 3 || w(sl.fs) == 0) continue;
      const ad::Jet<3> a = arm_jet(master_dofs(u, sl.m).tail<3>(), sl.r, sl.comp);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
          const int fi = asm_.dofs.free_of(sl.m, 3 + i), fj = asm_.dofs.free_of(sl.m, 3 + j);
          if (fi >= 0 && fj >= 0)
            t.emplace_back(red_of_free[static_cast<std::size_t>(fi)], red_of_free[static_cast<std::size_t>(fj)],
                           w(sl.fs) * a.h[static_cast<std::size_t>(3 * i + j)]);
        }
    }
  }

  // ---- the reduced system (q): every driver works on it ----
  VecX internal_force(const VecX& q, double* elem_scale = nullptr) const {
    if (!linked()) return xinternal_force(q, elem_scale);
    const VecX u = full(q);
    return jacobian(u).transpose() * xinternal_force(u, elem_scale);
  }
  VecX external_force(const VecX& q, double lam) const {
    if (!linked()) return xexternal_force(q, lam);
    const VecX u = full(q);
    return jacobian(u).transpose() * xexternal_force(u, lam);
  }
  SpMat tangent(const VecX& q) const {
    if (!linked()) return xtangent(q);
    const VecX u = full(q);
    const SpMat J = jacobian(u);
    SpMat K = (J.transpose() * xtangent(u) * J).pruned();
    std::vector<Triplet> t;
    add_constraint_hessian(u, xinternal_force(u), t);
    SpMat H(n_red, n_red);
    H.setFromTriplets(t.begin(), t.end());
    return K + H;
  }
  SpMat load_stiffness(const VecX& q) const {
    if (!linked()) return xload_stiffness(q);
    const VecX u = full(q);
    const SpMat J = jacobian(u);
    SpMat K = (J.transpose() * xload_stiffness(u) * J).pruned();
    std::vector<Triplet> t;
    add_constraint_hessian(u, xexternal_force(u, 1.0), t);
    SpMat H(n_red, n_red);
    H.setFromTriplets(t.begin(), t.end());
    return K + H;
  }
  std::vector<Reaction> reactions(const VecX& q, double lam) const { return xreactions(full(q), lam); }
  void recanonicalize(VecX& q) const {
    if (!linked()) return xrecanonicalize(q);
    VecX u = full(q);
    xrecanonicalize(u);
    q = restrict(u);
  }

  Eigen::Matrix<double, 18, 1> gather(const VecX& u, const CorotShell& cs) const {
    Eigen::Matrix<double, 18, 1> q = Eigen::Matrix<double, 18, 1>::Zero();
    for (int k = 0; k < 3; ++k)
      for (int c = 0; c < 6; ++c) {
        const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(cs.nodes[static_cast<std::size_t>(k)], c))];
        if (f >= 0) q(6 * k + c) = u(f);
      }
    return q;
  }
  int shell_dof(const CorotShell& cs, int a) const {
    return asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(cs.nodes[static_cast<std::size_t>(a / 6)], a % 6))];
  }

  Eigen::Matrix<double, 12, 1> gather(const VecX& u, const CorotElem& ce) const {
    Eigen::Matrix<double, 12, 1> q = Eigen::Matrix<double, 12, 1>::Zero();
    for (int c = 0; c < 6; ++c) {
      int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.ni, c))];
      int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.nj, c))];
      if (f1 >= 0) q(c) = u(f1);
      if (f2 >= 0) q(6 + c) = u(f2);
    }
    return q;
  }

  template <class F>
  void add_element(VecX& out, const CorotElem& ce, F&& value) const {
    for (int c = 0; c < 6; ++c) {
      int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.ni, c))];
      int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.nj, c))];
      if (f1 >= 0) out(f1) += value(c);
      if (f2 >= 0) out(f2) += value(6 + c);
    }
  }

  bool has_load_stiffness() const { return !ext_moments.empty() || !bar_loads.empty() || slave_loaded; }

  // d(external force)/du per unit load factor (non-symmetric): the
  // rotation-vector chart makes a spatially fixed moment state-dependent
  // (J(theta)^T M) and a span load follows the corotated frame. Exact, from
  // one first-order AD pass per moment / span load.
  SpMat xload_stiffness(const VecX& u) const {
    std::vector<Triplet> trips;
    for (const auto& [ni, M] : ext_moments) {
      int fx[3];
      V3t<ad::Grad<3>> th;
      for (int c = 0; c < 3; ++c) {
        fx[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, 3 + c))];
        th[static_cast<std::size_t>(c)] = ad::Grad<3>::seed(fx[c] >= 0 ? u(fx[c]) : 0.0, c);
      }
      const V3t<double> m{M(0), M(1), M(2)};
      const V3t<ad::Grad<3>> g = jt_times(th, m);
      for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
          if (fx[r] >= 0 && fx[c] >= 0) trips.emplace_back(fx[r], fx[c], g[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)]);
    }
    for (const CorotBarLoad& bl : bar_loads) {
      const CorotElem& ce = elems[bl.ei];
      const Eigen::Matrix<double, 12, 1> q = gather(u, ce);
      std::array<ad::Grad<12>, 12> x;
      for (int i = 0; i < 12; ++i) x[static_cast<std::size_t>(i)] = ad::Grad<12>::seed(q(i), i);
      const std::array<ad::Grad<12>, 12> f = corot::bar_load_vector(ce, bl, x);
      int g[12];
      for (int c = 0; c < 6; ++c) {
        g[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.ni, c))];
        g[6 + c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.nj, c))];
      }
      for (int r = 0; r < 12; ++r)
        for (int c = 0; c < 12; ++c)
          if (g[r] >= 0 && g[c] >= 0) trips.emplace_back(g[r], g[c], f[static_cast<std::size_t>(r)].g[static_cast<std::size_t>(c)]);
    }
    SpMat K(asm_.dofs.n_free, asm_.dofs.n_free);
    K.setFromTriplets(trips.begin(), trips.end());
    return K;
  }

  VecX xexternal_force(const VecX& u, double lam) const {
    VecX F = lam * Fext;
    for (const auto& [ni, M] : ext_moments) {
      Vec3 th = Vec3::Zero();
      int fx[3];
      for (int c = 0; c < 3; ++c) {
        fx[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, 3 + c))];
        if (fx[c] >= 0) th(c) = u(fx[c]);
      }
      Vec3 gm = rot_left_jacobian(th).transpose() * (lam * M);
      for (int c = 0; c < 3; ++c)
        if (fx[c] >= 0) F(fx[c]) += gm(c);
    }
    // Span loads: lower on the CURRENT corotated frame, push equivalent nodal
    // forces back to global; moments go through J(theta)^T like ext_moments.
    for (const CorotBarLoad& bl : bar_loads) {
      const CorotElem& ce = elems[bl.ei];
      const std::array<double, 12> f = corot::bar_load_vector(ce, bl, CorotElem::arr(gather(u, ce)));
      add_element(F, ce, [&](int k) { return lam * f[static_cast<std::size_t>(k)]; });
    }
    return F;
  }

  // `elem_scale` (optional) receives max |f_e| over every element end force,
  // supported DOFs included: the stress level of the structure, which stays
  // finite when the assembled free-DOF force vanishes (an unloaded deformed
  // equilibrium, lambda crossing 0 on a snap-through path).
  VecX xinternal_force(const VecX& u, double* elem_scale = nullptr) const {
    VecX f = VecX::Zero(u.size());
    if (elem_scale) *elem_scale = 0;
    for (const CorotElem& ce : elems) {
      auto q = gather(u, ce);
      auto g = ce.grad(q);
      if (elem_scale) *elem_scale = std::max(*elem_scale, g.template lpNorm<Eigen::Infinity>());
      for (int c = 0; c < 6; ++c) {
        int f1 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.ni, c))];
        int f2 = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.nj, c))];
        if (f1 >= 0) f(f1) += g(c);
        if (f2 >= 0) f(f2) += g(6 + c);
      }
    }
    for (const CorotShell& cs : shells) {
      const auto g = cs.grad(gather(u, cs));
      if (elem_scale) *elem_scale = std::max(*elem_scale, g.template lpNorm<Eigen::Infinity>());
      for (int a = 0; a < 18; ++a) {
        const int fa = shell_dof(cs, a);
        if (fa >= 0) f(fa) += g(a);
      }
    }
    for (const CorotCable& c : cables) {
      const CableEval ev = cable_eval(u, c);
      if (elem_scale) *elem_scale = std::max(*elem_scale, ev.T_max);
      for (int k = 0; k < 3; ++k) {
        const int f1 = asm_.dofs.free_of(c.ni, k), f2 = asm_.dofs.free_of(c.nj, k);
        if (f1 >= 0) f(f1) -= ev.f1(k);
        if (f2 >= 0) f(f2) -= ev.f2(k);
      }
    }
    for (const CorotLink& cl : clinks) {
      const ad::Jet<12> e = cl.jet(gather2(u, cl.ni, cl.nj));
      if (elem_scale)
        for (double g : e.g) *elem_scale = std::max(*elem_scale, std::fabs(g));
      for (int a = 0; a < 12; ++a) {
        const int fa = asm_.dofs.free_of(a < 6 ? cl.ni : cl.nj, a % 6);
        if (fa >= 0) f(fa) += e.g[static_cast<std::size_t>(a)];
      }
    }
    f += spring.cwiseProduct(u);
    return f;
  }
  CableEval cable_eval(const VecX& u, const CorotCable& c) const {
    Vec3 p1 = c.X1, p2 = c.X2;
    for (int k = 0; k < 3; ++k) {
      const int f1 = asm_.dofs.free_of(c.ni, k), f2 = asm_.dofs.free_of(c.nj, k);
      if (f1 >= 0) p1(k) += u(f1);
      if (f2 >= 0) p2(k) += u(f2);
    }
    return eval_cable(p1, p2, c.L0, c.EA, c.w);
  }
  Eigen::Matrix<double, 12, 1> gather2(const VecX& u, int ni, int nj) const {
    Eigen::Matrix<double, 12, 1> q = Eigen::Matrix<double, 12, 1>::Zero();
    for (int c = 0; c < 6; ++c) {
      const int f1 = asm_.dofs.free_of(ni, c), f2 = asm_.dofs.free_of(nj, c);
      if (f1 >= 0) q(c) = u(f1);
      if (f2 >= 0) q(6 + c) = u(f2);
    }
    return q;
  }

  SpMat xtangent(const VecX& u) const {
    std::vector<Triplet> trips;
    for (const CorotElem& ce : elems) {
      auto q = gather(u, ce);
      auto K = ce.hess(q);
      int g[12];
      for (int c = 0; c < 6; ++c) {
        g[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.ni, c))];
        g[6 + c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ce.nj, c))];
      }
      for (int a = 0; a < 12; ++a)
        for (int bb = 0; bb < 12; ++bb)
          if (g[a] >= 0 && g[bb] >= 0) trips.emplace_back(g[a], g[bb], 0.5 * (K(a, bb) + K(bb, a)));
    }
    for (const CorotShell& cs : shells) {
      const auto K = cs.hess(gather(u, cs));
      int g[18];
      for (int a = 0; a < 18; ++a) g[a] = shell_dof(cs, a);
      for (int a = 0; a < 18; ++a)
        for (int bb = 0; bb < 18; ++bb)
          if (g[a] >= 0 && g[bb] >= 0) trips.emplace_back(g[a], g[bb], 0.5 * (K(a, bb) + K(bb, a)));
    }
    for (const CorotCable& c : cables) {
      const Mat3 k3 = cable_eval(u, c).k3;
      int g[6];
      for (int k = 0; k < 3; ++k) {
        g[k] = asm_.dofs.free_of(c.ni, k);
        g[3 + k] = asm_.dofs.free_of(c.nj, k);
      }
      for (int a = 0; a < 3; ++a)
        for (int bb = 0; bb < 3; ++bb) {
          const double v = k3(a, bb);
          if (g[a] >= 0 && g[bb] >= 0) trips.emplace_back(g[a], g[bb], v);
          if (g[3 + a] >= 0 && g[3 + bb] >= 0) trips.emplace_back(g[3 + a], g[3 + bb], v);
          if (g[a] >= 0 && g[3 + bb] >= 0) trips.emplace_back(g[a], g[3 + bb], -v);
          if (g[3 + a] >= 0 && g[bb] >= 0) trips.emplace_back(g[3 + a], g[bb], -v);
        }
    }
    for (const CorotLink& cl : clinks) {
      const ad::Jet<12> e = cl.jet(gather2(u, cl.ni, cl.nj));
      int g[12];
      for (int a = 0; a < 12; ++a) g[a] = asm_.dofs.free_of(a < 6 ? cl.ni : cl.nj, a % 6);
      for (int a = 0; a < 12; ++a)
        for (int bb = 0; bb < 12; ++bb)
          if (g[a] >= 0 && g[bb] >= 0) trips.emplace_back(g[a], g[bb], e.h[static_cast<std::size_t>(12 * a + bb)]);
    }
    for (Eigen::Index i = 0; i < spring.size(); ++i)
      if (spring(i) != 0) trips.emplace_back(static_cast<int>(i), static_cast<int>(i), spring(i));
    SpMat Kt(asm_.dofs.n_free, asm_.dofs.n_free);
    Kt.setFromTriplets(trips.begin(), trips.end());
    return Kt;
  }

  // Support reactions at the state u, load level lam: on every node the
  // element forces minus the external actions, over ALL its DOFs (blocked
  // ones included), kept where a support blocks a DOF or carries a spring —
  // the force the support exerts on the structure, as in linear statics.
  std::vector<Reaction> xreactions(const VecX& u, double lam) const {
    const Model& m = *model;
    std::vector<Vec6> r = node_imbalance(u, lam);
    // a slave's imbalance (its constraint force) travels to its master
    // through the constraint's derivatives (blocked master DOFs included)
    for (const CorotSlave& sl : slaves) {
      const double rho = r[static_cast<std::size_t>(sl.sn)](sl.comp);
      Vec6& rm = r[static_cast<std::size_t>(sl.m)];
      rm(sl.comp) += rho;
      if (sl.comp < 3) {
        const ad::Jet<3> a = arm_jet(master_dofs(u, sl.m).tail<3>(), sl.r, sl.comp);
        for (int k = 0; k < 3; ++k) rm(3 + k) += rho * a.g[static_cast<std::size_t>(k)];
      }
    }
    for (const CorotLinSlave& ls : lin_slaves) {
      if (ls.all_terms.empty()) continue;
      const double rho = r[static_cast<std::size_t>(ls.sn)](ls.comp);
      for (const auto& [ni, c, coef] : ls.all_terms) r[static_cast<std::size_t>(ni)](c) += coef * rho;
      r[static_cast<std::size_t>(ls.sn)](ls.comp) -= rho;
    }
    std::vector<Reaction> out;
    for (const Support& sp : m.supports) {
      if (!sp.any()) continue;
      Reaction re;
      re.node = sp.node;
      const Vec6& v = r[static_cast<std::size_t>(m.node_index(sp.node))];
      for (int c = 0; c < 6; ++c)
        if (sp.fixed[static_cast<std::size_t>(c)] || sp.spring[static_cast<std::size_t>(c)] != 0) re.R(c) = v(c);
      out.push_back(re);
    }
    return out;
  }
  // per node: its elements' forces minus the external actions (all DOFs)
  std::vector<Vec6> node_imbalance(const VecX& u, double lam) const {
    const Model& m = *model;
    std::vector<Vec6> r(m.nodes.size(), Vec6::Zero());
    auto scatter = [&](const CorotElem& ce, const auto& v) {
      for (int c = 0; c < 6; ++c) {
        r[static_cast<std::size_t>(ce.ni)](c) += v(c);
        r[static_cast<std::size_t>(ce.nj)](c) += v(6 + c);
      }
    };
    for (const CorotElem& ce : elems) scatter(ce, ce.grad(gather(u, ce)));
    for (const CorotShell& cs : shells) {
      const auto g = cs.grad(gather(u, cs));
      for (int a = 0; a < 18; ++a) r[static_cast<std::size_t>(cs.nodes[static_cast<std::size_t>(a / 6)])](a % 6) += g(a);
    }
    for (const CorotBarLoad& bl : bar_loads) {
      const CorotElem& ce = elems[bl.ei];
      const std::array<double, 12> f = corot::bar_load_vector(ce, bl, CorotElem::arr(gather(u, ce)));
      Eigen::Matrix<double, 12, 1> v;
      for (int k = 0; k < 12; ++k) v(k) = -lam * f[static_cast<std::size_t>(k)];
      scatter(ce, v);
    }
    for (std::size_t ni = 0; ni < m.nodes.size(); ++ni) r[ni].head<3>() -= lam * nodal_forces[ni];
    for (const auto& [ni, M] : ext_moments) {
      Vec3 th = Vec3::Zero();
      for (int c = 0; c < 3; ++c) {
        const int f = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, 3 + c))];
        if (f >= 0) th(c) = u(f);
      }
      r[static_cast<std::size_t>(ni)].tail<3>() -= rot_left_jacobian(th).transpose() * (lam * M);
    }
    for (const CorotLink& cl : clinks) {
      const ad::Jet<12> e = cl.jet(gather2(u, cl.ni, cl.nj));
      for (int a = 0; a < 12; ++a) r[static_cast<std::size_t>(a < 6 ? cl.ni : cl.nj)](a % 6) += e.g[static_cast<std::size_t>(a)];
    }
    for (const CorotCable& c : cables) {
      const CableEval ev = cable_eval(u, c);
      r[static_cast<std::size_t>(c.ni)].head<3>() -= ev.f1;
      r[static_cast<std::size_t>(c.nj)].head<3>() -= ev.f2;
    }
    return r;
  }

  void xrecanonicalize(VecX& u) const {
    for (int ni = 0; ni < asm_.dofs.n_nodes; ++ni) {
      Vec3 th = Vec3::Zero();
      int fx[3];
      for (int c = 0; c < 3; ++c) {
        fx[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, 3 + c))];
        if (fx[c] >= 0) th(c) = u(fx[c]);
      }
      double t = th.norm();
      if (t > 0.9 * kPi) {
        Vec3 th2 = (1.0 - 2 * kPi / t) * th;
        for (int c = 0; c < 3; ++c)
          if (fx[c] >= 0) u(fx[c]) = th2(c);
      }
    }
  }
};

// One Newton solve at load level lam*Fext from a warm start. Returns iterations.
int newton_solve(const CorotSystem& sys, double lam, VecX& u, const CorotOptions& opt) {
  int iters = 0;
  double last_step = std::numeric_limits<double>::infinity();
  for (int it = 0;; ++it) {
    hooks_check(opt.hooks);
    require(it < opt.max_iterations, "corotational: Newton did not converge");
    VecX F = sys.external_force(u, lam);
    VecX fint = sys.internal_force(u);
    VecX r = F - fint;
    double scale = std::max({F.lpNorm<Eigen::Infinity>(), fint.lpNorm<Eigen::Infinity>(), 1e-30});
    // the elements' round-off floor: the deformational vector is a
    // difference of positions (cancellation ~ eps |X|), so a residual below
    // 64 eps max ||K|| |X| is noise, not imbalance — accepted only once the
    // last Newton correction is negligible (a mechanism or an unloaded start
    // below the floor is never taken for an equilibrium)
    const double rn = r.lpNorm<Eigen::Infinity>();
    if (rn <= opt.tol_rel * scale || (rn <= sys.force_floor && last_step <= 1e-8 * u.lpNorm<Eigen::Infinity>())) break;
    ++iters;  // Newton steps (tangent solves)
    VecX du;
    if (sys.has_load_stiffness()) {
      // Consistent tangent K_int - lam K_load: non-symmetric, sparse LU.
      SpMat Kt = sys.tangent(u) - lam * sys.load_stiffness(u);
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(Kt);
      require(lu.info() == Eigen::Success, "corotational: singular tangent");
      du = lu.solve(r);
      require(du.allFinite(), "corotational: singular tangent");
    } else {
      SpMat Kt = sys.tangent(u);
      FactorCache fac;
      fac.factorize(Kt);
      du = fac.solve(r);
    }
    // Line search (Crisfield, Non-linear FE Analysis of Solids and
    // Structures vol. 1, sec. 9.3): accept a trial step once the directional
    // derivative of the potential along du, du^T r, has dropped below 0.8 of
    // its start value — the energy measure the exact Newton step minimizes —
    // or once ||r||_inf has dropped. The inf-norm alone rejected the full
    // Newton step on slender members (its first iterate carries a large but
    // harmless axial residual), cutting the step to crumbs: a tip-loaded
    // cantilever at PL^2/EI = 2 took ~700 iterations over 20 substeps; it now
    // takes 8 in one.
    if (rn <= sys.force_floor) {
      // residual within the round-off: no descent to measure, the full
      // (noise-sized) correction is taken and its size decides the stop
      last_step = du.lpNorm<Eigen::Infinity>();
      u += du;
      sys.recanonicalize(u);
      continue;
    }
    double r0 = r.lpNorm<Eigen::Infinity>();
    const double g0 = std::fabs(du.dot(r));
    double s = 1.0;
    bool accepted = false;
    for (int back = 0; back <= opt.max_backtracks; ++back) {
      VecX ut = u + s * du;
      try {
        VecX rt = sys.external_force(ut, lam) - sys.internal_force(ut);
        if (rt.allFinite() && (std::fabs(du.dot(rt)) <= 0.8 * g0 || rt.lpNorm<Eigen::Infinity>() < r0)) {
          last_step = s * du.lpNorm<Eigen::Infinity>();
          u = ut;
          accepted = true;
          break;
        }
      } catch (const Error&) {
      }
      s *= 0.5;
    }
    require(accepted, "corotational: line search failed");
    sys.recanonicalize(u);
  }
  return iters;
}

}  // namespace

CorotResults solve_corotational(const Model& model, Id case_id, const CorotOptions& opt) {
  require(model.continua.empty(), "solve_corotational: continuum elements are not supported (linear statics, modal and dynamics take them)");
  CorotSystem sys = CorotSystem::build(model, case_id, opt.measure, true);
  VecX u = VecX::Zero(sys.n());
  CorotResults out;
  // Secant predictor between equal load substeps: start substep k + 1 from
  // u_k + (u_k - u_(k-1)) — the previous increment extrapolated — instead of
  // u_k, whose linearized first step overshoots the axial strain on slender
  // members; fall back to u_k if Newton fails from the prediction.
  VecX u_prev = u;
  for (int sub = 1; sub <= opt.substeps; ++sub) {
    const double lam = static_cast<double>(sub) / opt.substeps;
    if (opt.hooks)
      opt.hooks->report(static_cast<double>(sub - 1) / opt.substeps,
                        "substep " + std::to_string(sub) + "/" + std::to_string(opt.substeps));
    VecX trial = sub > 1 ? VecX(2 * u - u_prev) : u;
    u_prev = u;
    try {
      out.iterations += newton_solve(sys, lam, trial, opt);
      u = trial;
    } catch (const Cancelled&) {
      throw;
    } catch (const Error&) {
      if (sub == 1) throw;
      out.iterations += newton_solve(sys, lam, u, opt);
    }
  }

  hooks_report(opt.hooks, 1.0, "done");
  out.reactions = sys.reactions(u, 1.0);
  u = sys.full(u);  // the free DOFs (slaves included)
  if (!model.skew_supports.empty()) {
    const std::vector<Vec6> imb = sys.node_imbalance(u, 1.0);
    for (const SkewSupport& sk : model.skew_supports) {
      CaseResult::SkewReaction sr;
      sr.node = sk.node;
      const Vec6& v = imb[static_cast<std::size_t>(model.node_index(sk.node))];
      for (int k = 0; k < 6; ++k)
        if (sk.fixed[static_cast<std::size_t>(k)]) {
          const Vec3 a = sk.axes.row(k % 3).transpose();
          sr.local(k) = a.dot(v.segment<3>(k < 3 ? 0 : 3));
          sr.global.segment<3>(k < 3 ? 0 : 3) += sr.local(k) * a;
        }
      out.skew_reactions.push_back(sr);
    }
  }
  out.u_free = u;
  out.node_disp.assign(static_cast<std::size_t>(sys.asm_.dofs.n_nodes), Vec6::Zero());
  for (int ni = 0; ni < sys.asm_.dofs.n_nodes; ++ni)
    for (int c = 0; c < 6; ++c) {
      int fi = sys.asm_.dofs.free_ix[static_cast<std::size_t>(sys.asm_.dofs.full(ni, c))];
      if (fi >= 0) out.node_disp[static_cast<std::size_t>(ni)](c) = u(fi);
    }
  // Member end forces from the deformational state.
  out.bar_forces.assign(model.bars.size(), Vec12::Zero());
  for (const auto& ce : sys.elems) {
    auto q = sys.gather(u, ce);
    auto d = ce.deform(q);
    Eigen::Matrix<double, 7, 1> fl = ce.Kbar * d;
    Vec12& f = out.bar_forces[ce.bar_ix];
    f(6) = fl(0);
    f(0) = -fl(0);
    f.segment<3>(3) = fl.segment<3>(1);
    f.segment<3>(9) = fl.segment<3>(4);
    // shears from end-moment equilibrium about the corotated chord
    Vec3 x21 = ce.chord(q);
    double ln = x21.norm();
    f(1) = -(fl(3) + fl(6)) / ln;   // Vy from mz moments
    f(7) = -f(1);
    f(2) = (fl(2) + fl(5)) / ln;    // Vz from my moments
    f(8) = -f(2);
  }
  // cables: the end tensions of the deformed catenaries
  out.cable_T.assign(model.bars.size(), {0.0, 0.0});
  out.cable_H.assign(model.bars.size(), 0.0);
  for (const CorotCable& c : sys.cables) {
    const CableEval ev = sys.cable_eval(u, c);
    out.cable_T[c.bar_ix] = {ev.f1.norm(), ev.f2.norm()};
    out.cable_H[c.bar_ix] = ev.H;
  }
  return out;
}

// ---------------- second order (P-Delta) ----------------

SecondOrderResults solve_second_order(const Model& model, Id case_id,
                                      const SecondOrderOptions& opt) {
  require(model.continua.empty(), "solve_second_order: continuum elements are not supported (linear statics, modal and dynamics take them)");
  hooks_check(opt.hooks);
  model.validate();
  Assembled asm_ = assemble_elements(model);
  SpMat K = assemble_stiffness(model, asm_);
  CaseLoads cl = lower_case(model, asm_, model.load_case(case_id));
  require(cl.prescribed.empty(), "second order: settlements are not supported");
  StaticOptions sopt;
  sopt.stations_per_bar = opt.stations_per_bar;

  FactorCache fac;
  fac.factorize(K);
  fac.require_spd();
  VecX u = fac.solve(cl.rhs);
  std::vector<double> N(model.bars.size(), 0.0);
  std::vector<Vec3> Sn(model.shells.size(), Vec3::Zero());

  // exact_stability: the span loads' fixed-end forces under the current N
  // (beam_column_f0) replace the linear ones in the load vector.
  std::vector<Vec12> f0x;
  auto exact_rhs = [&](const std::vector<double>& Nf) {
    f0x.assign(asm_.elements.size(), Vec12::Zero());
    VecX rhs = cl.rhs;
    for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
      const Bar& b = model.bars[e];
      if (b.kind == BarKind::Cable) continue;
      const BarElement& el = asm_.elements[e];
      f0x[e] = beam_column_f0(model, el, cl.span[e], cl.f0[e], Nf[e]);
      const Vec12 d = el.T().transpose() *
                      (el.condense_f0_with(el.k_exact_local(Nf[e]), f0x[e]) - el.condense_f0(cl.f0[e]));
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      for (int c = 0; c < 6; ++c) {
        int gi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
        int gj = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nj, c))];
        if (gi >= 0) rhs(gi) += d(c);
        if (gj >= 0) rhs(gj) += d(6 + c);
      }
    }
    return rhs;
  };
  VecX rhs_it = cl.rhs;
  bool linear_exit = false;

  int iterations = 0;
  double du_prev = -1.0;
  StaticOptions bar_ends;
  bar_ends.stations_per_bar = 2;
  for (int it = 0;; ++it) {
    hooks_check(opt.hooks);
    if (opt.hooks)
      opt.hooks->report(static_cast<double>(it) / std::max(opt.max_iterations, 1),
                        "P-Delta iteration " + std::to_string(it + 1));
    require(it < opt.max_iterations, "second order: fixed point did not converge");
    ++iterations;
    // axial + shell membrane fields from the current iterate (linear recovery)
    CaseResult tmp = recover_case(model, asm_, cl, u, bar_ends, 0, "");
    std::vector<double> N1 = axial_field(tmp);
    std::vector<Vec3> Sn1 = shell_membrane_field(tmp);
    bool zero = true;
    for (double v : N1) zero &= (v == 0.0);
    for (const Vec3& v : Sn1) zero &= (v.squaredNorm() == 0.0);
    if (zero) {  // bit-identical to the linear solve
      linear_exit = true;
      break;
    }

    SpMat Kt;
    SpMat Kg;
    if (opt.exact_stability) {
      // exact member tangents replace K + Kg wholesale
      std::vector<Triplet> trips;
      for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
        const Bar& b = model.bars[e];
        if (b.kind == BarKind::Cable) continue;
        Mat12 ke = asm_.elements[e].ke_exact(N1[e]);
        int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
        int g[12];
        for (int c = 0; c < 6; ++c) {
          g[c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
          g[6 + c] = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(nj, c))];
        }
        for (int a = 0; a < 12; ++a)
          for (int bb = 0; bb < 12; ++bb)
            if (g[a] >= 0 && g[bb] >= 0 && ke(a, bb) != 0) trips.emplace_back(g[a], g[bb], ke(a, bb));
      }
      for (std::size_t f = 0; f < asm_.spring.size(); ++f)
        if (asm_.spring[f] != 0) {
          int fi = asm_.dofs.free_ix[f];
          if (fi >= 0) trips.emplace_back(fi, fi, asm_.spring[f]);
        }
      scatter_full_springs(asm_, trips);  // ST1 APPUI EL
      // Facet shells: ke + consistent membrane Kg (no exact stability-function
      // analogue exists for a plate — documented in shell.h).
      for (const ShellElem& se : asm_.shells) {
        Mat18 kes = se.el.ke();
        const Vec3& s_n = Sn1[se.index];
        if (s_n.squaredNorm() != 0) kes += se.el.kg(s_n(0), s_n(1), s_n(2));
        int gd[18];
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 6; ++c)
            gd[6 * k + c] = asm_.dofs.free_ix[static_cast<std::size_t>(
                asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c))];
        for (int a = 0; a < 18; ++a)
          for (int bb = 0; bb < 18; ++bb)
            if (gd[a] >= 0 && gd[bb] >= 0 && kes(a, bb) != 0)
              trips.emplace_back(gd[a], gd[bb], kes(a, bb));
      }
      Kt = SpMat(asm_.dofs.n_free, asm_.dofs.n_free);
      Kt.setFromTriplets(trips.begin(), trips.end());
      Kg = (Kt - K).eval();  // for the residual/reaction correction terms
      rhs_it = exact_rhs(N1);
    } else {
      Kg = assemble_geometric(model, asm_, N1, &Sn1);
      Kt = (K + Kg).eval();
    }
    FactorCache ft;
    ft.factorize(Kt);
    require(ft.inertia().n_neg == 0,
            "second order: load exceeds the critical load of the current state");
    VecX u1 = ft.solve(rhs_it);

    double du = (u1 - u).lpNorm<Eigen::Infinity>();
    double un = u1.lpNorm<Eigen::Infinity>();
    VecX res = rhs_it - K * u1 - Kg * u1;
    double scale = std::max({rhs_it.lpNorm<Eigen::Infinity>(), (K * u1).lpNorm<Eigen::Infinity>(),
                             (Kg * u1).lpNorm<Eigen::Infinity>(), 1e-30});
    u = u1;
    N = N1;
    Sn = Sn1;
    const bool res_ok = res.lpNorm<Eigen::Infinity>() <= opt.tol_rel_res * scale;
    // Round-off floor (exact_stability): the load vector f0(N) is re-evaluated
    // from N every iteration, so once the contraction reaches the round-off of
    // f0 (amplified by the conditioning of Kt) the change stops decreasing
    // instead of vanishing. A NON-decreasing change already below 1e-9 |u| is
    // that floor — a contracting fixed point decreases strictly above it.
    const bool floor = opt.exact_stability && du_prev >= 0 && du >= du_prev &&
                       du <= 1e-9 * std::max(un, 1e-300);
    du_prev = du;
    if ((du <= opt.tol_rel_du * std::max(un, 1e-300) || floor) && res_ok)
      break;
  }

  SecondOrderResults out;
  out.result = recover_case(model, asm_, cl, u, sopt, case_id, model.load_case(case_id).name);
  out.result.non_superposable = true;
  out.iterations = iterations;
  if (opt.exact_stability && !linear_exit) {
    // Exact beam-column members: end forces k_exact(N) u - f0(N), stations
    // M = M1 + N w on the exact deflected shape; reactions corrected by
    // (ke_exact - ke) u minus the change of the equivalent span loads.
    VecX ufull = VecX::Zero(kDofPerNode * asm_.dofs.n_nodes);
    for (int f = 0; f < ufull.size(); ++f) {
      int fi = asm_.dofs.free_ix[static_cast<std::size_t>(f)];
      if (fi >= 0) ufull(f) = u(fi);
    }
    VecX fg = VecX::Zero(ufull.size());
    for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
      const Bar& b = model.bars[e];
      if (b.kind == BarKind::Cable) continue;
      const BarElement& el = asm_.elements[e];
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        ue(c) = ufull(asm_.dofs.full(ni, c));
        ue(6 + c) = ufull(asm_.dofs.full(nj, c));
      }
      recover_beam_column(model, el, cl.span[e], f0x[e], N[e], ue, sopt, out.result.bars[e]);
      const Vec12 fe = (el.ke_exact(N[e]) - bar_ke_assembled(model, el)) * ue -
                       el.T().transpose() * (el.condense_f0_with(el.k_exact_local(N[e]), f0x[e]) -
                                             el.condense_f0(cl.f0[e]));
      for (int c = 0; c < 6; ++c) {
        fg(asm_.dofs.full(ni, c)) += fe(c);
        fg(asm_.dofs.full(nj, c)) += fe(6 + c);
      }
    }
    for (const ShellElem& se : asm_.shells) {
      const Vec3& s_n = Sn[se.index];
      if (s_n.squaredNorm() == 0.0) continue;
      Vec18 ue;
      for (int k = 0; k < 3; ++k)
        for (int c = 0; c < 6; ++c)
          ue(6 * k + c) = ufull(asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c));
      Vec18 fe = se.el.kg(s_n(0), s_n(1), s_n(2)) * ue;
      for (int k = 0; k < 3; ++k)
        for (int c = 0; c < 6; ++c)
          fg(asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c)) += fe(6 * k + c);
    }
    for (Reaction& re : out.result.reactions) {
      int ni = model.node_index(re.node);
      for (int c = 0; c < 6; ++c) re.R(c) += fg(asm_.dofs.full(ni, c));
    }
    hooks_report(opt.hooks, 1.0, "done");
    return out;
  }
  // Geometric reaction correction: R += sum_e kg(N_e) u_e at every node.
  {
    VecX ufull = VecX::Zero(kDofPerNode * asm_.dofs.n_nodes);
    for (int f = 0; f < ufull.size(); ++f) {
      int fi = asm_.dofs.free_ix[static_cast<std::size_t>(f)];
      if (fi >= 0) ufull(f) = u(fi);
    }
    VecX fg = VecX::Zero(ufull.size());
    for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
      if (model.bars[e].kind == BarKind::Cable || N[e] == 0.0) continue;
      const Bar& b = model.bars[e];
      int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
      Vec12 ue;
      for (int c = 0; c < 6; ++c) {
        ue(c) = ufull(asm_.dofs.full(ni, c));
        ue(6 + c) = ufull(asm_.dofs.full(nj, c));
      }
      Vec12 fe = asm_.elements[e].kg(N[e]) * ue;
      for (int c = 0; c < 6; ++c) {
        fg(asm_.dofs.full(ni, c)) += fe(c);
        fg(asm_.dofs.full(nj, c)) += fe(6 + c);
      }
    }
    // Same geometric correction for the facet shells: R += kg(n_e) u_e.
    for (const ShellElem& se : asm_.shells) {
      const Vec3& s_n = Sn[se.index];
      if (s_n.squaredNorm() == 0.0) continue;
      Vec18 ue;
      for (int k = 0; k < 3; ++k)
        for (int c = 0; c < 6; ++c)
          ue(6 * k + c) = ufull(asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c));
      Vec18 fe = se.el.kg(s_n(0), s_n(1), s_n(2)) * ue;
      for (int k = 0; k < 3; ++k)
        for (int c = 0; c < 6; ++c)
          fg(asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c)) += fe(6 * k + c);
    }
    for (Reaction& re : out.result.reactions) {
      int ni = model.node_index(re.node);
      for (int c = 0; c < 6; ++c) re.R(c) += fg(asm_.dofs.full(ni, c));
    }
  }
  hooks_report(opt.hooks, 1.0, "done");
  return out;
}

// ---------------- arc length ----------------

ArcTrace arc_length_trace(const Model& model, Id case_id, const ArcLengthOptions& opt) {
  require(model.continua.empty(), "arc_length_trace: continuum elements are not supported (linear statics, modal and dynamics take them)");
  CorotSystem sys = CorotSystem::build(model, case_id, opt.measure);
  require(sys.ext_moments.empty(),
          "arc length: nodal moments are not supported (force loading only)");
  require(sys.bar_loads.empty(),
          "arc length: bar span loads are not supported (their equivalent nodal "
          "vector is state-dependent; the continuation assumes a constant "
          "reference load — use nodal forces)");
  require(!sys.slave_loaded, "arc length: a nodal force on a rigid link's slave makes the reference load "
                             "state-dependent (load the master)");
  int n = sys.n();
  VecX u = VecX::Zero(n);
  double lambda = 0;
  double dl = opt.dl0;
  VecX du_prev = VecX::Zero(n);

  ArcTrace trace;
  {
    ArcPoint p0;
    p0.lambda = 0;
    p0.u = sys.full(u);
    trace.points.push_back(p0);
  }
  int prev_neg = 0;

  for (int step = 0; step < opt.max_steps; ++step) {
    hooks_check(opt.hooks);
    if (opt.hooks)
      opt.hooks->report(static_cast<double>(step) / std::max(opt.max_steps, 1),
                        "step " + std::to_string(step + 1) + "/" + std::to_string(opt.max_steps));
    bool converged = false;
    VecX u_new;
    double lam_new = 0;
    VecX Du;
    int iters_used = 0;  // corrector iterations of the accepted attempt
    for (int cut = 0; cut <= opt.max_cuts && !converged; ++cut) {
      try {
        SpMat Kt = sys.tangent(u);
        FactorCache fac;
        fac.factorize(Kt);
        VecX duF = fac.solve(sys.Fref);
        double dlam = dl / duF.norm();
        if (step > 0 && du_prev.dot(duF) < 0) dlam = -dlam;
        double lam = lambda + dlam;
        VecX uu = u + dlam * duF;
        Du = uu - u;
        bool ok = false;
        for (int it = 0; it < opt.max_iterations; ++it) {
          double elem_scale = 0;
          VecX fint = sys.internal_force(uu, &elem_scale);
          VecX r = lam * sys.Fref - fint;
          // The scale never collapses: at lambda ~ 0 on a snap-through path
          // both |lambda F| and the assembled |f_int| vanish — while the
          // members stay highly stressed (element force level), or at an
          // exactly unstrained mirror state (reference load |F|) — and a
          // tolerance relative to ~0 sits below the round-off of the element
          // forces, cutting the step to failure.
          double scale = std::max({std::fabs(lam) * sys.Fref.lpNorm<Eigen::Infinity>(),
                                   fint.lpNorm<Eigen::Infinity>(), elem_scale,
                                   sys.Fref.lpNorm<Eigen::Infinity>(), 1e-30});
          if (r.lpNorm<Eigen::Infinity>() <= std::max(opt.tol_rel * scale, sys.force_floor)) {
            ok = true;
            iters_used = it;
            break;
          }
          SpMat Kt2 = sys.tangent(uu);
          FactorCache f2;
          f2.factorize(Kt2);
          VecX dur = f2.solve(r);
          VecX duF2 = f2.solve(sys.Fref);
          VecX t = Du + dur;
          double a = duF2.dot(duF2);
          double b = 2 * duF2.dot(t);
          double c = t.dot(t) - dl * dl;
          double disc = b * b - 4 * a * c;
          if (disc < 0) break;  // step cut
          double r1 = (-b + std::sqrt(disc)) / (2 * a);
          double r2 = (-b - std::sqrt(disc)) / (2 * a);
          // Crisfield root selection: maximize dot(Du_new, Du).
          VecX c1 = t + r1 * duF2, c2 = t + r2 * duF2;
          double pick = c1.dot(Du) >= c2.dot(Du) ? r1 : r2;
          uu += dur + pick * duF2;
          lam += pick;
          Du = uu - u;
        }
        if (ok) {
          u_new = uu;
          lam_new = lam;
          converged = true;
        } else {
          dl /= 2;
        }
      } catch (const Error&) {
        dl /= 2;
      }
    }
    require(converged, "arc length: step failed after cuts");

    du_prev = u_new - u;
    u = u_new;
    lambda = lam_new;
    sys.recanonicalize(u);

    ArcPoint p;
    p.lambda = lambda;
    p.u = sys.full(u);
    try {
      SpMat Kt = sys.tangent(u);
      FactorCache fac;
      fac.factorize(Kt);
      p.n_neg = fac.inertia().n_neg;
    } catch (const Error&) {
      p.n_neg = -1;  // singular at the point
    }
    p.crossing = (p.n_neg != prev_neg);
    prev_neg = p.n_neg;
    trace.points.push_back(p);

    if (opt.desired_iterations > 0) {
      // Adaptive radius: scale by sqrt(desired / used), clamped to [1/2, 2]
      // per step and to [dl0/16, 4 dl0] overall.
      double f = std::clamp(std::sqrt(static_cast<double>(opt.desired_iterations) / std::max(1, iters_used)),
                            0.5, 2.0);
      dl = std::clamp(dl * f, opt.dl0 / 16, 4 * opt.dl0);
    }
    if (opt.stop_lambda > 0 && std::fabs(lambda) >= opt.stop_lambda) break;
    if (opt.stop_disp_inf > 0 && u.lpNorm<Eigen::Infinity>() >= opt.stop_disp_inf) break;
  }
  hooks_report(opt.hooks, 1.0, "done");
  return trace;
}

// ---------------- ruin load ----------------

RuinResult ruin_load(const Model& model, Id case_id, const RuinOptions& opt) {
  require(model.continua.empty(), "ruin_load: continuum elements are not supported (linear statics, modal and dynamics take them)");
  CorotSystem sys = CorotSystem::build(model, case_id, opt.corot.measure);
  VecX u = VecX::Zero(sys.n());

  // A quasi-static ruin under load control shows up as a Newton stall, an
  // indefinite/singular converged tangent, or a SNAP: Newton lands on a far
  // stable branch after a dynamic jump. The snap is certified by the tangent
  // inertia at the MIDPOINT between the warm state and the converged state —
  // a smooth quasi-static continuation stays in the stable region, a snap
  // crosses the unstable one.
  auto inertia_neg = [&](const VecX& state) {
    SpMat Kt = sys.tangent(state);
    FactorCache fac;
    fac.factorize(Kt);
    return fac.inertia().n_neg;
  };
  auto stable_at = [&](double lam, VecX& state, std::string& oracle) {
    VecX trial = state;
    try {
      newton_solve(sys, lam, trial, opt.corot);
    } catch (const Cancelled&) {
      throw;
    } catch (const Error&) {
      oracle = "stall";
      return false;
    }
    try {
      if (inertia_neg(trial) > 0) {
        oracle = "indefinite";
        return false;
      }
      VecX mid = 0.5 * (state + trial);
      if (inertia_neg(mid) > 0) {
        oracle = "snap";
        return false;
      }
    } catch (const Error&) {
      oracle = "singular";
      return false;
    }
    state = trial;
    return true;
  };

  RuinResult out;
  std::string oracle;
  double lo = 0, hi = opt.lambda_init;
  VecX ulo = u;
  const SolveHooks* hk = opt.corot.hooks;
  int probe = 0;
  while (true) {
    hooks_check(hk);
    hooks_report(hk, 0.0, "bracketing probe " + std::to_string(++probe));
    VecX state = ulo;
    if (!stable_at(hi, state, oracle)) break;
    lo = hi;
    ulo = state;
    hi *= opt.growth;
    require(hi <= opt.lambda_max, "ruin: no failure below lambda cap");
  }
  // Bisection progress: the halvings done over the halvings needed to reach
  // the tolerance from the bracket, kept non-decreasing.
  const double w0 = hi - lo;
  double frac = 0.5;
  while (hi - lo > opt.tol_rel * hi) {
    hooks_check(hk);
    if (hk) {
      const double need = std::log2(std::max(w0 / (opt.tol_rel * hi), 2.0));
      frac = std::max(frac, 0.5 + 0.5 * std::min(1.0, std::log2(w0 / (hi - lo)) / need));
      hk->report(frac, "bisection probe " + std::to_string(++probe));
    }
    double mid = 0.5 * (lo + hi);
    VecX state = ulo;
    std::string o2;
    if (stable_at(mid, state, o2)) {
      lo = mid;
      ulo = state;
    } else {
      hi = mid;
      oracle = o2;
    }
  }
  out.lambda_lower = lo;
  out.lambda_upper = hi;
  out.oracle = oracle;
  hooks_report(hk, 1.0, "done");
  return out;
}

CorotElementState corot_element_state(const Model& model, Id bar, const Vec12& q,
                                      CorotRotationMeasure measure) {
  const Assembled asm_ = assemble_elements(model);
  const std::size_t e = static_cast<std::size_t>(model.bar_index(bar));
  require(model.bars[e].kind != BarKind::Cable, "corotational: cables belong to solve_cable_static");
  const CorotElem ce = make_corot_elem(model, asm_, e, measure);
  CorotElementState st;
  st.energy = ce.energy(q);
  st.force = ce.grad(q);
  st.tangent = ce.hess(q);
  return st;
}

CorotStagedResults solve_corotational_staged(const Model& model, const std::vector<CorotStage>& stages,
                                             const CorotOptions& opt) {
  model.validate();
  require(model.continua.empty(), "solve_corotational_staged: continuum elements are not supported");
  require(!stages.empty(), "staged corotational: no stage");
  const std::size_t nn = model.nodes.size();
  std::vector<char> active_bar(model.bars.size(), 0), placed(nn, 0);
  std::vector<Vec6> D(nn, Vec6::Zero());
  std::vector<Load> loads;
  std::map<Id, double> l0;
  CorotStagedResults out;
  for (std::size_t k = 0; k < stages.size(); ++k) {
    const CorotStage& st = stages[k];
    std::vector<std::size_t> fresh;
    for (Id b : st.activate) {
      const std::size_t e = static_cast<std::size_t>(model.bar_index(b));
      require(!active_bar[e], "staged corotational: bar " + std::to_string(b) + " activated twice");
      active_bar[e] = 1;
      fresh.push_back(e);
    }
    for (const auto& [c, L] : st.cable_l0) {
      const Bar& cb = model.bar(c);
      require(cb.kind == BarKind::Cable && std::isfinite(L) && L > 0, "staged corotational: cable_l0 needs a cable and L0 > 0");
      l0[c] = L;
    }
    if (st.case_id != 0)
      for (const Load& l : model.load_case(st.case_id).loads) loads.push_back(l);
    // place the new nodes: rigid continuation through the new bars
    bool grew = true;
    std::vector<char> todo(fresh.size(), 1);
    while (grew) {
      grew = false;
      for (std::size_t q = 0; q < fresh.size(); ++q) {
        if (!todo[q]) continue;
        const Bar& b = model.bars[fresh[q]];
        const std::size_t i = static_cast<std::size_t>(model.node_index(b.ni)), j = static_cast<std::size_t>(model.node_index(b.nj));
        const bool pi = placed[i] || k == 0, pj = placed[j] || k == 0;
        if (pi && pj) {
          placed[i] = placed[j] = 1;
          todo[q] = 0;
          continue;
        }
        if (!pi && !pj) continue;
        const std::size_t from = pi ? i : j, to = pi ? j : i;
        const Mat3 R = rot_exp(D[from].tail<3>());
        const Vec3 X = model.nodes[to].x - model.nodes[from].x;
        D[to].head<3>() = D[from].head<3>() + (R - Mat3::Identity()) * X;
        D[to].tail<3>() = D[from].tail<3>();
        placed[to] = 1;
        todo[q] = 0;
        grew = true;
      }
    }
    for (std::size_t q = 0; q < fresh.size(); ++q)
      require(!todo[q], "staged corotational: a bar activated with neither end in place (activate it with its neighbours)");
    // the stage's structure
    Model mk = model;
    mk.bars.clear();
    for (std::size_t e = 0; e < model.bars.size(); ++e)
      if (active_bar[e]) {
        Bar b = model.bars[e];
        if (b.kind == BarKind::Cable && l0.count(b.id)) b.cable->unstressed_length = l0[b.id];
        mk.bars.push_back(b);
      }
    for (const Load& l : loads)
      if (l.bar != 0 && l.kind != Load::Kind::NodalForce && l.kind != Load::Kind::SelfWeight)
        require(active_bar[static_cast<std::size_t>(model.bar_index(l.bar))], "staged corotational: a load on an inactive bar");
    mk.cases.clear();
    mk.add_case("stage", loads);
    CorotSystem sys = CorotSystem::build(mk, mk.cases.front().id, opt.measure, true);
    require(!sys.linked(), "staged corotational: rigid links, constraints and skew supports are not supported");
    const DofMap& d = sys.asm_.dofs;
    VecX u = VecX::Zero(sys.n());
    for (int ni = 0; ni < d.n_nodes; ++ni)
      for (int c = 0; c < 6; ++c) {
        const int f = d.free_of(ni, c);
        if (f >= 0) u(f) = D[static_cast<std::size_t>(ni)](c);
      }
    out.iterations += newton_solve(sys, 1.0, u, opt);
    for (int ni = 0; ni < d.n_nodes; ++ni)
      for (int c = 0; c < 6; ++c) {
        const int f = d.free_of(ni, c);
        if (f >= 0) D[static_cast<std::size_t>(ni)](c) = u(f);
      }
    std::vector<Vec6> snap(nn, Vec6::Zero());
    for (std::size_t ni = 0; ni < nn; ++ni)
      if (placed[ni] || k == 0) snap[ni] = D[ni];
    out.node_disp.push_back(snap);
  }
  return out;
}

CorotDynResult corot_newmark(const Model& model, Id case_id, const std::function<double(double)>& lam,
                             const CorotDynOptions& opt) {
  require(model.continua.empty(), "corot_newmark: continuum elements are not supported");
  require(opt.dt > 0 && opt.steps > 0 && opt.beta > 0 && opt.gamma >= 0.5, "corot_newmark: dt, steps, beta > 0, gamma >= 1/2");
  require(static_cast<bool>(lam), "corot_newmark: a load function is required");
  CorotSystem sys = CorotSystem::build(model, case_id, opt.measure, true);
  require(!sys.linked(), "corot_newmark: rigid links, constraints and skew supports are not supported");
  const int n = sys.n();
  const SpMat M = assemble_mass(model, sys.asm_, opt.lumped_mass);
  require(M.rows() == n, "corot_newmark: mass size");
  VecX u = VecX::Zero(n), v = VecX::Zero(n), a = VecX::Zero(n);
  if (!opt.initial_displacement.empty()) {
    require(opt.initial_displacement.size() == model.nodes.size(),
            "corot_newmark: initial_displacement must hold one entry per model node");
    for (std::size_t ni = 0; ni < model.nodes.size(); ++ni)
      for (int c = 0; c < 6; ++c) {
        const int f = sys.asm_.dofs.free_of(static_cast<int>(ni), c);
        if (f >= 0) u(f) = opt.initial_displacement[ni](c);
        else
          require(opt.initial_displacement[ni](c) == 0.0,
                  "corot_newmark: initial displacement on a blocked DOF of node " + std::to_string(model.nodes[ni].id));
      }
  }
  const SpMat K0 = sys.tangent(u);
  const SpMat C = (opt.rayleigh_a * M + opt.rayleigh_b * K0).pruned();
  const bool galpha = opt.rho_inf >= 0;
  require(!galpha || opt.rho_inf <= 1, "corot_newmark: rho_inf in [0, 1] (or -1: Newmark)");
  const double am = galpha ? (2 * opt.rho_inf - 1) / (opt.rho_inf + 1) : 0.0;
  const double af = galpha ? opt.rho_inf / (opt.rho_inf + 1) : 0.0;
  const double dt = opt.dt, be = galpha ? 0.25 * (1 - am + af) * (1 - am + af) : opt.beta,
               ga = galpha ? 0.5 - am + af : opt.gamma;
  for (const auto& [node, vel] : opt.initial_velocity)
    for (int c = 0; c < 6; ++c) {
      const int f = sys.asm_.dofs.free_of(model.node_index(node), c);
      if (f >= 0) v(f) = vel(c);
    }
  {
    const VecX r0 = sys.external_force(u, lam(0.0)) - sys.internal_force(u) - C * v;
    const VecX md = M.diagonal();
    bool diag = true;
    for (int k = 0; k < M.outerSize() && diag; ++k)
      for (SpMat::InnerIterator it(M, k); it; ++it)
        if (it.row() != it.col() && it.value() != 0) diag = false;
    if (diag) {
      for (int i = 0; i < n; ++i) a(i) = md(i) > 0 ? r0(i) / md(i) : 0.0;
    } else {
      Eigen::SimplicialLDLT<SpMat> ld(M);
      if (ld.info() == Eigen::Success) a = ld.solve(r0);
    }
  }
  std::vector<int> obs;
  for (const auto& [node, c] : opt.observe) obs.push_back(sys.asm_.dofs.free_of(model.node_index(node), c));
  CorotDynResult out;
  auto record = [&](double t) {
    out.t.push_back(t);
    std::vector<double> row;
    for (int f : obs) row.push_back(f >= 0 ? u(f) : 0.0);
    out.u.push_back(row);
  };
  record(0.0);
  VecX Gn = sys.external_force(u, lam(0.0)) - sys.internal_force(u);  // F - f at the start of the step
  for (int s = 1; s <= opt.steps; ++s) {
    const double t = s * dt, l = lam(t);
    const VecX un = u, vn = v, an = a;
    VecX ut = un + dt * vn + dt * dt * (0.5 - be) * an;  // predictor (a = 0 in the step)
    double r_prev = std::numeric_limits<double>::infinity();
    for (int it = 0;; ++it) {
      require(it < opt.max_iterations, "corot_newmark: Newton did not converge at step " + std::to_string(s));
      const VecX at = (ut - un - dt * vn) / (be * dt * dt) - (0.5 / be - 1) * an;
      const VecX vt = vn + dt * ((1 - ga) * an + ga * at);
      const VecX F = sys.external_force(ut, l), fi = sys.internal_force(ut);
      const VecX Ma = M * ((1 - am) * at + am * an), Cv = C * ((1 - af) * vt + af * vn);
      const VecX r = (1 - af) * (F - fi) + af * Gn - Ma - Cv;
      const double scale = std::max({F.lpNorm<Eigen::Infinity>(), fi.lpNorm<Eigen::Infinity>(), Ma.lpNorm<Eigen::Infinity>(), 1e-30});
      const double rn = r.lpNorm<Eigen::Infinity>();
      // converged, at the force floor, or stalled on round-off (no decrease
      // within a thousand times the tolerance)
      const bool stalled = it > 1 && rn >= 0.5 * r_prev && rn <= 1e3 * opt.tol_rel * scale;
      r_prev = rn;
      if (rn <= opt.tol_rel * scale || (it > 0 && rn <= sys.force_floor) || stalled) {
        u = ut;
        a = at;
        v = vt;
        Gn = F - fi;
        break;
      }
      ++out.iterations;
      SpMat Kt = sys.tangent(ut);
      if (sys.has_load_stiffness()) Kt = (Kt - l * sys.load_stiffness(ut)).eval();
      const SpMat Ke = ((1 - af) * Kt + ((1 - am) / (be * dt * dt)) * M + ((1 - af) * ga / (be * dt)) * C).eval();
      Eigen::SparseLU<SpMat, Eigen::COLAMDOrdering<int>> lu;
      lu.compute(Ke);
      require(lu.info() == Eigen::Success, "corot_newmark: singular effective tangent");
      const VecX du = lu.solve(r);
      require(du.allFinite(), "corot_newmark: singular effective tangent");
      ut += du;
    }
    record(t);
  }
  out.node_disp.assign(static_cast<std::size_t>(sys.asm_.dofs.n_nodes), Vec6::Zero());
  out.node_vel.assign(static_cast<std::size_t>(sys.asm_.dofs.n_nodes), Vec6::Zero());
  for (int ni = 0; ni < sys.asm_.dofs.n_nodes; ++ni)
    for (int c = 0; c < 6; ++c) {
      const int fi = sys.asm_.dofs.free_ix[static_cast<std::size_t>(sys.asm_.dofs.full(ni, c))];
      if (fi >= 0) {
        out.node_disp[static_cast<std::size_t>(ni)](c) = u(fi);
        out.node_vel[static_cast<std::size_t>(ni)](c) = v(fi);
      }
    }
  return out;
}

CorotTangentProbe corot_tangent_probe(const Model& model, Id case_id, const VecX& u, double lam,
                                      CorotRotationMeasure measure, bool allow_cables) {
  CorotSystem sys = CorotSystem::build(model, case_id, measure, allow_cables);
  require(!sys.linked(), "corot probe: rigid links are not supported (the probe works on the free DOFs)");
  require(u.size() == sys.asm_.dofs.n_free, "corot probe: u must hold the free DOFs");
  CorotTangentProbe p;
  p.residual = sys.external_force(u, lam) - sys.internal_force(u);
  SpMat K = sys.tangent(u);
  if (sys.has_load_stiffness()) K = K - lam * sys.load_stiffness(u);
  p.tangent = MatX(K);
  return p;
}

}  // namespace keta

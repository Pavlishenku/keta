#include "keta/assembly.h"

#include <algorithm>

#include <functional>
#include <map>

#include "keta/cable.h"
#include "keta/parallel.h"
#include "keta/shell_tendon.h"
#include "keta/springmass.h"

namespace keta {

Mat12 bar_ke_assembled(const Model& model, const BarElement& el) {
  if (model.tapered.empty()) return el.ke();  // fast, bit-identical default path
  const TaperedBar* tb = model.tapered_of(el.bar().id);
  if (!tb) return el.ke();
  // No release is allowed on a tapered bar (Model::validate), so `k_local()`
  // IS the raw prismatic member block: subtracting it and adding the tapered
  // one keeps the eccentricity transform and the Winkler bed untouched.
  return el.ke() + el.T().transpose() * (tb->ke - el.k_local()) * el.T();
}

namespace {
Mat3 skew(const Vec3& a) {
  Mat3 S;
  S << 0, -a.z(), a.y(), a.z(), 0, -a.x(), -a.y(), a.x(), 0;
  return S;
}
}  // namespace

LinkElem link_element(const Model& model, const ElasticLink& l) {
  LinkElem e;
  e.id = l.id;
  e.ni = model.node_index(l.ni);
  e.nj = model.node_index(l.nj);
  const Vec3 xi = model.node(l.ni).x, xj = model.node(l.nj).x;
  const Vec3 P = xi + l.s * (xj - xi);
  Mat3 R = Mat3::Identity();
  if (l.axes)
    R = *l.axes;
  else if ((xj - xi).norm() > 1e-12)
    R = bar_rotation(xi, xj, l.roll);
  // global relative motion at P: u + theta x a = u - [a]x theta
  Eigen::Matrix<double, 6, 12> Bg = Eigen::Matrix<double, 6, 12>::Zero();
  Bg.block<3, 3>(0, 0) = -Mat3::Identity();
  Bg.block<3, 3>(0, 3) = skew(P - xi);
  Bg.block<3, 3>(0, 6) = Mat3::Identity();
  Bg.block<3, 3>(0, 9) = -skew(P - xj);
  Bg.block<3, 3>(3, 3) = -Mat3::Identity();
  Bg.block<3, 3>(3, 9) = Mat3::Identity();
  Mat6 R6 = Mat6::Zero();
  R6.block<3, 3>(0, 0) = R;
  R6.block<3, 3>(3, 3) = R;
  e.B = R6 * Bg;
  e.k = 0.5 * (l.k + l.k.transpose());
  e.K = e.B.transpose() * e.k * e.B;
  return e;
}

bool links_have_arms(const Model& model) {
  for (const RigidLink& l : model.rigid_links)
    for (Id s : l.slaves)
      if ((model.node(s).x - model.node(l.master).x).norm() > 0) return true;
  for (const ElasticLink& l : model.links)
    if ((model.node(l.nj).x - model.node(l.ni).x).norm() > 0) return true;
  return false;
}

namespace {
// Resolved master-slave expressions of the rigid links over non-slave full
// DOFs (chains substituted, cycles and double constraints rejected).
std::vector<LinkReduction::Expr> rigid_link_expressions(const Model& model, const DofMap& d) {
  const std::size_t nfull = static_cast<std::size_t>(kDofPerNode) * static_cast<std::size_t>(d.n_nodes);
  std::vector<std::vector<std::pair<int, double>>> raw(nfull);
  std::vector<char> is_slave(nfull, 0);
  std::vector<Id> owner(nfull, 0);
  for (const RigidLink& l : model.rigid_links) {
    const int m = model.node_index(l.master);
    const Vec3 xm = model.node(l.master).x;
    for (Id sid : l.slaves) {
      const int sn = model.node_index(sid);
      const Vec3 r = model.node(sid).x - xm;
      // theta x r rows: x: ty rz - tz ry; y: tz rx - tx rz; z: tx ry - ty rx
      const Mat3 C = -skew(r);  // (theta x r) = -[r]x theta
      for (int c = 0; c < 6; ++c) {
        if (!l.dofs[static_cast<std::size_t>(c)]) continue;
        const std::size_t f = static_cast<std::size_t>(d.full(sn, c));
        require(!is_slave[f], "rigid link " + std::to_string(l.id) + ": node " + std::to_string(sid) +
                                  " component " + std::to_string(c) + " is already constrained (link " +
                                  std::to_string(owner[f]) + ")");
        require(!d.blocked[f], "rigid link " + std::to_string(l.id) + ": slave node " + std::to_string(sid) +
                                   " component " + std::to_string(c) +
                                   " is supported — support the master, or join the support by an elastic link");
        is_slave[f] = 1;
        owner[f] = l.id;
        auto& e = raw[f];
        e.emplace_back(d.full(m, c), 1.0);
        if (c < 3)
          for (int a = 0; a < 3; ++a)
            if (C(c, a) != 0) e.emplace_back(d.full(m, 3 + a), C(c, a));
      }
    }
  }
  // general constraints: the first term's DOF eliminated
  auto eliminate = [&](std::size_t f, const std::string& who) {
    require(!is_slave[f], who + ": its eliminated DOF already carries a constraint");
    require(!d.blocked[f], who + ": its eliminated DOF is supported");
    is_slave[f] = 1;
  };
  for (const Constraint& c : model.constraints) {
    const std::string who = "constraint " + std::to_string(c.id);
    const ConstraintTerm& t0 = c.terms.front();
    const std::size_t f = static_cast<std::size_t>(d.full(model.node_index(t0.node), t0.comp));
    eliminate(f, who);
    for (std::size_t k = 1; k < c.terms.size(); ++k) {
      const ConstraintTerm& t = c.terms[k];
      if (t.coef != 0) raw[f].emplace_back(d.full(model.node_index(t.node), t.comp), -t.coef / t0.coef);
    }
  }
  // skew supports: a_k . u = 0 per blocked local component. The blocked rows
  // of the translations (then of the rotations) are brought to reduced row
  // echelon form over the node's still-free global components, the pivot
  // of each row the largest one left (partial pivoting): each row then
  // eliminates its pivot component and holds none of the other pivots, so
  // two blocked rotated directions (a guide and a longitudinal stop in a
  // plan frame) make no cyclic chain. A single row per group is the plain
  // largest-component elimination.
  for (const SkewSupport& sk : model.skew_supports) {
    const std::string who = "skew support " + std::to_string(sk.id);
    const int ni = model.node_index(sk.node);
    for (const int base : {0, 3}) {
      std::vector<Vec3> rows;
      std::vector<int> local;
      for (int k = base; k < base + 3; ++k)
        if (sk.fixed[static_cast<std::size_t>(k)]) {
          rows.push_back(sk.axes.row(k % 3).transpose());
          local.push_back(k);
        }
      std::vector<int> pivot(rows.size(), -1);
      std::array<bool, 3> used{};
      for (std::size_t r = 0; r < rows.size(); ++r) {
        int best = -1;
        for (int c = 0; c < 3; ++c) {
          const std::size_t f = static_cast<std::size_t>(d.full(ni, base + c));
          if (used[static_cast<std::size_t>(c)] || is_slave[f] || d.blocked[f] || std::fabs(rows[r](c)) < 1e-12)
            continue;
          if (best < 0 || std::fabs(rows[r](c)) > std::fabs(rows[r](best))) best = c;
        }
        require(best >= 0, who + ": no free global component left to eliminate for local component " +
                               std::to_string(local[r]));
        pivot[r] = best;
        used[static_cast<std::size_t>(best)] = true;
        rows[r] /= rows[r](best);
        for (std::size_t o = 0; o < rows.size(); ++o)
          if (o != r && rows[o](best) != 0) rows[o] -= rows[o](best) * rows[r];
      }
      for (std::size_t r = 0; r < rows.size(); ++r) {
        const std::size_t f = static_cast<std::size_t>(d.full(ni, base + pivot[r]));
        eliminate(f, who);
        for (int c = 0; c < 3; ++c)
          if (c != pivot[r] && rows[r](c) != 0) raw[f].emplace_back(d.full(ni, base + c), -rows[r](c));
      }
    }
  }
  std::vector<char> state(nfull, 0);  // 0 raw, 1 resolving, 2 resolved
  std::vector<std::vector<std::pair<int, double>>> res(nfull);
  std::function<const std::vector<std::pair<int, double>>&(int)> resolve =
      [&](int f) -> const std::vector<std::pair<int, double>>& {
    const std::size_t uf = static_cast<std::size_t>(f);
    if (state[uf] == 2) return res[uf];
    require(state[uf] == 0, "rigid links: cyclic master-slave chain through node index " +
                                std::to_string(f / kDofPerNode));
    state[uf] = 1;
    std::map<int, double> acc;
    for (const auto& [k, c] : raw[uf]) {
      if (is_slave[static_cast<std::size_t>(k)]) {
        for (const auto& [k2, c2] : resolve(k)) acc[k2] += c * c2;
      } else {
        acc[k] += c;
      }
    }
    for (const auto& [k, c] : acc)
      if (c != 0) res[uf].emplace_back(k, c);
    state[uf] = 2;
    return res[uf];
  };
  std::vector<LinkReduction::Expr> out;
  for (std::size_t f = 0; f < nfull; ++f)
    if (is_slave[f]) out.push_back({static_cast<int>(f), resolve(static_cast<int>(f))});
  return out;
}
}  // namespace

SpMat LinkReduction::reduce(const SpMat& A) const {
  SpMat r = (T.transpose() * A * T).pruned();
  r.makeCompressed();
  return r;
}
VecX LinkReduction::reduce(const VecX& f) const { return T.transpose() * f; }
MatX LinkReduction::reduce(const MatX& f) const { return T.transpose() * f; }
VecX LinkReduction::expand(const VecX& ur) const { return T * ur; }

Assembled assemble_elements(const Model& model, bool links_supported, const std::vector<double>* prestress) {
  model.validate();
  require(model.composite.empty(),
          "composite bars (Model::composite) are handled by the phasing — lower them for this "
          "analysis with composite_lowered(model)");
  for (const ElasticLink& l : model.links)
    require(!l.has_law(), "elastic link " + std::to_string(l.id) +
                              ": a non-linear law is solved by solve_static_link_laws (keta/link_nl.h)");
  require(!model.has_links() || links_supported,
          "node-to-node links (rigid or elastic) are handled by the linear static, modal, spectrum, buckling and phasing "
          "analyses only — this analysis does not support them");
  Assembled out;
  out.elements.reserve(model.bars.size());
  for (const Bar& b : model.bars) out.elements.emplace_back(model, b);
  out.shells = shell_elems(model);  // empty in a bar-only model

  int nn = static_cast<int>(model.nodes.size());
  DofMap& d = out.dofs;
  d.n_nodes = nn;
  d.blocked.assign(static_cast<std::size_t>(kDofPerNode) * nn, false);
  out.spring.assign(static_cast<std::size_t>(kDofPerNode) * nn, 0.0);

  for (const Support& s : model.supports) {
    int ni = model.node_index(s.node);
    for (int c = 0; c < 6; ++c) {
      if (s.fixed[static_cast<std::size_t>(c)]) d.blocked[static_cast<std::size_t>(d.full(ni, c))] = true;
      if (s.spring[static_cast<std::size_t>(c)] != 0)
        out.spring[static_cast<std::size_t>(d.full(ni, c))] += s.spring[static_cast<std::size_t>(c)];
    }
    // ST1 APPUI EL: full 6x6 elastic support, rotated to the global frame.
    if (s.spring6 && s.spring6->squaredNorm() > 0)
      out.spring6.emplace_back(ni, support_stiffness_global(s));
  }

  // Zero-stiffness auto-fix: accumulate element diagonal contributions; a DOF
  // with an exactly-zero total (element diagonals are non-negative) and no
  // spring is eliminated as if supported, and reported.
  std::vector<double> diag(static_cast<std::size_t>(kDofPerNode) * nn, 0.0);
  // A full 6x6 elastic support is stiffness too: its diagonal keeps the
  // auto-fix pass from eliminating the DOFs it braces.
  for (const auto& [sni, k6] : out.spring6)
    for (int c = 0; c < 6; ++c) diag[static_cast<std::size_t>(d.full(sni, c))] += k6(c, c);
  for (std::size_t e = 0; e < out.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) {
      // A cable contributes translational stiffness at solve time; pre-mark
      // its anchors' translations as stiffness sources (rotations stay free
      // to be auto-fixed — a cable carries no moment).
      for (Id nid : {b.ni, b.nj}) {
        int ni = model.node_index(nid);
        for (int c = 0; c < 3; ++c) diag[static_cast<std::size_t>(d.full(ni, c))] += 1.0;
      }
      continue;
    }
    Mat12 k = bar_ke_assembled(model, out.elements[e]);
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (int c = 0; c < 6; ++c) {
      diag[static_cast<std::size_t>(d.full(ni, c))] += k(c, c);
      diag[static_cast<std::size_t>(d.full(nj, c))] += k(6 + c, 6 + c);
    }
    if (prestress && e < prestress->size() && (*prestress)[e] != 0.0)
      for (int c = 0; c < 3; ++c) {
        const double g = std::fabs((*prestress)[e]) / out.elements[e].L();
        diag[static_cast<std::size_t>(d.full(ni, c))] += g;
        diag[static_cast<std::size_t>(d.full(nj, c))] += g;
      }
  }
  // Facet shells contribute stiffness on all six DOFs of their three nodes
  // (the Allman membrane carries real rz stiffness), so a shell-only node
  // keeps every DOF.
  for (const ShellElem& se : out.shells) {
    const Mat18& k = se.el.ke();
    for (int c = 0; c < 3; ++c)
      for (int q = 0; q < 6; ++q)
        diag[static_cast<std::size_t>(d.full(se.nodes[static_cast<std::size_t>(c)], q))] +=
            k(6 * c + q, 6 * c + q);
  }
  // Continuum elements: stiffness on the translations they use
  for (const Continuum& c : model.continua) {
    out.continua.push_back(continuum_kernel(model, c));
    const ContinuumKernel& k = out.continua.back();
    for (std::size_t i = 0; i < k.nodes.size(); ++i)
      for (int a = 0; a < k.dim; ++a)
        diag[static_cast<std::size_t>(d.full(k.nodes[i], a))] += k.K(k.dim * static_cast<int>(i) + a, k.dim * static_cast<int>(i) + a);
  }
  for (const PointMass& pm : model.point_masses) {
    // A point mass is not stiffness; it does not prevent auto-fix.
    (void)model.node(pm.node);
  }
  // Elastic links are stiffness on their two nodes.
  for (const ElasticLink& l : model.links) {
    out.links.push_back(link_element(model, l));
    const LinkElem& e = out.links.back();
    for (int c = 0; c < 6; ++c) {
      diag[static_cast<std::size_t>(d.full(e.ni, c))] += e.K(c, c);
      diag[static_cast<std::size_t>(d.full(e.nj, c))] += e.K(6 + c, 6 + c);
    }
  }
  // Rigid links: the slave DOFs and the master DOFs their expressions use
  // are kept (their stiffness arrives through T^t K T); a reduced DOF left
  // without stiffness then fails the factorization (typed), never silently.
  std::vector<char> keep(static_cast<std::size_t>(kDofPerNode) * nn, 0);
  std::vector<LinkReduction::Expr> exprs;
  if (!model.rigid_links.empty() || !model.constraints.empty() || !model.skew_supports.empty()) {
    exprs = rigid_link_expressions(model, d);
    for (const auto& e : exprs) {
      keep[static_cast<std::size_t>(e.slave)] = 1;
      for (const auto& [k, c] : e.terms) keep[static_cast<std::size_t>(k)] = 1;
    }
  }

  d.free_ix.assign(static_cast<std::size_t>(kDofPerNode) * nn, -1);
  int nf = 0;
  for (int ni = 0; ni < nn; ++ni) {
    for (int c = 0; c < 6; ++c) {
      std::size_t f = static_cast<std::size_t>(d.full(ni, c));
      if (d.blocked[f]) continue;
      if (diag[f] == 0.0 && out.spring[f] == 0.0 && !keep[f]) {
        d.auto_fixed.push_back({model.nodes[static_cast<std::size_t>(ni)].id, c});
        d.blocked[f] = false;  // auto-fixed, distinct from user-blocked; both eliminated
        continue;
      }
      d.free_ix[f] = nf++;
    }
  }
  d.n_free = nf;  // nf == 0 is legal (e.g. a fully clamped single bar)
  if (!exprs.empty()) {
    LinkReduction& red = out.red;
    red.active = true;
    std::vector<char> slave_free(static_cast<std::size_t>(nf), 0);
    for (const auto& e : exprs) slave_free[static_cast<std::size_t>(d.free_ix[static_cast<std::size_t>(e.slave)])] = 1;
    std::vector<int> red_ix(static_cast<std::size_t>(nf), -1);
    int nr = 0;
    for (int i = 0; i < nf; ++i)
      if (!slave_free[static_cast<std::size_t>(i)]) red_ix[static_cast<std::size_t>(i)] = nr++;
    std::vector<Triplet> t;
    for (int i = 0; i < nf; ++i)
      if (red_ix[static_cast<std::size_t>(i)] >= 0) t.emplace_back(i, red_ix[static_cast<std::size_t>(i)], 1.0);
    for (const auto& e : exprs) {
      const int i = d.free_ix[static_cast<std::size_t>(e.slave)];
      for (const auto& [k, c] : e.terms) {
        const int fk = d.free_ix[static_cast<std::size_t>(k)];
        if (fk >= 0) t.emplace_back(i, red_ix[static_cast<std::size_t>(fk)], c);
      }
    }
    red.n_red = nr;
    red.T.resize(nf, nr);
    red.T.setFromTriplets(t.begin(), t.end());
    red.T.makeCompressed();
    red.slaves = std::move(exprs);
  }
  return out;
}

namespace {

template <class Emit>
void scatter12(const Model& model, const DofMap& d, const Bar& b, const Mat12& k, Emit&& emit) {
  int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
  int gdof[12];
  for (int c = 0; c < 6; ++c) {
    gdof[c] = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
    gdof[6 + c] = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
  }
  for (int a = 0; a < 12; ++a) {
    if (gdof[a] < 0) continue;
    for (int bb = 0; bb < 12; ++bb) {
      if (gdof[bb] < 0) continue;
      double v = k(a, bb);
      if (v != 0.0) emit(gdof[a], gdof[bb], v);
    }
  }
}

// The 18 full-DOF indices of a facet shell, element DOF order
// ([ux,uy,uz,rx,ry,rz] per node, nodes in winding order).
void shell_gdofs(const DofMap& d, const ShellElem& se, int gdof[18]) {
  for (int k = 0; k < 3; ++k)
    for (int c = 0; c < 6; ++c)
      gdof[6 * k + c] = d.full(se.nodes[static_cast<std::size_t>(k)], c);
}

template <class Emit>
void scatter18(const DofMap& d, const ShellElem& se, const Mat18& k, Emit&& emit) {
  int full[18];
  shell_gdofs(d, se, full);
  int free_ix[18];
  for (int a = 0; a < 18; ++a) free_ix[a] = d.free_ix[static_cast<std::size_t>(full[a])];
  for (int a = 0; a < 18; ++a) {
    if (free_ix[a] < 0) continue;
    for (int b = 0; b < 18; ++b) {
      if (free_ix[b] < 0) continue;
      double v = k(a, b);
      if (v != 0.0) emit(free_ix[a], free_ix[b], v);
    }
  }
}

}  // namespace

namespace {
void scatter_continuum(const DofMap& d, const ContinuumKernel& k, const MatX& A, std::vector<Triplet>& trips) {
  const int n = static_cast<int>(k.nodes.size());
  for (int i = 0; i < n; ++i)
    for (int a = 0; a < k.dim; ++a) {
      const int fi = d.free_of(k.nodes[static_cast<std::size_t>(i)], a);
      if (fi < 0) continue;
      for (int j = 0; j < n; ++j)
        for (int b = 0; b < k.dim; ++b) {
          const int fj = d.free_of(k.nodes[static_cast<std::size_t>(j)], b);
          const double v = A(k.dim * i + a, k.dim * j + b);
          if (fj >= 0 && v != 0) trips.emplace_back(fi, fj, v);
        }
    }
}
}  // namespace

SpMat assemble_stiffness(const Model& model, const Assembled& asm_, bool cables_added_by_caller) {
  std::vector<Triplet> trips;
  trips.reserve(asm_.elements.size() * 100 + 16);
  for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) {
      require(cables_added_by_caller,
              "bar " + std::to_string(b.id) +
                  " is a catenary cable, a nonlinear element a linear analysis cannot take: solve the cables with "
                  "solve_cable_static or solve_corotational, linearize about that state (tangent_of_cable_state) or "
                  "replace them by their Ernst trusses (ernst_substitute)");
      continue;  // the caller adds the cables
    }
    scatter12(model, asm_.dofs, b, bar_ke_assembled(model, asm_.elements[e]),
              [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
  }
  // Facet shells, AFTER the bar loop, in model.shells order.
  for (const ShellElem& se : asm_.shells)
    scatter18(asm_.dofs, se, se.el.ke(),
              [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
  for (const ContinuumKernel& k : asm_.continua) scatter_continuum(asm_.dofs, k, k.K, trips);
  for (std::size_t f = 0; f < asm_.spring.size(); ++f) {
    if (asm_.spring[f] != 0.0) {
      int fi = asm_.dofs.free_ix[f];
      if (fi >= 0) trips.emplace_back(fi, fi, asm_.spring[f]);
    }
  }
  scatter_full_springs(asm_, trips);  // ST1 APPUI EL
  for (const LinkElem& e : asm_.links) {  // elastic links
    int g[12];
    for (int c = 0; c < 6; ++c) {
      g[c] = asm_.dofs.free_of(e.ni, c);
      g[6 + c] = asm_.dofs.free_of(e.nj, c);
    }
    for (int a = 0; a < 12; ++a) {
      if (g[a] < 0) continue;
      for (int b = 0; b < 12; ++b)
        if (g[b] >= 0 && e.K(a, b) != 0.0) trips.emplace_back(g[a], g[b], e.K(a, b));
    }
  }
  SpMat K(asm_.dofs.n_free, asm_.dofs.n_free);
  K.setFromTriplets(trips.begin(), trips.end());
  return K;
}

SpMat assemble_mass(const Model& model, const Assembled& asm_, bool lumped) {
  std::vector<Triplet> trips;
  for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) continue;
    if (bar_mass_is_null(model, b.id)) continue;  // ST1 MASSE NUL
    if (model.bar_mass_selective && !bar_mass_is_auto(model, b.id)) continue;
    Mat12 me = asm_.elements[e].me(bar_mass_lumped(model, b.id, lumped));
    if (bar_mass_full_diagonal(model, b.id)) {
      // MASSE DI: rotary term of the HRZ lumping, M L^2 / 78 at each end, on
      // the rotation about local z (PLANE: ST1's cantilever mode 4 to
      // 0.00 %). About local y it lengthens spatial modes away from
      // ST1 (0.4 % -> 1.1 %), so it is left out there.
      const double L = asm_.elements[e].L(), M = 2 * me(1, 1);
      for (int k : {5, 11}) me(k, k) = M * L * L / 78;
    }
    scatter12(model, asm_.dofs, b, me,
              [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
  }
  for (const PointMass& pm : model.point_masses) {
    int ni = model.node_index(pm.node);
    for (int c = 0; c < 3; ++c) {
      int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
      if (fi >= 0 && pm.mass != 0) trips.emplace_back(fi, fi, pm.mass);
      int fr = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, 3 + c))];
      if (fr >= 0 && pm.inertia[c] != 0) trips.emplace_back(fr, fr, pm.inertia[c]);
    }
  }
  // Catenary cables: the mass of the unstrained cable, w L0 / g,
  // lumped half at each end on the three translations — the curved element
  // has no consistent mass here (mesh-convergent, like the beaded string).
  {
    bool any_cable = false;
    for (const Bar& b : model.bars) any_cable |= b.kind == BarKind::Cable;
    if (any_cable) {
      const std::vector<double> l0 = resolve_cable_l0(model);
      for (std::size_t e = 0; e < model.bars.size(); ++e) {
        const Bar& b = model.bars[e];
        if (b.kind != BarKind::Cable || bar_mass_is_null(model, b.id)) continue;
        const double half = 0.5 * cable_ea_w(model, b).second * l0[e] / kGravity;
        if (half == 0.0) continue;
        for (Id nid : {b.ni, b.nj}) {
          const int ni = model.node_index(nid);
          for (int c = 0; c < 3; ++c) {
            const int fi = asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))];
            if (fi >= 0) trips.emplace_back(fi, fi, half);
          }
        }
      }
    }
  }
  scatter_bar_line_masses(model, asm_, trips);  // ST1 MASSE UNI / LIN / PARABOL
  scatter_nodal_masses(model, asm_, trips);     // ST1 MASSE PROPRE NOEUD
  // Facet shells, AFTER the bar loop, in model.shells order (consistent
  // rho*t*A + rotary rho*t^3/12, no drilling mass — see shell.h).
  for (const ShellElem& se : asm_.shells) {
    if (se.rho == 0.0) continue;
    scatter18(asm_.dofs, se, se.el.me(se.rho, lumped),
              [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
  }
  for (const ContinuumKernel& k : asm_.continua) {
    if (k.M.squaredNorm() == 0) continue;
    if (!lumped) {
      scatter_continuum(asm_.dofs, k, k.M, trips);
      continue;
    }
    // HRZ lumping: the diagonal scaled to the element mass per direction
    MatX L = MatX::Zero(k.M.rows(), k.M.cols());
    for (int a = 0; a < k.dim; ++a) {
      double tot = 0, dg = 0;
      for (int i = a; i < k.M.rows(); i += k.dim) {
        dg += k.M(i, i);
        for (int j = a; j < k.M.cols(); j += k.dim) tot += k.M(i, j);
      }
      for (int i = a; i < k.M.rows(); i += k.dim) L(i, i) = k.M(i, i) * tot / dg;
    }
    scatter_continuum(asm_.dofs, k, L, trips);
  }
  SpMat M(asm_.dofs.n_free, asm_.dofs.n_free);
  M.setFromTriplets(trips.begin(), trips.end());
  return M;
}

SpMat assemble_geometric(const Model& model, const Assembled& asm_, const std::vector<double>& N,
                         const std::vector<Vec3>* shell_n) {
  require(N.size() == asm_.elements.size(), "assemble_geometric: axial field size mismatch");
  std::vector<Triplet> trips;
  for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) continue;
    if (N[e] == 0.0) continue;
    scatter12(model, asm_.dofs, b, asm_.elements[e].kg(N[e]),
              [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
  }
  if (shell_n) {
    require(shell_n->size() == asm_.shells.size(),
            "assemble_geometric: shell membrane field size mismatch");
    for (const ShellElem& se : asm_.shells) {
      const Vec3& n = (*shell_n)[se.index];
      if (n.squaredNorm() == 0.0) continue;
      scatter18(asm_.dofs, se, se.el.kg(n(0), n(1), n(2)),
                [&](int i, int j, double v) { trips.emplace_back(i, j, v); });
    }
  }
  SpMat G(asm_.dofs.n_free, asm_.dofs.n_free);
  G.setFromTriplets(trips.begin(), trips.end());
  return G;
}

// ---------------- load lowering ----------------

CaseLoads empty_case_loads(const Assembled& asm_) {
  const DofMap& d = asm_.dofs;
  CaseLoads cl;
  cl.rhs = VecX::Zero(d.n_free);
  cl.nodal_full = VecX::Zero(kDofPerNode * d.n_nodes);
  cl.f0.assign(asm_.elements.size(), Vec12::Zero());
  cl.span.assign(asm_.elements.size(), BarSpanLoads{});
  cl.shell_recovery.assign(asm_.shells.size(), ShellRecovery{});
  cl.shell_f.assign(asm_.shells.size(), Vec18::Zero());
  return cl;
}

namespace {
// A bar load lowered by the default branch of lower_load: its fixed-end
// forces f0 (local) and their global nodal equivalent, both pure functions of
// the load — computable ahead, in parallel, then applied in order.
bool plain_bar_load(const Load& l) {
  switch (l.kind) {
    case Load::Kind::NodalForce:
    case Load::Kind::SupportDisp:
    case Load::Kind::SelfWeight:
    case Load::Kind::SurfacePressure:
    case Load::Kind::SurfaceThermal:
    case Load::Kind::PolygonPressure:
    case Load::Kind::Hydrostatic:
    case Load::Kind::ShellPrestress:
    case Load::Kind::ShellImposed: return false;
    default: return true;
  }
}
struct PreLowered {
  Vec12 f0, eq;
};
void lower_load_impl(const Model& model, const Assembled& asm_, const std::string& case_name,
                     const Load& l, CaseLoads& cl, const PreLowered* pre);
}  // namespace

CaseLoads lower_case(const Model& model, const Assembled& asm_, const LoadCase& lc) {
  CaseLoads cl = empty_case_loads(asm_);
  for (const Load& l : lc.loads) lower_load(model, asm_, lc.name, l, cl);
  return cl;
}

CaseLoads lower_case(const Model& model, const Assembled& asm_, const LoadCase& lc, int threads) {
  if (threads == 1 || lc.loads.size() < 64) return lower_case(model, asm_, lc);
  CaseLoads cl = empty_case_loads(asm_);
  if (!model.bars.empty()) (void)model.bar_index(model.bars.front().id);  // index built before the workers
  std::vector<PreLowered> pre(lc.loads.size());
  std::vector<char> has(lc.loads.size(), 0);
  parallel_for(lc.loads.size(), threads, [&](std::size_t i) {
    const Load& l = lc.loads[i];
    if (!plain_bar_load(l)) return;
    const std::size_t e = static_cast<std::size_t>(model.bar_index(l.bar));
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) return;  // the sequential pass raises its error
    // prismatic bars only (raw_f0's path for them); a CARA VAR / bedded
    // bar keeps the sequential lowering
    if (model.tapered_of(b.id)) return;
    pre[i].f0 = asm_.elements[e].f0_local(l, model.material(b.material));
    pre[i].eq = asm_.elements[e].T().transpose() * asm_.elements[e].condense_f0(pre[i].f0);
    has[i] = 1;
  });
  for (std::size_t i = 0; i < lc.loads.size(); ++i)
    lower_load_impl(model, asm_, lc.name, lc.loads[i], cl, has[i] ? &pre[i] : nullptr);
  return cl;
}

void lower_load(const Model& model, const Assembled& asm_, const std::string& case_name,
                const Load& l, CaseLoads& cl) {
  lower_load_impl(model, asm_, case_name, l, cl, nullptr);
}

namespace {
void lower_load_impl(const Model& model, const Assembled& asm_, const std::string& case_name,
                     const Load& l, CaseLoads& cl, const PreLowered* pre) {
  const DofMap& d = asm_.dofs;

  // Scatter a facet shell's 18 consistent nodal forces (element DOF order).
  auto scatter_shell = [&](const ShellElem& se, const Vec18& f18) {
    int gdof[18];
    shell_gdofs(d, se, gdof);
    if (se.index < cl.shell_f.size()) cl.shell_f[se.index] += f18;
    for (int a = 0; a < 18; ++a) {
      cl.nodal_full(gdof[a]) += f18(a);
      int fi = d.free_ix[static_cast<std::size_t>(gdof[a])];
      if (fi >= 0) cl.rhs(fi) += f18(a);
    }
  };

  auto scatter_equiv = [&](std::size_t e, const Vec12& f_nodal) {
    const Bar& b = model.bars[e];
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    for (int c = 0; c < 6; ++c) {
      int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
      if (fi >= 0) cl.rhs(fi) += f_nodal(c);
      int fj = d.free_ix[static_cast<std::size_t>(d.full(nj, c))];
      if (fj >= 0) cl.rhs(fj) += f_nodal(6 + c);
      cl.nodal_full(d.full(ni, c)) += f_nodal(c);
      cl.nodal_full(d.full(nj, c)) += f_nodal(6 + c);
    }
  };

  auto record_into = [&](BarSpanLoads& sp, std::size_t e, const Load& l, const BarElement& el) {
    const Bar& b = model.bars[e];
    const Material& mat = model.material(b.material);
    const SectionProps& sec = model.section(b.section).props;
    double L = el.L();
    auto to_local = [&](const Vec3& v, LoadFrame fr) {
      return fr == LoadFrame::Local ? v : Vec3(el.R() * v);
    };
    auto proj_scale = [&](LoadFrame fr) -> Vec3 {
      if (fr != LoadFrame::GlobalProjected) return Vec3::Ones();
      Vec3 ex = el.R().row(0);
      return Vec3(std::sqrt(std::max(0.0, 1 - ex.x() * ex.x())),
                  std::sqrt(std::max(0.0, 1 - ex.y() * ex.y())),
                  std::sqrt(std::max(0.0, 1 - ex.z() * ex.z())));
    };
    switch (l.kind) {
      case Load::Kind::BarUdl:
      case Load::Kind::BarLinear: {
        double a = l.relative ? l.a * L : l.a;
        double bb = l.b < 0 ? L : (l.relative ? l.b * L : l.b);
        Vec3 s = proj_scale(l.frame);
        Vec3 qa = to_local(l.F.cwiseProduct(s), l.frame);
        Vec3 qb = to_local((l.kind == Load::Kind::BarUdl ? l.F : l.F2).cwiseProduct(s), l.frame);
        DistSeg sg{std::clamp(a, 0.0, L), std::clamp(bb, 0.0, L), qa, qb};
        sg.ma = to_local(l.M, l.frame);
        sg.mb = to_local(l.kind == Load::Kind::BarUdl ? l.M : l.M2, l.frame);
        sp.segs.push_back(sg);
        break;
      }
      case Load::Kind::BarParabolic: {
        double a = l.relative ? l.a * L : l.a;
        double bb = l.b < 0 ? L : (l.relative ? l.b * L : l.b);
        Vec3 s = proj_scale(l.frame);
        DistSeg sg{std::clamp(a, 0.0, L), std::clamp(bb, 0.0, L), to_local(l.F.cwiseProduct(s), l.frame),
                   to_local(l.F2.cwiseProduct(s), l.frame)};
        sg.qc = to_local(l.Fmid.cwiseProduct(s), l.frame) - 0.5 * (sg.qa + sg.qb);
        sp.segs.push_back(sg);
        break;
      }
      case Load::Kind::BarPoint: {
        double a = l.relative ? l.a * L : l.a;
        sp.points.push_back({std::clamp(a, 0.0, L), to_local(l.F, l.frame), to_local(l.M, l.frame)});
        break;
      }
      case Load::Kind::BarThermal:
        sp.eps0 += mat.alpha * l.dt;
        sp.kappa_z += mat.alpha * l.gy;  // gradient along y bends about z
        sp.kappa_y -= mat.alpha * l.gz;  // theta_y = -dw/dx (see element.cpp)
        break;
      case Load::Kind::Prestrain:
        sp.eps0 += l.strain;
        sp.kappa_y += l.kappa_y;
        sp.kappa_z += l.kappa_z;
        sp.kappa_x += l.kappa_x;
        break;
      case Load::Kind::PrestrainPoint: {
        const Load d = l.frame == LoadFrame::Global ? dislocation_local(l, el.R()) : l;
        double a = d.relative ? d.a * L : d.a;
        sp.kinks.push_back({std::clamp(a, 0.0, L), d.strain, d.kappa_y, d.kappa_z, d.dvy, d.dvz,
                            d.dthx});
        break;
      }
      case Load::Kind::PrestrainLinear:
        sp.eps0 += l.strain;
        sp.kappa_y += l.kappa_y;
        sp.kappa_z += l.kappa_z;
        sp.kappa_x += l.kappa_x;
        if (L > 0) {
          sp.deps0 += (l.strain2 - l.strain) / L;
          sp.dkappa_y += (l.kappa_y2 - l.kappa_y) / L;
          sp.dkappa_z += (l.kappa_z2 - l.kappa_z) / L;
          sp.dkappa_x += (l.kappa_x2 - l.kappa_x) / L;
        }
        break;
      case Load::Kind::PrestrainParabolic: {
        // f(x) = f_i + (f_j - f_i + 4c) x/L - 4c x^2/L^2, c = f_m - (f_i + f_j)/2
        require(L > 0, "parabolic prestrain on a zero-length bar " + std::to_string(l.bar));
        auto add = [&](double fi, double fm, double fj, double& v0, double& d1, double& d2) {
          const double c = fm - 0.5 * (fi + fj);
          v0 += fi;
          d1 += (fj - fi + 4.0 * c) / L;
          d2 += -4.0 * c / (L * L);
        };
        add(l.strain, l.strain_m, l.strain2, sp.eps0, sp.deps0, sp.c2eps0);
        add(l.kappa_y, l.kappa_y_m, l.kappa_y2, sp.kappa_y, sp.dkappa_y, sp.c2kappa_y);
        add(l.kappa_z, l.kappa_z_m, l.kappa_z2, sp.kappa_z, sp.dkappa_z, sp.c2kappa_z);
        add(l.kappa_x, l.kappa_x_m, l.kappa_x2, sp.kappa_x, sp.dkappa_x, sp.c2kappa_x);
        break;
      }
      case Load::Kind::SelfWeight: {
        // A variable section (CARA VAR) weighs rho A(x): 64 linear pieces
        // (exact for a linear A, below any printed digit for a parabolic
        // one), as its equivalent loads (tapered f0) take it.
        const TaperedBar* tb = model.tapered_of(el.bar().id);
        if (tb && tb->props) {
          constexpr int kPieces = 64;
          for (int k = 0; k < kPieces; ++k) {
            const double u0 = static_cast<double>(k) / kPieces, u1 = static_cast<double>(k + 1) / kPieces;
            Vec3 qa = el.R() * (mat.rho * tb->props(u0).A * l.g);
            Vec3 qb = el.R() * (mat.rho * tb->props(u1).A * l.g);
            sp.segs.push_back({u0 * L, u1 * L, qa, qb});
          }
          break;
        }
        Vec3 ql = el.R() * (mat.rho * sec.A * l.g);
        sp.segs.push_back({0.0, L, ql, ql});
        break;
      }
      default: break;
    }
  };
  auto record_span = [&](std::size_t e, const Load& l, const BarElement& el) {
    record_into(cl.span[e], e, l, el);
  };

  // ST1 CARA VAR: a tapered bar lowers its span loads through the exact
  // flexibility integrals of `tapered.h` instead of the prismatic shape
  // functions (see `TaperedBar` in model.h).
  auto raw_f0 = [&](std::size_t e, const Load& l, const Material& mat) -> Vec12 {
    const TaperedBar* tb = model.tapered_of(model.bars[e].id);
    if (tb && tb->f0_span) {  // custom block (bef.h): from the load's span record
      BarSpanLoads one;
      record_into(one, e, l, asm_.elements[e]);
      return tb->f0_span(one);
    }
    if (tb) return tb->f0(l);
    return asm_.elements[e].f0_local(l, mat);
  };


  {
    switch (l.kind) {
      case Load::Kind::NodalForce: {
        int ni = model.node_index(l.node);
        for (int c = 0; c < 3; ++c) {
          cl.nodal_full(d.full(ni, c)) += l.F(c);
          cl.nodal_full(d.full(ni, 3 + c)) += l.M(c);
          int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
          if (fi >= 0) cl.rhs(fi) += l.F(c);
          int fr = d.free_ix[static_cast<std::size_t>(d.full(ni, 3 + c))];
          if (fr >= 0) cl.rhs(fr) += l.M(c);
        }
        break;
      }
      case Load::Kind::SupportDisp: {
        int ni = model.node_index(l.node);
        for (int c = 0; c < 6; ++c) {
          if (!l.disp[static_cast<std::size_t>(c)]) continue;
          int f = d.full(ni, c);
          require(d.blocked[static_cast<std::size_t>(f)],
                  "case '" + case_name + "': prescribed displacement at node " +
                      std::to_string(l.node) + " DOF " + std::to_string(c) +
                      " which is not rigidly blocked by a support");
          for (const auto& [pf, pv] : cl.prescribed)
            require(pf != f, "case '" + case_name + "': DOF prescribed twice");
          cl.prescribed.emplace_back(f, *l.disp[static_cast<std::size_t>(c)]);
        }
        break;
      }
      case Load::Kind::SelfWeight: {
        for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
          const Bar& b = model.bars[e];
          if (b.kind == BarKind::Cable) continue;
          const Material& mat = model.material(b.material);
          Vec12 f0 = raw_f0(e, l, mat);
          cl.f0[e] += f0;
          scatter_equiv(e, asm_.elements[e].T().transpose() * asm_.elements[e].condense_f0(f0));
          record_span(e, l, asm_.elements[e]);
        }
        // Facet shells: rho*t*g force/area, tributary A/3 per node.
        for (const ShellElem& se : asm_.shells)
          scatter_shell(se, se.el.self_weight_load(se.rho, l.g));
        // Continuum elements: the consistent body force int rho N g = M (g per node)
        for (const ContinuumKernel& k : asm_.continua) {
          const int n = static_cast<int>(k.nodes.size());
          VecX gv(k.dim * n);
          for (int i = 0; i < n; ++i)
            for (int a = 0; a < k.dim; ++a) gv(k.dim * i + a) = l.g(a);
          const VecX f = k.M * gv;
          for (int i = 0; i < n; ++i)
            for (int a = 0; a < k.dim; ++a) {
              const int g = d.full(k.nodes[static_cast<std::size_t>(i)], a);
              cl.nodal_full(g) += f(k.dim * i + a);
              const int fi = d.free_ix[static_cast<std::size_t>(g)];
              if (fi >= 0) cl.rhs(fi) += f(k.dim * i + a);
            }
        }
        for (const PointMass& pm : model.point_masses) {
          int ni = model.node_index(pm.node);
          Vec3 F = pm.mass * l.g;
          for (int c = 0; c < 3; ++c) {
            cl.nodal_full(d.full(ni, c)) += F(c);
            int fi = d.free_ix[static_cast<std::size_t>(d.full(ni, c))];
            if (fi >= 0) cl.rhs(fi) += F(c);
          }
        }
        break;
      }
      case Load::Kind::SurfacePressure: {
        // Applied to every ShellTri carrying this slab id.
        for (const ShellElem& se : asm_.shells) {
          if (se.slab != l.slab) continue;
          switch (l.sframe) {
            case SurfaceFrame::Global: scatter_shell(se, se.el.pressure_load_global(l.F)); break;
            case SurfaceFrame::LocalNormal:
              scatter_shell(se, se.el.pressure_load_normal(l.F.z()));
              break;
            case SurfaceFrame::Projected:
              scatter_shell(se, se.el.pressure_load_projected(l.F));
              break;
          }
        }
        break;
      }
      case Load::Kind::Hydrostatic: {
        for (const ShellElem& se : asm_.shells)
          if (l.slab == 0 || se.slab == l.slab)
            scatter_shell(se, se.el.hydrostatic_load(l.hydro_gamma, l.hydro_point, l.g, l.hydro_side));
        break;
      }
      case Load::Kind::ShellPrestress: {
        const ShellTendon* tn = nullptr;
        for (const ShellTendon& t : model.shell_tendons)
          if (t.id == l.node) tn = &t;
        require(tn != nullptr, "shell prestress: no tendon " + std::to_string(l.node));
        for (const auto& [k, f] : shell_tendon_loads(model, asm_.shells, *tn, l.strain)) scatter_shell(asm_.shells[k], f);
        break;
      }
      case Load::Kind::ShellImposed: {
        Vec18 d;
        for (int k = 0; k < 6; ++k) d.segment<3>(3 * k) = l.poly[static_cast<std::size_t>(k)];
        bool found = false;
        for (const ShellElem& se : asm_.shells)
          if (se.id == l.node) {
            scatter_shell(se, se.el.ke() * d);
            cl.shell_recovery[se.index].dstar += d;
            found = true;
          }
        require(found, "imposed shell field: no shell " + std::to_string(l.node));
        break;
      }
      case Load::Kind::PolygonPressure: {
        // the plane normal to the projection, an orthonormal (e1, e2) in it
        const Vec3 n = l.proj.normalized();
        const Vec3 e1 = (std::fabs(n.x()) < 0.9 ? Vec3::UnitX() : Vec3::UnitY()).cross(n).normalized();
        const Vec3 e2 = n.cross(e1);
        using P2 = Eigen::Vector2d;
        auto to2 = [&](const Vec3& p) { return P2(p.dot(e1), p.dot(e2)); };
        auto ring2 = [&](const std::vector<Vec3>& r) {
          std::vector<P2> o;
          for (const Vec3& p : r) o.push_back(to2(p));
          return o;
        };
        const std::vector<P2> outer = ring2(l.poly);
        std::vector<std::vector<P2>> holes;
        for (const auto& h : l.holes) holes.push_back(ring2(h));
        // area and first moments of (subject polygon) clipped by a convex
        // CCW triangle (Sutherland-Hodgman: exact area integrals even for a
        // non-convex subject)
        auto clip_moments = [](const std::vector<P2>& subj, const std::array<P2, 3>& tri, double& A, P2& S) {
          std::vector<P2> poly = subj;
          for (int k = 0; k < 3 && !poly.empty(); ++k) {
            const P2 a = tri[static_cast<std::size_t>(k)], b = tri[static_cast<std::size_t>((k + 1) % 3)];
            auto side = [&](const P2& p) { return (b.x() - a.x()) * (p.y() - a.y()) - (b.y() - a.y()) * (p.x() - a.x()); };
            std::vector<P2> out;
            for (std::size_t i = 0; i < poly.size(); ++i) {
              const P2 &p = poly[i], &q = poly[(i + 1) % poly.size()];
              const double sp = side(p), sq = side(q);
              if (sp >= 0) out.push_back(p);
              if ((sp >= 0) != (sq >= 0)) out.push_back(p + (q - p) * (sp / (sp - sq)));
            }
            poly.swap(out);
          }
          A = 0;
          S = P2::Zero();
          for (std::size_t i = 0; i < poly.size(); ++i) {
            const P2 &p = poly[i], &q = poly[(i + 1) % poly.size()];
            const double c = p.x() * q.y() - q.x() * p.y();
            A += c / 2;
            S += (p + q) * c / 6;
          }
        };
        auto signed_area = [](const std::vector<P2>& r) {
          double A = 0;
          for (std::size_t i = 0; i < r.size(); ++i)
            A += r[i].x() * r[(i + 1) % r.size()].y() - r[(i + 1) % r.size()].x() * r[i].y();
          return A / 2;
        };
        // orientation-independent: every ring taken CCW
        auto ccw = [&](std::vector<P2> r) {
          if (signed_area(r) < 0) std::reverse(r.begin(), r.end());
          return r;
        };
        const std::vector<P2> outer_ccw = ccw(outer);
        std::vector<std::vector<P2>> holes_ccw;
        for (const auto& h : holes) holes_ccw.push_back(ccw(h));
        for (const ShellElem& se : asm_.shells) {
          if (l.slab != 0 && se.slab != l.slab) continue;
          std::array<P2, 3> tri;
          for (int k = 0; k < 3; ++k)
            tri[static_cast<std::size_t>(k)] = to2(model.nodes[static_cast<std::size_t>(se.nodes[static_cast<std::size_t>(k)])].x);
          double At = (tri[1].x() - tri[0].x()) * (tri[2].y() - tri[0].y()) - (tri[2].x() - tri[0].x()) * (tri[1].y() - tri[0].y());
          if (std::fabs(At) <= 1e-14 * (tri[1] - tri[0]).squaredNorm()) continue;  // edge-on to the projection
          std::array<P2, 3> t = tri;
          if (At < 0) std::swap(t[1], t[2]);
          double A, Ah;
          P2 S, Sh;
          clip_moments(outer_ccw, t, A, S);
          for (const auto& h : holes_ccw) {
            clip_moments(h, t, Ah, Sh);
            A -= Ah;
            S -= Sh;
          }
          if (A <= 0) continue;
          const P2 c = S / A;
          // barycentric weights of the overlap's centroid in the facet
          // (node order of the facet, not of the swapped t)
          const double det = (tri[1].x() - tri[0].x()) * (tri[2].y() - tri[0].y()) -
                             (tri[2].x() - tri[0].x()) * (tri[1].y() - tri[0].y());
          const double w1 = ((c.x() - tri[0].x()) * (tri[2].y() - tri[0].y()) -
                             (tri[2].x() - tri[0].x()) * (c.y() - tri[0].y())) / det;
          const double w2 = ((tri[1].x() - tri[0].x()) * (c.y() - tri[0].y()) -
                             (c.x() - tri[0].x()) * (tri[1].y() - tri[0].y())) / det;
          const double w[3] = {1 - w1 - w2, w1, w2};
          Vec18 f = Vec18::Zero();
          for (int k = 0; k < 3; ++k) f.segment<3>(6 * k) = l.F * (A * w[k]);
          scatter_shell(se, f);
        }
        break;
      }
      case Load::Kind::SurfaceThermal: {
        // Initial-strain method: eps0 = alpha*dt*[1,1,0],
        // kappa0 = alpha*grad*[1,1,0]. The (eps0, kappa0) pair is ALSO recorded
        // per shell so recovery subtracts it (a free expansion must read zero).
        for (const ShellElem& se : asm_.shells) {
          if (se.slab != l.slab) continue;
          Vec3 eps0(se.alpha * l.dt, se.alpha * l.dt, 0.0);
          Vec3 kappa0(se.alpha * l.grad, se.alpha * l.grad, 0.0);
          scatter_shell(se, se.el.thermal_load(eps0, kappa0));
          cl.shell_recovery[se.index].eps0 += eps0;
          cl.shell_recovery[se.index].kappa0 += kappa0;
        }
        break;
      }
      default: {  // bar loads
        std::size_t e = static_cast<std::size_t>(model.bar_index(l.bar));
        const Bar& b = model.bars[e];
        require(b.kind != BarKind::Cable, "bar load on a cable");
        if (pre) {  // computed ahead (lower_case with threads): the same values
          cl.f0[e] += pre->f0;
          scatter_equiv(e, pre->eq);
          record_span(e, l, asm_.elements[e]);
          break;
        }
        const Material& mat = model.material(b.material);
        Vec12 f0 = raw_f0(e, l, mat);
        cl.f0[e] += f0;
        scatter_equiv(e, asm_.elements[e].T().transpose() * asm_.elements[e].condense_f0(f0));
        record_span(e, l, asm_.elements[e]);
        break;
      }
    }
  }
}
}  // namespace

void fold_link_offsets(const Model& model, const Assembled& asm_, const SpMat& K, CaseLoads& cl) {
  (void)model;
  const LinkReduction& red = asm_.red;
  cl.link_u0.resize(0);
  if (!red.active || cl.prescribed.empty()) return;
  std::map<int, double> pv(cl.prescribed.begin(), cl.prescribed.end());
  VecX u0 = VecX::Zero(asm_.dofs.n_free);
  bool any = false;
  for (const auto& e : red.slaves) {
    double v = 0;
    for (const auto& [k, c] : e.terms) {
      auto it = pv.find(k);
      if (it != pv.end()) v += c * it->second;
    }
    if (v != 0) {
      u0(asm_.dofs.free_ix[static_cast<std::size_t>(e.slave)]) = v;
      any = true;
    }
  }
  if (!any) return;
  cl.rhs -= K * u0;
  cl.link_u0 = std::move(u0);
}

void fold_prescribed(const Model& model, const Assembled& asm_, CaseLoads& cl) {
  if (cl.prescribed.empty()) return;
  const DofMap& d = asm_.dofs;
  // Walk each element touching a prescribed DOF and subtract K[free, p]*u_p.
  for (std::size_t e = 0; e < asm_.elements.size(); ++e) {
    const Bar& b = model.bars[e];
    if (b.kind == BarKind::Cable) continue;
    int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
    int gdof[12];
    for (int c = 0; c < 6; ++c) {
      gdof[c] = d.full(ni, c);
      gdof[6 + c] = d.full(nj, c);
    }
    Mat12 k = bar_ke_assembled(model, asm_.elements[e]);
    for (const auto& [pf, pv] : cl.prescribed) {
      for (int a = 0; a < 12; ++a) {
        if (gdof[a] != pf) continue;
        for (int r = 0; r < 12; ++r) {
          int fr = d.free_ix[static_cast<std::size_t>(gdof[r])];
          if (fr >= 0) cl.rhs(fr) -= k(r, a) * pv;
        }
      }
    }
  }
  // Same element-implicit fold for the facet shells.
  for (const ShellElem& se : asm_.shells) {
    int gdof[18];
    shell_gdofs(d, se, gdof);
    const Mat18& k = se.el.ke();
    for (const auto& [pf, pv] : cl.prescribed) {
      for (int a = 0; a < 18; ++a) {
        if (gdof[a] != pf) continue;
        for (int r = 0; r < 18; ++r) {
          int fr = d.free_ix[static_cast<std::size_t>(gdof[r])];
          if (fr >= 0) cl.rhs(fr) -= k(r, a) * pv;
        }
      }
    }
  }
  // And for the continuum elements (their translations).
  for (const ContinuumKernel& k : asm_.continua) {
    const int n = static_cast<int>(k.nodes.size()), nd = k.dim * n;
    std::vector<int> gdof(static_cast<std::size_t>(nd));
    for (int i = 0; i < n; ++i)
      for (int a = 0; a < k.dim; ++a) gdof[static_cast<std::size_t>(k.dim * i + a)] = d.full(k.nodes[static_cast<std::size_t>(i)], a);
    for (const auto& [pf, pv] : cl.prescribed)
      for (int a = 0; a < nd; ++a) {
        if (gdof[static_cast<std::size_t>(a)] != pf) continue;
        for (int r = 0; r < nd; ++r) {
          const int fr = d.free_ix[static_cast<std::size_t>(gdof[static_cast<std::size_t>(r)])];
          if (fr >= 0) cl.rhs(fr) -= k.K(r, a) * pv;
        }
      }
  }
}

bool has_parabolic(const std::vector<DistSeg>& segs) {
  return std::any_of(segs.begin(), segs.end(), [](const DistSeg& g) { return !g.qc.isZero(); });
}

void require_linear_span(const std::vector<DistSeg>& segs, const std::string& who) {
  require(!has_parabolic(segs), who + ": parabolic bar loads are not supported here (linear span loads only)");
}

void require_no_imposed_twist(const BarSpanLoads& sp, const std::string& who) {
  require(sp.kappa_x == 0 && sp.dkappa_x == 0 && sp.c2kappa_x == 0,
          who + ": an imposed twist rate (kappa_x) is not supported here");
}

bool has_parabolic_strain(const BarSpanLoads& sp) {
  return sp.c2eps0 != 0 || sp.c2kappa_y != 0 || sp.c2kappa_z != 0 || sp.c2kappa_x != 0;
}

void require_linear_strain(const BarSpanLoads& sp, const std::string& who) {
  require(!has_parabolic_strain(sp),
          who + ": a parabolic imposed strain is not supported here (linear imposed fields only)");
}

}  // namespace keta

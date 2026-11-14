#include "keta/substructure.h"

#include <map>
#include <set>

#include "keta/assembly.h"

namespace keta {

namespace {
// P^t A Q for index lists (dense)
MatX block(const SpMat& A, const std::vector<int>& rows, const std::vector<int>& cols) {
  std::vector<int> cpos(static_cast<std::size_t>(A.cols()), -1), rpos(static_cast<std::size_t>(A.rows()), -1);
  for (std::size_t j = 0; j < cols.size(); ++j) cpos[static_cast<std::size_t>(cols[j])] = static_cast<int>(j);
  for (std::size_t i = 0; i < rows.size(); ++i) rpos[static_cast<std::size_t>(rows[i])] = static_cast<int>(i);
  MatX out = MatX::Zero(static_cast<Eigen::Index>(rows.size()), static_cast<Eigen::Index>(cols.size()));
  for (int k = 0; k < A.outerSize(); ++k)
    for (SpMat::InnerIterator it(A, k); it; ++it) {
      const int r = rpos[static_cast<std::size_t>(it.row())], c = cpos[static_cast<std::size_t>(it.col())];
      if (r >= 0 && c >= 0) out(r, c) += it.value();
    }
  return out;
}
SpMat sparse_block(const SpMat& A, const std::vector<int>& idx) {
  std::vector<int> pos(static_cast<std::size_t>(A.rows()), -1);
  for (std::size_t i = 0; i < idx.size(); ++i) pos[static_cast<std::size_t>(idx[i])] = static_cast<int>(i);
  std::vector<Triplet> t;
  for (int k = 0; k < A.outerSize(); ++k)
    for (SpMat::InnerIterator it(A, k); it; ++it) {
      const int r = pos[static_cast<std::size_t>(it.row())], c = pos[static_cast<std::size_t>(it.col())];
      if (r >= 0 && c >= 0) t.emplace_back(r, c, it.value());
    }
  SpMat out(static_cast<Eigen::Index>(idx.size()), static_cast<Eigen::Index>(idx.size()));
  out.setFromTriplets(t.begin(), t.end());
  return out;
}
}  // namespace

SubstructuredModal solve_modal_substructured(const Model& model, const std::vector<Substructure>& subs,
                                             const std::vector<Id>& masters, const ModalOptions& opt) {
  model.validate();
  require(!subs.empty(), "substructuring: no substructure");
  require(opt.axial_state == nullptr, "substructuring: a prestressed modal analysis is not substructured");
  Assembled asm_ = assemble_elements(model, true);
  require(!asm_.red.active,
          "substructuring: rigid links, constraints and skew supports (master-slave reductions) are not "
          "substructured — model them in the residual structure of a full analysis");
  const DofMap& d = asm_.dofs;
  require(d.n_free > 0, "substructuring: no free DOFs");
  const SpMat K = assemble_stiffness(model, asm_);
  const SpMat M = assemble_mass(model, asm_, opt.lumped_mass);
  require(M.nonZeros() > 0, "substructuring: model has no mass");

  // owner of each node: -1 untouched, s >= 0 one substructure, -2 shared /
  // residual / master
  const int nn = static_cast<int>(model.nodes.size());
  std::vector<int> owner(static_cast<std::size_t>(nn), -1);
  auto touch = [&](Id node, int s) {
    int& o = owner[static_cast<std::size_t>(model.node_index(node))];
    o = (o == -1 || o == s) ? s : -2;
  };
  std::vector<int> bar_sub(model.bars.size(), -1), shell_sub(model.shells.size(), -1);
  std::map<Id, std::size_t> bar_ix, shell_ix;
  for (std::size_t i = 0; i < model.bars.size(); ++i) bar_ix[model.bars[i].id] = i;
  for (std::size_t i = 0; i < model.shells.size(); ++i) shell_ix[model.shells[i].id] = i;
  for (std::size_t s = 0; s < subs.size(); ++s) {
    require(subs[s].n_modes >= 0, "substructuring: n_modes must be >= 0");
    for (Id b : subs[s].bars) {
      auto it = bar_ix.find(b);
      require(it != bar_ix.end(), "substructuring: unknown bar " + std::to_string(b));
      require(bar_sub[it->second] == -1, "substructuring: bar " + std::to_string(b) + " in two substructures");
      bar_sub[it->second] = static_cast<int>(s);
    }
    for (Id e : subs[s].shells) {
      auto it = shell_ix.find(e);
      require(it != shell_ix.end(), "substructuring: unknown shell " + std::to_string(e));
      require(shell_sub[it->second] == -1, "substructuring: shell " + std::to_string(e) + " in two substructures");
      shell_sub[it->second] = static_cast<int>(s);
    }
  }
  for (std::size_t i = 0; i < model.bars.size(); ++i) {
    touch(model.bars[i].ni, bar_sub[i] < 0 ? -2 : bar_sub[i]);
    touch(model.bars[i].nj, bar_sub[i] < 0 ? -2 : bar_sub[i]);
  }
  for (std::size_t i = 0; i < model.shells.size(); ++i)
    for (Id n : model.shells[i].n) touch(n, shell_sub[i] < 0 ? -2 : shell_sub[i]);
  for (const Continuum& c : model.continua)  // continua stay in the residual structure
    for (Id n : c.n) touch(n, -2);
  for (const ElasticLink& l : model.links) {
    touch(l.ni, -2);
    touch(l.nj, -2);
  }
  for (Id n : masters) touch(n, -2);

  std::vector<std::vector<int>> interior(subs.size());
  std::vector<int> boundary;
  std::vector<char> is_interior(static_cast<std::size_t>(d.n_free), 0);
  for (int i = 0; i < nn; ++i)
    for (int c = 0; c < kDofPerNode; ++c) {
      const int f = d.free_of(i, c);
      if (f < 0) continue;
      const int o = owner[static_cast<std::size_t>(i)];
      if (o >= 0) {
        interior[static_cast<std::size_t>(o)].push_back(f);
        is_interior[static_cast<std::size_t>(f)] = 1;
      }
    }
  for (int f = 0; f < d.n_free; ++f)
    if (!is_interior[static_cast<std::size_t>(f)]) boundary.push_back(f);
  for (std::size_t s = 0; s < subs.size(); ++s) std::sort(interior[s].begin(), interior[s].end());

  SubstructuredModal out;
  out.n_boundary = static_cast<int>(boundary.size());
  out.boundary_free = boundary;
  int n_red = out.n_boundary;
  std::vector<int> q0(subs.size());
  for (std::size_t s = 0; s < subs.size(); ++s) {
    require(!interior[s].empty(), "substructuring: substructure " + std::to_string(s) + " has no interior DOF");
    require(subs[s].n_modes <= static_cast<int>(interior[s].size()),
            "substructuring: more fixed-interface modes than interior DOFs");
    q0[s] = n_red;
    n_red += subs[s].n_modes;
  }
  out.n_reduced = n_red;

  std::vector<Triplet> tt;
  for (std::size_t j = 0; j < boundary.size(); ++j) tt.emplace_back(boundary[j], static_cast<int>(j), 1.0);
  for (std::size_t s = 0; s < subs.size(); ++s) {
    const std::vector<int>& I = interior[s];
    SubstructureInfo info;
    info.n_interior = static_cast<int>(I.size());
    const SpMat Kii = sparse_block(K, I), Mii = sparse_block(M, I);
    // the boundary DOFs this substructure couples to
    std::set<int> bset;
    std::vector<char> inI(static_cast<std::size_t>(d.n_free), 0);
    for (int f : I) inI[static_cast<std::size_t>(f)] = 1;
    for (int k = 0; k < K.outerSize(); ++k)
      for (SpMat::InnerIterator it(K, k); it; ++it)
        if (inI[static_cast<std::size_t>(it.row())] && !inI[static_cast<std::size_t>(it.col())])
          bset.insert(static_cast<int>(it.col()));
    const std::vector<int> bl(bset.begin(), bset.end());
    FactorCache fac;
    fac.factorize(Kii);
    fac.require_spd();
    if (!bl.empty()) {
      const MatX Psi = -fac.solve(MatX(block(K, I, bl)));
      std::vector<int> bpos(static_cast<std::size_t>(d.n_free), -1);
      for (std::size_t j = 0; j < boundary.size(); ++j) bpos[static_cast<std::size_t>(boundary[j])] = static_cast<int>(j);
      for (std::size_t j = 0; j < bl.size(); ++j)
        for (std::size_t i = 0; i < I.size(); ++i) {
          const double v = Psi(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
          if (v != 0) tt.emplace_back(I[i], bpos[static_cast<std::size_t>(bl[j])], v);
        }
    }
    if (subs[s].n_modes > 0) {
      std::vector<PencilPair> pairs = static_cast<int>(I.size()) <= opt.dense_threshold
                                          ? dense_pencil(Kii, Mii)
                                          : k_inner_lanczos(fac, Kii, Mii, subs[s].n_modes);
      int kept = 0;
      for (const PencilPair& p : pairs) {
        if (kept == subs[s].n_modes) break;
        require(p.mu > 0, "substructuring: a fixed-interface mode without mass");
        VecX phi = p.x / std::sqrt(p.x.dot(Mii * p.x));
        for (std::size_t i = 0; i < I.size(); ++i)
          if (phi(static_cast<Eigen::Index>(i)) != 0) tt.emplace_back(I[i], q0[s] + kept, phi(static_cast<Eigen::Index>(i)));
        info.fixed_interface_hz.push_back(1 / std::sqrt(p.mu) / (2 * kPi));
        ++kept;
      }
      require(kept == subs[s].n_modes, "substructuring: fixed-interface modes not converged");
    }
    out.subs.push_back(std::move(info));
  }
  out.T.resize(d.n_free, n_red);
  out.T.setFromTriplets(tt.begin(), tt.end());
  out.K = MatX(out.T.transpose() * K * out.T);
  out.M = MatX(out.T.transpose() * M * out.T);
  out.K = 0.5 * (out.K + out.K.transpose()).eval();
  out.M = 0.5 * (out.M + out.M.transpose()).eval();
  std::vector<PencilPair> pairs = dense_pencil(out.K.sparseView(), out.M.sparseView());
  for (PencilPair& p : pairs) p.x = out.T * p.x;
  out.modal = modal_from_pairs(model, asm_, M, pairs, opt.n_modes, opt.support_coupling && !opt.lumped_mass);
  return out;
}

}  // namespace keta

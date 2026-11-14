#include "keta/eigsolve.h"

#include <algorithm>

#include <Eigen/Eigenvalues>
#include <Eigen/SparseCholesky>

#include "keta/springmass.h"
#include "keta/supernodal.h"

namespace keta {

double pencil_residual(const SpMat& K, const SpMat& B, double mu, const VecX& x) {
  VecX bx = B * x, kx = K * x;
  double num = (bx - mu * kx).norm();
  double den = bx.norm() + std::fabs(mu) * kx.norm();
  return den > 0 ? num / den : 0.0;
}

std::vector<PencilPair> dense_pencil(const SpMat& K, const SpMat& B) {
  Eigen::Index n = K.rows();
  MatX Kd = MatX(K), Bd = MatX(B);
  Eigen::LLT<MatX> llt(Kd);
  require(llt.info() == Eigen::Success, "eig: stiffness not positive definite");
  // Reduce B y = mu K y with K = L L^t: A = L^-1 B L^-t, A u = mu u, y = L^-t u.
  MatX Linv_B = llt.matrixL().solve(Bd);
  MatX A = llt.matrixL().solve(Linv_B.transpose()).transpose();
  A = 0.5 * (A + A.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<MatX> es(A);
  require(es.info() == Eigen::Success, "eig: dense eigensolve failed");
  std::vector<PencilPair> out;
  for (Eigen::Index i = n - 1; i >= 0; --i) {  // descending mu
    PencilPair p;
    p.mu = es.eigenvalues()(i);
    p.x = llt.matrixU().solve(es.eigenvectors().col(i));
    p.residual = pencil_residual(K, B, p.mu, p.x);
    out.push_back(std::move(p));
  }
  return out;
}

namespace {

// splitmix64 (Steele, Lea & Flood 2014): a reproducible stream for the
// pseudo-random start columns.
std::uint64_t splitmix64(std::uint64_t& state) {
  std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
  return z ^ (z >> 31);
}

// Start column `column` of a Lanczos block. Column 0 of seed 0 is the
// historical deterministic 1, 1/2, 1/3, ... (single-vector runs unchanged);
// every other column is uniform in [-1, 1), reproducible from (seed, column).
VecX start_column(Eigen::Index n, int column, std::uint64_t seed) {
  VecX v(n);
  if (column == 0 && seed == 0) {
    for (Eigen::Index i = 0; i < n; ++i) v(i) = 1.0 / static_cast<double>(i + 1);
    return v;
  }
  std::uint64_t state = (seed << 32) ^ static_cast<std::uint64_t>(column + 1);
  for (Eigen::Index i = 0; i < n; ++i)
    v(i) = static_cast<double>(splitmix64(state) >> 11) * 0x1.0p-52 - 1.0;
  return v;
}

// A K-orthonormal family with its K-images (a^T K b = a.dot(K b)).
struct KFamily {
  std::vector<VecX> x, kx;

  std::size_t size() const { return x.size(); }
  void push(VecX v, VecX kv) {
    x.push_back(std::move(v));
    kx.push_back(std::move(kv));
  }
  // One classical Gram-Schmidt pass of w against the family.
  void sweep(VecX& w) const {
    for (std::size_t i = 0; i < x.size(); ++i) {
      double c = w.dot(kx[i]);
      w -= c * x[i];
    }
  }
};

// Full reorthogonalization: two Gram-Schmidt passes against the deflated
// (locked) vectors and the Krylov basis.
void reorthogonalize(VecX& w, const KFamily& locked, const KFamily& basis) {
  for (int pass = 0; pass < 2; ++pass) {
    locked.sweep(w);
    basis.sweep(w);
  }
}

// The projected matrix T = Q^T B Q of the first `n_blocks` blocks: the
// block-tridiagonal Lanczos matrix with diagonal blocks A_j and sub-diagonal
// blocks R_{j+1} (the K-QR factors of the next block).
MatX lanczos_matrix(const std::vector<MatX>& A, const std::vector<MatX>& R, std::size_t n_blocks) {
  Eigen::Index m = 0;
  for (std::size_t j = 0; j < n_blocks; ++j) m += A[j].rows();
  MatX T = MatX::Zero(m, m);
  Eigen::Index at = 0;
  for (std::size_t j = 0; j < n_blocks; ++j) {
    const Eigen::Index b = A[j].rows();
    T.block(at, at, b, b) = A[j];
    if (j + 1 < n_blocks) {
      const MatX& r = R[j];  // rows: block j+1, cols: block j
      T.block(at + b, at, r.rows(), r.cols()) = r;
      T.block(at, at + b, r.cols(), r.rows()) = r.transpose();
    }
    at += b;
  }
  return T;
}

// Error bound of the Ritz pair y of T (the K-norm of A x - mu x for x = Q y):
// |R_next y_last|, y_last the rows of y on the last block.
double ritz_residual(const MatX& R_next, const VecX& y_last) {
  if (R_next.size() == 1) return std::fabs(R_next(0, 0) * y_last(0));
  return (R_next * y_last).norm();
}

}  // namespace

std::vector<PencilPair> k_inner_block_lanczos(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                              const BlockLanczosOptions& opt,
                                              const std::vector<VecX>& locked_x,
                                              const std::vector<VecX>& k_locked) {
  require(locked_x.size() == k_locked.size(), "lanczos: locked vectors without their K-images");
  require(opt.n_wanted > 0 && opt.block_size > 0, "lanczos: n_wanted and block_size must be positive");
  const Eigen::Index n = K.rows();
  KFamily locked{locked_x, k_locked};
  const int room = static_cast<int>(n) - static_cast<int>(locked.size());
  int m_max = opt.max_columns > 0 ? opt.max_columns
                                  : std::min(room, 8 * opt.n_wanted + 60 + (opt.block_size - 1));
  m_max = std::min(m_max, room);
  if (m_max <= 0) return {};

  // Start block: K-orthonormal, deflated against the locked vectors.
  KFamily block;
  for (int c = 0; c < std::min(opt.block_size, m_max); ++c) {
    VecX v = start_column(n, c, opt.seed);
    const double raw = v.norm();
    reorthogonalize(v, locked, block);
    if (c > 0 && v.norm() <= 1e-10 * raw) continue;  // dependent column
    VecX kv = K * v;
    double nrm = std::sqrt(std::max(0.0, v.dot(kv)));
    if (c == 0 && locked.size() == 0) require(nrm > 0, "lanczos: degenerate start vector");
    if (!(nrm > 0)) continue;
    v /= nrm;
    kv /= nrm;
    block.push(std::move(v), std::move(kv));
  }
  if (block.size() == 0) return {};

  KFamily basis;                // K-orthonormal Lanczos vectors (Q) and K Q
  std::vector<MatX> A, R;       // diagonal and sub-diagonal blocks of T
  double mu_scale = 0;          // largest |Rayleigh quotient| seen: the scale of the operator
  KFamily previous;             // the block before the current one

  while (!block.x.empty() && static_cast<int>(basis.size()) < m_max) {
    // Keep the basis within its budget: drop the trailing columns of the block.
    const std::size_t fits = static_cast<std::size_t>(m_max) - basis.size();
    if (block.size() > fits) {
      block.x.resize(fits);
      block.kx.resize(fits);
      R.back().conservativeResize(static_cast<Eigen::Index>(fits), R.back().cols());
    }
    const Eigen::Index b = static_cast<Eigen::Index>(block.size());
    for (Eigen::Index c = 0; c < b; ++c) basis.push(block.x[c], block.kx[c]);

    // W = K^-1 B V, A_j = V^T B V, then the three-term recurrence.
    std::vector<VecX> W(static_cast<std::size_t>(b));
    for (Eigen::Index c = 0; c < b; ++c) W[c] = Kfac.solve(VecX(B * block.x[c]));
    MatX Aj(b, b);
    for (Eigen::Index c = 0; c < b; ++c)
      for (Eigen::Index d = 0; d < b; ++d) Aj(c, d) = W[c].dot(block.kx[d]);
    Aj = (0.5 * (Aj + Aj.transpose())).eval();
    for (Eigen::Index c = 0; c < b; ++c) mu_scale = std::max(mu_scale, std::fabs(Aj(c, c)));
    for (Eigen::Index c = 0; c < b; ++c) {
      for (Eigen::Index d = 0; d < b; ++d) W[c] -= Aj(d, c) * block.x[d];
      if (!R.empty())  // V_{j-1} R_j^T: R_j has the rows of the current block
        for (Eigen::Index d = 0; d < R.back().cols(); ++d) W[c] -= R.back()(c, d) * previous.x[d];
      reorthogonalize(W[c], locked, basis);
    }
    A.push_back(std::move(Aj));

    // K-QR of the residual block (modified Gram-Schmidt). A column whose norm
    // vanishes has exhausted the invariant subspace along it: it leaves.
    KFamily next;
    MatX Rn = MatX::Zero(b, b);
    for (Eigen::Index c = 0; c < b; ++c) {
      VecX kw = K * W[c];
      for (std::size_t d = 0; d < next.size(); ++d) {
        double r = W[c].dot(next.kx[d]);
        W[c] -= r * next.x[d];
        kw -= r * next.kx[d];
        Rn(static_cast<Eigen::Index>(d), c) = r;
      }
      double nrm = std::sqrt(std::max(0.0, W[c].dot(kw)));
      if (nrm < 1e-14 || nrm < 1e-12 * mu_scale) continue;
      Rn(static_cast<Eigen::Index>(next.size()), c) = nrm;
      next.push(W[c] / nrm, kw / nrm);
    }
    Rn.conservativeResize(static_cast<Eigen::Index>(next.size()), b);
    R.push_back(Rn);

    // Early stop: every wanted Ritz value of the largest mu has converged.
    const int k = static_cast<int>(basis.size());
    if (opt.stop && next.size() > 0 && k >= opt.n_wanted + 10 &&
        static_cast<int>(A.size()) % opt.stop->check_every == 0 && k < m_max) {
      Eigen::SelfAdjointEigenSolver<MatX> ek(lanczos_matrix(A, R, A.size()));
      if (ek.info() == Eigen::Success) {
        const double top = std::fabs(ek.eigenvalues()(k - 1));
        bool done = true;
        for (int i = k - 1; i >= std::max(0, k - opt.n_wanted) && done; --i) {
          const VecX y_last = ek.eigenvectors().col(i).tail(b);
          done = ritz_residual(Rn, y_last) <= opt.stop->tol * top;
        }
        if (done) break;
      }
    }
    previous = std::move(block);
    block = std::move(next);
  }

  // Rayleigh-Ritz on the block-tridiagonal T.
  const MatX T = lanczos_matrix(A, R, A.size());
  const int m = static_cast<int>(T.rows());
  Eigen::SelfAdjointEigenSolver<MatX> es(T);
  require(es.info() == Eigen::Success, "lanczos: tridiagonal eigensolve failed");

  std::vector<PencilPair> out;
  const double mu_top = m > 0 ? std::fabs(es.eigenvalues()(m - 1)) : 0.0;
  int accepted = 0;
  for (int i = m - 1; i >= 0; --i) {  // descending mu
    PencilPair p;
    p.mu = es.eigenvalues()(i);
    if (opt.stop && opt.stop->lazy_vectors &&
        (accepted >= opt.n_wanted || !(p.mu > 1e-12 * std::max(1.0, mu_top))))
      break;
    VecX y = es.eigenvectors().col(i);
    VecX x = VecX::Zero(n);
    for (int k = 0; k < m; ++k) x += y(k) * basis.x[static_cast<std::size_t>(k)];
    p.x = std::move(x);
    p.residual = pencil_residual(K, B, p.mu, p.x);
    if (p.residual <= (opt.stop ? opt.stop->accept : 0.0)) ++accepted;
    out.push_back(std::move(p));
  }
  return out;
}

std::vector<PencilPair> k_inner_lanczos(const FactorCache& Kfac, const SpMat& K, const SpMat& B,
                                        int n_wanted, int max_iters, const LanczosStop* stop) {
  BlockLanczosOptions o;
  o.n_wanted = n_wanted;
  o.block_size = 1;
  o.max_columns = max_iters;
  o.stop = stop;
  return k_inner_block_lanczos(Kfac, K, B, o);
}

// Rigid-body influence vector for direction d over free DOFs (rotations about
// the global origin: translations e_axis x p, plus 1 on the parallel rotation).
VecX rigid_influence_vector(const Model& model, const DofMap& dofs, int d) {
  VecX r = VecX::Zero(dofs.n_free);
  for (int ni = 0; ni < dofs.n_nodes; ++ni) {
    const Vec3& p = model.nodes[static_cast<std::size_t>(ni)].x;
    Vec6 rn = Vec6::Zero();
    if (d < 3) {
      rn(d) = 1.0;
    } else {
      Vec3 e = Vec3::Zero();
      e(d - 3) = 1.0;
      Vec3 t = e.cross(p);
      rn.segment<3>(0) = t;
      rn(d) = 1.0;
    }
    for (int c = 0; c < 6; ++c) {
      int fi = dofs.free_ix[static_cast<std::size_t>(dofs.full(ni, c))];
      if (fi >= 0) r(fi) += rn(c);
    }
  }
  return r;
}

VecX ground_coupling(const Model& model, const Assembled& asm_, bool lumped, int d) {
  const DofMap& dm = asm_.dofs;
  VecX c = VecX::Zero(dm.n_free);
  if (lumped) return c;  // a lumped mass couples no dof to another
  // The mass on every dof that follows the ground (free or blocked): the
  // same assembly on a DOF map where the blocked dofs are numbered too.
  Assembled ext = asm_;
  int n = 0;
  bool any_blocked = false;
  for (std::size_t f = 0; f < dm.free_ix.size(); ++f) {
    const bool blocked = dm.free_ix[f] < 0 && f < dm.blocked.size() && dm.blocked[f];
    any_blocked = any_blocked || blocked;
    ext.dofs.free_ix[f] = dm.free_ix[f] >= 0 || blocked ? n++ : -1;
  }
  if (!any_blocked) return c;
  ext.dofs.n_free = n;
  // r on the blocked dofs only: M_fs r_s, the rows of the free dofs
  VecX rs = rigid_influence_vector(model, ext.dofs, d);
  for (std::size_t f = 0; f < dm.free_ix.size(); ++f)
    if (dm.free_ix[f] >= 0) rs(ext.dofs.free_ix[f]) = 0.0;
  const VecX g = assemble_mass(model, ext, false) * rs;
  for (std::size_t f = 0; f < dm.free_ix.size(); ++f)
    if (dm.free_ix[f] >= 0) c(dm.free_ix[f]) = g(ext.dofs.free_ix[f]);
  return c;
}

std::vector<Vec6> expand_unit_crest(const DofMap& dofs, const VecX& x) {
  // Scale so the largest-magnitude free component is exactly +1 (signed).
  double best = 0;
  for (Eigen::Index i = 0; i < x.size(); ++i)
    if (std::fabs(x(i)) > std::fabs(best)) best = x(i);
  double s = best != 0 ? 1.0 / best : 1.0;
  std::vector<Vec6> out(static_cast<std::size_t>(dofs.n_nodes), Vec6::Zero());
  for (int ni = 0; ni < dofs.n_nodes; ++ni)
    for (int c = 0; c < 6; ++c) {
      int fi = dofs.free_ix[static_cast<std::size_t>(dofs.full(ni, c))];
      if (fi >= 0) out[static_cast<std::size_t>(ni)](c) = s * x(fi);
    }
  return out;
}

namespace {

// omega^2 of a pair of the shifted pencil (K - sigma M) x = (1/mu) M x.
double omega_squared(const PencilPair& p, double shift) { return 1.0 / p.mu + shift; }

// A rigid mode of a free body: omega^2 = 0 to the accuracy of the shifted solve.
bool is_rigid_mode(const PencilPair& p, double shift) {
  return shift < 0 && p.mu > 0 && std::fabs(omega_squared(p, shift)) <= 1e-6 * std::fabs(shift);
}

}  // namespace

ModalResults modal_from_pairs(const Model& model, const Assembled& asm_, const SpMat& M,
                              const std::vector<PencilPair>& pairs, int n_modes,
                              bool support_coupling, double shift) {
  ModalResults r;
  r.shift = shift;
  // Massless directions have mu -> 0; the floor is relative to the largest mu
  // of an elastic mode (a rigid mode's mu = 1/|sigma| is no scale).
  double mu_max = 0.0;
  for (const PencilPair& p : pairs)
    if (!is_rigid_mode(p, shift)) mu_max = std::max(mu_max, std::fabs(p.mu));
  if (mu_max == 0.0 && !pairs.empty()) mu_max = std::fabs(pairs.front().mu);
  int rejected = 0;
  for (const PencilPair& p : pairs) {
    if (static_cast<int>(r.omegas.size()) >= n_modes) break;
    if (!(p.mu > 1e-12 * std::max(1.0, mu_max))) break;  // massless directions: mu -> 0
    if (p.residual > 1e-6) {                              // unconverged Ritz artifact
      ++rejected;
      continue;
    }
    double om = 0.0;  // a rigid mode
    if (!is_rigid_mode(p, shift)) {
      const double om2 = omega_squared(p, shift);
      require(om2 > 0, "modal: negative eigenvalue omega^2 = " + std::to_string(om2) +
                           " rad^2/s^2 — the stiffness is not positive semi-definite (beyond-critical "
                           "axial state)");
      om = std::sqrt(om2);
    }
    // M-normalize.
    VecX phi = p.x;
    double mm = phi.dot(M * phi);
    require(mm > 0, "modal: non-positive modal mass");
    phi /= std::sqrt(mm);
    // Deterministic sign: the largest-magnitude component positive, the
    // crest convention of the unit-crest shapes, so participation factors and
    // signed spectral fields do not depend on the eigensolver's sign choice.
    Eigen::Index crest = 0;
    for (Eigen::Index i = 1; i < phi.size(); ++i)
      if (std::fabs(phi(i)) > std::fabs(phi(crest))) crest = i;
    if (phi.size() > 0 && phi(crest) < 0) phi = -phi;
    r.omegas.push_back(om);
    r.frequencies_hz.push_back(om / (2 * kPi));
    r.modes_free.push_back(phi);
    r.modes.push_back(expand_unit_crest(asm_.dofs, phi));
    r.residuals.push_back(p.residual);
  }
  if (rejected > 0)  // a rejected pair lies below modes returned: they are not the lowest ones
    r.warning = "modal: " + std::to_string(rejected) +
                " eigenpair(s) rejected for a pencil residual above 1e-6 (ill-conditioned stiffness, e.g. "
                "very soft springs): the modes returned are not the lowest ones";

  for (int d = 0; d < 6; ++d) {
    VecX rd = rigid_influence_vector(model, asm_.dofs, d);
    r.total_mass[static_cast<std::size_t>(d)] = rd.dot(M * rd);
  }
  // Gamma = phi^T (M_ff r_f + M_fs r_s): the moving supports' consistent
  // mass coupling included (ground_coupling; zero for a lumped mass)
  std::array<VecX, 6> cpl;
  std::array<bool, 6> has_cpl{};
  if (support_coupling && !r.modes_free.empty())
    for (int d = 0; d < 6; ++d) {
      cpl[static_cast<std::size_t>(d)] = ground_coupling(model, asm_, false, d);
      has_cpl[static_cast<std::size_t>(d)] = cpl[static_cast<std::size_t>(d)].squaredNorm() > 0;
    }
  // With the coupling, the mass the whole spectrum mobilises is c^T M_ff^-1 c
  // for c = M_ff r_f + M_fs r_s (the sum of Gamma^2 over EVERY mode of the
  // pencil: the effective-mass ratios then sum to 1 over the full spectrum
  // and never exceed it on a truncated one), not r_f^T M_ff r_f. Solved on
  // the dofs carrying mass (M is PSD with zero rows on the massless dofs,
  // where c is zero too; on the others it is definite). Without coupling
  // the two agree and the quadratic form is kept (bit-identical).
  // (one factorization of the massive block serves the six directions)
  std::vector<int> keep;
  Eigen::SimplicialLDLT<SpMat, Eigen::Lower, KetaOrdering> ldlt;
  SupernodalLDLT sn;  // large models (keta/supernodal.h)
  bool use_sn = false;
  if (std::any_of(has_cpl.begin(), has_cpl.end(), [](bool b) { return b; })) {
    std::vector<int> at(static_cast<std::size_t>(M.rows()), -1);
    for (int i = 0; i < M.rows(); ++i)
      if (M.coeff(i, i) > 0) {
        at[static_cast<std::size_t>(i)] = static_cast<int>(keep.size());
        keep.push_back(i);
      }
    std::vector<Eigen::Triplet<double>> tr;
    tr.reserve(static_cast<std::size_t>(M.nonZeros()));
    for (int k = 0; k < M.outerSize(); ++k)
      for (SpMat::InnerIterator it(M, k); it; ++it) {
        const int a = at[static_cast<std::size_t>(it.row())], b = at[static_cast<std::size_t>(it.col())];
        if (a >= 0 && b >= 0) tr.emplace_back(a, b, it.value());
      }
    SpMat Mn(static_cast<Eigen::Index>(keep.size()), static_cast<Eigen::Index>(keep.size()));
    Mn.setFromTriplets(tr.begin(), tr.end());
    if (Mn.rows() >= kSupernodalMin) {
      sn.analyze(Mn);
      require(sn.factorize(Mn), "modal: the mass matrix of the massive dofs is not definite");
      use_sn = true;
    } else {
      ldlt.compute(Mn);
      require(ldlt.info() == Eigen::Success, "modal: the mass matrix of the massive dofs is not definite");
    }
  }
  for (int d = 0; d < 6; ++d) {
    if (!has_cpl[static_cast<std::size_t>(d)]) continue;
    const VecX c = M * rigid_influence_vector(model, asm_.dofs, d) + cpl[static_cast<std::size_t>(d)];
    VecX cn(static_cast<Eigen::Index>(keep.size()));
    for (std::size_t k = 0; k < keep.size(); ++k) cn(static_cast<Eigen::Index>(k)) = c(keep[k]);
    r.total_mass[static_cast<std::size_t>(d)] = cn.dot(use_sn ? sn.solve(cn) : VecX(ldlt.solve(cn)));
  }
  for (const VecX& phi : r.modes_free) {
    std::array<double, 6> gam{}, eff{};
    VecX Mphi = M * phi;
    for (int d = 0; d < 6; ++d) {
      VecX rd = rigid_influence_vector(model, asm_.dofs, d);
      gam[static_cast<std::size_t>(d)] = Mphi.dot(rd);
      if (has_cpl[static_cast<std::size_t>(d)])
        gam[static_cast<std::size_t>(d)] += phi.dot(cpl[static_cast<std::size_t>(d)]);
      eff[static_cast<std::size_t>(d)] =
          gam[static_cast<std::size_t>(d)] * gam[static_cast<std::size_t>(d)];
    }
    r.participation.push_back(gam);
    r.effective_mass.push_back(eff);
  }
  return r;
}

namespace {

// The pencil of a modal analysis: the free-dof stiffness (plus the geometric
// stiffness of a prestress state) and mass, and the pair actually solved —
// reduced by the rigid links (T^t K T, T^t M T) when the model has some.
struct ModalPencil {
  Assembled asm_;
  SpMat K, M;    // free dofs; M also serves the participation factors
  SpMat Kr, Mr;  // reduced by the rigid links (empty without links)

  bool reduced() const { return asm_.red.active; }
  const SpMat& Ks() const { return reduced() ? Kr : K; }
  const SpMat& Ms() const { return reduced() ? Mr : M; }
};

ModalPencil build_modal_pencil(const Model& model, bool lumped_mass,
                               const std::vector<double>* axial_state, const SpMat* stiffness = nullptr) {
  ModalPencil p;
  p.asm_ = assemble_elements(model, true, axial_state);
  require(p.asm_.dofs.n_free > 0, "modal: no free DOFs");
  if (stiffness) {
    require(stiffness->rows() == p.asm_.dofs.n_free && stiffness->cols() == p.asm_.dofs.n_free,
            "modal: the given stiffness has " + std::to_string(stiffness->rows()) + " rows for " +
                std::to_string(p.asm_.dofs.n_free) + " free DOFs");
    p.K = *stiffness;
  } else {
    p.K = assemble_stiffness(model, p.asm_);
  }
  if (axial_state && !stiffness) {
    require(model.continua.empty(), "modal: a prestressed modal analysis does not take continuum elements (no Kg)");
    require(!links_have_arms(model),
            "modal: a prestressed modal analysis does not take links with arms (their second-order "
            "term is not modelled)");
    SpMat Kg = assemble_geometric(model, p.asm_, *axial_state);
    p.K = (p.K + Kg).eval();  // vibration about the loaded state
  }
  p.M = assemble_mass(model, p.asm_, lumped_mass);
  require(p.M.nonZeros() > 0, "modal: model has no mass (set rho or add point masses)");
  if (p.reduced()) {
    p.Kr = p.asm_.red.reduce(p.K);
    p.Mr = p.asm_.red.reduce(p.M);
  }
  return p;
}

SpMat shifted(const SpMat& K, const SpMat& M, double sigma) {
  if (sigma == 0.0) return K;
  SpMat A = (K - sigma * M).eval();
  A.makeCompressed();
  return A;
}

// K is positive definite: no negative or zero pivot, and no pivot so small
// against the largest that the matrix is singular to rounding (a free body
// shows its rigid modes as such pivots, of either sign).
bool positive_definite(const FactorCache& f) {
  const Inertia in = f.inertia();
  return in.n_neg == 0 && in.n_zero == 0 && f.min_abs_pivot() > 1e-13 * f.max_abs_pivot();
}

// The n lowest eigenpairs of (Ks - sigma Ms, Ms) by the block Lanczos,
// certified complete by the Sturm count; deflated restarts while pairs are
// missing below the cut.
struct CertifiedPairs {
  std::vector<PencilPair> pairs;  // descending mu (ascending omega^2)
  int sturm_count = -1;
  std::string warning;
};

CertifiedPairs certified_lanczos(const FactorCache& fac, const SpMat& Ks_shifted, const SpMat& Ks,
                                 const SpMat& Ms, double sigma, int n_modes, int block_size) {
  constexpr int kMaxRestarts = 10;
  const LanczosStop stop;  // converged modes only, vectors for those consumed
  const auto omega2 = [&](const PencilPair& p) { return 1.0 / p.mu + sigma; };
  // Upper bound of the finite eigenvalues: the dofs that carry mass.
  int n_massive = 0;
  for (Eigen::Index i = 0; i < Ms.rows(); ++i)
    if (Ms.coeff(i, i) > 0) ++n_massive;
  const int n_wanted = std::min(n_modes, n_massive);

  CertifiedPairs out;
  std::vector<VecX> locked, k_locked;
  int missing = n_wanted;
  int runs = 0;
  for (int attempt = 0; attempt <= kMaxRestarts && missing > 0; ++attempt) {
    ++runs;
    BlockLanczosOptions o;
    o.n_wanted = missing;
    o.block_size = block_size;
    o.seed = static_cast<std::uint64_t>(attempt);
    o.stop = &stop;
    std::size_t before = out.pairs.size();
    for (PencilPair& p : k_inner_block_lanczos(fac, Ks_shifted, Ms, o, locked, k_locked)) {
      if (!(p.mu > 0) || p.residual > stop.accept) continue;
      VecX kx = Ks_shifted * p.x;
      const double nrm = std::sqrt(std::max(0.0, p.x.dot(kx)));
      if (!(nrm > 0)) continue;
      locked.push_back(p.x / nrm);
      k_locked.push_back(kx / nrm);
      out.pairs.push_back(std::move(p));
    }
    std::sort(out.pairs.begin(), out.pairs.end(),
              [](const PencilPair& a, const PencilPair& b) { return a.mu > b.mu; });
    if (out.pairs.empty()) break;
    if (out.pairs.size() == before && attempt > 0) break;  // no progress: stop restarting

    // Sturm certificate: every eigenvalue below a cut just above the highest
    // mode kept must be among the pairs found.
    const std::size_t kept = std::min<std::size_t>(static_cast<std::size_t>(n_wanted), out.pairs.size());
    const double top = omega2(out.pairs[kept - 1]);
    const double cut = top + 1e-6 * std::max(std::fabs(top), std::fabs(sigma)) + 1e-300;
    out.sturm_count = eigenvalue_count_below(Ks, Ms, cut);
    const int found = static_cast<int>(std::count_if(out.pairs.begin(), out.pairs.end(),
                                                     [&](const PencilPair& p) { return omega2(p) < cut; }));
    missing = std::max(out.sturm_count - found, n_wanted - static_cast<int>(out.pairs.size()));
  }
  if (static_cast<int>(out.pairs.size()) > n_wanted) out.pairs.resize(static_cast<std::size_t>(n_wanted));
  if (missing > 0)
    out.warning = "modal: the Sturm count reports " + std::to_string(out.sturm_count) +
                  " eigenvalue(s) up to the highest mode found but the Lanczos certified only " +
                  std::to_string(out.pairs.size()) + " of the " + std::to_string(n_modes) +
                  " mode(s) asked, after " + std::to_string(runs) + " Lanczos run(s)";
  return out;
}

}  // namespace

double automatic_shift(const SpMat& K, const SpMat& M) {
  double ratio = 0.0;
  for (Eigen::Index i = 0; i < K.rows(); ++i) {
    const double m = M.coeff(i, i);
    if (m > 0) ratio = std::max(ratio, K.coeff(i, i) / m);
  }
  require(ratio > 0, "modal: no dof carries both stiffness and mass");
  return -1e-8 * ratio;
}

int eigenvalue_count_below(const SpMat& K, const SpMat& M, double c) {
  for (int attempt = 0; attempt < 4; ++attempt) {
    FactorCache f;
    try {
      f.factorize(shifted(K, M, c));
      return f.inertia().n_neg;
    } catch (const Error&) {  // c on an eigenvalue: a zero pivot
      c += 1e-9 * std::max(std::fabs(c), 1e-300);
    }
  }
  fail("modal: Sturm count — K - c M stays singular near c = " + std::to_string(c));
}

namespace {

// One solve of the pencil (Ks - sigma Ms, Ms). `factor_of_ks`: the factor of
// the unshifted Ks when the caller has it (sigma == 0, Lanczos path).
ModalResults solve_shifted(const Model& model, const ModalPencil& pen, const ModalOptions& opt,
                           double sigma, const FactorCache* factor_of_ks) {
  const SpMat& Ks = pen.Ks();
  const SpMat& Ms = pen.Ms();
  const SpMat K_shifted = shifted(Ks, Ms, sigma);
  std::vector<PencilPair> pairs;
  CertifiedPairs certified;
  if (Ks.rows() <= opt.dense_threshold) {
    pairs = dense_pencil(K_shifted, Ms);
  } else {
    FactorCache own;
    const FactorCache* fac = factor_of_ks;
    if (!fac) {
      own.factorize(K_shifted);
      fac = &own;
    }
    require(positive_definite(*fac),
            "modal: K - sigma M is not positive definite (sigma = " + std::to_string(sigma) +
                ") — the stiffness has negative eigenvalues (beyond-critical axial state) or a massless "
                "mechanism");
    certified = certified_lanczos(*fac, K_shifted, Ks, Ms, sigma, opt.n_modes, opt.block_size);
    pairs = std::move(certified.pairs);
  }
  if (pen.reduced())
    for (PencilPair& p : pairs) p.x = pen.asm_.red.expand(p.x);
  ModalResults r = modal_from_pairs(model, pen.asm_, pen.M, pairs, opt.n_modes,
                                    opt.support_coupling && !opt.lumped_mass, sigma);
  r.sturm_count = certified.sturm_count;
  if (!certified.warning.empty()) r.warning = certified.warning;
  return r;
}

// K positive definite (the unshifted pencil applies). On the Lanczos path
// the factor is kept in `fac` for the solve.
bool stiffness_positive_definite(const ModalPencil& pen, const ModalOptions& opt, FactorCache& fac) {
  const SpMat& Ks = pen.Ks();
  if (Ks.rows() <= opt.dense_threshold) return Eigen::LLT<MatX>(MatX(Ks)).info() == Eigen::Success;
  try {
    fac.factorize(Ks);
  } catch (const Error&) {  // a zero pivot: singular
    return false;
  }
  return positive_definite(fac);
}

}  // namespace

ModalResults solve_modal(const Model& model, const ModalOptions& opt) {
  hooks_check(opt.hooks);
  model.validate();
  hooks_report(opt.hooks, 0.0, "assembly");
  require(opt.block_size > 0, "modal: block_size must be positive");
  const ModalPencil pen = build_modal_pencil(model, opt.lumped_mass, opt.axial_state, opt.stiffness);
  hooks_check(opt.hooks);
  hooks_report(opt.hooks, 0.5, pen.Ks().rows() <= opt.dense_threshold ? "eigensolve" : "Lanczos");

  ModalResults r;
  if (!std::isnan(opt.shift)) {
    r = solve_shifted(model, pen, opt, opt.shift, nullptr);
  } else {
    // The unshifted pencil when K is positive definite and it delivers a
    // certified spectrum; otherwise (a free body, or a stiffness so
    // ill-conditioned that pairs fail their residual test) the automatic
    // shift, which bounds the conditioning of K - sigma M.
    FactorCache fac;
    const bool unshifted = stiffness_positive_definite(pen, opt, fac);
    if (unshifted) r = solve_shifted(model, pen, opt, 0.0, pen.Ks().rows() <= opt.dense_threshold ? nullptr : &fac);
    if (!unshifted || !r.warning.empty())
      r = solve_shifted(model, pen, opt, automatic_shift(pen.Ks(), pen.Ms()), nullptr);
  }
  hooks_report(opt.hooks, 1.0, "done");
  return r;
}

BucklingResults solve_buckling(const Model& model, const BucklingOptions& opt) {
  require(model.continua.empty(), "solve_buckling: continuum elements are not supported (linear statics, modal and dynamics take them)");
  hooks_check(opt.hooks);
  model.validate();
  require(!model.cases.empty(), "buckling: no load cases");
  const LoadCase& lc = opt.case_id == 0 ? model.cases.front() : model.load_case(opt.case_id);

  hooks_report(opt.hooks, 0.0, "pre-stress");
  StaticContext ctx = StaticContext::build(model, true);
  for (const ElasticLink& l : model.links)
    require((model.node(l.nj).x - model.node(l.ni).x).norm() == 0,
            "buckling: an elastic link between distinct nodes carries a second-order term that is not "
            "modelled; coincident-node elastic links and rigid links (with arms) are");
  ctx.cache.require_spd();
  StaticOptions sopt;
  sopt.stations_per_bar = 2;
  CaseResult ref = ctx.solve_case(lc, sopt);
  std::vector<double> N = axial_field(ref);
  std::vector<Vec3> sn = shell_membrane_field(ref);

  bool any = false;
  for (double v : N)
    if (v != 0) any = true;
  for (const Vec3& v : sn)
    if (v.squaredNorm() != 0) any = true;
  BucklingResults out;
  if (!any) return out;  // zero prestress field: nothing to buckle

  SpMat Kg = assemble_geometric(model, ctx.asm_, N, &sn);
  // Rigid arms of the links: a slave at r from its master moves by
  // theta x r + 1/2 theta x (theta x r); the force P the slave hands to the
  // rigid arm (its loads minus its elements' forces, the constraint force of
  // the reference state) then adds to the master's rotations the second
  // variation of -P . u2 (dead P):  (P . r) I - (P r' + r P') / 2
  // (an inverted pendulum: -W h on the two horizontal rotations).
  if (!model.rigid_links.empty()) {
    CaseLoads cl = lower_case(model, ctx.asm_, lc);
    fold_prescribed(model, ctx.asm_, cl);
    fold_link_offsets(model, ctx.asm_, ctx.K, cl);
    const VecX v = ctx.cache.solve(cl.rhs);
    const VecX R = cl.rhs - ctx.K * v;
    const DofMap& d = ctx.asm_.dofs;
    std::vector<Triplet> t;
    for (int k = 0; k < Kg.outerSize(); ++k)
      for (SpMat::InnerIterator it(Kg, k); it; ++it) t.emplace_back(it.row(), it.col(), it.value());
    for (const RigidLink& l : model.rigid_links) {
      const int m = model.node_index(l.master);
      for (Id sid : l.slaves) {
        const int sn = model.node_index(sid);
        const Vec3 r = model.node(sid).x - model.node(l.master).x;
        if (r.squaredNorm() == 0) continue;
        Vec3 P = Vec3::Zero();
        for (int c = 0; c < 3; ++c) {
          const int fi = d.free_of(sn, c);
          if (l.dofs[static_cast<std::size_t>(c)] && fi >= 0) P(c) = R(fi);
        }
        const Mat3 B = P.dot(r) * Mat3::Identity() - 0.5 * (P * r.transpose() + r * P.transpose());
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) {
            const int fa = d.free_of(m, 3 + a), fb = d.free_of(m, 3 + b);
            if (fa >= 0 && fb >= 0 && B(a, b) != 0) t.emplace_back(fa, fb, B(a, b));
          }
      }
    }
    Kg.setFromTriplets(t.begin(), t.end());
  }
  std::vector<PencilPair> pairs;
  hooks_check(opt.hooks);
  hooks_report(opt.hooks, 0.5, "eigensolve");
  if (ctx.asm_.red.active) {
    // rigid links: the reduced pencil, buckling modes expanded by T
    const SpMat Kr = ctx.asm_.red.reduce(ctx.K), Kgr = ctx.asm_.red.reduce(Kg);
    if (Kr.rows() <= opt.dense_threshold) {
      pairs = dense_pencil(Kr, Kgr);
    } else {
      FactorCache fr;
      fr.factorize(Kr);
      pairs = k_inner_lanczos(fr, Kr, Kgr, opt.n_modes);
    }
    for (PencilPair& p : pairs) p.x = ctx.asm_.red.expand(p.x);
  } else if (ctx.asm_.dofs.n_free <= opt.dense_threshold) {
    pairs = dense_pencil(ctx.K, Kg);
  } else {
    pairs = k_inner_lanczos(ctx.cache, ctx.K, Kg, opt.n_modes);
  }
  // K + lambda Kg singular  <=>  Kg x = -(1/lambda) K x. Smallest positive
  // lambda = most negative mu. Collect ascending lambda.
  std::vector<std::pair<double, const PencilPair*>> lam;
  for (const PencilPair& p : pairs)
    if (p.mu < -1e-14 && p.residual < 1e-6) lam.emplace_back(-1.0 / p.mu, &p);
  std::sort(lam.begin(), lam.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  for (const auto& [l, p] : lam) {
    if (static_cast<int>(out.load_factors.size()) >= opt.n_modes) break;
    out.load_factors.push_back(l);
    out.modes.push_back(expand_unit_crest(ctx.asm_.dofs, p->x));
    out.residuals.push_back(p->residual);
    // the unit-crest scaling of expand_unit_crest, then x^T K x
    double best = 0;
    for (Eigen::Index i = 0; i < p->x.size(); ++i)
      if (std::fabs(p->x(i)) > std::fabs(best)) best = p->x(i);
    const VecX xs = best != 0 ? VecX(p->x / best) : p->x;
    const VecX kgx = -(Kg * xs);
    const double g = xs.dot(kgx);
    out.gen_stiffness.push_back(g);
    out.participation.push_back(g > 0 ? kgx.sum() / std::sqrt(g) : 0.0);
  }
  {
    const VecX ones = VecX::Ones(Kg.rows());
    out.total_kg = -ones.dot(Kg * ones);
  }
  hooks_report(opt.hooks, 1.0, "done");
  return out;
}

double buckling_factor_bisection(const Model& model, Id case_id, double lambda_max,
                                 double tol_rel, bool exact_stability) {
  require(model.continua.empty(), "buckling_factor_bisection: continuum elements are not supported (linear statics, modal and dynamics take them)");
  model.validate();
  require(!model.cases.empty(), "buckling: no load cases");
  const LoadCase& lc = case_id == 0 ? model.cases.front() : model.load_case(case_id);
  StaticContext ctx = StaticContext::build(model);
  StaticOptions sopt;
  sopt.stations_per_bar = 2;
  CaseResult ref = ctx.solve_case(lc, sopt);
  std::vector<double> N = axial_field(ref);
  std::vector<Vec3> sn = shell_membrane_field(ref);
  SpMat Kg = assemble_geometric(model, ctx.asm_, N, &sn);

  auto flipped = [&](double lambda) {
    SpMat Kt;
    if (exact_stability) {
      // exact stability-function tangent at the scaled axial state: the first
      // inertia flip is the EXACT critical multiplier even with one element
      // per member (ST1-style closed-form beam-column).
      std::vector<Triplet> trips;
      for (std::size_t e = 0; e < ctx.asm_.elements.size(); ++e) {
        const Bar& b = model.bars[e];
        if (b.kind == BarKind::Cable) continue;
        Mat12 ke;
        try {
          ke = ctx.asm_.elements[e].ke_exact(lambda * N[e]);
        } catch (const Error&) {
          return true;  // a member driven past its own trig range is beyond
        }
        int ni = model.node_index(b.ni), nj = model.node_index(b.nj);
        int g[12];
        for (int c = 0; c < 6; ++c) {
          g[c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(ctx.asm_.dofs.full(ni, c))];
          g[6 + c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(ctx.asm_.dofs.full(nj, c))];
        }
        for (int a = 0; a < 12; ++a)
          for (int bb = 0; bb < 12; ++bb)
            if (g[a] >= 0 && g[bb] >= 0 && ke(a, bb) != 0)
              trips.emplace_back(g[a], g[bb], ke(a, bb));
      }
      for (std::size_t f = 0; f < ctx.asm_.spring.size(); ++f)
        if (ctx.asm_.spring[f] != 0) {
          int fi = ctx.asm_.dofs.free_ix[f];
          if (fi >= 0) trips.emplace_back(fi, fi, ctx.asm_.spring[f]);
        }
      scatter_full_springs(ctx.asm_, trips);  // ST1 APPUI EL
      // Facet shells stay LINEAR under the exact bar tangents: ke plus the
      // lambda-scaled membrane geometric stiffness (same consistent-Kg form
      // the non-exact branch uses — no stability-function analogue exists for
      // a plate, and the ST1 exact path is a BAR feature).
      for (const ShellElem& se : ctx.asm_.shells) {
        const Vec3& s_n = sn[se.index];
        Mat18 kes = se.el.ke();
        if (s_n.squaredNorm() != 0)
          kes += lambda * se.el.kg(s_n(0), s_n(1), s_n(2));
        int gd[18];
        for (int k = 0; k < 3; ++k)
          for (int c = 0; c < 6; ++c)
            gd[6 * k + c] = ctx.asm_.dofs.free_ix[static_cast<std::size_t>(
                ctx.asm_.dofs.full(se.nodes[static_cast<std::size_t>(k)], c))];
        for (int a = 0; a < 18; ++a)
          for (int bb = 0; bb < 18; ++bb)
            if (gd[a] >= 0 && gd[bb] >= 0 && kes(a, bb) != 0)
              trips.emplace_back(gd[a], gd[bb], kes(a, bb));
      }
      Kt = SpMat(ctx.asm_.dofs.n_free, ctx.asm_.dofs.n_free);
      Kt.setFromTriplets(trips.begin(), trips.end());
    } else {
      Kt = ctx.K + lambda * Kg;
    }
    FactorCache f;
    try {
      f.factorize(Kt);
    } catch (const Error&) {
      return true;  // exactly singular: counted as flipped
    }
    return f.inertia().n_neg > 0;
  };

  require(!flipped(0.0), "buckling bisection: base state already beyond critical");
  double lo = 0.0, hi = 1.0;
  while (!flipped(hi)) {
    lo = hi;
    hi *= 2;
    require(hi <= lambda_max, "buckling bisection: no critical factor below cap");
  }
  while (hi - lo > tol_rel * hi) {
    double mid = 0.5 * (lo + hi);
    (flipped(mid) ? hi : lo) = mid;
  }
  return 0.5 * (lo + hi);
}

int mode_count_in_band(const Model& model, double f_lo_hz, double f_hi_hz, bool lumped_mass) {
  const ModalPencil pen = build_modal_pencil(model, lumped_mass, nullptr);
  const SpMat& Ks = pen.Ks();
  const SpMat& Ms = pen.Ms();
  const auto omega2 = [](double f_hz) { return std::pow(2 * kPi * f_hz, 2); };
  // [0, f]: from just below zero, so that the rigid modes of a free body
  // (omega = 0, a singular K) are counted in the band.
  const double lower = f_lo_hz > 0 ? omega2(f_lo_hz) : automatic_shift(Ks, Ms);
  return eigenvalue_count_below(Ks, Ms, omega2(f_hi_hz)) - eigenvalue_count_below(Ks, Ms, lower);
}

Model apply_imperfection(const Model& model, const std::vector<Vec6>& shape, double scale) {
  require(shape.size() == model.nodes.size(), "imperfection: shape/node count mismatch");
  require(std::isfinite(scale), "imperfection: non-finite scale");
  Model out = model;
  for (std::size_t i = 0; i < out.nodes.size(); ++i)
    out.nodes[i].x += scale * shape[i].segment<3>(0);
  return out;
}

}  // namespace keta

namespace keta::ltb {

double critical_uniform_moment(double L, double E, double Iz, double G, double J, double Iw) {
  keta::require(L > 0 && E > 0 && Iz > 0 && G > 0 && J > 0 && Iw >= 0, "ltb: bad input");
  double k = keta::kPi / L;
  return k * std::sqrt(E * Iz * G * J) * std::sqrt(1.0 + k * k * E * Iw / (G * J));
}

double critical_uniform_moment_monosym(double L, double E, double Iz, double G, double J,
                                       double Iw, double beta_x) {
  keta::require(L > 0 && E > 0 && Iz > 0 && G > 0 && J > 0, "ltb: bad input");
  double k = keta::kPi / L;
  double k2EIz = k * k * E * Iz;
  double half = beta_x / 2.0;
  return k2EIz * (half + std::sqrt(half * half + (E * Iw * k * k + G * J) / k2EIz));
}

}  // namespace keta::ltb

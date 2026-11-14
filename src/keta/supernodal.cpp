// Supernodal multifrontal LDL^t — see supernodal.h.
#include "keta/supernodal.h"

#include <algorithm>
#include <cmath>

#include "keta/ordering.h"

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace keta {

namespace {
std::atomic<int> g_threads{1};
// run fn(0 .. nb - 1) on up to T threads (block b on thread b mod T); the
// blocks are independent, so the result does not depend on T
template <class Fn>
void parallel_blocks(int nb, int T, const Fn& fn) {
  const int nt = std::min(T, nb);
  if (nt <= 1) {
    for (int b = 0; b < nb; ++b) fn(b);
    return;
  }
  std::vector<std::thread> th;
  for (int t = 0; t < nt; ++t)
    th.emplace_back([&, t] {
      for (int b = t; b < nb; b += nt) fn(b);
    });
  for (std::thread& x : th) x.join();
}
}  // namespace

void set_solver_threads(int n) {
  require(n >= 0, "solver threads: n >= 0 (0 = all hardware threads)");
  g_threads = n == 0 ? std::max(1, static_cast<int>(std::thread::hardware_concurrency())) : n;
}
int solver_threads() { return g_threads; }

void SupernodalLDLT::analyze(const SpMat& Ain) {
  n_ = Ain.rows();
  require(Ain.cols() == n_, "supernodal: a square matrix");
  // the full symmetric pattern (from the lower triangle)
  SpMat L = Ain.triangularView<Eigen::Lower>();
  SpMat A = L.selfadjointView<Eigen::Lower>();
  A.makeCompressed();
  KetaOrdering::PermutationType P;
  KetaOrdering()(A, P);
  perm_ = P.indices();
  iperm_.resize(n_);
  for (Eigen::Index k = 0; k < n_; ++k) iperm_(perm_(k)) = static_cast<int>(k);
  // elimination tree and column counts of the permuted matrix (up-looking)
  std::vector<int> parent(static_cast<std::size_t>(n_), -1), tag(static_cast<std::size_t>(n_)), cnt(static_cast<std::size_t>(n_), 0);
  for (Eigen::Index k = 0; k < n_; ++k) {
    tag[static_cast<std::size_t>(k)] = static_cast<int>(k);
    for (SpMat::InnerIterator it(A, perm_(k)); it; ++it) {
      int i = iperm_(static_cast<Eigen::Index>(it.row()));
      if (i >= k) continue;
      while (tag[static_cast<std::size_t>(i)] != k) {
        if (parent[static_cast<std::size_t>(i)] == -1) parent[static_cast<std::size_t>(i)] = static_cast<int>(k);
        ++cnt[static_cast<std::size_t>(i)];
        tag[static_cast<std::size_t>(i)] = static_cast<int>(k);
        i = parent[static_cast<std::size_t>(i)];
      }
    }
  }
  // fundamental supernodes: j + 1 = parent(j), j its only child, same structure below
  std::vector<int> nchild(static_cast<std::size_t>(n_), 0);
  for (Eigen::Index j = 0; j < n_; ++j)
    if (parent[static_cast<std::size_t>(j)] >= 0) ++nchild[static_cast<std::size_t>(parent[static_cast<std::size_t>(j)])];
  sn_start_.clear();
  std::vector<int> sn_of(static_cast<std::size_t>(n_));
  for (Eigen::Index j = 0; j < n_; ++j) {
    const bool merge = j > 0 && parent[static_cast<std::size_t>(j - 1)] == j && nchild[static_cast<std::size_t>(j)] == 1 &&
                       cnt[static_cast<std::size_t>(j - 1)] == cnt[static_cast<std::size_t>(j)] + 1 &&
                       static_cast<int>(j) - sn_start_.back() < 512;
    if (!merge) sn_start_.push_back(static_cast<int>(j));
    sn_of[static_cast<std::size_t>(j)] = static_cast<int>(sn_start_.size()) - 1;
  }
  // relaxed amalgamation: a supernode joins the next one when that one holds
  // its parent column and the explicit zeros stay few (or both are small):
  // larger dense blocks for the BLAS-3 kernels
  {
    sn_start_.push_back(static_cast<int>(n_));
    std::vector<int> merged{sn_start_[0]};
    for (std::size_t q = 1; q + 1 < sn_start_.size(); ++q) {
      const int f0 = merged.back(), b = sn_start_[q], e = sn_start_[q + 1];
      const bool parent_in = parent[static_cast<std::size_t>(b - 1)] == b;
      bool join = false;
      if (parent_in && e - f0 <= 256) {
        // zeros added to the columns [f0, b): rows of the merged block below each
        double nz = 0, zeros = 0;
        const double below_e = cnt[static_cast<std::size_t>(e - 1)];
        for (int c = f0; c < e; ++c) {
          const double full = (e - 1 - c) + below_e;
          nz += full;
          zeros += full - cnt[static_cast<std::size_t>(c)];
        }
        join = (e - f0 <= 16) || zeros <= 0.05 * nz;
      }
      if (!join) merged.push_back(b);
    }
    sn_start_ = merged;
    for (std::size_t q = 0; q < sn_start_.size(); ++q) {
      const int e = q + 1 < sn_start_.size() ? sn_start_[q + 1] : static_cast<int>(n_);
      for (int c = sn_start_[q]; c < e; ++c) sn_of[static_cast<std::size_t>(c)] = static_cast<int>(q);
    }
  }
  const int ns = static_cast<int>(sn_start_.size());
  sn_start_.push_back(static_cast<int>(n_));
  sn_parent_.assign(static_cast<std::size_t>(ns), -1);
  for (int s = 0; s < ns; ++s) {
    const int last = sn_start_[static_cast<std::size_t>(s) + 1] - 1;
    const int p = parent[static_cast<std::size_t>(last)];
    if (p >= 0) sn_parent_[static_cast<std::size_t>(s)] = sn_of[static_cast<std::size_t>(p)];
  }
  // the permuted lower entries per column, and their original positions
  origpos_.clear();
  colmap_.assign(static_cast<std::size_t>(n_), {});
  for (int c = 0; c < L.outerSize(); ++c)
    for (SpMat::InnerIterator it(L, c); it; ++it) {
      const int r = static_cast<int>(it.row());
      int pr = iperm_(r), pc = iperm_(c);
      if (pr < pc) std::swap(pr, pc);
      colmap_[static_cast<std::size_t>(pc)].push_back({pr, static_cast<int>(origpos_.size())});
      origpos_.push_back({r, c});
    }
  // row structures: own pattern below the pivots, plus the children's
  rows_.assign(static_cast<std::size_t>(ns), {});
  std::vector<std::vector<int>> kids(static_cast<std::size_t>(ns));
  for (int s = 0; s < ns; ++s)
    if (sn_parent_[static_cast<std::size_t>(s)] >= 0) kids[static_cast<std::size_t>(sn_parent_[static_cast<std::size_t>(s)])].push_back(s);
  std::vector<int> mark(static_cast<std::size_t>(n_), -1);
  for (int s = 0; s < ns; ++s) {
    const int f = sn_start_[static_cast<std::size_t>(s)], e = sn_start_[static_cast<std::size_t>(s) + 1];
    std::vector<int>& R = rows_[static_cast<std::size_t>(s)];
    for (int c = f; c < e; ++c) {
      R.push_back(c);
      mark[static_cast<std::size_t>(c)] = s;
    }
    std::vector<int> below;
    auto add = [&](int r) {
      if (r >= e && mark[static_cast<std::size_t>(r)] != s) {
        mark[static_cast<std::size_t>(r)] = s;
        below.push_back(r);
      }
    };
    for (int c = f; c < e; ++c)
      for (const auto& [r, pos] : colmap_[static_cast<std::size_t>(c)]) add(r);
    for (int k : kids[static_cast<std::size_t>(s)]) {
      const std::vector<int>& Rk = rows_[static_cast<std::size_t>(k)];
      const int wk = sn_start_[static_cast<std::size_t>(k) + 1] - sn_start_[static_cast<std::size_t>(k)];
      for (std::size_t q = static_cast<std::size_t>(wk); q < Rk.size(); ++q) add(Rk[q]);
    }
    std::sort(below.begin(), below.end());
    R.insert(R.end(), below.begin(), below.end());
    require(static_cast<int>(below.size()) == cnt[static_cast<std::size_t>(e - 1)],
            "supernodal: inconsistent symbolic structure");
  }
  L_.assign(static_cast<std::size_t>(ns), MatX());
  Dp_.resize(n_);
}

bool SupernodalLDLT::factorize(const SpMat& Ain) {
  require(Ain.rows() == n_, "supernodal: factorize after analyze, same size");
  // the original values by position (the analysed lower pattern)
  SpMat L = Ain.triangularView<Eigen::Lower>();
  std::vector<double> val(origpos_.size(), 0.0);
  {
    std::size_t q = 0;
    for (int c = 0; c < L.outerSize(); ++c)
      for (SpMat::InnerIterator it(L, c); it; ++it) {
        require(q < origpos_.size() && origpos_[q].first == it.row() && origpos_[q].second == c,
                "supernodal: the pattern changed since analyze");
        val[q++] = it.value();
      }
  }
  const int ns = static_cast<int>(sn_start_.size()) - 1;
  std::vector<MatX> upd(static_cast<std::size_t>(ns));  // update matrices waiting for their parent
  std::vector<std::vector<int>> kids(static_cast<std::size_t>(ns));
  for (int s = 0; s < ns; ++s)
    if (sn_parent_[static_cast<std::size_t>(s)] >= 0) kids[static_cast<std::size_t>(sn_parent_[static_cast<std::size_t>(s)])].push_back(s);
  const int T = solver_threads();
  // the work of one supernode (one thread; the blocks of a large front may
  // fan out, on a FIXED partition: the arithmetic never depends on T)
  auto process = [&](int s, std::vector<int>& pos) {
    const int f = sn_start_[static_cast<std::size_t>(s)], w = sn_start_[static_cast<std::size_t>(s) + 1] - f;
    const std::vector<int>& R = rows_[static_cast<std::size_t>(s)];
    const int m = static_cast<int>(R.size());
    for (int i = 0; i < m; ++i) pos[static_cast<std::size_t>(R[static_cast<std::size_t>(i)])] = i;
    MatX F = MatX::Zero(m, m);  // lower triangle used
    for (int c = f; c < f + w; ++c)
      for (const auto& [r, p] : colmap_[static_cast<std::size_t>(c)])
        F(pos[static_cast<std::size_t>(r)], c - f) += val[static_cast<std::size_t>(p)];
    for (int k : kids[static_cast<std::size_t>(s)]) {  // extend-add, children in index order
      const std::vector<int>& Rk = rows_[static_cast<std::size_t>(k)];
      const int wk = sn_start_[static_cast<std::size_t>(k) + 1] - sn_start_[static_cast<std::size_t>(k)];
      const MatX& U = upd[static_cast<std::size_t>(k)];
      const int mk = static_cast<int>(Rk.size()) - wk;
      std::vector<int> map(static_cast<std::size_t>(mk));
      for (int i = 0; i < mk; ++i) map[static_cast<std::size_t>(i)] = pos[static_cast<std::size_t>(Rk[static_cast<std::size_t>(wk + i)])];
      for (int jj = 0; jj < mk; ++jj) {
        const int cj = map[static_cast<std::size_t>(jj)];
        for (int i = jj; i < mk; ++i) F(map[static_cast<std::size_t>(i)], cj) += U(i, jj);
      }
      upd[static_cast<std::size_t>(k)] = MatX();
    }
    // unpivoted LDL^t of the pivot block, BLOCKED: panels of kPanel columns
    // factored column by column (on the panel and the rows below it inside
    // the block), then the trailing block updated by one GEMM, L21 D1 L21'
    // (lower triangle) — a fixed partition, the same arithmetic always
    constexpr int kPanel = 48;
    for (int k0 = 0; k0 < w; k0 += kPanel) {
      const int kb = std::min(kPanel, w - k0), k1 = k0 + kb;
      for (int k = k0; k < k1; ++k) {
        const double d = F(k, k);
        if (d == 0.0 || !std::isfinite(d)) return false;
        Dp_(f + k) = d;
        if (k + 1 < w) {
          const VecX l = F.col(k).segment(k + 1, w - k - 1) / d;
          for (int jj = k + 1; jj < k1; ++jj) F.col(jj).segment(jj, w - jj).noalias() -= (d * l(jj - k - 1)) * l.segment(jj - k - 1, w - jj);
          F.col(k).segment(k + 1, w - k - 1) = l;
        }
        F(k, k) = 1.0;
      }
      const int nt = w - k1;  // the trailing block
      if (nt > 0) {
        MatX LD = F.block(k1, k0, nt, kb);
        for (int c = 0; c < kb; ++c) LD.col(c) *= Dp_(f + k0 + c);
        const MatX Lp = F.block(k1, k0, nt, kb);
        F.block(k1, k1, nt, nt).triangularView<Eigen::Lower>() -= Lp * LD.transpose();
      }
    }
    const int mb = m - w;
    if (mb > 0) {
      // W = F21 L11^-T (row blocks), L21 = W D^-1, U = F22 - W L21' (column
      // panels, lower part): blocks of kBlock, independent of each other
      constexpr int kBlock = 256;
      const MatX L11 = F.topLeftCorner(w, w);
      MatX W = F.bottomLeftCorner(mb, w);
      const int nb = (mb + kBlock - 1) / kBlock;
      parallel_blocks(nb, T, [&](int b) {
        const int r0 = b * kBlock, nr = std::min(kBlock, mb - r0);
        MatX Wb = W.middleRows(r0, nr);
        L11.triangularView<Eigen::UnitLower>().transpose().solveInPlace<Eigen::OnTheRight>(Wb);  // Wb L11' = F21 rows
        W.middleRows(r0, nr) = Wb;
      });
      MatX L21 = W;
      for (int k = 0; k < w; ++k) L21.col(k) /= Dp_(f + k);
      MatX U = F.bottomRightCorner(mb, mb);
      parallel_blocks(nb, T, [&](int b) {
        const int c0 = b * kBlock, nc = std::min(kBlock, mb - c0);
        U.block(c0, c0, mb - c0, nc).noalias() -= W.middleRows(c0, mb - c0) * L21.middleRows(c0, nc).transpose();
      });
      F.bottomLeftCorner(mb, w) = L21;
      if (sn_parent_[static_cast<std::size_t>(s)] >= 0) upd[static_cast<std::size_t>(s)] = std::move(U);
    }
    if (!F.allFinite()) return false;
    L_[static_cast<std::size_t>(s)] = F.leftCols(w);
    return true;
  };
  if (T <= 1) {
    std::vector<int> pos(static_cast<std::size_t>(n_), -1);
    for (int s = 0; s < ns; ++s)
      if (!process(s, pos)) return false;
    return true;
  }
  // the supernodal tree in parallel: a supernode is ready once its children
  // are done; each is computed whole by one worker (results independent of
  // the schedule: the same operations on the same data)
  std::vector<int> pending(static_cast<std::size_t>(ns), 0);
  std::vector<int> ready;
  for (int s = 0; s < ns; ++s) {
    pending[static_cast<std::size_t>(s)] = static_cast<int>(kids[static_cast<std::size_t>(s)].size());
    if (pending[static_cast<std::size_t>(s)] == 0) ready.push_back(s);
  }
  std::reverse(ready.begin(), ready.end());
  std::mutex mu;
  std::condition_variable cv;
  int done = 0;
  bool failed = false;
  auto worker = [&] {
    std::vector<int> pos(static_cast<std::size_t>(n_), -1);
    for (;;) {
      int s;
      {
        std::unique_lock<std::mutex> lk(mu);
        cv.wait(lk, [&] { return !ready.empty() || done == ns || failed; });
        if (done == ns || failed) return;
        s = ready.back();
        ready.pop_back();
      }
      const bool ok = process(s, pos);
      {
        std::lock_guard<std::mutex> lk(mu);
        if (!ok) failed = true;
        ++done;
        const int p = sn_parent_[static_cast<std::size_t>(s)];
        if (ok && p >= 0 && --pending[static_cast<std::size_t>(p)] == 0) ready.push_back(p);
      }
      cv.notify_all();
    }
  };
  std::vector<std::thread> th;
  for (int t = 0; t < T; ++t) th.emplace_back(worker);
  for (std::thread& t : th) t.join();
  return !failed;
}

VecX SupernodalLDLT::solve(const VecX& b) const {
  MatX B = b;
  return solve(B).col(0);
}

MatX SupernodalLDLT::solve(const MatX& B) const {
  require(B.rows() == n_, "supernodal: right-hand side size");
  MatX Y(n_, B.cols());
  for (Eigen::Index k = 0; k < n_; ++k) Y.row(k) = B.row(perm_(k));
  const int ns = static_cast<int>(sn_start_.size()) - 1;
  for (int s = 0; s < ns; ++s) {  // forward: L
    const int f = sn_start_[static_cast<std::size_t>(s)], w = sn_start_[static_cast<std::size_t>(s) + 1] - f;
    const MatX& Ls = L_[static_cast<std::size_t>(s)];
    const std::vector<int>& R = rows_[static_cast<std::size_t>(s)];
    MatX X = Y.middleRows(f, w);
    Ls.topRows(w).triangularView<Eigen::UnitLower>().solveInPlace(X);
    Y.middleRows(f, w) = X;
    const int mb = static_cast<int>(R.size()) - w;
    if (mb > 0) {
      const MatX T = Ls.bottomRows(mb) * X;
      for (int i = 0; i < mb; ++i) Y.row(R[static_cast<std::size_t>(w + i)]) -= T.row(i);
    }
  }
  for (Eigen::Index k = 0; k < n_; ++k) Y.row(k) /= Dp_(k);
  for (int s = ns - 1; s >= 0; --s) {  // backward: L'
    const int f = sn_start_[static_cast<std::size_t>(s)], w = sn_start_[static_cast<std::size_t>(s) + 1] - f;
    const MatX& Ls = L_[static_cast<std::size_t>(s)];
    const std::vector<int>& R = rows_[static_cast<std::size_t>(s)];
    const int mb = static_cast<int>(R.size()) - w;
    MatX X = Y.middleRows(f, w);
    if (mb > 0) {
      MatX G(mb, Y.cols());
      for (int i = 0; i < mb; ++i) G.row(i) = Y.row(R[static_cast<std::size_t>(w + i)]);
      X.noalias() -= Ls.bottomRows(mb).transpose() * G;
    }
    Ls.topRows(w).triangularView<Eigen::UnitLower>().transpose().solveInPlace(X);
    Y.middleRows(f, w) = X;
  }
  MatX Xo(n_, B.cols());
  for (Eigen::Index k = 0; k < n_; ++k) Xo.row(perm_(k)) = Y.row(k);
  return Xo;
}

}  // namespace keta

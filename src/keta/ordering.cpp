// Nested dissection and symbolic fill — see ordering.h.
#include "keta/ordering.h"

#include <algorithm>

namespace keta {

namespace {
struct Dissector {
  const SpMat& A;
  std::vector<int> mark, order;
  int stamp = 0;
  explicit Dissector(const SpMat& a) : A(a), mark(static_cast<std::size_t>(a.rows()), 0) {}
  // BFS levels from s within the vertices marked st
  std::vector<std::vector<int>> levels(int s, int st) {
    std::vector<std::vector<int>> L{{s}};
    mark[static_cast<std::size_t>(s)] = -st;
    for (;;) {
      std::vector<int> nx;
      for (int v : L.back())
        for (SpMat::InnerIterator it(A, v); it; ++it) {
          const int w = static_cast<int>(it.row());
          if (mark[static_cast<std::size_t>(w)] == st) {
            mark[static_cast<std::size_t>(w)] = -st;
            nx.push_back(w);
          }
        }
      if (nx.empty()) break;
      std::sort(nx.begin(), nx.end());
      L.push_back(std::move(nx));
    }
    for (const auto& l : L)
      for (int v : l) mark[static_cast<std::size_t>(v)] = st;
    return L;
  }
  void run(std::vector<int> nodes) {
    while (!nodes.empty()) {
      if (nodes.size() <= 64) {
        std::sort(nodes.begin(), nodes.end());
        order.insert(order.end(), nodes.begin(), nodes.end());
        return;
      }
      const int st = ++stamp;
      for (int v : nodes) mark[static_cast<std::size_t>(v)] = st;
      std::vector<std::vector<int>> L = levels(nodes.front(), st);
      std::size_t reached = 0;
      for (const auto& l : L) reached += l.size();
      if (reached < nodes.size()) {  // disconnected: this component, then the rest
        std::vector<int> comp, rest;
        for (const auto& l : L)
          for (int v : l) {
            comp.push_back(v);
            mark[static_cast<std::size_t>(v)] = 0;
          }
        for (int v : nodes)
          if (mark[static_cast<std::size_t>(v)] == st) rest.push_back(v);
        for (int v : rest) mark[static_cast<std::size_t>(v)] = 0;
        run(std::move(comp));
        nodes = std::move(rest);
        continue;
      }
      // pseudo-peripheral start: restart from the far end while the depth grows
      for (int k = 0; k < 4; ++k) {
        std::vector<std::vector<int>> L2 = levels(L.back().front(), st);
        if (L2.size() <= L.size()) break;
        L = std::move(L2);
      }
      for (int v : nodes) mark[static_cast<std::size_t>(v)] = 0;
      const std::size_t n = nodes.size();
      std::size_t below = 0, best = L.size(), best_size = n + 1;
      for (std::size_t i = 0; i < L.size(); ++i) {
        const std::size_t above = n - below - L[i].size();
        if (below >= 0.35 * n && above >= 0.35 * n && L[i].size() < best_size) {
          best = i;
          best_size = L[i].size();
        }
        below += L[i].size();
      }
      if (best == L.size()) {  // no balanced level: the median one
        below = 0;
        for (best = 0; best + 1 < L.size() && below + L[best].size() < n / 2; ++best) below += L[best].size();
      }
      std::vector<int> a, b;
      for (std::size_t i = 0; i < L.size(); ++i)
        if (i < best) a.insert(a.end(), L[i].begin(), L[i].end());
        else if (i > best) b.insert(b.end(), L[i].begin(), L[i].end());
      std::vector<int> sep = L[best];
      if (a.empty() || b.empty()) {  // a path-like remainder: no split
        std::sort(nodes.begin(), nodes.end());
        order.insert(order.end(), nodes.begin(), nodes.end());
        return;
      }
      run(std::move(a));
      run(std::move(b));
      std::sort(sep.begin(), sep.end());
      order.insert(order.end(), sep.begin(), sep.end());
      return;
    }
  }
};
}  // namespace

Eigen::VectorXi nested_dissection_order(const SpMat& A) {
  Dissector d(A);
  std::vector<int> all(static_cast<std::size_t>(A.rows()));
  for (std::size_t i = 0; i < all.size(); ++i) all[i] = static_cast<int>(i);
  d.run(std::move(all));
  Eigen::VectorXi o(A.rows());
  for (std::size_t k = 0; k < d.order.size(); ++k) o(static_cast<Eigen::Index>(k)) = d.order[k];
  return o;
}

SymbolicFill symbolic_fill(const SpMat& A, const Eigen::VectorXi& order) {
  // up-looking symbolic elimination (Eigen's analyzePattern_preordered):
  // row k of L reaches, from every i < k of A's pattern, up the elimination
  // tree until a vertex already tagged k
  const Eigen::Index n = A.rows();
  Eigen::VectorXi pos(n), parent = Eigen::VectorXi::Constant(n, -1), tag(n), cnt = Eigen::VectorXi::Zero(n);
  for (Eigen::Index k = 0; k < n; ++k) pos(order(k)) = static_cast<int>(k);
  for (Eigen::Index k = 0; k < n; ++k) {
    tag(k) = static_cast<int>(k);
    for (SpMat::InnerIterator it(A, order(k)); it; ++it) {
      int i = pos(static_cast<Eigen::Index>(it.row()));
      if (i >= k) continue;
      while (tag(i) != k) {
        if (parent(i) == -1) parent(i) = static_cast<int>(k);
        ++cnt(i);
        tag(i) = static_cast<int>(k);
        i = parent(i);
      }
    }
  }
  SymbolicFill f;
  for (Eigen::Index j = 0; j < n; ++j) {
    f.nnz += cnt(j);
    f.work += static_cast<double>(cnt(j)) * cnt(j);
  }
  return f;
}

}  // namespace keta

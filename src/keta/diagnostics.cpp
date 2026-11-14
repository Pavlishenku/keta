#include "keta/diagnostics.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <unordered_map>

#include "keta/eigsolve.h"

namespace keta {

std::vector<Diagnostic> diagnose(const Model& model, bool mechanism_probe, double tol) {
  std::vector<Diagnostic> out;
  auto add = [&](DiagnosticKind k, Severity sev, std::vector<Id> nodes, std::vector<Id> bars,
                 int dof, std::string msg, std::string sug) {
    out.push_back({k, sev, std::move(nodes), std::move(bars), dof, std::move(msg), std::move(sug)});
  };

  // coincident nodes: the pairs (i < j) closer than tol, reported in (i, j)
  // order. A uniform grid of cell tol: two such nodes lie in the same or in
  // adjacent cells, so each node meets only its 27-cell neighbourhood
  // (O(N) for spread-out nodes) instead of every other node (O(N^2)).
  std::vector<std::pair<std::size_t, std::size_t>> close;
  {
    const std::size_t N = model.nodes.size();
    bool gridable = tol > 0 && std::isfinite(tol);
    for (const Node& nd : model.nodes)
      for (int c = 0; c < 3 && gridable; ++c)
        gridable = std::isfinite(nd.x(c)) && std::fabs(nd.x(c) / tol) < 1e15;
    if (gridable) {
      using Key = std::array<long long, 3>;
      struct KeyHash {
        std::size_t operator()(const Key& k) const {
          std::size_t h = 1469598103934665603ULL;
          for (long long v : k) h = (h ^ static_cast<std::size_t>(v)) * 1099511628211ULL;
          return h;
        }
      };
      std::unordered_map<Key, std::vector<std::size_t>, KeyHash> cells;
      std::vector<Key> key(N);
      for (std::size_t i = 0; i < N; ++i) {
        for (int c = 0; c < 3; ++c)
          key[i][static_cast<std::size_t>(c)] = static_cast<long long>(std::floor(model.nodes[i].x(c) / tol));
        cells[key[i]].push_back(i);
      }
      for (std::size_t i = 0; i < N; ++i)
        for (long long dx = -1; dx <= 1; ++dx)
          for (long long dy = -1; dy <= 1; ++dy)
            for (long long dz = -1; dz <= 1; ++dz) {
              auto it = cells.find({key[i][0] + dx, key[i][1] + dy, key[i][2] + dz});
              if (it == cells.end()) continue;
              for (std::size_t j : it->second)
                if (j > i && (model.nodes[i].x - model.nodes[j].x).norm() < tol) close.emplace_back(i, j);
            }
      std::sort(close.begin(), close.end());
    } else {
      // non-finite coordinates or tolerance, or a grid too fine for the
      // extent: the direct pair test
      for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = i + 1; j < N; ++j)
          if ((model.nodes[i].x - model.nodes[j].x).norm() < tol) close.emplace_back(i, j);
    }
  }
  for (const auto& [i, j] : close)
    add(DiagnosticKind::CoincidentNodes, Severity::Warning,
        {model.nodes[i].id, model.nodes[j].id}, {}, -1,
        "nodes " + std::to_string(model.nodes[i].id) + " and " +
            std::to_string(model.nodes[j].id) + " coincide within tolerance",
        "weld them or increase the spacing — an unwelded seam carries no force");

  // zero-length / degenerate bars, degenerate properties
  for (const Bar& b : model.bars) {
    double L = 0;
    bool ok = true;
    try {
      L = model.bar_chord_length(b);
    } catch (const Error&) {
      ok = false;
    }
    if (!ok) continue;
    if (L < tol)
      add(DiagnosticKind::ZeroLengthBar, Severity::Error, {b.ni, b.nj}, {b.id}, -1,
          "bar " + std::to_string(b.id) + " has near-zero length",
          "remove the bar or move its nodes apart");
    try {
      const SectionProps& p = model.section(b.section).props;
      const Material& mat = model.material(b.material);
      if (!(p.A > 0) || !(mat.E > 0))
        add(DiagnosticKind::DegenerateProperty, Severity::Error, {}, {b.id}, -1,
            "bar " + std::to_string(b.id) + " has zero axial rigidity",
            "give the section a positive area and the material a positive modulus");
    } catch (const Error&) {
    }
  }

  // isolated nodes and floating islands (union-find over bars)
  std::vector<int> parent(model.nodes.size());
  for (std::size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
  std::function<int(int)> find = [&](int a) {
    while (parent[static_cast<std::size_t>(a)] != a) a = parent[static_cast<std::size_t>(a)] =
                                          parent[static_cast<std::size_t>(parent[static_cast<std::size_t>(a)])];
    return a;
  };
  std::vector<bool> touched(model.nodes.size(), false);
  for (const Bar& b : model.bars) {
    int i = model.node_index(b.ni), j = model.node_index(b.nj);
    touched[static_cast<std::size_t>(i)] = touched[static_cast<std::size_t>(j)] = true;
    parent[static_cast<std::size_t>(find(i))] = find(j);
  }
  for (std::size_t i = 0; i < model.nodes.size(); ++i)
    if (!touched[i])
      add(DiagnosticKind::IsolatedNode, Severity::Warning, {model.nodes[i].id}, {}, -1,
          "node " + std::to_string(model.nodes[i].id) + " is connected to nothing",
          "remove it or connect a bar");
  // islands: a component with no support
  std::unordered_map<int, bool> comp_supported;
  for (std::size_t i = 0; i < model.nodes.size(); ++i)
    if (touched[i]) comp_supported.emplace(find(static_cast<int>(i)), false);
  for (const Support& s : model.supports) {
    int r = find(model.node_index(s.node));
    if (comp_supported.count(r)) comp_supported[r] = true;
  }
  for (const auto& [root, supported] : comp_supported)
    if (!supported) {
      std::vector<Id> nodes;
      for (std::size_t i = 0; i < model.nodes.size(); ++i)
        if (touched[i] && find(static_cast<int>(i)) == root) nodes.push_back(model.nodes[i].id);
      add(DiagnosticKind::FloatingIsland, Severity::Error, std::move(nodes), {}, -1,
          "a connected group of bars has no support at all",
          "support the island or connect it to the supported structure");
    }

  // mechanism probe: factorization verdict first, localization only on failure
  if (mechanism_probe) {
    try {
      Assembled asm_ = assemble_elements(model);
      if (asm_.dofs.n_free > 0) {
        SpMat K = assemble_stiffness(model, asm_);
        FactorCache fac;
        bool singular = false;
        int n_neg = 0;
        try {
          fac.factorize(K);
          // Signed pivot classification: |d| <= tol*max is a NULL pivot (an
          // effective mechanism — a rounding-level null may come out slightly
          // negative and must not masquerade as indefiniteness); only
          // d < -tol*max counts as genuinely negative.
          double dmax = fac.max_abs_pivot();
          if (fac.min_abs_pivot() <= 1e-10 * dmax) singular = true;
          if (!singular && fac.inertia().n_neg > 0) n_neg = fac.inertia().n_neg;
        } catch (const Error&) {
          singular = true;
        }
        if (n_neg > 0) {
          add(DiagnosticKind::IndefiniteStiffness, Severity::Error, {}, {}, -1,
              "the stiffness matrix is indefinite (" + std::to_string(n_neg) +
                  " negative pivots)",
              "check spring signs and geometric assumptions");
        } else if (singular) {
          // localize with the lowest pencil mode of (K + shift, I)
          Diagnostic d{DiagnosticKind::Mechanism, Severity::Error, {}, {}, -1,
                       "the model contains a mechanism (singular stiffness)",
                       "restrain the reported DOF or connect the reported nodes"};
          try {
            SpMat I(asm_.dofs.n_free, asm_.dofs.n_free);
            I.setIdentity();
            SpMat Ks = (K + 1e-6 * I).eval();
            auto pairs = dense_pencil(Ks, I);
            if (!pairs.empty()) {
              const VecX& mode = pairs.front().x;  // largest mu of I vs K = softest direction
              // find the largest components
              double best = 0;
              int best_dof = -1;
              for (Eigen::Index i = 0; i < mode.size(); ++i)
                if (std::fabs(mode(i)) > best) {
                  best = std::fabs(mode(i));
                  best_dof = static_cast<int>(i);
                }
              for (int ni = 0; ni < asm_.dofs.n_nodes && best_dof >= 0; ++ni)
                for (int c = 0; c < 6; ++c)
                  if (asm_.dofs.free_ix[static_cast<std::size_t>(asm_.dofs.full(ni, c))] ==
                      best_dof) {
                    d.nodes.push_back(model.nodes[static_cast<std::size_t>(ni)].id);
                    d.dof = c;
                  }
            }
          } catch (const Error&) {
            d.message += " (localization probe failed)";
          }
          out.push_back(std::move(d));
        }
      }
    } catch (const Error& e) {
      add(DiagnosticKind::Mechanism, Severity::Warning, {}, {}, -1,
          std::string("mechanism probe could not run: ") + e.what(),
          "fix the reported modeling errors first");
    }
  }
  return out;
}

}  // namespace keta

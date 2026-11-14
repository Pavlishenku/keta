#include "keta/combination.h"

#include <cmath>
#include <map>
#include <sstream>

namespace keta {

namespace {
struct Choice {
  double contrib = 0;  // signed adverse measure s * f * e
  double factor = 0;
  int action = -1;
};
// the admissible factors of an action in a role
std::vector<double> factors(const CombAction& a, const CombinationRules& r, LimitState ls, int variant, bool lead) {
  const bool uls = ls == LimitState::ULS;
  if (a.kind == ActionKind::Permanent) {
    const double sup = a.gamma_sup >= 0 ? a.gamma_sup : r.gamma_G_sup, inf = a.gamma_inf >= 0 ? a.gamma_inf : r.gamma_G_inf;
    if (!uls) return {1.0};
    if (variant == 2) return {r.xi * sup, inf};  // 6.10b
    return {sup, inf};
  }
  if (a.kind == ActionKind::Settlement) {
    const double sup = a.gamma_sup >= 0 ? a.gamma_sup : r.gamma_set_sup, inf = a.gamma_inf >= 0 ? a.gamma_inf : r.gamma_set_inf;
    if (!uls) return {1.0, 0.0};
    return {sup, inf};
  }
  // variable
  switch (ls) {
    case LimitState::ULS:
      if (variant == 1) return {a.gamma_sup * a.psi0, 0.0};  // 6.10a: no leading action
      return {lead ? a.gamma_sup : a.gamma_sup * a.psi0, 0.0};
    case LimitState::SLS_characteristic: return {lead ? 1.0 : a.psi0, 0.0};
    case LimitState::SLS_frequent: return {lead ? a.psi1 : a.psi2, 0.0};
    case LimitState::SLS_quasi_permanent: return {a.psi2, 0.0};
  }
  return {0.0};
}

CombinationExtreme solve(const std::vector<double>& e, const std::vector<CombAction>& A, const CombinationRules& r,
                         LimitState ls, bool maximise) {
  require(e.size() == A.size(), "combination: one effect per action");
  const double s = maximise ? 1.0 : -1.0;
  // groups: exclusive tags, singletons otherwise
  std::map<std::string, std::vector<int>> tagged;
  std::vector<std::vector<int>> groups;
  for (std::size_t a = 0; a < A.size(); ++a) {
    require(std::isfinite(e[a]), "combination: non-finite effect");
    if (A[a].exclusive.empty())
      groups.push_back({static_cast<int>(a)});
    else
      tagged[A[a].exclusive].push_back(static_cast<int>(a));
  }
  for (auto& [t, g] : tagged) {
    bool perm = A[static_cast<std::size_t>(g.front())].kind == ActionKind::Permanent;
    for (int a : g)
      require((A[static_cast<std::size_t>(a)].kind == ActionKind::Permanent) == perm,
              "combination: an exclusive group mixes permanent and other actions");
    groups.push_back(g);
  }
  // variants: ULS 6.10ab -> 1 (6.10a) and 2 (6.10b); otherwise 0
  std::vector<int> variants{0};
  if (ls == LimitState::ULS && r.eq610ab) variants = {1, 2};
  CombinationExtreme best;
  double best_score = -1e300;
  for (int variant : variants) {
    const bool leads = !(ls == LimitState::SLS_quasi_permanent || (ls == LimitState::ULS && variant == 1));
    std::vector<int> lead_options{-1};
    if (leads)
      for (std::size_t a = 0; a < A.size(); ++a)
        if (A[a].kind == ActionKind::Variable && A[a].may_lead) lead_options.push_back(static_cast<int>(a));
    for (int L : lead_options) {
      double score = 0;
      std::vector<Choice> picks;
      for (const std::vector<int>& g : groups) {
        const bool perm = A[static_cast<std::size_t>(g.front())].kind == ActionKind::Permanent;
        bool hasL = false;
        for (int a : g) hasL = hasL || a == L;
        Choice c;
        c.contrib = perm ? -1e300 : 0.0;  // a variable / settlement group may stay empty
        for (int a : g) {
          if (hasL && a != L) continue;
          for (double f : factors(A[static_cast<std::size_t>(a)], r, ls, variant, a == L)) {
            const double v = s * f * e[static_cast<std::size_t>(a)];
            if (v > c.contrib) c = {v, f, a};
          }
        }
        score += c.contrib;
        if (c.action >= 0 && c.factor != 0) picks.push_back(c);
      }
      if (score > best_score) {
        best_score = score;
        best.value = s * score;
        best.terms.clear();
        for (const Choice& c : picks) best.terms.push_back({A[static_cast<std::size_t>(c.action)].case_id, c.factor});
        std::ostringstream lab;
        lab << (ls == LimitState::ULS ? (variant == 1 ? "ULS 6.10a" : variant == 2 ? "ULS 6.10b" : "ULS 6.10")
                : ls == LimitState::SLS_characteristic ? "SLS characteristic"
                : ls == LimitState::SLS_frequent       ? "SLS frequent"
                                                       : "SLS quasi-permanent");
        if (L >= 0) lab << ", leading: " << (A[static_cast<std::size_t>(L)].name.empty() ? std::to_string(A[static_cast<std::size_t>(L)].case_id) : A[static_cast<std::size_t>(L)].name);
        best.label = lab.str();
      }
    }
  }
  return best;
}
}  // namespace

CombinationExtreme combination_extreme(const std::vector<double>& effects, const std::vector<CombAction>& actions,
                                       const CombinationRules& rules, LimitState ls, bool maximise) {
  return solve(effects, actions, rules, ls, maximise);
}

CombinationEnvelope combination_envelope(const std::vector<const CaseResult*>& R, const std::vector<CombAction>& A,
                                         const CombinationRules& rules, LimitState ls) {
  require(!R.empty() && R.size() == A.size(), "combination envelope: one result per action");
  for (const CaseResult* r : R)
    require(r && r->node_disp.size() == R[0]->node_disp.size() && r->reactions.size() == R[0]->reactions.size() &&
                r->bars.size() == R[0]->bars.size() && r->shells.size() == R[0]->shells.size(),
            "combination envelope: results of different shapes");
  CombinationEnvelope out;
  out.lo.value = *R[0];
  out.hi.value = *R[0];
  std::map<std::string, unsigned> index;
  auto intern = [&](const CombinationExtreme& c) {
    std::ostringstream k;
    k << c.label;
    for (const CombTerm& t : c.terms) k << '|' << t.case_id << ':' << t.factor;
    auto it = index.find(k.str());
    if (it != index.end()) return it->second;
    const unsigned i = static_cast<unsigned>(out.combinations.size());
    out.combinations.push_back(c);
    out.combinations.back().value = 0;
    index[k.str()] = i;
    return i;
  };
  std::vector<double> e(R.size());
  auto field = [&](auto&& get, auto&& set, unsigned* slo, unsigned* shi) {
    for (std::size_t a = 0; a < R.size(); ++a) e[a] = get(*R[a]);
    const CombinationExtreme hi = solve(e, A, rules, ls, true), lo = solve(e, A, rules, ls, false);
    set(out.hi.value, hi.value);
    set(out.lo.value, lo.value);
    *shi = intern(hi);
    *slo = intern(lo);
  };
  const CaseResult& a = *R[0];
  out.lo.disp_src.assign(a.node_disp.size(), {});
  out.hi.disp_src.assign(a.node_disp.size(), {});
  for (std::size_t ni = 0; ni < a.node_disp.size(); ++ni)
    for (int c = 0; c < 6; ++c)
      field([&](const CaseResult& r) { return r.node_disp[ni](c); }, [&](CaseResult& r, double v) { r.node_disp[ni](c) = v; },
            &out.lo.disp_src[ni][static_cast<std::size_t>(c)], &out.hi.disp_src[ni][static_cast<std::size_t>(c)]);
  out.lo.reaction_src.assign(a.reactions.size(), {});
  out.hi.reaction_src.assign(a.reactions.size(), {});
  for (std::size_t ri = 0; ri < a.reactions.size(); ++ri)
    for (int c = 0; c < 6; ++c)
      field([&](const CaseResult& r) { return r.reactions[ri].R(c); }, [&](CaseResult& r, double v) { r.reactions[ri].R(c) = v; },
            &out.lo.reaction_src[ri][static_cast<std::size_t>(c)], &out.hi.reaction_src[ri][static_cast<std::size_t>(c)]);
  out.lo.bar_src.resize(a.bars.size());
  out.hi.bar_src.resize(a.bars.size());
  auto comp = [](const Station& s, int c) {
    const double v[6] = {s.N, s.Vy, s.Vz, s.T, s.My, s.Mz};
    return v[c];
  };
  auto set_comp = [](Station& s, int c, double x) {
    double* p[6] = {&s.N, &s.Vy, &s.Vz, &s.T, &s.My, &s.Mz};
    *p[c] = x;
  };
  for (std::size_t b = 0; b < a.bars.size(); ++b) {
    out.lo.bar_src[b].stations.resize(a.bars[b].stations.size());
    out.hi.bar_src[b].stations.resize(a.bars[b].stations.size());
    for (int c = 0; c < 12; ++c)
      field([&](const CaseResult& r) { return r.bars[b].end_forces_local(c); },
            [&](CaseResult& r, double v) { r.bars[b].end_forces_local(c) = v; }, &out.lo.bar_src[b].end[static_cast<std::size_t>(c)],
            &out.hi.bar_src[b].end[static_cast<std::size_t>(c)]);
    for (std::size_t k = 0; k < a.bars[b].stations.size(); ++k) {
      for (int c = 0; c < 6; ++c)
        field([&](const CaseResult& r) { return comp(r.bars[b].stations[k], c); },
              [&](CaseResult& r, double v) { set_comp(r.bars[b].stations[k], c, v); },
              &out.lo.bar_src[b].stations[k].forces[static_cast<std::size_t>(c)],
              &out.hi.bar_src[b].stations[k].forces[static_cast<std::size_t>(c)]);
      for (int c = 0; c < 3; ++c) {
        field([&](const CaseResult& r) { return r.bars[b].stations[k].defl(c); },
              [&](CaseResult& r, double v) { r.bars[b].stations[k].defl(c) = v; },
              &out.lo.bar_src[b].stations[k].defl[static_cast<std::size_t>(c)],
              &out.hi.bar_src[b].stations[k].defl[static_cast<std::size_t>(c)]);
        field([&](const CaseResult& r) { return r.bars[b].stations[k].rot(c); },
              [&](CaseResult& r, double v) { r.bars[b].stations[k].rot(c) = v; },
              &out.lo.bar_src[b].stations[k].rot[static_cast<std::size_t>(c)],
              &out.hi.bar_src[b].stations[k].rot[static_cast<std::size_t>(c)]);
      }
    }
  }
  out.lo.shell_src.assign(a.shells.size(), ShellSrc{});
  out.hi.shell_src.assign(a.shells.size(), ShellSrc{});
  for (std::size_t sh = 0; sh < a.shells.size(); ++sh) {
    auto get = [](const ShellResultants& v, int c) {
      const double x[12] = {v.nx, v.ny, v.nxy, v.mx, v.my, v.mxy, v.qx, v.qy, v.m1, v.m2, v.n1, v.n2};
      return x[c];
    };
    auto put = [](ShellResultants& v, int c, double x) {
      double* p[12] = {&v.nx, &v.ny, &v.nxy, &v.mx, &v.my, &v.mxy, &v.qx, &v.qy, &v.m1, &v.m2, &v.n1, &v.n2};
      *p[c] = x;
    };
    for (int c = 0; c < 12; ++c)
      field([&](const CaseResult& r) { return get(r.shells[sh].r, c); }, [&](CaseResult& r, double v) { put(r.shells[sh].r, c, v); },
            &out.lo.shell_src[sh][static_cast<std::size_t>(c)], &out.hi.shell_src[sh][static_cast<std::size_t>(c)]);
    out.lo.value.shells[sh].r.m_theta = out.lo.value.shells[sh].r.n_theta = 0;
    out.hi.value.shells[sh].r.m_theta = out.hi.value.shells[sh].r.n_theta = 0;
  }
  return out;
}

}  // namespace keta

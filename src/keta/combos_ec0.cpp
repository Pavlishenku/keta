#include "keta/combos_ec0.h"

namespace keta {

VariableAction traffic_gr1a(Id id) { return {id, "gr1a", 0.75, 0.75, 0.0, 1.35}; }
VariableAction wind_action(Id id) { return {id, "wind", 0.6, 0.2, 0.0, 1.5}; }
VariableAction thermal_action(Id id) { return {id, "thermal", 0.6, 0.6, 0.5, 1.5}; }
VariableAction snow_action(Id id) { return {id, "snow", 0.8, 0.0, 0.0, 1.5}; }

namespace {

std::vector<Id> generate(Model& model, const std::vector<PermanentAction>& G,
                         const std::vector<VariableAction>& Q, const std::string& tag,
                         double gq_scale, bool lead_full, double psi_lead(const VariableAction&),
                         double psi_others(const VariableAction&), bool gamma_pairs,
                         const Ec0Options* uls) {
  std::vector<Id> out;
  auto emit = [&](const std::string& name, double gG, double xi, int lead) {
    std::vector<std::pair<Id, double>> terms;
    for (const auto& g : G) terms.emplace_back(g.case_id, gG * xi);
    for (std::size_t q = 0; q < Q.size(); ++q) {
      const VariableAction& a = Q[q];
      double gq = gq_scale > 0 ? a.gamma_q * gq_scale : 1.0;
      double factor;
      if (static_cast<int>(q) == lead)
        factor = gq * (lead_full ? 1.0 : psi_lead(a));
      else
        factor = gq * psi_others(a);
      if (factor != 0) terms.emplace_back(a.case_id, factor);
    }
    out.push_back(model.add_combo(name, terms));
  };

  int n_lead = Q.empty() ? 0 : static_cast<int>(Q.size());
  // xi (6.10b) reduces the unfavourable permanent actions only: the
  // favourable ones stay at gamma_inf (EN 1990 Table A2.4(B), note to 6.10b)
  auto gammas = [&](int lead, const std::string& base, double xi) {
    if (gamma_pairs) {
      emit(base + " Gsup", G.empty() ? 1.35 : G.front().gamma_sup, xi, lead);
      emit(base + " Ginf", G.empty() ? 1.0 : G.front().gamma_inf, 1.0, lead);
    } else {
      emit(base, 1.0, 1.0, lead);
    }
  };

  if (n_lead == 0) {
    gammas(-1, tag + " (G only)", 1.0);
  } else {
    // 6.10a has no leading action (every variable action at psi0): one pair,
    // permanent actions unfavourable then favourable
    if (uls && uls->use_610ab) gammas(-2, tag + " 6.10a", 1.0);
    for (int lead = 0; lead < n_lead; ++lead) {
      std::string base = tag + " lead " + Q[static_cast<std::size_t>(lead)].name;
      if (uls && uls->use_610ab) {
        // 6.10b: the lead in full, xi on the unfavourable permanent actions
        gammas(lead, base + " 6.10b", uls->xi);
      } else {
        gammas(lead, base, 1.0);
      }
    }
  }
  return out;
}

}  // namespace

std::vector<Id> generate_uls(Model& model, const std::vector<PermanentAction>& G,
                             const std::vector<VariableAction>& Q, const Ec0Options& opt) {
  return generate(
      model, G, Q, "ULS", 1.0, true, [](const VariableAction& a) { return a.psi0; },
      [](const VariableAction& a) { return a.psi0; }, true, &opt);
}
std::vector<Id> generate_sls_characteristic(Model& model, const std::vector<PermanentAction>& G,
                                            const std::vector<VariableAction>& Q) {
  return generate(
      model, G, Q, "SLS char", -1.0, true, [](const VariableAction& a) { return a.psi0; },
      [](const VariableAction& a) { return a.psi0; }, false, nullptr);
}
std::vector<Id> generate_sls_frequent(Model& model, const std::vector<PermanentAction>& G,
                                      const std::vector<VariableAction>& Q) {
  return generate(
      model, G, Q, "SLS freq", -1.0, false, [](const VariableAction& a) { return a.psi1; },
      [](const VariableAction& a) { return a.psi2; }, false, nullptr);
}
std::vector<Id> generate_sls_quasi_permanent(Model& model, const std::vector<PermanentAction>& G,
                                             const std::vector<VariableAction>& Q) {
  std::vector<Id> out;
  std::vector<std::pair<Id, double>> terms;
  for (const auto& g : G) terms.emplace_back(g.case_id, 1.0);
  for (const auto& a : Q)
    if (a.psi2 != 0) terms.emplace_back(a.case_id, a.psi2);
  out.push_back(model.add_combo("SLS qp", terms));
  return out;
}

}  // namespace keta

// EN 1990 combination generator (French National Annex defaults):
// ULS STR eq. 6.10 (and the 6.10a/6.10b pair), SLS characteristic / frequent /
// quasi-permanent — every variable action taking its turn as leading action,
// the others reduced by psi0/psi1/psi2; permanent actions at gamma_sup/inf.
#pragma once
#include "keta/model.h"

namespace keta {

struct VariableAction {
  Id case_id = 0;
  std::string name;
  double psi0 = 0.7, psi1 = 0.5, psi2 = 0.3;
  double gamma_q = 1.5;   // traffic per French NA: 1.35
};
struct PermanentAction {
  Id case_id = 0;
  double gamma_sup = 1.35, gamma_inf = 1.0;
};

struct Ec0Options {
  bool use_610ab = false;   // false = single 6.10; true = 6.10a + 6.10b pair
  double xi = 0.85;         // 6.10b reduction on permanent actions
};

// Common psi sets (EN 1990 A2, French NA):
VariableAction traffic_gr1a(Id case_id);   // TS/UDL group: psi0 varies; use 0.75/0.75/0
VariableAction wind_action(Id case_id);    // 0.6 / 0.2 / 0
VariableAction thermal_action(Id case_id); // 0.6 / 0.6 / 0.5
VariableAction snow_action(Id case_id);    // 0.8 / 0 / 0 (snow on a bridge, Table A2.1)

// Generate combinations and append them to the model. Returns the created ids.
std::vector<Id> generate_uls(Model& model, const std::vector<PermanentAction>& G,
                             const std::vector<VariableAction>& Q, const Ec0Options& opt = {});
std::vector<Id> generate_sls_characteristic(Model& model, const std::vector<PermanentAction>& G,
                                            const std::vector<VariableAction>& Q);
std::vector<Id> generate_sls_frequent(Model& model, const std::vector<PermanentAction>& G,
                                      const std::vector<VariableAction>& Q);
std::vector<Id> generate_sls_quasi_permanent(Model& model, const std::vector<PermanentAction>& G,
                                             const std::vector<VariableAction>& Q);

}  // namespace keta

// Combination Family and Envelope components.
#include <algorithm>
#include <cctype>
#include <set>

#include "keta/param/components.h"

namespace keta::param {

namespace {

// Plain defaults applied when the matching pin is unwired — ordinary inputs,
// no design code is built in.
constexpr double kPermUnfav = 1.35;
constexpr double kPermFav = 1.0;
constexpr double kVarFactor = 1.5;
constexpr double kPsi = 0.7;

// Integer ids of a pin, duplicates dropped (first occurrence kept).
std::vector<std::int64_t> ids_dedup(const Tree& t) {
  std::vector<std::int64_t> out;
  std::set<std::int64_t> seen;
  for (const Value& v : t.non_nulls())
    if (auto i = v.as_int(); i && seen.insert(*i).second) out.push_back(*i);
  return out;
}

// Station quantity selector: a case-insensitive name or an int 0..5, in the
// Station force order N, Vy, Vz, Mt, My, Mz.
int station_quantity(const Value& v) {
  if (const std::string* s = v.as_text()) {
    std::string k;
    for (char c : *s)
      if (!std::isspace(static_cast<unsigned char>(c))) k.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    const char* names[] = {"n", "vy", "vz", "mt", "my", "mz"};
    for (int i = 0; i < 6; ++i)
      if (k == names[i]) return i;
    throw ComponentError("unknown station quantity \"" + *s + "\" (expected N, Vy, Vz, Mt, My, Mz or an int 0 to 5)");
  }
  if (auto i = v.as_int()) {
    if (*i >= 0 && *i <= 5) return static_cast<int>(*i);
    throw ComponentError("station quantity index " + std::to_string(*i) + " out of range (expected 0 to 5)");
  }
  throw ComponentError("station quantity must be a name (N, Vy, Vz, Mt, My, Mz) or an int 0 to 5, got a " +
                       v.type_tag());
}

double station_force(const Station& st, int k) {
  switch (k) {
    case 0: return st.N;
    case 1: return st.Vy;
    case 2: return st.Vz;
    case 3: return st.T;
    case 4: return st.My;
    default: return st.Mz;
  }
}

// The dominant-source label: its own name, else C{id} / LC{id}.
std::string source_label(const EnvSource& s) {
  if (!s.name.empty()) return s.name;
  return (s.is_combo ? "C" : "LC") + std::to_string(s.id);
}

}  // namespace

// ============================ Combination Family ============================

ComponentInfo ComboFamilyNode::info() const {
  return {"keta.combo_family", "Combination Family", "Analyze", "Generatively expand a family of load combinations"};
}
std::vector<PinSpec> ComboFamilyNode::inputs() const {
  return {PinSpec::list("P", "int", "Permanent case ids"),
          PinSpec::list("V", "int", "Variable case ids").opt(),
          PinSpec::item("GS", "number", "Permanent factor, unfavourable").opt(),
          PinSpec::item("GI", "number", "Permanent factor, favourable").opt(),
          PinSpec::item("GV", "number", "Variable factor").opt(),
          PinSpec::list("PS", "number", "Accompaniment factor psi per variable case (paired with V; last repeats)").opt(),
          PinSpec::item("N", "text", "Name prefix").opt()};
}
std::vector<PinSpec> ComboFamilyNode::outputs() const {
  return {PinSpec::list("CB", "combo", "The expanded combination family")};
}
Outputs<Value> ComboFamilyNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::vector<std::int64_t> perm = ids_dedup(input_or_empty(in, 0));
  if (perm.empty()) throw ComponentError::missing("P");
  const std::vector<std::int64_t> vars = ids_dedup(input_or_empty(in, 1));
  const double gs = first_number(input_or_empty(in, 2)).value_or(kPermUnfav);
  const double gi = first_number(input_or_empty(in, 3)).value_or(kPermFav);
  const double gv = first_number(input_or_empty(in, 4)).value_or(kVarFactor);
  std::vector<double> psi;
  for (const Value& v : input_or_empty(in, 5).non_nulls())
    if (auto x = v.as_number()) psi.push_back(*x);
  auto psi_for = [&](std::size_t j) { return psi.empty() ? kPsi : (j < psi.size() ? psi[j] : psi.back()); };
  const std::string prefix = first_text(input_or_empty(in, 6)).value_or("ULS");

  Tree cb;
  cb.ensure_path(Path::single(0));
  int n = 0;
  auto emit = [&](std::string label, std::vector<std::pair<double, Id>> terms) {
    ComboSpec c;
    c.name = prefix + std::to_string(++n) + ": " + label;
    c.terms = std::move(terms);
    cb.append(Path::single(0), Value::of(std::move(c)));
  };
  const std::pair<double, const char*> signs[] = {{gs, "+"}, {gi, "-"}};
  // Two permanent-only combinations: unfavourable (G+) then favourable (G-).
  for (const auto& [g, tag] : signs) {
    std::vector<std::pair<double, Id>> terms;
    for (std::int64_t p : perm) terms.emplace_back(g, p);
    emit(std::string("G") + tag, std::move(terms));
  }
  // Each leading variable case: G+ then G-, leader at GV, the others at GV*psi.
  for (std::size_t i = 0; i < vars.size(); ++i)
    for (const auto& [g, tag] : signs) {
      std::vector<std::pair<double, Id>> terms;
      for (std::int64_t p : perm) terms.emplace_back(g, p);
      terms.emplace_back(gv, vars[i]);
      for (std::size_t j = 0; j < vars.size(); ++j)
        if (j != i) terms.emplace_back(gv * psi_for(j), vars[j]);
      emit(std::string("G") + tag + " lead LC" + std::to_string(vars[i]), std::move(terms));
    }
  return {cb};
}

// ============================ Envelope ============================

ComponentInfo EnvelopeNode::info() const {
  return {"keta.envelope", "Envelope", "Analyze", "Min/max envelope with the dominant source behind each extremum"};
}
std::vector<PinSpec> EnvelopeNode::inputs() const {
  return {PinSpec::item("R", "results", "Analyzed model"),
          PinSpec::item("Q", "text", "Station quantity for MIN/MAX/DMIN/DMAX (N, Vy, Vz, Mt, My, Mz)").opt(),
          PinSpec::item("N", "text", "Envelope name").opt()};
}
std::vector<PinSpec> EnvelopeNode::outputs() const {
  return {PinSpec::item("E", "envelope", "The envelope"),
          PinSpec::tree("MIN", "number", "Min of Q: one branch per bar, one item per station"),
          PinSpec::tree("MAX", "number", "Max of Q, same shape"),
          PinSpec::tree("DMIN", "text", "Dominant source behind each MIN item, same shape"),
          PinSpec::tree("DMAX", "text", "Dominant source behind each MAX item, same shape")};
}
Outputs<Value> EnvelopeNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  ResultsRef analyzed;
  for (const Value& v : input_or_empty(in, 0).non_nulls())
    if (const auto* r = v.get<ResultsRef>(); r && *r) {
      analyzed = *r;
      break;
    }
  if (!analyzed) throw ComponentError::missing("R");
  int q = 4;  // My
  {
    const auto vals = input_or_empty(in, 1).non_nulls();
    if (!vals.empty()) q = station_quantity(vals.front());
  }

  // Over the combinations when there are any, else over the load cases.
  const StaticResults& res = analyzed->results;
  const bool is_combo = !res.combos.empty();
  const std::vector<CaseResult>& chosen = is_combo ? res.combos : res.cases;
  auto out = std::make_shared<EnvelopeModel>();
  out->model = analyzed->model;
  out->name = first_text(input_or_empty(in, 2));
  std::vector<const CaseResult*> sources;
  for (const CaseResult& c : chosen) {
    sources.push_back(&c);
    out->sources.push_back({c.id, c.name, is_combo});
  }
  try {
    out->env = envelope(sources);
  } catch (const std::exception& e) {
    throw ComponentError(std::string("envelope failed: ") + e.what());
  }
  out->analyzed = analyzed;  // the section-stress envelope, on demand

  Tree mn, mx, dmn, dmx;
  const auto& lo = out->env.lo;
  const auto& hi = out->env.hi;
  for (std::size_t b = 0; b < lo.value.bars.size(); ++b) {
    const Path p = Path::single(static_cast<int>(b));
    for (std::size_t s = 0; s < lo.value.bars[b].stations.size(); ++s) {
      mn.append(p, Value::number(station_force(lo.value.bars[b].stations[s], q)));
      mx.append(p, Value::number(station_force(hi.value.bars[b].stations[s], q)));
      dmn.append(p, Value::text(source_label(out->sources[lo.bar_src[b].stations[s].forces[static_cast<std::size_t>(q)]])));
      dmx.append(p, Value::text(source_label(out->sources[hi.bar_src[b].stations[s].forces[static_cast<std::size_t>(q)]])));
    }
  }
  return {single(Value::of(EnvelopeRef(std::move(out)))), mn, mx, dmn, dmx};
}

}  // namespace keta::param

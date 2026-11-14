#include "keta/etude.h"

#include <algorithm>
#include <cmath>

namespace keta {

// ===================== 6.1 — study sections (ETUDE) =====================

double bar_axis_length(const Model& model, const Bar& bar) {
  // Same axis the element (and therefore the station recovery) uses: the chord
  // between the ECCENTRIC ends, see element.cpp.
  Vec3 pi = model.node(bar.ni).x + bar.ecc[0];
  Vec3 pj = model.node(bar.nj).x + bar.ecc[1];
  double L = (pj - pi).norm();
  require(L > 1e-12, "ETUDE: bar " + std::to_string(bar.id) + " has zero length");
  return L;
}

bool etude_step_is_legal(double lpas) {
  if (!finite(lpas) || lpas <= 0 || lpas > 1) return false;
  double inv = 1.0 / lpas;
  long n = std::lround(inv);
  if (n < 1) return false;
  if (std::fabs(inv - static_cast<double>(n)) > 1e-9 * static_cast<double>(n)) return false;
  // Terminating-decimal denominator: n = 2^a * 5^b (1/3, 1/7 rejected, p.97).
  while (n % 2 == 0) n /= 2;
  while (n % 5 == 0) n /= 5;
  return n == 1;
}

void EtudePlan::insert(Id bar, double L, double x, bool effort, bool depla) {
  std::vector<StudySection>& v = bars_[bar];
  double tol = 1e-9 * std::max(L, 1.0);
  for (StudySection& s : v)
    if (std::fabs(s.x - x) <= tol) {  // same section: merge the families
      s.effort = s.effort || effort;
      s.depla = s.depla || depla;
      return;
    }
  StudySection s;
  s.x = x;
  s.effort = effort;
  s.depla = depla;
  v.insert(std::upper_bound(v.begin(), v.end(), s,
                            [](const StudySection& a, const StudySection& b) { return a.x < b.x; }),
           s);
}

void EtudePlan::add_default(const Model& model) {
  model.index();
  for (const Bar& b : model.bars) {
    double L = bar_axis_length(model, b);
    insert(b.id, L, 0.0, true, false);
    insert(b.id, L, L, true, false);
  }
}

void EtudePlan::add(const Model& model, const EtudeSpec& spec) {
  model.index();
  bool effort = spec.effort, depla = spec.depla;
  if (!effort && !depla) effort = depla = true;

  switch (spec.mode) {
    case EtudeMode::Default:  // "addition aux sections définies PAR DEFAUT"
      bars_.clear();
      add_default(model);
      break;
    case EtudeMode::Seul:  // "sans prise en compte des sections par défaut"
      bars_.clear();
      break;
    case EtudeMode::Sup:  // "addition aux sections DEJA definies"
      break;
  }

  for (const EtudeLine& line : spec.lines) {
    std::vector<Id> ids;
    if (line.all) {
      // TOUT over a model with no active bar (a deck that DESACTIVER TOUT
      // before its PHASAGE) reaches nothing — not an input error.
      for (const Bar& b : model.bars) ids.push_back(b.id);
      if (ids.empty()) continue;
    } else {
      ids = line.bars;
    }
    require(!ids.empty(), "ETUDE: empty bar list (use TOUT to reach every bar)");

    for (Id bid : ids) {
      const Bar& b = model.bar(bid);  // typed error on an unknown bar
      double L = bar_axis_length(model, b);
      double tol = 1e-9 * std::max(L, 1.0);
      std::vector<double> xs;  // in the line's own unit (relative or absolute)

      if (line.step > 0) {
        double lo = 0, hi = line.absolute ? L : 1.0;
        if (line.abscissas.size() == 2) {
          lo = line.abscissas[0];
          hi = line.abscissas[1];
        } else {
          require(line.abscissas.empty(),
                  "ETUDE PAS: the abscissa list must be a range (two values) or empty");
        }
        require(finite(lo) && finite(hi) && hi > lo, "ETUDE PAS: empty or non-finite range");
        double span = hi - lo;
        if (!line.absolute) {
          require(etude_step_is_legal(line.step),
                  "ETUDE PAS " + std::to_string(line.step) +
                      ": a relative step must be an exact fraction 1/n of the bar "
                      "(1/2, 1/5, 1/10 legal; 1/3, 1/7 not) [ST1 v24 p.97]");
          double cnt = span / line.step;
          long n = std::lround(cnt);
          require(n >= 1 && std::fabs(cnt - static_cast<double>(n)) <= 1e-9 * static_cast<double>(n),
                  "ETUDE PAS: the relative step does not divide the requested range evenly");
          for (long k = 0; k <= n; ++k)
            xs.push_back(lo + span * static_cast<double>(k) / static_cast<double>(n));
        } else {
          require(line.step <= span + tol, "ETUDE PAS: absolute step longer than the range");
          long n = static_cast<long>(std::floor(span / line.step + 1e-9));
          for (long k = 0; k <= n; ++k) xs.push_back(lo + line.step * static_cast<double>(k));
          if (std::fabs(xs.back() - hi) > tol) xs.push_back(hi);  // the range end always studied
        }
      } else {
        require(!line.abscissas.empty(), "ETUDE: no abscissa given (SE liste2)");
        xs = line.abscissas;
      }

      for (double v : xs) {
        require(finite(v), "ETUDE: non-finite abscissa");
        double x = line.absolute ? v : v * L;
        require(x >= -tol && x <= L + tol,
                "ETUDE: abscissa " + std::to_string(x) + " m outside bar " + std::to_string(bid) +
                    " (length " + std::to_string(L) + " m)");
        insert(bid, L, std::clamp(x, 0.0, L), effort, depla);
      }
    }
  }
}

const std::vector<StudySection>* EtudePlan::sections(Id bar) const {
  auto it = bars_.find(bar);
  return it == bars_.end() ? nullptr : &it->second;
}

std::map<Id, std::vector<double>> EtudePlan::stations_abs(StudyFamily which) const {
  std::map<Id, std::vector<double>> out;
  for (const auto& [bid, v] : bars_) {
    std::vector<double> xs;
    for (const StudySection& s : v) {
      if (which == StudyFamily::Effort && !s.effort) continue;
      if (which == StudyFamily::Depla && !s.depla) continue;
      xs.push_back(s.x);
    }
    if (!xs.empty()) out.emplace(bid, std::move(xs));
  }
  return out;
}

// ================ 6.2 — generalized stress (CONTRG) ================

namespace {

void set_rule(ContrgDef& def, const std::vector<Id>& bar_ids, const ContrgRule& rule) {
  require(!bar_ids.empty(), "CONTRG: empty bar list");
  for (double c : {rule.or_.cfn, rule.or_.cfty, rule.or_.cftz, rule.or_.cfmx, rule.or_.cfmy,
                   rule.or_.cfmz, rule.or_.cst, rule.mi.cfn, rule.mi.cfty, rule.mi.cftz,
                   rule.mi.cfmx, rule.mi.cfmy, rule.mi.cfmz, rule.mi.cst, rule.ex.cfn,
                   rule.ex.cfty, rule.ex.cftz, rule.ex.cfmx, rule.ex.cfmy, rule.ex.cfmz,
                   rule.ex.cst})
    require(finite(c), "CONTRG: non-finite coefficient");
  for (Id b : bar_ids) def.bars[b] = rule;
}

double lerp_coeff(double a, double b, double xx) { return a + (b - a) * xx; }

double para_coeff(double cor, double cmi, double cex, double xx) {
  // Interpolating parabola through (0, cor), (1/2, cmi), (1, cex) — p.101
  // (whose printed constant term "c = cfex" is a typo, see etude.h).
  double a = 2.0 * (cor - 2.0 * cmi + cex);
  double b = -3.0 * cor + 4.0 * cmi - cex;
  return a * xx * xx + b * xx + cor;
}

}  // namespace

void ContrgDef::set_cte(const std::vector<Id>& bar_ids, const ContrgCoeffs& c) {
  ContrgRule r;
  r.kind = ContrgKind::Cte;
  r.or_ = c;
  set_rule(*this, bar_ids, r);
}

void ContrgDef::set_var_lin(const std::vector<Id>& bar_ids, const ContrgCoeffs& or_,
                            const ContrgCoeffs& ex) {
  ContrgRule r;
  r.kind = ContrgKind::VarLin;
  r.or_ = or_;
  r.ex = ex;
  set_rule(*this, bar_ids, r);
}

void ContrgDef::set_var_para(const std::vector<Id>& bar_ids, const ContrgCoeffs& or_,
                             const ContrgCoeffs& mi, const ContrgCoeffs& ex) {
  ContrgRule r;
  r.kind = ContrgKind::VarPara;
  r.or_ = or_;
  r.mi = mi;
  r.ex = ex;
  set_rule(*this, bar_ids, r);
}

Station st1_signs(const Station& st) {
  // every effort as ST1 lists it (the ST1 edition negates all six, in
  // SPATIALE too; a CONTRG coefficient on TY reads the listed Ty)
  Station s = st;
  s.N = -st.N;
  s.Vy = -st.Vy;
  s.Vz = -st.Vz;
  s.T = -st.T;
  s.My = -st.My;
  s.Mz = -st.Mz;
  return s;
}

ContrgCoeffs contrg_coeffs(const ContrgDef& def, Id bar, double xx) {
  auto it = def.bars.find(bar);
  require(it != def.bars.end(), "CONTRG " + std::to_string(def.id) + ": not defined on bar " +
                                    std::to_string(bar));
  require(finite(xx) && xx >= -1e-9 && xx <= 1.0 + 1e-9,
          "CONTRG: fractional abscissa outside [0, 1]");
  double t = std::clamp(xx, 0.0, 1.0);
  const ContrgRule& r = it->second;
  ContrgCoeffs c;
  switch (r.kind) {
    case ContrgKind::Cte:
      c = r.or_;
      break;
    case ContrgKind::VarLin:
      c.cfn = lerp_coeff(r.or_.cfn, r.ex.cfn, t);
      c.cfty = lerp_coeff(r.or_.cfty, r.ex.cfty, t);
      c.cftz = lerp_coeff(r.or_.cftz, r.ex.cftz, t);
      c.cfmx = lerp_coeff(r.or_.cfmx, r.ex.cfmx, t);
      c.cfmy = lerp_coeff(r.or_.cfmy, r.ex.cfmy, t);
      c.cfmz = lerp_coeff(r.or_.cfmz, r.ex.cfmz, t);
      c.cst = lerp_coeff(r.or_.cst, r.ex.cst, t);
      break;
    case ContrgKind::VarPara:
      c.cfn = para_coeff(r.or_.cfn, r.mi.cfn, r.ex.cfn, t);
      c.cfty = para_coeff(r.or_.cfty, r.mi.cfty, r.ex.cfty, t);
      c.cftz = para_coeff(r.or_.cftz, r.mi.cftz, r.ex.cftz, t);
      c.cfmx = para_coeff(r.or_.cfmx, r.mi.cfmx, r.ex.cfmx, t);
      c.cfmy = para_coeff(r.or_.cfmy, r.mi.cfmy, r.ex.cfmy, t);
      c.cfmz = para_coeff(r.or_.cfmz, r.mi.cfmz, r.ex.cfmz, t);
      c.cst = para_coeff(r.or_.cst, r.mi.cst, r.ex.cst, t);
      break;
  }
  return c;
}

double contrg_at(const ContrgDef& def, Id bar, double xx, const Station& st, bool include_cst) {
  ContrgCoeffs c = contrg_coeffs(def, bar, xx);
  const Station s = def.convention == ForceConvention::St1 ? st1_signs(st) : st;
  double v = c.cfn * s.N + c.cfty * s.Vy + c.cftz * s.Vz + c.cfmx * s.T + c.cfmy * s.My +
             c.cfmz * s.Mz;
  if (include_cst) v += c.cst;
  require(finite(v), "CONTRG " + std::to_string(def.id) + ": non-finite generalized stress");
  return v;
}

// ================ 13.2 / 13.4 — ENV and ENV COMB ================

SectionFibers section_fibers(const Section& sec) {
  SectionFibers f;
  if (!sec.stress_points.empty()) {
    double ymin = 0, ymax = 0, zmin = 0, zmax = 0;
    bool first = true;
    for (const StressPoint& p : sec.stress_points) {
      if (first) {
        ymin = ymax = p.y;
        zmin = zmax = p.z;
        first = false;
      }
      ymin = std::min(ymin, p.y);
      ymax = std::max(ymax, p.y);
      zmin = std::min(zmin, p.z);
      zmax = std::max(zmax, p.z);
    }
    f.vy = std::max(0.0, ymax);
    f.wy = std::max(0.0, -ymin);
    f.vz = std::max(0.0, zmax);
    f.wz = std::max(0.0, -zmin);
    return f;
  }
  // Fallback: symmetric fibers from the elastic moduli (Wely = Iy / z_max).
  const SectionProps& p = sec.props;
  require(p.Wely > 0 && p.Welz > 0,
          "CONTR: section " + std::to_string(sec.id) +
              " has neither stress points nor elastic moduli — extreme fibers "
              "vy/wy/vz/wz are undefined (ST1 CARA p.66)");
  f.vz = f.wz = p.Iy / p.Wely;
  f.vy = f.wy = p.Iz / p.Welz;
  return f;
}

double fiber_sigma(const Model& model, const Bar& bar, const Station& st, ContrFiber fiber,
                   bool st1_sign) {
  const Section& sec = model.section(bar.section);
  SectionFibers f = section_fibers(sec);
  double y = 0, z = 0;
  switch (fiber) {
    case ContrFiber::VY: y = f.vy; break;
    case ContrFiber::WY: y = -f.wy; break;
    case ContrFiber::VZ: z = f.vz; break;
    case ContrFiber::WZ: z = -f.wz; break;
    case ContrFiber::VYVZ: y = f.vy; z = f.vz; break;
    case ContrFiber::VYWZ: y = f.vy; z = -f.wz; break;
    case ContrFiber::WYVZ: y = -f.wy; z = f.vz; break;
    case ContrFiber::WYWZ: y = -f.wy; z = -f.wz; break;
  }
  const SectionProps& p = sec.props;
  require(p.A > 0, "CONTR: section " + std::to_string(sec.id) + " has zero area");
  double s = st.N / p.A;
  if (p.Iy > 0) s += st.My * z / p.Iy;
  if (p.Iz > 0) s -= st.Mz * y / p.Iz;
  return st1_sign ? -s : s;
}

namespace {

// Results share the model's bar order, so the caller's index is a hit in one
// compare; the linear scan is only the fallback (mixed-shape sources).
const BarResult& bar_result(const CaseResult& r, Id bar, std::size_t hint = 0) {
  if (hint < r.bars.size() && r.bars[hint].bar == bar) return r.bars[hint];
  for (const BarResult& br : r.bars)
    if (br.bar == bar) return br;
  fail("ENV: bar " + std::to_string(bar) + " not in the stored results");
}

const Station& station_of(const CaseResult& r, Id bar, std::size_t k, std::size_t hint = 0) {
  const BarResult& br = bar_result(r, bar, hint);
  require(k < br.stations.size(), "ENV: station index out of range on bar " + std::to_string(bar));
  return br.stations[k];
}

double station_value(const Model& model, const EnvSpec& spec, const CaseResult& r, Id bar,
                     std::size_t k, double xx, std::size_t hint = 0) {
  const BarResult& br = bar_result(r, bar, hint);
  require(k < br.stations.size(), "ENV: station index out of range on bar " + std::to_string(bar));
  const Station& st = br.stations[k];
  switch (spec.what) {
    case EnvWhat::Effort:
      require(spec.comp >= 0 && spec.comp < 6, "ENV EFFORT: component out of range (0..5)");
      switch (spec.comp) {
        case 0: return st.N;
        case 1: return st.Vy;
        case 2: return st.Vz;
        case 3: return st.T;
        case 4: return st.My;
        default: return st.Mz;
      }
    case EnvWhat::Contr:
      return fiber_sigma(model, model.bar(bar), st, spec.fiber, spec.sigma_st1_sign);
    case EnvWhat::Contrg:
      require(spec.contrg != nullptr, "ENV CONTRG: no generalized stress given");
      return contrg_at(*spec.contrg, bar, xx, st, spec.contrg_include_cst);
    case EnvWhat::Press: {
      require(spec.comp >= 0 && spec.comp < 4, "ENV PRESS: component out of range (0..3)");
      require(k < br.soil_p.size() && k < br.soil_m.size(),
              "ENV PRESS: bar " + std::to_string(bar) + " carries no Winkler bed");
      return spec.comp < 3 ? br.soil_p[k](spec.comp) : br.soil_m[k];
    }
    case EnvWhat::Depla:
      fail("ENV: DEPLA is a node quantity, not a section quantity");
  }
  fail("ENV: unsupported privileged quantity");
}

// One term of the box: the extreme contribution of a source, choosing the
// favourable coefficient (cmin or cmax) and — for an ENV-of-ENV source — the
// favourable alternative among the wrapped bound states.
template <class F>
double pick_corner(const EnvSource& s, F&& value_of, bool maximize, double& coef,
                   std::size_t& alt, bool skip_null = false) {
  const CaseResult* one[1] = {s.r};
  const CaseResult* const* items = s.alts.empty() ? one : s.alts.data();
  std::size_t count = s.alts.empty() ? 1 : s.alts.size();
  require(finite(s.cmin) && finite(s.cmax), "ENV: non-finite source coefficient");
  // skip_null (EnvSource): an alternative whose value is null here does not
  // compete unless they all are
  const bool by_raw = skip_null && count > 1 && s.raw.size() == count;
  std::vector<double> vals(count), raws(count);
  double scale = 0;
  for (std::size_t a = 0; a < count; ++a) {
    require(items[a] != nullptr, "ENV: source without a result");
    vals[a] = value_of(*items[a]);
    require(finite(vals[a]), "ENV: non-finite source value");
    raws[a] = by_raw ? value_of(*s.raw[a]) : vals[a];
    scale = std::max(scale, std::fabs(raws[a]));
  }
  bool first = true;
  double best = 0;
  for (std::size_t a = 0; a < count; ++a) {
    const double v = vals[a];
    if (skip_null && count > 1 && std::fabs(raws[a]) <= 1e-9 * scale) continue;
    for (double c : {s.cmin, s.cmax}) {
      double t = c * v;
      if (first || (maximize ? t > best : t < best)) {
        best = t;
        coef = c;
        alt = a;
        first = false;
      }
    }
  }
  return best;
}

const CaseResult* source_state(const EnvSource& s, std::size_t alt) {
  return s.alts.empty() ? s.r : s.alts[alt];
}

Vec6 reaction_of(const CaseResult& r, Id node) {
  for (const Reaction& re : r.reactions)
    if (re.node == node) return re.R;
  fail("ENV REAC: node " + std::to_string(node) + " carries no reaction");
}

}  // namespace

EnvSource env_source_of_env(const Envelope& e, double cmin, double cmax) {
  EnvSource s;
  s.cmin = cmin;
  s.cmax = cmax;
  s.alts = {&e.lo.value, &e.hi.value};
  return s;
}

double env_value(const Model& model, const EnvSpec& spec, const CaseResult& r) {
  if (spec.what == EnvWhat::Depla) {
    require(spec.comp >= 0 && spec.comp < 6, "ENV DEPLA: component out of range (0..5)");
    int ni = model.node_index(spec.node);
    require(static_cast<std::size_t>(ni) < r.node_disp.size(),
            "ENV DEPLA: node outside the stored results");
    return r.node_disp[static_cast<std::size_t>(ni)](spec.comp);
  }
  double xx = 0;
  if (spec.what == EnvWhat::Contrg) {
    const Station& st = station_of(r, spec.bar, spec.station);
    xx = st.x / bar_axis_length(model, model.bar(spec.bar));
  }
  return station_value(model, spec, r, spec.bar, spec.station, xx);
}

EnvRun env_run(const Model& model, const EnvSpec& spec, const std::vector<EnvSource>& sources) {
  require(!sources.empty(), "ENV: no sources");
  model.index();
  EnvRun out;
  std::size_t n = sources.size();
  out.coef_max.assign(n, 0.0);
  out.coef_min.assign(n, 0.0);
  out.alt_max.assign(n, 0);
  out.alt_min.assign(n, 0);

  auto value_of = [&](const CaseResult& r) { return env_value(model, spec, r); };
  for (std::size_t i = 0; i < n; ++i) {
    const bool skip = sources[i].skip_null && spec.what != EnvWhat::Depla;
    out.max += pick_corner(sources[i], value_of, true, out.coef_max[i], out.alt_max[i], skip);
    out.min += pick_corner(sources[i], value_of, false, out.coef_min[i], out.alt_min[i], skip);
  }

  // Companion states: the WHOLE state assembled with the very coefficients the
  // extreme was taken at (consistent state, never mixed per-field extrema).
  std::vector<const CaseResult*> smax, smin;
  for (std::size_t i = 0; i < n; ++i) {
    smax.push_back(source_state(sources[i], out.alt_max[i]));
    smin.push_back(source_state(sources[i], out.alt_min[i]));
  }
  out.state_max = combine_results(smax, out.coef_max, 0, "env max");
  out.state_min = combine_results(smin, out.coef_min, 0, "env min");

  if (spec.has_reac) {
    out.has_reac = true;
    out.reac_coef_max.assign(n, 0.0);
    out.reac_coef_min.assign(n, 0.0);
    std::vector<std::size_t> amax(n, 0), amin(n, 0);
    require(spec.reac_comp >= 0 && spec.reac_comp < 6, "ENV REAC: component out of range (0..5)");
    auto reac_val = [&](const CaseResult& r) { return reaction_of(r, spec.reac_node)(spec.reac_comp); };
    for (std::size_t i = 0; i < n; ++i) {
      out.reac_max += pick_corner(sources[i], reac_val, true, out.reac_coef_max[i], amax[i]);
      out.reac_min += pick_corner(sources[i], reac_val, false, out.reac_coef_min[i], amin[i]);
      out.reac_state_max +=
          out.reac_coef_max[i] * reaction_of(*source_state(sources[i], amax[i]), spec.reac_node);
      out.reac_state_min +=
          out.reac_coef_min[i] * reaction_of(*source_state(sources[i], amin[i]), spec.reac_node);
    }
  }
  return out;
}

std::vector<EnvPoint> env_run_sections(const Model& model, const EnvSpec& spec,
                                       const std::vector<EnvSource>& sources) {
  require(!sources.empty(), "ENV: no sources");
  require(spec.what != EnvWhat::Depla,
          "ENV: DEPLA is a node quantity — use env_run() for node displacements");
  model.index();
  std::size_t n = sources.size();
  std::vector<EnvPoint> out;

  // Every source shares the shape of the first one (same model, same stations).
  const CaseResult& shape = *source_state(sources.front(), 0);

  for (std::size_t bi = 0; bi < shape.bars.size(); ++bi) {
    const BarResult& br = shape.bars[bi];
    if (br.stations.empty()) continue;
    if (spec.what == EnvWhat::Contrg) {
      require(spec.contrg != nullptr, "ENV CONTRG: no generalized stress given");
      if (!spec.contrg->defined_for(br.bar)) continue;  // p.101: partial definition is legal
    }
    if (spec.what == EnvWhat::Press && br.soil_p.empty()) continue;  // no Winkler bed here
    double L = bar_axis_length(model, model.bar(br.bar));

    for (std::size_t k = 0; k < br.stations.size(); ++k) {
      EnvPoint p;
      p.bar = br.bar;
      p.station = k;
      p.x = br.stations[k].x;
      p.coef_max.assign(n, 0.0);
      p.coef_min.assign(n, 0.0);
      p.alt_max.assign(n, 0);
      p.alt_min.assign(n, 0);
      double xx = L > 0 ? p.x / L : 0.0;
      auto value_of = [&](const CaseResult& r) {
        return station_value(model, spec, r, br.bar, k, xx, bi);
      };
      for (std::size_t i = 0; i < n; ++i) {
        p.max += pick_corner(sources[i], value_of, true, p.coef_max[i], p.alt_max[i],
                             sources[i].skip_null);
        p.min += pick_corner(sources[i], value_of, false, p.coef_min[i], p.alt_min[i],
                             sources[i].skip_null);
      }
      // Concomitant station state, assembled at the same corner.
      p.st_max.x = p.st_min.x = p.x;
      for (std::size_t i = 0; i < n; ++i) {
        const Station& a = station_of(*source_state(sources[i], p.alt_max[i]), br.bar, k, bi);
        const Station& b = station_of(*source_state(sources[i], p.alt_min[i]), br.bar, k, bi);
        double ca = p.coef_max[i], cb = p.coef_min[i];
        p.st_max.N += ca * a.N;
        p.st_max.Vy += ca * a.Vy;
        p.st_max.Vz += ca * a.Vz;
        p.st_max.T += ca * a.T;
        p.st_max.My += ca * a.My;
        p.st_max.Mz += ca * a.Mz;
        p.st_max.defl += ca * a.defl;
        p.st_max.rot += ca * a.rot;
        p.st_min.N += cb * b.N;
        p.st_min.Vy += cb * b.Vy;
        p.st_min.Vz += cb * b.Vz;
        p.st_min.T += cb * b.T;
        p.st_min.My += cb * b.My;
        p.st_min.Mz += cb * b.Mz;
        p.st_min.defl += cb * b.defl;
        p.st_min.rot += cb * b.rot;
      }
      out.push_back(std::move(p));
    }
  }
  return out;
}

}  // namespace keta

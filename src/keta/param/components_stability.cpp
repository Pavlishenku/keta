// Second Order, Buckling, Modal, Imperfection, Nonlinear, Arc Length, Cable
// Static, Harmonic, PSD and Spectrum components.
//
// Each node drives one keta engine solve (nonlinear.h, eigsolve.h, cable.h,
// dynamics.h) on the model
// wired into M. keta's drivers have no cooperative-cancel hooks yet, so a
// cancel is honoured before the solve starts (never mid-solve) and progress is
// reported once, on completion.
#include <algorithm>
#include <cmath>
#include <cstdint>

#include "keta/cable.h"
#include "keta/dynamics.h"
#include "keta/assembly.h"
#include "keta/eigsolve.h"
#include "keta/nonlinear.h"
#include "keta/param/components.h"

namespace keta::param {

namespace {

ModelRef model_input(const Tree& t) {
  for (const auto& kv : t.branches())
    for (const auto& o : kv.second)
      if (o)
        if (const auto* m = o->get<ModelRef>(); m && *m) return *m;
  throw ComponentError::missing("M");
}

// N pin override, clamped to [0, 100000].
int modes_input(const Tree& t, int dflt) {
  if (auto n = first_int(t)) return static_cast<int>(std::clamp<std::int64_t>(*n, 0, 100'000));
  return dflt;
}

// Reference case: the C pin when wired, else the model's first case.
Id case_input(const Tree& t, const Model& m, const char* what) {
  if (auto c = first_int(t)) return static_cast<Id>(*c);
  if (m.cases.empty()) throw ComponentError(std::string(what) + ": the model has no load cases");
  return m.cases.front().id;
}

template <class F>
auto run_solver(const SolveCtx& ctx, const char* context, F&& f) -> decltype(f()) {
  if (ctx.cancelled()) throw ComponentError::cancelled();
  try {
    return f();
  } catch (const ComponentError&) {
    throw;
  } catch (const Cancelled&) {
    throw ComponentError::cancelled();
  } catch (const std::exception& e) {
    throw ComponentError(std::string(context) + ": " + e.what());
  }
}

int int_param(const json::Value& p, const char* key, int cur, int min) {
  // Non-negative integers only: anything else keeps the current value.
  if (const auto* v = p.get(key); v && v->is_number()) {
    const double x = v->as_number();
    if (x == std::trunc(x) && x >= min && x <= 1e9) return static_cast<int>(x);
  }
  return cur;
}

}  // namespace

// ============================ Second Order ============================

ComponentInfo SecondOrderNode::info() const {
  return {"keta.second_order", "Second Order", "Analyze", "Second-order elastic (P-Delta) static analysis"};
}
std::vector<PinSpec> SecondOrderNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"), PinSpec::item("C", "int", "Load-case id (first case if unset)").opt()};
}
std::vector<PinSpec> SecondOrderNode::outputs() const {
  return {PinSpec::item("R", "results", "Model + second-order results"),
          PinSpec::item("I", "number", "P-Delta iterations used (0 = linear)")};
}
Outputs<Value> SecondOrderNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const Id case_id = case_input(input_or_empty(in, 1), *model, "second-order analysis failed");
  SecondOrderOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.max_iterations = max_iterations;
  SecondOrderResults r =
      run_solver(ctx, "second-order analysis failed", [&] { return solve_second_order(*model, case_id, opt); });
  // The engine counts loop passes; I reports tangent factorizations. They
  // differ only on a zero axial state, where the engine stops after one pass
  // without factorizing: that is the linear solution, reported as I = 0.
  bool zero_state = true;
  for (double n : axial_field(r.result)) zero_state &= (n == 0.0);
  for (const Vec3& s : shell_membrane_field(r.result)) zero_state &= (s.squaredNorm() == 0.0);
  const int iterations = zero_state ? 0 : r.iterations;

  auto out = std::make_shared<AnalyzedModel>();
  out->model = model;
  out->results.cases.push_back(std::move(r.result));
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(ResultsRef(std::move(out)))), single(Value::number(iterations))};
}
json::Value SecondOrderNode::params() const {
  json::Value p = json::Value::object();
  p.set("max_iterations", max_iterations);
  return p;
}
void SecondOrderNode::set_params(const json::Value& p) { max_iterations = int_param(p, "max_iterations", max_iterations, 0); }

// ============================ Buckling ============================

ComponentInfo BucklingNode::info() const {
  return {"keta.buckling", "Buckling", "Analyze", "Linear buckling: λ·(reference case) critical loads"};
}
std::vector<PinSpec> BucklingNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::item("N", "int", "Number of modes (default from the node)").opt(),
          PinSpec::item("C", "int", "Reference load-case id (first case if unset)").opt()};
}
std::vector<PinSpec> BucklingNode::outputs() const {
  return {PinSpec::list("F", "number", "Critical load factors λ (ascending)"),
          PinSpec::item("B", "buckling", "Model + buckling results")};
}
Outputs<Value> BucklingNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  BucklingOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.n_modes = modes_input(input_or_empty(in, 1), n_modes);
  opt.case_id = case_input(input_or_empty(in, 2), *model, "buckling failed");
  BucklingResults r;
  if (opt.n_modes > 0)
    r = run_solver(ctx, "buckling failed", [&] { return solve_buckling(*model, opt); });
  Tree factors;
  factors.ensure_path(Path::single(0));
  for (double f : r.load_factors) factors.append(Path::single(0), Value::number(f));
  auto out = std::make_shared<BuckledModel>();
  out->model = model;
  out->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {factors, single(Value::of(BucklingRef(std::move(out))))};
}
json::Value BucklingNode::params() const {
  json::Value p = json::Value::object();
  p.set("n_modes", n_modes);
  return p;
}
void BucklingNode::set_params(const json::Value& p) { n_modes = int_param(p, "n_modes", n_modes, 0); }

// ============================ Modal ============================

ComponentInfo ModalNode::info() const {
  return {"keta.modal", "Modal", "Analyze", "Modal analysis: natural frequencies & mode shapes"};
}
std::vector<PinSpec> ModalNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"), PinSpec::item("N", "int", "Number of modes (default from the node)").opt()};
}
std::vector<PinSpec> ModalNode::outputs() const {
  return {PinSpec::list("F", "number", "Natural frequencies (Hz, ascending)"),
          PinSpec::item("R", "modal", "Model + modal results"),
          PinSpec::tree("P", "number",
                        "Participation factors Γ: one branch per mode {i} holding 6 numbers [x, y, z, rot-x, rot-y, "
                        "rot-z]. Parallel to F. Zero-filled for a fully-restrained model (no free DOFs)."),
          PinSpec::tree("EM", "number",
                        "Effective modal mass (= participation²): one branch per mode {i} holding 6 numbers [x, y, z, "
                        "rot-x, rot-y, rot-z]. Parallel to F. Zero-filled for a fully-restrained model (no free DOFs).")};
}
Outputs<Value> ModalNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  ModalOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.n_modes = modes_input(input_or_empty(in, 1), n_modes);
  ModalResults r;
  if (opt.n_modes > 0) r = run_solver(ctx, "modal analysis failed", [&] { return solve_modal(*model, opt); });
  Tree freqs, participation, effective_mass;
  freqs.ensure_path(Path::single(0));
  for (double f : r.frequencies_hz) freqs.append(Path::single(0), Value::number(f));
  for (std::size_t i = 0; i < r.frequencies_hz.size(); ++i) {
    const Path p = Path::single(static_cast<int>(i));
    for (int d = 0; d < 6; ++d) {
      participation.append(p, Value::number(i < r.participation.size() ? r.participation[i][d] : 0.0));
      effective_mass.append(p, Value::number(i < r.effective_mass.size() ? r.effective_mass[i][d] : 0.0));
    }
  }
  auto out = std::make_shared<ModalModel>();
  out->model = model;
  out->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {freqs, single(Value::of(ModalRef(std::move(out)))), participation, effective_mass};
}
json::Value ModalNode::params() const {
  json::Value p = json::Value::object();
  p.set("n_modes", n_modes);
  return p;
}
void ModalNode::set_params(const json::Value& p) { n_modes = int_param(p, "n_modes", n_modes, 0); }

// ============================ Arc Length ============================

ComponentInfo ArcLengthNode::info() const {
  return {"keta.arc_length", "Arc Length", "Analyze", "Arc-length continuation: the equilibrium path through limit points"};
}
std::vector<PinSpec> ArcLengthNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::item("C", "int", "Load-case id (first case if unset)").opt(),
          PinSpec::item("L", "number", "Arc length (displacement increment per step)").opt(),
          PinSpec::item("N", "int", "Max accepted steps").opt()};
}
std::vector<PinSpec> ArcLengthNode::outputs() const {
  return {PinSpec::item("T", "arctrace", "Model + traced equilibrium path"),
          PinSpec::list("F", "number", "Load multiplier λ of every converged point")};
}
Outputs<Value> ArcLengthNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const Id case_id = case_input(input_or_empty(in, 1), *model, "arc-length trace failed");
  ArcLengthOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.measure = CorotRotationMeasure::Sin;  // sine-form EICR rotations
  opt.dl0 = first_number(input_or_empty(in, 2)).value_or(0.02);
  opt.max_steps = static_cast<int>(std::clamp<std::int64_t>(first_int(input_or_empty(in, 3)).value_or(400), 1, 1'000'000));
  // The node's settings (the engine's own defaults are looser).
  opt.max_iterations = 20;
  opt.max_cuts = 8;
  opt.tol_rel = 1e-10;
  opt.desired_iterations = desired_iterations;
  opt.stop_lambda = stop_lambda > 0 ? stop_lambda : 0.0;
  opt.stop_disp_inf = 0.0;
  if (!(std::isfinite(opt.dl0) && opt.dl0 > 0))
    throw ComponentError("arc-length trace failed: arc length L must be finite and > 0");
  ArcTrace trace = run_solver(ctx, "arc-length trace failed", [&] { return arc_length_trace(*model, case_id, opt); });
  Tree lambdas;
  lambdas.ensure_path(Path::single(0));
  for (const ArcPoint& p : trace.points) lambdas.append(Path::single(0), Value::number(p.lambda));
  auto out = std::make_shared<ArcTraceModel>();
  out->model = model;
  out->trace = std::move(trace);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(ArcTraceRef(std::move(out)))), lambdas};
}
json::Value ArcLengthNode::params() const {
  json::Value p = json::Value::object();
  p.set("desired_iterations", desired_iterations);
  p.set("stop_lambda", stop_lambda);
  return p;
}
void ArcLengthNode::set_params(const json::Value& p) {
  desired_iterations = int_param(p, "desired_iterations", desired_iterations, 0);
  if (const auto* v = p.get("stop_lambda"); v && v->is_number()) stop_lambda = v->as_number();
}

// ============================ Cable Static ============================

ComponentInfo CableStaticNode::info() const {
  return {"keta.cable_static", "Cable Static", "Analyze", "Elastic-catenary cable + linear-frame static analysis"};
}
std::vector<PinSpec> CableStaticNode::inputs() const {
  return {PinSpec::item("M", "model", "The model (must carry cables)"),
          PinSpec::item("C", "int", "Load-case id (first case if unset)").opt(),
          PinSpec::item("S", "int", "Load substeps (default 1)").opt()};
}
std::vector<PinSpec> CableStaticNode::outputs() const {
  return {PinSpec::item("R", "results", "Model + frame results (deformed shape)"),
          PinSpec::item("K", "cablestatic", "Model + cable-static state (tensions, polylines)"),
          PinSpec::item("I", "number", "Total Newton iterations")};
}
Outputs<Value> CableStaticNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const Id case_id = case_input(input_or_empty(in, 1), *model, "cable-static analysis failed");
  CableStaticOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.stations_per_bar = stations;
  opt.substeps = static_cast<int>(std::clamp<std::int64_t>(first_int(input_or_empty(in, 2)).value_or(1), 1, 100'000));
  CableStaticResults r =
      run_solver(ctx, "cable-static analysis failed", [&] { return solve_cable_static(*model, case_id, opt); });
  // The converged frame case repackaged as ordinary results (one case), so
  // every downstream results consumer reads it unchanged.
  auto frame = std::make_shared<AnalyzedModel>();
  frame->model = model;
  frame->results.cases.push_back(r.result);
  const double iterations = static_cast<double>(r.iterations);
  auto cable = std::make_shared<CableStaticModel>();
  cable->model = model;
  cable->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(ResultsRef(std::move(frame)))), single(Value::of(CableStaticRef(std::move(cable)))),
          single(Value::number(iterations))};
}
json::Value CableStaticNode::params() const {
  json::Value p = json::Value::object();
  p.set("stations", stations);
  return p;
}
void CableStaticNode::set_params(const json::Value& p) { stations = int_param(p, "stations", stations, 0); }

// ============================ Harmonic ============================

ComponentInfo HarmonicNode::info() const {
  return {"keta.harmonic", "Harmonic", "Analyze", "Harmonic forced-vibration (frequency response) by modal superposition"};
}
std::vector<PinSpec> HarmonicNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::item("C", "int", "Load-case index (0-based; default 0)").opt(),
          PinSpec::list("F", "number", "Frequency sweep (Hz) — wire a Series node"),
          PinSpec::item("Z", "number", "Modal damping ratio zeta (default 0.02)").opt(),
          PinSpec::item("N", "int", "Modes to superpose (0 = all finite; default 0)").opt()};
}
std::vector<PinSpec> HarmonicNode::outputs() const {
  return {PinSpec::item("H", "harmonic", "Model + harmonic frequency-response"),
          PinSpec::item("F", "number", "Peak-response frequency (Hz)"),
          PinSpec::item("A", "number", "Peak displacement amplitude (m)"),
          PinSpec::item("R", "results", "Deformed envelope at the peak frequency")};
}
Outputs<Value> HarmonicNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const std::size_t case_index =
      static_cast<std::size_t>(std::max<std::int64_t>(first_int(input_or_empty(in, 1)).value_or(0), 0));
  std::vector<double> freqs;
  for (const Value& v : input_or_empty(in, 2).non_nulls())
    if (auto x = v.as_number()) freqs.push_back(*x);
  if (freqs.empty()) throw ComponentError("harmonic: F needs at least one frequency (Hz) — wire a Series node");
  DampingSpec damping;
  damping.kind = DampingSpec::Kind::Modal;
  damping.zeta = first_number(input_or_empty(in, 3)).value_or(0.02);
  const int n_modes =
      static_cast<int>(std::clamp<std::int64_t>(first_int(input_or_empty(in, 4)).value_or(0), 0, 100'000));
  if (case_index >= model->cases.size())
    throw ComponentError("harmonic analysis failed: load-case index " + std::to_string(case_index) +
                         " is out of range (the model has " + std::to_string(model->cases.size()) + " cases)");
  const LoadCase& lc = model->cases[case_index];
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  HarmonicResults r = run_solver(ctx, "harmonic analysis failed", [&] {
    return solve_harmonic(*model, lc.id, freqs, damping, n_modes, stations, &hooks);
  });
  // Peak translational amplitude over the whole sweep (first maximum wins).
  double peak_amp = 0;
  std::size_t peak_fi = 0;
  for (std::size_t fi = 0; fi < r.freqs_hz.size(); ++fi)
    for (const Vec6& a : r.disp_amp[fi])
      for (int c = 0; c < 3; ++c)
        if (a(c) > peak_amp) {
          peak_amp = a(c);
          peak_fi = fi;
        }
  const double peak_freq = r.freqs_hz.empty() ? 0.0 : r.freqs_hz[peak_fi];
  // Deformed-envelope snapshot: |u| signed by the real part, per component.
  CaseResult snap;
  snap.id = lc.id;
  snap.name = lc.name;
  if (!r.disp_re.empty())
    for (std::size_t ni = 0; ni < r.disp_re[peak_fi].size(); ++ni) {
      Vec6 d;
      for (int c = 0; c < 6; ++c) {
        const double re = r.disp_re[peak_fi][ni](c), im = r.disp_im[peak_fi][ni](c);
        d(c) = std::copysign(std::hypot(re, im), re);
      }
      snap.node_disp.push_back(d);
    }
  auto env = std::make_shared<AnalyzedModel>();
  env->model = model;
  env->results.cases.push_back(std::move(snap));
  auto h = std::make_shared<HarmonicModel>();
  h->model = model;
  h->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(HarmonicRef(std::move(h)))), single(Value::number(peak_freq)),
          single(Value::number(peak_amp)), single(Value::of(ResultsRef(std::move(env))))};
}
json::Value HarmonicNode::params() const {
  json::Value p = json::Value::object();
  p.set("stations", stations);
  return p;
}
void HarmonicNode::set_params(const json::Value& p) { stations = int_param(p, "stations", stations, 0); }

// ============================ Imperfection ============================

ComponentInfo ImperfectionNode::info() const {
  return {"keta.imperfection", "Imperfection", "Analyze", "Seed a buckling mode as an initial geometric imperfection"};
}
std::vector<PinSpec> ImperfectionNode::inputs() const {
  return {PinSpec::item("B", "buckling", "Buckling result (model + modes)"),
          PinSpec::item("M", "int", "Mode index (0 = lowest; default from the node)").opt(),
          PinSpec::item("S", "number", "Imperfection amplitude / scale (m)").opt()};
}
std::vector<PinSpec> ImperfectionNode::outputs() const { return {PinSpec::item("M", "model", "The imperfect model")}; }
Outputs<Value> ImperfectionNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const BucklingRef* b = first_of<BucklingRef>(input_or_empty(in, 0));
  if (!b || !*b) throw ComponentError::missing("B");
  const std::size_t m =
      static_cast<std::size_t>(std::max<std::int64_t>(first_int(input_or_empty(in, 1)).value_or(mode), 0));
  const double s = first_number(input_or_empty(in, 2)).value_or(scale);
  const BucklingResults& r = (*b)->results;
  if (m >= r.modes.size())
    throw ComponentError("imperfection failed: imperfection: buckling mode index " + std::to_string(m) +
                         " out of range (" + std::to_string(r.modes.size()) + " mode(s) available)");
  if (!std::isfinite(s))
    throw ComponentError("imperfection failed: imperfection: non-finite imperfection scale " +
                         Value::number(s).describe());
  return {single(Value::of(ModelRef(std::make_shared<const Model>(apply_imperfection(*(*b)->model, r.modes[m], s)))))};
}
json::Value ImperfectionNode::params() const {
  json::Value p = json::Value::object();
  p.set("mode", mode);
  p.set("scale", scale);
  return p;
}
void ImperfectionNode::set_params(const json::Value& p) {
  mode = int_param(p, "mode", mode, 0);
  if (const auto* v = p.get("scale"); v && v->is_number()) scale = v->as_number();
}

// ============================ Nonlinear ============================

ComponentInfo NonlinearNode::info() const {
  return {"keta.nonlinear", "Nonlinear", "Analyze", "Geometrically-nonlinear (co-rotational) large-displacement analysis"};
}
std::vector<PinSpec> NonlinearNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::item("C", "int", "Load-case id (first case if unset)").opt(),
          PinSpec::item("S", "int", "Load substeps (default 1)").opt()};
}
std::vector<PinSpec> NonlinearNode::outputs() const {
  return {PinSpec::item("D", "nonlinear", "Model + converged large-displacement state"),
          PinSpec::item("I", "number", "Total Newton iterations")};
}
Outputs<Value> NonlinearNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const Id case_id = case_input(input_or_empty(in, 1), *model, "nonlinear analysis failed");
  CorotOptions opt;
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  opt.hooks = &hooks;
  opt.measure = CorotRotationMeasure::Sin;  // sine-form EICR rotations
  opt.tol_rel = 1e-10;
  opt.max_iterations = max_iterations;
  opt.substeps = static_cast<int>(std::clamp<std::int64_t>(first_int(input_or_empty(in, 2)).value_or(1), 1, 100'000));
  CorotResults r = run_solver(ctx, "nonlinear analysis failed", [&] { return solve_corotational(*model, case_id, opt); });
  const double iterations = static_cast<double>(r.iterations);
  auto out = std::make_shared<NonlinearModel>();
  out->model = model;
  out->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(NonlinearRef(std::move(out)))), single(Value::number(iterations))};
}
json::Value NonlinearNode::params() const {
  json::Value p = json::Value::object();
  p.set("max_iterations", max_iterations);
  return p;
}
void NonlinearNode::set_params(const json::Value& p) { max_iterations = int_param(p, "max_iterations", max_iterations, 0); }

// ============================ PSD ============================

ComponentInfo PsdNode::info() const {
  return {"keta.psd", "PSD", "Analyze", "Random (PSD) response to a white-noise-scaled load pattern"};
}
std::vector<PinSpec> PsdNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::item("C", "int", "Load-case index (0-based; default 0)").opt(),
          PinSpec::list("F", "number", "Frequency sweep (Hz) — wire a Series node"),
          PinSpec::item("Z", "number", "Modal damping ratio zeta (default 0.02)").opt(),
          PinSpec::item("S", "number", "White input PSD level S_alpha (default 1.0)").opt()};
}
std::vector<PinSpec> PsdNode::outputs() const {
  return {PinSpec::item("P", "psd", "Model + random-response PSD"),
          PinSpec::item("V", "number", "Peak displacement-PSD value (m^2/Hz)"),
          PinSpec::item("F", "number", "Frequency of the peak PSD (Hz)")};
}
Outputs<Value> PsdNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  const std::size_t case_index =
      static_cast<std::size_t>(std::max<std::int64_t>(first_int(input_or_empty(in, 1)).value_or(0), 0));
  std::vector<double> freqs;
  for (const Value& v : input_or_empty(in, 2).non_nulls())
    if (auto x = v.as_number()) freqs.push_back(*x);
  if (freqs.empty()) throw ComponentError("psd: F needs at least one frequency (Hz) — wire a Series node");
  DampingSpec damping;
  damping.kind = DampingSpec::Kind::Modal;
  damping.zeta = first_number(input_or_empty(in, 3)).value_or(0.02);
  const double s0 = first_number(input_or_empty(in, 4)).value_or(1.0);
  if (case_index >= model->cases.size())
    throw ComponentError("PSD analysis failed: load-case index " + std::to_string(case_index) +
                         " is out of range (the model has " + std::to_string(model->cases.size()) + " cases)");
  const Id case_id = model->cases[case_index].id;
  const std::vector<double> s_alpha(freqs.size(), s0);
  const SolveHooks hooks = ctx_hooks(ctx, info().name);
  PsdResults r = run_solver(ctx, "PSD analysis failed", [&] {
    return solve_psd(*model, case_id, freqs, s_alpha, damping, 0, stations, &hooks);
  });
  double peak = 0;
  std::size_t peak_fi = 0;
  for (std::size_t fi = 0; fi < r.freqs_hz.size(); ++fi)
    for (const Vec6& d : r.disp_psd[fi])
      for (int c = 0; c < 3; ++c)
        if (d(c) > peak) {
          peak = d(c);
          peak_fi = fi;
        }
  const double peak_freq = r.freqs_hz.empty() ? 0.0 : r.freqs_hz[peak_fi];
  auto out = std::make_shared<PsdModel>();
  out->model = model;
  out->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(PsdRef(std::move(out)))), single(Value::number(peak)), single(Value::number(peak_freq))};
}
json::Value PsdNode::params() const {
  json::Value p = json::Value::object();
  p.set("stations", stations);
  return p;
}
void PsdNode::set_params(const json::Value& p) { stations = int_param(p, "stations", stations, 0); }

// ============================ Spectrum ============================

ComponentInfo SpectrumNode::info() const {
  return {"keta.spectrum", "Spectrum", "Analyze", "Response-spectrum (spectral) seismic analysis by modal superposition"};
}
std::vector<PinSpec> SpectrumNode::inputs() const {
  return {PinSpec::item("M", "model", "The model"),
          PinSpec::list("T", "number", "Spectrum periods (s), ascending"),
          PinSpec::list("S", "number", "Spectral pseudo-acceleration Sa (m/s^2)"),
          PinSpec::item("D", "int", "Excitation direction 0..5 (default 0)").opt(),
          PinSpec::item("Z", "number", "Modal damping ratio zeta (default 0.05)").opt(),
          PinSpec::item("R", "text", "Combination rule: srss | cqc (default srss)").opt()};
}
std::vector<PinSpec> SpectrumNode::outputs() const {
  return {PinSpec::item("SP", "spectrum", "Model + response-spectrum results"),
          PinSpec::item("V", "number", "Combined base shear (N)"),
          PinSpec::item("U", "results", "Combined peak displacement field")};
}
Outputs<Value> SpectrumNode::solve(const SolveCtx& ctx, const std::vector<Tree>& in) const {
  ModelRef model = model_input(input_or_empty(in, 0));
  std::vector<double> periods, sas;
  for (const Value& v : input_or_empty(in, 1).non_nulls())
    if (auto x = v.as_number()) periods.push_back(*x);
  for (const Value& v : input_or_empty(in, 2).non_nulls())
    if (auto x = v.as_number()) sas.push_back(*x);
  if (periods.empty() || sas.empty())
    throw ComponentError("spectrum: T (periods) and S (Sa) must each have at least one value");
  if (periods.size() != sas.size())
    throw ComponentError("spectrum: T has " + std::to_string(periods.size()) + " period(s) but S has " +
                         std::to_string(sas.size()) + " Sa value(s) — they must match");
  const int dir = static_cast<int>(std::max<std::int64_t>(first_int(input_or_empty(in, 3)).value_or(0), 0));
  DampingSpec damping;
  damping.kind = DampingSpec::Kind::Modal;
  damping.zeta = first_number(input_or_empty(in, 4)).value_or(0.05);
  std::string rule_txt = "srss";
  for (const Value& v : input_or_empty(in, 5).non_nulls())
    if (const std::string* t = v.as_text()) {
      rule_txt = *t;
      break;
    }
  std::string rule = rule_txt;
  rule.erase(0, rule.find_first_not_of(" \t\r\n"));
  rule.erase(rule.find_last_not_of(" \t\r\n") + 1);
  std::transform(rule.begin(), rule.end(), rule.begin(), [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  SpectrumOptions opt;
  if (rule == "srss")
    opt.rule = ComboRule::Srss;
  else if (rule == "cqc")
    opt.rule = ComboRule::Cqc;
  else
    throw ComponentError("spectrum: rule must be \"srss\" or \"cqc\", got \"" + rule + "\"");
  // Table validation.
  const std::string bad = "spectrum table invalid: bad analysis options: response spectrum: ";
  for (std::size_t i = 0; i < periods.size(); ++i) {
    if (!std::isfinite(periods[i]) || periods[i] < 0)
      throw ComponentError(bad + "every period must be finite and non-negative");
    if (!std::isfinite(sas[i]) || sas[i] < 0) throw ComponentError(bad + "every Sa ordinate must be finite and non-negative");
  }
  for (std::size_t i = 1; i < periods.size(); ++i)
    if (periods[i] <= periods[i - 1]) throw ComponentError(bad + "periods must be strictly increasing");
  if (dir > 5)
    throw ComponentError("response-spectrum analysis failed: excitation direction " + std::to_string(dir) +
                         " out of range (0..5)");
  ResponseSpectrum spectrum;
  spectrum.kind = ResponseSpectrum::Kind::Table;
  for (std::size_t i = 0; i < periods.size(); ++i) spectrum.table.emplace_back(periods[i], sas[i]);
  opt.dir = dir;
  SpectrumResults r = run_solver(ctx, "response-spectrum analysis failed", [&] {
    opt.n_modes = assemble_elements(*model).dofs.n_free;  // all modes
    return solve_spectrum(*model, spectrum, damping, opt);
  });
  const double base_shear = r.base_shear;
  CaseResult snap;
  snap.id = 1;
  snap.name = "spectrum dir " + std::to_string(dir);
  snap.node_disp = r.combined_disp;
  auto env = std::make_shared<AnalyzedModel>();
  env->model = model;
  env->results.cases.push_back(std::move(snap));
  auto out = std::make_shared<SpectrumModel>();
  out->model = model;
  out->results = std::move(r);
  if (ctx.progress) ctx.progress->report(info().name, "done", 1.0);
  return {single(Value::of(SpectrumRef(std::move(out)))), single(Value::number(base_shear)),
          single(Value::of(ResultsRef(std::move(env))))};
}

}  // namespace keta::param

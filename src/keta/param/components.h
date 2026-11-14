// keta::param — the first-party component set.
//
// FORMAT CONTRACT. Type ids ("keta.slider", ...), pin names, pin ORDER and
// persisted parameter keys are part of the document format: append-only.
// Documents serialize wires by (node, pin index): new pins are ALWAYS appended
// and read defensively (`input_or_empty`), never inserted or reordered.
//
// The set: Params (Slider, Int, Panel, Data (Internalized), Cluster / Cluster
// In / Cluster Out), Tree (with the Data Dam), Maths, Lists, Sets, Geometry
// (Point, Vector, Line, Series, Frame Grid, Grid, Polar Grid, Connect),
// Topology (Weld, Intersect & Node, Divide, Shatter, Cull Short Lines, End
// Points), Model (CSV Points, CSV Bars, Material, Catalog Section, Rect
// Section, Section (Exact / Editor / Composite), Beam, Cable, Slab, Surface
// Load, Surface Thermal, Support, Point Mass, Point Load, Line Load, Gravity,
// Linear Load, Point Load (Bar), Thermal Load, Settlement, Foundation, Case
// Label) and Analyze (Assemble, Analyze, Combination, Combination Family,
// Envelope, Second Order, Buckling, Modal, Imperfection, Nonlinear, Arc
// Length, Cable Static, Harmonic, PSD, Spectrum), plus Tendon, Phase, Phased
// Analysis and the moving-load family (Load Path, Influence Line, Convoy,
// Moving Load, IL Extremes, Influence Surface). Assemble declares its full
// pin list (slabs included) so wire indices stay stable.
// A cable is a Cable-kind bar (ids after the beams) carrying EA, w and its L0
// route on CableProps, with one shared inert placeholder material/section.
//
// Conventions (documented, tested):
//   * Gravity uses the standard g = 9.80665 m/s^2;
//   * a bar load authored on a line whose direction is OPPOSITE to the welded
//     bar is mirrored onto the bar (linear-load ends swapped, positions taken
//     from the other end);
//   * a support that blocks no DOF and carries no spring is omitted by
//     Assemble (keta's model validation rejects it; it has no mechanical
//     effect);
//   * Section (Exact) rejects a polygon whose axes are not principal (keta's
//     frame element requires them); catalog profiles carry their root fillets.
#pragma once
#include <atomic>
#include <chrono>
#include <memory>
#include <mutex>
#include <string>

#include "keta/linsolve.h"
#include "keta/param/helpers.h"

namespace keta::param {

// ---------------- Params ----------------

struct NumberSlider : ComponentBase<NumberSlider, Value> {
  double value = 1.0, min = 0.0, max = 10.0;
  std::string label;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct IntSlider : ComponentBase<IntSlider, Value> {
  std::int64_t value = 3, min = 1, max = 30;
  std::string label;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct Panel : ComponentBase<Panel, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// ---------------- Geometry ----------------

struct PointXyz : ComponentBase<PointXyz, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct VectorXyz : ComponentBase<VectorXyz, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct LineAB : ComponentBase<LineAB, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct Series : ComponentBase<Series, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// A whole multi-storey, multi-bay building frame from six numbers: column,
// X-beam and Y-beam lines plus the base and windward-face (x = 0) points, each
// one deterministic single-branch list.
struct FrameGrid : ComponentBase<FrameGrid, Value> {
  static constexpr std::int64_t kMaxDim = 60;
  static constexpr std::size_t kMaxBars = 100'000;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Geometry generators: a point lattice, then a rule that wires it. A BRANCH
// is one vertical fiber, an ITEM index the level.
// Grid: branch {i;j} = the fiber at column (i, j), NZ points at
// O + (i DX, j DY, k DZ). Polar Grid: branch {i} = angular fiber i, NZ points at
// O + (r cos t, r sin t, k H) with r = R[min(k, len-1)] and t = A i / N (full
// +-360 sweep, no duplicate seam) or A i / (N - 1) (partial, both ends) plus
// k TW. Connect: one line per source point to the point found by offsetting its
// tree address — DK on the item, DJ on the last path index (positional within
// the sibling group, W wraps it), DI on the second-to-last (positional, never
// wraps); a missing target is skipped, never an error.
struct GridNode : ComponentBase<GridNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct PolarGridNode : ComponentBase<PolarGridNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct ConnectNode : ComponentBase<ConnectNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// ---------------- Model ----------------

struct MaterialNode : ComponentBase<MaterialNode, Value> {
  std::string preset = "steel";  // "steel" | "concrete" | "custom"
  double e = 210e9, nu = 0.3, rho = 7850.0;
  Material material() const;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct SectionRect : ComponentBase<SectionRect, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct BeamNode : ComponentBase<BeamNode, Value> {
  bool truss = false;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Line -> elastic-catenary cable. EA [N], w [N/m] and the route value V (L0
// in m for "length", H in N for "h", the end-2 tension T in N for
// "end_tension") come from the pins when wired, else from the node params.
struct CableNode : ComponentBase<CableNode, Value> {
  std::string route = "h";
  double ea = 1.0e8, w = 0.0, value = 1.0e4;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Plates & shells. Slab: outline ring (>= 3 points, not closed) + one hole
// ring per branch of H + thickness T + material M + target mesh size S
// (default 0.5 m) + optional result-axes override LX -> a SlabSpec that
// Assemble's PL pin meshes. Surface Load: a uniform pressure W [N/m^2] on a
// slab in case C, frame param global | localnormal | projected. Surface
// Thermal: dt [K] + gradient dT/dz [K/m] on a slab in case C. Both loads carry
// the slab's outline; Assemble resolves the meshed slab by it.
struct SlabNode : ComponentBase<SlabNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct SurfaceLoadNode : ComponentBase<SurfaceLoadNode, Value> {
  std::string frame = "global";
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};
struct SurfaceThermalNode : ComponentBase<SurfaceThermalNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct SupportNode : ComponentBase<SupportNode, Value> {
  std::array<bool, 6> blocked{true, true, true, true, true, true};
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct PointMassNode : ComponentBase<PointMassNode, Value> {
  double mass = 100.0;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Shared shape of the case-carrying load nodes (case id = node state).
template <class Derived>
struct CaseNode : ComponentBase<Derived, Value> {
  Id case_id = 1;
  json::Value params() const override {
    json::Value p = json::Value::object();
    p.set("case", case_id);
    return p;
  }
  void set_params(const json::Value& p) override {
    if (const auto* c = p.get("case"); c && c->is_number()) case_id = c->as_int();
  }
};

// "global" | "local" | "projected" -> LoadFrame (unknown strings read as global).
LoadFrame parse_frame(const std::string& s);

template <class Derived>
struct FramedCaseNode : CaseNode<Derived> {
  std::string frame = "global";
  json::Value params() const override {
    json::Value p = CaseNode<Derived>::params();
    p.set("frame", frame);
    return p;
  }
  void set_params(const json::Value& p) override {
    CaseNode<Derived>::set_params(p);
    if (const auto* f = p.get("frame"); f && f->is_string()) frame = f->as_string();
  }
};

struct PointLoadNode : CaseNode<PointLoadNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct LineLoadNode : FramedCaseNode<LineLoadNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct GravityNode : CaseNode<GravityNode> {
  double factor = 1.0;  // multiplier on g
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct LinearLoadNode : FramedCaseNode<LinearLoadNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct BarPointLoadNode : FramedCaseNode<BarPointLoadNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct ThermalLoadNode : CaseNode<ThermalLoadNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct SettlementNode : CaseNode<SettlementNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct FoundationNode : ComponentBase<FoundationNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct CaseLabelNode : CaseNode<CaseLabelNode> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// ---------------- Analyze ----------------

// Beams + supports + loads (+ masses, combinations, foundations) -> a validated
// keta::Model. Coincident endpoints weld into shared nodes within `tolerance`
// (two members transfer force iff they share a node); sections and materials
// deduplicate by value; loads group into cases by case id (sorted).
struct AssembleNode : ComponentBase<AssembleNode, Value> {
  double tolerance = 1e-3;  // node-merge tolerance [m]
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// A factorization persisting across recomputes, SHARED by copies of the node
// that owns it (an async worker evaluates a copy of the graph: a cold cache per
// dispatch would refactor every time). FactorCache exact-compares K before any
// reuse, so the cache is a pure performance hint.
struct SharedFactorCache {
  std::mutex mu;
  FactorCache cache;
  FactorCache::Stats stats() {
    std::lock_guard<std::mutex> lk(mu);
    return cache.stats();
  }
};

// Linear static analysis (factorization kept in a SharedFactorCache), with
// the section-point stresses of every station on demand (exact Tri6 solves
// kept in a SharedSectionCache, see LazySectionStress).
struct AnalyzeNode : ComponentBase<AnalyzeNode, Value> {
  int stations = 11;
  using SharedCache = SharedFactorCache;
  std::shared_ptr<SharedCache> cache = std::make_shared<SharedCache>();
  std::shared_ptr<SharedSectionCache> stress_cache = std::make_shared<SharedSectionCache>();
  FactorCache::Stats cache_stats() const;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct CombinationNode : ComponentBase<CombinationNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Second-order elastic (P-Delta) statics of one load case (C, first case when
// unwired): the fixed point of [K + Kg(N(u))] u = F, repackaged as ordinary
// one-case Results so every downstream consumer reads it like Analyze's. I is
// the number of tangent factorizations (0 = zero axial state, i.e. the linear
// solution).
struct SecondOrderNode : ComponentBase<SecondOrderNode, Value> {
  int max_iterations = 30;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Linear (eigenvalue) buckling of the reference case C (first when unwired):
// F = the critical multipliers lambda of that case's loads, ascending; an
// all-tension / zero-axial reference yields an empty F (a result, not an error).
struct BucklingNode : ComponentBase<BucklingNode, Value> {
  int n_modes = 6;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Free vibration K phi = om^2 M phi (consistent mass): F = frequencies [Hz]
// ascending, P / EM = participation factors / effective modal masses, one
// branch {i} of six numbers [x, y, z, rx, ry, rz] per mode.
struct ModalNode : ComponentBase<ModalNode, Value> {
  int n_modes = 6;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Arc-length continuation (Crisfield cylindrical constraint) of case C (first
// when unwired) through limit points: T = the traced path, F = the load
// multiplier of every converged point (the zero-load start included). L =
// arc radius (displacement-norm increment, default 0.02), N = accepted-step
// cap (default 400). desired_iterations = 0 keeps L fixed; stop_lambda <= 0
// disables the |lambda| stop.
struct ArcLengthNode : ComponentBase<ArcLengthNode, Value> {
  int desired_iterations = 0;
  double stop_lambda = 0.0;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Seed buckling mode M (0 = lowest) of a Buckling result as an initial
// geometric imperfection: every node offset by S x the mode's translations
// (unit-crest shape, so S is the peak offset in m). Emits a plain model.
struct ImperfectionNode : ComponentBase<ImperfectionNode, Value> {
  int mode = 0;
  double scale = 0.01;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Large-displacement corotational statics of case C (first when unwired) in S
// load substeps, with the sin-form rotation measure: D = the converged state,
// I = Newton steps. Span loads and self-weight are accepted (keta's
// corotational driver lowers them on the current frame, with their load
// stiffness in the tangent).
struct NonlinearNode : ComponentBase<NonlinearNode, Value> {
  int max_iterations = 30;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Random response to the 0-based case C's pattern scaled by a white-noise
// multiplier of level S over the sweep F, one modal damping ratio Z, all
// modes: P = the response PSD, V / F = the peak translational displacement PSD
// and its frequency. keta's RMS integration needs >= 2 frequencies.
struct PsdNode : ComponentBase<PsdNode, Value> {
  int stations = 11;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Response spectrum: table (T periods ascending, S pseudo-accelerations),
// direction D 0..5, one modal damping ratio Z, rule R srss | cqc, all modes:
// SP = the results, V = the combined base shear, U = the combined peak
// displacement field (componentwise >= 0; no stations, no reactions).
struct SpectrumNode : ComponentBase<SpectrumNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Coupled elastic-catenary cables + linear frame, case C (first when unwired),
// S load substeps (cable self-weight always full): R = the converged frame
// case as ordinary results, K = the cable state (per-cable H / T max), I =
// Newton steps.
struct CableStaticNode : ComponentBase<CableStaticNode, Value> {
  int stations = 11;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Steady-state harmonic response by modal superposition: case C is a 0-BASED
// index (not an id), F the frequency sweep [Hz], Z one modal damping ratio, N
// the modes superposed (0 = all). F / A out = the peak translational amplitude
// over the sweep and its frequency; R = the deformed envelope there (each
// component |u| signed by its real part; no stations, no reactions — an
// amplitude envelope has no quasi-static diagram).
struct HarmonicNode : ComponentBase<HarmonicNode, Value> {
  int stations = 11;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// ---------------- Maths ----------------
//
// Item-wise over map_items (longest list). A NaN / inf is ALWAYS a typed error:
// domain guards reject the knowable cases (division / modulo by zero, sqrt of a
// negative, invalid powers) before computing and every result is
// finiteness-checked, failing the whole node — never a silent value seeded
// into a structural model. Add / Subtract / Multiply / Divide carry the
// Grasshopper mixed overloads (vector +- vector, point +- vector, point - point
// = vector, vector * number, vector / number).

enum class UnaryOp { Negate, Abs, Sqrt, Sin, Cos, Tan, Floor, Ceil, Round };
struct UnaryMathNode : ComponentBase<UnaryMathNode, Value> {
  UnaryOp op = UnaryOp::Negate;
  explicit UnaryMathNode(UnaryOp o = UnaryOp::Negate) : op(o) {}
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

enum class BinaryOp { Add, Subtract, Multiply, Divide, Power, Modulo, Larger, Smaller, And, Or };
struct BinaryMathNode : ComponentBase<BinaryMathNode, Value> {
  BinaryOp op = BinaryOp::Add;
  explicit BinaryMathNode(BinaryOp o = BinaryOp::Add) : op(o) {}
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct NotNode : ComponentBase<NotNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct EqualsNode : ComponentBase<EqualsNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct RangeNode : ComponentBase<RangeNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct RemapNode : ComponentBase<RemapNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct MassAdditionNode : ComponentBase<MassAdditionNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// A textual formula (see keta/param/expr.h) over the FIXED variable pins
// x y z u v w (fixed so wires never dangle when the formula changes; an
// unwired variable is 0). Parsed once per solve; a parse error names the
// source position.
struct ExpressionNode : ComponentBase<ExpressionNode, Value> {
  std::string formula = "x+y";
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// ---------------- Lists, Trees, Sets ----------------
//
// Lists operate PER BRANCH (paths preserved, nothing flattened across
// branches); secondary pattern / key / index inputs align branch-for-branch,
// repeating their last branch. Nulls are first-class members. Trees are pure
// value-agnostic restructurings delegating to DataTree. Sets select by
// geometric predicate and carry named FLAT-index selections (SetSpec). An empty
// input flows out empty — an operation on nothing is nothing, never an error.

#define KETA_PARAM_STATELESS(Name)                                                         \
  struct Name : ComponentBase<Name, Value> {                                               \
    ComponentInfo info() const override;                                                   \
    std::vector<PinSpec> inputs() const override;                                          \
    std::vector<PinSpec> outputs() const override;                                         \
    Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;     \
  }

// Lists
KETA_PARAM_STATELESS(ListItemNode);
KETA_PARAM_STATELESS(ListLengthNode);
KETA_PARAM_STATELESS(SubListNode);
KETA_PARAM_STATELESS(ShiftListNode);
KETA_PARAM_STATELESS(ListReverseNode);
KETA_PARAM_STATELESS(ListInsertNode);
KETA_PARAM_STATELESS(SortListNode);
KETA_PARAM_STATELESS(CullPatternNode);
KETA_PARAM_STATELESS(DispatchNode);
// Trees
KETA_PARAM_STATELESS(GraftNode);
KETA_PARAM_STATELESS(FlattenNode);
KETA_PARAM_STATELESS(SimplifyNode);
KETA_PARAM_STATELESS(TrimTreeNode);
KETA_PARAM_STATELESS(FlipMatrixNode);
KETA_PARAM_STATELESS(PruneTreeNode);
KETA_PARAM_STATELESS(ShiftPathsNode);
KETA_PARAM_STATELESS(MergeTreesNode);
KETA_PARAM_STATELESS(EntwineNode);
KETA_PARAM_STATELESS(CleanTreeNode);
KETA_PARAM_STATELESS(ParamViewerNode);
// Sets
KETA_PARAM_STATELESS(CullIndexNode);
KETA_PARAM_STATELESS(CullDuplicatesNode);
KETA_PARAM_STATELESS(SortPointsNode);
KETA_PARAM_STATELESS(SelectBoxNode);
KETA_PARAM_STATELESS(SelectLevelNode);
KETA_PARAM_STATELESS(SelectLengthNode);
KETA_PARAM_STATELESS(SelectOrientationNode);
KETA_PARAM_STATELESS(SelectProximityNode);
KETA_PARAM_STATELESS(SelectBySectionNode);
KETA_PARAM_STATELESS(SelectByMaterialNode);
KETA_PARAM_STATELESS(GroupByNode);
KETA_PARAM_STATELESS(CreateSetNode);
KETA_PARAM_STATELESS(SetMembersNode);
KETA_PARAM_STATELESS(FilterBySetNode);

#undef KETA_PARAM_STATELESS

// Compact synopsis of a tree (Param Viewer, pin tooltips): "{n} branches,
// {m} items (types)[, k null]" then one "{path} (N=k)" line per branch, capped
// at 30 branches.
std::string tree_synopsis(const Tree& t);

// Union / intersection / difference / symmetric difference of two sets.
struct SetOpNode : ComponentBase<SetOpNode, Value> {
  std::string op = "union";  // "union" | "intersection" | "difference" | "symmetric difference"
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// A named set from a 3D viewport pick: model user ids of `kind` ("nodes" |
// "bars"). With the element list E wired, bar id k maps to the flat slot of the
// k-th NON-NULL element (Assemble numbers bars over the non-null beams);
// otherwise index = id - 1. Ids are positional: an upstream reorder re-targets
// a stored pick (documented limitation).
struct SelectionNode : ComponentBase<SelectionNode, Value> {
  std::string kind = "bars";
  std::vector<std::int64_t> ids;
  std::string name = "Selection";
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// A European rolled profile by designation (keta/catalog.h). The optional D
// pin maps one section per D SLOT, nulls included (a null takes the param
// designation) so the output stays aligned with the lines it is paired with.
struct SectionCatalogNode : ComponentBase<SectionCatalogNode, Value> {
  std::string designation = "IPE300";
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Exact section of an arbitrary polygon: outer ring P plus hole rings H (one
// branch each), points read as (y, z) = (x, y). The closed-form polygon
// properties (A, Iy, Iz, Wel) come from Green's theorem; the Tri6 warping FE
// (section_exact.h) replaces J, Avy, Avz and yields Iw and the shear centre.
// `fineness` scales the ~3000-element standard mesh (N overrides it).
struct SectionExactNode : ComponentBase<SectionExactNode, Value> {
  double fineness = 1.0;
  double nu = 0.3;  // Poisson ratio of the shear-area integrals
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Exact section whose polygon lives in the params (edited interactively in a
// GUI), no input pins; outputs identical in order and tags to Section (Exact)
// and produced by the same core. Params in SI metres: outer [[y, z], ...]
// (>= 3 points, a shorter ring is ignored on load), holes, nu, fineness,
// grid_step and snap (editor-only, carried verbatim).
struct SectionEditorNode : ComponentBase<SectionEditorNode, Value> {
  std::vector<Eigen::Vector2d> outer{{-0.15, -0.25}, {0.15, -0.25}, {0.15, 0.25}, {-0.15, 0.25}};
  std::vector<std::vector<Eigen::Vector2d>> holes;
  double nu = 0.3;
  double fineness = 1.0;
  double grid_step = 0.05;  // editor grid / snap step [m]
  bool snap = true;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Modular-ratio (transformed) multi-material section: region rings R (one
// branch each) matched with materials M, the first material being the
// reference. Overlapping rings resolve last-wins; the analysed domain is the
// union of the rings. Outputs the transformed section, the reference material
// and the effective J_eff = (GJ)_eff / G_ref and Iw.
struct SectionCompositeNode : ComponentBase<SectionCompositeNode, Value> {
  double fineness = 1.0;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Generative combination family from permanent (P) and variable (V) case ids:
// G+ and G- permanent-only combinations, then for each leading variable case a
// G+ and a G- combination with GV on the leader and GV*psi on the others. Every
// factor is an input with a plain default (1.35 / 1.0 / 1.5 / 0.7): no design
// code is built in.
struct ComboFamilyNode : ComponentBase<ComboFamilyNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Min/max envelope of an analyzed model over its combinations (else its
// cases) with the dominant source behind every extremum; MIN / MAX / DMIN /
// DMAX pull one station quantity Q (N, Vy, Vz, Mt, My, Mz; default My) out as
// one branch per bar, one item per station.
struct EnvelopeNode : ComponentBase<EnvelopeNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// ---------------- Moving loads (influence lines & convoys) ----------------
//
// The canvas surface over keta's influence engine (influence.h,
// influence_surface.h): a travel path of contiguous beams, a response target,
// convoy patterns, their extremes and the per-position load cases of a full
// re-analysis sweep. Ordinates are response per unit force along the load
// direction D (default (0, 0, -1); its magnitude scales them).

struct LoadPathNode : ComponentBase<LoadPathNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Q: N, Vy, Vz, Mt, My, Mz (bar force at the bar section nearest T), Rx..Rmz
// (reaction) or dx..rz (displacement) at the node nearest T; exact case first
// ("Rz" reaction vs "rz" rotation), then case-insensitive. S samples per bar.
// The factorization persists across recomputes (SharedFactorCache): moving the
// target, the path or the quantity on an unchanged model costs one solve.
struct InfluenceLineNode : ComponentBase<InfluenceLineNode, Value> {
  std::shared_ptr<SharedFactorCache> cache = std::make_shared<SharedFactorCache>();
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

struct ConvoyNode : ComponentBase<ConvoyNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Sweeps a convoy along a path: one load set per reference position (case
// B + k), and one ULS combination per loaded position when G is wired.
struct MovingLoadNode : ComponentBase<MovingLoadNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Worst response of an influence line under an optional convoy and / or an
// optional UDL on the sign-matched zones (the load-absent state is always a
// candidate: MAX >= 0 >= MIN).
struct IlExtremesNode : ComponentBase<IlExtremesNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Influence surface over a deck node grid (one branch of G per row, each
// point snapped to the nearest model node). Factorization kept as for
// InfluenceLineNode.
struct InfluenceSurfaceNode : ComponentBase<InfluenceSurfaceNode, Value> {
  std::shared_ptr<SharedFactorCache> cache = std::make_shared<SharedFactorCache>();
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// ---------------- Topology & import ----------------

// Weld every line endpoint into a shared node cloud within T (default 1 mm),
// through Assemble's own welder: L (endpoints snapped, collapsed lines null),
// N (unique nodes, first-appearance order), A / B (0-based node index of each
// line end, aligned with the input).
struct WeldNode : ComponentBase<WeldNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// The rest of the Topology family, same tolerance doctrine (T defaults to
// 1 mm, Assemble's weld tolerance):
//   Intersect & Node — split crossing / near-missing / T-junction lines at ONE
//     shared point per meeting (L segments in the source branch, S = per-branch
//     source index, X = the meeting points); parallel overlaps and
//     endpoint-to-endpoint touches are left to Weld.
//   Divide — N equal segments per line, child branch {p;n} per line (P points,
//     S segments).
//   Shatter — split lines at the projection of the points of P within T.
//   Cull Short Lines — drop lines shorter than T (C count, M mask).
//   End Points — each line's start A and end B, aligned with the input.
struct IntersectNode : ComponentBase<IntersectNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct DivideLineNode : ComponentBase<DivideLineNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct ShatterNode : ComponentBase<ShatterNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct CullShortNode : ComponentBase<CullShortNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};
struct EndPointsNode : ComponentBase<EndPointsNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Points from a CSV file: 0-based column mapping, `unit` (m / cm / mm / in /
// ft) scaled to metres at the parse edge, one item per data row (null +
// warning on a bad row), ids from a column or the 1-based row number.
struct CsvPointsNode : ComponentBase<CsvPointsNode, Value> {
  std::string path;
  std::string delimiter = "auto";  // "auto" | "," | ";" | "\t"
  bool has_header = true;
  std::string unit = "m";
  std::optional<std::size_t> map_id, map_x = 0, map_y = 1, map_z = 2;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Bars from a CSV file whose end columns reference the ids of a wired
// points + ids source (CSV Points), plus section / material tag columns.
struct CsvBarsNode : ComponentBase<CsvBarsNode, Value> {
  std::string path;
  std::string delimiter = "auto";
  bool has_header = true;
  std::optional<std::size_t> map_id;
  std::size_t map_a = 0, map_b = 1;
  std::optional<std::size_t> map_section, map_material;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// ---------------- Prestress & construction phasing ----------------

// A prestressing tendon: one tendon per branch of L (chain lines or beams in
// travel order), its eccentricity profile from the passage points P (x = the
// normalized station, y / z = e_y / e_z in m) or, in Polynomial mode, the
// xi-basis coefficients Y / Z, and its force profile F over the normalized
// stations S. `mode`: Polyline / Parabola / Spline / Polynomial; `kind`:
// Post / Pre (data only). EP / AP (appended) feed the staged K = 1 losses.
struct TendonNode : ComponentBase<TendonNode, Value> {
  std::string mode = "Parabola";
  std::string kind = "Post";
  std::string name;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// One construction stage: bars / supports switched on and off, cases applied,
// tendons jacked (simultaneous) to a target fraction of their declared force.
struct PhaseNode : ComponentBase<PhaseNode, Value> {
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
};

// Staged construction: the model and its ordered phases resolved (bars by
// endpoints, supports by position and blocked pattern, tendons structurally)
// into a keta::PhasePlan and solved by incremental superposition.
// The per-substructure factor cache persists across recomputes and is shared
// by copies of the node, like Analyze's (a pure performance hint: every slot
// exact-compares K before reuse).
struct SharedPhasedFactorCache {
  std::mutex mu;
  PhasedFactorCache cache;
  FactorCache::Stats stats() {
    std::lock_guard<std::mutex> lk(mu);
    return cache.stats();
  }
};
struct PhasedAnalyzeNode : ComponentBase<PhasedAnalyzeNode, Value> {
  int stations = 11;
  double tolerance = 1e-3;  // geometry resolution against the model [m]
  std::shared_ptr<SharedPhasedFactorCache> cache = std::make_shared<SharedPhasedFactorCache>();
  std::shared_ptr<SharedSectionCache> stress_cache = std::make_shared<SharedSectionCache>();
  FactorCache::Stats cache_stats() const { return cache->stats(); }
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// The model tendon record Assemble lowers a TendonSpec into (no bar chain:
// that comes from Assemble's welding). Shared so Phased Analysis matches a
// phase's tendon against Model::tendons with exactly Assemble's fields.
Tendon tendon_record(const TendonSpec& spec);

// ---------------- Internalized data and Data Dam ----------------

// The static pin-tag universe: a persisted tag outside it (a future or
// malformed tag) maps to "any", a universal wire.
const std::vector<std::string>& all_type_tags();
std::string known_type_tag(const std::string& tag);

// Data (Internalized): a source node carrying a frozen DataTree in its params
// (a long-lived document surface):
//   {"branches": [{"path": [0], "items": [{"t": "num", "v": 1.5}, null]}]}
// Item tags num / int / bool / text / point / vector / line are append-only;
// an unknown or malformed item is kept VERBATIM (re-saved unchanged) and is a
// typed error only when solved. The output pin's tag is the shared tag of the
// non-null items, else "any". JSON numbers are doubles in keta: an int item
// is exact up to 2^53 and an integral "v" is accepted as an int.
struct BakedDataNode : ComponentBase<BakedDataNode, Value> {
  struct Item {
    std::optional<Value> value;        // a bakeable value
    std::optional<json::Value> raw;    // an unknown item, verbatim
  };  // both empty = a null slot
  struct BranchRec {
    std::vector<std::int64_t> path;
    std::vector<Item> items;
  };
  std::vector<BranchRec> branches;
  std::string dominant = "any";

  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
  // Freeze a solved tree (throws BakeError on an unbakeable kind or a
  // non-finite float).
  static BakedDataNode from_tree(const Tree& t);
};

// Why a bake was refused; the wording of what() is contractual.
class BakeError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// Bake node `id` in place: its solved single-output tree becomes a Data
// (Internalized) node at the same position, the output's consumers are
// re-wired onto it and the source node is removed (its upstream cone stays).
// Throws BakeError, leaving the graph untouched, when the node is missing, has
// not solved, has several outputs or carries unbakeable data.
NodeId bake_node(Graph<Value>& g, NodeId id);

// Data Dam: a value-agnostic pass-through that HOLDS the last data it let
// through until released. The first solve (or the first after a document
// load) captures, holds and passes the input through; afterwards it returns
// the held trees unchanged (the engine's early cutoff keeps the downstream
// cone cached) and only records that fresh input is waiting (stale, with the
// time it first diverged). request_release() queues one release for the next
// solve. The working state is NOT persisted and is SHARED between clones (the
// async worker solves a copy of the graph, the UI observes the same cell).
// Mode: manual, or auto-release after N s of staleness (timing driven by the
// GUI; the component only records stamps).
struct DataDamNode : ComponentBase<DataDamNode, Value> {
  struct State {
    std::mutex m;
    std::optional<std::vector<Tree>> held;
    std::atomic<bool> release_pending{false};
    std::atomic<bool> stale{false};
    std::optional<std::chrono::steady_clock::time_point> stale_since;
  };
  std::shared_ptr<State> state = std::make_shared<State>();
  std::optional<float> after_secs;  // nullopt = manual

  void request_release() const { state->release_pending = true; }
  bool is_stale() const { return state->stale; }
  std::optional<std::chrono::steady_clock::time_point> stale_since() const;

  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// ---------------- Clusters ----------------

// A subgraph packaged as one node. The cluster's params carry the inner
// graph as a graph-only document STRING (echoed verbatim on save, so a
// load -> save is byte-stable) plus its boundary pins [{name, tag, dir}]. The
// live inner graph is installed by the document layer (`hydrate_clusters`,
// document.h) after a load or paste — the solve context carries no registry
// — and is shared, immutable, between clones. A solve copies it, injects the
// k-th input into the k-th cluster_in marker and harvests the k-th
// cluster_out (POSITIONAL, in node order — the fold creates the markers in
// pin order), forwarding the cancel / progress context; an inner failure is
// reported naming the inner node.
inline constexpr const char* kClusterTypeId = "keta.cluster";
inline constexpr const char* kClusterInTypeId = "keta.cluster_in";
inline constexpr const char* kClusterOutTypeId = "keta.cluster_out";

// Inner-graph input marker: emits the tree its cluster injects (empty when
// standalone).
struct ClusterInNode : ComponentBase<ClusterInNode, Value> {
  std::string name, tag;
  std::optional<Tree> injected;  // runtime only
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

// Inner-graph output marker: passes its input through for the cluster to
// harvest.
struct ClusterOutNode : ComponentBase<ClusterOutNode, Value> {
  std::string name, tag;
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;
};

struct ClusterPin {
  std::string name, tag;
  bool out = false;  // "in" / "out"; an unknown dir reads as "in"
  friend bool operator==(const ClusterPin&, const ClusterPin&) = default;
};

struct ClusterNode : ComponentBase<ClusterNode, Value> {
  std::string name;
  std::string doc;                // the inner graph-only document (source of truth)
  std::vector<ClusterPin> pins;   // declaration order, inputs and outputs interleaved
  std::shared_ptr<const Graph<Value>> inner;  // hydrated template; null until hydrated
  ComponentInfo info() const override;
  std::vector<PinSpec> inputs() const override;
  std::vector<PinSpec> outputs() const override;
  Outputs<Value> solve(const SolveCtx&, const std::vector<Tree>&) const override;
  json::Value params() const override;
  void set_params(const json::Value&) override;  // drops any hydration
};

// ---------------- Registry ----------------

struct ComponentFactory {
  ComponentInfo info;
  std::function<std::unique_ptr<Comp>()> make;
};

class Registry {
 public:
  void push(ComponentFactory f) { factories_.push_back(std::move(f)); }
  const ComponentFactory* by_type_id(const std::string& id) const;
  // Follows at most ONE alias hop for a retired id (documents keep loading).
  const ComponentFactory* resolve_type_id(const std::string& id) const;
  std::unique_ptr<Comp> make(const std::string& id) const;  // nullptr when unknown
  std::vector<std::string> categories() const;              // display order
  const std::vector<ComponentFactory>& factories() const { return factories_; }

 private:
  std::vector<ComponentFactory> factories_;
};

// The standard first-party component set, palette-ordered.
Registry standard_registry();

}  // namespace keta::param

// keta::param — the universal wire value of the parametric layer and the
// geometric specifications the model components exchange before `Assemble`
// welds them into a keta::Model.
//
// Specs are GEOMETRIC, never id-based: a beam carries its axis line, a support
// its point, a load its lines — Assemble resolves them onto welded nodes and
// bars within its tolerance. All quantities are SI (m, N, Pa, kg, rad).
//
// EQUALITY drives the engine's early cutoff, so every alternative must compare
// by content, cheaply:
//   * light values (numbers, points, specs) compare member-wise;
//   * keta::Material / keta::Section have no operator== in the kernel, so the
//     comparison is written out here (materials_equal / sections_equal) — a
//     field added to either struct must be added there too;
//   * heavy shared payloads (Model, Results, Buckling, Modal, Envelope, Phased,
//     Influence, CableStatic, Harmonic, ArcTrace, Nonlinear, Psd, Spectrum)
//     compare by POINTER IDENTITY. This
//     is conservative: a re-solve always reads as "changed" and the downstream
//     cone recomputes — never a stale result, only a missed cutoff. Deep model
//     equality is a later refinement.
//
// Serialization (documents, baked data) is by type tag, never by variant
// position, so the alternative order below carries no compatibility promise.
#pragma once
#include <array>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <variant>
#include <vector>

#include "keta/cable.h"
#include "keta/dynamics.h"
#include "keta/eigsolve.h"
#include "keta/envelope.h"
#include "keta/influence.h"
#include "keta/model.h"
#include "keta/nonlinear.h"
#include "keta/phasing.h"
#include "keta/section_stress.h"
#include "keta/statics.h"

namespace keta::param {

// ---------------- geometry ----------------

struct Point {
  Vec3 p = Vec3::Zero();
  friend bool operator==(const Point& a, const Point& b) { return a.p == b.p; }
};
struct Vector {
  Vec3 v = Vec3::Zero();
  friend bool operator==(const Vector& a, const Vector& b) { return a.v == b.v; }
};
struct Line {
  Vec3 a = Vec3::Zero(), b = Vec3::Zero();
  double length() const { return (b - a).norm(); }
  friend bool operator==(const Line& x, const Line& y) { return x.a == y.a && x.b == y.b; }
};

// ---------------- model specs ----------------

bool materials_equal(const Material& a, const Material& b);
bool sections_equal(const Section& a, const Section& b);

struct MaterialSpec {
  Material m;  // id ignored (Assemble renumbers)
  friend bool operator==(const MaterialSpec& a, const MaterialSpec& b) { return materials_equal(a.m, b.m); }
};

// A finalized section (props computed); `name` is its display / dedup name.
struct SectionSpec {
  Section s;
  friend bool operator==(const SectionSpec& a, const SectionSpec& b) { return sections_equal(a.s, b.s); }
};

struct BeamSpec {
  Line line;
  SectionSpec section;
  MaterialSpec material;
  bool truss = false;
  Releases releases{};
  double roll = 0.0;          // rotation of local y/z about the axis [rad]
  std::optional<Vec3> ecc;    // rigid-arm offset of the section axis (global, both ends)
  friend bool operator==(const BeamSpec& a, const BeamSpec& b) {
    return a.line == b.line && a.section == b.section && a.material == b.material &&
           a.truss == b.truss && a.releases == b.releases && a.roll == b.roll && a.ecc == b.ecc;
  }
};

struct SupportSpec {
  Vec3 point = Vec3::Zero();
  std::array<bool, 6> blocked{true, true, true, true, true, true};
  std::array<std::optional<double>, 6> springs{};  // elastic stiffness per DOF
  std::optional<Line> line;                        // line support (every node on the segment)
  friend bool operator==(const SupportSpec& a, const SupportSpec& b) {
    return a.point == b.point && a.blocked == b.blocked && a.springs == b.springs && a.line == b.line;
  }
};

struct PointMassSpec {
  Vec3 point = Vec3::Zero();
  double mass = 0.0;               // [kg]
  Vec3 inertia = Vec3::Zero();     // rotary inertia about global axes [kg m^2]
  friend bool operator==(const PointMassSpec& a, const PointMassSpec& b) {
    return a.point == b.point && a.mass == b.mass && a.inertia == b.inertia;
  }
};

// Winkler bed on beams (bar-local moduli), attached by Assemble as a property.
struct FoundationSpec {
  std::vector<Line> lines;
  double kfx = 0, kfy = 0, kfz = 0, kmx = 0;
  friend bool operator==(const FoundationSpec&, const FoundationSpec&) = default;
};

// Loads: one alternative per kind; `case_id` groups them into load cases.
struct NodalLoad {
  Vec3 point = Vec3::Zero();
  Vec3 F = Vec3::Zero(), M = Vec3::Zero();
  friend bool operator==(const NodalLoad& a, const NodalLoad& b) {
    return a.point == b.point && a.F == b.F && a.M == b.M;
  }
};
struct UdlLoad {
  std::vector<Line> lines;
  Vec3 w = Vec3::Zero();  // [N/m]
  LoadFrame frame = LoadFrame::Global;
  friend bool operator==(const UdlLoad& a, const UdlLoad& b) {
    return a.lines == b.lines && a.w == b.w && a.frame == b.frame;
  }
};
struct GravityLoad {
  Vec3 g = Vec3(0, 0, -kGravity);
  friend bool operator==(const GravityLoad& a, const GravityLoad& b) { return a.g == b.g; }
};
struct LinearLoad {
  std::vector<Line> lines;
  Vec3 w1 = Vec3::Zero(), w2 = Vec3::Zero();  // at x1 / x2
  LoadFrame frame = LoadFrame::Global;
  double x1 = 0.0, x2 = 1.0;  // normalized range along the authored line
  friend bool operator==(const LinearLoad& a, const LinearLoad& b) {
    return a.lines == b.lines && a.w1 == b.w1 && a.w2 == b.w2 && a.frame == b.frame &&
           a.x1 == b.x1 && a.x2 == b.x2;
  }
};
struct BarPointLoad {
  std::vector<Line> lines;
  double at = 0.5;  // normalized position along the authored line
  Vec3 F = Vec3::Zero();
  LoadFrame frame = LoadFrame::Global;
  friend bool operator==(const BarPointLoad& a, const BarPointLoad& b) {
    return a.lines == b.lines && a.at == b.at && a.F == b.F && a.frame == b.frame;
  }
};
struct ThermalLoad {
  std::vector<Line> lines;
  double dt = 0, gy = 0, gz = 0;
  friend bool operator==(const ThermalLoad&, const ThermalLoad&) = default;
};
struct SettlementLoad {
  Vec3 point = Vec3::Zero();
  std::array<std::optional<double>, 6> d{};
  friend bool operator==(const SettlementLoad& a, const SettlementLoad& b) {
    return a.point == b.point && a.d == b.d;
  }
};
struct CaseLabel {
  std::string name;
  friend bool operator==(const CaseLabel&, const CaseLabel&) = default;
};

// Surface loads on a slab, identified by the slab's OUTLINE (the same clone
// Assemble's PL pin receives, so the match is exact).
struct SurfacePressureLoad {
  std::vector<Vec3> outline;
  Vec3 w = Vec3::Zero();  // [N/m^2]; LocalNormal uses w.z
  SurfaceFrame frame = SurfaceFrame::Global;
  friend bool operator==(const SurfacePressureLoad& a, const SurfacePressureLoad& b) {
    return a.outline == b.outline && a.w == b.w && a.frame == b.frame;
  }
};
struct SurfaceThermalLoad {
  std::vector<Vec3> outline;
  double dt = 0, grad = 0;  // [K], dT/dz [K/m]
  friend bool operator==(const SurfaceThermalLoad&, const SurfaceThermalLoad&) = default;
};

struct LoadSpec {
  using Kind = std::variant<NodalLoad, UdlLoad, GravityLoad, LinearLoad, BarPointLoad, ThermalLoad,
                            SettlementLoad, CaseLabel, SurfacePressureLoad, SurfaceThermalLoad>;
  Kind kind;
  Id case_id = 1;
  friend bool operator==(const LoadSpec&, const LoadSpec&) = default;
};

struct ComboSpec {
  std::optional<std::string> name;
  std::vector<std::pair<double, Id>> terms;  // (factor, case id)
  friend bool operator==(const ComboSpec&, const ComboSpec&) = default;
};

// A named selection: FLAT positions into a tree (every slot counted, nulls
// included, branches in path order — the DataTree::flatten order). Invariant:
// `indices` sorted ascending and deduplicated (use make()).
struct SetSpec {
  std::string name;
  std::vector<std::int64_t> indices;
  static SetSpec make(std::string name, std::vector<std::int64_t> idx);
  friend bool operator==(const SetSpec&, const SetSpec&) = default;
};

// A prestressing tendon before Assemble resolves its chain: the chain lines
// in travel order, the eccentricity profile and the force profile over the
// absolute chain abscissa s in [0, S], S = the sum of the line lengths. The
// eccentricity is in each bar's own section frame.
struct TendonSpec {
  std::vector<Line> lines;
  std::vector<PolySegment> profile;              // piecewise cubic, s-basis
  std::vector<std::pair<double, double>> force;  // (s, P) breakpoints, P >= 0
  bool pretension = false;                       // display / data only
  std::optional<double> ep, ap;                  // staged-tensioning K = 1 losses
  Id case_id = 1;                                // case carrying the equivalent loads
  std::optional<std::string> name;
  friend bool operator==(const TendonSpec& a, const TendonSpec& b) {
    if (a.profile.size() != b.profile.size()) return false;
    for (std::size_t i = 0; i < a.profile.size(); ++i) {
      const PolySegment &x = a.profile[i], &y = b.profile[i];
      if (x.s0 != y.s0 || x.s1 != y.s1 || x.cy != y.cy || x.cz != y.cz) return false;
      if (x.exact != y.exact || x.R != y.R || x.th != y.th || !(x.C == y.C) || !(x.U == y.U) ||
          !(x.D == y.D))
        return false;
      for (std::size_t k = 0; k < 4; ++k)
        if (!(x.p[k] == y.p[k])) return false;
    }
    return a.lines == b.lines && a.force == b.force && a.pretension == b.pretension && a.ep == b.ep &&
           a.ap == b.ap && a.case_id == b.case_id && a.name == b.name;
  }
};

// One construction stage before resolution: bar and support geometry to
// switch, case ids to apply, tendons to jack with their target fraction of
// the declared force.
struct PhaseSpec {
  std::optional<std::string> name;
  std::vector<Line> activate, deactivate;
  std::vector<SupportSpec> add_supports, remove_supports;
  std::vector<Id> cases;
  std::vector<std::pair<TendonSpec, double>> tension;
  friend bool operator==(const PhaseSpec&, const PhaseSpec&) = default;
};

// An elastic-catenary cable before assembly: anchor points, axial rigidity EA
// [N], weight per unstrained length w [N/m], and how the unstrained length is
// given — directly (Length, m) or from a target horizontal tension H / end-2
// tension magnitude T at the undeformed chord (N).
enum class CableRoute { Length, HorizontalTension, EndTension };
struct CableSpec {
  Vec3 a = Vec3::Zero(), b = Vec3::Zero();
  double ea = 0, w = 0;
  CableRoute route = CableRoute::HorizontalTension;
  double value = 0;
  friend bool operator==(const CableSpec&, const CableSpec&) = default;
};

// A meshable slab before Assemble: outline ring (not closed) + hole rings,
// constant thickness [m], material, target mesh edge length [m] and an
// optional result-axes override (a direction).
struct SlabSpec {
  std::vector<Vec3> outline;
  std::vector<std::vector<Vec3>> holes;
  double thickness = 0;
  MaterialSpec material;
  double mesh_size = 0.5;
  std::optional<Vec3> local_x;
  friend bool operator==(const SlabSpec& a, const SlabSpec& b) {
    return a.outline == b.outline && a.holes == b.holes && a.thickness == b.thickness && a.material == b.material &&
           a.mesh_size == b.mesh_size && a.local_x == b.local_x;
  }
};

// A moving-load travel path: ordered, contiguous beam lines (Load Path).
struct PathSpec {
  std::string name;
  std::vector<Line> lines;
  friend bool operator==(const PathSpec&, const PathSpec&) = default;
};

// A convoy load pattern: point axles (offset [m], load [N]) and uniform
// segments (start offset [m], length [m], TOTAL load [N]).
struct ConvoySpec {
  std::string name;
  std::vector<std::pair<double, double>> axles;
  std::vector<std::array<double, 3>> segments;
  friend bool operator==(const ConvoySpec&, const ConvoySpec&) = default;
};

// ---------------- heavy payloads (pointer-identity equality) ----------------

// Section-point stresses riding with static results: the exact per-section
// solves (mesh + unit fields, for a renderer's nodal field) and the
// per-station extremes, parallel to the results.
struct SectionStressBundle {
  std::shared_ptr<const SectionStressEngine> engine;
  StressResults results;
};

// Solved sections persisting across recomputes, shared by copies of the node
// that owns them (content-keyed: a load-only recompute re-solves nothing).
struct SharedSectionCache {
  std::mutex mu;
  SectionEngineCache cache;
};

// The section stresses of static results, built through a shared cache.
// `hooks` (nullable) is polled between sections.
std::shared_ptr<const SectionStressBundle> section_stress_bundle(const Model& model,
                                                                 const StaticResults& results,
                                                                 SharedSectionCache& cache,
                                                                 const SolveHooks* hooks = nullptr);

// Section stresses computed on FIRST ACCESS, once, thread-safe (a throw
// leaves it unset: the next access retries). Lazy because keta is headless:
// a moving-load sweep of a hundred cases would otherwise pay the recovery
// on every recompute nobody reads.
class LazySectionStress {
 public:
  explicit LazySectionStress(std::shared_ptr<SharedSectionCache> cache) : cache_(std::move(cache)) {}
  const SectionStressBundle& get(const Model& model, const StaticResults& results) const {
    std::call_once(once_, [&] { value_ = section_stress_bundle(model, results, *cache_); });
    return *value_;
  }

 private:
  std::shared_ptr<SharedSectionCache> cache_;
  mutable std::once_flag once_;
  mutable std::shared_ptr<const SectionStressBundle> value_;
};

struct AnalyzedModel {
  std::shared_ptr<const Model> model;
  StaticResults results;
  // Appended: set by Analyze and Phased Analysis (null elsewhere — a driver
  // whose results are not linear-static carries no section stresses).
  std::shared_ptr<const LazySectionStress> stress_source;
  // The section-point stresses (computed on first call), or null.
  const SectionStressBundle* stress() const {
    return stress_source ? &stress_source->get(*model, results) : nullptr;
  }
};

// A model paired with its linear-buckling factors and unit-crest modes (the
// Buckling node's B output, the Imperfection node's input).
struct BuckledModel {
  std::shared_ptr<const Model> model;
  BucklingResults results;
};

// A model paired with its natural frequencies and mode shapes (Modal's R).
struct ModalModel {
  std::shared_ptr<const Model> model;
  ModalResults results;
};

// One enveloped source: a load case or a combination, by user id and name
// (empty name = none: labelled LC{id} / C{id}).
struct EnvSource {
  Id id = 0;
  std::string name;
  bool is_combo = false;
};

// An analyzed model's min/max envelope over its combinations (else its cases)
// with the dominant source index behind every extremum (Envelope's E output).
struct EnvelopeModel {
  std::shared_ptr<const Model> model;
  std::optional<std::string> name;
  std::vector<EnvSource> sources;  // indexed by the envelope's *_src entries
  Envelope env;
  // Appended: the analyzed model the envelope was taken from.
  std::shared_ptr<const AnalyzedModel> analyzed;
  // The section-stress envelope over the same sources (computed on first
  // call, with the section stresses themselves), or nullopt when the
  // analyzed model carries none. Throws on a defect (shape mismatch,
  // non-finite stress).
  std::optional<StressEnvelope> stress_env() const {
    const SectionStressBundle* b = analyzed ? analyzed->stress() : nullptr;
    if (!b) return std::nullopt;
    return stress_envelope(b->results, name.value_or(""));
  }
};

using ModelRef = std::shared_ptr<const Model>;
using ResultsRef = std::shared_ptr<const AnalyzedModel>;

// A staged analysis (Phased Analysis's P output): the per-phase engine record
// (tendon jack / loss / effective forces) and its two result views, cases =
// phases (id = 1-based phase index, name = phase name).
struct PhasedModel {
  std::shared_ptr<const Model> model;
  PhasedResults phased;
  ResultsRef cumulated, incremental;
};
using PhasedRef = std::shared_ptr<const PhasedModel>;

// The dimensional family of an influence ordinate (set from the response, so
// a display edge never re-parses the label).
enum class IlDim { Force, Moment, DispTranslation, DispRotation };

// An influence line with its model, response label, load direction and
// ordinate dimension (Influence Line's IL output).
struct InfluenceModel {
  std::shared_ptr<const Model> model;
  InfluenceLine line;
  std::string label;
  Vec3 dir = Vec3(0, 0, -1);
  IlDim dim = IlDim::Moment;
};
using InfluenceRef = std::shared_ptr<const InfluenceModel>;

// A cable model's converged coupled equilibrium (Cable Static's K output).
struct CableStaticModel {
  std::shared_ptr<const Model> model;
  CableStaticResults results;
};
using CableStaticRef = std::shared_ptr<const CableStaticModel>;

// A frequency sweep's steady-state harmonic response (Harmonic's H output).
struct HarmonicModel {
  std::shared_ptr<const Model> model;
  HarmonicResults results;
};
using HarmonicRef = std::shared_ptr<const HarmonicModel>;

// A traced equilibrium path through limit points (Arc Length's T output).
struct ArcTraceModel {
  std::shared_ptr<const Model> model;
  ArcTrace trace;
};
using ArcTraceRef = std::shared_ptr<const ArcTraceModel>;

// A converged large-displacement (corotational) state (Nonlinear's D output).
struct NonlinearModel {
  std::shared_ptr<const Model> model;
  CorotResults results;
};
using NonlinearRef = std::shared_ptr<const NonlinearModel>;

// A random-vibration response (PSD's P output).
struct PsdModel {
  std::shared_ptr<const Model> model;
  PsdResults results;
};
using PsdRef = std::shared_ptr<const PsdModel>;

// A response-spectrum analysis (Spectrum's SP output).
struct SpectrumModel {
  std::shared_ptr<const Model> model;
  SpectrumResults results;
};
using SpectrumRef = std::shared_ptr<const SpectrumModel>;
using BucklingRef = std::shared_ptr<const BuckledModel>;
using ModalRef = std::shared_ptr<const ModalModel>;
using EnvelopeRef = std::shared_ptr<const EnvelopeModel>;

// ---------------- the wire value ----------------

class Value {
 public:
  using Data = std::variant<double, std::int64_t, bool, std::string, Point, Vector, Line, MaterialSpec,
                            SectionSpec, BeamSpec, SupportSpec, PointMassSpec, LoadSpec, ComboSpec,
                            ModelRef, ResultsRef, FoundationSpec, SetSpec, BucklingRef, ModalRef,
                            EnvelopeRef, TendonSpec, PhaseSpec, PhasedRef, PathSpec,
                            ConvoySpec, InfluenceRef, CableSpec, CableStaticRef, HarmonicRef,
                            ArcTraceRef, NonlinearRef, PsdRef, SpectrumRef, SlabSpec>;

  Value() : d_(0.0) {}
  static Value number(double x) { return Value(Data(std::in_place_type<double>, x)); }
  static Value integer(std::int64_t i) { return Value(Data(std::in_place_type<std::int64_t>, i)); }
  static Value boolean(bool b) { return Value(Data(std::in_place_type<bool>, b)); }
  static Value text(std::string s) { return Value(Data(std::in_place_type<std::string>, std::move(s))); }
  template <class T>
  static Value of(T v) {
    return Value(Data(std::in_place_type<T>, std::move(v)));
  }

  // Borrowing accessors (data, get, as_line, as_text, as_set) point INTO the
  // value: lvalues only, a call on a temporary value does not compile.
  const Data& data() const& { return d_; }
  const Data& data() const&& = delete;
  template <class T>
  const T* get() const& {
    return std::get_if<T>(&d_);
  }
  template <class T>
  const T* get() const&& = delete;

  // Coercing accessors: numbers accept ints; ints accept integral numbers;
  // bools accept numbers (non-zero); points and vectors read each other.
  std::optional<double> as_number() const;
  std::optional<std::int64_t> as_int() const;
  std::optional<bool> as_bool() const;
  std::optional<Vec3> as_point() const;
  std::optional<Vec3> as_vector() const { return as_point(); }
  const Line* as_line() const& { return get<Line>(); }
  const Line* as_line() const&& = delete;
  const std::string* as_text() const& { return get<std::string>(); }
  const std::string* as_text() const&& = delete;
  const SetSpec* as_set() const& { return get<SetSpec>(); }
  const SetSpec* as_set() const&& = delete;

  // Wire type tag ("number", "point", "beam", "model", ...).
  std::string type_tag() const;
  // Short human-readable rendering (Panel display, error messages).
  std::string describe() const;

  friend bool operator==(const Value& a, const Value& b) { return a.d_ == b.d_; }

 private:
  explicit Value(Data d) : d_(std::move(d)) {}
  Data d_;
};

}  // namespace keta::param

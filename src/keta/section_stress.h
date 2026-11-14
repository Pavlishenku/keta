// Section-point stresses at every study station: the per-station frame
// internal forces combined with the exact Tri6 warping engine
// (keta/section_exact.h) into fibre-extreme normal, resultant shear and von
// Mises stresses.
//
// Physics. A station carries N, Vy, Vz, T, My, Mz in the bar local frame. The
// pointwise section stress is
//   sigma(y, z) = N/A + (My (Iz z - Iyz y) - Mz (Iy y - Iyz z)) / (Iy Iz - Iyz^2)
//   tau         = St-Venant torsion + transverse-shear flow, recovered by the
//                 warping solve (omega, psi_y, psi_z)
//   sigma_vm    = sqrt(sigma^2 + 3 (tau_y^2 + tau_z^2))
// (tension-positive). It is LINEAR in the six actions, so per distinct
// (section geometry, nu) the warping problem is solved ONCE and six unit-action
// nodal fields are stored; a station's field is their combination, O(nodes)
// with no re-solve. Nodal values are the average of the incident elements'
// values at that node (the transverse-shear field is element-discontinuous).
//
// Candidate set (per section, built once): every boundary node (corners and
// midsides of the edges owned by one element) — sigma is affine in (y, z) so
// its extremes are boundary-bound and sig_min / sig_max are EXACT —, the
// governing nodes of each unit action (sigma min/max, |tau| and vm max, which
// keeps interior transverse-shear peaks) and the |tau| / vm governing node of
// every signed pair of unit actions. tau_max / vm_max of a combined field can
// still sit elsewhere: they are maxima over the candidate set, not bounds.
// Extremes are over mesh NODES.
//
// Geometry. A section's rings are its outline polygons, centroid at the
// origin, the outer one CCW and the holes CW: the closed-form shapes are
// polygonized — circles with 48 points — and a rolled I profile keeps its
// four root fillets (8 chords per quarter circle): keta's catalog props carry
// them, and a fillet-free outline would read the fibre stress of a section
// ~4 % less stiff than the one the analysis used (IPE 300: Iy 7.998e-5 vs
// 8.356e-5 m^4). A ShapePolygon is used as given, centred. The mesh targets a
// number of elements (Standard = 3000), material E = 1 with the bar
// material's nu. A section carrying an attached COMPOSITE exact backend
// (several materials, Section (Composite)) is used as attached: its materials
// keep their own nu.
//
// Fallback. A section without geometry (ShapeGeneric) or whose exact solve
// failed degrades to the closed-form stress points (sigma only, no shear),
// listed with its reason in `StressResults::fallback_sections`; a section
// without stress points degrades further to the uniform N/A.
//
// Units: SI (Pa, m); locations in the section (y, z) centroid frame.
#pragma once
#include <array>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "keta/hooks.h"
#include "keta/section_exact.h"
#include "keta/statics.h"

namespace keta {

// Mesh density of the exact solves: target element count per section.
enum class MeshResolution { Coarse, Standard, Fine };
int target_elements(MeshResolution r);  // 1200 / 3000 / 9000

// Outline rings of a section, centroid at the origin (outer CCW, holes CW);
// empty when the section carries no geometry (ShapeGeneric).
std::vector<std::vector<Vec2>> section_rings(const Section& s);

struct ShearExtreme {
  double tau_max = 0;  // max |tau| = hypot(tau_y, tau_z) [Pa]
  Vec2 tau_max_at = Vec2::Zero();
  double vm_max = 0;   // max von Mises [Pa]
  Vec2 vm_max_at = Vec2::Zero();
};

struct StationStress {
  double sig_min = 0, sig_max = 0;  // algebraic fibre extremes [Pa]
  Vec2 sig_min_at = Vec2::Zero(), sig_max_at = Vec2::Zero();
  std::optional<ShearExtreme> shear;  // exact path only
  bool has_shear() const { return shear.has_value(); }
};

struct BarStress {
  Id bar = 0;
  std::vector<StationStress> stations;  // parallel to BarResult::stations
};

struct CaseStress {
  Id id = 0;
  std::string name;
  std::vector<BarStress> bars;  // parallel to CaseResult::bars
};

struct StressFallback {
  enum class Reason { NoGeometry, SolveFailed };
  Id section = 0;
  Reason reason = Reason::NoGeometry;
  std::string message;  // the exact solve's error (SolveFailed only)
  friend bool operator==(const StressFallback&, const StressFallback&) = default;
};

struct StressResults {
  std::vector<CaseStress> cases;   // parallel to StaticResults::cases
  std::vector<CaseStress> combos;  // parallel to StaticResults::combos
  std::vector<StressFallback> fallback_sections;  // sorted by id
};

// One solved section: mesh, unit-action nodal fields and candidate nodes.
struct SectionStressEntry {
  std::vector<Vec2> nodes;                    // (y, z) centroid frame [m]
  std::vector<std::array<int, 6>> elements;   // Tri6 connectivity
  // Per unit action [N, Vy, Vz, My, Mz, T], per node (sigma, tau_y, tau_z)
  // [Pa per unit SI action].
  std::array<std::vector<Vec3>, 6> unit;
  std::vector<int> candidates;                // sorted node indices
  std::shared_ptr<const ExactSection> solved; // for arbitrary-point sampling
  Vec3 combine(int k, const std::array<double, 6>& a) const;

  // The candidates' coefficients gathered contiguously (structure of arrays)
  // for the per-station scan. The exact engine's sigma depends on N, My, Mz
  // only and its shear on Vy, Vz, T only (verified exactly at build:
  // `separable`), so a candidate costs 9 multiply-adds instead of 18 — and
  // the sum order of the non-zero terms is the general combine's, so the
  // values are bit-identical to it.
  struct CandidateTable {
    std::vector<Vec2> loc;
    std::vector<double> sN, sMy, sMz;           // sigma per unit N, My, Mz
    std::vector<double> tyVy, tyVz, tyT;        // tau_y per unit Vy, Vz, T
    std::vector<double> tzVy, tzVz, tzT;        // tau_z per unit Vy, Vz, T
  } table;
  bool separable = false;
};

// Content-keyed cache of solved sections persisting across recomputes (the
// Analyze node holds one). After a complete build it holds exactly the
// sections that build touched (a dimension slider does not accumulate dead
// solves); an interrupted build evicts nothing.
class SectionEngineCache {
 public:
  struct Stats {
    int builds = 0, reuses = 0;
  };
  const Stats& stats() const { return stats_; }
  std::size_t size() const { return map_.size(); }

 private:
  friend class SectionStressEngine;
  // (ring coordinates, nu, target elements) or the attached backend's address.
  using Ring = std::vector<std::pair<double, double>>;
  using Key = std::tuple<std::vector<Ring>, double, int, const void*>;
  std::map<Key, std::shared_ptr<const SectionStressEntry>> map_;
  Stats stats_;
};

class SectionStressEngine {
 public:
  // Solve every distinct (section, nu) a bar references. `hooks` (nullable)
  // is polled between sections: a cancellation throws Cancelled.
  static SectionStressEngine build(const Model& model, MeshResolution res,
                                   SectionEngineCache& cache, const SolveHooks* hooks = nullptr);
  static SectionStressEngine build(const Model& model, MeshResolution res = MeshResolution::Standard);

  const std::vector<StressFallback>& fallback_sections() const { return fallback_; }
  // Null when the section fell back.
  const SectionStressEntry* entry(Id section, double nu) const;
  // Full nodal field (sigma, tau_y, tau_z, vm) for actions [N, Vy, Vz, My,
  // Mz, T]; empty when the section fell back.
  std::vector<Eigen::Vector4d> nodal_field(Id section, double nu, const std::array<double, 6>& a) const;
  // Exact stress at a section point (centroid frame); nullopt when the section
  // fell back or the point is outside it.
  std::optional<StressState> point_stress(Id section, double nu, const std::array<double, 6>& a,
                                          const Vec2& yz) const;

 private:
  std::map<std::pair<Id, double>, std::shared_ptr<const SectionStressEntry>> entries_;
  std::vector<StressFallback> fallback_;
};

// [N, Vy, Vz, My, Mz, T] of a station.
std::array<double, 6> station_actions(const Station& st);

// `threads`: 0 = std::thread::hardware_concurrency(), 1 = sequential. Every
// (source, bar) is written into its own pre-sized slot, so the result does not
// depend on the thread count. Small recoveries stay sequential.
StressResults recover_stress_results(const Model& model, const StaticResults& results,
                                     const SectionStressEngine& engine, int threads = 0);

// Per-scalar min/max with the dominant source index (combinations when there
// are any, else cases). Throws on an empty source set, a shape mismatch or a
// non-finite stress.
struct StationStressEnv {
  double sig_min = 0, sig_max = 0, tau_max = 0, vm_max = 0;
  int sig_min_src = 0, sig_max_src = 0, tau_max_src = 0, vm_max_src = 0;
  bool has_shear = false;
};
struct BarStressEnv {
  Id bar = 0;
  std::vector<StationStressEnv> stations;
};
struct StressEnvelope {
  std::string name;
  bool sources_are_combos = false;
  std::vector<Id> source_ids;
  std::vector<BarStressEnv> bars;
};
StressEnvelope stress_envelope(const StressResults& stress, const std::string& name = "");

}  // namespace keta

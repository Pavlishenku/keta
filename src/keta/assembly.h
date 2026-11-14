// DOF numbering (support elimination + zero-stiffness auto-fix), sparse
// assembly of K / M / Kg, and per-case load lowering.
#pragma once
#include "keta/element.h"
#include "keta/shell.h"
#include "keta/continuum.h"

namespace keta {

struct CaseLoads;

struct AutoFixedDof {
  Id node = 0;
  int comp = 0;  // 0..5
};

// Free-DOF numbering. Rigidly blocked DOFs (supports) are eliminated; DOFs
// that receive exactly zero stiffness from every incident element and spring
// are auto-fixed (classic case: rotations of pure-truss nodes).
struct DofMap {
  int n_nodes = 0;
  int n_free = 0;
  std::vector<int> free_ix;    // full dof (6*node_ix+c) -> free index, or -1
  std::vector<bool> blocked;   // rigidly blocked by a support
  std::vector<AutoFixedDof> auto_fixed;

  int full(int node_ix, int comp) const { return kDofPerNode * node_ix + comp; }
  int free_of(int node_ix, int comp) const { return free_ix[static_cast<std::size_t>(full(node_ix, comp))]; }
};

// Exact master-slave elimination of the rigid links (model.h RigidLink). The
// assembled matrices and load vectors stay in the EXTENDED free space (slave
// DOFs keep their free index, so every scatter and recovery is unchanged);
// the solve runs on the reduced space u_ext = T u_r + u0:
//   K_r = T^t K T,  f_r = T^t (f - K u0),
// u0 carrying the slave offsets of prescribed (supported, moved) masters.
struct LinkReduction {
  bool active = false;
  int n_red = 0;
  SpMat T;  // n_free x n_red
  // per slave full DOF: its expression over NON-slave full DOFs (blocked
  // masters included — their prescribed values give u0, their reactions
  // collect the slave's constraint force)
  struct Expr {
    int slave = 0;                             // full dof
    std::vector<std::pair<int, double>> terms; // (full dof, coefficient)
  };
  std::vector<Expr> slaves;
  SpMat reduce(const SpMat& A) const;  // T^t A T
  VecX reduce(const VecX& f) const;    // T^t f
  MatX reduce(const MatX& f) const;
  VecX expand(const VecX& ur) const;   // T ur
};

// Elastic link lowered for assembly: node indices, the 12x12 global
// stiffness, and the map from the 12 nodal DOFs to the link-frame relative
// displacement (d, dtheta) at the spring point.
struct LinkElem {
  Id id = 0;
  int ni = 0, nj = 0;
  Eigen::Matrix<double, 6, 12> B = Eigen::Matrix<double, 6, 12>::Zero();
  Mat6 k = Mat6::Zero();
  Mat12 K = Mat12::Zero();  // B^t k B
};

// Cached per-model element set + DOF map — built once, shared by every driver.
struct Assembled {
  std::vector<BarElement> elements;  // parallel to model.bars
  // Flat facet shells, parallel to model.shells (empty in a bar-only model).
  // Their 18x18 stiffness is scattered AFTER the bar loop, in model.shells
  // order — determinism is load-bearing for the factorization cache.
  std::vector<ShellElem> shells;
  DofMap dofs;
  std::vector<double> spring;        // per full-DOF spring stiffness (from supports)
  // ST1 APPUI EL full elastic supports: (node index, 6x6 GLOBAL-frame stiffness),
  // already rotated by the support's ROT alpha. Independent of, and additive to,
  // the per-DOF `spring` diagonal above. See springmass.h.
  std::vector<std::pair<int, Mat6>> spring6;
  // node-to-node links (empty / inactive without links)
  std::vector<LinkElem> links;
  LinkReduction red;
  // continuum elements (G15), model.continua order
  std::vector<ContinuumKernel> continua;
};

// `links_supported`: the caller handles node-to-node links (it scatters the
// elastic links through assemble_stiffness and solves on the reduced space of
// `red`, see FactorCache::set_reduction). A model carrying links is REJECTED
// otherwise — no analysis silently ignores a link.
// `prestress` (per bar axial force, optional): a bar under a non-zero axial
// force is a stiffness source on its nodes' translations for the zero-stiffness
// auto-fix — the geometric stiffness of a prestressed truss or cable net holds
// the transverse DOFs its linear stiffness does not (prestressed modal).
Assembled assemble_elements(const Model& model, bool links_supported = false,
                            const std::vector<double>* prestress = nullptr);
// Elastic links: 12x12 global stiffness of one link (also used by recovery).
LinkElem link_element(const Model& model, const ElasticLink& l);
// Rigid links carrying a second-order term in a geometric stiffness: a link
// with a non-zero arm (master and slave apart), or an elastic link of
// non-zero arm. Linear buckling / prestressed modal analyses reject them.
bool links_have_arms(const Model& model);
// u0 of a lowered case (slave offsets of moved masters) in the extended free
// space, and the matching rhs correction f -= K u0. Zero when no prescribed
// master DOF has slaves.
void fold_link_offsets(const Model& model, const Assembled& asm_, const SpMat& K, CaseLoads& cl);

// Nodal-frame stiffness ACTUALLY assembled for one element: the element's own
// `ke()` for a prismatic bar, and for a bar carrying an exact tapered kernel
// (ST1 `CARA VAR`, see `TaperedBar` in model.h) the same quantity with the
// prismatic member block swapped for the tapered one,
//     T^t (k_tapered + k_soil) T = ke() + T^t (k_tapered - k_local()) T,
// so the rigid-arm transform and the Winkler bed the element already folded in
// are preserved. Bit-identical to `el.ke()` when the model carries no tapered
// bar.
Mat12 bar_ke_assembled(const Model& model, const BarElement& el);

// Reduced (free x free) matrices. Deterministic: serial scatter in bar order,
// facet shells scattered AFTER the bar loop in model.shells order (stiffness
// AND mass — a shell-bearing model has a full shell mass matrix, see shell.h).
// Catenary cables are nonlinear: a cable bar is a typed error unless the
// caller adds the cables itself (`cables_added_by_caller`: solve_cable_static,
// the tangent of a cable state) — a linear analysis never drops them silently.
SpMat assemble_stiffness(const Model& model, const Assembled& asm_, bool cables_added_by_caller = false);
SpMat assemble_mass(const Model& model, const Assembled& asm_, bool lumped = false);
// Geometric stiffness for a per-bar axial force field N (tension-positive) and
// an OPTIONAL per-shell membrane prestress field (nx, ny, nxy per element in
// model.shells order — `shell_membrane_field` in statics.h extracts it from a
// CaseResult). Convention: tangent = K + Kg; compression softens. Passing
// nullptr on a shell-bearing model means "no shell prestress": the shells
// contribute stiffness and mass but no destabilising term (the honest choice
// for the per-bar `_fields` drivers whose signature cannot carry a shell
// state); a non-null field must match model.shells in size.
SpMat assemble_geometric(const Model& model, const Assembled& asm_, const std::vector<double>& N,
                         const std::vector<Vec3>* shell_n = nullptr);

// ---------------- per-case load lowering ----------------

// Span-load description kept per bar for exact station recovery.
struct DistSeg {
  double a = 0, b = 0;  // [m] along the bar
  Vec3 qa = Vec3::Zero(), qb = Vec3::Zero();  // local force intensities [N/m]
  Vec3 ma = Vec3::Zero(), mb = Vec3::Zero();  // local moment intensities [N·m/m]
  // Parabolic part (BarParabolic), the excess at mid-range over the chord:
  // q(u) = qa + (qb - qa) u / l + 4 qc u (l - u) / l^2, u = x - a, l = b - a.
  Vec3 qc = Vec3::Zero();

  // q at u from a (local force intensity).
  Vec3 q_at(double u) const {
    const double l = b - a;
    if (!(l > 0)) return qa;
    const double t = u / l;
    return qa + (qb - qa) * t + 4.0 * t * (1.0 - t) * qc;
  }
};
// A span record with a parabolic segment; the consumers that carry linear
// span loads only refuse it by name (require_linear_span) instead of
// dropping the parabolic part.
bool has_parabolic(const std::vector<DistSeg>& segs);
void require_linear_span(const std::vector<DistSeg>& segs, const std::string& who);
struct PointAt {
  double a = 0;
  Vec3 F = Vec3::Zero(), M = Vec3::Zero();  // local
};
// Concentrated dislocation (ST1 DEFOR CON) at abscissa `a`: an axial gap and
// rotation kinks the deflection recovery must see as jumps of u / v' / w'.
struct KinkAt {
  double a = 0;                      // [m] along the bar
  double dux = 0;                    // axial gap [m]
  double theta_y = 0, theta_z = 0;   // rotation kinks [rad], local axes
  double dvy = 0, dvz = 0, dthx = 0;   // deflection jumps [m] / twist jump [rad]
};
struct BarSpanLoads {
  std::vector<DistSeg> segs;
  std::vector<PointAt> points;
  double eps0 = 0;                  // imposed axial strain at x=0 (thermal + prestrain)
  double kappa_y = 0, kappa_z = 0;  // imposed curvatures at x=0 (thermal gradients)
  // Linear variation of the imposed fields (ST1 DEFOR LIN): value(x) =
  // value_0 + slope * x. Zero slopes keep every uniform path bit-identical.
  double deps0 = 0;                       // d(eps0)/dx [1/m]
  double dkappa_y = 0, dkappa_z = 0;      // d(kappa)/dx [1/m^2]
  std::vector<KinkAt> kinks;              // concentrated dislocations
  // Imposed twist rate d(theta_x)/dx at x = 0 and its slope (appended).
  double kappa_x = 0, dkappa_x = 0;       // [rad/m], [rad/m^2]
  // Quadratic term of the imposed fields (appended; PrestrainParabolic):
  // value(x) = value_0 + slope x + c2 x^2. Zero keeps the linear paths.
  double c2eps0 = 0, c2kappa_y = 0, c2kappa_z = 0, c2kappa_x = 0;
};
// A span record with an imposed twist rate; the consumers that do not carry
// it refuse it by name.
void require_no_imposed_twist(const BarSpanLoads& sp, const std::string& who);
// A span record with a parabolic imposed field (a quadratic term); the
// consumers that carry linear imposed fields only refuse it by name.
bool has_parabolic_strain(const BarSpanLoads& sp);
void require_linear_strain(const BarSpanLoads& sp, const std::string& who);

struct CaseLoads {
  VecX rhs;                                        // free-DOF load vector
  std::vector<std::pair<int, double>> prescribed;  // (full dof, value)
  std::vector<Vec12> f0;                           // per-bar local FEF totals (raw)
  std::vector<BarSpanLoads> span;                  // per-bar span-load records
  VecX nodal_full;                                 // applied nodal loads on ALL dofs (for reactions)
  // Per-shell thermal initial strain / curvature (model.shells order) that
  // recovery must SUBTRACT — otherwise a free thermal expansion reports a
  // spurious resultant. Accumulated as SurfaceThermal loads are lowered.
  std::vector<ShellRecovery> shell_recovery;
  // per-shell consistent nodal loads (element DOF order, summed over the
  // case's shell loads): the element's own share, k_e u - shell_f is its
  // internal nodal force (section cuts, keta/section_cut.h)
  std::vector<Vec18> shell_f;
  // rigid links: slave offsets u0 of prescribed masters (extended free
  // space; empty = none) — added to the solved displacements by recovery
  VecX link_u0;
};

// Lower one load case. Throws on a prescribed DOF that is not rigidly blocked.
CaseLoads lower_case(const Model& model, const Assembled& asm_, const LoadCase& lc);
// The same, the equivalent nodal forces of the bar loads computed on
// `threads` workers (keta/parallel.h; 0 = hardware concurrency) and applied
// in load order: bit-identical to the sequential lowering.
CaseLoads lower_case(const Model& model, const Assembled& asm_, const LoadCase& lc, int threads);
// The building blocks of lower_case: an all-zero CaseLoads sized for `asm_`,
// and the accumulation of ONE load into it (`case_name` only labels errors).
// A bar load touches only its bar's f0 / span record and its 12 DOFs, so a
// caller sweeping many single-load cases can reuse one CaseLoads and reset
// those entries instead of re-lowering a whole case per position.
CaseLoads empty_case_loads(const Assembled& asm_);
void lower_load(const Model& model, const Assembled& asm_, const std::string& case_name,
                const Load& l, CaseLoads& cl);

// Element-implicit fold of prescribed displacements: rhs -= K_fp * u_p.
void fold_prescribed(const Model& model, const Assembled& asm_, CaseLoads& cl);

}  // namespace keta

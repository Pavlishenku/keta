// Static solve of a model whose elastic links carry NON-LINEAR laws
// (model.h LinkLaw: tension or compression only, gaps, elastic-perfectly
// plastic springs, any monotone piecewise-linear curve) — bearings that lift,
// buffers with a clearance, fuse springs.
//
// With monotone laws the equilibrium is the minimum of the convex potential
//     Pi(u) = 1/2 u' K u - f' u + sum_j Psi_j(b_j' u),  Psi_j = int F_j,
// K the linear structure (rigid links eliminated exactly, keta/assembly.h),
// b_j the link-frame component of each law. Newton on Pi with a backtracking
// line search on the energy: piecewise quadratic and convex, so once every
// law sits on its final segment the step is exact — finite termination, the
// equilibrium to round-off, no tolerance on the laws (no penalty, no
// regularisation). Linear analyses reject these links (assemble_elements).
// A law with F(0) != 0 (at-rest earth pressure, preloaded spring) acts from
// the start: at d = 0 its force is already applied to the structure.
// One load case at a time: the response is not superposable.
#pragma once
#include <limits>
#include <memory>

#include "keta/hinge_law.h"
#include "keta/statics.h"

namespace keta {

using HingeLawFamilyRef = std::shared_ptr<const HingeLawFamily>;

struct LinkLawOptions {
  int stations_per_bar = 11;
  int max_iter = 200;
};
struct LinkLawResult {
  CaseResult result;  // links: the law components carry F(d)
  int iterations = 0;
};
LinkLawResult solve_static_link_laws(const Model& model, Id case_id, const LinkLawOptions& opt = {});

// The engine behind solve_static_link_laws, reusable by drivers that change
// the laws between solves (keta/link_history.h): the linear part of the
// model (every law component taken out of its link), its assembly, and the
// law components' deformation operators b_j. `extra` lists further link
// components (link index in model.links, component) handled as law
// components although the model gives them no law; `laws` passed to solve()
// and recover() follow `comps()` order (the model's laws first, in link
// then component order, then `extra`).
class LinkLawCore {
 public:
  struct Comp {
    std::size_t link = 0;
    int c = 0;
    std::vector<std::pair<int, double>> bfree;  // free dof, coefficient
    std::vector<std::pair<int, double>> bfull;  // full dof, coefficient
  };
  // `allow_softening`: laws whose force decreases are accepted (the convex
  // solve() then does not apply: the path-following driver below uses the
  // core); otherwise they are a typed rejection
  LinkLawCore(const Model& model, const std::vector<std::pair<std::size_t, int>>& extra = {},
              bool allow_softening = false);
  const std::vector<Comp>& comps() const { return comps_; }
  int n_free() const { return asm_.dofs.n_free; }
  // a case lowered on the linear part, prescribed displacements and rigid
  // link offsets folded
  CaseLoads lower(const LoadCase& lc) const;
  // the part of each component's deformation set by prescribed DOFs and
  // rigid-link offsets (linear in the loads)
  std::vector<double> d0(const CaseLoads& cl) const;
  // Exact Newton on the convex energy for the right-hand side `rhs`, the
  // offsets `d0` and the laws; `v` is the starting point and the result.
  // Returns the iterations.
  int solve(const VecX& rhs, const std::vector<double>& d0, const std::vector<const LinkLaw*>& laws, VecX& v,
            int max_iter = 200) const;
  double deformation(std::size_t j, const VecX& v, const std::vector<double>& d0) const;
  // the CaseResult of the state v (the laws' forces F_j as nodal actions,
  // reported on their link components)
  CaseResult recover(const CaseLoads& cl, const VecX& v, const std::vector<double>& F, int stations, Id id,
                     const std::string& name) const;
  const SpMat& K() const { return K_; }
  const Assembled& assembled() const { return asm_; }

 private:
  Model w_;
  Assembled asm_;
  SpMat K_;
  std::vector<Comp> comps_;
};

// ---------------- path following through softening ----------------
// Static path of a model whose link laws may SOFTEN (a decreasing branch:
// rotational hinges past their peak, FEMA / ASCE 41 backbones)
// and whose rotational components may follow an M-theta law family of the
// link's own axial force (keta/hinge_law.h): displacement control of one
// DOF (bordered Newton on (u, lambda): r = f_dead + lambda f - K u - sum
// b_j F_j(b_j' u) = 0, u_ctrl = target), exact consistent tangent
// K + sum k_j b_j b_j' + dF/dN k_a b_j b_a' (non-symmetric with a family:
// sparse LU), step halving on failure. A dead case (permanent loads) is
// applied first at lambda = 0. Piecewise-linear laws: each step converges
// once the segments are right, and every point of the path is an exact
// equilibrium (no regularisation of the softening).
struct LinkHingeFamily {
  Id link = 0;
  int comp = 5;        // the rotational component following the family
  HingeLawFamilyRef family;
  int axial_comp = 0;  // the component whose force is N (a law component of the same link)
};
struct LinkPathOptions {
  Id ctrl_node = 0;
  int ctrl_comp = 0;
  std::vector<double> targets;  // piecewise path 0 -> t1 -> t2 ...
  int steps_per_leg = 20;
  int max_iter = 50;
  int max_cuts = 10;
  double tol = 1e-10;
  Id dead_case = 0;  // permanent loads applied before the path (0: none)
  std::vector<LinkHingeFamily> families;
  // Appended: ARC-LENGTH path following instead of displacement control
  // (when arc_steps > 0; `targets` then unused), for SNAP-BACKS that
  // displacement control cannot pass. The path is followed CELL by CELL — a
  // cell fixes each law's segment and, for an M-theta(N) family, its N
  // interval: inside it the forces are smooth (linear for a plain law;
  // for a family whose points move with N, rational in N, with dF/dN in
  // closed form). Each step: the tangent (null vector of [Kt -f], unit in
  // |du|^2 + arc_psi^2 |f|^2 dlam^2), the first boundary of the cell met
  // (linear estimate), and a bordered Newton corrector on the arc's normal
  // plane — or on that boundary's own equation, which lands EXACTLY on the
  // event (a breakpoint, a moving breakpoint theta_i(N), an N grid value);
  // a boundary crossed on the way retargets the step. At an event the cell
  // whose direction is consistent (every boundary at zero opening into it)
  // is kept, the cells continuing the motion first — so a turn sharper than
  // 90 degrees (a snap-back) is followed. Plain laws: one corrector
  // iteration (the path is piecewise linear). A point is recorded every
  // arc_length and at every event (arc_steps points at most); the run stops
  // once u_ctrl passes arc_stop_ctrl (if finite). The first step heads
  // towards increasing u_ctrl (arc_first_sign = -1: decreasing). A
  // family's N leaving its table ends the run (typed).
  int arc_steps = 0;
  double arc_length = 0;
  double arc_psi = 0;
  double arc_stop_ctrl = std::numeric_limits<double>::quiet_NaN();
  int arc_first_sign = 1;
};
struct LinkPathPoint {
  double lambda = 0, u_ctrl = 0;
  std::vector<double> d, F;  // per law component (core order: the model's laws, then the families')
};
struct LinkPathResult {
  std::vector<LinkPathPoint> trace;
  std::vector<Vec6> node_disp;  // final state
  int iterations = 0;
};
LinkPathResult solve_link_laws_path(const Model& model, Id case_id, const LinkPathOptions& opt);

}  // namespace keta

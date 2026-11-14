// Non-linear time integration of small discrete systems with hysteretic
// springs: Newmark-beta (average acceleration by default) with a full
// Newton-Raphson on the TOTAL dynamic equilibrium of every step,
//     M a(n+1) + C v(n+1) + K u(n+1) + f_s(u(n+1)) = p(t(n+1)),
// a and v expressed by the Newmark relations, the consistent tangent of the
// springs in the Jacobian, iterated to a residual of tol x (|p| + |forces|)
// — the converged discrete solution, not a fixed number of corrections
// (Chopra, Dynamics of Structures, §5.7). Springs connect two DOFs (or one
// DOF and the ground) and carry a uniaxial hysteresis law with history:
//   * ElastoPlastic: bilinear, kinematic hardening ratio alpha (0 = EPP),
//     exact return mapping;
//   * Takeda: modified Takeda (Otani 1974) moment-rotation hysteresis on a
//     bilinear skeleton, unloading stiffness k0 (phi_y / phi_max)^alpha,
//     reloading toward the previous peak (beta = 0);
//   * triple friction pendulum (Fenz & Constantinou 2008): only its static
//     backbone under a constant vertical load (tfpb_backbone); no hysteresis
//     law of it is provided.
// Dense matrices: for small SDOF/MDOF systems (typed rejection above
// 2000 DOFs).
#pragma once
#include <functional>
#include <memory>
#include <vector>

#include "keta/core.h"

namespace keta {

class Hysteresis {
 public:
  virtual ~Hysteresis() = default;
  // Force and tangent at the trial deformation d (relative to the committed
  // state's history); no history change until commit().
  virtual void trial(double d, double& f, double& kt) const = 0;
  virtual void commit(double d) = 0;
  virtual std::unique_ptr<Hysteresis> clone() const = 0;
};

class ElastoPlastic : public Hysteresis {
 public:
  ElastoPlastic(double k, double fy, double alpha = 0.0);
  void trial(double d, double& f, double& kt) const override;
  void commit(double d) override;
  std::unique_ptr<Hysteresis> clone() const override { return std::make_unique<ElastoPlastic>(*this); }

 private:
  double k_, fy_, alpha_;
  double dp_ = 0;     // plastic deformation
  double back_ = 0;   // back force (kinematic hardening)
  void state(double d, double& f, double& kt, double& dp, double& back) const;
};

// Modified Takeda (Otani 1974), beta = 0: bilinear skeleton (k0 = My /
// phi_y, kh = (Mu - My) / (phi_u - phi_y), symmetric), unloading stiffness
// k0 (phi_y / phi_m)^alpha (phi_m the largest deformation of the history),
// reloading in straight lines toward the apex of the latest inner loop on
// that side (a reversal point off the skeleton), then toward the outer peak
// of that side (or its yield point if never yielded), then along the skeleton.
// The response is walked segment by segment from the committed state.
class Takeda : public Hysteresis {
 public:
  Takeda(double My, double phi_y, double Mu, double phi_u, double alpha);
  void trial(double d, double& f, double& kt) const override;
  void commit(double d) override;
  std::unique_ptr<Hysteresis> clone() const override { return std::make_unique<Takeda>(*this); }

 private:
  struct Pt {
    double d, f;
  };
  struct State {
    double d = 0, f = 0;
    int last_dir = 0;            // +1 / -1 last motion, 0 virgin
    bool on_skeleton = true;
    double dmax = 0;             // largest |d| so far
    Pt peak[2] = {{0, 0}, {0, 0}};  // outer peaks [0] positive, [1] negative
    bool yielded[2] = {false, false};
    std::vector<Pt> apex[2];     // inner-loop apexes per side, latest last
  };
  double My_, phiy_, Mu_, phiu_, alpha_, k0_, kh_;
  State st_;
  void walk(State& s, double d, double& f, double& kt) const;
  double skeleton(double d, double& kt) const;
};

// Quasi-static force control of one hysteresis law: the deformation at which
// the law (from its committed state) reaches the force f_target, found on
// its monotone piecewise-linear path by bracketing + bisection to round-off;
// committed. Returns the deformation.
double hysteresis_to_force(Hysteresis& h, double d_committed, double f_target);

// Triple friction pendulum bearing (Fenz & Constantinou 2008), backbone of
// the five sliding regimes under the vertical load N, normalised force
// F / N: the inflection points
//   F1 = mu1 (u = 0), F2 = mu2, F3 = mu3, F4 = mu2 + d2 / Re2,
//   F5 = mu3 + d3 / Re3, F6 = mu1 + d1 / Re1 + d3 / Re3 (capacity),
// with the regime stiffnesses 1/(2 Re1), 1/(Re1 + Re2), 1/(Re2 + Re3),
// 1/(Re1 + Re3), 1/(2 Re1). Validity checked (mu1 < mu2 < mu3, Re1 < Re2 <=
// Re3, d1 > (mu3 - mu1) Re1). `tfpb_displacement` is the monotone backbone
// u(F/N) (typed rejection above the capacity F6).
struct TfpbParams {
  double Re1 = 0, Re2 = 0, Re3 = 0, mu1 = 0, mu2 = 0, mu3 = 0, d1 = 0, d2 = 0, d3 = 0;
};
struct TfpbBackbone {
  double F[6] = {0, 0, 0, 0, 0, 0};  // F1..F6 (normalised)
  double u[6] = {0, 0, 0, 0, 0, 0};  // displacements at those points [m]
};
TfpbBackbone tfpb_backbone(const TfpbParams& p);
double tfpb_displacement(const TfpbParams& p, double F_over_N);

struct NlSpring {
  int i = 0, j = -1;  // DOF indices; j = -1: to the ground (force on i = -f)
  std::shared_ptr<Hysteresis> law;
};

struct NlSystem {
  MatX M, C, K;  // linear parts (n x n)
  std::vector<NlSpring> springs;
};

struct NewmarkNlOptions {
  double gamma = 0.5, beta = 0.25;
  double dt = 0.01;
  int steps = 100;
  double tol = 1e-12;
  int max_iterations = 50;
};

struct NlHistory {
  std::vector<double> t;
  std::vector<VecX> u, v, a;
  std::vector<std::vector<double>> spring_force;  // per step, per spring
  int newton_iterations = 0;
};

NlHistory newmark_nonlinear(const NlSystem& sys, const std::function<VecX(double)>& p,
                            const NewmarkNlOptions& opt, const VecX& u0 = VecX(), const VecX& v0 = VecX());

}  // namespace keta

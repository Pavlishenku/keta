// Nonlinear time history of the ASSEMBLED model with devices on its links:
// seismic isolators (lead-rubber bearings, friction pendulums), yielding
// fuses, viscous dampers. Sparse counterpart of keta/nldyn.h.
//
// --- Devices ---
// On one component of an elastic link (keta/model.h ElasticLink, the
// component's linear stiffness replaced): a hysteresis law of the
// component deformation d = b'u (keta/nldyn.h: ElastoPlastic, Takeda, ...)
// and / or a viscous damper F = C |v|^alpha sign(v) of its velocity
// (alpha = 1: linear), continued linearly below |v| < v_lin (the
// power law has an infinite slope at rest). Helpers:
//  * lead-rubber bearing: bilinear, elastic k1, characteristic strength Qd,
//    post-yield k2 — Fy = Qd + k2 dy, dy = Qd / (k1 - k2);
//  * friction pendulum: F = (W / R) d + mu W sign(v), the friction rigid up
//    to a small elastic slip dy: bilinear k1 = mu W / dy + W / R, post-yield
//    W / R.
//
// --- Scheme ---
// Newmark (average acceleration by default) on
//     M a + C v + K u + sum_j b_j F_j(b_j'u, b_j'v) = p(t),
// C = a M + b K (Rayleigh, on the linear structure), Newton on every step to
// tol x (|p| + |forces|) with the consistent tangent (dF/dd + gamma / (beta
// dt) dF/dv on each device), sparse LDLT. Loads: load cases times time
// functions, and a ground acceleration along a global direction (relative
// motion: p = -M r a_g, r the unit translations along it on the free DOFs).
// The discrete energy balance of the average-acceleration scheme holds to
// the solver tolerance and is reported: input work = kinetic + strain
// energy + Rayleigh dissipation + device work.
// Rigid links (master-slave reduction) are a typed rejection. The model's
// link laws (keta/model.h LinkLaw: gaps and contact, compression / tension
// only, multilinear elastic, force at rest) take part as conservative
// devices — the energy balance then holds with their stored energy in the
// device work.
#pragma once
#include <functional>
#include <memory>

#include "keta/nldyn.h"
#include "keta/statics.h"

namespace keta {

struct LinkDevice {
  Id link = 0;
  int comp = 0;
  std::shared_ptr<Hysteresis> law;  // null: none
  double C = 0, alpha = 1, v_lin = 1e-4;
};
std::shared_ptr<Hysteresis> lrb_law(double k1, double Qd, double k2);
std::shared_ptr<Hysteresis> fps_law(double W, double R, double mu, double dy = 1e-3);
// a link law (gap, compression / tension only, multilinear elastic, force at
// rest) as a conservative device; the model's own link laws are added
// automatically by newmark_model (after the caller's devices, in link and
// component order)
std::shared_ptr<Hysteresis> elastic_law(const LinkLaw& l);

struct NlExcitation {
  std::vector<std::pair<Id, std::function<double(double)>>> cases;  // load case x f(t)
  int ground_dir = -1;                                               // 0..2, -1: none
  std::function<double(double)> ground_accel;
};

struct ModelNlDynOptions {
  double dt = 0.01;
  int steps = 100;
  double gamma = 0.5, beta = 0.25;
  double tol = 1e-10;
  int max_iterations = 50;
  double rayleigh_a = 0, rayleigh_b = 0;
  bool lumped_mass = false;
  std::vector<std::pair<Id, int>> observe;  // (node, dof)
};

struct ModelNlDynResult {
  std::vector<double> t;
  std::vector<std::vector<double>> u;              // per step: observed DOFs
  std::vector<std::vector<double>> device_force;   // per step: per device
  std::vector<std::vector<double>> device_deform;  // per step: per device
  std::vector<double> device_work;                 // per device, sum of 1/2 (F_n + F_n+1) (d_n+1 - d_n)
  double input_work = 0, kinetic = 0, strain = 0, rayleigh_work = 0;  // final energy terms
  int newton_iterations = 0;
};

ModelNlDynResult newmark_model(const Model& model, const std::vector<LinkDevice>& devices, const NlExcitation& load,
                               const ModelNlDynOptions& opt);

}  // namespace keta

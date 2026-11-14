// Railway bridge dynamics (EN 1991-2 §6.4.6): the HSLM-A universal trains,
// and the speed sweep of a train crossing the bridge — modal superposition,
// envelopes of displacements and accelerations per speed, the static
// envelope of the same train, and the dynamic factor.
//
// --- HSLM-A (EN 1991-2 Figure 6.12, Table 6.3) ---
// Train k = 1..10: N intermediate coaches of length D, bogie axle spacing d,
// axle load P (Table 6.3). Axles, from the head: a power car 0 / 3 / 14 /
// 17 m, then 3.525 m to the end-coach bogie, N + 3 bogies (axle pairs d
// apart) at the coach pitch D — articulated coaches share their bogies — then
// 3.525 m to the mirrored power car: 2N + 14 axles. The axle offsets are
// <= 0 (the train behind its head), loads in N.
//
// --- Speed sweep (rail_speed_sweep) ---
// The train enters the path (a chain of bars, its head at the path start at
// t = 0) and runs at constant speed v until its last axle has left, plus a
// free-vibration time. Modal superposition on the given modes (M-normalised,
// keta/eigsolve.h), modal damping (DampingSpec). EXACT ingredients:
//  * modal forces: for a unit force at abscissa x of a bar, the consistent
//    nodal load is the element's shape functions at x, so its projection on
//    a mode is a cubic in x (Hermite bending, linear axial): four samples per
//    bar and their Lagrange interpolant give it exactly;
//  * each modal oscillator is advanced by the exact recurrence of a damped
//    linear oscillator under a load linear over the step (Nigam & Jennings;
//    Chopra, Dynamics of Structures, Table 5.2.1); the time step only
//    samples the moving load (dt = T_min / steps_per_period, and at most
//    max_travel / v).
// Results per speed and observed DOF: max / min displacement, max |accel|.
// Static envelope: the influence line of each observed displacement under
// the same train (keta/influence.h). Dynamic factor = max over speeds of
// max |u_dyn| / max |u_stat| (EN 1991-2 §6.4.6.5 phi'_dyn).
#pragma once
#include "keta/dynamics.h"
#include "keta/influence.h"

namespace keta {

ConvoyPattern rail_hslm_a(int train);

struct RailSweepOptions {
  std::vector<double> speeds;               // [m/s]
  DampingSpec damping;                      // modal
  std::vector<std::pair<Id, int>> observe;  // (node, dof 0..5)
  double free_time = 0;                     // after the last axle leaves [s]
  int steps_per_period = 20;                // on the highest mode used
  double max_travel = 0.5;                  // [m] of train travel per step at most
};

struct RailSweepPoint {
  double speed = 0;
  std::vector<double> u_max, u_min, a_max;  // per observed DOF
};

struct RailSweepResult {
  std::vector<RailSweepPoint> points;
  std::vector<double> static_max, static_min;  // per observed DOF
  std::vector<double> dynamic_factor;          // max |u_dyn| / max |u_stat|
  std::vector<double> critical_speed;          // the speed of that maximum
};

RailSweepResult rail_speed_sweep(const Model& model, const ModalResults& modal, const std::vector<Id>& path_bars,
                                 const ConvoyPattern& train, const Vec3& dir, const RailSweepOptions& opt);

// --- Vehicle-bridge interaction (rail_speed_sweep_vbi) ---
// Each axle carries an unsprung mass m_w riding on the deck (no
// separation) and a sprung mass m_s on a suspension (k, c) (a sprung-mass
// system). With w(s, t) = sum phi_n(s) q_n(t) the deck
// displacement along `dir` and the axle at s = v t + offset, the wheel
// follows the deck with the full kinematics of a moving point:
//   y_w = phi' q,  y_w' = phi' q' + v dphi' q,
//   y_w'' = phi' q'' + 2 v dphi' q' + v^2 ddphi' q
// (phi, its first and second derivatives along the path from the same
// exact cubic per bar), and the contact force on the deck is
//   F = (m_w + m_s) g + k (y_s - y_w) + c (y_s' - y_w') - m_w y_w''.
// Bridge modes and sprung masses form one linear time-varying system,
// advanced by the average-acceleration Newmark scheme (second order in dt;
// it tends to the moving-force limit as the sprung masses vanish). The
// static envelope is the trains' static weights.
struct SprungAxle {
  double offset = 0;   // <= 0, behind the head
  double m_wheel = 0;  // unsprung [kg]
  double m_body = 0;   // sprung [kg]
  double k = 0, c = 0; // suspension
};
struct VehicleTrain {
  std::vector<SprungAxle> axles;
  double g = 9.81;
};
RailSweepResult rail_speed_sweep_vbi(const Model& model, const ModalResults& modal, const std::vector<Id>& path_bars,
                                     const VehicleTrain& train, const Vec3& dir, const RailSweepOptions& opt);

// (oscillator_step, the exact step of a damped linear oscillator, is in
// keta/dynamics.h.)

}  // namespace keta

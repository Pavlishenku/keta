// Footbridges under pedestrian loading — the Sétra / AFGC 2006 methodology
// ("Passerelles piétonnes : évaluation du comportement vibratoire sous
// l'action des piétons", guide méthodologique, Sétra/AFGC, 2006; English
// edition "Footbridges: assessment of vibrational behaviour of footbridges
// under pedestrian loading", Sétra, October 2006) and two dynamics helpers
// that go with it (a periodic Fourier load factor, a high-frequency modal
// damping filter).
//
// The guide's chain, step by step (chapter 2):
//   * footbridge class I to IV (traffic) -> crowd density d of each load case
//     (table of 2.4.1: III 0.5, II 0.8, I 1.0 pedestrian/m2);
//   * natural frequency -> frequency range 1 to 4 (tables 2.3 / 2.4) ->
//     load cases to run (table 2.5);
//   * per load case and mode: a harmonic SURFACE load
//       d x F0 x cos(2 pi f t) x n_eq / N x psi(f)          [N/m2]
//     with F0 = 280 / 140 / 35 N (vertical / longitudinal / lateral, first
//     harmonic), 70 / 35 / 7 N (second harmonic), N = d S the number of
//     pedestrians, n_eq = 10.8 sqrt(xi N) (random frequencies and phases,
//     cases 1 and class II case 3) or 1.85 sqrt(N) (random phases, case 2 and
//     class I case 3), psi(f) the trapezoidal reduction of figures 2.3 / 2.4;
//   * the load follows the sign of the mode shape (every contribution
//     unfavourable, 2.4.1.1) and drives the mode at resonance;
//   * the peak acceleration falls in a comfort range 1 to 4 (tables 2.1 /
//     2.2); laterally it is also kept below 0.10 m/s2 (no lock-in).
// Computation only: the ranges are returned as numbers, never as a verdict.
#pragma once
#include <functional>
#include <vector>

#include "keta/eigsolve.h"

namespace keta {

enum class FootbridgeClass { I = 1, II = 2, III = 3, IV = 4 };
// Direction of the vibration (and of the pedestrian force).
enum class PedestrianDirection { Vertical, Longitudinal, Lateral };

// Equivalent number of pedestrians acting in phase at the natural frequency:
// 10.8 sqrt(xi n) for random frequencies and phases (sparse or dense crowd,
// d < 1 P/m2), 1.85 sqrt(n) for random phases only (very dense crowd).
double setra_equivalent_pedestrians(double n, double xi, bool very_dense);

// Frequency range 1 (maximum risk) to 4 (negligible) of a natural frequency
// [Hz]. Vertical and longitudinal: 1 = [1.7, 2.1], 2 = [1, 1.7) U (2.1, 2.6],
// 3 = (2.6, 5], 4 = [0, 1) U (5, inf). Lateral: 1 = [0.5, 1.1],
// 2 = [0.3, 0.5) U (1.1, 1.3], 3 = (1.3, 2.5], 4 = the rest. A frequency on a
// boundary takes the range of higher risk.
int setra_frequency_range(double f_hz, PedestrianDirection dir);

// Comfort range 1 (maximum comfort) to 4 (unacceptable) of a peak
// acceleration [m/s2]. Vertical: <= 0.5, <= 1, <= 2.5, above. Horizontal
// (longitudinal or lateral): <= 0.15, <= 0.3, <= 0.8, above.
int setra_comfort_range(double accel, PedestrianDirection dir);
// Lateral acceleration above which forced synchronisation (lock-in) may
// start: kept below 0.10 m/s2 whatever the comfort level (table 2.2).
inline constexpr double kSetraLockInAccel = 0.10;

// Reduction factor psi(f) of the guide (figure 2.3 for the first harmonic,
// figure 2.4 for the second), a trapezoid in f [Hz]:
//   harmonic 1, vertical / longitudinal: 0 -> 1 on [1, 1.7], 1 on [1.7, 2.1],
//     1 -> 0 on [2.1, 2.6]; lateral: [0.3, 0.5], [0.5, 1.1], [1.1, 1.3];
//   harmonic 2, vertical / longitudinal: [2.6, 3.4], [3.4, 4.2], [4.2, 5];
//     lateral: [1.3, 1.7], [1.7, 2.1], [2.1, 2.5].
// harmonic is 1 or 2.
double setra_psi(double f_hz, PedestrianDirection dir, int harmonic = 1);

// Amplitude F0 [N] of the force of one pedestrian: first harmonic 280 / 140 /
// 35 N, second harmonic 70 / 35 / 7 N (vertical / longitudinal / lateral).
double setra_pedestrian_force(PedestrianDirection dir, int harmonic = 1);

// Crowd density d [pedestrians/m2] of load case 1, 2 or 3 for a class
// (case 1: III 0.5, II 0.8; case 2: I 1.0; case 3: II 0.8, I 1.0). Throws
// when the case does not apply to the class (e.g. case 2 for class II).
double setra_crowd_density(FootbridgeClass cls, int load_case);

// Load cases to run for a mode in frequency range `range` (table 2.5):
// class III {1} in range 1; class II {1} in ranges 1-2, {3} in range 3;
// class I {2} in ranges 1-2, {3} in range 3; empty otherwise (class IV, range 4).
std::vector<int> setra_load_cases(FootbridgeClass cls, int range);

// Mass of the pedestrians [kg/m2] to add to the deck: 70 kg per pedestrian
// (700 N), d x 70 (guide 2.3; the frequencies are bracketed between the empty
// deck and d = 1).
inline double setra_pedestrian_mass(double density) { return 70.0 * density; }

// Amplitude [N/m2] of the harmonic surface load of a load case:
//   d F0 n_eq / N psi(f),  N = d S,
// n_eq of 10.8 sqrt(xi N) (case 1, class II case 3) or 1.85 sqrt(N) (case 2,
// class I case 3); case 3 takes the second-harmonic force and psi. `area`
// [m2] is the loaded deck area S, `xi` the structural damping ratio.
struct SetraSurfaceLoad {
  double density = 0;        // d [P/m2]
  double pedestrians = 0;    // N = d S
  double equivalent = 0;     // n_eq
  double force = 0;          // F0 [N]
  double psi = 0;            // reduction factor
  double amplitude = 0;      // d F0 n_eq / N psi [N/m2]
  double load_frequency = 0; // frequency of the load = the natural frequency f [Hz]
};
SetraSurfaceLoad setra_surface_load(FootbridgeClass cls, int load_case, PedestrianDirection dir, double f_hz,
                                    double xi, double area);

// ---------------- the load that follows a mode shape ----------------

// A harmonic line load on deck bars whose sign follows the mode shape.
struct PedestrianModeLoad {
  std::vector<Id> bars;        // loaded bars (empty = every Frame bar)
  double q = 0;                // amplitude [N/m] = surface load x loaded width
  int comp = 2;                // global direction 0, 1, 2 (x, y, z)
  // false (the guide): q sign(phi) on every bar, every contribution
  // unfavourable — the maximal modal force. true: q only where phi > 0, zero
  // elsewhere (a constant-direction partial load).
  bool sign_cancelled = false;
};

// The load case of the pedestrian load on mode `mode` of `modal` (solved on
// `model`): uniform global loads q s(x) e_comp, s(x) the sign of the
// component `comp` of the mode shape along each bar — the element field of
// the unit-crest shape (ModalResults::modes): axial linear, transverse
// Hermite cubics of the end displacements and rotations, releases
// re-injected. A bar is cut at the zeros of that cubic (found exactly, one
// bisection per monotonic piece), so a mode node inside a bar is honoured.
// Feed it to solve_harmonic (steady state at the resonance), to
// time_history_modal / time_integrate through pedestrian_modal_force, or to
// statics.
LoadCase pedestrian_mode_load_case(const Model& model, const ModalResults& modal, int mode,
                                   const PedestrianModeLoad& load, Id case_id = 0);

// F(step) = F_case cos(2 pi f step dt), free-DOF vector of the load case above
// (lowered on `model`), the time function of time_integrate /
// time_history_modal.
std::function<VecX(int)> pedestrian_modal_force(const Model& model, const ModalResults& modal, int mode,
                                                const PedestrianModeLoad& load, double f_hz, double dt);

// ---------------- helpers ----------------

// Periodic load factor a0 + sum_i [a_i cos(2 pi i t / T) + b_i sin(2 pi i t / T)],
// a_cos[0] / b_sin[0] the first harmonic (i = 1).
double fourier_factor(double t, double a0, const std::vector<double>& a_cos, const std::vector<double>& b_sin,
                      double period);
// The factor sampled at t_n = n dt, n = 0..nsteps.
std::vector<double> fourier_load_factors(double a0, const std::vector<double>& a_cos,
                                         const std::vector<double>& b_sin, double period, double dt, int nsteps);

// High-frequency modal damping filter: xi_i for f_i <= f_c, xi_i + d_c above
// (over-damps the modes beyond the cut frequency). For the analyses that take
// per-mode damping (modal superposition, spectra, harmonic response); the
// direct integration takes Rayleigh damping only. freqs and xi of one length.
std::vector<double> high_frequency_filtered_damping(const std::vector<double>& freqs_hz,
                                                    const std::vector<double>& xi, double f_c, double d_c);

}  // namespace keta

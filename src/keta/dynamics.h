// Dynamics: direct time integration (Newmark-beta / HHT-alpha), modal
// superposition time history, steady-state harmonic response (complex modal
// FRF), random-vibration PSD response, and seismic response spectrum
// (EN 1998-1 elastic/design + tabulated) with SRSS / CQC (Der Kiureghian,
// Rosenblueth-Elishakoff for unequal damping) / ten-percent / double-sum
// combination rules, rigid-mode gating and missing-mass correction.
#pragma once
#include <map>

#include "keta/eigsolve.h"

namespace keta {

// ---------------- damping ----------------

struct DampingSpec {
  enum class Kind { Modal, PerMode, Rayleigh };
  Kind kind = Kind::Modal;
  double zeta = 0.05;                // Modal: one ratio for all modes
  std::vector<double> per_mode;      // PerMode: lowest-frequency first
  double alpha = 0, beta = 0;        // Rayleigh: C = alpha M + beta K

  double zeta_for(double omega, std::size_t mode_index) const;
  static DampingSpec rayleigh_two_mode(double om1, double zeta1, double om2, double zeta2);
};

// ---------------- direct time integration ----------------

struct IntegrationScheme {
  double gamma = 0.5, beta = 0.25;  // Newmark average acceleration
  double hht_alpha = 0.0;           // HHT-alpha in [-1/3, 0]; 0 = plain Newmark
  static IntegrationScheme newmark_average() { return {}; }
  static IntegrationScheme newmark_linear() { return {0.5, 1.0 / 6.0, 0.0}; }
  static IntegrationScheme hht(double alpha);
};

struct TimeHistoryResults {
  double dt = 0;
  std::vector<VecX> u, v, a;    // free-DOF fields per step (including t=0)
  std::vector<double> energy;   // 1/2 v'Mv + 1/2 u'Ku
};

// force(step) returns the free-DOF load vector at time step*dt (step 0..n).
TimeHistoryResults time_integrate(const Model& model, int n_steps, double dt,
                                  const std::function<VecX(int)>& force,
                                  const DampingSpec& damping = {},
                                  const IntegrationScheme& scheme = {}, bool lumped_mass = false);

// The same from given initial conditions u0, v0 on the free DOFs (extended
// numbering of assemble_elements(model, true); rigid links refused): the
// initial acceleration from M a0 = F(0) - C v0 - K u0. time_integrate is
// this from rest.
TimeHistoryResults time_integrate_from(const Model& model, int n_steps, double dt,
                                       const std::function<VecX(int)>& force, const VecX& u0, const VecX& v0,
                                       const DampingSpec& damping = {}, const IntegrationScheme& scheme = {},
                                       bool lumped_mass = false);

// Modal-superposition time history (same scheme reduced to scalar oscillators).
TimeHistoryResults time_history_modal(const Model& model, const ModalResults& modal, int n_steps,
                                      double dt, const std::function<VecX(int)>& force,
                                      const DampingSpec& damping = {},
                                      const IntegrationScheme& scheme = {});

// Excitation builders.
// Uniform base acceleration: F(t) = -M r_dir * accel[step].
// support_coupling: the consistent mass coupling to the moving supports is
// added to M r (ground_coupling, ModalOptions::support_coupling; on by
// default, the exact load of a consistent mass).
std::function<VecX(int)> seismic_force(const Model& model, int dir,
                                       const std::vector<double>& accel, bool lumped_mass = false,
                                       bool support_coupling = true);
// F(t) = P * sin(omega t) on one nodal DOF.
std::function<VecX(int)> harmonic_force(const Model& model, Id node, int comp, double amplitude,
                                        double omega, double dt);

// ---------------- harmonic steady state ----------------

struct HarmonicResults {
  std::vector<double> freqs_hz;
  // [freq][node] complex displacement, split re/im.
  std::vector<std::vector<Vec6>> disp_re, disp_im;
  std::vector<std::vector<Vec6>> disp_amp;  // per-component amplitude
  // COMPLEX station diagrams (appended; filled only when `stations > 0`):
  // [freq][bar] internal forces + deflections recovered from the real and
  // imaginary displacement parts by the exact STATIC station machinery — the
  // recovery operator is linear, the applied span loads ride the REAL
  // (in-phase) part. `bars_amp` holds the per-component amplitude
  // hypot(re, im); the phase is atan2(im, re) per component. Documented
  // approximation: the element's DISTRIBUTED INERTIA is not re-spread along
  // the span (exactly the consistent-mass discretization already made by the
  // modal basis) — station equilibrium uses the static span loads, so the
  // diagrams are exact up to the mesh, not beyond it.
  std::vector<std::vector<BarResult>> bars_re, bars_im, bars_amp;
  // The superposed modes (appended): circular frequencies [rad/s] and the
  // damping ratio each one took, lowest first.
  std::vector<double> omegas, zetas;
};

// F(t) = Re(P e^{i omega t}) with P from the load case; exact modal FRF
// q_i = (phi_i' P) / (om_i^2 - om^2 + 2 i zeta_i om_i om).
// `stations` > 0 additionally recovers complex station diagrams on every bar
// (that many stations per bar). `hooks` (keta/hooks.h, nullable): polled
// before the modal solve and at every sweep frequency ("frequency i/n").
HarmonicResults solve_harmonic(const Model& model, Id case_id, const std::vector<double>& freqs_hz,
                               const DampingSpec& damping = {}, int n_modes = 0,
                               int stations = 0, const SolveHooks* hooks = nullptr);

// ---------------- random vibration (PSD) ----------------

struct PsdResults {
  std::vector<double> freqs_hz;
  std::vector<std::vector<Vec6>> disp_psd;  // response PSD per node/component
  std::vector<Vec6> rms;                    // per node
  // Station-level RMS (appended; filled only when `stations > 0`): per bar,
  // stations whose N/Vy/Vz/T/My/Mz/defl hold the RMS of that quantity —
  // sqrt(int |H_q(f)|^2 S_alpha df) by the same linearity as the nodal PSD
  // (H_q the complex station transfer from the harmonic sweep).
  std::vector<BarResult> rms_bars;
};

// Fully-correlated load pattern F(t) = alpha(t) P with input PSD S_alpha(f)
// [unit^2/Hz]: S_y(f) = |H_y(f)|^2 S_alpha(f), RMS = sqrt(int S_y df).
// `stations` > 0 additionally integrates station-force RMS on every bar.
PsdResults solve_psd(const Model& model, Id case_id, const std::vector<double>& freqs_hz,
                     const std::vector<double>& s_alpha, const DampingSpec& damping = {},
                     int n_modes = 0, int stations = 0, const SolveHooks* hooks = nullptr);

// ---------------- response spectrum ----------------

struct ResponseSpectrum {
  enum class Kind { Table, Ec8Elastic, Ec8Design };
  Kind kind = Kind::Table;
  std::vector<std::pair<double, double>> table;  // (period s, Sa m/s^2), clamped ends
  // EC8 parameters:
  double ag = 0;      // design ground acceleration [m/s^2]
  double S = 1.0;     // soil factor
  double TB = 0.15, TC = 0.4, TD = 2.0;
  double q = 1.0;     // behaviour factor (design spectrum)
  double beta_floor = 0.2;  // lower bound factor (design spectrum)
  // Appended: the plateau amplification (2.5; 3.0 for an EN 1998-1 vertical
  // spectrum) and a FIXED damping correction (> 0; 0 = eta(zeta) =
  // sqrt(10 / (5 + 100 zeta)) >= 0.55, the EN 1998-1 default).
  double plateau = 2.5;
  double eta_fixed = 0.0;
  // Appended: a TABLE given at the damping table_zeta (> 0) is corrected to
  // the modal damping by eta(zeta) / eta(table_zeta) (EN 1998-1 eta), for
  // periods >= table_correction_from (0: everywhere; the zero-period
  // acceleration is commonly left uncorrected).
  double table_zeta = 0.0;
  double table_correction_from = 0.0;

  double sa(double period, double zeta) const;
  static ResponseSpectrum ec8_elastic(double ag, double S, double TB, double TC, double TD);
  static ResponseSpectrum ec8_design(double ag, double S, double TB, double TC, double TD,
                                     double q);
};

// Modal combination rules: SRSS, CQC (Der Kiureghian), the 10 % rule
// (US NRC RG 1.92), the double sum (Rosenblueth); appended: ABS (sum of the
// absolute modal values, an upper bound) and NRL (Naval Research Laboratory:
// the largest |value| plus the SRSS of the others).
enum class ComboRule { Srss, Cqc, TenPercent, DoubleSum, Abs, Nrl };

struct SpectrumOptions {
  int dir = 2;  // 0..5 (2 = vertical Z, 0/1 horizontal)
  ComboRule rule = ComboRule::Cqc;
  int n_modes = 10;
  bool missing_mass = false;  // add residual static pseudo-mode at the ZPA
  double duration = 20.0;     // DoubleSum strong-motion duration t_d [s]
  // Appended: > 0 recovers per-mode station diagrams on every bar (that many
  // stations per bar) and combines them with the same rule.
  int stations = 0;
  // Appended: the per-bar recovery abscissae [m] (the study sections), taking
  // precedence over `stations` where given — results then combine with the
  // static ones station by station.
  const std::map<Id, std::vector<double>>* stations_abs = nullptr;
  // Appended: participation with the consistent-mass coupling to the moving
  // supports (ModalOptions::support_coupling; on by default, the exact
  // earthquake load of a consistent mass).
  bool support_coupling = true;
  // Appended: the lumped (diagonal) mass for the modes (ModalOptions::lumped_mass)
  bool lumped_mass = false;
  // Appended: excitation along a unit vector of translations (any
  // horizontal direction, an inclined one); zero = `dir`. Gamma_v = v .
  // (Gamma_x, Gamma_y, Gamma_z), exactly (the participation is linear in the
  // influence vector); missing mass sum v_c^2 m_c - sum Gamma_v^2.
  Vec6 dir_vector = Vec6::Zero();
};

struct SpectrumResults {
  std::vector<double> frequencies_hz;
  std::vector<double> zetas;
  std::vector<double> participation;         // Gamma_i for the chosen direction
  std::vector<std::vector<Vec6>> modal_disp; // signed per-mode peak fields
  std::vector<Vec6> combined_disp;           // componentwise >= 0
  std::vector<double> modal_base_shear;
  double base_shear = 0;
  double missing_mass_shear = 0;
  // Station diagrams (appended; filled only when opt.stations > 0):
  // `modal_bars[i]` are the SIGNED per-mode peak diagrams recovered from the
  // signed modal field Gamma_i Sd_i phi_i (load-free static recovery), and
  // `combined_bars` applies the combination rule per station and per
  // component, componentwise >= 0 — exactly the convention of the nodal
  // `modal_disp` / `combined_disp` pair (the sign information lives in the
  // per-mode diagrams, the dominant mode's being the physical shape). The
  // missing-mass correction applies to the base shear only, not to the
  // station diagrams (documented).
  std::vector<std::vector<BarResult>> modal_bars;
  std::vector<BarResult> combined_bars;
  // Support reactions of the same load-free recoveries (appended, filled with
  // the station diagrams): signed per mode, combined componentwise >= 0.
  std::vector<std::vector<Reaction>> modal_reactions;
  std::vector<Reaction> combined_reactions;
};

SpectrumResults solve_spectrum(const Model& model, const ResponseSpectrum& spectrum,
                               const DampingSpec& damping = {}, const SpectrumOptions& opt = {});

// Directional combination: SRSS or the EN 1998-1 30% rule over per-direction
// combined displacement fields.
std::vector<Vec6> combine_directional_srss(const std::vector<std::vector<Vec6>>& dirs);
std::vector<Vec6> combine_directional_30pct(const std::vector<std::vector<Vec6>>& dirs);

// The exact step of a damped linear oscillator u'' + 2 zeta w u' + w^2 u =
// p(t), p linear from p0 to p1 over dt (closed-form particular + homogeneous
// solution; w = 0: a free mass).
void oscillator_step(double w, double zeta, double dt, double p0, double p1, double& u, double& v);

// Pseudo-acceleration response spectrum Sa(T) = w^2 max |u| of an
// acceleration history (linear between samples), each oscillator advanced
// by the exact step: exact for the piecewise-linear history.
std::vector<double> response_spectrum_of(const std::vector<double>& accel, double dt, const std::vector<double>& periods,
                                         double zeta);

// Floor response spectrum: the ABSOLUTE acceleration history of (node,
// comp) under a ground acceleration along global translation `dir`, by
// modal superposition on the given modes (each modal oscillator advanced
// exactly for the piecewise-linear ground motion; the floor acceleration is
// sampled at the same dt), then its response spectrum at zeta_floor.
struct FloorSpectrum {
  std::vector<double> accel;  // absolute floor acceleration per step
  std::vector<double> periods, sa;
};
FloorSpectrum floor_response_spectrum(const Model& model, const ModalResults& modal, Id node, int comp, int dir,
                                      const std::vector<double>& ground_accel, double dt, const DampingSpec& damping,
                                      const std::vector<double>& periods, double zeta_floor);

// Composite modal damping of a structure of several materials (the
// strain-energy weighting used for non-proportional damping, e.g. steel and
// concrete, EN 1998-2): zeta_i = sum_e zeta_e
// E_e,i / sum_e E_e,i, E_e,i = 1/2 phi_i' K_e phi_i over bars, shells and
// elastic links (zeta_links). A per-mode damping spec.
DampingSpec strain_energy_damping(const Model& model, const ModalResults& modal, const std::map<Id, double>& zeta_of_material,
                                  double zeta_links = 0.0);

// Cross-correlation coefficients (exposed for direct use).
double cqc_rho(double om_i, double om_j, double zi, double zj);
double double_sum_rho(double om_i, double om_j, double zi, double zj, double td);

}  // namespace keta

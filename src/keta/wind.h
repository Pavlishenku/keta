// Wind engineering:
//  * EN 1991-1-4 mean-wind profile (log law, terrain roughness, turbulence
//    intensity, peak velocity pressure) and per-bar static drag/lift/moment
//    loads emitted as an ordinary LoadCase;
//  * vortex shedding: critical speed n b / St, Scruton number, and the
//    EN 1991-1-4 Annex E lock-in amplitude estimate;
//  * spectral buffeting: modal quasi-steady response with the Kaimal or von
//    Karman u-spectrum, Davenport coherence joint acceptance, mechanical
//    admittance |H|^2, RMS + Davenport peak factor;
//  * classical 2-DOF (heave-pitch) flutter of a flat-plate deck section by
//    the p-k method with full Theodorsen aerodynamics (Jones' approximation
//    of C(k)), general elastic-axis offset a and static unbalance x_alpha;
//    the critical speed is the first zero crossing of the branch damping.
#pragma once
#include "keta/dynamics.h"

namespace keta {

// ---------------- mean wind (EN 1991-1-4) ----------------

struct WindProfile {
  double vb = 26.0;    // basic wind speed [m/s]
  double z0 = 0.05;    // roughness length [m] (terrain II = 0.05)
  double zmin = 2.0;
  double rho = 1.225;  // air density [kg/m^3]
  double co = 1.0;     // orography factor
  double kI = 1.0;

  double kr() const { return 0.19 * std::pow(z0 / 0.05, 0.07); }
  double vm(double z) const;   // mean speed
  double iv(double z) const;   // turbulence intensity
  double qp(double z) const;   // peak velocity pressure [Pa]
};

struct AeroSection {
  double B = 1.0;               // reference width [m]
  double cd = 1.0, cl = 0.0, cm = 0.0;   // static coefficients
  double dcl = 0.0, dcm = 0.0;  // incidence slopes [1/rad] (aeroelastic terms)
};

// Static wind load case: drag along `wind_dir` (horizontal unit vector),
// lift vertical, per-bar from the peak pressure at the bar mid height.
LoadCase mean_wind_case(const Model& model, const WindProfile& wind, const AeroSection& aero,
                        const Vec3& wind_dir, Id case_id, const std::string& name);

// ---------------- vortex shedding ----------------

double vortex_critical_speed(double natural_freq, double width, double strouhal = 0.18);
double scruton_number(double mass_per_length, double damping_ratio, double rho, double width);
// EN 1991-1-4 Annex E amplitude estimate: y_max / b.
double vortex_amplitude_ratio(double strouhal, double scruton, double K = 0.13, double Kw = 0.6,
                              double clat = 0.7);

// ---------------- spectral buffeting ----------------

struct BuffetingOptions {
  double U = 30.0;       // mean speed at deck level [m/s]
  double Iu = 0.15;      // along-wind turbulence intensity
  double Lu = 100.0;     // integral length scale [m]
  double rho = 1.225;
  double coherence_c = 10.0;  // Davenport decay coefficient
  bool von_karman = false;    // false = EN 1991-1-4 Annex B (Kaimal type), true = von Karman
  int n_freq = 200;           // log grid points
  double gust_duration = 600.0;
};

// One-sided u-spectrum S_u(n) [ (m/s)^2 / Hz ], sigma_u = Iu U, f_L = n Lu / U:
// EN 1991-1-4 (B.2) 6.8 f_L / (1 + 10.2 f_L)^(5/3), or von Karman
// 4 f_L / (1 + 70.8 f_L^2)^(5/6); both integrate to sigma_u^2 over n.
double u_spectrum(double n, const BuffetingOptions& o);

struct BuffetingModeResult {
  double sigma_q = 0;   // modal RMS (unit modal mass)
  double peak_factor = 0;
  double nu = 0;        // up-crossing rate
};
struct BuffetingResults {
  std::vector<BuffetingModeResult> modes;
  std::vector<Vec6> rms_disp;   // nodal RMS field
  std::vector<Vec6> peak_disp;  // g * rms
};

// Quasi-steady drag buffeting of the given modes: per-bar loading strip with
// the drag direction from `wind_dir`, joint acceptance with Davenport
// coherence over the bar-midpoint distances.
BuffetingResults spectral_buffeting(const Model& model, const ModalResults& modal,
                                    const AeroSection& aero, const Vec3& wind_dir,
                                    const DampingSpec& damping, const BuffetingOptions& opt);

// ---------------- flutter ----------------

struct FlutterSection {
  double m = 0;       // mass per unit span [kg/m]
  double I = 0;       // pitch inertia per unit span about the elastic axis [kg m]
  double S = 0;       // static unbalance m*x_alpha*b [kg] (per span)
  double b = 0;       // HALF chord [m]
  double a = 0;       // elastic-axis offset aft of midchord, in half-chords
  double omega_h = 0; // uncoupled heave frequency [rad/s]
  double omega_a = 0; // uncoupled pitch frequency [rad/s]
  double zeta_h = 0, zeta_a = 0;
  double rho = 1.225;
};

struct FlutterResult {
  double U_flutter = 0;      // critical speed [m/s]
  double omega_flutter = 0;  // flutter frequency [rad/s]
  bool found = false;
};

// Theodorsen circulation function (Jones approximation), k = omega b / U.
std::complex<double> theodorsen_c(double k);

// p-k sweep over U in [u_min, u_max]; first damping zero crossing bisected.
FlutterResult flutter_speed(const FlutterSection& s, double u_min, double u_max, int n_sweep = 60);

}  // namespace keta

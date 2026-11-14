// Transient 2-D heat conduction in a bar cross-section exposed to fire.
//
// Physics
// -------
// rho(T) c(T) dT/dt = div( lambda(T) grad T )                       (EN 1992-1-2 §3.3)
// solved on a uniform Cartesian grid of square cells of side MAILLE by a
// cell-centred finite-volume scheme (identical to the classic 5-point FDM
// Laplacian in the interior) integrated in time by BACKWARD EULER.
//
// The material properties are LAGGED: lambda, rho and c are evaluated at the
// temperature of the PREVIOUS step. Two consequences, both deliberate:
//   * the linear system stays symmetric positive definite and its SPARSITY
//     PATTERN is constant for the whole run, so Eigen's SimplicialLDLT is
//     `analyzePattern`-ed once and only re-`factorize`-d each step;
//   * the scheme is first-order in time in the nonlinear properties. Backward
//     Euler is unconditionally stable, so unlike ST1's explicit scheme
//     (manual annex B5 p381, which is CFL-bound and therefore forced into very
//     short steps on fine meshes) the step size is an accuracy choice only.
//     Default dt = 30 s.
//
// Surface (Robin) boundary condition, EN 1991-1-2 §3.1:
//   q_net = alpha_c (T_g - T_s) + eps sigma ((T_g+273.15)^4 - (T_s+273.15)^4)
// The radiative part is written exactly as an equivalent linear coefficient
//   h_rad = eps sigma (Tg_K^2 + Ts_K^2)(Tg_K + Ts_K)   so that
//   q_net = (alpha_c + h_rad)(T_g - T_s) = h_surf (T_g - T_s),
// and the wall temperature T_s is then the root of the quartic obtained by
// writing the same flux through the series network
//   surface film  ->  passive protection (PROTECTION res, m^2K/W)  ->  half cell
// exactly as ST1 does (manual annex B5 p382, eq. 5 + 6). The COMPOSITE SURFACE
// CONDUCTANCE that reaches the first cell centre is therefore
//
//   U_bnd = 1 / ( 1/h_surf + res + (MAILLE/2)/lambda_cell )      [W/m^2K]
//
// with T_s recovered from T_s = T_g - q/h_surf and fixed-point-iterated (on the
// previous step's cell temperature) before the implicit solve. `res` is treated
// as constant, as ST1 does and documents (manual p216).
//
// Boundary kinds (ST1 `EXPOSITION FEU`, manual p216-217):
//   CN / HC / HCM   standard, hydrocarbon and French "majoree" gas curves,
//                   delegated to keta's ExposureCurve (fire.h)
//   ADIABATIQUE     zero flux
//   TEMP t          constant gas temperature
//   LOI_GAZ  table  piecewise-linear gas temperature (minutes, degC)
//   LOI_PAROI table piecewise-linear WALL temperature -> Dirichlet through the
//                   protection + half-cell resistance
//
// Geometries (ST1 `SECTION`, manual p218-222):
//   DALLE      1-D through the thickness, faces SUP(=VY) / INF(=WY)
//   RECTANGLE  2-D, four faces VY/WY (+z/-z) and VZ/WZ (+y/-y)
//   PERIMETRE  arbitrary polygon rasterised by point-in-polygon; each boundary
//              face takes the exposure of the NEAREST polygon segment
//
// Units are SI and METRES ARE MANDATORY (manual p214: "L'unite de longueur de
// reference est le metre : il est indispensable de travailler en metre pour
// faire les calculs au feu"). Temperatures are in degrees Celsius everywhere in
// the API; the +273.15 offset is applied only inside the radiation term.
//
// Known limitations (never silent):
//   * no concrete SPALLING model — as ST1 (manual p219), the user subtracts the
//     spalled thickness from the section geometry;
//   * no moisture transport: the water content only enters through the EN
//     1992-1-2 §3.3.2 specific-heat peak;
//   * the section is assumed prismatic and the conduction along the bar axis is
//     neglected (2-D section problem), which is the EN 1992-1-2 annex B setting.
#pragma once
#include <functional>

#include "keta/fire.h"

namespace keta {

// ---------------------------------------------------------------------------
// Thermal material — ST1 `MAT i ... E FEU EC2 ...` (manual p214-215)
// ---------------------------------------------------------------------------
struct ThermalMaterial {
  enum class Conductivity { Constant, Ec2 };  // ST1 CONDUCTIVITE conduc | EC2
  enum class DensityLaw { Ec2, Constant };    // ST1 RHO_FEU rho <EC2, CTE>

  Aggregate aggregate = Aggregate::Siliceous;  // GRANULAT SILICE / CALCAIRE
  double uu = 0.0;         // UU: water content, weight of water / weight of concrete (0..0.03)
  Conductivity cond = Conductivity::Ec2;
  double lambda_const = 1.6;  // CONDUCTIVITE conduc [W/mK] (simplified analysis)
  // CONDUCTIVITE EC2 alpha: lambda = alpha*lambda_lower + (1-alpha)*lambda_upper,
  // the two limits of EN 1992-1-2 §3.3.3. ST1 default alpha = 1 (lower limit).
  // ST1's named shortcuts (manual p214-215) map onto this ratio as:
  //   AN NORMAL (traditional concretes, French NA §3.3.3) -> alpha = 1 (LOWER)
  //   AN BHP    (high-performance concretes, NA §6.3)     -> alpha = 0 (UPPER)
  double lambda_alpha = 1.0;
  double emissivity = 0.7;  // EMISSIVITE (EN 1992-1-2 §2.2 default 0.7)
  double rho20 = 2300.0;    // RHO_FEU [kg/m^3] — fire density, NOT the weight density
  DensityLaw rho_law = DensityLaw::Ec2;
  double c_const = 0.0;  // CHALEUR_SPECIFIQUE_MASSIQUE [J/kgK]; > 0 => constant

  // EN 1992-1-2 §3.3.3 thermal conductivity [W/mK].
  double lambda(double theta_c) const;
  // EN 1992-1-2 §3.3.2 specific heat [J/kgK], with the moisture peak between
  // 100 and 115 C: c_peak = 900 (u = 0), 1470 (u = 1.5 %), 2020 (u = 3 %),
  // linear in between, then linear back onto the dry curve at 200 C.
  double cp(double theta_c) const;
  // EN 1992-1-2 §3.3.2(3) density ratio applied to rho20 [kg/m^3].
  double rho(double theta_c) const;
  double rho_cp(double theta_c) const { return rho(theta_c) * cp(theta_c); }
  // Appended — user tables (degC, value), piecewise linear, held constant
  // beyond their ends: conductivity [W/mK] and specific heat [J/kgK]
  // replacing the laws above when given; and the VOLUMETRIC ENTHALPY H(T)
  // [J/m^3] (latent heat of a phase change included), which replaces
  // rho c: each step then solves H(T_new) - H(T_old) = dt div(lambda grad T)
  // by Newton (tangent capacity dH/dT) safeguarded by an exact line search
  // (the residual is the gradient of a convex functional) — energy conserved
  // exactly, whatever the step crosses (beyond its ends H continues with its
  // end slopes).
  std::vector<std::pair<double, double>> lambda_table, cp_table, enthalpy_table;
  double enthalpy(double theta_c) const;
  // the capacity of a step from t0 to t1 [J/m^3K]: the enthalpy chord, or
  // rho c at t0
  double step_capacity(double t0, double t1) const;

  static ThermalMaterial ec2_an_normal();  // CONDUCTIVITE EC2 AN NORMAL (alpha = 1)
  static ThermalMaterial ec2_an_bhp();     // CONDUCTIVITE EC2 AN BHP    (alpha = 0)
  // Constant-property material: makes the module usable for any homogeneous
  // solid (steel plates agree with EN 1993-1-2).
  static ThermalMaterial constant_props(double lambda, double rho, double cp,
                                        double emissivity = 0.7);
};

// ---------------------------------------------------------------------------
// Exposure — ST1 `EXPOSITION FEU i ... FIN` (manual p216-217)
// ---------------------------------------------------------------------------
struct Exposure {
  enum class Kind {
    Standard,     // CN  — ISO 834
    Hydrocarbon,  // HC
    Hcm,          // HCM — hydrocarbure majoree (CETU)
    Adiabatic,    // ADIABATIQUE
    ConstantGas,  // TEMP temp
    GasTable,     // LOI_GAZ   (minutes, degC)
    WallTable     // LOI_PAROI (minutes, degC) — prescribed WALL temperature
  };
  Kind kind = Kind::Standard;
  double temp = 20.0;                            // TEMP
  std::vector<std::pair<double, double>> table;  // LOI_GAZ / LOI_PAROI
  double protection = 0.0;   // PROTECTION res [m^2K/W], series surface resistance
  bool user_convection = false;  // CONVECTION alphac (else CONVECTION EC1)
  double alpha_c_user = 25.0;

  bool adiabatic() const { return kind == Kind::Adiabatic; }
  bool dirichlet() const { return kind == Kind::WallTable; }
  // Gas (or, for LOI_PAROI, wall) temperature [degC] at t minutes. Tables are
  // held constant beyond their extreme dates (manual p216).
  double gas_temp(double t_minutes) const;
  // Convection coefficient [W/m^2K]. CONVECTION EC1 defaults follow
  // EN 1991-1-2 §3.2 (25 for the standard curve, 50 for hydrocarbon curves);
  // for TEMP, ST1 switches to the unexposed-surface value 4 W/m^2K while the
  // wall stays below 50 C (manual p217).
  double alpha_c(double t_wall_c) const;
  // Appended: a convection coefficient [W/m^2K] as a function of the wall
  // temperature [degC] (piecewise linear, constant beyond the ends),
  // replacing the rules above when given (a film coefficient h(T)).
  std::vector<std::pair<double, double>> convection_table;
};

// ---------------------------------------------------------------------------
// Section geometry — ST1 `SECTION is ... FIN` (manual p218-222)
// ---------------------------------------------------------------------------
// Face naming, keta's reading of manual p218: the DALLE faces named SUP(VY) and
// INF(WY) are the upper/lower ones, i.e. the faces whose outward normal is
// +z / -z in the section local frame (the thickness EPAISSEUR runs along z, the
// width LARGEUR along y). By extension the RECTANGLE's lateral faces VZ / WZ are
// the +y / -y ones.
struct FireSection {
  enum class Kind { Slab, Rectangle, Perimeter };
  Kind kind = Kind::Slab;

  double thickness = 0.0;  // EPAISSEUR [m] — along local z
  double width = 1.0;      // LARGEUR   [m] — along local y (DALLE default 1 m)

  // Exposure indices into the `exposures` vector passed to solve_heat2d;
  // a negative index means "unexposed" and is treated as adiabatic.
  int expo_vy = -1;  // +z (SUP)
  int expo_wy = -1;  // -z (INF)
  int expo_vz = -1;  // +y
  int expo_wz = -1;  // -y

  // PERIMETRE: polygon vertices stored as (y, z, 0) in the section local frame,
  // direct or indirect order (manual p221), and one exposure index per segment
  // [i, i+1] (the last closing [n-1, 0]).
  std::vector<Vec3> polygon;
  std::vector<int> seg_exposure;

  double maille = 0.01;  // MAILLE [m] — cell size; ST1 default 0.01 m

  // Smallest dimension the mesh must resolve (manual p389: the mesh size must
  // be at least three times smaller than it).
  double min_dimension() const;
  void validate() const;  // typed rejects: units, degenerate geometry, coarse mesh
};

// ---------------------------------------------------------------------------
// Grid + results
// ---------------------------------------------------------------------------
struct HeatCell {
  double y = 0, z = 0;  // cell centre in the section local frame [m]
  double area = 0;      // cell area [m^2] (per unit length along the bar)
};

struct HeatGrid {
  // Cell pitch. hz is the MAILLE snapped so an integer number of cells spans the
  // section exactly; hy equals hz except for a DALLE, which is meshed 1-D
  // through its thickness and where hy is the full LARGEUR (one cell wide).
  double hy = 0, hz = 0;
  int ny = 0, nz = 0;           // bounding-box cell counts
  double y0 = 0, z0 = 0;        // centre of cell (0, 0)
  std::vector<int> index;       // iz*ny + iy -> active cell index, or -1
  std::vector<HeatCell> cells;  // active cells, row-major in (iz, iy)

  int at(int iy, int iz) const;
  // Index of the active cell containing (y, z); -1 when outside the section.
  int locate(double y, double z) const;
};

struct HeatResult {
  HeatGrid grid;
  std::vector<double> times;         // output instants [minutes], ascending
  std::vector<std::vector<double>> T;  // [snapshot][cell] temperature [degC]

  int snap(int i) const;  // resolves -1 to the last snapshot; throws if out of range
  const std::vector<double>& field(int i = -1) const;
  // Temperature at a point of the section [degC]. Throws when the point falls
  // outside the meshed section.
  double T_at(double y, double z, int i = -1) const;
  double mean_temp(int i = -1) const;  // area-weighted
  double max_temp(int i = -1) const;
};

struct HeatOptions {
  double duration = 60.0;         // FEU DUREE [minutes]
  double t_init = 20.0;           // FEU TEMP  [degC] — initial structure temperature
  double dt = 30.0;               // time step [s]
  std::vector<double> snapshots;  // extra output instants [minutes]; `duration` always output
  int max_cells = 2000000;        // typed reject above this
  // Appended: an initial temperature FIELD (y, z) -> degC replacing t_init.
  std::function<double(double, double)> t_initial;
};

// Solve the transient section heat problem. `exposures` is the ST1 exposure
// table; `section` indexes into it.
HeatResult solve_heat2d(const FireSection& section, const ThermalMaterial& mat,
                        const std::vector<Exposure>& exposures, const HeatOptions& opt = {});

}  // namespace keta

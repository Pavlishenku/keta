#include "keta/heat2d.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace keta {
namespace {

constexpr double kSigmaSB = 5.67e-8;  // Stefan-Boltzmann [W/m^2K^4]
constexpr double kAbsZero = 273.15;   // degC -> K offset used in the radiation term

double clamp_theta(double t) { return std::min(std::max(t, 20.0), 1200.0); }

// Linear interpolation of the EN 1992-1-2 §3.3.2(2) specific-heat peak values
// c_peak(u) = 900 (u = 0), 1470 (u = 1.5 %), 2020 (u = 3 %).
double cp_peak(double u) {
  double x = std::min(std::max(u, 0.0), 0.03);
  if (x <= 0.015) return 900.0 + (1470.0 - 900.0) * (x / 0.015);
  return 1470.0 + (2020.0 - 1470.0) * ((x - 0.015) / 0.015);
}

}  // namespace

// ---------------------------------------------------------------------------
// ThermalMaterial
// ---------------------------------------------------------------------------

namespace {
// piecewise-linear table, held constant beyond its ends
double table_at(const std::vector<std::pair<double, double>>& t, double x) {
  if (x <= t.front().first) return t.front().second;
  if (x >= t.back().first) return t.back().second;
  for (std::size_t i = 1; i < t.size(); ++i)
    if (x <= t[i].first)
      return t[i - 1].second + (t[i].second - t[i - 1].second) * (x - t[i - 1].first) / (t[i].first - t[i - 1].first);
  return t.back().second;
}
void check_table(const std::vector<std::pair<double, double>>& t, const char* what) {
  for (std::size_t i = 0; i < t.size(); ++i) {
    require(std::isfinite(t[i].first) && std::isfinite(t[i].second), std::string("heat2d: non-finite ") + what + " table");
    if (i > 0) require(t[i].first > t[i - 1].first, std::string("heat2d: ") + what + " table: temperatures must increase");
  }
}
}  // namespace

double ThermalMaterial::enthalpy(double x) const {
  const auto& t = enthalpy_table;
  require(t.size() >= 2, "heat2d: the enthalpy table needs two points");
  auto seg = [&](std::size_t i) { return (t[i + 1].second - t[i].second) / (t[i + 1].first - t[i].first); };
  if (x <= t.front().first) return t.front().second + seg(0) * (x - t.front().first);
  if (x >= t.back().first) return t.back().second + seg(t.size() - 2) * (x - t.back().first);
  for (std::size_t i = 1; i < t.size(); ++i)
    if (x <= t[i].first) return t[i - 1].second + seg(i - 1) * (x - t[i - 1].first);
  return t.back().second;
}

double ThermalMaterial::step_capacity(double t0, double t1) const {
  if (enthalpy_table.empty()) return rho_cp(t0);
  if (std::fabs(t1 - t0) > 1e-12 * (1 + std::fabs(t0))) return (enthalpy(t1) - enthalpy(t0)) / (t1 - t0);
  const double d = 1e-6 * (1 + std::fabs(t0));  // the slope at t0
  return (enthalpy(t0 + d) - enthalpy(t0 - d)) / (2 * d);
}

double ThermalMaterial::lambda(double theta) const {
  if (!lambda_table.empty()) return table_at(lambda_table, theta);
  if (cond == Conductivity::Constant) {
    require(lambda_const > 0, "heat2d: constant conductivity must be > 0");
    return lambda_const;
  }
  double x = clamp_theta(theta) / 100.0;
  double upper = 2.0 - 0.2451 * x + 0.0107 * x * x;   // EN 1992-1-2 §3.3.3 eq. (3.11)
  double lower = 1.36 - 0.136 * x + 0.0057 * x * x;   // EN 1992-1-2 §3.3.3 eq. (3.12)
  double a = std::min(std::max(lambda_alpha, 0.0), 1.0);
  return a * lower + (1.0 - a) * upper;
}

double ThermalMaterial::cp(double theta) const {
  if (!cp_table.empty()) return table_at(cp_table, theta);
  if (c_const > 0) return c_const;
  double t = clamp_theta(theta);
  // Dry curve, EN 1992-1-2 §3.3.2(1).
  double dry;
  if (t <= 100.0)
    dry = 900.0;
  else if (t <= 200.0)
    dry = 900.0 + (t - 100.0);
  else if (t <= 400.0)
    dry = 1000.0 + (t - 200.0) / 2.0;
  else
    dry = 1100.0;
  double peak = cp_peak(uu);
  if (peak <= 900.0 || t <= 100.0 || t > 200.0) return dry;
  if (t <= 115.0) return peak;                       // moisture plateau 100..115 C
  return peak + (1000.0 - peak) * (t - 115.0) / 85.0;  // back onto the dry curve at 200 C
}

double ThermalMaterial::rho(double theta) const {
  require(rho20 > 0, "heat2d: fire density RHO_FEU must be > 0");
  if (rho_law == DensityLaw::Constant) return rho20;
  double t = clamp_theta(theta);  // EN 1992-1-2 §3.3.2(3)
  if (t <= 115.0) return rho20;
  if (t <= 200.0) return rho20 * (1.0 - 0.02 * (t - 115.0) / 85.0);
  if (t <= 400.0) return rho20 * (0.98 - 0.03 * (t - 200.0) / 200.0);
  return rho20 * (0.95 - 0.07 * (t - 400.0) / 800.0);
}

ThermalMaterial ThermalMaterial::ec2_an_normal() {
  ThermalMaterial m;
  m.cond = Conductivity::Ec2;
  m.lambda_alpha = 1.0;
  return m;
}

ThermalMaterial ThermalMaterial::ec2_an_bhp() {
  ThermalMaterial m;
  m.cond = Conductivity::Ec2;
  m.lambda_alpha = 0.0;
  return m;
}

ThermalMaterial ThermalMaterial::constant_props(double lam, double rho_, double c, double em) {
  require(lam > 0 && rho_ > 0 && c > 0, "heat2d: constant properties must be > 0");
  ThermalMaterial m;
  m.cond = Conductivity::Constant;
  m.lambda_const = lam;
  m.rho_law = DensityLaw::Constant;
  m.rho20 = rho_;
  m.c_const = c;
  m.emissivity = em;
  return m;
}

// ---------------------------------------------------------------------------
// Exposure
// ---------------------------------------------------------------------------

double Exposure::gas_temp(double t) const {
  ExposureCurve c;
  switch (kind) {
    case Kind::Standard: c.kind = ExposureCurve::Kind::Iso834; return c.gas_temp(t);
    case Kind::Hydrocarbon: c.kind = ExposureCurve::Kind::Hydrocarbon; return c.gas_temp(t);
    case Kind::Hcm: c.kind = ExposureCurve::Kind::Hcm; return c.gas_temp(t);
    case Kind::GasTable:
    case Kind::WallTable:
      require(!table.empty(), "heat2d: LOI_GAZ / LOI_PAROI needs at least one (t, T) pair");
      c.kind = ExposureCurve::Kind::Table;
      c.table = table;
      return c.gas_temp(t);
    case Kind::ConstantGas:
    case Kind::Adiabatic: return temp;
  }
  return temp;
}

double Exposure::alpha_c(double t_wall) const {
  if (!convection_table.empty()) return table_at(convection_table, t_wall);
  if (user_convection) {
    require(alpha_c_user >= 0, "heat2d: CONVECTION alphac must be >= 0");
    return alpha_c_user;
  }
  switch (kind) {
    case Kind::Standard: return 25.0;  // EN 1991-1-2 §3.2.1
    // EN 1991-1-2 §3.2.2 gives 50 W/m^2K for the hydrocarbon curve; the ST1
    // manual's table on p217 lists 25 for HC, which contradicts both EN 1991-1-2
    // and ST1's own HCM entry. keta follows the standard (and keta's own
    // ExposureCurve::alpha_c in fire.cpp): 50 for both HC and HCM.
    case Kind::Hydrocarbon:
    case Kind::Hcm: return 50.0;
    case Kind::GasTable:
    case Kind::WallTable: return 50.0;  // manual p217: "programme par extension"
    case Kind::ConstantGas: return t_wall < 50.0 ? 4.0 : 25.0;  // manual p217
    case Kind::Adiabatic: return 0.0;
  }
  return 25.0;
}

// ---------------------------------------------------------------------------
// FireSection
// ---------------------------------------------------------------------------

namespace {

// Squared distance from (py, pz) to the segment [a, b] in the (y, z) plane.
double dist2_to_segment(double py, double pz, const Vec3& a, const Vec3& b) {
  double ay = a.x(), az = a.y(), by = b.x(), bz = b.y();
  double vy = by - ay, vz = bz - az;
  double len2 = vy * vy + vz * vz;
  double t = len2 > 0 ? ((py - ay) * vy + (pz - az) * vz) / len2 : 0.0;
  t = std::min(std::max(t, 0.0), 1.0);
  double dy = py - (ay + t * vy), dz = pz - (az + t * vz);
  return dy * dy + dz * dz;
}

// Ray-casting point-in-polygon in the (y, z) plane; orientation-agnostic
// (ST1 allows the perimeter in either direction, manual p221).
bool point_in_polygon(const std::vector<Vec3>& poly, double py, double pz) {
  bool inside = false;
  std::size_t n = poly.size();
  for (std::size_t i = 0, j = n - 1; i < n; j = i++) {
    double yi = poly[i].x(), zi = poly[i].y();
    double yj = poly[j].x(), zj = poly[j].y();
    if ((zi > pz) != (zj > pz)) {
      double yx = yi + (pz - zi) * (yj - yi) / (zj - zi);
      if (py < yx) inside = !inside;
    }
  }
  return inside;
}

}  // namespace

double FireSection::min_dimension() const {
  switch (kind) {
    // A DALLE is meshed through its thickness only (1-D), so the width does not
    // constrain the mesh.
    case Kind::Slab: return thickness;
    case Kind::Rectangle: return std::min(thickness, width);
    case Kind::Perimeter: {
      // The manual's rule (p389) bears on the smallest DRAWING dimension —
      // the thinnest plate of the outline — not on the shortest polyline
      // segment: official decks approximate fillets with millimetre-scale
      // segments (exemple10's IPN) that must not drive the mesh guard. The
      // bounding box is the robust proxy keta can compute for an arbitrary
      // outline; sub-maille segments are simply absorbed by rasterization.
      double xmin = std::numeric_limits<double>::infinity(), xmax = -xmin;
      double ymin = xmin, ymax = -xmin;
      for (const Vec3& p : polygon) {
        xmin = std::min(xmin, p.x());
        xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y());
        ymax = std::max(ymax, p.y());
      }
      return std::min(xmax - xmin, ymax - ymin);
    }
  }
  return 0.0;
}

void FireSection::validate() const {
  require(maille > 0 && finite(maille), "fire section: MAILLE must be a finite positive length");
  // Metres are mandatory (manual p214). A "mesh" of 10 mm entered as 10 (i.e. in
  // millimetres) or a 300 mm web entered as 300 lands far outside any real
  // cross-section, so guard the plausible metric range instead of silently
  // meshing a 300 m beam.
  require(maille <= 1.0, "fire section: MAILLE " + std::to_string(maille) +
                             " m is implausible — fire input must be in METRES (ST1 manual p214)");
  if (kind == Kind::Perimeter) {
    require(polygon.size() >= 3, "fire section: PERIMETRE needs at least 3 vertices");
    require(seg_exposure.size() == polygon.size(),
            "fire section: PERIMETRE needs one EXPOSITION per segment (" +
                std::to_string(polygon.size()) + " expected, " +
                std::to_string(seg_exposure.size()) + " given)");
    for (const Vec3& p : polygon)
      require(finite(p), "fire section: non-finite PERIMETRE vertex");
  } else {
    require(thickness > 0 && finite(thickness), "fire section: EPAISSEUR must be > 0");
    require(width > 0 && finite(width), "fire section: LARGEUR must be > 0");
    require(thickness <= 100.0 && width <= 100.0,
            "fire section: dimensions look non-metric — fire input must be in METRES "
            "(ST1 manual p214)");
  }
  double dmin = min_dimension();
  require(dmin > 0 && finite(dmin), "fire section: degenerate geometry (zero smallest dimension)");
  // Manual p389: "il faut une taille de maille trois fois inferieure a la
  // dimension la plus petite du dessin".
  require(maille <= dmin / 3.0,
          "fire section: MAILLE " + std::to_string(maille) +
              " m is too coarse — it must be at least 3x smaller than the smallest "
              "section dimension (" + std::to_string(dmin) + " m, i.e. MAILLE <= " +
              std::to_string(dmin / 3.0) + " m; ST1 manual p389)");
}

// ---------------------------------------------------------------------------
// HeatGrid / HeatResult
// ---------------------------------------------------------------------------

int HeatGrid::at(int iy, int iz) const {
  if (iy < 0 || iz < 0 || iy >= ny || iz >= nz) return -1;
  return index[static_cast<std::size_t>(iz) * static_cast<std::size_t>(ny) +
               static_cast<std::size_t>(iy)];
}

int HeatGrid::locate(double y, double z) const {
  require(hy > 0 && hz > 0, "heat2d: empty grid");
  int iy = static_cast<int>(std::lround((y - y0) / hy));
  int iz = static_cast<int>(std::lround((z - z0) / hz));
  int c = at(iy, iz);
  if (c >= 0) return c;
  // Fall back on the nearest active cell in a growing ring: rebars sitting a
  // hair outside a rasterised outline must still find their concrete cell.
  for (int r = 1; r <= 3; ++r) {
    int best = -1;
    double bd = 0;
    for (int dz = -r; dz <= r; ++dz)
      for (int dy = -r; dy <= r; ++dy) {
        if (std::abs(dz) != r && std::abs(dy) != r) continue;
        int k = at(iy + dy, iz + dz);
        if (k < 0) continue;
        double d = (cells[static_cast<std::size_t>(k)].y - y) *
                       (cells[static_cast<std::size_t>(k)].y - y) +
                   (cells[static_cast<std::size_t>(k)].z - z) *
                       (cells[static_cast<std::size_t>(k)].z - z);
        if (best < 0 || d < bd) {
          best = k;
          bd = d;
        }
      }
    if (best >= 0) return best;
  }
  return -1;
}

int HeatResult::snap(int i) const {
  require(!T.empty(), "heat2d: no snapshot recorded");
  if (i < 0) return static_cast<int>(T.size()) - 1;
  require(i < static_cast<int>(T.size()),
          "heat2d: snapshot index " + std::to_string(i) + " out of range");
  return i;
}

const std::vector<double>& HeatResult::field(int i) const {
  return T[static_cast<std::size_t>(snap(i))];
}

double HeatResult::T_at(double y, double z, int i) const {
  int c = grid.locate(y, z);
  require(c >= 0, "heat2d: point (" + std::to_string(y) + ", " + std::to_string(z) +
                      ") lies outside the meshed section");
  return field(i)[static_cast<std::size_t>(c)];
}

double HeatResult::mean_temp(int i) const {
  const std::vector<double>& f = field(i);
  double num = 0, den = 0;
  for (std::size_t k = 0; k < grid.cells.size(); ++k) {
    num += grid.cells[k].area * f[k];
    den += grid.cells[k].area;
  }
  require(den > 0, "heat2d: empty section");
  return num / den;
}

double HeatResult::max_temp(int i) const {
  const std::vector<double>& f = field(i);
  double m = -1e300;
  for (double v : f) m = std::max(m, v);
  return m;
}

// ---------------------------------------------------------------------------
// Mesh construction
// ---------------------------------------------------------------------------

namespace {

// A boundary face of an active cell: exposure, area, half-cell distance.
struct BoundaryFace {
  int cell = 0;
  int expo = -1;     // index into the exposure table; -1 = unexposed (adiabatic)
  double area = 0;   // face area per unit bar length [m]
  double dist = 0;   // cell centre -> face distance [m]
  double ts = 20.0;  // running surface temperature [degC]
};

struct InteriorFace {
  int a = 0, b = 0;
  double area = 0;  // [m]
  double dist = 0;  // centre-to-centre distance [m]
};

struct Mesh {
  HeatGrid grid;
  std::vector<InteriorFace> interior;
  std::vector<BoundaryFace> boundary;
};

Mesh build_mesh(const FireSection& s, const std::vector<Exposure>& expo, int max_cells) {
  s.validate();
  Mesh m;
  HeatGrid& g = m.grid;

  double ymin = 0, zmin = 0;
  int nyc = 0, nzc = 0;
  // Face exposures of the bounding box, in the order (+y, -y, +z, -z).
  int box_expo[4] = {-1, -1, -1, -1};

  if (s.kind == FireSection::Kind::Slab) {
    // 1-D through the thickness: one column of cells, lateral faces adiabatic.
    nyc = 1;
    nzc = std::max(1, static_cast<int>(std::llround(s.thickness / s.maille)));
    g.hz = s.thickness / nzc;  // snap the pitch so the thickness is spanned exactly
    g.hy = s.width;
    ymin = -s.width / 2;
    zmin = -s.thickness / 2;
    box_expo[2] = s.expo_vy;
    box_expo[3] = s.expo_wy;
  } else if (s.kind == FireSection::Kind::Rectangle) {
    nyc = std::max(1, static_cast<int>(std::llround(s.width / s.maille)));
    nzc = std::max(1, static_cast<int>(std::llround(s.thickness / s.maille)));
    g.hy = s.width / nyc;
    g.hz = s.thickness / nzc;
    ymin = -s.width / 2;
    zmin = -s.thickness / 2;
    box_expo[0] = s.expo_vz;
    box_expo[1] = s.expo_wz;
    box_expo[2] = s.expo_vy;
    box_expo[3] = s.expo_wy;
  } else {
    double y_lo = 1e300, y_hi = -1e300, z_lo = 1e300, z_hi = -1e300;
    for (const Vec3& p : s.polygon) {
      y_lo = std::min(y_lo, p.x());
      y_hi = std::max(y_hi, p.x());
      z_lo = std::min(z_lo, p.y());
      z_hi = std::max(z_hi, p.y());
    }
    require(y_hi > y_lo && z_hi > z_lo, "fire section: degenerate PERIMETRE outline");
    nyc = std::max(1, static_cast<int>(std::ceil((y_hi - y_lo) / s.maille)));
    nzc = std::max(1, static_cast<int>(std::ceil((z_hi - z_lo) / s.maille)));
    g.hy = g.hz = s.maille;
    // Centre the raster on the outline so a section symmetric about y = 0 or
    // z = 0 keeps that symmetry in the mesh (manual p221).
    ymin = 0.5 * (y_lo + y_hi) - 0.5 * nyc * s.maille;
    zmin = 0.5 * (z_lo + z_hi) - 0.5 * nzc * s.maille;
  }

  require(static_cast<long long>(nyc) * nzc <= max_cells,
          "heat2d: MAILLE " + std::to_string(s.maille) + " m would produce " +
              std::to_string(static_cast<long long>(nyc) * nzc) +
              " cells (limit " + std::to_string(max_cells) + ")");

  g.ny = nyc;
  g.nz = nzc;
  g.y0 = ymin + 0.5 * g.hy;
  g.z0 = zmin + 0.5 * g.hz;

  g.index.assign(static_cast<std::size_t>(nyc) * static_cast<std::size_t>(nzc), -1);
  for (int iz = 0; iz < nzc; ++iz)
    for (int iy = 0; iy < nyc; ++iy) {
      double cy = g.y0 + iy * g.hy;
      double cz = g.z0 + iz * g.hz;
      if (s.kind == FireSection::Kind::Perimeter && !point_in_polygon(s.polygon, cy, cz))
        continue;
      g.index[static_cast<std::size_t>(iz) * static_cast<std::size_t>(nyc) +
              static_cast<std::size_t>(iy)] = static_cast<int>(g.cells.size());
      g.cells.push_back({cy, cz, g.hy * g.hz});
    }
  require(!g.cells.empty(), "heat2d: the mesh is empty — check MAILLE against the outline");

  // Faces. Directions: (+y, -y, +z, -z).
  const int dy[4] = {1, -1, 0, 0};
  const int dz[4] = {0, 0, 1, -1};
  for (int iz = 0; iz < nzc; ++iz)
    for (int iy = 0; iy < nyc; ++iy) {
      int c = g.at(iy, iz);
      if (c < 0) continue;
      for (int d = 0; d < 4; ++d) {
        int nb = g.at(iy + dy[d], iz + dz[d]);
        // A face normal to y has height hz and a centre-to-centre distance hy;
        // a face normal to z has width hy and distance hz.
        double area = (d < 2) ? g.hz : g.hy;
        double pitch = (d < 2) ? g.hy : g.hz;
        if (nb >= 0) {
          if (nb > c) m.interior.push_back({c, nb, area, pitch});
          continue;
        }
        BoundaryFace bf;
        bf.cell = c;
        bf.area = area;
        bf.dist = 0.5 * pitch;
        if (s.kind == FireSection::Kind::Perimeter) {
          // Nearest polygon segment to the face midpoint (manual p220: the
          // exposure is carried by the segment, the mesh by the raster).
          double fy = g.cells[static_cast<std::size_t>(c)].y + 0.5 * dy[d] * g.hy;
          double fz = g.cells[static_cast<std::size_t>(c)].z + 0.5 * dz[d] * g.hz;
          std::size_t n = s.polygon.size(), best = 0;
          double bd = 1e300;
          for (std::size_t k = 0; k < n; ++k) {
            double dd = dist2_to_segment(fy, fz, s.polygon[k], s.polygon[(k + 1) % n]);
            if (dd < bd) {
              bd = dd;
              best = k;
            }
          }
          bf.expo = s.seg_exposure[best];
        } else {
          bf.expo = box_expo[d];
        }
        if (bf.expo >= 0) {
          require(bf.expo < static_cast<int>(expo.size()),
                  "heat2d: exposure index " + std::to_string(bf.expo) + " is out of range");
          if (expo[static_cast<std::size_t>(bf.expo)].adiabatic()) continue;
        } else {
          continue;  // unexposed face = adiabatic
        }
        m.boundary.push_back(bf);
      }
    }
  return m;
}

int csc_index(const SpMat& A, int r, int c) {
  for (int p = A.outerIndexPtr()[c]; p < A.outerIndexPtr()[c + 1]; ++p)
    if (A.innerIndexPtr()[p] == r) return p;
  fail("heat2d: internal — missing sparsity entry");
}

}  // namespace

// ---------------------------------------------------------------------------
// Solver
// ---------------------------------------------------------------------------

HeatResult solve_heat2d(const FireSection& section, const ThermalMaterial& mat,
                        const std::vector<Exposure>& exposures, const HeatOptions& opt) {
  require(opt.dt > 0 && finite(opt.dt), "heat2d: time step must be > 0");
  require(opt.duration >= 0 && finite(opt.duration), "heat2d: FEU DUREE must be >= 0 minutes");
  require(finite(opt.t_init), "heat2d: FEU TEMP must be finite");

  Mesh mesh = build_mesh(section, exposures, opt.max_cells);
  const HeatGrid& g = mesh.grid;
  int n = static_cast<int>(g.cells.size());

  // Output instants: the requested snapshots plus the fire duration, sorted.
  std::vector<double> outs = opt.snapshots;
  outs.push_back(opt.duration);
  std::sort(outs.begin(), outs.end());
  outs.erase(std::unique(outs.begin(), outs.end()), outs.end());
  for (double t : outs)
    require(t >= 0 && t <= opt.duration + 1e-12,
            "heat2d: snapshot at " + std::to_string(t) +
                " min is outside the fire duration " + std::to_string(opt.duration) + " min");

  HeatResult res;
  res.grid = g;

  check_table(mat.lambda_table, "conductivity");
  check_table(mat.cp_table, "specific heat");
  check_table(mat.enthalpy_table, "enthalpy");
  for (std::size_t k = 1; k < mat.enthalpy_table.size(); ++k)
    require(mat.enthalpy_table[k].second >= mat.enthalpy_table[k - 1].second, "heat2d: the enthalpy must not decrease");
  for (const Exposure& e : exposures) check_table(e.convection_table, "convection");
  std::vector<double> T(static_cast<std::size_t>(n), opt.t_init);
  if (opt.t_initial)
    for (int i = 0; i < n; ++i)
      T[static_cast<std::size_t>(i)] = opt.t_initial(g.cells[static_cast<std::size_t>(i)].y, g.cells[static_cast<std::size_t>(i)].z);
  std::vector<double> lam(static_cast<std::size_t>(n), 0.0);
  for (BoundaryFace& b : mesh.boundary) b.ts = T[static_cast<std::size_t>(b.cell)];

  // --- constant sparsity pattern (lower triangle) ---
  std::vector<Triplet> pat;
  pat.reserve(static_cast<std::size_t>(n) + mesh.interior.size());
  for (int i = 0; i < n; ++i) pat.emplace_back(i, i, 1.0);
  for (const InteriorFace& f : mesh.interior) {
    int r = std::max(f.a, f.b), c = std::min(f.a, f.b);
    pat.emplace_back(r, c, 0.0);
  }
  SpMat A(n, n);
  A.setFromTriplets(pat.begin(), pat.end());
  A.makeCompressed();
  std::vector<int> diag_ix(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) diag_ix[static_cast<std::size_t>(i)] = csc_index(A, i, i);
  std::vector<int> face_ix(mesh.interior.size());
  for (std::size_t k = 0; k < mesh.interior.size(); ++k) {
    const InteriorFace& f = mesh.interior[k];
    face_ix[k] = csc_index(A, std::max(f.a, f.b), std::min(f.a, f.b));
  }

  Eigen::SimplicialLDLT<SpMat> ldlt;
  ldlt.analyzePattern(A);

  VecX rhs(n), x(n);
  double* vals = A.valuePtr();
  int nnz = static_cast<int>(A.nonZeros());

  auto record = [&](double t_min) {
    res.times.push_back(t_min);
    res.T.push_back(T);
  };

  std::size_t next_out = 0;
  if (outs[0] <= 0) {
    record(0.0);
    ++next_out;
  }

  double t = 0.0;                        // [s]
  double t_end = opt.duration * 60.0;    // [s]
  int guard = 0;
  const int max_steps = 10000000;
  while (t < t_end - 1e-9 && next_out < outs.size()) {
    require(++guard < max_steps, "heat2d: too many time steps — increase dt");
    double step = opt.dt;
    double t_next_out = outs[next_out] * 60.0;
    if (t + step > t_next_out) step = t_next_out - t;
    double t_new_min = (t + step) / 60.0;

    // Lagged properties (previous step's temperatures).
    for (int i = 0; i < n; ++i) lam[static_cast<std::size_t>(i)] = mat.lambda(T[static_cast<std::size_t>(i)]);
    const std::vector<double> T_old = T;
    std::fill(vals, vals + nnz, 0.0);
    rhs.setZero();
    for (std::size_t k = 0; k < mesh.interior.size(); ++k) {
      const InteriorFace& f = mesh.interior[k];
      double la = lam[static_cast<std::size_t>(f.a)], lb = lam[static_cast<std::size_t>(f.b)];
      // Series of two half-cells: 1 / (d/(2 la) + d/(2 lb)).
      double u = (la + lb > 0) ? 2.0 * la * lb / (f.dist * (la + lb)) : 0.0;
      u *= f.area;
      vals[diag_ix[static_cast<std::size_t>(f.a)]] += u;
      vals[diag_ix[static_cast<std::size_t>(f.b)]] += u;
      vals[face_ix[k]] -= u;
    }
    for (BoundaryFace& b : mesh.boundary) {
      const Exposure& e = exposures[static_cast<std::size_t>(b.expo)];
      double tg = e.gas_temp(t_new_min);
      double lc = lam[static_cast<std::size_t>(b.cell)];
      double r_half = b.dist / lc;  // half-cell conduction resistance [m^2K/W]
      double u_area = 0.0, t_drive = tg;
      if (e.dirichlet()) {
        // LOI_PAROI: the wall temperature itself is prescribed; only the
        // protection and the half cell stand between it and the cell centre.
        double rt = e.protection + r_half;
        u_area = rt > 0 ? 1.0 / rt : 0.0;
        b.ts = tg;
      } else {
        // Fixed point on the wall temperature (root of the quartic obtained by
        // equating the surface flux with the flux through protection + cell).
        double ti = T[static_cast<std::size_t>(b.cell)];
        double ts = b.ts;
        double h_surf = 0;
        for (int it = 0; it < 25; ++it) {
          double tgk = tg + kAbsZero, tsk = ts + kAbsZero;
          double h_rad = mat.emissivity * kSigmaSB * (tgk * tgk + tsk * tsk) * (tgk + tsk);
          h_surf = e.alpha_c(ts) + h_rad;
          double rt = (h_surf > 0 ? 1.0 / h_surf : 0.0) + e.protection + r_half;
          double q = rt > 0 ? (tg - ti) / rt : 0.0;
          double ts_new = h_surf > 0 ? tg - q / h_surf : tg;
          double d = ts_new - ts;
          ts = ts_new;
          if (std::fabs(d) < 1e-8) break;
        }
        b.ts = ts;
        double rt = (h_surf > 0 ? 1.0 / h_surf : 0.0) + e.protection + r_half;
        u_area = rt > 0 ? 1.0 / rt : 0.0;
      }
      double u = u_area * b.area;
      vals[diag_ix[static_cast<std::size_t>(b.cell)]] += u;
      rhs(b.cell) += u * t_drive;
    }

    if (mat.enthalpy_table.empty()) {
      // lagged capacity: one linear solve
      for (int i = 0; i < n; ++i) {
        const double t0 = T_old[static_cast<std::size_t>(i)];
        const double cap = mat.rho_cp(t0) * g.cells[static_cast<std::size_t>(i)].area / step;
        vals[diag_ix[static_cast<std::size_t>(i)]] += cap;
        rhs(i) += cap * t0;
      }
      ldlt.factorize(A);
      require(ldlt.info() == Eigen::Success, "heat2d: conduction matrix factorization failed");
      x = ldlt.solve(rhs);
      require(x.allFinite(), "heat2d: non-finite temperature field");
      for (int i = 0; i < n; ++i) T[static_cast<std::size_t>(i)] = x(i);
    } else {
      // enthalpy: R(T) = a/dt (H(T) - H(T_old)) + K T - f = 0 is the gradient
      // of a strictly convex functional (H non-decreasing, K SPD): Newton
      // with the tangent capacity, safeguarded by an exact line search on
      // the monotone directional derivative — convergent whatever the
      // latent heat or the step, and H(T_new) - H(T_old) balances the
      // conducted heat exactly at convergence
      const SpMat K = A;  // conduction + boundary, without capacity
      const std::vector<double> base(vals, vals + nnz);
      const VecX f = rhs;
      VecX ca(n), Hold(n), Tk(n);
      for (int i = 0; i < n; ++i) {
        ca(i) = g.cells[static_cast<std::size_t>(i)].area / step;
        Hold(i) = mat.enthalpy(T_old[static_cast<std::size_t>(i)]);
        Tk(i) = T_old[static_cast<std::size_t>(i)];
      }
      auto resid = [&](const VecX& tv) {
        VecX r = K.selfadjointView<Eigen::Lower>() * tv - f;
        for (int i = 0; i < n; ++i) r(i) += ca(i) * (mat.enthalpy(tv(i)) - Hold(i));
        return r;
      };
      VecX R = resid(Tk);
      for (int iter = 0;; ++iter) {
        require(iter < 200, "heat2d: the enthalpy iteration did not converge");
        std::copy(base.begin(), base.end(), vals);
        for (int i = 0; i < n; ++i) vals[diag_ix[static_cast<std::size_t>(i)]] += ca(i) * mat.step_capacity(Tk(i), Tk(i));
        ldlt.factorize(A);
        require(ldlt.info() == Eigen::Success, "heat2d: conduction matrix factorization failed");
        const VecX d = -ldlt.solve(R);
        require(d.allFinite(), "heat2d: non-finite temperature field");
        const double g0 = d.dot(R);  // < 0: a descent direction
        double alpha = 1.0;
        VecX Rn = resid(Tk + d);
        double g1 = d.dot(Rn);
        if (g0 < 0 && g1 > 0) {  // overshoot: the root of phi'(alpha) on [0, 1] (Illinois)
          double lo = 0, glo = g0, hi = 1, ghi = g1;
          for (int ls = 0; ls < 60; ++ls) {
            alpha = (lo * ghi - hi * glo) / (ghi - glo);
            Rn = resid(Tk + alpha * d);
            const double ga = d.dot(Rn);
            if (std::fabs(ga) <= 1e-3 * std::fabs(g0)) break;
            if (ga < 0) {
              lo = alpha;
              glo = ga;
              ghi *= 0.5;
            } else {
              hi = alpha;
              ghi = ga;
              glo *= 0.5;
            }
          }
        }
        Tk += alpha * d;
        R = Rn;
        const double change = alpha * d.cwiseAbs().maxCoeff(), scale = std::max(1.0, Tk.cwiseAbs().maxCoeff());
        if (change <= 1e-10 * scale) break;
      }
      for (int i = 0; i < n; ++i) T[static_cast<std::size_t>(i)] = Tk(i);
    }

    t += step;
    if (t >= t_next_out - 1e-9) {
      record(outs[next_out]);
      ++next_out;
    }
  }
  // Zero-duration (or all snapshots already recorded) runs still report a field.
  while (next_out < outs.size()) {
    record(outs[next_out]);
    ++next_out;
  }
  return res;
}

}  // namespace keta

// Cross-sections: parametric shapes with closed-form properties (A, Iy, Iz, J,
// shear areas, elastic/plastic moduli) plus generic polygons (with holes) whose
// in-plane properties come from exact Green-theorem integrals. Stress recovery
// points are generated per shape for sigma = N/A + My*z/Iy - Mz*y/Iz.
//
// Section local axes: y and z are the principal axes used by the frame element;
// the section lies in the local y-z plane. For polygon input the coordinates are
// re-centered to the centroid; principal rotation is NOT applied automatically
// (a warning-grade error is raised when products of inertia are significant).
#pragma once
#include <limits>
#include <memory>
#include <variant>

#include "keta/core.h"

namespace keta {

// Result of the exact Tri6 warping/torsion engine (`keta/section_exact.h`).
// Declared here only so a `Section` can carry one; this header stays free of the
// finite-element machinery.
struct ExactSection;

struct SectionProps {
  double A = 0;     // area [m^2]
  double Avy = 0;   // shear area along local y (0 = rigid in shear -> Euler-Bernoulli)
  double Avz = 0;   // shear area along local z
  double Iy = 0;    // second moment about local y [m^4] (bending in x-z plane)
  double Iz = 0;    // second moment about local z (bending in x-y plane)
  double J = 0;     // St-Venant torsion constant [m^4]
  double Wely = 0;  // elastic modulus |Iy / z_max|
  double Welz = 0;
  double Wply = 0;  // plastic modulus (0 = unknown)
  double Wplz = 0;
  double Ip = 0;    // polar moment for the torsional mass inertia (0 = Iy + Iz)
  // --- appended ---
  // Warping constant C_M [m^6] about the shear centre (0 = Saint-Venant
  // torsion only), read by the Vlasov drivers (keta/warping.h) when their
  // spec gives none for the section; the I shape computes tf b^3 (h - tf)^2 / 24.
  double Iw = 0;
  // Wagner monosymmetry coefficient beta_x [m] (0 for a doubly symmetric
  // section). Carried for the section data; the Vlasov drivers cover the
  // doubly symmetric case and refuse a non-zero value.
  double beta_x = 0;
};

// A named point of the section outline where stresses are recovered.
struct StressPoint {
  std::string name;
  double y = 0, z = 0;  // section coordinates [m]
  double tau_factor_y = 0;  // fraction of Vy/Avy contributing to shear stress here
  double tau_factor_z = 0;  // fraction of Vz/Avz contributing here
};

// ---- parametric shapes (dimensions in m) ----
struct ShapeRect { double b = 0, h = 0; };                       // b along y, h along z
struct ShapeTube { double d = 0, t = 0; };                       // outer diameter, wall
struct ShapeCircle { double d = 0; };
// Doubly-symmetric I/H. `r` = root fillet radius of a hot-rolled profile (0 =
// welded / no fillet, the historical closed forms, bit-identical). With r > 0
// the four web-flange fillets enter A, Iy, Iz, Wpl and Wel exactly (circular
// spandrels), J follows El Darwish & Johnston (1965) and Avz the EN 1993-1-1
// §6.2.6(3)a rolled-section rule; Avy and the stress points stay the
// fillet-free ones. The fiber mesh and a hand-fed exact engine ignore r; the
// section-point stress rings (section_stress.h) keep the fillets.
struct ShapeI { double h = 0, b = 0, tw = 0, tf = 0, r = 0; };
struct ShapeBox { double h = 0, b = 0, tw = 0, tf = 0; };        // rectangular hollow
struct ShapeT { double h = 0, b = 0, tw = 0, tf = 0; };          // flange on top (+z)
struct ShapePolygon {
  std::vector<Vec3> outer;               // (y, z, 0) vertices, CCW
  std::vector<std::vector<Vec3>> holes;  // CW or CCW, orientation normalized internally
};
struct ShapeGeneric {};                   // properties given directly, no geometry

using Shape =
    std::variant<ShapeGeneric, ShapeRect, ShapeCircle, ShapeTube, ShapeI, ShapeBox, ShapeT,
                 ShapePolygon>;

struct Section {
  Id id = 0;
  std::string name;
  Shape shape = ShapeGeneric{};
  SectionProps props;                    // computed by finalize() for parametric shapes
  std::vector<StressPoint> stress_points;

  // Optional exact backend: null unless `attach_exact()` (section_exact.h) ran.
  // When present it carries the warping-FE J / Iw / shear centre / shear areas
  // and the stress evaluators; `props.J`, `props.Avy` and `props.Avz` are then
  // the exact ones. Shared, so every bar using this section holds one solve.
  // NOT serialised (io.cpp): it is a derived quantity, re-attached on demand.
  std::shared_ptr<const ExactSection> exact;

  // --- appended ---
  // Shear areas of a parametric shape: the energy-consistent ones (Saint-
  // Venant flexure) for the Poisson ratio of the material of its bars, set
  // when a bar is added (Model::add_bar; a section shared by two ratios is
  // duplicated). `shear_nu`: the ratio they were computed for (NaN: not yet);
  // `avy_shape` / `avz_shape`: the closed-form placeholders of finalize(), so
  // shear areas set by hand (0 = Euler-Bernoulli) are kept.
  double shear_nu = std::numeric_limits<double>::quiet_NaN();
  double avy_shape = -1, avz_shape = -1;
  Id shear_base = 0;  // the section a ratio-specific copy was made from

  // Compute props + stress points from the shape (no-op for ShapeGeneric where
  // the caller supplied props directly). Throws on degenerate dimensions.
  void finalize();
};

// Shear areas of the parametric shapes: Energy (default) = the energy-
// consistent ones of Saint-Venant flexure for the material's Poisson ratio
// (closed forms for the rectangle, circle and tube, the warping solve
// otherwise; Model::resolve_shear); ShapeFactor = the classical factors 5/6
// (rectangle), 0.9 (circle).
enum class ShearAreas { Energy, ShapeFactor };
void set_default_shear_areas(ShearAreas s);
ShearAreas default_shear_areas();

// Closed-form property computation per shape (also used by the fiber mesher).
SectionProps props_rect(const ShapeRect&);
SectionProps props_circle(const ShapeCircle&);
SectionProps props_tube(const ShapeTube&);
SectionProps props_i(const ShapeI&);
SectionProps props_box(const ShapeBox&);
SectionProps props_t(const ShapeT&);
// Polygon: exact area/centroid/inertia by Green's theorem; J and shear areas by
// thin-wall/compact heuristics (documented approximation), Wel from bounding box.
SectionProps props_polygon(const ShapePolygon&);

}  // namespace keta

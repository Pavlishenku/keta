// Outlines of hot-rolled profiles beyond the I/H catalog (catalog.h):
// angles (EN 10056-1), taper-flange I beams INP (DIN 1025-1) and parallel-
// flange channels UPE (DIN 1026-2 / EN 10279), built as polygons with their
// root and toe fillets as circular arcs tangent to the adjacent faces
// (polygonised with `segments` chords per fillet: the area deficit of the
// chords is r^2 (theta - sin theta)/2 per chord, < 1e-8 relative at 256),
// and the full polygon inertia — A, centroid, Iy, Iz AND the product Iyz,
// with the principal angle — by Green's theorem (exact on the polygon).
// Units: metres. Axes: y horizontal, z vertical, origin at the outline's
// lower-left reference (angle heel / I and U web axis at mid-height... see
// each builder).
#pragma once
#include <optional>
#include <string>
#include <vector>

#include "keta/core.h"

namespace keta {

// A simple polygon (CCW corners) whose corner k is rounded by radius r[k]
// (0 = sharp): the arc tangent to both adjacent edges, convex or concave.
std::vector<Vec3> filleted_polygon(const std::vector<Vec3>& corners, const std::vector<double>& radii,
                                   int segments = 256);

struct PolygonInertia {
  double A = 0, yc = 0, zc = 0;
  double Iy = 0, Iz = 0, Iyz = 0;  // centroidal, Iy = int (z - zc)^2 dA, Iyz = int (y - yc)(z - zc) dA
  double I1 = 0, I2 = 0, alpha = 0;  // principal values and the angle of axis 1 from y [rad]
};
PolygonInertia polygon_inertia(const std::vector<Vec3>& ring);

// Equal-leg angle L h x t (heel at the origin, legs along +y and +z), root
// radius r1, toe radius r2.
std::vector<Vec3> outline_angle(double h, double t, double r1, double r2, int segments = 256);
// Taper-flange I beam INP: depth h, width b, web tw, flange thickness tf at
// b/4 from the toe, inner flange slope `slope` (0.14 for DIN 1025-1 up to
// INP 300), root radius r1 (= tw), toe radius r2; symmetric about y = 0 and
// z = h/2 (origin at the bottom flange's mid-width).
std::vector<Vec3> outline_inp(double h, double b, double tw, double tf, double r1, double r2,
                              double slope = 0.14, int segments = 256);
// Parallel-flange channel UPE: depth h, width b, web tw, flanges tf, root
// radius r (web back at y = 0).
std::vector<Vec3> outline_upe(double h, double b, double tw, double tf, double r, int segments = 256);

// Tabulated outlines by designation ("L100x10", "I100", "UPE100"; tolerant
// of case and spaces); nullopt when unknown. Only the sizes whose dimensions
// are in the table below.
std::optional<std::vector<Vec3>> profile_outline(const std::string& designation);

}  // namespace keta

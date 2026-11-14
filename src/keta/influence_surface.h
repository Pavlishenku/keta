// Influence SURFACES on a node grid — the 2-D analogue of `influence.h`.
//
// ---------------------------------------------------------------------------
// Physics
// ---------------------------------------------------------------------------
// Muller-Breslau: eta(x_i, y_j) is the value of ONE chosen response (a
// support reaction component, a nodal displacement component, or a member
// force component at a station) produced by a UNIT load standing on grid
// node (i, j). It is computed by the adjoint sweep of the influence LINE
// (unit_load_sweep, keta/influence.h): one factorization of K, one solve
// K w = g for the response, then O(1) per load position — the same numbers
// as a direct multi-RHS sweep over the nx*ny unit loads, up to round-off. A
// surface is a line whose load positions happen to be laid out on a grid.
//
// The response is LINEAR in the load, so a load pattern p(x, y) gives
//     E = integral over the deck of p(x,y) eta(x,y) dx dy
// for a distributed pattern, and E = sum_k P_k eta(x_k, y_k) for wheels. This
// is exactly the 2-D version of the convoy / UDL algebra in `influence.h`.
//
// ---------------------------------------------------------------------------
// Interpolation and the two integrals (documented numerical choices)
// ---------------------------------------------------------------------------
// Between the nodes eta is BILINEAR on each grid cell (`at`). Two integrals are
// offered and they use the same interpolant in two different ways:
//
//  * `integral_box` — EXACT integral of the bilinear field over an arbitrary
//    clipped rectangle, from the separable 1-D shape-function integrals
//    (int N_a dx)(int N_b dy). This is what a rolling PATCH load uses: the
//    patch is applied whole, favourable or not, like a truck contact area.
//
//  * `integral_pos` / `integral_neg` — the FAVOURABLE-REGION integral of a free
//    UDL that may load any number of disjoint zones. The zero set of a bilinear
//    function is a hyperbola, which has no tidy closed-form clip, so each cell
//    is split into two triangles along BOTH diagonals and the two splittings
//    are averaged. On a triangle the interpolant is linear, its zero set is a
//    straight line, and the sign split is EXACT (Sutherland-Hodgman clip +
//    fan integration). Averaging the two diagonals reproduces the bilinear cell
//    integral exactly on a sign-uniform cell:
//        (2a+b+c+2d)/6 averaged with (a+2b+2c+d)/6  =  (a+b+c+d)/4.
//    Cells straddling the zero contour carry an O(h^2) split error, named here
//    rather than hidden.
//
// ---------------------------------------------------------------------------
// Units, conventions, limitations
// ---------------------------------------------------------------------------
// SI: m, N, Pa. `dir` is the direction of the unit load (e.g. (0,0,-1) for a
// downward wheel); eta is then "response per newton". The grid must be
// rectangular in PLAN (global X and Y), with strictly increasing gridlines; the
// node elevations may vary (a cambered deck grillage is fine) since the surface
// is parameterised by plan coordinates. Cable bars are excluded from the
// reaction accumulation exactly as in the influence line — an influence surface
// is a LINEAR quantity and a cable structure is not.
#pragma once
#include "keta/influence.h"

namespace keta {

// A rectangular arrangement of existing model nodes, row-major:
// index(i, j) = j * nx + i, i along the grid's X axis, j along Y.
struct NodeGrid {
  int nx = 0, ny = 0;
  std::vector<Id> nodes;
  std::vector<double> x, y;  // gridlines [m], strictly increasing
  Id at(int i, int j) const;
};

// Validate that `nodes` (row-major, nx*ny of them) really form a plan-
// rectangular grid and extract the gridlines. `tol` is relative to the grid
// extent. Throws a typed `Error` on a ragged / non-monotone / duplicated grid.
NodeGrid build_node_grid(const Model& model, const std::vector<Id>& nodes, int nx, int ny,
                         double tol = 1e-9);

struct InfluenceSurface {
  int nx = 0, ny = 0;
  std::vector<double> x, y;
  std::vector<double> eta;  // row-major, response per unit load

  double node_eta(int i, int j) const;
  double at(double px, double py) const;  // bilinear; 0 outside the grid
  double integral_box(double x0, double x1, double y0, double y1) const;  // exact, clipped
  double integral() const;                                                // whole grid
  double integral_pos() const;                                            // eta > 0 zones
  double integral_neg() const;                                            // eta < 0 zones
  std::pair<double, double> point_extreme() const;  // (max, min) of a single unit load
};

// One factorization, nx*ny unit loads of direction `dir`, one response.
InfluenceSurface influence_surface(const Model& model, const NodeGrid& grid, const Vec3& dir,
                                   const ResponseKind& response);
// Same sweep against a caller-owned factor cache (see influence_line), with
// optional cancel / progress hooks (see unit_load_sweep).
InfluenceSurface influence_surface(const Model& model, const NodeGrid& grid, const Vec3& dir,
                                   const ResponseKind& response, FactorCache& cache,
                                   const SolveHooks* hooks = nullptr);

// ---------------- surface algebra (the 2-D convoy) ----------------

struct SurfaceExtremes {
  double max = 0, min = 0;
  double rx_max = 0, ry_max = 0;  // reference position producing the maximum
  double rx_min = 0, ry_min = 0;
};

// Rectangular patch a (along x) by b (along y) carrying the pressure `p`
// [N/m^2], stepped over the grid; the reference position is the patch CENTRE
// and the part of the patch hanging off the grid contributes zero. The patch is
// applied WHOLE (favourable and unfavourable parts alike).
SurfaceExtremes patch_extreme(const InfluenceSurface& s, double a, double b, double p,
                              int steps_x = 60, int steps_y = 60);

// Free UDL of pressure `p` that may load ANY number of zones: the favourable-
// region integral. max = p * integral_pos, min = p * integral_neg.
SurfaceExtremes udl_surface_extreme(const InfluenceSurface& s, double p);

// A wheel set: offsets from the moving reference point and their loads.
struct WheelSet {
  std::vector<std::array<double, 3>> wheels;  // (dx, dy, P>0 factor)
  WheelSet mirrored_x() const;                // travel in the other direction
};

// Slide the wheel set over the grid; wheels off the deck contribute zero.
SurfaceExtremes wheel_set_extreme(const InfluenceSurface& s, const WheelSet& ws, int steps_x = 60,
                                  int steps_y = 60);

// Value of a wheel set at one explicit reference position (the sweep kernel,
// exposed so a hand placement can be checked against the swept extreme).
double wheel_set_value(const InfluenceSurface& s, const WheelSet& ws, double rx, double ry);

}  // namespace keta

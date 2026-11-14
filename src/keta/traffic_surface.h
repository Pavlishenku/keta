// EN 1991-2 Load Model 1 placed on an influence SURFACE (keta/influence_
// surface.h): the notional lanes, their numbering, the tandem systems and
// the UDL arranged for the most adverse effect, and the concomitant effect
// of that arrangement on other surfaces.
//
// --- Rules (EN 1991-2 §4.2.4, §4.3.2) ---
// Carriageway [y0, y1] across, [x0, x1] along (plan, the surface's axes).
// n notional lanes of width w (EN division, §4.2.3), anywhere across the
// carriageway, not overlapping; lane of rank r = 1..3 carries a tandem
// system of two axles alpha_Qr Q_rk 1.2 m apart, each of two wheels 2 m
// apart, centred in the lane (general effects), and the UDL alpha_qr q_rk on
// the UNFAVOURABLE parts of the lane only, longitudinally and transversally;
// the remaining area carries alpha_q q_rk on its unfavourable parts. Ranks
// above 3 carry the remaining-area UDL and no tandem. A tandem that is
// favourable wherever it stands in its lane is left off (loads act only
// where unfavourable).
//
// --- Exact ingredients ---
// The surface is bilinear per grid cell. For a lane centred at yc, the
// tandem's effect is piecewise linear in its position x (each wheel inside
// one cell), so its extreme is attained where a wheel stands on a gridline
// x: those positions are all evaluated. The favourable integral over a
// rectangle (a lane strip, the carriageway) uses the surface's own
// favourable-region rule (each cell split along both diagonals, linear
// triangles clipped by the rectangle and by the zero line, averaged —
// keta/influence_surface.h), exactly. Lane positions are searched on a
// transverse step (step_y) together with every position where a wheel or a
// lane edge meets a gridline y (the breakpoints of the effect); the best
// non-overlapping lanes and their ranks are found exhaustively.
//
// --- Concomitant effects ---
// lm1_surface_value: the effect on `target` of the arrangement found for
// `driver` — same lanes, same tandem positions, UDL on the zones where the
// driver is unfavourable.
#pragma once
#include "keta/influence_surface.h"
#include "keta/traffic.h"

namespace keta {

// Integral of `value` over the rectangle [x0, x1] x [y0, y1] restricted to
// where sign * driver > 0 (sign 0: no restriction), both surfaces on one grid,
// by the two-diagonal linear triangles (exact for that rule).
double surface_integral_where(const InfluenceSurface& value, const InfluenceSurface& driver, int sign, double x0,
                              double x1, double y0, double y1);

struct Lm1SurfaceOptions {
  double y0 = 0, y1 = 0;  // carriageway across
  double x0 = -1e300, x1 = 1e300;  // along (clipped to the surface)
  double lane_width = 3.0;
  int n_lanes = 0;  // 0: EN 1991-2 §4.2.3 division of the carriageway width
  Lm1Factors factors;
  double step_y = 0.05;
  bool minimise = false;  // the most negative effect instead of the most positive
  // Appended: boxes (x0, x1, y0, y1) where no LM1 load stands (the zone
  // around a special vehicle, EN 1991-2 Annex A): a tandem may not overlap
  // one, the UDLs (lanes and remaining area) leave them out — exactly (the
  // favourable integrals of the boxes are subtracted, the tandem searched on
  // the admissible intervals and their ends).
  std::vector<std::array<double, 4>> exclude;
};

struct Lm1Lane {
  int rank = 0;
  double y_start = 0;  // lane edge (lane = [y_start, y_start + w])
  bool tandem = true;  // false: the tandem is favourable everywhere, left off
  double ts_x = 0;     // tandem's first axle
  double ts = 0, udl = 0;  // the lane's two contributions
};

struct Lm1SurfaceResult {
  double value = 0;
  std::vector<Lm1Lane> lanes;  // ranks 1..min(n, 3)
  double remaining = 0;        // the remaining area's contribution (and lanes of rank > 3)
};

Lm1SurfaceResult lm1_surface_extreme(const InfluenceSurface& s, const Lm1SurfaceOptions& o);
double lm1_surface_value(const InfluenceSurface& target, const InfluenceSurface& driver, const Lm1SurfaceResult& r,
                         const Lm1SurfaceOptions& o);

}  // namespace keta

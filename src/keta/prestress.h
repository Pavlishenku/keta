// Prestress tendons by the equivalent-load method.
//  * Profile e_y(s), e_z(s): piecewise cubic over the bar chain (section
//    frame, centroid origin — e doubles as the moment arm).
//  * Force profile P(s): jacking force with EN 1992-1-1 §5.10.5.2 friction
//    P(s) = P0 exp(-mu (theta(s) + k s)) (theta = cumulated angular deviation,
//    the angle swept by the 3-D tangent), anchorage draw-in by the mirror rule,
//    one active anchor at either end or double-end jacking (max of both
//    curves, the two draw-ins released one after the other).
//  * Lowering per bar: anchor thrust +-P, tangent kicks +-(P e)', anchor
//    moments (My = +P e_z, Mz = -P e_y), distributed deviation w = (P e)''
//    and axial drag P' — a self-equilibrated load system reproducing the
//    isostatic identity N = -P, My = -P e_z, Mz = +P e_y exactly; on a
//    hyperstatic system the excess is the EN 1992-1-1 §5.10.7 secondary
//    (parasitic) moment, linear between supports.
#pragma once
#include "keta/statics.h"

namespace keta {

// Evaluate the tendon profile at chain abscissa s.
double tendon_ey(const Tendon& t, double s);
double tendon_ez(const Tendon& t, double s);
double tendon_dey(const Tendon& t, double s);
double tendon_dez(const Tendon& t, double s);
double tendon_d2ey(const Tendon& t, double s);
double tendon_d2ez(const Tendon& t, double s);

// Piecewise-linear force profile after friction (+ optional draw-in).
struct ForceProfile {
  std::vector<double> s;  // ascending abscissae
  std::vector<double> P;
  double at(double x) const;
  double slope(double x) const;
};
ForceProfile tendon_force_profile(const Model& model, const Tendon& t, int samples_per_meter = 20);

// Equivalent loads of the tendon at `fraction` of its jack force. Loads are
// expressed in the LOCAL frame of each chain bar.
std::vector<Load> tendon_loads(const Model& model, const Tendon& t, double fraction = 1.0,
                               int n_sub = 8);
// The same, knowing the force is nil outside the chain range [lo, hi]: the
// bars outside it get the same loads with zero intensities, unevaluated (the
// loads of a unit loss at one node of a tendon's grid).
std::vector<Load> tendon_loads(const Model& model, const Tendon& t, double fraction, int n_sub,
                               double lo, double hi);

// Convenience: chain length and per-bar spans.
double tendon_chain_length(const Model& model, const Tendon& t);

}  // namespace keta

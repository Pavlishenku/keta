// Railway fatigue traffic — EN 1991-2 Annex D.
//
//  * The 12 service trains of figures D.1 to D.12: axle positions from the
//    dimension chains of the figures, axle loads as the figures give them in
//    kN (225 kN for a 22.5 t axle: the Annex counts 10 kN per tonne, so
//    sum(Q) / 10 kN is the train mass of the mix tables exactly).
//  * The three traffic mixes, Tables D.1 to D.3 (trains per day), each about
//    25e6 t per year and track.
//  * The fatigue dynamic enhancement of a real train (§6.4.6.6 / Annex D):
//    1 + 0.5 (phi' + 0.5 phi''), phi' and phi'' from Annex C — K = v / (2 L
//    n0), phi' = K / (1 - K + K^4) (1.325 for K >= 0.76), phi'' = alpha/100
//    [56 exp(-(L/10)^2) + 50 (L n0 / 80 - 1) exp(-(L/20)^2)] >= 0, alpha =
//    min(v / 22, 1). n0 is the bridge's first frequency (the lower limit of
//    EN 1991-2 Figure 6.10 by `rail_n0_lower`).
#pragma once
#include <vector>

#include "keta/influence.h"

namespace keta {

struct ServiceTrain {
  int number = 0;           // 1..12
  double mass_t = 0;        // sum(Q) of the figure / 10 kN per tonne [t]
  double v_max_kmh = 0;     // maximum speed V of the figure [km/h]
  ConvoyPattern pattern;    // axles (offset from the first axle [m], load [N])
  double total_load() const;
  double length() const;    // first to last axle [m]
};
ServiceTrain annex_d_train(int number);

enum class RailFatigueMix { Standard, Heavy250, Light };
struct RailMixEntry {
  int train = 0;
  int per_day = 0;
  double tonnage_mt = 0;  // per_day x mass x 365 / 1e6 [1e6 t / year]
};
std::vector<RailMixEntry> annex_d_mix(RailFatigueMix mix);
double annex_d_mix_tonnage(RailFatigueMix mix);  // [1e6 t / year]

// Lower limit of the first natural frequency, EN 1991-2 Figure 6.10:
// 80 / L (4 <= L <= 20 m), 23.58 L^-0.592 (20 < L <= 100 m) [Hz].
double rail_n0_lower(double L_phi);
// phi', phi'' of a real train (Annex C) and the fatigue factor 1 + 0.5 (phi' + 0.5 phi'').
std::pair<double, double> rail_real_train_increments(double v_ms, double L_phi, double n0);
double rail_fatigue_dynamic_factor(double v_ms, double L_phi, double n0);

}  // namespace keta

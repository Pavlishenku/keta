// Résal effect: the share of the shear force carried by the inclined chords
// of a variable-depth member (EN 1992-1-1 §6.2.1(1), V_ccd and V_td; BPEL 91
// art. 7.1.2 "effet Résal").
//
// A chord (flange) whose centroid line runs along (1, dy/dx, dz/dx) carries
// its axial force N_f = sigma_f A_f along that line: on the face of outward
// normal +x the force N_f (1, dy/dx, dz/dx). The shear force of the section
// (keta/statics.h Station: resultant of the stresses on that face, dMy/dx =
// Vz) is therefore the web shear plus the transverse components of the chords:
//
//   Vz = Vz_web + sum_f N_f dz_f/dx ,   Vy = Vy_web + sum_f N_f dy_f/dx .
//
// The chord stress is the linear stress of the section at the chord centroid,
// sigma_f = N / A + My z_f / Iy - Mz y_f / Iz (principal centroidal axes), so
// the relief is linear in the forces of the section, with factors that depend
// on the section alone:
//
//   Vz - Vz_web = kz . (N, My, Mz) ,  kz = sum_f A_f dz_f/dx (1 / A, z_f / Iy, -y_f / Iz) .
//
// Two equal chords h apart with no web, the lower one inclined at alpha under
// a moment M alone: N_f = -+ M / h and the web keeps V - (M / h) tan alpha.
// The module returns raw values; whether the relief helps (a haunch deepening
// toward the larger moment) or loads the web follows from the signs.
#pragma once
#include <vector>

#include "keta/core.h"

namespace keta {

struct ResalChord {
  double area = 0;             // [m^2]
  double y = 0, z = 0;         // centroid, from the centroid of the section [m]
  double dy_dx = 0, dz_dx = 0; // slopes of the centroid line along the bar
};

struct ResalSection {
  double A = 0, Iy = 0, Iz = 0;  // whole section, principal centroidal axes
  std::vector<ResalChord> chords;
};

// Shear carried by the chords per unit of (N, My, Mz): relief_y = y . f,
// relief_z = z . f with f = (N, My, Mz). Typed errors: non-positive A, Iy, Iz,
// a chord of non-positive area or non-finite data.
struct ResalFactors {
  Vec3 y = Vec3::Zero(), z = Vec3::Zero();
};
ResalFactors resal_factors(const ResalSection& s);

struct WebShear {
  double Vy = 0, Vz = 0;              // shear left to the web
  double relief_y = 0, relief_z = 0;  // shear carried by the chords (V - V_web)
};
WebShear resal_web_shear(const ResalFactors& k, double N, double Vy, double Vz, double My, double Mz);

}  // namespace keta

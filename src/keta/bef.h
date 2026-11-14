// Exact beam on a zoned elastic (Winkler) foundation — ST1 "BARRE SUR SOL
// ELASTIQUE" (note RDM ST1 §6, after PSH2 / A. L. Millan, Annales ITBTP n°408).
//
// ---------------------------------------------------------------------------
// Physics
// ---------------------------------------------------------------------------
// Euler-Bernoulli member (shear deformation neglected, as the note and the
// manual p.60 require: no SY/SZ on CARA PSE) resting on distributed springs
// of constant modulus per zone:
//
//   bending (each plane)   EI v'''' + k v = q(x)          (moduli kfy, kfz)
//   axial                  EA u''   - k u = -q_x(x) (+ EA eps0')   (kfx)
//   torsion                GJ th''  - k th = 0                       (kmx)
//
// The note solves them with the transfer functions
//   f_p(x) = sum_n (-k/EI)^n x^(p+4n) / (p+4n)!      (bending, f_p' = f_(p-1))
//   g_p(x) = sum_n ( k/EA)^n x^(p+2n) / (p+2n)!      (axial / torsion)
// whose short series "fait disparaître toute discontinuité pour K = 0"; f_4,
// f_5 (g_2, g_3) are the Duhamel particular solutions of the linear loads.
// keta uses exactly these functions.
//
// Conditioning. Transfer matrices over a long, stiffly bedded bar grow like
// e^(beta L) — the manual itself warns (p.60) that "si l'on définit la poutre
// par l'intermédiaire d'une seule barre, ST1 peut fournir des valeurs
// erronées". keta cuts every zone into pieces with r l^4 <= 1 (bending) and
// r l^2 <= 1 (axial / torsion), builds each piece's EXACT stiffness and
// fixed-end forces, and condenses the interior piece nodes statically: the
// result is still the exact solution of the ODE (no discretization error),
// and stays well conditioned for any beta L.
//
// The member block replaces keta's Hermite consistent bed for the ST1
// `CARA PSE` zones; it is installed on `Model::tapered` (custom member block:
// stiffness, equivalent loads from the span record, station recovery), so the
// statics, the influence lines (SURCH) and every driver built on the stiffness
// see it. Unlike the Hermite bed, the station INTERNAL FORCES include the soil
// reaction (the transfer solution is evaluated at the station), as in ST1.
//
// Inputs: local axes, SI units; zone moduli [N/m^2] (force per unit length per
// unit displacement), kmx [N m / m / rad]. Loads come as keta span records
// (local distributed segments, point forces / moments, imposed strain and
// curvature, uniform or linear). A concentrated dislocation (DEFOR CON) on a
// bedded bar is a typed rejection.
#pragma once
#include "keta/statics.h"

namespace keta {

struct BefZone {
  double x1 = 0;  // end abscissa of the zone [m] (the last zone ends at L)
  double kx = 0, ky = 0, kz = 0, kt = 0;
};

struct BefBar {
  double L = 0;
  double EA = 0, EIy = 0, EIz = 0, GJ = 0;
  std::vector<BefZone> zones;  // contiguous from x = 0, ascending x1
};

// Exact local 12x12 member block (keta DOF order u v w thx thy thz, i then j).
Mat12 bef_ke(const BefBar& bar);
// Equivalent nodal loads (local 12) of a span-load record: minus the
// fixed-end forces of the exact clamped solution.
Vec12 bef_f0(const BefBar& bar, const BarSpanLoads& sp);

// Installs the exact block on `model.tapered[bar_id]` for a bar of the model
// (no release, no Bar::soil allowed: typed rejections; an eccentric bar keeps
// its rigid arms, the bed lies along the offset axis).
void install_bef(Model& model, Id bar_id, const std::vector<BefZone>& zones);

// A bar whose uniform Winkler bed (Bar::soil) the exact block can carry: a
// frame bar without release, CARA VAR kernel or shear deformation (the
// exact block is Euler-Bernoulli; a Timoshenko bar keeps the Hermite bed
// rather than lose its shear deformation).
bool exact_bed_applies(const Model& model, const Bar& bar);
// Replaces, in `model`, the Hermite bed of every bar for which
// exact_bed_applies by the exact member block (one zone over the bar).
// Returns the number of bars converted. Linear statics run on such a copy:
// exact displacements, and internal forces at the stations that include
// the soil reaction (the Hermite bed has neither).
int install_exact_beds(Model& model);

}  // namespace keta

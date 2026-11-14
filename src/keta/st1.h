// SETRA/CEREMA ST1 deck importer (pragmatic subset of the v24 language):
//   OPTION PLANE|GRILL|SPATIALE, TITRE, NOEUD / GENER NOEUD, BARRE / GENER
//   BARRE (id n1 n2 or id DE n1 A n2), APPUI (rigid DX..RZ, DECOL POS/NEG,
//   EL DI diagonal springs), EXC (rigid end offsets OR/EX, global frame),
//   ART (end releases OR/EX RY/RZ/RX), BETA (local roll), CARA <ids|TOUT>
//   SX/SY/SZ/IX/IY/IZ/VY/WY/VZ/WZ, CONS <ids|TOUT> E NU RO TEMP, CHARG blocks
//   (POIDS PROPRE, NOEUD forces, BARRE CON/UNI/LIN with XL + REL/ABS +
//   GLO/LOC/PRO, APPUI settlements DX..RZ, TEMP UNI/GY/GZ, DEFOR UNI/CON/LIN),
//   COMB result combinations, ETUDE PAS, EXEC CHARG/MODE/FLAMBEMENT,
//   integer lists with "i a j [pas p]" ranges, # comments, ; separators.
// Unknown statements are collected as warnings, never silently dropped.
// ST1 units are user-consistent; `gravity` scales POIDS PROPRE (use 1.0 for
// tonne-force decks whose RO is already a force density).
#pragma once
#include "keta/model.h"

namespace keta {

struct St1Import {
  Model model;
  std::string title;
  std::vector<std::string> warnings;
  bool exec_requested = false;   // EXEC CHARG
  bool exec_modal = false;       // EXEC MODE
  bool exec_buckling = false;    // EXEC FLAMBEMENT
  double etude_pas = 0;          // ETUDE PAS (study-section step hint; 0 = unset)
};

struct St1Options {
  // Units multiplier of POIDS PROPRE: ST1's RO is a WEIGHT density, so a
  // self-weight weighs gravity x RO x A per unit length (1 = ST1 semantics;
  // e.g. 9.80665 for a deck whose RO is a mass density in kg/m3 and whose
  // forces are in N). The mass density is RO / G_DYN regardless.
  double gravity = 1.0;
  // RESU edition: false = keta's own tables ('== CHARG 1 ...', N compression-
  // positive, keta station signs); true = the tables as ST1 prints them
  // (keta/st1_edit.h: block headers, Fortran formats, ST1 sign conventions),
  // so a result file compares line for line with an ST1 one. `keta session`
  // turns it on.
  bool st1_edition = false;
  // SURCH convoys (Bc, Bt, Mc, TS, LM ...) moved on the PAS_CONV grid ONLY,
  // as ST1 does (manual p.147), instead of the grid PLUS the exact breakpoints
  // of the influence line (keta's default: the exact extremum, always equal or
  // above ST1's). Off by default; `keta session --pas-conv`.
  bool pas_conv_grid_only = false;
  // Creep of a PHASAGE under a varying stress: true (default) = the exact
  // linear superposition with the law itself (keta/phasing.h
  // CreepEngine::Exact); false = ST1's own equivalent-time method
  // (CreepEngine::St1, identified on ST1's own results, exact only
  // under a constant stress). `keta session --fluage-st1`.
  bool creep_superposition = true;
  // Deferred losses of the tendons in a PHASAGE: true (default) = solved by
  // the structure (keta/phasing.h PhasingOptions::exact_losses: the free
  // strain of the concrete at the tendon, the relaxation on the local stress,
  // the loss of each tendon fed back on all of them until it converges);
  // false = EN 1992-1-1 eq. 5.46 restarted on each interval, as ST1.
  // `keta session --pertes-st1`.
  bool exact_losses = true;
  // Workers of the PHASAGE station recovery (keta/parallel.h): 0 (default) =
  // hardware concurrency, 1 = sequential. The results are bit-identical for
  // any count (each bar writes its own record, reductions run in bar order).
  // `keta session --threads N`.
  int threads = 0;
};

St1Import import_st1(const std::string& text, const St1Options& opt = {});
St1Import import_st1_file(const std::string& path, const St1Options& opt = {});

}  // namespace keta

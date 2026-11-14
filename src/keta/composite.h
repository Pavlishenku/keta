// Composite bars made of section PARTS (steel girder + concrete slab, a box
// cast in several pours...): each part has its own material, section and
// position, and is cast (activated) in its own phase.
//
// The member is EXACT for full interaction: plane sections stay plane across
// the parts, so the active parts behave as ONE beam whose section is their
// sum transformed to a reference modulus — about its elastic centroid
//     A* = sum n_i A_i,  z_c = sum n_i A_i z_i / A*,
//     Iy* = sum n_i (Iy_i + A_i (z_i - z_c)^2),  Iz* = sum n_i Iz_i,
//     J* = sum (G_i / G_ref) J_i,  Av* = sum (G_i / G_ref) Av_i,
// n_i = E_i / E_ref (E_ref: the first active part) — carried by a prismatic
// element on rigid arms from the bar's node line to that centroid
// (Bar::ecc). With such an element the member is nodally exact (no stacking
// of eccentric bars, whose axial and bending fields are incompatible inside
// an element). The parts lie on the bar's local z axis (y_i = 0) with their
// principal axes parallel to the bar's, so the section stays principal.
//
// Phasing (keta/phasing.h): a part is cast stress-free in its phase; the
// stress of each part is kept per station from the generalized strains of
// every increment, and the creep and shrinkage of its material act on it
// (exact Boltzmann superposition per part, the restraint of the other parts
// giving the eigenstresses). Results of a composite bar are reported about
// its NODE LINE (the centroid moves when a part is cast). Outside the
// phasing, `composite_lowered` gives the model with every part active.
#pragma once
#include "keta/model.h"

namespace keta {

struct CompositeState {
  double E_ref = 0, nu_ref = 0, rho = 0;
  SectionProps props;           // transformed, about the centroid
  double zc = 0;                // centroid, local z from the node line
  std::vector<double> E;        // modulus of each part (0 = inactive)
};

// The transformed section of the active parts of composite bar `bar`, with
// the moduli E_parts (one per part; <= 0 = inactive).
CompositeState composite_state(const Model& model, Id bar, const std::vector<double>& E_parts);

// Replace, in `w`, bar `bar` by its transformed member (a new section, a new
// material carrying E_ref / nu_ref / the mass of the active parts, the arms
// to the centroid). `w` must hold the bar and the materials.
void lower_composite_bar(Model& w, Id bar, const CompositeState& s);

// The model with every composite bar lowered with ALL its parts active at
// their materials' moduli (Model::composite emptied): the finished member for
// the analyses outside the phasing.
Model composite_lowered(const Model& model);

}  // namespace keta

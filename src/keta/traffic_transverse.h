// Transverse distribution of traffic loads.
//
//  * Courbon's method (rigid cross-beams, girders without torsional
//    stiffness, all of the same span): a load P at transverse position y is
//    shared as R_i = P w_i [1/W + (y - y_c)(y_i - y_c) / J], W = sum w_j, the
//    stiffness centre y_c = sum w_j y_j / W, J = sum w_j (y_j - y_c)^2 — the
//    rigid-body (translation + rotation) of the deck on springs w_j (the
//    girders' bending stiffness EI_j / L^3 at mid-span). Equal girders: K_i =
//    1/n + e x_i / sum x_j^2, the textbook form. K_i is LINEAR in y, so its mean
//    over a lane is its value at the lane centre.
//  * A separable influence surface eta(x, y) = IL(x) K(y), K a table of the
//    transverse coefficient (any n values, linear between): every surface
//    placement (patches, wheels, UDL zones) then runs on a longitudinal line
//    and a transverse table. On the grid of the line samples x and the table
//    abscissae y the bilinear surface IS the product of the two piecewise-
//    linear interpolants (tensor-product interpolation), exact.
#pragma once
#include <vector>

#include "keta/influence_surface.h"

namespace keta {

// Courbon coefficients K_i(y) of all girders (positions `girder_y`, relative
// stiffnesses `weights`, empty = equal) for a unit load at y.
std::vector<double> courbon_coefficients(const std::vector<double>& girder_y, const std::vector<double>& weights,
                                         double y);

// eta(x_i, y_j) = il.eta[i] * k[j] on the grid (il.s) x (ky), ky strictly increasing.
InfluenceSurface separable_surface(const InfluenceLine& il, const std::vector<double>& ky,
                                   const std::vector<double>& k);

}  // namespace keta

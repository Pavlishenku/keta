// Scalar field conduction on continuum elements (keta/continuum.h): heat
// (temperature) or groundwater seepage (hydraulic head) — the same
// equation, div(k grad u) + s = c du/dt.
//
// Elements: plane (thickness t), axisymmetric (per radian: r dr dz) and
// solid continua of the model, the same shape functions and Gauss rules as
// the mechanics. Conductivity per material, orthotropic along the global
// axes and linear in the field: k_i(u) = k_i + k1 u;
// volumetric capacity c for the transient. Boundary data: fixed values at
// nodes, flux q (input per unit area) and convection h (u - u_env) on
// element faces (2-D: edge a -> a + 1; solid: the hexahedron faces of
// continuum_face_load), point sources, and conductance LINKS between two
// nodes or from a node to an environment value (thermal gap links).
// Steady: Newton on R(u) = K(u) u - f with the exact tangent (k(u) term
// included) to round-off; transient: backward Euler with the consistent
// capacity matrix, Newton at each step. The "reactions" are the flows the
// fixed nodes supply (positive: into the body) — e.g. the discharge of a
// seepage problem.
//
// continuum_thermal_load: the consistent nodal forces of a temperature
// field (nodal values, interpolated) on an element, int B' D eps0 dV,
// eps0 = alpha (T - T0) (plane strain: the (1 + nu) of the restrained
// thickness; axisymmetric: the hoop component too) — the coupling of a
// conduction result to the mechanics.
#pragma once
#include <map>
#include <vector>

#include "keta/model.h"

namespace keta {

struct ConductionMaterial {
  double kx = 1, ky = 1, kz = 1;  // conductivity along the global axes
  double k1 = 0;                  // k_i(u) = k_i + k1 u
  double capacity = 0;            // rho c (transient)
};
struct ConductionProblem {
  std::map<Id, ConductionMaterial> materials;  // material id -> its conduction
  std::map<Id, double> fixed;                   // node -> imposed value
  struct Face {
    Id element = 0;
    int face = 0;
    double q = 0;               // flux INTO the body per unit area
    double h = 0, u_env = 0;    // convection h (u - u_env) out of the body
  };
  std::vector<Face> faces;
  std::map<Id, double> sources;  // node -> point input
  struct Link {
    Id a = 0, b = 0;  // b = 0: to the environment value u_env
    double k = 0, u_env = 0;
  };
  std::vector<Link> links;
  double u_initial = 0;                // transient: the initial field
  std::map<Id, double> initial;        // ... or per node
};
struct ConductionOptions {
  std::vector<double> times;  // empty: steady state; else the output times (from 0), backward Euler
  int substeps = 1;           // steps between two output times
  int max_iter = 50;
  double tol = 1e-12;
};
struct ConductionResult {
  std::vector<std::map<Id, double>> u;          // per output (steady: one), node -> value
  std::vector<std::map<Id, double>> reactions;  // the flows the fixed nodes supply
  int iterations = 0;
};
ConductionResult solve_conduction(const Model& model, const ConductionProblem& p, const ConductionOptions& opt = {});

std::vector<Load> continuum_thermal_load(const Model& model, Id element, const std::map<Id, double>& T, double T0);

}  // namespace keta

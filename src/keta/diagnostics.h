// Model diagnostics: localizes what validate() cannot reject — coincident
// nodes, zero-length bars, isolated nodes, floating islands (union-find over
// bars and supports), degenerate properties, and the factorization-gated
// mechanism probe: the reduced K is factored with the solver's own LDLT; the
// EXISTENCE verdict is the factorization's (zero/negative pivots), and only
// then a low-eigenmode probe localizes the mechanism to nodes/DOFs.
// Never throws on a broken model — a diagnosis IS the product.
#pragma once
#include "keta/model.h"

namespace keta {

enum class DiagnosticKind {
  CoincidentNodes,
  ZeroLengthBar,
  IsolatedNode,
  FloatingIsland,
  DegenerateProperty,
  Mechanism,
  IndefiniteStiffness,
};
enum class Severity { Error, Warning };

struct Diagnostic {
  DiagnosticKind kind;
  Severity severity;
  std::vector<Id> nodes;
  std::vector<Id> bars;
  int dof = -1;
  std::string message;
  std::string suggestion;
};

std::vector<Diagnostic> diagnose(const Model& model, bool mechanism_probe = true,
                                 double coincident_tol = 1e-6);

}  // namespace keta

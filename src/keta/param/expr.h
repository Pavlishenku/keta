// keta::param — the formula language of the Expression component: a
// recursive-descent Pratt parser to a small AST, evaluated over six variable
// slots x y z u v w.
//
// Grammar, lowest to highest precedence:
//   comparisons  < > <= >= == !=   -> 1.0 / 0.0   (left, binding power 1)
//   + -                                          (left, 3)
//   * / %                                        (left, 5)
//   unary - +                                    (prefix, 7)
//   ^ power                                      (RIGHT-associative, 8)
//   atoms: numbers (1, .5, 2e-3), (expr), constants pi e, the six variables,
//          calls f(...): sin cos tan asin acos atan sqrt abs floor ceil round
//          ln log(base 10) exp (1 argument); atan2 min max pow (2 arguments).
//
// NaN / inf are NEVER silent: knowable domain errors are rejected before
// computing (division / modulo by zero, sqrt of a negative, ln/log of a
// non-positive, asin/acos outside [-1, 1], a non-integer power of a negative
// base, zero to a negative power) and every operator / function result is
// finiteness-checked; failures come back as ExprError naming the operation.
#pragma once
#include <memory>
#include <stdexcept>
#include <string>

namespace keta::param::expr {

// Parse failure, with the byte position it was detected at.
struct ParseError : std::runtime_error {
  std::size_t pos;
  ParseError(std::size_t p, const std::string& m) : std::runtime_error(m), pos(p) {}
};

// Evaluation failure (domain or non-finite result).
struct EvalError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct Node;  // AST node (opaque)

class Expr {
 public:
  static Expr parse(const std::string& source);  // throws ParseError
  double eval(const double vars[6]) const;       // throws EvalError
 private:
  std::shared_ptr<const Node> root_;
};

// Shared domain guards and finiteness check (also used by the Maths nodes, so
// every message string has one source). Each returns nullptr when admissible.
const char* guard_sqrt(double x);
const char* guard_pow(double a, double b);
const char* guard_div(double divisor);
const char* guard_mod(double divisor);
// x when finite, else throws EvalError("<op> produced an undefined result
// (NaN)" / "numeric overflow in <op>").
double finite_or_throw(double x, const std::string& op);

}  // namespace keta::param::expr

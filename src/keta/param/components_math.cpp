// Maths components.
#include <algorithm>
#include <cctype>
#include <cmath>

#include "keta/param/components.h"
#include "keta/param/expr.h"

namespace keta::param {

namespace {

const char* const kCategory = "Maths";

double num(const Value& v) {
  if (auto x = v.as_number()) return *x;
  throw ComponentError("expected a number, got " + v.type_tag());
}
bool boolean(const Value& v) {
  if (auto b = v.as_bool()) return *b;
  throw ComponentError("expected a boolean, got " + v.type_tag());
}

// The node-side adapters of the shared guards / finiteness check.
double finite(double x, const std::string& op) {
  try {
    return expr::finite_or_throw(x, op);
  } catch (const expr::EvalError& e) {
    throw ComponentError(e.what());
  }
}
Vec3 finite_vec3(const Vec3& v, const std::string& op) {
  finite(v.x(), op + " (x component)");
  finite(v.y(), op + " (y component)");
  finite(v.z(), op + " (z component)");
  return v;
}
void guard(const char* why) {
  if (why) throw ComponentError(why);
}

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

struct UnaryInfo {
  const char* id;
  const char* name;
  const char* desc;
};
UnaryInfo unary_info(UnaryOp op) {
  switch (op) {
    case UnaryOp::Negate: return {"keta.negate", "Negate", "Negate a number (-x)"};
    case UnaryOp::Abs: return {"keta.abs", "Absolute", "Absolute value |x|"};
    case UnaryOp::Sqrt: return {"keta.sqrt", "Square Root", "Square root of x; a negative input errors the node"};
    case UnaryOp::Sin: return {"keta.sin", "Sine", "Sine of x (radians)"};
    case UnaryOp::Cos: return {"keta.cos", "Cosine", "Cosine of x (radians)"};
    case UnaryOp::Tan: return {"keta.tan", "Tangent", "Tangent of x (radians)"};
    case UnaryOp::Floor: return {"keta.floor", "Floor", "Round x down to the nearest integer"};
    case UnaryOp::Ceil: return {"keta.ceil", "Ceiling", "Round x up to the nearest integer"};
    case UnaryOp::Round: return {"keta.round", "Round", "Round x to the nearest integer"};
  }
  return {"", "", ""};
}

double apply_unary(UnaryOp op, double x) {
  switch (op) {
    case UnaryOp::Negate: return -x;
    case UnaryOp::Abs: return std::fabs(x);
    case UnaryOp::Sqrt: guard(expr::guard_sqrt(x)); return std::sqrt(x);
    case UnaryOp::Sin: return std::sin(x);
    case UnaryOp::Cos: return std::cos(x);
    case UnaryOp::Tan: return std::tan(x);
    case UnaryOp::Floor: return std::floor(x);
    case UnaryOp::Ceil: return std::ceil(x);
    case UnaryOp::Round: return std::round(x);  // half away from zero
  }
  return x;
}

struct BinaryInfo {
  const char* id;
  const char* name;
  const char* desc;
  const char* operand_tag;  // "number" | "bool"
  const char* result_tag;
};
BinaryInfo binary_info(BinaryOp op) {
  switch (op) {
    case BinaryOp::Add: return {"keta.add", "Add", "A + B (numbers, vectors, or point + vector)", "number", "number"};
    case BinaryOp::Subtract:
      return {"keta.subtract", "Subtract", "A - B (numbers, vectors, point - vector, or point - point)", "number", "number"};
    case BinaryOp::Multiply: return {"keta.multiply", "Multiply", "A x B (numbers, or vector x number)", "number", "number"};
    case BinaryOp::Divide:
      return {"keta.divide", "Divide", "A / B (numbers, or vector / number); dividing by zero errors the node", "number", "number"};
    case BinaryOp::Power:
      return {"keta.power", "Power",
              "A raised to the power B; a non-integer power of a negative base, or a zero base with a negative "
              "exponent, errors the node",
              "number", "number"};
    case BinaryOp::Modulo: return {"keta.modulo", "Modulo", "Remainder of A divided by B; a zero divisor errors the node", "number", "number"};
    case BinaryOp::Larger: return {"keta.larger", "Larger", "True when A > B", "number", "bool"};
    case BinaryOp::Smaller: return {"keta.smaller", "Smaller", "True when A < B", "number", "bool"};
    case BinaryOp::And: return {"keta.gate_and", "And", "Logical AND of A and B", "bool", "bool"};
    case BinaryOp::Or: return {"keta.gate_or", "Or", "Logical OR of A and B", "bool", "bool"};
  }
  return {"", "", "", "", ""};
}

Value combine(BinaryOp op, const Value& a, const Value& b) {
  auto x = a.as_number(), y = b.as_number();
  const auto* pa = a.get<Point>();
  const auto* pb = b.get<Point>();
  const auto* va = a.get<Vector>();
  const auto* vb = b.get<Vector>();
  switch (op) {
    case BinaryOp::Add:
      if (x && y) return Value::number(finite(*x + *y, "add"));
      if (va && vb) return Value::of(Vector{finite_vec3(va->v + vb->v, "add")});
      if (pa && vb) return Value::of(Point{finite_vec3(pa->p + vb->v, "add")});
      if (va && pb) return Value::of(Point{finite_vec3(pb->p + va->v, "add")});
      if (pa && pb) return Value::of(Point{finite_vec3(pa->p + pb->p, "add")});
      throw ComponentError("cannot add " + a.type_tag() + " and " + b.type_tag());
    case BinaryOp::Subtract:
      if (x && y) return Value::number(finite(*x - *y, "subtract"));
      if (va && vb) return Value::of(Vector{finite_vec3(va->v - vb->v, "subtract")});
      if (pa && vb) return Value::of(Point{finite_vec3(pa->p - vb->v, "subtract")});
      if (pa && pb) return Value::of(Vector{finite_vec3(pa->p - pb->p, "subtract")});  // connecting vector
      throw ComponentError("cannot subtract " + b.type_tag() + " from " + a.type_tag());
    case BinaryOp::Multiply:
      if (x && y) return Value::number(finite(*x * *y, "multiply"));
      if (va && y) return Value::of(Vector{finite_vec3(va->v * *y, "multiply")});
      if (x && vb) return Value::of(Vector{finite_vec3(vb->v * *x, "multiply")});
      throw ComponentError("cannot multiply " + a.type_tag() + " and " + b.type_tag());
    case BinaryOp::Divide:
      if (x && y) {
        guard(expr::guard_div(*y));
        return Value::number(finite(*x / *y, "divide"));
      }
      if (va && y) {
        // True component-wise division (not v * (1/s): a subnormal s would
        // give inf and a false NaN).
        guard(expr::guard_div(*y));
        return Value::of(Vector{finite_vec3(Vec3(va->v.x() / *y, va->v.y() / *y, va->v.z() / *y), "divide")});
      }
      throw ComponentError("cannot divide " + a.type_tag() + " by " + b.type_tag());
    case BinaryOp::Power: {
      double p = num(a), q = num(b);
      guard(expr::guard_pow(p, q));
      return Value::number(finite(std::pow(p, q), "power"));
    }
    case BinaryOp::Modulo: {
      double p = num(a), q = num(b);
      guard(expr::guard_mod(q));
      return Value::number(finite(std::fmod(p, q), "modulo"));
    }
    case BinaryOp::Larger: return Value::boolean(num(a) > num(b));
    case BinaryOp::Smaller: return Value::boolean(num(a) < num(b));
    case BinaryOp::And: return Value::boolean(boolean(a) && boolean(b));
    case BinaryOp::Or: return Value::boolean(boolean(a) || boolean(b));
  }
  return Value();
}

}  // namespace

// ---------------- unary ----------------

ComponentInfo UnaryMathNode::info() const {
  auto u = unary_info(op);
  return {u.id, u.name, kCategory, u.desc};
}
std::vector<PinSpec> UnaryMathNode::inputs() const { return {PinSpec::item("x", "number", "Input value")}; }
std::vector<PinSpec> UnaryMathNode::outputs() const { return {PinSpec::item("R", "number", "Result")}; }
Outputs<Value> UnaryMathNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  const std::string opname = lower(unary_info(op).name);
  return {map_items({{"x", &in[0]}}, [&](const std::vector<const Value*>& it) {
    return Value::number(finite(apply_unary(op, num(*it[0])), opname));
  })};
}

// ---------------- binary ----------------

ComponentInfo BinaryMathNode::info() const {
  auto b = binary_info(op);
  return {b.id, b.name, kCategory, b.desc};
}
std::vector<PinSpec> BinaryMathNode::inputs() const {
  auto b = binary_info(op);
  return {PinSpec::item("A", b.operand_tag, "First operand"), PinSpec::item("B", b.operand_tag, "Second operand")};
}
std::vector<PinSpec> BinaryMathNode::outputs() const { return {PinSpec::item("R", binary_info(op).result_tag, "Result")}; }
Outputs<Value> BinaryMathNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {map_items({{"A", &in[0]}, {"B", &in[1]}},
                    [&](const std::vector<const Value*>& it) { return combine(op, *it[0], *it[1]); })};
}

// ---------------- Not / Equals / Range / Remap / Mass Addition ----------------

ComponentInfo NotNode::info() const { return {"keta.gate_not", "Not", kCategory, "Logical NOT of a boolean"}; }
std::vector<PinSpec> NotNode::inputs() const { return {PinSpec::item("x", "bool", "Input value")}; }
std::vector<PinSpec> NotNode::outputs() const { return {PinSpec::item("R", "bool", "Result")}; }
Outputs<Value> NotNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  return {map_items({{"x", &in[0]}}, [](const std::vector<const Value*>& it) { return Value::boolean(!boolean(*it[0])); })};
}

ComponentInfo EqualsNode::info() const {
  return {"keta.equals", "Equals", kCategory, "True when |A - B| <= E (E optional, default 0)"};
}
std::vector<PinSpec> EqualsNode::inputs() const {
  return {PinSpec::item("A", "number", "First operand"), PinSpec::item("B", "number", "Second operand"),
          PinSpec::item("E", "number", "Tolerance (default 0)").opt()};
}
std::vector<PinSpec> EqualsNode::outputs() const { return {PinSpec::item("R", "bool", "Result")}; }
Outputs<Value> EqualsNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree e = number_or(input_or_empty(in, 2), 0.0);
  return {map_items({{"A", &in[0]}, {"B", &in[1]}, {"E", &e}}, [](const std::vector<const Value*>& it) {
    return Value::boolean(std::fabs(num(*it[0]) - num(*it[1])) <= num(*it[2]));
  })};
}

ComponentInfo RangeNode::info() const {
  return {"keta.range", "Range", kCategory, "N+1 evenly spaced numbers over the domain [D0, D1]"};
}
std::vector<PinSpec> RangeNode::inputs() const {
  return {PinSpec::item("D0", "number", "Domain start (default 0)").opt(),
          PinSpec::item("D1", "number", "Domain end (default 1)").opt(),
          PinSpec::item("N", "int", "Number of steps (default 10)").opt()};
}
std::vector<PinSpec> RangeNode::outputs() const { return {PinSpec::list("R", "number", "The evenly spaced numbers")}; }
Outputs<Value> RangeNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  double d0 = first_number(in[0]).value_or(0.0), d1 = first_number(in[1]).value_or(1.0);
  std::int64_t n = first_int(in[2]).value_or(10);
  if (n < 0 || n > 1'000'000) throw ComponentError("N " + std::to_string(n) + " out of range 0..=1e6");
  double step = n == 0 ? 0.0 : finite((d1 - d0) / static_cast<double>(n), "range");
  Tree t;
  for (std::int64_t i = 0; i <= n; ++i)
    t.append(Path::single(0), Value::number(finite(d0 + step * static_cast<double>(i), "range")));
  return {t};
}

ComponentInfo RemapNode::info() const {
  return {"keta.remap", "Remap", kCategory, "Remap V from the source domain [S0, S1] to the target domain [T0, T1]"};
}
std::vector<PinSpec> RemapNode::inputs() const {
  return {PinSpec::item("V", "number", "Value to remap"),
          PinSpec::item("S0", "number", "Source domain start (default 0)").opt(),
          PinSpec::item("S1", "number", "Source domain end (default 1)").opt(),
          PinSpec::item("T0", "number", "Target domain start (default 0)").opt(),
          PinSpec::item("T1", "number", "Target domain end (default 1)").opt()};
}
std::vector<PinSpec> RemapNode::outputs() const { return {PinSpec::item("R", "number", "Remapped value")}; }
Outputs<Value> RemapNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree s0 = number_or(input_or_empty(in, 1), 0.0), s1 = number_or(input_or_empty(in, 2), 1.0);
  Tree t0 = number_or(input_or_empty(in, 3), 0.0), t1 = number_or(input_or_empty(in, 4), 1.0);
  return {map_items({{"V", &in[0]}, {"S0", &s0}, {"S1", &s1}, {"T0", &t0}, {"T1", &t1}},
                    [](const std::vector<const Value*>& it) {
                      double v = num(*it[0]), a0 = num(*it[1]), a1 = num(*it[2]), b0 = num(*it[3]), b1 = num(*it[4]);
                      double den = a1 - a0;
                      double r = std::fabs(den) < 1e-15 ? b0 : b0 + (v - a0) / den * (b1 - b0);
                      return Value::number(finite(r, "remap"));
                    })};
}

ComponentInfo MassAdditionNode::info() const {
  return {"keta.mass_addition", "Mass Addition", kCategory,
          "Sum every branch of numbers (S = total, Pr = running partial sums)"};
}
std::vector<PinSpec> MassAdditionNode::inputs() const { return {PinSpec::list("L", "number", "List of numbers to sum").opt()}; }
std::vector<PinSpec> MassAdditionNode::outputs() const {
  return {PinSpec::item("S", "number", "Sum of each branch"), PinSpec::list("Pr", "number", "Running partial sums of each branch")};
}
Outputs<Value> MassAdditionNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  Tree sums, partials;
  for (const auto& [path, branch] : in[0].branches()) {
    double sum = 0;
    auto& pr = partials.ensure_path(path);
    for (const auto& o : branch) {
      if (!o) continue;
      sum = finite(sum + num(*o), "mass addition");
      pr.emplace_back(Value::number(sum));
    }
    sums.append(path, Value::number(sum));
  }
  return {sums, partials};
}

// ---------------- Expression ----------------

ComponentInfo ExpressionNode::info() const {
  return {"keta.expression", "Expression", kCategory, "Evaluate a textual formula over the inputs x, y, z, u, v, w"};
}
std::vector<PinSpec> ExpressionNode::inputs() const {
  std::vector<PinSpec> pins;
  for (const char* v : {"x", "y", "z", "u", "v", "w"}) pins.push_back(PinSpec::item(v, "number", "Variable (default 0)").opt());
  return pins;
}
std::vector<PinSpec> ExpressionNode::outputs() const { return {PinSpec::item("R", "number", "The evaluated result")}; }
Outputs<Value> ExpressionNode::solve(const SolveCtx&, const std::vector<Tree>& in) const {
  expr::Expr ast;
  try {
    ast = expr::Expr::parse(formula);
  } catch (const expr::ParseError& e) {
    throw ComponentError("expression parse error at " + std::to_string(e.pos) + ": " + e.what());
  }
  static const char* const names[6] = {"x", "y", "z", "u", "v", "w"};
  std::vector<Tree> defaulted;
  defaulted.reserve(6);
  for (std::size_t k = 0; k < 6; ++k) defaulted.push_back(number_or(input_or_empty(in, k), 0.0));
  std::vector<std::pair<std::string, const Tree*>> pins;
  for (std::size_t k = 0; k < 6; ++k) pins.emplace_back(names[k], &defaulted[k]);
  return {map_items(pins, [&](const std::vector<const Value*>& it) {
    double vars[6];
    for (std::size_t k = 0; k < 6; ++k) vars[k] = num(*it[k]);
    try {
      // Operators finiteness-check themselves; a bare literal (1e400) does not,
      // so the top-level result is checked too.
      return Value::number(expr::finite_or_throw(ast.eval(vars), "expression"));
    } catch (const expr::EvalError& e) {
      throw ComponentError("expression '" + formula + "': " + e.what());
    }
  })};
}
json::Value ExpressionNode::params() const {
  json::Value p = json::Value::object();
  p.set("formula", formula);
  return p;
}
void ExpressionNode::set_params(const json::Value& p) {
  if (const auto* f = p.get("formula"); f && f->is_string()) formula = f->as_string();
}

}  // namespace keta::param

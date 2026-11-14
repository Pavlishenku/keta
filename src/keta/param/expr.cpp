#include "keta/param/expr.h"

#include <charconv>
#include <cctype>
#include <cmath>
#include <optional>
#include <vector>

namespace keta::param::expr {

const char* guard_sqrt(double x) { return x < 0.0 ? "square root of a negative number" : nullptr; }
const char* guard_pow(double a, double b) {
  if (a == 0.0 && b < 0.0) return "zero raised to a negative power";
  if (std::isinf(b)) return nullptr;  // IEEE pow defines it; overflow is caught after
  if (a < 0.0 && std::trunc(b) != b) return "invalid power (non-integer exponent of a negative base)";
  return nullptr;
}
const char* guard_div(double d) { return d == 0.0 ? "division by zero" : nullptr; }
const char* guard_mod(double d) { return d == 0.0 ? "modulo by zero" : nullptr; }

double finite_or_throw(double x, const std::string& op) {
  if (std::isfinite(x)) return x;
  if (std::isnan(x)) throw EvalError(op + " produced an undefined result (NaN)");
  throw EvalError("numeric overflow in " + op);
}

namespace {

void check(const char* why) {
  if (why) throw EvalError(why);
}

enum class Op { Add, Sub, Mul, Div, Mod, Pow, Lt, Gt, Le, Ge, Eq, Ne };

}  // namespace

struct Node {
  enum class Kind { Num, Var, Neg, Bin, Fn1, Fn2 } kind = Kind::Num;
  double num = 0;
  int var = 0;
  Op op = Op::Add;
  std::string fname;
  double (*f1)(double) = nullptr;
  double (*f2)(double, double) = nullptr;
  std::shared_ptr<const Node> a, b;

  double eval(const double v[6]) const {
    switch (kind) {
      case Kind::Num:
        return num;
      case Kind::Var:
        return v[var];
      case Kind::Neg:
        return -a->eval(v);
      case Kind::Fn1: {
        double x = a->eval(v);
        if (fname == "sqrt") check(guard_sqrt(x));
        if ((fname == "asin" || fname == "acos") && !(x >= -1.0 && x <= 1.0))
          throw EvalError(fname + ": argument outside [-1, 1]");
        if ((fname == "ln" || fname == "log") && x <= 0.0) throw EvalError("logarithm of a non-positive number");
        return finite_or_throw(f1(x), fname);
      }
      case Kind::Fn2: {
        double x = a->eval(v), y = b->eval(v);
        if (fname == "pow") check(guard_pow(x, y));
        return finite_or_throw(f2(x, y), fname);
      }
      case Kind::Bin: {
        double x = a->eval(v), y = b->eval(v);
        switch (op) {
          case Op::Add: return finite_or_throw(x + y, "add");
          case Op::Sub: return finite_or_throw(x - y, "subtract");
          case Op::Mul: return finite_or_throw(x * y, "multiply");
          case Op::Div: check(guard_div(y)); return finite_or_throw(x / y, "divide");
          case Op::Mod: check(guard_mod(y)); return finite_or_throw(std::fmod(x, y), "modulo");
          case Op::Pow: check(guard_pow(x, y)); return finite_or_throw(std::pow(x, y), "power");
          case Op::Lt: return x < y ? 1.0 : 0.0;
          case Op::Gt: return x > y ? 1.0 : 0.0;
          case Op::Le: return x <= y ? 1.0 : 0.0;
          case Op::Ge: return x >= y ? 1.0 : 0.0;
          case Op::Eq: return x == y ? 1.0 : 0.0;
          case Op::Ne: return x != y ? 1.0 : 0.0;
        }
      }
    }
    return 0.0;
  }
};

namespace {

using NodeP = std::shared_ptr<const Node>;

enum class T { Num, Ident, Plus, Minus, Star, Slash, Percent, Caret, LParen, RParen, Comma, Lt, Gt, Le, Ge, EqEq, Ne };

struct Tok {
  T t;
  std::size_t pos;
  double num = 0;
  std::string text;
};

const char* const kVars[6] = {"x", "y", "z", "u", "v", "w"};

std::string describe(const Tok& k) {
  switch (k.t) {
    case T::Num: return "number";
    case T::Ident: return "'" + k.text + "'";
    case T::LParen: return "'('";
    case T::RParen: return "')'";
    case T::Comma: return "','";
    default: return "operator";
  }
}

std::vector<Tok> tokenize(const std::string& s) {
  std::vector<Tok> out;
  std::size_t i = 0;
  while (i < s.size()) {
    char c = s[i];
    if (std::isspace(static_cast<unsigned char>(c))) {
      ++i;
      continue;
    }
    std::size_t start = i;
    char next = i + 1 < s.size() ? s[i + 1] : '\0';
    auto one = [&](T t) {
      out.push_back({t, start, 0, {}});
      ++i;
    };
    switch (c) {
      case '+': one(T::Plus); continue;
      case '-': one(T::Minus); continue;
      case '*': one(T::Star); continue;
      case '/': one(T::Slash); continue;
      case '%': one(T::Percent); continue;
      case '^': one(T::Caret); continue;
      case '(': one(T::LParen); continue;
      case ')': one(T::RParen); continue;
      case ',': one(T::Comma); continue;
      case '<':
        out.push_back({next == '=' ? T::Le : T::Lt, start, 0, {}});
        i += next == '=' ? 2 : 1;
        continue;
      case '>':
        out.push_back({next == '=' ? T::Ge : T::Gt, start, 0, {}});
        i += next == '=' ? 2 : 1;
        continue;
      case '=':
        if (next != '=') throw ParseError(start, "expected '==' (single '=' is not valid)");
        out.push_back({T::EqEq, start, 0, {}});
        i += 2;
        continue;
      case '!':
        if (next != '=') throw ParseError(start, "expected '!=' after '!'");
        out.push_back({T::Ne, start, 0, {}});
        i += 2;
        continue;
      default:
        break;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      std::size_t j = i;
      while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
      if (j < s.size() && (s[j] == 'e' || s[j] == 'E')) {
        ++j;
        if (j < s.size() && (s[j] == '+' || s[j] == '-')) ++j;
        while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
      }
      std::string text = s.substr(i, j - i);
      double v = 0;
      auto [p, ec] = std::from_chars(text.data(), text.data() + text.size(), v);
      if (ec != std::errc() || p != text.data() + text.size()) {
        // from_chars reports out_of_range for huge literals; read them as +inf
        // so the finiteness check names them ("1e400" -> overflow).
        if (ec == std::errc::result_out_of_range && p == text.data() + text.size())
          v = std::strtod(text.c_str(), nullptr);
        else
          throw ParseError(start, "invalid number '" + text + "'");
      }
      out.push_back({T::Num, start, v, {}});
      i = j;
      continue;
    }
    if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
      std::size_t j = i;
      while (j < s.size() && (std::isalnum(static_cast<unsigned char>(s[j])) || s[j] == '_')) ++j;
      out.push_back({T::Ident, start, 0, s.substr(i, j - i)});
      i = j;
      continue;
    }
    throw ParseError(start, std::string("unexpected character '") + c + "'");
  }
  return out;
}

std::optional<std::pair<int, int>> infix_bp(T t) {
  switch (t) {
    case T::EqEq: case T::Ne: case T::Lt: case T::Gt: case T::Le: case T::Ge: return std::pair{1, 2};
    case T::Plus: case T::Minus: return std::pair{3, 4};
    case T::Star: case T::Slash: case T::Percent: return std::pair{5, 6};
    case T::Caret: return std::pair{8, 7};  // right-associative
    default: return std::nullopt;
  }
}
constexpr int kPrefixBp = 7;

Op to_op(T t) {
  switch (t) {
    case T::Plus: return Op::Add;
    case T::Minus: return Op::Sub;
    case T::Star: return Op::Mul;
    case T::Slash: return Op::Div;
    case T::Percent: return Op::Mod;
    case T::Caret: return Op::Pow;
    case T::Lt: return Op::Lt;
    case T::Gt: return Op::Gt;
    case T::Le: return Op::Le;
    case T::Ge: return Op::Ge;
    case T::EqEq: return Op::Eq;
    default: return Op::Ne;
  }
}

NodeP leaf(double v) {
  auto n = std::make_shared<Node>();
  n->num = v;
  return n;
}

double c_round(double x) { return std::round(x); }  // half away from zero

struct Parser {
  std::vector<Tok> toks;
  std::size_t i = 0, eof = 0;

  const Tok* peek() const { return i < toks.size() ? &toks[i] : nullptr; }

  void expect(T want, const char* what) {
    if (i >= toks.size()) throw ParseError(eof, std::string("expected ") + what + ", found end of input");
    if (toks[i].t != want)
      throw ParseError(toks[i].pos, std::string("expected ") + what + ", found " + describe(toks[i]));
    ++i;
  }

  NodeP expr_bp(int min_bp) {
    if (i >= toks.size()) throw ParseError(eof, "unexpected end of input");
    Tok tok = toks[i++];
    NodeP lhs;
    switch (tok.t) {
      case T::Num: lhs = leaf(tok.num); break;
      case T::Minus: {
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::Neg;
        n->a = expr_bp(kPrefixBp);
        lhs = n;
        break;
      }
      case T::Plus: lhs = expr_bp(kPrefixBp); break;
      case T::LParen:
        lhs = expr_bp(0);
        expect(T::RParen, "')'");
        break;
      case T::Ident: lhs = ident_or_call(tok); break;
      default: throw ParseError(tok.pos, "unexpected token " + describe(tok));
    }
    while (const Tok* p = peek()) {
      auto bp = infix_bp(p->t);
      if (!bp || bp->first < min_bp) break;
      T op = p->t;
      ++i;
      auto n = std::make_shared<Node>();
      n->kind = Node::Kind::Bin;
      n->op = to_op(op);
      n->a = lhs;
      n->b = expr_bp(bp->second);
      lhs = n;
    }
    return lhs;
  }

  NodeP ident_or_call(const Tok& tok) {
    const std::string& name = tok.text;
    if (const Tok* p = peek(); p && p->t == T::LParen) {
      ++i;
      std::vector<NodeP> args;
      if (const Tok* q = peek(); !q || q->t != T::RParen) {
        for (;;) {
          args.push_back(expr_bp(0));
          if (const Tok* c = peek(); c && c->t == T::Comma) {
            ++i;
            continue;
          }
          break;
        }
      }
      expect(T::RParen, "')'");
      return call(name, std::move(args), tok.pos);
    }
    if (name == "pi") return leaf(3.14159265358979323846);
    if (name == "e") return leaf(2.71828182845904523536);
    for (int k = 0; k < 6; ++k)
      if (name == kVars[k]) {
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::Var;
        n->var = k;
        return n;
      }
    throw ParseError(tok.pos, "unknown variable '" + name + "'");
  }

  static NodeP call(const std::string& name, std::vector<NodeP> args, std::size_t pos) {
    struct F1 {
      const char* n;
      double (*f)(double);
    };
    static const F1 f1s[] = {{"sin", [](double x) { return std::sin(x); }},   {"cos", [](double x) { return std::cos(x); }},
                             {"tan", [](double x) { return std::tan(x); }},   {"asin", [](double x) { return std::asin(x); }},
                             {"acos", [](double x) { return std::acos(x); }}, {"atan", [](double x) { return std::atan(x); }},
                             {"sqrt", [](double x) { return std::sqrt(x); }}, {"abs", [](double x) { return std::fabs(x); }},
                             {"floor", [](double x) { return std::floor(x); }}, {"ceil", [](double x) { return std::ceil(x); }},
                             {"round", c_round},                              {"exp", [](double x) { return std::exp(x); }},
                             {"ln", [](double x) { return std::log(x); }},    {"log", [](double x) { return std::log10(x); }}};
    for (const F1& f : f1s)
      if (name == f.n) {
        if (args.size() != 1)
          throw ParseError(pos, name + "() takes 1 argument, got " + std::to_string(args.size()));
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::Fn1;
        n->fname = f.n;
        n->f1 = f.f;
        n->a = args[0];
        return n;
      }
    struct F2 {
      const char* n;
      double (*f)(double, double);
    };
    static const F2 f2s[] = {{"atan2", [](double a, double b) { return std::atan2(a, b); }},
                             {"min", [](double a, double b) { return std::fmin(a, b); }},
                             {"max", [](double a, double b) { return std::fmax(a, b); }},
                             {"pow", [](double a, double b) { return std::pow(a, b); }}};
    for (const F2& f : f2s)
      if (name == f.n) {
        if (args.size() != 2)
          throw ParseError(pos, name + "() takes 2 arguments, got " + std::to_string(args.size()));
        auto n = std::make_shared<Node>();
        n->kind = Node::Kind::Fn2;
        n->fname = f.n;
        n->f2 = f.f;
        n->a = args[0];
        n->b = args[1];
        return n;
      }
    throw ParseError(pos, "unknown function '" + name + "'");
  }
};

}  // namespace

Expr Expr::parse(const std::string& source) {
  Parser p{tokenize(source), 0, source.size()};
  NodeP root = p.expr_bp(0);
  if (p.i < p.toks.size()) throw ParseError(p.toks[p.i].pos, "unexpected trailing token " + describe(p.toks[p.i]));
  Expr e;
  e.root_ = std::move(root);
  return e;
}

double Expr::eval(const double vars[6]) const { return root_->eval(vars); }

}  // namespace keta::param::expr

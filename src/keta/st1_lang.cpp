#include "keta/st1_lang.h"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <set>
#include <sstream>

namespace keta::st1 {
namespace {

std::string lower(std::string s) {
  for (char& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
std::string upper(std::string s) {
  for (char& c : s) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

std::string fmt_num(double v) {
  char buf[40];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

// Reserved words that must never be shadowed by a user variable (ST1 keywords
// are reserved, manual p.30). Only the short/ambiguous ones need listing —
// longer command words simply never collide because substitution only fires
// for identifiers DEFINED in the environment, and defining these is refused.
const std::set<std::string>& keyword_blacklist() {
  static const std::set<std::string> k = {
      "a",    "pas",  "de",   "or",   "ex",   "mi",  "x",   "y",    "z",    "e",
      "nu",   "ro",   "sx",   "sy",   "sz",   "ix",  "iy",  "iz",   "vy",   "vz",
      "wy",   "wz",   "n",    "ty",   "tz",   "mx",  "my",  "mz",   "dx",   "dy",
      "dz",   "rx",   "ry",   "rz",   "fx",   "fy",  "fz",  "tout", "toutes",
      "uni",  "lin",  "con",  "glo",  "loc",  "pro", "rel", "abs",  "el",   "di",
      "id",   "se",   "sup",  "inf",  "pos",  "neg", "cte", "var",  "si",   "pour",
      "sinon", "et",  "ou",   "non",  "vide", "dim", "get", "fin",  "max",  "min"};
  return k;
}

}  // namespace

bool Env::has_scalar(const std::string& n) const { return scalars.count(n) > 0; }
double Env::scalar(const std::string& n) const {
  auto it = scalars.find(n);
  require(it != scalars.end(), "ST1: variable inconnue '" + n + "'");
  return it->second;
}

// ---------------------------------------------------------------- numbers --
bool is_st1_number(const std::string& s) {
  if (s.empty()) return false;
  std::size_t i = (s[0] == '+' || s[0] == '-') ? 1 : 0;
  bool digits = false, dot = false, exp = false, exp_digits = false;
  for (; i < s.size(); ++i) {
    char c = s[i];
    if (std::isdigit(static_cast<unsigned char>(c))) {
      (exp ? exp_digits : digits) = true;
    } else if (c == '.' && !dot && !exp) {
      dot = true;
    } else if ((c == 'e' || c == 'E' || c == 'd' || c == 'D') && digits && !exp) {
      exp = true;
      if (i + 1 < s.size() && (s[i + 1] == '+' || s[i + 1] == '-')) ++i;
    } else {
      return false;
    }
  }
  return digits && (!exp || exp_digits);
}

double st1_num(const std::string& s) {
  std::string t = s;
  for (char& c : t)
    if (c == 'd' || c == 'D') c = 'e';
  return std::stod(t);
}

// ------------------------------------------------------------------- lexer --
std::vector<Statement> lex(const std::string& text) {
  // Physical lines -> merged logical lines (trailing-comma continuation).
  std::vector<std::string> phys;
  {
    std::stringstream ss(text);
    std::string line;
    while (std::getline(ss, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      phys.push_back(line);
    }
  }
  // Strip comments (respecting quotes), then merge continuations.
  auto strip_comment = [](const std::string& l) {
    std::string out;
    char q = 0;  // the open quote, ' or " (a string never holds the other)
    for (std::size_t i = 0; i < l.size(); ++i) {
      char c = l[i];
      if (c == '\'' || c == '"') {
        if (!q) q = c;
        else if (q == c) q = 0;
      }
      if (c == '#' && !q) break;
      out += c;
    }
    return out;
  };
  std::vector<std::string> logical;
  std::vector<int> logical_line;  // physical line (1-based) each logical line starts on
  std::string pending;
  int pending_line = 0;
  for (std::size_t k = 0; k < phys.size(); ++k) {
    std::string l = strip_comment(phys[k]);
    std::size_t end = l.find_last_not_of(" \t");
    bool cont = end != std::string::npos && l[end] == ',';
    if (pending.empty()) pending_line = static_cast<int>(k) + 1;
    pending += l;
    if (cont) {
      pending += ' ';
      continue;
    }
    logical.push_back(pending);
    logical_line.push_back(pending_line);
    pending.clear();
  }
  if (!pending.empty()) {
    logical.push_back(pending);
    logical_line.push_back(pending_line);
  }

  // Tokenize each logical line, split on ';' into statements, isolate the
  // << >> delimiters and ',' separators as standalone tokens.
  std::vector<Statement> out;
  for (std::size_t li = 0; li < logical.size(); ++li) {
    const std::string& line = logical[li];
    const int line_no = logical_line[li];
    Statement cur;
    auto push = [&](Statement st) {
      for (Tok& t : st) t.line = line_no;
      out.push_back(std::move(st));
    };
    auto flush = [&]() {
      if (!cur.empty()) push(cur);
      cur.clear();
    };
    std::size_t i = 0;
    while (i < line.size()) {
      char c = line[i];
      if (c == ' ' || c == '\t') {
        ++i;
        continue;
      }
      if (c == ';') {
        flush();
        ++i;
        continue;
      }
      if (c == '\'' || c == '"') {  // quoted string ('…' or "…"), a doubled quote collapses
        const char q = c;
        std::string s;
        ++i;
        while (i < line.size()) {
          if (line[i] == q && i + 1 < line.size() && line[i + 1] == q) {
            s += q;
            i += 2;
          } else if (line[i] == q) {
            ++i;
            break;
          } else {
            s += line[i++];
          }
        }
        cur.push_back({s, true});
        continue;
      }
      if (c == '<' && i + 1 < line.size() && line[i + 1] == '<') {
        flush();
        push({{"<<", false}});
        i += 2;
        continue;
      }
      if (c == '>' && i + 1 < line.size() && line[i + 1] == '>') {
        flush();
        push({{">>", false}});
        i += 2;
        continue;
      }
      if (c == ',') {
        cur.push_back({",", false});
        ++i;
        continue;
      }
      // Plain token: run to the next blank / ';' / ',' / quote / << / >>.
      std::string s;
      while (i < line.size()) {
        char d = line[i];
        if (d == ' ' || d == '\t' || d == ';' || d == ',' || d == '\'' || d == '"') break;
        if ((d == '<' && i + 1 < line.size() && line[i + 1] == '<') ||
            (d == '>' && i + 1 < line.size() && line[i + 1] == '>'))
          break;
        s += d;
        ++i;
      }
      cur.push_back({s, false});
    }
    flush();
  }
  return out;
}

// -------------------------------------------------------------- expression --
namespace expr {

struct Lexer {
  const std::string& s;
  const Env& env;
  std::size_t i = 0;
  void skip() {
    while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) ++i;
  }
  bool eat(const char* op) {
    skip();
    std::size_t n = std::strlen(op);
    if (s.compare(i, n, op) == 0) {
      i += n;
      return true;
    }
    return false;
  }
  bool peek_char(char c) {
    skip();
    return i < s.size() && s[i] == c;
  }
  bool done() {
    skip();
    return i >= s.size();
  }
};

double parse_ou(Lexer& lx);

double call_fn(const std::string& name, double x) {
  if (name == "abs") return std::fabs(x);
  if (name == "int") return std::trunc(x);  // truncation toward zero (documented)
  if (name == "sqrt") {
    require(x >= 0, "ST1: sqrt d'un negatif");
    return std::sqrt(x);
  }
  if (name == "sin") return std::sin(x);
  if (name == "cos") return std::cos(x);
  if (name == "tan") return std::tan(x);
  if (name == "asin") return std::asin(x);
  if (name == "acos") return std::acos(x);
  if (name == "atan") return std::atan(x);
  if (name == "exp") return std::exp(x);
  if (name == "log") {
    require(x > 0, "ST1: log d'un non-positif");
    return std::log(x);
  }
  if (name == "sinh") return std::sinh(x);
  if (name == "cosh") return std::cosh(x);
  if (name == "tanh") return std::tanh(x);
  fail("ST1: fonction inconnue '" + name + "'");
}

bool is_fn(const std::string& name) {
  static const std::set<std::string> f = {"abs",  "int",  "sqrt", "sin",  "cos",
                                          "tan",  "asin", "acos", "atan", "exp",
                                          "log",  "sinh", "cosh", "tanh"};
  return f.count(name) > 0;
}

double parse_primary(Lexer& lx) {
  lx.skip();
  require(lx.i < lx.s.size(), "ST1: expression tronquee '" + lx.s + "'");
  char c = lx.s[lx.i];
  if (c == '(') {
    ++lx.i;
    double v = parse_ou(lx);
    require(lx.eat(")"), "ST1: parenthese fermante manquante dans '" + lx.s + "'");
    return v;
  }
  if (c == '-') {
    ++lx.i;
    return -parse_primary(lx);
  }
  if (c == '+') {
    ++lx.i;
    return parse_primary(lx);
  }
  if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
    std::size_t j = lx.i;
    bool exp = false;
    while (j < lx.s.size()) {
      char d = lx.s[j];
      if (std::isdigit(static_cast<unsigned char>(d)) || d == '.') {
        ++j;
      } else if ((d == 'e' || d == 'E' || d == 'd' || d == 'D') && !exp &&
                 j + 1 < lx.s.size() &&
                 (std::isdigit(static_cast<unsigned char>(lx.s[j + 1])) ||
                  lx.s[j + 1] == '+' || lx.s[j + 1] == '-')) {
        exp = true;
        j += (lx.s[j + 1] == '+' || lx.s[j + 1] == '-') ? 2 : 1;
      } else {
        break;
      }
    }
    std::string num = lx.s.substr(lx.i, j - lx.i);
    lx.i = j;
    require(is_st1_number(num), "ST1: nombre invalide '" + num + "'");
    return st1_num(num);
  }
  if (std::isalpha(static_cast<unsigned char>(c)) || c == '_' || c == '$') {
    std::size_t j = lx.i;
    while (j < lx.s.size() &&
           (std::isalnum(static_cast<unsigned char>(lx.s[j])) || lx.s[j] == '_' ||
            lx.s[j] == '$'))
      ++j;
    std::string name = lower(lx.s.substr(lx.i, j - lx.i));
    lx.i = j;
    if (name == "pi") return std::acos(-1.0);
    if (name == "non") return parse_primary(lx) != 0.0 ? 0.0 : 1.0;
    if (lx.peek_char('(')) {
      ++lx.i;  // '('
      double arg = parse_ou(lx);
      require(lx.eat(")"), "ST1: parenthese fermante manquante dans '" + lx.s + "'");
      if (is_fn(name)) return call_fn(name, arg);
      auto it = lx.env.arrays.find(name);
      require(it != lx.env.arrays.end(), "ST1: tableau ou fonction inconnu '" + name + "'");
      long idx = std::lround(arg);
      require(idx >= 1 && idx <= static_cast<long>(it->second.size()),
              "ST1: indice hors bornes pour '" + name + "'");
      return it->second[static_cast<std::size_t>(idx - 1)];
    }
    if (lx.env.has_scalar(name)) return lx.env.scalar(name);
    auto lt = lx.env.lists.find(name);
    if (lt != lx.env.lists.end() && lt->second.size() == 1) return lt->second[0];
    fail("ST1: variable inconnue '" + name + "' dans '" + lx.s + "'");
  }
  fail("ST1: expression invalide '" + lx.s + "'");
}

double parse_pow(Lexer& lx) {
  double base = parse_primary(lx);
  if (lx.eat("**")) {
    double e = parse_pow(lx);  // right-associative
    return std::pow(base, e);
  }
  return base;
}

double parse_mul(Lexer& lx) {
  double v = parse_pow(lx);
  for (;;) {
    lx.skip();
    // '/=' is the not-equal comparator, not a division.
    if (lx.peek_char('/') && !(lx.i + 1 < lx.s.size() && lx.s[lx.i + 1] == '=')) {
      ++lx.i;
      double d = parse_pow(lx);
      require(d != 0.0, "ST1: division par zero dans '" + lx.s + "'");
      v /= d;
    } else if (lx.peek_char('*') && !(lx.i + 1 < lx.s.size() && lx.s[lx.i + 1] == '*')) {
      ++lx.i;
      v *= parse_pow(lx);
    } else {
      return v;
    }
  }
}

double parse_add(Lexer& lx) {
  double v = parse_mul(lx);
  for (;;) {
    lx.skip();
    if (lx.peek_char('+')) {
      ++lx.i;
      v += parse_mul(lx);
    } else if (lx.peek_char('-')) {
      ++lx.i;
      v -= parse_mul(lx);
    } else {
      return v;
    }
  }
}

double parse_cmp(Lexer& lx) {
  double v = parse_add(lx);
  lx.skip();
  if (lx.eat("==") || (lx.peek_char('=') && (++lx.i, true)))
    return v == parse_add(lx) ? 1.0 : 0.0;
  if (lx.eat("/=")) return v != parse_add(lx) ? 1.0 : 0.0;
  if (lx.eat(">=")) return v >= parse_add(lx) ? 1.0 : 0.0;
  if (lx.eat("<=")) return v <= parse_add(lx) ? 1.0 : 0.0;
  if (lx.eat(">")) return v > parse_add(lx) ? 1.0 : 0.0;
  if (lx.eat("<")) return v < parse_add(lx) ? 1.0 : 0.0;
  return v;
}

// Word operator scan: 'et' / 'ou' as standalone identifiers.
bool eat_word(Lexer& lx, const char* w) {
  lx.skip();
  std::size_t n = std::strlen(w), i = lx.i;
  if (i + n > lx.s.size()) return false;
  for (std::size_t k = 0; k < n; ++k)
    if (std::tolower(static_cast<unsigned char>(lx.s[i + k])) != w[k]) return false;
  // must not be part of a longer identifier
  if (i + n < lx.s.size() &&
      (std::isalnum(static_cast<unsigned char>(lx.s[i + n])) || lx.s[i + n] == '_'))
    return false;
  lx.i = i + n;
  return true;
}

double parse_et(Lexer& lx) {
  double v = parse_cmp(lx);
  while (eat_word(lx, "et")) {
    double r = parse_cmp(lx);
    v = (v != 0.0 && r != 0.0) ? 1.0 : 0.0;
  }
  return v;
}

double parse_ou(Lexer& lx) {
  double v = parse_et(lx);
  while (eat_word(lx, "ou")) {
    double r = parse_et(lx);
    v = (v != 0.0 || r != 0.0) ? 1.0 : 0.0;
  }
  return v;
}

}  // namespace expr

double eval_expr(const std::string& text, const Env& env) {
  expr::Lexer lx{text, env};
  // Leading 'non'
  double v = expr::parse_ou(lx);
  require(lx.done(), "ST1: fin d'expression inattendue dans '" + text + "'");
  return v;
}

// ------------------------------------------------------------- interpreter --
struct Interp::Node {
  enum class Kind { Stmt, For, If } kind = Kind::Stmt;
  Statement stmt;                 // Stmt
  std::string var;                // For: loop variable (lower)
  Statement spec;                 // For: list spec tokens
  bool faire = false;             // For: FAIRE bounds semantics (start,end[,step])
  std::string test;               // If: expression text
  std::vector<Node> body, els;
  int line = 0;                   // line of the statement that opened the node
};

Interp::Interp(Env& env, Sink sink) : env_(env), sink_(std::move(sink)) {
  loader_ = [](const std::string& path) {
    std::ifstream f(path);
    require(f.good(), "ST1: LIRE — fichier introuvable '" + path + "'");
    std::stringstream ss;
    ss << f.rdbuf();
    return ss.str();
  };
}

namespace {
bool is_delim(const Statement& st, const char* d) {
  return st.size() == 1 && !st[0].str && st[0].s == d;
}
}  // namespace

// Parse a run of statements into a node tree.
std::vector<Interp::Node> Interp::parse(const std::vector<Statement>& sts, std::size_t& i,
                                        Until until) const {
  std::vector<Node> out;
  while (i < sts.size()) {
    if (until == Until::One && !out.empty()) return out;
    const Statement& st = sts[i];
    parse_line_ = st.empty() ? parse_line_ : st[0].line;
    if (is_delim(st, ">>")) {
      require(until == Until::Block, "ST1: '>>' sans '<<'");
      ++i;
      return out;
    }
    std::string head = st.empty() || st[0].str ? "" : upper(st[0].s);
    if (head == "POUR" || head == "FAIRE") {
      // POUR var = liste   |   FAIRE var = debut,fin[,pas]  (Fortran DO style,
      // attested in real decks: FAIRE i=1,ncoug).
      require(st.size() >= 2 && !st[1].str, "ST1: POUR/FAIRE invalide");
      Node n;
      n.kind = Node::Kind::For;
      n.line = st.empty() ? 0 : st[0].line;
      n.faire = head == "FAIRE";
      // The variable and '=' may be glued (no-blank rule: 'pour i=1 a n').
      std::size_t k = 2;
      std::string first = st[1].s;
      auto eq = first.find('=');
      if (eq != std::string::npos) {
        n.var = lower(first.substr(0, eq));
        if (eq + 1 < first.size()) {
          n.spec.push_back({first.substr(eq + 1), false});
        }
      } else {
        n.var = lower(first);
        require(k < st.size() && st[k].s == "=", "ST1: POUR attend 'var = liste'");
        ++k;
      }
      require(!n.var.empty(), "ST1: POUR/FAIRE invalide");
      require(!keyword_blacklist().count(n.var),
              "ST1: nom de variable reserve '" + n.var + "'");
      for (; k < st.size(); ++k) n.spec.push_back(st[k]);
      require(!n.spec.empty(), "ST1: POUR attend 'var = liste'");
      ++i;
      // Body: '<<' block or the single next node.
      if (i < sts.size() && is_delim(sts[i], "<<")) {
        ++i;
        n.body = parse(sts, i, Until::Block);
      } else {
        require(i < sts.size(), "ST1: POUR sans corps");
        n.body = parse(sts, i, Until::One);
      }
      out.push_back(std::move(n));
      continue;
    }
    if (head == "SI" || head.rfind("SI(", 0) == 0) {
      // Test = concatenation of tokens until parens balance. The no-blank
      // rule may glue the head: 'SI( kxal > 1e-7 )' tokenizes as 'SI(' ...
      Node n;
      n.kind = Node::Kind::If;
      n.line = st.empty() ? 0 : st[0].line;
      Statement expanded;
      if (head != "SI") {
        expanded.push_back({"SI", false});
        expanded.push_back({st[0].s.substr(2), false});
        for (std::size_t k2 = 1; k2 < st.size(); ++k2) expanded.push_back(st[k2]);
      }
      const Statement& sv = head == "SI" ? st : expanded;
      std::size_t k = 1;
      int depth = 0;
      std::string test;
      for (; k < sv.size(); ++k) {
        const std::string& t = sv[k].s;
        for (char c : t) {
          if (c == '(') ++depth;
          if (c == ')') --depth;
        }
        if (!test.empty()) test += ' ';
        test += t;
        if (depth == 0) {
          ++k;
          break;
        }
      }
      require(depth == 0 && !test.empty(), "ST1: test SI mal parenthese");
      n.test = test;
      if (k < sv.size()) {
        // Inline instruction on the SI line.
        Node inner;
        inner.kind = Node::Kind::Stmt;
        for (; k < sv.size(); ++k) inner.stmt.push_back(sv[k]);
        n.body.push_back(std::move(inner));
        ++i;
      } else {
        ++i;
        if (i < sts.size() && is_delim(sts[i], "<<")) {
          ++i;
          n.body = parse(sts, i, Until::Block);
        } else {
          require(i < sts.size(), "ST1: SI sans corps");
          n.body = parse(sts, i, Until::One);
        }
      }
      // Optional SINON
      if (i < sts.size() && !sts[i].empty() && !sts[i][0].str &&
          upper(sts[i][0].s) == "SINON") {
        const Statement& sn = sts[i];
        if (sn.size() > 1) {
          Node inner;
          inner.kind = Node::Kind::Stmt;
          for (std::size_t k2 = 1; k2 < sn.size(); ++k2) inner.stmt.push_back(sn[k2]);
          n.els.push_back(std::move(inner));
          ++i;
        } else {
          ++i;
          if (i < sts.size() && is_delim(sts[i], "<<")) {
            ++i;
            n.els = parse(sts, i, Until::Block);
          } else {
            require(i < sts.size(), "ST1: SINON sans corps");
            n.els = parse(sts, i, Until::One);
          }
        }
      }
      out.push_back(std::move(n));
      continue;
    }
    Node n;
    n.kind = Node::Kind::Stmt;
    n.line = st.empty() ? 0 : st[0].line;
    n.stmt = st;
    out.push_back(std::move(n));
    ++i;
  }
  require(until != Until::Block, "ST1: '<<' sans '>>'");
  return out;
}

void Interp::fail_here() {
  if (error_.empty()) error_.push_back(where());
}

Interp::Flow Interp::exec(const std::vector<Node>& nodes) {
  for (const Node& n : nodes) {
    if (!stack_.empty()) stack_.back().line = n.line;
    try {
      const Flow f = exec_node(n);
      if (f != Flow::Normal) return f;
    } catch (...) {
      fail_here();
      throw;
    }
  }
  return Flow::Normal;
}

Interp::Flow Interp::exec_node(const Node& n) {
  {
    switch (n.kind) {
      case Node::Kind::For: {
        // Evaluate the list spec: numbers/exprs with ',' separators,
        // 'a'-ranges with optional 'pas', or a named list. FAIRE bounds are
        // Fortran DO semantics: start,end[,step].
        std::vector<double> values;
        const Statement sub = substitute(n.spec);
        std::size_t k = 0;
        auto val = [&](const Tok& t) { return eval_expr(t.s, env_); };
        if (n.faire) {
          std::vector<double> b;
          for (const Tok& t : sub)
            if (t.s != ",") b.push_back(val(t));
          require(b.size() >= 2 && b.size() <= 3, "ST1: FAIRE attend debut,fin[,pas]");
          double step = b.size() == 3 ? b[2] : 1.0;
          require(step != 0, "ST1: FAIRE — pas nul");
          if (step > 0)
            for (double v = b[0]; v <= b[1] + 1e-12; v += step) values.push_back(v);
          else
            for (double v = b[0]; v >= b[1] - 1e-12; v += step) values.push_back(v);
          bool broke2 = false;
          for (double v : values) {
            env_.scalars[n.var] = v;
            Flow f = exec(n.body);
            if (f == Flow::Break) {
              broke2 = true;
              break;
            }
            if (f == Flow::Return) return Flow::Return;
          }
          (void)broke2;
          break;
        }
        while (k < sub.size()) {
          if (sub[k].s == ",") {
            ++k;
            continue;
          }
          std::string low = lower(sub[k].s);
          if (!sub[k].str && env_.lists.count(low)) {
            for (double v : env_.lists.at(low)) values.push_back(v);
            ++k;
            continue;
          }
          double start = val(sub[k]);
          ++k;
          if (k < sub.size() && lower(sub[k].s) == "a") {
            ++k;
            require(k < sub.size(), "ST1: plage 'a' tronquee");
            double end = val(sub[k]);
            ++k;
            double step = 1;
            if (k + 1 < sub.size() + 1 && k < sub.size() && lower(sub[k].s) == "pas") {
              ++k;
              require(k < sub.size(), "ST1: PAS tronque");
              step = val(sub[k]);
              ++k;
            }
            require(step != 0, "ST1: PAS nul");
            if (step > 0)
              for (double v = start; v <= end + 1e-12; v += step) values.push_back(v);
            else
              for (double v = start; v >= end - 1e-12; v += step) values.push_back(v);
          } else {
            values.push_back(start);
          }
        }
        bool broke = false;
        for (double v : values) {
          env_.scalars[n.var] = v;
          Flow f = exec(n.body);
          if (f == Flow::Break) {
            broke = true;
            break;
          }
          if (f == Flow::Return) return Flow::Return;
          // Continue falls through to the next iteration.
        }
        (void)broke;  // loop variable keeps its last value (manual p.259)
        break;
      }
      case Node::Kind::If: {
        // Substitute variables inside the test text before evaluation is
        // unnecessary: eval_expr resolves variables itself.
        double t = eval_expr(n.test, env_);
        Flow f = (t != 0.0) ? exec(n.body) : exec(n.els);
        if (f != Flow::Normal) return f;
        break;
      }
      case Node::Kind::Stmt: {
        Flow f = exec_stmt(n.stmt);
        if (f != Flow::Normal) return f;
        break;
      }
    }
  }
  return Flow::Normal;
}

Interp::Flow Interp::exec_stmt(const Statement& st) {
  if (st.empty()) return Flow::Normal;
  std::string head = st[0].str ? "" : upper(st[0].s);

  if (head == "PASSER") return Flow::Continue;
  if (head == "BRISER") return Flow::Break;
  if (head == "RETOUR") return Flow::Return;
  if (head == "PAUSE") return Flow::Normal;  // interactive pause: no-op headless

  if (head == "DIM") {
    // DIM nom(n)
    require(st.size() >= 2 && !st[1].str, "ST1: DIM invalide");
    const std::string& t = st[1].s;
    auto p = t.find('(');
    require(p != std::string::npos && t.back() == ')', "ST1: DIM attend nom(n)");
    std::string name = lower(t.substr(0, p));
    require(!keyword_blacklist().count(name), "ST1: nom reserve '" + name + "'");
    double n = eval_expr(t.substr(p + 1, t.size() - p - 2), env_);
    long sz = std::lround(n);
    require(sz >= 1, "ST1: dimension de tableau invalide");
    env_.arrays[name].assign(static_cast<std::size_t>(sz), 0.0);
    return Flow::Normal;
  }

  if (head == "ECRIRE" || head == "EC") {  // 'ec' abbreviation attested in decks
    std::string outl;
    const Statement sub = substitute(st);
    for (std::size_t k = 1; k < sub.size(); ++k) {
      if (!outl.empty()) outl += ' ';
      outl += sub[k].s;
    }
    if (echo_) echo_(outl);
    return Flow::Normal;
  }

  if (head == "LIRE") {
    require(st.size() >= 2, "ST1: LIRE invalide");
    if (!st[1].str && upper(st[1].s) == "CONSOLE") return Flow::Normal;  // headless
    require(st[1].str, "ST1: LIRE attend un nom de fichier entre quotes");
    require(lire_depth_ < 8, "ST1: LIRE imbrique trop profondement");
    std::string text = loader_(st[1].s);
    ++lire_depth_;
    Flow f;
    try {
      f = run_nodes(text, st[1].s);
    } catch (...) {
      --lire_depth_;
      throw;
    }
    --lire_depth_;
    // RETOUR inside a file stops that file only.
    (void)f;
    return Flow::Normal;
  }

  // Assignment forms: `name = expr`, `name = list`, `name(i) = expr`,
  // or a single token `name=expr` (no-blank rule).
  {
    std::string lhs;
    Statement rhs;
    if (st.size() >= 3 && !st[0].str && !st[1].str && st[1].s == "=") {
      lhs = st[0].s;
      for (std::size_t k = 2; k < st.size(); ++k) rhs.push_back(st[k]);
    } else if (st.size() >= 2 && !st[0].str && !st[1].str && st[1].s.size() >= 2 &&
               st[1].s[0] == '=' && st[1].s[1] != '=' &&
               (std::isalpha(static_cast<unsigned char>(st[0].s[0])) || st[0].s[0] == '_')) {
      // 'lg =0' — the '=' glued to the VALUE side (no-blank rule cuts both ways).
      bool ident = true;
      for (char c : st[0].s)
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '(' ||
              c == ')'))
          ident = false;
      if (ident) {
        lhs = st[0].s;
        rhs.push_back({st[1].s.substr(1), false});
        for (std::size_t k = 2; k < st.size(); ++k) rhs.push_back(st[k]);
      }
    }
    if (lhs.empty() && st.size() >= 1 && !st[0].str) {
      const std::string& t = st[0].s;
      auto p = t.find('=');
      if (p != std::string::npos && p > 0 && t.find("==") == std::string::npos &&
          t[p - 1] != '<' && t[p - 1] != '>' && t[p - 1] != '/') {
        // token like a=1+2 — only treat as assignment when the lhs is a
        // plain or indexed identifier.
        std::string cand = t.substr(0, p);
        bool ident = std::isalpha(static_cast<unsigned char>(cand[0])) || cand[0] == '_';
        for (char c : cand)
          if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '(' ||
                c == ')' || c == '$'))
            ident = false;
        if (ident) {
          lhs = cand;
          // 'lg0= 10/nb': the '=' may end the token — never push an empty
          // expression token.
          if (p + 1 < t.size()) rhs.push_back({t.substr(p + 1), false});
          for (std::size_t k = 1; k < st.size(); ++k) rhs.push_back(st[k]);
        }
      }
    }
    if (!lhs.empty()) {
      do_assign(lhs, rhs);
      return Flow::Normal;
    }
  }

  // Engine statement: substitute and hand to the sink.
  sink_(substitute(st));
  return Flow::Normal;
}

void Interp::do_assign(const std::string& lhs_raw, const Statement& rhs) {
  std::string lhs = lower(lhs_raw);
  auto p = lhs.find('(');
  if (p != std::string::npos) {
    // Indexed array element.
    require(lhs.back() == ')', "ST1: affectation indicee invalide");
    std::string name = lhs.substr(0, p);
    auto it = env_.arrays.find(name);
    require(it != env_.arrays.end(), "ST1: tableau non declare '" + name + "' (DIM)");
    long idx = std::lround(eval_expr(lhs.substr(p + 1, lhs.size() - p - 2), env_));
    require(idx >= 1 && idx <= static_cast<long>(it->second.size()),
            "ST1: indice hors bornes pour '" + name + "'");
    require(rhs.size() == 1 && !rhs[0].str, "ST1: affectation indicee attend une valeur");
    it->second[static_cast<std::size_t>(idx - 1)] = eval_expr(rhs[0].s, env_);
    return;
  }
  require(!keyword_blacklist().count(lhs), "ST1: nom de variable reserve '" + lhs + "'");
  // List literal? (contains ',' tokens, an 'a' range, or the word 'vide')
  bool is_list = false;
  for (const Tok& t : rhs) {
    std::string low = lower(t.s);
    if (t.s == "," || low == "a" || low == "vide" || env_.lists.count(low)) is_list = true;
  }
  if (rhs.size() > 1) {
    // Multiple bare values also form a list.
    bool all_val = true;
    for (const Tok& t : rhs) {
      if (t.s == "," || lower(t.s) == "a" || lower(t.s) == "pas") continue;
      if (t.str) all_val = false;
    }
    if (all_val) is_list = true;
  }
  if (is_list) {
    std::vector<double> vals;
    std::size_t k = 0;
    const Statement sub = substitute(rhs);
    while (k < sub.size()) {
      if (sub[k].s == ",") {
        ++k;
        continue;
      }
      if (lower(sub[k].s) == "vide") {
        ++k;
        continue;
      }
      std::string low = lower(sub[k].s);
      if (!sub[k].str && env_.lists.count(low)) {
        for (double v : env_.lists.at(low)) vals.push_back(v);
        ++k;
        continue;
      }
      double start = eval_expr(sub[k].s, env_);
      ++k;
      if (k < sub.size() && lower(sub[k].s) == "a") {
        ++k;
        require(k < sub.size(), "ST1: plage 'a' tronquee");
        double end = eval_expr(sub[k].s, env_);
        ++k;
        double step = 1;
        if (k < sub.size() && lower(sub[k].s) == "pas") {
          ++k;
          require(k < sub.size(), "ST1: PAS tronque");
          step = eval_expr(sub[k].s, env_);
          ++k;
        }
        require(step != 0, "ST1: PAS nul");
        if (step > 0)
          for (double v = start; v <= end + 1e-12; v += step) vals.push_back(v);
        else
          for (double v = start; v >= end - 1e-12; v += step) vals.push_back(v);
      } else {
        vals.push_back(start);
      }
    }
    env_.lists[lhs] = std::move(vals);
    env_.scalars.erase(lhs);
    return;
  }
  require(rhs.size() == 1 && !rhs[0].str, "ST1: affectation invalide pour '" + lhs + "'");
  env_.scalars[lhs] = eval_expr(rhs[0].s, env_);
  env_.lists.erase(lhs);
}

Statement Interp::substitute(const Statement& st) const {
  Statement out;
  for (const Tok& t : st) {
    if (t.str || t.s == "," || t.s == "<<" || t.s == ">>") {
      out.push_back(t);
      continue;
    }
    if (is_st1_number(t.s)) {
      // Normalize Fortran d-exponents so downstream parsers see e-notation.
      std::string s = t.s;
      bool d = false;
      for (char& c : s)
        if (c == 'd' || c == 'D') {
          c = 'e';
          d = true;
        }
      out.push_back({d ? s : t.s, false});
      continue;
    }
    std::string low = lower(t.s);
    bool pure_ident = std::isalpha(static_cast<unsigned char>(t.s[0])) || t.s[0] == '_' ||
                      t.s[0] == '$';
    for (char c : t.s)
      if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '$'))
        pure_ident = false;
    if (pure_ident) {
      if (keyword_blacklist().count(low) && t.s[0] != '$') {
        out.push_back(t);
        continue;
      }
      if (env_.has_scalar(low)) {
        out.push_back({fmt_num(env_.scalar(low)), false});
        continue;
      }
      auto lt = env_.lists.find(low);
      if (lt != env_.lists.end()) {
        // A list variable stands for its comma-separated list (`civil =
        // civil,1` then `SURCH civil 1.2` in Exemple1): without the commas
        // a trailing element would read as a coefficient.
        for (std::size_t k = 0; k < lt->second.size(); ++k) {
          if (k > 0) out.push_back({",", false});
          out.push_back({fmt_num(lt->second[k]), false});
        }
        continue;
      }
      out.push_back(t);  // unknown identifier: engine keyword
      continue;
    }
    // Token with expression characters: try to evaluate; keep verbatim when
    // it is not an expression (e.g. W800_600, file paths).
    bool has_expr_char = false;
    for (char c : t.s)
      if (c == '+' || c == '-' || c == '*' || c == '/' || c == '(' || c == '$') has_expr_char = true;
    if (has_expr_char) {
      try {
        double v = eval_expr(t.s, env_);
        out.push_back({fmt_num(v), false});
        continue;
      } catch (const Error&) {
        // fall through: keyword-like token
      }
    }
    out.push_back(t);
  }
  return out;
}

Interp::Flow Interp::run_nodes(const std::string& text, const std::string& file) {
  stack_.push_back({file, 0});
  try {
    std::vector<Statement> sts = lex(text);
    std::size_t i = 0;
    parse_line_ = 0;
    std::vector<Node> nodes;
    try {
      nodes = parse(sts, i, Until::End);
    } catch (...) {
      stack_.back().line = parse_line_;
      fail_here();
      throw;
    }
    const Flow f = exec(nodes);
    stack_.pop_back();
    return f;
  } catch (...) {
    stack_.pop_back();
    throw;
  }
}

void Interp::run(const std::string& text, const std::string& file) {
  error_.clear();
  run_nodes(text, file);
}

}  // namespace keta::st1

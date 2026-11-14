// ST1 pseudo-programming interpreter (manual p.30-31, p.259-277).
//
// The ST1 language layers a small scripting facility over the engine
// commands: scalar variables and expressions (`nom = 1+2`, operators
// + - * / ** with ** right-associative, functions abs int sqrt sin cos tan
// asin acos atan exp log sinh cosh tanh, constant pi = acos(-1)); named
// numeric lists (`ls = 1,2,5 a 8`, `vide` = empty); DIM arrays with 1-based
// indexing `nom(i)`; POUR loops over lists with PASSER / BRISER; SI/SINON
// tests with comparators == (or =) /= > < >= <= and word operators
// NON ET OU; ECRIRE output; PAUSE; LIRE file inclusion with RETOUR;
// GET result-buffer variables prefixed `$` (case-insensitive, set by the
// runtime through the shared environment).
//
// Lexical rules implemented here (p.30-31):
//   * `#` comments to end of line (outside quoted strings);
//   * `;` separates statements on one line;
//   * a line whose last non-blank character is `,` continues on the next;
//   * quoted strings with doubled inner quotes (`'Poussage d''un pont'`);
//   * Fortran-style `d` exponents (.2d-3) are accepted anywhere a number is;
//   * `<<` / `>>` block delimiters are split into standalone statements;
//   * `,` inside token runs separates list items and is kept as a token.
//
// The interpreter walks the statement tree and hands every ENGINE statement
// (anything that is not pure control flow) to a sink callback, with every
// expression token already substituted by its numeric value and every
// named-list variable expanded in place. The sink executes engine commands
// immediately, so a GET performed by the runtime is visible to the very
// next SI/POUR condition (co-routine style, not a pre-expansion).
#pragma once
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "keta/core.h"

namespace keta::st1 {

// One lexed token. `str` marks a quoted string (quotes stripped, inner
// doubled quotes collapsed) so 'TOUT' the string never collides with TOUT
// the keyword.
struct Tok {
  std::string s;
  bool str = false;
  int line = 0;  // 1-based physical line of the statement in its file (0: unknown)
};
using Statement = std::vector<Tok>;

// Shared variable environment. Keys are lower-cased; `$` prefixes are kept
// (GET buffer variables live here too, written by the runtime).
struct Env {
  std::map<std::string, double> scalars;
  std::map<std::string, std::vector<double>> lists;
  std::map<std::string, std::vector<double>> arrays;  // DIM name(n), 1-based

  bool has_scalar(const std::string& lower_name) const;
  double scalar(const std::string& lower_name) const;  // fails loud if absent
};

// Expression evaluation over a self-contained expression string (its own
// character-level lexer: identifiers, numbers with e/d exponents, operators,
// parentheses, single-argument functions, ET/OU/NON word operators,
// comparators). Comparisons and logic yield 1.0 / 0.0.
double eval_expr(const std::string& text, const Env& env);

// Numeric literal check/parse shared with the command layer: accepts
// 1e6, .2d-3, -4.5, +3. Returns false for anything else.
bool is_st1_number(const std::string& s);
double st1_num(const std::string& s);

// Lex a full text into statements (comments, continuations, `;`, quotes,
// `<<`/`>>` isolation, comma tokens).
std::vector<Statement> lex(const std::string& text);

class Interp {
 public:
  using Sink = std::function<void(const Statement&)>;
  using Echo = std::function<void(const std::string&)>;   // ECRIRE lines
  using FileLoader = std::function<std::string(const std::string&)>;

  // A place in the deck: the file (as named to run() or to LIRE) and the
  // 1-based line of the statement (0 when unknown).
  struct Where {
    std::string file;
    int line = 0;
  };

  Interp(Env& env, Sink sink);
  void set_echo(Echo e) { echo_ = std::move(e); }
  void set_file_loader(FileLoader f) { loader_ = std::move(f); }

  // Interpret a text: control flow handled here, engine statements
  // substituted and passed to the sink in execution order. `file` names the
  // text in where() and error_where() (erreur.txt, progress reports).
  void run(const std::string& text, const std::string& file = "");

  // The statement being executed (LIRE files included), for the sink.
  Where where() const { return stack_.empty() ? Where{} : stack_.back(); }
  // Where the last run() failed: the innermost statement that threw, set
  // before the exception leaves run(). Empty after a run that succeeded.
  const std::vector<Where>& error_where() const { return error_; }

 private:
  struct Node;
  // Parse statements into a node tree. Until::End consumes everything,
  // Until::Block stops at '>>', Until::One parses exactly one node (the
  // body of a block-less POUR/SI).
  enum class Until { End, Block, One };
  std::vector<Node> parse(const std::vector<Statement>& sts, std::size_t& i,
                          Until until) const;
  // Execution status threaded through nested blocks.
  enum class Flow { Normal, Break, Continue, Return };
  Flow exec(const std::vector<Node>& nodes);
  Flow exec_node(const Node& n);
  // Lex, parse and execute one text (run() or a LIRE file) under `file`.
  Flow run_nodes(const std::string& text, const std::string& file);
  Flow exec_stmt(const Statement& st);
  void do_assign(const std::string& lhs_raw, const Statement& rhs);
  Statement substitute(const Statement& st) const;

  Env& env_;
  Sink sink_;
  Echo echo_;
  FileLoader loader_;
  int lire_depth_ = 0;
  std::vector<Where> stack_;   // run() text, then the open LIRE files
  std::vector<Where> error_;   // 0 or 1 entry: the failing statement
  mutable int parse_line_ = 0; // statement parse() is on (its errors)
  void fail_here();            // record where() in error_ once (innermost wins)
};

}  // namespace keta::st1

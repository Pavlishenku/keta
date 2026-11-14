#include "keta/json.h"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "keta/core.h"

namespace keta::json {

bool Value::as_bool() const {
  if (type_ != Type::Bool) fail("json: not a bool");
  return bool_;
}
double Value::as_number() const {
  if (type_ != Type::Number) fail("json: not a number");
  return num_;
}
std::int64_t Value::as_int() const {
  double n = as_number();
  double r = std::nearbyint(n);
  if (std::fabs(n - r) > 1e-9) fail("json: not an integer");
  return static_cast<std::int64_t>(r);
}
const std::string& Value::as_string() const {
  if (type_ != Type::String) fail("json: not a string");
  return str_;
}
const Array& Value::as_array() const {
  if (type_ != Type::Array) fail("json: not an array");
  return arr_;
}
Array& Value::as_array() {
  if (type_ != Type::Array) fail("json: not an array");
  return arr_;
}
const Members& Value::members() const {
  if (type_ != Type::Object) fail("json: not an object");
  return obj_;
}
const Value* Value::get(const std::string& key) const {
  if (type_ != Type::Object) return nullptr;
  for (const auto& [k, v] : obj_)
    if (k == key) return &v;
  return nullptr;
}
const Value& Value::at(const std::string& key) const {
  const Value* v = get(key);
  if (!v) fail("json: missing key '" + key + "'");
  return *v;
}
Value& Value::set(const std::string& key, Value v) {
  if (type_ != Type::Object) fail("json: not an object");
  for (auto& [k, existing] : obj_)
    if (k == key) {
      existing = std::move(v);
      return existing;
    }
  obj_.emplace_back(key, std::move(v));
  return obj_.back().second;
}
double Value::num_or(const std::string& key, double dflt) const {
  const Value* v = get(key);
  return v && v->is_number() ? v->as_number() : dflt;
}
bool Value::bool_or(const std::string& key, bool dflt) const {
  const Value* v = get(key);
  return v && v->is_bool() ? v->as_bool() : dflt;
}
std::string Value::str_or(const std::string& key, const std::string& dflt) const {
  const Value* v = get(key);
  return v && v->is_string() ? v->as_string() : dflt;
}
void Value::push_back(Value v) {
  if (type_ != Type::Array) fail("json: not an array");
  arr_.push_back(std::move(v));
}

// ---------------- serializer ----------------

namespace {

void dump_string(const std::string& s, std::string& out) {
  out += '"';
  for (unsigned char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (c < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof buf, "\\u%04x", c);
          out += buf;
        } else {
          out += static_cast<char>(c);
        }
    }
  }
  out += '"';
}

void dump_number(double n, std::string& out) {
  if (!std::isfinite(n)) fail("json: cannot serialize non-finite number");
  double r = std::nearbyint(n);
  char buf[32];
  if (n == r && std::fabs(n) < 1e15) {
    std::snprintf(buf, sizeof buf, "%.0f", n);
  } else {
    std::snprintf(buf, sizeof buf, "%.17g", n);
  }
  out += buf;
}

void dump_value(const Value& v, std::string& out, int indent, int depth) {
  auto newline = [&](int d) {
    if (indent > 0) {
      out += '\n';
      out.append(static_cast<std::size_t>(indent) * d, ' ');
    }
  };
  switch (v.type()) {
    case Value::Type::Null: out += "null"; break;
    case Value::Type::Bool: out += v.as_bool() ? "true" : "false"; break;
    case Value::Type::Number: dump_number(v.as_number(), out); break;
    case Value::Type::String: dump_string(v.as_string(), out); break;
    case Value::Type::Array: {
      const Array& a = v.as_array();
      if (a.empty()) {
        out += "[]";
        break;
      }
      out += '[';
      for (std::size_t i = 0; i < a.size(); ++i) {
        if (i) out += ',';
        newline(depth + 1);
        dump_value(a[i], out, indent, depth + 1);
      }
      newline(depth);
      out += ']';
      break;
    }
    case Value::Type::Object: {
      const Members& m = v.members();
      if (m.empty()) {
        out += "{}";
        break;
      }
      out += '{';
      bool first = true;
      for (const auto& [k, val] : m) {
        if (!first) out += ',';
        first = false;
        newline(depth + 1);
        dump_string(k, out);
        out += ':';
        if (indent > 0) out += ' ';
        dump_value(val, out, indent, depth + 1);
      }
      newline(depth);
      out += '}';
      break;
    }
  }
}

// ---------------- parser ----------------

struct Parser {
  const char* p;
  const char* end;

  void ws() {
    while (p < end && (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r')) ++p;
  }
  [[noreturn]] void err(const std::string& msg) { fail("json parse: " + msg); }
  char peek() {
    if (p >= end) err("unexpected end of input");
    return *p;
  }
  void expect(char c) {
    if (p >= end || *p != c) err(std::string("expected '") + c + "'");
    ++p;
  }

  Value parse_value() {
    ws();
    char c = peek();
    switch (c) {
      case '{': return parse_object();
      case '[': return parse_array();
      case '"': return Value(parse_string());
      case 't':
        literal("true");
        return Value(true);
      case 'f':
        literal("false");
        return Value(false);
      case 'n':
        literal("null");
        return Value(nullptr);
      default: return parse_number();
    }
  }

  void literal(const char* lit) {
    std::size_t n = std::strlen(lit);
    if (static_cast<std::size_t>(end - p) < n || std::strncmp(p, lit, n) != 0)
      err(std::string("invalid literal, expected ") + lit);
    p += n;
  }

  std::string parse_string() {
    expect('"');
    std::string s;
    while (true) {
      if (p >= end) err("unterminated string");
      unsigned char c = static_cast<unsigned char>(*p++);
      if (c == '"') break;
      if (c == '\\') {
        if (p >= end) err("bad escape");
        char e = *p++;
        switch (e) {
          case '"': s += '"'; break;
          case '\\': s += '\\'; break;
          case '/': s += '/'; break;
          case 'b': s += '\b'; break;
          case 'f': s += '\f'; break;
          case 'n': s += '\n'; break;
          case 'r': s += '\r'; break;
          case 't': s += '\t'; break;
          case 'u': {
            if (end - p < 4) err("bad \\u escape");
            unsigned code = 0;
            for (int i = 0; i < 4; ++i) {
              char h = *p++;
              code <<= 4;
              if (h >= '0' && h <= '9') code |= static_cast<unsigned>(h - '0');
              else if (h >= 'a' && h <= 'f') code |= static_cast<unsigned>(h - 'a' + 10);
              else if (h >= 'A' && h <= 'F') code |= static_cast<unsigned>(h - 'A' + 10);
              else err("bad hex digit in \\u escape");
            }
            // UTF-8 encode (surrogate pairs folded to replacement — model files are ASCII/UTF-8 BMP).
            if (code < 0x80) {
              s += static_cast<char>(code);
            } else if (code < 0x800) {
              s += static_cast<char>(0xC0 | (code >> 6));
              s += static_cast<char>(0x80 | (code & 0x3F));
            } else {
              s += static_cast<char>(0xE0 | (code >> 12));
              s += static_cast<char>(0x80 | ((code >> 6) & 0x3F));
              s += static_cast<char>(0x80 | (code & 0x3F));
            }
            break;
          }
          default: err("bad escape");
        }
      } else {
        s += static_cast<char>(c);
      }
    }
    return s;
  }

  Value parse_number() {
    const char* start = p;
    if (p < end && (*p == '-' || *p == '+')) ++p;
    while (p < end && ((*p >= '0' && *p <= '9') || *p == '.' || *p == 'e' || *p == 'E' ||
                       *p == '+' || *p == '-'))
      ++p;
    if (p == start) err("invalid number");
    std::string tok(start, p);
    char* parsed_end = nullptr;
    double v = std::strtod(tok.c_str(), &parsed_end);
    if (parsed_end != tok.c_str() + tok.size()) err("invalid number '" + tok + "'");
    return Value(v);
  }

  Value parse_array() {
    expect('[');
    Value a = Value::array();
    ws();
    if (peek() == ']') {
      ++p;
      return a;
    }
    while (true) {
      a.push_back(parse_value());
      ws();
      char c = peek();
      if (c == ',') {
        ++p;
        continue;
      }
      if (c == ']') {
        ++p;
        return a;
      }
      err("expected ',' or ']'");
    }
  }

  Value parse_object() {
    expect('{');
    Value o = Value::object();
    ws();
    if (peek() == '}') {
      ++p;
      return o;
    }
    while (true) {
      ws();
      std::string key = parse_string();
      ws();
      expect(':');
      o.set(key, parse_value());
      ws();
      char c = peek();
      if (c == ',') {
        ++p;
        continue;
      }
      if (c == '}') {
        ++p;
        return o;
      }
      err("expected ',' or '}'");
    }
  }
};

}  // namespace

std::string Value::dump(int indent) const {
  std::string out;
  dump_value(*this, out, indent, 0);
  return out;
}

Value Value::parse(const std::string& text) {
  Parser parser{text.data(), text.data() + text.size()};
  Value v = parser.parse_value();
  parser.ws();
  if (parser.p != parser.end) fail("json parse: trailing content");
  return v;
}

}  // namespace keta::json

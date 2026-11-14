// Minimal self-contained JSON DOM (parse + serialize). No external dependency.
// Numbers are double; object member order is preserved (vector of pairs) so
// serialized models diff cleanly and results are byte-deterministic.
#pragma once
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace keta::json {

class Value;
using Array = std::vector<Value>;
using Members = std::vector<std::pair<std::string, Value>>;

class Value {
 public:
  enum class Type { Null, Bool, Number, String, Array, Object };

  Value() : type_(Type::Null) {}
  Value(std::nullptr_t) : type_(Type::Null) {}
  Value(bool b) : type_(Type::Bool), bool_(b) {}
  Value(double n) : type_(Type::Number), num_(n) {}
  Value(int n) : type_(Type::Number), num_(n) {}
  Value(std::int64_t n) : type_(Type::Number), num_(static_cast<double>(n)) {}
  Value(std::size_t n) : type_(Type::Number), num_(static_cast<double>(n)) {}
  Value(const char* s) : type_(Type::String), str_(s) {}
  Value(std::string s) : type_(Type::String), str_(std::move(s)) {}
  Value(Array a) : type_(Type::Array), arr_(std::move(a)) {}
  Value(Members m) : type_(Type::Object), obj_(std::move(m)) {}

  static Value object() { return Value(Members{}); }
  static Value array() { return Value(Array{}); }

  Type type() const { return type_; }
  bool is_null() const { return type_ == Type::Null; }
  bool is_bool() const { return type_ == Type::Bool; }
  bool is_number() const { return type_ == Type::Number; }
  bool is_string() const { return type_ == Type::String; }
  bool is_array() const { return type_ == Type::Array; }
  bool is_object() const { return type_ == Type::Object; }

  bool as_bool() const;
  double as_number() const;
  std::int64_t as_int() const;
  const std::string& as_string() const;
  const Array& as_array() const;
  Array& as_array();
  const Members& members() const;

  // Object access. get() returns nullptr when absent.
  const Value* get(const std::string& key) const;
  const Value& at(const std::string& key) const;  // throws keta::Error when absent
  Value& set(const std::string& key, Value v);    // insert or replace
  bool has(const std::string& key) const { return get(key) != nullptr; }

  // Typed convenience getters with defaults (for optional model fields).
  double num_or(const std::string& key, double dflt) const;
  bool bool_or(const std::string& key, bool dflt) const;
  std::string str_or(const std::string& key, const std::string& dflt) const;

  void push_back(Value v);

  std::string dump(int indent = 0) const;
  static Value parse(const std::string& text);  // throws keta::Error on malformed input

 private:
  Type type_;
  bool bool_ = false;
  double num_ = 0.0;
  std::string str_;
  Array arr_;
  Members obj_;
};

}  // namespace keta::json

#pragma once

#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

namespace json_deserializer {

using Number = std::variant<long, double>;
using String = std::string;
using Bool = bool;
using Null = std::nullptr_t;
struct Value;
using Object = std::map<String, Value>;
using Array = std::vector<Value>;

struct Value : public std::variant<Number, String, Bool, Null, Array, Object> {
  using base_type = std::variant<Number, String, Bool, Null, Array, Object>;
  using base_type::base_type;
};

inline void
skip_spaces(std::istream& is) {
  while (std::isspace(is.peek()))
    is.get();

  if (not is.good())
    throw std::invalid_argument("unexpected end of stream");
}

bool
deserialize(std::istream& is, Number& number) {
  using int_type = std::istream::int_type;
  bool is_float = false;
  bool has_exp = false;
  bool has_negative_exp = false;
  std::stringstream ss;

  skip_spaces(is);

  if (int_type c = is.peek(); c == '-') {
    ss.put(static_cast<char>(is.get()));
  } else if (not is.good())
    throw std::invalid_argument("unexpected end of stream");

  for (;;) {
    int_type c_int = is.peek();
    if (is.good()) {
      char c = static_cast<char>(c_int);

      if (c == '.') {
        if (is_float)
          throw std::invalid_argument("more than one '.'");
        else
          is_float = true;
      } else if (std::tolower(c) == 'e') {
        if (has_exp)
          throw std::invalid_argument("more than one 'e'");
        is_float = true;
        has_exp = true;
      } else if (c == '-') {
        if (has_negative_exp)
          throw std::invalid_argument("more than one '-' in exponent value");
        has_negative_exp = true;
      } else if (not std::isdigit(c)) {
        break;
      }

      ss.put(static_cast<char>(is.get()));
    } else if (is.eof())
      break;
    else
      throw std::invalid_argument("unexpected end of stream");
  }

  if (is_float) {
    number = std::stod(ss.str());
  } else {
    try {
      number = std::stol(ss.str());
    } catch (...) {
      return false;
    }
  }

  return true;
}

bool
deserialize(std::istream& is, String& out) {
  using int_type = std::istream::int_type;

  std::stringstream ss;
  skip_spaces(is);

  if (int_type c = is.peek(); c != '"')
    return false;
  is.get();

  for (;;) {
    int_type c = is.get();
    if (not is.good())
      throw std::invalid_argument("unexpected end of stream");

    if (c == '"')
      break;
    else
      ss.put(static_cast<char>(c));
  }

  out = ss.str();

  return true;
}

bool
deserialize(std::istream& is, Bool& out) {
  std::stringstream ss;
  skip_spaces(is);

  constexpr std::size_t bufsize = 4;
  std::istream::char_type buf[bufsize];
  is.read(buf, bufsize);
  auto read_bytes = is.gcount();

  if (not is.good() or read_bytes != 4) {
    while (read_bytes > 0) {
      --read_bytes;
      is.putback(buf[read_bytes]);
    }
    return false;
  }

  std::transform(buf, buf + read_bytes, buf,
                 [](char c) { return std::tolower(c); });
  if (std::equal(buf, buf + read_bytes, "true")) {
    out = true;
  } else if (std::equal(buf, buf + bufsize, "fals") and
             std::tolower(is.peek()) == 'e') {
    is.get();
    out = false;
  } else {
    while (read_bytes > 0) {
      --read_bytes;
      is.putback(buf[read_bytes]);
    }

    return false;
  }

  return true;
}

bool
deserialize(std::istream& is, Null) {
  skip_spaces(is);

  constexpr std::size_t bufsize = 4;
  std::istream::char_type buf[bufsize];
  is.read(buf, bufsize);
  auto read_bytes = is.gcount();

  if (not is.good() or read_bytes != 4) {
    while (read_bytes > 0) {
      --read_bytes;
      is.putback(buf[read_bytes]);
    }
    return false;
  }

  if (std::equal(buf, buf + bufsize, "null")) {
    return true;
  } else {
    while (read_bytes > 0) {
      --read_bytes;
      is.putback(buf[read_bytes]);
    }
    return false;
  }
}

bool deserialize(std::istream& is, Value& value);

bool
deserialize(std::istream& is, Object& object) {
  using int_type = std::istream::int_type;
  object.clear();
  skip_spaces(is);

  if (int_type c = is.peek(); c != '{')
    return false;
  is.get();

  skip_spaces(is);

  if (int_type c = is.peek(); c == '}') {
    is.get();
    return true;
  }

  for (;;) {
    String name;
    if (not deserialize(is, name))
      throw std::invalid_argument("expected name");

    skip_spaces(is);

    if (int_type c = is.get(); c != ':')
      throw std::invalid_argument("expected colon");

    skip_spaces(is);

    Value value;
    if (not deserialize(is, value))
      throw std::invalid_argument("expected value");

    object.emplace(std::move(name), std::move(value));

    skip_spaces(is);

    if (int_type c = is.get(); c == '}') {
      break;
    } else if (c != ',')
      throw std::invalid_argument("expected comma");
  }

  return true;
}

bool
deserialize(std::istream& is, Array& array) {
  using int_type = std::istream::int_type;
  array.clear();
  skip_spaces(is);

  if (int_type c = is.peek(); c != '[')
    return false;
  is.get();

  if (int_type c = is.peek(); c == ']') {
    is.get();
    return true;
  }

  for (;;) {
    skip_spaces(is);
    Value value;
    if (not deserialize(is, value))
      throw std::invalid_argument("expected value");

    array.emplace_back(std::move(value));

    skip_spaces(is);

    if (int_type c = is.get(); c == ']') {
      break;
    } else if (c != ',')
      throw std::invalid_argument("expected comma");
  }

  return true;
}

bool
deserialize(std::istream& is, Value& value) {
  if (Array array; deserialize(is, array))
    value = std::move(array);
  else if (Object object; deserialize(is, object))
    value = std::move(object);
  else if (Bool b; deserialize(is, b))
    value = b;
  else if (deserialize(is, Null{}))
    value = Null{};
  else if (String string; deserialize(is, string))
    value = std::move(string);
  else if (Number number; deserialize(is, number))
    value = number;
  else
    return false;

  return true;
}

} // namespace json_deserializer

#pragma once

#include "jsonify_base.hpp"

#include <algorithm>

namespace results {

template <typename CharT, typename Traits, typename T>
std::enable_if_t<detail::is_jsonificable_v<std::decay_t<T>>,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, T &&t) {
  namespace rng = std::ranges;
  using t_nocvref = std::remove_cv_t<std::remove_reference_t<T>>;

  if constexpr (std::is_same_v<t_nocvref, bool>) {
    if (t) {
      return os << "true";
    } else {
      return os << "false";
    }
  } else if constexpr (std::is_arithmetic_v<t_nocvref>) {
    if constexpr (sizeof(t_nocvref) == sizeof(char)) {
      return os << static_cast<std::conditional_t<std::is_signed_v<t_nocvref>,
                                                  int, unsigned>>(t);
    } else {
      IosSaver iosSaver(os);
      return os << std::fixed << t;
    }
  } else if constexpr (detail::is_optional_v<t_nocvref>) {
    if (t)
      return jsonify(os, *t);
    else
      return os << "null";
  } else if constexpr (detail::is_stringlike_v<t_nocvref>) {
    os << '"';
    return detail::streamEscapedString(
               os, std::basic_string_view<CharT, Traits>(std::forward<T>(t)))
           << '"';
  } else if constexpr (detail::is_maplike_v<t_nocvref>) {
    static_assert(detail::is_stringlike_v<typename t_nocvref::key_type>);
    os << '{';
    const auto last = rng::prev(rng::end(t));

    auto streamElement = [&](const auto &element) {
      const auto &[key, value] = element;
      jsonify(os, key, value);
    };

    auto iter = rng::begin(t);
    if (iter != rng::end(t)) {
      streamElement(*iter++);
      rng::for_each(iter, rng::end(t), [&](const auto &element) {
        os << ',';
        streamElement(element);
      });
    }

    os << '}';
    return os;
  } else if constexpr (detail::is_iterable_v<t_nocvref>) {
    os << '[';
    auto iter = rng::begin(t);
    if (iter != rng::end(t)) {
      jsonify(os, *iter++);
      rng::for_each(iter, rng::end(t), [&](const auto &element) {
        os << ',';
        jsonify(os, element);
      });
    }

    os << ']';
    return os;
  } else
    static_assert(detail::false_var_v<T>,
                  "T is not neither an arithmetic value, a string-like object, "
                  "a map-like object or an iterable object");
}

template <typename CharT, typename Traits, typename Key, typename Value>
std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os, Key &&key, Value &&value) {
  static_assert(detail::is_stringlike_v<std::decay_t<Key>>);

  jsonify(os, std::forward<Key>(key)) << ':';
  jsonify(os, std::forward<Value>(value));

  return os;
}

namespace detail {

template <typename CharT, typename Traits, typename Arg1, typename Arg2,
          typename... Args>
std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os, Arg1 &&arg1, Arg2 &&arg2,
        Args &&...args) {
  static_assert(detail::is_stringlike_v<std::decay_t<Arg1>>);

  results::jsonify(os, std::forward<Arg1>(arg1), std::forward<Arg2>(arg2));
  if constexpr (sizeof...(args) > 0) {
    os << ',';
    return jsonify(os, std::forward<Args>(args)...);
  } else
    return os;
}

} // namespace detail

template <typename CharT, typename Traits, typename... Args>
std::enable_if_t<(sizeof...(Args) > 3) and (sizeof...(Args) % 2) == 0,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, Args &&...args) {
  return detail::jsonify(os, std::forward<Args>(args)...);
}

} // namespace results

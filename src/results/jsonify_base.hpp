#pragma once

#include "../ios_saver.hpp"

#include <optional>
#include <ostream>
#include <range/v3/iterator/concepts.hpp>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace results {

namespace detail {

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &
streamEscapedString(std::basic_ostream<CharT, Traits> &os,
                    const std::basic_string_view<CharT, Traits> str) {
  using str_t = decltype(str);

  typename str_t::size_type pos = 0;
  for (; pos < str.size();) {
    auto endPos = str.find('"', pos);
    if (endPos == str_t::npos) {
      os << str.substr(pos);
      return os;
    }

    if (endPos != pos)
      os << str.substr(pos, endPos - pos);
    os << "\\\"";
    pos = endPos + 1;
  }

  return os;
}

template <typename> struct is_optional : std::false_type {};

template <typename T> struct is_optional<std::optional<T>> : std::true_type {};

template <typename T> constexpr bool is_optional_v = is_optional<T>::value;

template <typename, typename = void>
struct is_maplike : std::bool_constant<false> {};

template <typename T>
struct is_maplike<
    T, std::void_t<
           typename T::key_type, typename T::mapped_type,
           std::enable_if_t<std::is_same_v<
               typename T::value_type, std::pair<const typename T::key_type,
                                                 typename T::mapped_type>>>>>
    : std::bool_constant<true> {};

template <typename T> constexpr bool is_maplike_v = is_maplike<T>::value;

template <typename, typename = void>
struct is_basic_string : std::bool_constant<false> {};

template <typename T>
struct is_basic_string<
    T, std::enable_if_t<std::is_same_v<
           T, std::basic_string<typename T::value_type, typename T::traits_type,
                                typename T::allocator_type>>>>
    : std::bool_constant<true> {};

template <typename, typename = void>
struct is_basic_string_view : std::bool_constant<false> {};

template <typename T>
struct is_basic_string_view<
    T, std::enable_if_t<
           std::is_same_v<T, std::basic_string_view<typename T::value_type,
                                                    typename T::traits_type>>>>
    : std::bool_constant<true> {};

template <typename, typename = void>
struct is_stringlike : std::bool_constant<false> {};

template <typename T>
struct is_stringlike<T, std::enable_if_t<std::disjunction_v<
                            std::is_convertible<T, const char *>,
                            is_basic_string<T>, is_basic_string_view<T>>>>
    : std::bool_constant<true> {};

template <typename T> constexpr bool is_stringlike_v = is_stringlike<T>::value;

static_assert(is_stringlike_v<char[5]>);
static_assert(is_stringlike_v<std::string>);
static_assert(is_stringlike_v<std::string_view>);

template <typename, typename = void>
struct is_iterable : std::bool_constant<false> {};

template <typename T>
struct is_iterable<T, std::enable_if_t<ranges::Iterator<typename T::iterator>>>
    : std::bool_constant<true> {};

template <typename T>
struct is_iterable<T,
                   std::enable_if_t<std::is_pointer_v<T> or std::is_array_v<T>>>
    : std::bool_constant<true> {};

template <typename T> constexpr bool is_iterable_v = is_iterable<T>::value;

struct nonesuch {};

template <typename, typename = void> struct key_type_if_exists {
  using type = nonesuch;
};

template <typename T>
struct key_type_if_exists<T, std::void_t<typename T::key_type>> {
  using type = typename T::key_type;
};
template <typename T>
using key_type_if_exists_t = typename key_type_if_exists<T>::type;

template <typename, typename = void> struct value_type_if_exists {
  using type = nonesuch;
};

template <typename T>
struct value_type_if_exists<T, std::void_t<typename T::value_type>> {
  using type = typename T::value_type;
};
template <typename T>
using value_type_if_exists_t = typename value_type_if_exists<T>::type;

template <typename, typename = void>
struct is_jsonificable : std::bool_constant<false> {};

template <typename T>
struct is_jsonificable<
    T,
    std::enable_if_t<std::is_arithmetic_v<T> or is_optional_v<T> or
                     is_stringlike_v<T> or is_maplike_v<T> or is_iterable_v<T>>>
    : std::bool_constant<true> {};

template <typename T>
constexpr bool is_jsonificable_v = is_jsonificable<T>::value;

static_assert(is_jsonificable_v<char[5]>);

template <typename...> struct false_var : std::false_type {};

template <typename... T> constexpr bool false_var_v = false_var<T...>::value;

} // namespace detail

template <typename CharT, typename Traits, typename T>
std::enable_if_t<detail::is_jsonificable_v<std::decay_t<T>>,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, T &&t);

template <typename CharT, typename Traits, typename Key, typename Value>
std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os, Key &&key, Value &&value);

template <typename CharT, typename Traits, typename... Args>
std::enable_if_t<(sizeof...(Args) > 3) and (sizeof...(Args) % 2) == 0,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, Args &&...args);

} // namespace results

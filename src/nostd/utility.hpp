#pragma once

#include "type_traits.hpp"

#include <utility>
#include <vector>

namespace nostd {

namespace detail {

template <typename Seq> struct make_integer_sequence_rev;

template <typename T, T... Ns>
struct make_integer_sequence_rev<std::integer_sequence<T, Ns...>>
    : std::integer_sequence<T, (sizeof...(Ns) - Ns - 1)...> {};

} // namespace detail

template <typename T, T N>
using make_integer_sequence_rev = typename detail::make_integer_sequence_rev<
    std::make_integer_sequence<T, N>>;

template <std::size_t N>
using make_index_sequence_rev = make_integer_sequence_rev<std::size_t, N>;

/*
template <std::size_t Index, typename T>
struct Enumerated {
  static constexpr std::size_t index = Index;
  using type = T;

  template <typename U>
  explicit Enumerated(U&& u) noexcept(noexcept(T(std::forward<U>(u))))
      : value(std::forward<U>(u)) {}

  T value;
};

template <std::size_t Index, typename T>
auto
make_enumerated(T&& t) {
  using type = std::conditional_t<std::is_rvalue_reference_v<T>,
                                  std::remove_reference_t<T>, T>;
  return Enumerated<Index, type>(std::forward<T>(t));
}

static_assert(std::is_same_v<decltype(make_enumerated<0>(3).value), int>);
static_assert(std::is_same_v<
              decltype(make_enumerated<0>(std::declval<int&>()).value), int&>);
static_assert(std::is_same_v<
              decltype(make_enumerated<0>(std::declval<int&&>()).value), int>);

namespace detail {

template <typename Fn, std::size_t... Idx, typename... Ts>
auto
call_enumerated_args(Fn fn, std::index_sequence<Idx...>, Ts&&... ts) {
  std::tuple<Ts...> args(std::forward<Ts>(ts)...);

  return fn(make_enumerated<Idx>(
      std::forward<std::tuple_element_t<Idx, decltype(args)>>(args))...);
}

} // namespace detail

template <typename Fn, typename... Ts>
auto
call_enumerated_args(Fn fn, Ts&&... ts) {
  return detail::call_enumerated_args(std::move(fn),
                                      std::make_index_sequence<sizeof...(Ts)>(),
                                      std::forward<Ts>(ts)...);
}
*/

template <typename T>
constexpr std::make_signed_t<T> as_signed(const T &t) noexcept {
  using signed_type = std::make_signed_t<T>;
  assert(t <= std::numeric_limits<signed_type>::max());
  return static_cast<signed_type>(t);
}

template <typename T>
constexpr std::make_unsigned_t<T> as_unsigned(const T &t) noexcept {
  assert(t >= 0);
  return static_cast<std::make_unsigned_t<T>>(t);
}

template <typename T>
constexpr make_signed_upcast_t<T> as_signed_upcast(const T &t) noexcept {
  using signed_type = make_signed_upcast_t<T>;
  return static_cast<signed_type>(t);
}

template <typename T>
constexpr make_unsigned_downcast_t<T>
as_unsigned_downcast(const T &t) noexcept {
  using unsigned_type = make_unsigned_downcast_t<T>;
  assert(t >= 0);
  assert(t <= std::numeric_limits<unsigned_type>::max());
  return static_cast<unsigned_type>(t);
}

} // namespace nostd

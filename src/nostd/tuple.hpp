#pragma once

#include "type_traits.hpp"

#include <cstddef>
#include <tuple>
#include <utility>

namespace nostd {

namespace detail {

template <typename Tuple, typename Indices> struct tuple_until_impl;

template <typename Tuple, std::size_t... Idx>
struct tuple_until_impl<Tuple, std::index_sequence<Idx...>> {
  using type = std::tuple<std::tuple_element_t<Idx, Tuple>...>;
};

template <typename Fun, typename Args, std::size_t... Idx>
auto apply_until(Fun &&fun, Args &&args, std::index_sequence<Idx...>) noexcept(
    noexcept(fun(std::get<Idx>(std::forward<Args>(args))...))) {
  return fun(std::get<Idx>(std::forward<Args>(args))...);
}

} // namespace detail

template <std::size_t Index, typename... Ts>
struct tuple_from_args_until
    : detail::tuple_until_impl<std::tuple<Ts...>,
                               std::make_index_sequence<Index>> {};

template <std::size_t Index, typename... Ts>
using tuple_from_args_until_t =
    typename tuple_from_args_until<Index, Ts...>::type;

template <std::size_t Index, typename Fun, typename Args>
auto apply_until(Fun &&fun, Args &&args) noexcept(noexcept(
    detail::apply_until(std::forward<Fun>(fun), std::forward<Args>(args),
                        std::make_index_sequence<Index>()))) {
  return detail::apply_until(std::forward<Fun>(fun), std::forward<Args>(args),
                             std::make_index_sequence<Index>());
}

} // namespace nostd

#pragma once

#include "vec2d_element_builder_ct_parts_traits.hpp"

#include <tuple>

namespace het::detail {

struct vec2d_element_empty {};

template <auto Size>
struct vec2d_element_fixed_size {
  using size_type = decltype(Size);
  static constexpr size_type size = Size;
};

struct vec2d_element_dynamic_size {};

template <typename Fun>
struct vec2d_element_dynamic_size_from_callable : Fun {
  using Fun::operator();

  template <typename _Fun>
  constexpr vec2d_element_dynamic_size_from_callable(_Fun&& fun)
      : Fun(std::forward<_Fun>(fun)) {}
};

template <typename SizeType>
struct checked_vec2d_size_type : nostd::identity_type<SizeType> {
  static_assert(is_vec2d_size_type_v<SizeType>,
                "SizeType must be vec2d_element_static_size or "
                "vec2d_element_dynamic_size");
};

template <typename SizeType>
using checked_vec2d_size_type_t =
    typename checked_vec2d_size_type<SizeType>::type;

template <typename SizeType>
struct vec2d_element_default_uninitialized
    : checked_vec2d_size_type_t<SizeType> {

  template <typename _SizeType>
  constexpr vec2d_element_default_uninitialized(_SizeType&& sizeType)
      : checked_vec2d_size_type_t<SizeType>(std::forward<_SizeType>(sizeType)) {
  }
};

template <typename SizeType>
struct vec2d_element_default_construction_default
    : checked_vec2d_size_type_t<SizeType> {

  template <typename _SizeType>
  constexpr vec2d_element_default_construction_default(_SizeType&& sizeType)
      : checked_vec2d_size_type_t<SizeType>(std::forward<_SizeType>(sizeType)) {
  }
};

template <typename SizeType, typename... Values>
struct vec2d_element_default_construction
    : checked_vec2d_size_type_t<SizeType> {
  using values_type = std::tuple<Values...>;

  template <std::size_t Index>
  using value_type = std::tuple_element_t<Index, values_type>;

  constexpr vec2d_element_default_construction() = default;
  constexpr vec2d_element_default_construction(
      const vec2d_element_default_construction&) = default;
  constexpr vec2d_element_default_construction(
      vec2d_element_default_construction&&) = default;

  template <typename _SizeType, typename... Us,
            typename = std::enable_if_t<not std::is_same_v<
                std::decay_t<_SizeType>, vec2d_element_default_construction>>>
  constexpr vec2d_element_default_construction(_SizeType&& sizeType,
                                               Us&&... values)
      : checked_vec2d_size_type_t<SizeType>(std::forward<_SizeType>(sizeType)),
        init_values(std::forward<Us>(values)...) {}

  std::tuple<Values...> init_values;
};

template <typename SizeType, typename... Values>
inline auto
make_vec2d_element_default_construction(Values&&... values) noexcept(
    noexcept(std::tuple(std::forward<Values>(values)...))) {
  // TODO: handle reference_wrapper
  return vec2d_element_default_construction<SizeType, std::decay_t<Values>...>(
      std::forward<Values>(values)...);
}

template <template <typename...> typename Parent, typename... Ts>
struct vec2d_build_parts_common : std::tuple<Ts...> {
  template <std::size_t Idx>
  using type = std::tuple_element_t<Idx, std::tuple<Ts...>>;
  static constexpr std::size_t arity = sizeof...(Ts);
  using last_type = type<arity - 1>;

  constexpr vec2d_build_parts_common() = default;
  constexpr vec2d_build_parts_common(const vec2d_build_parts_common&) = default;
  constexpr vec2d_build_parts_common(vec2d_build_parts_common&&) = default;

  template <typename... Us>
  explicit constexpr vec2d_build_parts_common(Us&&... us)
      : std::tuple<Ts...>(std::forward<Us>(us)...) {
    static_assert(sizeof...(Us) == sizeof...(Ts));
  }

  template <typename T>
  constexpr auto
  replace_last_with(T&& t) {
    static_assert(sizeof...(Ts) > 0);
    return replace_last_with_impl(
        std::forward<T>(t), std::make_index_sequence<sizeof...(Ts) - 1>());
  }

  constexpr decltype(auto)
  take_last() {
    if constexpr (std::is_reference_v<last_type>)
      return std::get<sizeof...(Ts) - 1>(*this);
    else
      return std::move(std::get<sizeof...(Ts) - 1>(*this));
  }

  template <typename T>
  constexpr auto
  append(T&& t) {
    return append_impl(std::forward<T>(t), std::index_sequence_for<Ts...>());
  }

private:
  template <typename Indices, typename T>
  struct replace_last_with_impl_t;

public:
  template <typename T>
  using replace_last_with_t = typename replace_last_with_impl_t<
      std::make_index_sequence<sizeof...(Ts) - 1>, T>::type;

  template <typename T>
  using append_t = Parent<Ts..., T>;

private:
  template <typename T, std::size_t... Idx>
  constexpr auto
  replace_last_with_impl(T&& t, std::index_sequence<Idx...>) {
    static_assert(sizeof...(Idx) == sizeof...(Ts) - 1);

    return Parent<type<Idx>..., std::decay_t<T>>(
        std::get<Idx>(std::move(*this))..., std::forward<T>(t));
  }

  template <typename T, std::size_t... Idx>
  constexpr append_t<T>
  append_impl(T&& t, std::index_sequence<Idx...>) {
    return append_t<std::decay_t<T>>(std::get<Idx>(std::move(*this))...,
                                     std::forward<T>(t));
  }

  template <typename T, std::size_t... Idx>
  struct replace_last_with_impl_t<std::index_sequence<Idx...>, T> {
    static_assert(sizeof...(Idx) == sizeof...(Ts) - 1);

    using type = Parent<std::tuple_element_t<Idx, std::tuple<Ts...>>..., T>;
  };
};

template <template <typename...> typename Parent>
struct vec2d_build_parts_common<Parent> : std::tuple<> {
  template <std::size_t>
  using type = void;
  static constexpr std::size_t arity = 0;
  using last_type = void;

  constexpr vec2d_build_parts_common() = default;

  template <typename T>
  constexpr auto
  append(T&& t) {
    return Parent<std::decay_t<T>>(std::forward<T>(t));
  }

  template <typename T>
  using append_t = Parent<T>;
};

template <typename... Ts>
struct vec2d_build_parts : vec2d_build_parts_common<vec2d_build_parts, Ts...> {
  using base_type = vec2d_build_parts_common<vec2d_build_parts, Ts...>;
  template <std::size_t Idx>
  using type = std::tuple_element_t<Idx, std::tuple<Ts...>>;
  static constexpr std::size_t arity = sizeof...(Ts);
  using last_type = type<arity - 1>;
  using base_type::base_type;
};

} // namespace het::detail

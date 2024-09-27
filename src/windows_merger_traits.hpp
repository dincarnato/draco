#pragma once

#include "tiny_fraction.hpp"

#include "nostd/type_traits.hpp"

#include <cstdint>

namespace windows_merger {

template <typename T> struct WindowsMergerResize {
  using type = T;
  T size;

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 std::is_convertible_v<std::decay_t<U>, T> and
                                 std::is_scalar_v<T>,
                             int> = 0>
  constexpr WindowsMergerResize(U &&u) : size(static_cast<T>(u)) {}

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 std::is_convertible_v<std::decay_t<U>, T> and
                                 not std::is_scalar_v<T>,
                             int> = 0>
  constexpr WindowsMergerResize(U &&u) : size(std::forward<U>(u)) {}

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 not std::is_convertible_v<std::decay_t<U>, T>,
                             int> = 0>
  constexpr explicit WindowsMergerResize(U &&u) : size(std::forward<U>(u)) {}
};

template <typename T>
WindowsMergerResize(T &&size) -> WindowsMergerResize<std::decay_t<T>>;

template <typename T> struct WindowsMergerReserve {
  using type = T;
  T size;

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 std::is_convertible_v<std::decay_t<U>, T> and
                                 std::is_scalar_v<T>,
                             int> = 0>
  constexpr WindowsMergerReserve(U &&u) : size(static_cast<T>(u)) {}

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 std::is_convertible_v<std::decay_t<U>, T> and
                                 not std::is_scalar_v<T>,
                             int> = 0>
  constexpr WindowsMergerReserve(U &&u) : size(std::forward<U>(u)) {}

  template <typename U,
            std::enable_if_t<std::is_convertible_v<T, std::decay_t<U>> and
                                 not std::is_convertible_v<std::decay_t<U>, T>,
                             int> = 0>
  constexpr explicit WindowsMergerReserve(U &&u) : size(std::forward<U>(u)) {}
};

template <typename T>
WindowsMergerReserve(T &&size) -> WindowsMergerReserve<std::decay_t<T>>;

struct WindowsMergerWindow;
struct WindowsMergerWindowBase;

template <typename> struct WindowsMergerWindowAccessor;

template <typename> struct WindowsMergerWindowBaseAccessor;

struct WindowsMergerTraits {
  using clusters_size_type = std::uint8_t;
  using bases_size_type = std::uint16_t;
  using windows_size_type = std::size_t;
  using weight_type = TinyFraction;
  using coverage_type = unsigned;

  template <typename T>
  using weight_getter_t =
      decltype(std::declval<const T &>().weight(clusters_size_type()));

  template <typename T>
  using coverage_getter_t = decltype(std::declval<const T &>().coverage());

  template <typename T>
  using clusters_size_getter_t =
      decltype(std::declval<const T &>().clusters_size());

  template <typename T>
  static constexpr bool is_window_baselike_v =
      nostd::is_detected_convertible_v<const weight_type &, weight_getter_t,
                                       T> and
      nostd::is_detected_convertible_v<const coverage_type &, coverage_getter_t,
                                       T> and
      nostd::is_detected_convertible_v<clusters_size_type,
                                       clusters_size_getter_t, T>;

  template <typename T, typename = void>
  struct is_lvalue_window_baselike_reference : std::false_type {};

  template <typename T>
  static constexpr bool is_lvalue_window_baselike_reference_v =
      is_lvalue_window_baselike_reference<T>::value;

  template <typename T, typename = void>
  struct is_rvalue_window_baselike_reference : std::false_type {};

  template <typename T>
  static constexpr bool is_rvalue_window_baselike_reference_v =
      is_rvalue_window_baselike_reference<T>::value;

  template <typename T>
  using begin_index_getter_t =
      decltype(std::declval<const T &>().begin_index());

  template <typename T>
  using end_index_getter_t = decltype(std::declval<const T &>().end_index());

  template <typename T>
  using size_getter_t = decltype(std::declval<const T &>().size());

  template <typename T>
  using coverages_getter_t = decltype(std::declval<const T &>().coverages());

  template <typename T>
  using begin_index_setter_t = decltype(std::declval<T &>().set_begin_index(
      std::declval<bases_size_type>()));

  template <typename T>
  static constexpr bool is_windowlike_v =
      nostd::is_detected_exact_v<clusters_size_type, clusters_size_getter_t,
                                 T> and
      nostd::is_detected_exact_v<bases_size_type, begin_index_getter_t, T> and
      nostd::is_detected_exact_v<bases_size_type, end_index_getter_t, T> and
      nostd::is_detected_exact_v<bases_size_type, size_getter_t, T> and
      nostd::is_detected_v<coverages_getter_t, T> and
      nostd::is_detected_v<begin_index_setter_t, T>;

  template <typename T> struct is_window : std::false_type {};

  template <typename T> static constexpr bool is_window_v = is_window<T>::value;

  template <typename T>
  static constexpr bool is_cache_window_baselike_v =
      nostd::is_detected_convertible_v<const weight_type &, weight_getter_t,
                                       T> and
      not nostd::is_detected_convertible_v<const coverage_type &,
                                           coverage_getter_t, T> and
      nostd::is_detected_convertible_v<clusters_size_type,
                                       clusters_size_getter_t, T>;

  template <typename T> struct is_reshaper_arg : std::false_type {};

  template <typename T>
  struct is_reshaper_arg<WindowsMergerResize<T>> : std::true_type {};

  template <typename T>
  struct is_reshaper_arg<WindowsMergerReserve<T>> : std::true_type {};

  template <typename T>
  static constexpr bool is_reshaper_arg_v = is_reshaper_arg<T>::value;

  template <typename T> struct is_reserve_reshaper : std::false_type {};

  template <typename T>
  struct is_reserve_reshaper<WindowsMergerReserve<T>> : std::true_type {};

  template <typename T>
  static constexpr bool is_reserve_reshaper_v = is_reserve_reshaper<T>::value;

  template <typename T> struct is_resize_reshaper : std::false_type {};

  template <typename T>
  struct is_resize_reshaper<WindowsMergerResize<T>> : std::true_type {};

  template <typename T>
  static constexpr bool is_resize_reshaper_v = is_resize_reshaper<T>::value;
};

template <>
struct WindowsMergerTraits::is_lvalue_window_baselike_reference<
    WindowsMergerWindowBase &, void> : std::true_type {};

template <typename T>
struct WindowsMergerTraits::is_lvalue_window_baselike_reference<
    T,
    std::enable_if_t<std::is_same_v<
        std::decay_t<T>, WindowsMergerWindowBaseAccessor<WindowsMergerWindow>>>>
    : std::true_type {};

template <typename T>
struct WindowsMergerTraits::is_lvalue_window_baselike_reference<
    T, std::enable_if_t<std::is_same_v<
           std::decay_t<T>,
           WindowsMergerWindowBaseAccessor<WindowsMergerWindow const>>>>
    : std::true_type {};

template <>
struct WindowsMergerTraits::is_rvalue_window_baselike_reference<
    WindowsMergerWindowBase &&, void> : std::true_type {};

template <typename T>
struct WindowsMergerTraits::is_rvalue_window_baselike_reference<
    T, std::enable_if_t<std::is_same_v<
           std::decay_t<T>,
           WindowsMergerWindowBaseAccessor<WindowsMergerWindow &&>>>>
    : std::true_type {};

template <>
struct WindowsMergerTraits::is_window<WindowsMergerWindow> : std::true_type {};

} // namespace windows_merger

#pragma once

#include <iterator>
#include <tuple>
#include <type_traits>

namespace nostd {

template <typename T> struct is_reference_wrapper : std::false_type {};

template <typename T>
struct is_reference_wrapper<std::reference_wrapper<T>> : std::true_type {};

template <typename T, typename U, typename = void>
struct is_reference_wrapper_to : std::false_type {};

template <typename T, typename U>
struct is_reference_wrapper_to<
    std::reference_wrapper<T>, U,
    std::enable_if_t<std::is_same<std::reference_wrapper<T>, U>::value>>
    : std::true_type {};

template <typename T, typename = void>
struct is_equality_comparable : std::false_type {};

template <typename T>
struct is_equality_comparable<
    T, std::enable_if_t<std::is_convertible_v<
           decltype(std::declval<T>() == std::declval<T>()), bool>>>
    : std::true_type {};

template <typename Iter>
constexpr bool is_equality_comparable_v = is_equality_comparable<Iter>::value;

template <typename T, typename U> struct copy_const {
  using type = std::conditional_t<std::is_const_v<T>, std::add_const_t<U>,
                                  std::remove_const_t<U>>;
};

template <typename T, typename U>
using copy_const_t = typename copy_const<T, U>::type;

// Credits to cppreference.com
#ifndef DETECTOR_IDIOM
#define DETECTOR_IDIOM

struct nonesuch {
  nonesuch() = delete;
  ~nonesuch() = delete;
  nonesuch(nonesuch const &) = delete;
  void operator=(nonesuch const &) = delete;
};

namespace detail {
template <class Default, class AlwaysVoid, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
  // Note that std::void_t is a C++17 feature
  using value_t = std::true_type;
  using type = Op<Args...>;
};

} // namespace detail

template <template <class...> class Op, class... Args>
using is_detected =
    typename detail::detector<nonesuch, void, Op, Args...>::value_t;

template <template <class...> class Op, class... Args>
using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;

template <class Default, template <class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

template <template <class...> class Op, class... Args>
constexpr bool is_detected_v = is_detected<Op, Args...>::value;

template <class Default, template <class...> class Op, class... Args>
using detected_or_t = typename detected_or<Default, Op, Args...>::type;

template <class Expected, template <class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

template <class Expected, template <class...> class Op, class... Args>
constexpr bool is_detected_exact_v =
    is_detected_exact<Expected, Op, Args...>::value;

template <class To, template <class...> class Op, class... Args>
using is_detected_convertible =
    std::is_convertible<detected_t<Op, Args...>, To>;

template <class To, template <class...> class Op, class... Args>
constexpr bool is_detected_convertible_v =
    is_detected_convertible<To, Op, Args...>::value;

#endif

template <typename T, std::size_t = 0> struct identity_type { using type = T; };

template <typename T, std::size_t Index = 0>
using identity_type_t = typename identity_type<T, Index>::type;

template <typename T> struct is_integral_constant : std::false_type {};

template <typename T, T Value>
struct is_integral_constant<std::integral_constant<T, Value>> : std::true_type {
};

template <typename T>
constexpr bool is_integral_constant_v = is_integral_constant<T>::value;

namespace detail {

template <typename Fun, typename Arg, typename Seq>
struct is_applicable_array_impl;

template <typename Fun, typename Arg, std::size_t... Idx>
struct is_applicable_array_impl<Fun, Arg, std::index_sequence<Idx...>>
    : std::is_invocable<Fun, identity_type_t<Arg, Idx>...> {};

template <typename R, typename Fun, typename Arg, typename Seq>
struct is_applicable_r_array_impl;

template <typename R, typename Fun, typename Arg, std::size_t... Idx>
struct is_applicable_r_array_impl<R, Fun, Arg, std::index_sequence<Idx...>>
    : std::is_invocable_r<R, Fun, identity_type_t<Arg, Idx>...> {};

} // namespace detail

template <typename Fun, typename Args>
struct is_applicable : std::false_type {};

template <typename Fun, typename... Ts>
struct is_applicable<Fun, std::tuple<Ts...>> : std::is_invocable<Fun, Ts...> {};

template <typename Fun, typename T, std::size_t Size>
struct is_applicable<Fun, std::array<T, Size>>
    : detail::is_applicable_array_impl<Fun, T, std::make_index_sequence<Size>> {
};

template <typename Fun, typename Args>
constexpr bool is_applicable_v = is_applicable<Fun, Args>::value;

template <typename R, typename Fun, typename Args>
struct is_applicable_r : std::false_type {};

template <typename R, typename Fun, typename... Ts>
struct is_applicable_r<R, Fun, std::tuple<Ts...>>
    : std::is_invocable_r<R, Fun, Ts...> {};

template <typename R, typename Fun, typename T, std::size_t Size>
struct is_applicable_r<R, Fun, std::array<T, Size>>
    : detail::is_applicable_r_array_impl<R, Fun, T,
                                         std::make_index_sequence<Size>> {};

template <typename R, typename Fun, typename Args>
constexpr bool is_applicable_r_v = is_applicable_r<R, Fun, Args>::value;

template <typename Fun, typename Args> struct apply_result;

template <typename Fun, typename... Ts>
struct apply_result<Fun, std::tuple<Ts...>> : std::invoke_result<Fun, Ts...> {};

template <typename Fun, typename Args>
using apply_result_t = typename apply_result<Fun, Args>::type;

template <typename T> struct is_tuple : std::false_type {};

template <typename... Ts>
struct is_tuple<std::tuple<Ts...>> : std::true_type {};

template <typename T> constexpr bool is_tuple_v = is_tuple<T>::value;

namespace detail {
template <typename T, typename, typename... Args>
struct is_constructible_impl : std::false_type {};

template <typename T, typename... Args>
struct is_constructible_impl<
    T, std::void_t<decltype(T(std::declval<Args>()...))>, Args...>
    : std::true_type {};

} // namespace detail

template <typename T, typename... Args>
struct is_constructible : detail::is_constructible_impl<T, void, Args...> {};

template <typename T, typename... Args>
constexpr bool is_constructible_v = is_constructible<T, Args...>::value;

template <typename T, typename Tuple>
struct is_constructible_from_tuple : std::false_type {};

// TODO: Args must be wrap-deferenced
template <typename T, typename... Args>
struct is_constructible_from_tuple<T, std::tuple<Args...>>
    : is_constructible<T, Args...> {};

template <typename T, typename Tuple>
constexpr bool is_constructible_from_tuple_v =
    is_constructible_from_tuple<T, Tuple>::value;

template <typename T> struct is_array_class : std::false_type {};

template <typename T, std::size_t N>
struct is_array_class<std::array<T, N>> : std::true_type {};

template <typename T>
constexpr bool is_array_class_v = is_array_class<T>::value;

template <typename T, typename = void> struct make_signed_upcast {
  using type = T;
};

template <typename T>
struct make_signed_upcast<T, std::enable_if_t<std::is_same_v<T, std::size_t>>> {
  using type = std::ptrdiff_t;
};

template <typename T>
struct make_signed_upcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint8_t> and
                        not std::is_same_v<T, std::size_t>>> {
  using type = std::int16_t;
};

template <typename T>
struct make_signed_upcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint16_t> and
                        not std::is_same_v<T, std::size_t>>> {
  using type = std::int32_t;
};

template <typename T>
struct make_signed_upcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint32_t> and
                        not std::is_same_v<T, std::size_t>>> {
  using type = std::int64_t;
};

template <typename T>
using make_signed_upcast_t = typename make_signed_upcast<T>::type;

template <typename T, typename = void> struct make_unsigned_downcast {
  using type = T;
};

template <typename T>
struct make_unsigned_downcast<
    T, std::enable_if_t<std::is_same_v<T, std::ptrdiff_t>>> {
  using type = std::size_t;
};

template <typename T>
struct make_unsigned_downcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint64_t> and
                        not std::is_same_v<T, std::ptrdiff_t>>> {
  using type = std::uint32_t;
};

template <typename T>
struct make_unsigned_downcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint32_t> and
                        not std::is_same_v<T, std::ptrdiff_t>>> {
  using type = std::uint16_t;
};

template <typename T>
struct make_unsigned_downcast<
    T, std::enable_if_t<std::is_same_v<T, std::uint16_t> and
                        not std::is_same_v<T, std::ptrdiff_t>>> {
  using type = std::uint8_t;
};

template <typename T>
using make_unsigned_downcast_t = typename make_unsigned_downcast<T>::type;

} // namespace nostd

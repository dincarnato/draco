#pragma once

#include "../../nostd/type_traits.hpp"

namespace het::detail {

// Forward declarations

struct vec2d_element_empty;
template <auto>
struct vec2d_element_fixed_size;
struct vec2d_element_dynamic_size;
template <typename Fun>
struct vec2d_element_dynamic_size_from_callable;
template <typename T>
struct vec2d_element_default_uninitialized;
template <typename T>
struct vec2d_element_default_construction_default;
template <typename T, typename... Values>
struct vec2d_element_default_construction;
template <typename... Ts>
struct vec2d_build_parts;

// End of forward declarations

template <typename T>
struct is_vec2d_empty : std::false_type {};

template <>
struct is_vec2d_empty<vec2d_element_empty> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_empty_v = is_vec2d_empty<T>::value;

template <typename T>
struct is_vec2d_fixed_size : std::false_type {};

template <std::size_t Size>
struct is_vec2d_fixed_size<vec2d_element_fixed_size<Size>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_fixed_size_v = is_vec2d_fixed_size<T>::value;

template <typename T>
struct is_vec2d_dynamic_size : std::false_type {};

template <>
struct is_vec2d_dynamic_size<vec2d_element_dynamic_size> : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size<vec2d_element_dynamic_size_from_callable<Fun>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_dynamic_size_v = is_vec2d_dynamic_size<T>::value;

template <typename T>
struct is_vec2d_size_type : std::false_type {};

template <std::size_t Size>
struct is_vec2d_size_type<vec2d_element_fixed_size<Size>> : std::true_type {};

template <>
struct is_vec2d_size_type<vec2d_element_dynamic_size> : std::true_type {};

template <typename Fun>
struct is_vec2d_size_type<vec2d_element_dynamic_size_from_callable<Fun>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_size_type_v = is_vec2d_size_type<T>::value;

template <typename T>
struct is_vec2d_element_default_uninitialized : std::false_type {};

template <typename SizeType>
struct is_vec2d_element_default_uninitialized<
    vec2d_element_default_uninitialized<SizeType>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_element_default_uninitialized_v =
    is_vec2d_element_default_uninitialized<T>::value;

template <typename T>
struct is_vec2d_element_default_construction_default : std::false_type {};

template <typename SizeType>
struct is_vec2d_element_default_construction_default<
    vec2d_element_default_construction_default<SizeType>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_element_default_construction_default_v =
    is_vec2d_element_default_construction_default<T>::value;

template <typename T>
struct is_vec2d_element_default_construction : std::false_type {};

template <typename SizeType, typename... Values>
struct is_vec2d_element_default_construction<
    vec2d_element_default_construction<SizeType, Values...>> : std::true_type {
};

template <typename T>
constexpr bool is_vec2d_element_default_construction_v =
    is_vec2d_element_default_construction<T>::value;

template <typename T>
struct is_vec2d_element_with_defaults : std::false_type {};

template <typename SizeType>
struct is_vec2d_element_with_defaults<
    vec2d_element_default_uninitialized<SizeType>> : std::true_type {};

template <typename SizeType>
struct is_vec2d_element_with_defaults<
    vec2d_element_default_construction_default<SizeType>> : std::true_type {};

template <typename SizeType, typename... Values>
struct is_vec2d_element_with_defaults<
    vec2d_element_default_construction<SizeType, Values...>> : std::true_type {
};

template <typename T>
constexpr bool is_vec2d_element_with_defaults_v =
    is_vec2d_element_with_defaults<T>::value;

template <typename T>
struct is_vec2d_fixed_size_part : std::false_type {};

template <std::size_t Size>
struct is_vec2d_fixed_size_part<vec2d_element_fixed_size<Size>>
    : std::true_type {};

template <std::size_t Size>
struct is_vec2d_fixed_size_part<
    vec2d_element_default_uninitialized<vec2d_element_fixed_size<Size>>>
    : std::true_type {};

template <std::size_t Size>
struct is_vec2d_fixed_size_part<
    vec2d_element_default_construction_default<vec2d_element_fixed_size<Size>>>
    : std::true_type {};

template <std::size_t Size, typename Value>
struct is_vec2d_fixed_size_part<
    vec2d_element_default_construction<vec2d_element_fixed_size<Size>, Value>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_fixed_size_part_v = is_vec2d_fixed_size_part<T>::value;

template <typename T>
struct is_vec2d_dynamic_size_part : std::false_type {};

template <>
struct is_vec2d_dynamic_size_part<vec2d_element_dynamic_size> : std::true_type {
};

template <>
struct is_vec2d_dynamic_size_part<
    vec2d_element_default_uninitialized<vec2d_element_dynamic_size>>
    : std::true_type {};

template <>
struct is_vec2d_dynamic_size_part<
    vec2d_element_default_construction_default<vec2d_element_dynamic_size>>
    : std::true_type {};

template <typename Value>
struct is_vec2d_dynamic_size_part<
    vec2d_element_default_construction<vec2d_element_dynamic_size, Value>>
    : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_part<vec2d_element_dynamic_size_from_callable<Fun>>
    : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_part<vec2d_element_default_uninitialized<
    vec2d_element_dynamic_size_from_callable<Fun>>> : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_part<vec2d_element_default_construction_default<
    vec2d_element_dynamic_size_from_callable<Fun>>> : std::true_type {};

template <typename Value, typename Fun>
struct is_vec2d_dynamic_size_part<vec2d_element_default_construction<
    vec2d_element_dynamic_size_from_callable<Fun>, Value>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_dynamic_size_part_v =
    is_vec2d_dynamic_size_part<T>::value;

template <typename T>
struct is_vec2d_size_part
    : std::bool_constant<is_vec2d_dynamic_size_part_v<T> or
                         is_vec2d_fixed_size_part_v<T>> {};

template <typename T>
constexpr bool is_vec2d_size_part_v = is_vec2d_size_part<T>::value;

template <typename T>
struct is_vec2d_dynamic_size_from_callable_part : std::false_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_from_callable_part<
    vec2d_element_dynamic_size_from_callable<Fun>> : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_from_callable_part<
    vec2d_element_default_uninitialized<
        vec2d_element_dynamic_size_from_callable<Fun>>> : std::true_type {};

template <typename Fun>
struct is_vec2d_dynamic_size_from_callable_part<
    vec2d_element_default_construction_default<
        vec2d_element_dynamic_size_from_callable<Fun>>> : std::true_type {};

template <typename Value, typename Fun>
struct is_vec2d_dynamic_size_from_callable_part<
    vec2d_element_default_construction<
        vec2d_element_dynamic_size_from_callable<Fun>, Value>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_dynamic_size_from_callable_part_v =
    is_vec2d_dynamic_size_from_callable_part<T>::value;

template <typename T>
struct is_vec2d_build_parts : std::false_type {};

template <>
struct is_vec2d_build_parts<void> : std::true_type {};

template <typename... Ts>
struct is_vec2d_build_parts<vec2d_build_parts<Ts...>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_build_parts_v = is_vec2d_build_parts<T>::value;

} // namespace het::detail

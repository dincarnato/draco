#pragma once

#include <tuple>
#include <type_traits>

namespace het::detail {

// Forward delarations
//
template <typename... Ts> struct vec2d_init_parts;

struct vec2d_rt_part_empty;
struct vec2d_rt_part_skip;

template <typename> struct vec2d_rt_part_with_lines;

template <typename> struct vec2d_rt_part_skip_first_lines;

template <typename> struct vec2d_rt_part_with_max_size;

template <typename> struct vec2d_rt_part_with_size;

template <typename, typename> struct vec2d_rt_part_with_fun_size;

template <typename, typename> struct vec2d_rt_part_transform;

struct vec2d_rt_part_construction_placeholder;

template <typename T> struct vec2d_rt_part_construct_default;

template <typename T> struct vec2d_rt_part_uninitialized;

template <typename T, typename Value> struct vec2d_rt_part_construct_copies;

template <typename T, typename... ArgsTuples> struct vec2d_rt_part_construct;

template <typename T, typename Iter> struct vec2d_rt_part_construct_from_iter;

template <typename T, typename Iter>
struct vec2d_rt_part_construct_from_args_iter;

template <typename SizeT, typename T> struct vec2d_rt_part_from_address;

struct vec2d_rt_part_unwinder;

// End of forward declarations

template <typename T> struct is_vec2d_init_parts : std::false_type {};

template <typename... Ts>
struct is_vec2d_init_parts<vec2d_init_parts<Ts...>> : std::true_type {};

template <typename... Ts>
constexpr bool is_vec2d_init_parts_v = is_vec2d_init_parts<Ts...>::value;

template <typename T> struct is_vec2d_rt_part_empty : std::false_type {};

template <>
struct is_vec2d_rt_part_empty<vec2d_rt_part_empty> : std::true_type {};

template <typename Fun>
struct is_vec2d_rt_part_empty<vec2d_rt_part_transform<vec2d_rt_part_empty, Fun>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_empty_v = is_vec2d_rt_part_empty<T>::value;

template <typename T> struct is_vec2d_rt_part_skip : std::false_type {};

template <>
struct is_vec2d_rt_part_skip<vec2d_rt_part_skip> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_skip_v = is_vec2d_rt_part_skip<T>::value;

template <typename T> struct is_vec2d_rt_part_with_lines : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_with_lines<vec2d_rt_part_with_lines<SizeT>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_with_lines_v =
    is_vec2d_rt_part_with_lines<T>::value;

template <typename T>
struct is_vec2d_rt_part_skip_first_lines : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_skip_first_lines<vec2d_rt_part_skip_first_lines<SizeT>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_skip_first_lines_v =
    is_vec2d_rt_part_skip_first_lines<T>::value;

template <typename T>
struct is_vec2d_rt_part_with_max_size : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_with_max_size<vec2d_rt_part_with_max_size<SizeT>>
    : std::true_type {};

template <typename SizeT, typename Fun>
struct is_vec2d_rt_part_with_max_size<
    vec2d_rt_part_transform<vec2d_rt_part_with_max_size<SizeT>, Fun>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_with_max_size_v =
    is_vec2d_rt_part_with_max_size<T>::value;

template <typename T> struct is_vec2d_rt_part_with_size : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_with_size<vec2d_rt_part_with_size<SizeT>>
    : std::true_type {};

template <typename SizeT, typename Fun>
struct is_vec2d_rt_part_with_size<
    vec2d_rt_part_transform<vec2d_rt_part_with_size<SizeT>, Fun>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_with_size_v =
    is_vec2d_rt_part_with_size<T>::value;

template <typename T>
struct is_vec2d_rt_part_with_fun_size : std::false_type {};

template <typename Fun, typename SizeT>
struct is_vec2d_rt_part_with_fun_size<vec2d_rt_part_with_fun_size<Fun, SizeT>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_with_fun_size_v =
    is_vec2d_rt_part_with_fun_size<T>::value;

template <typename T> struct is_vec2d_rt_part_transform : std::false_type {};

template <typename Fun, typename Base>
struct is_vec2d_rt_part_transform<vec2d_rt_part_transform<Fun, Base>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_transform_v =
    is_vec2d_rt_part_transform<T>::value;

template <typename T>
struct is_vec2d_rt_part_construction_placeholder : std::false_type {};

template <>
struct is_vec2d_rt_part_construction_placeholder<
    vec2d_rt_part_construction_placeholder> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construction_placeholder_v =
    is_vec2d_rt_part_construction_placeholder<T>::value;

template <typename T>
struct is_vec2d_rt_part_construct_default : std::false_type {};

template <typename T>
struct is_vec2d_rt_part_construct_default<vec2d_rt_part_construct_default<T>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construct_default_v =
    is_vec2d_rt_part_construct_default<T>::value;

template <typename T>
struct is_vec2d_rt_part_uninitialized : std::false_type {};

template <typename T>
struct is_vec2d_rt_part_uninitialized<vec2d_rt_part_uninitialized<T>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_uninitialized_v =
    is_vec2d_rt_part_uninitialized<T>::value;

template <typename T>
struct is_vec2d_rt_part_construct_copies : std::false_type {};

template <typename T, typename Value>
struct is_vec2d_rt_part_construct_copies<
    vec2d_rt_part_construct_copies<T, Value>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construct_copies_v =
    is_vec2d_rt_part_construct_copies<T>::value;

template <typename T> struct is_vec2d_rt_part_construct : std::false_type {};

template <typename T, typename... ArgsTuples>
struct is_vec2d_rt_part_construct<vec2d_rt_part_construct<T, ArgsTuples...>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construct_v =
    is_vec2d_rt_part_construct<T>::value;

template <typename T>
struct is_vec2d_rt_part_construct_from_iter : std::false_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_construct_from_iter<
    vec2d_rt_part_construct_from_iter<T, Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construct_from_iter_v =
    is_vec2d_rt_part_construct_from_iter<T>::value;

template <typename T>
struct is_vec2d_rt_part_construct_from_args_iter : std::false_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_construct_from_args_iter<
    vec2d_rt_part_construct_from_args_iter<T, Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_construct_from_args_iter_v =
    is_vec2d_rt_part_construct_from_args_iter<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_construct_base : std::false_type {};

template <>
struct is_vec2d_rt_part_generic_construct_base<vec2d_rt_part_empty>
    : std::true_type {};

template <>
struct is_vec2d_rt_part_generic_construct_base<
    vec2d_rt_part_construction_placeholder> : std::true_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_base<
    vec2d_rt_part_with_max_size<SizeT>> : std::true_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_base<vec2d_rt_part_with_size<SizeT>>
    : std::true_type {};

template <typename Fun, typename SizeT>
struct is_vec2d_rt_part_generic_construct_base<
    vec2d_rt_part_with_fun_size<Fun, SizeT>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_construct_base_v =
    is_vec2d_rt_part_generic_construct_base<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_construct : std::false_type {};

template <typename T>
struct is_vec2d_rt_part_generic_construct<vec2d_rt_part_construct_default<T>>
    : std::true_type {};

template <typename T>
struct is_vec2d_rt_part_generic_construct<vec2d_rt_part_uninitialized<T>>
    : std::true_type {};

template <typename T, typename Value>
struct is_vec2d_rt_part_generic_construct<
    vec2d_rt_part_construct_copies<T, Value>> : std::true_type {};

template <typename T, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_construct<
    vec2d_rt_part_construct<T, ArgsTuples...>> : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_construct<
    vec2d_rt_part_construct_from_iter<T, Iter>> : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_construct<
    vec2d_rt_part_construct_from_args_iter<T, Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_construct_v =
    is_vec2d_rt_part_generic_construct<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_construct_with_size : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_size<SizeT>>>
    : std::true_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_size<SizeT>>>
    : std::true_type {};

template <typename SizeT, typename Value>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_size<SizeT>, Value>>
    : std::true_type {};

template <typename SizeT, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_construct<vec2d_rt_part_with_size<SizeT>, ArgsTuples...>>
    : std::true_type {};

template <typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_construct_from_iter<vec2d_rt_part_with_size<SizeT>, Iter>>
    : std::true_type {};

template <typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_size<
    vec2d_rt_part_construct_from_args_iter<vec2d_rt_part_with_size<SizeT>,
                                           Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_construct_with_size_v =
    is_vec2d_rt_part_generic_construct_with_size<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_construct_with_fun_size : std::false_type {};

template <typename Fun, typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_fun_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_fun_size<Fun, SizeT>>>
    : std::true_type {};

template <typename Fun, typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_fun_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_fun_size<Fun, SizeT>>>
    : std::true_type {};

template <typename Fun, typename SizeT, typename Value>
struct is_vec2d_rt_part_generic_construct_with_fun_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_fun_size<Fun, SizeT>,
                                   Value>> : std::true_type {};

template <typename Fun, typename SizeT, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_construct_with_fun_size<vec2d_rt_part_construct<
    vec2d_rt_part_with_fun_size<Fun, SizeT>, ArgsTuples...>> : std::true_type {
};

template <typename Fun, typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_fun_size<
    vec2d_rt_part_construct_from_iter<vec2d_rt_part_with_fun_size<Fun, SizeT>,
                                      Iter>> : std::true_type {};

template <typename Fun, typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_fun_size<
    vec2d_rt_part_construct_from_args_iter<
        vec2d_rt_part_with_fun_size<Fun, SizeT>, Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_construct_with_fun_size_v =
    is_vec2d_rt_part_generic_construct_with_fun_size<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_construct_with_max_size : std::false_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_max_size<SizeT>>>
    : std::true_type {};

template <typename SizeT>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_max_size<SizeT>>>
    : std::true_type {};

template <typename SizeT, typename Value>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_max_size<SizeT>, Value>>
    : std::true_type {};

template <typename SizeT, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_construct<vec2d_rt_part_with_max_size<SizeT>, ArgsTuples...>>
    : std::true_type {};

template <typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_construct_from_iter<vec2d_rt_part_with_max_size<SizeT>, Iter>>
    : std::true_type {};

template <typename SizeT, typename Iter>
struct is_vec2d_rt_part_generic_construct_with_max_size<
    vec2d_rt_part_construct_from_args_iter<vec2d_rt_part_with_max_size<SizeT>,
                                           Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_construct_with_max_size_v =
    is_vec2d_rt_part_generic_construct_with_max_size<T>::value;

template <typename T> struct is_vec2d_rt_part_from_address : std::false_type {};

template <typename SizeT, typename T>
struct is_vec2d_rt_part_from_address<vec2d_rt_part_from_address<SizeT, T>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_from_address_v =
    is_vec2d_rt_part_from_address<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_with_size : std::false_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_size<vec2d_rt_part_with_size<T>>
    : std::true_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_size<T>>>
    : std::true_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_size<T>>> : std::true_type {
};

template <typename T, typename Value>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_size<T>, Value>>
    : std::true_type {};

template <typename T, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_construct<vec2d_rt_part_with_size<T>, ArgsTuples...>>
    : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_construct_from_iter<vec2d_rt_part_with_size<T>, Iter>>
    : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_size<
    vec2d_rt_part_construct_from_args_iter<vec2d_rt_part_with_size<T>, Iter>>
    : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_with_size_v =
    is_vec2d_rt_part_generic_with_size<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_with_fun_size : std::false_type {};

template <typename Fun, typename T>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_with_fun_size<Fun, T>> : std::true_type {};

template <typename Fun, typename T>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_fun_size<Fun, T>>>
    : std::true_type {};

template <typename Fun, typename T>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_fun_size<Fun, T>>>
    : std::true_type {};

template <typename Fun, typename T, typename Value>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_fun_size<Fun, T>, Value>>
    : std::true_type {};

template <typename Fun, typename T, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_construct<vec2d_rt_part_with_fun_size<Fun, T>, ArgsTuples...>>
    : std::true_type {};

template <typename Fun, typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_fun_size<vec2d_rt_part_construct_from_iter<
    vec2d_rt_part_with_fun_size<Fun, T>, Iter>> : std::true_type {};

template <typename Fun, typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_fun_size<
    vec2d_rt_part_construct_from_args_iter<vec2d_rt_part_with_fun_size<Fun, T>,
                                           Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_with_fun_size_v =
    is_vec2d_rt_part_generic_with_fun_size<T>::value;

template <typename T>
struct is_vec2d_rt_part_generic_with_max_size : std::false_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_max_size<vec2d_rt_part_with_max_size<T>>
    : std::true_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_construct_default<vec2d_rt_part_with_max_size<T>>>
    : std::true_type {};

template <typename T>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_uninitialized<vec2d_rt_part_with_max_size<T>>>
    : std::true_type {};

template <typename T, typename Value>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_construct_copies<vec2d_rt_part_with_max_size<T>, Value>>
    : std::true_type {};

template <typename T, typename... ArgsTuples>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_construct<vec2d_rt_part_with_max_size<T>, ArgsTuples...>>
    : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_construct_from_iter<vec2d_rt_part_with_max_size<T>, Iter>>
    : std::true_type {};

template <typename T, typename Iter>
struct is_vec2d_rt_part_generic_with_max_size<
    vec2d_rt_part_construct_from_args_iter<vec2d_rt_part_with_max_size<T>,
                                           Iter>> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_generic_with_max_size_v =
    is_vec2d_rt_part_generic_with_max_size<T>::value;

template <typename T> struct is_vec2d_rt_part_unwinder : std::false_type {};

template <>
struct is_vec2d_rt_part_unwinder<vec2d_rt_part_unwinder> : std::true_type {};

template <typename T>
constexpr bool is_vec2d_rt_part_unwinder_v =
    is_vec2d_rt_part_unwinder<T>::value;

} // namespace het::detail

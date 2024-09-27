#pragma once

#include "../../nostd/type_traits.hpp"
#include "../allocator_traits.hpp"
#include "vec2d_element_builder_rt_parts_traits.hpp"

#include <iterator>
#include <tuple>

namespace het::detail {

// Forward declarations

struct vec2d_rt_part_empty;

template <template <typename...> typename Parent, typename... Ts>
struct vec2d_build_parts_common;

template <template <typename, typename, typename> typename TemplatedBuilder,
          typename Alloc, typename BuildParts, typename InitParts>
struct vec2d_element_builder_rt_common;

template <typename Alloc, typename BuildParts,
          typename InitParts = vec2d_init_parts<vec2d_rt_part_empty>>
struct vec2d_element_builder_rt_build_at;

template <typename Alloc, typename BuildParts,
          typename InitParts = vec2d_init_parts<vec2d_rt_part_empty>>
struct vec2d_element_builder_rt_move_to;

template <typename Alloc, typename BuildParts,
          typename InitParts = vec2d_init_parts<vec2d_rt_part_empty>>
struct vec2d_element_builder_rt_destroy;

// End of forward declarations

template <typename... Ts>
struct vec2d_init_parts : vec2d_build_parts_common<vec2d_init_parts, Ts...> {
  using base_type = vec2d_build_parts_common<vec2d_init_parts, Ts...>;
  template <std::size_t Idx>
  using type = std::tuple_element_t<Idx, std::tuple<Ts...>>;
  static constexpr std::size_t arity = sizeof...(Ts);
  using last_type = type<arity - 1>;

  template <typename T>
  using replace_last_with_t =
      typename base_type::template replace_last_with_t<T>;

  template <typename T> using append_t = vec2d_init_parts<Ts..., T>;

  using base_type::base_type;
};

struct vec2d_rt_part_empty {};
struct vec2d_rt_part_skip {};

template <typename SizeT> struct vec2d_rt_part_with_lines {
  using size_type = SizeT;

  inline vec2d_rt_part_with_lines(const vec2d_rt_part_with_lines &) = default;
  inline vec2d_rt_part_with_lines(vec2d_rt_part_with_lines &&) = default;

  SizeT lines;
};

template <typename SizeT> struct vec2d_rt_part_skip_first_lines {
  using size_type = SizeT;

  template <typename _SizeT>
  inline vec2d_rt_part_skip_first_lines(
      const vec2d_rt_part_with_lines<SizeT> &lines, _SizeT &&skip) noexcept
      : lines(lines.lines), skip_lines(std::forward<_SizeT>(skip)) {}

  template <typename _SizeT>
  inline vec2d_rt_part_skip_first_lines(vec2d_rt_part_with_lines<SizeT> &&lines,
                                        _SizeT skip) noexcept
      : lines(std::move(lines.lines)), skip_lines(std::forward<_SizeT>(skip)) {}

  SizeT lines;
  SizeT skip_lines;
};

template <typename SizeT> struct vec2d_rt_part_with_max_size {
  using size_type = SizeT;

  inline vec2d_rt_part_with_max_size(const vec2d_rt_part_with_max_size &) =
      default;
  inline vec2d_rt_part_with_max_size(vec2d_rt_part_with_max_size &&) = default;

  SizeT max_size;
};

template <typename SizeT> struct vec2d_rt_part_with_size {
  using size_type = SizeT;

  inline vec2d_rt_part_with_size(const vec2d_rt_part_with_size &) = default;
  inline vec2d_rt_part_with_size(vec2d_rt_part_with_size &&) = default;

  SizeT size;
  SizeT max_size;
};

template <typename Fun, typename SizeT>
struct vec2d_rt_part_with_fun_size : Fun {
  using function_type = Fun;

  inline vec2d_rt_part_with_fun_size(const vec2d_rt_part_with_fun_size &) =
      default;
  inline vec2d_rt_part_with_fun_size(vec2d_rt_part_with_fun_size &&) = default;

  SizeT max_size;
};

template <typename Base, typename Fun>
struct vec2d_rt_part_transform : Base, Fun {
  using base_type = Base;
  using function_type = Fun;

  template <typename _Base, typename _Fun>
  inline vec2d_rt_part_transform(_Base &&base, _Fun &&fun) noexcept(
      noexcept(Fun(std::forward<_Fun>(fun))) and
      noexcept(Base(std::forward<_Base>(base))))
      : base_type(std::forward<_Base>(base)),
        function_type(std::forward<_Fun>(fun)) {}
};

template <typename Base, typename Fun>
vec2d_rt_part_transform(Base &&, Fun &&)
    -> vec2d_rt_part_transform<std::decay_t<Base>, std::decay_t<Fun>>;

struct vec2d_rt_part_construction_placeholder {};

template <typename T> struct vec2d_rt_part_construct_default : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);
  template <typename _T>
  inline explicit vec2d_rt_part_construct_default(_T &&t) noexcept(
      noexcept(T(std::forward<_T>(t))))
      : T(std::forward<_T>(t)) {}
};

template <typename T>
vec2d_rt_part_construct_default(T &&)
    -> vec2d_rt_part_construct_default<std::decay_t<T>>;

template <typename T> struct vec2d_rt_part_uninitialized : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);

  template <typename _T>
  inline explicit vec2d_rt_part_uninitialized(_T &&t) noexcept(
      noexcept(T(std::forward<_T>(t))))
      : T(std::forward<_T>(t)) {}
};

template <typename T>
vec2d_rt_part_uninitialized(T &&)
    -> vec2d_rt_part_uninitialized<std::decay_t<T>>;

template <typename T, typename Value>
struct vec2d_rt_part_construct_copies : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);

  template <typename _Value>
  inline vec2d_rt_part_construct_copies(T &&t, _Value &&value) noexcept(
      noexcept(T(std::move(t))) and
      noexcept(Value(std::forward<_Value>(value))))
      : T(std::move(t)), value(std::forward<_Value>(value)) {}

  Value value;
};

template <typename T, typename Value>
vec2d_rt_part_construct_copies(T &&, Value &&)
    -> vec2d_rt_part_construct_copies<std::decay_t<T>, std::decay_t<Value>>;

template <typename T, typename... ArgsTuples>
struct vec2d_rt_part_construct : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);
  static_assert((nostd::is_tuple_v<ArgsTuples> && ...));
  constexpr static std::size_t args_tuples_arity = sizeof...(ArgsTuples);

  template <typename... _ArgsTuples>
  inline vec2d_rt_part_construct(T &&t, _ArgsTuples &&...args_tuples) noexcept(
      noexcept(T(std::move(t))) and
      noexcept(
          std::tuple<ArgsTuples...>(std::forward<_ArgsTuples>(args_tuples)...)))
      : T(std::move(t)),
        args_tuples(std::forward<_ArgsTuples>(args_tuples)...) {}

  template <typename... Args>
  inline auto new_appended(Args &&...args) noexcept(
      noexcept(new_appended_impl(std::make_tuple(std::forward<Args>(args)...),
                                 std::index_sequence_for<ArgsTuples...>()))) {
    return new_appended_impl(std::make_tuple(std::forward<Args>(args)...),
                             std::index_sequence_for<ArgsTuples...>());
  }

  template <typename... NewArgsTuples>
  inline auto new_appending_tuples(NewArgsTuples &&...new_args_tuples) noexcept(
      noexcept(new_appending_tuples_impl(
          std::index_sequence_for<ArgsTuples...>(),
          std::forward<NewArgsTuples>(new_args_tuples)...))) {
    static_assert((nostd::is_tuple_v<std::decay_t<NewArgsTuples>> && ...));

    return new_appending_tuples_impl(
        std::index_sequence_for<ArgsTuples...>(),
        std::forward<NewArgsTuples>(new_args_tuples)...);
  }

  std::tuple<ArgsTuples...> args_tuples;

private:
  template <typename NewTuple, std::size_t... Idx>
  inline auto new_appended_impl(
      NewTuple &&newTupleArgs,
      std::index_sequence<
          Idx...>) noexcept(std::is_nothrow_move_constructible_v<T> and
                            (std::is_nothrow_move_constructible_v<ArgsTuples> &&
                             ...) and
                            std::is_nothrow_move_constructible_v<
                                std::decay_t<NewTuple>>) {
    return vec2d_rt_part_construct<T, ArgsTuples..., std::decay_t<NewTuple>>(
        std::move(static_cast<T &>(*this)),
        std::move(std::get<Idx>(args_tuples))..., std::move(newTupleArgs));
  }

  template <typename... NewArgsTuples, std::size_t... Idx>
  inline auto new_appending_tuples_impl(
      std::index_sequence<Idx...>,
      NewArgsTuples
          &&...new_args_tuples) noexcept(std::
                                             is_nothrow_move_constructible_v<
                                                 T> and
                                         (std::is_nothrow_move_constructible_v<
                                              ArgsTuples> &&
                                          ...) and
                                         (std::is_nothrow_move_constructible_v<
                                              NewArgsTuples> &&
                                          ...)) {
    return vec2d_rt_part_construct<T, ArgsTuples...,
                                   std::decay_t<NewArgsTuples>...>(
        std::move(static_cast<T &>(*this)),
        std::move(std::get<Idx>(args_tuples))...,
        std::forward<NewArgsTuples>(new_args_tuples)...);
  }
};

template <typename T, typename... ArgsTuples>
vec2d_rt_part_construct(T &&, ArgsTuples &&...)
    -> vec2d_rt_part_construct<std::decay_t<T>, std::decay_t<ArgsTuples>...>;

template <typename T, typename Iter>
struct vec2d_rt_part_construct_from_iter : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);

  template <typename _Iter>
  inline vec2d_rt_part_construct_from_iter(T &&t, _Iter &&iter) noexcept(
      noexcept(T(std::move(t))) and noexcept(Iter(std::forward<_Iter>(iter))))
      : T(std::move(t)), iter(std::forward<_Iter>(iter)) {}

  Iter iter;
};

template <typename T, typename Iter>
struct vec2d_rt_part_construct_from_args_iter : T {
  static_assert(is_vec2d_rt_part_generic_construct_base_v<T>);
  static_assert(
      nostd::is_tuple_v<typename std::iterator_traits<Iter>::value_type>);

  template <typename _Iter>
  inline vec2d_rt_part_construct_from_args_iter(T &&t, _Iter &&iter) noexcept(
      noexcept(T(std::move(t))) and noexcept(std::forward<_Iter>(iter)))
      : T(std::move(t)), iter(std::forward<_Iter>(iter)) {}

  Iter iter;
};

template <typename LinesT, typename T>
struct vec2d_rt_part_from_address : LinesT {
  static_assert(is_vec2d_rt_part_with_lines_v<LinesT> or
                is_vec2d_rt_part_skip_first_lines_v<LinesT>);
  using base_type = LinesT;
  using size_type = typename base_type::size_type;

  template <typename _LinesT>
  inline vec2d_rt_part_from_address(_LinesT &&lines, T *from_address) noexcept(
      std::is_nothrow_move_constructible_v<size_type>)
      : base_type(std::forward<_LinesT>(lines)), from_address(from_address) {}

  T *from_address;
};

template <typename LinesT, typename T>
vec2d_rt_part_from_address(LinesT &&lines, T *from_address)
    -> vec2d_rt_part_from_address<std::decay_t<LinesT>, std::decay_t<T>>;

struct vec2d_rt_part_unwinder {};

} // namespace het::detail

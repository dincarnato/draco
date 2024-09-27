#pragma once

#include "../allocator_traits.hpp"
#include "vec2d_element_builder_rt_common.hpp"
#include "vec2d_element_builder_rt_construction_base.hpp"
#include "vec2d_element_builder_rt_parts.hpp"
#include "vec2d_element_builder_rt_parts_traits.hpp"

namespace het::detail {

template <typename Alloc, typename BuildParts>
vec2d_element_builder_rt_build_at(
    const BuildParts &build_parts,
    typename allocator_traits<Alloc>::first_pointer,
    Alloc &) -> vec2d_element_builder_rt_build_at<std::decay_t<Alloc>,
                                                  std::decay_t<BuildParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
vec2d_element_builder_rt_build_at(
    const BuildParts &build_parts,
    typename allocator_traits<Alloc>::first_pointer, Alloc &, InitParts &&)
    -> vec2d_element_builder_rt_build_at<
        std::decay_t<Alloc>, std::decay_t<BuildParts>, std::decay_t<InitParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
struct vec2d_element_builder_rt_build_at
    : vec2d_element_builder_rt_common<vec2d_element_builder_rt_build_at, Alloc,
                                      BuildParts, InitParts>,
      vec2d_element_builder_rt_construction_base<
          vec2d_element_builder_rt_build_at, Alloc, BuildParts, InitParts> {
  using base_type =
      vec2d_element_builder_rt_common<vec2d_element_builder_rt_build_at, Alloc,
                                      BuildParts, InitParts>;
  using construction_base_type = vec2d_element_builder_rt_construction_base<
      vec2d_element_builder_rt_build_at, Alloc, BuildParts, InitParts>;
  friend vec2d_element_builder_rt_construction_base<
      vec2d_element_builder_rt_build_at, Alloc, BuildParts, InitParts>;

  using base_type::base_type;
  using init_part_type = typename base_type::init_parts_type;
  using allocator_type = typename base_type::allocator_type;
  using build_parts_type = typename base_type::build_parts_type;
  using init_parts_type = typename base_type::init_parts_type;
  using base_address_type = typename base_type::base_address_type;

  using base_type::address;
  using base_type::alloc;
  using base_type::build_parts;
  using base_type::init_parts;

  using base_type::destroy_element;
  using base_type::new_appended;
  using base_type::new_replaced_last;
  using base_type::unwind_destroy_line;
  using base_type::unwind_destroy_line_partial;

  template <std::size_t Index>
  static constexpr std::size_t get_element_index_for_init() noexcept {
    static_assert(Index >= 1);
    return Index - 1;
  }

  static constexpr bool is_lines_init() noexcept {
    return init_parts_type::arity == 1;
  }

  [[nodiscard]] inline auto
  next() noexcept(noexcept(this->new_appended(vec2d_rt_part_empty()))) {
    static_assert(init_parts_type::arity >= 2,
                  "next() must be called after construct_default(), "
                  "uninitialized(), construct_copies(), construct_from_iter(), "
                  "construct_from_args_iter() or, only in case of parts built "
                  "with fixed_size(), construct() or construct_all()");

    using last_init_part_type = typename init_parts_type::last_type;
    static_assert(is_vec2d_rt_part_generic_construct_v<last_init_part_type>,
                  "next() must be called after construct_default(), "
                  "uninitialized(), construct_copies(), construct_from_iter(), "
                  "construct_from_args_iter() or, only in case of parts built "
                  "with fixed_size(), construct() or construct_all()");

    if constexpr (is_vec2d_rt_part_construct_v<last_init_part_type>) {
      static constexpr std::size_t element_index = init_parts_type::arity - 2;
      using last_build_part_type =
          typename build_parts_type::template type<element_index>;

      static_assert(is_vec2d_fixed_size_part_v<last_build_part_type>);
      static_assert(
          last_build_part_type::size == last_init_part_type::args_tuples_arity,
          "next() is called after construct(), but not all the elements "
          "specified by fixed_size().set_size() have been constructed. Use "
          "more construct() calls or use construct_all() to specify the "
          "construction of the remaining elements");
    }

    static_assert(
        init_parts_type::arity <= build_parts_type::arity,
        "next() cannot be called on the last element. Use done() instead");

    return new_appended(vec2d_rt_part_empty());
  }

  inline auto done() {
    static_assert(init_parts_type::arity >= 2);

    using last_init_part_type = typename init_parts_type::last_type;
    static_assert(is_vec2d_rt_part_generic_construct_v<last_init_part_type>,
                  "done() must be called after construct_default(), "
                  "uninitialized(), construct_copies(), construct_from_iter(), "
                  "construct_from_args_iter() or, only in case of parts built "
                  "with fixed_size(), construct() or construct_all()");

    if constexpr (is_vec2d_rt_part_construct_v<last_init_part_type>) {
      static constexpr std::size_t element_index = init_parts_type::arity - 2;
      using last_build_part_type =
          typename build_parts_type::template type<element_index>;

      static_assert(is_vec2d_fixed_size_part_v<last_build_part_type>);
      static_assert(
          last_build_part_type::size == last_init_part_type::args_tuples_arity,
          "done() is called after construct(), but not all the elements "
          "specified by fixed_size().set_size() have been constructed. Use "
          "more construct() calls or use construct_all() to specify the "
          "construction of the remaining elements");
    }

    static_assert(
        init_parts_type::arity == build_parts_type::arity + 1,
        "done() can be called only on the last element. Use next() instead");

    perform_build(std::make_index_sequence<build_parts_type::arity>());
    return new_appended(vec2d_rt_part_unwinder());
  }

  inline void unwind() noexcept {
    static_assert(
        is_vec2d_rt_part_unwinder_v<typename init_parts_type::last_type>,
        "unwind() can be called only after the last done()");

    unwind_impl(std::make_index_sequence<build_parts_type::arity>());
  }

private:
  template <std::size_t... Idx>
  inline void perform_build(std::index_sequence<Idx...>) const {
    using init_type = typename init_parts_type::template type<0>;
    static_assert(is_vec2d_rt_part_with_lines_v<init_type> or
                  is_vec2d_rt_part_skip_first_lines_v<init_type>);
    using lines_type = typename init_type::size_type;
    auto &&lines = std::get<0>(init_parts).lines;

    const auto element_max_sizes =
        std::tuple(base_type::template get_init_part_max_size<Idx + 1>()...);

    const auto zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<init_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<init_type>);
        return std::get<0>(init_parts).skip_lines;
      }
    }();

    for (lines_type line_index = zero_line; line_index < lines; ++line_index) {
      std::size_t constructed_elements = 0;
      try {
        (std::apply(
             [&](auto &&...sizes) {
               construction_base_type::
                   template perform_all_elements_build_from_init<Idx + 1>(
                       line_index, std::forward<decltype(sizes)>(sizes)...);
               ++constructed_elements;
             },
             element_max_sizes),
         ...);
      } catch (...) {
        std::apply(
            [&](auto &&...sizes) {
              unwind_destroy_line_partial(
                  line_index, constructed_elements,
                  nostd::make_index_sequence_rev<sizeof...(Idx)>(),
                  std::forward<decltype(sizes)>(sizes)...);
            },
            element_max_sizes);

        while (line_index > zero_line) {
          std::apply(
              [&](auto &&...sizes) {
                unwind_destroy_line(
                    --line_index,
                    nostd::make_index_sequence_rev<sizeof...(Idx)>(),
                    std::forward<decltype(sizes)>(sizes)...);
              },
              element_max_sizes);
        }

        throw;
      }
    }
  }

  template <std::size_t... Idx>
  inline void unwind_impl(std::index_sequence<Idx...>) noexcept {
    using init_type = typename init_parts_type::template type<0>;
    using lines_type = typename init_type::size_type;
    const auto element_max_sizes =
        std::tuple(base_type::template get_init_part_max_size<Idx + 1>()...);

    const auto zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<init_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<init_type>);
        return std::get<0>(init_parts).skip_lines;
      }
    }();

    lines_type line_index = std::get<0>(init_parts).lines;
    while (line_index > zero_line) {
      std::apply(
          [&](auto &&...sizes) {
            unwind_destroy_line(
                --line_index,
                nostd::make_index_sequence_rev<build_parts_type::arity>(),
                std::forward<decltype(sizes)>(sizes)...);
          },
          element_max_sizes);
    }
  }
};

} // namespace het::detail

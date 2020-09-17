#pragma once

#include "vec2d_element_builder_rt_common.hpp"

namespace het::detail {

template <typename Alloc, typename BuildParts>
vec2d_element_builder_rt_destroy(
    const BuildParts& build_parts,
    typename allocator_traits<Alloc>::first_pointer, Alloc&)
    ->vec2d_element_builder_rt_destroy<std::decay_t<Alloc>,
                                       std::decay_t<BuildParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
vec2d_element_builder_rt_destroy(
    const BuildParts& build_parts,
    typename allocator_traits<Alloc>::first_pointer, Alloc&, InitParts &&)
    ->vec2d_element_builder_rt_destroy<
        std::decay_t<Alloc>, std::decay_t<BuildParts>, std::decay_t<InitParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
struct [[nodiscard]] vec2d_element_builder_rt_destroy
    : vec2d_element_builder_rt_common<vec2d_element_builder_rt_destroy, Alloc,
                                      BuildParts, InitParts> {
  using base_type =
      vec2d_element_builder_rt_common<vec2d_element_builder_rt_destroy, Alloc,
                                      BuildParts, InitParts>;
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
  using base_type::unwind_destroy_element;
  using base_type::unwind_destroy_line;

  template <std::size_t Index>
  static constexpr std::size_t get_element_index_for_init() noexcept {
    static_assert(Index >= 1);
    return Index - 1;
  }

  static constexpr bool is_lines_init() noexcept {
    return init_parts_type::arity == 1;
  }

  inline auto next() noexcept(
      noexcept(this->new_appended(vec2d_rt_part_empty()))) {
    static_assert(init_parts_type::arity >= 2);
    static constexpr std::size_t element_index = init_parts_type::arity - 2;

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(is_vec2d_rt_part_with_max_size_v<
                          typename init_parts_type::last_type> or
                          is_vec2d_rt_part_with_fun_size_v<
                              typename init_parts_type::last_type>,
                      "next() must be called after with_max_size() for parts "
                      "built with dynamic_size() but without "
                      "size_from_callable() or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_with_size_v<typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "next() must be called after with_size() for parts built with "
            "dynamic_size() and size_from_callable() or after with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "next() must be called after from() or next() for "
          "static sized parts");
    }

    static_assert(
        init_parts_type::arity <= build_parts_type::arity,
        "next() cannot be called on the last element. Use done() instead");

    return new_appended(vec2d_rt_part_empty());
  }

  inline auto done() {
    static_assert(init_parts_type::arity >= 2);
    static constexpr std::size_t element_index = init_parts_type::arity - 2;

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(is_vec2d_rt_part_with_max_size_v<
                          typename init_parts_type::last_type> or
                          is_vec2d_rt_part_with_fun_size_v<
                              typename init_parts_type::last_type>,
                      "done() must be called after with_max_size() for parts "
                      "built with dynamic_size() but without "
                      "size_from_callable() or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_with_size_v<typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "done() must be called after with_size() for parts built with "
            "dynamic_size() and size_from_callable() or after with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "done() must be called after from() or next() for "
          "static sized parts");
    }

    static_assert(
        init_parts_type::arity == build_parts_type::arity + 1,
        "done() can be called only on the last element. Use next() instead");

    perform_destroy(std::make_index_sequence<build_parts_type::arity>());
  }

private:
  template <std::size_t... Idx>
  inline void perform_destroy(std::index_sequence<Idx...>) noexcept {
    using init_type = typename init_parts_type::template type<0>;
    static_assert(is_vec2d_rt_part_with_lines_v<init_type> or
                  is_vec2d_rt_part_skip_first_lines_v<init_type>);
    using lines_type = typename init_type::size_type;

    const auto element_max_sizes =
        std::tuple(base_type::template get_init_part_max_size<Idx + 1>()...);

    auto line_index = std::get<0>(init_parts).lines;
    const auto zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<init_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<init_type>);
        return std::get<0>(init_parts).skip_lines;
      }
    }();

    while (line_index > zero_line) {
      std::apply(
          [&, this](auto&&... sizes) {
            unwind_destroy_line(
                --line_index, nostd::make_index_sequence_rev<sizeof...(Idx)>(),
                std::forward<decltype(sizes)>(sizes)...);
          },
          element_max_sizes);
    }
  }
};

} // namespace het::detail

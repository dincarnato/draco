#pragma once

#include "../allocator.hpp"
#include "../allocator_traits.hpp"
#include "vec2d_element_builder_ct_parts.hpp"
#include "vec2d_element_builder_ct_parts_traits.hpp"
#include "vec2d_element_builder_rt_parts.hpp"

#include "../../nostd/tuple.hpp"
#include "../../nostd/type_traits.hpp"
#include "../../nostd/utility.hpp"

#include <boost/callable_traits.hpp>

#include <cstdlib>
#include <iostream>
#include <optional>
#include <tuple>
#include <variant>

namespace het::detail {

template <typename Alloc, typename BuildParts, typename Idx>
struct vec2d_element_builder_ct_base {
  using allocator_type = Alloc;
  using build_parts_type = BuildParts;

  constexpr vec2d_element_builder_ct_base() = default;
  constexpr explicit vec2d_element_builder_ct_base(
      const build_parts_type& build_parts)
      : build_parts(build_parts) {}
  constexpr explicit vec2d_element_builder_ct_base(
      build_parts_type&& build_parts)
      : build_parts(std::move(build_parts)) {}

protected:
  build_parts_type build_parts;
};

template <typename Alloc, typename BuildParts, typename Idx>
struct vec2d_element_builder_ct;

template <typename Alloc, typename BuildParts, typename Idx>
struct vec2d_element_builder_rt;

template <>
struct vec2d_element_builder_ct<void, void,
                                std::integral_constant<std::size_t, 0>> {
  using build_parts_type = void;

  template <typename Alloc>
  constexpr auto
  set_allocator_type() const {
    return vec2d_element_builder_ct<Alloc,
                                    vec2d_build_parts<vec2d_element_empty>,
                                    std::integral_constant<std::size_t, 0>>{};
  }
};

template <typename Alloc, typename BuildParts, typename Idx>
struct vec2d_element_builder_ct
    : vec2d_element_builder_ct_base<Alloc, BuildParts, Idx> {
  using build_parts_type = BuildParts;
  using base_type = vec2d_element_builder_ct_base<Alloc, BuildParts, Idx>;
  static constexpr std::size_t index_value = Idx::value;

  using base_type::base_type;
  using base_type::build_parts;

  struct incomplete_fixed_size {
    constexpr incomplete_fixed_size(vec2d_element_builder_ct builder)
        : builder(std::move(builder)) {}

    template <std::size_t Sz>
    constexpr auto
    set_size() {
      static_assert(
          is_vec2d_empty_v<typename build_parts_type::last_type>,
          "fixed_size must be called after set_allocator_type() or next()");

      using new_build_parts_type =
          typename build_parts_type::template replace_last_with_t<
              vec2d_element_fixed_size<Sz>>;
      using new_builder_type =
          vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

      return new_builder_type(builder.build_parts.replace_last_with(
          vec2d_element_fixed_size<Sz>()));
    }

  private:
    vec2d_element_builder_ct builder;
  };

  constexpr auto
  fixed_size() {
    static_assert(
        is_vec2d_empty_v<typename build_parts_type::last_type>,
        "fixed_size must be called after set_allocator_type() or next()");
    return incomplete_fixed_size(std::move(*this));
  }

  constexpr auto
  dynamic_size() {
    static_assert(
        is_vec2d_empty_v<typename build_parts_type::last_type>,
        "dynamic_size must be called after set_allocator_type() or next()");

    using new_build_parts_type =
        typename build_parts_type::template replace_last_with_t<
            vec2d_element_dynamic_size>;
    using new_builder_type =
        vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

    return new_builder_type(
        build_parts.replace_last_with(vec2d_element_dynamic_size()));
  }

  template <typename Fun>
  constexpr auto
  size_from_callable(Fun&& fun) {
    static_assert(
        is_vec2d_dynamic_size_v<typename build_parts_type::last_type> and
            not is_vec2d_dynamic_size_from_callable_part_v<
                typename build_parts_type::last_type>,
        "size_from_callable() must be called after dynamic_size()");
    namespace ct = boost::callable_traits;
    static_assert(
        (std::tuple_size_v<ct::args_t<std::decay_t<Fun>>> <= index_value),
        "size_from_callable() must be called with an invocable "
        "object (function, lambda or struct with operator()), that "
        "accept a number of pointers less or equal than the number "
        "of previous constructed elements");

    using new_part_type =
        vec2d_element_dynamic_size_from_callable<std::decay_t<Fun>>;
    using new_build_parts_type =
        typename build_parts_type::template replace_last_with_t<new_part_type>;
    using new_builder_type =
        vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

    return new_builder_type(
        build_parts.replace_last_with(new_part_type(std::forward<Fun>(fun))));
  }

  constexpr auto
  default_construction_default() {
    static_assert(is_vec2d_size_part_v<typename build_parts_type::last_type>,
                  "default_construction_default() must be called after "
                  "fixed_size(), dynamic_size() or size_from_callable()");

    using new_part_type = vec2d_element_default_construction_default<
        typename build_parts_type::last_type>;
    using new_build_parts_type =
        typename build_parts_type::template replace_last_with_t<new_part_type>;
    using new_builder_type =
        vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

    return new_builder_type(build_parts.replace_last_with(
        new_part_type(std::move(std::get<Idx::value>(build_parts)))));
  }

  template <typename... Args>
  constexpr auto
  default_construction(Args&&... args) {
    static_assert(is_vec2d_size_part_v<typename build_parts_type::last_type>,
                  "default_construction() must be called after "
                  "fixed_size(), dynamic_size() or size_from_callable()");

    using new_part_type =
        vec2d_element_default_construction<typename build_parts_type::last_type,
                                           std::decay_t<Args>...>;
    using new_build_parts_type =
        typename build_parts_type::template replace_last_with_t<new_part_type>;
    using new_builder_type =
        vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

    return new_builder_type(build_parts.replace_last_with(
        new_part_type(build_parts.take_last(), std::forward<Args>(args)...)));
  }

  constexpr auto
  default_uninitialized() {
    static_assert(is_vec2d_size_part_v<typename build_parts_type::last_type>,
                  "default_uninitialized() must be called after "
                  "fixed_size(), dynamic_size() or size_from_callable()");
    static_assert(
        std::is_trivially_destructible_v<
            typename allocator_traits<Alloc>::template value_type<Idx::value>>,
        "the type must be trivially destructible in order to leave "
        "it uninitialized");

    using new_part_type = vec2d_element_default_uninitialized<
        typename build_parts_type::last_type>;
    using new_build_parts_type =
        typename build_parts_type::template replace_last_with_t<new_part_type>;
    using new_builder_type =
        vec2d_element_builder_ct<Alloc, new_build_parts_type, Idx>;

    return new_builder_type(
        build_parts.replace_last_with(new_part_type(build_parts.take_last())));
  }

  constexpr auto
  next() {
    static_assert(
        is_vec2d_element_with_defaults_v<typename build_parts_type::last_type>,
        "next() must be called after default_uninitialized(), "
        "default_construction_default() or default_construction()");
    static_assert(
        Idx::value < allocator_traits<Alloc>::args_arity::value - 1,
        "next() cannot be called on the last element. Use done() instead");

    using new_build_parts_type =
        typename build_parts_type::template append_t<vec2d_element_empty>;
    using new_builder_type = vec2d_element_builder_ct<
        Alloc, new_build_parts_type,
        std::integral_constant<std::size_t, Idx::value + 1>>;

    return new_builder_type(build_parts.append(vec2d_element_empty()));
  }

  constexpr auto
  done() {
    static_assert(
        is_vec2d_element_with_defaults_v<typename build_parts_type::last_type>,
        "done() must be called after default_uninitialized(), "
        "default_construction_default() or default_construction()");
    static_assert(
        Idx::value == allocator_traits<Alloc>::args_arity::value - 1,
        "done() can be called only on the last element. Use next() instead");

    using new_builder_type = vec2d_element_builder_ct<
        Alloc, build_parts_type,
        std::integral_constant<std::size_t, Idx::value + 1>>;

    return new_builder_type(std::move(build_parts));
  }
};

template <typename Alloc, typename BuildParts>
struct vec2d_element_builder_ct<
    Alloc, BuildParts,
    std::integral_constant<std::size_t,
                           allocator_traits<Alloc>::args_arity::value>>
    : vec2d_element_builder_ct_base<
          Alloc, BuildParts,
          std::integral_constant<std::size_t,
                                 allocator_traits<Alloc>::args_arity::value>> {
  using base_type = vec2d_element_builder_ct_base<
      Alloc, BuildParts,
      std::integral_constant<std::size_t,
                             allocator_traits<Alloc>::args_arity::value>>;

  using base_type::base_type;
  using build_parts_type = BuildParts;
  static constexpr std::size_t arity =
      allocator_traits<Alloc>::args_arity::value;
  template <std::size_t Idx>
  using build_part_type = typename build_parts_type::template type<Idx>;

  inline auto
  build_at(typename allocator_traits<Alloc>::first_pointer address,
           Alloc& alloc) const
      noexcept(noexcept(vec2d_element_builder_rt_build_at(
          base_type::build_parts, address, alloc))) {
    return vec2d_element_builder_rt_build_at(base_type::build_parts, address,
                                             alloc);
  }

  inline auto
  move_to(typename allocator_traits<Alloc>::first_pointer address,
          Alloc& alloc) const
      noexcept(noexcept(vec2d_element_builder_rt_move_to(base_type::build_parts,
                                                         address, alloc))) {
    return vec2d_element_builder_rt_move_to(base_type::build_parts, address,
                                            alloc);
  }

  inline auto
  destroy(typename allocator_traits<Alloc>::first_pointer address,
          Alloc& alloc) const
      noexcept(noexcept(vec2d_element_builder_rt_destroy(base_type::build_parts,
                                                         address, alloc))) {
    return vec2d_element_builder_rt_destroy(base_type::build_parts, address,
                                            alloc);
  }
};

} // namespace het::detail


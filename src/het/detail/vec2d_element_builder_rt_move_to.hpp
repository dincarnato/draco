#pragma once

#include "vec2d_element_builder_rt_common.hpp"
#include "vec2d_element_builder_rt_construction_base.hpp"
#include "../../nostd/utility.hpp"

#include <boost/callable_traits.hpp>

namespace het::detail {

template <typename Alloc, typename BuildParts>
vec2d_element_builder_rt_move_to(
    const BuildParts& build_parts,
    typename allocator_traits<Alloc>::first_pointer, Alloc&)
    ->vec2d_element_builder_rt_move_to<std::decay_t<Alloc>,
                                       std::decay_t<BuildParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
vec2d_element_builder_rt_move_to(
    const BuildParts& build_parts,
    typename allocator_traits<Alloc>::first_pointer, Alloc&, InitParts &&)
    ->vec2d_element_builder_rt_move_to<
        std::decay_t<Alloc>, std::decay_t<BuildParts>, std::decay_t<InitParts>>;

template <typename Alloc, typename BuildParts, typename InitParts>
struct vec2d_element_builder_rt_move_to
    : vec2d_element_builder_rt_common<vec2d_element_builder_rt_move_to, Alloc,
                                      BuildParts, InitParts>,
      vec2d_element_builder_rt_construction_base<
          vec2d_element_builder_rt_move_to, Alloc, BuildParts, InitParts> {
  using base_type =
      vec2d_element_builder_rt_common<vec2d_element_builder_rt_move_to, Alloc,
                                      BuildParts, InitParts>;
  using construction_base_type = vec2d_element_builder_rt_construction_base<
      vec2d_element_builder_rt_move_to, Alloc, BuildParts, InitParts>;
  friend vec2d_element_builder_rt_construction_base<
      vec2d_element_builder_rt_move_to, Alloc, BuildParts, InitParts>;

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
  using base_type::get_element_index;
  using base_type::new_appended;
  using base_type::new_replaced_last;
  using base_type::unwind_destroy_line;
  using base_type::unwind_destroy_line_partial;

  template <std::size_t Index>
  static constexpr std::size_t
  get_element_index_for_init() noexcept {
    static_assert(Index >= 1);

    if constexpr (Index <= build_parts_type::arity)
      return Index - 1;
    else if constexpr (Index <= build_parts_type::arity * 2)
      return Index - build_parts_type::arity - 1;
    else {
      static_assert(Index > build_parts_type::arity * 2 + 1 and
                    Index <= build_parts_type::arity * 3 + 1);
      return Index - build_parts_type::arity * 2 - 2;
    }
  }

  static constexpr bool
  is_lines_init() noexcept {
    return init_parts_type::arity == 1 or
           init_parts_type::arity == build_parts_type::arity * 2 + 2;
  }

  [[nodiscard]] inline auto
  from(typename allocator_traits<Alloc>::template pointer<0> address) noexcept(
      noexcept(this->new_replaced_last(vec2d_rt_part_from_address(
          std::declval<vec2d_rt_part_with_lines<
              typename allocator_traits<Alloc>::lines_size_type>&&>(),
          address)))) {

    static_assert(init_parts_type::arity == build_parts_type::arity * 2 + 2 and
                      (is_vec2d_rt_part_with_lines_v<
                           typename init_parts_type::last_type> or
                       is_vec2d_rt_part_skip_first_lines_v<
                           typename init_parts_type::last_type>),
                  "from() must be called after done().with_lines()");

    return new_replaced_last(
        vec2d_rt_part_from_address(init_parts.take_last(), address));
  }

  [[nodiscard]] inline auto
  remaining() noexcept(
      noexcept(this->new_appended(vec2d_rt_part_construction_placeholder()))) {
    static_assert(init_parts_type::arity >= 2);
    static_assert(init_parts_type::arity == build_parts_type::arity + 1,
                  "remaining() must be called on the last element after "
                  "[copy|move]_to() and before from()");
    static constexpr std::size_t element_index = get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(is_vec2d_rt_part_with_max_size_v<
                          typename init_parts_type::last_type> or
                          is_vec2d_rt_part_with_fun_size_v<
                              typename init_parts_type::last_type>,
                      "remaining() must be called after with_max_size() for "
                      "parts built with dynamic_size() but without "
                      "size_from_callable() or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_with_size_v<typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "remaining() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "remaining() must be called after next() for "
          "static sized parts");
    }

    return new_appended(vec2d_rt_part_construction_placeholder());
  }

  [[nodiscard]] inline auto
  skip_remaining() noexcept(noexcept(this->skip_remaning_impl(
      std::make_index_sequence<init_parts_type::arity>(),
      std::make_index_sequence<build_parts_type::arity>()))) {

    static_assert(init_parts_type::arity >= 2);
    static_assert(init_parts_type::arity == build_parts_type::arity + 1,
                  "skip_remaining() must be called on the last element after "
                  "[copy|move]_to() and before from()");
    static constexpr std::size_t element_index = get_element_index();

    // TODO: check that not construction is needed using the actual size values

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_with_max_size_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "skip_remaining() must be called after with_max_size() for "
            "parts built with dynamic_size() but without size_from_callable() "
            "or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_with_size_v<typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "skip_remaining() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "skip_remaining() must be called after next() for "
          "static sized parts");
    }

    return skip_remaning_impl(
        std::make_index_sequence<init_parts_type::arity>(),
        std::make_index_sequence<build_parts_type::arity>());
  }

  template <typename Fun>
  [[nodiscard]] inline auto
  transform(Fun&& fun) noexcept(noexcept(this->new_replaced_last(
      vec2d_rt_part_transform(init_parts.take_last(),
                              std::forward<Fun>(fun))))) {
    static_assert(init_parts_type::arity >= 2);
    static_assert(
        init_parts_type::arity > build_parts_type::arity * 2 + 2,
        "transform() must be called on elements after from().start(), after "
        "the declaration of the size and before next() or done()");

    static constexpr std::size_t element_index = get_element_index();

    namespace ct = boost::callable_traits;
    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(is_vec2d_rt_part_with_max_size_v<
                          typename init_parts_type::last_type>,
                      "transform() must be called after with_max_size() for "
                      "parts built with dynamic_size() but without "
                      "size_from_callable()");
      else {
        static_assert(
            is_vec2d_rt_part_with_size_v<typename init_parts_type::last_type>,
            "transform() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable()");
      }

      static_assert(
          std::tuple_size_v<ct::args_t<std::decay_t<Fun>>> == 1,
          "transform() for parts built with dynamic_size() must be a callable "
          "object taking one const reference to the relative element");

      static_assert(
          std::is_convertible_v<
              std::tuple_element_t<0, ct::args_t<std::decay_t<Fun>>>,
              typename allocator_traits<
                  allocator_type>::template const_reference<element_index>>,
          "transform() for parts built with dynamic_size() must be a callable "
          "object taking a const reference to the relative element");

      static_assert(
          std::is_convertible_v<ct::return_type_t<std::decay_t<Fun>>,
                                typename allocator_traits<allocator_type>::
                                    template value_type<element_index>>,
          "transform() for parts built with dynamic_size() must be a callable "
          "object returning and instance of the relative element");
    } else {
      static_assert(
          is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "transform() must be called after start() or next() for static sized "
          "parts");

      static_assert(std::tuple_size_v<ct::args_t<std::decay_t<Fun>>> ==
                        last_build_part_type::size,
                    "transform() for parts built with static_size() must be a "
                    "callable object taking a number of const references equal "
                    "to the size specified by set_size()");

      if constexpr (last_build_part_type::size == 1) {
        static_assert(
            nostd::is_applicable_r_v<
                typename allocator_traits<allocator_type>::template value_type<
                    element_index>,
                std::decay_t<Fun>,
                std::array<typename allocator_traits<allocator_type>::
                               template value_type<element_index> const&,
                           1>> or
                nostd::is_applicable_r_v<
                    std::array<typename allocator_traits<allocator_type>::
                                   template value_type<element_index>,
                               1>,
                    std::decay_t<Fun>,
                    std::array<typename allocator_traits<allocator_type>::
                                   template value_type<element_index> const&,
                               1>>,
            "transform() for parts built with static_size() and set_size<1>() "
            "must be a callable object taking 1 const reference to the "
            "relative element. Furthermore, the callable must return an "
            "instance of T or an array of T with size 1, where T is the type "
            "of the relative element.");
      } else {
        static_assert(
            nostd::is_applicable_r_v<
                std::array<typename allocator_traits<allocator_type>::
                               template value_type<element_index>,
                           last_build_part_type::size>,
                std::decay_t<Fun>,
                std::array<typename allocator_traits<allocator_type>::
                               template value_type<element_index> const&,
                           last_build_part_type::size>>,
            "transform() for parts built with static_size() must be a "
            "callable object taking N const references to the relative "
            "element, where N is the size specified by set_size(). "
            "Furthermore, the callable must return an array of values, "
            "with the type of the relative element and the size equal "
            "to the value specified by set_size()");
      }
    }

    return new_replaced_last(vec2d_rt_part_transform(init_parts.take_last(),
                                                     std::forward<Fun>(fun)));
  }

  [[nodiscard]] inline auto
  next() noexcept(noexcept(this->new_appended(vec2d_rt_part_empty()))) {
    static_assert(init_parts_type::arity >= 2);
    static_assert(init_parts_type::arity != build_parts_type::arity + 1,
                  "next() cannot be called on the last element after "
                  "[copy|move]_to(). Use remaining() instead in order to "
                  "define the construction behaviour for the empty space");

    static constexpr std::size_t element_index = get_element_index();

    if constexpr (init_parts_type::arity <= build_parts_type::arity + 1 or
                  init_parts_type::arity > build_parts_type::arity * 2 + 2) {
      using last_build_part_type =
          typename build_parts_type::template type<element_index>;
      if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
        if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                          last_build_part_type>)
          static_assert(is_vec2d_rt_part_with_max_size_v<
                            typename init_parts_type::last_type> or
                            is_vec2d_rt_part_with_fun_size_v<
                                typename init_parts_type::last_type>,
                        "next() must be called after with_max_size() for "
                        "parts built with dynamic_size() but without "
                        "size_from_callable() or after with_fun_size()");
        else {
          static_assert(is_vec2d_rt_part_with_size_v<
                            typename init_parts_type::last_type> or
                            is_vec2d_rt_part_with_fun_size_v<
                                typename init_parts_type::last_type>,
                        "next() must be called after with_size() for parts "
                        "built with dynamic_size() and size_from_callable() or "
                        "after with_fun_size()");
        }
      } else {
        static_assert(
            is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
            "next() must be called after start() or next() for "
            "static sized parts");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_generic_construct_v<
              typename init_parts_type::last_type>,
          "next() must be called after construct_default(), "
          "uninitialized(), construct_copies(), construct_from_iter(), "
          "construct_from_args_iter() or, only in case of parts built "
          "with fixed_size(), construct() or construct_all()");

      using last_build_part_type =
          typename build_parts_type::template type<element_index>;
      if constexpr (is_vec2d_rt_part_construct_v<last_build_part_type>) {
        static_assert(is_vec2d_fixed_size_part_v<build_parts_type>);

        static_assert(
            last_build_part_type::args_tuples_arity == init_parts_type::size,
            "next() is called after construct(), but not all the elements "
            "specified by fixed_size().set_size() have been constructed. Use "
            "more construct() calls or use construct_all() to specify the "
            "construction of the remaining elements");
      }
    }

    static_assert(
        element_index < build_parts_type::arity,
        "next() cannot be called on the last element. Use done() instead");

    return new_appended(vec2d_rt_part_empty());
  }

  inline auto
  done() {
    static_assert(init_parts_type::arity >= 2);
    static_assert(init_parts_type::arity != build_parts_type::arity + 1,
                  "done() cannot be called on the last element after "
                  "[copy|move]_to(). Use remaining() instead in order to "
                  "define the construction behaviour for the empty space");

    static constexpr std::size_t element_index = get_element_index();

    if constexpr (init_parts_type::arity <= build_parts_type::arity + 1 or
                  init_parts_type::arity > build_parts_type::arity * 2 + 2) {
      using last_build_part_type =
          typename build_parts_type::template type<element_index>;
      if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
        if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                          last_build_part_type>)
          static_assert(is_vec2d_rt_part_with_max_size_v<
                            typename init_parts_type::last_type> or
                            is_vec2d_rt_part_with_fun_size_v<
                                typename init_parts_type::last_type>,
                        "done() must be called after with_max_size() for "
                        "parts built with dynamic_size() but without "
                        "size_from_callable() or after with_fun_size()");
        else {
          static_assert(is_vec2d_rt_part_with_size_v<
                            typename init_parts_type::last_type> or
                            is_vec2d_rt_part_with_fun_size_v<
                                typename init_parts_type::last_type>,
                        "done() must be called after with_size() for parts "
                        "built with dynamic_size() and size_from_callable() or "
                        "after with_fun_size()");
        }
      } else {
        static_assert(
            is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
            "done() must be called after start() or next() for "
            "static sized parts");
      }
    } else {
      using last_build_part_type =
          typename build_parts_type::template type<element_index>;
      if constexpr (is_vec2d_rt_part_construct_v<last_build_part_type>) {
        static_assert(is_vec2d_fixed_size_part_v<build_parts_type>);

        static_assert(
            last_build_part_type::args_tuples_arity == init_parts_type::size,
            "done() is called after construct(), but not all the elements "
            "specified by fixed_size().set_size() have been constructed. Use "
            "more construct() calls or use construct_all() to specify the "
            "construction of the remaining elements");
      }
    }

    static_assert(
        element_index == build_parts_type::arity - 1,
        "done() can be called only on the last element. Use next() instead");

    if constexpr (init_parts_type::arity == build_parts_type::arity * 3 + 2) {
      perform_move(std::make_index_sequence<build_parts_type::arity>());

      return new_appended(vec2d_rt_part_unwinder());
    } else {
      return new_appended(vec2d_rt_part_empty());
    }
  }

  inline void
  unwind() noexcept {
    static_assert(
        is_vec2d_rt_part_unwinder_v<typename init_parts_type::last_type>,
        "unwind() can be called only after the last done()");

    return unwind_impl(std::make_index_sequence<build_parts_type::arity>());
  }

private:
  // TODO noexcept specifier
  template <std::size_t... Idx>
  inline void
  perform_move(std::index_sequence<Idx...>) {

    constexpr std::size_t init_part_from_index =
        build_parts_type::arity * 2 + 1;
    using init_part_from =
        typename init_parts_type::template type<init_part_from_index>;
    using init_part_to = typename init_parts_type::template type<0>;

    static_assert(is_vec2d_rt_part_with_lines_v<init_part_to> or
                  is_vec2d_rt_part_skip_first_lines_v<init_part_to>);
    static_assert(is_vec2d_rt_part_from_address_v<init_part_from>);

    const auto from_address =
        std::get<init_part_from_index>(init_parts).from_address;
    const auto to_element_max_sizes =
        std::tuple(base_type::template get_init_part_max_size<Idx + 1>()...);
    const auto from_element_max_sizes = std::tuple(
        base_type::template get_init_part_max_size<Idx + init_part_from_index +
                                                   1>()...);
    static_assert(std::tuple_size_v<decltype(to_element_max_sizes)> ==
                  std::tuple_size_v<decltype(from_element_max_sizes)>);

    using lines_type = typename init_part_from::size_type;
    static_assert(std::is_same_v<lines_type, typename init_part_to::size_type>);
    auto&& from_lines = std::get<init_part_from_index>(init_parts).lines;
    auto&& to_lines = std::get<0>(init_parts).lines;

    const auto from_zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<
                        typename init_part_from::base_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<
                      typename init_part_from::base_type>);
        return std::get<init_part_from_index>(init_parts).skip_lines;
      }
    }();

    const auto to_zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<init_part_to>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<init_part_to>);
        return std::get<0>(init_parts).skip_lines;
      }
    }();

    lines_type to_line_index = to_zero_line;
    lines_type from_line_index = from_zero_line;
    for (; from_line_index < from_lines and to_line_index < to_lines;
         ++from_line_index, ++to_line_index) {
      std::size_t moved_elements = 0;
      try {
        (
            [&] {
              perform_element_move<Idx>(from_address, from_line_index,
                                        to_line_index, from_element_max_sizes,
                                        to_element_max_sizes);
              ++moved_elements;
            }(),
            ...);
      } catch (...) {
        unwind_line_partial(
            from_address, from_line_index, to_line_index, moved_elements,
            from_element_max_sizes, to_element_max_sizes,
            nostd::make_index_sequence_rev<build_parts_type::arity>());
        while (from_line_index > from_zero_line and
               to_line_index > to_zero_line)
          unwind_line(
              from_address, --from_line_index, --to_line_index,
              from_element_max_sizes, to_element_max_sizes,
              nostd::make_index_sequence_rev<build_parts_type::arity>());

        throw;
      }
    }

    if constexpr (not is_vec2d_rt_part_skip_v<
                      typename init_parts_type::template type<
                          build_parts_type::arity + 1>>) {
      const lines_type remaining_lines = to_lines - to_line_index;
      const lines_type remaining_line_zero(0);
      for (lines_type remaining_line_index = remaining_line_zero;
           remaining_line_index < remaining_lines;
           ++remaining_line_index, ++to_line_index) {

        std::size_t constructed_elements = 0;
        try {
          (std::apply(
               [&, this](auto&&... sizes) {
                 constexpr std::size_t remaining_part_index =
                     build_parts_type::arity + 1 + Idx;

                 construction_base_type::
                     template perform_all_elements_build_from_init<
                         remaining_part_index>(
                         to_line_index,
                         std::forward<decltype(sizes)>(sizes)...);
                 ++constructed_elements;
               },
               to_element_max_sizes),
           ...);
        } catch (...) {
          std::apply(
              [&, this](auto&&... sizes) {
                unwind_destroy_line_partial(
                    to_line_index, constructed_elements,
                    nostd::make_index_sequence_rev<sizeof...(Idx)>(),
                    std::forward<decltype(sizes)>(sizes)...);
              },
              to_element_max_sizes);

          while (remaining_line_index > remaining_line_zero) {
            std::apply(
                [&, this](auto&&... sizes) {
                  unwind_destroy_line(
                      --to_line_index,
                      nostd::make_index_sequence_rev<sizeof...(Idx)>(),
                      std::forward<decltype(sizes)>(sizes)...);
                  --remaining_line_index;
                },
                to_element_max_sizes);
          }

          while (from_line_index > from_zero_line and
                 to_line_index > to_zero_line)
            unwind_line(
                from_address, --from_line_index, --to_line_index,
                from_element_max_sizes, to_element_max_sizes,
                nostd::make_index_sequence_rev<build_parts_type::arity>());

          throw;
        }
      }
    }
  }

  // TODO noexcept specifier
  template <std::size_t Index, typename FromAddress, typename FromLineSize,
            typename ToLineSize, typename FromElementSizes,
            typename ToElementSizes>
  inline void
  perform_element_move(FromAddress&& from_address,
                       FromLineSize&& from_line_index,
                       ToLineSize&& to_line_index,
                       FromElementSizes&& from_element_sizes,
                       ToElementSizes&& to_element_sizes) {
    const auto first_pointer_from = std::apply(
        [&, this](auto&&... sizes) {
          return allocator_traits<Alloc>::template first_pointer_of<Index>(
              *alloc, from_address, from_line_index,
              std::forward<decltype(sizes)>(sizes)...);
        },
        from_element_sizes);

    const auto first_pointer_to = std::apply(
        [&, this](auto&&... sizes) {
          return allocator_traits<Alloc>::template first_pointer_of<Index>(
              *alloc, address, to_line_index,
              std::forward<decltype(sizes)>(sizes)...);
        },
        to_element_sizes);

    using element_type =
        typename allocator_traits<allocator_type>::template value_type<Index>;

    constexpr std::size_t init_parts_from_index =
        build_parts_type::arity * 2 + 1;
    constexpr std::size_t from_init_part_index =
        init_parts_from_index + 1 + Index;

    auto&& n_from_elements = std::apply(
        [&, this](auto&&... sizes) {
          return base_type::template get_init_part_n_elements_at<
              from_init_part_index>(from_address,
                                    std::forward<FromLineSize>(from_line_index),
                                    std::forward<decltype(sizes)>(sizes)...);
        },
        from_element_sizes);

    auto&& n_to_elements = std::apply(
        [&, this](auto&&... sizes) {
          return base_type::template get_init_part_n_elements<Index + 1>(
              std::forward<ToLineSize>(to_line_index),
              std::forward<decltype(sizes)>(sizes)...);
        },
        to_element_sizes);

    using NElementsType = std::decay_t<decltype(n_to_elements)>;
    static_assert(
        std::is_same_v<NElementsType, std::decay_t<decltype(n_from_elements)>>);
    auto&& moving_n_elements = std::min(n_from_elements, n_to_elements);
    auto remaining_n_elements = [&] {
      if (n_to_elements > n_from_elements)
        return n_to_elements - n_from_elements;
      else
        return NElementsType(0);
    }();

    using from_init_part_type =
        typename init_parts_type::template type<from_init_part_index>;

    auto from_data = [&]() {
      if constexpr (is_vec2d_rt_part_transform_v<from_init_part_type>) {
        return base_type::template call_transform_on_fixed<
            from_init_part_index>(first_pointer_from);
      } else {
        return []() {};
      }
    }();

    auto data_from_ptr = [&] {
      if constexpr (is_vec2d_rt_part_transform_v<from_init_part_type>) {
        if constexpr (nostd::is_array_class_v<
                          std::decay_t<decltype(from_data)>>)
          return from_data.data();
        else
          return &from_data;
      } else
        return first_pointer_from;
    }();
    auto data_to_ptr = first_pointer_to;

    for (NElementsType element_index(0); element_index < moving_n_elements;
         ++element_index, ++data_from_ptr, ++data_to_ptr) {
      try {
        if constexpr (is_vec2d_rt_part_transform_v<from_init_part_type> or
                      (std::is_nothrow_move_constructible_v<element_type> and
                       std::is_nothrow_move_assignable_v<element_type>)) {
          allocator_traits<allocator_type>::construct(
              *alloc, data_to_ptr, std::move(*data_from_ptr));
        } else {
          allocator_traits<allocator_type>::construct(*alloc, data_to_ptr,
                                                      *data_from_ptr);
        }
      } catch (...) {
        unwind_element_move<Index>(data_from_ptr, data_to_ptr, element_index,
                                   NElementsType(0));

        throw;
      }
    }

    if constexpr (not is_vec2d_rt_part_skip_v<
                      typename init_parts_type::template type<
                          build_parts_type::arity + 1>>) {
      try {
        constexpr std::size_t remaining_init_part_index =
            build_parts_type::arity + 1 + Index;
        construction_base_type::template perform_elements_build_from_init<
            remaining_init_part_index>(first_pointer_to, remaining_n_elements,
                                       moving_n_elements);
      } catch (...) {
        unwind_element_move<Index>(data_from_ptr, data_to_ptr,
                                   moving_n_elements, NElementsType(0));
        throw;
      }
    }
  }

  template <std::size_t Index, typename T, typename ElementSizeT_A,
            typename ElementSizeT_B>
  inline void
  unwind_element_move(T* from_element_exception_ptr [[maybe_unused]],
                      T* to_element_exception_ptr,
                      ElementSizeT_A exception_index,
                      ElementSizeT_B&& first_element_index) noexcept {
    static_assert(std::is_same_v<std::decay_t<ElementSizeT_A>,
                                 std::decay_t<ElementSizeT_B>>);

    constexpr std::size_t init_parts_from_index =
        build_parts_type::arity * 2 + 1;
    constexpr std::size_t from_init_part_index =
        init_parts_from_index + 1 + Index;
    using from_init_part_type =
        typename init_parts_type::template type<from_init_part_index>;
    using element_type =
        typename allocator_traits<allocator_type>::template value_type<Index>;

    while (exception_index > first_element_index) {
      --to_element_exception_ptr;
      if constexpr (not is_vec2d_rt_part_transform_v<from_init_part_type> and
                    std::is_nothrow_move_constructible_v<element_type> and
                    std::is_nothrow_move_assignable_v<element_type>)
        *--from_element_exception_ptr = std::move(*to_element_exception_ptr);

      allocator_traits<allocator_type>::destroy(*alloc,
                                                to_element_exception_ptr);
      --exception_index;
    }
  }

  template <typename FromAddress, typename FromLineSize, typename ToLineSize,
            typename FromElementSizes, typename ToElementSizes,
            std::size_t... Idx>
  inline void
  unwind_line_partial(FromAddress&& from_address,
                      FromLineSize&& from_line_index,
                      ToLineSize&& to_line_index,
                      std::size_t element_exception_index,
                      FromElementSizes&& from_element_sizes,
                      ToElementSizes&& to_element_sizes,
                      std::index_sequence<Idx...>) const noexcept {
    static_assert(sizeof...(Idx) < 2 or std::get<0>(std::array{Idx...}) >
                                            std::get<1>(std::array{Idx...}));
    (
        [&] {
          if (Idx < element_exception_index) {
            unwind_element<Idx>(from_address, from_line_index, to_line_index,
                                from_element_sizes, to_element_sizes);
          }
        }(),
        ...);
  }

  template <typename FromAddress, typename FromLineSize, typename ToLineSize,
            typename FromElementSizes, typename ToElementSizes,
            std::size_t... Idx>
  inline void
  unwind_line(FromAddress&& from_address, FromLineSize&& from_line_index,
              ToLineSize&& to_line_index, FromElementSizes&& from_element_sizes,
              ToElementSizes&& to_element_sizes,
              std::index_sequence<Idx...>) const noexcept {
    static_assert(sizeof...(Idx) < 2 or std::get<0>(std::array{Idx...}) >
                                            std::get<1>(std::array{Idx...}));
    (
        [&] {
          unwind_element<Idx>(from_address, from_line_index, to_line_index,
                              from_element_sizes, to_element_sizes);
        }(),
        ...);
  }

  template <std::size_t Index, typename FromAddress, typename FromLineSize,
            typename ToLineSize, typename FromElementSizes,
            typename ToElementSizes>
  inline void
  unwind_element(FromAddress&& from_address, FromLineSize&& from_line_index,
                 ToLineSize&& to_line_index,
                 FromElementSizes&& from_element_sizes,
                 ToElementSizes&& to_element_sizes) const noexcept {

    using element_type = typename build_parts_type::template type<Index>;
    auto n_elements = std::apply(
        [&, this](auto&&... sizes) {
          return base_type::template get_init_part_n_elements<Index + 1>(
              to_line_index, sizes...);
        },
        to_element_sizes);

    auto data_ptr_from =
        std::apply(
            [&, this](auto&&... sizes) {
              return allocator_traits<Alloc>::template first_pointer_of<Index>(
                  *alloc, from_address, from_line_index, sizes...);
            },
            std::forward<FromElementSizes>(from_element_sizes)) +
        n_elements;

    auto data_ptr_to =
        std::apply(
            [&, this](auto&&... sizes) {
              return allocator_traits<Alloc>::template first_pointer_of<Index>(
                  *alloc, address, to_line_index, sizes...);
            },
            std::forward<ToElementSizes>(to_element_sizes)) +
        n_elements;

    const std::decay_t<decltype(n_elements)> element_zero(0);
    while (n_elements > element_zero) {
      if constexpr (std::is_nothrow_move_constructible_v<element_type> and
                    std::is_nothrow_move_assignable_v<element_type>) {
        *--data_ptr_from = std::move(*--data_ptr_to);
        destroy_element(data_ptr_to);
      } else
        destroy_element(--data_ptr_to);

      --n_elements;
    }
  }

  template <std::size_t... Idx>
  inline void
  unwind_impl(std::index_sequence<Idx...>) noexcept {
    constexpr std::size_t init_parts_from_index =
        build_parts_type::arity * 2 + 1;
    using to_init_type = typename init_parts_type::template type<0>;
    using lines_type = typename to_init_type::size_type;
    const auto from_address =
        std::get<init_parts_from_index>(init_parts).from_address;
    const auto to_element_max_sizes =
        std::tuple(base_type::template get_init_part_max_size<Idx + 1>()...);
    const auto from_element_max_sizes = std::tuple(
        base_type::template get_init_part_max_size<Idx + init_parts_from_index +
                                                   1>()...);
    static_assert(std::tuple_size_v<decltype(to_element_max_sizes)> ==
                  std::tuple_size_v<decltype(from_element_max_sizes)>);

    auto from_line_index = std::get<init_parts_from_index>(init_parts).lines;
    auto to_line_index = std::get<0>(init_parts).lines;

    const auto to_zero_line = [&] {
      if constexpr (is_vec2d_rt_part_with_lines_v<to_init_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<to_init_type>);
        return std::get<0>(init_parts).skip_lines;
      }
    }();

    const auto from_zero_line = [&] {
      using from_init_type =
          typename init_parts_type::template type<init_parts_from_index>;
      if constexpr (is_vec2d_rt_part_with_lines_v<
                        typename from_init_type::base_type>)
        return lines_type(0);
      else {
        static_assert(is_vec2d_rt_part_skip_first_lines_v<
                      typename from_init_type::base_type>);
        return std::get<init_parts_from_index>(init_parts).skip_lines;
      }
    }();

    while (from_line_index > from_zero_line and to_line_index > to_zero_line)
      unwind_line(from_address, --from_line_index, --to_line_index,
                  from_element_max_sizes, to_element_max_sizes,
                  nostd::make_index_sequence_rev<build_parts_type::arity>());
  }

  template <std::size_t... InitsIdx, std::size_t... SkipIdx>
  inline auto
  skip_remaning_impl(
      std::index_sequence<InitsIdx...>,
      std::index_sequence<
          SkipIdx...>) noexcept((std::
                                     is_nothrow_move_constructible_v<
                                         typename init_parts_type::
                                             template type<InitsIdx>> &&
                                 ...)) {
    using new_init_parts_type = vec2d_init_parts<
        typename init_parts_type::template type<InitsIdx>...,
        nostd::identity_type_t<vec2d_rt_part_skip, SkipIdx>...>;

    return vec2d_element_builder_rt_move_to<Alloc, BuildParts,
                                            new_init_parts_type>(
        *base_type::build_parts, base_type::address, *base_type::alloc,
        new_init_parts_type(
            std::move(std::get<InitsIdx>(base_type::init_parts))...,
            nostd::identity_type_t<vec2d_rt_part_skip, SkipIdx>()...));
  }
};

} // namespace het::detail

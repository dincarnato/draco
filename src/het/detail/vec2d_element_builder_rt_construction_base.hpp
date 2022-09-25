#pragma once

#include "../allocator_traits.hpp"
#include "vec2d_element_builder_ct_parts_traits.hpp"
#include "vec2d_element_builder_rt_parts.hpp"
#include "vec2d_element_builder_rt_parts_traits.hpp"

#include <cassert>

namespace het::detail {

template <template <typename, typename, typename> typename TemplatedBuilder,
          typename Alloc, typename BuildParts, typename InitParts>
struct vec2d_element_builder_rt_construction_base {
  template <typename _Alloc, typename _BuildParts, typename _InitParts>
  using templated_builder_type =
      TemplatedBuilder<_Alloc, _BuildParts, _InitParts>;
  using builder_type = templated_builder_type<Alloc, BuildParts, InitParts>;

  using allocator_type = Alloc;
  using build_parts_type = BuildParts;
  using init_parts_type = InitParts;
  using base_address_type = typename allocator_traits<Alloc>::first_pointer;

  // TODO noexcept specifications
  [[nodiscard]] inline auto construct_default() {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_max_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_default() must be called after with_max_size() for "
            "parts built with dynamic_size() but without size_from_callable() "
            "or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_default() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_construction_placeholder_v<
              typename init_parts_type::last_type> or
              is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "construct_default() must be called after with_lines() or next() for "
          "static sized parts");
    }

    return static_cast<builder_type &>(*this).new_replaced_last(
        vec2d_rt_part_construct_default(
            static_cast<builder_type &>(*this).init_parts.take_last()));
  }

  // TODO noexcept specifications
  [[nodiscard]] inline auto uninitialized() {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using element_type =
        typename allocator_traits<Alloc>::template value_type<element_index>;
    static_assert(
        std::is_trivially_destructible_v<element_type>,
        "uninitialized() can be used only for trivially destructible types");

    using last_build_part_type =
        typename build_parts_type::template type<init_parts_type::arity - 2>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_max_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "uninitialized() must be called after with_max_size() for "
            "parts built with dynamic_size() but without size_from_callable() "
            "or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "uninitialized() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_construction_placeholder_v<
              typename init_parts_type::last_type> or
              is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "uninitialized() must be called after with_lines() or next() for "
          "static sized parts");
    }

    static_assert(
        not is_vec2d_element_default_uninitialized_v<last_build_part_type>,
        "the type is already uninitialized by default. Use construct_default() "
        "instead");

    return static_cast<builder_type &>(*this).new_replaced_last(
        vec2d_rt_part_uninitialized(
            static_cast<builder_type &>(*this).init_parts.take_last()));
  }

  template <typename Value>
  [[nodiscard]] inline auto construct_copies(Value &&value) noexcept(
      noexcept(static_cast<builder_type &>(*this).new_replaced_last(
          vec2d_rt_part_construct_copies(
              static_cast<builder_type &>(*this).init_parts.take_last(),
              std::forward<Value>(value))))) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_max_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_copies() must be called after with_max_size() for "
            "parts built with dynamic_size() but without size_from_callable() "
            "or after with_size_fun()");
      else {
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_copies() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_size_fun()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_construction_placeholder_v<
              typename init_parts_type::last_type> or
              is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "construct_copies() must be called after with_lines() or next() for "
          "static sized parts");
    }

    return static_cast<builder_type &>(*this).new_replaced_last(
        vec2d_rt_part_construct_copies(
            static_cast<builder_type &>(*this).init_parts.take_last(),
            std::forward<Value>(value)));
  }

  // TODO: noexcept specifier
  template <typename... Args>
  [[nodiscard]] inline auto construct(Args &&...args) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    static_assert(
        is_vec2d_fixed_size_part_v<last_build_part_type>,
        "construct() must be used with parts built with fixed_size(), "
        "otherwise the number of elements to construct is unknown at "
        "compile-time. Use construct_from_iter() or construct_from_args_iter() "
        "instead or consider using fixed_size() to build the part");

    static_assert(
        is_vec2d_rt_part_construction_placeholder_v<
            typename init_parts_type::last_type> or
            is_vec2d_rt_part_empty_v<typename init_parts_type::last_type> or
            is_vec2d_rt_part_construct_v<typename init_parts_type::last_type>,
        "construct() must be called after with_lines(), next() or another "
        "construct()");

    static_assert(
        nostd::is_constructible_v<typename allocator_traits<Alloc>::
                                      template value_type<element_index>,
                                  Args...>,
        "the underline element is not constructable using the passed "
        "arguments");

    if constexpr (is_vec2d_rt_part_construct_v<
                      typename init_parts_type::last_type>) {
      using construct_type = typename init_parts_type::last_type;
      static_assert(construct_type::args_tuples_arity <
                        last_build_part_type::size,
                    "construct() cannot be called more than the size set using "
                    "fixed_size().set_size()");
      return static_cast<builder_type &>(*this).new_replaced_last(
          static_cast<builder_type &>(*this)
              .init_parts.take_last()
              .new_appended(std::forward<Args>(args)...));
    } else {
      return static_cast<builder_type &>(*this).new_replaced_last(
          vec2d_rt_part_construct(
              static_cast<builder_type &>(*this).init_parts.take_last(),
              std::tuple(std::forward<Args>(args)...)));
    }
  }

  // TODO: noexcept specifier
  template <typename... TupleArgs>
  [[nodiscard]] inline auto construct_all(TupleArgs &&...tuple_args) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();
    static_assert((nostd::is_tuple_v<std::decay_t<TupleArgs>> && ...),
                  "all the argument of construct_all() must be tuples");

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    static_assert(
        is_vec2d_fixed_size_part_v<last_build_part_type>,
        "construct_all() must be used with parts built with fixed_size(), "
        "otherwise the number of elements to construct is unknown at "
        "compile-time. Use construct_from_iter() or construct_from_args_iter() "
        "instead or consider using fixed_size() to build the part");

    static_assert(
        is_vec2d_rt_part_construction_placeholder_v<
            typename init_parts_type::last_type> or
            is_vec2d_rt_part_empty_v<typename init_parts_type::last_type> or
            is_vec2d_rt_part_construct_v<typename init_parts_type::last_type>,
        "construct_all() must be called after with_lines(), next() or another "
        "construct()");

    static_assert((nostd::is_constructible_from_tuple_v<
                       typename allocator_traits<Alloc>::template value_type<
                           element_index>,
                       std::decay_t<TupleArgs>> &&
                   ...),
                  "all the tuple must contain sets of values that can be used "
                  "to construct the underlying element");

    if constexpr (is_vec2d_rt_part_construct_v<
                      typename init_parts_type::last_type>) {
      using construct_type = typename init_parts_type::last_type;

      static_assert(
          construct_type::args_tuples_arity + sizeof...(TupleArgs) ==
              last_build_part_type::size,
          "construct_all() must be used to initialize all the remaining "
          "elements that are not constructed using construct(). Use a number "
          "of arguments according to the value specified by "
          "fixed_size().set_size() and the number of construct() calls");
      return static_cast<builder_type &>(*this).new_replaced_last(
          static_cast<builder_type &>(*this)
              .init_parts.take_last()
              .new_appending_tuples(std::forward<TupleArgs>(tuple_args)...));
    } else {
      static_assert(
          sizeof...(TupleArgs) == last_build_part_type::size,
          "construct_all() must be used to initialize all the elements. The "
          "number of arguments must be the same as the value used for "
          "fixed_size().set_size()");
      return static_cast<builder_type &>(*this).new_replaced_last(
          vec2d_rt_part_construct(
              static_cast<builder_type &>(*this).init_parts.take_last(),
              std::forward<TupleArgs>(tuple_args)...));
    }
  }

  template <typename Iter>
  [[nodiscard]] inline auto construct_from_iter(Iter &&iter) noexcept(
      noexcept(static_cast<builder_type &>(*this).new_appended(
          vec2d_rt_part_construct_from_iter(std::forward<Iter>(iter))))) {

    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_max_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_from_iter() must be called after with_max_size() for "
            "parts built with dynamic_size() but without size_from_callable() "
            "or after with_fun_size()");
      else {
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_from_iter() must be called after with_size() for parts "
            "built with dynamic_size() and size_from_callable() or after "
            "with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_construction_placeholder_v<
              typename init_parts_type::last_type> or
              is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "construct_from_iter() must be called after with_lines() or next() "
          "for static sized parts");
    }

    static_assert(
        nostd::is_constructible_v<typename allocator_traits<Alloc>::
                                      template value_type<element_index>,
                                  typename Iter::reference>,
        "the passed iterator's reference type must be usable as argument to be "
        "passed to the constructor of the underlying element");

    return static_cast<builder_type &>(*this).new_appended(
        vec2d_rt_part_construct_from_iter(
            static_cast<builder_type &>(*this).init_parts.take_last(),
            std::forward<Iter>(iter)));
  }

  template <typename Iter>
  [[nodiscard]] inline auto construct_from_args_iter(Iter &&iter) noexcept(
      noexcept(static_cast<builder_type &>(*this).new_appended(
          vec2d_rt_part_construct_from_args_iter(std::forward<Iter>(iter))))) {

    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using last_build_part_type =
        typename build_parts_type::template type<element_index>;
    if constexpr (is_vec2d_dynamic_size_part_v<last_build_part_type>) {
      if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                        last_build_part_type>)
        static_assert(
            is_vec2d_rt_part_construction_placeholder_v<
                typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_max_size_v<
                    typename init_parts_type::last_type> or
                is_vec2d_rt_part_with_fun_size_v<
                    typename init_parts_type::last_type>,
            "construct_from_args_iter() must be called after with_max_size() "
            "for parts built with dynamic_size() but without "
            "size_from_callable() or after with_fun_size()");
      else {
        static_assert(is_vec2d_rt_part_construction_placeholder_v<
                          typename init_parts_type::last_type> or
                          is_vec2d_rt_part_with_size_v<
                              typename init_parts_type::last_type> or
                          is_vec2d_rt_part_with_fun_size_v<
                              typename init_parts_type::last_type>,
                      "construct_from_args_iter() must be called after "
                      "with_size() for parts built with dynamic_size() and "
                      "size_from_callable() or after with_fun_size()");
      }
    } else {
      static_assert(
          is_vec2d_rt_part_construction_placeholder_v<
              typename init_parts_type::last_type> or
              is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
          "construct_from_args_iter() must be called after with_lines() or "
          "next() for static sized parts");
    }

    static_assert(
        nostd::is_constructible_from_tuple_v<
            typename allocator_traits<Alloc>::template value_type<
                element_index>,
            typename Iter::value_type>,
        "the passed iterator's value type must be a tuple of arguments that "
        "can be usable for the construction of the underlying element");

    return static_cast<builder_type &>(*this).new_appended(
        vec2d_rt_part_construct_from_all_iter(
            static_cast<builder_type &>(*this).init_parts.take_last(),
            std::forward<Iter>(iter)));
  }

protected:
  template <std::size_t Index, typename LineSizeT, typename... Sizes>
  inline void perform_all_elements_build_from_init(LineSizeT &&line_index,
                                                   Sizes &&...sizes) const {
    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t element_index =
        builder_type::template get_element_index_for_init<Index>();

    auto first_pointer =
        allocator_traits<allocator_type>::template first_pointer_of<
            element_index>(*builder.alloc, builder.address, line_index,
                           sizes...);

    const auto &&n_elements = builder.template get_init_part_n_elements<Index>(
        std::forward<LineSizeT>(line_index), std::forward<Sizes>(sizes)...);

    perform_elements_build_from_init<Index>(first_pointer, n_elements);
  }

  template <std::size_t Index, typename T, typename ElementSizeT>
  inline void
  perform_elements_build_from_init(T *first_pointer, ElementSizeT &&n_elements,
                                   std::size_t skipped = std::size_t(0)) const {
    constexpr std::size_t element_index =
        builder_type::template get_element_index_for_init<Index>();

    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<element_index>>);

    using init_part_type = typename init_parts_type::template type<Index>;

    if constexpr (is_vec2d_rt_part_construct_default_v<init_part_type>) {
      perform_element_build_from_init_construct_default<Index>(
          first_pointer + skipped, std::forward<ElementSizeT>(n_elements));
    } else if constexpr (is_vec2d_rt_part_construct_copies_v<init_part_type>) {
      perform_element_build_from_init_construct_copies<Index>(
          first_pointer + skipped, std::forward<ElementSizeT>(n_elements));
    } else if constexpr (is_vec2d_rt_part_construct_v<init_part_type>) {
      assert(n_elements + skipped == init_part_type::args_tuples_arity);
      perform_element_build_from_init_construct<Index>(
          first_pointer, skipped,
          std::make_index_sequence<init_part_type::args_tuples_arity>());
    } else if constexpr (is_vec2d_rt_part_construct_from_iter_v<
                             init_part_type>) {
      perform_element_build_from_init_construct_from_iter<Index>(
          first_pointer + skipped, std::forward<ElementSizeT>(n_elements));
    } else if constexpr (is_vec2d_rt_part_construct_from_args_iter_v<
                             init_part_type>) {
      perform_element_build_from_init_construct_from_args_iter<Index>(
          first_pointer + skipped, std::forward<ElementSizeT>(n_elements));
    } else {
      static_assert(is_vec2d_rt_part_uninitialized_v<init_part_type>);
    }
  }

private:
  template <std::size_t Index, typename T, typename NElements>
  inline void perform_element_build_from_init_construct_default(
      T *first_pointer, NElements &&n_elements) const {

    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t element_index =
        builder_type::template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<element_index>>);
    using build_part_type =
        typename build_parts_type::template type<element_index>;

    const auto &build_part = std::get<element_index>(*builder.build_parts);

    auto data_ptr = first_pointer;
    for (std::decay_t<NElements> element_index(0); element_index < n_elements;
         ++element_index, ++data_ptr) {
      try {
        if constexpr (is_vec2d_element_default_construction_default_v<
                          build_part_type>)
          allocator_traits<Alloc>::construct(*builder.alloc, data_ptr);
        else if constexpr (is_vec2d_element_default_construction_v<
                               build_part_type>) {
          std::apply(
              [&, data_ptr](auto &&...inits) {
                allocator_traits<allocator_type>::construct(
                    *builder.alloc, data_ptr,
                    std::forward<decltype(inits)>(inits)...);
              },
              build_part.init_values);
        } else {
          static_assert(
              is_vec2d_element_default_uninitialized_v<build_part_type>);
        }
      } catch (...) {
        const std::decay_t<NElements> element_zero(0);
        while (element_index > element_zero) {
          builder.destroy_element(--data_ptr);
          --element_index;
        }

        throw;
      }
    }
  }

  template <std::size_t Index, typename T, typename NElements>
  inline void perform_element_build_from_init_construct_copies(
      T *first_pointer, NElements &&n_elements) const {

    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t element_index =
        builder_type::template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<element_index>>);

    const auto &init_part = std::get<Index>(builder.init_parts);
    auto data_ptr = first_pointer;
    for (std::decay_t<NElements> element_index(0); element_index < n_elements;
         ++element_index, ++data_ptr) {
      try {
        allocator_traits<Alloc>::construct(*builder.alloc, data_ptr,
                                           init_part.value);
      } catch (...) {
        const std::decay_t<NElements> element_zero(0);
        while (element_index > element_zero) {
          builder.destroy_element(--data_ptr);
          --element_index;
        }

        throw;
      }
    }
  }

  template <std::size_t Index, typename T, std::size_t... Idx>
  inline void
  perform_element_build_from_init_construct(T *first_pointer,
                                            std::size_t skipped,
                                            std::index_sequence<Idx...>) const {
    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t build_part_index =
        builder_type::template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<build_part_index>>);

    const auto &init_part = std::get<Index>(builder.init_parts);
    std::size_t constructed_elements = 0;
    auto data_ptr = first_pointer;
    try {
      (std::apply(
           [skipped, &data_ptr, &constructed_elements,
            &builder](auto &&...args) {
             if (Idx >= skipped) {
               allocator_traits<Alloc>::construct(
                   *builder.alloc, data_ptr,
                   std::forward<decltype(args)>(args)...);

               ++constructed_elements;
             }
             ++data_ptr;
           },
           std::get<Idx>(init_part.args_tuples)),
       ...);
    } catch (...) {
      while (constructed_elements > 0) {
        builder.destroy_element(--data_ptr);
        --constructed_elements;
      }

      throw;
    }
  }

  template <std::size_t Index, typename T, typename NElements>
  inline void perform_element_build_from_init_construct_from_iter(
      T *first_pointer, NElements &&n_elements) const {

    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t build_part_index =
        builder.template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<build_part_index>>);

    const auto &init_part = std::get<Index>(builder.init_parts);

    std::decay_t<NElements> element_index(0);
    auto data_ptr = first_pointer;
    for (auto iter = init_part.iter; element_index < n_elements;
         ++iter, ++element_index, ++data_ptr) {
      try {
        allocator_traits<allocator_type>::construct(*builder.alloc, data_ptr,
                                                    *iter);
      } catch (...) {
        const std::decay_t<NElements> element_zero(0);
        while (element_index > element_zero) {
          builder.destroy_element(--data_ptr);
          --element_index;
        }

        throw;
      }
    }
  }

  template <std::size_t Index, typename T, typename NElements>
  inline void perform_element_build_from_init_construct_from_args_iter(
      T *first_pointer, NElements &&n_elements) const {

    const auto &builder = static_cast<const builder_type &>(*this);

    constexpr std::size_t build_part_index =
        builder_type::template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<allocator_type>::
                                template pointer<build_part_index>>);

    const auto &init_part = std::get<Index>(builder.init_parts);

    std::decay_t<NElements> element_index(0);
    auto data_ptr = first_pointer;
    for (auto iter = init_part.iter; element_index < n_elements;
         ++iter, ++element_index, ++data_ptr) {
      try {
        std::apply(
            [this, data_ptr, &builder](auto &&...args) {
              allocator_traits<allocator_type>::construct(
                  *builder.alloc, data_ptr,
                  std::forward<decltype(args)>(args)...);
            },
            *iter);
      } catch (...) {
        const std::decay_t<NElements> element_zero(0);
        while (element_index > element_zero) {
          builder.destroy_element(--data_ptr);
          --element_index;
        }

        throw;
      }
    }
  }
};

} // namespace het::detail

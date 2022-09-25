#pragma once

#include "../allocator_traits.hpp"
#include "vec2d_element_builder_ct_parts.hpp"
#include "vec2d_element_builder_rt_parts.hpp"

#include <boost/callable_traits.hpp>

namespace het::detail {

template <template <typename, typename, typename> typename TemplatedBuilder,
          typename Alloc, typename BuildParts, typename InitParts>
struct vec2d_element_builder_rt_common {
  template <typename _Alloc, typename _BuildParts, typename _InitParts>
  using templated_builder_type =
      TemplatedBuilder<_Alloc, _BuildParts, _InitParts>;
  using builder_type = templated_builder_type<Alloc, BuildParts, InitParts>;

  using allocator_type = Alloc;
  using build_parts_type = BuildParts;
  using init_parts_type = InitParts;
  using base_address_type = typename allocator_traits<Alloc>::first_pointer;

  template <template <typename, typename, typename> typename, typename,
            typename, typename>
  friend struct vec2d_element_builder_rt_common;

protected:
  const build_parts_type *build_parts;
  base_address_type address;
  allocator_type *alloc;
  init_parts_type init_parts;

  template <typename Part>
  using new_appended_t =
      templated_builder_type<allocator_type, build_parts_type,
                             typename init_parts_type::template append_t<Part>>;

  template <typename Part>
  using new_replaced_last_t = templated_builder_type<
      allocator_type, build_parts_type,
      typename init_parts_type::template replace_last_with_t<Part>>;

  template <typename Part>
  [[nodiscard]] inline auto new_appended(Part &&part) noexcept(
      noexcept(new_appended_t<std::decay_t<Part>>(
          *build_parts, address, *alloc,
          init_parts.append(std::forward<Part>(part))))) {

    return new_appended_t<std::decay_t<Part>>(
        *build_parts, address, *alloc,
        init_parts.append(std::forward<Part>(part)));
  }

  template <typename Part>
  [[nodiscard]] inline auto new_replaced_last(Part &&part) noexcept(
      noexcept(new_replaced_last_t<std::decay_t<Part>>(
          *build_parts, address, *alloc,
          init_parts.replace_last_with(std::forward<Part>(part))))) {

    return new_replaced_last_t<std::decay_t<Part>>(
        *build_parts, address, *alloc,
        init_parts.replace_last_with(std::forward<Part>(part)));
  }

private:
  // TODO: noexcept specifier
  template <std::size_t Index, typename T, std::size_t... Idx>
  inline auto call_transform_on_fixed_impl(T *first_pointer,
                                           std::index_sequence<Idx...>) {
    using init_part_type = typename init_parts_type::template type<Index>;
    auto &init_part = std::get<Index>(init_parts);
    return static_cast<typename init_part_type::function_type &>(init_part)(
        first_pointer[Idx]...);
  }

public:
  inline vec2d_element_builder_rt_common() = default;
  inline vec2d_element_builder_rt_common(vec2d_element_builder_rt_common &&) =
      default;
  inline vec2d_element_builder_rt_common &
  operator=(vec2d_element_builder_rt_common &&) = default;
  inline vec2d_element_builder_rt_common(
      const build_parts_type &build_parts, base_address_type address,
      Alloc &alloc) noexcept(std::
                                 is_nothrow_copy_constructible_v<
                                     base_address_type>)
      : build_parts(&build_parts), address(address), alloc(&alloc) {
    static_assert(init_parts_type::arity > 0);
    static_assert(
        std::is_same_v<InitParts, vec2d_init_parts<vec2d_rt_part_empty>>,
        "Initial constructor can be used only when there are not init parts");
  }
  template <typename _InitParts>
  inline vec2d_element_builder_rt_common(
      const build_parts_type &build_parts, base_address_type address,
      Alloc &alloc,
      _InitParts &&init_parts) noexcept(std::
                                            is_nothrow_copy_constructible_v<
                                                base_address_type>
                                                and noexcept(init_parts_type(
                                                    std::forward<_InitParts>(
                                                        init_parts))))
      : build_parts(&build_parts), address(address), alloc(&alloc),
        init_parts(std::forward<_InitParts>(init_parts)) {}

  static constexpr std::size_t get_element_index() noexcept {
    static_assert(init_parts_type::arity > 0);
    return builder_type::template get_element_index_for_init<
        init_parts_type::arity - 1>();
  }

  template <typename SizeT>
  [[nodiscard]] inline auto
  with_lines(SizeT &&n) noexcept(noexcept(new_replaced_last(
      vec2d_rt_part_with_lines<std::decay_t<SizeT>>{std::forward<SizeT>(n)}))) {
    static_assert(builder_type::is_lines_init());
    static_assert(
        is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
        "with_lines() must be called after build_at() or after done()");

    return new_replaced_last(
        vec2d_rt_part_with_lines<std::decay_t<SizeT>>{std::forward<SizeT>(n)});
  }

  template <typename SizeT>
  [[nodiscard]] inline auto skip_lines(SizeT &&n) noexcept(noexcept(
      new_replaced_last(vec2d_rt_part_skip_first_lines<std::decay_t<SizeT>>{
          init_parts.take_last(), std::forward<SizeT>(n)}))) {
    static_assert(builder_type::is_lines_init());
    static_assert(
        is_vec2d_rt_part_with_lines_v<typename init_parts_type::last_type>,
        "skip_lines() must be called after with_lines()");

    return new_replaced_last(
        vec2d_rt_part_skip_first_lines<std::decay_t<SizeT>>{
            init_parts.take_last(), std::forward<SizeT>(n)});
  }

  [[nodiscard]] inline auto
  start() noexcept(noexcept(new_appended(vec2d_rt_part_empty()))) {
    static_assert(builder_type::is_lines_init());
    static_assert(
        not is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
        "start() must be called after with_lines() or from() (in case of "
        "copy_to() or move_to())");

    return new_appended(vec2d_rt_part_empty());
  }

  // TODO noexcept specifications
  template <typename SizeT> [[nodiscard]] inline auto with_max_size(SizeT &&n) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using element_size_type =
        typename allocator_traits<Alloc>::template size_type<element_index>;
    static_assert(std::is_convertible_v<SizeT, element_size_type>,
                  "SizeT must be convertible to the relative size_type in the "
                  "allocator_trait");
    static_assert(is_vec2d_rt_part_empty_v<typename init_parts_type::last_type>,
                  "with_max_size() must be called after with_lines(), "
                  "from()(if available) or next()");
    static_assert(
        is_vec2d_dynamic_size_part_v<
            typename build_parts_type::template type<element_index>>,
        "with_max_size() can only be called if the relative part is created "
        "with dynamic_size()");

    return new_replaced_last(
        vec2d_rt_part_with_max_size<SizeT>{std::forward<SizeT>(n)});
  }

  // TODO noexcept specifications
  template <typename SizeT> [[nodiscard]] inline auto with_size(SizeT &&n) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    using element_size_type =
        typename allocator_traits<Alloc>::template size_type<element_index>;
    static_assert(std::is_convertible_v<SizeT, element_size_type>,
                  "SizeT must be convertible to the relative size_type in the "
                  "allocator_trait");

    static_assert(
        is_vec2d_rt_part_with_max_size_v<typename init_parts_type::last_type>,
        "with_size() must be called after with_max_size()");
    static_assert(
        std::is_same_v<std::decay_t<SizeT>,
                       typename init_parts_type::last_type::size_type>,
        "SizeT must be the same type as used with with_max_size()");
    using part_type = typename build_parts_type::template type<element_index>;
    static_assert(
        is_vec2d_dynamic_size_part_v<part_type> and
            not is_vec2d_dynamic_size_from_callable_part_v<part_type>,
        "with_size() can only be called if the relative part is created "
        "with dynamic_size() and without size_from_callable()");

    return new_replaced_last(vec2d_rt_part_with_size<SizeT>{
        std::forward<SizeT>(n), init_parts.take_last().max_size});
  }

  // TODO noexcept specifications
  template <typename Fun> [[nodiscard]] inline auto with_fun_size(Fun &&fun) {
    static constexpr std::size_t element_index =
        builder_type::get_element_index();

    namespace ct = boost::callable_traits;
    using fun_return_type = ct::return_type_t<std::decay_t<Fun>>;

    using element_size_type =
        typename allocator_traits<Alloc>::template size_type<element_index>;
    static_assert(std::is_convertible_v<fun_return_type, element_size_type>,
                  "the callable object must return a type convertible to the "
                  "relative size_type in the allocator_trait");

    static_assert(
        is_vec2d_rt_part_with_max_size_v<typename init_parts_type::last_type>,
        "with_fun_size() must be called after with_max_size()");
    using part_type = typename build_parts_type::template type<element_index>;
    static_assert(is_vec2d_dynamic_size_part_v<part_type>,
                  "with_fun_size() can only be called if the relative part is "
                  "created with dynamic_size()");

    return new_replaced_last(
        vec2d_rt_part_with_fun_size<
            std::decay_t<Fun>, typename init_parts_type::last_type::size_type>{
            std::forward<Fun>(fun), init_parts.take_last().max_size});
  }

protected:
  // TODO noexcept specifier
  template <std::size_t Index, typename LineSizeT, typename... Sizes>
  [[nodiscard]] inline auto
  get_init_part_n_elements_at(base_address_type address, LineSizeT &&line_index,
                              Sizes &&...sizes) const noexcept {
    using init_part_type = typename init_parts_type::template type<Index>;
    static constexpr std::size_t build_part_index =
        builder_type::template get_element_index_for_init<Index>();
    using build_part_type =
        typename build_parts_type::template type<build_part_index>;

    if constexpr (is_vec2d_fixed_size_part_v<build_part_type>)
      return build_part_type::size;
    else {
      static_assert(is_vec2d_dynamic_size_part_v<build_part_type>);
      using element_size_type = typename allocator_traits<
          Alloc>::template size_type<build_part_index>;
      const element_size_type max_size = std::get<Index>(init_parts).max_size;

      if constexpr (is_vec2d_rt_part_generic_with_fun_size_v<init_part_type>) {
        namespace ct = boost::callable_traits;

        const element_size_type size = call_dynamic_size_from_callable_at(
            address, std::forward<LineSizeT>(line_index),
            std::get<Index>(init_parts),
            std::make_index_sequence<
                std::tuple_size_v<ct::args_t<init_part_type>>>(),
            std::forward<Sizes>(sizes)...);
        return std::min(size, max_size);

      } else if constexpr (is_vec2d_dynamic_size_from_callable_part_v<
                               build_part_type>) {
        namespace ct = boost::callable_traits;

        const auto &build_part = std::get<build_part_index>(*build_parts);
        const element_size_type size = call_dynamic_size_from_callable_at(
            address, std::forward<LineSizeT>(line_index), build_part,
            std::make_index_sequence<
                std::tuple_size_v<ct::args_t<build_part_type>>>(),
            std::forward<Sizes>(sizes)...);
        return std::min(size, max_size);

      } else {
        static_assert(is_vec2d_rt_part_with_size_v<init_part_type>);
        const element_size_type size = std::get<Index>(init_parts).size;
        return std::min(size, max_size);
      }
    }
  }

  template <std::size_t Index, typename LineSizeT, typename... Sizes>
  [[nodiscard]] inline auto get_init_part_n_elements(LineSizeT &&line_index,
                                                     Sizes &&...sizes) const
      noexcept(noexcept(get_init_part_n_elements_at<Index>(
          address, std::forward<LineSizeT>(line_index),
          std::forward<Sizes>(sizes)...))) {

    return get_init_part_n_elements_at<Index>(
        address, std::forward<LineSizeT>(line_index),
        std::forward<Sizes>(sizes)...);
  }

  template <std::size_t Index>
  [[nodiscard]] inline auto get_init_part_max_size() const noexcept {
    const auto &init_part = std::get<Index>(init_parts);

    using build_part_type = typename build_parts_type::template type<
        builder_type::template get_element_index_for_init<Index>()>;
    using init_part_type = std::decay_t<decltype(init_part)>;

    if constexpr (is_vec2d_fixed_size_part_v<build_part_type>)
      return build_part_type::size;
    else {
      static_assert(is_vec2d_dynamic_size_part_v<build_part_type>);
      static_assert(is_vec2d_rt_part_generic_with_size_v<init_part_type> or
                    is_vec2d_rt_part_generic_with_fun_size_v<init_part_type> or
                    is_vec2d_rt_part_generic_with_max_size_v<init_part_type>);
      return init_part.max_size;
    }
  }

  // TODO noexcept specifier
  template <typename LineIndexT, typename Fun, std::size_t... Idx,
            typename... Sizes>
  [[nodiscard]] inline auto
  call_dynamic_size_from_callable(LineIndexT &&line_index, Fun &&fun,
                                  std::index_sequence<Idx...> idx,
                                  Sizes &&...sizes) const {

    return call_dynamic_size_from_callable_at(
        address, std::forward<LineIndexT>(line_index), std::forward<Fun>(fun),
        std::move(idx), std::forward<Sizes>(sizes)...);
  }

  template <typename LineIndexT, typename Fun, std::size_t... Idx,
            typename... Sizes>
  [[nodiscard]] inline auto call_dynamic_size_from_callable_at(
      base_address_type address, LineIndexT &&line_index, Fun &&fun,
      std::index_sequence<Idx...>, Sizes &&...sizes) const {

    return fun(allocator_traits<Alloc>::template first_pointer_of<Idx>(
        *alloc, address, line_index, sizes...)...);
  }

  template <typename T> inline void destroy_element(T *ptr) const noexcept {
    allocator_traits<Alloc>::destroy(*alloc, ptr);
  }

  template <typename LineIndexT, std::size_t... Idx, typename... Sizes>
  inline void unwind_destroy_line(LineIndexT &&line_index,
                                  std::index_sequence<Idx...>,
                                  Sizes &&...sizes) const noexcept {
    static_assert(sizeof...(Idx) < 2 or std::get<0>(std::array{Idx...}) >
                                            std::get<1>(std::array{Idx...}));
    ([&] { unwind_destroy_element<Idx>(line_index, sizes...); }(), ...);
  }

  template <typename LineIndexT, std::size_t... Idx, typename... Sizes>
  inline void unwind_destroy_line_partial(LineIndexT &&line_index,
                                          std::size_t element_exception_index,
                                          std::index_sequence<Idx...>,
                                          Sizes &&...sizes) const noexcept {
    static_assert(sizeof...(Idx) < 2 or std::get<0>(std::array{Idx...}) >
                                            std::get<1>(std::array{Idx...}));
    (
        [&] {
          if (Idx < element_exception_index) {
            unwind_destroy_element<Idx>(line_index, sizes...);
          }
        }(),
        ...);
  }

  template <std::size_t Index, typename LineIndexT, typename... Sizes>
  inline void unwind_destroy_element(LineIndexT &&line_index,
                                     Sizes &&...sizes) const noexcept {

    auto n_elements = get_init_part_n_elements<Index + 1>(line_index, sizes...);
    auto data_ptr = allocator_traits<Alloc>::template first_pointer_of<Index>(
                        *alloc, address, std::forward<LineIndexT>(line_index),
                        std::forward<Sizes>(sizes)...) +
                    n_elements;

    const std::decay_t<decltype(n_elements)> element_zero(0);
    while (n_elements > element_zero) {
      destroy_element(--data_ptr);
      --n_elements;
    }
  }

  template <std::size_t Index, typename T>
  inline auto call_transform_on_fixed(T *first_pointer) noexcept(
      noexcept(call_transform_on_fixed_impl<Index>(
          std::move(first_pointer),
          std::make_index_sequence<build_parts_type::template type<
              builder_type ::template get_element_index_for_init<Index>()>::
                                       size>()))) {
    constexpr std::size_t build_part_index =
        builder_type::template get_element_index_for_init<Index>();
    static_assert(
        std::is_same_v<T *, typename allocator_traits<Alloc>::template pointer<
                                build_part_index>>);

    using build_part_type =
        typename build_parts_type::template type<build_part_index>;
    static_assert(is_vec2d_fixed_size_part_v<build_part_type>);

    return call_transform_on_fixed_impl<Index>(
        first_pointer, std::make_index_sequence<build_part_type::size>());
  }
};

} // namespace het::detail

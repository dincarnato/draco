#pragma once

#include "../nostd/utility.hpp"

#include <tuple>

namespace het {

template <typename... Ts> struct allocator {
  using allocator_type = allocator;

  template <typename... _Ts> using generic_allocator_type = allocator<_Ts...>;
  using args_arity = std::integral_constant<std::size_t, sizeof...(Ts)>;
  using is_packed = std::false_type;

  template <std::size_t Idx>
  using value_type = std::tuple_element_t<Idx, std::tuple<Ts...>>;

  template <std::size_t> using difference_type = std::ptrdiff_t;
  using lines_difference_type = std::ptrdiff_t;

  template <std::size_t> using size_type = std::size_t;
  using lines_size_type = std::size_t;

  using propagate_on_container_copy_assignment = std::true_type;
  using propagate_on_container_move_assignment = std::true_type;
  using propagate_on_container_swap = std::true_type;
  using is_always_equal = std::true_type;

  template <typename... Sizes>
  constexpr lines_size_type size_of_line(Sizes... sizes) const noexcept {
    static_assert(sizeof...(Sizes) == sizeof...(Ts),
                  "the number of passed Sizes must be the same as the "
                  "templated arguments");

    if constexpr (sizeof...(Sizes) == 0)
      return 0;
    else
      return size_of_line_impl(std::index_sequence_for<Sizes...>(), sizes...);
  }

  template <typename... Sizes>
  [[nodiscard]] inline value_type<0> *allocate(lines_size_type count,
                                               Sizes... sizes) {
    static_assert(sizeof...(Sizes) == args_arity::value,
                  "the number of passed Sizes must be the same as the "
                  "templated arguments");

    return reinterpret_cast<value_type<0> *>(
        ::operator new(size_of_line(sizes...) * count,
                       std::align_val_t(alignof(value_type<0>))));
  }

  template <typename T, typename... Sizes>
  inline void deallocate(T *p, lines_size_type, Sizes...) {
    static_assert(sizeof...(Sizes) == args_arity::value,
                  "the number of passed Sizes must be the same as the "
                  "templated arguments");

    ::operator delete(p, std::align_val_t(alignof(value_type<0>)));
  }

  template <std::size_t Idx, typename... Sizes>
  constexpr lines_difference_type offset_of(Sizes... sizes) const noexcept {
    static_assert(Idx < sizeof...(Ts),
                  "Idx must be less than the size of templated arguments");
    static_assert(Idx <= sizeof...(Sizes),
                  "At least Idx templated arguments are needed");

    if constexpr (sizeof...(Sizes) == 0)
      return 0;
    else
      return offset_of<Idx>(std::index_sequence_for<Sizes...>(), sizes...);
  }

  template <std::size_t Idx, typename... Sizes>
  constexpr value_type<Idx> *first_pointer_of(value_type<0> *first,
                                              lines_size_type line_index,
                                              Sizes... sizes) const noexcept {
    return reinterpret_cast<value_type<Idx> *>(
        reinterpret_cast<char *>(first) + size_of_line(sizes...) * line_index +
        offset_of<Idx>(sizes...));
  }

  template <typename... Rhs>
  constexpr bool operator==(const allocator<Rhs...> &) const noexcept {
    return true;
  }

  template <typename... Rhs>
  constexpr bool operator!=(const allocator<Rhs...> &) noexcept {
    return false;
  }

private:
  template <std::size_t Idx, std::size_t... Indices, typename... Sizes>
  constexpr lines_difference_type offset_of(std::index_sequence<Indices...>,
                                            Sizes... sizes) const noexcept {
    lines_size_type size = 0;
    ((void)(size =
                size +
                ((Indices <= Idx)
                     ? ((size % alignof(std::tuple_element_t<
                                        Indices, std::tuple<Ts...>>) ==
                                 0
                             ? 0
                             : (sizeof(std::tuple_element_t<
                                       Indices, std::tuple<Ts...>>) -
                                (size % alignof(std::tuple_element_t<
                                                Indices, std::tuple<Ts...>>)))))
                     : 0) +
                ((Indices < Idx)
                     ? static_cast<lines_size_type>(
                           nostd::as_signed_upcast(sizeof(
                               std::tuple_element_t<Indices,
                                                    std::tuple<Ts...>>)) *
                           sizes)
                     : 0)),
     ...);
    return static_cast<lines_difference_type>(size);
  }

  template <std::size_t... Idx, typename... Sizes>
  constexpr lines_size_type size_of_line_impl(std::index_sequence<Idx...>,
                                              Sizes... sizes) const noexcept {
    lines_size_type size = 0;
    ((void)(size =
                size +
                ((size % alignof(
                             std::tuple_element_t<Idx, std::tuple<Ts...>>) ==
                          0
                      ? 0
                      : (sizeof(std::tuple_element_t<Idx, std::tuple<Ts...>>) -
                         (size % alignof(std::tuple_element_t<
                                         Idx, std::tuple<Ts...>>))))) +
                static_cast<lines_size_type>(
                    nostd::as_signed_upcast(
                        sizeof(std::tuple_element_t<Idx, std::tuple<Ts...>>)) *
                    sizes)),
     ...);

    if (auto remainder =
            size % alignof(std::tuple_element_t<0, std::tuple<Ts...>>);
        remainder != 0)
      size += alignof(std::tuple_element_t<0, std::tuple<Ts...>>) - remainder;
    return size;
  }
};

static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .size_of_line(1, 2, 2, 1) == 34);
static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .offset_of<0>(1, 2, 2, 1) == 0);
static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .offset_of<0>() == 0);
static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .offset_of<1>(1, 2, 2, 1) == 4);
static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .offset_of<2>(1, 2, 2, 1) == 16);
static_assert(
    allocator<std::uint16_t, std::uint32_t, std::uint64_t, std::uint8_t>{}
        .offset_of<3>(1, 2, 2, 1) == 32);

} // namespace het

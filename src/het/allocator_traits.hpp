#pragma once

#include "allocator.hpp"
#include "../nostd/type_traits.hpp"

#include <memory>
#include <tuple>

namespace het {

template <typename T>
struct is_het_allocator : std::false_type {};

template <typename... Ts>
struct is_het_allocator<allocator<Ts...>> : std::true_type {};

template <typename T>
constexpr bool is_het_allocator_v = is_het_allocator<T>::value;

namespace detail {

  template <typename T>
  using has_pointer = typename T::pointer;
  template <typename T>
  using has_const_pointer = typename T::const_pointer;
  template <typename T>
  using has_void_pointer = typename T::void_pointer;
  template <typename T>
  using has_const_void_pointer = typename T::const_void_pointer;
  template <typename T>
  using has_difference_type = typename T::difference_type;
  template <typename T>
  using has_lines_difference_type = typename T::lines_difference_type;
  template <typename T>
  using has_size_type = typename T::size_type;
  template <typename T>
  using has_lines_size_type = typename T::lines_size_type;
  template <typename T>
  using has_propagate_on_container_copy_assignment =
      typename T::propagate_on_container_copy_assignment;
  template <typename T>
  using has_propagate_on_container_move_assignment =
      typename T::propagate_on_container_move_assignment;
  template <typename T>
  using has_propagate_on_container_swap =
      typename T::propagate_on_container_swap;
  template <typename T>
  using has_is_always_equal = typename T::is_always_equal;
  template <typename T>
  using has_is_packed = std::bool_constant<T::is_packed>;
  template <typename T>
  using has_construct = typename T::construct;
  template <typename T>
  using has_destroy = typename T::destroy;
  template <typename T>
  using has_max_size = typename T::max_size;
  template <typename T>
  using has_select_on_container_copy_construction =
      typename T::select_on_container_copy_construction;

} // namespace detail

template <typename HetAlloc>
struct allocator_traits {
  static_assert(is_het_allocator_v<HetAlloc>);

  using allocator_type = HetAlloc;
  template <typename... Ts>
  using generic_allocator_type =
      typename HetAlloc::template generic_allocator_type<Ts...>;

  using args_arity = typename HetAlloc::args_arity;
  static_assert(args_arity::value > 0,
                "at least one argument type is required");
  using is_packed =
      nostd::detected_or_t<std::false_type, detail::has_is_packed, HetAlloc>;

  template <std::size_t Idx>
  using value_type = typename HetAlloc::template value_type<Idx>;

  template <std::size_t Idx>
  using pointer =
      nostd::detected_or_t<value_type<Idx>*, detail::has_pointer, HetAlloc>;

  using first_pointer = pointer<0>;

  template <std::size_t Idx>
  using const_pointer =
      nostd::detected_or_t<typename std::pointer_traits<pointer<Idx>>::
                               template rebind<const value_type<Idx>>,
                           detail::has_const_pointer, HetAlloc>;

  using first_const_pointer = const_pointer<0>;

  template <std::size_t Idx>
  using void_pointer = nostd::detected_or_t<
      typename std::pointer_traits<pointer<Idx>>::template rebind<void>,
      detail::has_void_pointer, HetAlloc>;

  template <std::size_t Idx>
  using const_void_pointer = nostd::detected_or_t<
      typename std::pointer_traits<pointer<Idx>>::template rebind<const void>,
      detail::has_const_void_pointer, HetAlloc>;

  template <std::size_t Idx>
  using difference_type = nostd::detected_or_t<
      typename std::pointer_traits<pointer<Idx>>::difference_type,
      detail::has_difference_type, HetAlloc>;

  using lines_difference_type = nostd::detected_or_t<
      typename std::pointer_traits<first_pointer>::difference_type,
      detail::has_lines_difference_type, HetAlloc>;

  template <std::size_t Idx>
  using size_type =
      nostd::detected_or_t<std::make_unsigned_t<difference_type<Idx>>,
                           detail::has_size_type, HetAlloc>;

  using lines_size_type =
      nostd::detected_or_t<std::make_unsigned_t<lines_difference_type>,
                           detail::has_lines_size_type, HetAlloc>;

  using propagate_on_container_copy_assignment =
      nostd::detected_or_t<std::false_type,
                           detail::has_propagate_on_container_copy_assignment,
                           HetAlloc>;

  using propagate_on_container_move_assignment =
      nostd::detected_or_t<std::false_type,
                           detail::has_propagate_on_container_move_assignment,
                           HetAlloc>;

  using propagate_on_container_swap =
      nostd::detected_or_t<std::false_type,
                           detail::has_propagate_on_container_swap, HetAlloc>;

  using is_always_equal =
      nostd::detected_or_t<typename std::is_empty<HetAlloc>::type,
                           detail::has_is_always_equal, HetAlloc>;

  template <typename... Sizes>
  [[nodiscard]] static inline first_pointer
  allocate(HetAlloc& alloc, lines_size_type count, Sizes... sizes) {
    return check_and_run(std::index_sequence_for<Sizes...>(),
                         [&] { return alloc.allocate(count, sizes...); },
                         sizes...);
  }

  template <typename... Sizes>
  [[nodiscard]] static inline first_pointer
  allocate(HetAlloc& alloc, Sizes... sizes, lines_size_type count,
           const_void_pointer<0> hint) {
    return check_and_run(std::index_sequence_for<Sizes...>(),
                         [&] { return alloc.allocate(sizes..., count, hint); },
                         sizes...);
  }

  template <typename... Sizes>
  static inline void
  deallocate(HetAlloc& alloc, first_pointer p, lines_size_type count,
             Sizes... sizes) {
    check_and_run(std::index_sequence_for<Sizes...>(),
                  [&] { alloc.deallocate(p, count, sizes...); }, sizes...);
  }

  template <typename T, typename... Args>
  static inline void
  construct(HetAlloc& alloc, T* p, Args&&... args) {
    if constexpr (nostd::is_detected_v<detail::has_construct, HetAlloc>)
      alloc.construct(p, std::forward<Args>(args)...);
    else
      ::new (p) T(std::forward<Args>(args)...);
  }

  template <typename T, typename... Args>
  static inline void
  destroy(HetAlloc& alloc, T* p) {
    if constexpr (nostd::is_detected_v<detail::has_destroy, HetAlloc>)
      alloc.template destroy(p);
    else
      p->~T();
  }

  static inline HetAlloc
  select_on_container_copy_construction(const HetAlloc& a) {
    if constexpr (nostd::is_detected_v<
                      detail::has_select_on_container_copy_construction,
                      HetAlloc>)
      return a.select_on_container_copy_construction();
    else
      return a;
  }

  template <typename... Sizes>
  static constexpr lines_size_type
  size_of_line(const HetAlloc& alloc, Sizes... sizes) noexcept {
    return check_and_run(std::index_sequence_for<Sizes...>(),
                         [&] { return alloc.size_of_line(sizes...); },
                         sizes...);
  }

  template <std::size_t Idx, typename... Sizes>
  static constexpr lines_difference_type
  offset_of(const HetAlloc& alloc, Sizes... sizes) noexcept {
    return alloc.template offset_of<Idx>(sizes...);
  }

  template <std::size_t Idx, typename... Sizes>
  static constexpr value_type<Idx>*
  first_pointer_of(const HetAlloc& alloc, first_pointer first,
                   lines_size_type line_index, Sizes... sizes) noexcept {
    return check_and_run(std::index_sequence_for<Sizes...>(),
                         [&] {
                           return alloc.template first_pointer_of<Idx>(
                               first, line_index, sizes...);
                         },
                         sizes...);
  }

private:
  template <typename Fun, typename std::size_t... Idx, typename... Sizes>
  static inline decltype(auto)
  check_and_run(std::index_sequence<Idx...>, Fun fun, Sizes...) {
    static_assert(sizeof...(Sizes) == args_arity::value,
                  "the size for every type must be specified");
    static_assert(
        (std::is_convertible_v<
             std::tuple_element_t<Idx, std::tuple<std::decay_t<Sizes>...>>,
             size_type<Idx>> &&
         ...),
        "all Sizes type must correspond to the allocator size_type-s");
    return fun();
  }
};

} // namespace het

#pragma once

#include "windows_merger_windows.hpp"

#include "windows_merger_exceptions.hpp"
#include "windows_merger_window_accessor.hpp"
#include "windows_merger_window_accessor_coverages_accessor.hpp"
#include "windows_merger_window_base.hpp"
#include "windows_merger_window_base_accessor.hpp"
#include "windows_merger_window_coverages_accessor.hpp"
#include "windows_merger_window_weights_accessor.hpp"
#include "windows_merger_windows_iterator.hpp"

#include <algorithm>
#include <cassert>
#include <range/v3/iterator/concepts.hpp>

namespace windows_merger {

inline WindowsMergerWindows::WindowsMergerWindows(
    clusters_size_type n_clusters) noexcept
    : n_clusters(n_clusters) {}

inline WindowsMergerWindows::WindowsMergerWindows(
    clusters_size_type n_clusters, bases_size_type n_bases) noexcept
    : n_clusters(n_clusters), n_bases_capacity(n_bases) {}

inline WindowsMergerWindows::WindowsMergerWindows(
    WindowsMergerWindows&& other) noexcept
    : n_clusters(other.n_clusters),
      n_bases_capacity(std::exchange(other.n_bases_capacity, 0)),
      n_windows_capacity(std::exchange(other.n_windows_capacity, 0)),
      n_windows(std::exchange(other.n_windows, 0)),
      windows_data(std::exchange(other.windows_data, nullptr)),
      coverages_data(std::exchange(other.coverages_data, nullptr)) {}

inline WindowsMergerWindows&
WindowsMergerWindows::operator=(WindowsMergerWindows&& other) noexcept {
  destroy_and_deallocate();

  n_clusters = other.n_clusters;
  n_bases_capacity = std::exchange(other.n_bases_capacity, 0);
  n_windows_capacity = std::exchange(other.n_windows_capacity, 0);
  n_windows = std::exchange(other.n_windows, 0);
  windows_data = std::exchange(other.windows_data, nullptr);
  coverages_data = std::exchange(other.coverages_data, nullptr);

  return *this;
}

inline WindowsMergerWindows::WindowsMergerWindows(clusters_size_type n_clusters,
                                                  bases_size_type n_bases,
                                                  windows_size_type n_windows)
    : n_clusters(n_clusters), n_bases_capacity(n_bases),
      n_windows_capacity(n_windows), n_windows(n_windows) {
  if (n_clusters > 0 and n_bases > 0 and n_windows_capacity > 0)
    allocate_impl();
}

inline void
WindowsMergerWindows::allocate_impl() noexcept(false) {
  assert(n_clusters > 0);
  assert(n_bases_capacity > 0);

  const auto n_clusters_bases = n_clusters * n_bases_capacity;
  windows_data = windows_allocator_traits::allocate(*this, n_windows_capacity,
                                                    2, n_clusters_bases);
  auto windows_unwinder = [&] {
    try {
      return windows_builder.build_at(windows_data, *this)
          .with_lines(n_windows_capacity)
          .start()
          .construct_default()
          .next()
          .with_max_size(static_cast<bases_size_type>(n_clusters_bases))
          .with_fun_size([n_clusters =
                              n_clusters](const bases_size_type* start_end) {
            assert(start_end[0] <= start_end[1]);
            return static_cast<bases_size_type>((start_end[1] - start_end[0]) *
                                                n_clusters);
          })
          .construct_default()
          .done();
    } catch (...) {
      windows_allocator_traits::deallocate(
          *this, windows_data, n_windows_capacity, 2, n_clusters_bases);
      throw;
    }
  }();

  coverages_data = [&] {
    try {
      return coverages_allocator_traits::allocate(*this, n_windows_capacity, 1,
                                                  n_bases_capacity);
    } catch (...) {
      windows_unwinder.unwind();
      windows_allocator_traits::deallocate(
          *this, windows_data, n_windows_capacity, 2, n_clusters_bases);
      throw;
    }
  }();

  try {
    coverages_builder.build_at(coverages_data, *this)
        .with_lines(n_windows_capacity)
        .start()
        .construct_default()
        .next()
        .with_max_size(n_bases_capacity)
        .construct_default()
        .done();
  } catch (...) {
    coverages_allocator_traits::deallocate(
        *this, coverages_data, n_windows_capacity, 1, n_bases_capacity);
    windows_unwinder.unwind();
    windows_allocator_traits::deallocate(
        *this, windows_data, n_windows_capacity, 2, n_clusters_bases);
    throw;
  }
}

inline WindowsMergerWindows::~WindowsMergerWindows() noexcept {
  destroy_and_deallocate();
}

inline void
WindowsMergerWindows::destroy_and_deallocate() noexcept {
  static_assert(std::is_trivially_destructible_v<weight_type>);
  if (n_clusters > 0 and n_bases_capacity > 0 and n_windows_capacity > 0) {
    windows_builder.destroy(windows_data, *this)
        .with_lines(n_windows)
        .start()
        .next()
        .with_max_size(
            static_cast<bases_size_type>(n_bases_capacity * n_clusters))
        .with_fun_size([n_clusters =
                            n_clusters](const bases_size_type* start_end) {
          assert(start_end[0] <= start_end[1]);
          return static_cast<bases_size_type>((start_end[1] - start_end[0]) *
                                              n_clusters);
        })
        .done();
    windows_allocator_traits::deallocate(*this, windows_data,
                                         n_windows_capacity, 2,
                                         n_clusters * n_bases_capacity);

    coverages_builder.destroy(coverages_data, *this)
        .with_lines(n_windows)
        .start()
        .next()
        .with_max_size(n_bases_capacity)
        .done();
    coverages_allocator_traits::deallocate(
        *this, coverages_data, n_windows_capacity, 1, n_bases_capacity);
  }
}

inline bool
WindowsMergerWindows::operator==(const WindowsMergerWindows& other) const
    noexcept {
  return static_cast<const windows_allocator&>(*this) ==
             static_cast<const windows_allocator&>(other) &&
         static_cast<const coverages_allocator&>(*this) ==
             static_cast<const coverages_allocator&>(other) &&
         n_clusters == other.n_clusters && n_bases_capacity &&
         other.n_bases_capacity &&
         n_windows_capacity == other.n_windows_capacity &&
         n_windows == other.n_windows && windows_data == other.windows_data &&
         coverages_data == other.coverages_data;
}

inline bool
WindowsMergerWindows::operator!=(const WindowsMergerWindows& other) const
    noexcept {
  return !operator==(other);
}

inline auto
WindowsMergerWindows::clusters_size() const noexcept -> clusters_size_type {
  return n_clusters;
}

inline auto
WindowsMergerWindows::bases_capacity() const noexcept -> bases_size_type {
  return n_bases_capacity;
}

inline auto
WindowsMergerWindows::windows_size() const noexcept -> windows_size_type {
  return n_windows;
}

inline auto
WindowsMergerWindows::windows_capacity() const noexcept -> windows_size_type {
  return n_windows_capacity;
}

inline auto WindowsMergerWindows::operator[](windows_size_type index) &
    noexcept -> window_accessor {
  return window_accessor{*this, index};
}

inline auto WindowsMergerWindows::operator[](windows_size_type index) const &
    noexcept -> const_window_accessor {
  return const_window_accessor{*this, index};
}

inline auto WindowsMergerWindows::operator[](windows_size_type index) &&
    noexcept -> move_window_accessor {
  return move_window_accessor{*this, index};
}

template <std::size_t Index>
inline auto
WindowsMergerWindows::get_window_first_pointer(
    windows_size_type window_index) const noexcept {
  static_assert(Index < 2);
  return windows_allocator_traits::first_pointer_of<Index>(
      *this, windows_data, window_index, 2, n_clusters * n_bases_capacity);
}

template <std::size_t Index>
inline auto
WindowsMergerWindows::get_window_first_pointer(
    windows_size_type window_index) noexcept {
  return const_cast<typename windows_allocator_traits::template pointer<Index>>(
      const_cast<const WindowsMergerWindows&>(*this)
          .get_window_first_pointer<Index>(window_index));
}

template <std::size_t Index>
inline auto
WindowsMergerWindows::get_coverage_first_pointer(
    windows_size_type window_index) const noexcept {
  return coverages_allocator_traits::first_pointer_of<Index>(
      *this, coverages_data, window_index, 1, n_bases_capacity);
}

template <std::size_t Index>
inline auto
WindowsMergerWindows::get_coverage_first_pointer(
    windows_size_type window_index) noexcept {
  return const_cast<
      typename coverages_allocator_traits::template pointer<Index>>(
      const_cast<const WindowsMergerWindows&>(*this)
          .template get_coverage_first_pointer<Index>(window_index));
}

template <typename Base>
auto
WindowsMergerWindows::base_emplace_back(windows_size_type window_index,
                                        Base&& base) noexcept(false)
    -> base_accessor {
  using new_base_type = std::decay_t<Base>;
  static_assert(traits_type::template is_window_baselike_v<new_base_type>);

  assert(n_clusters > 0);
  assert(n_windows_capacity > 0);

  assert(window_index < n_windows);
  if (base.clusters_size() != n_clusters)
    throw exceptions::InvalidClustersSize(
        "base object has a number of clusters different from the windows");

  if (n_bases_capacity == 0) {
    n_bases_capacity = 1;
    allocate_impl();
  }

  const auto window_start_end = get_window_first_pointer<0>(window_index);
  assert(window_start_end[0] <= window_start_end[1]);
  const auto window_size =
      static_cast<bases_size_type>(window_start_end[1] - window_start_end[0]);

  if (window_size == n_bases_capacity)
    reshape(n_bases_capacity * 2, WindowsMergerReserve(n_windows_capacity));

  constexpr bool weights_build_is_noexcept =
      not std::is_lvalue_reference_v<Base> and
      noexcept(weight_type(base.weight(std::declval<clusters_size_type>()))) and
      noexcept(std::declval<new_base_type&>().weight(
                   std::declval<clusters_size_type>()) =
                   std::declval<weight_type&&>());

  auto new_base_ptr =
      get_window_first_pointer<1>(window_index) + window_size * n_clusters;
  for (clusters_size_type cluster_index = 0; cluster_index < n_clusters;
       ++cluster_index) {
    try {
      if constexpr (weights_build_is_noexcept) {
        windows_allocator_traits::construct(
            *this, new_base_ptr + cluster_index,
            std::move(base.weight(cluster_index)));
      } else {
        windows_allocator_traits::construct(*this, new_base_ptr + cluster_index,
                                            base.weight(cluster_index));
      }
    } catch (...) {
      while (cluster_index > 0) {
        --cluster_index;

        if constexpr (weights_build_is_noexcept)
          base.weight(cluster_index) = std::move(new_base_ptr[cluster_index]);

        windows_allocator_traits::destroy(*this, new_base_ptr + cluster_index);
      }

      throw;
    }
  }

  constexpr bool coverages_build_is_noexcept =
      not std::is_lvalue_reference_v<Base> and
      noexcept(coverage_type(base.coverage())) and
      noexcept(std::declval<new_base_type&>().coverage() =
                   std::declval<coverage_type&&>());

  auto new_coverage_ptr =
      get_coverage_first_pointer<1>(window_index) + window_size;
  try {
    if constexpr (coverages_build_is_noexcept)
      coverages_allocator_traits::construct(*this, new_coverage_ptr,
                                            std::move(base.coverage()));
    else
      coverages_allocator_traits::construct(*this, new_coverage_ptr,
                                            base.coverage());
  } catch (...) {
    if constexpr (coverages_build_is_noexcept)
      base.coverage() = std::move(*new_coverage_ptr);

    for (clusters_size_type cluster_index = n_clusters; cluster_index > 0;) {
      --cluster_index;

      if constexpr (weights_build_is_noexcept)
        base.weight(cluster_index) = std::move(new_base_ptr[cluster_index]);

      windows_allocator_traits::destroy(*this, new_base_ptr + cluster_index);
    }

    throw;
  }

  ++get_window_first_pointer<0>(window_index)[1];
  ++*get_coverage_first_pointer<0>(window_index);

  return base_accessor(*this, window_index, window_size);
}

template <template <typename> typename WindowsReshaper>
void
WindowsMergerWindows::reshape(
    bases_size_type bases_capacity,
    WindowsReshaper<windows_size_type>&& windows_reshaper) noexcept(false) {
  using windows_reshaper_type =
      std::decay_t<WindowsReshaper<windows_size_type>>;
  static_assert(traits_type::template is_reshaper_arg_v<windows_reshaper_type>,
                "windows_reshaper must be a WindowsMergerResize<T> or a "
                "WindowsMergerReserve<T> type");
  static_assert(std::is_convertible_v<typename windows_reshaper_type::type,
                                      windows_size_type>,
                "The type stored in the windows_reshaper must be convertible "
                "to a windows_size_type");
  assert(n_clusters > 0);

  auto new_windows_size = [&] {
    if constexpr (traits_type::template is_reserve_reshaper_v<
                      windows_reshaper_type>) {
      return n_windows;
    } else {
      return static_cast<windows_size_type>(windows_reshaper.size);
    }
  }();

  auto new_windows_capacity = [&] {
    if constexpr (traits_type::template is_reserve_reshaper_v<
                      windows_reshaper_type>) {
      return std::max({static_cast<windows_size_type>(windows_reshaper.size),
                       n_windows_capacity, new_windows_size});
    } else {
      return std::max(n_windows_capacity, new_windows_size);
    }
  }();

  if (n_windows_capacity == 0 or n_bases_capacity == 0) {
    if (new_windows_capacity != 0)
      n_windows_capacity = new_windows_capacity;

    if (bases_capacity != 0)
      n_bases_capacity = bases_capacity;

    if (n_windows_capacity != 0 and n_bases_capacity != 0)
      allocate_impl();
    else
      return;
  }

  auto new_bases_capacity =
      std::max(static_cast<bases_size_type>(bases_capacity), n_bases_capacity);

  if (new_bases_capacity == n_bases_capacity and
      new_windows_capacity == n_windows_capacity and
      new_windows_size == n_windows)
    return;

  if (new_bases_capacity != n_bases_capacity or
      new_windows_capacity > n_windows_capacity) {
    const auto new_clusters_bases =
        static_cast<bases_size_type>(new_bases_capacity * n_clusters);
    const auto old_clusters_bases =
        static_cast<bases_size_type>(n_bases_capacity * n_clusters);

    auto new_windows_data = windows_allocator_traits::allocate(
        *this, new_windows_capacity, 2, new_clusters_bases);

    auto windows_unwinder = [&] {
      try {
        return windows_builder.move_to(new_windows_data, *this)
            .with_lines(new_windows_size)
            .start()
            .next()
            .with_max_size(new_clusters_bases)
            .with_fun_size(
                [n_clusters = n_clusters](const bases_size_type* start_end) {
                  assert(start_end[0] <= start_end[1]);
                  return static_cast<bases_size_type>(
                      (start_end[1] - start_end[0]) * n_clusters);
                })
            .remaining()
            .construct_default()
            .next()
            .with_max_size(new_clusters_bases)
            .with_fun_size(
                [n_clusters = n_clusters](const bases_size_type* start_end) {
                  assert(start_end[0] <= start_end[1]);
                  return static_cast<bases_size_type>(
                      (start_end[1] - start_end[0]) * n_clusters);
                })
            .construct_default()
            .done()
            .with_lines(n_windows)
            .from(windows_data)
            .start()
            .next()
            .with_max_size(old_clusters_bases)
            .with_fun_size(
                [n_clusters = n_clusters](const bases_size_type* start_end) {
                  assert(start_end[0] <= start_end[1]);
                  return static_cast<bases_size_type>(
                      (start_end[1] - start_end[0]) * n_clusters);
                })
            .done();
      } catch (...) {
        windows_allocator_traits::deallocate(*this, new_windows_data,
                                             new_windows_capacity, 2,
                                             new_clusters_bases);
        throw;
      }
    }();

    auto new_coverages_data = [&] {
      try {
        return coverages_allocator_traits::allocate(*this, new_windows_capacity,
                                                    1, new_bases_capacity);
      } catch (...) {
        windows_unwinder.unwind();
        windows_allocator_traits::deallocate(*this, new_windows_data,
                                             new_windows_capacity, 2,
                                             new_clusters_bases);
        throw;
      }
    }();

    try {
      coverages_builder.move_to(new_coverages_data, *this)
          .with_lines(new_windows_size)
          .start()
          .next()
          .with_max_size(new_bases_capacity)
          .remaining()
          .construct_default()
          .next()
          .with_max_size(new_bases_capacity)
          .construct_default()
          .done()
          .with_lines(n_windows)
          .from(coverages_data)
          .start()
          .next()
          .with_max_size(n_bases_capacity)
          .done();
    } catch (...) {
      coverages_allocator_traits::deallocate(*this, new_coverages_data,
                                             new_windows_capacity, 1,
                                             new_bases_capacity);
      windows_unwinder.unwind();
      windows_allocator_traits::deallocate(
          *this, new_windows_data, new_windows_capacity, 2, new_clusters_bases);
      throw;
    }

    windows_builder.destroy(windows_data, *this)
        .with_lines(n_windows)
        .start()
        .next()
        .with_max_size(old_clusters_bases)
        .with_fun_size([n_clusters =
                            n_clusters](const bases_size_type* start_end) {
          assert(start_end[0] <= start_end[1]);
          return static_cast<bases_size_type>((start_end[1] - start_end[0]) *
                                              n_clusters);
        })
        .done();
    windows_allocator_traits::deallocate(
        *this, windows_data, n_windows_capacity, 2, old_clusters_bases);

    coverages_builder.destroy(coverages_data, *this)
        .with_lines(n_windows)
        .start()
        .next()
        .with_max_size(n_bases_capacity)
        .done();
    coverages_allocator_traits::deallocate(
        *this, coverages_data, n_windows_capacity, 1, n_bases_capacity);

    windows_data = new_windows_data;
    coverages_data = new_coverages_data;
  } else if constexpr (WindowsMergerTraits::template is_resize_reshaper_v<
                           windows_reshaper_type>) {
    static_assert(not WindowsMergerTraits::template is_reserve_reshaper_v<
                  windows_reshaper_type>);
    const auto clusters_bases =
        static_cast<bases_size_type>(n_bases_capacity * n_clusters);

    if (new_windows_size > n_windows) {
      auto windows_unwinder =
          windows_builder.build_at(windows_data, *this)
              .with_lines(new_windows_size)
              .skip_lines(n_windows)
              .start()
              .construct_default()
              .next()
              .with_max_size(clusters_bases)
              .with_fun_size(
                  [n_clusters = n_clusters](const bases_size_type* start_end) {
                    assert(start_end[0] <= start_end[1]);
                    return static_cast<bases_size_type>(
                        (start_end[1] - start_end[0]) * n_clusters);
                  })
              .construct_default()
              .done();

      try {
        coverages_builder.build_at(coverages_data, *this)
            .with_lines(new_windows_size)
            .skip_lines(n_windows)
            .start()
            .construct_default()
            .next()
            .with_max_size(n_bases_capacity)
            .construct_default()
            .done();
      } catch (...) {
        windows_unwinder.unwind();
        throw;
      }
    } else if (new_windows_size < n_windows) {
      windows_builder.destroy(windows_data, *this)
          .with_lines(n_windows)
          .skip_lines(new_windows_size)
          .start()
          .next()
          .with_max_size(clusters_bases)
          .with_fun_size([n_clusters =
                              n_clusters](const bases_size_type* start_end) {
            assert(start_end[0] <= start_end[1]);
            return static_cast<bases_size_type>((start_end[1] - start_end[0]) *
                                                n_clusters);
          })
          .done();

      coverages_builder.destroy(coverages_data, *this)
          .with_lines(n_windows)
          .skip_lines(new_windows_size)
          .start()
          .next()
          .with_max_size(n_bases_capacity)
          .done();
    }
  }

  n_bases_capacity = new_bases_capacity;
  n_windows_capacity = new_windows_capacity;
  n_windows = new_windows_size;
}

inline auto
WindowsMergerWindows::emplace_back() noexcept(false) -> window_accessor {
  const auto new_window_index = n_windows;
  if (n_bases_capacity > 0)
    reshape(n_bases_capacity, resizer_type(n_windows + 1));
  else {
    if (n_windows_capacity == n_windows) {
      if (n_windows_capacity == 0)
        n_windows_capacity = 1;
      else
        n_windows_capacity *= 2;
    }
    ++n_windows;
  }
  return window_accessor(*this, new_window_index);
}

inline void
WindowsMergerWindows::set_begin_index(
    windows_size_type window_index,
    bases_size_type new_begin_index) noexcept(false) {
  assert(window_index < n_windows);

  if (n_bases_capacity == 0)
    reshape(1, reserver_type(n_windows_capacity));

  const auto start_end = get_window_first_pointer<0>(window_index);
  const bases_size_type new_end_index =
      new_begin_index +
      static_cast<bases_size_type>(start_end[1] - start_end[0]);

  start_end[0] = new_begin_index;
  start_end[1] = new_end_index;
}

inline bool
WindowsMergerWindows::is_allocated() const noexcept {
  return n_bases_capacity > 0 and n_windows_capacity > 0;
}

inline auto
WindowsMergerWindows::front() noexcept -> window_accessor {
  assert(n_windows > 0);
  return window_accessor(*this, 0);
}

inline auto
WindowsMergerWindows::front() const noexcept -> const_window_accessor {
  assert(n_windows > 0);
  return const_window_accessor(*this, 0);
}

inline auto
WindowsMergerWindows::back() noexcept -> window_accessor {
  assert(n_windows > 0);
  return window_accessor(*this, n_windows - 1);
}

inline auto
WindowsMergerWindows::back() const noexcept -> const_window_accessor {
  assert(n_windows > 0);
  return const_window_accessor(*this, n_windows - 1);
}

inline void
WindowsMergerWindows::destroy_range(windows_size_type window_index,
                                    bases_size_type first,
                                    bases_size_type last) noexcept {
  {
    auto weights_pointer = get_window_first_pointer<1>(window_index);
    const auto weights_ptr_end = weights_pointer + last;
    for (auto weight_ptr = weights_pointer + first;
         weight_ptr < weights_ptr_end; ++weight_ptr) {
      windows_allocator_traits::destroy(*this, weight_ptr);
    }
  }

  {
    auto coverages_pointer = get_coverage_first_pointer<1>(window_index);
    const auto coverages_ptr_end = coverages_pointer + last;
    for (auto coverage_ptr = coverages_pointer + first;
         coverage_ptr < coverages_ptr_end; ++coverage_ptr) {
      coverages_allocator_traits::destroy(*this, coverage_ptr);
    }
  }
}

inline void
WindowsMergerWindows::copy_construct_from_range(
    const WindowsMergerWindows& from_windows,
    windows_size_type from_window_index, bases_size_type from_first,
    bases_size_type from_last, windows_size_type to_window_index,
    bases_size_type
        to_first) noexcept(noexcept(std::
                                        is_nothrow_copy_constructible_v<
                                            weight_type>) and
                           noexcept(std::is_nothrow_copy_constructible_v<
                                    coverage_type>)) {
  assert(to_window_index < n_windows);
  assert(from_window_index < from_windows.n_windows);

  assert(from_first <= from_windows[from_window_index].size());
  assert(from_last <= from_windows[from_window_index].size());
  assert(to_first + (from_last - from_first) <= n_bases_capacity);

  {
    auto from_weights_pointer =
        from_windows.template get_window_first_pointer<1>(from_window_index);
    auto to_weights_pointer = get_window_first_pointer<1>(to_window_index);

    auto from_first_weight = from_first * n_clusters;
    auto from_last_weight = from_last * n_clusters;
    auto to_first_weight = to_first * n_clusters;
    for (; from_first_weight < from_last_weight;
         ++from_first_weight, ++to_first_weight)
      windows_allocator_traits::construct(
          *this, to_weights_pointer + to_first_weight,
          from_weights_pointer[from_first_weight]);
  }

  {
    auto from_coverages_pointer =
        from_windows.template get_coverage_first_pointer<1>(from_window_index);
    auto to_coverages_pointer = get_coverage_first_pointer<1>(to_window_index);
    for (; from_first < from_last; ++from_first, ++to_first)
      windows_allocator_traits::construct(*this,
                                          to_coverages_pointer + to_first,
                                          from_coverages_pointer[from_first]);
  }
}

inline void
WindowsMergerWindows::copy_construct_from_range(
    const WindowsMergerWindow& from_windows, bases_size_type from_first,
    bases_size_type from_last, windows_size_type to_window_index,
    bases_size_type
        to_first) noexcept(noexcept(std::
                                        is_nothrow_copy_constructible_v<
                                            weight_type>) and
                           noexcept(std::is_nothrow_copy_constructible_v<
                                    coverage_type>)) {
  assert(to_window_index < n_windows);

  assert(from_first <= from_windows.size());
  assert(from_last <= from_windows.size());
  assert(to_first + (from_last - from_first) <= n_bases_capacity);

  {
    auto&& linear_weights = from_windows.linear_weights();
    auto to_weights_pointer = get_window_first_pointer<1>(to_window_index);

    assert(from_first <=
           std::numeric_limits<bases_size_type>::max() / n_clusters);
    assert(from_last <=
           std::numeric_limits<bases_size_type>::max() / n_clusters);

    auto from_first_weight =
        static_cast<bases_size_type>(from_first * n_clusters);
    auto from_last_weight =
        static_cast<bases_size_type>(from_last * n_clusters);
    auto to_first_weight = to_first * n_clusters;
    for (; from_first_weight < from_last_weight;
         ++from_first_weight, ++to_first_weight)
      windows_allocator_traits::construct(*this,
                                          to_weights_pointer + to_first_weight,
                                          linear_weights[from_first_weight]);
  }

  {
    auto&& coverages = from_windows.coverages();
    auto to_coverages_pointer = get_coverage_first_pointer<1>(to_window_index);
    for (; from_first < from_last; ++from_first, ++to_first)
      windows_allocator_traits::construct(
          *this, to_coverages_pointer + to_first, coverages[from_first]);
  }
}

inline void
WindowsMergerWindows::move_construct_from_range(
    const WindowsMergerWindows& from_windows,
    windows_size_type from_window_index, bases_size_type from_first,
    bases_size_type from_last, windows_size_type to_window_index,
    bases_size_type
        to_first) noexcept(noexcept(std::
                                        is_nothrow_move_constructible_v<
                                            weight_type>) and
                           noexcept(std::is_nothrow_move_constructible_v<
                                    coverage_type>)) {
  assert(to_window_index < n_windows);
  assert(from_window_index < from_windows.n_windows);

  assert(from_first <= from_windows[from_window_index].size());
  assert(from_last <= from_windows[from_window_index].size());
  assert(to_first + (from_last - from_first) <= n_bases_capacity);

  {
    auto from_weights_pointer =
        from_windows.template get_window_first_pointer<1>(from_window_index);
    auto to_weights_pointer = get_window_first_pointer<1>(to_window_index);

    auto from_first_weight = from_first * n_clusters;
    auto from_last_weight = from_last * n_clusters;
    auto to_first_weight = to_first * n_clusters;
    for (; from_first_weight < from_last_weight;
         ++from_first_weight, ++to_first_weight)
      windows_allocator_traits::construct(
          *this, to_weights_pointer + to_first_weight,
          std::move(from_weights_pointer[from_first_weight]));
  }

  {
    auto from_coverages_pointer =
        from_windows.template get_coverage_first_pointer<1>(from_window_index);
    auto to_coverages_pointer = get_coverage_first_pointer<1>(to_window_index);
    for (; from_first < from_last; ++from_first, ++to_first)
      windows_allocator_traits::construct(
          *this, to_coverages_pointer + to_first,
          std::move(from_coverages_pointer[from_first]));
  }
}

inline void
WindowsMergerWindows::move_construct_from_range(
    WindowsMergerWindow&& from_windows, bases_size_type from_first,
    bases_size_type from_last, windows_size_type to_window_index,
    bases_size_type
        to_first) noexcept(noexcept(std::
                                        is_nothrow_move_constructible_v<
                                            weight_type>) and
                           noexcept(std::is_nothrow_move_constructible_v<
                                    coverage_type>)) {
  assert(to_window_index < n_windows);
  assert(to_first + (from_last - from_first) <= n_bases_capacity);

  {
    auto&& linear_weights = from_windows.linear_weights();
    auto to_weights_pointer = get_window_first_pointer<1>(to_window_index);

    assert(from_first <=
           std::numeric_limits<bases_size_type>::max() / n_clusters);
    assert(from_last <=
           std::numeric_limits<bases_size_type>::max() / n_clusters);

    auto from_first_weight =
        static_cast<bases_size_type>(from_first * n_clusters);
    auto from_last_weight =
        static_cast<bases_size_type>(from_last * n_clusters);
    auto to_first_weight = to_first * n_clusters;
    for (; from_first_weight < from_last_weight;
         ++from_first_weight, ++to_first_weight)
      windows_allocator_traits::construct(
          *this, to_weights_pointer + to_first_weight,
          std::move(linear_weights[from_first_weight]));
  }

  {
    auto&& coverages = from_windows.coverages();
    auto to_coverages_pointer = get_coverage_first_pointer<1>(to_window_index);
    for (; from_first < from_last; ++from_first, ++to_first)
      windows_allocator_traits::construct(*this,
                                          to_coverages_pointer + to_first,
                                          std::move(coverages[from_first]));
  }
}

inline auto
    WindowsMergerWindows::begin() &
    noexcept -> iterator {
  return iterator(*this);
}

inline auto
    WindowsMergerWindows::end() &
    noexcept -> iterator {
  return iterator(*this, n_windows);
}

inline auto
    WindowsMergerWindows::begin() const
    & noexcept -> const_iterator {
  return const_iterator(*this);
}

inline auto
    WindowsMergerWindows::end() const
    & noexcept -> const_iterator {
  return const_iterator(*this, n_windows);
}

inline auto
    WindowsMergerWindows::begin() &&
    noexcept -> move_iterator {
  return move_iterator(*this);
}

inline auto
    WindowsMergerWindows::end() &&
    noexcept -> move_iterator {
  return move_iterator(*this, n_windows);
}

inline auto
    WindowsMergerWindows::rbegin() &
    noexcept -> reverse_iterator {
  return reverse_iterator(iterator(*this, n_windows));
}

inline auto
    WindowsMergerWindows::rend() &
    noexcept -> reverse_iterator {
  return reverse_iterator(iterator(*this));
}

inline auto
    WindowsMergerWindows::rbegin() const
    & noexcept -> const_reverse_iterator {
  return const_reverse_iterator(const_iterator(*this, n_windows));
}

inline auto
    WindowsMergerWindows::rend() const
    & noexcept -> const_reverse_iterator {
  return const_reverse_iterator(const_iterator(*this));
}

inline auto
    WindowsMergerWindows::rbegin() &&
    noexcept -> move_reverse_iterator {
  return move_reverse_iterator(move_iterator(*this, n_windows));
}

inline auto
    WindowsMergerWindows::rend() &&
    noexcept -> move_reverse_iterator {
  return move_reverse_iterator(move_iterator(*this));
}

template <typename T, typename U>
auto
operator==(T&& t, U&& u) noexcept -> std::enable_if_t<
    WindowsMergerTraits::template is_window_baselike_v<std::decay_t<T>> and
        WindowsMergerTraits::template is_window_baselike_v<std::decay_t<U>>,
    bool> {

  return t.coverage() == u.coverage() and
         t.clusters_size() == u.clusters_size() and
         [t_weights = t.weights(), u_weights = u.weights()] {
           return ranges::equal(t_weights, u_weights);
         }();
}

template <typename T, typename U>
auto
operator!=(T&& t, U&& u) noexcept -> std::enable_if_t<
    WindowsMergerTraits::template is_window_baselike_v<std::decay_t<T>> and
        WindowsMergerTraits::template is_window_baselike_v<std::decay_t<U>>,
    bool> {
  return not(t == u);
}

template <typename T, typename U>
auto
operator==(T&& t, U&& u) noexcept -> std::enable_if_t<
    WindowsMergerTraits::template is_windowlike_v<std::decay_t<T>> and
        WindowsMergerTraits::template is_windowlike_v<std::decay_t<U>>,
    bool> {

  return t.clusters_size() == u.clusters_size() and
         t.begin_index() == u.begin_index() and t.size() == u.size() and
         ranges::equal(t, u);
}

template <typename T, typename U>
auto
operator!=(T&& t, U&& u) noexcept -> std::enable_if_t<
    WindowsMergerTraits::template is_windowlike_v<std::decay_t<T>> and
        WindowsMergerTraits::template is_windowlike_v<std::decay_t<U>>,
    bool> {
  return not(t == u);
}

static_assert(::ranges::RandomAccessIterator<WindowsMergerWindows::iterator>);
static_assert(
    ::ranges::RandomAccessIterator<WindowsMergerWindows::move_iterator>);
static_assert(
    ::ranges::RandomAccessIterator<WindowsMergerWindows::const_iterator>);
static_assert(
    ::ranges::RandomAccessIterator<WindowsMergerWindows::reverse_iterator>);
static_assert(::ranges::RandomAccessIterator<
              WindowsMergerWindows::move_reverse_iterator>);
static_assert(::ranges::RandomAccessIterator<
              WindowsMergerWindows::const_reverse_iterator>);

} // namespace windows_merger

#pragma once

#include "windows_merger_exceptions.hpp"
#include "windows_merger_window.hpp"
#include "windows_merger_window_accessor.hpp"
#include "windows_merger_window_accessor_coverages_accessor.hpp"
#include "windows_merger_window_base.hpp"
#include "windows_merger_window_base_iterator.hpp"

#include <cassert>
#include <utility>

namespace windows_merger {

template <typename Merger>
WindowsMergerWindowAccessor<Merger>::WindowsMergerWindowAccessor(
    Merger &merger, index_type window_index) noexcept
    : merger(&merger), _index(window_index) {
  assert(window_index < merger.windows_size());
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::index() const noexcept -> index_type {
  return _index;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::begin_index() const noexcept
    -> bases_size_type {
  if (merger->is_allocated())
    return std::as_const(*merger).template get_window_first_pointer<0>(
        _index)[0];
  else
    return bases_size_type(0);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::end_index() const noexcept
    -> bases_size_type {
  if (merger->is_allocated())
    return std::as_const(*merger).template get_window_first_pointer<0>(
        _index)[1];
  else
    return bases_size_type(0);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::size() const noexcept
    -> bases_size_type {
  if (merger->is_allocated()) {
    const auto start_end =
        std::as_const(*merger).template get_window_first_pointer<0>(_index);
    assert(start_end[0] <= start_end[1]);
    return start_end[1] - start_end[0];
  } else
    return bases_size_type(0);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::front() const noexcept
    -> base_accessor {
  assert(size() > 0);
  return base_accessor(*merger, _index, 0);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::back() const noexcept
    -> base_accessor {
  const auto window_size = size();
  assert(window_size > 0);
  return base_accessor(*merger, _index, window_size - 1);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator[](
    bases_size_type index) const noexcept -> base_accessor {
  assert(index < size());
  return base_accessor(*merger, _index, index);
}

template <typename Merger>
template <typename... Args>
auto WindowsMergerWindowAccessor<Merger>::emplace_back(Args &&...args) const
    noexcept(false) -> base_accessor {
  static_assert(not std::is_const_v<Merger>);
  return merger->base_emplace_back(
      _index, WindowsMergerWindowBase(std::forward<Args>(args)...));
}

template <typename Merger>
template <typename Base>
void WindowsMergerWindowAccessor<Merger>::push_back(Base &&base) const
    noexcept(false) {
  static_assert(not std::is_const_v<Merger>);
  static_assert(std::is_same_v<std::decay_t<Base>, WindowsMergerWindowBase>);
  merger->base_emplace_back(_index, std::forward<Base>(base));
}

template <typename Merger>
void WindowsMergerWindowAccessor<Merger>::set_begin_index(
    bases_size_type begin_index) const noexcept(false) {
  merger->set_begin_index(_index, begin_index);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::coverages() const noexcept
    -> coverages_accessor {
  return coverages_accessor(*merger, _index);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::clusters_size() const noexcept
    -> clusters_size_type {
  return merger->n_clusters;
}

template <typename Merger>
bool WindowsMergerWindowAccessor<Merger>::empty() const noexcept {
  return size() == 0;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::begin() const noexcept -> iterator {
  return iterator(*merger, _index);
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::end() const noexcept -> iterator {
  return iterator(*merger, _index, size());
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::rbegin() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*merger, _index, size()));
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::rend() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*merger, _index));
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows> const &rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows const> const &rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows const> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows &&> const &rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindowAccessor<WindowsMergerWindows &&> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindow const &rhs) const noexcept(false) -> self const & {
  assign_from_window(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::operator=(
    WindowsMergerWindow &&rhs) const noexcept(false) -> self const & {
  assign_from_window(std::move(rhs));
  return *this;
}

template <typename Merger>
template <typename Accessor>
void WindowsMergerWindowAccessor<Merger>::assign_from_accessor(
    Accessor &&rhs) const noexcept(false) {
  const auto n_clusters = merger->n_clusters;
  if (n_clusters != rhs.merger->n_clusters)
    throw exceptions::InvalidClustersSize(
        "cannot perform (copy|move)-assignment from a window with a different "
        "number of clusters");

  const auto rhs_start_end =
      rhs.merger->template get_window_first_pointer<0>(rhs._index);
  const auto new_bases_size =
      static_cast<bases_size_type>(rhs_start_end[1] - rhs_start_end[0]);
  if (merger->n_bases_capacity < new_bases_size)
    merger->reshape(new_bases_size,
                    typename decayed_merger_type::reserver_type(0));

  const auto start_end = merger->template get_window_first_pointer<0>(_index);
  const auto old_size =
      static_cast<bases_size_type>(start_end[1] - start_end[0]);
  const auto bases_to_assign = std::min(new_bases_size, old_size);

  {
    const auto weights_pointer =
        merger->template get_window_first_pointer<1>(_index);
    const auto rhs_weights_pointer =
        rhs.merger->template get_window_first_pointer<1>(rhs._index);
    if constexpr (std::is_rvalue_reference_v<
                      typename std::decay_t<Accessor>::merger_type>) {
      std::move(rhs_weights_pointer,
                rhs_weights_pointer + bases_to_assign * n_clusters,
                weights_pointer);
    } else {
      std::copy(rhs_weights_pointer,
                rhs_weights_pointer + bases_to_assign * n_clusters,
                weights_pointer);
    }
  }

  {
    const auto coverages_pointer =
        merger->template get_coverage_first_pointer<1>(_index);
    const auto rhs_coverages_pointer =
        rhs.merger->template get_coverage_first_pointer<1>(rhs._index);
    if constexpr (std::is_rvalue_reference_v<
                      typename std::decay_t<Accessor>::merger_type>) {
      std::move(rhs_coverages_pointer, rhs_coverages_pointer + bases_to_assign,
                coverages_pointer);
    } else {
      std::copy(rhs_coverages_pointer, rhs_coverages_pointer + bases_to_assign,
                coverages_pointer);
    }
  }

  if (old_size > new_bases_size) {
    merger->destroy_range(_index, new_bases_size, old_size);
  } else if (old_size < new_bases_size) {
    if constexpr (std::is_rvalue_reference_v<
                      typename std::decay_t<Accessor>::merger_type>) {
      merger->move_construct_from_range(*rhs.merger, rhs._index, old_size,
                                        new_bases_size, _index, old_size);
    } else {
      merger->copy_construct_from_range(*rhs.merger, rhs._index, old_size,
                                        new_bases_size, _index, old_size);
    }
  }

  if constexpr (std::is_rvalue_reference_v<
                    typename std::decay_t<Accessor>::merger_type>) {
    start_end[0] = std::move(rhs_start_end[0]);
    start_end[1] = std::move(rhs_start_end[1]);
  } else {
    start_end[0] = rhs_start_end[0];
    start_end[1] = rhs_start_end[1];
  }

  {
    const auto coverages_size =
        merger->template get_coverage_first_pointer<0>(_index);
    const auto rhs_coverages_size =
        rhs.merger->template get_coverage_first_pointer<0>(rhs._index);

    if constexpr (std::is_rvalue_reference_v<
                      typename std::decay_t<Accessor>::merger_type>) {
      *coverages_size = std::move(*rhs_coverages_size);
    } else {
      *coverages_size = *rhs_coverages_size;
    }
  }
}

template <typename Merger>
template <typename Window>
void WindowsMergerWindowAccessor<Merger>::assign_from_window(
    Window &&window) const noexcept(false) {
  const auto n_clusters = merger->n_clusters;
  if (n_clusters != window.clusters_size())
    throw exceptions::InvalidClustersSize(
        "cannot perform (copy|move)-assignment from a window with a different "
        "number of clusters");

  const auto new_bases_size = window.size();
  if (merger->n_bases_capacity < new_bases_size)
    merger->reshape(new_bases_size,
                    typename decayed_merger_type::reserver_type(0));

  const auto start_end = merger->template get_window_first_pointer<0>(_index);
  const auto old_size =
      static_cast<bases_size_type>(start_end[1] - start_end[0]);
  const auto bases_to_assign = std::min(new_bases_size, old_size);

  {
    const auto weights_pointer =
        merger->template get_window_first_pointer<1>(_index);
    auto &&weights = window.linear_weights();
    if constexpr (std::is_lvalue_reference_v<Window>) {
      std::ranges::copy(std::ranges::begin(weights),
                        std::ranges::next(std::ranges::begin(weights),
                                          bases_to_assign * n_clusters),
                        weights_pointer);
    } else {
      std::ranges::move(std::ranges::begin(weights),
                        std::ranges::next(std::ranges::begin(weights),
                                          bases_to_assign * n_clusters),
                        weights_pointer);
    }
  }

  {
    const auto coverages_pointer =
        merger->template get_coverage_first_pointer<1>(_index);
    auto &&coverages = window.coverages();
    if constexpr (std::is_lvalue_reference_v<Window>) {
      std::ranges::copy(
          std::ranges::begin(coverages),
          std::ranges::next(std::ranges::begin(coverages), bases_to_assign),
          coverages_pointer);
    } else {
      std::ranges::move(
          std::ranges::begin(coverages),
          std::ranges::next(std::ranges::begin(coverages), bases_to_assign),
          coverages_pointer);
    }
  }

  if (old_size > new_bases_size) {
    merger->destroy_range(_index, new_bases_size, old_size);
  } else if (old_size < new_bases_size) {
    if constexpr (std::is_lvalue_reference_v<Window>) {
      merger->copy_construct_from_range(window, old_size, new_bases_size,
                                        _index, old_size);
    } else {
      merger->move_construct_from_range(std::move(window), old_size,
                                        new_bases_size, _index, old_size);
    }
  }

  start_end[0] = window.begin_index();
  start_end[1] = window.end_index();

  const auto coverages_size =
      merger->template get_coverage_first_pointer<0>(_index);

  *coverages_size = window.size();
}

template <typename Merger>
bool WindowsMergerWindowAccessor<Merger>::is_merger_allocated() const noexcept {
  return merger->is_allocated();
}

template <typename Merger>
auto WindowsMergerWindowAccessor<Merger>::merger_get_bases_capacity()
    const noexcept -> bases_size_type {
  return merger->n_bases_capacity;
}

template <typename Merger>
template <std::size_t Index>
auto WindowsMergerWindowAccessor<Merger>::merger_get_window_first_pointer(
    windows_size_type window_index) const noexcept {
  return merger->template get_window_first_pointer<Index>(window_index);
}

template <typename Merger>
template <std::size_t Index>
auto WindowsMergerWindowAccessor<Merger>::merger_get_coverage_first_pointer(
    windows_size_type window_index) const noexcept {
  return merger->template get_coverage_first_pointer<Index>(window_index);
}

template <typename Merger>
template <typename RhsMerger>
void WindowsMergerWindowAccessor<Merger>::merger_move_construct_from_range(
    RhsMerger &from_merger, windows_size_type from_window_index,
    bases_size_type from_first_base, bases_size_type from_last_base,
    windows_size_type to_window_index, bases_size_type to_first_base) const
    noexcept(std::is_nothrow_move_constructible_v<weight_type> and
             std::is_nothrow_move_constructible_v<coverage_type>) {
  merger->move_construct_from_range(from_merger, from_window_index,
                                    from_first_base, from_last_base,
                                    to_window_index, to_first_base);
}

template <typename Merger>
void WindowsMergerWindowAccessor<Merger>::merger_destroy_range(
    windows_size_type window_index, bases_size_type first,
    bases_size_type last) noexcept {
  return merger->destroy_range(window_index, first, last);
}

template <typename Merger>
void WindowsMergerWindowAccessor<Merger>::resize(bases_size_type value) const
    noexcept(false) {
  if (value > merger->n_bases_capacity)
    merger->reshape(value, typename decayed_merger_type::reserver_type(0));

  const auto weights_start_end =
      merger->template get_window_first_pointer<0>(_index);
  const auto old_size =
      static_cast<bases_size_type>(weights_start_end[1] - weights_start_end[0]);
  const auto n_clusters = merger->n_clusters;

  if (value < old_size) {
    merger->destroy_range(_index, value, old_size);
  } else if (value > old_size) {
    const auto weights = merger->template get_window_first_pointer<1>(_index);
    const auto weights_start = weights + old_size * n_clusters;
    const auto weights_end = weights + value * n_clusters;
    for (weight_type *weight_ptr = weights_start; weight_ptr < weights_end;
         ++weight_ptr) {
      try {
        windows_allocator_traits::construct(*merger, weight_ptr);
      } catch (...) {
        while (weight_ptr > weights_start) {
          --weight_ptr;
          windows_allocator_traits::destroy(*merger, weight_ptr);
        }

        throw;
      }
    }

    const auto coverages =
        merger->template get_coverage_first_pointer<1>(_index);
    for (bases_size_type coverage_index = old_size; coverage_index < value;
         ++coverage_index) {
      try {
        windows_allocator_traits::construct(*merger,
                                            coverages + coverage_index);
      } catch (...) {
        while (coverage_index > old_size) {
          --coverage_index;
          windows_allocator_traits::destroy(*merger,
                                            coverages + coverage_index);
        }

        for (weight_type *weight_ptr = weights_end;
             weight_ptr > weights_start;) {
          --weight_ptr;
          windows_allocator_traits::destroy(*merger, weight_ptr);
        }

        throw;
      }
    }
  }

  weights_start_end[1] = weights_start_end[0] + value;
  const auto coverages_size =
      merger->template get_coverage_first_pointer<0>(_index);
  *coverages_size = value;
}

template <typename Merger>
void WindowsMergerWindowAccessor<Merger>::clear() const noexcept {
  {
    auto start_end = merger_get_window_first_pointer<0>(_index);
    merger->destroy_range(_index, 0, start_end[1] - start_end[0]);
    start_end[1] = start_end[0];
  }

  static_assert(std::is_trivially_destructible_v<coverage_type>);
  auto start_end = merger_get_coverage_first_pointer<0>(_index);
  start_end[0] = coverage_type(0);
}

template <typename Merger>
void WindowsMergerWindowAccessor<Merger>::reorder_clusters(
    std::vector<clusters_size_type> clusters_indices) noexcept(false) {
  using pair_type = std::array<clusters_size_type, 2>;
  thread_local std::vector<pair_type> swap_indices;

  auto const n_clusters = merger->n_clusters;
  assert(n_clusters == clusters_indices.size());
  if (n_clusters == 0)
    return;

  auto const n_bases = size();
  if (n_bases == 0)
    return;

  swap_indices.clear();
  constexpr auto temp_data_index =
      std::numeric_limits<clusters_size_type>::max();

  {
    auto position = clusters_size_type(0);
    for (auto clusters_indices_iter = std::begin(clusters_indices),
              clusters_indices_end = std::end(clusters_indices);
         clusters_indices_iter < clusters_indices_end;
         ++clusters_indices_iter, ++position) {

      auto cluster_index = *clusters_indices_iter;
      if (cluster_index != position) {
        if (auto &next_cluster_index = clusters_indices[cluster_index];
            next_cluster_index == position) {

          swap_indices.push_back(pair_type{position, cluster_index});
          std::swap(clusters_indices[position], next_cluster_index);
        } else {
          swap_indices.push_back(pair_type{position, temp_data_index});
          auto current_pos = position;
          while (cluster_index != position) {
            swap_indices.push_back(pair_type{current_pos, cluster_index});
            current_pos = cluster_index;
            cluster_index =
                std::exchange(clusters_indices[cluster_index], cluster_index);
          }

          swap_indices.push_back(pair_type{current_pos, temp_data_index});
        }
      }
    }
  }

  auto const bases_start = merger->template get_window_first_pointer<1>(_index);
  auto const get_base = [&](bases_size_type base_index) {
    return bases_start + base_index * n_clusters;
  };

  weight_type temp_weight;
  for (bases_size_type base_index = 0; base_index < n_bases; ++base_index) {
    auto const base_ptr = get_base(base_index);
    for (auto &&[cluster_index_a, cluster_index_b] : swap_indices) {
      if (cluster_index_b == temp_data_index) {
        std::swap(base_ptr[cluster_index_a], temp_weight);
      } else {
        std::swap(base_ptr[cluster_index_a], base_ptr[cluster_index_b]);
      }
    }
  }
}

template <typename AccessorL, typename AccessorR>
std::enable_if_t<
    std::is_same_v<std::decay_t<AccessorL>,
                   WindowsMergerWindowAccessor<
                       typename std::decay_t<AccessorL>::merger_type>> and
    std::is_same_v<std::decay_t<AccessorR>,
                   WindowsMergerWindowAccessor<
                       typename std::decay_t<AccessorR>::merger_type>>>
swap(AccessorL &&lhs, AccessorR &&rhs) noexcept(false) {
  if (lhs.merger == rhs.merger and lhs._index == rhs._index)
    return;

  assert(lhs.is_merger_allocated());
  assert(rhs.is_merger_allocated());

  using accessor_l_t = std::decay_t<AccessorL>;
  using accessor_r_t = std::decay_t<AccessorR>;
  using lhs_merger_t = typename accessor_l_t::merger_type;
  using rhs_merger_t = typename accessor_r_t::merger_type;

  const auto n_clusters = lhs.clusters_size();
  if (n_clusters != rhs.clusters_size())
    throw exceptions::InvalidClustersSize(
        "cannot perform swap between windows with a different number of "
        "clusters");

  auto lhs_start_end =
      lhs.template merger_get_window_first_pointer<0>(lhs._index);
  auto rhs_start_end =
      rhs.template merger_get_window_first_pointer<0>(rhs._index);
  const auto lhs_bases_size =
      static_cast<typename lhs_merger_t::bases_size_type>(lhs_start_end[1] -
                                                          lhs_start_end[0]);
  const auto rhs_bases_size =
      static_cast<typename rhs_merger_t::bases_size_type>(rhs_start_end[1] -
                                                          rhs_start_end[0]);
  assert(lhs_bases_size <= lhs.merger_get_bases_capacity());
  assert(rhs_bases_size <= rhs.merger_get_bases_capacity());
  const auto bases_to_swap = std::min(lhs_bases_size, rhs_bases_size);

  if (lhs.merger_get_bases_capacity() < rhs_bases_size) {
    lhs.merger->reshape(rhs_bases_size,
                        typename lhs_merger_t::reserver_type(0));
    lhs_start_end = lhs.template merger_get_window_first_pointer<0>(lhs._index);
  }

  if (rhs.merger_get_bases_capacity() < lhs_bases_size) {
    rhs.merger->reshape(lhs_bases_size,
                        typename rhs_merger_t::reserver_type(0));
    rhs_start_end = rhs.template merger_get_window_first_pointer<0>(rhs._index);
  }

  {
    const auto lhs_weights =
        lhs.template merger_get_window_first_pointer<1>(lhs._index);
    const auto rhs_weights =
        rhs.template merger_get_window_first_pointer<1>(rhs._index);
    std::swap_ranges(lhs_weights, lhs_weights + bases_to_swap * n_clusters,
                     rhs_weights);
  }

  {
    const auto lhs_coverages =
        lhs.template merger_get_coverage_first_pointer<1>(lhs._index);
    const auto rhs_coverages =
        rhs.template merger_get_coverage_first_pointer<1>(rhs._index);
    std::swap_ranges(lhs_coverages, lhs_coverages + bases_to_swap,
                     rhs_coverages);
  }

  if (lhs_bases_size < rhs_bases_size) {
    lhs.merger_move_construct_from_range(*rhs.merger, rhs._index,
                                         lhs_bases_size, rhs_bases_size,
                                         lhs._index, lhs_bases_size);
    rhs.merger_destroy_range(rhs._index, lhs_bases_size, rhs_bases_size);
  } else if (rhs_bases_size < lhs_bases_size) {
    rhs.merger_move_construct_from_range(*lhs.merger, lhs._index,
                                         rhs_bases_size, lhs_bases_size,
                                         rhs._index, rhs_bases_size);
    lhs.merger_destroy_range(lhs._index, rhs_bases_size, lhs_bases_size);
  }

  std::swap(lhs_start_end[0], rhs_start_end[0]);
  std::swap(lhs_start_end[1], rhs_start_end[1]);
  const auto lhs_coverage_size =
      lhs.template merger_get_coverage_first_pointer<0>(lhs._index);
  const auto rhs_coverage_size =
      rhs.template merger_get_coverage_first_pointer<0>(rhs._index);
  std::swap(*lhs_coverage_size, *rhs_coverage_size);
}

static_assert(std::random_access_iterator<
              WindowsMergerWindowAccessor<WindowsMergerWindows>::iterator>);
static_assert(
    std::random_access_iterator<
        WindowsMergerWindowAccessor<const WindowsMergerWindows>::iterator>);
static_assert(std::random_access_iterator<
              WindowsMergerWindowAccessor<WindowsMergerWindows &&>::iterator>);
static_assert(
    std::random_access_iterator<
        WindowsMergerWindowAccessor<WindowsMergerWindows>::reverse_iterator>);
static_assert(std::random_access_iterator<WindowsMergerWindowAccessor<
                  const WindowsMergerWindows>::reverse_iterator>);
static_assert(std::random_access_iterator<WindowsMergerWindowAccessor<
                  WindowsMergerWindows &&>::reverse_iterator>);

} // namespace windows_merger

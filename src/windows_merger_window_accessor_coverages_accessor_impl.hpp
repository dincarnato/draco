#pragma once

#include "windows_merger_window_accessor_coverages_accessor.hpp"

namespace windows_merger {

template <typename Merger>
WindowsMergerWindowAccessorCoveragesAccessor<Merger>::
    WindowsMergerWindowAccessorCoveragesAccessor(
        Merger &merger, index_type window_index) noexcept
    : merger(&merger), window_index(window_index) {}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::index()
    const noexcept -> index_type {
  return window_index;
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::size() const noexcept
    -> bases_size_type {
  if (merger->is_allocated())
    return *merger->template get_coverage_first_pointer<0>(window_index);
  else
    return bases_size_type(0);
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::front()
    const noexcept -> reference {
  assert(size() > 0);
  return merger->template get_coverage_first_pointer<1>(window_index)[0];
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::back() const noexcept
    -> reference {
  const auto bases_size = size();
  assert(bases_size > 0);
  return merger->template get_coverage_first_pointer<1>(
      window_index)[bases_size - 1];
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::operator[](
    bases_size_type index) const noexcept -> reference {
  assert(index < size());
  return merger->template get_coverage_first_pointer<1>(window_index)[index];
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::begin()
    const noexcept -> iterator {
  if (merger->is_allocated())
    return merger->template get_coverage_first_pointer<1>(window_index);
  else
    return nullptr;
}

template <typename Merger>
auto WindowsMergerWindowAccessorCoveragesAccessor<Merger>::end() const noexcept
    -> iterator {
  if (merger->is_allocated())
    return merger->template get_coverage_first_pointer<1>(window_index) +
           size();
  else
    return nullptr;
}

} // namespace windows_merger

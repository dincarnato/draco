#pragma once

#include "windows_merger_window_base_iterator.hpp"

#include "windows_merger_window_base.hpp"
#include "windows_merger_window_base_accessor.hpp"

namespace windows_merger {

template <typename Merger>
WindowsMergerWindowBaseIterator<Merger>::WindowsMergerWindowBaseIterator(
    Merger& merger, windows_size_type window_index,
    signed_bases_size_type base_index) noexcept
    : merger(&merger), window_index(window_index), base_index(base_index) {}

template <typename Merger>
auto WindowsMergerWindowBaseIterator<Merger>::operator*() const noexcept
    -> reference {
  assert(base_index < std::numeric_limits<bases_size_type>::max());
  return base_accessor(*merger, window_index,
                       static_cast<bases_size_type>(base_index));
}

template <typename Merger>
auto WindowsMergerWindowBaseIterator<Merger>::
operator[](difference_type offset) const noexcept -> reference {
  return base_accessor(*merger, window_index,
                       static_cast<bases_size_type>(base_index + offset));
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::operator++() noexcept -> self& {
  ++base_index;
  return *this;
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::operator++(int) noexcept -> self {
  auto new_iter = *this;
  ++base_index;
  return new_iter;
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::operator--() noexcept -> self& {
  --base_index;
  return *this;
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::operator--(int) noexcept -> self {
  auto new_iter = *this;
  --base_index;
  return new_iter;
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::
operator+=(difference_type offset) noexcept -> self& {
  base_index += offset;
  return *this;
}

template <typename Merger>
auto
WindowsMergerWindowBaseIterator<Merger>::
operator-=(difference_type offset) noexcept -> self& {
  base_index -= offset;
  return *this;
}

} // namespace windows_merger

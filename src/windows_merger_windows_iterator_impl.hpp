#pragma once

#include "windows_merger_window_accessor.hpp"
#include "windows_merger_windows_iterator.hpp"

namespace windows_merger {

template <typename Merger>
WindowsMergerWindowsIterator<Merger>::WindowsMergerWindowsIterator(
    Merger &merger, signed_windows_size_type window_index) noexcept
    : merger(&merger), window_index(window_index) {}

template <typename Merger>
WindowsMergerWindowsIterator<Merger>::WindowsMergerWindowsIterator(
    Merger &merger, windows_size_type window_index) noexcept
    : merger(&merger),
      window_index(static_cast<signed_windows_size_type>(window_index)) {}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator*() const noexcept
    -> reference {
  assert(static_cast<std::make_unsigned_t<signed_windows_size_type>>(
             window_index) < std::numeric_limits<windows_size_type>::max());
  return accessor(*merger, static_cast<windows_size_type>(window_index));
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator[](
    difference_type offset) const noexcept -> reference {
  return accessor(*merger,
                  static_cast<windows_size_type>(window_index + offset));
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator++() noexcept -> self & {
  ++window_index;
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator++(int) noexcept -> self {
  auto new_iter = *this;
  ++window_index;
  return new_iter;
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator--() noexcept -> self & {
  --window_index;
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator--(int) noexcept -> self {
  auto new_iter = *this;
  --window_index;
  return new_iter;
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator+=(
    difference_type offset) noexcept -> self & {
  window_index += offset;
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowsIterator<Merger>::operator-=(
    difference_type offset) noexcept -> self & {
  window_index -= offset;
  return *this;
}

} // namespace windows_merger

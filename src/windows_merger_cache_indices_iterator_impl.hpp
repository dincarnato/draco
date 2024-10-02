#pragma once

#include "windows_merger_cache_indices_accessor.hpp"
#include "windows_merger_cache_indices_iterator.hpp"
#include "windows_merger_cache_indices_line.hpp"
#include <iterator>

namespace windows_merger {

template <typename Merger>
WindowsMergerCacheIndicesIterator<Merger>::WindowsMergerCacheIndicesIterator(
    Merger &merger, signed_windows_size_type line_index) noexcept
    : merger(&merger), line_index(line_index) {}

template <typename Merger>
WindowsMergerCacheIndicesIterator<Merger>::WindowsMergerCacheIndicesIterator(
    Merger &merger, windows_size_type line_index) noexcept
    : merger(&merger),
      line_index(static_cast<signed_windows_size_type>(line_index)) {}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator*() const noexcept
    -> reference {
  assert(static_cast<std::make_unsigned_t<signed_windows_size_type>>(
             line_index) < std::numeric_limits<windows_size_type>::max());
  return accessor(*merger, static_cast<windows_size_type>(line_index));
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator[](
    difference_type offset) const noexcept -> reference {
  return accessor(*merger, static_cast<windows_size_type>(line_index + offset));
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator++() noexcept
    -> self & {
  ++line_index;
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator++(int) noexcept
    -> self {
  auto new_iter = *this;
  ++line_index;
  return new_iter;
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator--() noexcept
    -> self & {
  --line_index;
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator--(int) noexcept
    -> self {
  auto new_iter = *this;
  --line_index;
  return new_iter;
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator+=(
    difference_type offset) noexcept -> self & {
  line_index += offset;
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesIterator<Merger>::operator-=(
    difference_type offset) noexcept -> self & {
  line_index -= offset;
  return *this;
}

} // namespace windows_merger

#include "windows_merger_windows.hpp"

namespace windows_merger {

static_assert(std::random_access_iterator<
              WindowsMergerCacheIndicesIterator<WindowsMergerWindows>>);
static_assert(std::random_access_iterator<
              WindowsMergerCacheIndicesIterator<const WindowsMergerWindows>>);

} // namespace windows_merger

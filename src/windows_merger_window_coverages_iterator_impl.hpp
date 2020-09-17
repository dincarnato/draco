#pragma once

#include "windows_merger_window_coverages_iterator.hpp"

namespace windows_merger {

template <typename Window>
WindowsMergerWindowCoveragesIterator<
    Window>::WindowsMergerWindowCoveragesIterator(Window& window,
                                                  signed_bases_size_type
                                                      base_index) noexcept
    : window(&window), base_index(base_index) {}

template <typename Window>
auto WindowsMergerWindowCoveragesIterator<Window>::operator*() const noexcept
    -> reference {
  assert(base_index >= 0);
  assert(base_index < std::numeric_limits<bases_size_type>::max());
  return window->_bases[static_cast<bases_size_type>(base_index)]._coverage;
}

template <typename Window>
auto WindowsMergerWindowCoveragesIterator<Window>::
operator[](difference_type offset) const noexcept -> reference {
  assert(base_index <=
         std::numeric_limits<signed_bases_size_type>::max() - offset);
  const auto new_index = base_index + offset;

  assert(new_index >= 0);
  assert(new_index < std::numeric_limits<bases_size_type>::max());
  return window->_bases[static_cast<bases_size_type>(new_index)]._coverage;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::operator++() noexcept -> self& {
  ++base_index;
  return *this;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::operator++(int) noexcept -> self {
  auto new_iter = *this;
  ++base_index;
  return new_iter;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::operator--() noexcept -> self& {
  --base_index;
  return *this;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::operator--(int) noexcept -> self {
  auto new_iter = *this;
  --base_index;
  return new_iter;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::
operator+=(difference_type offset) noexcept -> self& {
  base_index += offset;
  return *this;
}

template <typename Window>
auto
WindowsMergerWindowCoveragesIterator<Window>::
operator-=(difference_type offset) noexcept -> self& {
  base_index -= offset;
  return *this;
}

} // namespace windows_merger

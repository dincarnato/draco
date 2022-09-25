#pragma once

#include "windows_merger_window_coverages_accessor.hpp"
#include "windows_merger_window_coverages_iterator.hpp"

namespace windows_merger {

template <typename Window>
WindowsMergerWindowCoveragesAccessor<
    Window>::WindowsMergerWindowCoveragesAccessor(Window &window) noexcept
    : _window(&window) {}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::size() const noexcept
    -> bases_size_type {
  return _window->size();
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::front() const noexcept
    -> reference {
  return _window->_bases.front()._coverage;
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::back() const noexcept
    -> reference {
  return _window->_bases.back()._coverage;
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::operator[](
    bases_size_type index) const noexcept -> reference {
  return _window->_bases[index]._coverage;
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::begin() const noexcept
    -> iterator {
  return iterator(*_window);
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::end() const noexcept
    -> iterator {
  return iterator(*_window, _window->size());
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::rbegin() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*_window, _window->size()));
}

template <typename Window>
auto WindowsMergerWindowCoveragesAccessor<Window>::rend() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*_window));
}

} // namespace windows_merger

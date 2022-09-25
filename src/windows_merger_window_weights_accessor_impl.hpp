#pragma once

#include "windows_merger_window_weights_accessor.hpp"
#include "windows_merger_window_weights_iterator.hpp"

#include <range/v3/iterator/concepts.hpp>

namespace windows_merger {

template <typename Window>
WindowsMergerWindowWeightsAccessor<Window>::WindowsMergerWindowWeightsAccessor(
    Window &window) noexcept
    : _window(&window) {}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::size() const noexcept
    -> bases_size_type {
  return _window->size() * _window->clusters_size();
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::front() const noexcept
    -> reference {
  assert(not _window->_bases.empty());
  auto &&front_base = _window->_bases.front();
  assert(not front_base._weights.empty());
  return front_base._weights.front();
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::back() const noexcept
    -> reference {
  assert(not _window->_bases.empty());
  auto &&back_base = _window->_bases.back();
  assert(not back_base._weights.empty());
  return back_base._weights.back();
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::operator[](
    bases_size_type index) const noexcept -> reference {
  const auto clusters_size = _window->clusters_size();

  assert(index < _window->_bases.size() * clusters_size);
  const auto base_index = static_cast<bases_size_type>(index / clusters_size);
  auto &&base = _window->_bases[base_index];
  assert(base.clusters_size() == clusters_size);

  return base._weights[static_cast<bases_size_type>(index % clusters_size)];
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::begin() const noexcept
    -> iterator {
  return iterator(*_window);
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::end() const noexcept
    -> iterator {
  return iterator(*_window, size());
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::rbegin() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*_window, size()));
}

template <typename Window>
auto WindowsMergerWindowWeightsAccessor<Window>::rend() const noexcept
    -> reverse_iterator {
  return reverse_iterator(iterator(*_window));
}

} // namespace windows_merger

#include "windows_merger_windows.hpp"

namespace windows_merger {

static_assert(
    ::ranges::RandomAccessIterator<
        WindowsMergerWindowWeightsAccessor<WindowsMergerWindow>::iterator>);
static_assert(::ranges::RandomAccessIterator<WindowsMergerWindowWeightsAccessor<
                  const WindowsMergerWindow>::iterator>);
static_assert(
    ::ranges::RandomAccessIterator<
        WindowsMergerWindowWeightsAccessor<WindowsMergerWindow &&>::iterator>);

} // namespace windows_merger

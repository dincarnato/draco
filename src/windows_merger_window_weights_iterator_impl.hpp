#pragma once

#include "windows_merger_window_weights_iterator.hpp"

namespace windows_merger {

template <typename Window>
WindowsMergerWindowWeightsIterator<Window>::WindowsMergerWindowWeightsIterator(
    Window &window, signed_bases_size_type weight_index) noexcept
    : window(&window), weight_index(weight_index),
      n_clusters(window.clusters_size()) {}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator*() const noexcept
    -> reference {
  assert(weight_index < std::numeric_limits<bases_size_type>::max());
  return window->_bases[static_cast<bases_size_type>(weight_index / n_clusters)]
      ._weights[static_cast<clusters_size_type>(weight_index % n_clusters)];
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator[](
    difference_type offset) const noexcept -> reference {
  auto new_index = static_cast<bases_size_type>(weight_index + offset);
  return window->_bases[static_cast<bases_size_type>(new_index / n_clusters)]
      ._weights[static_cast<clusters_size_type>(new_index % n_clusters)];
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator++() noexcept
    -> self & {
  ++weight_index;
  return *this;
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator++(int) noexcept
    -> self {
  auto new_iter = *this;
  ++weight_index;
  return new_iter;
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator--() noexcept
    -> self & {
  --weight_index;
  return *this;
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator--(int) noexcept
    -> self {
  auto new_iter = *this;
  --weight_index;
  return new_iter;
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator+=(
    difference_type offset) noexcept -> self & {
  weight_index += offset;
  return *this;
}

template <typename Window>
auto WindowsMergerWindowWeightsIterator<Window>::operator-=(
    difference_type offset) noexcept -> self & {
  weight_index -= offset;
  return *this;
}

} // namespace windows_merger

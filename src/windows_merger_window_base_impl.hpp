#pragma once

#include "windows_merger_window_base.hpp"

#include "windows_merger_window_base_accessor.hpp"
#include "windows_merger_window_base_weights_accessor.hpp"
#include "windows_merger_windows.hpp"

#include <cassert>

namespace windows_merger {

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    clusters_size_type n_clusters) noexcept(false)
    : _weights(n_clusters), _coverage(0) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> const
        &rhs) noexcept(false)
    : _weights(std::begin(rhs.weights()), std::end(rhs.weights())),
      _coverage(rhs.coverage()) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> &&rhs) noexcept(false)
    : WindowsMergerWindowBase(rhs) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> const
        &rhs) noexcept(false)
    : _weights(std::begin(rhs.weights()), std::end(rhs.weights())),
      _coverage(rhs.coverage()) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const>
        &&rhs) noexcept(false)
    : WindowsMergerWindowBase(rhs) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> const
        &rhs) noexcept(false)
    : _weights(std::move_iterator(std::begin(rhs.weights())),
               std::move_iterator(std::end(rhs.weights()))),
      _coverage(std::move(rhs.coverage())) {}

inline WindowsMergerWindowBase::WindowsMergerWindowBase(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&>
        &&rhs) noexcept(false)
    : WindowsMergerWindowBase(rhs) {}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> const
        &rhs) noexcept(false) -> self & {
  assign_from_accessor(rhs);
  return *this;
}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> &&rhs) noexcept(false)
    -> self & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> const
        &rhs) noexcept(false) -> self & {
  assign_from_accessor(rhs);
  return *this;
}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const>
        &&rhs) noexcept(false) -> self & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> const
        &rhs) noexcept(false) -> self & {
  assign_from_accessor(rhs);
  return *this;
}

inline auto WindowsMergerWindowBase::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&>
        &&rhs) noexcept(false) -> self & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Accessor>
void WindowsMergerWindowBase::assign_from_accessor(
    Accessor &&accessor) noexcept(false) {
  _weights.resize(accessor.clusters_size());
  auto &&accessor_weights = accessor.weights();

  using merger_type = typename std::decay_t<Accessor>::merger_type;

  if constexpr (std::is_rvalue_reference_v<merger_type>) {
    std::move(std::begin(accessor_weights), std::end(accessor_weights),
              std::begin(_weights));
    _coverage = std::move(accessor.coverage());
  } else {
    std::copy(std::begin(accessor_weights), std::end(accessor_weights),
              std::begin(_weights));
    _coverage = accessor.coverage();
  }
}

inline auto WindowsMergerWindowBase::weight(clusters_size_type index) noexcept
    -> weight_type & {
  assert(index < _weights.size());
  return _weights[index];
}

inline auto WindowsMergerWindowBase::weight(
    clusters_size_type index) const noexcept -> const weight_type & {
  assert(index < _weights.size());
  return _weights[index];
}

inline auto WindowsMergerWindowBase::coverage() noexcept -> coverage_type & {
  return _coverage;
}

inline auto
WindowsMergerWindowBase::coverage() const noexcept -> const coverage_type & {
  return _coverage;
}

inline auto
WindowsMergerWindowBase::clusters_size() const noexcept -> clusters_size_type {
  return static_cast<clusters_size_type>(_weights.size());
}

inline auto
WindowsMergerWindowBase::weights() noexcept -> weights_accessor_type {
  return weights_accessor_type(_weights.data(),
                               _weights.data() + std::size(_weights));
}

inline auto WindowsMergerWindowBase::weights() const noexcept
    -> const_weights_accessor_type {
  return const_weights_accessor_type(_weights.data(),
                                     _weights.data() + std::size(_weights));
}

} // namespace windows_merger

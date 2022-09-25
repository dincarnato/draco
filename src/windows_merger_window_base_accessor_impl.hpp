#pragma once

#include "windows_merger_window_base.hpp"
#include "windows_merger_window_base_accessor.hpp"

namespace windows_merger {

template <typename Merger>
WindowsMergerWindowBaseAccessor<Merger>::WindowsMergerWindowBaseAccessor(
    Merger &merger, windows_size_type window_index,
    bases_size_type base_index) noexcept
    : merger(&merger), _window_index(window_index), _base_index(base_index) {}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> const &rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> const &rhs)
    const noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> const &rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBase const &rhs) const noexcept(false) -> self const & {
  assign_from_base(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::operator=(
    WindowsMergerWindowBase &&rhs) const noexcept(false) -> self const & {
  assign_from_base(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::weight(
    clusters_size_type index) const noexcept -> weight_type & {
  assert(index < clusters_size());
  const auto weights_base =
      merger->template get_window_first_pointer<1>(_window_index);
  return weights_base[_base_index * merger->n_clusters + index];
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::coverage() const noexcept
    -> coverage_type & {
  const auto coverages_base =
      merger->template get_coverage_first_pointer<1>(_window_index);
  return coverages_base[_base_index];
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::clusters_size() const noexcept
    -> clusters_size_type {
  return merger->n_clusters;
}

template <typename Merger>
auto WindowsMergerWindowBaseAccessor<Merger>::weights() const noexcept
    -> weights_accessor_type {
  const auto weights_base =
      merger->template get_window_first_pointer<1>(_window_index);
  return weights_accessor_type(weights_base + _base_index * merger->n_clusters,
                               weights_base +
                                   (_base_index + 1) * merger->n_clusters);
}

template <typename Merger>
void WindowsMergerWindowBaseAccessor<Merger>::perform_checks() const noexcept {
  assert(merger->is_allocated());
  assert(_window_index < merger->n_windows);

#ifndef NDEBUG
  auto window_start_end =
      merger->template get_window_first_pointer<0>(_window_index);
  assert(_base_index < window_start_end[1] - window_start_end[0]);
#endif
}

template <typename Merger>
template <typename Accessor>
void WindowsMergerWindowBaseAccessor<Merger>::assign_from_accessor(
    Accessor &&rhs) const noexcept {
  perform_checks();
  rhs.perform_checks();
  assert(clusters_size() == rhs.clusters_size());

  using rhs_merger = typename std::decay_t<Accessor>::merger_type;

  const auto lhs_weights = weights();
  const auto rhs_weights = rhs.weights();

  if constexpr (std::is_rvalue_reference_v<rhs_merger>) {
    std::move(std::begin(rhs_weights), std::end(rhs_weights),
              std::begin(lhs_weights));
    coverage() = std::move(rhs.coverage());
  } else {
    std::copy(std::begin(rhs_weights), std::end(rhs_weights),
              std::begin(lhs_weights));
    coverage() = rhs.coverage();
  }
}

template <typename Merger>
template <typename Base>
void WindowsMergerWindowBaseAccessor<Merger>::assign_from_base(
    Base &&rhs) const noexcept {
  perform_checks();
  assert(clusters_size() == rhs.clusters_size());

  const auto lhs_weights = weights();
  const auto rhs_weights = rhs.weights();

  if constexpr (std::is_rvalue_reference_v<Base>) {
    std::move(std::begin(rhs_weights), std::end(rhs_weights),
              std::begin(lhs_weights));
    coverage() = std::move(rhs.coverage());
  } else {
    std::copy(std::begin(rhs_weights), std::end(rhs_weights),
              std::begin(lhs_weights));
    coverage() = rhs.coverage();
  }
}

} // namespace windows_merger

#pragma once

#include "windows_merger_traits.hpp"

#include "tiny_fraction.hpp"

#include <vector>

namespace windows_merger {

struct WindowsMergerWindows;

template <typename>
struct WindowsMergerWindowBaseAccessor;

template <typename>
struct WindowsMergerWindowBaseWeightsAccessor;

struct WindowsMergerWindowBase {
  template <typename>
  friend struct WindowsMergerWindowWeightsAccessor;
  template <typename>
  friend struct WindowsMergerWindowWeightsIterator;
  template <typename>
  friend struct WindowsMergerWindowCoveragesAccessor;
  template <typename>
  friend struct WindowsMergerWindowCoveragesIterator;

  using traits_type = WindowsMergerTraits;
  using weight_type = typename traits_type::weight_type;
  using coverage_type = typename traits_type::coverage_type;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using weights_accessor_type =
      WindowsMergerWindowBaseWeightsAccessor<weight_type>;
  using const_weights_accessor_type =
      WindowsMergerWindowBaseWeightsAccessor<const weight_type>;
  using self = WindowsMergerWindowBase;

  WindowsMergerWindowBase() = default;
  explicit WindowsMergerWindowBase(clusters_size_type n_clusters) noexcept(
      false);

  WindowsMergerWindowBase(const self&) = default;
  WindowsMergerWindowBase(self&&) = default;
  WindowsMergerWindowBase(WindowsMergerWindowBaseAccessor<
                          WindowsMergerWindows> const&) noexcept(false);
  WindowsMergerWindowBase(
      WindowsMergerWindowBaseAccessor<WindowsMergerWindows>&&) noexcept(false);
  WindowsMergerWindowBase(WindowsMergerWindowBaseAccessor<
                          WindowsMergerWindows const> const&) noexcept(false);
  WindowsMergerWindowBase(WindowsMergerWindowBaseAccessor<
                          WindowsMergerWindows const>&&) noexcept(false);
  WindowsMergerWindowBase(WindowsMergerWindowBaseAccessor<
                          WindowsMergerWindows&&> const&) noexcept(false);
  WindowsMergerWindowBase(WindowsMergerWindowBaseAccessor<
                          WindowsMergerWindows&&>&&) noexcept(false);

  self& operator=(self const&) = default;
  self& operator=(self&&) = default;

  self& operator=(WindowsMergerWindowBaseAccessor<
                  WindowsMergerWindows> const&) noexcept(false);
  self& operator=(
      WindowsMergerWindowBaseAccessor<WindowsMergerWindows>&&) noexcept(false);
  self& operator=(WindowsMergerWindowBaseAccessor<
                  WindowsMergerWindows const> const&) noexcept(false);
  self& operator=(WindowsMergerWindowBaseAccessor<
                  WindowsMergerWindows const>&&) noexcept(false);
  self& operator=(WindowsMergerWindowBaseAccessor<
                  WindowsMergerWindows&&> const&) noexcept(false);
  self&
  operator=(WindowsMergerWindowBaseAccessor<WindowsMergerWindows&&>&&) noexcept(
      false);

  weight_type& weight(clusters_size_type index) noexcept;
  const weight_type& weight(clusters_size_type index) const noexcept;

  clusters_size_type clusters_size() const noexcept;

  coverage_type& coverage() noexcept;
  const coverage_type& coverage() const noexcept;

  weights_accessor_type weights() noexcept;
  const_weights_accessor_type weights() const noexcept;

private:
  template <typename Accessor>
  void assign_from_accessor(Accessor&& accessor) noexcept(false);

  std::vector<weight_type> _weights;
  coverage_type _coverage;
};

} // namespace windows_merger

#include "windows_merger_window_base_impl.hpp"

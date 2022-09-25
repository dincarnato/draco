#pragma once

#include "windows_merger_traits.hpp"

#include "nostd/type_traits.hpp"

namespace windows_merger {

template <typename> struct WindowsMergerWindowAccessor;

template <typename> struct WindowsMergerWindowBaseWeightsAccessor;

struct WindowsMergerWindows;

template <typename Merger> struct WindowsMergerWindowBaseAccessor {
  template <typename> friend struct WindowsMergerWindowBaseAccessor;

  using merger_type = std::decay_t<Merger>;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using window_accessor = WindowsMergerWindowAccessor<Merger>;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using windows_size_type = typename traits_type::windows_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using coverage_type =
      nostd::copy_const_t<Merger, typename traits_type::coverage_type>;
  using weight_type =
      nostd::copy_const_t<Merger, typename traits_type::weight_type>;
  using weights_accessor_type =
      WindowsMergerWindowBaseWeightsAccessor<weight_type>;
  using self = WindowsMergerWindowBaseAccessor;

  WindowsMergerWindowBaseAccessor(Merger &merger,
                                  windows_size_type window_index,
                                  bases_size_type base_index) noexcept;
  WindowsMergerWindowBaseAccessor(const WindowsMergerWindowBaseAccessor &) =
      default;
  WindowsMergerWindowBaseAccessor(WindowsMergerWindowBaseAccessor &&) = default;

  self const &
  operator=(WindowsMergerWindowBaseAccessor<WindowsMergerWindows> const &) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowBaseAccessor<WindowsMergerWindows> &&) const
      noexcept(false);
  self const &operator=(
      WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> const &) const
      noexcept(false);
  self const &operator=(
      WindowsMergerWindowBaseAccessor<WindowsMergerWindows const> &&) const
      noexcept(false);
  self const &operator=(
      WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> const &) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowBaseAccessor<WindowsMergerWindows &&> &&) const
      noexcept(false);

  self const &operator=(WindowsMergerWindowBase const &) const noexcept(false);
  self const &operator=(WindowsMergerWindowBase &&) const noexcept(false);

  weight_type &weight(clusters_size_type index) const noexcept;
  coverage_type &coverage() const noexcept;
  clusters_size_type clusters_size() const noexcept;
  weights_accessor_type weights() const noexcept;

private:
  void perform_checks() const noexcept;

  template <typename Accessor>
  void assign_from_accessor(Accessor &&rhs) const noexcept;
  template <typename Base> void assign_from_base(Base &&rhs) const noexcept;

  merger_pointer_type merger;
  windows_size_type _window_index;
  bases_size_type _base_index;
};

} // namespace windows_merger

#include "windows_merger_window_base_accessor_impl.hpp"

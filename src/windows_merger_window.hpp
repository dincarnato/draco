#pragma once

#include "windows_merger_traits.hpp"

#include <range/v3/core.hpp>

#include <vector>

class WeightedClusters;

namespace windows_merger {

struct WindowsMergerWindowBase;
struct WindowsMergerWindow;
struct WindowsMergerWindows;
template <typename> struct WindowsMergerWindowAccessor;
template <typename> struct WindowsMergerWindowCoveragesAccessor;
template <typename> struct WindowsMergerWindowWeightsAccessor;

struct WindowsMergerWindow {
  template <typename> friend struct WindowsMergerWindowCoveragesAccessor;
  template <typename> friend struct WindowsMergerWindowCoveragesIterator;
  template <typename> friend struct WindowsMergerWindowWeightsAccessor;
  template <typename> friend struct WindowsMergerWindowWeightsIterator;

  using traits_type = WindowsMergerTraits;
  using index_type = typename traits_type::windows_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using weight_type = typename traits_type::weight_type;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using coverage_type = typename traits_type::coverage_type;
  using base_type = WindowsMergerWindowBase;
  using bases_type = std::vector<base_type>;
  using coverages_accessor =
      WindowsMergerWindowCoveragesAccessor<WindowsMergerWindow>;
  using const_coverages_accessor =
      WindowsMergerWindowCoveragesAccessor<const WindowsMergerWindow>;
  using linear_weights_accessor =
      WindowsMergerWindowWeightsAccessor<WindowsMergerWindow>;
  using const_linear_weights_accessor =
      WindowsMergerWindowWeightsAccessor<const WindowsMergerWindow>;
  using iterator = typename bases_type::iterator;
  using const_iterator = typename bases_type::const_iterator;
  using move_iterator = ranges::move_iterator<typename bases_type::iterator>;
  using reverse_iterator = typename bases_type::reverse_iterator;
  using const_reverse_iterator = typename bases_type::const_reverse_iterator;
  using move_reverse_iterator =
      ranges::move_iterator<typename bases_type::reverse_iterator>;
  using self = WindowsMergerWindow;

  WindowsMergerWindow() = default;
  explicit WindowsMergerWindow(clusters_size_type n_clusters) noexcept(
      std::is_nothrow_default_constructible_v<bases_type>);

  WindowsMergerWindow(WindowsMergerWindow const &) = default;
  WindowsMergerWindow(WindowsMergerWindow &&) = default;

  WindowsMergerWindow(WindowsMergerWindowAccessor<WindowsMergerWindows> const
                          &) noexcept(false);
  WindowsMergerWindow(
      WindowsMergerWindowAccessor<WindowsMergerWindows> &&) noexcept(false);
  WindowsMergerWindow(
      WindowsMergerWindowAccessor<WindowsMergerWindows const> const
          &) noexcept(false);
  WindowsMergerWindow(WindowsMergerWindowAccessor<WindowsMergerWindows const>
                          &&) noexcept(false);
  WindowsMergerWindow(WindowsMergerWindowAccessor<WindowsMergerWindows &&> const
                          &) noexcept(false);
  WindowsMergerWindow(
      WindowsMergerWindowAccessor<WindowsMergerWindows &&> &&) noexcept(false);

  self &operator=(WindowsMergerWindow const &) = default;
  self &operator=(WindowsMergerWindow &&) = default;

  self &operator=(WindowsMergerWindowAccessor<WindowsMergerWindows> const
                      &) noexcept(false);
  self &operator=(
      WindowsMergerWindowAccessor<WindowsMergerWindows> &&) noexcept(false);
  self &operator=(WindowsMergerWindowAccessor<WindowsMergerWindows const> const
                      &) noexcept(false);
  self &operator=(WindowsMergerWindowAccessor<WindowsMergerWindows const>
                      &&) noexcept(false);
  self &operator=(WindowsMergerWindowAccessor<WindowsMergerWindows &&> const
                      &) noexcept(false);
  self &operator=(
      WindowsMergerWindowAccessor<WindowsMergerWindows &&> &&) noexcept(false);

  bases_size_type begin_index() const noexcept;
  bases_size_type end_index() const noexcept;
  bases_size_type size() const noexcept;
  bool empty() const noexcept;
  clusters_size_type clusters_size() const noexcept;
  coverages_accessor coverages() noexcept;
  const_coverages_accessor coverages() const noexcept;
  linear_weights_accessor linear_weights() noexcept;
  const_linear_weights_accessor linear_weights() const noexcept;
  void set_begin_index(bases_size_type begin_index) noexcept;

  base_type const &front() const noexcept;
  base_type const &back() const noexcept;

  base_type &front() noexcept;
  base_type &back() noexcept;

  iterator begin() & noexcept;
  iterator end() & noexcept;
  const_iterator begin() const & noexcept;
  const_iterator end() const & noexcept;
  move_iterator begin() && noexcept;
  move_iterator end() && noexcept;
  reverse_iterator rbegin() & noexcept;
  reverse_iterator rend() & noexcept;
  const_reverse_iterator rbegin() const & noexcept;
  const_reverse_iterator rend() const & noexcept;
  move_reverse_iterator rbegin() && noexcept;
  move_reverse_iterator rend() && noexcept;

  base_type const &operator[](bases_size_type index) const noexcept;
  base_type &operator[](bases_size_type index) noexcept;

  template <typename... Args>
  base_type &emplace_back(Args &&...args) noexcept(false);

  template <typename Base> void push_back(Base &&base) noexcept(false);

private:
  template <typename Accessor> void assign_from_accessor(Accessor &&accessor);

  bases_type _bases;
  bases_size_type _begin_index = 0;
  clusters_size_type n_clusters = 0;
};

} // namespace windows_merger

#include "windows_merger_window_impl.hpp"

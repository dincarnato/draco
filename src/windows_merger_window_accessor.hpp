#pragma once

#include "windows_merger_traits.hpp"
#include "windows_merger_windows.hpp"

#include <type_traits>
#include <utility>

namespace windows_merger {

template <typename> struct WindowsMergerWindowAccessor;
template <typename> struct WindowsMergerWindowBaseAccessor;
template <typename> struct WindowsMergerWindowAccessorCoveragesAccessor;
template <typename> struct WindowsMergerWindowBaseIterator;

template <typename Merger> struct WindowsMergerWindowAccessor {
  template <typename>
  friend struct ::windows_merger::WindowsMergerWindowAccessor;

  using merger_type = Merger;
  using decayed_merger_type = std::decay_t<Merger>;
  using merger_pointer_type = std::remove_reference_t<Merger> *;
  using traits_type = WindowsMergerTraits;
  using index_type = typename traits_type::windows_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using windows_size_type = typename traits_type::windows_size_type;
  using coverage_type = typename traits_type::coverage_type;
  using windows_allocator_traits =
      typename std::decay_t<Merger>::windows_allocator_traits;
  using weight_type = typename traits_type::weight_type;
  using base_accessor = WindowsMergerWindowBaseAccessor<Merger>;
  using coverages_accessor =
      WindowsMergerWindowAccessorCoveragesAccessor<Merger>;
  using iterator = WindowsMergerWindowBaseIterator<Merger>;
  using reverse_iterator =
      ranges::reverse_iterator<WindowsMergerWindowBaseIterator<Merger>>;
  using self = WindowsMergerWindowAccessor;

  explicit WindowsMergerWindowAccessor(Merger &merger,
                                       index_type window_index) noexcept;
  WindowsMergerWindowAccessor(const WindowsMergerWindowAccessor &) = default;
  WindowsMergerWindowAccessor(WindowsMergerWindowAccessor &&) = default;

  self const &
  operator=(WindowsMergerWindowAccessor<WindowsMergerWindows> const &) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowAccessor<WindowsMergerWindows> &&) const
      noexcept(false);
  self const &operator=(
      WindowsMergerWindowAccessor<WindowsMergerWindows const> const &) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowAccessor<WindowsMergerWindows const> &&) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowAccessor<WindowsMergerWindows &&> const &) const
      noexcept(false);
  self const &
  operator=(WindowsMergerWindowAccessor<WindowsMergerWindows &&> &&) const
      noexcept(false);

  self const &operator=(WindowsMergerWindow const &) const noexcept(false);
  self const &operator=(WindowsMergerWindow &&) const noexcept(false);

  template <typename AccessorL, typename AccessorR>
  friend std::enable_if_t<
      std::is_same_v<std::decay_t<AccessorL>,
                     WindowsMergerWindowAccessor<
                         typename std::decay_t<AccessorL>::merger_type>> and
      std::is_same_v<std::decay_t<AccessorR>,
                     WindowsMergerWindowAccessor<
                         typename std::decay_t<AccessorR>::merger_type>>>
  swap(AccessorL &&lhs, AccessorR &&rhs) noexcept(false);

  index_type index() const noexcept;
  bases_size_type begin_index() const noexcept;
  bases_size_type end_index() const noexcept;
  bases_size_type size() const noexcept;
  coverages_accessor coverages() const noexcept;
  clusters_size_type clusters_size() const noexcept;
  bool empty() const noexcept;
  void set_begin_index(bases_size_type begin_index) const noexcept(false);

  iterator begin() const noexcept;
  iterator end() const noexcept;

  reverse_iterator rbegin() const noexcept;
  reverse_iterator rend() const noexcept;

  base_accessor front() const noexcept;
  base_accessor back() const noexcept;

  base_accessor operator[](bases_size_type index) const noexcept;

  template <typename... Args>
  base_accessor emplace_back(Args &&...args) const noexcept(false);

  template <typename Base> void push_back(Base &&base) const noexcept(false);

  void resize(bases_size_type value) const noexcept(false);
  void clear() const noexcept;

  void reorder_clusters(
      std::vector<clusters_size_type> clusters_indices) noexcept(false);

private:
  template <typename Accessor>
  void assign_from_accessor(Accessor &&accessor) const noexcept(false);
  template <typename Window>
  void assign_from_window(Window &&window) const noexcept(false);

  inline bool is_merger_allocated() const noexcept;
  inline bases_size_type merger_get_bases_capacity() const noexcept;

  template <std::size_t Index>
  auto merger_get_window_first_pointer(
      windows_size_type window_index) const noexcept;

  template <std::size_t Index>
  auto merger_get_coverage_first_pointer(
      windows_size_type window_index) const noexcept;

  template <typename RhsMerger>
  void merger_move_construct_from_range(RhsMerger &from_merger,
                                        windows_size_type from_window_index,
                                        bases_size_type from_first_base,
                                        bases_size_type from_last_base,
                                        windows_size_type to_window_index,
                                        bases_size_type to_first_base) const
      noexcept(std::is_nothrow_move_constructible_v<weight_type>
                   and std::is_nothrow_move_constructible_v<coverage_type>);

  inline void merger_destroy_range(windows_size_type window_index,
                                   bases_size_type first,
                                   bases_size_type last) noexcept;

  merger_pointer_type merger;
  index_type _index;
};

} // namespace windows_merger

#include "windows_merger_window_accessor_impl.hpp"

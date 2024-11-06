#pragma once

#include "windows_merger_traits.hpp"

#include "het/allocator.hpp"
#include "het/vec2d_builder.hpp"

#include "tiny_fraction.hpp"

namespace windows_merger {

template <typename> struct WindowsMergerWindowsIterator;
template <typename> struct WindowsMergerWindowAccessorCoveragesAccessor;

struct WindowsMergerWindow;
template <typename> struct WindowsMergerWindowAccessor;
template <typename> struct WindowsMergerWindowBaseAccessor;

struct WindowsMergerWindows
    : het::allocator<typename WindowsMergerTraits::bases_size_type,
                     typename WindowsMergerTraits::weight_type>,
      het::allocator<typename WindowsMergerTraits::bases_size_type,
                     typename WindowsMergerTraits::coverage_type> {
  template <typename> friend struct WindowsMergerWindowAccessor;
  template <typename> friend struct WindowsMergerWindowBaseAccessor;
  template <typename>
  friend struct WindowsMergerWindowAccessorCoveragesAccessor;

  using traits_type = WindowsMergerTraits;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using windows_size_type = typename traits_type::windows_size_type;
  using weight_type = typename traits_type::weight_type;
  using coverage_type = typename traits_type::coverage_type;
  using windows_allocator = het::allocator<bases_size_type, weight_type>;
  using windows_allocator_traits = het::allocator_traits<windows_allocator>;
  using coverages_allocator = het::allocator<bases_size_type, coverage_type>;
  using coverages_allocator_traits = het::allocator_traits<coverages_allocator>;
  using window_accessor = WindowsMergerWindowAccessor<WindowsMergerWindows>;
  using const_window_accessor =
      WindowsMergerWindowAccessor<const WindowsMergerWindows>;
  using move_window_accessor =
      WindowsMergerWindowAccessor<WindowsMergerWindows &&>;
  using base_accessor = WindowsMergerWindowBaseAccessor<WindowsMergerWindows>;
  using const_base_accessor =
      WindowsMergerWindowBaseAccessor<const WindowsMergerWindows>;
  using reserver_type = WindowsMergerReserve<windows_size_type>;
  using resizer_type = WindowsMergerResize<windows_size_type>;
  using iterator = WindowsMergerWindowsIterator<WindowsMergerWindows>;
  using move_iterator = WindowsMergerWindowsIterator<WindowsMergerWindows &&>;
  using const_iterator =
      WindowsMergerWindowsIterator<const WindowsMergerWindows>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using move_reverse_iterator = std::reverse_iterator<move_iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  WindowsMergerWindows() = default;
  explicit WindowsMergerWindows(clusters_size_type n_clusters) noexcept;
  WindowsMergerWindows(clusters_size_type n_clusters,
                       bases_size_type n_bases) noexcept;
  WindowsMergerWindows(clusters_size_type n_clusters, bases_size_type n_bases,
                       windows_size_type n_windows);
  ~WindowsMergerWindows() noexcept;

  WindowsMergerWindows(const WindowsMergerWindows &) = delete;
  WindowsMergerWindows &operator=(const WindowsMergerWindows &) = delete;

  WindowsMergerWindows(WindowsMergerWindows &&) noexcept;
  WindowsMergerWindows &operator=(WindowsMergerWindows &&) noexcept;

  bool operator==(const WindowsMergerWindows &other) const noexcept;
  bool operator!=(const WindowsMergerWindows &other) const noexcept;

  clusters_size_type clusters_size() const noexcept;
  bases_size_type bases_capacity() const noexcept;
  windows_size_type windows_size() const noexcept;
  windows_size_type windows_capacity() const noexcept;

  window_accessor operator[](windows_size_type index) & noexcept;
  const_window_accessor operator[](windows_size_type index) const & noexcept;
  move_window_accessor operator[](windows_size_type index) && noexcept;
  window_accessor front() noexcept;
  const_window_accessor front() const noexcept;
  window_accessor back() noexcept;
  const_window_accessor back() const noexcept;

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

  template <template <typename> typename WindowsReshaper>
  void reshape(
      bases_size_type bases_reshaper,
      WindowsReshaper<windows_size_type> &&windows_reshaper) noexcept(false);

  window_accessor emplace_back() noexcept(false);

private:
  using windows_data_type =
      typename windows_allocator_traits::template pointer<0>;
  using coverages_data_type =
      typename windows_allocator_traits::template pointer<0>;

  static constexpr auto windows_builder =
      het::vec2d_builder::build()
          .template set_allocator_type<windows_allocator>()
          .fixed_size()
          .template set_size<2>()
          .default_construction(bases_size_type(0))
          .next()
          .dynamic_size()
          .default_uninitialized()
          .done();

  static constexpr auto coverages_builder =
      het::vec2d_builder::build()
          .template set_allocator_type<coverages_allocator>()
          .fixed_size()
          .template set_size<1>()
          .default_construction(bases_size_type(0))
          .next()
          .dynamic_size()
          .size_from_callable(
              [](const bases_size_type *bases) { return *bases; })
          .default_uninitialized()
          .done();

  template <std::size_t Index>
  auto get_window_first_pointer(windows_size_type window_index) noexcept;

  template <std::size_t Index>
  auto get_window_first_pointer(windows_size_type window_index) const noexcept;

  template <std::size_t Index>
  auto get_coverage_first_pointer(windows_size_type window_index) noexcept;

  template <std::size_t Index>
  auto
  get_coverage_first_pointer(windows_size_type window_index) const noexcept;

  template <typename Base>
  base_accessor base_emplace_back(windows_size_type window_index,
                                  Base &&base) noexcept(false);

  void allocate_impl() noexcept(false);
  inline bool is_allocated() const noexcept;
  void set_begin_index(windows_size_type window_index,
                       bases_size_type new_begin_index) noexcept(false);
  void destroy_range(windows_size_type window_index, bases_size_type first,
                     bases_size_type last) noexcept;
  void copy_construct_from_range(
      const WindowsMergerWindows &from_windows,
      windows_size_type from_window_index, bases_size_type from_first,
      bases_size_type from_last, windows_size_type to_window_index,
      bases_size_type
          to_first) noexcept(noexcept(std::
                                          is_nothrow_copy_constructible_v<
                                              weight_type>) and
                             noexcept(std::is_nothrow_copy_constructible_v<
                                      coverage_type>));
  void copy_construct_from_range(
      const WindowsMergerWindow &from_window, bases_size_type from_first,
      bases_size_type from_last, windows_size_type to_window_index,
      bases_size_type
          to_first) noexcept(noexcept(std::
                                          is_nothrow_copy_constructible_v<
                                              weight_type>) and
                             noexcept(std::is_nothrow_copy_constructible_v<
                                      coverage_type>));
  void move_construct_from_range(
      const WindowsMergerWindows &from_windows,
      windows_size_type from_window_index, bases_size_type from_first,
      bases_size_type from_last, windows_size_type to_window_index,
      bases_size_type
          to_first) noexcept(noexcept(std::
                                          is_nothrow_move_constructible_v<
                                              weight_type>) and
                             noexcept(std::is_nothrow_move_constructible_v<
                                      coverage_type>));
  void move_construct_from_range(
      WindowsMergerWindow &&from_windows, bases_size_type from_first,
      bases_size_type from_last, windows_size_type to_window_index,
      bases_size_type
          to_first) noexcept(noexcept(std::
                                          is_nothrow_move_constructible_v<
                                              weight_type>) and
                             noexcept(std::is_nothrow_move_constructible_v<
                                      coverage_type>));
  void destroy_and_deallocate() noexcept;

  clusters_size_type n_clusters = 0;
  bases_size_type n_bases_capacity = 0;
  windows_size_type n_windows_capacity = 0;
  windows_size_type n_windows = 0;
  windows_data_type windows_data{};
  coverages_data_type coverages_data{};
};

} // namespace windows_merger

#include "windows_merger_windows_impl.hpp"

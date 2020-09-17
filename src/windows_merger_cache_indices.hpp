#pragma once

#include "windows_merger_cache_indices_accessor.hpp"
#include "windows_merger_traits.hpp"

#include "het/allocator.hpp"
#include "het/allocator_traits.hpp"
#include "het/vec2d_builder.hpp"

#include <range/v3/core.hpp>

namespace windows_merger {

template <typename>
struct WindowsMergerCacheIndicesIterator;

struct WindowsMergerCacheIndices : het::allocator<WindowsMergerTraits::windows_size_type, WindowsMergerTraits::windows_size_type> {
  template <typename>
  friend struct WindowsMergerCacheIndicesAccessor;

  template <typename>
  friend struct WindowsMergerCacheIndicesIterator;

  using traits_type = WindowsMergerTraits;
  using windows_size_type = typename traits_type::windows_size_type;
  using allocator = het::allocator<windows_size_type, windows_size_type>;
  using allocator_traits = het::allocator_traits<allocator>;
  using indices_accessor =
      WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices>;
  using const_indices_accessor =
      WindowsMergerCacheIndicesAccessor<const WindowsMergerCacheIndices>;
  using reserver_type = WindowsMergerReserve<windows_size_type>;
  using resizer_type = WindowsMergerResize<windows_size_type>;
  using iterator = WindowsMergerCacheIndicesIterator<WindowsMergerCacheIndices>;
  using const_iterator =
      WindowsMergerCacheIndicesIterator<const WindowsMergerCacheIndices>;
  using reverse_iterator = ranges::reverse_iterator<
      WindowsMergerCacheIndicesIterator<WindowsMergerCacheIndices>>;
  using const_reverse_iterator = ranges::reverse_iterator<
      WindowsMergerCacheIndicesIterator<const WindowsMergerCacheIndices>>;

  WindowsMergerCacheIndices() = default;
  explicit WindowsMergerCacheIndices(windows_size_type line_capacity) noexcept;
  explicit WindowsMergerCacheIndices(windows_size_type line_capacity,
                                     windows_size_type lines) noexcept(false);
  ~WindowsMergerCacheIndices() noexcept;

  WindowsMergerCacheIndices(WindowsMergerCacheIndices const&) = delete;
  WindowsMergerCacheIndices&
  operator=(WindowsMergerCacheIndices const&) = delete;

  WindowsMergerCacheIndices(WindowsMergerCacheIndices&&) noexcept;
  WindowsMergerCacheIndices& operator=(WindowsMergerCacheIndices&&) noexcept;

  windows_size_type size() const noexcept;
  windows_size_type capacity() const noexcept;
  windows_size_type line_capacity() const noexcept;

  indices_accessor operator[](windows_size_type line_index) noexcept;
  const_indices_accessor operator[](windows_size_type line_index) const
      noexcept;

  template <template <typename> typename LinesReshaper>
  void
  reshape(windows_size_type line_capacity,
          LinesReshaper<windows_size_type>&& lines_reshaper) noexcept(false);

  iterator begin() noexcept;
  iterator end() noexcept;
  reverse_iterator rbegin() noexcept;
  reverse_iterator rend() noexcept;

  const_iterator begin() const noexcept;
  const_iterator end() const noexcept;
  const_reverse_iterator rbegin() const noexcept;
  const_reverse_iterator rend() const noexcept;

private:
  using data_type = typename allocator_traits::first_pointer;
  static constexpr auto builder =
      het::vec2d_builder::build()
          .template set_allocator_type<allocator>()
          .fixed_size()
          .template set_size<1>()
          .default_construction(windows_size_type(0))
          .next()
          .dynamic_size()
          .size_from_callable(
              [](windows_size_type const* size) { return *size; })
          .default_uninitialized()
          .done();

  void allocate_impl() noexcept(false);

  template <std::size_t Index>
  auto get_indices_first_pointer(windows_size_type line_index) noexcept;

  template <std::size_t Index>
  auto get_indices_first_pointer(windows_size_type line_index) const noexcept;

  windows_size_type& index_emplace_back(windows_size_type line_index,
                                        windows_size_type index) noexcept;
  void deallocate_and_destroy() noexcept;

  windows_size_type _lines_size = 0;
  windows_size_type _lines_capacity = 0;
  windows_size_type _line_capacity = 0;
  data_type data;
};

} // namespace windows_merger

#include "windows_merger_cache_indices_impl.hpp"

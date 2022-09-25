#pragma once

#include "windows_merger_cache_indices.hpp"
#include "windows_merger_cache_indices_iterator.hpp"

#include <algorithm>
#include <cassert>

namespace windows_merger {

inline WindowsMergerCacheIndices::WindowsMergerCacheIndices(
    windows_size_type line_capacity) noexcept
    : _line_capacity(line_capacity) {}

inline WindowsMergerCacheIndices::WindowsMergerCacheIndices(
    windows_size_type line_capacity, windows_size_type lines) noexcept(false)
    : _lines_size(lines), _lines_capacity(lines),
      _line_capacity(line_capacity) {
  if (_lines_capacity > 0 and _line_capacity > 0)
    allocate_impl();
}

inline WindowsMergerCacheIndices::WindowsMergerCacheIndices(
    WindowsMergerCacheIndices &&other) noexcept
    : _lines_size(std::exchange(other._lines_size, 0)),
      _lines_capacity(std::exchange(other._lines_capacity, 0)),
      _line_capacity(std::exchange(other._line_capacity, 0)),
      data(std::exchange(other.data, nullptr)) {}

inline WindowsMergerCacheIndices &WindowsMergerCacheIndices::operator=(
    WindowsMergerCacheIndices &&other) noexcept {
  deallocate_and_destroy();

  _lines_size = std::exchange(other._lines_size, 0);
  _lines_capacity = std::exchange(other._lines_capacity, 0);
  _line_capacity = std::exchange(other._line_capacity, 0);
  data = std::exchange(other.data, nullptr);

  return *this;
}

inline WindowsMergerCacheIndices::~WindowsMergerCacheIndices() noexcept {
  deallocate_and_destroy();
}

inline void WindowsMergerCacheIndices::deallocate_and_destroy() noexcept {
  builder.destroy(data, *this)
      .with_lines(_lines_size)
      .start()
      .next()
      .with_max_size(_line_capacity)
      .done();

  if (_lines_capacity > 0 and _line_capacity > 0)
    allocator_traits::deallocate(*this, data, _lines_capacity, 1,
                                 _line_capacity);
}

inline void WindowsMergerCacheIndices::allocate_impl() noexcept(false) {
  assert(_lines_capacity > 0);
  assert(_line_capacity > 0);

  data = allocator_traits::allocate(*this, _lines_capacity, 1, _line_capacity);
  try {
    builder.build_at(data, *this)
        .with_lines(_lines_size)
        .start()
        .construct_default()
        .next()
        .with_max_size(_line_capacity)
        .construct_default()
        .done();
  } catch (...) {
    allocator_traits::deallocate(*this, data, _lines_capacity, 1,
                                 _line_capacity);
    throw;
  }
}

template <template <typename> typename LinesReshaper>
void WindowsMergerCacheIndices::reshape(
    windows_size_type line_capacity,
    LinesReshaper<windows_size_type> &&lines_reshaper) noexcept(false) {
  using lines_reshaper_type = std::decay_t<LinesReshaper<windows_size_type>>;
  static_assert(traits_type::template is_reshaper_arg_v<lines_reshaper_type>,
                "lines_reshaper must be a WindowsMergerResize<T> or a "
                "WindowsMergerReserve<T> type");
  static_assert(std::is_convertible_v<typename lines_reshaper_type::type,
                                      windows_size_type>,
                "The type stored in the lines_reshaper must be convertible "
                "to a windows_size_type");
  auto new_lines_size = [&] {
    if constexpr (traits_type::template is_reserve_reshaper_v<
                      lines_reshaper_type>) {
      return _lines_size;
    } else {
      return static_cast<windows_size_type>(lines_reshaper.size);
    }
  }();

  auto new_lines_capacity = [&] {
    if constexpr (traits_type::template is_reserve_reshaper_v<
                      lines_reshaper_type>) {
      return std::max({static_cast<windows_size_type>(lines_reshaper.size),
                       _lines_capacity, new_lines_size});
    } else {
      return std::max(_lines_capacity, new_lines_size);
    }
  }();

  if (_lines_capacity == 0 or _line_capacity == 0) {
    if (new_lines_capacity != 0)
      _lines_capacity = new_lines_capacity;

    if (line_capacity != 0)
      _line_capacity = line_capacity;

    if (_lines_capacity != 0 and _line_capacity != 0)
      allocate_impl();
    else
      return;
  }

  auto new_line_capacity =
      std::max(static_cast<windows_size_type>(line_capacity), _line_capacity);

  if (new_line_capacity == _line_capacity and
      new_lines_capacity == _lines_capacity and new_lines_size == _lines_size)
    return;

  if (new_line_capacity != _line_capacity or
      new_lines_capacity > _lines_capacity) {
    auto new_data = allocator_traits::allocate(*this, new_lines_capacity, 1,
                                               new_line_capacity);
    try {
      builder.move_to(new_data, *this)
          .with_lines(new_lines_size)
          .start()
          .next()
          .with_max_size(new_line_capacity)
          .skip_remaining()
          .done()
          .with_lines(_lines_size)
          .from(data)
          .start()
          .next()
          .with_max_size(_line_capacity)
          .done();
    } catch (...) {
      allocator_traits::deallocate(*this, data, new_lines_capacity, 1,
                                   new_line_capacity);
      throw;
    }

    builder.destroy(data, *this)
        .with_lines(_lines_size)
        .start()
        .next()
        .with_max_size(_line_capacity)
        .done();
    allocator_traits::deallocate(*this, data, _lines_capacity, 1,
                                 _line_capacity);

    data = new_data;
  } else if constexpr (WindowsMergerTraits::template is_resize_reshaper_v<
                           lines_reshaper_type>) {
    if (new_lines_size > _lines_size) {
      builder.build_at(data, *this)
          .with_lines(new_lines_size)
          .skip_lines(_lines_size)
          .start()
          .construct_default()
          .next()
          .with_max_size(new_line_capacity)
          .construct_default()
          .done();
    } else if (new_lines_size < _lines_size) {
      builder.destroy(data, *this)
          .with_lines(_lines_size)
          .skip_lines(new_lines_size)
          .start()
          .next()
          .with_max_size(new_line_capacity)
          .done();
    }
  }

  _line_capacity = new_line_capacity;
  _lines_capacity = new_lines_capacity;
  _lines_size = new_lines_size;
}

template <std::size_t Index>
inline auto WindowsMergerCacheIndices::get_indices_first_pointer(
    windows_size_type line_index) const noexcept {
  static_assert(Index < 2);
  return allocator_traits::first_pointer_of<Index>(*this, data, line_index, 1,
                                                   _line_capacity);
}

template <std::size_t Index>
inline auto WindowsMergerCacheIndices::get_indices_first_pointer(
    windows_size_type line_index) noexcept {
  return const_cast<typename allocator_traits::template pointer<Index>>(
      const_cast<const WindowsMergerCacheIndices &>(*this)
          .get_indices_first_pointer<Index>(line_index));
}

inline auto
WindowsMergerCacheIndices::index_emplace_back(windows_size_type line_index,
                                              windows_size_type index) noexcept
    -> windows_size_type & {
  assert(_lines_capacity > 0);

  assert(line_index < _lines_size);

  if (_line_capacity == 0) {
    _line_capacity = 1;
    allocate_impl();
  }

  auto line_size = get_indices_first_pointer<0>(line_index);
  if (*line_size == _line_capacity) {
    reshape(_line_capacity * 2, reserver_type(_lines_capacity));
    line_size = get_indices_first_pointer<0>(line_index);
  }

  const auto new_index = get_indices_first_pointer<1>(line_index) + *line_size;
  allocator_traits::construct(*this, new_index, index);

  ++*line_size;
  return *new_index;
}

inline auto WindowsMergerCacheIndices::size() const noexcept
    -> windows_size_type {
  return _lines_size;
}

inline auto WindowsMergerCacheIndices::capacity() const noexcept
    -> windows_size_type {
  return _lines_capacity;
}

inline auto WindowsMergerCacheIndices::line_capacity() const noexcept
    -> windows_size_type {
  return _line_capacity;
}

inline auto
WindowsMergerCacheIndices::operator[](windows_size_type line_index) noexcept
    -> indices_accessor {
  return indices_accessor(*this, line_index);
}

inline auto WindowsMergerCacheIndices::operator[](
    windows_size_type line_index) const noexcept -> const_indices_accessor {
  return const_indices_accessor(*this, line_index);
}

inline auto WindowsMergerCacheIndices::begin() noexcept -> iterator {
  return iterator{*this};
}

inline auto WindowsMergerCacheIndices::end() noexcept -> iterator {
  return iterator{*this, _lines_size};
}

inline auto WindowsMergerCacheIndices::rbegin() noexcept -> reverse_iterator {
  return reverse_iterator{iterator{*this, _lines_size}};
}

inline auto WindowsMergerCacheIndices::rend() noexcept -> reverse_iterator {
  return reverse_iterator{iterator{*this}};
}

inline auto WindowsMergerCacheIndices::begin() const noexcept
    -> const_iterator {
  return const_iterator{*this};
}

inline auto WindowsMergerCacheIndices::end() const noexcept -> const_iterator {
  return const_iterator{*this, _lines_size};
}

inline auto WindowsMergerCacheIndices::rbegin() const noexcept
    -> const_reverse_iterator {
  return const_reverse_iterator{const_iterator{*this, _lines_size}};
}

inline auto WindowsMergerCacheIndices::rend() const noexcept
    -> const_reverse_iterator {
  return const_reverse_iterator{const_iterator{*this}};
}

} // namespace windows_merger

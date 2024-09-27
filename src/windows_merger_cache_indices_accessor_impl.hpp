#pragma once

#include "windows_merger_cache_indices_accessor.hpp"
#include "windows_merger_cache_indices_line.hpp"

namespace windows_merger {

template <typename Merger>
WindowsMergerCacheIndicesAccessor<Merger>::WindowsMergerCacheIndicesAccessor(
    Merger &merger, windows_size_type line_index) noexcept
    : merger(&merger), line_index(line_index) {}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::index() const noexcept
    -> windows_size_type {
  return line_index;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::size() const noexcept
    -> windows_size_type {
  assert(line_index < merger->size());
  return *merger->template get_indices_first_pointer<0>(line_index);
}

template <typename Merger>
bool WindowsMergerCacheIndicesAccessor<Merger>::empty() const noexcept {
  return size() == 0;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::front() const noexcept
    -> reference {
  assert(size() > 0);
  return *merger->template get_indices_first_pointer<1>(line_index);
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::back() const noexcept
    -> reference {
  assert(size() > 0);
  return merger->template get_indices_first_pointer<1>(line_index)[size() - 1];
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::begin() const noexcept
    -> iterator {
  return merger->template get_indices_first_pointer<1>(line_index);
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::end() const noexcept
    -> iterator {
  return merger->template get_indices_first_pointer<1>(line_index) + size();
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::rbegin() const noexcept
    -> reverse_iterator {
  return ranges::reverse_iterator<iterator>(end());
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::rend() const noexcept
    -> reverse_iterator {
  return ranges::reverse_iterator<iterator>(begin());
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator[](
    windows_size_type index) const noexcept -> reference {
  assert(index < size());
  return merger->template get_indices_first_pointer<1>(line_index)[index];
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::emplace_back(
    windows_size_type index) const noexcept(false) -> windows_size_type & {
  return merger->index_emplace_back(line_index, index);
}

template <typename Merger>
void WindowsMergerCacheIndicesAccessor<Merger>::push_back(
    windows_size_type index) const noexcept(false) {
  merger->index_emplace_back(line_index, index);
}

template <typename Merger>
void WindowsMergerCacheIndicesAccessor<Merger>::resize(
    windows_size_type value) const noexcept(false) {
  if (value > merger->_line_capacity)
    merger->reshape(value, typename decayed_merger_type::reserver_type(0));

  static_assert(std::is_trivially_constructible_v<windows_size_type>);
  static_assert(std::is_trivially_destructible_v<windows_size_type>);
  *merger->template get_indices_first_pointer<0>(line_index) = value;
}

template <typename Merger>
void WindowsMergerCacheIndicesAccessor<Merger>::clear() const noexcept {
  static_assert(std::is_trivially_constructible_v<windows_size_type>);
  static_assert(std::is_trivially_destructible_v<windows_size_type>);
  *merger->template get_indices_first_pointer<0>(line_index) = 0;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices> const &rhs)
    const noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices const> const
        &rhs) const noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices &&> const &rhs)
    const noexcept(false) -> self const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices const> &&rhs)
    const noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesAccessor<WindowsMergerCacheIndices &&> &&rhs) const
    noexcept(false) -> self const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesLine const &line) const
    noexcept(false) -> self const & {
  assign_from_line(line);
  return *this;
}

template <typename Merger>
auto WindowsMergerCacheIndicesAccessor<Merger>::operator=(
    WindowsMergerCacheIndicesLine &&line) const
    noexcept(false) -> self const & {
  assign_from_line(std::move(line));
  return *this;
}

template <typename Merger>
WindowsMergerCacheIndicesAccessor<
    Merger>::operator WindowsMergerCacheIndicesLine() const noexcept(false) {
  return WindowsMergerCacheIndicesLine(*this);
}

template <typename Merger>
template <typename Accessor>
void WindowsMergerCacheIndicesAccessor<Merger>::assign_from_accessor(
    Accessor &&accessor) const noexcept(false) {
  assert(line_index < merger->_lines_size);
  assert(accessor.line_index < accessor.merger->_lines_size);

  const auto rhs_size = *accessor.merger->template get_indices_first_pointer<0>(
      accessor.line_index);

  if (merger->_line_capacity < rhs_size)
    merger->reshape(rhs_size, typename decayed_merger_type::reserver_type(0));

  const auto lhs_size =
      merger->template get_indices_first_pointer<0>(line_index);
  static_assert(std::is_trivially_destructible_v<windows_size_type>);

  const auto rhs_indices =
      accessor.merger->template get_indices_first_pointer<1>(
          accessor.line_index);
  const auto lhs_indices =
      merger->template get_indices_first_pointer<1>(line_index);
  std::copy(rhs_indices, rhs_indices + rhs_size, lhs_indices);
  *lhs_size = rhs_size;
}

template <typename Merger>
template <typename Line>
void WindowsMergerCacheIndicesAccessor<Merger>::assign_from_line(
    Line &&line) const noexcept(false) {
  assert(line_index < merger->_lines_size);

  const auto rhs_size = line.size();

  if (merger->_line_capacity < rhs_size)
    merger->reshape(rhs_size, typename decayed_merger_type::reserver_type(0));

  const auto lhs_size =
      merger->template get_indices_first_pointer<0>(line_index);
  static_assert(std::is_trivially_destructible_v<windows_size_type>);

  const auto lhs_indices =
      merger->template get_indices_first_pointer<1>(line_index);
  ranges::copy(line, lhs_indices);
  *lhs_size = rhs_size;
}

} // namespace windows_merger

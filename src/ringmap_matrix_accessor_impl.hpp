#pragma once

#include "ringmap_matrix_accessor.hpp"

#include <algorithm>
#include <cassert>
#include <functional>
#include <tuple>

template <typename Matrix>
template <typename Row, typename>
RingmapMatrixAccessor<Matrix>::RingmapMatrixAccessor(Row row,
                                                     unsigned col) noexcept
    : row(row), col(col) {}

template <typename Matrix>
RingmapMatrixAccessor<Matrix>::operator value_type() const noexcept {
  return static_cast<value_type>(find().second);
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(value_type value) const
    noexcept(false) -> RingmapMatrixAccessor const & {

  if (auto &&[iter, found] = find(); not value and found) {
    auto &&[begin, end] =
        std::ranges::remove(iter, std::ranges::end(*row), col);
    row->erase(begin, end);
  } else if (value and not found) {
    row->emplace_back(col);
    std::ranges::sort(*row);
  }

  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<RingmapMatrix> const &rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<const RingmapMatrix> const &rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<RingmapMatrix &&> const &rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<RingmapMatrix> &&rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<const RingmapMatrix> &&rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::operator=(
    RingmapMatrixAccessor<RingmapMatrix &&> &&rhs) const noexcept
    -> RingmapMatrixAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixAccessor<Matrix>::find() const noexcept
    -> std::pair<row_iterator, bool> {
  assert(std::ranges::is_sorted(*row));
  auto iter = std::ranges::upper_bound(*row, col, std::less_equal<unsigned>());
  bool found = (iter != std::ranges::end(*row) and *iter == col);
  return {std::move(iter), found};
}

template <typename Matrix>
template <typename Accessor>
void RingmapMatrixAccessor<Matrix>::copy_from_accessor(
    Accessor &&accessor) const noexcept {
  *this = static_cast<value_type>(accessor);
}

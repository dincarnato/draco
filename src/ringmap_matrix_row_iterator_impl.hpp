#pragma once

#include "ringmap_matrix_row_accessor.hpp"
#include "ringmap_matrix_row_iterator.hpp"

#include <cassert>
#include <numeric>

template <typename Matrix>
RingmapMatrixRowIterator<Matrix>::RingmapMatrixRowIterator(
    Matrix& matrix, row_type* const row) noexcept
    : matrix(&matrix), row(row) {}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator++() noexcept -> self& {
  ++row;
  return *this;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator++(int) noexcept -> self {
  auto other = *this;
  ++row;
  return other;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator--() noexcept -> self& {
  --row;
  return *this;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator--(int) noexcept -> self {
  auto other = *this;
  --row;
  return other;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator+=(difference_type offset) noexcept
    -> self& {
  row += offset;
  return *this;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator+(difference_type offset) const
    noexcept -> self {
  return self{*matrix, row + offset};
}

template <typename Matrix>
RingmapMatrixRowIterator<Matrix>
operator+(typename RingmapMatrixRowIterator<Matrix>::difference_type offset,
          RingmapMatrixRowIterator<Matrix> const& iter) noexcept {
  return iter + offset;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator-=(difference_type offset) noexcept
    -> self& {
  row -= offset;
  return *this;
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator-(difference_type offset) const
    noexcept -> self {
  return self{*matrix, row - offset};
}

template <typename Matrix>
auto
RingmapMatrixRowIterator<Matrix>::operator-(self const& rhs) const noexcept
    -> difference_type {
  assert(matrix == rhs.matrix);
  assert(matrix);
  return row - rhs.row;
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator==(const RingmapMatrixRowIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row == other.row;
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator!=(const RingmapMatrixRowIterator& other) const noexcept {
  return not operator==(other);
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator<(const RingmapMatrixRowIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row < other.row;
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator<=(const RingmapMatrixRowIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row <= other.row;
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator>=(const RingmapMatrixRowIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row >= other.row;
}

template <typename Matrix>
bool
RingmapMatrixRowIterator<Matrix>::
operator>(const RingmapMatrixRowIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row > other.row;
}

template <typename Matrix>
auto RingmapMatrixRowIterator<Matrix>::operator*() const noexcept -> reference {
  return {*matrix, row};
}

template <typename Matrix>
auto RingmapMatrixRowIterator<Matrix>::operator[](difference_type offset) const
    noexcept -> reference {
  return {*matrix, row + offset};
}

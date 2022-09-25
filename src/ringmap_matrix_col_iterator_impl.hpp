#pragma once

#include "ringmap_matrix_col_accessor.hpp"
#include "ringmap_matrix_col_iterator.hpp"

#include <cassert>

template <typename Matrix>
RingmapMatrixColIterator<Matrix>::RingmapMatrixColIterator(
    Matrix &matrix, unsigned col) noexcept
    : matrix(&matrix), col(col) {}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator++() noexcept -> self & {
  ++col;
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator++(int) noexcept -> self {
  auto other = *this;
  ++col;
  return other;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator--() noexcept -> self & {
  --col;
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator--(int) noexcept -> self {
  auto other = *this;
  --col;
  return other;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator+=(
    difference_type offset) noexcept -> self & {
  col += offset;
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator+(
    difference_type offset) const noexcept -> self {
  return self{*matrix, col + offset};
}

template <typename Matrix>
RingmapMatrixColIterator<Matrix>
operator+(typename RingmapMatrixColIterator<Matrix>::difference_type offset,
          RingmapMatrixColIterator<Matrix> const &iter) noexcept {
  return iter + offset;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator-=(
    difference_type offset) noexcept -> self & {
  col -= offset;
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator-(
    difference_type offset) const noexcept -> self {
  return self{*matrix, col - offset};
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator-(self const &rhs) const noexcept
    -> difference_type {
  assert(matrix == rhs.matrix);
  assert(matrix);
  return col - rhs.col;
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator==(
    const RingmapMatrixColIterator &other) const noexcept {
  assert(matrix == other.matrix);
  return col == other.col;
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator!=(
    const RingmapMatrixColIterator &other) const noexcept {
  return not operator==(other);
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator<(
    const RingmapMatrixColIterator &other) const noexcept {
  assert(matrix == other.matrix);
  return col < other.col;
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator<=(
    const RingmapMatrixColIterator &other) const noexcept {
  assert(matrix == other.matrix);
  return col <= other.col;
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator>=(
    const RingmapMatrixColIterator &other) const noexcept {
  assert(matrix == other.matrix);
  return col >= other.col;
}

template <typename Matrix>
bool RingmapMatrixColIterator<Matrix>::operator>(
    const RingmapMatrixColIterator &other) const noexcept {
  assert(matrix == other.matrix);
  return col > other.col;
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator*() const noexcept -> reference {
  return {*matrix, col};
}

template <typename Matrix>
auto RingmapMatrixColIterator<Matrix>::operator[](
    difference_type offset) const noexcept -> reference {
  return {*matrix, col + offset};
}

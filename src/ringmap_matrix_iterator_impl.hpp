#pragma once

#include "ringmap_matrix_iterator.hpp"

#include <cassert>

template <typename Matrix, RingmapMatrixIteratorType Type>
RingmapMatrixIterator<Matrix, Type>::RingmapMatrixIterator(
    Matrix& matrix, row_type* const row, unsigned col) noexcept
    : matrix(&matrix), row(row), col(col) {}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator++() noexcept -> self& {
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    ++col;
  else if constexpr (Type == RingmapMatrixIteratorType::rows)
    ++row;
  else {
    ++col;
    if (col >= matrix->bases) {
      col = 0;
      ++row;
    }
  }

  return *this;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator++(int) noexcept -> self {
  auto other = *this;
  ++*this;

  return other;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator--() noexcept -> self& {
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    --col;
  else if constexpr (Type == RingmapMatrixIteratorType::rows)
    --row;
  else {
    if (col == 0) {
      assert(matrix->bases > 0);
      --row;
      col = matrix->bases - 1;
    } else
      --col;
  }

  return *this;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator--(int) noexcept -> self {
  auto other = *this;
  --*this;

  return other;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator+=(difference_type offset) noexcept
    -> self& {
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    col += offset;
  else if constexpr (Type == RingmapMatrixIteratorType::rows)
    row += offset;
  else {
    col += offset;
    row += col / matrix->bases;
    col %= matrix->bases;
  }

  return *this;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator+(difference_type offset) const
    noexcept -> self {
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    return {*matrix, row, col + offset};
  else if constexpr (Type == RingmapMatrixIteratorType::cols)
    return {*matrix, row + offset, col};
  else {
    row_type* otherRow = row;
    unsigned otherCol = col + offset;

    otherRow += otherCol / matrix->bases;
    otherCol %= matrix->bases;

    return {*matrix, otherRow, otherCol};
  }
}

template <typename Matrix, RingmapMatrixIteratorType Type>
RingmapMatrixIterator<Matrix, Type>
operator+(typename RingmapMatrixIterator<Matrix, Type>::difference_type offset,
          RingmapMatrixIterator<Matrix, Type> const& rhs) noexcept {
  return rhs + offset;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator-=(difference_type offset) noexcept
    -> self& {
  return *this += -offset;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator-(difference_type offset) const
    noexcept -> self {
  return *this + -offset;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto
RingmapMatrixIterator<Matrix, Type>::operator-(self const& rhs) const noexcept
    -> difference_type {
  assert(matrix == rhs.matrix);
  if constexpr (Type == RingmapMatrixIteratorType::cols) {
    return static_cast<difference_type>(col) -
           static_cast<difference_type>(rhs.col);
  } else if constexpr (Type == RingmapMatrixIteratorType::rows) {
    return row - rhs.row;
  } else {
    return (row - rhs.row) * matrix->bases +
           (static_cast<difference_type>(col) -
            static_cast<difference_type>(rhs.col));
  }
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator==(const RingmapMatrixIterator& other) const noexcept {
  assert(matrix == other.matrix);
  return row == other.row and col == other.col;
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator!=(const RingmapMatrixIterator& other) const noexcept {
  return not operator==(other);
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator<(const RingmapMatrixIterator& other) const noexcept {
  assert(matrix == other.matrix);
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    return col < other.col;
  else if constexpr (Type == RingmapMatrixIteratorType::rows)
    return row < other.row;
  else {
    if (row == other.row)
      return col < other.col;
    return row < other.row;
  }
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator<=(const RingmapMatrixIterator& other) const noexcept {
  assert(matrix == other.matrix);
  if constexpr (Type == RingmapMatrixIteratorType::cols)
    return col <= other.col;
  else if constexpr (Type == RingmapMatrixIteratorType::rows)
    return row <= other.row;
  else {
    if (row == other.row)
      return col <= other.col;
    return row <= other.row;
  }
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator>=(const RingmapMatrixIterator& other) const noexcept {
  return not(*this < other);
}

template <typename Matrix, RingmapMatrixIteratorType Type>
bool
RingmapMatrixIterator<Matrix, Type>::
operator>(const RingmapMatrixIterator& other) const noexcept {
  return not(*this <= other);
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto RingmapMatrixIterator<Matrix, Type>::operator*() const noexcept
    -> reference {
  return {row, col};
}

template <typename Matrix, RingmapMatrixIteratorType Type>
auto RingmapMatrixIterator<Matrix, Type>::
operator[](difference_type offset) const noexcept -> reference {
  auto temp = *this;
  temp += offset;
  return *temp;
}

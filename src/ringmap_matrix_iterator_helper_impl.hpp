#pragma once

#include "ringmap_matrix_col_iterator.hpp"
#include "ringmap_matrix_iterator_helper.hpp"
#include "ringmap_matrix_row_iterator.hpp"

namespace detail {

template <typename Matrix>
RingmapMatrixRowIteratorHelper<Matrix>::RingmapMatrixRowIteratorHelper(
    Matrix& matrix) noexcept
    : matrix(&matrix) {}

template <typename Matrix>
auto
RingmapMatrixRowIteratorHelper<Matrix>::begin() const noexcept -> iterator {
  return {*matrix, matrix->data.data()};
}

template <typename Matrix>
auto
RingmapMatrixRowIteratorHelper<Matrix>::end() const noexcept -> iterator {
  return {*matrix, matrix->data.data() + matrix->readsCount};
}

template <typename Matrix>
RingmapMatrixColIteratorHelper<Matrix>::RingmapMatrixColIteratorHelper(
    Matrix& matrix) noexcept
    : matrix(&matrix) {}

template <typename Matrix>
auto
RingmapMatrixColIteratorHelper<Matrix>::begin() const noexcept -> iterator {
  return {*matrix, 0};
}

template <typename Matrix>
auto
RingmapMatrixColIteratorHelper<Matrix>::end() const noexcept -> iterator {
  return {*matrix, matrix->bases};
}

} /* namespace detail */

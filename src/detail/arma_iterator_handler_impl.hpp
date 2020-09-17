#pragma once

#include "arma_iterator_handler.hpp"
#include "arma_iterator_helper.hpp"

namespace detail {

template <typename Mat, ArmaIteratorDirection direction>
ArmaIteratorHandler<Mat, direction>::ArmaIteratorHandler(Mat* matrix) noexcept
    : matrix(matrix) {}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaIteratorHandler<Mat, direction>::begin() const noexcept -> iterator {
  return {matrix};
}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaIteratorHandler<Mat, direction>::end() const noexcept -> iterator {
  return {matrix, 0};
}

template <typename Mat, ArmaIteratorDirection direction>
std::size_t
ArmaIteratorHandler<Mat, direction>::size() const noexcept {
  if constexpr (direction == ArmaIteratorDirection::rows)
    return matrix->n_rows;
  else {
    static_assert(direction == ArmaIteratorDirection::cols);
    return matrix->n_cols;
  }
}

} // namespace detail

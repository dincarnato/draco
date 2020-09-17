#pragma once

class RingmapMatrix;
template <typename>
class RingmapMatrixRowIterator;
template <typename>
class RingmapMatrixColIterator;

namespace detail {

template <typename Matrix>
class RingmapMatrixRowIteratorHelper {
public:
  RingmapMatrixRowIteratorHelper() = default;
  RingmapMatrixRowIteratorHelper(Matrix& matrix) noexcept;

  using pointer = std::remove_reference_t<Matrix>*;
  using iterator = RingmapMatrixRowIterator<Matrix>;

  iterator begin() const noexcept;
  iterator end() const noexcept;

private:
  pointer matrix;
};

template <typename Matrix>
class RingmapMatrixColIteratorHelper {
public:
  RingmapMatrixColIteratorHelper() = default;
  RingmapMatrixColIteratorHelper(Matrix& matrix) noexcept;

  using pointer = std::remove_reference_t<Matrix>*;
  using iterator = RingmapMatrixColIterator<Matrix>;

  iterator begin() const noexcept;
  iterator end() const noexcept;

private:
  pointer matrix;
};

} /* namespace detail */

#include "ringmap_matrix_iterator_helper_impl.hpp"

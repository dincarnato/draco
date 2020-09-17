#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_traits.hpp"

#include <iterator>
#include <type_traits>

template <typename>
class RingmapMatrixColAccessor;
namespace detail {
template <typename>
class RingmapMatrixColIteratorHelper;
}

template <typename Matrix>
class RingmapMatrixColIterator {
  friend class RingmapMatrix;
  template <typename>
  friend class detail::RingmapMatrixColIteratorHelper;

public:
  using pointer = std::remove_reference_t<Matrix>*;
  using value_type = std::vector<bool>;
  using reference = RingmapMatrixColAccessor<Matrix>;
  using row_type = std::conditional_t<std::is_const<Matrix>::value,
                                      const ringmap_matrix::row_type,
                                      typename ringmap_matrix::row_type>;
  using difference_type = std::ptrdiff_t;
  using iterator_category = ranges::random_access_iterator_tag;
  using self = RingmapMatrixColIterator;

  RingmapMatrixColIterator() = default;
  RingmapMatrixColIterator(Matrix& matrix, unsigned col) noexcept;

  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;
  self& operator+=(difference_type offset) noexcept;
  self operator+(difference_type offset) const noexcept;

  template <typename _Matrix>
  friend RingmapMatrixColIterator<_Matrix>
  operator+(typename RingmapMatrixColIterator<_Matrix>::difference_type offset,
            RingmapMatrixColIterator<_Matrix> const& iter) noexcept;

  self& operator-=(difference_type offset) noexcept;
  self operator-(difference_type offset) const noexcept;
  difference_type operator-(self const& rhs) const noexcept;

  bool operator==(const self& other) const noexcept;
  bool operator!=(const self& other) const noexcept;
  bool operator<(const self& other) const noexcept;
  bool operator<=(const self& other) const noexcept;
  bool operator>=(const self& other) const noexcept;
  bool operator>(const self& other) const noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type offset) const noexcept;

private:
  Matrix* matrix = nullptr;
  unsigned col;
};

#include "ringmap_matrix_col_iterator_impl.hpp"

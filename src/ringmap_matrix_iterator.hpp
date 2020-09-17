#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_traits.hpp"

#include <range/v3/core.hpp>

template <typename>
class RingmapMatrixAccessor;

template <typename Matrix, RingmapMatrixIteratorType Type>
class RingmapMatrixIterator {
  template <typename>
  friend class RingmapMatrixRowAccessor;
  template <typename>
  friend class RingmapMatrixColAccessor;
  friend class RingmapMatrix;

public:
  static constexpr RingmapMatrixIteratorType iterator_type = Type;

  using value_type = ringmap_matrix::value_type;
  using pointer = std::remove_reference_t<Matrix>*;
  using difference_type = std::ptrdiff_t;
  using reference = RingmapMatrixAccessor<Matrix>;
  using iterator_category = ranges::random_access_iterator_tag;
  using row_type = std::conditional_t<std::is_const<Matrix>::value,
                                      const ringmap_matrix::row_type,
                                      ringmap_matrix::row_type>;
  using self = RingmapMatrixIterator;

  RingmapMatrixIterator() = default;
  RingmapMatrixIterator(Matrix& matrix, row_type* const row,
                        unsigned col) noexcept;

  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;
  self& operator+=(difference_type offset) noexcept;
  self operator+(difference_type offset) const noexcept;

  template <typename _Matrix, RingmapMatrixIteratorType _Type>
  friend RingmapMatrixIterator<_Matrix, _Type> operator+(
      typename RingmapMatrixIterator<_Matrix, _Type>::difference_type offset,
      RingmapMatrixIterator<_Matrix, _Type> const& iter) noexcept;

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
  pointer matrix = nullptr;
  row_type* row = nullptr;
  unsigned col;
};

#include "ringmap_matrix_iterator_impl.hpp"

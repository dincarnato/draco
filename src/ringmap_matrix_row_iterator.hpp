#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_traits.hpp"

#include <iterator>
#include <range/v3/core.hpp>
#include <type_traits>
#include <vector>

template <typename>
class RingmapMatrixRowAccessor;
namespace detail {
template <typename>
class RingmapMatrixRowIteratorHelper;
}

template <typename Matrix>
class RingmapMatrixRowIterator {
  friend class RingmapMatrix;
  template <typename>
  friend class detail::RingmapMatrixRowIteratorHelper;

public:
  using pointer = std::remove_reference_t<Matrix>*;
  using value_type = std::vector<bool>;
  using reference = RingmapMatrixRowAccessor<Matrix>;
  using row_type = std::conditional_t<std::is_const<Matrix>::value,
                                      const ringmap_matrix::row_type,
                                      typename ringmap_matrix::row_type>;
  using difference_type = std::ptrdiff_t;
  using iterator_category = ranges::random_access_iterator_tag;
  using self = RingmapMatrixRowIterator;

  RingmapMatrixRowIterator() = default;
  RingmapMatrixRowIterator(Matrix& matrix, row_type* const row) noexcept;

  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;
  self& operator+=(difference_type offset) noexcept;
  self operator+(difference_type offset) const noexcept;

  template <typename _Matrix>
  friend RingmapMatrixRowIterator<_Matrix>
  operator+(typename RingmapMatrixRowIterator<_Matrix>::difference_type offset,
            RingmapMatrixRowIterator<_Matrix> const& iter) noexcept;

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
};

#include "ringmap_matrix_row_iterator_impl.hpp"

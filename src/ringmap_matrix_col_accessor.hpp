#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_traits.hpp"

#include <type_traits>

template <typename, RingmapMatrixIteratorType> class RingmapMatrixIterator;
class RingmapMatrix;

template <typename Matrix> class RingmapMatrixColAccessor {
  template <typename> friend class RingmapMatrixColIterator;
  friend class RingmapMatrix;

public:
  using matrix_type = Matrix;
  using pointer = std::remove_reference_t<Matrix> *;
  using value_type = std::vector<bool>;
  using iterator =
      RingmapMatrixIterator<Matrix, RingmapMatrixIteratorType::rows>;
  using row_type = std::conditional_t<std::is_const<Matrix>::value,
                                      const ringmap_matrix::row_type,
                                      ringmap_matrix::row_type>;

  RingmapMatrixColAccessor() = default;
  RingmapMatrixColAccessor(Matrix &matrix, unsigned col) noexcept;
  RingmapMatrixColAccessor(RingmapMatrixColAccessor const &) = default;
  RingmapMatrixColAccessor(RingmapMatrixColAccessor &&) = default;

  RingmapMatrixColAccessor const &operator=(value_type const &column) const
      noexcept(false);
  operator value_type() const noexcept(false);

  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<RingmapMatrix> const &rhs) const
      noexcept(false);
  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<const RingmapMatrix> const &rhs) const
      noexcept(false);
  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<RingmapMatrix &&> const &rhs) const
      noexcept(false);
  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<RingmapMatrix> &&rhs) const
      noexcept(false);
  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<const RingmapMatrix> &&rhs) const
      noexcept(false);
  RingmapMatrixColAccessor const &
  operator=(RingmapMatrixColAccessor<RingmapMatrix &&> &&rhs) const
      noexcept(false);

  iterator begin() const noexcept;
  iterator end() const noexcept;

  unsigned size() const noexcept;
  std::size_t sum() const noexcept;
  double mean() const noexcept;
  template <typename URBG> void shuffle(URBG &&g) const noexcept(false);

  RingmapMatrixAccessor<Matrix> operator[](std::size_t index) const noexcept;

private:
  template <typename Accessor>
  void copy_from_accessor(Accessor &&accessor) const noexcept(false);

  pointer matrix = nullptr;
  unsigned col;
};

#include "ringmap_matrix_col_accessor_impl.hpp"

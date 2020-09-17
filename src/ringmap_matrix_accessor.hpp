#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_row.hpp"
#include "ringmap_matrix_traits.hpp"

#include <type_traits>
#include <utility>

class RingmapMatrix;

template <typename Matrix>
class RingmapMatrixAccessor {
  friend class RingmapMatrix;
  template <typename>
  friend class RingmapMatrixRowIterator;
  template <typename>
  friend class RingmapMatrixColIterator;
  template <typename>
  friend class RingmapMatrixRowAccessor;
  template <typename>
  friend class RingmapMatrixColAccessor;
  template <typename, RingmapMatrixIteratorType>
  friend class RingmapMatrixIterator;

public:
  using value_type = ringmap_matrix::value_type;
  using matrix_type = Matrix;
  using row_type = std::conditional_t<std::is_const_v<Matrix>,
                                      const ringmap_matrix::row_type,
                                      ringmap_matrix::row_type>;
  using row_iterator = std::conditional_t<std::is_const_v<Matrix>,
                                          typename row_type::const_iterator,
                                          typename row_type::iterator>;

  RingmapMatrixAccessor() = default;
  template <typename Row, typename = std::enable_if_t<std::is_same<
                              Row, std::add_pointer_t<row_type>>::value>>
  RingmapMatrixAccessor(Row row, unsigned col) noexcept;
  RingmapMatrixAccessor(RingmapMatrixAccessor const&) = default;
  RingmapMatrixAccessor(RingmapMatrixAccessor&&) = default;

  std::pair<row_iterator, bool> find() const noexcept;

  RingmapMatrixAccessor const& operator=(value_type value) const
      noexcept(false);

  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<RingmapMatrix> const& rhs) const noexcept;
  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<const RingmapMatrix> const& rhs) const
      noexcept;
  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<RingmapMatrix&&> const& rhs) const noexcept;
  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<RingmapMatrix>&& rhs) const noexcept;
  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<const RingmapMatrix>&& rhs) const noexcept;
  RingmapMatrixAccessor const&
  operator=(RingmapMatrixAccessor<RingmapMatrix&&>&& rhs) const noexcept;

  operator value_type() const noexcept;

private:
  template <typename Accessor>
  void copy_from_accessor(Accessor&& accessor) const noexcept;

  row_type* row;
  unsigned col;
};

#include "ringmap_matrix_accessor_impl.hpp"

#pragma once

#include "ringmap_matrix_iterator_type.hpp"
#include "ringmap_matrix_traits.hpp"

#include <type_traits>

#include <armadillo>

class RingmapMatrix;
struct RingmapMatrixRow;
template <typename, RingmapMatrixIteratorType> class RingmapMatrixIterator;
template <typename> class RingmapMatrixAccessor;

template <typename Matrix> class RingmapMatrixRowAccessor {
  template <typename> friend class RingmapMatrixRowIterator;
  template <typename> friend class RingmapMatrixRowAccessor;
  friend class RingmapMatrix;

public:
  using matrix_type = Matrix;
  using pointer = std::remove_reference_t<Matrix> *;
  using value_type = std::vector<bool>;
  using iterator =
      RingmapMatrixIterator<Matrix, RingmapMatrixIteratorType::cols>;
  using row_type = std::conditional_t<std::is_const<Matrix>::value,
                                      const ringmap_matrix::row_type,
                                      ringmap_matrix::row_type>;
  using reference = std::conditional_t<std::is_rvalue_reference_v<Matrix>,
                                       ringmap_matrix::row_type &&, row_type &>;

  RingmapMatrixRowAccessor() noexcept = default;
  RingmapMatrixRowAccessor(RingmapMatrixRowAccessor const &) = default;
  RingmapMatrixRowAccessor(RingmapMatrixRowAccessor &&) = default;

  template <typename M = Matrix>
  RingmapMatrixRowAccessor(
      Matrix &matrix,
      std::enable_if_t<std::is_const<M>::value, row_type *const> row) noexcept;

  template <typename M = Matrix>
  RingmapMatrixRowAccessor(
      Matrix &matrix,
      std::enable_if_t<not std::is_const<M>::value, row_type *const>
          row) noexcept;

  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<RingmapMatrix> const &rhs) const
      noexcept(false);
  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<const RingmapMatrix> const &rhs) const
      noexcept(false);
  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<RingmapMatrix &&> const &rhs) const
      noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>);
  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<RingmapMatrix> &&rhs) const
      noexcept(false);
  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<const RingmapMatrix> &&rhs) const
      noexcept(false);
  RingmapMatrixRowAccessor const &
  operator=(RingmapMatrixRowAccessor<RingmapMatrix &&> &&rhs) const
      noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>);

  RingmapMatrixRowAccessor const &operator=(value_type const &row) const
      noexcept(false);
  RingmapMatrixRowAccessor const &operator=(value_type &&row) const
      noexcept(false);
  operator value_type() const noexcept(false);

  constexpr ringmap_matrix::base_index_type begin_index() const noexcept;
  constexpr ringmap_matrix::base_index_type end_index() const noexcept;

  void set_begin_index(ringmap_matrix::base_index_type value) noexcept;
  void set_end_index(ringmap_matrix::base_index_type value) noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

  std::size_t sum() const noexcept;
  double mean() const noexcept;
  template <typename URBG>
  void shuffle(URBG &&g) const
      noexcept(std::is_nothrow_swappable_v<ringmap_matrix::value_type>);

  template <typename T> T convTo() const noexcept(false);
  arma::rowvec operator*(const arma::mat &eigenVectors) const noexcept(false);

  constexpr reference modifiedIndices() const noexcept;

  RingmapMatrixAccessor<Matrix> operator[](unsigned index) const noexcept;
  template <typename _Matrix>
  bool
  operator==(const RingmapMatrixRowAccessor<_Matrix> &other) const noexcept;

  template <typename ML, typename MR>
  friend bool operator<(RingmapMatrixRowAccessor<ML> const &lhs,
                        RingmapMatrixRowAccessor<MR> const &rhs) noexcept;

  template <typename M>
  friend bool operator<(RingmapMatrixRowAccessor<M> const &lhs,
                        ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator<(ringmap_matrix::row_type const &lhs,
                        RingmapMatrixRowAccessor<M> const &rhs) noexcept;

  template <typename ML, typename MR>
  friend bool operator<=(RingmapMatrixRowAccessor<ML> const &lhs,
                         RingmapMatrixRowAccessor<MR> const &rhs) noexcept;

  template <typename M>
  friend bool operator<=(RingmapMatrixRowAccessor<M> const &lhs,
                         ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator<=(ringmap_matrix::row_type const &lhs,
                         RingmapMatrixRowAccessor<M> const &rhs) noexcept;

  template <typename ML, typename MR>
  friend bool operator>=(RingmapMatrixRowAccessor<ML> const &lhs,
                         RingmapMatrixRowAccessor<MR> const &rhs) noexcept;

  template <typename M>
  friend bool operator>=(RingmapMatrixRowAccessor<M> const &lhs,
                         ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator>=(ringmap_matrix::row_type const &lhs,
                         RingmapMatrixRowAccessor<M> const &rhs) noexcept;

  template <typename ML, typename MR>
  friend bool operator>(RingmapMatrixRowAccessor<ML> const &lhs,
                        RingmapMatrixRowAccessor<MR> const &rhs) noexcept;

  template <typename M>
  friend bool operator>(RingmapMatrixRowAccessor<M> const &lhs,
                        ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator>(ringmap_matrix::row_type const &lhs,
                        RingmapMatrixRowAccessor<M> const &rhs) noexcept;

  template <typename M>
  friend bool operator==(RingmapMatrixRowAccessor<M> const &lhs,
                         ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator==(ringmap_matrix::row_type const &lhs,
                         RingmapMatrixRowAccessor<M> const &rhs) noexcept;

  template <typename M>
  friend bool operator!=(RingmapMatrixRowAccessor<M> const &lhs,
                         ringmap_matrix::row_type const &rhs) noexcept;

  template <typename M>
  friend bool operator!=(ringmap_matrix::row_type const &lhs,
                         RingmapMatrixRowAccessor<M> const &rhs) noexcept;

private:
  template <typename Accessor>
  void assign_from_accessor(Accessor &&accessor) const noexcept(false);
  template <typename Accessor>
  void move_from_accessor(Accessor &&accessor) const
      noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>);
  template <typename Row> void assign_from_row(Row &&row) const noexcept(false);

  pointer matrix = nullptr;
  row_type *row = nullptr;
};

#include "ringmap_matrix_row_accessor_impl.hpp"

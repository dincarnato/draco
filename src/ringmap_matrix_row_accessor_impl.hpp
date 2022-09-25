#pragma once

#include "ringmap_matrix_iterator.hpp"
#include "ringmap_matrix_row_accessor.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

template <typename Matrix>
template <typename M>
RingmapMatrixRowAccessor<Matrix>::RingmapMatrixRowAccessor(
    Matrix &matrix,
    std::enable_if_t<std::is_const<M>::value, row_type *const> row) noexcept
    : matrix(&matrix), row(row) {}

template <typename Matrix>
template <typename M>
RingmapMatrixRowAccessor<Matrix>::RingmapMatrixRowAccessor(
    Matrix &matrix,
    std::enable_if_t<not std::is_const<M>::value, row_type *const> row) noexcept
    : matrix(&matrix), row(row) {
  matrix.readsCount = std::max(
      static_cast<unsigned>(std::distance(matrix.data.data(), row)) + 1,
      matrix.readsCount);
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<RingmapMatrix> const &rhs) const noexcept(false)
    -> RingmapMatrixRowAccessor const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<const RingmapMatrix> const &rhs) const
    noexcept(false) -> RingmapMatrixRowAccessor const & {
  assign_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<RingmapMatrix &&> const &rhs) const
    noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>)
        -> RingmapMatrixRowAccessor const & {
  move_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<RingmapMatrix> &&rhs) const noexcept(false)
    -> RingmapMatrixRowAccessor const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<const RingmapMatrix> &&rhs) const noexcept(false)
    -> RingmapMatrixRowAccessor const & {
  assign_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(
    RingmapMatrixRowAccessor<RingmapMatrix &&> &&rhs) const
    noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>)
        -> RingmapMatrixRowAccessor const & {
  move_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::begin() const noexcept -> iterator {
  assert(matrix);
  return {*matrix, row, 0};
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::end() const noexcept -> iterator {
  assert(matrix);
  return {*matrix, row, matrix->bases};
}

template <typename Matrix>
std::size_t RingmapMatrixRowAccessor<Matrix>::sum() const noexcept {
  return ranges::accumulate(begin(), end(), 0u);
}

template <typename Matrix>
double RingmapMatrixRowAccessor<Matrix>::mean() const noexcept {
  assert(matrix);
  return static_cast<double>(sum()) / matrix->bases;
}

template <typename Matrix>
arma::rowvec
RingmapMatrixRowAccessor<Matrix>::operator*(const arma::mat &eigenVectors) const
    noexcept(false) {
  assert(matrix);
  assert(eigenVectors.n_rows == matrix->cols_size());
  arma::rowvec proj(eigenVectors.n_cols, arma::fill::zeros);
  for (unsigned index : *row) {
    for (unsigned col = 0; col < eigenVectors.n_cols; ++col)
      proj[col] += eigenVectors(index, col);
  }

  return proj;
}

template <typename Matrix>
template <typename T>
T RingmapMatrixRowAccessor<Matrix>::convTo() const noexcept(false) {
  assert(matrix);
  T out(matrix->bases, arma::fill::zeros);
  for (unsigned index : *row)
    out[index] = 1;
  return out;
}

template <typename Matrix>
template <typename URBG>
void RingmapMatrixRowAccessor<Matrix>::shuffle(URBG &&g) const
    noexcept(std::is_nothrow_swappable_v<ringmap_matrix::value_type>) {
  std::uniform_int_distribution<unsigned> distribution(0, matrix->bases - 1);
  auto iter = ranges::begin(*row);
  for (unsigned &index : *row) {
    do {
      index = distribution(g);
    } while (ranges::find(ranges::begin(*row), iter, index) != iter);

    ++iter;
  }
  ranges::sort(*row);
}

template <typename Matrix>
RingmapMatrixAccessor<Matrix>
RingmapMatrixRowAccessor<Matrix>::operator[](unsigned index) const noexcept {
  assert(matrix);
  return {row, index};
}

template <typename Matrix>
template <typename _Matrix>
bool RingmapMatrixRowAccessor<Matrix>::operator==(
    const RingmapMatrixRowAccessor<_Matrix> &other) const noexcept {
  static_assert(std::is_same_v<std::decay_t<_Matrix>, RingmapMatrix>);
  assert(matrix);
  return *row == *other.row;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::modifiedIndices() const noexcept
    -> reference {
  assert(matrix);
  return *row;
}

template <typename Matrix>
template <typename Accessor>
void RingmapMatrixRowAccessor<Matrix>::assign_from_accessor(
    Accessor &&accessor) const noexcept(false) {

  assert(matrix);
  assert(matrix->bases == accessor.matrix->bases);
  *row = *accessor.row;
}

template <typename Matrix>
template <typename Accessor>
void RingmapMatrixRowAccessor<Matrix>::move_from_accessor(
    Accessor &&accessor) const
    noexcept(std::is_nothrow_move_assignable_v<ringmap_matrix::row_type>) {

  assert(matrix);
  assert(matrix->bases == accessor.matrix->bases);
  *row = std::move(*accessor.row);
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(value_type const &row) const
    noexcept(false) -> RingmapMatrixRowAccessor const & {
  assign_from_row(row);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixRowAccessor<Matrix>::operator=(value_type &&row) const
    noexcept(false) -> RingmapMatrixRowAccessor const & {
  assign_from_row(std::move(row));
  return *this;
}

template <typename Matrix>
RingmapMatrixRowAccessor<Matrix>::operator value_type() const noexcept(false) {
  value_type out(matrix->bases, false);
  ranges::view::zip(*this, out) |
      ranges::for_each([](auto &&accessor, bool &value) {
        if (accessor.get())
          value = true;
      });

  return out;
}

template <typename Matrix>
template <typename Row>
void RingmapMatrixRowAccessor<Matrix>::assign_from_row(Row &&row) const
    noexcept(false) {
  ranges::view::zip(row, *this) |
      ranges::for_each([](bool value, auto &&accessor) {
        if (value)
          accessor.set();
        else
          accessor.clear();
      });
}

template <typename Matrix>
ringmap_matrix::base_index_type
RingmapMatrixRowAccessor<Matrix>::begin_index() const noexcept {
  assert(matrix);
  return row->begin_index;
}

template <typename Matrix>
ringmap_matrix::base_index_type
RingmapMatrixRowAccessor<Matrix>::end_index() const noexcept {
  assert(matrix);
  return row->end_index;
}

template <typename Matrix>
void RingmapMatrixRowAccessor<Matrix>::set_begin_index(
    ringmap_matrix::base_index_type value) noexcept {
  row->begin_index = value;
}

template <typename Matrix>
void RingmapMatrixRowAccessor<Matrix>::set_end_index(
    ringmap_matrix::base_index_type value) noexcept {
  row->end_index = value;
}

template <typename ML, typename MR>
bool operator<(RingmapMatrixRowAccessor<ML> const &lhs,
               RingmapMatrixRowAccessor<MR> const &rhs) noexcept {
  return *lhs.row < *rhs.row;
}

template <typename M>
bool operator<(RingmapMatrixRowAccessor<M> const &lhs,
               ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row < rhs;
}

template <typename M>
bool operator<(ringmap_matrix::row_type const &lhs,
               RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs < *rhs.row;
}

template <typename ML, typename MR>
bool operator<=(RingmapMatrixRowAccessor<ML> const &lhs,
                RingmapMatrixRowAccessor<MR> const &rhs) noexcept {
  return *lhs.row <= *rhs.row;
}

template <typename M>
bool operator<=(RingmapMatrixRowAccessor<M> const &lhs,
                ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row <= rhs;
}

template <typename M>
bool operator<=(ringmap_matrix::row_type const &lhs,
                RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs <= *rhs.row;
}

template <typename ML, typename MR>
bool operator>=(RingmapMatrixRowAccessor<ML> const &lhs,
                RingmapMatrixRowAccessor<MR> const &rhs) noexcept {
  return *lhs.row >= *rhs.row;
}

template <typename M>
bool operator>=(RingmapMatrixRowAccessor<M> const &lhs,
                ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row >= rhs;
}

template <typename M>
bool operator>=(ringmap_matrix::row_type const &lhs,
                RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs >= *rhs.row;
}

template <typename ML, typename MR>
bool operator>(RingmapMatrixRowAccessor<ML> const &lhs,
               RingmapMatrixRowAccessor<MR> const &rhs) noexcept {
  return *lhs.row > *rhs.row;
}

template <typename M>
bool operator>(RingmapMatrixRowAccessor<M> const &lhs,
               ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row > rhs;
}

template <typename M>
bool operator>(ringmap_matrix::row_type const &lhs,
               RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs > *rhs.row;
}

template <typename M>
bool operator==(RingmapMatrixRowAccessor<M> const &lhs,
                ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row == rhs;
}

template <typename M>
bool operator==(ringmap_matrix::row_type const &lhs,
                RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs == *rhs.row;
}

template <typename M>
bool operator!=(RingmapMatrixRowAccessor<M> const &lhs,
                ringmap_matrix::row_type const &rhs) noexcept {
  return *lhs.row != rhs;
}

template <typename M>
bool operator!=(ringmap_matrix::row_type const &lhs,
                RingmapMatrixRowAccessor<M> const &rhs) noexcept {
  return lhs != *rhs.row;
}

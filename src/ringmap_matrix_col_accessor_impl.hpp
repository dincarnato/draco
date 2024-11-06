#pragma once

#include "ringmap_matrix_col_accessor.hpp"
#include "ringmap_matrix_iterator.hpp"

#include <algorithm>
#include <numeric>
#include <random>
#include <stdexcept>
#include <utility>

template <typename Matrix>
RingmapMatrixColAccessor<Matrix>::RingmapMatrixColAccessor(
    Matrix &matrix, unsigned col) noexcept
    : matrix(&matrix), col(col) {}

template <typename Matrix>
unsigned RingmapMatrixColAccessor<Matrix>::size() const noexcept {
  return matrix->readsCount;
}

template <typename Matrix>
std::size_t RingmapMatrixColAccessor<Matrix>::sum() const noexcept {
  return std::accumulate(std::ranges::begin(*this), std::ranges::end(*this),
                         0u);
}

template <typename Matrix>
double RingmapMatrixColAccessor<Matrix>::mean() const noexcept {
  return static_cast<double>(sum()) / matrix->readsCount;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::begin() const noexcept -> iterator {
  return {*matrix, matrix->data.data(), col};
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::end() const noexcept -> iterator {
  return {*matrix, matrix->data.data() + matrix->readsCount, col};
}

template <typename Matrix>
template <typename URBG>
void RingmapMatrixColAccessor<Matrix>::shuffle(URBG &&g) const noexcept(false) {
  std::uniform_int_distribution<unsigned> distribution(0,
                                                       matrix->readsCount - 1);
  for (auto &&row : matrix->data) {
    auto colIter = std::ranges::lower_bound(row, col);
    if (colIter == std::ranges::end(row) or *colIter != col)
      continue;

    auto &otherRow = matrix->data[distribution(g)];
    if (std::ranges::binary_search(otherRow, col))
      continue;

    row.erase(colIter);
    otherRow.insert(std::ranges::upper_bound(otherRow, col), col);
  }
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<RingmapMatrix> const &rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<const RingmapMatrix> const &rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<RingmapMatrix &&> const &rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(rhs);
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<RingmapMatrix> &&rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<const RingmapMatrix> &&rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(
    RingmapMatrixColAccessor<RingmapMatrix &&> &&rhs) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  copy_from_accessor(std::move(rhs));
  return *this;
}

template <typename Matrix>
RingmapMatrixAccessor<Matrix>
RingmapMatrixColAccessor<Matrix>::operator[](std::size_t index) const noexcept {
  return {matrix->data.data() + index, col};
}

template <typename Matrix>
template <typename Accessor>
void RingmapMatrixColAccessor<Matrix>::copy_from_accessor(Accessor &&rhs) const
    noexcept(false) {
  if (matrix->data.size() < rhs.matrix->data.size())
    matrix->data.resize(rhs.matrix->data.size());

  if (matrix->readsCount < rhs.matrix->readsCount)
    matrix->readsCount = rhs.matrix->readsCount;

  auto rowIter = std::ranges::begin(matrix->data);
  auto const rowEnd = std::ranges::end(matrix->data);
  auto rhsRowIter = std::ranges::begin(std::as_const(rhs.matrix->data));
  for (; rowIter < rowEnd; ++rowIter, ++rhsRowIter) {
    auto &&rhsRow = *rhsRowIter;
    auto &row = *rowIter;
    if (std::ranges::find(rhsRow, rhs.col) == std::ranges::end(rhsRow)) {
      auto &&[begin, end] = std::ranges::remove(row, col);
      row.erase(begin, end);
    } else if (std::ranges::find(row, col) == std::ranges::end(row)) {
      row.emplace_back(col);
      std::ranges::sort(row);
    }
  }
}

template <typename Matrix>
auto RingmapMatrixColAccessor<Matrix>::operator=(value_type const &column) const
    noexcept(false) -> RingmapMatrixColAccessor const & {
  assert(matrix);
  if (column.size() != matrix->readsCount)
    throw std::runtime_error("invalid column size");

  auto col_iter = std::ranges::begin(column);
  auto this_iter = std::ranges::begin(this);
  for (; col_iter < std::ranges::end(column) and
         this_iter < std::ranges::end(*this);
       ++col_iter, ++this_iter) {
    bool element = *col_iter;
    auto &&accessor = *this_iter;

    if (element)
      accessor.set();
    else
      accessor.clear();
  }

  return *this;
}

template <typename Matrix>
RingmapMatrixColAccessor<Matrix>::operator value_type() const noexcept(false) {
  value_type out(matrix->readsCount, false);

  auto this_iter = std::ranges::begin(this);
  auto out_iter = std::ranges::begin(out);
  for (;
       this_iter < std::ranges::end(*this) and out_iter < std::ranges::end(out);
       ++this_iter, ++out_iter) {
    auto &&accessor = *this_iter;

    if (accessor.get())
      *out_iter = true;
  }

  return out;
}

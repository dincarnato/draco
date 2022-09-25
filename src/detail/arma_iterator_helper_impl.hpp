#pragma once

#include "arma_accessor.hpp"
#include "arma_iterator_helper.hpp"

namespace detail {

template <typename Mat, ArmaIteratorDirection direction>
ArmaIteratorHelper<Mat, direction>::ArmaIteratorHelper(Mat *matrix) noexcept
    : accessor(std::in_place_t{}, matrix, 0) {}

template <typename Mat, ArmaIteratorDirection direction>
ArmaIteratorHelper<Mat, direction>::ArmaIteratorHelper(Mat *matrix,
                                                       int) noexcept
    : accessor([&] {
        if constexpr (direction == ArmaIteratorDirection::rows) {
          return reference(std::in_place_t{}, matrix, matrix->n_rows);
        } else {
          return reference(std::in_place_t{}, matrix, matrix->n_cols);
        }
      }()) {}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator++() noexcept -> self & {
  assert(accessor);
  ++*accessor;
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator++(int) noexcept -> self {
  assert(accessor);
  auto other = *this;
  ++*accessor;
  return other;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator--() noexcept -> self & {
  assert(accessor);
  --*accessor;
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator--(int) noexcept -> self {
  assert(accessor);
  auto other = *this;
  --*accessor;
  return other;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator+=(
    difference_type offset) noexcept -> self & {
  assert(accessor);
  *accessor += offset;
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator+(
    difference_type offset) const noexcept -> self {
  auto other = *this;
  other += offset;
  return other;
}

template <typename Mat, ArmaIteratorDirection direction>
ArmaIteratorHelper<Mat, direction>
operator+(typename ArmaIteratorHelper<Mat, direction>::difference_type offset,
          ArmaIteratorHelper<Mat, direction> const &iter) noexcept {
  return iter + offset;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator-=(
    difference_type offset) noexcept -> self & {
  assert(accessor);
  *accessor -= offset;
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator-(
    difference_type offset) const noexcept -> self {
  auto other = *this;
  other -= offset;
  return other;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator-(
    ArmaIteratorHelper const &rhs) const noexcept -> difference_type {
  assert(accessor);
  assert(rhs.accessor);
  return *accessor - *rhs.accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
auto ArmaIteratorHelper<Mat, direction>::operator*() const noexcept
    -> reference {
  assert(accessor);
  return *accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator==(
    const self &other) const noexcept {
  return accessor == other.accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator!=(
    const self &other) const noexcept {
  return accessor != other.accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator<(
    const self &other) const noexcept {
  assert(accessor);
  assert(other.accessor);
  return *accessor < *other.accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator<=(
    const self &other) const noexcept {
  assert(accessor);
  assert(other.accessor);
  return *accessor <= *other.accessor;
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator>=(
    const self &other) const noexcept {
  return not(*this < other);
}

template <typename Mat, ArmaIteratorDirection direction>
bool ArmaIteratorHelper<Mat, direction>::operator>(
    const self &other) const noexcept {
  return not(*this <= other);
}
} // namespace detail

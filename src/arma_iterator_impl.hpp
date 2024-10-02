#pragma once

#include "arma_iterator.hpp"
#include <concepts>
#include <iterator>

template <typename Mat>
ArmaIterator<Mat>::ArmaIterator(Mat &matrix) noexcept : matrix(&matrix) {}

template <typename Mat>
auto ArmaIterator<Mat>::rows() const noexcept -> rows_handler {
  return rows_handler(matrix);
}

template <typename Mat>
auto ArmaIterator<Mat>::cols() const noexcept -> cols_handler {
  return cols_handler(matrix);
}

static_assert(std::copyable<typename detail::ArmaIteratorHandler<
                  arma::mat, detail::ArmaIteratorDirection::rows>::iterator>);
static_assert(
    std::input_or_output_iterator<typename detail::ArmaIteratorHandler<
        arma::mat, detail::ArmaIteratorDirection::rows>::iterator>);

static_assert(std::random_access_iterator<typename detail::ArmaIteratorHandler<
                  arma::mat, detail::ArmaIteratorDirection::rows>::iterator>);
static_assert(
    std::random_access_iterator<typename detail::ArmaIteratorHandler<
        arma::mat const, detail::ArmaIteratorDirection::rows>::iterator>);

static_assert(std::random_access_iterator<typename detail::ArmaIteratorHandler<
                  arma::mat, detail::ArmaIteratorDirection::cols>::iterator>);
static_assert(
    std::random_access_iterator<typename detail::ArmaIteratorHandler<
        arma::mat const, detail::ArmaIteratorDirection::cols>::iterator>);

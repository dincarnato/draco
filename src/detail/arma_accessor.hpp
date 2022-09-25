#pragma once

#include "arma_iterator_direction.hpp"

#include <armadillo>
#include <type_traits>

namespace detail {

template <typename Mat, ArmaIteratorDirection direction>
class ArmaAccessor
    : public std::conditional_t<
          direction == ArmaIteratorDirection::rows,
          arma::subview_row<typename std::remove_const_t<Mat>::elem_type>,
          arma::subview_col<typename std::remove_const_t<Mat>::elem_type>> {
  using matrix_type = Mat;
  using base_matrix_type = std::remove_const_t<Mat>;
  static constexpr ArmaIteratorDirection direction_value = direction;
  using matrix_value_type = typename base_matrix_type::elem_type;
  using value_type = matrix_value_type;
  using base_type = std::conditional_t<direction == ArmaIteratorDirection::rows,
                                       arma::subview_row<matrix_value_type>,
                                       arma::subview_col<matrix_value_type>>;
  using subview_base_type = arma::subview<matrix_value_type>;
  using concrete_type =
      std::conditional_t<direction == ArmaIteratorDirection::rows, arma::rowvec,
                         arma::colvec>;
  using self = ArmaAccessor;

public:
  using base_type::base_type;
  using base_type::operator();

  ArmaAccessor() = default;
  ArmaAccessor(self const &) = default;
  ArmaAccessor(self &&) = default;
  ArmaAccessor(Mat *matrix, std::size_t index) noexcept;

  self const &operator=(self const &accessor) const noexcept;
  self const &operator=(self &&accessor) const noexcept;

  template <typename _Mat>
  self const &
  operator=(ArmaAccessor<_Mat, direction> const &accessor) const noexcept;

  template <typename _Mat>
  self const &
  operator=(ArmaAccessor<_Mat, direction> &&accessor) const noexcept;

  ArmaAccessor const &operator=(concrete_type const &concrete) const noexcept;
  ArmaAccessor const &operator=(concrete_type &&concrete) const noexcept;

  std::size_t size() const noexcept;

private:
  Mat *matrix = nullptr;
};

} // namespace detail

#include "arma_accessor_impl.hpp"

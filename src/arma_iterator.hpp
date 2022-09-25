#pragma once

#include "detail/arma_accessor.hpp"
#include "detail/arma_iterator_direction.hpp"
#include "detail/arma_iterator_handler.hpp"
#include "detail/arma_iterator_helper.hpp"

#include <armadillo>
#include <memory>
#include <type_traits>

namespace arma {

template <typename Mat>
struct is_subview_row<
    detail::ArmaAccessor<Mat, detail::ArmaIteratorDirection::rows>>
    : std::enable_if_t<is_arma_type<Mat>::value, std::true_type> {};

template <typename Mat>
struct is_subview_col<
    detail::ArmaAccessor<Mat, detail::ArmaIteratorDirection::cols>>
    : std::enable_if_t<is_arma_type<Mat>::value, std::true_type> {};

template <typename Mat>
struct Proxy<detail::ArmaAccessor<Mat, detail::ArmaIteratorDirection::rows>>
    : std::enable_if_t<is_arma_type<Mat>::value,
                       Proxy<subview_row<typename Mat::elem_type>>> {
  using Proxy<subview_row<typename Mat::elem_type>>::Proxy;
};

template <typename Mat>
struct Proxy<detail::ArmaAccessor<Mat, detail::ArmaIteratorDirection::cols>>
    : std::enable_if_t<is_arma_type<Mat>::value,
                       Proxy<subview_col<typename Mat::elem_type>>> {
  using Proxy<subview_col<typename Mat::elem_type>>::Proxy;
};

} /* namespace arma */

template <typename Mat> class ArmaIterator {
public:
  static_assert(not std::is_reference_v<Mat>);
  using matrix_type = Mat;
  using reference = Mat &;
  using pointer = Mat *;

  using rows_handler =
      detail::ArmaIteratorHandler<Mat, detail::ArmaIteratorDirection::rows>;
  using cols_handler =
      detail::ArmaIteratorHandler<Mat, detail::ArmaIteratorDirection::cols>;

  ArmaIterator(Mat &matrix) noexcept;

  rows_handler rows() const noexcept;
  cols_handler cols() const noexcept;

private:
  pointer matrix;
};

#include "arma_iterator_impl.hpp"

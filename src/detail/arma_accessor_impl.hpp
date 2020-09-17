#pragma once

#include "arma_accessor.hpp"

namespace detail {

template <typename Mat, ArmaIteratorDirection direction>
ArmaAccessor<Mat, direction>::ArmaAccessor(Mat* matrix,
                                           std::size_t index) noexcept
    : base_type([&] {
        if constexpr (direction == ArmaIteratorDirection::rows) {
          assert(index < matrix->n_rows);
          return matrix->row(index);
        } else {
          assert(index < matrix->n_cols);
          return matrix->col(index);
        }
      }()),
      matrix(matrix) {}

template <typename Mat, ArmaIteratorDirection direction>
std::size_t
ArmaAccessor<Mat, direction>::size() const noexcept {
  assert(matrix);
  if constexpr (direction == ArmaIteratorDirection::rows)
    return matrix->n_cols;
  else
    return matrix->n_rows;
}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaAccessor<Mat, direction>::operator=(self const& accessor) const noexcept
    -> self const& {
  assert(matrix);
  assert(accessor.matrix);

  base_type::operator=(static_cast<base_type const&>(accessor));
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaAccessor<Mat, direction>::operator=(self&& accessor) const noexcept
    -> self const& {
  assert(matrix);
  assert(accessor.matrix);

  base_type::operator=(static_cast<base_type&&>(std::move(accessor)));
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
template <typename _Mat>
auto
ArmaAccessor<Mat, direction>::
operator=(ArmaAccessor<_Mat, direction> const& accessor) const noexcept
    -> self const& {
  using rhs_mat = std::decay_t<_Mat>;
  static_assert(arma::is_Mat<rhs_mat>::value);
  assert(matrix);
  assert(accessor.matrix);

  base_type::operator=(
      static_cast<typename rhs_mat::base_type const&>(accessor));
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
template <typename _Mat>
auto
ArmaAccessor<Mat, direction>::
operator=(ArmaAccessor<_Mat, direction>&& accessor) const noexcept
    -> self const& {
  using rhs_mat = std::decay_t<_Mat>;
  static_assert(arma::is_Mat<rhs_mat>::value);
  assert(matrix);
  assert(accessor.matrix);

  base_type::operator=(
      static_cast<typename rhs_mat::base_type&&>(std::move(accessor)));
  return *this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaAccessor<Mat, direction>::operator=(concrete_type const& concrete) const
    noexcept -> self const& {
  assert(matrix);
  base_type::operator=(concrete);
  return this;
}

template <typename Mat, ArmaIteratorDirection direction>
auto
ArmaAccessor<Mat, direction>::operator=(concrete_type&& concrete) const noexcept
    -> self const& {
  assert(matrix);
  base_type::operator=(std::move(concrete));
  return this;
}

} // namespace detail

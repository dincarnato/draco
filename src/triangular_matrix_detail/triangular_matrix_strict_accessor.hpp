#pragma once

#include "triangular_matrix_strict_reference_wrapper.hpp"

#include <type_traits>

template <typename, typename> struct TriangularMatrixStrict;

namespace detail {

template <typename TrMat> struct TriangularMatrixStrictAccessor {
  using triangular_matrix_type = std::decay_t<TrMat>;
  using value_type = typename triangular_matrix_type::value_type;
  using size_type = typename triangular_matrix_type::size_type;
  using matrix_pointer = TrMat *;
  using reference =
      std::conditional_t<std::is_const_v<TrMat>, const value_type &,
                         TriangularMatrixStrictReferenceWrapper<TrMat>>;

  TriangularMatrixStrictAccessor(TrMat &tr_mat, size_type index) noexcept
      : tr_mat(&tr_mat), _index(index) {}

  reference operator[](size_type index) const noexcept {
    if constexpr (std::is_const_v<TrMat>) {
      constexpr static value_type zero(0);

      if (_index == index)
        return zero;

      const auto data_index = [&] {
        if (_index < index)
          return tr_mat->get_offset(_index, index);
        else
          return tr_mat->get_offset(index, _index);
      }();
      return tr_mat->_data[data_index];
    } else {
      return TriangularMatrixStrictReferenceWrapper(*tr_mat, _index, index);
    }
  }

private:
  matrix_pointer tr_mat;
  size_type _index;
};

} // namespace detail

#pragma once

#include "arma_iterator_direction.hpp"

#include <armadillo>
#include <optional>
#include <range/v3/core.hpp>

namespace detail {

template <typename, ArmaIteratorDirection>
class ArmaAccessor;

template <typename Mat, ArmaIteratorDirection direction>
class ArmaIteratorHelper {
  template <typename, ArmaIteratorDirection>
  friend class ArmaIteratorHandler;

public:
  using value_type =
      std::conditional_t<direction == ArmaIteratorDirection::rows, arma::rowvec,
                         arma::colvec>;
  using reference = ArmaAccessor<Mat, direction>;
  using difference_type = std::ptrdiff_t;
  using iterator_category = ranges::random_access_iterator_tag;
  using self = ArmaIteratorHelper;

  ArmaIteratorHelper() = default;
  ArmaIteratorHelper(Mat* matrix) noexcept;
  ArmaIteratorHelper(Mat* matrix, int) noexcept;

  self& operator++() noexcept;
  self operator++(int) noexcept;
  self& operator--() noexcept;
  self operator--(int) noexcept;
  self& operator+=(difference_type offset) noexcept;
  self operator+(difference_type offset) const noexcept;

  template <typename _Mat, ArmaIteratorDirection _direction>
  friend ArmaIteratorHelper<_Mat, _direction> operator+(
      typename ArmaIteratorHelper<_Mat, _direction>::difference_type offset,
      ArmaIteratorHelper<_Mat, _direction> const& iter) noexcept;

  self& operator-=(difference_type offset) noexcept;
  self operator-(difference_type offset) const noexcept;
  difference_type operator-(self const& rhs) const noexcept;

  bool operator==(const self& other) const noexcept;
  bool operator!=(const self& other) const noexcept;
  bool operator<(const self& other) const noexcept;
  bool operator<=(const self& other) const noexcept;
  bool operator>=(const self& other) const noexcept;
  bool operator>(const self& other) const noexcept;

  reference operator*() const noexcept;
  reference operator[](difference_type offset) const noexcept;

private:
  std::optional<reference> accessor;
};

} // namespace detail

#include "arma_iterator_helper_impl.hpp"

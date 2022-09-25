#pragma once

#include "arma_iterator_direction.hpp"

#include <type_traits>

namespace detail {

template <typename, ArmaIteratorDirection> class ArmaIteratorHelper;

template <typename Mat, ArmaIteratorDirection direction>
class ArmaIteratorHandler {
public:
  using iterator = ArmaIteratorHelper<Mat, direction>;
  using value_type = typename iterator::value_type;
  using reference = typename iterator::reference;

  explicit ArmaIteratorHandler(Mat *matrix) noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;

  std::size_t size() const noexcept;

private:
  Mat *matrix;
};

} // namespace detail

#include "arma_iterator_handler_impl.hpp"

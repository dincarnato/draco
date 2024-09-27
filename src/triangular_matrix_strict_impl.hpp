#pragma once

#include "triangular_matrix_strict.hpp"

#include <cmath>

template <typename T, typename Alloc>
TriangularMatrixStrict<T, Alloc>::TriangularMatrixStrict(
    const Alloc &alloc) noexcept(noexcept(Alloc()))
    : _n_elements(0), _data(alloc) {}

template <typename T, typename Alloc>
TriangularMatrixStrict<T, Alloc>::TriangularMatrixStrict(size_type elements,
                                                         const Alloc &alloc)
    : _n_elements(elements), _data(get_size_from_elements(elements), alloc) {}

template <typename T, typename Alloc>
TriangularMatrixStrict<T, Alloc>::TriangularMatrixStrict(size_type elements,
                                                         const T &value,
                                                         const Alloc &alloc)
    : _n_elements(elements),
      _data(get_size_from_elements(elements), value, alloc) {}

template <typename T, typename Alloc>
constexpr inline auto TriangularMatrixStrict<T, Alloc>::get_size_from_elements(
    size_type elements) noexcept -> size_type {
  if (elements < 2)
    return 0;
  else
    return (elements - 1) * elements / 2;
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::get_offset(
    size_type row, size_type col) const noexcept -> size_type {
  return get_offset_with_elements(row, col, _n_elements);
}

template <typename T, typename Alloc>
constexpr inline auto
TriangularMatrixStrict<T, Alloc>::get_offset_with_elements(
    size_type row, size_type col, size_type n_elements) noexcept -> size_type {
  assert(row < col);
  return (col - row - 1) + row * (n_elements * 2 - row - 1) / 2;
}

template <typename T, typename Alloc>
TriangularMatrixStrict<T, Alloc>::TriangularMatrixStrict(
    std::initializer_list<std::initializer_list<T>> data, const Alloc &alloc)
    : _n_elements(std::size(data)),
      _data(get_size_from_elements(_n_elements), alloc) {

  const auto elements_size = std::size(data);
  auto write_iter = std::begin(_data);

  size_type row_index = 0;
  const auto end_data = std::end(data);
  for (auto row_iter = std::begin(data); row_iter < end_data;
       ++row_iter, ++row_index) {
    const auto &row = *row_iter;

    if (std::size(row) != elements_size)
      throw std::invalid_argument(
          "expected a triangular-matrix-like as argument; the number of "
          "elements per row is inconsistent");

    size_type col_index = 0;
    const auto end_row = std::end(row);
    for (auto col_iter = std::begin(row); col_iter < end_row;
         ++col_iter, ++col_index) {
      if (*col_iter != std::begin(std::begin(data)[col_index])[row_index])
        throw std::invalid_argument("expected a triangular-matrix-like as "
                                    "argument; the matrix is not symmetric");
    }

    if (std::begin(row)[row_index] != T(0))
      throw std::invalid_argument("expected a triangular-matrix-like as "
                                  "argument; diagonal must be zero");

    for (auto col_iter = std::next(
             std::begin(row),
             static_cast<
                 typename std::iterator_traits<iterator>::difference_type>(
                 row_index + size_type(1)));
         col_iter < end_row; ++col_iter) {
      *write_iter++ = *col_iter;
    }
  }
}

template <typename T, typename Alloc>
inline void
TriangularMatrixStrict<T, Alloc>::resize(size_type new_size) noexcept(false) {
  const auto new_storage_size = get_size_from_elements(new_size);
  storage_type new_storage(new_storage_size);

  const auto max_elements = std::min(_n_elements, new_size);

  for (size_type row = 0; row < max_elements; ++row) {
    for (size_type col = row + 1; col < max_elements; ++col) {
      auto new_storage_offset = get_offset_with_elements(row, col, new_size);
      new_storage[new_storage_offset] = _data[get_offset(row, col)];
    }
  }

  _data = std::move(new_storage);
  _n_elements = new_size;
}

template <typename T, typename Alloc>
inline void TriangularMatrixStrict<T, Alloc>::resize(
    size_type new_size, value_type const &value) noexcept(false) {
  const auto new_storage_size = get_size_from_elements(new_size);
  storage_type new_storage(new_storage_size);

  const auto max_elements = std::min(_n_elements, new_size);

  for (size_type row = 0; row < max_elements; ++row) {
    for (size_type col = row + 1; col < max_elements; ++col) {
      auto new_storage_offset = get_offset_with_elements(row, col, new_size);
      new_storage[new_storage_offset] = _data[get_offset(row, col)];
    }

    for (size_type col = max_elements; col < new_size; ++col) {
      auto new_storage_offset = get_offset_with_elements(row, col, new_size);
      new_storage[new_storage_offset] = value;
    }
  }

  for (size_type row = max_elements; row < new_size; ++row) {
    for (size_type col = row + 1; col < new_size; ++col) {
      auto new_storage_offset = get_offset_with_elements(row, col, new_size);
      new_storage[new_storage_offset] = value;
    }
  }

  _data = std::move(new_storage);
  _n_elements = new_size;
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::begin() noexcept -> iterator {
  return _data.begin();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::begin() const noexcept -> const_iterator {
  return _data.begin();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::cbegin() const noexcept -> const_iterator {
  return _data.cbegin();
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::end() noexcept -> iterator {
  return _data.end();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::end() const noexcept -> const_iterator {
  return _data.end();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::cend() const noexcept -> const_iterator {
  return _data.cend();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::rbegin() noexcept -> reverse_iterator {
  return _data.rbegin();
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::rbegin() const noexcept
    -> const_reverse_iterator {
  return _data.rbegin();
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::crbegin() const noexcept
    -> const_reverse_iterator {
  return _data.crbegin();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::rend() noexcept -> reverse_iterator {
  return _data.rend();
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::rend() const noexcept
    -> const_reverse_iterator {
  return _data.rend();
}

template <typename T, typename Alloc>
inline auto TriangularMatrixStrict<T, Alloc>::crend() const noexcept
    -> const_reverse_iterator {
  return _data.crend();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::size() const noexcept -> size_type {
  return _data.size();
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::rows_size() const noexcept -> size_type {
  return _n_elements;
}

template <typename T, typename Alloc>
inline auto
TriangularMatrixStrict<T, Alloc>::cols_size() const noexcept -> size_type {
  return _n_elements;
}

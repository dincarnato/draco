#pragma once

#include <memory>
#include <vector>

#include "triangular_matrix_detail/triangular_matrix_strict_accessor.hpp"
#include "triangular_matrix_detail/triangular_matrix_strict_reference_wrapper.hpp"

template <typename T, typename Alloc = std::allocator<T>>
struct TriangularMatrixStrict {
  using value_type = T;
  using allocator = Alloc;
  using storage_type = std::vector<T, Alloc>;
  using size_type = typename storage_type::size_type;
  using iterator = typename storage_type::iterator;
  using const_iterator = typename storage_type::const_iterator;
  using reverse_iterator = typename storage_type::reverse_iterator;
  using const_reverse_iterator = typename storage_type::const_reverse_iterator;
  using accessor =
      detail::TriangularMatrixStrictAccessor<TriangularMatrixStrict>;
  using const_accessor =
      detail::TriangularMatrixStrictAccessor<const TriangularMatrixStrict>;

  template <typename>
  friend struct detail::TriangularMatrixStrictAccessor;
  template <typename>
  friend struct detail::TriangularMatrixStrictReferenceWrapper;

  TriangularMatrixStrict() = default;
  explicit TriangularMatrixStrict(const Alloc& alloc) noexcept(
      noexcept(Alloc()));
  explicit TriangularMatrixStrict(size_type elements,
                                  const Alloc& alloc = Alloc());
  explicit TriangularMatrixStrict(size_type elements, const T& value,
                                  const Alloc& alloc = Alloc());
  explicit TriangularMatrixStrict(
      std::initializer_list<std::initializer_list<T>> data,
      const Alloc& alloc = Alloc());

  TriangularMatrixStrict(const TriangularMatrixStrict&) = default;
  TriangularMatrixStrict(TriangularMatrixStrict&&) = default;
  TriangularMatrixStrict(const TriangularMatrixStrict& other,
                         const Alloc& alloc);
  TriangularMatrixStrict(TriangularMatrixStrict&& other, const Alloc& alloc);

  TriangularMatrixStrict& operator=(const TriangularMatrixStrict&) = default;
  TriangularMatrixStrict& operator=(TriangularMatrixStrict&&) = default;

  inline void resize(size_type new_size) noexcept(false);
  inline void resize(size_type new_size,
                     value_type const& value) noexcept(false);

  inline iterator begin() noexcept;
  inline const_iterator begin() const noexcept;
  inline const_iterator cbegin() const noexcept;
  inline iterator end() noexcept;
  inline const_iterator end() const noexcept;
  inline const_iterator cend() const noexcept;

  inline reverse_iterator rbegin() noexcept;
  inline const_reverse_iterator rbegin() const noexcept;
  inline const_reverse_iterator crbegin() const noexcept;
  inline reverse_iterator rend() noexcept;
  inline const_reverse_iterator rend() const noexcept;
  inline const_reverse_iterator crend() const noexcept;

  inline size_type size() const noexcept;
  inline size_type rows_size() const noexcept;
  inline size_type cols_size() const noexcept;

  accessor operator[](size_type index) noexcept {
    return accessor(*this, index);
  }

  const_accessor operator[](size_type index) const noexcept {
    return const_accessor(*this, index);
  }

private:
  size_type _n_elements{0};
  storage_type _data;

  static constexpr inline size_type
  get_size_from_elements(size_type elements) noexcept;

  static constexpr inline size_type
  get_offset_with_elements(size_type row, size_type col,
                           size_type n_elements) noexcept;

  inline size_type get_offset(size_type row, size_type col) const noexcept;
};

#include "triangular_matrix_strict_impl.hpp"

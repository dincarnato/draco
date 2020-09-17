#pragma once

#include <type_traits>

template <typename, typename>
struct TriangularMatrixStrict;

namespace detail {

template <typename TrMat>
struct TriangularMatrixStrictReferenceWrapper {
  static_assert(not std::is_const_v<TrMat>);

  using triangular_matrix_type = TrMat;
  using value_type = typename triangular_matrix_type::value_type;
  using size_type = typename triangular_matrix_type::size_type;
  using matrix_pointer = TrMat*;
  using reference = value_type&;
  using pointer = value_type*;

  TriangularMatrixStrictReferenceWrapper(TrMat& tr_mat, size_type indexA,
                                         size_type indexB) noexcept
      : value([&]() -> pointer {
          if (indexA == indexB)
            return nullptr;
          else {
            const auto index = [&] {
              if (indexA < indexB)
                return tr_mat.get_offset(indexA, indexB);
              else
                return tr_mat.get_offset(indexB, indexA);
            }();

            return &tr_mat._data[index];
          }
        }()) {}

  inline void
  check_validity() const noexcept(
#ifndef NDEBUG
      false
#else
      true
#endif
  ) {
#ifndef NDEBUG
    if (value == nullptr)
      throw std::runtime_error(
          "cannot set the diagonal of a strictly triangular matrix");
#endif
  }

  inline operator value_type&() noexcept(noexcept(check_validity())) {
    return get();
  }

  inline value_type&
  get() noexcept(noexcept(check_validity())) {
    check_validity();
    return *value;
  }

  inline const value_type&
  get_const() const noexcept {
    if (value == nullptr)
      return zero;
    else
      return *value;
  }

  template <typename T>
  inline value_type&
  operator=(T&& rhs) noexcept(noexcept(check_validity()) and
                              noexcept(*this->value = std::forward<T>(rhs))) {
    check_validity();
    return *value = std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator+=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value += std::forward<T>(rhs))) {
    check_validity();
    return *value += std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator-=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value -= std::forward<T>(rhs))) {
    check_validity();
    return *value -= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator*=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value *= std::forward<T>(rhs))) {
    check_validity();
    return *value *= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator/=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value /= std::forward<T>(rhs))) {
    check_validity();
    return *value /= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator%=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value %= std::forward<T>(rhs))) {
    check_validity();
    return *value %= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator^=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value ^= std::forward<T>(rhs))) {
    check_validity();
    return *value ^= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator&=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value &= std::forward<T>(rhs))) {
    check_validity();
    return *value &= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator|=(T&& rhs) noexcept(noexcept(check_validity()) and
                               noexcept(*this->value |= std::forward<T>(rhs))) {
    check_validity();
    return *value |= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator<<=(T&& rhs) noexcept(noexcept(check_validity()) and
                                noexcept(*this->value <<=
                                         std::forward<T>(rhs))) {
    check_validity();
    return *value <<= std::forward<T>(rhs);
  }

  template <typename T>
  inline value_type&
  operator>>=(T&& rhs) noexcept(noexcept(check_validity()) and
                                noexcept(*this->value >>=
                                         std::forward<T>(rhs))) {
    check_validity();
    return *value >>= std::forward<T>(rhs);
  }

  inline value_type&
  operator++() noexcept(noexcept(check_validity()) and
                        noexcept(++*this->value)) {
    check_validity();
    return ++*value;
  }

  inline value_type
  operator++(int) noexcept(noexcept(check_validity()) and noexcept((*this->value)++)) {
    check_validity();
    return (*value)++;
  }

  inline value_type&
  operator--() noexcept(noexcept(check_validity()) and noexcept(--*this->value)) {
    check_validity();
    return --*value;
  }

  inline value_type
  operator--(int) noexcept(noexcept(check_validity()) and noexcept((*this->value)--)) {
    check_validity();
    return (*value)--;
  }

private: 
  pointer value;
  static constexpr value_type zero{0};
};

} // namespace detail

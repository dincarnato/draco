#pragma once

#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace detail {

template <typename T> struct TinyFractionBase {
  template <typename> friend struct TinyFractionBase;

  using value_type = T;
  using higher_value_type = std::conditional_t<
      std::is_signed_v<T>,
      std::conditional_t<std::is_same_v<T, std::int8_t>, std::int16_t,
                         std::conditional_t<std::is_same_v<T, std::int16_t>,
                                            std::int32_t, std::int64_t>>,
      std::conditional_t<std::is_same_v<T, std::uint8_t>, std::uint16_t,
                         std::conditional_t<std::is_same_v<T, std::uint16_t>,
                                            std::uint32_t, std::uint64_t>>>;
  using lower_value_type = std::conditional_t<
      std::is_signed_v<T>,
      std::conditional_t<std::is_same_v<T, std::int64_t>, std::int32_t,
                         std::conditional_t<std::is_same_v<T, std::int32_t>,
                                            std::int16_t, std::int8_t>>,
      std::conditional_t<std::is_same_v<T, std::uint64_t>, std::uint32_t,
                         std::conditional_t<std::is_same_v<T, std::uint32_t>,
                                            std::uint16_t, std::uint8_t>>>;
  static constexpr auto normalizer_value =
      std::is_signed_v<value_type> ? std::numeric_limits<std::int8_t>::max()
                                   : std::numeric_limits<std::uint8_t>::max();
  static constexpr auto max_value = std::numeric_limits<value_type>::max();
  static constexpr auto min_value = std::numeric_limits<value_type>::lowest();

  constexpr TinyFractionBase() = default;
  constexpr TinyFractionBase(const TinyFractionBase &) = default;
  constexpr TinyFractionBase(TinyFractionBase &&) = default;

  constexpr explicit TinyFractionBase(float value)
      : _value([&] {
          if (value < 0.f or value > 1.f)
            throw std::out_of_range("value must be between 0 and 1");

          return static_cast<value_type>(round(value * normalizer_value));
        }()) {}
  constexpr explicit TinyFractionBase(double value)
      : _value([&] {
          if (value < 0. or value > 1.)
            throw std::out_of_range("value must be between 0 and 1");

          return static_cast<value_type>(round(value * normalizer_value));
        }()) {}
  constexpr TinyFractionBase(TinyFractionBase<higher_value_type> higher)
      : _value([&] {
          if (higher._value < min_value or higher._value > max_value)
            throw std::out_of_range(
                "value cannot be represented from higher kind");

          return static_cast<value_type>(higher._value);
        }()) {}

  constexpr TinyFractionBase &operator=(const TinyFractionBase &) = default;
  constexpr TinyFractionBase &operator=(TinyFractionBase &&) = default;

  constexpr TinyFractionBase &
  operator=(TinyFractionBase<higher_value_type> higher) {
    if (higher._value < min_value or higher._value > max_value)
      throw std::out_of_range("value cannot be represented from higher kind");

    _value = static_cast<value_type>(higher._value);
    return *this;
  }

  constexpr TinyFractionBase &operator=(float value) {
    if (value < 0.f or value > 1.f)
      throw std::out_of_range("value must be between 0 and 1");

    _value = static_cast<value_type>(value * normalizer_value);
    return *this;
  }

  constexpr TinyFractionBase &operator=(double value) {
    if (value < 0. or value > 1.)
      throw std::out_of_range("value must be between 0 and 1");

    _value = static_cast<value_type>(value * normalizer_value);
    return *this;
  }

  constexpr operator float() const noexcept {
    return static_cast<float>(_value) / normalizer_value;
  }

  constexpr operator double() const noexcept {
    return static_cast<double>(_value) / normalizer_value;
  }

  constexpr auto operator+(const TinyFractionBase &other) const noexcept {
    if constexpr (std::is_same_v<value_type, std::uint8_t> or
                  std::is_same_v<value_type, std::int8_t>)
      return TinyFractionBase<higher_value_type>{static_cast<higher_value_type>(
          static_cast<higher_value_type>(_value) +
          static_cast<higher_value_type>(other._value))};
    else {
      assert(max_value - other._value <= _value);
      return TinyFractionBase{static_cast<value_type>(_value + other._value)};
    }
  }

  template <typename U>
  constexpr std::enable_if_t<
      std::is_same_v<U, lower_value_type> and
          not std::is_same_v<value_type, std::uint8_t> and
          not std::is_same_v<value_type, std::int8_t>,
      TinyFractionBase>
  operator+(const TinyFractionBase<U> &other) const noexcept {
    assert(max_value - other._value >= _value);
    return TinyFractionBase{static_cast<value_type>(
        _value + static_cast<value_type>(other._value))};
  }

  constexpr auto operator-(const TinyFractionBase &other) const noexcept {
    assert(min_value + other._value <= _value);
    return TinyFractionBase{static_cast<value_type>(_value - other._value)};
  }

  constexpr auto operator*(const TinyFractionBase &other) const noexcept {
    return TinyFractionBase<higher_value_type>{static_cast<higher_value_type>(
        static_cast<higher_value_type>(_value) *
        static_cast<higher_value_type>(other._value) / normalizer_value)};
  }

  template <typename U>
  constexpr std::enable_if_t<std::is_arithmetic_v<std::decay_t<U>> and
                                 not std::is_signed_v<std::decay_t<U>>,
                             TinyFractionBase<higher_value_type>>
  operator*(U rhs) const noexcept {
    return TinyFractionBase<higher_value_type>{static_cast<higher_value_type>(
        static_cast<higher_value_type>(_value) * rhs)};
  }

  template <typename U>
  constexpr std::enable_if_t<
      std::is_arithmetic_v<std::decay_t<U>> and
          std::is_signed_v<std::decay_t<U>>,
      std::conditional_t<
          std::is_signed_v<value_type>, TinyFractionBase,
          TinyFractionBase<std::make_signed_t<higher_value_type>>>>
  operator*(U rhs) const noexcept {
    if constexpr (std::is_signed_v<value_type>)
      return TinyFractionBase{static_cast<value_type>(_value * rhs)};
    else {
      using signed_type = std::make_signed_t<higher_value_type>;
      return TinyFractionBase<signed_type>{
          static_cast<signed_type>(static_cast<signed_type>(_value) / 2 * rhs)};
    }
  }

  template <typename U>
  constexpr std::enable_if_t<std::is_arithmetic_v<std::decay_t<U>> and
                                 not std::is_signed_v<std::decay_t<U>>,
                             TinyFractionBase>
  operator/(U rhs) const noexcept {
    assert(rhs != 0);
    return TinyFractionBase{static_cast<value_type>(_value / rhs)};
  }

  template <typename U>
  constexpr std::enable_if_t<
      std::is_arithmetic_v<std::decay_t<U>> and
          std::is_signed_v<std::decay_t<U>>,
      std::conditional_t<std::is_signed_v<value_type>, TinyFractionBase,
                         TinyFractionBase<std::make_signed_t<value_type>>>>
  operator/(U rhs) const noexcept {
    assert(rhs != 0);
    if constexpr (std::is_signed_v<value_type>)
      return TinyFractionBase{static_cast<value_type>(_value / rhs)};
    else {
      using signed_type = std::make_signed_t<value_type>;
      return TinyFractionBase<signed_type>{
          static_cast<signed_type>(static_cast<signed_type>(_value / 2) / rhs)};
    }
  }

  constexpr auto upcast() const noexcept {
    static_assert(not std::is_same_v<value_type, std::uint64_t> and
                      not std::is_same_v<value_type, std::int64_t>,
                  "cannot explicitly upcast type anymore");
    return TinyFractionBase<higher_value_type>{
        static_cast<higher_value_type>(_value)};
  }

  constexpr auto downcast() const noexcept {
    static_assert(not std::is_same_v<value_type, std::uint8_t> and
                      not std::is_same_v<value_type, std::int8_t>,
                  "cannot explicitly downcast type anymore");
    return TinyFractionBase<lower_value_type>{
        static_cast<lower_value_type>(_value)};
  }

  constexpr bool operator==(const TinyFractionBase &rhs) const noexcept {
    return _value == rhs._value;
  }

  constexpr bool operator!=(const TinyFractionBase &rhs) const noexcept {
    return _value != rhs._value;
  }

  constexpr bool operator<(const TinyFractionBase &rhs) const noexcept {
    return _value < rhs._value;
  }

  constexpr bool operator<=(const TinyFractionBase &rhs) const noexcept {
    return _value <= rhs._value;
  }

  constexpr bool operator>=(const TinyFractionBase &rhs) const noexcept {
    return _value >= rhs._value;
  }

  constexpr bool operator>(const TinyFractionBase &rhs) const noexcept {
    return _value > rhs._value;
  }

protected:
  value_type _value;

  constexpr explicit TinyFractionBase(value_type value) noexcept
      : _value(value) {}

  template <typename U> constexpr static std::decay_t<U> round(U u) noexcept {
    using u_type = std::decay_t<U>;
    static_assert(std::is_floating_point_v<u_type>);
    auto const trunc = static_cast<u_type>(static_cast<long>(u));
    auto const inc = static_cast<u_type>(static_cast<long>(u + 0.5));
    assert((u < 0.5 and u >= -0.5) or
           inc != static_cast<u_type>(static_cast<long>(u - 0.5)));
    if (inc > trunc)
      return inc;
    else
      return trunc;
  }
};

} // namespace detail

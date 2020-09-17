#pragma once

#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace cte {

namespace detail {

template <typename>
struct FloatReprProperties;

template <>
struct FloatReprProperties<std::uint16_t> {
  constexpr static std::size_t exponent_bits = 5;
  constexpr static std::size_t mantissa_bits = 10;
};

template <>
struct FloatReprProperties<std::uint32_t> {
  constexpr static std::size_t exponent_bits = 8;
  constexpr static std::size_t mantissa_bits = 23;
};

template <>
struct FloatReprProperties<std::uint64_t> {
  constexpr static std::size_t exponent_bits = 11;
  constexpr static std::size_t mantissa_bits = 52;
};

} // namespace detail

template <typename T>
struct FloatingPoint {
  static_assert(std::is_floating_point_v<T>,
                "T must be a valid floating point type");
  static_assert(
      std::numeric_limits<T>::is_iec559,
      "Sorry, your floating point type is not conforming to IEEE 754 standard "
      "representation, which is the only one supported in this library");

  using repr_type = std::conditional_t<
      sizeof(T) == 2, std::uint16_t,
      std::conditional_t<
          sizeof(T) == 4, std::uint32_t,
          std::conditional_t<sizeof(T) == 8, std::uint64_t, void>>>;
  using signed_repr_type = std::make_signed_t<repr_type>;

  static_assert(not std::is_same_v<repr_type, void>,
                "Unsupported floating type (too big or too small?)");

  constexpr static std::size_t exponent_bits =
      detail::FloatReprProperties<repr_type>::exponent_bits;
  constexpr static std::size_t mantissa_bits =
      detail::FloatReprProperties<repr_type>::mantissa_bits;
  constexpr static repr_type exponent_bias =
      (repr_type(1) << (exponent_bits - 1)) - 1;
  constexpr static repr_type denorm_exp = ((repr_type(1) << exponent_bits) - 1);
  constexpr static repr_type infinity = denorm_exp << mantissa_bits;
  constexpr static repr_type nan = infinity | repr_type(1);

  constexpr FloatingPoint() noexcept : repr(get_representation(T(0))) {}
  constexpr FloatingPoint(T value) noexcept : repr(get_representation(value)) {}

private:
  explicit constexpr FloatingPoint(repr_type value) noexcept : repr(value) {}

public:
  constexpr static FloatingPoint
  from_representation(repr_type value) noexcept {
    return FloatingPoint(value);
  }

  constexpr bool
  is_positive() const noexcept {
    return (repr & (repr_type(1) << (exponent_bits + mantissa_bits))) == 0;
  }

  constexpr repr_type
  mantissa() const noexcept {
    return repr & ((repr_type(1) << mantissa_bits) - 1);
  }

  constexpr repr_type
  exponent() const noexcept {
    return (repr >> mantissa_bits) & ((repr_type(1) << exponent_bits) - 1);
  }

  constexpr signed_repr_type
  exponent_unbiased() const noexcept {
    return static_cast<signed_repr_type>(exponent()) -
           static_cast<signed_repr_type>(exponent_bias);
  }

  constexpr repr_type
  representation() const noexcept {
    return repr;
  }
  constexpr operator repr_type() const noexcept { return repr; }

  constexpr operator T() const noexcept {
    auto const raw_exponent = this->exponent();
    if (raw_exponent == 0) {
      if (is_positive())
        return T(0.);
      else
        return T(-0.);
    } else if (raw_exponent == denorm_exp) {
      auto mantissa = this->mantissa();
      if (mantissa == 0) {
        if (is_positive())
          return std::numeric_limits<T>::infinity();
        else
          return -std::numeric_limits<T>::infinity();
      } else {
        if (is_positive())
          return std::numeric_limits<T>::quiet_NaN();
        else
          return -std::numeric_limits<T>::quiet_NaN();
      }
    }

    T out = 1;

    {
      auto mantissa = this->mantissa();
      T cur = 0.5;
      for (repr_type mantissa_bit = 0;
           mantissa_bit < mantissa_bits and mantissa != 0; ++mantissa_bit) {
        if (mantissa & (repr_type(1) << (mantissa_bits - 1))) {
          out += cur;
        }
        mantissa <<= 1;
        cur /= T(2);
      }
    }

    auto exponent = exponent_unbiased();
    if (exponent >= 0) {
      for (; exponent > 0; --exponent) {
        out *= T(2);
      }
    } else {
      for (; exponent < 0; ++exponent) {
        out /= T(2);
      }
    }

    if (not is_positive())
      out = -out;

    return out;
  }

  constexpr bool
  is_inf() const noexcept {
    return (repr & ((repr_type(1) << (mantissa_bits + exponent_bits)) - 1)) ==
           infinity;
  }
  constexpr bool
  is_nan() const noexcept {
    return (repr & ((repr_type(1) << (mantissa_bits + exponent_bits)) - 1)) ==
           nan;
  }

  constexpr FloatingPoint
  truncate() const noexcept {
    if (is_inf() or is_nan())
      return *this;

    auto exp = exponent_unbiased();
    auto const mantissa = [&] {
      if (exp < 0) {
        exp = -static_cast<signed_repr_type>(exponent_bias);
        return repr_type(0);
      } else if (exp >= static_cast<signed_repr_type>(mantissa_bits)) {
        return this->mantissa();
      } else {
        return this->mantissa() &
               ((~repr_type(0)) ^
                ((repr_type(1)
                  << (mantissa_bits - static_cast<std::size_t>(exp))) -
                 1));
      }
    }();

    return from_elements(
        !is_positive(),
        static_cast<repr_type>(exp +
                               static_cast<signed_repr_type>(exponent_bias)),
        mantissa);
  }

private:
  constexpr static repr_type
  get_representation(T x) noexcept {
    auto const sign_repr = [&] {
      if (x < T(0)) {
        x = -x;
        return repr_type(1) << (sizeof(repr_type) *
                                    std::numeric_limits<unsigned char>::digits -
                                1);
      } else {
        return repr_type(0);
      }
    }();

    auto const exp_mantissa = [&] {
      if (x == T(0)) {
        return repr_type(0);
      } else if (x == std::numeric_limits<T>::infinity()) {
        return infinity;
      } else if (x != x) {
        return nan;
      } else {
        signed_repr_type exponent = 0;
        repr_type mantissa = 0;

        auto raw_mantissa = x;
        if (x >= T(2)) {
          while (raw_mantissa >= T(2)) {
            raw_mantissa /= T(2);
            ++exponent;
          }
        } else {
          while (raw_mantissa < T(1)) {
            raw_mantissa *= T(2);
            --exponent;
          }
        }

        raw_mantissa -= T(1);
        for (std::size_t mantissa_bit = 0; mantissa_bit < mantissa_bits;
             ++mantissa_bit) {
          mantissa <<= 1;
          raw_mantissa *= T(2);
          if (raw_mantissa >= T(1)) {
            mantissa |= repr_type(1);
            raw_mantissa -= T(1);
          }
        }

        return (((repr_type(1) << (exponent_bits - 1)) - 1 +
                 static_cast<repr_type>(exponent))
                << mantissa_bits) |
               mantissa;
      }
    }();

    // Exponent and mantissa value must not have the sign bit set
    assert((exp_mantissa &
            (repr_type(1) << (sizeof(repr_type) *
                                  std::numeric_limits<unsigned char>::digits -
                              1))) == repr_type(0));
    return sign_repr | exp_mantissa;
  }

  constexpr static FloatingPoint
  from_elements(repr_type sign, repr_type exponent,
                repr_type mantissa) noexcept {
    assert((sign & ~repr_type(1)) == 0);
    assert((exponent & ~((repr_type(1) << exponent_bits) - 1)) == 0);
    assert((mantissa & ~((repr_type(1) << mantissa_bits) - 1)) == 0);

    return FloatingPoint(mantissa | (exponent << mantissa_bits) |
                         (sign << (mantissa_bits + exponent_bits)));
  }

  repr_type repr;
};

using Float = FloatingPoint<float>;
using Double = FloatingPoint<double>;

} // namespace cte

static_assert(cte::Float(0.f).is_positive());
static_assert(cte::Double(0.).is_positive());
static_assert(cte::Float(123.f).is_positive());
static_assert(cte::Double(123.).is_positive());
static_assert(cte::Float(.0005f).is_positive());
static_assert(cte::Double(.0005).is_positive());
static_assert(cte::Float(64.00001f).is_positive());
static_assert(not cte::Float(-123.f).is_positive());
static_assert(not cte::Double(-123.).is_positive());
static_assert(not cte::Float(-0.0005f).is_positive());
static_assert(not cte::Double(-0.0005).is_positive());
static_assert(not cte::Float(-64.00001f).is_positive());

static_assert(cte::Float(0.f).exponent() == 0u);
static_assert(cte::Double(0.).exponent() == 0ul);
static_assert(cte::Float(123.f).exponent() == 133u);
static_assert(cte::Double(123.).exponent() == 1029ul);
static_assert(cte::Float(0.0005f).exponent() == 116u);
static_assert(cte::Double(0.0005).exponent() == 1012ul);
static_assert(cte::Float(64.00001f).exponent() == 133u);

static_assert(cte::Float(0.f).mantissa() == 0u);
static_assert(cte::Double(0.).mantissa() == 0ul);
static_assert(cte::Float(123.f).mantissa() == 7733248u);
static_assert(cte::Double(123.).mantissa() == 4151755906482176ul);
static_assert(cte::Float(0.0005f).mantissa() == 201327u);
static_assert(cte::Double(0.0005).mantissa() == 108086391056892ul);
static_assert(cte::Float(64.00001f).mantissa() == 1u);

static_assert(cte::Float(std::numeric_limits<float>::infinity()).is_inf());
static_assert(cte::Float(-std::numeric_limits<float>::infinity()).is_inf());
static_assert(cte::Float(std::numeric_limits<float>::quiet_NaN()).is_nan());
static_assert(cte::Double(std::numeric_limits<double>::infinity()).is_inf());
static_assert(cte::Double(-std::numeric_limits<double>::infinity()).is_inf());
static_assert(cte::Double(std::numeric_limits<double>::quiet_NaN()).is_nan());

// GCC bug #88173 does not allow to compile the following tests
#ifndef __GNUC__
static_assert(cte::Float(-std::numeric_limits<float>::quiet_NaN()).is_nan());
static_assert(cte::Double(-std::numeric_limits<double>::quiet_NaN()).is_nan());
#endif

static_assert(static_cast<float>(cte::Float(0.f)) == 0.f);
static_assert(static_cast<float>(cte::Float(123.f)) == 123.f);
static_assert(static_cast<float>(cte::Float(-123.f)) == -123.f);
static_assert(static_cast<float>(cte::Float(.0005f)) == .0005f);
static_assert(static_cast<float>(cte::Float(-.0005f)) == -.0005f);
static_assert(static_cast<float>(cte::Float(64.00001f)) == 64.00001f);
static_assert(static_cast<float>(cte::Float(-64.00001f)) == -64.00001f);
static_assert(
    static_cast<float>(cte::Float(std::numeric_limits<float>::infinity())) ==
    std::numeric_limits<float>::infinity());
static_assert(
    static_cast<float>(cte::Float(-std::numeric_limits<float>::infinity())) ==
    -std::numeric_limits<float>::infinity());

static_assert(static_cast<double>(cte::Double(0.f)) == 0.f);
static_assert(static_cast<double>(cte::Double(123.f)) == 123.f);
static_assert(static_cast<double>(cte::Double(-123.f)) == -123.f);
static_assert(static_cast<double>(cte::Double(.0005f)) == .0005f);
static_assert(static_cast<double>(cte::Double(-.0005f)) == -.0005f);
static_assert(static_cast<double>(cte::Double(64.00001f)) == 64.00001f);
static_assert(static_cast<double>(cte::Double(-64.00001f)) == -64.00001f);
static_assert(
    static_cast<double>(cte::Double(std::numeric_limits<double>::infinity())) ==
    std::numeric_limits<double>::infinity());
static_assert(static_cast<double>(
                  cte::Double(-std::numeric_limits<double>::infinity())) ==
              -std::numeric_limits<double>::infinity());

static_assert(static_cast<float>(cte::Float(1.12345f).truncate()) == 1.f);
static_assert(static_cast<float>(cte::Float(12.78698f).truncate()) == 12.f);
static_assert(static_cast<float>(cte::Float(0.12345f).truncate()) == 0.f);
static_assert(
    static_cast<float>(cte::Float(1234567890123456789.f).truncate()) ==
    1234567890123456789.f);
static_assert(static_cast<float>(cte::Float(123.456789f).truncate()) == 123.f);
static_assert(static_cast<float>(cte::Float(1.12345f).truncate()) == 1.f);
static_assert(static_cast<float>(cte::Float(0.12345f).truncate()) == 0.f);
static_assert(
    static_cast<float>(cte::Float(1234567890123456789.f).truncate()) ==
    1234567890123456789.f);
static_assert(static_cast<float>(cte::Float(123.456789f).truncate()) == 123.f);

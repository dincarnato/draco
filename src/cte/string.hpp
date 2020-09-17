#pragma once

#include <array>
#include <cassert>
#include <string_view>
#include <type_traits>

#include "floating_point.hpp"

namespace cte {

enum class StringArgType {
  Literal,
  Float,
  Double,
};

namespace detail {

template <auto n, StringArgType arg_type>
constexpr std::size_t
get_float_integral_part_size() noexcept {
  static_assert(arg_type != StringArgType::Literal);

  using value_type =
      std::conditional_t<arg_type == StringArgType::Float, float, double>;
  using raw_value_type =
      std::conditional_t<arg_type == StringArgType::Float, Float, Double>;

  auto const raw_value = raw_value_type::from_representation(n);

  std::size_t integral_part_size = 0;

  auto tmp = static_cast<value_type>(raw_value);
  if (tmp < value_type(0)) {
    tmp = -tmp;
    ++integral_part_size;
  }

  if (tmp < 1) {
    ++integral_part_size;
  } else {
    while (tmp >= value_type(1)) {
      tmp /= value_type(10);
      ++integral_part_size;
    }
  }

  return integral_part_size;
}

struct AdjustIntegralPartSizeResult {
  std::size_t decimal_part_size;
  std::size_t exponent_part;
};

template <auto n, StringArgType arg_type>
constexpr AdjustIntegralPartSizeResult
adjust_integral_part_size(std::size_t& integral_part_size) noexcept {
  constexpr std::size_t integral_part_max_size = 8;

  AdjustIntegralPartSizeResult result{0, 0};

  if (integral_part_size > integral_part_max_size) {
    result.decimal_part_size = integral_part_size - 1;
    result.exponent_part = integral_part_size - 1;
    integral_part_size = 1;
  }

  return result;
}

template <auto n, StringArgType arg_type>
constexpr std::size_t
get_float_decimal_part_size(std::size_t initial_decimal_part_size) noexcept {
  using raw_value_type =
      std::conditional_t<arg_type == StringArgType::Float, Float, Double>;

  constexpr std::size_t decimal_part_max_size = [] {
    if constexpr (arg_type == StringArgType::Float) {
      return std::size_t(6);
    } else {
      return std::size_t(10);
    }
  }();

  auto const raw_value = raw_value_type::from_representation(n);
  std::size_t decimal_part_size = 0;

  auto const exp = raw_value.exponent_unbiased();
  if (exp < static_cast<typename raw_value_type::signed_repr_type>(
                raw_value_type::mantissa_bits)) {
    auto mantissa =
        raw_value.mantissa() &
        ((typename raw_value_type::repr_type(1)
          << (raw_value_type::mantissa_bits - static_cast<std::size_t>(exp))) -
         1);
    for (std::size_t mantissa_bit_index = 0;
         mantissa_bit_index < raw_value_type::mantissa_bits;
         ++mantissa_bit_index) {
      if ((mantissa & typename raw_value_type::repr_type(1)) == 1) {
        decimal_part_size = std::min(
            raw_value_type::mantissa_bits - mantissa_bit_index -
                static_cast<std::size_t>(exp) + initial_decimal_part_size,
            decimal_part_max_size);
        break;
      }
      mantissa >>= 1;
    }

    if (decimal_part_size != 0)
      ++decimal_part_size;
  }

  return decimal_part_size;
}

constexpr std::size_t
get_float_exponent_part_size(std::size_t exponent_part) noexcept {
  if (exponent_part == 0) {
    return std::size_t(0);
  } else {
    std::size_t exponent_part_size = 1;
    while (exponent_part != 0) {
      exponent_part /= 10;
      ++exponent_part_size;
    }

    return exponent_part_size;
  }
}

template <auto n, StringArgType arg_type>
constexpr std::size_t
string_representation_size() noexcept {
  if constexpr (arg_type == StringArgType::Literal) {
    using n_type = std::decay_t<decltype(n)>;
    auto repr_size = [&] {
      if constexpr (not std::is_same_v<n_type, bool>) {
        std::size_t repr_size = 0;
        auto x = n;
        do {
          ++repr_size;
          x /= 10;
        } while (x != 0);

        return repr_size;
      } else {
        return std::size_t(1);
      }
    }();

    if constexpr (std::is_signed_v<n_type>)
      ++repr_size;

    return repr_size;
  } else {
    std::size_t integral_part_size =
        get_float_integral_part_size<n, arg_type>();
    std::size_t decimal_part_size = 0;
    std::size_t exponent_part = 0;
    {
      auto result = adjust_integral_part_size<n, arg_type>(integral_part_size);
      decimal_part_size = result.decimal_part_size;
      exponent_part = result.exponent_part;
    }

    decimal_part_size =
        get_float_decimal_part_size<n, arg_type>(decimal_part_size);

    auto const exponent_part_size = get_float_exponent_part_size(exponent_part);

    return integral_part_size + decimal_part_size + exponent_part_size;
  }
}

} // namespace detail

template <std::size_t N>
struct string {
  template <std::size_t>
  friend struct string;

  using iterator = char*;
  using const_iterator = char const*;
  using reverse_iterator = iterator;
  using const_reverse_iterator = const_iterator;
  static constexpr std::size_t max_size = N;

  constexpr string() noexcept : data() {}
  constexpr string(char const (&s)[N + 1]) noexcept : data() {
    char const* s_iter = &s[0];
    for (auto iter = std::begin(data), end = std::end(data); iter < end;
         ++iter, ++s_iter) {
      *iter = *s_iter;
    }
  }

  constexpr string(std::array<char, N + 1> const& data) noexcept : data(data) {
    assert(data.back() == '\0');
  }

  constexpr string(std::array<char, N + 1>&& data) noexcept
      : data(std::move(data)) {
    assert(data.back() == '\0');
  }

  constexpr operator std::string_view() const noexcept {
    return std::string_view(data.data(), N);
  }

  constexpr char const*
  c_str() const noexcept {
    return data.data();
  }

  constexpr string<N>
  replace(char old_char, char new_char) const noexcept {
    string<N> new_string;
    auto old_iter = std::begin(data);
    auto const old_iter_end = std::prev(std::end(data));
    auto new_iter = std::begin(new_string.data);
    for (; old_iter < old_iter_end; ++old_iter, ++new_iter) {
      if (auto const c = *old_iter; c == old_char) {
        *new_iter = new_char;
      } else {
        *new_iter = c;
      }
    }
    new_string.data.back() = '\0';

    return new_string;
  }

  template <std::size_t M>
  constexpr string<M + N - 1>
  append(char const (&other)[M]) const noexcept {
    // literal strings includes the null terminator
    return append(string<M - 1>(other));
  }

  template <std::size_t M>
  constexpr string<M + N>
  append(string<M> const& other) const noexcept {
    string<M + N> new_string;

    auto new_iter = std::begin(new_string.data);
    for (auto iter = std::begin(data), iter_end = std::prev(std::end(data));
         iter < iter_end and *iter != '\0'; ++iter, ++new_iter) {
      *new_iter = *iter;
    }

    for (auto iter = std::begin(other.data),
              iter_end = std::prev(std::end(other.data));
         iter < iter_end and *iter != '\0'; ++iter, ++new_iter) {
      *new_iter = *iter;
    }

    auto const new_string_end = std::end(new_string);
    for (; new_iter < new_string_end; ++new_iter) {
      *new_iter = '\0';
    }
    return new_string;
  }

  constexpr iterator
  begin() noexcept {
    return data.begin();
  }
  constexpr const_iterator
  begin() const noexcept {
    return data.begin();
  }

  constexpr iterator
  end() noexcept {
    return data.end();
  }

  constexpr const_iterator
  end() const noexcept {
    return data.end();
  }

  constexpr std::size_t
  size() const noexcept {
    for (std::size_t index = 0; index < N; ++index) {
      if (data[index] == '\0')
        return index;
    }

    return N;
  }

  constexpr bool
  is_empty() const noexcept {
    if constexpr (N == 0)
      return true;
    else
      return data[0] == '\0';
  }

private:
  std::array<char, N + 1> data;
};

// literal strings includes the null terminator
template <std::size_t N>
string(char const (&)[N])->string<N - 1>;

template <auto t, StringArgType arg_type = StringArgType::Literal>
constexpr auto
into_string() {
  constexpr std::size_t string_size =
      detail::string_representation_size<t, arg_type>();
  std::array<char, string_size + 1> data = {};

  if constexpr (arg_type == StringArgType::Literal) {
    using integral_type = std::decay_t<decltype(t)>;
    integral_type integral = t;

    if constexpr (std::is_same_v<integral_type, bool>) {
      data[0] = '0' + static_cast<char>(integral);
    } else {
      auto iter = std::next(std::rbegin(data));
      if constexpr (std::is_signed_v<integral_type>) {
        if (integral >= 0) {
          *iter++ = '\0';
        } else {
          data[0] = '-';
          integral = -1;
        }
      }

      do {
        *iter++ = '0' + static_cast<std::int8_t>(integral % 10);
        integral /= 10;
      } while (integral != 0);
    }
    data.back() = '\0';
  } else {
    using raw_value_type =
        std::conditional_t<arg_type == StringArgType::Float, Float, Double>;
    using value_type =
        std::conditional_t<arg_type == StringArgType::Float, float, double>;

    auto const raw_value = raw_value_type::from_representation(t);
    std::size_t char_index = 0;

    if (not raw_value.is_positive()) {
      data[0] = '-';
      char_index = 1;
    }

    auto integral_part_size =
        detail::get_float_integral_part_size<t, arg_type>();
    std::size_t decimal_part_size = 0;
    std::size_t exponent_part = 0;
    {
      auto result =
          detail::adjust_integral_part_size<t, arg_type>(integral_part_size);
      decimal_part_size = result.decimal_part_size;
      exponent_part = result.exponent_part;
    }

    // Scientific notation is still unsupported
    assert(exponent_part == 0);

    auto integral_part = [&] {
      if (exponent_part == 0) {
        if (raw_value.is_positive()) {
          return static_cast<std::uint64_t>(
              static_cast<value_type>(raw_value.truncate()));
        } else {
          --integral_part_size;
          assert(char_index == 1);
          return static_cast<std::uint64_t>(
              -static_cast<value_type>(raw_value.truncate()));
        }
      } else {
        auto value = static_cast<value_type>(raw_value);
        for (std::size_t index = 0; index < exponent_part; ++index) {
          value /= value_type(10);
        }
        if (value < 0) {
          --integral_part_size;
          assert(char_index == 1);
          value *= value_type(-1);
        }

        return static_cast<std::uint64_t>(
            static_cast<value_type>(raw_value_type(value).truncate()));
      }
    }();
    assert(integral_part_size > 0);

    for (std::size_t index = char_index + integral_part_size;
         index > char_index;) {
      --index;
      assert(integral_part_size == 1 or integral_part != 0);
      data[index] = '0' + static_cast<char>(integral_part % 10);
      integral_part /= 10;
    }

    decimal_part_size =
        detail::get_float_decimal_part_size<t, arg_type>(decimal_part_size);
    char_index += integral_part_size;
    if (decimal_part_size != 0) {
      --decimal_part_size;
      data[char_index++] = '.';

      auto const starting_index = char_index;
      auto tmp = static_cast<value_type>(raw_value);
      if (tmp < 0)
        tmp *= -1;

      tmp -= static_cast<value_type>(raw_value_type(tmp).truncate());
      for (std::size_t index = 0; index < decimal_part_size;
           ++index, ++char_index) {
        tmp *= value_type(10);
        auto const truncated_value =
            static_cast<value_type>(raw_value_type(tmp).truncate());
        assert(truncated_value >= 0 and truncated_value <= 9);
        data[char_index] = '0' + static_cast<char>(truncated_value);
        tmp -= truncated_value;
      }

      for (std::size_t index = char_index - 1; index > starting_index;
           --index) {
        if (data[index] == '0')
          data[index] = '\0';
        else
          break;
      }
    }
    assert(char_index < data.size());

    while (char_index < data.size())
      data[char_index++] = '\0';
  }

  return string<string_size>(std::move(data));
}

template <std::size_t N, std::size_t M>
constexpr bool
operator==(string<N> const& lhs, string<M> const& rhs) noexcept {
  auto lhs_iter = std::begin(lhs);
  auto const lhs_end_iter = std::end(lhs);
  auto rhs_iter = std::begin(rhs);
  auto const rhs_end_iter = std::end(rhs);

  for (; lhs_iter < lhs_end_iter and rhs_iter < rhs_end_iter;
       ++lhs_iter, ++rhs_iter) {
    auto const lhs_char = *lhs_iter;
    auto const rhs_char = *lhs_iter;
    if (not(lhs_char == rhs_char))
      return false;
    else if (lhs_char == '\0' and rhs_char == '\0')
      return true;
  }

  return false;
}

template <std::size_t N, std::size_t M>
constexpr bool
operator!=(string<N> const& lhs, string<M> const& rhs) noexcept {
  return not(lhs == rhs);
}

namespace detail {

template <typename T, typename = void>
struct into_string_helper;

template <typename T>
struct into_string_helper<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  using type = cte::FloatingPoint<T>;
  static constexpr StringArgType arg_type =
      std::is_same_v<T, float> ? StringArgType::Float : StringArgType::Double;
  static_assert(arg_type == StringArgType::Float or std::is_same_v<T, double>,
                "only float and double types are supported for now");

  static constexpr typename type::repr_type
  representation(T t) noexcept {
    return type(t).representation();
  }
};

template <typename T>
struct into_string_helper<T, std::enable_if_t<std::is_literal_type_v<T>>> {
  using type = T;
  static constexpr StringArgType arg_type = StringArgType::Literal;

  static constexpr T
  representation(T t) noexcept {
    return t;
  }
};

} // namespace detail
} // namespace cte

#define CTE_INTO_STRING(x)                                                     \
  cte::into_string<                                                            \
      cte::detail::into_string_helper<decltype(1)>::representation(1),         \
      cte::detail::into_string_helper<decltype(1)>::arg_type>()

static_assert(cte::into_string<true>() == cte::string("1"));
static_assert(cte::into_string<false>() == cte::string("0"));
static_assert(cte::into_string<0>() == cte::string("0"));
static_assert(std::is_same_v<decltype(cte::into_string<0>()), cte::string<2>>);
static_assert(cte::into_string<9>() == cte::string("9"));
static_assert(cte::into_string<10>() == cte::string("10"));
static_assert(cte::into_string<12345678>() == cte::string("12345678"));
static_assert(cte::into_string<-12345678>() == cte::string("-12345678"));
static_assert(cte::into_string<char(-42)>() == cte::string("-42"));
static_assert(cte::into_string<char(42)>() == cte::string("42"));
static_assert(cte::into_string<static_cast<unsigned char>(42)>() ==
              cte::string("42"));
static_assert(cte::into_string<42u>() == cte::string("42"));
static_assert(cte::into_string<-42l>() == cte::string("-42"));
static_assert(cte::into_string<42l>() == cte::string("42"));
static_assert(cte::into_string<42ul>() == cte::string("42"));
static_assert(cte::into_string<-42ll>() == cte::string("-42"));
static_assert(cte::into_string<42ll>() == cte::string("42"));
static_assert(cte::into_string<42ull>() == cte::string("42"));

static_assert(
    cte::detail::string_representation_size<cte::Float(12.f).representation(),
                                            cte::StringArgType::Float>() == 2);
static_assert(
    cte::detail::string_representation_size<cte::Float(12.5f).representation(),
                                            cte::StringArgType::Float>() == 4);
static_assert(
    cte::detail::string_representation_size<cte::Float(12.75f).representation(),
                                            cte::StringArgType::Float>() == 5);
static_assert(
    cte::detail::string_representation_size<
        cte::Float(12.875f).representation(), cte::StringArgType::Float>() ==
    6);
static_assert(
    cte::detail::string_representation_size<
        cte::Float(12.9375f).representation(), cte::StringArgType::Float>() ==
    7);
static_assert(
    cte::detail::string_representation_size<
        cte::Float(12.78125f).representation(), cte::StringArgType::Float>() ==
    8);
static_assert(
    cte::detail::string_representation_size<
        cte::Float(-12.78125f).representation(), cte::StringArgType::Float>() ==
    9);

static_assert(
    cte::detail::string_representation_size<cte::Float(-12.1f).representation(),
                                            cte::StringArgType::Float>() == 10);

static_assert(cte::into_string<cte::Float(12.f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12"));
static_assert(cte::into_string<cte::Float(12.5f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12.5"));
static_assert(cte::into_string<cte::Float(12.75f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12.75"));
static_assert(cte::into_string<cte::Float(12.875f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12.875"));
static_assert(cte::into_string<cte::Float(12.9375f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12.9375"));
static_assert(cte::into_string<cte::Float(12.78125f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("12.78125"));
static_assert(cte::into_string<cte::Float(-12.78125f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("-12.78125"));
static_assert(cte::into_string<cte::Float(0.5f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("0.5"));
static_assert(cte::into_string<cte::Float(-0.5f).representation(),
                               cte::StringArgType::Float>() ==
              cte::string("-0.5"));

static_assert(cte::into_string<cte::Double(12.).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12"));
static_assert(cte::into_string<cte::Double(12.5).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12.5"));
static_assert(cte::into_string<cte::Double(12.75).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12.75"));
static_assert(cte::into_string<cte::Double(12.875).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12.875"));
static_assert(cte::into_string<cte::Double(12.9375).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12.9375"));
static_assert(cte::into_string<cte::Double(12.78125).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("12.78125"));
static_assert(cte::into_string<cte::Double(-12.78125).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("-12.78125"));
static_assert(cte::into_string<cte::Double(0.5).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("0.5"));
static_assert(cte::into_string<cte::Double(-0.5).representation(),
                               cte::StringArgType::Double>() ==
              cte::string("-0.5"));

static_assert(CTE_INTO_STRING(1) == cte::string("1"));
static_assert(CTE_INTO_STRING(1.f) == cte::string("1"));
static_assert(CTE_INTO_STRING(1.) == cte::string("1"));
static_assert(CTE_INTO_STRING(1.5f) == cte::string("1.5"));
static_assert(CTE_INTO_STRING(1.5) == cte::string("1.5"));

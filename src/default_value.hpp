#pragma once

#include <cassert>
#include <type_traits>

#include "cte/floating_point.hpp"
#include "cte/string.hpp"

enum class DefaultValueType {
  Integral,
  Float,
  Double,
  String,
  None,
};

template <auto value, DefaultValueType type>
struct DefaultValue;

template <auto _value>
struct DefaultValue<_value, DefaultValueType::Integral> {
  using value_type = decltype(_value);
  static_assert(std::is_integral_v<value_type>);
  constexpr static bool is_available = true;

  constexpr DefaultValue() = default;
  template <typename... Ts>
  constexpr DefaultValue(Ts...) noexcept {}

  constexpr value_type
  value() const noexcept {
    return _value;
  }

  constexpr auto
  into_string() const noexcept {
    return cte::into_string<_value>();
  }
};

template <auto repr>
struct DefaultValue<repr, DefaultValueType::Float> {
  using value_type = float;
  using repr_type = typename cte::FloatingPoint<float>::repr_type;
  static_assert(std::is_same_v<repr_type, decltype(repr)>);
  constexpr static bool is_available = true;

  constexpr DefaultValue() = default;
  template <typename... Ts>
  constexpr DefaultValue(Ts...) noexcept {}

  constexpr value_type
  value() const noexcept {
    return static_cast<value_type>(
        cte::FloatingPoint<float>::from_representation(repr));
  }

  constexpr auto
  into_string() const noexcept {
    return cte::into_string<repr, cte::StringArgType::Float>();
  }
};

template <auto repr>
struct DefaultValue<repr, DefaultValueType::Double> {
  using value_type = double;
  using repr_type = typename cte::FloatingPoint<double>::repr_type;
  static_assert(std::is_same_v<repr_type, decltype(repr)>);
  constexpr static bool is_available = true;

  constexpr DefaultValue() = default;
  template <typename... Ts>
  constexpr DefaultValue(Ts...) noexcept {}

  constexpr value_type
  value() const noexcept {
    return static_cast<value_type>(
        cte::FloatingPoint<double>::from_representation(repr));
  }

  constexpr auto
  into_string() const noexcept {
    return cte::into_string<repr, cte::StringArgType::Double>();
  }
};

template <std::size_t size>
struct DefaultValue<size, DefaultValueType::String> {
  using value_type = cte::string<size>;
  constexpr static bool is_available = true;

  constexpr DefaultValue(cte::string<size> value) noexcept
      : _value(std::move(value)) {}

  constexpr DefaultValue(char const (&value)[size + 1]) noexcept
      : _value(cte::string(value)) {}

  constexpr value_type
  value() const noexcept {
    return _value;
  }

  constexpr auto
  into_string() const noexcept {
    return _value;
  }

private:
  cte::string<size> _value;
};

template <>
struct DefaultValue<0, DefaultValueType::None> {
  constexpr static bool is_available = false;
  constexpr DefaultValue() = default;
};

template <typename T>
struct is_default_value : std::false_type {};

template <auto value, DefaultValueType value_type>
struct is_default_value<DefaultValue<value, value_type>> : std::true_type {};

template <typename T>
constexpr bool is_default_value_v = is_default_value<T>::value;

template <typename T, typename D, typename = void>
struct is_convertible_from_default_value : std::false_type {};

template <typename T, typename D>
struct is_convertible_from_default_value<
    T, D,
    std::enable_if_t<is_default_value_v<D> and
                     std::is_convertible_v<typename D::value_type, T>>>
    : std::true_type {};

template <typename T, std::size_t size>
struct is_convertible_from_default_value<
    T, DefaultValue<size, DefaultValueType::String>,
    std::enable_if_t<std::is_convertible_v<char const*, T>>> : std::true_type {
};

template <typename T, typename D>
constexpr bool is_convertible_from_default_value_v =
    is_convertible_from_default_value<T, D>::value;

namespace detail {

template <typename T, typename = void>
struct DefaultValueArgValue;

template <typename T>
struct DefaultValueArgValue<
    T, std::enable_if_t<std::is_integral_v<T> and
                        not std::is_same_v<std::decay_t<T>, char const*>>> {
  using value_type = T;
  constexpr static DefaultValueType type = DefaultValueType::Integral;

  constexpr DefaultValueArgValue(T t) noexcept : _value(t) {}

  constexpr value_type
  value() const noexcept {
    return _value;
  }

private:
  value_type _value;
};

template <typename T>
struct DefaultValueArgValue<T, std::enable_if_t<std::is_floating_point_v<T>>> {
  using floating_point_type = cte::FloatingPoint<T>;
  using value_type = typename floating_point_type::repr_type;

  constexpr static DefaultValueType type = std::is_same_v<T, float>
                                               ? DefaultValueType::Float
                                               : DefaultValueType::Double;

  constexpr DefaultValueArgValue(T t) noexcept
      : _value(floating_point_type(t).representation()) {}

  constexpr value_type
  value() const noexcept {
    return _value;
  }

private:
  value_type _value;
};

template <std::size_t N>
struct DefaultValueArgValue<char const (&)[N], void> {
  using value_type = std::size_t;
  constexpr static DefaultValueType type = DefaultValueType::String;

  constexpr DefaultValueArgValue([[maybe_unused]] char const (&s)[N]) noexcept {
    assert(s[N - 1] == '\0');
  }

  constexpr value_type
  value() const noexcept {
    return N - 1;
  }
};

template <std::size_t N>
struct DefaultValueArgValue<cte::string<N>, void> {
  using value_type = std::size_t;
  constexpr static DefaultValueType type = DefaultValueType::String;

  constexpr DefaultValueArgValue(cte::string<N>) noexcept {}

  constexpr value_type
  value() const noexcept {
    return N;
  }
};

} // namespace detail

#define MAKE_DEFAULT_VALUE(x)                                                  \
  ::DefaultValue<::detail::DefaultValueArgValue<decltype(x)>(x).value(),       \
                 ::detail::DefaultValueArgValue<decltype(x)>::type>(x)

#define MAKE_NO_DEFAULT_VALUE ::DefaultValue<0, ::DefaultValueType::None>()

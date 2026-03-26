#pragma once

#include "cte/string.hpp"
#include "default_value.hpp"
#include "results/jsonify_base.hpp"

#include <exception>
#include <ostream>
#include <stdexcept>
#include <string>

namespace args {

enum class WeightsInitialization {
  Kmeans,
  Uniform,
  Random,
};

} // namespace args

template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &
operator>>(std::basic_istream<CharT, Traits> &in,
           ::args::WeightsInitialization &weights_initialization) {
  std::string tmp;
  in >> tmp;
  if (tmp == "uniform") {
    weights_initialization = ::args::WeightsInitialization::Uniform;
  } else if (tmp == "kmeans") {
    weights_initialization = ::args::WeightsInitialization::Kmeans;
  } else if (tmp == "random") {
    weights_initialization = ::args::WeightsInitialization::Random;
  } else {
    throw std::runtime_error("invalid weights initialization");
  }

  return in;
}

namespace results {

template <typename CharT, typename Traits, typename T>
std::enable_if_t<detail::is_jsonificable_v<std::decay_t<T>>,
                 std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, T &&t);

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os,
        ::args::WeightsInitialization const &weights_initialization) {
  switch (weights_initialization) {
  case ::args::WeightsInitialization::Uniform:
    return jsonify(os, "uniform");
  case ::args::WeightsInitialization::Kmeans:
    return jsonify(os, "kmeans");
  case ::args::WeightsInitialization::Random:
    return jsonify(os, "random");
  default:
    throw std::runtime_error("unreachable");
  }
}

} // namespace results

namespace detail {
template <> struct DefaultValueArgValue<::args::WeightsInitialization> {
  using value_type = ::args::WeightsInitialization;
  constexpr static DefaultValueType type = DefaultValueType::String;

  constexpr DefaultValueArgValue(::args::WeightsInitialization value) noexcept
      : value_(value) {}

  constexpr value_type value() const noexcept { return value_; }

private:
  value_type value_;
};

} // namespace detail

template <::args::WeightsInitialization weights>
struct DefaultValue<weights, DefaultValueType::String> {
  using value_type = ::args::WeightsInitialization;
  using repr_type = ::args::WeightsInitialization;
  constexpr static bool is_available = true;

  constexpr DefaultValue() = default;
  template <typename... Ts>
  constexpr DefaultValue(Ts... args) noexcept
      : value_(std::forward<Ts...>(args...)) {}

  constexpr value_type value() const noexcept { return value_; }

  constexpr auto into_string() const noexcept {
    switch (value_) {
    case ::args::WeightsInitialization::Uniform:
      return cte::string("uniform");
    case ::args::WeightsInitialization::Kmeans:
      return cte::string("kmeans\0");
    case ::args::WeightsInitialization::Random:
      return cte::string("random\0");
    default:
      std::terminate();
    }
  }

private:
  value_type value_;
};

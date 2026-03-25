#pragma once

#include "cte/string.hpp"
#include "default_value.hpp"
#include "results/jsonify_base.hpp"

#include <exception>
#include <ostream>
#include <stdexcept>
#include <string>

namespace args {

enum class PriorsInitialization {
  Uniform,
  Weights,
  Random,
};

} // namespace args

template <typename CharT, typename Traits>
std::basic_istream<CharT, Traits> &
operator>>(std::basic_istream<CharT, Traits> &in,
           ::args::PriorsInitialization &priors_initialization) {
  std::string tmp;
  in >> tmp;
  if (tmp == "uniform") {
    priors_initialization = ::args::PriorsInitialization::Uniform;
  } else if (tmp == "weights") {
    priors_initialization = ::args::PriorsInitialization::Weights;
  } else if (tmp == "random") {
    priors_initialization = ::args::PriorsInitialization::Random;
  } else {
    throw std::runtime_error("invalid priors initialization");
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
        ::args::PriorsInitialization const &priors_initialization) {
  switch (priors_initialization) {
  case ::args::PriorsInitialization::Uniform:
    return jsonify(os, "uniform");
  case ::args::PriorsInitialization::Weights:
    return jsonify(os, "weights");
  case ::args::PriorsInitialization::Random:
    return jsonify(os, "random");
  default:
    throw std::runtime_error("unreachable");
  }
}

} // namespace results

namespace detail {
template <> struct DefaultValueArgValue<::args::PriorsInitialization> {
  using value_type = ::args::PriorsInitialization;
  constexpr static DefaultValueType type = DefaultValueType::String;

  constexpr DefaultValueArgValue(::args::PriorsInitialization value) noexcept
      : value_(value) {}

  constexpr value_type value() const noexcept { return value_; }

private:
  value_type value_;
};

} // namespace detail

template <::args::PriorsInitialization priors>
struct DefaultValue<priors, DefaultValueType::String> {
  using value_type = ::args::PriorsInitialization;
  using repr_type = ::args::PriorsInitialization;
  constexpr static bool is_available = true;

  constexpr DefaultValue() = default;
  template <typename... Ts>
  constexpr DefaultValue(Ts... args) noexcept
      : value_(std::forward<Ts...>(args...)) {}

  constexpr value_type value() const noexcept { return value_; }

  constexpr auto into_string() const noexcept {
    switch (value_) {
    case ::args::PriorsInitialization::Uniform:
      return cte::string("uniform");
    case ::args::PriorsInitialization::Weights:
      return cte::string("weights");
    case ::args::PriorsInitialization::Random:
      return cte::string("random\0");
    default:
      std::terminate();
    }
  }

private:
  value_type value_;
};

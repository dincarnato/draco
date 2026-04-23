#pragma once

#include "weibull_fitter.hpp"

#include <cmath>
#include <numeric>

#include <dlib/optimization.h>

inline WeibullFitter::WeibullFitter(
    PerturbedEigengap const &perturbed_eigengap) noexcept
    : perturbed_eigengap(perturbed_eigengap) {}

inline WeibullParams WeibullFitter::fit() noexcept(false) {
  WeibullParams fit_params;
  auto const &perturbed_eigengap = this->perturbed_eigengap.get();

  auto const size = perturbed_eigengap.size();
  if (size == 0) {
    fit_params.shape = 0.;
    fit_params.scale = 0.;
    return fit_params;
  }

  bool updated_params = false;
  {
    if (last_k == 0.) {
      last_k = 1.2;
      log_k = 0.;
      updated_params = true;
    }

    if (last_l == 0.) {
      auto data = perturbed_eigengap;
      auto const data_begin = std::begin(data);
      auto median_iter =
          std::next(data_begin, static_cast<std::ptrdiff_t>(size / 2));
      std::nth_element(data_begin, median_iter, std::end(data));
      last_l = *median_iter / std::log(2.);
      log_l = std::log(last_l);
      updated_params = true;
    }

    if (updated_params)
      l_k_inv = 1. / std::pow(last_l, last_k);
  }

  if (size < last_size) {
    last_size = 0;
    log_x_cache.clear();
    cum_log_x = 0.;
    cum_x_k = 0.;
    cum_x_k_log_x = 0.;
  }

  if (last_size != size) {
    log_x_cache.resize(size);
    auto data_iter = std::next(std::begin(perturbed_eigengap),
                               static_cast<std::ptrdiff_t>(last_size));
    auto const data_end = std::end(perturbed_eigengap);
    auto log_iter = std::next(std::begin(log_x_cache),
                              static_cast<std::ptrdiff_t>(last_size));

    for (; data_iter < data_end; ++data_iter, ++log_iter) {
      double const x = *data_iter;
      assert(x > 0.);
      double const log_x = std::log(x);
      double const x_k = std::pow(x, last_k);

      *log_iter = log_x;
      cum_log_x += log_x;
      cum_x_k += x_k;
      cum_x_k_log_x += x_k * log_x;
    }
  }

  if (not updated_params and size == last_size and last_size != 0) {
    fit_params.shape = last_k;
    fit_params.scale = last_l;
    return fit_params;
  }

  last_size = size;
  params_type params{last_k, last_l};
  auto current_min_params = min_params;
  auto current_max_params = max_params;
  for (;;) {
    try {
      dlib::find_max_box_constrained(
          dlib::lbfgs_search_strategy(search_points),
          dlib::objective_delta_stop_strategy(delta_stop),
          [this](auto &&params) { return pdf_log_likelihood(params); },
          [this](auto &&params) { return pdf_log_likelihood_deriv(params); },
          params, current_min_params, current_max_params);
    } catch (std::exception const &) {
      break;
    }

    bool updated_limits = false;
    for (unsigned param_index = 0; param_index < 2; ++param_index) {
      double const param = params(param_index);
      auto &&current_max_param = current_max_params(param_index);
      auto &&current_min_param = current_min_params(param_index);

      if (std::abs(param - current_max_param) <=
          current_max_params(param_index) * 0.05) {
        current_max_param = param * 2.;
        updated_limits = true;
      }

      if (std::abs(param - current_min_param) <=
              current_min_params(param_index) * 0.05 and
          current_min_params(param_index) > 1.) {
        current_min_param = param / 2.;
        updated_limits = true;
      }
    }

    if (not updated_limits)
      break;

    current_min_params(0) = std::max(current_min_params(0), 1.);
    current_max_params(0) = std::max(current_max_params(0), 1.);
  }

  fit_params.shape = params(0);
  fit_params.scale = params(1);
  return fit_params;
}

inline void WeibullFitter::update_k_expressions(std::size_t start_index,
                                                KExpr initial_data) noexcept {
  auto const &perturbed_eigengap = this->perturbed_eigengap.get();
  auto data_iter = std::next(std::begin(perturbed_eigengap),
                             static_cast<std::ptrdiff_t>(start_index));
  auto const data_end = std::end(perturbed_eigengap);
  auto log_iter = std::next(std::begin(log_x_cache),
                            static_cast<std::ptrdiff_t>(start_index));

  for (; data_iter < data_end; ++data_iter, ++log_iter) {
    double const x_k = std::pow(*data_iter, last_k);

    initial_data.x_k += x_k;
    initial_data.x_k_log_x += x_k * *log_iter;
  }

  cum_x_k = initial_data.x_k;
  cum_x_k_log_x = initial_data.x_k_log_x;
}

inline double
WeibullFitter::pdf_log_likelihood_slow(params_type const &params) noexcept {
  assert(this->perturbed_eigengap.get().size() == last_size);

  double const k = params(0);
  double const l = params(1);

  assert(k > 0.);
  assert(l > 0.);
  auto const l_k = std::pow(l, k);
  double k_over_l_k;
  if (l_k == 0)
    return std::numeric_limits<double>::infinity();
  else
    k_over_l_k = k / l_k;

  auto &&perturbed_eigengap = this->perturbed_eigengap.get();
  auto const result = std::accumulate(
      std::begin(perturbed_eigengap), std::end(perturbed_eigengap), 0.,
      [&](double acc, double x) {
        auto const y = std::log(k_over_l_k * std::pow(x, k - 1.) *
                                std::exp(-std::pow(x / l, k)));
        if (std::isinf(y)) {
          return std::numeric_limits<double>::lowest() /
                 static_cast<double>(perturbed_eigengap.size() + 1);
        } else {
          return acc + y;
        }
      });
  return result;
}

inline double
WeibullFitter::pdf_log_likelihood(params_type const &params) noexcept {
  assert(this->perturbed_eigengap.get().size() == last_size);

  double const k = params(0);
  double const l = params(1);
  auto const n = static_cast<double>(last_size);

  update_params(k, l);
  assert(k > 0.);
  assert(l > 0.);

  auto const result =
      n * (log_k - k * log_l) + (k - 1.) * cum_log_x - l_k_inv * cum_x_k;
  return result;
}

inline auto
WeibullFitter::pdf_log_likelihood_deriv(params_type const &params) noexcept
    -> params_type {
  assert(this->perturbed_eigengap.get().size() == last_size);

  double const k = params(0);
  double const l = params(1);
  auto const n = static_cast<double>(last_size);

  update_params(k, l);
  assert(k > 0.);
  assert(l > 0.);

  double const k_deriv = n * (1. / k - log_l) + cum_log_x +
                         l_k_inv * (log_l * cum_x_k - cum_x_k_log_x);
  double const l_deriv = k / l * (l_k_inv * cum_x_k - n);
  return params_type{k_deriv, l_deriv};
}

inline void WeibullFitter::update_params(double k, double l) noexcept {
  bool changed_params = false;

  if (k != last_k) {
    assert(k > 0.);
    last_k = k;
    log_k = std::log(k);
    update_k_expressions(0, KExpr{0., 0.});
    changed_params = true;
  }

  if (l != last_l) {
    assert(l > 0.);
    last_l = l;
    log_l = std::log(l);
    changed_params = true;
  }

  if (changed_params) {
    auto const l_k = std::pow(l, k);
    if (l_k == 0)
      l_k_inv = std::numeric_limits<double>::infinity();
    else
      l_k_inv = 1. / l_k;
  }
}

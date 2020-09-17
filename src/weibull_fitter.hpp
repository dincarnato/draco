#pragma once

#include "ptba_types.hpp"

#include <functional>

#include <dlib/matrix.h>

struct WeibullParams {
  double shape;
  double scale;
};

namespace test {
struct WeibullFitter;
} // namespace test

struct WeibullFitter {
  friend struct test::WeibullFitter;

  using params_type = dlib::matrix<double, 0, 1>;
  using data_iterator = PerturbedEigengap::const_iterator;

  WeibullFitter(PerturbedEigengap const& perturbed_eigengap) noexcept;
  WeibullParams fit() noexcept(false);

  static constexpr std::size_t search_points = 30;
  static constexpr double delta_stop = 1e-7;
  inline static const params_type min_params{1., 1e-4};
  inline static const params_type max_params{3., 2.};

private:
  struct KExpr {
    double x_k;
    double x_k_log_x;
  };

  double pdf_log_likelihood_slow(params_type const& params) noexcept;
  double pdf_log_likelihood(params_type const& params) noexcept;
  params_type pdf_log_likelihood_deriv(params_type const& params) noexcept;

  void update_k_expressions(std::size_t start_index,
                            KExpr initial_data) noexcept;
  void update_params(double k, double l) noexcept;

  std::reference_wrapper<const PerturbedEigengap> perturbed_eigengap;
  std::size_t last_size = 0;
  std::vector<double> log_x_cache;
  double last_k = 0.;
  double last_l = 0.;
  double log_k = 0.;
  double log_l = 0.;
  double l_k_inv = 0.;
  double cum_log_x = 0.;
  double cum_x_k = 0.;
  double cum_x_k_log_x = 0.;
};

#include "weibull_fitter_impl.hpp"

#include "test/weibull_fitter.hpp"

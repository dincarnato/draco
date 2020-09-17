#pragma once

#include "../weibull_fitter.hpp"

namespace test {

struct WeibullFitter {
  using params_type = ::WeibullFitter::params_type;

  static double
  pdf_log_likelihood(::WeibullFitter& weibull_fitter,
                     params_type const& params) noexcept {
    return weibull_fitter.pdf_log_likelihood(params);
  }

  static params_type
  pdf_log_likelihood_deriv(::WeibullFitter& weibull_fitter,
                           params_type const& params) noexcept {
    return weibull_fitter.pdf_log_likelihood_deriv(params);
  }
};

} // namespace test

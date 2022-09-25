#include "weibull_fitter.hpp"
#include "ptba_types.hpp"
#include "test/weibull_fitter.hpp"

#include <algorithm>
#include <cassert>
#include <random>

#include <boost/math/distributions/weibull.hpp>

static void test_pdf_log_likelihood(double shape, double scale) {
  constexpr std::size_t n_data = 2000;

  thread_local std::mt19937 random_gen(std::random_device{}());
  std::weibull_distribution<double> weibull_dist(shape, scale);

  PerturbedEigengap data(n_data);
  std::generate(std::begin(data), std::end(data),
                [&] { return weibull_dist(random_gen); });

  auto const log_likelihood = [&] {
    boost::math::weibull_distribution<double> dist(shape, scale);
    return std::accumulate(std::begin(data), std::end(data), 0.,
                           [&](double acc, double x) {
                             return acc + std::log(boost::math::pdf(dist, x));
                           });
  }();

  WeibullFitter fitter(data);
  fitter.fit();
  auto const fitted_log_likelihood = test::WeibullFitter::pdf_log_likelihood(
      fitter, WeibullFitter::params_type{shape, scale});
  assert(std::abs(log_likelihood - fitted_log_likelihood) < 0.1);
}

static void test_with_params(double shape, double scale) {
  constexpr std::size_t n_data = 2000;

  thread_local std::mt19937 random_gen(std::random_device{}());
  std::weibull_distribution<double> weibull_dist(shape, scale);

  PerturbedEigengap data(n_data);
  std::generate(std::begin(data), std::end(data),
                [&] { return weibull_dist(random_gen); });

  auto get_log_likelihood = [&] {
    boost::math::weibull_distribution<double> dist(shape, scale);
    return std::accumulate(std::begin(data), std::end(data), 0.,
                           [&](double acc, double x) {
                             return acc + std::log(boost::math::pdf(dist, x));
                           });
  };

  WeibullFitter fitter(data);

  auto check_result = [&] {
    auto result = fitter.fit();
    auto const fitted_log_likelihood = test::WeibullFitter::pdf_log_likelihood(
        fitter, WeibullFitter::params_type{result.shape, result.scale});
    double const log_likelihood = get_log_likelihood();
    assert(std::abs(log_likelihood - fitted_log_likelihood) <
           std::max(std::abs(log_likelihood) * 0.05, 5.));
  };

  check_result();

  data.resize(n_data * 2);
  std::copy(std::begin(data), std::next(std::begin(data), n_data),
            std::next(std::begin(data), n_data));
  check_result();

  data.resize(n_data / 2);
  check_result();
}

int main() {
  test_pdf_log_likelihood(1.42, 0.8);
  test_pdf_log_likelihood(3.5, 0.8);
  test_pdf_log_likelihood(1.2, 2.5);
  test_pdf_log_likelihood(3.5, 2.5);
  test_pdf_log_likelihood(1.2, 0.00008);
  test_pdf_log_likelihood(3.5, 0.00008);

  test_with_params(1.42, 0.8);
  test_with_params(3.5, 0.8);
  test_with_params(1.2, 2.5);
  test_with_params(3.5, 2.5);
  test_with_params(1.2, 0.00008);
  test_with_params(3.5, 0.00008);
}

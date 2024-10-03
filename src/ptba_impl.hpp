#pragma once

#include "kolmogorov_smirnov.hpp"
#include "ptba.hpp"

#include <algorithm>
#include <cstddef>

#include <boost/math/distributions.hpp>
#include <limits>

template <typename Distribution>
bool Ptba::is_distribution(
    Distribution &&distribution,
    PerturbedEigengap const &perturbed_data) noexcept(false) {

  assert(not perturbed_data.empty());
  std::size_t const data_size = perturbed_data.size();
  auto data = perturbed_data;
  std::sort(std::begin(data), std::end(data));

  std::vector<double> cumulative_data(data_size);
  std::copy(std::begin(data), std::end(data), std::begin(cumulative_data));
  std::partial_sum(std::begin(cumulative_data), std::end(cumulative_data),
                   std::begin(cumulative_data));
  std::transform(
      std::begin(cumulative_data), std::end(cumulative_data),
      std::begin(cumulative_data),
      [max = cumulative_data.back()](auto value) { return value / max; });

  std::vector<double> diffs(data_size);
  {
    auto data_iter = std::begin(data);
    auto const data_end = std::end(data);
    auto cumulative_data_iter = std::begin(cumulative_data);
    auto diffs_iter = std::begin(diffs);
    for (; data_iter < data_end;
         ++data_iter, ++diffs_iter, ++cumulative_data_iter) {
      auto value = *data_iter;
      auto cum_value = *cumulative_data_iter;
      auto &diff = *diffs_iter;

      auto cdf_value = boost::math::cdf(distribution, value);
      diff = std::abs(cdf_value - cum_value);
    }
  }

  auto const threshold =
      kolmogorov_smirnov_critical_value(static_cast<unsigned>(data_size), 0.10);
  auto const after_high_diff_index =
      static_cast<std::ptrdiff_t>(static_cast<double>(data_size) * 0.90);
  assert(after_high_diff_index > 0);
  assert(data_size < std::numeric_limits<std::ptrdiff_t>::max() and
         after_high_diff_index <= static_cast<std::ptrdiff_t>(data_size));
  auto const high_diff_iter =
      std::next(std::begin(diffs), after_high_diff_index - 1);
  std::nth_element(std::begin(diffs), high_diff_iter, std::end(diffs));
  return *high_diff_iter < threshold;
}

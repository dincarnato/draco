#pragma once

#include "../expectation_maximization.hpp"
#include "args.hpp"
#include "weighted_clusters.hpp"

#include <random>
#include <utility>

namespace test {
namespace expectation_maximization {

struct ExpectationMaximization : ::ExpectationMaximization {
  template <typename G>
    requires std::uniform_random_bit_generator<std::remove_cvref_t<G>>
  ExpectationMaximization(RingmapMatrix const &ringmap,
                          WeightedClusters &weights, Args const &args, G &&rng)
      : ::ExpectationMaximization(ringmap, weights, args,
                                  std::forward<G>(rng)) {}

  double expectation() noexcept {
    return ::ExpectationMaximization::expectation();
  }

  void maximization() noexcept { ::ExpectationMaximization::maximization(); }

  constexpr auto responsibilities(this auto &self) -> decltype(auto) {
    return (self.responsibilities_);
  }

  constexpr auto log_likelihoods(this auto &self) -> decltype(auto) {
    return (self.log_likelihoods_);
  }

  constexpr auto weights_buffer(this auto &self) -> decltype(auto) {
    return (self.weights_buffer_);
  }

  constexpr auto priors(this auto &self) -> decltype(auto) {
    return (self.priors_);
  }
};

} // namespace expectation_maximization

using ExpectationMaximization =
    expectation_maximization::ExpectationMaximization;

} // namespace test

#pragma once

#include "expectation_maximization/responsibilities.hpp"
#include "ringmap_matrix_row_accessor.hpp"
#include "weighted_clusters.hpp"

#include <random>
#include <type_traits>
#include <variant>
#include <vector>

struct CompactRingmap;
struct CompactRingmapRow;
class RingmapMatrix;
struct Args;

namespace expectation_maximization {

struct Converged {
  std::uint16_t after_iterations;
};

struct MaxIterations {};

using Convergence = std::variant<Converged, MaxIterations>;

struct Result {
  double log_likelihood;
  Convergence convergence;
};

struct ExpectationMaximization {
  template <typename G>
    requires std::uniform_random_bit_generator<std::remove_cvref_t<G>>
  ExpectationMaximization(RingmapMatrix const &ringmap,
                          WeightedClusters &weights, Args const &args, G &&rng);

  Result run() noexcept;

  /**
   * Writes to `assignments` the number of reads assigned to each cluster based
   * on the weights and the priors.
   *
   * `assignments` represents the output for the assignments. Its length must
   * be equal to the number of clusters.
   *
   * `buffer` is just a temporary buffer for calculations. Its length must be
   * equal to the number of clusters.
   */
  template <typename R>
  void read_assignment(RingmapMatrixRowAccessor<R> ringmap_row,
                       std::span<std::uint32_t> assignments,
                       std::span<double> buffer, std::mt19937 &rng) const;

  constexpr std::span<const double> priors() const noexcept { return priors_; }

protected:
  double expectation() noexcept;
  void maximization() noexcept;

  RingmapMatrix const *ringmap_;
  WeightedClusters *weights_;
  Args const *args_;
  std::vector<double> priors_;
  Responsibilities responsibilities_;
  Responsibilities weights_buffer_;
};

void weighted_priors_initialization(std::span<double> priors,
                                    RingmapMatrix const &ringmap,
                                    WeightedClusters &weights) noexcept;

} // namespace expectation_maximization

using ExpectationMaximization =
    expectation_maximization::ExpectationMaximization;

#include "args.hpp"
#include "logger.hpp"
#include "ringmap_matrix.hpp"
#include "span_formatter.hpp"

#include <algorithm>
#include <functional>
#include <span>

namespace expectation_maximization {

template <typename G>
  requires std::uniform_random_bit_generator<std::remove_cvref_t<G>>
ExpectationMaximization::ExpectationMaximization(RingmapMatrix const &ringmap,
                                                 WeightedClusters &weights,
                                                 Args const &args, G &&rng)
    : ringmap_(&ringmap), weights_(&weights), args_(&args),
      priors_(weights.getClustersSize()),
      responsibilities_(ringmap.rows_size(),
                        static_cast<std::uint8_t>(weights.getClustersSize())),
      weights_buffer_(static_cast<std::uint32_t>(weights.getElementsSize()),
                      static_cast<std::uint8_t>(weights.getClustersSize())) {
  if (weights.getElementsSize() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error(
        "too many bases to initialize expectation-maximization");
  }
  if (weights.getClustersSize() > std::numeric_limits<std::uint8_t>::max()) {
    throw std::runtime_error(
        "too many clusters to initialize expectation-maximization");
  }

  switch (args.expectation_maximization_priors_initialization()) {
  case args::PriorsInitialization::Uniform:
    std::ranges::fill(priors_,
                      1. / static_cast<double>(weights.getClustersSize()));
    break;
  case args::PriorsInitialization::Weights:
    weighted_priors_initialization(priors_, *ringmap_, *weights_);
    logger::trace("Initialization of priors by weights: {}",
                  SpanFormatter(priors_));
    break;
  case args::PriorsInitialization::Random:
    auto random_gen = std::uniform_real_distribution<float>(1e-6f, 1.f);
    std::ranges::generate(priors_, [&] { return random_gen(rng); });
    auto priors_sum = std::ranges::fold_left(priors_, 0., std::plus<>{});
    for (auto &prior : priors_) {
      prior /= priors_sum;
    }
    logger::trace("Initialization of random priors: {}",
                  SpanFormatter(priors_));
    break;
  }
}

template <typename R>
void ExpectationMaximization::read_assignment(
    RingmapMatrixRowAccessor<R> ringmap_row,
    std::span<std::uint32_t> assignments, std::span<double> buffer,
    std::mt19937 &rng) const {
  assert(std::size(assignments) == weights_->getClustersSize());

  double probability_sum = 0.;
  for (auto &&[probability, cluster_weights, cluster_prior] :
       std::views::zip(buffer, weights_->clusters(), priors_)) {
    probability =
        std::ranges::fold_left(ringmap_row.modifiedIndices() |
                                   std::views::transform([&](auto base_index) {
                                     return cluster_weights[base_index];
                                   }),
                               cluster_prior, std::multiplies{});
    probability_sum += probability;
  }

  for (auto &&[assignment, probability] :
       std::views::zip(assignments, buffer)) {
    assignment = static_cast<std::uint32_t>(std::round(
        static_cast<double>(std::size(ringmap_row.modifiedIndices())) *
        probability / probability_sum));
  }

  auto total_assignments = std::ranges::fold_left(
      assignments, static_cast<std::uint32_t>(0), std::plus{});
  auto assignments_difference = static_cast<std::int32_t>(
      static_cast<std::int64_t>(std::size(ringmap_row.modifiedIndices())) -
      static_cast<std::int64_t>(total_assignments));
  if (assignments_difference != 0) {
    std::uniform_int_distribution<std::uint8_t> chooser(
        static_cast<std::uint8_t>(0),
        static_cast<std::uint8_t>(weights_->getClustersSize() - 1));

    for (;;) {
      auto &assignment = assignments[chooser(rng)];
      auto new_assignment =
          static_cast<std::int32_t>(assignment) + assignments_difference;
      if (new_assignment < 0) {
        continue;
      }
      assignment = static_cast<std::uint32_t>(new_assignment);
      break;
    }
  }

  assert(std::ranges::fold_left(assignments, static_cast<std::uint32_t>(0),
                                std::plus{}) ==
         std::size(ringmap_row.modifiedIndices()));
}

} // namespace expectation_maximization

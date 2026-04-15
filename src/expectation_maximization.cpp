#include "expectation_maximization.hpp"
#include "args.hpp"
#include "weighted_clusters.hpp"
#include <algorithm>
#include <armadillo>
#include <cmath>
#include <functional>
#include <limits>
#include <ranges>

namespace expectation_maximization {

void weighted_priors_initialization(std::span<double> priors,
                                    RingmapMatrix const &ringmap,
                                    WeightedClusters &weights) noexcept {
  for (auto ringmap_rows = ringmap.rows(); auto &&ringmap_row : ringmap_rows) {
    auto row_counts =
        static_cast<double>(std::size(ringmap_row.modifiedIndices()));
    for (auto &&[prior, weights_cluster] :
         std::views::zip(priors, weights.clusters())) {
      prior = std::ranges::fold_left(
          ringmap_row.modifiedIndices() |
              std::views::transform([&](auto base_index) {
                return weights_cluster[base_index] * row_counts;
              }),
          prior, std::plus<>{});
    }
  }
  auto priors_sum = std::ranges::fold_left(priors, 0., std::plus<>{});
  for (auto &prior : priors) {
    prior /= priors_sum;
  }
}

double ExpectationMaximization::expectation() noexcept {
  auto log_likelihood = 0.;
  auto max_log_likelihood = -std::numeric_limits<double>::infinity();
  for (auto &&[row, responsibilities_row, row_index] : std::views::zip(
           ringmap_->rows(), responsibilities_.rows(), std::views::iota(0))) {
    for (auto &&[prior, cluster_weights, responsibility, cluster_index] :
         std::views::zip(priors_, weights_->clusters(), responsibilities_row,
                         std::views::iota(0))) {
      auto indices_iter = std::ranges::begin(row.modifiedIndices());
      responsibility = std::ranges::fold_left(
          std::views::zip(std::views::iota(0uz), cluster_weights) |
              std::views::drop(row.begin_index()) |
              std::views::take(row.end_index() - row.begin_index()) |
              std::views::transform([&](auto &&tuple) {
                auto &&[base_index, base_weight] = tuple;
                if (indices_iter != std::ranges::end(row.modifiedIndices()) and
                    base_index == *indices_iter) {
                  ++indices_iter;
                  return std::log(base_weight + 1e-10);
                } else {
                  return std::log((1. - base_weight) + 1e-10);
                }
              }),
          std::log(prior + 1e-10), std::plus{});

      max_log_likelihood = std::max(max_log_likelihood, responsibility);
    }

    double responsibilities_row_sum = 0.;
    for (auto &responsibility : responsibilities_row) {
      responsibility = std::exp(responsibility - max_log_likelihood);
      responsibilities_row_sum += responsibility;
    }
    double normalizer = 1. / responsibilities_row_sum;
    for (auto &responsibility : responsibilities_row) {
      responsibility *= normalizer;
    }

    log_likelihood += max_log_likelihood + std::log(responsibilities_row_sum);
  }

  return log_likelihood;
}

void ExpectationMaximization::maximization() noexcept {
  std::ranges::fill(priors_, 0.);
  std::ranges::fill(weights_buffer_.raw_data(), 0.);
  std::ranges::for_each(
      std::views::zip(ringmap_->rows(), responsibilities_.rows()),
      [&](auto &&tuple) {
        auto &&[ringmap_row, responsibilities_row] = tuple;

        for (auto &&[prior, responsibility, weights_row] : std::views::zip(
                 priors_, responsibilities_row, weights_buffer_.clusters())) {
          prior += responsibility;

          for (auto indices = ringmap_row.modifiedIndices();
               auto modification_index : indices) {
            weights_row[modification_index] += responsibility;
          }
        }
      });

  for (auto weights_buffer_clusters = weights_buffer_.clusters();
       auto &&[weights_buffer_cluster, weights_cluster,
               cluster_responsibility_count] :
       std::views::zip(weights_buffer_clusters, weights_->clusters(),
                       priors_)) {
    if (cluster_responsibility_count < 1e-12) {
      std::ranges::fill(weights_cluster, 0.5f);
    } else {
      auto denominator = cluster_responsibility_count + 2e-6;
      for (auto &&[tmp_weight, weight] :
           std::views::zip(*weights_buffer_cluster, weights_cluster)) {
        auto new_weight = (tmp_weight + 1e-6) / denominator;
        weight = std::clamp(static_cast<float>(new_weight), 1e-6f, 1.f - 1e-6f);
      }
    }
  }

  auto n_rows = static_cast<double>(ringmap_->rows_size());
  std::ranges::transform(priors_, std::ranges::begin(priors_),
                         [&](auto prior) { return prior / n_rows; });
}

Result ExpectationMaximization::run() noexcept {
  auto const max_iterations = args_->expectation_maximization_max_iterations();
  auto const tolerance = args_->expectation_maximization_tolerance();

  if (max_iterations == 0) {
    return Result{
        .log_likelihood = -std::numeric_limits<double>::infinity(),
        .convergence = MaxIterations{},
    };
  }

  double previous_log_likelihood = expectation();
  maximization();
  for (std::uint16_t iteration = 1; iteration < max_iterations; ++iteration) {
    auto log_likelihood = expectation();
    if (std::abs(log_likelihood - previous_log_likelihood) < tolerance) {
      return Result{
          .log_likelihood = log_likelihood,
          .convergence = Converged{
              .after_iterations = static_cast<std::uint16_t>(iteration + 1),
          }};
    }

    maximization();
    previous_log_likelihood = log_likelihood;
  }

  return Result{
      .log_likelihood = previous_log_likelihood,
      .convergence = MaxIterations{},
  };
}

} // namespace expectation_maximization

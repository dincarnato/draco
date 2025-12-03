#include "clusters_replicates.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace clusters_replicates {

void reorder_best_permutation(
    std::vector<WeightedClusters> &replicates_clusters) {
  std::size_t const n_replicates = std::size(replicates_clusters);
  if (n_replicates <= 1) {
    return;
  }

  PermutationsDistances permutation_distances(replicates_clusters);
  ReplicatesClustersPermutations replicates_clusters_permutations(
      permutation_distances.replicates(), permutation_distances.clusters());

  double best_clusters_distance = std::numeric_limits<double>::infinity();
  auto best_replicates_clusters_permutations = replicates_clusters_permutations;
  for (;;) {
    // TODO: maybe we should cache part of the calculation
    double total_clusters_distance = 0.;
    for (std::uint8_t replicate_1 = 0;
         replicate_1 < permutation_distances.replicates(); ++replicate_1) {
      auto replicate_1_clusters_permutations =
          replicates_clusters_permutations[replicate_1];

      for (std::uint8_t replicate_2 = replicate_1 + 1;
           replicate_2 < permutation_distances.replicates(); ++replicate_2) {
        auto replicate_2_clusters_permutations =
            replicates_clusters_permutations[replicate_2];
        auto replicates_distances =
            permutation_distances.replicates_pair_distances(replicate_1,
                                                            replicate_2);

        auto current_clusters_distance = replicates_distances.clusters_distance(
            replicate_1_clusters_permutations,
            replicate_2_clusters_permutations);
        total_clusters_distance += current_clusters_distance;
      }
    }

    if (total_clusters_distance < best_clusters_distance) {
      best_clusters_distance = total_clusters_distance;
      std::ranges::copy(
          replicates_clusters_permutations.indices(),
          std::ranges::begin(best_replicates_clusters_permutations.indices()));
    }

    std::uint8_t replicate_index = permutation_distances.replicates() - 1;
    for (;;) {
      if (std::ranges::next_permutation(
              replicates_clusters_permutations[replicate_index])
              .found) {
        break;
      }

      // We don't need to create different combinations for the first replicate
      if (replicate_index <= 1) {
        goto finished;
      } else {
        --replicate_index;
      }
    }
  }
finished:

  assert(not std::isinf(best_clusters_distance));
  std::ranges::for_each(
      std::views::zip(best_replicates_clusters_permutations,
                      replicates_clusters) |
          std::views::drop(1),
      [&](auto &&tuple) {
        auto &&[replicate_permutations, replicate_weighted_clusters] = tuple;

        for (std::uint8_t cluster_index = 0;
             cluster_index < permutation_distances.clusters() - 1;
             ++cluster_index) {

          auto &current_cluster_index = replicate_permutations[cluster_index];
          if (cluster_index != current_cluster_index) {
            auto swapping_cluster_iter = std::ranges::find(
                replicate_permutations | std::views::drop(cluster_index + 1),
                cluster_index);
            assert(swapping_cluster_iter !=
                   std::ranges::end(replicate_permutations));
            replicate_weighted_clusters.swap_clusters(cluster_index,
                                                      current_cluster_index);
            std::swap(current_cluster_index, *swapping_cluster_iter);
          }
        }
        assert(std::ranges::is_sorted(replicate_permutations));
      });
}

PermutationsDistances::PermutationsDistances(
    std::vector<WeightedClusters> &replicates_clusters)
    : replicates_(([&] {
        auto replicates = std::size(replicates_clusters);
        if (replicates > std::numeric_limits<std::uint8_t>::max()) {
          throw std::runtime_error("number of replicates is too high");
        }
        return static_cast<std::uint8_t>(replicates);
      })()),
      clusters_(([&] {
        if (replicates_clusters.empty()) {
          return static_cast<std::uint8_t>(0);
        }

        auto clusters = replicates_clusters[0].getClustersSize();
        if (clusters > std::numeric_limits<std::uint8_t>::max()) {
          throw std::runtime_error("number of clusters is too high");
        }
        return static_cast<std::uint8_t>(clusters);
      })()),
      replicates_combinations_(
          static_cast<std::uint16_t>(replicates_ * (replicates_ - 1) / 2)),
      distances_(distances_size(replicates_, clusters_)) {
  std::size_t distance_index = 0;
  for (std::uint8_t replicate_a_index = 0; replicate_a_index < replicates_;
       ++replicate_a_index) {
    auto const &replicate_a = replicates_clusters[replicate_a_index];

    for (std::uint8_t replicate_b_index = replicate_a_index + 1;
         replicate_b_index < replicates_; ++replicate_b_index) {
      auto const &replicate_b = replicates_clusters[replicate_b_index];

      for (auto replicate_a_clusters = replicate_a.clusters();
           auto &&cluster_a : replicate_a_clusters) {
        for (auto replicate_b_clusters = replicate_b.clusters();
             auto &&cluster_b : replicate_b_clusters) {
          auto distance = std::ranges::fold_left(
              std::views::zip(cluster_a, cluster_b) |
                  std::views::transform([](auto &&weights) {
                    return std::abs(std::get<0>(weights) -
                                    std::get<1>(weights));
                  }),
              0., std::plus{});
          distances_[distance_index++] = distance;
        }
      }
    }
  }
}

} // namespace clusters_replicates

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
      replicates_combinations(
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

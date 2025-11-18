#pragma once

#include "weighted_clusters.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>

namespace clusters_replicates {

constexpr std::size_t distances_size(std::uint8_t replicates,
                                     std::uint8_t clusters) noexcept {
  if (replicates <= 1) {
    return 0;
  }

  std::size_t replicates_combinations =
      static_cast<std::size_t>(replicates) *
      static_cast<std::size_t>(replicates - 1) / 2;
  return replicates_combinations * clusters * clusters;
}

struct PermutationsDistances {
  PermutationsDistances(std::vector<WeightedClusters> &replicates_clusters);

private:
  std::uint8_t replicates_;
  std::uint8_t clusters_;
  std::uint16_t replicates_combinations;
  std::vector<double> distances_;
};

} // namespace clusters_replicates

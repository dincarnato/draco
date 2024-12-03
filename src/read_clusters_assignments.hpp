#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

struct ReadClustersAssignments {
  ReadClustersAssignments(std::size_t clusters) : clusters_(clusters) {};

  std::vector<std::uint32_t> &cluster(std::uint8_t cluster_index) noexcept {
    return clusters_[cluster_index];
  }

  std::vector<std::uint32_t> const &
  cluster(std::uint8_t cluster_index) const noexcept {
    return clusters_[cluster_index];
  }

  std::vector<std::vector<std::uint32_t>> clusters() const noexcept {
    return clusters_;
  }

protected:
  std::vector<std::vector<std::uint32_t>> clusters_;
};

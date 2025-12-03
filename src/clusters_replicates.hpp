#pragma once

#include "results/transcript.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <numeric>
#include <ranges>
#include <span>
#include <vector>

namespace clusters_replicates {

void reorder_best_permutation(
    std::vector<WeightedClusters> &replicates_clusters,
    results::Transcript const &transcript, unsigned window_index);

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

struct ReplicatesPairDistances {
  constexpr ReplicatesPairDistances(std::span<double const> distances,
                                    std::uint8_t clusters) noexcept
      : distances_(distances), clusters_(clusters) {
    assert(clusters * clusters == std::size(distances));
  }

  constexpr double clusters_distance(
      std::span<const std::uint8_t> replicate_a_coupling_indices,
      std::span<const std::uint8_t> replicate_b_coupling_indices) noexcept {
    assert(std::size(replicate_a_coupling_indices) ==
           std::size(replicate_b_coupling_indices));
    assert(std::size(replicate_a_coupling_indices) == clusters_);
    assert(std::ranges::all_of(replicate_a_coupling_indices,
                               [&](auto index) { return index < clusters_; }));
    assert(std::ranges::all_of(replicate_b_coupling_indices,
                               [&](auto index) { return index < clusters_; }));

    return std::ranges::fold_left(
        std::views::zip(replicate_a_coupling_indices,
                        replicate_b_coupling_indices),
        0., [&](double acc, auto &&tuple) {
          auto &&[cluster_a_index, cluster_b_index] = tuple;
          auto linear_index = static_cast<std::size_t>(cluster_a_index) *
                                  static_cast<std::size_t>(clusters_) +
                              static_cast<std::size_t>(cluster_b_index);
          assert(linear_index < std::size(distances_));
          return acc + distances_[linear_index];
        });
  }

private:
  std::span<double const> distances_;
  std::uint8_t clusters_;
};

struct PermutationsDistances {
  PermutationsDistances(std::vector<WeightedClusters> &replicates_clusters);

  constexpr std::span<double const> distances() const noexcept {
    return distances_;
  }
  constexpr std::uint8_t replicates() const noexcept { return replicates_; }
  constexpr std::uint16_t replicates_combinationns() const noexcept {
    return replicates_combinations_;
  }
  constexpr std::uint8_t clusters() const noexcept { return clusters_; }

  constexpr auto
  replicate_linear_index(std::uint8_t replicate_1_index,
                         std::uint8_t replicate_2_index) const noexcept {
    return static_cast<std::size_t>(
        replicates_combinations_ -
        (static_cast<std::uint16_t>(replicates_ - replicate_1_index) *
         static_cast<std::uint16_t>(replicates_ - replicate_1_index - 1) / 2) +
        static_cast<std::uint16_t>(replicate_2_index - replicate_1_index - 1));
  }

  constexpr auto
  replicates_pair_distances(std::uint8_t replicate_1_index,
                            std::uint8_t replicate_2_index) const {
    auto clusters_squared =
        static_cast<std::size_t>(static_cast<std::uint16_t>(clusters_) *
                                 static_cast<std::uint16_t>(clusters_));
    auto begin_index =
        replicate_linear_index(replicate_1_index, replicate_2_index) *
        clusters_squared;
    assert(begin_index < std::size(distances_));
    assert(begin_index + clusters_squared <= std::size(distances_));

    auto replicates_distances =
        std::span(std::ranges::next(std::ranges::begin(distances_),
                                    static_cast<std::ptrdiff_t>(begin_index)),
                  clusters_squared);
    return ReplicatesPairDistances(replicates_distances, clusters_);
  }

private:
  std::uint8_t replicates_;
  std::uint8_t clusters_;
  std::uint16_t replicates_combinations_;
  std::vector<double> distances_;
};

struct ReplicatesClustersPermutations {
  constexpr ReplicatesClustersPermutations(std::uint8_t replicates,
                                           std::uint8_t clusters)
      : indices_(static_cast<std::size_t>(replicates) *
                 static_cast<std::size_t>(clusters)),
        clusters_(clusters) {
    auto chunks = indices_ | std::views::chunk(clusters);
    for (auto &&chunk : chunks) {
      std::ranges::iota(chunk, static_cast<std::uint8_t>(0));
    }
  }

  constexpr std::span<std::uint8_t>
  operator[](std::size_t replicate_index) noexcept {
    return std::span(
        std::ranges::next(
            std::ranges::begin(indices_),
            static_cast<std::ptrdiff_t>(replicate_index *
                                        static_cast<std::size_t>(clusters_))),
        static_cast<std::size_t>(clusters_));
  }

  constexpr std::span<std::uint8_t const> indices() const noexcept {
    return indices_;
  }
  constexpr std::span<std::uint8_t> indices() noexcept { return indices_; }

  constexpr auto replicates(this auto &self) noexcept {
    return self.indices_ | std::views::chunk(self.clusters_);
  }

  constexpr auto clusters(this auto &self) noexcept {
    return std::views::iota(static_cast<std::uint8_t>(0), self.clusters_) |
           std::views::transform([&](auto cluster_index) {
             return self.indices_ | std::views::drop(cluster_index) |
                    std::views::stride(self.clusters_);
           });
  }

  constexpr auto begin() noexcept { return std::ranges::begin(replicates()); }
  constexpr auto end() noexcept { return std::ranges::end(replicates()); }

private:
  std::vector<std::uint8_t> indices_;
  std::uint8_t clusters_;
};

} // namespace clusters_replicates

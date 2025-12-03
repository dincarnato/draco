#include "clusters_replicates.hpp"
#include "results/transcript.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <numeric>
#include <ranges>
#include <vector>

using namespace clusters_replicates;

static void test_distances_size() {
  assert(distances_size(1, 5) == 0);
  assert(distances_size(2, 5) == 25);
  assert(distances_size(4, 5) == 25 * 6);
  assert(distances_size(6, 4) == 16 * 15);
}

static auto make_coupling_indices(auto... index) {
  return std::array{static_cast<std::uint8_t>(index)...};
};

static double calc_clusters_distance(auto const &cluster_a,
                                     auto const &cluster_b) {
  return std::ranges::fold_left(std::views::zip(cluster_a, cluster_b) |
                                    std::views::transform([](auto &&weights) {
                                      return std::abs(std::get<0>(weights) -
                                                      std::get<1>(weights));
                                    }),
                                0., std::plus{});
}

static double calc_replicates_clusters_distance(
    WeightedClusters const &replicate_a,
    std::span<const std::uint8_t> replicate_a_permutations,
    WeightedClusters const &replicate_b,
    std::span<const std::uint8_t> replicate_b_permutations) noexcept {
  return std::ranges::fold_left(
      std::views::zip(replicate_a_permutations, replicate_b_permutations) |
          std::views::transform([&](auto &&tuple) {
            auto [cluster_a_index, cluster_b_index] = tuple;
            auto &&cluster_a = replicate_a.cluster(cluster_a_index);
            auto &&cluster_b = replicate_b.cluster(cluster_b_index);

            return calc_clusters_distance(cluster_a, cluster_b);
          }),
      0., std::plus{});
}

static void test_permutations_distances_constructor() {
  std::array<std::initializer_list<float>, 4> replicate_1{{
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.4f, 0.2f, 0.4f, 0.2f, 0.1f},
      {0.3f, 0.5f, 0.1f, 0.4f, 0.2f},
      {0.2f, 0.1f, 0.3f, 0.1f, 0.3f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_2 = {{
      {0.3f, 0.4f, 0.1f, 0.2f, 0.2f},
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.4f, 0.3f, 0.6f, 0.2f, 0.1f},
      {0.2f, 0.1f, 0.1f, 0.3f, 0.3f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_3 = {{
      {0.5f, 0.3f, 0.3f, 0.2f, 0.3f},
      {0.2f, 0.4f, 0.1f, 0.3f, 0.2f},
      {0.2f, 0.1f, 0.3f, 0.3f, 0.1f},
      {0.1f, 0.2f, 0.3f, 0.2f, 0.4f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_4 = {{
      {0.5f, 0.3f, 0.3f, 0.2f, 0.3f},
      {0.2f, 0.1f, 0.3f, 0.3f, 0.1f},
      {0.1f, 0.2f, 0.3f, 0.2f, 0.4f},
      {0.2f, 0.4f, 0.1f, 0.3f, 0.2f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_5{{
      {0.2f, 0.1f, 0.3f, 0.1f, 0.3f},
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.3f, 0.5f, 0.1f, 0.4f, 0.2f},
      {0.4f, 0.2f, 0.4f, 0.2f, 0.1f},
  }};

  std::vector<WeightedClusters> replicates_clusters{
      WeightedClusters{
          replicate_1[0],
          replicate_1[1],
          replicate_1[2],
          replicate_1[3],
      },
      WeightedClusters{
          replicate_2[0],
          replicate_2[1],
          replicate_2[2],
          replicate_2[3],
      },
      WeightedClusters{
          replicate_3[0],
          replicate_3[1],
          replicate_3[2],
          replicate_3[3],
      },
      WeightedClusters{
          replicate_4[0],
          replicate_4[1],
          replicate_4[2],
          replicate_4[3],
      },
      WeightedClusters{
          replicate_5[0],
          replicate_5[1],
          replicate_5[2],
          replicate_5[3],
      },
  };

  PermutationsDistances permutations_distances(replicates_clusters);
  assert(permutations_distances.clusters() == 4);
  assert(permutations_distances.replicates() == 5);
  auto distances = permutations_distances.distances();
  assert(std::size(distances) == distances_size(5, 4));

  std::size_t distance_index = 0;
  for (std::size_t replicate_a_index = 0; replicate_a_index < 5;
       ++replicate_a_index) {
    auto &replicate_a = replicates_clusters[replicate_a_index];
    for (std::size_t replicate_b_index = replicate_a_index + 1;
         replicate_b_index < 5; ++replicate_b_index) {
      auto &replicate_b = replicates_clusters[replicate_b_index];
      for (std::size_t cluster_a_index = 0; cluster_a_index < 4;
           ++cluster_a_index) {
        for (std::size_t cluster_b_index = 0; cluster_b_index < 4;
             ++cluster_b_index, ++distance_index) {
          assert(distances[distance_index] ==
                 calc_clusters_distance(replicate_a.cluster(cluster_a_index),
                                        replicate_b.cluster(cluster_b_index)));
        }
      }
    }
  }
}

static void test_replicates_pair_distances() {
  std::array<std::initializer_list<float>, 4> replicate_1{{
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.4f, 0.2f, 0.4f, 0.2f, 0.1f},
      {0.3f, 0.5f, 0.1f, 0.4f, 0.2f},
      {0.2f, 0.1f, 0.3f, 0.1f, 0.3f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_2 = {{
      {0.3f, 0.4f, 0.1f, 0.2f, 0.2f},
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.4f, 0.3f, 0.6f, 0.2f, 0.1f},
      {0.2f, 0.1f, 0.1f, 0.3f, 0.3f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_3 = {{
      {0.5f, 0.3f, 0.3f, 0.2f, 0.3f},
      {0.2f, 0.4f, 0.1f, 0.3f, 0.2f},
      {0.2f, 0.1f, 0.3f, 0.3f, 0.1f},
      {0.1f, 0.2f, 0.3f, 0.2f, 0.4f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_4 = {{
      {0.5f, 0.3f, 0.3f, 0.2f, 0.3f},
      {0.2f, 0.1f, 0.3f, 0.3f, 0.1f},
      {0.1f, 0.2f, 0.3f, 0.2f, 0.4f},
      {0.2f, 0.4f, 0.1f, 0.3f, 0.2f},
  }};

  std::array<std::initializer_list<float>, 4> replicate_5{{
      {0.2f, 0.1f, 0.3f, 0.1f, 0.3f},
      {0.1f, 0.2f, 0.2f, 0.3f, 0.4f},
      {0.3f, 0.5f, 0.1f, 0.4f, 0.2f},
      {0.4f, 0.2f, 0.4f, 0.2f, 0.1f},
  }};

  std::vector<WeightedClusters> replicates_clusters{
      WeightedClusters{
          replicate_1[0],
          replicate_1[1],
          replicate_1[2],
          replicate_1[3],
      },
      WeightedClusters{
          replicate_2[0],
          replicate_2[1],
          replicate_2[2],
          replicate_2[3],
      },
      WeightedClusters{
          replicate_3[0],
          replicate_3[1],
          replicate_3[2],
          replicate_3[3],
      },
      WeightedClusters{
          replicate_4[0],
          replicate_4[1],
          replicate_4[2],
          replicate_4[3],
      },
      WeightedClusters{
          replicate_5[0],
          replicate_5[1],
          replicate_5[2],
          replicate_5[3],
      },
  };

  PermutationsDistances permutations_distances(replicates_clusters);

  auto assert_replicate_distances = [&](std::uint8_t replicate_a_index,
                                        std::uint8_t replicate_b_index) {
    auto replicates_distances =
        permutations_distances.replicates_pair_distances(replicate_a_index,
                                                         replicate_b_index);

    std::array<std::uint8_t, 4> replicate_a_permutations;
    std::array<std::uint8_t, 4> replicate_b_permutations;

    std::ranges::iota(replicate_a_permutations, static_cast<std::uint8_t>(0));
    do {
      std::ranges::iota(replicate_b_permutations, static_cast<std::uint8_t>(0));
      do {
        assert(replicates_distances.clusters_distance(
                   replicate_a_permutations, replicate_b_permutations) ==
               calc_replicates_clusters_distance(
                   replicates_clusters[replicate_a_index],
                   replicate_a_permutations,
                   replicates_clusters[replicate_b_index],
                   replicate_b_permutations));
      } while (std::ranges::next_permutation(replicate_b_permutations).found);
    } while (std::ranges::next_permutation(replicate_a_permutations).found);
  };

  assert_replicate_distances(0, 1);
  assert_replicate_distances(0, 3);
  assert_replicate_distances(1, 3);
}

static void test_reorder_best_permutation_more_iterations() {
  std::array<std::initializer_list<float>, 5> replicate_1{{
      {0.10f, 0.10f, 0.25f, 0.35f, 0.30f},
      {0.10f, 0.30f, 0.35f, 0.20f, 0.25f},
      {0.30f, 0.40f, 0.10f, 0.35f, 0.15f},
      {0.15f, 0.05f, 0.20f, 0.05f, 0.10f},
      {0.05f, 0.15f, 0.10f, 0.05f, 0.20f},
  }};

  std::array<std::initializer_list<float>, 5> replicate_2{{
      {0.30f, 0.40f, 0.10f, 0.35f, 0.15f},
      {0.10f, 0.30f, 0.35f, 0.20f, 0.25f},
      {0.05f, 0.15f, 0.10f, 0.05f, 0.20f},
      {0.15f, 0.05f, 0.20f, 0.05f, 0.10f},
      {0.10f, 0.10f, 0.25f, 0.35f, 0.30f},
  }};

  std::array<std::initializer_list<float>, 5> replicate_3{{
      {0.15f, 0.05f, 0.20f, 0.05f, 0.10f},
      {0.05f, 0.15f, 0.10f, 0.05f, 0.20f},
      {0.30f, 0.40f, 0.10f, 0.35f, 0.15f},
      {0.10f, 0.10f, 0.25f, 0.35f, 0.30f},
      {0.10f, 0.30f, 0.35f, 0.20f, 0.25f},
  }};

  std::array<std::initializer_list<float>, 5> replicate_4{{
      {0.10f, 0.10f, 0.25f, 0.35f, 0.30f},
      {0.30f, 0.40f, 0.10f, 0.35f, 0.15f},
      {0.10f, 0.30f, 0.35f, 0.20f, 0.25f},
      {0.15f, 0.05f, 0.20f, 0.05f, 0.10f},
      {0.05f, 0.15f, 0.10f, 0.05f, 0.20f},
  }};

  std::vector<WeightedClusters> replicates_weighted_clusters{
      WeightedClusters{
          replicate_1[0],
          replicate_1[1],
          replicate_1[2],
          replicate_1[3],
          replicate_1[4],
      },
      WeightedClusters{
          replicate_2[0],
          replicate_2[1],
          replicate_2[2],
          replicate_2[3],
          replicate_2[4],
      },
      WeightedClusters{
          replicate_3[0],
          replicate_3[1],
          replicate_3[2],
          replicate_3[3],
          replicate_3[4],
      },
      WeightedClusters{
          replicate_4[0],
          replicate_4[1],
          replicate_4[2],
          replicate_4[3],
          replicate_4[4],
      },
  };

  results::Transcript transcript(5);
  reorder_best_permutation(replicates_weighted_clusters, transcript, 0);
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(0),
                            replicate_1[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(1),
                            replicate_1[1]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(2),
                            replicate_1[2]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(3),
                            replicate_1[3]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(4),
                            replicate_1[4]));

  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(0),
                            replicate_2[4]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(1),
                            replicate_2[1]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(2),
                            replicate_2[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(3),
                            replicate_2[3]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(4),
                            replicate_2[2]));

  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(0),
                            replicate_3[3]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(1),
                            replicate_3[4]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(2),
                            replicate_3[2]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(3),
                            replicate_3[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(4),
                            replicate_3[1]));

  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(0),
                            replicate_4[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(1),
                            replicate_4[2]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(2),
                            replicate_4[1]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(3),
                            replicate_4[3]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(4),
                            replicate_4[4]));
}

static void test_reorder_best_permutation() {
  std::array<std::initializer_list<float>, 3> replicate_1{{
      {0.10f, 0.20f, 0.30f, 0.40f, 0.50f},
      {0.60f, 0.30f, 0.60f, 0.20f, 0.30f},
      {0.30f, 0.50f, 0.10f, 0.40f, 0.20f},
  }};

  std::array<std::initializer_list<float>, 3> replicate_2 = {
      {{0.33f, 0.50f, 0.13f, 0.40f, 0.20f},
       {0.10f, 0.20f, 0.27f, 0.40f, 0.47f},
       {0.57f, 0.30f, 0.60f, 0.20f, 0.34f}}};

  std::array<std::initializer_list<float>, 3> replicate_3 = {
      {{0.60f, 0.30f, 0.60f, 0.24f, 0.30f},
       {0.27f, 0.50f, 0.10f, 0.37f, 0.20f},
       {0.13f, 0.20f, 0.30f, 0.40f, 0.50f}}};

  std::vector<WeightedClusters> replicates_weighted_clusters{
      WeightedClusters{
          replicate_1[0],
          replicate_1[1],
          replicate_1[2],
      },
      WeightedClusters{
          replicate_2[0],
          replicate_2[1],
          replicate_2[2],
      },
      WeightedClusters{
          replicate_3[0],
          replicate_3[1],
          replicate_3[2],
      },
  };

  results::Transcript transcript(3);
  reorder_best_permutation(replicates_weighted_clusters, transcript, 3);
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(0),
                            replicate_1[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(1),
                            replicate_1[1]));
  assert(std::ranges::equal(replicates_weighted_clusters[0].cluster(2),
                            replicate_1[2]));

  assert(std::ranges::equal(replicates_weighted_clusters[1].cluster(0),
                            replicate_2[1]));
  assert(std::ranges::equal(replicates_weighted_clusters[1].cluster(1),
                            replicate_2[2]));
  assert(std::ranges::equal(replicates_weighted_clusters[1].cluster(2),
                            replicate_2[0]));

  assert(std::ranges::equal(replicates_weighted_clusters[2].cluster(0),
                            replicate_3[2]));
  assert(std::ranges::equal(replicates_weighted_clusters[2].cluster(1),
                            replicate_3[0]));
  assert(std::ranges::equal(replicates_weighted_clusters[2].cluster(2),
                            replicate_3[1]));
}

int main() {
  test_distances_size();
  test_permutations_distances_constructor();
  test_replicates_pair_distances();
  test_reorder_best_permutation();
  test_reorder_best_permutation_more_iterations();
}

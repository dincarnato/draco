#include "clusters_replicates.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
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

int main() {
  test_distances_size();
  test_permutations_distances_constructor();
}

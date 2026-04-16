#include "test/expectation_maximization.hpp"
#include "compact_ringmap.hpp"
#include "expectation_maximization.hpp"
#include "mutation_map_transcript_read.hpp"
#include "ringmap_matrix.hpp"
#include "test/args.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <random>
#include <span>
#include <string_view>
#include <vector>

static CompactRingmap
make_compact_ringmap(std::size_t n_bases,
                     std::span<MutationMapTranscriptRead> const &reads) {
  RingmapMatrix matrix(static_cast<std::uint32_t>(reads.size()),
                       static_cast<unsigned>(n_bases));
  for (auto const &read : reads)
    matrix.addRead(read);
  return CompactRingmap(matrix);
}

static void test_single_cluster_prior_is_one() {
  constexpr std::size_t n_bases = 4;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 4u, .indices = {1, 3}},
      {.begin = 0u, .end = 4u, .indices = {0}},
      {.begin = 0u, .end = 4u, .indices = {1, 3}},
  };
  auto ringmap = make_compact_ringmap(n_bases, reads);

  WeightedClusters weights({{0.3f, 0.9f, 0.5f, 0.7f}});

  std::array<double, 1> priors{0.0};
  expectation_maximization::weighted_priors_initialization(priors, ringmap,
                                                           weights);

  assert(std::abs(priors[0] - 1.0) < 1e-12);
}

static void test_symmetric_weights_produce_uniform_priors() {
  constexpr std::size_t n_bases = 4;
  constexpr std::size_t n_clusters = 3;
  constexpr float each_base_weight = 1.f / static_cast<float>(n_clusters);

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 4u, .indices = {0, 2}},
      {.begin = 0u, .end = 4u, .indices = {1}},
      {.begin = 0u, .end = 4u, .indices = {0, 2}},
  };
  auto ringmap = make_compact_ringmap(n_bases, reads);

  // All clusters, all bases: same weight w.
  WeightedClusters weights(n_bases, n_clusters, false);
  for (auto &&weights_element : weights) {
    std::ranges::fill(weights_element, each_base_weight);
  }

  std::vector<double> priors(n_clusters, 0.);
  expectation_maximization::weighted_priors_initialization(priors, ringmap,
                                                           weights);

  assert(std::ranges::all_of(priors, [](auto prior) {
    return std::abs(prior - 1.0 / static_cast<double>(n_clusters)) < 1e-12;
  }));
}

static void test_exact_values_two_clusters() {
  constexpr std::size_t n_bases = 2;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 2u, .indices = {0}},
      {.begin = 0u, .end = 2u, .indices = {0}},
      {.begin = 0u, .end = 2u, .indices = {1}},
  };
  auto ringmap = make_compact_ringmap(n_bases, reads);

  WeightedClusters weights({
      {0.8f, 0.2f},
      {0.2f, 0.8f},
  });

  std::array<double, 2> priors{0., 0.};
  expectation_maximization::weighted_priors_initialization(priors, ringmap,
                                                           weights);

  assert(std::abs(std::ranges::fold_left(priors, 0., std::plus{}) - 1.) <
         1e-12);
  assert(std::abs(priors[0] - 0.6) < 1e-6);
  assert(std::abs(priors[1] - 0.4) < 1e-6);
}

static void test_exact_values_multiple_patterns() {
  constexpr std::size_t n_bases = 4;
  constexpr std::array<float, 2> expected_unnormalized_priors{5.3f, 5.7f};
  constexpr float expected_unnormalized_priors_sum =
      std::ranges::fold_left(expected_unnormalized_priors, 0.f, std::plus<>{});

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 4u, .indices = {0, 1}},
      {.begin = 0u, .end = 4u, .indices = {0, 1}},
      {.begin = 0u, .end = 4u, .indices = {2, 3}},
      {.begin = 0u, .end = 4u, .indices = {2, 3}},
      {.begin = 0u, .end = 4u, .indices = {2, 3}},
      {.begin = 0u, .end = 4u, .indices = {1}},
  };
  auto ringmap = make_compact_ringmap(n_bases, reads);

  WeightedClusters weights({
      {0.7f, 0.3f, 0.6f, 0.4f},
      {0.3f, 0.7f, 0.4f, 0.6f},
  });

  std::array<double, 2> priors{0., 0.};
  expectation_maximization::weighted_priors_initialization(priors, ringmap,
                                                           weights);

  assert(std::abs(std::ranges::fold_left(priors, 0., std::plus{}) - 1.) <
         1e-12);
  assert(std::ranges::all_of(
      std::views::zip(priors, expected_unnormalized_priors), [&](auto &&tuple) {
        auto [prior, expected_unnormalized_prior] = tuple;
        return (prior - expected_unnormalized_prior /
                            expected_unnormalized_priors_sum) < 1e-6;
      }));
}

static void test_mixed_empty_and_modified_reads() {
  constexpr std::size_t n_bases = 2;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 2u, .indices = {0}},
      {.begin = 0u, .end = 2u, .indices = {0}},
      {.begin = 0u, .end = 2u, .indices = {}},
      {.begin = 0u, .end = 2u, .indices = {}},
      {.begin = 0u, .end = 2u, .indices = {}},
  };
  auto ringmap = make_compact_ringmap(n_bases, reads);

  WeightedClusters weights({
      {0.8f, 0.2f},
      {0.2f, 0.8f},
  });

  std::array<double, 2> priors{0., 0.};
  expectation_maximization::weighted_priors_initialization(priors, ringmap,
                                                           weights);

  assert(std::abs(priors[0] - 0.8) < 1e-6);
  assert(std::abs(priors[1] - 0.2) < 1e-6);
}

static std::array<MutationMapTranscriptRead, 13> reads{
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {1, 3, 4}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {1, 3, 4, 5}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {2, 3}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {4}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {0, 2}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {0, 2}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {1, 2}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {6, 7}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {6, 7}},
    MutationMapTranscriptRead{.begin = 0u, .end = 8u, .indices = {0, 6}},
};

void test_steps() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGATTTATATGGACCATA";
  RingmapMatrix matrix(static_cast<std::uint32_t>(std::size(reads)),
                       std::size(sequence));
  for (auto const &read : reads) {
    matrix.addRead(read);
  }

  auto ringmap = CompactRingmap(matrix);
  WeightedClusters weights({{
      {0.1f, 0.6f, 0.4f, 0.5f, 0.4f, 0.3f, 0.1f, 0.1f},
      {0.6f, 0.3f, 0.5f, 0.4f, 0.4f, 0.3f, 0.3f, 0.1f},
      {0.3f, 0.1f, 0.1f, 0.1f, 0.2f, 0.4f, 0.6f, 0.8f},
  }});

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;
  std::mt19937 rng(0);
  test::ExpectationMaximization expectation_maximization(ringmap, weights, args,
                                                         rng);
  {
    auto log_likelihood = expectation_maximization.expectation();
    assert(std::abs(log_likelihood - (-53.01284286)) < 1e-6);
    auto &&responsibilities = expectation_maximization.responsibilities();
    assert(std::ranges::equal(
        responsibilities.raw_data(),
        std::array{
            0.4677, 0.2829, 0.2494, 0.5540, 0.3352, 0.1108, 0.0736,
            0.9012, 0.0252, 0.0166, 0.5226, 0.4608, 0.7900, 0.2048,
            0.0052, 0.6193, 0.3746, 0.0061, 0.0038, 0.0089, 0.9873,
            0.8954, 0.1032, 0.0015, 0.8946, 0.1031, 0.0023,
        },
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));

    expectation_maximization.maximization();
    assert(std::ranges::equal(
        expectation_maximization.priors(),
        std::array{0.40981685, 0.33172291, 0.25846024},
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));
    assert(std::ranges::equal(
        weights.raw(),
        std::array{3.0732050e-02, 5.3915546e-01, 1.5214959e-01, 4.8426310e-01,
                   9.5321760e-02, 2.6688100e-03, 2.9213605e-01, 5.5232527e-01,
                   1.8382650e-02, 4.5221861e-01, 1.3469805e-01, 2.9410300e-03,
                   4.3997478e-01, 1.2554759e-01, 3.4099170e-02, 1.6792221e-01,
                   2.3904400e-02, 6.8208000e-04, 4.5443700e-03, 1.2531589e-01,
                   7.2481871e-01, 1.4301100e-03, 4.1221900e-03, 5.8768337e-01},
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));
  }

  // Second iteration
  {
    auto log_likelihood = expectation_maximization.expectation();
    assert(std::abs(log_likelihood - (-41.86657741)) < 1e-6);
    auto &&responsibilities = expectation_maximization.responsibilities();
    assert(std::ranges::equal(
        responsibilities.raw_data(),
        std::array{
            3.6775281e-01, 3.9811205e-01, 2.3413514e-01, 8.1536598e-01,
            1.6130727e-01, 2.3326750e-02, 8.2933600e-03, 9.9035058e-01,
            1.3560500e-03, 2.9997000e-04, 3.7604271e-01, 6.2365732e-01,
            7.3354844e-01, 2.6639117e-01, 6.0390000e-05, 6.2098666e-01,
            3.7894924e-01, 6.4100000e-05, 2.7300000e-06, 2.6852000e-04,
            9.9972875e-01, 9.9583117e-01, 4.1685400e-03, 2.9000000e-07,
            9.9949228e-01, 5.0771000e-04, 0.0000000e+00,
        },
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));

    expectation_maximization.maximization();
    assert(std::ranges::equal(
        expectation_maximization.priors(),
        std::array{0.40656732, 0.33561085, 0.25782184},
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));
    assert(std::ranges::equal(
        weights.raw(),
        std::array{3.1951700e-03, 5.4017320e-01, 1.8688220e-01, 5.1630618e-01,
                   6.2129720e-02, 1.8400000e-05, 2.5941816e-01, 6.0189706e-01,
                   8.4662000e-04, 4.9500935e-01, 8.7928380e-02, 1.9510000e-05,
                   5.3178616e-01, 3.8044160e-02, 6.9600900e-03, 1.8910539e-01,
                   1.1660000e-04, 1.0000000e-06, 5.7980000e-05, 8.6313470e-02,
                   7.8262579e-01, 1.2200000e-06, 1.2332000e-04, 5.9655307e-01},
        [](auto a, auto b) { return std::abs(a - b) < 1e-3; }));
  }

  // Let's finish running the whole algorithm
  auto result = expectation_maximization.run();
  auto const &converged =
      std::get<expectation_maximization::Converged>(result.convergence);
  assert(converged.after_iterations == 21);
  assert(std::abs(result.log_likelihood - (-37.38889328)) < 1e-6);
  assert(std::ranges::equal(
      expectation_maximization.priors(),
      std::array{0.2306184, 0.53861242, 0.23076918},
      [](auto a, auto b) { return std::abs(a - b) < 1e-6; }));
  assert(std::ranges::equal(
      expectation_maximization.responsibilities().raw_data(),
      std::array{
          1.4000000e-07, 9.9999944e-01, 4.2000000e-07, 9.9800977e-01,
          1.9902300e-03, 0.0000000e+00, 0.0000000e+00, 1.0000000e+00,
          0.0000000e+00, 0.0000000e+00, 1.8900000e-06, 9.9999811e-01,
          0.0000000e+00, 1.0000000e+00, 0.0000000e+00, 0.0000000e+00,
          1.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00,
          1.0000000e+00, 9.9998621e-01, 1.3790000e-05, 0.0000000e+00,
          1.0000000e+00, 0.0000000e+00, 0.0000000e+00,
      },
      [](auto a, auto b) { return std::abs(a - b) < 1e-6; }));
  assert(std::ranges::equal(
      weights.raw(),
      std::array{1.0000000e-06, 2.8563458e-01, 3.3333288e-01, 6.6709808e-01,
                 1.4281915e-01, 1.0000000e-06, 1.0000000e-06, 5.7126848e-01,
                 1.0000000e-06, 6.6709808e-01, 1.4281915e-01, 1.0000000e-06,
                 9.9999900e-01, 2.8024000e-04, 1.0000000e-06, 3.3355146e-01,
                 1.0000000e-06, 1.0000000e-06, 1.0000000e-06, 1.0000000e-06,
                 9.9999900e-01, 1.0000000e-06, 1.0000000e-06, 6.6666669e-01},
      [](auto a, auto b) { return std::abs(a - b) < 1e-6; }));
}

void test_read_assignment() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGATTTATATGGACCATA";
  RingmapMatrix matrix(static_cast<std::uint32_t>(std::size(reads)),
                       std::size(sequence));
  for (auto const &read : reads) {
    for (std::uint8_t index = 0; index < 20; ++index) {
      matrix.addRead(read);
    }
  }

  auto ringmap = CompactRingmap(matrix);
  WeightedClusters weights({{
      {0.1f, 0.7f, 0.4f, 0.5f, 0.4f, 0.3f, 0.1f, 0.1f},
      {0.6f, 0.2f, 0.5f, 0.4f, 0.4f, 0.3f, 0.3f, 0.1f},
      {0.3f, 0.1f, 0.1f, 0.1f, 0.2f, 0.4f, 0.6f, 0.8f},
  }});

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;
  std::mt19937 rng(0);
  test::ExpectationMaximization expectation_maximization(ringmap, weights, args,
                                                         rng);
  std::vector<std::uint32_t> assignments(3);
  std::vector<double> buffer(3);

  auto row = ringmap.row(0);
  assert(row.indices().empty());
  expectation_maximization.read_assignment(row, assignments, buffer, rng);
  assert(std::ranges::equal(assignments, std::array{23u, 21u, 16u}));

  row = ringmap.row(2);
  assert(row.count() == 40);
  expectation_maximization.read_assignment(row, assignments, buffer, rng);
  assert(std::ranges::equal(assignments, std::array{2u, 37u, 1u}));

  row = ringmap.row(6);
  assert(row.count() == 40);
  expectation_maximization.read_assignment(row, assignments, buffer, rng);
  assert(std::ranges::equal(assignments, std::array{0u, 1u, 39u}));

  row = ringmap.row(7);
  assert(row.count() == 20);
  expectation_maximization.read_assignment(row, assignments, buffer, rng);
  assert(std::ranges::equal(assignments, std::array{19u, 1u, 0u}));
}

int main() {
  test_single_cluster_prior_is_one();
  test_symmetric_weights_produce_uniform_priors();
  test_exact_values_two_clusters();
  test_exact_values_multiple_patterns();
  test_mixed_empty_and_modified_reads();
  test_steps();
  test_read_assignment();
}

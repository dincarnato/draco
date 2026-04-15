#include "reassignment.hpp"
#include "draco.hpp"
#include "mutation_map_transcript_read.hpp"
#include "priors_initialization.hpp"
#include "results/transcript.hpp"
#include "results/window.hpp"
#include "ringmap_data.hpp"
#include "ringmap_matrix.hpp"
#include "test/args.hpp"
#include "test/reassignment.hpp"
#include "weighted_clusters.hpp"

#include <cassert>
#include <cstdint>
#include <optional>
#include <random>
#include <vector>

static RingmapData
make_ringmap(std::string_view sequence,
             std::span<MutationMapTranscriptRead const> reads,
             Args const &args) {
  RingmapMatrix matrix(static_cast<std::uint32_t>(reads.size()),
                       static_cast<unsigned>(sequence.size()));
  for (auto const &read : reads) {
    matrix.addRead(read);
  }
  return RingmapData(sequence, std::move(matrix), 0u,
                     static_cast<unsigned>(sequence.size()), args);
}

static results::Window make_window(unsigned short begin_index,
                                   WeightedClusters const &weights,
                                   unsigned window_size) {
  return results::Window(begin_index, weights,
                         std::vector<unsigned>(window_size, 1u));
}

static void test_single_cluster() {
  constexpr std::string_view sequence = "AGCTAATT";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads(
      5, {.begin = 0u, .end = 8u, .indices = {3}});

  auto ringmap = make_ringmap(sequence, reads, args);

  WeightedClusters weights(window_size, 1u);
  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(42);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_fractions_sum_to_one_two_clusters() {
  constexpr std::string_view sequence = "AGCTAATTCCGCCGAT";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 16u, .indices = {1, 3, 4}},
      {.begin = 0u, .end = 16u, .indices = {}},
      {.begin = 0u, .end = 16u, .indices = {1, 3, 4}},
      {.begin = 0u, .end = 16u, .indices = {}},
      {.begin = 0u, .end = 16u, .indices = {10, 11, 12}},
      {.begin = 0u, .end = 16u, .indices = {10, 11, 12}},
      {.begin = 0u, .end = 16u, .indices = {1, 3, 4}},
      {.begin = 0u, .end = 16u, .indices = {10, 11, 12}},
  };

  auto ringmap = make_ringmap(sequence, reads, args);

  WeightedClusters weights({
      {0.85f, 0.875f, 0.5f, 0.875f, 0.875f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
       0.125f, 0.125f, 0.125f, 0.5f, 0.5f, 0.5f},
      {0.15f, 0.125f, 0.5f, 0.125f, 0.125f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f,
       0.875f, 0.875f, 0.875f, 0.5f, 0.5f, 0.5f},
  });

  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(0);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_well_separated_clusters_converge_correctly() {
  constexpr std::string_view sequence = "AGCTAATT";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());
  constexpr auto reads_per_cluster = 10u;

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads;
  for (unsigned read_index = 0; read_index < reads_per_cluster; ++read_index) {
    reads.push_back({.begin = 0u, .end = 8u, .indices = {0, 1}});
  }
  for (unsigned read_index = 0; read_index < reads_per_cluster; ++read_index) {
    reads.push_back({.begin = 0u, .end = 8u, .indices = {6, 7}});
  }

  auto ringmap = make_ringmap(sequence, reads, args);

  WeightedClusters weights({
      {0.9f, 0.9f, 0.5f, 0.5f, 0.5f, 0.5f, 0.1f, 0.1f},
      {0.1f, 0.1f, 0.5f, 0.5f, 0.5f, 0.5f, 0.9f, 0.9f},
  });

  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(12345);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_output_shapes_match_weighted_clusters() {
  constexpr std::string_view sequence = "AGCTAATTCCG";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());
  constexpr auto n_clusters = 3uz;

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 11u, .indices = {0, 1}},
      {.begin = 0u, .end = 11u, .indices = {4, 5}},
      {.begin = 0u, .end = 11u, .indices = {8, 9}},
      {.begin = 0u, .end = 11u, .indices = {0, 1}},
      {.begin = 0u, .end = 11u, .indices = {4, 5}},
      {.begin = 0u, .end = 11u, .indices = {8, 9}},
  };

  auto ringmap = make_ringmap(sequence, reads, args);
  WeightedClusters weights(window_size, n_clusters);
  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(99);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_reads_shorter_than_window_are_skipped_gracefully() {
  constexpr std::string_view sequence = "AGCTAATT";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 8u, .indices = {1, 2}},
      {.begin = 0u, .end = 8u, .indices = {1, 2}},
      {.begin = 0u, .end = 8u, .indices = {5, 6}},
      {.begin = 0u, .end = 8u, .indices = {5, 6}},
      {.begin = 1u, .end = 5u, .indices = {2}},
      {.begin = 2u, .end = 6u, .indices = {3}},
  };

  auto ringmap = make_ringmap(sequence, reads, args);

  WeightedClusters weights({
      {0.875f, 0.875f, 0.5f, 0.5f, 0.5f, 0.125f, 0.125f, 0.5f},
      {0.125f, 0.125f, 0.5f, 0.5f, 0.5f, 0.875f, 0.875f, 0.5f},
  });

  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(7);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_base_filtered_ringmap_create_reduced_path() {
  constexpr std::string_view sequence = "AGCTAATT"; // 8 bases
  constexpr auto window_size = static_cast<unsigned>(sequence.size());
  constexpr auto reads_per_cluster = 1000u;

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads;
  for (unsigned i = 0; i < reads_per_cluster; ++i) {
    reads.push_back({.begin = 0u, .end = 8u, .indices = {2}});
  }
  for (unsigned i = 0; i < reads_per_cluster; ++i) {
    reads.push_back({.begin = 0u, .end = 8u, .indices = {5}});
  }

  auto ringmap = make_ringmap(sequence, reads, args);

  std::vector<RingmapData> filtered{ringmap};
  RingmapData::filter_bases_on_replicates_for_assignments(filtered);

  assert(filtered[0].bases_filtered());
  assert(filtered[0].data().cols_size() == 4);

  // Cluster 0 is biased towards base 2 (filtered col 0);
  // cluster 1 is biased towards base 5 (filtered col 1).
  WeightedClusters weights({
      {0.5f, 0.5f, 0.9f, 0.5f, 0.5f, 0.1f, 0.5f, 0.5f},
      {0.5f, 0.5f, 0.1f, 0.5f, 0.5f, 0.9f, 0.5f, 0.5f},
  });

  results::Window window = make_window(0, weights, window_size);

  std::mt19937 rng(55);
  std::vector<std::uint32_t> assignments_per_cluster;
  std::vector<double> buffer;
  std::vector<std::uint32_t> mapped_rows;
  std::vector<std::uint32_t> clusters_reads_count;

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};

  test::Reassignment({
                         .replicates_splitted_ringmaps = replicates_splitted,
                         .filtered_ringmaps = filtered,
                         .ptba_on_replicate_results = {},
                         .windows_max_clusters_constraints = constraints,
                         .transcript_result = &transcript,
                         .stop = &stop,
                         .args = &args,
                         .window_index = 0uz,
                         .window_size = window_size,
                     })
      .reweight_with_expectation_maximization_iteration(
          reassignment::ReweightAndReassignWithExpectationMaximizationIteration{
              .filtered_ringmap = &filtered[0],
              .window = &window,
              .rng = &rng,
              .assignments_per_cluster = &assignments_per_cluster,
              .buffer = &buffer,
              .mapped_rows = &mapped_rows,
              .clusters_reads_count = &clusters_reads_count,
              .replicate_index = 0uz,
          });
}

static void test_buffer_reuse_across_successive_calls() {
  constexpr std::string_view sequence = "AGCTAATT";
  constexpr auto window_size = static_cast<unsigned>(sequence.size());

  test::Args args;
  args.expectation_maximization_priors_initialization() =
      args::PriorsInitialization::Uniform;

  std::vector<MutationMapTranscriptRead> reads{
      {.begin = 0u, .end = 8u, .indices = {0, 1}},
      {.begin = 0u, .end = 8u, .indices = {0, 1}},
      {.begin = 0u, .end = 8u, .indices = {6, 7}},
      {.begin = 0u, .end = 8u, .indices = {6, 7}},
  };

  auto ringmap = make_ringmap(sequence, reads, args);

  std::mt19937 rng(0);
  // Garbage, the function must reinitialise these on every call.
  std::vector<std::uint32_t> assignments_per_cluster{999u, 999u};
  std::vector<double> buffer{-1.0, -1.0};
  std::vector<std::uint32_t> mapped_rows{42u, 42u, 42u};
  std::vector<std::uint32_t> clusters_reads_count{888u, 888u};

  results::Transcript transcript(1);
  bool stop = true;
  std::vector<std::optional<unsigned>> constraints(1, std::nullopt);
  std::vector<std::vector<RingmapData>> replicates_splitted{{ringmap}};
  std::vector<RingmapData> filtered{ringmap};

  test::Reassignment reassignment({
      .replicates_splitted_ringmaps = replicates_splitted,
      .filtered_ringmaps = filtered,
      .ptba_on_replicate_results = {},
      .windows_max_clusters_constraints = constraints,
      .transcript_result = &transcript,
      .stop = &stop,
      .args = &args,
      .window_index = 0uz,
      .window_size = window_size,
  });

  reassignment::ReweightAndReassignWithExpectationMaximizationIteration
      iteration_args{
          .filtered_ringmap = &filtered[0],
          .window = nullptr,
          .rng = &rng,
          .assignments_per_cluster = &assignments_per_cluster,
          .buffer = &buffer,
          .mapped_rows = &mapped_rows,
          .clusters_reads_count = &clusters_reads_count,
          .replicate_index = 0uz,
      };

  for (std::uint8_t call = 0; call < 3; ++call) {
    WeightedClusters weights({
        {0.9f, 0.9f, 0.5f, 0.5f, 0.5f, 0.5f, 0.1f, 0.1f},
        {0.1f, 0.1f, 0.5f, 0.5f, 0.5f, 0.5f, 0.9f, 0.9f},
    });
    results::Window window = make_window(0, weights, window_size);
    iteration_args.window = &window;

    reassignment.reweight_with_expectation_maximization_iteration(
        iteration_args);
  }
}

int main() {
  test_single_cluster();
  test_fractions_sum_to_one_two_clusters();
  test_well_separated_clusters_converge_correctly();
  test_output_shapes_match_weighted_clusters();
  test_reads_shorter_than_window_are_skipped_gracefully();
  test_base_filtered_ringmap_create_reduced_path();
  test_buffer_reuse_across_successive_calls();
}

#include "draco.hpp"
#include "args.hpp"
#include "results/transcript.hpp"
#include "ringmap_matrix.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <ranges>
#include <vector>

void test_merge_windows_and_add_window_results_not_merging() {
  Args args;
  constexpr std::size_t SEQUENCE_SIZE = 40uz;
  constexpr std::size_t WINDOW_SIZE = 20uz;
  std::vector<Window> windows{
      Window{.start_base = 0,
             .weights = WeightedClusters(WINDOW_SIZE, 1),
             .coverages = std::vector<unsigned>(WINDOW_SIZE)},
      Window{.start_base = 2,
             .weights = WeightedClusters(WINDOW_SIZE, 2),
             .coverages = std::vector<unsigned>(WINDOW_SIZE)},
      Window{.start_base = 4,
             .weights = WeightedClusters(WINDOW_SIZE, 1),
             .coverages = std::vector<unsigned>(WINDOW_SIZE)},
      Window{.start_base = 6,
             .weights = WeightedClusters(WINDOW_SIZE, 2),
             .coverages = std::vector<unsigned>(WINDOW_SIZE)},
  };
  constexpr unsigned WINDOWS = 4u;
  constexpr unsigned READS_PER_WINDOW = 10u;
  auto windows_reads_indices =
      std::views::iota(0u, WINDOWS) |
      std::views::transform([&](auto window_index) {
        return std::views::iota(0u, READS_PER_WINDOW) |
               std::views::transform([&](auto window_read_index) {
                 return window_index * READS_PER_WINDOW + window_read_index;
               }) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();
  auto sequence = std::views::repeat('A') | std::views::take(SEQUENCE_SIZE) |
                  std::ranges::to<std::string>();
  RingmapMatrix ringmap_matrix(WINDOWS * READS_PER_WINDOW, SEQUENCE_SIZE);
  RingmapData ringmap_data(sequence, std::move(ringmap_matrix), 0,
                           SEQUENCE_SIZE, args);

  results::Transcript transcript_results;
  merge_windows_and_add_window_results(windows, windows_reads_indices,
                                       ringmap_data, transcript_results, args);

  assert(transcript_results.windows.has_value());
  assert(std::size(*transcript_results.windows) == WINDOWS);
  assert(std::ranges::equal(
      *transcript_results.windows |
          std::views::transform([](auto const &window) {
            return window.weighted_clusters.getClustersSize();
          }),
      std::array{1, 2, 1, 2}));
}

int main() { test_merge_windows_and_add_window_results_not_merging(); }

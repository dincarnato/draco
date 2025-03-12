#include "draco.hpp"
#include "args.hpp"
#include "results/transcript.hpp"
#include "ringmap_matrix.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
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

  results::Transcript transcript_results(1);
  merge_windows_and_add_window_results(windows, windows_reads_indices,
                                       ringmap_data, transcript_results,
                                       WINDOW_SIZE, 0, args);

  assert(transcript_results.windows[0].has_value());
  assert(std::size(*transcript_results.windows[0]) == WINDOWS);
  assert(std::ranges::equal(
      *transcript_results.windows[0] |
          std::views::transform([](auto const &window) {
            return window.weighted_clusters.getClustersSize();
          }),
      std::array{1, 2, 1, 2}));
}

void test_make_windows_and_reads_indices_range_same_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(30)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, n_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 30u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();

  auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
      std::ranges::cbegin(windows), std::ranges::cend(windows),
      std::ranges::cbegin(windows_reads_indices),
      std::ranges::cend(windows_reads_indices), 0);

  assert(std::ranges::begin(windows_and_reads_indices_range) !=
         std::ranges::end(windows_and_reads_indices_range));
  assert(std::addressof(std::get<0>(
             *std::ranges::begin(windows_and_reads_indices_range))) ==
         std::addressof(*std::begin(windows)));
  assert(std::ranges::distance(windows_and_reads_indices_range) ==
         std::ranges::distance(windows));
}

void test_make_windows_and_reads_indices_range_separated_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(90)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        std::size_t actual_clusters;
        if (start_base < 30 || start_base >= 60) {
          actual_clusters = n_clusters;
        } else {
          actual_clusters = n_clusters + 1;
        }

        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, actual_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 90u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        std::ranges::cbegin(windows), std::ranges::cend(windows),
        std::ranges::cbegin(windows_reads_indices),
        std::ranges::cend(windows_reads_indices), 0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<0>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*std::begin(windows)));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 11);
  }

  {
    auto windows_iter = std::ranges::next(std::ranges::cbegin(windows), 10);
    auto windows_reads_indices_iter =
        std::ranges::next(std::ranges::cbegin(windows_reads_indices), 10);
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<0>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*windows_iter));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 11);
  }

  {
    auto windows_iter = std::ranges::next(std::ranges::cbegin(windows), 20);
    auto windows_reads_indices_iter =
        std::ranges::next(std::ranges::cbegin(windows_reads_indices), 20);
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<0>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*windows_iter));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 10);
  }
}

void test_make_windows_and_reads_indices_range_overlapping_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(90)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        std::size_t actual_clusters;
        if (start_base < 30 || (start_base > 39 && start_base <= 51)) {
          actual_clusters = n_clusters;
        } else {
          actual_clusters = n_clusters + 1;
        }

        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, actual_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 90u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        std::ranges::cbegin(windows), std::ranges::cend(windows),
        std::ranges::cbegin(windows_reads_indices),
        std::ranges::cend(windows_reads_indices), 2);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<0>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*std::begin(windows)));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 19);
  }
}

void test_update_iters_and_region_same_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(30)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, n_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 30u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();

  std::vector<std::optional<std::uint16_t>> previous_region_ends(n_clusters);
  auto windows_iter = std::ranges::cbegin(windows);
  auto windows_reads_indices_iter = std::ranges::cbegin(windows_reads_indices);

  auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
      windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
      std::ranges::cend(windows_reads_indices), 0);

  update_iters_and_region(windows_iter, windows_reads_indices_iter,
                          windows_and_reads_indices_range, previous_region_ends,
                          0);

  assert(windows_iter == std::ranges::cend(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::cend(windows_reads_indices));
  assert(std::size(previous_region_ends) == n_clusters);
  assert(
      std::ranges::all_of(previous_region_ends, [](auto &&previous_region_end) {
        return not previous_region_end.has_value();
      }));
}

void test_update_iters_and_region_separated_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(90)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        std::size_t actual_clusters;
        if (start_base < 30 || start_base >= 60) {
          actual_clusters = n_clusters;
        } else {
          actual_clusters = n_clusters + 1;
        }

        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, actual_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 90u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();

  std::vector<std::optional<std::uint16_t>> previous_region_ends(n_clusters +
                                                                 1);
  auto windows_iter = std::ranges::cbegin(windows);
  auto windows_reads_indices_iter = std::ranges::cbegin(windows_reads_indices);

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, 0);
  }

  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 10));
  assert(windows_iter->weights.getClustersSize() == n_clusters + 1);
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::cbegin(windows_reads_indices), 10));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 27 + window_size);
  assert(not previous_region_ends[2].has_value());

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, 0);
  }

  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 20));
  assert(windows_iter->weights.getClustersSize() == n_clusters);
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::cbegin(windows_reads_indices), 20));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 27 + window_size);
  assert(previous_region_ends[2].has_value());
  assert(*previous_region_ends[2] == 57 + window_size);

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, 0);
  }

  assert(windows_iter == std::ranges::cend(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::cend(windows_reads_indices));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 27 + window_size);
  assert(previous_region_ends[2].has_value());
  assert(*previous_region_ends[2] == 57 + window_size);
}

void test_update_iters_and_region_overlapping_clusters() {
  constexpr std::size_t n_clusters = 2;
  constexpr std::size_t window_size = 17;
  constexpr std::size_t window_stride = 3;

  auto windows =
      std::views::iota(static_cast<std::uint16_t>(0),
                       static_cast<std::uint16_t>(60)) |
      std::views::filter(
          [&](auto index) { return index % window_stride == 0; }) |
      std::views::transform([](auto start_base) {
        std::size_t actual_clusters;
        if (start_base < 30 || (start_base > 39 && start_base <= 51)) {
          actual_clusters = n_clusters;
        } else {
          actual_clusters = n_clusters + 1;
        }

        return Window{.start_base = start_base,
                      .weights = WeightedClusters(window_size, actual_clusters),
                      .coverages = {}};
      }) |
      std::ranges::to<std::vector>();
  auto windows_reads_indices =
      std::views::iota(0u, 60u) | std::views::filter([&](auto index) {
        return index % window_stride == 0;
      }) |
      std::views::transform([](auto start_base) {
        return std::views::iota(start_base, start_base + window_size) |
               std::ranges::to<std::vector>();
      }) |
      std::ranges::to<std::vector>();
  std::vector<std::optional<std::uint16_t>> previous_region_ends(n_clusters +
                                                                 1);
  auto windows_iter = std::ranges::cbegin(windows);
  auto windows_reads_indices_iter = std::ranges::cbegin(windows_reads_indices);

  constexpr auto min_bases_overlap = 8u;
  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), min_bases_overlap);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, min_bases_overlap);
  }

  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 10));
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::cbegin(windows_reads_indices), 10));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 27 + window_size - min_bases_overlap);
  assert(not previous_region_ends[2].has_value());

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), min_bases_overlap);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, min_bases_overlap);
  }

  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 14));
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::cbegin(windows_reads_indices), 14));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 27 + window_size - min_bases_overlap);
  assert(previous_region_ends[2].has_value());
  assert(*previous_region_ends[2] == 39 + window_size - min_bases_overlap);

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), min_bases_overlap);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, min_bases_overlap);
  }

  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 18));
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::cbegin(windows_reads_indices), 18));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 51 + window_size - min_bases_overlap);
  assert(previous_region_ends[2].has_value());
  assert(*previous_region_ends[2] == 39 + window_size - min_bases_overlap);

  {
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), 0);

    update_iters_and_region(windows_iter, windows_reads_indices_iter,
                            windows_and_reads_indices_range,
                            previous_region_ends, 0);
  }

  assert(windows_iter == std::ranges::cend(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::cend(windows_reads_indices));
  assert(not previous_region_ends[0].has_value());
  assert(previous_region_ends[1].has_value());
  assert(*previous_region_ends[1] == 51 + window_size - min_bases_overlap);
  assert(previous_region_ends[2].has_value());
  assert(*previous_region_ends[2] == 39 + window_size - min_bases_overlap);
}

int main() {
  test_merge_windows_and_add_window_results_not_merging();
  test_make_windows_and_reads_indices_range_same_clusters();
  test_make_windows_and_reads_indices_range_separated_clusters();
  test_make_windows_and_reads_indices_range_overlapping_clusters();
  test_update_iters_and_region_same_clusters();
  test_update_iters_and_region_separated_clusters();
  test_update_iters_and_region_overlapping_clusters();
}

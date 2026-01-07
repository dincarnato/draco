#include "draco.hpp"
#include "args.hpp"
#include "mutation_map.hpp"
#include "mutation_map_transcript.hpp"
#include "results/analysis.hpp"
#include "results/transcript.hpp"
#include "ringmap_matrix.hpp"
#include "weighted_clusters.hpp"
#include "window_clusters_with_confidence.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string_view>
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
      std::ranges::cbegin(windows), std::ranges::cbegin(windows),
      std::ranges::cend(windows), std::ranges::cbegin(windows_reads_indices),
      std::ranges::cend(windows_reads_indices), 0);

  assert(std::ranges::begin(windows_and_reads_indices_range) !=
         std::ranges::end(windows_and_reads_indices_range));
  assert(std::addressof(std::get<1>(
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
        std::ranges::cbegin(windows), std::ranges::cbegin(windows),
        std::ranges::cend(windows), std::ranges::cbegin(windows_reads_indices),
        std::ranges::cend(windows_reads_indices), 0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<1>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*std::begin(windows)));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 11);
  }

  {
    auto windows_iter = std::ranges::next(std::ranges::cbegin(windows), 10);
    auto windows_reads_indices_iter =
        std::ranges::next(std::ranges::cbegin(windows_reads_indices), 10);
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<1>(
               *std::ranges::begin(windows_and_reads_indices_range))) ==
           std::addressof(*windows_iter));
    assert(std::ranges::distance(windows_and_reads_indices_range) == 11);
  }

  {
    auto windows_iter = std::ranges::next(std::ranges::cbegin(windows), 20);
    auto windows_reads_indices_iter =
        std::ranges::next(std::ranges::cbegin(windows_reads_indices), 20);
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<1>(
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
        std::ranges::cbegin(windows), std::ranges::cbegin(windows),
        std::ranges::cend(windows), std::ranges::cbegin(windows_reads_indices),
        std::ranges::cend(windows_reads_indices), 2);

    assert(std::ranges::begin(windows_and_reads_indices_range) !=
           std::ranges::end(windows_and_reads_indices_range));
    assert(std::addressof(std::get<1>(
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
      std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
      windows_reads_indices_iter, std::ranges::cend(windows_reads_indices), 0);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        min_bases_overlap);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        min_bases_overlap);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        min_bases_overlap);

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
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        0);

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

void test_get_best_pre_collapsing_clusters_one_window() {
  constexpr std::size_t window_size = 10;
  constexpr std::size_t window_offset = 1;

  auto const make_window = [&](std::uint8_t n_clusters,
                               unsigned short start_base) {
    return Window{
        .start_base = start_base,
        .weights = WeightedClusters(window_size, n_clusters),
        .coverages = std::vector<unsigned>(),
    };
  };

  std::vector ptba_on_replicate_results{
      std::optional(
          PtbaOnReplicate{.pre_collapsing_clusters = std::vector{3u, 5u},
                          .windows =
                              std::vector{
                                  make_window(3, 0),
                                  make_window(5, 5),
                              },
                          .window_size = window_size,
                          .window_offset = window_offset}),
      std::optional(PtbaOnReplicate{
          .pre_collapsing_clusters = std::vector{1u, 1u},
          .windows =
              std::vector{
                  make_window(1, 0),
                  make_window(1, 5),
              },
          .window_size = window_size,
          .window_offset = window_offset,
      }),
      std::optional(PtbaOnReplicate{
          .pre_collapsing_clusters = std::vector{4u, 2u},
          .windows =
              std::vector{
                  make_window(4, 0),
                  make_window(2, 5),
              },
          .window_size = window_size,
          .window_offset = window_offset,
      }),
      std::optional(PtbaOnReplicate{
          .pre_collapsing_clusters = std::vector{2u, 1u},
          .windows =
              std::vector{
                  make_window(2, 0),
                  make_window(1, 5),
              },
          .window_size = window_size,
          .window_offset = window_offset,
      }),
  };

  {
    auto results = get_best_pre_collapsing_clusters(ptba_on_replicate_results,
                                                    "transcript");
    assert(std::size(results) == 2);
    assert(results[0].n_clusters == 2);
    assert(results[0].confidence == 1. / 4.);

    assert(results[0].n_clusters == 1);
    assert(results[0].confidence == 1. / 2.);
  }

  ptba_on_replicate_results.emplace_back(PtbaOnReplicate{
      .pre_collapsing_clusters = std::vector{6u, 6u},
      .windows =
          std::vector{
              make_window(6, 0),
              make_window(6, 5),
          },
      .window_size = window_size,
      .window_offset = window_offset,
  });
  {
    auto results = get_best_pre_collapsing_clusters(ptba_on_replicate_results,
                                                    "transcript");
    assert(std::size(results) == 2);
    assert(results[0].n_clusters == 3);
    assert(results[0].confidence == 1. / 5.);
    assert(results[1].n_clusters == 2);
    assert(results[0].confidence == 1. / 5.);
  }
}

void test_window_info_get_n_windows_and_precise_offset() {
  {
    auto n_windows_and_precise_offset =
        get_n_windows_and_precise_offset(29, 7, 3);
    assert(n_windows_and_precise_offset.n_windows == 8);
    assert(std::abs(n_windows_and_precise_offset.window_precise_offset -
                    3.142857143) < 0.0001);
  }

  {
    auto n_windows_and_precise_offset =
        get_n_windows_and_precise_offset(28, 7, 3);
    assert(n_windows_and_precise_offset.n_windows == 8);
    assert(n_windows_and_precise_offset.window_precise_offset == 3);
  }

  {
    auto n_windows_and_precise_offset =
        get_n_windows_and_precise_offset(5, 4, 3);
    assert(n_windows_and_precise_offset.n_windows == 1);
    assert(n_windows_and_precise_offset.window_precise_offset == 0.);
  }

  {
    auto n_windows_and_precise_offset =
        get_n_windows_and_precise_offset(5, 5, 3);
    assert(n_windows_and_precise_offset.n_windows == 1);
    assert(n_windows_and_precise_offset.window_precise_offset == 0.);
  }

  {
    auto n_windows_and_precise_offset =
        get_n_windows_and_precise_offset(10, 9, 1);
    assert(n_windows_and_precise_offset.n_windows == 2);
    assert(n_windows_and_precise_offset.window_precise_offset == 1.);
  }
}

void test_window_info_get_start_base() {
  {
    auto windows_info = WindowsInfo::from_size_and_offset(29, 7, 3);
    assert(windows_info.get_start_base(0) == 0);
    assert(windows_info.get_start_base(1) == 3);
    assert(windows_info.get_start_base(4) == 13);
    assert(windows_info.get_start_base(7) == 22);
  }

  {
    auto windows_info = WindowsInfo::from_size_and_offset(28, 7, 3);
    assert(windows_info.get_start_base(0) == 0);
    assert(windows_info.get_start_base(1) == 3);
    assert(windows_info.get_start_base(4) == 12);
    assert(windows_info.get_start_base(7) == 21);
  }
}

void test_add_detected_clusters_with_confidence() {
  auto pre_collapsing_clusters =
      std::views::repeat(
          PreCollapsingClusters{.n_clusters = 999, .confidence = 1.}) |
      std::views::take(10) | std::ranges::to<std::vector>();
  std::ranges::for_each(std::views::repeat(PreCollapsingClusters{
                            .n_clusters = 1, .confidence = 1}) |
                            std::views::take(5) | std::views::as_rvalue,
                        [&](auto &&value) {
                          pre_collapsing_clusters.push_back(std::move(value));
                        });
  std::ranges::for_each(std::views::repeat(PreCollapsingClusters{
                            .n_clusters = 1, .confidence = 0.5}) |
                            std::views::take(3) | std::views::as_rvalue,
                        [&](auto &&value) {
                          pre_collapsing_clusters.push_back(std::move(value));
                        });
  std::ranges::for_each(std::views::repeat(PreCollapsingClusters{
                            .n_clusters = 2, .confidence = 0.5}) |
                            std::views::take(4) | std::views::as_rvalue,
                        [&](auto &&value) {
                          pre_collapsing_clusters.push_back(std::move(value));
                        });
  std::ranges::for_each(std::views::repeat(PreCollapsingClusters{
                            .n_clusters = 2, .confidence = 1}) |
                            std::views::take(2) | std::views::as_rvalue,
                        [&](auto &&value) {
                          pre_collapsing_clusters.push_back(std::move(value));
                        });
  std::ranges::for_each(std::views::repeat(PreCollapsingClusters{
                            .n_clusters = 999, .confidence = 1}) |
                            std::views::take(5) | std::views::as_rvalue,
                        [&](auto &&value) {
                          pre_collapsing_clusters.push_back(std::move(value));
                        });

  std::vector<WindowClustersWithConfidence> detected_clusters_with_confidence;
  std::size_t window_size = 17;
  std::size_t window_offset = 5;
  auto windows_info =
      WindowsInfo::from_size_and_offset(165, static_cast<unsigned>(window_size),
                                        static_cast<unsigned>(window_offset));
  add_detected_clusters_with_confidence(
      detected_clusters_with_confidence,
      results::WindowRange{.window_index_begin = 10, .window_index_end = 24},
      pre_collapsing_clusters, windows_info);

  assert(std::size(detected_clusters_with_confidence) == 4);

  assert(detected_clusters_with_confidence[0].n_clusters == 1);
  assert(detected_clusters_with_confidence[0].confidence == 1.);
  assert(detected_clusters_with_confidence[0].start_base ==
         windows_info.get_start_base(10));
  assert(detected_clusters_with_confidence[0].end_base ==
         windows_info.get_start_base(14) + window_size);

  assert(detected_clusters_with_confidence[1].n_clusters == 1);
  assert(detected_clusters_with_confidence[1].confidence == 0.5);
  assert(detected_clusters_with_confidence[1].start_base ==
         windows_info.get_start_base(15));
  assert(detected_clusters_with_confidence[1].end_base ==
         windows_info.get_start_base(17) + window_size);

  assert(detected_clusters_with_confidence[2].n_clusters == 2);
  assert(detected_clusters_with_confidence[2].confidence == 0.5);
  assert(detected_clusters_with_confidence[2].start_base ==
         windows_info.get_start_base(18));
  assert(detected_clusters_with_confidence[2].end_base ==
         windows_info.get_start_base(21) + window_size);

  assert(detected_clusters_with_confidence[3].n_clusters == 2);
  assert(detected_clusters_with_confidence[3].confidence == 1.);
  assert(detected_clusters_with_confidence[3].start_base ==
         windows_info.get_start_base(22));
  assert(detected_clusters_with_confidence[3].end_base ==
         windows_info.get_start_base(23) + window_size);
}

void test_handle_transcripts_clusters_confidences() {
  constexpr std::string_view sequence = "ACACCCAAACACAAACAAACCCACAAACACAACACAC";
  constexpr unsigned n_reads = 1000;
  constexpr unsigned window_size = 10;
  constexpr unsigned window_offset = 3;
  constexpr unsigned n_windows =
      (std::size(sequence) - window_size) / window_offset;
  assert((std::size(sequence) - window_size) % window_offset == 0);

  std::vector<MutationMapTranscript> owned_transcripts{
      MutationMapTranscript(MutationMap(), "test_1", 0),
      MutationMapTranscript(MutationMap(), "test_2", 0),
      MutationMapTranscript(MutationMap(), "test_3", 0),
  };

  for (auto &transcript : owned_transcripts) {
    transcript.setSequence(sequence);
    transcript.setReads(n_reads);
  }

  auto transcripts = owned_transcripts |
                     std::views::transform(
                         [](auto const &transcript) { return &transcript; }) |
                     std::ranges::to<std::vector>();

  Args args;

  std::vector<RingmapData> owned_ringmap_data{
      RingmapData(sequence, RingmapMatrix(n_reads, std::size(sequence)), 0,
                  std::size(sequence), args),
      RingmapData(sequence, RingmapMatrix(n_reads, std::size(sequence)), 0,
                  std::size(sequence), args),
      RingmapData(sequence, RingmapMatrix(n_reads, std::size(sequence)), 0,
                  std::size(sequence), args),
  };
  auto ringmap_data =
      owned_ringmap_data |
      std::views::transform([](auto &ringmap_data) { return &ringmap_data; }) |
      std::ranges::to<std::vector>();

  results::Analysis analysis_result;

  std::optional<std::ofstream> raw_n_clusters_stream = std::nullopt;
  std::mutex raw_n_clusters_stream_mutex;

  auto const make_window = [&](std::uint8_t n_clusters,
                               unsigned short start_base) {
    return Window{
        .start_base = start_base,
        .weights = WeightedClusters(window_size, n_clusters),
        .coverages = std::vector<unsigned>(),
    };
  };

  HandleTranscripts{.transcripts = transcripts,
                    .ringmaps_data = ringmap_data,
                    .analysis_result = analysis_result,
                    .args = args,
                    .raw_n_clusters_stream = raw_n_clusters_stream,
                    .raw_n_clusters_stream_mutex = raw_n_clusters_stream_mutex,
                    .use_logger = false,
                    .allow_empty_patterns = true}(
      [&](auto replicate_index, auto const &, auto &) {
        std::vector<unsigned> pre_collapsing_clusters(n_windows, 1);
        switch (replicate_index) {
        case 0:
          break;
        case 1:
          std::ranges::fill_n(
              std::next(std::ranges::begin(pre_collapsing_clusters),
                        n_windows / 3),
              n_windows / 3, 2);
          break;
        case 2:
          std::ranges::fill(
              pre_collapsing_clusters | std::views::drop(n_windows / 2), 2);
          break;
        default:
          throw new std::runtime_error("unexpected replicate index");
        }

        auto windows = std::views::iota(0u, n_windows) |
                       std::views::transform([&](auto window_index) {
                         auto start_base = static_cast<unsigned short>(
                             window_index * window_offset);
                         return make_window(1, start_base);
                       }) |
                       std::ranges::to<std::vector>();

        PtbaOnReplicate ptba_on_replicate{
            .pre_collapsing_clusters = pre_collapsing_clusters,
            .windows = windows,
            .window_size = window_size,
            .window_offset = window_offset,
        };
        return std::optional(std::move(ptba_on_replicate));
      },
      [&](unsigned short start_base, unsigned short end_base,
          std::uint8_t n_clusters, std::vector<arma::mat> const &,
          results::Transcript const &) {
        return WeightedClusters(end_base - start_base, n_clusters);
      });

  auto &analysis_transcripts = analysis_result.transcripts();
  assert(std::size(analysis_transcripts) == 1);
  auto const &transcript = analysis_transcripts.front();
  assert(std::size(transcript.detected_clusters_with_confidence) == 4);

  assert(transcript.detected_clusters_with_confidence[0].start_base == 0);
  assert(transcript.detected_clusters_with_confidence[0].end_base ==
         (n_windows / 3 - 1) * window_offset + window_size);
  assert(transcript.detected_clusters_with_confidence[0].confidence == 1.);
  assert(transcript.detected_clusters_with_confidence[0].n_clusters == 1);

  assert(transcript.detected_clusters_with_confidence[1].start_base ==
         n_windows / 3 * window_offset);
  assert(transcript.detected_clusters_with_confidence[1].end_base ==
         (n_windows / 2 - 1) * window_offset + window_size);
  assert(std::abs(transcript.detected_clusters_with_confidence[1].confidence -
                  2. / 3.) < 0.0001);
  assert(transcript.detected_clusters_with_confidence[1].n_clusters == 1);

  assert(transcript.detected_clusters_with_confidence[2].start_base ==
         n_windows / 2 * window_offset);
  assert(transcript.detected_clusters_with_confidence[2].end_base ==
         (n_windows * 2 / 3 - 1) * window_offset + window_size);
  assert(std::abs(transcript.detected_clusters_with_confidence[2].confidence -
                  2. / 3.) < 0.0001);
  assert(transcript.detected_clusters_with_confidence[2].n_clusters == 2);

  assert(transcript.detected_clusters_with_confidence[3].start_base ==
         n_windows * 2 / 3 * window_offset);
  assert(transcript.detected_clusters_with_confidence[3].end_base ==
         (n_windows - 1) * window_offset + window_size);
  assert(std::abs(transcript.detected_clusters_with_confidence[3].confidence -
                  2. / 3.) < 0.0001);
  assert(transcript.detected_clusters_with_confidence[3].n_clusters == 1);
}

static void test_handle_overlapping_regions_should_break() {
  std::vector<Window> windows{Window{.start_base = 23,
                                     .weights = WeightedClusters(100, 2),
                                     .coverages = std::vector<unsigned>()}};
  auto windows_iter = std::ranges::cbegin(windows);
  std::vector<std::vector<unsigned>> windows_reads_indices(
      2, std::vector<unsigned>{0, 1, 2});
  std::span windows_reads_indices_span(windows_reads_indices);
  auto windows_reads_indices_iter =
      std::ranges::begin(windows_reads_indices_span);

  std::optional<short unsigned> previous_overlapping_region_end;
  auto result = handle_overlapping_regions(
      previous_overlapping_region_end, windows_iter, std::ranges::cend(windows),
      windows_reads_indices_iter, 2, 5);
  assert(result == HandleOverlappingRegionsResult::Break);
  assert(not previous_overlapping_region_end.has_value());
  assert(windows_iter == std::ranges::cbegin(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::begin(windows_reads_indices_span));

  previous_overlapping_region_end = 15;
  result = handle_overlapping_regions(previous_overlapping_region_end,
                                      windows_iter, std::ranges::cend(windows),
                                      windows_reads_indices_iter, 2, 5);
  assert(result == HandleOverlappingRegionsResult::Break);
  assert(*previous_overlapping_region_end == 15);
  assert(windows_iter == std::ranges::cbegin(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::begin(windows_reads_indices_span));
}

static void test_handle_overlapping_regions_advance_iters_and_region_end() {
  std::vector<Window> windows{
      Window{.start_base = 23,
             .weights = WeightedClusters(100, 2),
             .coverages = std::vector<unsigned>()},
      Window{.start_base = 25,
             .weights = WeightedClusters(100, 2),
             .coverages = std::vector<unsigned>()},
      Window{.start_base = 27,
             .weights = WeightedClusters(100, 1),
             .coverages = std::vector<unsigned>()},
      Window{.start_base = 29,
             .weights = WeightedClusters(100, 1),
             .coverages = std::vector<unsigned>()},
      Window{.start_base = 31,
             .weights = WeightedClusters(100, 2),
             .coverages = std::vector<unsigned>()},
  };

  auto windows_iter = std::ranges::cbegin(windows);
  std::vector<std::vector<unsigned>> windows_reads_indices(std::size(windows));
  std::span windows_reads_indices_span(windows_reads_indices);
  auto windows_reads_indices_iter =
      std::ranges::begin(windows_reads_indices_span);

  std::vector<std::optional<short unsigned>> previous_overlapping_region_ends{
      31, 27};
  auto result = handle_overlapping_regions(
      previous_overlapping_region_ends[1], windows_iter,
      std::ranges::cend(windows), windows_reads_indices_iter, 2, 5);
  assert(result == HandleOverlappingRegionsResult::Continue);
  assert(previous_overlapping_region_ends[1] == 120);
  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 2));
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::begin(windows_reads_indices_span), 2));

  result = handle_overlapping_regions(previous_overlapping_region_ends[0],
                                      windows_iter, std::ranges::cend(windows),
                                      windows_reads_indices_iter, 1, 5);
  assert(result == HandleOverlappingRegionsResult::Continue);
  assert(previous_overlapping_region_ends[0] == 124);
  assert(windows_iter == std::ranges::next(std::ranges::cbegin(windows), 4));
  assert(windows_reads_indices_iter ==
         std::ranges::next(std::ranges::begin(windows_reads_indices_span), 4));

  result = handle_overlapping_regions(previous_overlapping_region_ends[1],
                                      windows_iter, std::ranges::cend(windows),
                                      windows_reads_indices_iter, 2, 5);
  assert(result == HandleOverlappingRegionsResult::Continue);
  assert(previous_overlapping_region_ends[1] == 126);
  assert(windows_iter == std::ranges::cend(windows));
  assert(windows_reads_indices_iter ==
         std::ranges::end(windows_reads_indices_span));
}

int main() {
  test_merge_windows_and_add_window_results_not_merging();
  test_make_windows_and_reads_indices_range_same_clusters();
  test_make_windows_and_reads_indices_range_separated_clusters();
  test_make_windows_and_reads_indices_range_overlapping_clusters();
  test_update_iters_and_region_same_clusters();
  test_update_iters_and_region_separated_clusters();
  test_update_iters_and_region_overlapping_clusters();
  test_window_info_get_n_windows_and_precise_offset();
  test_window_info_get_start_base();
  test_add_detected_clusters_with_confidence();
  test_handle_transcripts_clusters_confidences();
  test_handle_overlapping_regions_should_break();
  test_handle_overlapping_regions_advance_iters_and_region_end();
}

#include "draco.hpp"
#include "args.hpp"
#include "fmt/base.h"
#include "fmt/ostream.h"
#include "graph_cut.hpp"
#include "logger.hpp"
#include "mutation_map_transcript.hpp"
#include "mutation_map_writer.hpp"
#include "ptba.hpp"
#include "results/analysis.hpp"
#include "results/transcript.hpp"
#include "results/window.hpp"
#include "ringmap_data.hpp"
#include "to_vector.hpp"
#include "windows_merger.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <iterator>
#include <limits>
#include <mutex>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <optional>
#include <ranges>
#include <set>
#include <sstream>
#include <string_view>

#include <armadillo>
#include <tbb/global_control.h>

#include <omp.h>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace fs = std::filesystem;
namespace ranges = std::ranges;

struct WindowsSpan {
  std::size_t begin;
  std::size_t end;
  unsigned n_clusters;
  unsigned max_n_clusters;

  constexpr std::ptrdiff_t size() const noexcept {
    return static_cast<std::ptrdiff_t>(end) -
           static_cast<std::ptrdiff_t>(begin);
  }
};

constexpr static auto const invalid_n_clusters =
    std::numeric_limits<unsigned>::max();

static std::vector<WindowsSpan>
get_windows_spans(std::vector<unsigned> const &windows_n_clusters,
                  std::vector<std::optional<unsigned>> const
                      &windows_max_clusters_constraints) noexcept(false) {

  auto get_max_n_clusters = [&](WindowsSpan const &span) {
    if (span.end == span.begin) {
      return invalid_n_clusters - 1;
    } else {
      return ranges::min(windows_max_clusters_constraints |
                         std::views::drop(span.begin) |
                         std::views::take(span.end - span.begin) |
                         std::views::transform([](auto &&constraint) {
                           if (constraint)
                             return *constraint;
                           else {
                             return invalid_n_clusters - 1;
                           }
                         }));
    }
    unsigned max_n_clusters = [&] {
      if (auto &&constraint = windows_max_clusters_constraints[span.begin];
          constraint) {
        return *constraint;
      } else {
        return invalid_n_clusters - 1;
      }
    }();

    for (auto window_index = span.begin; window_index < span.end;
         ++window_index) {
      if (auto const &constraint =
              windows_max_clusters_constraints[window_index];
          constraint) {
        max_n_clusters = std::min(max_n_clusters, *constraint);
      }
    }

    return max_n_clusters;
  };

  auto const windows_size = windows_n_clusters.size();
  std::vector<WindowsSpan> windows_spans;
  {
    WindowsSpan span{std::size_t(0), std::size_t(0), invalid_n_clusters,
                     invalid_n_clusters - 1};
    for (std::size_t window_index = 0; window_index < windows_size;
         ++window_index) {
      auto const window_n_clusters = windows_n_clusters[window_index];
      if (window_n_clusters != span.n_clusters) {
        span.end = window_index;
        if (span.size() > 0) {
          span.max_n_clusters = get_max_n_clusters(span);
          windows_spans.emplace_back(span);
        }

        span.begin = window_index;
        span.n_clusters = window_n_clusters;
      }
    }

    span.end = windows_size;
    if (span.size() > 0) {
      span.max_n_clusters = get_max_n_clusters(span);
      windows_spans.emplace_back(span);
    }
  }

  return windows_spans;
}

bool merge_spans(std::vector<WindowsSpan> &windows_spans) noexcept {
  bool updated = false;

  ranges::sort(windows_spans, {},
               [](auto &&windows_span) { return windows_span.begin; });

  auto windows_spans_iter = std::begin(windows_spans);
  auto const windows_spans_end = std::end(windows_spans);

  for (; windows_spans_iter < windows_spans_end; ++windows_spans_iter) {
    auto &first_windows_span = *windows_spans_iter;

    auto const span_clusters = first_windows_span.n_clusters;
    auto span_max_clusters = first_windows_span.max_n_clusters;
    auto other_windows_spans_iter = std::next(windows_spans_iter);
    auto last_span_end = first_windows_span.end;
    for (; other_windows_spans_iter < windows_spans_end;
         ++other_windows_spans_iter) {
      auto &&other_windows_span = *other_windows_spans_iter;

      if (other_windows_span.n_clusters != span_clusters)
        break;
      else {
        last_span_end = other_windows_span.end;
        span_max_clusters =
            std::min(span_max_clusters, other_windows_span.max_n_clusters);
      }
    }

    if (first_windows_span.end != last_span_end) {
      first_windows_span.end = last_span_end;
      first_windows_span.max_n_clusters = span_max_clusters;
      std::for_each(std::next(windows_spans_iter), other_windows_spans_iter,
                    [&](auto &&windows_span) {
                      windows_span.n_clusters = invalid_n_clusters;
                    });
      updated = true;
    }
  }

  return updated;
}

enum class LoopAction {
  None,
  Continue,
  Break,
};

template <typename F1, typename F2, typename F3>
void windows_span_expander(std::vector<unsigned> &windows_n_clusters,
                           std::vector<std::optional<unsigned>> const
                               &windows_max_clusters_constraints,
                           F1 &&loop_start_check, F2 &&same_clusters_check,
                           F3 &&different_clusters_check) noexcept(false) {
  auto const windows_size = windows_n_clusters.size();
  auto windows_spans =
      get_windows_spans(windows_n_clusters, windows_max_clusters_constraints);

  for (;;) {
    ranges::sort(windows_spans, {},
                 [](auto &&window_span) { return window_span.size(); });

    assert(ranges::none_of(windows_spans, [](auto &&windows_span) {
      return windows_span.n_clusters == invalid_n_clusters;
    }));

    assert(ranges::all_of(windows_spans, [](auto &&windows_span) {
      return windows_span.n_clusters <= windows_span.max_n_clusters;
    }));

    assert(std::ranges::all_of(windows_spans, [&](auto &&windows_span) {
      return ranges::all_of(
          windows_max_clusters_constraints |
              std::views::drop(windows_span.begin) |
              std::views::take(windows_span.end - windows_span.begin) |
              ranges::views::filter(
                  [](auto &&constraint) { return constraint.has_value(); }),
          [&](auto &&constraint) {
            return windows_span.max_n_clusters <= *constraint;
          });
    }));

    bool updated = false;
    for (auto &&window_span : windows_spans) {
      if (auto const action = loop_start_check(window_span);
          action == LoopAction::Continue) {
        continue;
      } else if (action == LoopAction::Break) {
        break;
      }

      auto const span_max_clusters = window_span.max_n_clusters;
      assert(window_span.n_clusters <= span_max_clusters);

      auto left_windows_span_ptr = [&] {
        if (window_span.begin > 0) {
          auto const left_window_span_iter =
              ranges::find(windows_spans, window_span.begin,
                           [](auto &&window_span) { return window_span.end; });
          assert(left_window_span_iter != ranges::end(windows_spans));
          return &*left_window_span_iter;
        } else
          return static_cast<WindowsSpan *>(nullptr);
      }();

      auto right_windows_span_ptr = [&] {
        if (window_span.end < windows_size) {
          auto const right_window_span_iter = ranges::find(
              windows_spans, window_span.end,
              [](auto &&window_span) { return window_span.begin; });
          assert(right_window_span_iter != ranges::end(windows_spans));
          return &*right_window_span_iter;
        } else
          return static_cast<WindowsSpan *>(nullptr);
      }();

      if (left_windows_span_ptr) {
        auto &left_windows_span = *left_windows_span_ptr;

        if (right_windows_span_ptr) {
          auto &right_windows_span = *right_windows_span_ptr;

          if (left_windows_span.n_clusters == right_windows_span.n_clusters) {
            if (same_clusters_check(left_windows_span, right_windows_span)) {

              if (left_windows_span.n_clusters <= span_max_clusters and
                  right_windows_span.n_clusters <= span_max_clusters) {
                left_windows_span.end = right_windows_span.end;
                right_windows_span.n_clusters = invalid_n_clusters;
                window_span.n_clusters = invalid_n_clusters;
                left_windows_span.max_n_clusters =
                    std::min(std::min(left_windows_span.max_n_clusters,
                                      right_windows_span.max_n_clusters),
                             span_max_clusters);

                updated = true;
                break;
              } else if (window_span.n_clusters != span_max_clusters) {
                window_span.n_clusters = span_max_clusters;
                updated = true;
                // Keep looping
              }
            }
          } else if (left_windows_span.n_clusters == span_max_clusters and
                     different_clusters_check(left_windows_span)) {
            left_windows_span.end = window_span.end;
            window_span.n_clusters = invalid_n_clusters;
            left_windows_span.max_n_clusters =
                std::min(left_windows_span.max_n_clusters, span_max_clusters);

            updated = true;
            break;
          } else if (right_windows_span.n_clusters == span_max_clusters and
                     different_clusters_check(right_windows_span)) {
            right_windows_span.begin = window_span.begin;
            window_span.n_clusters = invalid_n_clusters;
            right_windows_span.max_n_clusters =
                std::min(right_windows_span.max_n_clusters, span_max_clusters);

            updated = true;
            break;
          } else {
            auto const left_windows_span_size = left_windows_span.size();
            auto const right_windows_span_size = right_windows_span.size();
            if (left_windows_span_size >= right_windows_span_size) {
              if (different_clusters_check(left_windows_span)) {
                if (left_windows_span.n_clusters <= span_max_clusters) {
                  left_windows_span.end = window_span.end;
                  window_span.n_clusters = invalid_n_clusters;
                  left_windows_span.max_n_clusters = std::min(
                      left_windows_span.max_n_clusters, span_max_clusters);

                  updated = true;
                  break;
                } else if (window_span.n_clusters != span_max_clusters) {
                  window_span.n_clusters = span_max_clusters;
                  updated = true;
                  // Keep looping
                }
              }
            } else {
              if (different_clusters_check(right_windows_span)) {
                if (right_windows_span.n_clusters <= span_max_clusters) {
                  right_windows_span.begin = window_span.begin;
                  window_span.n_clusters = invalid_n_clusters;
                  right_windows_span.max_n_clusters = std::min(
                      right_windows_span.max_n_clusters, span_max_clusters);

                  updated = true;
                  break;
                } else if (window_span.n_clusters != span_max_clusters) {
                  window_span.n_clusters = span_max_clusters;
                  updated = true;
                  // Keep looping
                }
              }
            }
          }
        } else {
          if (different_clusters_check(left_windows_span)) {
            if (left_windows_span.n_clusters <= span_max_clusters) {
              left_windows_span.end = window_span.end;
              window_span.n_clusters = invalid_n_clusters;
              left_windows_span.max_n_clusters =
                  std::min(left_windows_span.max_n_clusters, span_max_clusters);

              updated = true;
              break;
            } else if (window_span.n_clusters != span_max_clusters) {
              window_span.n_clusters = span_max_clusters;
              updated = true;
              // Keep looping
            }
          }
        }
      } else if (right_windows_span_ptr) {
        auto &right_windows_span = *right_windows_span_ptr;

        if (different_clusters_check(right_windows_span)) {
          if (right_windows_span.n_clusters <= span_max_clusters) {
            right_windows_span.begin = window_span.begin;
            window_span.n_clusters = invalid_n_clusters;
            right_windows_span.max_n_clusters =
                std::min(right_windows_span.max_n_clusters, span_max_clusters);

            updated = true;
            break;
          } else if (window_span.n_clusters != span_max_clusters) {
            window_span.n_clusters = span_max_clusters;
            updated = true;
            // Keep looping
          }
        }
      }
    }

    if (not updated)
      updated = merge_spans(windows_spans);

    if (not updated)
      break;

    windows_spans.erase(
        std::remove_if(std::begin(windows_spans), std::end(windows_spans),
                       [](auto &&windows_span) {
                         return windows_span.n_clusters == invalid_n_clusters;
                       }),
        std::end(windows_spans));
  }

  auto const windows_n_clusters_begin = std::begin(windows_n_clusters);
  for (auto &&windows_span : windows_spans) {
    std::fill(std::next(windows_n_clusters_begin,
                        static_cast<std::ptrdiff_t>(windows_span.begin)),
              std::next(windows_n_clusters_begin,
                        static_cast<std::ptrdiff_t>(windows_span.end)),
              windows_span.n_clusters);
  }

  assert(ranges::all_of(windows_spans, [](auto &&windows_span) {
    return windows_span.n_clusters <= windows_span.max_n_clusters;
  }));

  assert(ranges::all_of(windows_spans, [&](auto &&windows_span) {
    return ranges::all_of(
        windows_max_clusters_constraints |
            std::views::drop(windows_span.begin) |
            std::views::take(windows_span.end - windows_span.begin) |
            std::views::filter(
                [](auto &&constraint) { return constraint.has_value(); }),
        [&](auto &&constraint) {
          return windows_span.max_n_clusters <= *constraint;
        });
  }));
}

void collapse_outlayer_clusters(std::vector<unsigned> &windows_n_clusters,
                                std::vector<std::optional<unsigned>> const
                                    &windows_max_clusters_constraints,
                                Args const &args) noexcept(false) {
  auto const max_collapsing_windows = args.max_collapsing_windows();
  auto const min_surrounding_windows_size = args.min_surrounding_windows_size();
  windows_span_expander(
      windows_n_clusters, windows_max_clusters_constraints,
      [=](auto &&window_span) {
        if (window_span.size() > max_collapsing_windows) {
          return LoopAction::Break;
        } else {
          return LoopAction::None;
        }
      },
      [=](auto &&left_windows_span, auto &&right_windows_span) {
        return left_windows_span.size() + right_windows_span.size() >=
               min_surrounding_windows_size;
      },
      [=](auto &&windows_span) {
        return windows_span.size() >= min_surrounding_windows_size;
      });
}

void set_uninformative_clusters_to_surrounding(
    std::vector<unsigned> &windows_n_clusters,
    std::vector<std::optional<unsigned>> const
        &windows_max_clusters_constraints) noexcept(false) {
  windows_span_expander(
      windows_n_clusters, windows_max_clusters_constraints,
      [](auto &&window_span) {
        if (window_span.n_clusters != 0) {
          return LoopAction::Continue;
        } else {
          return LoopAction::None;
        }
      },
      [](auto &&, auto &&) { return true; }, [](auto &&) { return true; });
}

void assign_reads_to_clusters(
    results::Window &window,
    RingmapData::clusters_assignment_type &&clusters_assignment,
    RingmapData const &ringmap, RingmapData const &filteredRingmap) {
  if (not window.bases_coverages) {
    window.bases_coverages = std::vector<std::vector<unsigned>>{};
  }
  auto &bases_coverages = *window.bases_coverages;
  bases_coverages.resize(
      window.weighted_clusters.getClustersSize(),
      std::vector<unsigned>(window.end_index - window.begin_index, 0));

  auto &&original_data = ringmap.data();

  auto &&original_indices_map = filteredRingmap.getReadsMap();
  for (auto &&clusters_assignment_pair : clusters_assignment) {
    auto &&read_clusters_assignments = std::get<1>(clusters_assignment_pair);

    auto &&clusters_assignments = read_clusters_assignments.clusters();
    auto const clusters_size = clusters_assignments.size();
    assert(window.weighted_clusters.getClustersSize() == clusters_size);
    assert(clusters_size <= std::numeric_limits<std::int8_t>::max());

    for (std::uint8_t cluster_index = 0;
         cluster_index < static_cast<std::uint8_t>(clusters_size);
         ++cluster_index) {
      auto &&cluster_assignments = clusters_assignments[cluster_index];

      for (auto filtered_read_index : cluster_assignments) {
        auto original_read_index = original_indices_map[filtered_read_index];
        window.assignments[original_read_index] =
            static_cast<std::int8_t>(cluster_index);

        auto &&original_row = original_data.row(original_read_index);

        auto const begin_index =
            std::max(original_row.begin_index(),
                     static_cast<unsigned>(window.begin_index));
        auto const end_index = std::min(
            original_row.end_index(), static_cast<unsigned>(window.end_index));

        auto &&cluster_bases_coverages = bases_coverages[cluster_index];
        assert(cluster_bases_coverages.size() >= end_index - begin_index);
        ranges::for_each(
            cluster_bases_coverages |
                std::views::drop(begin_index - window.begin_index) |
                std::views::take(end_index - begin_index),
            [](auto &&coverage) { ++coverage; });
      }
    }
  }
}

void dump_assignments(results::Transcript const &transcript,
                      results::Window &window, RingmapData const &ringmap,
                      std::optional<std::size_t> replicate_index,
                      std::string_view assignments_dump_directory) {
  auto mm_basename = std::format("{}_{}-{}", transcript.name,
                                 window.begin_index, window.end_index - 1);
  auto mm_filename = ([&]() {
    if (replicate_index.has_value()) {
      return std::format("{}_r{}.mm", mm_basename, *replicate_index + 1);
    } else {
      return std::format("{}.mm", mm_basename);
    }
  })();
  auto mm_path =
      std::filesystem::path{assignments_dump_directory} / mm_filename;
  MutationMapWriter mm_writer{mm_path};

  auto ringmap_data = ringmap.data();
  auto const n_reads = static_cast<std::uint32_t>(ringmap_data.rows_size());
  assert(n_reads == window.assignments.size());

  auto const n_clusters = static_cast<std::uint8_t>(window.fractions.size());
  for (std::uint8_t cluster_index = 0; cluster_index < n_clusters;
       ++cluster_index) {
    auto mm_transcript = mm_writer.transcript(
        std::format("{}_c{}", mm_basename, cluster_index), transcript.sequence);
    for (std::uint32_t read_index = 0; read_index < n_reads; ++read_index) {
      if (window.assignments[read_index] != cluster_index) {
        continue;
      }

      auto read = ringmap_data.row(read_index);
      assert(read.is_valid());

      std::stringstream indices;
      for (auto index : read.modifiedIndices()) {
        indices << index << ',';
      }

      mm_transcript.add_read(read);
    }
  }
}

void output_raw_n_clusters(std::ofstream &raw_n_clusters_stream,
                           std::mutex &raw_n_clusters_stream_mutex,
                           unsigned int window_size,
                           std::vector<Window> const &windows,
                           std::vector<unsigned int> const &windows_n_clusters,
                           results::Transcript &transcript_result) {
  auto windows_iter = std::cbegin(windows);
  auto windows_end = std::cend(windows);
  auto windows_n_clusters_iter = std::cbegin(windows_n_clusters);
  assert(std::distance(windows_iter, windows_end) ==
         std::distance(windows_n_clusters_iter, std::cend(windows_n_clusters)));

  {
    std::lock_guard<std::mutex> lock(raw_n_clusters_stream_mutex);

    for (; windows_iter < windows_end;
         ++windows_iter, ++windows_n_clusters_iter) {
      auto &&window = *windows_iter;
      auto n_clusters = *windows_n_clusters_iter;
      fmt::println(raw_n_clusters_stream, "{}\t{}\t{}\t{}",
                   transcript_result.name, window.start_base,
                   window.start_base + window_size, n_clusters);
    }
  }
}

void merge_windows_and_add_window_results(
    std::vector<Window> const &windows,
    std::span<std::vector<unsigned>> const &windows_reads_indices,
    RingmapData const &ringmap_data, results::Transcript &transcript_result,
    unsigned window_size, std::size_t replicate_index, Args const &args) {

  std::vector<std::optional<std::uint16_t>> previous_overlapping_region_ends;

  auto min_windows_overlap = static_cast<std::uint32_t>(
      static_cast<double>(window_size) *
      static_cast<double>(args.min_windows_overlap()));

  auto windows_iter = std::cbegin(windows);
  auto windows_reads_indices_iter = std::cbegin(windows_reads_indices);
  while (windows_iter != std::cend(windows)) {
    auto const &window = *windows_iter;
    assert(window.weights.getClustersSize() <=
           std::numeric_limits<unsigned>::max());
    auto const n_clusters =
        static_cast<unsigned>(window.weights.getClustersSize());
    if (std::size(previous_overlapping_region_ends) < n_clusters) {
      previous_overlapping_region_ends.resize(n_clusters);
    }
    if (n_clusters > 0) {
      auto &previous_overlapping_region_end =
          previous_overlapping_region_ends[n_clusters - 1];
      auto result = handle_overlapping_regions(
          previous_overlapping_region_end, windows_iter, std::cend(windows),
          windows_reads_indices_iter, n_clusters, min_windows_overlap);

      if (result == HandleOverlappingRegionsResult::Continue) {
        continue;
      }
    }

    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        std::ranges::cbegin(windows), windows_iter, std::ranges::cend(windows),
        windows_reads_indices_iter, std::ranges::cend(windows_reads_indices),
        min_windows_overlap);

    auto const get_window_coverages = [windows_and_reads_indices_range,
                                       n_clusters,
                                       &ringmap_data](std::size_t begin_index,
                                                      std::size_t end_index) {
      auto const window_size = end_index - begin_index;
      std::vector<unsigned> coverages(window_size, 0u);
      {
        std::set<std::uint32_t> reads_indices;
        std::ranges::for_each(
            windows_and_reads_indices_range |
                std::views::filter([&](auto const &tuple) {
                  return std::get<1>(tuple).weights.getClustersSize() ==
                         n_clusters;
                }) |
                std::views::transform(
                    [](auto const &tuple) { return std::get<2>(tuple); }),
            [&](auto const &indices) {
              reads_indices.insert(std::begin(indices), std::end(indices));
            });

        auto &&data = ringmap_data.data();
        assert(std::all_of(std::begin(reads_indices), std::end(reads_indices),
                           [nrows = data.rows_size()](auto const read_index) {
                             return read_index < nrows;
                           }));

        for (auto &&read_index : reads_indices) {
          auto &&row = data.row(read_index);
          auto const row_begin =
              std::max(row.begin_index(), static_cast<unsigned>(begin_index));
          auto const row_end =
              std::min(row.end_index(), static_cast<unsigned>(end_index));
          auto const row_size = static_cast<unsigned>(std::max(
              static_cast<int>(row_end) - static_cast<int>(row_begin), 0));

          ranges::for_each(coverages |
                               std::views::drop(row_begin - begin_index) |
                               std::views::take(row_size),
                           [](auto &&coverage) { ++coverage; });
        }
      }
      return coverages;
    };

    auto const add_result_window =
        [&transcript_result, replicate_index](results::Window &&result_window,
                                              std::size_t window_index_begin,
                                              std::size_t window_index_end) {
          if (transcript_result.windows[replicate_index])
            transcript_result.windows[replicate_index]->emplace_back(
                std::move(result_window));
          else
            transcript_result.windows[replicate_index].emplace(
                {std::move(result_window)});

          if (replicate_index == 0) {
            auto windows_size = std::size(*transcript_result.windows[0]);
            transcript_result.window_ranges.resize(windows_size);
            transcript_result.window_ranges[windows_size - 1] =
                results::WindowRange{.window_index_begin = window_index_begin,
                                     .window_index_end = window_index_end};
          }
        };

    auto const update_iters_and_region_helper =
        [&, windows_and_reads_indices_range]() {
          update_iters_and_region(windows_iter, windows_reads_indices_iter,
                                  windows_and_reads_indices_range,
                                  previous_overlapping_region_ends,
                                  min_windows_overlap);
        };

    auto filtered_windows =
        windows_and_reads_indices_range |
        std::views::filter([&](auto const &pair) {
          return std::get<1>(pair).weights.getClustersSize() == n_clusters;
        });

    if (n_clusters == 0) {
      if (args.report_uninformative()) {
        std::ranges::for_each(filtered_windows, [&](auto &&pair) {
          std::size_t window_index = std::get<0>(pair);
          auto &&window = std::get<1>(pair);
          auto const coverages = get_window_coverages(
              window.start_base, window.start_base + window.coverages.size());
          auto result_window =
              results::Window(window.start_base, window.weights, coverages);

          add_result_window(std::move(result_window), window_index,
                            window_index + 1);
        });
      }

      update_iters_and_region_helper();
      continue;
    }

    assert(n_clusters <= std::numeric_limits<std::uint8_t>::max());
    windows_merger::WindowsMerger windows_merger(
        static_cast<std::uint8_t>(n_clusters));
    auto const window_range_begin = std::get<0>(filtered_windows.front());
    auto const window_range_end = std::get<0>(filtered_windows.back()) + 1;
    assert(static_cast<std::ptrdiff_t>(window_range_end) -
               static_cast<std::ptrdiff_t>(window_range_begin) >=
           std::ranges::distance(filtered_windows));
    std::ranges::for_each(filtered_windows, [&](auto &&pair) {
      auto &&window = std::get<1>(pair);
      windows_merger.add_window(window.start_base, window.weights,
                                window.coverages);
    });

    auto const merged_window = windows_merger.merge();
    auto const coverages = get_window_coverages(merged_window.begin_index(),
                                                merged_window.end_index());
    auto result_window = results::Window(merged_window, coverages);
    add_result_window(std::move(result_window), window_range_begin,
                      window_range_end);

    update_iters_and_region_helper();
  }
}

inline static std::optional<PtbaOnReplicate>
ptba_on_replicate(std::size_t replicate_index, RingmapData const &ringmap_data,
                  Args const &args,
                  std::optional<std::ofstream> &raw_n_clusters_stream,
                  std::mutex &raw_n_clusters_stream_mutex,
                  results::Transcript &transcript_result) {
  auto const median_read_size = [&] {
    auto reads_sizes = ringmap_data.data().rows() |
                       std::views::transform([](auto &&row) {
                         assert(row.end_index() >= row.begin_index());
                         return static_cast<std::uint64_t>(row.end_index() -
                                                           row.begin_index());
                       }) |
                       more_ranges::to_vector();

    auto median_iter =
        ranges::next(ranges::begin(reads_sizes),
                     static_cast<std::ptrdiff_t>(reads_sizes.size() / 2));
    ranges::nth_element(reads_sizes, median_iter);
    return *median_iter;
  }();

  std::size_t transcript_size = ringmap_data.data().cols_size();
  const auto window_size = [&] {
    auto window_size_fraction_transcript_size =
        args.window_size_fraction_transcript_size();
    auto window_size_maybe_fraction = args.window_size();
    unsigned window_size_absolute;

    if (window_size_fraction_transcript_size > 0.) {
      window_size_absolute = static_cast<unsigned>(
          std::min(window_size_fraction_transcript_size, 1.) *
          static_cast<double>(transcript_size));
    } else if (window_size_maybe_fraction <= 1.) {
      window_size_absolute = static_cast<unsigned>(
          static_cast<double>(median_read_size) * window_size_maybe_fraction);
    } else {
      window_size_absolute =
          static_cast<unsigned>(std::round(window_size_maybe_fraction));
    }

    return std::min(window_size_absolute,
                    static_cast<unsigned>(transcript_size));
  }();
  const auto window_offset = [&] {
    auto &&window_shift_maybe_fraction = args.window_shift();
    if (window_shift_maybe_fraction < 1.) {
      return std::max(1u,
                      static_cast<unsigned>(static_cast<double>(window_size) *
                                            window_shift_maybe_fraction));
    } else {
      return static_cast<unsigned>(std::round(window_shift_maybe_fraction));
    }
  }();

  assert(window_size <= transcript_size);

  auto windows_info = WindowsInfo{
      .transcript_size = transcript_size,
      .window_size = window_size,
      .window_offset = window_offset,
  };
  auto n_windows_and_precise_offset =
      windows_info.get_n_windows_and_precise_offset();

  auto const n_windows = n_windows_and_precise_offset.n_windows;
  assert(n_windows > 0);
  std::vector<Window> windows(n_windows);
  if (n_windows > 1) {
    for (std::size_t window_index = 0; window_index < n_windows;
         ++window_index) {
      auto start_base = windows_info.get_start_base(
          n_windows_and_precise_offset, window_index);

      windows[window_index].start_base =
          static_cast<unsigned short>(start_base);
    }
  } else {
    windows[0].start_base = 0;
  }

  std::vector<unsigned> windows_n_clusters(windows.size());
  {
    auto windows_iter = std::cbegin(windows);
    auto const windows_end = std::cend(windows);
    auto &&windows_n_clusters_iter = std::begin(windows_n_clusters);

    for (; windows_iter < windows_end;
         ++windows_iter, ++windows_n_clusters_iter) {
      auto &&window = *windows_iter;
      auto &&window_n_clusters = *windows_n_clusters_iter;

      auto window_ringmap_data = ringmap_data.get_new_range(
          window.start_base, window.start_base + window_size);
      Ptba ptba(window_ringmap_data, args);

      auto const result = ptba.run();
      logger::on_debug_level(print_log_data, result.log_data, window,
                             static_cast<std::size_t>(std::distance(
                                 std::cbegin(windows), windows_iter)),
                             window_size, transcript_result, replicate_index);

      assert(result.significantIndices.size() <=
             std::numeric_limits<
                 std::remove_reference_t<decltype(window_n_clusters)>>::max());
      window_n_clusters =
          static_cast<std::remove_reference_t<decltype(window_n_clusters)>>(
              result.significantIndices.size());

      if (args.create_eigengaps_plots()) {
        auto const [eigengaps_filename, perturbed_eigengaps_filename] = [&] {
          std::array<std::string, 2> filenames;
          auto const start_base = window.start_base + 1;
          auto const end_base = window.start_base + window_size;
          std::stringstream buf;
          buf << "window_" << start_base << '-' << end_base << "_eigengaps.txt";
          filenames[0] = buf.str();

          buf.str("");
          buf << "window_" << start_base << '-' << end_base
              << "_perturbed_eigengaps.txt";
          filenames[1] = buf.str();

          return filenames;
        }();

        auto const result_dir =
            fs::path(args.eigengaps_plots_root_dir()) /
            std::format("{}_{}", transcript_result.name, replicate_index + 1);
        fs::create_directory(result_dir);
        Ptba::dumpEigenGaps(result.eigenGaps,
                            (result_dir / eigengaps_filename).c_str());
        Ptba::dumpPerturbedEigenGaps(
            result.perturbedEigenGaps,
            (result_dir / perturbed_eigengaps_filename).c_str());
      }
    }
  }

  if (raw_n_clusters_stream) {
    output_raw_n_clusters(*raw_n_clusters_stream, raw_n_clusters_stream_mutex,
                          window_size, windows, windows_n_clusters,
                          transcript_result);
    return std::nullopt;
  }

  assert(std::size(windows_n_clusters) == std::size(windows));
  return PtbaOnReplicate{
      .pre_collapsing_clusters = std::move(windows_n_clusters),
      .windows = std::move(windows),
      .window_size = window_size,
      .window_offset = window_offset,
  };
}

std::vector<PreCollapsingClusters> get_best_pre_collapsing_clusters(
    std::span<std::optional<PtbaOnReplicate>> ptba_on_replicate_results,
    std::string_view transcript_name) {
  auto const windows_size = std::size(ptba_on_replicate_results[0]->windows);
  std::vector<unsigned> window_pre_collapsing_clusters(
      std::size(ptba_on_replicate_results));
  // Take the median number of clusters, in case of even number of
  // replicates take the lowest value of the two median values.
  auto median_index = (std::size(ptba_on_replicate_results) - 1) / 2;

  return std::views::iota(0uz, windows_size) |
         std::views::transform([&](std::size_t window_index) {
           std::ranges::copy(
               ptba_on_replicate_results |
                   std::views::transform(
                       [&](const auto &ptba_on_replicate_result) {
                         assert(ptba_on_replicate_result.has_value());
                         assert(std::size(ptba_on_replicate_result
                                              ->pre_collapsing_clusters) >
                                window_index);
                         return ptba_on_replicate_result
                             ->pre_collapsing_clusters[window_index];
                       }),
               std::ranges::begin(window_pre_collapsing_clusters));

           std::ranges::sort(window_pre_collapsing_clusters);
           unsigned n_clusters = window_pre_collapsing_clusters[median_index];

           logger::on_warn_level([&] {
             if (std::ranges::any_of(window_pre_collapsing_clusters,
                                     [&](auto replicate_n_clusters) {
                                       return replicate_n_clusters !=
                                              n_clusters;
                                     })) {
               std::string all_n_clusters =
                   std::format("[{}", window_pre_collapsing_clusters[0]);
               for (auto n_clusters :
                    window_pre_collapsing_clusters | std::views::drop(1)) {
                 std::format_to(std::back_inserter(all_n_clusters), ", {}",
                                n_clusters);
               }
               std::format_to(std::back_inserter(all_n_clusters), "]");
               logger::warn(
                   "On transcript {} window {} the number of "
                   "detected clusters are not the same across replicates ({}). "
                   "The number of clusters that is going to be used is {}.",
                   transcript_name, window_index, all_n_clusters, n_clusters);
             }
           });

           float confidence =
               static_cast<float>(std::ranges::count(
                   window_pre_collapsing_clusters, n_clusters)) /
               static_cast<float>(std::size(ptba_on_replicate_results));
           return PreCollapsingClusters{
               n_clusters,
               confidence,
           };
         }) |
         std::ranges::to<std::vector>();
}

void handle_transcripts(
    std::vector<MutationMapTranscript const *> const &transcripts,
    std::vector<RingmapData *> const &ringmaps_data,
    results::Analysis &analysis_result, Args const &args,
    std::optional<std::ofstream> &raw_n_clusters_stream,
    std::mutex &raw_n_clusters_stream_mutex) {
  HandleTranscripts{
      .transcripts = transcripts,
      .ringmaps_data = ringmaps_data,
      .analysis_result = analysis_result,
      .args = args,
      .raw_n_clusters_stream = raw_n_clusters_stream,
      .raw_n_clusters_stream_mutex = raw_n_clusters_stream_mutex,
      .use_stdout = true,
      .allow_empty_patterns = false,
  }(
      [&](auto replicate_index, auto const &ringmap_data,
          auto &transcript_result) {
        return ptba_on_replicate(
            replicate_index, ringmap_data, args, raw_n_clusters_stream,
            raw_n_clusters_stream_mutex, transcript_result);
      },
      [&](unsigned short start_base, unsigned short end_base,
          std::uint8_t n_clusters,
          std::vector<arma::mat> const &replicates_covariance,
          results::Transcript const &transcript_result) {
        GraphCut graphCut(replicates_covariance);
        return graphCut.run(n_clusters, args.soft_clustering_weight_module(),
                            args.soft_clustering_initializations(),
                            args.soft_clustering_iterations(), start_base,
                            end_base, transcript_result);
      });
}

void add_detected_clusters_with_confidence(
    std::vector<WindowClustersWithConfidence>
        &detected_clusters_with_confidence,
    results::WindowRange const &region_range,
    std::vector<PreCollapsingClusters> const &pre_collapsing_clusters,
    WindowsInfo const &windows_info) {

  assert(region_range.window_index_begin <= std::size(pre_collapsing_clusters));
  assert(region_range.window_index_end <= std::size(pre_collapsing_clusters));
  auto &&region_window_clusters = std::ranges::subrange(
      std::ranges::next(
          std::ranges::cbegin(pre_collapsing_clusters),
          static_cast<std::ptrdiff_t>(region_range.window_index_begin)),
      std::ranges::next(
          std::ranges::cbegin(pre_collapsing_clusters),
          static_cast<std::ptrdiff_t>(region_range.window_index_end)));
  if (not region_window_clusters.empty()) {
    auto n_windows_and_precise_offset =
        windows_info.get_n_windows_and_precise_offset();

    auto region_window_clusters_with_start_base =
        std::views::zip(region_window_clusters,
                        std::views::iota(region_range.window_index_begin,
                                         region_range.window_index_end) |
                            std::views::transform([&](auto window_index) {
                              return windows_info.get_start_base(
                                  n_windows_and_precise_offset, window_index);
                            }));

    auto region_end_base = ([&](std::size_t start_base) {
      if (region_range.window_index_end > 0 and
          region_range.window_index_end > region_range.window_index_begin) {
        return std::min(
            start_base + windows_info.window_size,
            windows_info.get_start_base(n_windows_and_precise_offset,
                                        region_range.window_index_end - 1) +
                windows_info.window_size);
      } else {
        return std::max(start_base, windows_info.get_start_base(
                                        n_windows_and_precise_offset,
                                        region_range.window_index_begin));
      }
    });
    auto last_window_clusters = ([&]() {
      auto [first_window_clusters, start_base] =
          region_window_clusters_with_start_base.front();

      return WindowClustersWithConfidence{
          .n_clusters = first_window_clusters.n_clusters,
          .confidence = first_window_clusters.confidence,
          .start_base = std::max(
              start_base,
              windows_info.get_start_base(n_windows_and_precise_offset,
                                          region_range.window_index_begin)),
          .end_base = region_end_base(start_base),
      };
    })();

    std::ranges::for_each(
        region_window_clusters_with_start_base | std::views::drop(1),
        [&](const auto &pair) {
          auto &&[window_clusters, start_base] = pair;
          if (last_window_clusters != window_clusters) {
            detected_clusters_with_confidence.push_back(last_window_clusters);
            last_window_clusters.n_clusters = window_clusters.n_clusters;
            last_window_clusters.confidence = window_clusters.confidence;
            last_window_clusters.start_base = start_base;
          }

          last_window_clusters.end_base = region_end_base(start_base);
        });
    detected_clusters_with_confidence.push_back(last_window_clusters);
  }
}

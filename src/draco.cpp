#include "draco.hpp"
#include "args.hpp"
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
      auto const &previous_overlapping_region_end =
          previous_overlapping_region_ends[n_clusters - 1];
      if (previous_overlapping_region_end.has_value() and
          window.start_base < *previous_overlapping_region_end) {

        auto next_useful_window = std::ranges::find_if(
            windows_iter, std::cend(windows), [&](auto const &window) {
              return window.weights.getClustersSize() != n_clusters or
                     window.start_base >= *previous_overlapping_region_end;
            });
        windows_reads_indices_iter +=
            std::ranges::distance(windows_iter, next_useful_window);
        windows_iter = next_useful_window;

        continue;
      }
    }

    auto min_windows_overlap = static_cast<std::uint32_t>(
        static_cast<double>(window_size) *
        static_cast<double>(args.min_windows_overlap() / 100.));
    auto windows_and_reads_indices_range = make_windows_and_reads_indices_range(
        windows_iter, std::ranges::cend(windows), windows_reads_indices_iter,
        std::ranges::cend(windows_reads_indices), min_windows_overlap);

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
                std::views::filter([&](auto const &pair) {
                  return std::get<0>(pair).weights.getClustersSize() ==
                         n_clusters;
                }) |
                std::views::transform(
                    [](auto const &pair) { return std::get<1>(pair); }),
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
        [&transcript_result, replicate_index](results::Window &&result_window) {
          if (transcript_result.windows[replicate_index])
            transcript_result.windows[replicate_index]->emplace_back(
                std::move(result_window));
          else
            transcript_result.windows[replicate_index].emplace(
                {std::move(result_window)});
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
        std::views::transform(
            [](auto const &pair) { return std::get<0>(pair); }) |
        std::views::filter([&](auto const &window) {
          return window.weights.getClustersSize() == n_clusters;
        });

    if (n_clusters == 0) {
      if (args.report_uninformative()) {
        std::ranges::for_each(filtered_windows, [&](auto &&window) {
          auto const coverages = get_window_coverages(
              window.start_base, window.start_base + window.coverages.size());
          auto result_window =
              results::Window(window.start_base, window.weights, coverages);

          add_result_window(std::move(result_window));
        });
      }

      update_iters_and_region_helper();
      continue;
    }

    assert(n_clusters <= std::numeric_limits<std::uint8_t>::max());
    windows_merger::WindowsMerger windows_merger(
        static_cast<std::uint8_t>(n_clusters));
    std::ranges::for_each(filtered_windows, [&](auto &&window) {
      windows_merger.add_window(window.start_base, window.weights,
                                window.coverages);
    });

    auto const merged_window = windows_merger.merge();
    auto const coverages = get_window_coverages(merged_window.begin_index(),
                                                merged_window.end_index());
    auto result_window = results::Window(merged_window, coverages);
    add_result_window(std::move(result_window));

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
    auto &&window_size = args.window_size();
    if (window_size <= 0) {
      window_size = static_cast<unsigned>(
          static_cast<double>(median_read_size) * args.window_size_fraction());
    }

    return std::min(window_size, static_cast<unsigned>(transcript_size));
  }();
  const auto window_offset = [&] {
    auto &&window_shift = args.window_shift();
    if (window_shift > 0) {
      return window_shift;
    } else {
      return static_cast<unsigned>(static_cast<double>(window_size) *
                                   args.window_shift_fraction());
    }
  }();

  assert(window_size <= transcript_size);

  std::size_t n_windows = (transcript_size - window_size) / window_offset + 1;
  if (n_windows * window_offset + window_size < transcript_size)
    ++n_windows;

  assert(n_windows > 0);
  std::vector<Window> windows(n_windows);
  if (n_windows > 1) {
    auto const window_precise_offset =
        static_cast<double>(transcript_size - window_size) /
        static_cast<double>(n_windows - 1);
    for (std::size_t window_index = 0; window_index < n_windows;
         ++window_index) {
      auto start_base = static_cast<std::size_t>(std::round(
          static_cast<double>(window_index) * window_precise_offset));
      if (start_base + window_size > transcript_size)
        start_base = transcript_size - window_size;

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
  };
}

std::vector<unsigned> get_best_pre_collapsing_clusters(
    std::span<std::optional<PtbaOnReplicate>> ptba_on_replicate_results) {
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
           return window_pre_collapsing_clusters[median_index];
         }) |
         std::ranges::to<std::vector>();
}

void handle_transcripts(
    std::vector<MutationMapTranscript const *> const &transcripts,
    std::vector<RingmapData *> const &ringmaps_data,
    results::Analysis &analysis_result, Args const &args,
    std::optional<std::ofstream> &raw_n_clusters_stream,
    std::mutex &raw_n_clusters_stream_mutex) {
  auto const &first_transcript = *transcripts[0];
  if (std::ranges::any_of(ringmaps_data, [](auto const &ringmap_data) {
        return ringmap_data->data().rows_size() == 0;
      })) {
    std::cout << "\x1b[2K\r[+] Skipping transcript " << first_transcript.getId()
              << " (no reads)" << std::endl;
    return;
  }
  std::cout << "\x1b[2K\r[+] Analyzing transcript " << first_transcript.getId()
            << std::flush;

  results::Transcript transcript_result(std::size(transcripts));
  transcript_result.name = first_transcript.getId();
  std::ranges::copy(transcripts |
                        std::views::transform([](auto const &transcript) {
                          return transcript->getReadsSize();
                        }),
                    std::ranges::begin(transcript_result.reads));
  transcript_result.sequence = first_transcript.getSequence();
  assert(not transcript_result.name.empty());

  auto ptba_on_replicate_results =
      std::vector<std::optional<PtbaOnReplicate>>(std::size(transcripts));
  tbb::parallel_for(
      0uz, std::size(transcripts), [&](std::size_t replicate_index) {
        auto const &ringmap_data = *ringmaps_data[replicate_index];
        ptba_on_replicate_results[replicate_index] = ptba_on_replicate(
            replicate_index, ringmap_data, args, raw_n_clusters_stream,
            raw_n_clusters_stream_mutex, transcript_result);
      });

  if (not std::ranges::all_of(ptba_on_replicate_results,
                              [](const auto &ptba_on_replicate_result) {
                                return ptba_on_replicate_result.has_value();
                              })) {
    return;
  }

  if (std::ranges::any_of(
          ptba_on_replicate_results | std::views::drop(1),
          [&](auto &ptba_on_replicate_result) {
            return std::size(ptba_on_replicate_result->windows) !=
                       std::size(ptba_on_replicate_results[0]->windows) or
                   ptba_on_replicate_result->window_size !=
                       ptba_on_replicate_results[0]->window_size;
          })) {
    throw std::runtime_error(
        std::format("the number of windows for transcript {} is incoherent",
                    first_transcript.getId()));
  };

  // auto &windows = ptba_on_replicate_result->windows;
  auto const window_size = ptba_on_replicate_results[0]->window_size;
  auto const windows_size = std::size(ptba_on_replicate_results[0]->windows);

  auto const pre_collapsing_clusters =
      get_best_pre_collapsing_clusters(ptba_on_replicate_results);

  std::vector<unsigned> windows_n_clusters;
  std::vector<std::optional<unsigned>> windows_max_clusters_constraints(
      windows_size, std::nullopt);

  for (bool stop = false; not stop;) {
    stop = true;
    std::ranges::fill(transcript_result.windows, std::nullopt);

    windows_n_clusters = pre_collapsing_clusters;
    ([&]() {
      auto windows_n_clusters_iter = std::ranges::begin(windows_n_clusters);
      const auto windows_n_clusters_end = std::ranges::end(windows_n_clusters);

      auto windows_max_clusters_constraints_iter =
          std::ranges::begin(windows_max_clusters_constraints);
      const auto windows_max_clusters_constraints_end =
          std::ranges::end(windows_max_clusters_constraints);

      for (; windows_n_clusters_iter < windows_n_clusters_end and
             windows_max_clusters_constraints_iter <
                 windows_max_clusters_constraints_end;
           ++windows_n_clusters_iter, ++windows_max_clusters_constraints_iter) {
        auto &&window_n_clusters = *windows_n_clusters_iter;
        auto &&constraint = *windows_max_clusters_constraints_iter;

        if (constraint) {
          window_n_clusters = std::min(window_n_clusters, *constraint);
        }
      }
    })();

    if (args.set_uninformative_clusters_to_surrounding()) {
      set_uninformative_clusters_to_surrounding(
          windows_n_clusters, windows_max_clusters_constraints);

      [[maybe_unused]] constexpr auto const zero_clusters =
          [](auto n_clusters) { return n_clusters == 0; };
      assert(ranges::none_of(windows_n_clusters, zero_clusters) or
             ranges::all_of(windows_n_clusters, zero_clusters));
    }

    if (args.max_collapsing_windows() > 0)
      collapse_outlayer_clusters(windows_n_clusters,
                                 windows_max_clusters_constraints, args);

    if (args.set_all_uninformative_to_one()) {
      if (ranges::all_of(windows_n_clusters,
                         [](auto n_clusters) { return n_clusters == 0; })) {
        ranges::fill(windows_n_clusters, 1u);
      }
    }

    auto const &first_result_windows = ptba_on_replicate_results[0]->windows;
    if (not std::ranges::all_of(
            std::views::iota(0uz, windows_size), [&](auto window_index) {
              auto const &first_result_window =
                  first_result_windows[window_index];
              return std::ranges::all_of(
                  ptba_on_replicate_results | std::views::drop(1),
                  [&](auto const &ptba_on_replicate_result) {
                    auto const &window =
                        ptba_on_replicate_result->windows[window_index];
                    return window.start_base ==
                               first_result_window.start_base &&
                           window.weights.getClustersSize() ==
                               first_result_window.weights.getClustersSize();
                  });
            })) {
      throw std::runtime_error(
          std::format("incoherent windows between replicates for transcript {}",
                      first_transcript.getId()));
    };

    std::vector<std::vector<unsigned>> replicates_windows_reads_indices_vec(
        windows_size * std::size(ringmaps_data));
    auto replicates_windows_reads_indices =
        [&](std::size_t replicate_index,
            std::size_t window_index) -> std::vector<unsigned> & {
      return replicates_windows_reads_indices_vec[replicate_index *
                                                      windows_size +
                                                  window_index];
    };
    auto replicates_windows_reads_indices_replicate =
        [&](std::size_t replicate_index) {
          return std::span(
              std::ranges::next(
                  std::ranges::begin(replicates_windows_reads_indices_vec),
                  static_cast<std::ptrdiff_t>(replicate_index * windows_size)),
              windows_size);
        };

    {
      auto windows_n_clusters_iter = std::cbegin(windows_n_clusters);

      for (auto window_index = 0uz; window_index < windows_size;
           ++window_index, ++windows_n_clusters_iter) {

        assert(*windows_n_clusters_iter <=
               std::numeric_limits<std::uint8_t>::max());
        auto n_clusters = static_cast<std::uint8_t>(*windows_n_clusters_iter);

        auto replicates_filtered_data =
            std::views::zip(ringmaps_data, ptba_on_replicate_results,
                            std::views::iota(0uz)) |
            std::views::transform([&](auto &&tuple) {
              auto &&[ringmap_data, ptba_on_replicate_result, replicate_index] =
                  std::move(tuple);
              auto &window = ptba_on_replicate_result->windows[window_index];
              std::vector<unsigned> window_reads_indices;
              auto window_ringmap_data = ringmap_data->get_new_range(
                  window.start_base, window.start_base + window_size,
                  &window_reads_indices);

              assert(window_reads_indices.size() ==
                     window_ringmap_data.data().rows_size());
              assert(std::is_sorted(std::begin(window_reads_indices),
                                    std::end(window_reads_indices)));
              assert(std::unique(std::begin(window_reads_indices),
                                 std::end(window_reads_indices)) ==
                     std::end(window_reads_indices));
              window.coverages = window_ringmap_data.getBaseCoverages();

              assert(window_reads_indices.size() ==
                     window_ringmap_data.data().rows_size());

              replicates_windows_reads_indices(replicate_index, window_index) =
                  std::move(window_reads_indices);

              return window_ringmap_data;
            }) |
            std::views::as_rvalue | std::ranges::to<std::vector>();

        RingmapData::filter_bases_on_replicates(replicates_filtered_data);
        for (auto &filtered_data : replicates_filtered_data) {
          filtered_data.filterReads();
        }
        RingmapData::filter_bases_on_replicates(replicates_filtered_data);

        typename RingmapData::clusters_pattern_type patterns;
        auto const &first_result_window = first_result_windows[window_index];
        for (;;) {
          if (n_clusters > 1 and
              std::ranges::any_of(replicates_filtered_data,
                                  [](auto const &filtered_data) {
                                    return filtered_data.data().rows_size() > 0;
                                  })) {
            auto replicates_covariance =
                replicates_filtered_data |
                std::views::filter([](const auto &filtered_data) {
                  return filtered_data.data().rows_size() > 0;
                }) |
                std::views::transform([](const auto &filtered_data) {
                  return filtered_data.data().covariance(
                      filtered_data.getBaseWeights());
                }) |
                std::views::as_rvalue | std::ranges::to<std::vector>();

            GraphCut graphCut(replicates_covariance);

            assert(first_result_window.start_base <=
                   std::numeric_limits<
                       decltype(first_result_window.start_base)>::max() -
                       window_size);
            auto graphCutResults = graphCut.run(
                n_clusters, args.soft_clustering_weight_module(),
                args.soft_clustering_initializations(),
                args.soft_clustering_iterations(),
                first_result_window.start_base,
                first_result_window.start_base +
                    static_cast<decltype(first_result_window.start_base)>(
                        window_size),
                transcript_result);

            std::ranges::for_each(
                std::views::zip(replicates_filtered_data,
                                ptba_on_replicate_results),
                [&](auto pair) {
                  auto &&[filtered_data, ptba_on_replicate_result] = pair;
                  auto clusters =
                      filtered_data.getUnfilteredWeights(graphCutResults);

                  assert(clusters.getElementsSize() == window_size);
                  ptba_on_replicate_result->windows[window_index].weights =
                      std::move(clusters);
                });

            break;
          } else {
            for (auto &ptba_on_replicate_result : ptba_on_replicate_results) {
              ptba_on_replicate_result->windows[window_index].weights =
                  WeightedClusters(window_size, n_clusters);
            }
            break;
          }
        }
      }
    }

    std::ranges::for_each(
        std::views::zip(std::views::iota(0uz), ptba_on_replicate_results,
                        ringmaps_data),
        [&](auto &&tuple) {
          auto &&[replicate_index, ptba_on_replicate_result, ringmap_data_ptr] =
              std::move(tuple);
          auto const &ringmap_data = *ringmap_data_ptr;
          auto const &windows = ptba_on_replicate_result->windows;
          auto const &replicate_windows_reads_indices =
              replicates_windows_reads_indices_replicate(replicate_index);

          assert(std::size(windows) ==
                 std::size(replicate_windows_reads_indices));
          merge_windows_and_add_window_results(
              windows, replicate_windows_reads_indices, ringmap_data,
              transcript_result, window_size, replicate_index, args);
        });

    if (std::ranges::any_of(transcript_result.windows | std::views::drop(1),
                            [&](auto const &window) {
                              return (window.has_value() !=
                                      transcript_result.windows[0].has_value());
                            })) {
      std::cerr << "\nInconsistency detected between merged windows"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }

    if (transcript_result.windows[0]) {
      auto const &first_replicate_result_windows =
          *transcript_result.windows[0];
      auto replicates_splitted_ringmaps =
          std::views::zip(ringmaps_data, transcript_result.windows) |
          std::views::transform([](auto &&pair) {
            auto &&[ringmap_data, result_windows] = std::move(pair);
            auto splitted_ringmaps =
                RingmapData(*ringmap_data).split_into_windows(*result_windows);
            assert(std::size(splitted_ringmaps) == std::size(*result_windows));
            return splitted_ringmaps;
          }) |
          std::views::as_rvalue | std::ranges::to<std::vector>();

      std::vector<RingmapData> filtered_ringmaps(std::size(ringmaps_data));
      for (auto window_index :
           std::views::iota(0uz, std::size(first_replicate_result_windows))) {
        if (std::ranges::any_of(
                transcript_result.windows | std::views::drop(1) |
                    std::views::transform(
                        [&](const auto &replicate_windows) -> decltype(auto) {
                          return (*replicate_windows)[window_index];
                        }),
                [&](auto const &window) {
                  return window.weighted_clusters.getClustersSize() !=
                         first_replicate_result_windows[window_index]
                             .weighted_clusters.getClustersSize();
                })) {
          std::cerr << "\nInconsistency between weighted clusters sizes across "
                       "replicates"
                    << std::endl;
          std::exit(EXIT_FAILURE);
        }

        if (first_replicate_result_windows[window_index]
                .weighted_clusters.getClustersSize() == 0) {
          continue;
        }

        for (auto &&[splitted_ringmaps, windows] : std::views::zip(
                 replicates_splitted_ringmaps, transcript_result.windows)) {
          auto &window = (*windows)[window_index];
          auto const &ringmap = splitted_ringmaps[window_index];
          window.assignments.resize(ringmap.data().rows_size());
          ranges::fill(window.assignments, std::int8_t(-1));
        }

        std::ranges::copy(
            replicates_splitted_ringmaps |
                std::views::transform(
                    [&](auto const &replicate_splitted_ringmaps)
                        -> decltype(auto) {
                      return replicate_splitted_ringmaps[window_index];
                    }),
            std::ranges::begin(filtered_ringmaps));
        RingmapData::filter_bases_on_replicates(filtered_ringmaps);
        for (auto &filtered_ringmap : filtered_ringmaps) {
          filtered_ringmap.filterReads();
        }
        RingmapData::filter_bases_on_replicates(filtered_ringmaps);

        std::ranges::for_each(
            std::views::zip(
                std::views::iota(0uz), filtered_ringmaps,
                replicates_splitted_ringmaps,
                transcript_result.windows |
                    std::views::transform(
                        [&](auto &replicate_windows) -> decltype(auto) {
                          return (*replicate_windows)[window_index];
                        })) |
                std::views::take_while([&](auto const &) { return stop; }),
            [&](auto &&tuple) {
              auto replicate_index = std::get<0>(tuple);
              auto const &filtered_ringmap = std::get<1>(tuple);
              auto const &ringmap = std::get<2>(tuple)[window_index];
              auto &window = std::get<3>(tuple);

              auto &&fractions_result = filtered_ringmap.fractionReadsByWeights(
                  window.weighted_clusters, window_size,
                  args.skip_ambiguous_assignments());
              std::tie(window.fractions, window.patterns, std::ignore) =
                  std::move(fractions_result);
              assert(window.fractions.size() > 1 or window.fractions.empty() or
                     window.fractions[0] >= 0.01);

              bool const redundand_patterns = [&] {
                auto patterns_iter = std::cbegin(*window.patterns);
                auto const patterns_end = std::cend(*window.patterns);

                for (; patterns_iter < patterns_end; ++patterns_iter) {
                  auto &&cur_pattern = *patterns_iter;
                  auto const begin_cur_pattern = std::cbegin(cur_pattern);
                  auto const end_cur_pattern = std::cend(cur_pattern);

                  if (std::any_of(std::next(patterns_iter), patterns_end,
                                  [&](auto &&next_pattern) {
                                    return std::equal(begin_cur_pattern,
                                                      end_cur_pattern,
                                                      std::cbegin(next_pattern),
                                                      std::cend(next_pattern));
                                  })) {
                    return true;
                  }
                }

                return false;
              }();

              if (redundand_patterns or
                  ranges::any_of(
                      window.fractions,
                      [min_cluster_fraction =
                           args.minimum_cluster_fraction()](auto &&fraction) {
                        return fraction < min_cluster_fraction;
                      })) {
                stop = false;
                auto const result_window_begin = window.begin_index;
                auto const result_window_end = window.end_index;
                logger::debug(
                    "Transcript {}, replicate {}, window {} (bases {}-{}), "
                    "reducing number of clusters from {} to {} because a "
                    "redundant weights pattern is found",
                    transcript_result.name, replicate_index + 1,
                    window_index + 1, result_window_begin + 1,
                    result_window_end, window.fractions.size(),
                    window.fractions.size() - 1);
                assert(window.fractions.size() > 1);
                auto const new_clusters_constraint =
                    static_cast<unsigned>(window.fractions.size() - 1);

                for (auto &replicate_result : ptba_on_replicate_results) {
                  for (auto &&[window, window_constraint] :
                       std::views::zip(replicate_result->windows,
                                       windows_max_clusters_constraints)) {
                    if (window.start_base >= result_window_begin and
                        window.start_base + window_size <= result_window_end) {

                      window_constraint = new_clusters_constraint;
                    }
                  }
                }
              }

              if (not stop)
                return;

              *window.patterns =
                  filtered_ringmap.remapPatterns(*window.patterns);

              assign_reads_to_clusters(window,
                                       std::move(std::get<2>(fractions_result)),
                                       ringmap, filtered_ringmap);

              if (not args.assignments_dump_directory().empty()) {
                std::optional<std::size_t> usable_replicate_index;
                if (std::size(ptba_on_replicate_results) > 1) {
                  usable_replicate_index.emplace(replicate_index);
                }
                dump_assignments(transcript_result, window, ringmap,
                                 usable_replicate_index,
                                 args.assignments_dump_directory());
              }

              if (std::all_of(std::cbegin(*window.patterns),
                              std::cend(*window.patterns), [](auto &&pattern) {
                                return std::all_of(
                                    std::cbegin(pattern), std::cend(pattern),
                                    [](auto &&value) { return value == 0; });
                              })) {
                window.patterns = std::nullopt;
                window.bases_coverages = std::nullopt;
              }
            });
      }
    }
  }

  analysis_result.addTranscript(std::move(transcript_result));
}

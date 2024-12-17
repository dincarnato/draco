#include "draco.hpp"
#include "args.hpp"
#include "graph_cut.hpp"
#include "logger.hpp"
#include "mutation_map.hpp"
#include "mutation_map_writer.hpp"
#include "parallel/blocking_queue.hpp"
#include "ptba.hpp"
#include "results/analysis.hpp"
#include "results/transcript.hpp"
#include "results/window.hpp"
#include "ringmap_data.hpp"
#include "to_vector.hpp"
#include "windows_merger.hpp"

#include <algorithm>
#include <charconv>
#include <format>
#include <iostream>
#include <iterator>
#include <mutex>
#include <oneapi/tbb/info.h>
#include <oneapi/tbb/parallel_pipeline.h>
#include <optional>
#include <ranges>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>

#include <armadillo>
#include <tbb/global_control.h>

#include <omp.h>

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#elif __has_include(<experimental/filesystem>)
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#error "Missing filesystem header"
#endif

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
get_windows_spans(std::vector<Window> const &windows,
                  std::vector<unsigned> const &windows_n_clusters,
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

  auto const windows_size = windows.size();
  assert(windows_size == windows_n_clusters.size());
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
void windows_span_expander(std::vector<Window> const &windows,
                           std::vector<unsigned> &windows_n_clusters,
                           std::vector<std::optional<unsigned>> const
                               &windows_max_clusters_constraints,
                           F1 &&loop_start_check, F2 &&same_clusters_check,
                           F3 &&different_clusters_check) noexcept(false) {
  auto const windows_size = windows.size();
  auto windows_spans = get_windows_spans(windows, windows_n_clusters,
                                         windows_max_clusters_constraints);

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
    std::fill(std::next(windows_n_clusters_begin, windows_span.begin),
              std::next(windows_n_clusters_begin, windows_span.end),
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

void collapse_outlayer_clusters(std::vector<Window> const &windows,
                                std::vector<unsigned> &windows_n_clusters,
                                std::vector<std::optional<unsigned>> const
                                    &windows_max_clusters_constraints,
                                Args const &args) noexcept(false) {
  auto const max_collapsing_windows = args.max_collapsing_windows();
  auto const min_surrounding_windows_size = args.min_surrounding_windows_size();
  windows_span_expander(
      windows, windows_n_clusters, windows_max_clusters_constraints,
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
    std::vector<Window> const &windows,
    std::vector<unsigned> &windows_n_clusters,
    std::vector<std::optional<unsigned>> const
        &windows_max_clusters_constraints) noexcept(false) {
  windows_span_expander(
      windows, windows_n_clusters, windows_max_clusters_constraints,
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

    for (std::size_t cluster_index = 0; cluster_index < clusters_size;
         ++cluster_index) {
      auto &&cluster_assignments = clusters_assignments[cluster_index];

      for (auto filtered_read_index : cluster_assignments) {
        auto original_read_index = original_indices_map[filtered_read_index];
        window.assignments[original_read_index] = cluster_index;

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
                      std::string_view assignments_dump_directory) {
  auto mm_basename = std::format("{}_{}-{}", transcript.name,
                                 window.begin_index, window.end_index - 1);
  auto mm_path = std::filesystem::path{assignments_dump_directory} /
                 std::format("{}.mm", mm_basename);
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

void handle_transcript(MutationMapTranscript const &transcript,
                       RingmapData &ringmapData,
                       results::Analysis &analysisResult, Args const &args,
                       std::optional<std::ofstream> &raw_n_clusters_stream,
                       std::mutex &raw_n_clusters_stream_mutex) {
  if (ringmapData.data().rows_size() == 0) {
    std::cout << "\x1b[2K\r[+] Skipping transcript " << transcript.getId()
              << " (no reads)" << std::endl;
    return;
  }
  std::cout << "\x1b[2K\r[+] Analyzing transcript " << transcript.getId()
            << std::flush;

  results::Transcript transcriptResult;
  transcriptResult.name = transcript.getId();
  transcriptResult.reads = transcript.getReadsSize();
  transcriptResult.sequence = transcript.getSequence();
  assert(not transcriptResult.name.empty());

  auto const median_read_size = [&] {
    auto reads_sizes = ringmapData.data().rows() |
                       std::views::transform([](auto &&row) {
                         assert(row.end_index() >= row.begin_index());
                         return static_cast<std::uint64_t>(row.end_index() -
                                                           row.begin_index());
                       }) |
                       more_ranges::to_vector();

    auto median_iter =
        ranges::next(ranges::begin(reads_sizes), reads_sizes.size() / 2);
    ranges::nth_element(reads_sizes, median_iter);
    return *median_iter;
  }();

  std::size_t transcript_size = ringmapData.data().cols_size();
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

      auto window_ringmap_data = ringmapData.get_new_range(
          window.start_base, window.start_base + window_size);
      Ptba ptba(window_ringmap_data, args);

      auto const result = ptba.run();
      logger::on_debug_level(print_log_data, result.log_data, window,
                             static_cast<std::size_t>(std::distance(
                                 std::cbegin(windows), windows_iter)),
                             window_size, transcriptResult);
      if (args.create_eigengaps_plots()) {
        window_n_clusters = result.significantIndices.size();
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
            fs::path(args.eigengaps_plots_root_dir()) / transcriptResult.name;
        fs::create_directory(result_dir);
        Ptba::dumpEigenGaps(result.eigenGaps,
                            (result_dir / eigengaps_filename).c_str());
        Ptba::dumpPerturbedEigenGaps(
            result.perturbedEigenGaps,
            (result_dir / perturbed_eigengaps_filename).c_str());
      } else {
        window_n_clusters = result.significantIndices.size();
      }
    }
  }

  if (raw_n_clusters_stream) {
    std::stringstream raw_n_clusters_data;

    auto windows_iter = std::cbegin(windows);
    auto windows_end = std::cend(windows);
    auto windows_n_clusters_iter = std::cbegin(windows_n_clusters);
    assert(
        std::distance(windows_iter, windows_end) ==
        std::distance(windows_n_clusters_iter, std::cend(windows_n_clusters)));

    for (; windows_iter < windows_end;
         ++windows_iter, ++windows_n_clusters_iter) {
      auto &&window = *windows_iter;
      auto n_clusters = *windows_n_clusters_iter;
      raw_n_clusters_data << transcriptResult.name << '\t' << window.start_base
                          << '\t' << window.start_base + window_size << '\t'
                          << n_clusters << '\n';
    }

    {
      std::lock_guard<std::mutex> lock(raw_n_clusters_stream_mutex);
      *raw_n_clusters_stream << raw_n_clusters_data.rdbuf();
    }
    analysisResult.addTranscript(std::move(transcriptResult));
    return;
  }

  auto const pre_collapsing_clusters = std::move(windows_n_clusters);
  windows_n_clusters.clear();
  std::vector<std::optional<unsigned>> windows_max_clusters_constraints(
      windows.size(), std::nullopt);

  for (bool stop = false; not stop;) {
    stop = true;
    transcriptResult.windows = std::nullopt;

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
          windows, windows_n_clusters, windows_max_clusters_constraints);

      [[maybe_unused]] constexpr auto const zero_clusters =
          [](auto n_clusters) { return n_clusters == 0; };
      assert(ranges::none_of(windows_n_clusters, zero_clusters) or
             ranges::all_of(windows_n_clusters, zero_clusters));
    }

    if (args.max_collapsing_windows() > 0)
      collapse_outlayer_clusters(windows, windows_n_clusters,
                                 windows_max_clusters_constraints, args);

    if (args.set_all_uninformative_to_one()) {
      if (ranges::all_of(windows_n_clusters,
                         [](auto n_clusters) { return n_clusters == 0; })) {
        ranges::fill(windows_n_clusters, 1u);
      }
    }

    std::vector<std::vector<unsigned>> windows_reads_indices;
    windows_reads_indices.reserve(windows.size());

    {
      auto windows_iter = std::begin(windows);
      auto const windows_end = std::end(windows);
      auto windows_n_clusters_iter = std::cbegin(windows_n_clusters);

      for (; windows_iter < windows_end;
           ++windows_iter, ++windows_n_clusters_iter) {

        auto &&window = *windows_iter;
        auto n_clusters = *windows_n_clusters_iter;

        std::vector<unsigned> window_reads_indices;
        auto window_ringmap_data = ringmapData.get_new_range(
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
        windows_reads_indices.emplace_back(std::move(window_reads_indices));

        auto filtered_data = window_ringmap_data;
        filtered_data.filterBases();
        filtered_data.filterReads();
        filtered_data.filterBases();

        typename RingmapData::clusters_pattern_type patterns;
        for (;;) {
          if (n_clusters > 1 and filtered_data.data().rows_size() > 0) {
            auto covariance =
                filtered_data.data().covariance(filtered_data.getBaseWeights());
            GraphCut graphCut(covariance);

            auto graphCutResults = graphCut.run(
                n_clusters, args.soft_clustering_weight_module(),
                args.soft_clustering_initializations(),
                args.soft_clustering_iterations(), window.start_base,
                window.start_base + window_size, transcriptResult);
            auto clusters =
                filtered_data.getUnfilteredWeights(std::move(graphCutResults));

            assert(clusters.getElementsSize() == window_size);
            window.weights = std::move(clusters);
            break;
          } else {
            window.weights = WeightedClusters(window_size, n_clusters);
            break;
          }
        }
      }
    }

    {
      auto window_iter = std::begin(windows);
      auto window_reads_indices_iter = std::begin(windows_reads_indices);
      // std::ofstream merge_stream("merge.txt");
      while (window_iter != std::end(windows)) {
        unsigned const n_clusters = window_iter->weights.getClustersSize();
        auto last_window = std::find_if(
            std::next(window_iter), std::end(windows),
            [n_clusters](auto &&window) {
              return window.weights.getClustersSize() != n_clusters;
            });

        auto const last_window_reads_indices = std::next(
            window_reads_indices_iter, std::distance(window_iter, last_window));

        auto const get_window_coverages = [&window_reads_indices_iter,
                                           &last_window_reads_indices,
                                           &ringmapData](
                                              std::size_t begin_index,
                                              std::size_t end_index) {
          auto const window_size = end_index - begin_index;
          std::vector<unsigned> coverages(window_size, 0u);
          {
            std::set<std::size_t> reads_indices;
            std::for_each(window_reads_indices_iter, last_window_reads_indices,
                          [&](auto const &indices) {
                            reads_indices.insert(std::begin(indices),
                                                 std::end(indices));
                          });

            auto &&data = ringmapData.data();
            assert(
                std::all_of(std::begin(reads_indices), std::end(reads_indices),
                            [nrows = data.rows_size()](auto const read_index) {
                              return read_index < nrows;
                            }));

            for (auto &&read_index : reads_indices) {
              auto &&row = data.row(read_index);
              auto const row_begin = std::max(
                  row.begin_index(), static_cast<unsigned>(begin_index));
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

        auto const add_result_window = [&transcriptResult](
                                           results::Window &&result_window) {
          if (transcriptResult.windows)
            transcriptResult.windows->emplace_back(std::move(result_window));
          else
            transcriptResult.windows.emplace({std::move(result_window)});
        };

        if (n_clusters == 0) {
          if (args.report_uninformative()) {
            std::for_each(window_iter, last_window, [&](auto &&window) {
              auto const coverages = get_window_coverages(
                  window.start_base,
                  window.start_base + window.coverages.size());
              auto result_window =
                  results::Window(window.start_base, window.weights, coverages);

              add_result_window(std::move(result_window));
            });
          }

          window_iter = last_window;
          window_reads_indices_iter = last_window_reads_indices;
          continue;
        }

        windows_merger::WindowsMerger windows_merger(n_clusters);
        std::for_each(window_iter, last_window, [&](auto &&window) {
          windows_merger.add_window(window.start_base, window.weights,
                                    window.coverages);
        });

        auto const merged_window = windows_merger.merge();
        auto const coverages = get_window_coverages(merged_window.begin_index(),
                                                    merged_window.end_index());
        auto result_window = results::Window(merged_window, coverages);
        add_result_window(std::move(result_window));

        window_iter = last_window;
        window_reads_indices_iter = last_window_reads_indices;
      }
    }

    if (transcriptResult.windows) {
      auto &result_windows = *transcriptResult.windows;
      auto splitted_ringmaps =
          RingmapData(ringmapData).split_into_windows(result_windows);

      assert(splitted_ringmaps.size() == result_windows.size());

      auto windows_iter = std::begin(result_windows);
      auto const windows_end = std::end(result_windows);
      auto splitted_ringmaps_iter = std::begin(splitted_ringmaps);
      for (; windows_iter < windows_end;
           ++windows_iter, ++splitted_ringmaps_iter) {
        auto &window = *windows_iter;
        if (window.weighted_clusters.getClustersSize() == 0) {
          continue;
        }

        auto &ringmap = *splitted_ringmaps_iter;

        window.assignments.resize(ringmap.data().rows_size());
        ranges::fill(window.assignments, std::int8_t(-1));

        auto filteredRingmap = ringmap;
        filteredRingmap.filterBases();
        filteredRingmap.filterReads();
        filteredRingmap.filterBases();

        auto &&fractions_result = filteredRingmap.fractionReadsByWeights(
            window.weighted_clusters, window_size,
            args.skip_ambiguous_assignments());
        std::tie(window.fractions, window.patterns, std::ignore) =
            std::move(fractions_result);
        assert(window.fractions.size() > 1 or window.fractions.empty() or
               window.fractions[0] >= 0.01);

        bool const redundandPatterns = [&] {
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

        if (redundandPatterns or
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
              "Transcript {}, window {} (bases {}-{}), reducing number "
              "of clusters from {} to {} because a redundant weights "
              "pattern is found",
              transcriptResult.name,
              std::distance(std::begin(result_windows), windows_iter) + 1,
              result_window_begin + 1, result_window_end,
              window.fractions.size(), window.fractions.size() - 1);
          assert(window.fractions.size() > 1);
          auto const new_clusters_constraint =
              static_cast<unsigned>(window.fractions.size() - 1);

          ([&]() {
            auto windows_iter = std::ranges::begin(windows);
            const auto windows_end = std::ranges::end(windows);

            auto windows_max_clusters_constraints_iter =
                std::ranges::begin(windows_max_clusters_constraints);
            const auto windows_max_clusters_constraints_end =
                std::ranges::end(windows_max_clusters_constraints);

            for (; windows_iter < windows_end and
                   windows_max_clusters_constraints_iter <
                       windows_max_clusters_constraints_end;
                 ++windows_iter, ++windows_max_clusters_constraints_iter) {
              auto &&window = *windows_iter;
              auto &&window_constraint = *windows_max_clusters_constraints_iter;

              if (window.start_base >= result_window_begin and
                  window.start_base + window_size <= result_window_end) {

                window_constraint = new_clusters_constraint;
              }
            }
          })();
        }

        if (not stop)
          continue;

        *window.patterns = filteredRingmap.remapPatterns(*window.patterns);

        assign_reads_to_clusters(window,
                                 std::move(std::get<2>(fractions_result)),
                                 ringmap, filteredRingmap);

        if (not args.assignments_dump_directory().empty()) {
          dump_assignments(transcriptResult, window, ringmap,
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
      }
    }
  }

  analysisResult.addTranscript(std::move(transcriptResult));
}

int main(int argc, char *argv[]) {
  auto const args = Args(argc, argv);

  auto &&log_level = args.log_level();
  if (not log_level.empty()) {
    auto level = logger::parse_level(log_level);
    if (level.has_value()) {
      logger::instance.set_level(*level);
    } else {
      logger::warn("Invalid debug level \"{}\"", log_level);
    }
  }

  if (not args.assignments_dump_directory().empty()) {
    logger::debug("Creating directory {}", args.assignments_dump_directory());
    std::filesystem::create_directories(args.assignments_dump_directory());
  }

  std::cout << "\n[+] Starting DRACO analysis. This might take a while...\n";

  // Disabling OMP, we parallelize a higher level
  omp_set_num_threads(1);

  MutationMap mutationMap(args.mm_filename());
  results::Analysis analysisResult(args.output_filename());

  analysisResult.filename = args.mm_filename();

  if (args.create_eigengaps_plots()) {
    fs::create_directory(fs::path(args.eigengaps_plots_root_dir()));
  }

  parallel::blocking_queue<std::pair<MutationMapTranscript, RingmapData>> queue(
      10);
  std::thread reader([&] {
    RingmapData::enqueueRingmapsFromMutationMap(mutationMap, queue, args);
  });
  /*
  std::thread reader([&] {
    auto transcriptIter = std::next(std::begin(mutationMap), 13);
    queue.emplace(*transcriptIter, RingmapData(*transcriptIter, args));
    ++transcriptIter;
    queue.emplace(*transcriptIter, RingmapData(*transcriptIter, args));
    queue.finish();
  });
  */

  const auto max_allowed_parallelism = [&] {
    auto n_processors = args.n_processors();
    if (n_processors == 0) {
      n_processors = tbb::info::default_concurrency();
    }

    return static_cast<std::size_t>(n_processors);
  }();

  tbb::global_control tbb_global_control(
      tbb::global_control::max_allowed_parallelism, max_allowed_parallelism);

  std::optional<std::ofstream> raw_n_clusters_stream;
  std::mutex raw_n_clusters_stream_mutex;
  if (!args.output_raw_n_clusters().empty()) {
    raw_n_clusters_stream =
        std::ofstream(args.output_raw_n_clusters(),
                      std::ios_base::out | std::ios_base::trunc);
  }

  tbb::parallel_pipeline(
      max_allowed_parallelism,
      tbb::make_filter<void, void>(
          tbb::filter_mode::parallel, [&](tbb::flow_control &flow_control) {
            auto poppedData = queue.pop();
            if (not poppedData) {
              flow_control.stop();
              return;
            }

            auto const &transcript = std::get<0>(*poppedData);
            auto &ringmapData = std::get<1>(*poppedData);

            handle_transcript(transcript, ringmapData, analysisResult, args,
                              raw_n_clusters_stream,
                              raw_n_clusters_stream_mutex);
          }));

  reader.join();
  std::cout << "\n[+] All done.\n\n";
}

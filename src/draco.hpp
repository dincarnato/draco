#pragma once

#include "concepts.hpp"
#include "logger.hpp"
#include "mutation_map_transcript.hpp"
#include "results/analysis.hpp"
#include "weighted_clusters.hpp"

#include <armadillo>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <mutex>
#include <oneapi/tbb/parallel_for.h>
#include <ranges>
#include <vector>

struct Window {
  unsigned short start_base;
  WeightedClusters weights;
  std::vector<unsigned> coverages;
};

class RingmapData;

namespace results {
struct Transcript;
} // namespace results

struct Args;

void handle_transcripts(
    std::vector<MutationMapTranscript const *> const &transcripts,
    std::vector<RingmapData *> const &ringmapsData,
    results::Analysis &analysisResult, Args const &args,
    std::optional<std::ofstream> &raw_n_clusters_stream,
    std::mutex &raw_n_clusters_stream_mutex);

void merge_windows_and_add_window_results(
    std::vector<Window> const &windows,
    std::span<std::vector<unsigned>> const &windows_reads_indices,
    RingmapData const &ringmap_data, results::Transcript &transcript_result,
    unsigned window_size, std::size_t replicate_index, Args const &args);

template <typename WRII, std::sentinel_for<WRII> S>
  requires std::input_or_output_iterator<WRII> and
           std::same_as<typename WRII::value_type, std::vector<std::uint32_t>>
inline auto make_windows_and_reads_indices_range(
    std::vector<Window>::const_iterator windows_iter,
    std::vector<Window>::const_iterator windows_iter_end,
    WRII windows_reads_indices_iter, S windows_reads_indices_iter_end,
    std::uint32_t min_windows_overlap) {
  auto const n_clusters = windows_iter->weights.getClustersSize();
  auto const window_size = windows_iter->weights.getElementsSize();

  auto zip =
      std::views::zip(std::ranges::subrange(windows_iter, windows_iter_end),
                      std::ranges::subrange(windows_reads_indices_iter,
                                            windows_reads_indices_iter_end));

  auto zip_maybe_end = std::ranges::end(zip);
  auto last_same_clusters_start_base = windows_iter->start_base;
  assert(std::ranges::begin(zip) != std::ranges::end(zip));
  for (auto zip_iter = std::ranges::next(std::ranges::begin(zip)),
            zip_end = std::ranges::end(zip);
       zip_iter != zip_end; ++zip_iter) {
    auto const &window = std::get<0>(*zip_iter);

    auto last_window_clusters = window.weights.getClustersSize();
    if (last_window_clusters == n_clusters) {
      last_same_clusters_start_base = window.start_base;
      zip_maybe_end = std::ranges::end(zip);
    } else {
      if (zip_maybe_end == std::ranges::end(zip)) {
        zip_maybe_end = std::next(zip_iter);
      }

      if (n_clusters == 0 or last_same_clusters_start_base + window_size <=
                                 window.start_base + min_windows_overlap) {
        break;
      }
    }
  }

  return std::ranges::subrange(std::ranges::begin(zip), zip_maybe_end);
}

template <typename WRII, typename R>
  requires std::input_or_output_iterator<WRII> and
           std::same_as<typename WRII::value_type,
                        std::vector<unsigned int>> and
           std::ranges::range<R> and
           (std::same_as<std::ranges::range_value_t<R>,
                         std::tuple<Window, std::vector<unsigned int>>> ||
            std::same_as<std::ranges::range_value_t<R>,
                         std::pair<Window, std::vector<unsigned int>>>)
void update_iters_and_region(
    typename std::vector<Window>::const_iterator &window_iter,
    WRII &window_reads_indices_iter, R windows_and_reads_indices_range,
    std::vector<std::optional<std::uint16_t>> &previous_overlapping_region_ends,
    unsigned min_windows_overlap) {
  auto const n_clusters = window_iter->weights.getClustersSize();
  auto const window_size = window_iter->weights.getElementsSize();

  auto const end_iter = std::ranges::end(windows_and_reads_indices_range);
  auto last_same_clusters_size_iter =
      std::ranges::begin(windows_and_reads_indices_range);
  assert(std::get<0>(*last_same_clusters_size_iter).weights.getClustersSize() ==
         n_clusters);

  std::ptrdiff_t region_size;
  decltype(last_same_clusters_size_iter) first_different_clusters_size_iter;
  if (last_same_clusters_size_iter == end_iter) {
    region_size = 0;
    first_different_clusters_size_iter = end_iter;
  } else {
    region_size = 1;
    first_different_clusters_size_iter =
        std::ranges::next(last_same_clusters_size_iter);
  }
  for (; first_different_clusters_size_iter != end_iter;
       ++last_same_clusters_size_iter, ++first_different_clusters_size_iter,
       ++region_size) {
    auto const &first_different_clusters_size =
        *first_different_clusters_size_iter;
    if (std::get<0>(first_different_clusters_size).weights.getClustersSize() !=
        n_clusters) {
      break;
    }
  }
  assert(std::get<0>(*last_same_clusters_size_iter).weights.getClustersSize() ==
         n_clusters);

  // There is no need to store the end of the region if there are
  // not other set of windows with different number of clusters.
  if (first_different_clusters_size_iter !=
      std::end(windows_and_reads_indices_range)) {
    auto const &last_with_same_n_clusters = *last_same_clusters_size_iter;

    if (n_clusters > 0) {
      auto const &last_window_with_same_n_clusters =
          std::get<0>(last_with_same_n_clusters);
      previous_overlapping_region_ends[n_clusters - 1] =
          last_window_with_same_n_clusters.start_base + window_size -
          min_windows_overlap;
    }
  }

  window_iter += region_size;
  window_reads_indices_iter += region_size;
};

struct PtbaOnReplicate {
  std::vector<unsigned> pre_collapsing_clusters;
  std::vector<Window> windows;
  unsigned int window_size;
};

std::vector<unsigned> get_best_pre_collapsing_clusters(
    std::span<std::optional<PtbaOnReplicate>> ptba_on_replicate_results);

void set_uninformative_clusters_to_surrounding(
    std::vector<unsigned> &windows_n_clusters,
    std::vector<std::optional<unsigned>> const
        &windows_max_clusters_constraints) noexcept(false);

void collapse_outlayer_clusters(std::vector<unsigned> &windows_n_clusters,
                                std::vector<std::optional<unsigned>> const
                                    &windows_max_clusters_constraints,
                                Args const &args) noexcept(false);

void assign_reads_to_clusters(
    results::Window &window,
    RingmapData::clusters_assignment_type &&clusters_assignment,
    RingmapData const &ringmap, RingmapData const &filteredRingmap);

void dump_assignments(results::Transcript const &transcript,
                      results::Window &window, RingmapData const &ringmap,
                      std::optional<std::size_t> replicate_index,
                      std::string_view assignments_dump_directory);

struct HandleTranscripts {
  std::vector<MutationMapTranscript const *> const &transcripts;
  std::vector<RingmapData *> const &ringmaps_data;
  results::Analysis &analysis_result;
  Args const &args;
  std::optional<std::ofstream> &raw_n_clusters_stream;
  std::mutex &raw_n_clusters_stream_mutex;
  bool use_stdout;

  void operator()(
      InvocableR<std::optional<PtbaOnReplicate>, std::size_t,
                 RingmapData const &, results::Transcript &> auto
          &&ptba_on_replicate,
      InvocableR<WeightedClusters, unsigned short, unsigned short, std::uint8_t,
                 std::vector<arma::mat> const &,
                 results::Transcript const &> auto &&get_weighted_clusters) {
    auto const &first_transcript = *transcripts[0];
    if (std::ranges::any_of(ringmaps_data, [](auto const &ringmap_data) {
          return ringmap_data->data().rows_size() == 0;
        })) {
      if (use_stdout) {
        std::cout << "\x1b[2K\r[+] Skipping transcript "
                  << first_transcript.getId() << " (no reads)" << std::endl;
      }
      return;
    }
    if (use_stdout) {
      std::cout << "\x1b[2K\r[+] Analyzing transcript "
                << first_transcript.getId() << std::flush;
    }

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
              replicate_index, ringmap_data, transcript_result);
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
        const auto windows_n_clusters_end =
            std::ranges::end(windows_n_clusters);

        auto windows_max_clusters_constraints_iter =
            std::ranges::begin(windows_max_clusters_constraints);
        const auto windows_max_clusters_constraints_end =
            std::ranges::end(windows_max_clusters_constraints);

        for (; windows_n_clusters_iter < windows_n_clusters_end and
               windows_max_clusters_constraints_iter <
                   windows_max_clusters_constraints_end;
             ++windows_n_clusters_iter,
             ++windows_max_clusters_constraints_iter) {
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
        assert(std::ranges::none_of(windows_n_clusters, zero_clusters) or
               std::ranges::all_of(windows_n_clusters, zero_clusters));
      }

      if (args.max_collapsing_windows() > 0)
        collapse_outlayer_clusters(windows_n_clusters,
                                   windows_max_clusters_constraints, args);

      if (args.set_all_uninformative_to_one()) {
        if (std::ranges::all_of(windows_n_clusters, [](auto n_clusters) {
              return n_clusters == 0;
            })) {
          std::ranges::fill(windows_n_clusters, 1u);
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
        throw std::runtime_error(std::format(
            "incoherent windows between replicates for transcript {}",
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
                    static_cast<std::ptrdiff_t>(replicate_index *
                                                windows_size)),
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
                auto &&[ringmap_data, ptba_on_replicate_result,
                        replicate_index] = std::move(tuple);
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

                replicates_windows_reads_indices(replicate_index,
                                                 window_index) =
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
                std::ranges::any_of(
                    replicates_filtered_data, [](auto const &filtered_data) {
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

              assert(first_result_window.start_base <=
                     std::numeric_limits<
                         decltype(first_result_window.start_base)>::max() -
                         window_size);
              auto graphCutResults = get_weighted_clusters(
                  first_result_window.start_base,
                  first_result_window.start_base +
                      static_cast<decltype(first_result_window.start_base)>(
                          window_size),
                  n_clusters, replicates_covariance, transcript_result);

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
            auto &&[replicate_index, ptba_on_replicate_result,
                    ringmap_data_ptr] = std::move(tuple);
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
                                return (
                                    window.has_value() !=
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
              auto splitted_ringmaps = RingmapData(*ringmap_data)
                                           .split_into_windows(*result_windows);
              assert(std::size(splitted_ringmaps) ==
                     std::size(*result_windows));
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
            std::cerr
                << "\nInconsistency between weighted clusters sizes across "
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
            std::ranges::fill(window.assignments, std::int8_t(-1));
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

                auto &&fractions_result =
                    filtered_ringmap.fractionReadsByWeights(
                        window.weighted_clusters, window_size,
                        args.skip_ambiguous_assignments());
                std::tie(window.fractions, window.patterns, std::ignore) =
                    std::move(fractions_result);
                assert(window.fractions.size() > 1 or
                       window.fractions.empty() or window.fractions[0] >= 0.01);

                bool const redundand_patterns = [&] {
                  auto patterns_iter = std::cbegin(*window.patterns);
                  auto const patterns_end = std::cend(*window.patterns);

                  for (; patterns_iter < patterns_end; ++patterns_iter) {
                    auto &&cur_pattern = *patterns_iter;
                    auto const begin_cur_pattern = std::cbegin(cur_pattern);
                    auto const end_cur_pattern = std::cend(cur_pattern);

                    if (std::any_of(std::next(patterns_iter), patterns_end,
                                    [&](auto &&next_pattern) {
                                      return std::equal(
                                          begin_cur_pattern, end_cur_pattern,
                                          std::cbegin(next_pattern),
                                          std::cend(next_pattern));
                                    })) {
                      return true;
                    }
                  }

                  return false;
                }();

                if (redundand_patterns or
                    std::ranges::any_of(
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
                          window.start_base + window_size <=
                              result_window_end) {

                        window_constraint = new_clusters_constraint;
                      }
                    }
                  }
                }

                if (not stop)
                  return;

                *window.patterns =
                    filtered_ringmap.remapPatterns(*window.patterns);

                assign_reads_to_clusters(
                    window, std::move(std::get<2>(fractions_result)), ringmap,
                    filtered_ringmap);

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
                                std::cend(*window.patterns),
                                [](auto &&pattern) {
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
};

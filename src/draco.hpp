#pragma once

#include "weighted_clusters.hpp"

#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <mutex>
#include <ranges>
#include <vector>

struct Window {
  unsigned short start_base;
  WeightedClusters weights;
  std::vector<unsigned> coverages;
};

class MutationMapTranscript;
class RingmapData;

namespace results {
struct Analysis;
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

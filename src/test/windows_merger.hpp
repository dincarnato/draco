#pragma once

#include "../windows_merger.hpp"

namespace windows_merger {
namespace test {

struct WindowsMerger {
  using NaiveWindowsMerger = windows_merger::WindowsMerger;
  using bases_size_type = typename NaiveWindowsMerger::bases_size_type;
  using input_weights_type = typename NaiveWindowsMerger::input_weights_type;
  using input_coverages_type =
      typename NaiveWindowsMerger::input_coverages_type;
  using windows_size_type = typename NaiveWindowsMerger::windows_size_type;
  using clusters_size_type = typename NaiveWindowsMerger::clusters_size_type;
  using distance_type = typename NaiveWindowsMerger::distance_type;

  inline static void
  transform_window_to_native(NaiveWindowsMerger& windows_merger,
                             bases_size_type start_offset,
                             const input_weights_type& weights,
                             const input_coverages_type& coverages) {
    windows_merger.transform_window_to_native(start_offset, weights, coverages);
  }

  inline static void
  process_queue(NaiveWindowsMerger& windows_merger) {
    windows_merger.process_queue();
  }

  inline static void
  prepare_cache(NaiveWindowsMerger& windows_merger) {
    windows_merger.prepare_cache();
  }

  inline static void
  prepare_indices(NaiveWindowsMerger& windows_merger) {
    windows_merger.prepare_indices();
  }

  inline static void
  prepare_high_coverages(NaiveWindowsMerger& windows_merger) {
    windows_merger.prepare_high_coverages();
  }

  inline static void
  prepare_distances(NaiveWindowsMerger& windows_merger) {
    windows_merger.prepare_distances();
  }

  inline static void
  update_distance_between(NaiveWindowsMerger& windows_merger,
                          windows_size_type cache_index_a,
                          windows_size_type cache_index_b) noexcept {
    windows_merger.update_distance_between(cache_index_a, cache_index_b);
  }

  inline static std::pair<double, std::vector<clusters_size_type> const&>
  get_best_distance_and_permutation_between(
      const NaiveWindowsMerger& windows_merger, windows_size_type cache_index_a,
      windows_size_type cache_index_b) noexcept {
    return windows_merger.get_best_distance_and_permutation_between(
        cache_index_a, cache_index_b);
  }

  inline static distance_type
  distance_between_bases(
      const NaiveWindowsMerger& windows_merger,
      typename WindowsMergerWindows::const_base_accessor const& base_a,
      typename WindowsMergerWindows::const_base_accessor const& base_b,
      std::vector<clusters_size_type> const&
          cluster_matching_indices) noexcept {
    return windows_merger.distance_between_bases(base_a, base_b,
                                                 cluster_matching_indices);
  }

  inline static std::array<windows_size_type, 2>
  find_best_cached_pair(const NaiveWindowsMerger& windows_merger) noexcept {
    return windows_merger.find_best_cached_pair();
  }

  static inline void
  merge_cached_windows_into_first(
      NaiveWindowsMerger& windows_merger, windows_size_type first_index_a,
      windows_size_type second_index_a) noexcept(false) {
    windows_merger.merge_cached_windows_into_first(first_index_a,
                                                   second_index_a);
  }

  static inline void
  resize_and_collapse_cache(NaiveWindowsMerger& windows_merger,
                            bases_size_type new_capacity) noexcept(false) {
    windows_merger.resize_and_collapse_cache(new_capacity);
  }

  static inline void
  update_indices_for_post_collapsing(const NaiveWindowsMerger& windows_merger,
                                     windows_size_type& index_a,
                                     windows_size_type& index_b) noexcept {
    windows_merger.update_indices_for_post_collapsing(index_a, index_b);
  }

  static inline auto&
  dequeuing(NaiveWindowsMerger& windows_merger) {
    return windows_merger.dequeueing;
  }

  static inline const auto&
  dequeuing(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.dequeueing;
  }

  static inline auto&
  queue(NaiveWindowsMerger& windows_merger) {
    return windows_merger.queue;
  }

  static inline const auto&
  queue(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.queue;
  }

  static inline auto&
  windows(NaiveWindowsMerger& windows_merger) {
    return windows_merger.windows;
  }

  static inline const auto&
  windows(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.windows;
  }

  static inline auto&
  cache(NaiveWindowsMerger& windows_merger) {
    return windows_merger.cache;
  }

  static inline const auto&
  cache(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.cache;
  }

  static inline auto&
  cache_indices(NaiveWindowsMerger& windows_merger) {
    return windows_merger.cache_indices;
  }

  static inline const auto&
  cache_indices(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.cache_indices;
  }

  static inline auto&
  windows_distances(NaiveWindowsMerger& windows_merger) {
    return windows_merger.windows_distances;
  }

  static inline const auto&
  windows_distances(const NaiveWindowsMerger& windows_merger) {
    return windows_merger.windows_distances;
  }
};

} // namespace test
} // namespace windows_merger

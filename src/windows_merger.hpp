#pragma once

#include "windows_merger_cache_indices.hpp"
#include "windows_merger_exceptions.hpp"
#include "windows_merger_traits.hpp"
#include "windows_merger_windows.hpp"

#include "triangular_matrix_strict.hpp"

#include "weighted_clusters.hpp"

#include "parallel/blocking_queue.hpp"

#include <atomic>
#include <cstdint>
#include <tuple>
#include <vector>

namespace windows_merger {

namespace test {

struct WindowsMerger;

} // namespace test

struct WindowsMerger {
  using traits_type = WindowsMergerTraits;
  using coverage_type = typename traits_type::coverage_type;
  using bases_size_type = typename traits_type::bases_size_type;
  using clusters_size_type = typename traits_type::clusters_size_type;
  using windows_size_type = typename traits_type::windows_size_type;
  using input_weights_type = WeightedClusters;
  using input_coverages_type = std::vector<unsigned>;
  using distance_type = double;
  using distances_type = TriangularMatrixStrict<distance_type>;

  friend test::WindowsMerger;

private:
  using queue_type = parallel::blocking_queue<
      std::tuple<bases_size_type, input_weights_type, input_coverages_type>>;

public:
  explicit WindowsMerger(clusters_size_type n_clusters) noexcept;

  WindowsMerger(const WindowsMerger &) = delete;
  WindowsMerger &operator=(const WindowsMerger &) = delete;
  WindowsMerger(WindowsMerger &&) noexcept(
      std::is_nothrow_move_constructible_v<queue_type>
          and std::is_nothrow_move_constructible_v<WindowsMergerWindows> and
              std::is_nothrow_move_constructible_v<WindowsMergerCacheIndices>
                  and std::is_nothrow_move_constructible_v<distances_type>);
  WindowsMerger &operator=(WindowsMerger &&) noexcept(
      std::is_nothrow_move_assignable_v<queue_type>
          and std::is_nothrow_move_assignable_v<WindowsMergerWindows>
              and std::is_nothrow_move_assignable_v<WindowsMergerCacheIndices>
                  and std::is_nothrow_move_assignable_v<distances_type>);

  template <typename Weights, typename Coverages>
  void add_window(bases_size_type start_offset, Weights &&weights,
                  Coverages &&coverages);

  void wait_queue();
  const WindowsMergerWindows &get_windows() const noexcept;
  WindowsMergerWindow merge() noexcept(false);

private:
  void transform_window_to_native(bases_size_type start_offset,
                                  const input_weights_type &weights,
                                  const input_coverages_type &coverages);
  void process_queue();
  void prepare_cache();
  void prepare_indices();
  void prepare_distances();
  void prepare_high_coverages();
  void update_distance_between(windows_size_type cache_index_a,
                               windows_size_type cache_index_b) noexcept;
  std::pair<double, std::vector<clusters_size_type> const &>
  get_best_distance_and_permutation_between(
      windows_size_type cache_index_a,
      windows_size_type cache_index_b) const noexcept;
  void
  update_high_coverage_cache(windows_size_type cache_window_index) noexcept;

  distance_type distance_between_bases(
      typename WindowsMergerWindows::const_base_accessor const &base_a,
      typename WindowsMergerWindows::const_base_accessor const &base_b,
      std::vector<clusters_size_type> const &cluster_matching_indices)
      const noexcept;
  std::array<windows_size_type, 2> find_best_cached_pair() const noexcept;
  void merge_cached_windows_into_first(
      windows_size_type first_index_a,
      windows_size_type second_index_a) noexcept(false);
  void resize_and_collapse_cache(bases_size_type new_capacity) noexcept(false);
  void
  update_indices_for_post_collapsing(windows_size_type &index_a,
                                     windows_size_type &index_b) const noexcept;

  std::atomic_bool dequeueing = false;
  queue_type queue;
  WindowsMergerWindows windows;
  WindowsMergerWindows cache;
  WindowsMergerCacheIndices cache_indices;
  distances_type windows_distances;
  std::vector<coverage_type> cache_high_coverages;

  double non_overlap_penalty = 0.5;
  static constexpr bases_size_type initial_cache_bases_capacity = 4;
};

} // namespace windows_merger

#include "windows_merger_impl.hpp"

#include "test/windows_merger.hpp"

#include "windows_merger.hpp"

#include <numeric>
#include <type_traits>

namespace ranges = std::ranges;

namespace windows_merger {

WindowsMerger::WindowsMerger(clusters_size_type n_clusters) noexcept
    : windows(n_clusters), cache(n_clusters) {}

WindowsMerger::WindowsMerger(WindowsMerger &&other) noexcept(
    std::is_nothrow_move_constructible_v<queue_type> and
    std::is_nothrow_move_constructible_v<WindowsMergerWindows> and
    std::is_nothrow_move_constructible_v<WindowsMergerCacheIndices> and
    std::is_nothrow_move_constructible_v<distances_type>)
    : dequeueing(other.dequeueing.load(std::memory_order_acquire)),
      queue(std::move(other.queue)), windows(std::move(other.windows)),
      cache(std::move(other.cache)),
      cache_indices(std::move(other.cache_indices)),
      windows_distances(std::move(other.windows_distances)),
      non_overlap_penalty(other.non_overlap_penalty) {}

WindowsMerger &WindowsMerger::operator=(WindowsMerger &&other) noexcept(
    std::is_nothrow_move_assignable_v<queue_type> and
    std::is_nothrow_move_assignable_v<WindowsMergerWindows> and
    std::is_nothrow_move_assignable_v<WindowsMergerCacheIndices> and
    std::is_nothrow_move_assignable_v<distances_type>) {
  dequeueing.store(other.dequeueing.load(std::memory_order_acquire),
                   std::memory_order_release);
  queue = std::move(other.queue);
  windows = std::move(other.windows);
  cache = std::move(other.cache);
  cache_indices = std::move(other.cache_indices);
  windows_distances = std::move(other.windows_distances);
  non_overlap_penalty = other.non_overlap_penalty;

  return *this;
}

void WindowsMerger::process_queue() {
  for (;;) {
    auto maybe_window = queue.try_pop();
    if (maybe_window) {
      std::apply(
          [this](auto &&...args) {
            transform_window_to_native(std::forward<decltype(args)>(args)...);
          },
          std::move(*maybe_window));
    } else {
      dequeueing.store(false, std::memory_order::release);
      break;
    }
  }
}

void WindowsMerger::wait_queue() {
  queue.finished();
  if (bool expected = false; dequeueing.compare_exchange_strong(
          expected, true, std::memory_order::acq_rel)) {
    process_queue();
  }

  queue.finish();
  ranges::sort(windows, ranges::less{},
               [](auto &&window) { return window.begin_index(); });
}

void WindowsMerger::transform_window_to_native(
    bases_size_type start_offset, const input_weights_type &weights,
    const input_coverages_type &coverages) {

  const auto bases_size =
      static_cast<typename WindowsMergerTraits::bases_size_type>(
          weights.getElementsSize());
  windows.reshape(bases_size, WindowsMergerWindows::reserver_type(0));

  auto new_window = windows.emplace_back();
  new_window.resize(bases_size);
  new_window.set_begin_index(start_offset);

  {
    auto weights_iter = ranges::begin(weights);
    const auto weights_end = ranges::end(weights);
    auto window_iter = ranges::begin(new_window);

    for (; weights_iter < weights_end; ++weights_iter, ++window_iter) {
      auto &&weights_span = *weights_iter;
      auto &&base_accessor = *window_iter;
      auto &&base_weights = base_accessor.weights();
      ranges::transform(weights_span, ranges::begin(base_weights),
                        [](auto &&weight) { return TinyFraction(weight); });
    }
  }

  auto &&window_coverages = new_window.coverages();
  ranges::copy(coverages, ranges::begin(window_coverages));
}

const WindowsMergerWindows &WindowsMerger::get_windows() const noexcept {
  return windows;
}

WindowsMergerWindow WindowsMerger::merge() noexcept(false) {
  wait_queue();
  prepare_cache();
  prepare_indices();
  prepare_high_coverages();
  prepare_distances();

  for (;;) {
    const auto best_pair_indices = find_best_cached_pair();
    if (best_pair_indices[1] == std::numeric_limits<windows_size_type>::max()) {
      auto last_window_index = best_pair_indices[0];
      assert(ranges::count_if(
                 cache, [](auto &&window) { return not window.empty(); }) == 1);
      assert(ranges::count_if(cache_indices, [](auto &&indices) {
               return not indices.empty();
             }) == 1);

      return cache[last_window_index];
    }

    merge_cached_windows_into_first(best_pair_indices[0], best_pair_indices[1]);
  }
}

void WindowsMerger::prepare_cache() {
  cache.reshape(windows.bases_capacity(),
                WindowsMergerWindows::resizer_type(windows.windows_size()));

  ranges::copy(windows, ranges::begin(cache));
}

void WindowsMerger::prepare_indices() {
  const auto windows_size = cache.windows_size();
  cache_indices.reshape(initial_cache_bases_capacity,
                        WindowsMergerCacheIndices::resizer_type(windows_size));
  for (windows_size_type window_index = 0; window_index < windows_size;
       ++window_index) {
    auto &&accessor = cache_indices[window_index];
    accessor.emplace_back(window_index);
  }
}

void WindowsMerger::prepare_distances() {
  const auto windows_size = cache.windows_size();
  windows_distances = distances_type(
      windows_size, std::numeric_limits<distance_type>::infinity());

  for (windows_size_type window_a_index = 0; window_a_index < windows_size;
       ++window_a_index) {
    const auto window_a = std::as_const(cache)[window_a_index];
    const auto window_a_end = window_a.end_index();

    for (auto window_b_index = window_a_index + 1;
         window_b_index < windows_size; ++window_b_index) {
      const auto window_b = std::as_const(cache)[window_b_index];
      const auto window_b_begin = window_b.begin_index();

      if (window_b_begin >= window_a_end)
        break;

      update_distance_between(window_a_index, window_b_index);
    }
  }
}

void WindowsMerger::prepare_high_coverages() {
  const auto cache_size = cache.windows_size();
  cache_high_coverages.resize(cache_size);
  for (windows_size_type cache_window_index = 0;
       cache_window_index < cache_size; ++cache_window_index) {
    update_high_coverage_cache(cache_window_index);
  }
}

void WindowsMerger::update_distance_between(
    windows_size_type cache_index_a, windows_size_type cache_index_b) noexcept {
  windows_distances[cache_index_a][cache_index_b] = std::get<0>(
      get_best_distance_and_permutation_between(cache_index_a, cache_index_b));
}

auto WindowsMerger::get_best_distance_and_permutation_between(
    windows_size_type cache_index_a,
    windows_size_type cache_index_b) const noexcept
    -> std::pair<distance_type, std::vector<clusters_size_type> const &> {
  thread_local std::vector<clusters_size_type> permutation_indices;
  thread_local std::vector<clusters_size_type> best_permutation;

  auto best_result =
      std::make_pair(std::numeric_limits<distance_type>::infinity(),
                     std::cref(best_permutation));

  assert(not cache_indices[cache_index_a].empty());
  assert(not cache_indices[cache_index_b].empty());

  auto &&cache_window_a = std::as_const(cache)[cache_index_a];
  auto &&cache_window_b = std::as_const(cache)[cache_index_b];
  assert(not cache_indices[cache_index_a].empty());
  assert(not cache_indices[cache_index_b].empty());

  const auto clusters_size = cache.clusters_size();
  permutation_indices.resize(clusters_size);
  best_permutation.resize(clusters_size);
  std::iota(std::ranges::begin(permutation_indices),
            std::ranges::end(permutation_indices), clusters_size_type(0));
  std::ranges::copy(permutation_indices, ranges::begin(best_permutation));

  const auto cache_window_a_begin = cache_window_a.begin_index();
  const auto cache_window_a_end = cache_window_a.end_index();
  const auto cache_window_b_begin = cache_window_b.begin_index();
  const auto cache_window_b_end = cache_window_b.end_index();

  const auto first_base = std::max(cache_window_a_begin, cache_window_b_begin);
  const auto last_base = std::min(cache_window_a_end, cache_window_b_end);

  const auto window_a_high_coverage = cache_high_coverages[cache_index_a];
  const auto window_b_high_coverage = cache_high_coverages[cache_index_b];

  auto &best_distance = std::get<0>(best_result);

  const auto permutation_begin = ranges::begin(permutation_indices);
  const auto permutation_end = ranges::end(permutation_indices);

  do {
    auto current_distance = std::numeric_limits<distance_type>::infinity();

    for (bases_size_type base_index = first_base,
                         base_a_index = first_base - cache_window_a_begin,
                         base_b_index = first_base - cache_window_b_begin;
         base_index < last_base; ++base_index, ++base_a_index, ++base_b_index) {

      const auto cache_base_a_coverage =
          cache_window_a[base_a_index].coverage();
      const auto cache_base_b_coverage =
          cache_window_b[base_b_index].coverage();

      // TODO: bring this out of loop, result is independent from permutation
      const coverage_type cum_coverage =
          cache_base_a_coverage + cache_base_b_coverage;

      if (cum_coverage == 0)
        continue;

      const double base_a_normalizer =
          static_cast<double>(cache_base_a_coverage) / cum_coverage;
      const double base_b_normalizer =
          static_cast<double>(cache_base_b_coverage) / cum_coverage;

      auto &&cache_base_a = std::as_const(cache_window_a)[base_a_index];
      auto &&cache_base_b = std::as_const(cache_window_b)[base_b_index];

      for (clusters_size_type first_cluster_index = 0;
           first_cluster_index < clusters_size; ++first_cluster_index) {
        clusters_size_type second_cluster_index =
            permutation_indices[first_cluster_index];

        const auto cache_base_a_weight =
            static_cast<double>(cache_base_a.weight(first_cluster_index));
        const auto cache_base_b_weight =
            static_cast<double>(cache_base_b.weight(second_cluster_index));

        double mean = cache_base_a_weight * base_a_normalizer +
                      cache_base_b_weight * base_b_normalizer;

        bool intersecting = false;
        auto get_window_distance = [&](windows_size_type cache_window_index,
                                       clusters_size_type cluster_index,
                                       double weight) -> double {
          auto &&window_indices =
              std::as_const(cache_indices)[cache_window_index];
          auto &&retval = std::accumulate(
              std::ranges::begin(window_indices),
              std::ranges::end(window_indices), 0.,
              [&](double acc, windows_size_type window_index) {
                auto &&window = std::as_const(windows)[window_index];
                if (base_index < window.begin_index() or
                    base_index >= window.end_index())
                  return acc;
                else {
                  intersecting = true;
                  double value =
                      acc +
                      std::abs(
                          static_cast<double>(
                              window[base_index - window.begin_index()].weight(
                                  cluster_index)) -
                          mean) *
                          weight;
                  return value;
                }
              });
          return retval;
        };

        const auto base_a_coverage_normalizer =
            static_cast<double>(cache_base_a_coverage) /
            (window_a_high_coverage > 0 ? window_a_high_coverage : 1);
        const auto base_b_coverage_normalizer =
            static_cast<double>(cache_base_b_coverage) /
            (window_b_high_coverage > 0 ? window_b_high_coverage : 1);

        double fragment_distance =
            get_window_distance(cache_index_a, first_cluster_index,
                                base_a_coverage_normalizer) +
            get_window_distance(cache_index_b, second_cluster_index,
                                base_b_coverage_normalizer);
        if (intersecting) {
          if (current_distance ==
              std::numeric_limits<distance_type>::infinity())
            current_distance = 0.;

          current_distance += fragment_distance;
        }
      }
    }

    if (current_distance < best_distance) {
      best_distance = current_distance;
      ranges::copy(permutation_begin, permutation_end,
                   ranges::begin(best_permutation));
    }
  } while (std::next_permutation(permutation_begin, permutation_end));

  const double penalty =
      static_cast<double>(
          (first_base - std::min(cache_window_a_begin, cache_window_b_begin)) +
          (std::max(cache_window_a_end, cache_window_b_end) - last_base)) *
      non_overlap_penalty / 2.;
  best_distance += penalty;

  return best_result;
}

void WindowsMerger::update_high_coverage_cache(
    windows_size_type cache_window_index) noexcept {
  thread_local std::vector<coverage_type> temp_coverages;

  assert(cache_window_index < cache_high_coverages.size());

  auto &&cache_window = std::as_const(cache)[cache_window_index];
  auto &&coverages = cache_window.coverages();

  const auto coverages_size = coverages.size();
  if (coverages_size == 0) {
    cache_high_coverages[cache_window_index] = 0;
  } else {
    temp_coverages.resize(coverages_size);
    ranges::copy(coverages, ranges::begin(temp_coverages));

    auto percentile_90_index =
        static_cast<std::ptrdiff_t>(coverages_size * 0.9);
    if (percentile_90_index == coverages_size)
      --percentile_90_index;

    auto percentile_90 =
        ranges::next(ranges::begin(temp_coverages), percentile_90_index);

    assert(percentile_90 >= ranges::begin(temp_coverages));
    assert(percentile_90 < ranges::end(temp_coverages));

    ranges::nth_element(ranges::begin(temp_coverages), percentile_90,
                        ranges::end(temp_coverages));
    cache_high_coverages[cache_window_index] = *percentile_90;
  }
}

auto WindowsMerger::distance_between_bases(
    typename WindowsMergerWindows::const_base_accessor const &base_a,
    typename WindowsMergerWindows::const_base_accessor const &base_b,
    std::vector<clusters_size_type> const &cluster_matching_indices)
    const noexcept -> distance_type {

  distance_type distance = 0;
  const auto total_coverage =
      static_cast<double>(base_a.coverage() + base_b.coverage());
  const auto coverage_normalizer = [&] {
    if (total_coverage != 0.)
      return 1. / total_coverage;
    else
      return 1.;
  }();

  const auto base_a_normalizer =
      static_cast<double>(base_a.coverage()) * coverage_normalizer;
  const auto base_b_normalizer =
      static_cast<double>(base_b.coverage()) * coverage_normalizer;

  const auto clusters_size = cache.clusters_size();
  for (clusters_size_type first_cluster_index = 0;
       first_cluster_index < clusters_size; ++first_cluster_index) {
    clusters_size_type second_cluster_index =
        cluster_matching_indices[first_cluster_index];

    distance += std::abs(
        static_cast<distance_type>(base_a.weight(first_cluster_index)) *
            base_a_normalizer -
        static_cast<distance_type>(base_b.weight(second_cluster_index)) *
            base_b_normalizer);
  }

  assert(distance != 0);
  return distance;
}

auto WindowsMerger::find_best_cached_pair() const noexcept
    -> std::array<windows_size_type, 2> {
  auto best_pair = std::array<windows_size_type, 2>{
      std::numeric_limits<windows_size_type>::max(),
      std::numeric_limits<windows_size_type>::max()};
  auto best_distance = std::numeric_limits<distance_type>::infinity();
  std::size_t usable_windows = 0;

  const auto cache_size = cache.windows_size();
  for (windows_size_type first_window_index = 0;
       first_window_index < cache_size; ++first_window_index) {
    if (cache_indices[first_window_index].empty())
      continue;

    ++usable_windows;
    auto &&distances_line =
        std::as_const(windows_distances)[first_window_index];
    for (windows_size_type second_window_index = first_window_index + 1;
         second_window_index < cache_size; ++second_window_index) {
      if (cache_indices[second_window_index].empty())
        continue;

      auto distance = distances_line[second_window_index];
      assert(distance > 0);
      if (distance < best_distance) {
        best_pair[0] = first_window_index;
        best_pair[1] = second_window_index;

        best_distance = distance;
      }
    }

    if (best_distance == std::numeric_limits<distance_type>::infinity()) {
      best_pair[0] = first_window_index;
    }
  }

  if (usable_windows > 1 and
      best_distance == std::numeric_limits<distance_type>::infinity()) {

    struct BestPair {
      std::size_t first_index = std::numeric_limits<std::size_t>::max();
      std::size_t second_index = std::numeric_limits<std::size_t>::max();
      std::size_t new_size = std::numeric_limits<std::size_t>::max();
    };

    BestPair sized_best_pair;
    auto const n_windows = cache.windows_size();
    for (std::size_t first_window_index = 0; first_window_index < n_windows - 1;
         ++first_window_index) {
      auto &&first_window = cache[first_window_index];

      if (first_window.empty()) {
        continue;
      }

      for (std::size_t second_window_index = first_window_index + 1;
           second_window_index < n_windows; ++second_window_index) {
        auto &&second_window = cache[second_window_index];

        if (second_window.empty()) {
          continue;
        } else {
          if (first_window.end_index() < second_window.end_index()) {
            auto const new_size = static_cast<std::size_t>(
                second_window.end_index() - first_window.begin_index());
            if (new_size < sized_best_pair.new_size) {
              sized_best_pair.first_index = first_window_index;
              sized_best_pair.second_index = second_window_index;
              sized_best_pair.new_size = new_size;
            }

            break;
          } else {
            auto const new_size = first_window.size();
            if (new_size < sized_best_pair.new_size) {
              sized_best_pair.first_index = first_window_index;
              sized_best_pair.second_index = second_window_index;
              sized_best_pair.new_size = new_size;
            }
          }
        }
      }
    }

    assert(sized_best_pair.new_size != std::numeric_limits<std::size_t>::max());
    best_pair[0] = sized_best_pair.first_index;
    best_pair[1] = sized_best_pair.second_index;
  }

  return best_pair;
}

void WindowsMerger::merge_cached_windows_into_first(
    windows_size_type first_index,
    windows_size_type second_index) noexcept(false) {

  const auto new_window_size = [&] {
    auto &&first_cached_window = cache[first_index];
    auto &&second_cached_window = cache[second_index];

    assert(first_index < second_index);
    assert(first_cached_window.begin_index() <
           second_cached_window.begin_index());
    const auto new_window_size = [&] {
      const auto new_window_size = static_cast<bases_size_type>(
          second_cached_window.end_index() - first_cached_window.begin_index());
      const auto current_window_size = first_cached_window.size();
      return std::max(new_window_size, current_window_size);
    }();
    if (new_window_size > cache.bases_capacity()) {
      update_indices_for_post_collapsing(first_index, second_index);
      assert(new_window_size < std::numeric_limits<bases_size_type>::max() / 2);
      resize_and_collapse_cache(
          static_cast<bases_size_type>(new_window_size * 2));
    }

    return new_window_size;
  }();

  auto &&first_cached_window = cache[first_index];
  auto &&second_cached_window = cache[second_index];

  const auto first_cached_window_old_end_index =
      first_cached_window.end_index();
  auto permutation_indices =
      get_best_distance_and_permutation_between(first_index, second_index)
          .second;

  assert(new_window_size >= first_cached_window.size());
  first_cached_window.resize(new_window_size);

  auto first_win_base_index = static_cast<bases_size_type>(
      second_cached_window.begin_index() - first_cached_window.begin_index());
  bases_size_type second_win_base_index = 0;
  const auto intersection_size =
      static_cast<bases_size_type>(std::min(first_cached_window_old_end_index,
                                            second_cached_window.end_index()) -
                                   second_cached_window.begin_index());
  const auto clusters_size = cache.clusters_size();

  for (; second_win_base_index < intersection_size;
       ++first_win_base_index, ++second_win_base_index) {
    auto &&first_win_base = first_cached_window[first_win_base_index];
    auto &&second_win_base = second_cached_window[second_win_base_index];
    auto &&first_win_weights = first_win_base.weights();
    auto &&second_win_weights = second_win_base.weights();
    const auto first_win_base_coverage = first_win_base.coverage();
    const auto second_win_base_coverage = second_win_base.coverage();
    const auto coverages_sum =
        first_win_base_coverage + second_win_base_coverage;
    const auto first_win_base_coverage_weight = [&] {
      if (coverages_sum != 0)
        return static_cast<double>(first_win_base_coverage) / coverages_sum;
      else
        return 0.;
    }();
    const auto second_win_base_coverage_weight = [&] {
      if (coverages_sum != 0)
        return static_cast<double>(second_win_base_coverage) / coverages_sum;
      else
        return 0.;
    }();

    for (clusters_size_type first_cluster_index = 0;
         first_cluster_index < clusters_size; ++first_cluster_index) {

      const auto second_cluster_index =
          permutation_indices[first_cluster_index];
      first_win_weights[first_cluster_index] = static_cast<TinyFraction>(
          static_cast<double>(first_win_weights[first_cluster_index]) *
              first_win_base_coverage_weight +
          static_cast<double>(second_win_weights[second_cluster_index]) *
              second_win_base_coverage_weight);
    }

    first_win_base.coverage() += second_win_base_coverage;
  }

  for (auto second_cached_window_size = second_cached_window.size();
       second_win_base_index < second_cached_window_size;
       ++second_win_base_index, ++first_win_base_index) {
    auto &&first_win_base = first_cached_window[first_win_base_index];
    auto &&second_win_base = second_cached_window[second_win_base_index];
    auto &&first_win_weights = first_win_base.weights();
    auto &&second_win_weights = second_win_base.weights();

    for (clusters_size_type first_cluster_index = 0;
         first_cluster_index < clusters_size; ++first_cluster_index) {

      const auto second_cluster_index =
          permutation_indices[first_cluster_index];
      first_win_weights[first_cluster_index] =
          std::move(second_win_weights[second_cluster_index]);
    }
    first_win_base.coverage() = second_win_base.coverage();
  }

  second_cached_window.clear();

  for (windows_size_type window_index : cache_indices[second_index]) {
    auto &&window = windows[window_index];
    window.reorder_clusters(permutation_indices);
  }

  {
    const auto first_cache_indices_old_size = cache_indices[first_index].size();
    const auto new_line_size =
        first_cache_indices_old_size + cache_indices[second_index].size();
    cache_indices.reshape(
        0, typename WindowsMergerCacheIndices::reserver_type(new_line_size));

    auto &&first_cache_indices = cache_indices[first_index];
    auto &&second_cache_indices = cache_indices[second_index];
    first_cache_indices.resize(new_line_size);
    using diff_t =
        typename WindowsMergerCacheIndices::iterator::difference_type;
    ranges::move(
        second_cache_indices,
        ranges::next(ranges::begin(first_cache_indices),
                     static_cast<diff_t>(first_cache_indices_old_size)));
    second_cache_indices.clear();
  }

  update_high_coverage_cache(first_index);
  for (windows_size_type other_window_index = 0;
       other_window_index < first_index; ++other_window_index) {
    if (not cache_indices[other_window_index].empty())
      update_distance_between(other_window_index, first_index);
    else
      windows_distances[other_window_index][first_index] =
          std::numeric_limits<distance_type>::infinity();
  }
  for (windows_size_type other_window_index = first_index + 1,
                         windows_size = cache.windows_size();
       other_window_index < windows_size; ++other_window_index) {
    if (not cache_indices[other_window_index].empty())
      update_distance_between(first_index, other_window_index);
    else
      windows_distances[first_index][other_window_index] =
          std::numeric_limits<distance_type>::infinity();
  }
}

void WindowsMerger::update_indices_for_post_collapsing(
    windows_size_type &index_a, windows_size_type &index_b) const noexcept {
  bool updated_one_index = false;
  windows_size_type available_indices = 0;

  auto cache_iter = ranges::begin(cache);
  auto cache_indices_iter = ranges::begin(cache_indices);
  const auto windows_size = cache.windows_size();

  for (windows_size_type window_index = 0; window_index < windows_size;
       ++window_index, ++cache_iter, ++cache_indices_iter++) {

    if (window_index == index_a) {
      index_a = available_indices;
      if (updated_one_index)
        return;
      else
        updated_one_index = true;
    }

    if (window_index == index_b) {
      index_b = available_indices;
      if (updated_one_index)
        return;
      else
        updated_one_index = true;
    }

    if (not(*cache_indices_iter).empty())
      ++available_indices;
  }

  std::cerr << "invalid indices\n";
  std::terminate();
}

void WindowsMerger::resize_and_collapse_cache(
    bases_size_type new_capacity) noexcept(false) {
  assert(new_capacity > cache.bases_capacity());

  const auto new_windows_size = static_cast<windows_size_type>(ranges::count_if(
      cache_indices, [](auto &&indices) { return not indices.empty(); }));
  assert(new_windows_size ==
         static_cast<windows_size_type>(ranges::count_if(
             cache, [](auto &&window) { return not window.empty(); })));

  const auto n_clusters = cache.clusters_size();
  WindowsMergerWindows new_cache(n_clusters, new_capacity, new_windows_size);
  std::vector<coverage_type> new_cache_high_coverages(new_windows_size);
  WindowsMergerCacheIndices new_cache_indices(new_capacity, new_windows_size);
  std::vector<windows_size_type> valid_indices(new_windows_size);

  {
    auto cache_indices_iter = std::move(cache_indices).begin();
    const auto cache_indices_end = std::move(cache_indices).end();
    auto cache_iter = std::move(cache).begin();
    auto cache_high_coverages_iter = std::move(cache_high_coverages).begin();
    auto new_cache_iter = ranges::begin(new_cache);
    auto new_cache_indices_iter = ranges::begin(new_cache_indices);
    auto new_cache_high_coverages_iter =
        ranges::begin(new_cache_high_coverages);
    auto valid_indices_iter = ranges::begin(valid_indices);

    for (windows_size_type window_index = 0;
         cache_indices_iter < cache_indices_end;
         ++cache_indices_iter, ++cache_iter, ++cache_high_coverages_iter,
                           ++window_index) {

      auto &&cache_indices_line = *cache_indices_iter;
      if (not cache_indices_line.empty()) {

        *new_cache_iter++ = *cache_iter;
        *new_cache_indices_iter++ = cache_indices_line;

        *new_cache_high_coverages_iter++ = *cache_high_coverages_iter;
        *valid_indices_iter++ = window_index;
      }
    }

    assert(cache_iter == std::move(cache).end());
    assert(cache_high_coverages_iter == std::move(cache_high_coverages).end());

    assert(new_cache_iter == ranges::end(new_cache));
    assert(new_cache_indices_iter == ranges::end(new_cache_indices));
    assert(new_cache_high_coverages_iter ==
           ranges::end(new_cache_high_coverages));
    assert(valid_indices_iter == ranges::end(valid_indices));
  }

  distances_type new_windows_distances(new_windows_size);

  {
    auto const begin_valid_indices = ranges::begin(valid_indices);
    auto const end_valid_indices = ranges::end(valid_indices);

    auto new_windows_distances_iter = ranges::begin(new_windows_distances);
    for (auto first_index_iter = begin_valid_indices;
         first_index_iter < end_valid_indices; ++first_index_iter) {
      auto const first_index = *first_index_iter;
      for (auto second_index_iter = ranges::next(first_index_iter);
           second_index_iter < end_valid_indices; ++second_index_iter) {
        auto const second_index = *second_index_iter;
        *new_windows_distances_iter++ =
            windows_distances[first_index][second_index];
      }
    }

    assert(new_windows_distances_iter == ranges::end(new_windows_distances));
  }

  cache = std::move(new_cache);
  cache_indices = std::move(new_cache_indices);
  windows_distances = std::move(new_windows_distances);
  cache_high_coverages = std::move(new_cache_high_coverages);

  assert(ranges::none_of(cache_indices,
                         [](auto &&indices) { return indices.empty(); }));
  assert(ranges::none_of(cache, [](auto &&window) { return window.empty(); }));
}

} // namespace windows_merger

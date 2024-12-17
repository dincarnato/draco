#include "windows_merger.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <numeric>
#include <oneapi/tbb/parallel_for.h>
#include <random>
#include <string_view>

#include "windows_merger_deserializer.hpp"

using namespace windows_merger;
namespace ranges = std::ranges;

static constexpr std::size_t n_tests = 20;
static constexpr std::size_t window_size = 70;
static constexpr std::size_t sequence_length = 300;

static std::tuple<WeightedClusters, std::vector<unsigned>,
                  typename WindowsMerger::bases_size_type>
generate_random_windows(typename WindowsMerger::clusters_size_type n_clusters,
                        float coverage_mean = 120.f,
                        float coverage_stddev = 20.) {
  using coverage_type = WindowsMergerTraits::coverage_type;

  thread_local std::mt19937 random_gen(std::random_device{}());
  thread_local std::uniform_real_distribution<double> real_dist(0., 1.);
  thread_local std::uniform_int_distribution<
      typename WindowsMerger::bases_size_type>
      seq_length_dist(0, sequence_length - 10);
  std::normal_distribution<float> normal_dist(coverage_mean, coverage_stddev);

  WindowsMerger::bases_size_type start_index = seq_length_dist(random_gen);
  auto current_window_size = static_cast<WindowsMerger::bases_size_type>(
      std::min(window_size, sequence_length - start_index));
  auto out =
      std::make_tuple(WeightedClusters(current_window_size, n_clusters, false),
                      std::vector<unsigned>(current_window_size), start_index);
  auto &weighted_clusters = std::get<0>(out);
  auto &coverages = std::get<1>(out);

  for (auto &&span : weighted_clusters) {
    std::vector<double> weights(n_clusters);
    ranges::generate(weights, [&] { return real_dist(random_gen); });
    const double normalizer = std::accumulate(std::ranges::begin(weights),
                                              std::ranges::end(weights), 0.);
    ranges::copy(weights | std::views::transform([normalizer](double value) {
                   return weighted_clusters_weight_type(value / normalizer);
                 }),
                 ranges::begin(span));
  }

  ranges::generate(coverages, [&] {
    return static_cast<coverage_type>(normal_dist(random_gen));
  });

  return out;
}

static void test_add_windows_mt() {
  constexpr std::size_t n_clusters = 4;

  WindowsMerger merger(n_clusters);
  tbb::parallel_for(std::size_t(0), n_tests, [&](std::size_t) {
    auto [weighted_clusters, coverages, start_index] =
        generate_random_windows(n_clusters);
    merger.add_window(start_index, std::move(weighted_clusters),
                      std::move(coverages));
  });
  merger.wait_queue();
}

static void test_windows_sorted() {
  constexpr std::size_t n_clusters = 3;
  constexpr std::size_t n_windows = 100;

  std::vector<WeightedClusters> all_weighted_clusters;
  std::vector<std::vector<unsigned>> all_coverages;
  std::vector<typename WindowsMerger::bases_size_type> all_starts;
  all_weighted_clusters.reserve(n_windows);
  all_coverages.reserve(n_windows);
  all_starts.reserve(n_windows);

  for (std::size_t window_index = 0; window_index < n_windows; ++window_index) {
    auto [weighted_clusters, coverages, start_index] =
        generate_random_windows(n_clusters);
    all_weighted_clusters.emplace_back(std::move(weighted_clusters));
    all_coverages.emplace_back(std::move(coverages));
    all_starts.emplace_back(start_index);
  }

  WindowsMerger merger(n_clusters);
  for (std::size_t window_index = 0; window_index < n_windows; ++window_index)
    merger.add_window(all_starts[window_index],
                      all_weighted_clusters[window_index],
                      all_coverages[window_index]);
  merger.wait_queue();

  {
    std::vector<std::size_t> indices(n_windows);
    std::iota(ranges::begin(indices), ranges::end(indices), std::size_t(0));
    ranges::sort(indices, ranges::less{},
                 [&](std::size_t index) { return all_starts[index]; });

    decltype(all_weighted_clusters) new_all_weighted_clusters;
    decltype(all_coverages) new_all_coverages;

    new_all_weighted_clusters.reserve(n_windows);
    new_all_coverages.reserve(n_windows);

    for (std::size_t window_index = 0; window_index < n_windows;
         ++window_index) {
      new_all_weighted_clusters.emplace_back(
          all_weighted_clusters[indices[window_index]]);
      new_all_coverages.emplace_back(all_coverages[indices[window_index]]);
    }

    all_weighted_clusters = std::move(new_all_weighted_clusters);
    all_coverages = std::move(new_all_coverages);

    ranges::sort(all_starts);
  }

  const auto &merger_windows = merger.get_windows();
  for (std::size_t window_index = 0; window_index < n_windows; ++window_index) {
    auto &&window_accessor = merger_windows[window_index];
    auto &&weighted_clusters = all_weighted_clusters[window_index];
    auto &&coverages = all_coverages[window_index];

    assert(window_accessor.begin_index() == all_starts[window_index]);

    {
      auto window_accessor_iter = ranges::begin(window_accessor);
      auto const window_accessor_end = ranges::end(window_accessor);
      auto weighted_clusters_iter = ranges::begin(weighted_clusters);
      auto const weighted_clusters_end = ranges::end(weighted_clusters);

      for (; window_accessor_iter < window_accessor_end and
             weighted_clusters_iter < weighted_clusters_end;
           ++window_accessor_iter, ++weighted_clusters_iter) {
        auto &&base_accessor = *window_accessor_iter;
        auto &&base_spans = *weighted_clusters_iter;

        auto &&weights_accessor = base_accessor.weights();
        auto weights_accessor_iter = ranges::begin(weights_accessor);
        auto const weights_accessor_end = ranges::end(weights_accessor);
        auto base_spans_iter = ranges::begin(base_spans);
        auto const base_spans_end = ranges::end(base_spans);

        for (; weights_accessor_iter < weights_accessor_end and
               base_spans_iter < base_spans_end;
             ++weights_accessor_iter, ++base_spans_iter) {
          assert(*weights_accessor_iter == TinyFraction(*base_spans_iter));
        }

        assert(weights_accessor_iter == weights_accessor_end);
        assert(base_spans_iter == base_spans_end);
      }
      assert(window_accessor_iter == window_accessor_end);
      assert(weighted_clusters_iter == weighted_clusters_end);
    }

    auto &&window_coverages = window_accessor.coverages();
    assert(ranges::equal(window_coverages, coverages));
  }
}

static WindowsMerger
create_random_merger(WindowsMerger::windows_size_type n_windows = 300,
                     WindowsMerger::clusters_size_type n_clusters = 5) {
  WindowsMerger merger(n_clusters);
  for (std::size_t window_index = 0; window_index < n_windows; ++window_index) {
    auto [weighted_clusters, coverages, start_index] =
        generate_random_windows(n_clusters);
    merger.add_window(start_index, std::move(weighted_clusters),
                      std::move(coverages));
  }
  merger.wait_queue();
  return merger;
}

static void test_prepare_cache() {
  WindowsMerger merger = create_random_merger();
  test::WindowsMerger::prepare_cache(merger);
  auto &&windows = test::WindowsMerger::windows(std::as_const(merger));
  auto &&cache = test::WindowsMerger::cache(std::as_const(merger));

  assert(ranges::equal(windows, cache));

  // Explicit coverages check
  {
    auto windows_iter = ranges::begin(windows);
    auto const windows_end = ranges::end(windows);
    auto cache_iter = ranges::begin(cache);
    auto const cache_end = ranges::end(cache);
    for (; windows_iter < windows_end and cache_iter < cache_end;
         ++windows_iter, ++cache_iter) {
      auto &&window = *windows_iter;
      auto &&cached_window = *cache_iter;

      assert(ranges::equal(window.coverages(), cached_window.coverages()));
    }
  }
}

static void test_prepare_indices() {
  WindowsMerger merger = create_random_merger();
  test::WindowsMerger::prepare_cache(merger);
  test::WindowsMerger::prepare_indices(merger);

  auto &&cache_indices =
      test::WindowsMerger::cache_indices(std::as_const(merger));
  typename WindowsMerger::windows_size_type window_index = 0;
  for (auto &&window : cache_indices) {
    assert(window.size() == 1);
    assert(window[0] == window_index++);
  }
}

static void test_basic_prepare_distances() {
  using windows_size_type = typename WindowsMerger::windows_size_type;
  using distance_type = typename WindowsMerger::distance_type;

  WindowsMerger merger = create_random_merger(300, 3);
  test::WindowsMerger::prepare_cache(merger);
  test::WindowsMerger::prepare_indices(merger);
  test::WindowsMerger::prepare_high_coverages(merger);
  test::WindowsMerger::prepare_distances(merger);

  auto &&windows_distances =
      test::WindowsMerger::windows_distances(std::as_const(merger));
  auto &&cache = test::WindowsMerger::cache(std::as_const(merger));

  for (windows_size_type first_window_index = 0;
       first_window_index < cache.windows_size(); ++first_window_index) {
    auto &&first_cache_window = std::as_const(cache)[first_window_index];
    auto &&window_distances =
        std::as_const(windows_distances)[first_window_index];

    for (windows_size_type second_window_index = first_window_index + 1;
         second_window_index < cache.windows_size(); ++second_window_index) {
      auto &&second_cache_window = cache[second_window_index];
      auto &&distance = window_distances[second_window_index];

      if (first_cache_window.begin_index() >= second_cache_window.end_index() or
          second_cache_window.begin_index() >= first_cache_window.end_index())
        assert(distance == std::numeric_limits<distance_type>::infinity());
      else
        assert(distance >= 0);
    }
  }
}

bool operator==(
    typename WindowsMergerWindows::const_window_accessor merger_window,
    Window const &window) noexcept {
  using bases_size_type = WindowsMergerTraits::bases_size_type;
  using clusters_size_type = WindowsMergerTraits::clusters_size_type;

  if (window.clusters_weights.size() != merger_window.clusters_size() or
      window.begin != merger_window.begin_index() or
      window.end != merger_window.end_index())
    return false;

  for (bases_size_type base_index = 0; base_index < merger_window.size();
       ++base_index) {
    for (clusters_size_type cluster_index = 0;
         cluster_index < merger_window.clusters_size(); ++cluster_index) {

      auto const &cluster = window.clusters_weights[cluster_index];
      if (cluster.size() != merger_window.size())
        return false;
      TinyFraction fraction(cluster[base_index]);

      if (merger_window[base_index].weight(cluster_index) != fraction)
        return false;
    }
  }

  return ranges::equal(merger_window.coverages(), window.coverages);
}

bool operator!=(
    typename WindowsMergerWindows::const_window_accessor merger_window,
    Window const &window) noexcept {
  return not(merger_window == window);
}

bool operator==(WindowsMergerWindows const &merger_windows,
                std::vector<Window> const &windows) noexcept {
  if (merger_windows.windows_size() != windows.size())
    return false;

  auto merger_window_iter = merger_windows.begin();
  auto const merger_windows_end = merger_windows.end();
  auto window_iter = ranges::begin(windows);
  for (; merger_window_iter < merger_windows_end;
       ++merger_window_iter, ++window_iter) {
    if (*merger_window_iter != *window_iter)
      return false;
  }

  return true;
}

static void test_from_serialized_data(char const *filename) {
  using windows_size_type = WindowsMergerTraits::windows_size_type;
  using clusters_size_type = WindowsMergerTraits::clusters_size_type;
  namespace jsonde = json_deserializer;
  using namespace std::string_view_literals;

  std::ifstream data_file(filename);

  if (!data_file.good())
    throw std::runtime_error("cannot open serialized data file");

  clusters_size_type n_clusters;
  windows_size_type n_windows;
  unsigned short window_size;

  std::string line;

  {
    assert(std::getline(data_file, line));
    auto const first_space_pos = line.find(' ');
    assert(first_space_pos != std::string::npos);
    assert(first_space_pos > 0);
    assert(std::string_view(line.c_str() + first_space_pos) == " clusters");
    /*
    assert(std::from_chars(line.c_str(), line.c_str() + first_space_pos,
                           n_clusters)
               .ptr == line.c_str() + first_space_pos);
               */
    n_clusters = static_cast<clusters_size_type>(
        std::stoi(std::string(line.substr(0, first_space_pos))));
  }

  {
    assert(std::getline(data_file, line));
    auto const first_space_pos = line.find(' ');
    assert(first_space_pos != std::string::npos);
    assert(first_space_pos > 0);
    assert(std::string_view(line.c_str() + first_space_pos) == " windows");
    /*
    assert(
        std::from_chars(line.c_str(), line.c_str() + first_space_pos, n_windows)
            .ptr == line.c_str() + first_space_pos);
            */
    n_windows = static_cast<windows_size_type>(
        std::stoi(std::string(line.substr(0, first_space_pos))));
  }

  {
    assert(std::getline(data_file, line));
    auto const delimiter_pos = line.find('=');
    assert(delimiter_pos != std::string::npos);
    assert(std::string_view(line.c_str(), delimiter_pos) == "window_size");
    /*
    assert(std::from_chars(line.c_str() + delimiter_pos + 1,
                           line.c_str() + line.size(), window_size)
               .ptr == line.c_str() + line.size());
               */
    window_size = static_cast<unsigned short>(
        std::stoi(std::string(line.substr(delimiter_pos + 1))));
  }

  WindowsMerger windows_merger(n_clusters);

  for (unsigned window_index = 0; window_index < n_windows; ++window_index) {
    auto [weighted_clusters, coverages, start_base_index] =
        deserialize_initial_window(data_file, n_clusters, window_size);
    windows_merger.add_window(start_base_index, std::move(weighted_clusters),
                              std::move(coverages));
  }

  assert(std::getline(data_file, line));
  assert(line.empty());
  assert(std::getline(data_file, line));
  assert(line == "after sorting:");

  windows_merger.wait_queue();

  {
    auto windows = deserialize_windows(data_file, n_clusters);
    assert(test::WindowsMerger::windows(windows_merger) == windows);
  }
  assert(std::getline(data_file, line));
  assert(line.empty());

  test::WindowsMerger::prepare_cache(windows_merger);
  test::WindowsMerger::prepare_indices(windows_merger);
  test::WindowsMerger::prepare_high_coverages(windows_merger);
  test::WindowsMerger::prepare_distances(windows_merger);

  assert(std::getline(data_file, line));
  assert(line == "initial distances:");
  {
    auto const distances = [&] {
      jsonde::Array raw_distances;
      deserialize(data_file, raw_distances);

      std::vector<double> distances(raw_distances.size());
      ranges::transform(raw_distances, ranges::begin(distances),
                        [](jsonde::Value const &json_value) {
                          double value;
                          from_json_number_value(value, json_value);
                          return value;
                        });
      return distances;
    }();

    auto distance_iter = ranges::begin(distances);
    auto const distances_end = ranges::end(distances);
    for (windows_size_type window_index_a = 0; window_index_a < n_windows;
         ++window_index_a) {

      for (windows_size_type window_index_b = window_index_a + 1;
           window_index_b < n_windows; ++window_index_b, ++distance_iter) {
        assert(distance_iter < distances_end);

        if (std::isnan(*distance_iter)) {
          auto &&windows = std::as_const(windows_merger).get_windows();
          auto &&window_a = windows[window_index_a];
          auto &&window_b = windows[window_index_b];

          assert(window_a.end_index() <= window_b.begin_index());
        } else {
          double best_distance =
              test::WindowsMerger::get_best_distance_and_permutation_between(
                  windows_merger, window_index_a, window_index_b)
                  .first;
          assert(std::abs(best_distance - *distance_iter) < 1e-2);
        }
      }
    }

    assert(distance_iter == distances_end);
  }
  assert(std::getline(data_file, line));
  assert(line.empty());

  while (std::getline(data_file, line)) {
    if (line == "final window:")
      break;

    assert(std::string_view(line.c_str(), "best pair: "sv.size()) ==
           "best pair: ");
    std::array<windows_size_type, 2> best_pair_indices;
    {
      std::string_view raw_best_pairs(line.c_str() + "best pair: "sv.size());
      auto space_pos = raw_best_pairs.find(' ');
      assert(space_pos != std::string_view::npos);

      /*
      assert(std::from_chars(ranges::begin(raw_best_pairs),
                             std::next(ranges::begin(raw_best_pairs),
                                       static_cast<std::ptrdiff_t>(space_pos)),
                             best_pair_indices[0])
                 .ptr == std::next(ranges::begin(raw_best_pairs),
                                   static_cast<std::ptrdiff_t>(space_pos)));
      assert(
          std::from_chars(std::next(ranges::begin(raw_best_pairs),
                                    static_cast<std::ptrdiff_t>(space_pos + 1)),
                          ranges::end(raw_best_pairs), best_pair_indices[1])
              .ptr == ranges::end(raw_best_pairs));
      */
      best_pair_indices[0] = static_cast<windows_size_type>(
          std::stoi(std::string(raw_best_pairs.substr(0, space_pos))));
      best_pair_indices[1] = static_cast<windows_size_type>(
          std::stoi(std::string(raw_best_pairs.substr(space_pos + 1))));
    }

    assert(std::getline(data_file, line));
    assert(line == "merged window:");

    auto merged_window = deserialize_window(data_file, n_clusters);

    assert(std::getline(data_file, line));
    assert(line.empty());

    auto const &merger_cache = test::WindowsMerger::cache(windows_merger);
    assert(merger_cache.windows_size() > 1);
    auto const &merger_best_pair_indices =
        test::WindowsMerger::find_best_cached_pair(windows_merger);
    assert(merger_best_pair_indices == best_pair_indices);

    auto first_window_index = merger_best_pair_indices[0];
    auto const old_windows_size = merger_cache.windows_size();
    auto &&cache_indices =
        test::WindowsMerger::cache_indices(std::as_const(windows_merger));
    auto const empty_windows = ranges::count_if(
        ranges::begin(cache_indices),
        ranges::next(ranges::begin(cache_indices),
                     static_cast<std::ptrdiff_t>(first_window_index)),
        [](auto &&indices) { return indices.empty(); });
    test::WindowsMerger::merge_cached_windows_into_first(
        windows_merger, first_window_index, merger_best_pair_indices[1]);

    if (merger_cache.windows_size() != old_windows_size) {
      first_window_index -= static_cast<windows_size_type>(empty_windows);
    }
    assert(merger_cache[first_window_index] == merged_window);
  }

  auto final_window = deserialize_window(data_file, n_clusters);
  auto const &merger_best_pair_indices =
      test::WindowsMerger::find_best_cached_pair(windows_merger);
  assert(merger_best_pair_indices[1] ==
         std::numeric_limits<windows_size_type>::max());
  auto const last_window_index = merger_best_pair_indices[0];
  auto const &merger_cache = test::WindowsMerger::cache(windows_merger);
  assert(not merger_cache[last_window_index].empty());
  assert(ranges::count_if(merger_cache, [](auto &&window) {
           return not window.empty();
         }) == 1);
  assert(merger_cache[last_window_index] == final_window);

  WindowsMergerWindow const last_window_output =
      merger_cache[last_window_index];
  assert(last_window_output == final_window);
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " serialized_data.txt\n";
    return 2;
  }

  test_add_windows_mt();
  test_windows_sorted();
  test_prepare_cache();
  test_prepare_indices();
  test_basic_prepare_distances();
  test_from_serialized_data(argv[1]);
}

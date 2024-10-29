#pragma once

#include "../ios_saver.hpp"
#include "../ringmap_data.hpp"
#include "../ringmap_matrix_row.hpp"
#include "../weighted_clusters.hpp"

#include "jsonify_base.hpp"

#include <algorithm>
#include <iomanip>
#include <vector>

namespace windows_merger {
struct WindowsMergerWindow;
} // namespace windows_merger

namespace results {

struct Window {
  using clusters_fraction_type = RingmapData::clusters_fraction_type;
  using clusters_pattern_type = RingmapData::clusters_pattern_type;
  using clusters_assignment_type = RingmapData::clusters_assignment_type;

  Window() = default;
  Window(windows_merger::WindowsMergerWindow const &wmw,
         std::vector<unsigned> const &coverages) noexcept(false);
  Window(unsigned short begin_index, WeightedClusters const &weighted_clusters,
         std::vector<unsigned> const &coverages) noexcept(false);

  unsigned short begin_index, end_index;
  WeightedClusters weighted_clusters;
  clusters_fraction_type fractions;
  std::optional<clusters_pattern_type> patterns;
  std::optional<std::vector<std::vector<unsigned>>> bases_coverages;
  std::vector<std::int8_t> assignments;
  std::vector<unsigned> coverages;
};

} // namespace results

#include "../windows_merger_window.hpp"

namespace results {

inline Window::Window(windows_merger::WindowsMergerWindow const &wmw,
                      std::vector<unsigned> const &coverages) noexcept(false)
    : begin_index(wmw.begin_index()), end_index(wmw.end_index()),
      weighted_clusters(wmw.size(), wmw.clusters_size(), false),
      coverages(coverages) {
  auto wmw_iter = std::begin(wmw);
  auto const end_wmw = std::end(wmw);
  auto base_weights_iter = std::begin(weighted_clusters);

  for (; wmw_iter < end_wmw; ++wmw_iter, ++base_weights_iter) {
    auto &&wmw_base = *wmw_iter;
    auto &&wmw_weights = wmw_base.weights();
    auto &&base_weights = *base_weights_iter;
    [[maybe_unused]] auto end_base_weights =
        std::copy(std::begin(wmw_weights), std::end(wmw_weights),
                  std::begin(base_weights));
    assert(end_base_weights == std::end(base_weights));
  }
}

inline Window::Window(unsigned short begin_index,
                      WeightedClusters const &weighted_clusters,
                      std::vector<unsigned> const &coverages) noexcept(false)
    : begin_index(begin_index),
      end_index(begin_index + static_cast<unsigned short>(coverages.size())),
      weighted_clusters(weighted_clusters), coverages(coverages) {
  if (begin_index + coverages.size() >
      std::numeric_limits<unsigned short>::max()) {
    throw std::runtime_error(
        "end_index cannot be represented with an unsigned short");
  }
  auto const weighted_clusters_size = weighted_clusters.getElementsSize();
  if (weighted_clusters_size != 0 and
      weighted_clusters_size != coverages.size()) {
    throw std::runtime_error("Weights have a different size than coverages");
  }
}

template <typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os,
        WeightedClustersClusterWrapper<T> const &cluster_wrapper) {
  os << '[';
  auto iter = std::ranges::begin(cluster_wrapper);
  if (iter != std::ranges::end(cluster_wrapper)) {
    IosSaver iosSaver(os);
    os << std::fixed << std::setprecision(3);

    jsonify(os, *iter++);
    std::ranges::for_each(iter, std::ranges::end(cluster_wrapper),
                          [&](auto &&weight) {
                            os << ',';
                            jsonify(os, weight);
                          });
  }
  return os << ']';
}

template <typename CharT, typename Traits, typename T>
inline std::enable_if_t<std::is_convertible_v<std::decay_t<T>, Window>,
                        std::basic_ostream<CharT, Traits> &>
jsonify(std::basic_ostream<CharT, Traits> &os, T &&window) {
  os << '{';
  {
    IosSaver iosSaver(os);
    os << std::fixed << std::setprecision(3);

    jsonify(os, "start", window.begin_index, "end", window.end_index,
            "stoichiometries", window.fractions);

    if (window.patterns) {
      os << ',';
      jsonify(os, "counts", window.patterns);
    }

    if (window.bases_coverages) {
      os << ',';
      jsonify(os, "coverage", window.bases_coverages);
    }

    os << ',';
    jsonify(os, "preCoverage", window.coverages);
  }

  os << ",";
  jsonify(os, "weights") << ":[";
  auto &&clusters = window.weighted_clusters.clusters();
  auto iter = std::ranges::begin(clusters);
  if (iter != std::ranges::end(clusters)) {
    IosSaver iosSaver(os);
    os << std::fixed << std::setprecision(3);

    jsonify(os, *iter++);
    std::ranges::for_each(iter, std::ranges::end(clusters),
                          [&](const auto &cluster) {
                            os << ',';
                            jsonify(os, cluster);
                          });
  }

  return os << "]}";
}

template <typename CharT, typename Traits>
std::basic_ostream<CharT, Traits> &
jsonify(std::basic_ostream<CharT, Traits> &os,
        RingmapMatrixRow const &ringmap_row) {
  os << '{';
  jsonify(os, "start", ringmap_row.begin_index(), "end",
          ringmap_row.end_index(), "mutatedIndexes",
          static_cast<RingmapMatrixRow::base_type const &>(ringmap_row));
  return os << "}";
}

} // namespace results

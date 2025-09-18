#pragma once

#include <algorithm>
#include <armadillo>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iterator>
#include <random>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

namespace kmeans {

#define MAKE_DISTANCE_SQUARED(ty)                                              \
  inline ty distance_squared(ty a, ty b) noexcept {                            \
    auto diff = a - b;                                                         \
    return diff * diff;                                                        \
  }

MAKE_DISTANCE_SQUARED(double)

#define MAKE_DISTANCE_SQUARED_VEC(ty)                                          \
  template <typename T>                                                        \
  inline T distance_squared(ty const &a, ty const &b) noexcept {               \
    auto diff = a - b;                                                         \
    return arma::dot(diff, diff);                                              \
  }

MAKE_DISTANCE_SQUARED_VEC(arma::Col<T>);
MAKE_DISTANCE_SQUARED_VEC(arma::subview_row<T>);

template <typename T>
inline T distance_squared(arma::subview_row<T> const &a,
                          arma::Col<T> const &b) {
  auto diff = a.t() - b;
  return arma::dot(diff, diff);
}

constexpr double MAX_CENTROID_DIMENSION_DISTANCE = 1e-5;

template <typename T, typename Gen>
void kmeans_pp_centroids_initialization(arma::subview<T> const &points,
                                        std::uint8_t n_clusters,
                                        std::vector<arma::uword> &indices,
                                        Gen &&random_generator) {
  using distance_t = decltype(distance_squared(points[0], points[0]));

  indices.resize(points.n_rows);
  std::ranges::iota(indices, static_cast<arma::uword>(0));

  std::ranges::shuffle(indices, random_generator);

  if (std::size(indices) <= n_clusters) {
    return;
  }
  std::vector<distance_t> all_centroids_distances(std::size(indices) - 1);

  for (std::size_t centroid_index = 1; centroid_index < n_clusters;
       ++centroid_index) {
    std::span set_centroids_indices(std::ranges::begin(indices),
                                    centroid_index);
    std::span points_indices(indices | std::views::drop(centroid_index));
    std::span centroids_distance(all_centroids_distances |
                                 std::views::drop(centroid_index - 1));

    for (auto &&[point_index, centroid_distance] :
         std::views::zip(points_indices, centroids_distance)) {
      auto const &point = points.row(point_index);

      centroid_distance = std::ranges::min(
          set_centroids_indices |
          std::views::transform([&](auto centroid_index) {
            return distance_squared(point, points.row(centroid_index));
          }));
    }

    auto centroids_distances_sum =
        std::ranges::fold_left(centroids_distance, static_cast<distance_t>(0),
                               std::plus<distance_t>());
    std::uniform_real_distribution<distance_t> dist(0.,
                                                    centroids_distances_sum);
    auto threshold = dist(random_generator);
    auto cumulative_distance = static_cast<distance_t>(0);
    for (auto &&[point_index, centroid_distance] :
         std::views::zip(points_indices, all_centroids_distances)) {
      cumulative_distance += centroid_distance;
      if (cumulative_distance >= threshold) {
        std::swap(point_index, indices[centroid_index]);
        break;
      }
    }
  }

  indices.resize(n_clusters);
}

template <typename T, typename Gen>
std::vector<arma::uword>
kmeans_pp_centroids_initialization(arma::subview<T> const &points,
                                   std::uint8_t n_clusters,
                                   Gen &&random_generator) {
  std::vector<arma::uword> output;
  kmeans_pp_centroids_initialization(points, n_clusters, output,
                                     std::forward<Gen>(random_generator));
  return output;
}

template <typename T> struct Result {
  std::vector<std::vector<arma::uword>> clusters_indices;
  arma::Mat<T> centroids;
};

} // namespace kmeans

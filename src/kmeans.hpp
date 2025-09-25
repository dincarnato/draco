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

template <typename T, typename Gen>
void single_iteration(arma::subview<T> const &points, std::uint8_t n_clusters,
                      Result<T> &result,
                      std::vector<arma::uword> &centroids_indices,
                      Gen &&random_generator) {
  auto &centroids = result.centroids;
  auto &assignments = result.clusters_indices;

  kmeans_pp_centroids_initialization(points, n_clusters, centroids_indices,
                                     random_generator);
  centroids.resize(std::size(centroids_indices), points.n_cols);

  for (auto &&[row_index, centroid_index] : std::views::zip(
           std::views::iota(static_cast<arma::uword>(0)), centroids_indices)) {
    centroids.row(row_index) = points.row(centroid_index);
  }

  assignments.resize(n_clusters);

  // TODO: implement history check to avoid infinite iterations and/or
  // maximum number of iterations
  for (;;) {
    for (auto &assignment : assignments) {
      assignment.clear();
    }

    for (arma::uword point_index = 0; point_index < points.n_rows;
         ++point_index) {
      auto &&point = points.row(point_index);
      auto best_centroid_index = *std::ranges::min_element(
          std::views::iota(static_cast<arma::uword>(0),
                           static_cast<arma::uword>(n_clusters)),
          {}, [&](auto index) {
            auto centroid = centroids.row(index);
            return distance_squared(point, centroid);
          });
      assignments[static_cast<std::size_t>(best_centroid_index)].push_back(
          point_index);
    };

    bool finished = true;
    arma::Row<T> new_centroid;
    for (auto &&[index, assignment] : std::views::zip(
             std::views::iota(static_cast<arma::uword>(0)), assignments)) {
      auto centroid = centroids.row(index);
      new_centroid.resize(points.n_cols);
      new_centroid.fill(static_cast<T>(0));

      for (auto point_index : assignment) {
        auto &&point = points.row(point_index);
        new_centroid += point;
      }
      new_centroid *= static_cast<T>(1) / static_cast<T>(std::size(assignment));

      if (arma::any(arma::abs(new_centroid - centroid) >
                    static_cast<T>(MAX_CENTROID_DIMENSION_DISTANCE))) {
        centroid = new_centroid;
        finished = false;
      }
    }

    if (finished) {
      break;
    }
  }
}

template <typename T>
T result_cumulative_centroids_distance_squared(arma::subview<T> const &points,
                                               Result<T> const &result) {
  return std::ranges::fold_left(
      std::views::zip(result.clusters_indices,
                      std::views::iota(static_cast<arma::uword>(0)) |
                          std::views::transform([&](auto row_index) {
                            return result.centroids.row(row_index);
                          })) |
          std::views::transform([&](auto &&tuple) {
            auto &&[cluster_indices, centroid] = tuple;
            return std::ranges::fold_left(
                cluster_indices | std::views::transform([&](auto point_index) {
                  return distance_squared(points.row(point_index), centroid);
                }),
                static_cast<T>(0), std::plus<T>{});
          }),
      static_cast<T>(0), std::plus<T>{});
}

template <typename T, typename Gen>
Result<T> run(arma::subview<T> const &points, std::uint8_t n_clusters,
              std::uint16_t iterations, Gen &&random_generator) {
  iterations = std::max(iterations, static_cast<std::uint16_t>(1));

  if (n_clusters <= 1) {
    std::vector<std::vector<arma::uword>> clusters_indices{
        std::views::iota(static_cast<arma::uword>(0), points.n_rows) |
        std::ranges::to<std::vector>()};

    arma::Row<T> centroid(points.n_cols);
    for (arma::uword point_index = 0; point_index < points.n_rows;
         ++point_index) {
      auto &&point = points.row(point_index);
      centroid += point;
    }
    centroid *= static_cast<T>(1) / static_cast<T>(points.n_rows);
    return Result{.clusters_indices = clusters_indices,
                  .centroids = arma::mat{centroid}};
  }

  Result<T> best_result;
  std::vector<arma::uword> centroid_indices;
  single_iteration(points, n_clusters, best_result, centroid_indices,
                   random_generator);
  auto best_centroids_distance =
      result_cumulative_centroids_distance_squared(points, best_result);

  Result<T> result;
  for (std::uint16_t iteration = 1; iteration < iterations; ++iteration) {
    single_iteration(points, n_clusters, result, centroid_indices,
                     random_generator);
    auto centroids_distance =
        result_cumulative_centroids_distance_squared(points, result);
    if (centroids_distance < best_centroids_distance) {
      std::swap(best_result, result);
      best_centroids_distance = centroids_distance;
    }
  }

  return best_result;
}

} // namespace kmeans

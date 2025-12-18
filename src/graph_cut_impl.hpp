#pragma once

#include "clusters_replicates.hpp"
#include "graph_cut.hpp"
#include "kmeans.hpp"
#include "results/transcript.hpp"

#include <algorithm>
#include <armadillo>
#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <ranges>
#include <sstream>
#include <type_traits>
#include <utility>
#include <vector>

template <typename Fun, typename Ret = decltype(std::declval<Fun>()(
                            std::declval<const arma::mat &>()))>
static constexpr void checkGraphFunCallable() {
  static_assert(std::is_same<Ret, arma::mat>::value,
                "Fun must be of type arma::mat (const arma::mat&)");
}

template <typename Fun, typename Gen>
  requires requires(Fun fun) {
    { fun(std::declval<arma::mat const &>()) } -> std::same_as<arma::mat>;
  }
WeightedClusters GraphCut::partitionGraph(std::uint8_t nClusters,
                                          std::uint16_t kmeans_iterations,
                                          results::Transcript const &transcript,
                                          unsigned window_index, Fun graphFun,
                                          double distance_warning_threshold,
                                          Gen &&random_generator) const {
  auto const &firstAdjacency = adjacencies[0];
  assert(std::ranges::all_of(adjacencies | std::views::drop(1),
                             [&](auto const &adjacency) {
                               return adjacency.n_rows == firstAdjacency.n_rows;
                             }));

  checkGraphFunCallable<std::decay_t<Fun>>();

  if (nClusters < 2) {
    return WeightedClusters(firstAdjacency.n_rows, 1);
  }

  if (std::size(adjacencies) == 1) {
    return weighted_clusters_from_adjacency(nClusters, kmeans_iterations,
                                            graphFun, random_generator,
                                            adjacencies[0]);
  } else {
    auto all_weighted_clusters =
        adjacencies | std::views::transform([&](auto const &adjacency) {
          return weighted_clusters_from_adjacency(nClusters, kmeans_iterations,
                                                  graphFun, random_generator,
                                                  adjacency);
        }) |
        std::views::as_rvalue | std::ranges::to<std::vector>();

    clusters_replicates::reorder_best_permutation(all_weighted_clusters,
                                                  transcript, window_index,
                                                  distance_warning_threshold);
    return merge_weighted_clusters(std::move(all_weighted_clusters));
  }
}

template <typename Gen>
arma::mat eigenvectors_to_weighted_clusters(arma::mat const &eigenvectors,
                                            std::uint8_t n_clusters,
                                            std::uint16_t kmeans_iterations,
                                            Gen &&random_generator) {
  auto const n_bases = eigenvectors.n_rows;

  if (n_clusters == 2) {
    auto fiedler = eigenvectors.col(1);
    auto min = fiedler.min();
    auto max = fiedler.max();

    auto denominator = max - min;
    arma::mat weighted_clusters(n_clusters, n_bases);
    if (denominator < 1e-6) {
      weighted_clusters.fill(0.5);
    } else {
      weighted_clusters.row(0) = (fiedler.t() - min) * (1. / denominator);
      weighted_clusters.row(1) = 1. - weighted_clusters.row(0);
    }

    return weighted_clusters;
  }

  auto useful_eigenvecs =
      eigenvectors.submat(arma::span::all, arma::span(1, n_clusters - 1));
  useful_eigenvecs = arma::normalise(useful_eigenvecs, 2, 1);

  auto centroids = std::move(kmeans::run(useful_eigenvecs, n_clusters,
                                         kmeans_iterations, random_generator)
                                 .centroids);
  auto weights = pairwise_distances(centroids, useful_eigenvecs);
  weights = arma::exp(-weights);
  weights.clean(0.);
  weights.each_col([](arma::vec &col) { col /= arma::sum(col); });
  return weights;
}

template <typename Fun, typename Gen>
  requires requires(Fun fun) {
    { fun(std::declval<arma::mat const &>()) } -> std::same_as<arma::mat>;
  }
WeightedClusters weighted_clusters_from_adjacency(
    std::uint8_t n_clusters, std::uint16_t kmeans_iterations, Fun graph_fun,
    Gen &&random_generator, arma::mat const &adjacency) {
  auto graph = graph_fun(adjacency);
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, graph);

  // TODO: avoid allocating a matrix and then create a WeightedClusters.
  // Try to work directly with a WeightedClusters instead.
  auto weighted_clusters_mat = eigenvectors_to_weighted_clusters(
      eigenvectors, n_clusters, kmeans_iterations, random_generator);

  WeightedClusters weighted_clusters(weighted_clusters_mat.n_cols,
                                     weighted_clusters_mat.n_rows, false);
  std::ranges::for_each(
      std::views::zip(std::views::iota(static_cast<arma::uword>(0)) |
                          std::views::transform([&](auto cluster_index) {
                            return weighted_clusters_mat.row(cluster_index);
                          }),
                      weighted_clusters.clusters()),
      [](auto tuple) {
        auto &&[cluster_vec, cluster] = tuple;
        for (auto &&[element_index, weight] : std::views::zip(
                 std::views::iota(static_cast<arma::uword>(0)), cluster)) {
          weight = static_cast<float>(cluster_vec[element_index]);
        }
      });

  return weighted_clusters;
}

#pragma once

#include "clusters_replicates.hpp"
#include "clusters_traits.hpp"
#include "graph_cut.hpp"
#include "kmeans.hpp"
#include "logger.hpp"

#include <algorithm>
#include <armadillo>
#include <array>
#include <cassert>
#include <concepts>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <random>
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
                                          Fun graphFun,
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

    clusters_replicates::reorder_best_permutation(all_weighted_clusters);
    return merge_weighted_clusters(std::move(all_weighted_clusters));
  }
}

template <typename Fun>
HardClusters GraphCut::partitionGraphHard(std::uint8_t nClusters,
                                          Fun graphFun) const {
  checkGraphFunCallable<std::decay_t<Fun>>();
  auto const &firstAdjacency = adjacencies[0];
  assert(std::ranges::all_of(adjacencies | std::views::drop(1),
                             [&](auto const &adjacency) {
                               return adjacency.n_rows == firstAdjacency.n_rows;
                             }));

  auto graphs = adjacencies | std::views::transform([&](auto const &adjacency) {
                  return graphFun(adjacency);
                }) |
                std::ranges::to<std::vector>();

  HardClusters hardClusters = [&, this] {
    if (auto init = std::get_if<HardClusters>(&initialClusters))
      return *init;
    else {
      HardClusters hardClusters(firstAdjacency.n_rows, nClusters);
      std::mt19937 randomGen(std::random_device{}());
      std::uniform_int_distribution<std::uint8_t> clusterAssigner(0, nClusters -
                                                                         1);
      auto hardClustersWrapper = hardClusters.clusters();
      do {
        hardClusters.clear();

        for (unsigned baseIndex = 0; baseIndex < firstAdjacency.n_rows;
             ++baseIndex)
          hardClusters[baseIndex].set(clusterAssigner(randomGen));
      } while (std::ranges::any_of(hardClustersWrapper,
                                   [](const auto &cluster) {
                                     return std::ranges::count(cluster, true) <
                                            3;
                                   }) or
               std::ranges::any_of(graphs, [&](auto const &graph) {
                 return std::isinf(calculateCutScore(graph, hardClusters));
               }));

      return hardClusters;
    }
  }();

  auto clusters = hardClusters.clusters();

  for (const auto clustersEnd = std::ranges::end(clusters);;) {
    double bestScore = std::ranges::fold_left(
        graphs | std::views::transform([&](auto const &graph) {
          return calculateCutScore(graph, hardClusters);
        }),
        0., std::plus{});
    if (std::isnan(bestScore)) {
      std::stringstream graphs_content;

      auto graphs_to_dump =
          graphs | std::views::transform([&](auto const &graph) {
            auto score = calculateCutScore(graph, hardClusters);
            return std::pair(&graph, score);
          }) |
          std::views::filter(
              [](auto const &pair) { return std::isnan(std::get<1>(pair)); }) |
          std::views::transform([](auto &&pair) { return std::get<0>(pair); });
      (*std::ranges::begin(graphs_to_dump))->print(graphs_content);
      for (auto graph_to_dump : graphs_to_dump | std::views::drop(1)) {
        graphs_content << "\n\n";
        graph_to_dump->print(graphs_content);
      }

      logger::error("Invalid bestScore. Graphs content:\n{}",
                    graphs_content.str());
      std::abort();
    }

    auto bestFromClusterIter = clustersEnd;
    auto bestToClusterIter = bestFromClusterIter;
    // FIXME GCC warns about this that could be used uninitialized
    // (theoretically, it cannot happens...). Let's init it...
    auto bestBaseIndex = std::numeric_limits<std::size_t>::max();

    for (auto fromClusterIter = std::ranges::begin(clusters);
         fromClusterIter < clustersEnd; ++fromClusterIter) {
      auto fromCluster = *fromClusterIter;
      if (std::ranges::count(fromCluster, true) <= 3)
        continue;

      const auto fromClusterEnd = std::ranges::end(fromCluster);

      for (auto toClusterIter = std::ranges::begin(clusters);
           toClusterIter < clustersEnd; ++toClusterIter) {
        if (fromClusterIter == toClusterIter)
          continue;

        auto toCluster = *toClusterIter;

        auto fromClusterCurrentIter = std::ranges::begin(fromCluster);
        auto toClusterCurrentIter = std::ranges::begin(toCluster);
        std::size_t baseIndex = 0;
        for (; fromClusterCurrentIter < fromClusterEnd;
             ++fromClusterCurrentIter, ++toClusterCurrentIter, ++baseIndex) {
          auto fromClusterCurrent = *fromClusterCurrentIter;
          if (not fromClusterCurrent)
            continue;

          auto toClusterCurrent = *toClusterCurrentIter;
          if (toClusterCurrent)
            continue;

          fromClusterCurrent.swap(toClusterCurrent);

          double currentScore = std::ranges::fold_left(
              graphs | std::views::transform([&](auto const &graph) {
                return calculateCutScore(graph, hardClusters);
              }),
              0., std::plus{});
          assert(not std::isnan(currentScore));

          if (not std::isinf(currentScore) and currentScore < bestScore) {
            bestScore = currentScore;
            bestFromClusterIter = fromClusterIter;
            bestToClusterIter = toClusterIter;
            bestBaseIndex = baseIndex;
          }

          fromClusterCurrent.swap(toClusterCurrent);
        }
      }
    }

    if (bestFromClusterIter == clustersEnd)
      return hardClusters;

    {
      auto bestFromCluster = *bestFromClusterIter;
      auto bestToCluster = *bestToClusterIter;

      // FIXME this should not assert, but it is related to the suspicious
      // warning of GCC
      assert(bestBaseIndex != std::numeric_limits<std::size_t>::max());
      bestFromCluster[bestBaseIndex].swap(bestToCluster[bestBaseIndex]);
    }
  }
}

template <typename Clusters>
double GraphCut::calculateCutScore(const arma::mat &graph,
                                   const Clusters &clusters) const {
  static_assert(is_clusters_v<Clusters>);
  double score = 0;
  auto complement = clusters.complement();
  auto clustersWrapper = clusters.clusters();
  for (const auto &cluster : clustersWrapper) {
    double divisor = calculateCut(graph, cluster, cluster);
    if (divisor == 0.)
      return std::numeric_limits<double>::infinity();

    score += calculateCut(graph, cluster, complement.cluster(cluster.index())) /
             divisor;
    assert(not std::isnan(score));
  }

  return score;
}

template <typename ClusterA, typename ClusterB>
double GraphCut::calculateCut(const arma::mat &graph, const ClusterA &clusterA,
                              const ClusterB &clusterB) {
  static_assert(is_cluster_wrapper_v<ClusterA>);
  static_assert(is_cluster_wrapper_v<ClusterB>);
  double score = 0;
  using value_type = typename ClusterA::iterator::value_type;

  static constexpr bool valueIsIntegral = std::is_integral_v<value_type>;
  static_assert(valueIsIntegral ==
                std::is_integral_v<typename ClusterB::iterator::value_type>);

  std::size_t indexA = 0;
  for (auto &&valueA : clusterA) {
    const bool skip = [&] {
      if constexpr (valueIsIntegral)
        return not valueA;
      else
        return valueA == 0;
    }();

    if (skip) {
      ++indexA;
      continue;
    }

    const auto &colA = graph.col(indexA);
    std::size_t indexB = 0;
    for (auto &&valueB : clusterB) {
      if constexpr (valueIsIntegral) {
        if (valueB)
          score += colA(indexB);
      } else
        score += colA(indexB) * valueA * valueB;
      ++indexB;
    }

    ++indexA;
  }

  assert(not std::isnan(score));
  return score;
}

template <typename Clusters>
void GraphCut::setInitialClusters(Clusters &&clusters) {
  using clusters_type = std::decay_t<Clusters>;
  static_assert(std::is_same_v<clusters_type, HardClusters> or
                std::is_same_v<clusters_type, WeightedClusters>);

  initialClusters.emplace<clusters_type>(std::forward<Clusters>(clusters));
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

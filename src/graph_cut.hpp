#pragma once

#include "hard_clusters.hpp"
#include "results/transcript.hpp"
#include "weighted_clusters.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"

#include <armadillo>
#include <concepts>
#include <cstdint>
#include <utility>
#include <variant>
#include <vector>

class RnaSecondaryStructure;

class GraphCut {
  struct FuzzyCut {};
  struct HardCut {};

public:
  enum class Graph {
    symmetricLaplacian,
    randomWalkLaplacian,
    laplacian,
    probabilityTransition,
    modularity,
    adjacency
  };

  constexpr static FuzzyCut fuzzy{};
  constexpr static HardCut hard{};

  GraphCut() = default;
  GraphCut(std::vector<arma::mat> const &adjacencies,
           Graph type = Graph::symmetricLaplacian);

  WeightedClusters run(std::uint8_t nClusters,
                       std::uint16_t kmeans_iterations) const;
  HardClusters run(std::uint8_t nClusters, HardCut) const;

  double calculateClustersScore(
      const std::vector<std::vector<bool>> &rawClusters) const;
  double calculateClustersScore(const HardClusters &clusters) const;
  template <typename Clusters> void setInitialClusters(Clusters &&clusters);

private:
  inline arma::mat createGraph(const arma::mat &adjacency) const;
  inline arma::mat createSymmetricLaplacian(const arma::mat &adjacency) const;
  template <typename Fun, typename Gen>
    requires requires(Fun fun) {
      { fun(std::declval<arma::mat const &>()) } -> std::same_as<arma::mat>;
    }
  WeightedClusters partitionGraph(std::uint8_t nClusters,
                                  std::uint16_t kmeans_iterations, Fun graphFun,
                                  Gen &&random_generator) const;
  template <typename Fun>
  HardClusters partitionGraphHard(std::uint8_t nClusters, Fun graphFun) const;
  template <typename Clusters>
  double calculateCutScore(const arma::mat &graph,
                           const Clusters &clusters) const;
  template <typename ClusterA, typename ClusterB>
  static double calculateCut(const arma::mat &graph, const ClusterA &clusterA,
                             const ClusterB &clusterB);

  arma::mat getGraphWithNoLoops(const arma::mat &matrix) const;

  Graph graphType = Graph::symmetricLaplacian;
  std::vector<arma::mat> adjacencies;
  std::variant<std::monostate, HardClusters, WeightedClusters> initialClusters;
};

arma::mat pairwise_distances(arma::mat const &a, arma::subview<double> b);

template <typename Gen>
arma::mat eigenvectors_to_weighted_clusters(arma::mat const &eigenvectors,
                                            std::uint8_t n_clusters,
                                            std::uint16_t kmeans_iterations,
                                            Gen &&rundom_generator);

template <typename Fun, typename Gen>
  requires requires(Fun fun) {
    { fun(std::declval<arma::mat const &>()) } -> std::same_as<arma::mat>;
  }
WeightedClusters weighted_clusters_from_adjacency(
    std::uint8_t nClusters, std::uint16_t kmeans_iterations, Fun graphFun,
    Gen &&random_generator, arma::mat const &adjacency);

WeightedClusters
merge_weighted_clusters(std::vector<WeightedClusters> &&all_weighted_clusters);

#include "graph_cut_impl.hpp"

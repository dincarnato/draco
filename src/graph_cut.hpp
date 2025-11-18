#pragma once

#include "results/transcript.hpp"
#include "weighted_clusters.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"

#include <armadillo>
#include <vector>

class GraphCut {
public:
  enum class Graph {
    symmetricLaplacian,
    randomWalkLaplacian,
    laplacian,
    probabilityTransition,
    modularity,
    adjacency
  };

  GraphCut() = default;
  GraphCut(std::vector<arma::mat> const &adjacencies,
           Graph type = Graph::symmetricLaplacian);

  WeightedClusters run(std::uint8_t nClusters,
                       std::uint16_t kmeans_iterations) const;

  double calculateClustersScore(
      const std::vector<std::vector<bool>> &rawClusters) const;

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

  arma::mat getGraphWithNoLoops(const arma::mat &matrix) const;

  Graph graphType = Graph::symmetricLaplacian;
  std::vector<arma::mat> adjacencies;
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

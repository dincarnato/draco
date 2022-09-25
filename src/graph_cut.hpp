#pragma once

#include "hard_clusters.hpp"
#include "weighted_clusters.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"

#include <armadillo>
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
  GraphCut(const arma::mat &adjacency, Graph type = Graph::symmetricLaplacian);

  WeightedClusters run(std::uint8_t nClusters, std::uint16_t iterations,
                       FuzzyCut = fuzzy) const;
  HardClusters run(std::uint8_t nClusters, HardCut) const;

  double calculateClustersScore(
      const std::vector<std::vector<bool>> &rawClusters) const;
  double calculateClustersScore(const HardClusters &clusters) const;
  template <typename Clusters> void setInitialClusters(Clusters &&clusters);

private:
  inline arma::mat createGraph(const arma::mat &adjacency) const;
  inline arma::mat createSymmetricLaplacian(const arma::mat &adjacency) const;
  template <typename Fun>
  std::tuple<WeightedClusters, double> partitionGraph(std::uint8_t nClusters,
                                                      Fun graphFun) const;
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
  arma::mat adjacency;
  std::variant<std::monostate, HardClusters, WeightedClusters> initialClusters;
};

#include "graph_cut_impl.hpp"

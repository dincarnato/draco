#include "graph_cut.hpp"
#include "rna_secondary_structure.hpp"

#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>

GraphCut::GraphCut(const arma::mat &adjacency, Graph type)
    : graphType(type), adjacency(adjacency) {
  for (std::size_t row = 0; row < adjacency.n_rows; ++row) {
    if (arma::all(arma::abs(adjacency.row(row)) < 1e-6)) {
      std::cerr << std::setprecision(2);
      adjacency.print(std::cerr);
      throw std::runtime_error(
          "adjacency cannot have rows with all zero values");
    }
  }

  for (std::size_t col = 0; col < adjacency.n_cols; ++col) {
    if (arma::all(arma::abs(adjacency.col(col)) < 1e-6)) {
      std::cerr << std::setprecision(2);
      adjacency.print(std::cerr);
      throw std::runtime_error(
          "adjacency cannot have cols with all zero values");
    }
  }
}

arma::mat GraphCut::createGraph(const arma::mat &adjacency) const {
  switch (graphType) {
  case Graph::symmetricLaplacian: {
    arma::mat symLaplacian = createSymmetricLaplacian(adjacency);
    assert(not symLaplacian.has_nan());
    return symLaplacian;
  }
  case Graph::adjacency:
    return adjacency;
  default:
    throw std::runtime_error("graph type not implemeted, sorry :(");
  }
}

arma::mat GraphCut::createSymmetricLaplacian(const arma::mat &adjacency) const {
  arma::mat adjacency_zero_diag = arma::log(adjacency + 1.);
  adjacency_zero_diag.diag() = arma::zeros(adjacency_zero_diag.n_cols);
  arma::vec degree = arma::sum(adjacency_zero_diag, 1);
  arma::mat laplacian = arma::diagmat(degree) - adjacency_zero_diag;
  std::transform(std::begin(degree), std::end(degree), std::begin(degree),
                 [](double degree) {
                   if (std::abs(degree) < 1e-4)
                     return 1.;
                   else
                     return degree;
                 });
  arma::mat invSqrtDiagonal =
      arma::diagmat(arma::ones(degree.size()) / arma::sqrt(degree));

  return invSqrtDiagonal * laplacian * invSqrtDiagonal;
}

arma::mat GraphCut::getGraphWithNoLoops(const arma::mat &matrix) const {
  arma::mat symLaplacian = createGraph(matrix);
  symLaplacian.diag() = arma::zeros(symLaplacian.n_rows);
  return symLaplacian;
}

auto GraphCut::run(std::uint8_t nClusters, float weightModule,
                   std::uint16_t nTries, std::uint16_t iterations,
                   FuzzyCut) const -> WeightedClusters {
  if (nClusters < 2)
    throw std::logic_error("nClusters must be at least 2");

  auto createPartitionGraph = [this, nClusters, weightModule, nTries]() {
    return partitionGraph(
        nClusters, weightModule, nTries,
        [this](const auto &matrix) { return getGraphWithNoLoops(matrix); });
  };

  auto [clusters, score] = createPartitionGraph();
  for (std::uint16_t iteration = 1; iteration < iterations; ++iteration) {
    auto [newClusters, newScore] = createPartitionGraph();
    if (not std::isnan(newScore) and newScore < score) {
      clusters = std::move(newClusters);
      score = newScore;
    }
  }

  return clusters;
}

auto GraphCut::run(std::uint8_t nClusters, HardCut) const -> HardClusters {
  if (nClusters < 2)
    throw std::logic_error("nClusters must be at least 2");

  auto clusters = partitionGraphHard(nClusters, [this](const auto &matrix) {
    return getGraphWithNoLoops(matrix);
  });

  return clusters;
}

double GraphCut::calculateClustersScore(
    const std::vector<std::vector<bool>> &rawClusters) const {
  assert(not rawClusters.empty());
  const std::size_t nElements = adjacency.n_cols;

#ifndef NDEBUG
  for (const auto &cluster : rawClusters)
    assert(cluster.size() == nElements);
#endif

  HardClusters hardClusters(nElements,
                            static_cast<std::uint8_t>(rawClusters.size()));
  {
    auto &&clusters = hardClusters.clusters();
    auto clusterIter = std::ranges::begin(clusters);
    const auto clustersEnd = std::ranges::end(clusters);
    auto rawClusterIter = std::ranges::begin(rawClusters);

    for (; clusterIter < clustersEnd; ++clusterIter, ++rawClusterIter) {
      auto &&cluster = *clusterIter;
      auto &&rawCluster = *rawClusterIter;

      auto baseIter = std::ranges::begin(cluster);
      const auto clusterEnd = std::ranges::end(cluster);
      auto rawBaseIter = std::ranges::begin(rawCluster);

      for (; baseIter < clusterEnd; ++baseIter, ++rawBaseIter) {
        if (*rawBaseIter) {
          auto wrapper = *baseIter;
          wrapper.set();
        }
      }
    }
  }

  return calculateCutScore(getGraphWithNoLoops(adjacency), hardClusters);
}

double GraphCut::calculateClustersScore(const HardClusters &clusters) const {
  return calculateCutScore(getGraphWithNoLoops(adjacency), clusters);
}

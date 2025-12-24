#include "graph_cut.hpp"
#include "logger.hpp"
#include "rna_secondary_structure.hpp"

#include <algorithm>
#include <cassert>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_reduce.h>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <tuple>

GraphCut::GraphCut(std::vector<arma::mat> const &adjacencies, Graph type)
    : graphType(type), adjacencies(adjacencies) {
  for (const auto &adjacency : adjacencies) {
    for (std::size_t row = 0; row < adjacency.n_rows; ++row) {
      if (arma::all(arma::abs(adjacency.row(row)) < 1e-6)) {
        std::stringstream repr;
        repr << std::setprecision(2);
        adjacency.print(repr);

        logger::error("Adjacency cannot have rows with all zero values, "
                      "terminating. Matrix:\n{}",
                      repr.str());
        std::terminate();
      }
    }

    for (std::size_t col = 0; col < adjacency.n_cols; ++col) {
      if (arma::all(arma::abs(adjacency.col(col)) < 1e-6)) {
        std::stringstream repr;
        repr << std::setprecision(2);
        adjacency.print(repr);

        logger::error("Adjacency cannot have cols with all zero values, "
                      "terminating. Matrix:\n{}",
                      repr.str());
        std::terminate();
      }
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
                   std::uint16_t window_begin_index,
                   std::uint16_t window_end_index,
                   results::Transcript const &transcript, FuzzyCut) const
    -> WeightedClusters {
  if (nClusters < 2)
    throw std::logic_error("nClusters must be at least 2");

  auto createPartitionGraph = [this, nClusters, weightModule, nTries]() {
    return partitionGraph(
        nClusters, weightModule, nTries,
        [this](const auto &matrix) { return getGraphWithNoLoops(matrix); });
  };

  auto result = tbb::parallel_reduce(
      tbb::blocked_range<std::uint16_t>(
          0, std::max(iterations, static_cast<std::uint16_t>(1))),
      std::optional<std::tuple<WeightedClusters, double>>{},
      [&](const tbb::blocked_range<std::uint16_t> &range,
          std::optional<std::tuple<WeightedClusters, double>> &&best_result) {
        auto iter = std::begin(range);
        auto const end = std::end(range);
        if (not best_result.has_value() and iter != end) {
          logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                        "initial run on iteration {} ",
                        transcript.name, window_begin_index + 1,
                        window_end_index, iter + 1);
          best_result.emplace(createPartitionGraph());
          logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                        "initial best score on iteration {} is {}",
                        transcript.name, window_begin_index + 1,
                        window_end_index, iter + 1, std::get<1>(*best_result));
          ++iter;
        }

        auto [clusters, score] = *best_result;
        for (; iter != end; ++iter) {
          logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                        "running iteration {}",
                        transcript.name, window_begin_index + 1,
                        window_end_index, iter + 1);
          auto [newClusters, newScore] = createPartitionGraph();
          if (not std::isnan(newScore) and newScore < score) {
            logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                          "score on iteration {} ({}) "
                          "is better than the best score ({})",
                          transcript.name, window_begin_index + 1,
                          window_end_index, iter + 1, newScore, score);
            clusters = std::move(newClusters);
            score = newScore;
          } else {
            logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                          "score on iteration {} ({}) "
                          "is worse than the best score ({})",
                          transcript.name, window_begin_index + 1,
                          window_end_index, iter + 1, newScore, score);
          }
        }

        return best_result;
      },
      [window_begin_index, window_end_index, &transcript](
          std::optional<std::tuple<WeightedClusters, double>> &&result1,
          std::optional<std::tuple<WeightedClusters, double>> &&result2) {
        if (result1.has_value() && result2.has_value()) {
          auto const score1 = std::get<1>(*result1);
          auto const score2 = std::get<1>(*result2);
          if (std::isnan(score2) or
              (not std::isnan(score1) and score1 < score2)) {
            logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                          "using score {} that is better than {}",
                          transcript.name, window_begin_index + 1,
                          window_end_index, score1, score2);
            return result1;
          } else {
            logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                          "using score {} that is better than {}",
                          transcript.name, window_begin_index + 1,
                          window_end_index, score2, score1);
            return result2;
          }
        } else if (result2.has_value()) {
          logger::trace(
              "Graph partitioning, transcript {}, window {}-{}: using score {}",
              transcript.name, window_begin_index + 1, window_end_index,
              std::get<1>(*result2));
          return result2;
        } else {
          if (result1.has_value()) {
            logger::trace("Graph partitioning, transcript {}, window {}-{}: "
                          "using score {}",
                          transcript.name, window_begin_index + 1,
                          window_end_index, std::get<1>(*result1));
          }
          return result1;
        }
      });

  assert(result.has_value());
  logger::debug(
      "Graph partitioning, transcript {}, window {}-{}: best score is {}",
      transcript.name, window_begin_index + 1, window_end_index,
      std::get<1>(*result));
  return std::get<0>(std::move(*result));
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
  auto const &firstAdjacency = adjacencies[0];
  assert(std::ranges::all_of(adjacencies | std::views::drop(1),
                             [&](auto const &adjacency) {
                               return adjacency.n_cols == firstAdjacency.n_cols;
                             }));

  const std::size_t nElements = firstAdjacency.n_cols;

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

  return std::ranges::fold_left(
      adjacencies | std::views::transform([&](auto const &adjacency) {
        return calculateCutScore(getGraphWithNoLoops(adjacency), hardClusters);
      }),
      0., std::plus{});
}

double GraphCut::calculateClustersScore(const HardClusters &clusters) const {
  return std::ranges::fold_left(
      adjacencies | std::views::transform([&](auto const &adjacency) {
        return calculateCutScore(getGraphWithNoLoops(adjacency), clusters);
      }),
      0., std::plus{});
}

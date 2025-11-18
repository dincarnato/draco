#include "graph_cut.hpp"
#include "weighted_clusters.hpp"

#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cmath>
#include <exception>
#include <functional>
#include <iomanip>
#include <iostream>
#include <oneapi/tbb/blocked_range.h>
#include <oneapi/tbb/parallel_reduce.h>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <vector>

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

auto GraphCut::run(std::uint8_t nClusters,
                   std::uint16_t kmeans_iterations) const -> WeightedClusters {
  if (nClusters < 2)
    throw std::logic_error("nClusters must be at least 2");

  return partitionGraph(
      nClusters, kmeans_iterations,
      [this](const auto &matrix) { return getGraphWithNoLoops(matrix); },
      std::mt19937(std::random_device{}()));
}

arma::mat pairwise_distances(arma::mat const &a, arma::subview<double> b) {
  if (a.n_cols != b.n_cols) {
    throw new std::runtime_error("the matrices given to pairwise_distances "
                                 "must have the same number of columns");
  }

  arma::mat out(a.n_rows, b.n_rows);
  // Armadillo matrices are column-first, it's better to iterate over rows first
  // in order to improve caching
  for (arma::uword b_index = 0; b_index < b.n_rows; ++b_index) {
    auto b_row = b.row(b_index);
    for (arma::uword a_index = 0; a_index < a.n_rows; ++a_index) {
      out[a_index, b_index] =
          std::sqrt(arma::sum(arma::pow(a.row(a_index) - b_row, 2)));
    }
  }
  return out;
}

WeightedClusters
merge_weighted_clusters(std::vector<WeightedClusters> &&all_weighted_clusters) {
  auto n_replicates = static_cast<double>(std::size(all_weighted_clusters));

  WeightedClusters weights(std::move(all_weighted_clusters[0]));

  for (auto &&[base_index, base_weights] :
       std::views::zip(std::views::iota(0uz), weights)) {
    for (auto &&[cluster_index, weight] :
         std::views::zip(std::views::iota(0uz), base_weights)) {
      weight = static_cast<float>(
          std::ranges::fold_left(
              all_weighted_clusters | std::views::drop(1) |
                  std::views::transform([&](auto const &cluster_weights) {
                    return static_cast<double>(
                        cluster_weights[base_index][cluster_index]);
                  }),
              static_cast<double>(weight), std::plus{}) /
          n_replicates);
    }
  }

  return weights;
}

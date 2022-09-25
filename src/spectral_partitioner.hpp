#pragma once

#include <armadillo>
#include <array>
#include <vector>

struct SpectralPartitioner {
  static std::array<std::vector<std::size_t>, 2>
  bipartite(const arma::mat &adjacency);

  static arma::mat getNormalizedLaplacian(arma::mat adjacency);

  static std::pair<arma::mat, arma::vec>
  getDecomposition(const arma::mat &normalizedLaplacian);
};

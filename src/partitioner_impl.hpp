#pragma once

#include "arma_iterator.hpp"
#include "partitioner.hpp"
#include "utils.hpp"

#include <cmath>
#include <random>
#include <tuple>

#include <iomanip>

template <typename Read>
std::tuple<double, unsigned>
Partitioner::findNearestCentroid(Read&& read,
                                 const arma::mat& centroidsMatrix) {
  assert(centroidsMatrix.n_rows > 0);
  assert(centroidsMatrix.n_cols > 0);
  arma::mat coordinateWiseDifference =
      arma::repmat(read, centroidsMatrix.n_rows, 1) - centroidsMatrix;
  arma::vec distances = arma::sum(arma::pow(coordinateWiseDifference, 2), 1);

  assert(distances.size() == centroidsMatrix.n_rows);
  unsigned minDistanceIndex = arma::index_min(distances);
  return std::make_tuple(distances[minDistanceIndex], minDistanceIndex);
}

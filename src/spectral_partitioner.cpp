#include "spectral_partitioner.hpp"

#include <cassert>
#include <cmath>
#include <numeric>

arma::mat SpectralPartitioner::getNormalizedLaplacian(arma::mat adjacency) {
  assert(not adjacency.has_nan());
  assert(adjacency.n_rows == adjacency.n_cols);
  adjacency.diag() = arma::zeros(adjacency.n_rows);
  arma::vec degreeVec = arma::sum(adjacency, 1);
  arma::mat laplacian = (arma::diagmat(degreeVec) - adjacency);

  std::transform(std::begin(degreeVec), std::end(degreeVec),
                 std::begin(degreeVec), [](double degree) {
                   if (std::abs(degree) < 1e-4)
                     return 1.;
                   else
                     return degree;
                 });
  arma::mat inverseSqrtDegreeMat =
      arma::diagmat(arma::ones(adjacency.n_cols) / arma::sqrt(degreeVec));
  arma::mat normalizedLaplacian =
      inverseSqrtDegreeMat * laplacian * inverseSqrtDegreeMat;
  assert(not normalizedLaplacian.has_nan());
  return normalizedLaplacian;
}

std::pair<arma::mat, arma::vec>
SpectralPartitioner::getDecomposition(const arma::mat &normalizedLaplacian) {
  std::pair<arma::mat, arma::vec> eigenDecomposition;
  assert(not normalizedLaplacian.has_nan());
  arma::eig_sym(eigenDecomposition.second, eigenDecomposition.first,
                normalizedLaplacian);
  // We accept very small negative values
  if (arma::any(eigenDecomposition.second < -0.0001))
    throw std::logic_error("at least one eigenvalue is negative, the input "
                           "matrix is not a valid adjacency matrix");
  return eigenDecomposition;
}

std::array<std::vector<std::size_t>, 2>
SpectralPartitioner::bipartite(const arma::mat &adjacency) {
  constexpr unsigned maxIterations = 1000;

  arma::vec secondEigenVec;
  {
    auto eigenvecs = getDecomposition(getNormalizedLaplacian(adjacency)).first;
    assert(eigenvecs.n_cols >= 2);
    secondEigenVec = eigenvecs.col(1);
  }
  assert(secondEigenVec.size() > 1);
  assert(not secondEigenVec.has_nan());

  std::array<double, 2> centroids;
  std::array<std::vector<std::size_t>, 2> assignments;

  {
    std::vector<std::size_t> initialIndices(adjacency.n_rows);
    std::iota(std::begin(initialIndices), std::end(initialIndices),
              static_cast<std::size_t>(0));
    std::sort(std::begin(initialIndices), std::end(initialIndices),
              [&](const std::size_t indexA, const std::size_t indexB) {
                return secondEigenVec[indexA] < secondEigenVec[indexB];
              });
    // auto&& [minElement, maxElement] =
    // std::minmax_element(std::begin(secondEigenVec),
    // std::end(secondEigenVec)); if(minElement == maxElement)
    //    throw std::runtime_error("second eigenvector is constant");

    std::swap(initialIndices[0], initialIndices[secondEigenVec.size() / 4]);
    std::swap(initialIndices[1], initialIndices[secondEigenVec.size() * 4 / 3]);

    auto indexIter = std::begin(initialIndices);
    assignments[0].push_back(*indexIter);
    centroids[0] = secondEigenVec[*indexIter++];
    assignments[1].push_back(*indexIter);
    centroids[1] = secondEigenVec[*indexIter++];

    assert(std::none_of(
        std::begin(centroids), std::end(centroids),
        [](const double centroid) { return std::isnan(centroid); }));

    std::array<double, 2> distances;
    for (auto endIter = std::end(initialIndices); indexIter != endIter;
         ++indexIter) {
      const double value = secondEigenVec[*indexIter];
      for (unsigned clusterIndex = 0; clusterIndex < 2; ++clusterIndex)
        distances[clusterIndex] = std::pow(value - centroids[clusterIndex], 2);

      const auto bestCentroidIter =
          std::min_element(std::begin(distances), std::end(distances));
      const auto centroidIndex = static_cast<std::size_t>(
          std::distance(std::begin(distances), bestCentroidIter));
      auto &assignment = assignments[centroidIndex];
      std::size_t assignmentSize = assignment.size();
      *bestCentroidIter =
          (*bestCentroidIter * static_cast<double>(assignmentSize) + value) /
          static_cast<double>(assignmentSize + 1);
      assignment.push_back(*indexIter);
    }

    assert(std::all_of(
        std::begin(assignments), std::end(assignments),
        [](const auto &assignment) { return not assignment.empty(); }));
  }

  for (unsigned iterationIndex = 0; iterationIndex < maxIterations;
       ++iterationIndex) {
    std::array<std::vector<std::size_t>, 2> newAssignments;
    for (std::size_t elementIndex = 0; elementIndex < adjacency.n_rows;
         ++elementIndex) {
      double value = secondEigenVec[elementIndex];
      std::array<double, 2> distances;
      for (unsigned clusterIndex = 0; clusterIndex < 2; ++clusterIndex)
        distances[clusterIndex] = std::pow(value - centroids[clusterIndex], 2);

      newAssignments[static_cast<std::size_t>(
                         std::distance(std::begin(distances),
                                       std::min_element(std::begin(distances),
                                                        std::end(distances))))]
          .push_back(elementIndex);
    }
    assert(std::all_of(
        std::begin(newAssignments), std::end(newAssignments),
        [](const auto &assignment) { return not assignment.empty(); }));

    if (newAssignments == assignments)
      return assignments;

    assignments = std::move(newAssignments);
    for (unsigned clusterIndex = 0; clusterIndex < 2; ++clusterIndex) {
      const auto &assignment = assignments[clusterIndex];
      centroids[clusterIndex] =
          std::accumulate(std::begin(assignment), std::end(assignment), 0.,
                          [&](double accumulator, std::size_t index) {
                            return accumulator + secondEigenVec[index];
                          }) /
          static_cast<double>(assignment.size());
    }

    assert(std::none_of(
        std::begin(centroids), std::end(centroids),
        [](const double centroid) { return std::isnan(centroid); }));
  }

  assert(std::all_of(
      std::begin(assignments), std::end(assignments),
      [](const auto &assignment) { return not assignment.empty(); }));
  return assignments;
}

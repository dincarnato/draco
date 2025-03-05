#pragma once

#include "ringmap_data.hpp"
#include "utils.hpp"

#include <armadillo>
#include <functional>
#include <memory>
#include <tuple>
#include <vector>

class Partitioner {
public:
  Partitioner(const RingmapData &data);

  template <typename ArmaVec>
  std::vector<std::vector<unsigned>>
  buildAssignmentIndices(unsigned clusters,
                         ArmaVec &&eigenVectorsIndices) const;

  template <typename ArmaVec>
  std::vector<std::vector<unsigned>>
  buildAssignmentIndicesAndDumpScores(unsigned clusters,
                                      ArmaVec &&eigenVectorsIndices,
                                      const std::string &scoresFilename) const;

  void setReRunWithChanges(bool value = true);
  void setMinClusterAbundancy(double value);
  void setUseWeightedChanges(bool value = true);
  void setUseEigenGapWeights(bool value = true);
  void setMinClusters(unsigned value);

private:
  template <typename ArmaVec>
  std::vector<std::vector<unsigned>>
  buildAssignmentIndices(unsigned clusters, ArmaVec &&eigenVectorsIndices,
                         const std::string *scoresFilename) const;

  template <typename Read>
  static std::tuple<double, unsigned>
  findNearestCentroid(Read &&read, const arma::mat &centroidsMatrix);
  const RingmapData *ringmapData;

  double minClusterAbundancy = 0.1;
  unsigned minClusters = 0;
  bool reRunWithChanges = true;
  bool useWeightedChanges = false;
  bool useEigenGapWeights = true;
};

#include "partitioner_impl.hpp"

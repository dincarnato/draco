#include "partitioner.hpp"

Partitioner::Partitioner(const RingmapData &data) : ringmapData(&data) {}

void Partitioner::setReRunWithChanges(bool value) { reRunWithChanges = value; }

void Partitioner::setMinClusterAbundancy(double value) {
  minClusterAbundancy = value;
}

void Partitioner::setUseWeightedChanges(bool value) {
  useWeightedChanges = value;
}

void Partitioner::setUseEigenGapWeights(bool value) {
  useEigenGapWeights = value;
}

void Partitioner::setMinClusters(unsigned value) { minClusters = value; }

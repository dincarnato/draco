#include "partitioner.hpp"

Partitioner::Partitioner(const RingmapData& data) : ringmapData(&data) {}

std::tuple<std::unique_ptr<RingmapData>, const RingmapData*>
Partitioner::getFilteredData() const {
  const RingmapData* filteredData = ringmapData;
  std::unique_ptr<RingmapData> ownedFilteredData;
  if (filteredData->getModificationsFilter() < 2 or
      not filteredData->isFiltered()) {
    ownedFilteredData = std::make_unique<RingmapData>(*ringmapData);
    ownedFilteredData->filterReads();
    ownedFilteredData->filter();
    filteredData = ownedFilteredData.get();
  }

  return std::make_tuple(std::move(ownedFilteredData), filteredData);
}

void
Partitioner::setReRunWithChanges(bool value) {
  reRunWithChanges = value;
}

void
Partitioner::setMinClusterAbundancy(double value) {
  minClusterAbundancy = value;
}

void
Partitioner::setUseWeightedChanges(bool value) {
  useWeightedChanges = value;
}

void
Partitioner::setUseEigenGapWeights(bool value) {
  useEigenGapWeights = value;
}

void
Partitioner::setMinClusters(unsigned value) {
  minClusters = value;
}

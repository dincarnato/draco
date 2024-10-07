#pragma once

#include "weighted_clusters.hpp"

#include <vector>

struct Window {
  unsigned short start_base;
  WeightedClusters weights;
  std::vector<unsigned> coverages;
};

#pragma once

#include <cstddef>

struct WindowClustersWithConfidence {
  unsigned n_clusters;
  float confidence;
  std::size_t start_base;
  std::size_t end_base;
};

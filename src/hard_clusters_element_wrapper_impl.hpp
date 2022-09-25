#pragma once

#include "hard_clusters_element_wrapper.hpp"

#include <cassert>
#include <limits>

template <typename Cluster, bool complemented>
HardClustersElementWrapper<Cluster, complemented>::HardClustersElementWrapper(
    Cluster &cluster, std::size_t index) noexcept
    : _element(&cluster._clusters[index])
#ifndef NDEBUG
      ,
      _cluster(&cluster)
#endif
{
}

template <typename Cluster, bool complemented>
auto HardClustersElementWrapper<Cluster, complemented>::get() const noexcept
    -> index_type {
  static_assert(not complemented);
  return *_element;
}

template <typename Cluster, bool complemented>
void HardClustersElementWrapper<Cluster, complemented>::set(
    index_type clusterIndex) const noexcept {
  static_assert(not complemented);
  assert(clusterIndex < _cluster->nClusters);
  *_element = clusterIndex;
}

template <typename Cluster, bool complemented>
void HardClustersElementWrapper<Cluster, complemented>::clear() const noexcept {
  static_assert(not complemented);
  *_element = std::numeric_limits<index_type>::max();
}

template <typename Cluster, bool complemented>
bool HardClustersElementWrapper<Cluster, complemented>::isCluster(
    index_type index) const noexcept {
  if constexpr (not complemented)
    return *_element == index;
  else
    return *_element != index;
}

#pragma once

#include "hard_cluster_wrapper.hpp"
#include "hard_clusters_iterator.hpp"
#include "hard_clusters_wrapper.hpp"

#include <range/v3/algorithm.hpp>

template <typename Cluster, bool complemented>
HardClustersWrapper<Cluster, complemented>::HardClustersWrapper(
    Cluster& clusters) noexcept
    : _clusters(&clusters) {}

template <typename Cluster, bool complemented>
auto HardClustersWrapper<Cluster, complemented>::
operator[](index_type index) const noexcept -> cluster_wrapper {
  return {*_clusters, index};
}

template <typename Cluster, bool complemented>
auto
HardClustersWrapper<Cluster, complemented>::begin() const noexcept -> iterator {
  return iterator{*_clusters, 0};
}

template <typename Cluster, bool complemented>
auto
HardClustersWrapper<Cluster, complemented>::end() const noexcept -> iterator {
  return iterator{*_clusters, _clusters->nClusters};
}

template <typename Cluster, bool complemented>
std::size_t
HardClustersWrapper<Cluster, complemented>::size() const noexcept {
  assert(_clusters);
  return _clusters->getClustersSize();
}

template <typename Cluster, bool complemented>
std::size_t
HardClustersWrapper<Cluster, complemented>::elements_size() const noexcept {
  assert(_clusters);
  return _clusters->getElementsSize();
}

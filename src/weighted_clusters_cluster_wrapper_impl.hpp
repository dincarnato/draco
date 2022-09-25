#pragma once

#include "weighted_clusters_cluster.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"

#include <range/v3/core.hpp>

template <typename T>
WeightedClustersClusterWrapper<T>::WeightedClustersClusterWrapper(
    weighted_clusters_type &weightedClusters,
    std::ptrdiff_t clusterIndex) noexcept
    : weightedClusters(&weightedClusters), clusterIndex(clusterIndex) {
  assert(this->weightedClusters);
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<std::ptrdiff_t>(weightedClusters._clusters));
}

template <typename T>
auto WeightedClustersClusterWrapper<T>::begin() const noexcept -> iterator {
  assert(weightedClusters);
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<std::ptrdiff_t>(weightedClusters->_clusters));

  return iterator{*weightedClusters, clusterIndex, 0};
}

template <typename T>
auto WeightedClustersClusterWrapper<T>::end() const noexcept -> iterator {
  assert(weightedClusters);
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<std::ptrdiff_t>(weightedClusters->_clusters));

  return iterator{*weightedClusters, clusterIndex, iterator::endTag};
}

template <typename T>
auto WeightedClustersClusterWrapper<T>::operator[](
    std::size_t index) const noexcept -> reference {
  assert(weightedClusters);
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<std::ptrdiff_t>(weightedClusters->_clusters));

  return weightedClusters->weights[static_cast<std::size_t>(
      static_cast<std::ptrdiff_t>(index) *
          static_cast<std::ptrdiff_t>(weightedClusters->_clusters) +
      clusterIndex)];
}

template <typename T>
std::size_t WeightedClustersClusterWrapper<T>::index() const noexcept {
  assert(weightedClusters);
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<std::ptrdiff_t>(weightedClusters->_clusters));

  return static_cast<std::size_t>(clusterIndex);
}

template <typename T>
std::size_t WeightedClustersClusterWrapper<T>::size() const noexcept {
  assert(weightedClusters);
  return weightedClusters->elements;
}

template <typename T>
auto WeightedClustersClusterWrapper<T>::operator=(concrete_type const &cluster)
    const noexcept(std::is_nothrow_copy_assignable_v<std::decay_t<T>>)
        -> WeightedClustersClusterWrapper & {
  ranges::copy(cluster, ranges::begin(*this));
  return *this;
}

template <typename T>
auto WeightedClustersClusterWrapper<T>::operator=(concrete_type &&cluster) const
    noexcept(std::is_nothrow_move_assignable_v<std::decay_t<T>>)
        -> WeightedClustersClusterWrapper & {
  ranges::move(cluster, ranges::begin(*this));
  return *this;
}

template <typename T>
WeightedClustersClusterWrapper<T>::operator concrete_type() const
    noexcept(false) {
  return concrete_type(*this);
}

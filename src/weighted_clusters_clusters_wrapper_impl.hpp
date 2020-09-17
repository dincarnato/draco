#pragma once

#include "weighted_clusters_clusters_wrapper.hpp"

template <typename T>
WeightedClustersClustersWrapper<T>::WeightedClustersClustersWrapper(
    weighted_clusters_type& weightedClusters) noexcept
    : weightedClusters(&weightedClusters) {
  assert(this->weightedClusters);
}

template <typename T>
auto
WeightedClustersClustersWrapper<T>::
operator=(concrete_type const& concrete_clusters) const
    noexcept(std::is_nothrow_copy_assignable_v<T>)
        -> WeightedClustersClustersWrapper& {
  assert(concrete_clusters.size() == weightedClusters->_clusters);
  assert(concrete_clusters.elements_size() == weightedClusters->_elements);

  ranges::copy(concrete_clusters, begin());
  return *this;
}

template <typename T>
auto
WeightedClustersClustersWrapper<T>::
operator=(concrete_type&& concrete_clusters) const
    noexcept(std::is_nothrow_move_assignable_v<T>)
        -> WeightedClustersClustersWrapper& {
  assert(concrete_clusters.size() == weightedClusters->_clusters);
  assert(concrete_clusters.elements_size() == weightedClusters->_elements);

  ranges::move(concrete_clusters, begin());
  return *this;
}

template <typename T>
auto
WeightedClustersClustersWrapper<T>::begin() const noexcept -> iterator {
  assert(weightedClusters);
  return {*weightedClusters, 0};
}

template <typename T>
auto
WeightedClustersClustersWrapper<T>::end() const noexcept -> iterator {
  assert(weightedClusters);
  return iterator{*weightedClusters};
}

template <typename T>
std::size_t
WeightedClustersClustersWrapper<T>::size() const noexcept {
  assert(weightedClusters);
  return weightedClusters->getClustersSize();
}

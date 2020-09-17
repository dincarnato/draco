#pragma once

#include "weighted_clusters_cluster_iterator.hpp"

template <typename T>
WeightedClustersClusterIterator<T>::WeightedClustersClusterIterator(
    weighted_clusters_type& weightedClusters) noexcept
    : wrapper(weightedClusters,
              static_cast<difference_type>(weightedClusters._clusters)) {}

template <typename T>
WeightedClustersClusterIterator<T>::WeightedClustersClusterIterator(
    weighted_clusters_type& weightedClusters,
    difference_type clusterIndex) noexcept
    : wrapper(weightedClusters, clusterIndex) {
  assert(clusterIndex >= 0);
  assert(clusterIndex <=
         static_cast<difference_type>(weightedClusters._clusters));
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator++() noexcept
    -> WeightedClustersClusterIterator& {
  ++wrapper.clusterIndex;
  return *this;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator++(int) noexcept
    -> WeightedClustersClusterIterator {
  auto other = *this;
  ++wrapper.clusterIndex;

  return other;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator--() noexcept
    -> WeightedClustersClusterIterator& {
  --wrapper.clusterIndex;
  return *this;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator--(int) noexcept
    -> WeightedClustersClusterIterator {
  auto other = *this;
  --wrapper.clusterIndex;

  return other;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator+=(difference_type offset) noexcept
    -> WeightedClustersClusterIterator& {
  wrapper.clusterIndex += offset;
  return *this;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator-=(difference_type offset) noexcept
    -> WeightedClustersClusterIterator& {
  wrapper.clusterIndex -= offset;
  return *this;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator+(difference_type offset) const
    noexcept -> WeightedClustersClusterIterator {
  auto newIter = *this;
  newIter.wrapper.clusterIndex += offset;
  return newIter;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::operator-(difference_type offset) const
    noexcept -> WeightedClustersClusterIterator {
  auto newIter = *this;
  newIter.wrapper.clusterIndex -= offset;
  return newIter;
}

template <typename U>
WeightedClustersClusterIterator<U>
operator+(typename WeightedClustersClusterIterator<U>::difference_type offset,
          const WeightedClustersClusterIterator<U>& iter) noexcept {
  auto newIter = iter;
  newIter.wrapper.clusterIndex += offset;
  return newIter;
}

template <typename T>
auto WeightedClustersClusterIterator<T>::operator*() const noexcept
    -> reference {
  assert(wrapper.clusterIndex <
         static_cast<difference_type>(wrapper.weightedClusters->_clusters));
  return wrapper;
}

template <typename T>
auto WeightedClustersClusterIterator<T>::
operator[](difference_type offset) const noexcept -> reference {
  assert(wrapper.clusterIndex + offset <
         static_cast<difference_type>(wrapper.weightedClusters->_clusters));
  auto newWrapper = wrapper;
  newWrapper.clusterIndex += offset;
  return newWrapper;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator==(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex == other.wrapper.clusterIndex;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator!=(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex != other.wrapper.clusterIndex;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator<(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex < other.wrapper.clusterIndex;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator<=(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex <= other.wrapper.clusterIndex;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator>=(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex >= other.wrapper.clusterIndex;
}

template <typename T>
bool
WeightedClustersClusterIterator<T>::
operator>(const WeightedClustersClusterIterator& other) const noexcept {
  return wrapper.clusterIndex > other.wrapper.clusterIndex;
}

template <typename T>
auto
WeightedClustersClusterIterator<T>::
operator-(const WeightedClustersClusterIterator& other) const noexcept
    -> difference_type {
  return static_cast<difference_type>(wrapper.clusterIndex) -
         other.wrapper.clusterIndex;
}

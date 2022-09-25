#pragma once

#include "weighted_clusters_cluster.hpp"
#include "weighted_clusters_cluster_wrapper.hpp"

#include <range/v3/algorithm.hpp>

template <typename T>
template <typename U>
WeightedClustersCluster<T>::WeightedClustersCluster(
    WeightedClustersClusterWrapper<U> const &wrapper) noexcept(false) {
  init_from_wrapper(wrapper);
}

template <typename T>
template <typename U>
WeightedClustersCluster<T>::WeightedClustersCluster(
    WeightedClustersClusterWrapper<U> &&wrapper) noexcept(false) {
  init_from_wrapper(std::move(wrapper));
}

template <typename T>
template <typename Wrapper>
void WeightedClustersCluster<T>::init_from_wrapper(Wrapper &&wrapper) noexcept(
    false) {
  base_type::resize(wrapper.size());
  ranges::copy(wrapper, std::begin(*this));
}

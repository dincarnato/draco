#pragma once

#include <vector>

template <typename T>
class WeightedClustersClusterWrapper;

template <typename T>
struct WeightedClustersCluster : std::vector<T> {
  using base_type = std::vector<T>;

  WeightedClustersCluster() = default;
  template <typename U>
  WeightedClustersCluster(
      WeightedClustersClusterWrapper<U> const& wrapper) noexcept(false);

  template <typename U>
  WeightedClustersCluster(WeightedClustersClusterWrapper<U>&& wrapper) noexcept(
      false);

private:
  template <typename Wrapper>
  void init_from_wrapper(Wrapper&& wrapper) noexcept(false);
};

#include "weighted_clusters_cluster_impl.hpp"

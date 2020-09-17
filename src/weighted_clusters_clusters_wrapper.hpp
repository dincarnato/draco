#pragma once

#include "weighted_clusters_cluster_iterator.hpp"

#include <type_traits>

class WeightedClusters;
template <typename>
struct WeightedClustersClusters;

template <typename T>
class WeightedClustersClustersWrapper {
  using weighted_clusters_type =
      std::conditional_t<std::is_const<T>::value, const WeightedClusters,
                         WeightedClusters>;

public:
  using concrete_type = WeightedClustersClusters<std::decay_t<T>>;
  using iterator = WeightedClustersClusterIterator<T>;

  explicit WeightedClustersClustersWrapper(
      weighted_clusters_type& weightedClusters) noexcept;

  WeightedClustersClustersWrapper&
  operator=(concrete_type const& concrete_clusters) const
      noexcept(std::is_nothrow_copy_assignable_v<T>);
  WeightedClustersClustersWrapper&
  operator=(concrete_type&& concrete_clusters) const
      noexcept(std::is_nothrow_move_assignable_v<T>);

  operator concrete_type() const noexcept;

  iterator begin() const noexcept;
  iterator end() const noexcept;
  std::size_t size() const noexcept;

private:
  weighted_clusters_type* weightedClusters;
};

#include "weighted_clusters_clusters_wrapper_impl.hpp"
